#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include <valarray>
#include <sstream>
#include <filesystem>
#include <string>

using namespace std::string_literals;

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "miMaS/field.h"
#include "miMaS/complex_field.h"
#include "miMaS/weno.h"
#include "miMaS/fft.h"
#include "miMaS/array_view.h"
#include "miMaS/poisson.h"
#include "miMaS/rk.h"
#include "miMaS/config.h"
#include "miMaS/signal_handler.h"
#include "miMaS/iteration.h"

namespace math = boost::math::constants;
const std::complex<double> & I = std::complex<double>(0.,1.);

#define SQ(X) ((X)*(X))
#define Xi(i) (i*fh.step.dx+fh.range.x_min)
#define Vk(k) (k*fh.step.dv+fh.range.v_min)

auto
maxwellian ( double rho , double u , double T ) {
  return [=](double x,double v){ return rho/(std::sqrt(2.*math::pi<double>()*T))*std::exp( -0.5*SQ(v-u)/T ); };
}

double
error_H ( std::size_t Nx , std::size_t Nv , double dt )
{
  const double Tf = 5.;

  field<double,1> fh(boost::extents[Nv][Nx]);
  complex_field<double,1> hfh(boost::extents[Nv][Nx]);

  const double Kx = 0.5;
  fh.range.v_min = -8.; fh.range.v_max = 8.;
  fh.range.x_min =  0.; fh.range.x_max = 2./Kx*math::pi<double>();
  fh.compute_steps();

  ublas::vector<double> v(Nv,0.);
  std::generate( v.begin() , v.end() , [&,k=0]() mutable {return (k++)*fh.step.dv+fh.range.v_min;} );

  ublas::vector<double> kx(Nx); // beware, Nx need to be odd
  {
    double l = fh.range.len_x();
    for ( auto i=0 ; i<Nx/2 ; ++i ) { kx[i]    = 2.*math::pi<double>()*i/l; }
    for ( int i=-Nx/2 ; i<0 ; ++i ) { kx[Nx+i] = 2.*math::pi<double>()*i/l; }
  }

  const double alpha = 0.2 , ui = 2.;
  auto tb_M1 = maxwellian( 0.5*alpha , ui , 1. ) , tb_M2 = maxwellian( 0.5*alpha , -ui , 1. );
  for (field<double,2>::size_type k=0 ; k<fh.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<fh.size(1) ; ++i ) {
      fh[k][i] = ( tb_M1(Xi(i),Vk(k)) + tb_M2(Xi(i),Vk(k)) )*(1. + 0.01*std::cos(Kx*Xi(i)));
    }
    fft::fft(fh[k].begin(),fh[k].end(),hfh[k].begin());
  }

  std::size_t iter = 0;
  double current_time = 0.;

  ublas::vector<double> uc(Nx,0.);
  ublas::vector<double> E (Nx,0.);

  std::vector<double> H; H.reserve(int(std::ceil(Tf/dt))+1);

  const double rho_c = 1.-alpha;
  const double sqrt_rho_c = std::sqrt(rho_c);

  // init E with Poisson solver, init also ee, Emax, H and times
  {
    poisson<double> poisson_solver(Nx,fh.range.len_x());
    ublas::vector<double> rho(Nx,0.); rho = fh.density(); // compute density from init hot data
    for ( auto i=0 ; i<Nx ; ++i ) { rho[i] += (1.-alpha); } // add (1-alpha) for cold particules
    E = poisson_solver(rho);

    double total_energy = energy(fh,E);
    total_energy += 0.; // sum(rho_c*u_c*u_c) = 0 because u_c = 0 at time 0
    H.push_back( total_energy );
  }

  // initialize memory for all temporary variables
  ublas::vector<double> J(Nx,0.);
  fft::spectrum_ d(Nx);
  field<double,1> Edvf(tools::array_view<const std::size_t>(fh.shape(),2));
  ublas::vector<double> uc1(Nx) , uc2(Nx) , uc3(Nx) , uc4(Nx) , uc5(Nx) , uc6(Nx) , uc7(Nx),
                        E1 (Nx) , E2 (Nx) , E3 (Nx) , E4 (Nx) , E5 (Nx) , E6 (Nx) , E7 (Nx);
  complex_field<double,1> hfh1(boost::extents[Nv][Nx]) , hfh2(boost::extents[Nv][Nx]) ,
                          hfh3(boost::extents[Nv][Nx]) , hfh4(boost::extents[Nv][Nx]) ,
                          hfh5(boost::extents[Nv][Nx]) , hfh6(boost::extents[Nv][Nx]) ,
                          hfh7(boost::extents[Nv][Nx]) ;
  fh.write("init.dat");

  while (  current_time < Tf ) {

///////////////////////////////////////////////////////////////////////////////
// DP4(3) /////////////////////////////////////////////////////////////////////

    // STAGE 1
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh[k].begin(),hfh[k].end(),fh[k].begin()); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E);

      double c05 = std::cos(0.5*dt*sqrt_rho_c), s05 = std::sin(0.5*dt*sqrt_rho_c);
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc1[i] =  uc[i]*c05 + E[i]*s05/sqrt_rho_c - 0.5*dt*J[i]*s05/sqrt_rho_c;
        E1[i]  = -uc[i]*s05*sqrt_rho_c + E[i]*c05 - 0.5*dt*J[i]*c05;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh1[k][i] = hfh[k][i]*std::exp(-0.5*I*kx[i]*v[k]*dt) - 0.5*dt*d[i]*std::exp(-0.5*I*kx[i]*v[k]*dt);
        }
      }
    } // end stage 1

    // STAGE 2
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh1[k].begin(),hfh1[k].end(),fh[k].begin()); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E1);

      double c05 = std::cos(0.5*dt*sqrt_rho_c), s05 = std::sin(0.5*dt*sqrt_rho_c);
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc2[i] =  uc[i]*c05 + E[i]*s05/sqrt_rho_c;
        E2[i]  = -uc[i]*s05*sqrt_rho_c + E[i]*c05 - 0.5*dt*J[i];
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh2[k][i] = hfh[k][i]*std::exp(-0.5*I*kx[i]*v[k]*dt) - 0.5*dt*d[i];
        }
      }
    } // end stage 2

    // STAGE 3
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh2[k].begin(),hfh2[k].end(),fh[k].begin()); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E2);

      double c1  = std::cos(dt*sqrt_rho_c)     , s1  = std::sin(dt*sqrt_rho_c)    ,
             c05 = std::cos(0.5*dt*sqrt_rho_c) , s05 = std::sin(0.5*dt*sqrt_rho_c);
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc3[i] =  uc[i]*c1 + E[i]*s1/sqrt_rho_c - dt*J[i]*s05/sqrt_rho_c;
        E3[i]  = -uc[i]*s1*sqrt_rho_c + E[i]*c1 - dt*J[i]*c05;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh3[k][i] = hfh[k][i]*std::exp(-I*kx[i]*v[k]*dt) - dt*d[i]*std::exp(-0.5*I*kx[i]*v[k]*dt);
        }
      }
    } // end stage 3

    // STAGE 4
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh3[k].begin(),hfh3[k].end(),fh[k].begin()); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E3);

      double c1  = std::cos(dt*sqrt_rho_c)     , s1  = std::sin(dt*sqrt_rho_c)    ,
             c05 = std::cos(0.5*dt*sqrt_rho_c) , s05 = std::sin(0.5*dt*sqrt_rho_c);
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc4[i] = -(1./3.)*(  uc[i]*c1 + E[i]*s1/sqrt_rho_c ) + (1./3.)*(  uc1[i]*c05 + E1[i]*s05/sqrt_rho_c ) + (2./3.)*(  uc2[i]*c05 + E2[i]*s05/sqrt_rho_c ) + uc3[i]/3.;
        E4[i]  = -(1./3.)*( -uc[i]*s1*sqrt_rho_c + E[i]*c1 ) + (1./3.)*( -uc1[i]*s05*sqrt_rho_c + E1[i]*c05 ) + (2./3.)*( -uc2[i]*s05*sqrt_rho_c + E2[i]*c05 ) + E3[i]/3. - (1./6.)*dt*J[i];
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh4[k][i] = -(1./3.)*hfh[k][i]*std::exp(-I*kx[i]*v[k]*dt) + (1./3.)*hfh1[k][i]*std::exp(-0.5*I*kx[i]*v[k]*dt) + (2./3.)*hfh2[k][i]*std::exp(-0.5*I*kx[i]*v[k]*dt) + hfh3[k][i]/3. - (1./6.)*dt*d[i];
        }
      }
    } // end stage 4

    // STAGE 5
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh4[k].begin(),hfh4[k].end(),fh[k].begin()); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E4);

      double c1  = std::cos(dt*sqrt_rho_c)     , s1  = std::sin(dt*sqrt_rho_c)    ,
             c05 = std::cos(0.5*dt*sqrt_rho_c) , s05 = std::sin(0.5*dt*sqrt_rho_c);
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc5[i] = -(1./5.)*(  uc[i]*c1 + E[i]*s1/sqrt_rho_c ) + (1./5.)*(  uc1[i]*c05 + E1[i]*s05/sqrt_rho_c ) + (2./5.)*(  uc2[i]*c05 + E2[i]*s05/sqrt_rho_c ) + uc3[i]/5. + (2./5.)*uc4[i];
        E5[i]  = -(1./5.)*( -uc[i]*s1*sqrt_rho_c + E[i]*c1 ) + (1./5.)*( -uc1[i]*s05*sqrt_rho_c + E1[i]*c05 ) + (2./5.)*( -uc2[i]*s05*sqrt_rho_c + E2[i]*c05 ) + E3[i]/5. + (2./5.)*E4[i] - (1./10.)*dt*J[i];
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh5[k][i] = -(1./5.)*hfh[k][i]*std::exp(-I*kx[i]*v[k]*dt) + (1./5.)*hfh1[k][i]*std::exp(-0.5*I*kx[i]*v[k]*dt) + (2./5.)*hfh2[k][i]*std::exp(-0.5*I*kx[i]*v[k]*dt) + hfh3[k][i]/5. + (2./5.)*hfh4[k][i] - 0.1*dt*d[i];
        }
      }
    } // end stage 5

///////////////////////////////////////////////////////////////////////////////
// MONITORING /////////////////////////////////////////////////////////////////

    // SAVE TIME STEP
    std::copy(  uc4.begin() ,  uc4.end() ,  uc.begin() );
    std::copy(   E4.begin() ,   E4.end() ,   E.begin() );
    std::copy( hfh4.begin() , hfh4.end() , hfh.begin() );

    for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh[k].begin(),hfh[k].end(),fh[k].begin()); }
    double total_energy = energy(fh,E);
    {
      auto rhoh = fh.density();
      fft::spectrum_ hrhoh(Nx); hrhoh.fft(&rhoh[0]);
      fft::spectrum_ hE(Nx); hE.fft(&E[0]);
      fft::spectrum_ hrhoc(Nx);
      hrhoc[0] = I*kx[0]*hE[0] - hrhoh[0] + 1.;
      for ( auto i=1 ; i<Nx ; ++i )
        { hrhoc[i] = I*kx[i]*hE[i] - hrhoh[i]; }
      ublas::vector<double> rhoc (Nx,0.); hrhoc.ifft(rhoc.begin());

      for ( auto i=0 ; i<Nx ; ++i )
        { total_energy += rhoc[i]*uc[i]*uc[i]; }
    }
    H.push_back( total_energy );

    // increment time
    current_time += dt;

    ++iter;
    if ( current_time+dt > Tf ) { dt = Tf - current_time; }
  } // while (  current_time < Tf ) // end of time loop

  for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh[k].begin(),hfh[k].end(),fh[k].begin()); }
  fh.write("vp.dat");

  double h = std::abs(*std::max_element( H.begin() , H.end() , [&](double a,double b){return ( std::abs((a-H[0])/std::abs(H[0])) < std::abs((b-H[0])/std::abs(H[0])) );} ));
  return std::abs((h-H[0])/std::abs(H[0]));
}

int
main ( int argc , char const * argv[] )
{
  const std::size_t Nx = 75 , Nv = 1024;
  const double dt_max = 3.*16./Nv;

  for ( auto i=1 ; i<6 ; ++i ) {
    double dt = dt_max/double(i);
    std::cout << dt << " " << std::flush;
    double h = error_H(Nx,Nv,dt);
    std::cout << h << std::endl;
  }

  return 0;
}
