#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include <valarray>
#include <sstream>
#include <complex>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "miMaS/field.h"
#include "miMaS/complex_field.h"
#include "miMaS/weno.h"
#include "miMaS/fft.h"
#include "miMaS/array_view.h"
#include "miMaS/poisson.h"

namespace math = boost::math::constants;
const std::complex<double> & I = std::complex<double>(0.,1.);

#define SQ(X) ((X)*(X))
#define Xi(i) (i*fh.step.dx+fh.range.x_min)
#define Vk(k) (k*fh.step.dv+fh.range.v_min)

//#define ping(X) std::cout << __LINE__ << " " << #X << ":" << X << std::endl
//int debug = 0;

auto
maxwellian ( double rho , double u , double T ) {
  //std::cout << rho << " " << u << " " << T << std::endl;
  //std::cout << rho/(std::sqrt(2.*math::pi<double>()*T)) << std::endl;
  return [=](double x,double v){ return rho/(std::sqrt(2.*math::pi<double>()*T))*std::exp( -0.5*SQ(v-u)/T ); };
}

int main(int,char**)
{
  std::size_t Nx = 81, Nv = 128;

  // $(u_c,E,\hat{f}_h)$ and $f_h$
  ublas::vector<double> uc(Nx,0.);
  ublas::vector<double> E (Nx,0.);
  field<double,1> fh(boost::extents[Nv][Nx]);
  complex_field<double,1> hfh(boost::extents[Nv][Nx]);

  const double Kx = 0.5;
  // phase-space domain
  fh.range.v_min = -8.; fh.range.v_max = 8.;
  //fh.range.x_min =  0.; fh.range.x_max = 20.*math::pi<double>();
  fh.range.x_min =  0.; fh.range.x_max = 2./Kx*math::pi<double>();

  // compute dx, dv
  fh.step.dv = (fh.range.v_max-fh.range.v_min)/Nv;
  fh.step.dx = (fh.range.x_max-fh.range.x_min)/Nx;

  double dt = 1.5*fh.step.dv;
  double Tf = 20.;
  
  // velocity and frequency
  ublas::vector<double> v (Nv,0.); for ( std::size_t k=0 ; k<Nv ; ++k ) { v[k] = Vk(k); }
  const double l = fh.range.x_max-fh.range.x_min;
  ublas::vector<double> kx(Nx);
  for ( auto i=0 ; i<Nx/2 ; ++i ) { kx[i]    = 2.*math::pi<double>()*i/l; }
  for ( int i=-Nx/2 ; i<0 ; ++i ) { kx[Nx+i] = 2.*math::pi<double>()*i/l; }

  // initial condition
  double ui=2., alpha=0.2;
  auto tb_M1 = maxwellian(0.5*alpha,ui,1.) , tb_M2 = maxwellian(0.5*alpha,-ui,1.);
  auto v10_Mh = maxwellian(alpha,0.,1.);
  for (field<double,2>::size_type k=0 ; k<fh.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<fh.size(1) ; ++i ) {
      //fh[k][i] = ( 0.5*alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)) + 0.5*alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)+ui)) )*(1.+0.04*std::cos(0.3*Xi(i)));
      //fh[k][i] = ( 0.5*alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)) + 0.5*alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)+ui)) )*(1.+0.04*std::cos(Kx*Xi(i)));

      // tb
      fh[k][i] = ( tb_M1(Xi(i),Vk(k)) + tb_M2(Xi(i),Vk(k)) )*(1. + 0.01*std::cos(Kx*Xi(i)));
      // v10
      //fh[k][i] = ( std::pow(Vk(k),10)*v10_Mh(Xi(i),Vk(k))/945. )*(1. + 0.01*std::cos(Kx*Xi(i)));
    }
    fft::fft(&(fh[k][0]),&(fh[k][Nx-1])+1,&(hfh[k][0]));
  }
  fh.write("vphl/init.dat");

  std::cout << "Nx: " << Nx << "\n";
  std::cout << "Nv: " << Nv << "\n";
  std::cout << "v_min: " << fh.range.v_min << "\n";
  std::cout << "v_max: " << fh.range.v_max << "\n";
  std::cout << "x_min: " << fh.range.x_min << "\n";
  std::cout << "x_max: " << fh.range.x_max << "\n";
  std::cout << "dt: " << dt << "\n";
  std::cout << "dx: " << fh.step.dx << "\n";
  std::cout << "dv: " << fh.step.dv << "\n";
  std::cout << "Tf: " << Tf << "\n";
  std::cout << "f_0: " << "\"tb\"" << "\n";
  std::cout << std::endl;

  const double rho_c = 1.-alpha;
  const double sqrt_rho_c = std::sqrt(rho_c);
  // init E (electric field) with Poisson
  poisson<double> poisson_solver(Nx,l);
  ublas::vector<double> rho(Nx,0.);
  rho = fh.density(); // compute density from init data
  for ( auto i=0 ; i<Nx ; ++i ) { rho[i] += rho_c; } // add (1-alpha) for cold particules
  E = poisson_solver(rho);

  // monitoring data
  std::vector<double> ee;
  std::vector<double> Emax;
  std::vector<double> times;

  times.push_back(0);
  {
    double electric_energy = 0.;
    for ( const auto & ei : E ) { electric_energy += ei*ei*fh.step.dx; }
    ee.push_back( std::sqrt(electric_energy) );
  }
  Emax.push_back( std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );

  // initialize memory for all temporary variables
  ublas::vector<double> J(Nx);
  fft::spectrum_ d(Nx);
  field<double,1> Edvf(tools::array_view<const std::size_t>(fh.shape(),2));
  ublas::vector<double> uc1(Nx) , uc2(Nx) , uc3(Nx) , uc4(Nx) , uc5(Nx) , uc6(Nx),
                        E1(Nx)  , E2(Nx)  , E3(Nx)  , E4(Nx)  , E5(Nx)  , E6(Nx) ;
  complex_field<double,1> hfh1(boost::extents[Nv][Nx]),hfh2(boost::extents[Nv][Nx]),hfh3(boost::extents[Nv][Nx]),hfh4(boost::extents[Nv][Nx]),hfh5(boost::extents[Nv][Nx]),hfh6(boost::extents[Nv][Nx]);

  std::size_t i_t = 0;
  double current_time = 0.;
  while ( current_time < Tf ) {
    std::cout << " [" << std::setw(5) << i_t << "] " << i_t*dt << "\r" << std::flush;

/*
    // STAGE 1
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh[k][0]),&(hfh[k][0])+Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E);

      double c = std::cos(dt*sqrt_rho_c), s = std::sin(dt*sqrt_rho_c);
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc1[i] =  uc[i]*c + E[i]*s/sqrt_rho_c - dt*J[i]*s;
        E1[i]  = -uc[i]*s*sqrt_rho_c + E[i]*c - dt*J[i]*c;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh1[k][i] = hfh[k][i]*std::exp(-I*kx[i]*Vk(k)*dt) - dt*d[i]*std::exp(-I*kx[i]*Vk(k)*dt);
        }
      }
    } // end stage 1

    // STAGE 2
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh1[k][0]),&(hfh1[k][0])+Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E1);

      double c2 = std::cos(0.5*dt*sqrt_rho_c), s2 = std::sin(0.5*dt*sqrt_rho_c);
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc2[i] =  0.75*uc[i]*c2 + 0.75*E[i]*s2/sqrt_rho_c + 0.25*uc1[i]*c2 - 0.25*E1[i]*s2/sqrt_rho_c + 0.25*dt*J[i]*s2;
        E2[i]  = -0.75*uc[i]*s2*sqrt_rho_c + 0.75*E[i]*c2 + 0.25*uc1[i]*s2*sqrt_rho_c + 0.25*E1[i]*c2 - 0.25*dt*J[i]*c2;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh2[k][i] = 0.75*hfh[k][i]*std::exp(-I*kx[i]*Vk(k)*0.5*dt) + 0.25*hfh1[k][i]*std::exp(I*kx[i]*Vk(k)*0.5*dt) - 0.25*dt*d[i]*std::exp(I*kx[i]*Vk(k)*0.5*dt);
        }
      }
    } // end stage 2

    // STAGE 3
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh2[k][0]),&(hfh2[k][0])+Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E2);

      double c = std::cos(dt*sqrt_rho_c), s = std::sin(dt*sqrt_rho_c);
      double c2 = std::cos(0.5*dt*sqrt_rho_c), s2 = std::sin(0.5*dt*sqrt_rho_c);
      for ( auto i=0 ; i<Nx ; ++i ) {
        double tmp_uc =  (1./3.)*uc[i]*c + (1./3.)*E[i]*s/sqrt_rho_c + (2./3.)*uc2[i]*c2 + (2./3.)*E2[i]*s2/sqrt_rho_c - (2./3.)*dt*J[i]*s2;
        E[i]          = -(1./3.)*uc[i]*s*sqrt_rho_c + (1./3.)*E[i]*c - (2./3.)*uc2[i]*s2*sqrt_rho_c + (2./3.)*E2[i]*c2 - (2./3.)*dt*J[i]*c2;
        uc[i] = tmp_uc;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh[k][i] = (1./3.)*hfh[k][i]*std::exp(-I*kx[i]*Vk(k)*dt) + (2./3.)*hfh2[k][i]*std::exp(-I*kx[i]*Vk(k)*0.5*dt) - (2./3.)*dt*d[i]*std::exp(-I*kx[i]*Vk(k)*0.5*dt);
        }
      }
    } // end stage 3
*/
    // STAGE 1
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh[k][0]),&(hfh[k][0])+Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E);

      double c02 = std::cos(0.2*dt*sqrt_rho_c) , s02 = std::sin(0.2*dt*sqrt_rho_c);
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc1[i] =  uc[i]*c02 + E[i]*s02/sqrt_rho_c + 0.2*dt*J[i]*s02/sqrt_rho_c;
        E1[i]  = -uc[i]*s02*sqrt_rho_c + E[i]*c02 + 0.2*dt*J[i]*c02;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh1[k][i] = hfh[k][i]*std::exp(-0.2*I*kx[i]*Vk(k)*dt) + 0.2*dt*d[i]*std::exp(-0.2*I*kx[i]*Vk(k)*dt);
        }
      }
    } // end stage 1

    // STAGE 2
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh1[k][0]),&(hfh1[k][0])+Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E1);

      double c01 = std::cos(0.1*sqrt_rho_c*dt) , s01 = std::sin(0.1*sqrt_rho_c*dt) , c03 = std::cos(0.3*sqrt_rho_c*dt) , s03 = std::sin(0.3*sqrt_rho_c*dt);
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc2[i] = 0.625*(  uc[i]*c03 + E[i]*s03/sqrt_rho_c ) + 0.375*( uc1[i]*c01 + E1[i]*s01/sqrt_rho_c ) + 0.225*dt*J[i]*s01/sqrt_rho_c;
        E2[i]  = 0.625*( -uc[i]*s03*sqrt_rho_c + E[i]*c03 ) - 0.375*( sqrt_rho_c*uc1[i]*s01 + E1[i]*c01 ) + 0.225*dt*J[i]*c01;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh2[k][i] = 0.625*hfh[k][i]*std::exp(-0.3*I*kx[i]*Vk(k)*dt) + 0.375*hfh1[k][i]*std::exp(-0.1*I*kx[i]*Vk(k)*dt) + 0.225*dt*d[i]*std::exp(-0.1*I*kx[i]*Vk(k));
        }
      }
    } // end stage 2

    // STAGE 3
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh2[k][0]),&(hfh2[k][0])+Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E2);

      double c08 = std::cos(0.8*dt*sqrt_rho_c) , s08 = std::sin(0.8*dt*sqrt_rho_c) ,
             c06 = std::cos(0.6*dt*sqrt_rho_c) , s06 = std::sin(0.6*dt*sqrt_rho_c) ,
             c05 = std::cos(0.5*dt*sqrt_rho_c) , s05 = std::sin(0.5*dt*sqrt_rho_c) ;
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc3[i] = (175./27.)*(  uc[i]*c08 + E[i]*s08/sqrt_rho_c ) + (100./9.)*(  uc1[i]*c06 + E1[i]*s06/sqrt_rho_c ) - (448./27.)*(  uc2[i]*c05 + E2[i]*s05/sqrt_rho_c ) + (32./9.)*dt*J[i]*s05/sqrt_rho_c;
        E3[i]  = (175./27.)*( -uc[i]*s08*sqrt_rho_c + E[i]*c08 ) + (100./9.)*( -uc1[i]*s06*sqrt_rho_c + E1[i]*c06 ) - (448./27.)*( -uc2[i]*s05*sqrt_rho_c + E2[i]*c05 ) + (32./9.)*dt*J[i]*c05;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh3[k][i] = (175./27.)*hfh[k][i]*std::exp(-0.8*I*kx[i]*Vk(k)*dt) + (100./9.)*hfh1[k][i]*std::exp(-0.6*I*kx[i]*Vk(k)*dt) - (448./27.)*hfh2[k][i]*std::exp(-0.5*I*kx[i]*Vk(k)*dt) + (32./9.)*dt*d[i]*std::exp(-0.5*I*kx[i]*Vk(k));
        }
      }
    } // end stage 3

    // STAGE 4
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh3[k][0]),&(hfh3[k][0])+Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E3);

      double c89    = std::cos((8./9.)*sqrt_rho_c*dt)    , s89    = std::sin((8./9.)*sqrt_rho_c*dt)    ,
             c3145  = std::cos((31./45.)*sqrt_rho_c*dt)  , s3145  = std::sin((31./45.)*sqrt_rho_c*dt)  ,
             c5390  = std::cos((53./90.)*sqrt_rho_c*dt)  , s5390  = std::sin((53./90.)*sqrt_rho_c*dt)  ,
             c44827 = std::cos((448./27.)*sqrt_rho_c*dt) , s44827 = std::sin((448./27.)*sqrt_rho_c*dt) ,
             c445   = std::cos((4./45.)*sqrt_rho_c*dt)   , s445   = std::sin((4./45.)*sqrt_rho_c*dt)   ;
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc4[i] = (3551./6561.)*(  uc[i]*c89 + E[i]*s89/sqrt_rho_c ) + (7420./2187.)*(  uc1[i]*c3145 + E1[i]*s3145/sqrt_rho_c ) - (37376./6561.)*(  uc2[i]*c5390 + E2[i]*s5390/sqrt_rho_c ) + (2014./729.)*(  uc3[i]*c445 + E3[i]*s445/sqrt_rho_c ) - (212./729.)*dt*J[i]*s445/sqrt_rho_c;
        E4[i]  = (3551./6561.)*( -uc[i]*s89*sqrt_rho_c + E[i]*c89 ) + (7420./2187.)*( -uc1[i]*s3145*sqrt_rho_c + E1[i]*c3145 ) - (37376./6561.)*( -uc2[i]*s5390*sqrt_rho_c + E2[i]*c5390 ) + (2014./729.)*( -uc3[i]*s445*sqrt_rho_c + E3[i]*c445 ) - (212./729.)*dt*J[i]*c445;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh4[k][i] = (3551./6561.)*hfh[k][i]*std::exp(-(8./9.)*I*kx[i]*Vk(k)*dt) + (7420./2187.)*hfh1[k][i]*std::exp(-(31./45.)*I*kx[i]*Vk(k)*dt) - (37376./6561.)*hfh2[k][i]*std::exp(-(53./90.)*I*kx[i]*Vk(k)*dt) + (2014./729.)*hfh3[k][i]*std::exp(-(4./45.)*I*kx[i]*Vk(k)*dt) - (212./729.)*dt*d[i]*std::exp(-(4./45.)*kx[i]*Vk(k)*dt);
        }
      }
    } // end stage 4

    // STAGE 5
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh4[k][0]),&(hfh4[k][0])+Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E4);

      double c   = std::cos(sqrt_rho_c*dt)         , s   = std::sin(sqrt_rho_c*dt)         ,
             c08 = std::cos(0.8*sqrt_rho_c*dt)     , s08 = std::sin(0.8*sqrt_rho_c*dt)     ,
             c07 = std::cos(0.7*sqrt_rho_c*dt)     , s07 = std::sin(0.7*sqrt_rho_c*dt)     ,
             c02 = std::cos(0.2*sqrt_rho_c*dt)     , s02 = std::sin(0.2*sqrt_rho_c*dt)     ,
             c19 = std::cos((1./9.)*sqrt_rho_c*dt) , s19 = std::sin((1./9.)*sqrt_rho_c*dt) ;
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc5[i] = (313397./335808.)*(  uc[i]*c + E[i]*s/sqrt_rho_c ) + (424025./55968.)*(  uc1[i]*c08 + E1[i]*s08/sqrt_rho_c ) - (61400./5247.)*(  uc2[i]*c07 + E2[i]*s07/sqrt_rho_c ) + (96075./18656.)*(  uc3[i]*c02 + E3[i]*s02/sqrt_rho_c ) - (35721./37312.)*(  uc4[i]*c19 + E4[i]*s19/sqrt_rho_c ) - (5103./18656.)*dt*J[i]*s19/sqrt_rho_c;
        E5[i]  = (313397./335808.)*( -uc[i]*s*sqrt_rho_c + E[i]*c ) + (424025./55968.)*( -uc1[i]*s08*sqrt_rho_c + E1[i]*c08 ) - (61400./5247.)*( -uc2[i]*s07*sqrt_rho_c + E2[i]*c07 ) + (96075./18656.)*( -uc3[i]*s02*sqrt_rho_c + E3[i]*c02 ) - (35721./37312.)*( -uc4[i]*s19*sqrt_rho_c + E4[i]*c19 ) - (5103./18656.)*dt*J[i]*c19;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh5[k][i] = (313397./335808.)*hfh[k][i]*std::exp(-I*kx[i]*Vk(k)*dt) + (424025./55968.)*hfh1[k][i]*std::exp(-0.8*I*kx[i]*Vk(k)*dt) - (61400./5247.)*hfh2[k][i]*std::exp(-0.7*I*kx[i]*Vk(k)*dt) + (96075./18656.)*hfh3[k][i]*std::exp(-0.2*I*kx[i]*Vk(k)*dt) - (35721./37312.)*hfh4[k][i]*std::exp(-(1./9.)*I*kx[i]*Vk(k)*dt) - (5103./18656.)*dt*d[i]*std::exp(-(1./9.)*I*kx[i]*Vk(k)*dt);
        }
      }
    } // end stage 5

    // STAGE 6
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh5[k][0]),&(hfh5[k][0])+Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E5);

      double c   = std::cos(sqrt_rho_c*dt)         , s   = std::sin(sqrt_rho_c*dt)         ,
             c08 = std::cos(0.8*sqrt_rho_c*dt)     , s08 = std::sin(0.8*sqrt_rho_c*dt)     ,
             c07 = std::cos(0.7*sqrt_rho_c*dt)     , s07 = std::sin(0.7*sqrt_rho_c*dt)     ,
             c02 = std::cos(0.2*sqrt_rho_c*dt)     , s02 = std::sin(0.2*sqrt_rho_c*dt)     ,
             c19 = std::cos((1./9.)*sqrt_rho_c*dt) , s19 = std::sin((1./9.)*sqrt_rho_c*dt) ;
      for ( auto i=0 ; i<Nx ; ++i ) {
        uc6[i] = -(563./3456.)*(  uc[i]*c + E[i]*s/sqrt_rho_c ) - (575./252.)*(  uc1[i]*c08 + E1[i]*s08/sqrt_rho_c ) + (31400./10017.)*(  uc2[i]*c07 + E2[i]*s07/sqrt_rho_c ) + (325./1344.)*(  uc3[i]*c02 + E3[i]*s02/sqrt_rho_c ) - (7533./6784.)*(  uc4[i]*c19 + E4[i]*s19/sqrt_rho_c ) + (32./28.)*uc5[i];
        E6[i]  = -(563./3456.)*( -uc[i]*s*sqrt_rho_c + E[i]*c ) - (575./252.)*( -uc1[i]*s08*sqrt_rho_c + E1[i]*c08 ) + (31400./10017.)*( -uc2[i]*s07*sqrt_rho_c + E2[i]*c07 ) + (325./1344.)*( -uc3[i]*s02*sqrt_rho_c + E3[i]*c02 ) - (7533./6784.)*( -uc4[i]*s19*sqrt_rho_c + E4[i]*c19 ) + (32./28.)*E5[i] + (11./84.)*dt*J[i];
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh6[k][i] = -(563./3456.)*hfh[k][i]*std::exp(-I*kx[i]*Vk(k)*dt) - (575./252.)*hfh1[k][i]*std::exp(-0.8*I*kx[i]*Vk(k)*dt) + (31400./10017.)*hfh2[k][i]*std::exp(-0.7*I*kx[i]*Vk(k)*dt) + (325./1344.)*hfh3[k][i]*std::exp(-0.2*I*kx[i]*Vk(k)*dt) - (7533./6784.)*hfh4[k][i]*std::exp(-(1./9.)*I*kx[i]*Vk(k)*dt) + (32./28.)*hfh5[k][i] + (11./84.)*dt*d[i];
          hfh[i] = hfh6[i];
        }
      }
    } // end stage 6

    // STAGE 7
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh6[k][0]),&(hfh6[k][0])+Nx,&(fh[k][0])); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E6);

      double c   = std::cos(sqrt_rho_c*dt)         , s   = std::sin(sqrt_rho_c*dt)         ,
             c08 = std::cos(0.8*sqrt_rho_c*dt)     , s08 = std::sin(0.8*sqrt_rho_c*dt)     ,
             c07 = std::cos(0.7*sqrt_rho_c*dt)     , s07 = std::sin(0.7*sqrt_rho_c*dt)     ,
             c02 = std::cos(0.2*sqrt_rho_c*dt)     , s02 = std::sin(0.2*sqrt_rho_c*dt)     ,
             c19 = std::cos((1./9.)*sqrt_rho_c*dt) , s19 = std::sin((1./9.)*sqrt_rho_c*dt) ;
      for ( auto i=0 ; i<Nx ; ++i ) {
        double uc_tmp = (8813./172800.)*( uc[i]*c + E[i]*s/sqrt_rho_c ) + (41./180.)*(  uc1[i]*c08 + E1[i]*s08/sqrt_rho_c ) - (4294./10017.)*(  uc2[i]*c07 + E2[i]*s07/sqrt_rho_c ) + (263./384.)*(  uc3[i]*c02 + E3[i]*s02/sqrt_rho_c ) - (137781./339200.)*(  uc4[i]*c19 + E4[i]*s19/sqrt_rho_c ) + (803./4200)*uc5[i] + (17./25.)*uc6[i];
        E[i]          = (8813./172800.)*( uc[i]*c + E[i]*s/sqrt_rho_c ) + (41./180.)*( -uc1[i]*s08*sqrt_rho_c + E1[i]*c08 ) - (4294./10017.)*( -uc2[i]*s07*sqrt_rho_c + E2[i]*c07 ) + (263./384.)*( -uc3[i]*s02*sqrt_rho_c + E3[i]*c02 ) - (137781./339200.)*( -uc4[i]*s19*sqrt_rho_c + E4[i]*c19 ) + (803./4200)*E5[i] + (17./25.)*E6[i] + (1./40.)*dt*J[i];
        uc[i] = uc_tmp;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<Nx ; ++i ) {
          hfh[k][i] = (8813./172800.)*hfh[k][i]*std::exp(-I*kx[i]*Vk(k)*dt) + (41./180.)*hfh1[k][i]*std::exp(-0.8*I*kx[i]*Vk(k)*dt) - (4294./10017.)*hfh2[k][i]*std::exp(-0.7*I*kx[i]*Vk(k)*dt) + (263./384.)*hfh3[k][i]*std::exp(-0.2*I*kx[i]*Vk(k)*dt) - (137781./339200.)*hfh4[k][i]*std::exp(-(1./9.)*I*kx[i]*Vk(k)*dt) + (803./4200)*hfh5[k][i] + (17./25.)*hfh6[k][i] + (1./40.)*dt*d[i];
        }
      }
    } // end stage 7

    Emax.push_back( std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );

    double electric_energy = 0.;
    for ( const auto & ei : E ) { electric_energy += ei*ei*fh.step.dx; }
    ee.push_back( std::sqrt(electric_energy) );

    current_time += dt;
    ++i_t;
    times.push_back( current_time );
  } // while current_time < Tf
  std::cout<<" ["<<std::setw(5)<<i_t<<"] "<<i_t*dt<<std::endl;

  std::ofstream of;
  std::size_t count = 0;
  auto dt_y = [&,count=0](auto const& y) mutable { std::stringstream ss; ss<<times[count++]<<" "<<y; return ss.str(); };

  of.open("vphl/ee.dat");
  std::transform( ee.begin() , ee.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();

  of.open("vphl/Emax.dat");
  std::transform( Emax.begin() , Emax.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();

  for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh[k][0]),&(hfh[k][Nx-1])+1,&(fh[k][0])); }
  fh.write("vphl/vp.dat");


  return 0;
}
