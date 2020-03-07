#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include <valarray>
#include <sstream>
#include <filesystem>
#include <string>
#include <utility>

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
#include "miMaS/splitting.h"
#include "miMaS/lagrange5.h"
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

auto
error_H ( std::size_t Nx , std::size_t Nv , double dt )
{
  const double Tf = 6.5;

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

  const double alpha = 0.2 , ui = 3.4;
  auto tb_M1 = maxwellian( 0.5*alpha , ui , 1. ) , tb_M2 = maxwellian( 0.5*alpha , -ui , 1. );
  for (field<double,2>::size_type k=0 ; k<fh.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<fh.size(1) ; ++i ) {
      hfh[k][i] = 0.;
      fh[k][i] = ( tb_M1(Xi(i),Vk(k)) + tb_M2(Xi(i),Vk(k)) )*(1. + 0.01*std::cos(Kx*Xi(i)));
    }
    fft::fft(fh[k].begin(),fh[k].end(),hfh[k].begin());
  }

  std::size_t iter = 0;
  double current_time = 0.;

  ublas::vector<double> uc(Nx,0.);
  ublas::vector<double> E (Nx,0.);

  std::vector<double> H;
  std::vector<double> ee;
  std::vector<double> t;

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
  t.push_back(current_time);

  // initialize memory for all temporary variables
  splitting<double,1> Lie( fh , fh.range.len_x() , rho_c );

  ublas::vector<double> uc1(Nx) , uc2(Nx) , uc3(Nx) , uc4(Nx) , ucn(Nx),
                        E1 (Nx) , E2 (Nx) , E3 (Nx) , E4 (Nx) , En (Nx);
  complex_field<double,1> hfh1(boost::extents[Nv][Nx]) , hfh2(boost::extents[Nv][Nx]) ,
                          hfh3(boost::extents[Nv][Nx]) , hfh4(boost::extents[Nv][Nx]) ,
                          hfhn(boost::extents[Nv][Nx]) ;

  const double alpha1=1./(4.-std::cbrt(4.)), alpha2=alpha1, alpha3=1./(1.-SQ(std::cbrt(4.)));

  const double g1 = alpha1, g2 = alpha1+alpha2;
  const double w1 = (g2*(1.-g2))/(g1*(g1-1.)-g2*(g2-1.)) , w2 = 1.-w1 , w3 = w2 , w4 = w1;

  while (  current_time < Tf ) {
  /*std::cout << *std::max_element( hfh[100].begin() , hfh[100].end() , [](auto a,auto b){
    return std::abs(a) < std::abs(b);
  }) << " " << std::endl;*/

///////////////////////////////////////////////////////////////////////////////
// Suzuki /////////////////////////////////////////////////////////////////////

    /*    
    Lie.phi_a(0.5*dt,uc,E,hfh);
    Lie.phi_c(0.5*dt,uc,E,hfh);
    Lie.phi_b(dt,uc,E,hfh);
    Lie.phi_c(0.5*dt,uc,E,hfh);
    Lie.phi_a(0.5*dt,uc,E,hfh);
    */

    // Strang alpha1
    Lie.phi_a(alpha1*0.5*dt,uc,E,hfh);
    Lie.phi_c(alpha1*0.5*dt,uc,E,hfh);
    Lie.phi_b(alpha1*dt,uc,E,hfh);
    Lie.phi_c(alpha1*0.5*dt,uc,E,hfh);
    Lie.phi_a(alpha1*0.5*dt,uc,E,hfh);

    // Strang alpha2
    Lie.phi_a(alpha2*0.5*dt,uc,E,hfh);
    Lie.phi_c(alpha2*0.5*dt,uc,E,hfh);
    Lie.phi_b(alpha2*dt,uc,E,hfh);
    Lie.phi_c(alpha2*0.5*dt,uc,E,hfh);
    Lie.phi_a(alpha2*0.5*dt,uc,E,hfh);

    // Strang alpha3
    Lie.phi_a(alpha3*0.5*dt,uc,E,hfh);
    Lie.phi_c(alpha3*0.5*dt,uc,E,hfh);
    Lie.phi_b(alpha3*dt,uc,E,hfh);
    Lie.phi_c(alpha3*0.5*dt,uc,E,hfh);
    Lie.phi_a(alpha3*0.5*dt,uc,E,hfh);

    // Strang alpha2
    Lie.phi_a(alpha2*0.5*dt,uc,E,hfh);
    Lie.phi_c(alpha2*0.5*dt,uc,E,hfh);
    Lie.phi_b(alpha2*dt,uc,E,hfh);
    Lie.phi_c(alpha2*0.5*dt,uc,E,hfh);
    Lie.phi_a(alpha2*0.5*dt,uc,E,hfh);

    // Strang alpha1
    Lie.phi_a(alpha1*0.5*dt,uc,E,hfh);
    Lie.phi_c(alpha1*0.5*dt,uc,E,hfh);
    Lie.phi_b(alpha1*dt,uc,E,hfh);
    Lie.phi_c(alpha1*0.5*dt,uc,E,hfh);
    Lie.phi_a(alpha1*0.5*dt,uc,E,hfh);


///////////////////////////////////////////////////////////////////////////////
// MONITORING /////////////////////////////////////////////////////////////////

    for ( auto k=0 ; k<Nv ; ++k ) {
      fft::ifft(hfh[k].begin(),hfh[k].end(),fh[k].begin());
    }

    double total_energy = kinetic_energy(fh);
    for ( const auto ei : E ) {
      total_energy += ei*ei*fh.step.dx;
    }
    {
      auto rhoh = fh.density();
      fft::spectrum_ hrhoh(Nx); hrhoh.fft(&rhoh[0]);
      fft::spectrum_ hE(Nx); hE.fft(&E[0]);
      fft::spectrum_ hrhoc(Nx);
      hrhoc[0] = I*kx[0]*hE[0] - hrhoh[0] + 1.;
      for ( auto i=1 ; i<Nx ; ++i ) {
        hrhoc[i] = I*kx[i]*hE[i] - hrhoh[i];
      }
      ublas::vector<double> rhoc (Nx,0.); hrhoc.ifft(&rhoc[0]);

      for ( auto i=0 ; i<Nx ; ++i ) {
        total_energy += rhoc[i]*uc[i]*uc[i]*fh.step.dx;
      }
    }
    H.push_back( total_energy );

    double electric_energy = 0.;
    for ( const auto & ei : E ) { electric_energy += ei*ei*fh.step.dx; }
    ee.push_back( std::sqrt(electric_energy) );

    // increment time
    current_time += dt;
    t.push_back(current_time);


    ++iter;
    //if ( current_time+dt > Tf ) { dt = Tf - current_time; }
  } // while (  current_time < Tf ) // end of time loop

  std::ofstream of("H.dat");
  std::transform( H.cbegin() , H.cend() , std::ostream_iterator<std::string>(of,"\n") ,
    [&,count=0] ( const double h ) mutable {
      std::stringstream ss;
      ss << t[count++] << " " << h << " " << std::abs(H[0]-h)/std::abs(H[0]);
      return ss.str();
  });
  of.close();

  of.open("ee.dat");
  std::transform( ee.cbegin() , ee.cend() , std::ostream_iterator<std::string>(of,"\n") ,
    [&,count=0] ( const double e ) mutable {
      std::stringstream ss;
      ss << t[count++] << " " << e;
      return ss.str();
  });
  of.close();

  for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh[k].begin(),hfh[k].end(),fh[k].begin()); }
  fh.write("vp.dat");
  
  double h_max = std::abs(*std::max_element( H.begin() , H.end() ,
    [&] ( double a , double b ) {
      return std::abs(a-H[0])/std::abs(H[0]) < std::abs(b-H[0])/std::abs(H[0]);
    }
  ));
  
  double h_last = H.back();
  return std::make_pair(std::abs(H[0]-h_last)/std::abs(H[0]),std::abs(H[0]-h_max)/std::abs(H[0]));
}

int
main ( int argc , char const * argv[] )
{
  const std::size_t Nx = 243 ;//, Nv = 2048;
  //const double dt_ref = 1.;

  std::vector<double> h;
  /*
  for ( auto i=1 ; i<5 ; i++ ) {
    double dt = dt_ref/double(i);
    std::cout << dt << " ";
    h.push_back(error_H(Nx,Nv,dt));
    std::cout << h[i-1] << std::endl;
  }
  */
  
  std::ofstream of("err_h.dat");
  for ( auto dt : {2.,1.,0.5,0.25,0.1,0.05,0.01,0.005} ) {
  //for ( auto dt : {0.05} ) {
    std::cout << dt << " " << std::flush;
    double tmp = error_H(Nx,Nv,dt);
    std::cout << tmp << std::endl;
    of << dt << " " << tmp << std::endl;
    h.push_back(tmp);
  }
  of.close();
  /*
  std::ofstream of("err_h.dat");
  for ( auto Nv : {64,128,256,512,1024,2048} ) {
    // cfl = dv/dt = cste = 0.5
    double dt = 2.*16./Nv;
    std::cout << dt << " " << std::flush;
    auto tmp = error_H(Nx,Nv,dt);
    std::cout << tmp.first << " " << tmp.second << std::endl;
    of << dt << " " << tmp.first << " " << tmp.second << std::endl;
  }
  of.close();
  */
  return 0;
}
