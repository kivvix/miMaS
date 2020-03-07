#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include <valarray>
#include <sstream>
#include <complex>
#include <tuple>

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

struct iter_s {
  std::size_t iter;
  double dt;
  double current_time;
  double Lhfh;
  double LE;
};

#define save(data,dir,suffix,x_y) {\
  std::stringstream filename; filename << #data << "_" << suffix << ".dat"; \
  std::ofstream of( dir / filename.str() );\
  std::transform( data.begin() , data.end() , std::ostream_iterator<std::string>(of,"\n") , x_y );\
  of.close();\
}

namespace math = boost::math::constants;
const std::complex<double> & I = std::complex<double>(0.,1.);

#define SQ(X) ((X)*(X))
#define Xi(i) (i*fh.step.dx+fh.range.x_min)
#define Vk(k) (k*fh.step.dv+fh.range.v_min)

#define ping(X) std::cerr << __LINE__ << " " << #X << ":" << X << std::endl
int debug = 0;

auto
maxwellian ( double rho , double u , double T ) {
  //std::cout << rho << " " << u << " " << T << std::endl;
  //std::cout << rho/(std::sqrt(2.*math::pi<double>()*T)) << std::endl;
  return [=](double x,double v){ return rho/(std::sqrt(2.*math::pi<double>()*T))*std::exp( -0.5*SQ(v-u)/T ); };
}

int
main(int argc, char const *argv[])
{

  std::filesystem::path p("config.init");
  if ( argc > 1 ) {
    p = argv[1];
  }
  auto c = config(p);

  std::filesystem::create_directories(c.output_dir);
  std::cout << c << std::endl;
  std::ofstream oconfig( c.output_dir / "config.init" );
  oconfig << c << std::endl;
  oconfig.close();

/* --------------------------------------------------------------- */
  std::size_t Nx = c.Nx, Nv = c.Nv;

  // $(u_c,E,\hat{f}_h)$ and $f_h$
  ublas::vector<double> uc(Nx,0.);
  ublas::vector<double> E (Nx,0.);
  field<double,1> fh(boost::extents[Nv][Nx]);
  complex_field<double,1> hfh(boost::extents[Nv][Nx]);
  ublas::vector<double> uc1(Nx) , uc2(Nx) , uc3(Nx) , uc4(Nx) , ucn(Nx),
                        E1 (Nx) , E2 (Nx) , E3 (Nx) , E4 (Nx) , En (Nx);
  complex_field<double,1> hfh1(boost::extents[Nv][Nx]) , hfh2(boost::extents[Nv][Nx]) ,
                          hfh3(boost::extents[Nv][Nx]) , hfh4(boost::extents[Nv][Nx]) ,
                          hfhn(boost::extents[Nv][Nx]) ;

  const double Kx = 0.5;
  // phase-space domain
  fh.range.v_min = -8.; fh.range.v_max = 8.;
  fh.range.x_min =  0.; fh.range.x_max = 2./Kx*math::pi<double>();

  // compute dx, dv
  fh.step.dv = (fh.range.v_max-fh.range.v_min)/Nv;
  fh.step.dx = (fh.range.x_max-fh.range.x_min)/Nx;

  double dt = 0.5*fh.step.dv;
  
  // velocity and frequency
  ublas::vector<double> v (Nv,0.); for ( std::size_t k=0 ; k<Nv ; ++k ) { v[k] = Vk(k); }
  const double l = fh.range.x_max-fh.range.x_min;
  ublas::vector<double> kx(Nx);
  for ( auto i=0 ; i<Nx/2 ; ++i ) { kx[i]    = 2.*math::pi<double>()*i/l; }
  for ( int i=-Nx/2 ; i<0 ; ++i ) { kx[Nx+i] = 2.*math::pi<double>()*i/l; }

  // initial condition
  auto tb_M1 = maxwellian(0.5*c.alpha,c.ui,1.) , tb_M2 = maxwellian(0.5*c.alpha,-c.ui,1.);
  for (field<double,2>::size_type k=0 ; k<fh.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<fh.size(1) ; ++i ) {
      //fh[k][i] = ( 0.5*c.alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-c.ui)) + 0.5*c.alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)+c.ui)) )*(1.+0.04*std::cos(0.3*Xi(i)));
      //fh[k][i] = ( 0.5*c.alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-c.ui)) + 0.5*c.alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)+c.ui)) )*(1.+0.04*std::cos(Kx*Xi(i)));

      fh[k][i] = ( tb_M1(Xi(i),Vk(k)) + tb_M2(Xi(i),Vk(k)) )*(1. + 0.01*std::cos(Kx*Xi(i)));
    }
    fft::fft(&(fh[k][0]),&(fh[k][Nx-1])+1,&(hfh[k][0]));
  }
  fh.write( c.output_dir / "init_vhls.dat" );
/*
  std::cout << "Nx: " << Nx << "\n";
  std::cout << "Nv: " << Nv << "\n";
  std::cout << "v_min: " << fh.range.v_min << "\n";
  std::cout << "v_max: " << fh.range.v_max << "\n";
  std::cout << "x_min: " << fh.range.x_min << "\n";
  std::cout << "x_max: " << fh.range.x_max << "\n";
  std::cout << "dt: " << dt << "\n";
  std::cout << "dx: " << fh.step.dx << "\n";
  std::cout << "dv: " << fh.step.dv << "\n";
  std::cout << "Tf: " << c.Tf << "\n";
  std::cout << "f_0: " << "\"tb\"" << "\n";
  std::cout << std::endl;
*/
  const double rho_c = 1.-c.alpha;
  // init E (electric field) with Poisson
  {
    poisson<double> poisson_solver(Nx,l);
    ublas::vector<double> rho(Nx,0.);
    rho = fh.density(); // compute density from init data
    for ( auto i=0 ; i<Nx ; ++i ) { rho[i] += rho_c; } // add (1-alpha) for cold particules
    E = poisson_solver(rho);
  }

  // monitoring data
  std::vector<iter_s> iterations; iterations.reserve(int(std::ceil(c.Tf/dt))+1);
  std::vector<iter_s> success_iterations; success_iterations.reserve(int(std::ceil(c.Tf/dt))+1);
  std::vector<double> ee;
  std::vector<double> Emax;
  std::vector<double> H;
  std::vector<double> times;

  auto save_data = [&] ( std::string suffix ) {
    save(iterations,c.output_dir,suffix,[](auto const& it) { std::stringstream ss; ss << it.iter << " " << it.dt << " " << it.current_time << " " << it.Lhfh << " " << it.LE; return ss.str(); })
    save(success_iterations,c.output_dir,suffix,[](auto const& it) { std::stringstream ss; ss << it.iter << " " << it.dt << " " << it.current_time << " " << it.Lhfh << " " << it.LE; return ss.str(); })

    auto dt_y = [&,count=0](auto const& y) mutable { std::stringstream ss; ss<<times[count++]<<" "<<y; return ss.str(); };
    save(ee,c.output_dir,suffix,dt_y);
    save(Emax,c.output_dir,suffix,dt_y);
    save(H,c.output_dir,suffix,dt_y);
  };


  signals_handler::signals_handler<SIGINT,SIGILL>::handler( [&]( int signal ) -> void {
    std::cerr << "\n\033[41;97m ** End of execution after signal " << signal << " ** \033[0m\n";
    std::cerr << "\033[36msave data...\033[0m\n";

    save_data("vhls_SIGINT");
    for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh[k].begin(),hfh[k].end(),fh[k].begin()); }
    fh.write( c.output_dir / "vp_vhls_SIGINT.dat" );
  });

  /*
  signal_handler::signal_handler<SIGINT>::function_handler = [&](int signal) {
    std::cerr << "\n\033[41;97m ** End of execution after signal " << signal << " ** \033[0m\n";
    std::cerr << "\033[36msave data...\033[0m\n" ;

    save_data("vhls_SIGINT");
    for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh[k][0]),&(hfh[k][0])+c.Nx,&(fh[k][0])); }
    fh.write( c.output_dir / "vp_vhls_SIGINT.dat" );
  };
  signal_handler::signal_handler<SIGINT>::signal();
  */

  times.push_back(0);
  {
    double electric_energy = 0.;
    for ( const auto & ei : E ) { electric_energy += ei*ei*fh.step.dx; }
    ee.push_back( std::sqrt(electric_energy) );
  }
  Emax.push_back( std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );

  //U_type<double,1> U({uc,E,hfh};
  splitting<double,1> Lie( fh , l , rho_c );

  std::size_t i_t = 0;
  double current_time = 0.;
  const double alpha1=1./(4.-std::cbrt(4.)), alpha2=alpha1, alpha3=1./(1.-SQ(std::cbrt(4.)));

  const double g1 = alpha1, g2 = alpha1+alpha2;
  const double w1 = (g2*(1.-g2))/(g1*(g1-1.)-g2*(g2-1.)) , w2 = 1.-w1 , w3 = w2 , w4 = w1;

  while ( current_time < c.Tf ) {
    std::cout<<"\r ["<<std::setw(6)<<i_t<<"] "<< std::setw(8) << current_time << " (" << std::setw(9) << dt << ")"<<std::flush;
    /**
    // Strang classique
    Lie.phi_a(0.5*dt,uc,E,hfh);
    Lie.phi_b(0.5*dt,uc,E,hfh);
    Lie.phi_c(dt,uc,E,hfh);
    Lie.phi_b(0.5*dt,uc,E,hfh);
    Lie.phi_a(0.5*dt,uc,E,hfh);
    /**/

    std::copy( hfh.origin() , hfh.origin()+hfh.num_elements() , hfhn.origin() );
    std::copy(  E.begin() ,  E.end() ,  En.begin() );
    std::copy( uc.begin() , uc.end() , ucn.begin() );

    // Strang alpha1
    Lie.phi_a(alpha1*0.5*dt,uc,E,hfh);
    Lie.phi_c(alpha1*0.5*dt,uc,E,hfh);
    Lie.phi_b(alpha1*dt,uc,E,hfh);
    Lie.phi_c(alpha1*0.5*dt,uc,E,hfh);
    Lie.phi_a(alpha1*0.5*dt,uc,E,hfh);

    std::copy( hfh.origin() , hfh.origin()+hfh.num_elements() , hfh1.origin() );
    std::copy(  E.begin() ,  E.end() ,  E1.begin() );
    std::copy( uc.begin() , uc.end() , uc1.begin() );

    // Strang alpha2
    Lie.phi_a(alpha2*0.5*dt,uc,E,hfh);
    Lie.phi_c(alpha2*0.5*dt,uc,E,hfh);
    Lie.phi_b(alpha2*dt,uc,E,hfh);
    Lie.phi_c(alpha2*0.5*dt,uc,E,hfh);
    Lie.phi_a(alpha2*0.5*dt,uc,E,hfh);

    std::copy( hfh.origin() , hfh.origin()+hfh.num_elements() , hfh2.origin() );
    std::copy(  E.begin() ,  E.end() ,  E2.begin() );
    std::copy( uc.begin() , uc.end() , uc2.begin() );

    // Strang alpha3
    Lie.phi_a(alpha3*0.5*dt,uc,E,hfh);
    Lie.phi_c(alpha3*0.5*dt,uc,E,hfh);
    Lie.phi_b(alpha3*dt,uc,E,hfh);
    Lie.phi_c(alpha3*0.5*dt,uc,E,hfh);
    Lie.phi_a(alpha3*0.5*dt,uc,E,hfh);

    std::copy( hfh.origin() , hfh.origin()+hfh.num_elements() , hfh3.origin() );
    std::copy(  E.begin() ,  E.end() ,  E3.begin() );
    std::copy( uc.begin() , uc.end() , uc3.begin() );

    // Strang alpha2
    Lie.phi_a(alpha2*0.5*dt,uc,E,hfh);
    Lie.phi_c(alpha2*0.5*dt,uc,E,hfh);
    Lie.phi_b(alpha2*dt,uc,E,hfh);
    Lie.phi_c(alpha2*0.5*dt,uc,E,hfh);
    Lie.phi_a(alpha2*0.5*dt,uc,E,hfh);

    std::copy( hfh.origin() , hfh.origin()+hfh.num_elements() , hfh4.origin() );
    std::copy(  E.begin() ,  E.end() ,  E4.begin() );
    std::copy( uc.begin() , uc.end() , uc4.begin() );

    // Strang alpha1
    Lie.phi_a(alpha1*0.5*dt,uc,E,hfh);
    Lie.phi_c(alpha1*0.5*dt,uc,E,hfh);
    Lie.phi_b(alpha1*dt,uc,E,hfh);
    Lie.phi_c(alpha1*0.5*dt,uc,E,hfh);
    Lie.phi_a(alpha1*0.5*dt,uc,E,hfh);

    /**/

    double L_hfh = 0.;
    for ( auto k=0 ; k<Nv ; ++k ) {
      for ( auto i=0 ; i<Nx ; ++i ) {
        auto hfhtile_ik = -hfhn[k][i] + w1*(hfh1[k][i]+hfh4[k][i]) + w2*(hfh2[k][i]+hfh3[k][i]);
        L_hfh += SQ( std::abs( hfh[k][i]-hfhtile_ik ) )*fh.step.dx*fh.step.dv;
      }
    }
    L_hfh = std::sqrt(L_hfh);

    double L_E = 0.;
    for ( auto i=0 ; i<Nx ; ++i ) {
      double Etile_ik = -En[i] + w1*(E1[i]+E4[i]) + w2*(E2[i]+E3[i]);
      L_E += SQ( std::abs( E[i]-Etile_ik ) )*fh.step.dx;
    }
    L_E = std::sqrt(L_E);

    std::cout << " -- " << std::setw(10) << L_hfh << "        " << std::flush;

    if ( std::abs(L_hfh - c.tol) <= c.tol ) // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    {
      // SAVE TIME STEP

      Emax.push_back( std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );

      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh[k].begin(),hfh[k].end(),fh[k].begin()); }
      double electric_energy = 0.;
      for ( const auto & ei : E ) { electric_energy += ei*ei*fh.step.dx; }
      ee.push_back( std::sqrt(electric_energy) );
      double total_energy = energy(fh,E);
      {
        auto rhoh = fh.density();
        fft::spectrum_ hrhoh(c.Nx); hrhoh.fft(&rhoh[0]);
        fft::spectrum_ hE(c.Nx); hE.fft(&E[0]);
        fft::spectrum_ hrhoc(c.Nx);
        hrhoc[0] = I*kx[0]*hE[0] - hrhoh[0] + 1.;
        for ( auto i=1 ; i<c.Nx ; ++i ) {
          hrhoc[i] = I*kx[i]*hE[i] - hrhoh[i];
        }
        ublas::vector<double> rhoc (c.Nx,0.); hrhoc.ifft(&rhoc[0]);

        for ( auto i=0 ; i<c.Nx ; ++i ) {
          total_energy += rhoc[i]*uc[i]*uc[i];
        }
      }
      H.push_back( total_energy );

      current_time += dt;
      times.push_back( current_time );
      success_iterations.push_back( { i_t , dt , current_time , L_hfh } );
    }
    else {
      // REMAKE THE STEP
      std::copy( hfhn.origin() , hfhn.origin()+hfhn.num_elements() , hfh.origin() );
      std::copy( En.begin()  , En.end()  , E.begin() );
      std::copy( ucn.begin() , ucn.end() , uc.begin() );
    }

    iterations.push_back( { i_t , dt , current_time , L_hfh } );
    ++i_t;

    dt = std::pow( c.tol/L_hfh , 0.25 )*dt;
    if ( current_time+dt > c.Tf ) { dt = c.Tf - current_time; }
  } // while current_time < c.Tf
  std::cout<<" ["<<std::setw(5)<<i_t<<"] "<<i_t*dt<<std::endl;

  std::ofstream of;
  std::size_t count = 0;

  /*
  auto dt_y = [&,count=0](auto const& y) mutable { std::stringstream ss; ss<<times[count++]<<" "<<y; return ss.str(); };

  of.open( c.output_dir / "ee_vhls.dat" );
  std::transform( ee.begin() , ee.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();

  of.open( c.output_dir / "Emax_vhls.dat" );
  std::transform( Emax.begin() , Emax.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();
  
  of.open( c.output_dir / "H_vhls.dat" );
  std::transform( H.begin() , H.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();
  */
  save_data("vhls_suzuki");

  for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh[k][0]),&(hfh[k][Nx-1])+1,&(fh[k][0])); }
  fh.write( c.output_dir / "vp_vhls_suzuki.dat" );

  auto dx_y = [&,count=0](auto const& y) mutable { std::stringstream ss; ss<< fh.step.dx*(count++) <<" "<<y; return ss.str(); };
  save(E,c.output_dir,"vhls",dx_y);
  save(uc,c.output_dir,"vhls",dx_y);

  return 0;
}

