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
#include "miMaS/splitting.h"
#include "miMaS/lagrange5.h"
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

int
main ( int argc , char const * argv[] )
{
  std::filesystem::path p("config.init");
  if ( argc > 1 )
    { p = argv[1]; }
  auto c = config(p);
  c.name = "vhls";

  c.create_output_directory();
  std::ofstream ofconfig( c.output_dir / "config.init" );
  ofconfig << c << "\n";
  ofconfig.close();

/* ------------------------------------------------------------------------- */
  field<double,1> fh(boost::extents[c.Nv][c.Nx]);
  complex_field<double,1> hfh(boost::extents[c.Nv][c.Nx]);

  const double Kx = 0.5;
  fh.range.v_min = -8.; fh.range.v_max = 8.;
  fh.range.x_min =  0.; fh.range.x_max = 2./Kx*math::pi<double>();
  fh.compute_steps();

  ublas::vector<double> v(c.Nv,0.);
  std::generate( v.begin() , v.end() , [&,k=0]() mutable {return (k++)*fh.step.dv+fh.range.v_min;} );

  ublas::vector<double> kx(c.Nx); // beware, Nx need to be odd
  {
    double l = fh.range.len_x();
    for ( auto i=0 ; i<c.Nx/2 ; ++i ) { kx[i]      = 2.*math::pi<double>()*i/l; }
    for ( int i=-c.Nx/2 ; i<0 ; ++i ) { kx[c.Nx+i] = 2.*math::pi<double>()*i/l; }
  }

  auto tb_M1 = maxwellian( 0.5*c.alpha , c.ui , 1. ) , tb_M2 = maxwellian( 0.5*c.alpha , -c.ui , 1. );
  for (field<double,2>::size_type k=0 ; k<fh.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<fh.size(1) ; ++i ) {
      fh[k][i] = ( tb_M1(Xi(i),Vk(k)) + tb_M2(Xi(i),Vk(k)) )*(1. + 0.01*std::cos(Kx*Xi(i)));
    }
    fft::fft(fh[k].begin(),fh[k].end(),hfh[k].begin());
  }
  fh.write( c.output_dir / ("init_"+c.name+".dat") );

  iteration::iteration<double> iter;
  iter.iter = 0;
  iter.current_time = 0.;
  iter.dt = 0.5*fh.step.dv;

  ublas::vector<double> uc(c.Nx,0.);
  ublas::vector<double> E (c.Nx,0.);

  std::vector<double> ee;   ee.reserve(int(std::ceil(c.Tf/iter.dt))+1);
  std::vector<double> Emax; Emax.reserve(int(std::ceil(c.Tf/iter.dt))+1);
  std::vector<double> H;    H.reserve(int(std::ceil(c.Tf/iter.dt))+1);
  std::vector<iteration::iteration<double>> iterations; iterations.reserve(int(std::ceil(c.Tf/iter.dt))+1);
  std::vector<iteration::iteration<double>> success_iterations; success_iterations.reserve(int(std::ceil(c.Tf/iter.dt))+1);

  std::vector<double> times; times.reserve(int(std::ceil(c.Tf/iter.dt))+1);

  // lambda to save data, any time you want, where you want
  auto save_data = [&] ( std::string && suffix ) -> void {
    suffix = c.name + suffix;

    // save iterations informations
    auto writer_iter = [] ( auto const & it ) {
      std::stringstream ss; ss << it;
      return ss.str();
    };
    c << monitoring::data( "iterations_"+suffix+".dat"         , iterations         , writer_iter );
    c << monitoring::data( "success_iterations_"+suffix+".dat" , success_iterations , writer_iter );

    // save temporel data
    auto dt_y = [&,count=0] (auto const& y) mutable {
      std::stringstream ss; ss<<times[count++]<<" "<<y;
      return ss.str();
    };
    c << monitoring::data( "ee_"+suffix+".dat"   , ee   , dt_y );
    c << monitoring::data( "Emax_"+suffix+".dat" , Emax , dt_y );
    c << monitoring::data( "H_"+suffix+".dat"    , H    , dt_y );

    // save distribution function
    for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh[k].begin(),hfh[k].end(),fh[k].begin()); }
    fh.write( c.output_dir / ("vp_"+suffix+".dat") );
  };

  // to stop simulation at any time (and save data)
  signal_handler::signal_handler<SIGINT,SIGILL>::handler( [&]( int signal ) -> void {
    std::cerr << "\n\033[41;97m ** End of execution after signal " << signal << " ** \033[0m\n";
    std::cerr << "\033[38;5;202msave data...\033[0m\n";
    save_data("_SIGINT");
  });

  const double rho_c = 1.-c.alpha;

  // time 0
  // init E with Poisson solver, init also ee, Emax, H and times
  {
    poisson<double> poisson_solver(c.Nx,fh.range.len_x());
    ublas::vector<double> rho(c.Nx,0.);
    rho = fh.density(); // compute density from init data
    for ( auto i=0 ; i<c.Nx ; ++i ) { rho[i] += (1.-c.alpha); } // add (1-alpha) for cold particules
    E = poisson_solver(rho);

    Emax.push_back( std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );
    double electric_energy = std::sqrt(std::accumulate(
      E.begin() , E.end() , 0. ,
      [&] ( double partial_sum , double ei ) {
        return partial_sum + ei*ei*fh.step.dx;
      }
    ));
    ee.push_back( electric_energy );

    double total_energy = energy(fh,E);
    total_energy += 0.; // sum(rho_c*u_c*u_c) = 0 because u_c = 0 at time 0
    H.push_back( total_energy );
    times.push_back(0.);
  }

  // initialize memory for all temporary variables
  ublas::vector<double> J(c.Nx,0.);
  fft::spectrum_ d(c.Nx);
  field<double,1> Edvf(tools::array_view<const std::size_t>(fh.shape(),2));
  ublas::vector<double> ucn(c.Nx) , uctilde(c.Nx) , uc1(c.Nx) , uc2(c.Nx) , uc3(c.Nx) , uc4(c.Nx) ,
                        En (c.Nx) , Etilde (c.Nx) , E1 (c.Nx) , E2(c.Nx)  , E3 (c.Nx) , E4 (c.Nx) ;
  complex_field<double,1> hfhn(boost::extents[c.Nv][c.Nx]) , hfhtilde(boost::extents[c.Nv][c.Nx]) ,
                          hfh1(boost::extents[c.Nv][c.Nx]) , hfh2(boost::extents[c.Nv][c.Nx]) ,
                          hfh3(boost::extents[c.Nv][c.Nx]) , hfh4(boost::extents[c.Nv][c.Nx]) ;

  //U_type<double,1> U({uc,E,hfh};
  splitting<double,1> Lie( fh , fh.range.len_x() , rho_c );

  std::size_t i_t = 0;
  double current_time = 0.;
  const double alpha1=1./(4.-std::cbrt(4.)), alpha2=alpha1, alpha3=1./(1.-SQ(std::cbrt(4.)));

  const double g1 = alpha1, g2 = alpha1+alpha2;
  const double w1 = (g2*(1.-g2))/(g1*(g1-1.)-g2*(g2-1.)) , w2 = 1.-w1 , w3 = w2 , w4 = w1;

  while (  iter.current_time < c.Tf ) {
    std::cout << "\r" << iteration::time(iter) << std::flush;

    std::copy( hfh.origin() , hfh.origin()+hfh.num_elements() , hfhn.origin() );
    std::copy(  E.begin() ,  E.end() ,  En.begin() );
    std::copy( uc.begin() , uc.end() , ucn.begin() );

    // Strang alpha1
    Lie.phi_a(alpha1*0.5*iter.dt,uc,E,hfh);
    Lie.phi_c(alpha1*0.5*iter.dt,uc,E,hfh);
    Lie.phi_b(alpha1*iter.dt,uc,E,hfh);
    Lie.phi_c(alpha1*0.5*iter.dt,uc,E,hfh);
    Lie.phi_a(alpha1*0.5*iter.dt,uc,E,hfh);

    std::copy( hfh.origin() , hfh.origin()+hfh.num_elements() , hfh1.origin() );
    std::copy(  E.begin() ,  E.end() ,  E1.begin() );
    std::copy( uc.begin() , uc.end() , uc1.begin() );

    // Strang alpha2
    Lie.phi_a(alpha2*0.5*iter.dt,uc,E,hfh);
    Lie.phi_c(alpha2*0.5*iter.dt,uc,E,hfh);
    Lie.phi_b(alpha2*iter.dt,uc,E,hfh);
    Lie.phi_c(alpha2*0.5*iter.dt,uc,E,hfh);
    Lie.phi_a(alpha2*0.5*iter.dt,uc,E,hfh);

    std::copy( hfh.origin() , hfh.origin()+hfh.num_elements() , hfh2.origin() );
    std::copy(  E.begin() ,  E.end() ,  E2.begin() );
    std::copy( uc.begin() , uc.end() , uc2.begin() );

    // Strang alpha3
    Lie.phi_a(alpha3*0.5*iter.dt,uc,E,hfh);
    Lie.phi_c(alpha3*0.5*iter.dt,uc,E,hfh);
    Lie.phi_b(alpha3*iter.dt,uc,E,hfh);
    Lie.phi_c(alpha3*0.5*iter.dt,uc,E,hfh);
    Lie.phi_a(alpha3*0.5*iter.dt,uc,E,hfh);

    std::copy( hfh.origin() , hfh.origin()+hfh.num_elements() , hfh3.origin() );
    std::copy(  E.begin() ,  E.end() ,  E3.begin() );
    std::copy( uc.begin() , uc.end() , uc3.begin() );

    // Strang alpha2
    Lie.phi_a(alpha2*0.5*iter.dt,uc,E,hfh);
    Lie.phi_c(alpha2*0.5*iter.dt,uc,E,hfh);
    Lie.phi_b(alpha2*iter.dt,uc,E,hfh);
    Lie.phi_c(alpha2*0.5*iter.dt,uc,E,hfh);
    Lie.phi_a(alpha2*0.5*iter.dt,uc,E,hfh);

    std::copy( hfh.origin() , hfh.origin()+hfh.num_elements() , hfh4.origin() );
    std::copy(  E.begin() ,  E.end() ,  E4.begin() );
    std::copy( uc.begin() , uc.end() , uc4.begin() );

    // Strang alpha1
    Lie.phi_a(alpha1*0.5*iter.dt,uc,E,hfh);
    Lie.phi_c(alpha1*0.5*iter.dt,uc,E,hfh);
    Lie.phi_b(alpha1*iter.dt,uc,E,hfh);
    Lie.phi_c(alpha1*0.5*iter.dt,uc,E,hfh);
    Lie.phi_a(alpha1*0.5*iter.dt,uc,E,hfh);

    for ( auto i=0 ; i<c.Nx ; ++i ) {
       Etilde[i] = - En[i] + w1*( E1[i]+ E4[i]) + w2*( E2[i]+ E3[i]);
      uctilde[i] = -ucn[i] + w1*(uc1[i]+uc4[i]) + w2*(uc2[i]+uc3[i]);
    }
    for ( auto k=0 ; k<c.Nv ; ++k ) {
      for ( auto i=0 ; i<c.Nx ; ++i ) {
         hfhtilde[k][i] = - hfhn[k][i] + w1*( hfh1[k][i]+ hfh4[k][i]) + w2*( hfh2[k][i]+ hfh3[k][i]);
      }
    }

///////////////////////////////////////////////////////////////////////////////
// MONITORING /////////////////////////////////////////////////////////////////

    // CHECH local error for compute next time step
    iter.E_error(E,Etilde,fh.step.dx);
    iter.hfh_error(hfh,hfhtilde,fh.step.dx*fh.step.dv);
    iter.success = std::abs(iter.error() - c.tol) <= c.tol;

    std::cout << " -- " << iteration::error(iter) << std::flush;
    if ( iter.success )
    {
      Emax.push_back( std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );
      double electric_energy = std::sqrt(std::accumulate(
        E.begin() , E.end() , 0. ,
        [&] ( double partial_sum , double ei ) {
          return partial_sum + ei*ei*fh.step.dx;
        }
      ));
      ee.push_back( electric_energy );

      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh[k].begin(),hfh[k].end(),fh[k].begin()); }
      double total_energy = energy(fh,E);
      {
        auto rhoh = fh.density();
        fft::spectrum_ hrhoh(c.Nx); hrhoh.fft(rhoh.begin());
        fft::spectrum_ hE(c.Nx); hE.fft(E.begin());
        fft::spectrum_ hrhoc(c.Nx);
        hrhoc[0] = I*kx[0]*hE[0] - hrhoh[0] + 1.;
        for ( auto i=1 ; i<c.Nx ; ++i ) {
          hrhoc[i] = I*kx[i]*hE[i] - hrhoh[i];
        }
        ublas::vector<double> rhoc (c.Nx,0.); hrhoc.ifft(rhoc.begin());

        for ( auto i=0 ; i<c.Nx ; ++i ) {
          total_energy += rhoc[i]*uc[i]*uc[i];
        }
      }
      H.push_back( total_energy );

      // increment time
      iter.current_time += iter.dt;
      times.push_back( iter.current_time );
      success_iterations.push_back( iter );
    } else {
      // RESET THE STEP
      std::copy( hfhn.origin() , hfhn.origin()+hfhn.num_elements() , hfh.origin() );
      std::copy( En.begin()  , En.end()  , E.begin() );
      std::copy( ucn.begin() , ucn.end() , uc.begin() );
    }
    iterations.push_back( iter );


    ++iter.iter;
    iter.dt = std::pow( c.tol/iter.Lhfh , 0.25 )*iter.dt;
    if ( iter.current_time+iter.dt > c.Tf ) { iter.dt = c.Tf - iter.current_time; }
  } // while (  iter.current_time < c.Tf ) // end of time loop
  std::cout << "\r" << time(iter) << std::endl;

  save_data("suzuki");

  auto dx_y = [&,count=0](auto const& y) mutable {
    std::stringstream ss; ss<< fh.step.dx*(count++) <<" "<<y;
    return ss.str();
  };
  c << monitoring::data( "E_"+c.name+".dat" , E , dx_y );

  ublas::vector<double> rho (c.Nx,0.);
  ublas::vector<double> rhoc(c.Nx,0.);
  {
    ublas::vector<double> rhoh = fh.density();
    fft::spectrum_ hrhoh(c.Nx); hrhoh.fft(rhoh.begin());
    fft::spectrum_ hE(c.Nx);    hE.fft(E.begin());

    fft::spectrum_ hrho(c.Nx), hrhoc(c.Nx);

    hrho[0] = I*kx[0]*hE[0] + 1.;
    hrho[0] = I*kx[0]*hE[0] - hrhoh[0] + 1.;
    for ( auto i=1 ; i<c.Nx ; ++i ) {
      hrho[i]  = I*kx[i]*hE[i];
      hrhoc[i] = I*kx[i]*hE[i] - hrhoh[i];
    }
    hrho.ifft(rho.begin());
  }

  c << monitoring::data( "rho_"+c.name+".dat"  , rho  , dx_y );
  c << monitoring::data( "uc_"+c.name+".dat"   , uc   , dx_y );
  c << monitoring::data( "rhoc_"+c.name+".dat" , rhoc , dx_y );

  J = fh.courant();
  for ( auto i=0 ; i<c.Nx ; ++i ) {
    J[i] += rhoc[i]*uc[i];
  }
  c << monitoring::data( "J_"+c.name+".dat" , J , dx_y );

  return 0;
}