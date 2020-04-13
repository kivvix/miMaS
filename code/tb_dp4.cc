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

int
main ( int argc , char const * argv[] )
{
  std::filesystem::path p("config.init");
  if ( argc > 1 )
    { p = argv[1]; }
  auto c = config(p);
  c.name = "vhll";

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
  const double sqrt_rho_c = std::sqrt(rho_c);

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
  ublas::vector<double> uc1(c.Nx) , uc2(c.Nx) , uc3(c.Nx) , uc4(c.Nx) , uc5(c.Nx) , uc6(c.Nx) , uc7(c.Nx),
                        E1 (c.Nx) , E2(c.Nx)  , E3 (c.Nx) , E4 (c.Nx) , E5 (c.Nx) , E6 (c.Nx) , E7 (c.Nx);
  complex_field<double,1> hfh1(boost::extents[c.Nv][c.Nx]) , hfh2(boost::extents[c.Nv][c.Nx]) ,
                          hfh3(boost::extents[c.Nv][c.Nx]) , hfh4(boost::extents[c.Nv][c.Nx]) ,
                          hfh5(boost::extents[c.Nv][c.Nx]) , hfh6(boost::extents[c.Nv][c.Nx]) ,
                          hfh7(boost::extents[c.Nv][c.Nx]) ;

  while (  iter.current_time < c.Tf ) {
    std::cout << "\r" << iteration::time(iter) << std::flush;

///////////////////////////////////////////////////////////////////////////////
// DP4(3) /////////////////////////////////////////////////////////////////////

    // STAGE 1
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh[k].begin(),hfh[k].end(),fh[k].begin()); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E);

      double c05 = std::cos(0.5*iter.dt*sqrt_rho_c), s05 = std::sin(0.5*iter.dt*sqrt_rho_c);
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc1[i] =  uc[i]*c05 + E[i]*s05/sqrt_rho_c - 0.5*iter.dt*J[i]*s05/sqrt_rho_c;
        E1[i]  = -uc[i]*s05*sqrt_rho_c + E[i]*c05 - 0.5*iter.dt*J[i]*c05;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh1[k][i] = hfh[k][i]*std::exp(-0.5*I*kx[i]*v[k]*iter.dt) - 0.5*iter.dt*d[i]*std::exp(-0.5*I*kx[i]*v[k]*iter.dt);
        }
      }
    } // end stage 1

    // STAGE 2
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh1[k].begin(),hfh1[k].end(),fh[k].begin()); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E1);

      double c05 = std::cos(0.5*iter.dt*sqrt_rho_c), s05 = std::sin(0.5*iter.dt*sqrt_rho_c);
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc2[i] =  uc[i]*c05 + E[i]*s05/sqrt_rho_c;
        E2[i]  = -uc[i]*s05*sqrt_rho_c + E[i]*c05 - 0.5*iter.dt*J[i];
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh2[k][i] = hfh[k][i]*std::exp(-0.5*I*kx[i]*v[k]*iter.dt) - 0.5*iter.dt*d[i];
        }
      }
    } // end stage 2

    // STAGE 3
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh2[k].begin(),hfh2[k].end(),fh[k].begin()); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E2);

      double c1  = std::cos(iter.dt*sqrt_rho_c)     , s1  = std::sin(iter.dt*sqrt_rho_c)    ,
             c05 = std::cos(0.5*iter.dt*sqrt_rho_c) , s05 = std::sin(0.5*iter.dt*sqrt_rho_c);
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc3[i] =  uc[i]*c1 + E[i]*s1/sqrt_rho_c - iter.dt*J[i]*s05/sqrt_rho_c;
        E3[i]  = -uc[i]*s1*sqrt_rho_c + E[i]*c1 - iter.dt*J[i]*c05;
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh3[k][i] = hfh[k][i]*std::exp(-I*kx[i]*v[k]*iter.dt) - iter.dt*d[i]*std::exp(-0.5*I*kx[i]*v[k]*iter.dt);
        }
      }
    } // end stage 3

    // STAGE 4
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh3[k].begin(),hfh3[k].end(),fh[k].begin()); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E3);

      double c1  = std::cos(iter.dt*sqrt_rho_c)     , s1  = std::sin(iter.dt*sqrt_rho_c)    ,
             c05 = std::cos(0.5*iter.dt*sqrt_rho_c) , s05 = std::sin(0.5*iter.dt*sqrt_rho_c);
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc4[i] = -(1./3.)*(  uc[i]*c1 + E[i]*s1/sqrt_rho_c ) + (1./3.)*(  uc1[i]*c05 + E1[i]*s05/sqrt_rho_c ) + (2./3.)*(  uc2[i]*c05 + E2[i]*s05/sqrt_rho_c ) + uc3[i]/3.;
        E4[i]  = -(1./3.)*( -uc[i]*s1*sqrt_rho_c + E[i]*c1 ) + (1./3.)*( -uc1[i]*s05*sqrt_rho_c + E1[i]*c05 ) + (2./3.)*( -uc2[i]*s05*sqrt_rho_c + E2[i]*c05 ) + E3[i]/3. - (1./6.)*iter.dt*J[i];
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh4[k][i] = -(1./3.)*hfh[k][i]*std::exp(-I*kx[i]*v[k]*iter.dt) + (1./3.)*hfh1[k][i]*std::exp(-0.5*I*kx[i]*v[k]*iter.dt) + (2./3.)*hfh2[k][i]*std::exp(-0.5*I*kx[i]*v[k]*iter.dt) + hfh3[k][i]/3. - (1./6.)*iter.dt*d[i];
        }
      }
    } // end stage 4

    // STAGE 5
    {
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(hfh4[k].begin(),hfh4[k].end(),fh[k].begin()); }
      J = fh.courant();
      Edvf = weno::trp_v(fh,E4);

      double c1  = std::cos(iter.dt*sqrt_rho_c)     , s1  = std::sin(iter.dt*sqrt_rho_c)    ,
             c05 = std::cos(0.5*iter.dt*sqrt_rho_c) , s05 = std::sin(0.5*iter.dt*sqrt_rho_c);
      for ( auto i=0 ; i<c.Nx ; ++i ) {
        uc5[i] = -(1./5.)*(  uc[i]*c1 + E[i]*s1/sqrt_rho_c ) + (1./5.)*(  uc1[i]*c05 + E1[i]*s05/sqrt_rho_c ) + (2./5.)*(  uc2[i]*c05 + E2[i]*s05/sqrt_rho_c ) + uc3[i]/5. + (2./5.)*uc4[i];
        E5[i]  = -(1./5.)*( -uc[i]*s1*sqrt_rho_c + E[i]*c1 ) + (1./5.)*( -uc1[i]*s05*sqrt_rho_c + E1[i]*c05 ) + (2./5.)*( -uc2[i]*s05*sqrt_rho_c + E2[i]*c05 ) + E3[i]/5. + (2./5.)*E4[i] - (1./10.)*iter.dt*J[i];
      }
      for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) {
        d.fft(&(Edvf[k][0]));
        for ( auto i=0 ; i<c.Nx ; ++i ) {
          hfh5[k][i] = -(1./5.)*hfh[k][i]*std::exp(-I*kx[i]*v[k]*iter.dt) + (1./5.)*hfh1[k][i]*std::exp(-0.5*I*kx[i]*v[k]*iter.dt) + (2./5.)*hfh2[k][i]*std::exp(-0.5*I*kx[i]*v[k]*iter.dt) + hfh3[k][i]/5. + (2./5.)*hfh4[k][i] - 0.1*iter.dt*d[i];
        }
      }
    } // end stage 5

///////////////////////////////////////////////////////////////////////////////
// MONITORING /////////////////////////////////////////////////////////////////

    // CHECH local error for compute next time step
    iter.E_error(E5,E4,fh.step.dx);
    iter.hfh_error(hfh5,hfh4,fh.step.dx*fh.step.dv);
    iter.success = std::abs(iter.error() - c.tol) <= c.tol;

    std::cout << " -- " << iteration::error(iter) << std::flush;
    if ( iter.success )
    {
      // SAVE TIME STEP
      std::copy( uc4.begin()  , uc4.end()  , uc.begin()  );
      std::copy( E4.begin()   , E4.end()   , E.begin()   );
      std::copy( hfh4.begin() , hfh4.end() , hfh.begin() );

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
    }
    iterations.push_back( iter );


    ++iter.iter;
    iter.dt = std::pow( c.tol/iter.Lhfh , 0.25 )*iter.dt;
    if ( iter.current_time+iter.dt > c.Tf ) { iter.dt = c.Tf - iter.current_time; }
  } // while (  iter.current_time < c.Tf ) // end of time loop
  std::cout << "\r" << time(iter) << std::endl;

  save_data("dp4");

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
