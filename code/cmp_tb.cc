#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include <valarray>
#include <sstream>
#include <filesystem>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "miMaS/field.h"
#include "miMaS/weno.h"
#include "miMaS/fft.h"
#include "miMaS/array_view.h"
#include "miMaS/poisson.h"
#include "miMaS/rk.h"
#include "miMaS/config.h"

namespace o2 {
  template < typename _T , std::size_t NumDimsV >
  auto
  trp_v ( field<_T,NumDimsV> const & u , ublas::vector<_T> const& E )
  {
    field<_T,NumDimsV> trp(tools::array_view<const std::size_t>(u.shape(),NumDimsV+1));

    { auto k=0, km1=trp.size(0)-1;
      for ( auto i=0 ; i<trp.size(1) ; ++i ) {
        trp[k][i] = ( E(i)*(u[k+1][i]-u[km1][i])/(2.*u.step.dv) );
      }
    }
    for ( auto k=1 ; k<trp.size(0)-1 ; ++k ) {
      for ( auto i=0 ; i<trp.size(1) ; ++i ) {
        trp[k][i] = ( E(i)*(u[k+1][i]-u[k-1][i])/(2.*u.step.dv) );
      }
    }
    { auto k=trp.size(0)-1, kp1=0;
      for ( auto i=0 ; i<trp.size(1) ; ++i ) {
        trp[k][i] = ( E(i)*(u[kp1][i]-u[k-1][i])/(2.*u.step.dv) );
      }
    }

    return trp;
  }
}


namespace math = boost::math::constants;
const std::complex<double> & I = std::complex<double>(0.,1.);

#define SQ(X) ((X)*(X))
#define Xi(i) (i*f.step.dx+f.range.x_min)
#define Vk(k) (k*f.step.dv+f.range.v_min)


auto
maxwellian ( double rho , double u , double T ) {
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
  std::ofstream oconfig( c.output_dir / "config.init" );
  oconfig << c << std::endl;
  oconfig.close();

/* --------------------------------------------------------------- */
  field<double,1> f(boost::extents[c.Nv][c.Nx]);
  f.range.v_min = -8.; f.range.v_max = 8.;
  f.step.dv = (f.range.v_max-f.range.v_min)/c.Nv;

  double Kx = 0.5;
  f.range.x_min = 0.; f.range.x_max = 2./Kx*math::pi<double>();
  //f.range.x_min = 0.; f.range.x_max = 20.*math::pi<double>();
  f.step.dx = (f.range.x_max-f.range.x_min)/c.Nx;

  ublas::vector<double> v (c.Nv,0.);
  for ( std::size_t k=0 ; k<c.Nv ; ++k ) { v[k] = Vk(k); }

  ublas::vector<double> kx(c.Nx);
  {
    double l = f.range.len_x();
    for ( auto i=0 ; i<c.Nx/2 ; ++i ) { kx[i]      = 2.*math::pi<double>()*i/l; }
    for ( int i=-c.Nx/2 ; i<0 ; ++i ) { kx[c.Nx+i] = 2.*math::pi<double>()*i/l; }
  }

  double alpha = 0.2 , ui = 2.0;
  auto landau_M = maxwellian( 1.        ,  0. , 1.   );
  auto db_M1    = maxwellian( 0.5       , -ui , 1.   ) , db_M2  = maxwellian( 0.5       ,  ui , 1.  );
  //auto bot_M1   = maxwellian( 1.-alpha  ,  0. , 1.   ) , bot_M2 = maxwellian( alpha     ,  ui , 0.25);
  auto bot_M1   = maxwellian(0.9,0.,1.) , bot_M2 = maxwellian(0.2,4.5,0.25);
  auto tb_MC    = maxwellian( 1.-alpha  ,  0. , c.Tc ) ,
       tb_M1    = maxwellian( 0.5*alpha ,  ui , 1.   ) , tb_M2  = maxwellian( 0.5*alpha , -ui , 1.  );
  auto v10_MC   = maxwellian( 1.-alpha  ,  0. , c.Tc ) , v10_Mh = maxwellian( alpha     ,  0. , 1.  );

  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      //f[k][i] = ((1.-alpha)*M0c(Xi(i),Vk(k)) + 0.5*alpha*( Mpu(Xi(i),Vk(k)) + Mmu(Xi(i),Vk(k)) ))*(1.+0.01*std::cos(0.5*Xi(i)));

      //// landau damping : Kx=0.5
      //f[k][i] = landau_M(Xi(i),Vk(k))*(1.+0.001*std::cos(Kx*Xi(i)));
      //// double beam Kx=0.2, ui=2.4 ou Kx=0.2, ui=4.5
      //f[k][i] = (db_M1(Xi(i),Vk(k))+db_M2(Xi(i),Vk(k)))*(1.+0.001*std::cos(Kx*Xi(i)));
      //// bot Kx=0.5 , alpha=0.2 , ui=4.5
      //f[k][i] = (bot_M1(Xi(i),Vk(k)) + bot_M2(Xi(i),Vk(k)))*(1.+0.04*std::cos(0.3*Xi(i)));
      //// tb Kx=0.5 , ui=4. , alpha=0.2 , Tc=0.01
      f[k][i] = tb_MC(Xi(i),Vk(k)) + (tb_M1(Xi(i),Vk(k)) + tb_M2(Xi(i),Vk(k)) )*(1. + 0.01*std::cos(Kx*Xi(i)));
      //// v10 Kx=0.5 , alpha=0.2
      //f[k][i] = v10_MC(Xi(i),Vk(k)) + ( std::pow(Vk(k),10)*v10_Mh(Xi(i),Vk(k))/945. )*(1. + 0.01*std::cos(Kx*Xi(i)));
    }
  }
  f.write( c.output_dir / "init.dat" );

  unsigned int i_t = 0;
  double current_time = 0.;
  double dt = 1.433*f.step.dv/0.6;

  std::vector<double> ee;   ee.reserve(int(std::ceil(c.Tf/dt))+1);
  std::vector<double> Emax; Emax.reserve(int(std::ceil(c.Tf/dt))+1);
  std::vector<double> H;    H.reserve(int(std::ceil(c.Tf/dt))+1);
  std::vector<double> Ec;   Ec.reserve(int(std::ceil(c.Tf/dt))+1);

  std::vector<double> times; times.reserve(int(std::ceil(c.Tf/dt))+1);

  // space scheme
  auto weno = [&](field<double,1>const& f , ublas::vector<double> const& E )->field<double,1> { return weno::trp_v(f,E); };
  auto cd2  = [&](field<double,1>const& f , ublas::vector<double> const& E )->field<double,1> { return o2::trp_v(f,E); };

  // time scheme init
  lawson::RK33<poisson<double>> rk(c.Nx,c.Nv,f.range.len_x(),f.shape(),v,kx,weno);

  rk.E = rk.poisson_solver(f.density());
  {
    Emax.push_back( std::abs(*std::max_element( rk.E.begin() , rk.E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );
    double electric_energy = 0.;
    for ( const auto & ei : rk.E ) { electric_energy += ei*ei*f.step.dx; }
    ee.push_back( std::sqrt(electric_energy) );
    H.push_back( energy(f,rk.E) );
    Ec.push_back( kinetic_energy(f) );
    times.push_back( 0. );
  }


  while (  current_time < c.Tf ) {
    std::cout<<" ["<<std::setw(5)<<i_t<<"] "<< current_time <<"\r"<<std::flush;
    
    f = rk(f,dt);

    // end of time loop

    // MONITORING
    rk.E = rk.poisson_solver(f.density());
    Emax.push_back( std::abs(*std::max_element( rk.E.begin() , rk.E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );
    double electric_energy = 0.;
    for ( const auto & ei : rk.E ) { electric_energy += ei*ei*f.step.dx; }
    ee.push_back( std::sqrt(electric_energy) );
    H.push_back( energy(f,rk.E) );
    Ec.push_back( kinetic_energy(f) );

    //dt = std::min( 0.1 , SIGMA*f.step.dv/Emax[i_t] );

    // increment time
    ++i_t;
    current_time += dt;
    times.push_back( current_time );
  } // while (  i_t*dt < Tf )
  std::cout<<" ["<<std::setw(5)<<i_t<<"] "<<i_t*dt <<std::endl;

  f.write( c.output_dir / "vp.dat" );

  auto dt_y = [&,count=0](auto const& y) mutable { std::stringstream ss; ss<<times[count++]<<" "<<y; return ss.str(); };

  std::ofstream of;

  of.open( c.output_dir / "ee.dat" );
  std::transform( ee.begin() , ee.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();
  
  of.open( c.output_dir / "Emax.dat" );
  std::transform( Emax.begin() , Emax.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();
  
  of.open( c.output_dir / "H.dat" );
  std::transform( H.begin() , H.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();
  
  of.open( c.output_dir / "Ec.dat" );
  std::transform( Ec.begin() , Ec.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();

  auto dx_y = [&,count=0](auto const& y) mutable { std::stringstream ss; ss<< f.step.dx*(count++) <<" "<<y; return ss.str(); };

  of.open( c.output_dir / "E.dat" );
  std::transform( rk.E.begin() , rk.E.end() , std::ostream_iterator<std::string>(of,"\n") , dx_y );
  of.close();

  auto rho = f.density();
  of.open( c.output_dir / "rho.dat" );
  std::transform( rho.begin() , rho.end() , std::ostream_iterator<std::string>(of,"\n") , dx_y );
  of.close();

  auto J = f.courant();
  of.open( c.output_dir / "J.dat" );
  std::transform( J.begin() , J.end() , std::ostream_iterator<std::string>(of,"\n") , dx_y );
  of.close();

  of.open( c.output_dir / "energy.dat" );
  for ( auto i=0 ; i<times.size() ; ++i ) {
    of << times[i] << " " << ee[i] << " " << Ec[i] << " " << H[i] << " " << ee[i] + Ec[i] << "\n";
  }
  of.close();

  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      f[k][i] -= tb_MC(Xi(i),Vk(k)) + (tb_M1(Xi(i),Vk(k)) + tb_M2(Xi(i),Vk(k)) )*(1. + 0.01*std::cos(Kx*Xi(i)));
    }
  }
  f.write( c.output_dir / "diff.dat" );

  return 0;
}
