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
#include "miMaS/rk.h"
#include "miMaS/config.h"
#include "miMaS/signal_handler.h"
#include "miMaS/iteration.h"



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

struct surimi {
  ublas::vector<double> E;

  surimi ( std::size_t nx , double )
    : E(nx,1.)
  { ; }

  ~surimi ()
  { ; }

  ublas::vector<double>
  operator () ( ublas::vector<double> const& )
  {
    return E;
  }
};

auto
disque ( double x0 , double v0 , double R ) {
  return [=](double x,double v){ return ( SQ( x-x0 ) + SQ(v-v0) < R*R ) ? 1. : 0.; };
}

auto
test2 ( double vmin , double vmax , double kx ) {
  return [=](double x,double v){ return ( vmin < v && v < vmax ) ? std::cos(2.*math::pi<double>()*kx*x) : 0.; };
}

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
  c.name = "cm_e10m3";

  c.create_output_directory();
  std::ofstream ofconfig( c.output_dir / "config.init" );
  ofconfig << c << "\n";
  ofconfig.close();

/* ------------------------------------------------------------------------- */
  field<double,1> f(boost::extents[c.Nv][c.Nx]);
  field<double,1> f_sol(boost::extents[c.Nv][c.Nx]);

  complex_field<double,1> hf(boost::extents[c.Nv][c.Nx]);
  fft::spectrum_ d(c.Nx);

  std::vector<double> mass;

  const double Kx = 0.5;
  f.range.v_min = -3.; f.range.v_max = 3.;
  f.range.x_min = -3.; f.range.x_max = 3.;
  f.compute_steps();

  f_sol.range = f.range;
  f_sol.compute_steps();

  ublas::vector<double> v (c.Nv,1.); // velocity in x direction, transport at speed 1

  ublas::vector<double> kx(c.Nx);
  {
    double l = f.range.len_x();
    for ( auto i=0 ; i<c.Nx/2 ; ++i ) { kx[i]      = 2.*math::pi<double>()*i/l; }
    for ( int i=-c.Nx/2 ; i<0 ; ++i ) { kx[c.Nx+i] = 2.*math::pi<double>()*i/l; }
  }

  auto c0 = disque(0.,0.,1.);
  auto d0 = test2(-1.,1.,1./f.range.len_x());

  double m=0;
  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      f[k][i] = c0(Xi(i),Vk(k));
      f_sol[k][i] = c0(Xi(i)-c.Tf,Vk(k)-c.Tf);

      m += f[k][i]*f.step.dx*f.step.dv;
      //f[k][i] = d0(Xi(i),Vk(k));
      //f_sol[k][i] = d0(Xi(i)-c.Tf,Vk(k)-c.Tf);
    }
  }
  std::cout << m << std::endl;
  //d.fft(f[int(c.Nv/2.)].begin());
  //std::copy( d.begin() , d.end() , std::ostream_iterator<std::complex<double>>(std::cout," "));

  f.write( c.output_dir / ("init_"+c.name+".dat") );
  f_sol.write( c.output_dir / ("sol_"+c.name+".dat") );

  iteration::iteration<double> iter;
  iter.iter = 0;
  iter.current_time = 0.;
  //iter.dt = 0.050*f.step.dv;
  //iter.dt = 2.*std::sqrt(2.)*f.step.dv;
  iter.dt = 0.150*f.step.dv;

  std::vector<double> normL2; normL2.reserve(100);
  auto L2 = [&]( field<double,1>const& f )->double {
    return std::pow(
            std::accumulate( f.origin() , f.origin()+f.num_elements() , 0. ,
              []( double s , double fik )->double { return s + fik*fik; }
          ),2)*f.step.dx*f.step.dv;
  };

  // space scheme
  auto wenol = [&](field<double,1>const& f , ublas::vector<double> const& E )->field<double,1> { return wenolin::trp_v(f,E); };
  auto weno  = [&](field<double,1>const& f , ublas::vector<double> const& E )->field<double,1> { return weno::trp_v(f,E); };
  auto cd2   = [&](field<double,1>const& f , ublas::vector<double> const& E )->field<double,1> { return o2::trp_v(f,E); };

  // time scheme init
  //expRK::HochbruckOstermann<surimi> rk(c.Nx,c.Nv,f.range.len_x(),f.shape(),v,kx,cd2); // 0.250 , 0.501 , 1.702
  expRK::CoxMatthews<surimi> rk(c.Nx,c.Nv,f.range.len_x(),f.shape(),v,kx,cd2); // 0.150 , 0.450 , 1.351
  //expRK::RK22<surimi> rk(c.Nx,c.Nv,f.range.len_x(),f.shape(),v,kx,cd2);
  //lawson::RK33<surimi> rk(c.Nx,c.Nv,f.range.len_x(),f.shape(),v,kx,weno);

  //while (  iter.current_time < c.Tf ) {
  while ( iter.iter < 101 ) {
    std::cout << "\r" << iteration::time(iter) << std::flush;

    normL2.push_back(L2(f));
    f = rk(f,iter.dt);


    ++iter.iter;
    //if ( iter.current_time + iter.dt > c.Tf ) { iter.dt = c.Tf-iter.current_time; }
    iter.current_time += iter.dt;
  }
  std::cout << "\r" << time(iter) << std::endl;

  f.write( c.output_dir / ("vp_"+c.name+".dat") );

  /*
  field<double,1> f_diff(boost::extents[c.Nv][c.Nx]);
  f_diff.range = f.range;
  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      f_diff[k][i] = f[k][i]-d0(Xi(i)-c.Tf,Vk(k)-c.Tf);
    }
  }
  f_diff.write( c.output_dir / ("diff_"+c.name+".dat") );
  */

  std::ofstream of( c.output_dir / ("normL2"+c.name+".dat") );
  std::transform(
    normL2.begin() , normL2.end() ,
    std::ostream_iterator<std::string>(of,"\n") ,
    [&,count=0](double x) mutable {
      std::stringstream ss; ss << count << " " << iter.dt*(count++) << " " << x;
      return ss.str();
    }
  );
  of.close();

  return 0;
}
