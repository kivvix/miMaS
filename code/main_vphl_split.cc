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

int main(int,char**)
{
  std::size_t Nx = 135, Nv = 256;

  // $(u_c,E,\hat{f}_h)$ and $f_h$
  ublas::vector<double> uc(Nx,0.);
  ublas::vector<double> E (Nx,0.);
  field<double,1> fh(boost::extents[Nv][Nx]);
  complex_field<double,1> hfh(boost::extents[Nv][Nx]);

  const double Kx = 0.5;
  // phase-space domain
  fh.range.v_min = -8.; fh.range.v_max = 8.;
  fh.range.x_min =  0.; fh.range.x_max = 2./Kx*math::pi<double>();

  // compute dx, dv
  fh.step.dv = (fh.range.v_max-fh.range.v_min)/Nv;
  fh.step.dx = (fh.range.x_max-fh.range.x_min)/Nx;

  double dt = 0.05;//1.*fh.step.dv;
  double Tf = 200.;
  
  // velocity and frequency
  ublas::vector<double> v (Nv,0.); for ( std::size_t k=0 ; k<Nv ; ++k ) { v[k] = Vk(k); }
  const double l = fh.range.x_max-fh.range.x_min;
  ublas::vector<double> kx(Nx);
  for ( auto i=0 ; i<Nx/2 ; ++i ) { kx[i]    = 2.*math::pi<double>()*i/l; }
  for ( int i=-Nx/2 ; i<0 ; ++i ) { kx[Nx+i] = 2.*math::pi<double>()*i/l; }

  // initial condition
  double ui=2., alpha=0.2;
  auto tb_M1 = maxwellian(0.5*alpha,ui,1.) , tb_M2 = maxwellian(0.5*alpha,-ui,1.);
  for (field<double,2>::size_type k=0 ; k<fh.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<fh.size(1) ; ++i ) {
      //fh[k][i] = ( 0.5*alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)) + 0.5*alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)+ui)) )*(1.+0.04*std::cos(0.3*Xi(i)));
      //fh[k][i] = ( 0.5*alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)) + 0.5*alpha/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)+ui)) )*(1.+0.04*std::cos(Kx*Xi(i)));

      fh[k][i] = ( tb_M1(Xi(i),Vk(k)) + tb_M2(Xi(i),Vk(k)) )*(1. + 0.01*std::cos(Kx*Xi(i)));
    }
    fft::fft(&(fh[k][0]),&(fh[k][Nx-1])+1,&(hfh[k][0]));
  }
  fh.write("vphl/split/init.dat");

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
  // init E (electric field) with Poisson
  {
    poisson<double> poisson_solver(Nx,l);
    ublas::vector<double> rho(Nx,0.);
    rho = fh.density(); // compute density from init data
    for ( auto i=0 ; i<Nx ; ++i ) { rho[i] += rho_c; } // add (1-alpha) for cold particules
    E = poisson_solver(rho);
  }

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

  //U_type<double,1> U({uc,E,hfh};
  splitting<double,1> Lie( fh , l , rho_c );

  std::size_t i_t = 0;
  double current_time = 0.;
  while ( current_time < Tf ) {
    std::cout << " [" << std::setw(5) << i_t << "] " << i_t*dt << "\r" << std::flush;

    Lie.phi_a(dt,uc,E,hfh);
    Lie.phi_b(dt,uc,E,hfh);
    Lie.phi_c(dt,uc,E,hfh);
    //Lie.phi_b(0.5*dt,uc,E,hfh);
    //Lie.phi_a(0.5*dt,uc,E,hfh);

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

  of.open("vphl/split/ee.dat");
  std::transform( ee.begin() , ee.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();

  of.open("vphl/split/Emax.dat");
  std::transform( Emax.begin() , Emax.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();

  for ( auto k=0 ; k<hfh.shape()[0] ; ++k ) { fft::ifft(&(hfh[k][0]),&(hfh[k][Nx-1])+1,&(fh[k][0])); }
  fh.write("vphl/split/vp.dat");


  return 0;
}

