#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include <valarray>
#include <sstream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "miMaS/field.h"
#include "miMaS/weno.h"
#include "miMaS/fft.h"
#include "miMaS/array_view.h"
#include "miMaS/poisson.h"
#include "miMaS/rk.h"

#ifndef SIGMA
#define SIGMA (1.0)
#endif
#define E_MAX (2.0)

/*
template <unsigned int i>
struct phi
{
  static std::complex<double>
  operator () ( std::complex<double> const & z ) {
    static std::valarray<std::complex<double>> coeff(i);
    coeff[0] = 1.;

    for ( unsigned int k=1 ; k<coeff.size() ; ++k ) {
      coeff[k] = coeff[k-1] * z / (double(k));
    }

    return (std::exp(z) - std::accumulate( std::begin(coeff) , std::end(coeff) , std::complex<double>(0.,0.) ))/(std::pow(z,i));
  }
};
*/
template <unsigned int i>
std::complex<double>
phi ( std::complex<double> const & _z )
{
  std::valarray<std::complex<double>> coeff(i);
  coeff[0] = 1.;

  std::complex<double> z = _z;
  if ( _z == 0. ) { z = std::complex<double>(1.,0.); }

  for ( unsigned int k=1 ; k<coeff.size() ; ++k ) {
    coeff[k] = coeff[k-1] * z / (double(k));
  }
  //std::copy(std::begin(coeff),std::end(coeff),std::ostream_iterator<std::complex<double>>(std::cout," . "));
  //std::cout << std::endl;

  if ( z != 0. ) {
    return (std::exp(z) - std::accumulate( std::begin(coeff) , std::end(coeff) , std::complex<double>(0.,0.) ))/(std::pow(z,i));  
  }
  return coeff[i-1];
}

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
  //std::cout << rho << " " << u << " " << T << std::endl;
  //std::cout << rho/(std::sqrt(2.*math::pi<double>()*T)) << std::endl;
  return [=](double x,double v){ return rho/(std::sqrt(2.*math::pi<double>()*T))*std::exp( -0.5*SQ(v-u)/T ); };
}

int main(int,char**)
{
	std::size_t Nx = 135, Nv = 256 , Nb_iter=10;
	field<double,1> f(boost::extents[Nv][Nx]);

	f.range.v_min = -8.; f.range.v_max = 8.;
	f.step.dv = (f.range.v_max-f.range.v_min)/Nv;

  const double Kx = 0.5;
	//f.range.x_min = 0.; f.range.x_max = 20.*math::pi<double>();
  f.range.x_min = 0.; f.range.x_max = 2./Kx*math::pi<double>(); //10.0*math::pi<double>();
  //f.range.x_min = -8.; f.range.x_max = 8.;
	f.step.dx = (f.range.x_max-f.range.x_min)/Nx;

  //field<double,1> f_sol = f;
  //field<double,1> f_ini = f;
	
	ublas::vector<double> v (Nv,0.);
  ublas::vector<double> E (Nx,0.),rho(Nx);
  for ( std::size_t k=0 ; k<Nv ; ++k ) { v[k] = Vk(k); }
  //for ( std::size_t i=0 ; i<Nx ; ++i ) { E[i] = -Xi(i); }

	const double l = f.range.x_max-f.range.x_min;
	ublas::vector<double> kx(Nx);
	//for ( auto i=0 ; i<Nx/2+1 ; ++i )   { kx[i] = 2.*math::pi<double>()*i/l; }
	//for ( auto i=0 ; i<((Nx/2)) ; ++i ) { kx[i+Nx/2+1] = -kx[Nx/2-i]; }
  for ( auto i=0 ; i<Nx/2 ; ++i ) { kx[i]    = 2.*math::pi<double>()*i/l; }
  for ( int i=-Nx/2 ; i<0 ; ++i ) { kx[Nx+i] = 2.*math::pi<double>()*i/l; }
	
  double np = 0.9 , nb = 0.2 , ui = 4.5;
  double alpha = 0.1;
  double Tc = 0.01; //, u = 4.5;
  //auto Mmu = maxwellian(1.,-u,1.) , Mpu = maxwellian(1.,u,1.) , M0c = maxwellian(1.,0.,Tc);

  ui = 2.;
  alpha = 0.2;
  Tc = 0.01;
  auto landau_M = maxwellian(1.,0.,1.);
  auto db_M1 = maxwellian(0.5,-ui,1.) , db_M2 = maxwellian(0.5,ui,1.);
  auto bot_M1 = maxwellian(1.-alpha,0.,1.) , bot_M2 = maxwellian(alpha,ui,0.25);
  auto tb_MC = maxwellian(1.-alpha,0.,Tc) , tb_M1 = maxwellian(0.5*alpha,ui,1.) , tb_M2 = maxwellian(0.5*alpha,-ui,1.);
  auto v10_MC = maxwellian(1.-alpha,0,Tc) , v10_Mh = maxwellian(alpha,0.,1.);

  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      //f[k][i] = ((1.-alpha)*M0c(Xi(i),Vk(k)) + 0.5*alpha*( Mpu(Xi(i),Vk(k)) + Mmu(Xi(i),Vk(k)) ))*(1.+0.01*std::cos(0.5*Xi(i)));

      //// landau damping : Kx=0.5
      //f[k][i] = landau_M(Xi(i),Vk(k))*(1.+0.001*std::cos(Kx*Xi(i)));
      //// double beam Kx=0.2, ui=2.4 ou Kx=0.2, ui=4.5
      //f[k][i] = (db_M1(Xi(i),Vk(k))+db_M2(Xi(i),Vk(k)))*(1.+0.001*std::cos(Kx*Xi(i)));
      //// bot Kx=0.5 , alpha=0.2 , ui=4.5
      //f[k][i] = (bot_M1(Xi(i),Vk(k)) + bot_M2(Xi(i),Vk(k)))*(1.+0.01*std::cos(Kx*Xi(i)));
      //// tb Kx=0.5 , ui=4. , alpha=0.2 , Tc=0.01
      f[k][i] = tb_MC(Xi(i),Vk(k)) + (tb_M1(Xi(i),Vk(k)) + tb_M2(Xi(i),Vk(k)) )*(1. + 0.01*std::cos(Kx*Xi(i)));
      //// v10 Kx=0.5 , alpha=0.2
      //f[k][i] = v10_MC(Xi(i),Vk(k)) + ( std::pow(Vk(k),10)*v10_Mh(Xi(i),Vk(k))/945. )*(1. + 0.01*std::cos(Kx*Xi(i)));
    }
  }
  f.write("vphl/kin/init_tb.dat");

  poisson<double> poisson_solver(Nx,l);
  rho = f.density();
  E = poisson_solver(rho);
  

  //double Tf = 60.;//2*math::pi<double>();
  const double Tf = 200.;
  int i_t=0;

  double current_time = 0.;
  // SIGMA is the CFL number
  double dt = f.step.dv;//std::min( 0.1 , SIGMA*f.step.dv/std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );

  std::cout << "Nx: " << Nx << "\n";
  std::cout << "Nv: " << Nv << "\n";
  std::cout << "v_min: " << f.range.v_min << "\n";
  std::cout << "v_max: " << f.range.v_max << "\n";
  std::cout << "x_min: " << f.range.x_min << "\n";
  std::cout << "x_max: " << f.range.x_max << "\n";
  std::cout << "dt: " << dt << "\n";
  std::cout << "dx: " << f.step.dx << "\n";
  std::cout << "dv: " << f.step.dv << "\n";
  std::cout << "Tf: " << Tf << "\n";
  std::cout << "f_0: " << "\"tb\"" << "\n";
  std::cout << std::endl;

  std::ofstream info("info.yaml");

  info << "Nx: " << Nx << "\n";
  info << "Nv: " << Nv << "\n";
  info << "v_min: " << f.range.v_min << "\n";
  info << "v_max: " << f.range.v_max << "\n";
  info << "x_min: " << f.range.x_min << "\n";
  info << "x_max: " << f.range.x_max << "\n";
  info << "dt: " << dt << "\n";
  info << "dx: " << f.step.dx << "\n";
  info << "dv: " << f.step.dv << "\n";
  info << "Tf: " << Tf << "\n";
  info << "f_0: " << "\"tb\"" << "\n";
  info << std::endl;
  info.close();

  std::vector<double> ee;   ee.reserve(int(std::ceil(Tf/dt))+1);
  std::vector<double> Emax; Emax.reserve(int(std::ceil(Tf/dt))+1);
  std::vector<double> H;    H.reserve(int(std::ceil(Tf/dt))+1);

  std::vector<double> times; times.reserve(int(std::ceil(Tf/dt))+1);

  // space scheme
  auto weno = [&](field<double,1>const& f , ublas::vector<double> const& E )->field<double,1> { return weno::trp_v(f,E); };
  auto cd2 = [&](field<double,1>const& f , ublas::vector<double> const& E )->field<double,1> { return o2::trp_v(f,E); };

  // time scheme init
  lawson::RK33<poisson<double>> rk(Nx,Nv,l,f.shape(),v,kx,weno);

  rk.E = rk.poisson_solver(f.density());
  Emax.push_back( std::abs(*std::max_element( rk.E.begin() , rk.E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );
  double electric_energy = 0.;
  for ( const auto & ei : rk.E ) { electric_energy += ei*ei*f.step.dx; }
  ee.push_back( std::sqrt(electric_energy) );
  H.push_back( energy(f,rk.E) );
  times.push_back( 0. );

  while (  current_time < Tf ) {
    std::cout<<" ["<<std::setw(5)<<i_t<<"] "<< current_time <<"\r"<<std::flush;
    
    f = rk(f,dt);

    // end of time loop

    // MONITORING
    rk.E = rk.poisson_solver(f.density());
    Emax.push_back( std::abs(*std::max_element( rk.E.begin() , rk.E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} )) );
    electric_energy = 0.;
    for ( const auto & ei : rk.E ) { electric_energy += ei*ei*f.step.dx; }
    ee.push_back( std::sqrt(electric_energy) );
    H.push_back( energy(f,rk.E) );

    dt = std::min( 0.1 , SIGMA*f.step.dv/Emax[i_t] );

    // increment time
    ++i_t;
    current_time += dt;
    times.push_back( current_time );
	} // while (  i_t*dt < Tf )


//#define FOLDER "lukas/landau/"
#define FOLDER "vphl/kin/"
#define SPACE_SCHEME "weno"

  std::cout<<" ["<<std::setw(5)<<i_t<<"] "<<i_t*dt <<std::endl;

  std::stringstream ss; ss << FOLDER << "vp_" << rk.label << "_" << SPACE_SCHEME << ".dat";
  f.write("vphl/kin/vp_tb.dat");

  ss.str(std::string());
  std::ofstream of;
  std::size_t count = 0;
  auto dt_y = [&,count=0](auto const& y) mutable { std::stringstream ss; ss<<times[count++]<<" "<<y; return ss.str(); };
  //ss << FOLDER << "ee_" << rk.label << "_" << SPACE_SCHEME << ".dat";
  of.open("vphl/kin/ee_tb.dat"); ss.str(std::string());
  for ( auto i=0; i<ee.size() ; ++i ) {
    of << times[i] <<" " << ee[i] << "\n";
  }
  of.close();
  //ss << FOLDER << "H_" << rk.label << "_" << SPACE_SCHEME << ".dat";
  of.open("vphl/kin/H_tb.dat"); ss.str(std::string());
  for ( auto i=0; i<H.size() ; ++i ) {
    of << times[i] <<" " << (H[i]-H[0])/std::abs(H[0]) << "\n";
  }

  //double h = std::abs(*std::max_element( H.begin() , H.end() , [&](double a,double b){return ( std::abs((a-H[0])/std::abs(H[0])) < std::abs((b-H[0])/std::abs(H[0])) );} ));
  //std::cout << dt << " " << std::abs((h-H[0])/std::abs(H[0])) << "\n";

  of.close();
  //ss << FOLDER << "Emax_" << rk.label << "_" << SPACE_SCHEME << ".dat";
  of.open("vphl/kin/Emax_tb.dat"); ss.str(std::string());
  std::transform( Emax.begin() , Emax.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();

	return 0;
}
