#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "miMaS/field.h"
#include "miMaS/weno.h"
#include "miMaS/fft.h"
#include "miMaS/array_view.h"
#include "miMaS/poisson.h"

namespace math = boost::math::constants;
const std::complex<double> & I = std::complex<double>(0.,1.);

#define SQ(X) ((X)*(X))
#define Xi(i) (i*f.step.dx+f.range.x_min)
#define Vk(k) (k*f.step.dv+f.range.v_min)

int main(int,char**)
{
	std::size_t Nx = 81, Nv = 256 , Nb_iter=10;
	field<double,1> f(boost::extents[Nv][Nx]);

	f.range.v_min = -12.; f.range.v_max = 12;
	f.step.dv = (f.range.v_max-f.range.v_min)/Nv;
	f.range.x_min = 0.; f.range.x_max = 20.*math::pi<double>();
  //f.range.x_min = -8.; f.range.x_max = 8.;
	f.step.dx = (f.range.x_max-f.range.x_min)/Nx;

	const double dt = 0.05; //1.606*f.step.dv/0.6; //0.5*6.*math::pi<double>()/(Nv*f.range.v_max);

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
  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      //f[k][i] = ( std::exp(-0.5*SQ(Vk(k)))*np/std::sqrt(2.*math::pi<double>()) + nb/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)/0.25) )*(1.+0.04*std::cos(0.5*Xi(i)));
      f[k][i] = ( std::exp(-0.5*SQ(Vk(k)))*np/std::sqrt(2.*math::pi<double>()) + nb/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)/0.25) )*(1.+0.04*std::cos(0.3*Xi(i)));
      //f[k][i] = std::exp( -SQ(Xi(i)-6) );
      //f[k][i] = std::exp(-SQ(Xi(i)-3)/0.5 - SQ(Vk(k))/2.);
      //f_sol[k][i] = std::exp(-SQ(Xi(i)-3)/0.5 - SQ(Vk(k))/2.);
      //f[k][i] = SQ(v(k))*std::exp(-0.5*SQ(Vk(k)))/std::sqrt(2.*math::pi<double>())*(1.+0.01*std::cos(0.5*Xi(i)));
      //f_ini[k][i] = f[k][i];
      //f_sol[k][i] =  std::exp(-SQ(Xi(i)-3.-v[k]*Nb_iter*dt)/0.5 - SQ(Vk(k)-E[i]*Nb_iter*dt)/2.);
    }
  }
  f.write("init.dat");

  poisson<double> poisson_solver(Nx,l);
  rho = f.density();
  E = poisson_solver(rho);

  //double Tf = 60.;//2*math::pi<double>();
  double Tf = 0.5;
  int i_t=0;

  std::cout << "dt " << dt << std::endl;
  std::cout << "dx " << f.step.dx << std::endl;
  std::cout << "dv " << f.step.dv << std::endl;
  std::cout << "Tf " << Tf << std::endl;

  ublas::vector<double> ee(int(std::ceil(Tf/dt)),0.);
  ublas::vector<double> Emax(int(std::ceil(Tf/dt)),0.);
  ublas::vector<double> H(int(std::ceil(Tf/dt)),0.);
  //ublas::vector<double> ee(Nb_iter);
  //ublas::vector<double> Emax(Nb_iter);
  //ublas::vector<double> H(Nb_iter);

  field<double,1> f1(tools::array_view<const std::size_t>(f.shape(),2)),f2(tools::array_view<const std::size_t>(f.shape(),2)),f3(tools::array_view<const std::size_t>(f.shape(),2));
  fft::spectrum_ hf(Nx),hf1(Nx),hf2(Nx),hf3(Nx),hEdvf(Nx);

  while (  i_t*dt < Tf ) {
  //while (  i_t < Nb_iter ) {
  	//if (t%32==0) { std::cout<<"\r"<<t<<" "<<std::flush ; }
    std::cout<<" \r"<<i_t<<" "<<std::flush;
    Emax(i_t) = std::abs(*std::max_element( E.begin() , E.end() , [](double a,double b){return (std::abs(a) < std::abs(b));} ));


/**/
    rho = f.density();
    E = poisson_solver(rho);
    field<double,1> Edvf = weno::trp_v(f,E);
	  for ( auto k=0 ; k<f.size(0) ; ++k ) {
	  	hf.fft(&(f[k][0]));
	  	hEdvf.fft(&(Edvf[k][0]));

	  	for ( auto i=0 ; i<Nx ; ++i ) {
        hf1[i] = std::exp(-v(k)*I*kx[i]*dt)*( hf[i]-dt*hEdvf[i] );
	  	}
	  	hf1.ifft(&(f1[k][0]));
	  }

    rho = f1.density();
    E = poisson_solver(rho);
    Edvf = weno::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf2[i] = 0.75*std::exp(-0.5*v(k)*I*kx[i]*dt)*hf[i] + 0.25*std::exp(0.5*v(k)*I*kx[i]*dt)*( hf1[i]-dt*hEdvf[i] );
      }
      hf2.ifft(&(f2[k][0]));
    }

    rho = f2.density();
    E = poisson_solver(rho);
    Edvf = weno::trp_v(f2,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf2.fft(&(f2[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf[i] = (1./3.)*std::exp(-v(k)*I*kx[i]*dt)*hf[i] + (2./3.)*std::exp(-0.5*v(k)*I*kx[i]*dt)*( hf2[i]-dt*hEdvf[i] );
      }
      hf.ifft(&(f[k][0]));
    }

/**/

/*
    E = poisson_solver( f.density() );
    field<double,1> Edvf = weno::trp_v(f,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf1[i] = std::exp(-0.5*v(k)*I*kx[i]*dt)*( hf[i] - 0.5*dt*hEdvf[i] );
      }

      hf1.ifft(&(f1[k][0]));
    }

    E = poisson_solver( f1.density() );
    Edvf = weno::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf2[i] = std::exp(-0.5*v(k)*I*kx[i]*dt)*hf[i] - 0.5*dt*hEdvf[i] ;
      }
      
      hf2.ifft(&(f2[k][0]));
    }

    E = poisson_solver( f2.density() );
    Edvf = weno::trp_v(f2,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf3[i] = std::exp(-v(k)*I*kx[i]*dt)*hf[i] - dt*std::exp(-0.5*v(k)*I*kx[i]*dt)*hEdvf[i] ;
      }

      hf3.ifft(&(f3[k][0]));
    }

    E = poisson_solver( f3.density() );
    Edvf = weno::trp_v(f3,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hf2.fft(&(f2[k][0]));
      hf3.fft(&(f3[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        hf[i] = -(1./3.)*std::exp(-v(k)*I*kx[i]*dt)*hf[i] + (1./3.)*std::exp(-0.5*v(k)*I*kx[i]*dt)*hf1[i] + (2./3.)*std::exp(-0.5*v(k)*I*kx[i]*dt)*hf2[i] + (1./3.)*hf3[i] + (1./6.)*dt*hEdvf[i]  ;
      }

      hf.ifft(&(f[k][0]));
    }

*/

    rho = f.density();
    E = poisson_solver(rho);
    ee(i_t) = 0.;
    for ( auto i=0 ; i<Nx ; ++i ) {
      ee(i_t) += SQ(E(i))*f.step.dx;
    }
    H(i_t) = energy(f,E);
    ++i_t;
	}
  std::cout<<" \r"<<i_t<<std::endl;


  f.write("vp.dat");
  std::ofstream of;
  std::size_t count = 0;
  auto dt_y = [&,count=0](auto const& y) mutable { std::stringstream ss; ss<<count*dt<<" "<<y; return ss.str(); };
  of.open("ee.dat");
  for ( auto i=0; i<ee.size() ; ++i ) {
    of << i*dt <<" " << ee[i] << "\n";
  }
  of.close();
  of.open("H.dat");
  for ( auto i=0; i<H.size() ; ++i ) {
    of << i*dt <<" " << (H[i]-H[0])/std::abs(H[0]) << "\n";
  }

  double h = std::abs(*std::max_element( H.begin() , H.end() , [&](double a,double b){return ( std::abs((a-H[0])/std::abs(H[0])) < std::abs((b-H[0])/std::abs(H[0])) );} ));
  std::cout << dt << " " << std::abs((h-H[0])/std::abs(H[0])) << "\n";

  of.close();
  of.open("Emax.dat");
  std::transform( ee.begin() , ee.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();

	return 0;
}
