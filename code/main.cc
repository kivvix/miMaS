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

#define SQ(X) ((X)*(X))
#define Xi(i) (i*f.step.dx+f.range.x_min)
#define Vk(k) (k*f.step.dv+f.range.v_min)

int main(int,char**)
{
	std::size_t Nx = 65, Nv = 128 , Nb_iter=10;
	field<double,1> f(boost::extents[Nv][Nx]);

	f.range.v_min = -8.; f.range.v_max = 8;
	f.step.dv = (f.range.v_max-f.range.v_min)/Nv;
	f.range.x_min = 0.; f.range.x_max = 4.*math::pi<double>();
	f.step.dx = (f.range.x_max-f.range.x_min)/Nx;

	const double dt = 0.5*f.step.dv;

  field<double,1> f_sol = f;
  field<double,1> f_ini = f;
	
	ublas::vector<double> v (Nv,0.);
  ublas::vector<double> E (Nx,1.),rho(Nx);
  //for ( std::size_t k=0 ; k<Nv ; ++k ) { v[k] = Vk(k); }
  for ( std::size_t k=0 ; k<Nv ; ++k ) { v[k] = Vk(k); }

	const double l = f.range.x_max-f.range.x_min;
	ublas::vector<double> kx(Nx);
	for ( auto i=0 ; i<Nx/2+1 ; ++i )   { kx[i] = 2.*math::pi<double>()*i/l; }
	for ( auto i=0 ; i<((Nx/2)) ; ++i ) { kx[i+Nx/2+1] = -kx[Nx/2-i]; }
	
  double np = 0.9 , nb = 0.2 , ui = 4.5;
  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      //f[k][i] = ( std::exp(-0.5*SQ(Vk(k)))*np/std::sqrt(2.*math::pi<double>()) + nb/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)/0.25) )*(1.+0.04*std::cos(0.5*Xi(i)));
      f[k][i] = std::exp( -SQ(Xi(i)-6) );
      //f[k][i] = std::exp(-SQ(Xi(i)-3)/0.5 - SQ(Vk(k))/2.);
      //f_ini[k][i] = f[k][i];
      //f_sol[k][i] =  std::exp(-SQ(X(i)-0.9-v[k]*200*dt)/0.01 - SQ(V(k)-E[i]*200*dt)/0.05);
    }
  }
  f.write("init.dat");

  poisson<double> poisson_solver(Nx,l);
  //rho = f.density();
  //E = poisson_solver(rho);
  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      //f_sol[k][i] = ( std::exp(-0.5*SQ(Vk(k)-E(i)*Nb_iter*dt))*np/std::sqrt(2.*math::pi<double>()) + nb/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui-E(i)*Nb_iter*dt)/0.25) )*(1.+0.04*std::cos(0.5*(Xi(i)-v(k)*Nb_iter*dt)));
      f_sol[k][i] = std::exp( -SQ(Xi(i)-6-v(k)*Nb_iter*dt) );
      //f_sol[k][i] = std::exp(-SQ(Xi(i)-3-v(k)*Nb_iter*dt)/0.5 - SQ(Vk(k)-E(i)*Nb_iter*dt)/2.);
    }
  }

int BUG=0;
  for ( auto t=0 ; t<Nb_iter ; ++t ) {
  	//if (t%32==0) { std::cout<<"\r"<<t<<" "<<std::flush ; }
    //std::cout<<"\r"<<t<<" "<<std::flush;

std::cerr << "bug " << BUG++ << std::endl;
    field<double,1> Edvf = weno::trp_v(f,E);
std::cerr << "bug " << BUG++ << std::endl;
    field<double,1> f1=f,f2=f;
std::cerr << "bug " << BUG++ << std::endl;
	  fft::spectrum hf(Nx),hf1(Nx),hf2(Nx),hEdvf(Nx);
std::cerr << "bug " << BUG++ << std::endl;
    for ( auto i=0 ; i<Nx ; ++i ) {
std::cerr << i << std::endl;
      hf[i][fft::re] = hf[i][fft::im] = hf1[i][fft::re] = hf1[i][fft::im] = hf2[i][fft::re] = hf2[i][fft::im] = hEdvf[i][fft::re] = hEdvf[i][fft::im] = 0.;
    }
std::cerr << "bug " << BUG++ << std::endl;

	  for ( auto k=0 ; k<f.size(0) ; ++k ) {
	  	hf.fft(&(f[k][0]));
	  	hEdvf.fft(&(Edvf[k][0]));

	  	for ( auto i=0 ; i<Nx ; ++i ) {
	  		auto re = hf[i][fft::re] , im = hf[i][fft::im];
	  		hf1[i][fft::re] = std::cos(v(k)*kx[i]*dt)*(re-dt*hEdvf[i][fft::re]) + std::sin(v(k)*kx[i]*dt)*(im-dt*hEdvf[i][fft::im]);
	  		hf1[i][fft::im] = std::cos(v(k)*kx[i]*dt)*(im-dt*hEdvf[i][fft::im]) - std::sin(v(k)*kx[i]*dt)*(re-dt*hEdvf[i][fft::re]);
	  	}
	  	hf1.ifft(&(f1[k][0]));
	  }

    //rho = f1.density();
    //E = poisson_solver(rho);
    Edvf = weno::trp_v(f1,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf1.fft(&(f1[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        auto re = hf[i][fft::re] , im = hf[i][fft::im];
        hf2[i][fft::re] = 0.75*(std::cos(0.5*v(k)*kx[i]*dt)*re + std::sin(0.5*v(k)*kx[i]*dt)*im) + 0.25*( std::cos(0.5*v(k)*kx[i]*dt)*(hf1[i][fft::re]-dt*hEdvf[i][fft::re]) - std::sin(0.5*v(k)*kx[i]*dt)*(hf1[i][fft::im]-dt*hEdvf[i][fft::im]) );
        hf2[i][fft::im] = 0.75*(std::cos(0.5*v(k)*kx[i]*dt)*im - std::sin(0.5*v(k)*kx[i]*dt)*re) + 0.25*( std::cos(0.5*v(k)*kx[i]*dt)*(hf1[i][fft::im]-dt*hEdvf[i][fft::im]) + std::sin(0.5*v(k)*kx[i]*dt)*(hf1[i][fft::re]-dt*hEdvf[i][fft::re]) );
      }
      hf2.ifft(&(f2[k][0]));
    }

    //rho = f2.density();
    //E = poisson_solver(rho);
    Edvf = weno::trp_v(f2,E);
    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      hf.fft(&(f[k][0]));
      hf2.fft(&(f2[k][0]));
      hEdvf.fft(&(Edvf[k][0]));

      for ( auto i=0 ; i<Nx ; ++i ) {
        auto re = hf[i][fft::re] , im = hf[i][fft::im];
        hf[i][fft::re] = (1./3.)*(std::cos(v(k)*kx[i]*dt)*re + std::sin(v(k)*kx[i]*dt)*im) + (2./3.)*( std::cos(0.5*v(k)*kx[i]*dt)*(hf2[i][fft::re]-dt*hEdvf[i][fft::re]) + std::sin(0.5*v(k)*kx[i]*dt)*(hf2[i][fft::im]-dt*hEdvf[i][fft::im]) );
        hf[i][fft::im] = (1./3.)*(std::cos(v(k)*kx[i]*dt)*im - std::sin(v(k)*kx[i]*dt)*re) + (2./3.)*( std::cos(0.5*v(k)*kx[i]*dt)*(hf2[i][fft::im]-dt*hEdvf[i][fft::im]) - std::sin(0.5*v(k)*kx[i]*dt)*(hf2[i][fft::re]-dt*hEdvf[i][fft::re]) );
      }
      hf.ifft(&(f[k][0]));
    }

    //rho = f.density();
    //E = poisson_solver(rho);
    //std::copy(E.begin(),E.end(),std::ostream_iterator<double>(std::cout," ")); std::cout << std::endl;
    //std::cout << *std::max_element(E.begin(),E.end()) << std::endl;
	}


  f.write("vp.dat");
  f_sol.write("sol.dat");
  std::ofstream of("diff.dat");
  for ( auto k=0 ; k<Nv ; ++k ) {
    for (auto i=0 ; i<Nx ; ++i ) {
      of << Xi(i) << " " << Vk(k) << " " << f[k][i]-f_sol[k][i] << "\n";
    }
    of << std::endl;
  }
  of.close();
/*
  std::ofstream of("diff.dat");
  for ( auto k=0 ; k<Nv ; ++k ) {
    for (auto i=0 ; i<Nx ; ++i ) {
      of << X(i) << " " << V(k) << " " << f[k][i] << " " <<  f_sol[k][i] << "\n";
    }
    of << std::endl;
  }

  of.close();
*/
	return 0;
}
