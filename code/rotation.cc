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

struct err
{
  double infty=0.;
  double one=0.;
  double dt=0. , dx=0. , dv=0.;
  std::size_t n=16;
  std::size_t nb_iter=0;
};

std::ostream &
operator << ( std::ostream & os , err const& e ) {
  os << e.n << " " << " " << e.nb_iter << " " << e.dx << " " << e.dt << " " << e.one << " " << e.infty;
  return os;
}

err
rotation ( std::size_t N , int q=-2 )
{
	std::size_t Nx = N, Nv = N;
	field<double,1> f(boost::extents[Nv][Nx]);
  field<double,1> f_sol(boost::extents[Nv][Nx]);

	f.range.v_min = -10.; f.range.v_max = 10;
	f.step.dv = (f.range.v_max-f.range.v_min)/Nv;
	f.range.x_min = -10.; f.range.x_max = 10.;
	f.step.dx = (f.range.x_max-f.range.x_min)/Nx;


	const double dt = (math::pi<double>()*std::pow(2,q))*f.step.dv/f.range.v_max;
	

	ublas::vector<double> v (Nv,0.); for ( std::size_t k=0 ; k<Nv ; ++k ) { v[k] =  Vk(k); }
  ublas::vector<double> E (Nx,0.); for ( std::size_t i=0 ; i<Nx ; ++i ) { E[i] = -Xi(i); }
  
	const double l = f.range.x_max-f.range.x_min;
	ublas::vector<double> kx(Nx);
  for ( auto i=0 ; i<Nx/2 ; ++i ) { kx[i] = 2.*math::pi<double>()*i/l; }
  for ( int i=-Nx/2 ; i<0 ; ++i ) { kx[Nx+i] = 2.*math::pi<double>()*i/l; }
	
  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      f[k][i] = std::exp(-SQ(Xi(i)-3)/0.5 - SQ(Vk(k))/2.);
      f_sol[k][i] = f[k][i];
    }
  }

  double Tf = 2*math::pi<double>();
  int i_t=0;

  while ( i_t*dt < Tf ) {
    field<double,1> Edvf = weno::trp_v(f,E);
    field<double,1> f1=f,f2=f;
	  fft::spectrum hf(Nx),hf1(Nx),hf2(Nx),hEdvf(Nx);

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

    ++i_t;
	}

  err e;
  e.dt = dt;
  e.dx = f.step.dx;
  e.dv = f.step.dv;
  e.n = N;
  e.nb_iter = i_t;
  for ( auto k=0 ; k<Nv ; ++k ) {
    for (auto i=0 ; i<Nx ; ++i ) {
      e.one += std::abs(f[k][i]-f_sol[k][i])*f.step.dv*f.step.dx;
      e.infty = std::max(std::abs(f[k][i]-f_sol[k][i]),e.infty);
    }
  }

	return e;
}

int
main ( int argc , char** argv )
{
  for ( std::size_t exp=5 ; exp<11 ; ++exp ) {
    std::size_t n=1<<exp;
    err e = rotation(n);
    std::cout << e << std::endl;
  }
  
  return 0;
}

