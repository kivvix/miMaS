#include <iostream>
#include <algorithm>
#include <iterator>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "miMaS/field.h"
#include "miMaS/weno.h"
#include "miMaS/fft.h"
#include "miMaS/array_view.h"

namespace math = boost::math::constants;

int main(int,char**)
{
	std::size_t Nx = 33, Nv = 64;
	double dt = 0.01;
	field<double,1> f(boost::extents[Nv][Nx]);
	ublas::vector<double> v (Nv); std::generate(v.begin(),v.end(),[](){return -1.;});
	ublas::vector<double> E (Nx); std::generate(E.begin(),E.end(),[](){return 1.;});

	f.range.v_min = -18.; f.range.v_max = 18.;
	f.step.dv = (f.range.v_max-f.range.v_min)/Nv;
	f.range.x_min = 0.; f.range.x_max = 1.;
	f.step.dx = (f.range.x_max-f.range.x_min)/Nx;
	
	double l = f.range.x_max-f.range.x_min;
	
	ublas::vector<double> kx(Nx);
	std::generate(kx.begin(),kx.end(),[](){return 42.;});
	for ( auto i=0 ; i<Nx/2+1 ; ++i )   { kx[i] = i/l; }
	for ( auto i=0 ; i<((Nx/2)) ; ++i ) { kx[i+Nx/2+1] = -kx[Nx/2-i]; }
	std::copy(kx.begin(),kx.end(),std::ostream_iterator<double>(std::cout," ")); std::cout << std::endl;
	
  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      f[k][i] = cos(2.*math::pi<double>()*(f.step.dx*i+f.range.x_min))*cos((7*math::pi<double>())*(f.step.dv*k+f.range.v_min));
    }
  }
  f.write("init.dat");

  //fft::fft_x<double,1> fft_trpx( tools::array_view<const std::size_t>(f.shape(),2),Nx,f.range.x_max-f.range.x_min );

  for ( auto t=0 ; t<10 ; ++t ) {
  	std::cout << t << "\r";
	  field<double,1> Edvf = weno::trp_v(f,E);
	  fft::spectrum hf(Nx),hEdvf(Nx);
	  for ( auto k=0 ; k<f.size(0) ; ++k ) {
	  	hf.fft(&(f[k][0]));
	  	hEdvf.fft(&(Edvf[k][0]));

	  	for ( auto i=0 ; i<Nx ; ++i ) {
	  		auto re = hf[i][fft::re];
	  		auto im = hf[i][fft::im];
	  		hf[i][fft::re] = std::cos(v[k]*kx[i]*dt)*(re-dt*hEdvf[i][fft::re]) + std::sin(v[k]*kx[i]*dt)*(im-dt*hEdvf[i][fft::im]);
	  		hf[i][fft::im] = std::cos(v[k]*kx[i]*dt)*(im-dt*hEdvf[i][fft::im]) - std::sin(v[k]*kx[i]*dt)*(re-dt*hEdvf[i][fft::re]);
	  	}
	  	hf.ifft(&(f[k][0]));
	  }
	}

  f.write("vp.dat");

	return 0;
}
