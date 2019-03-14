#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include <valarray>
#include <complex>

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
#define Xi(i) (i*f0.step.dx+f0.range.x_min)
#define Vk(k) (k*f0.step.dv+f0.range.v_min)

int main(int,char**)
{
  const std::size_t NumDimV = 1;
	const int Nx = 64, Nv = 128, Nb_iter=100;

	field<double,NumDimV> f0( boost::extents[Nv][Nx] );

	f0.range.v_min = -10.; f0.range.v_max = 10.;
	f0.step.dv = (f0.range.v_max-f0.range.v_min)/Nv;
	f0.range.x_min = -10.; f0.range.x_max = 10.;
	f0.step.dx = (f0.range.x_max-f0.range.x_min)/Nx;

  double Tf = 2.*math::pi<double>();

	ublas::vector<double> v (Nv);
  ublas::vector<double> E (Nx);
  for ( std::size_t k=0 ; k<Nv ; ++k ) { v[k] = Vk(k); }
  for ( std::size_t i=0 ; i<Nx ; ++i ) { E[i] = -Xi(i); }

	const double lx = f0.range.x_max-f0.range.x_min;
  const double lv = f0.range.v_max-f0.range.v_min;
	ublas::vector<double> kx(Nx),kv(Nv);
  for ( int i=0  ; i<Nx/2 ; ++i ) { kx[i]    = 2.*math::pi<double>()*i/lx; }
  for ( int i=-Nx/2 ; i<0 ; ++i ) { kx[Nx+i] = 2.*math::pi<double>()*i/lx; }

  for ( int k=0  ; k<Nv/2 ; ++k ) { kv[k]    = 2.*math::pi<double>()*k/lv; }
  for ( int k=-Nv/2 ; k<0 ; ++k ) { kv[Nv+k] = 2.*math::pi<double>()*k/lv; }
	
  for (field<double,NumDimV>::size_type k=0 ; k<f0.size(0) ; ++k ) {
    for (field<double,NumDimV>::size_type i=0 ; i<f0.size(1) ; ++i ) {
      //f[k][i] = std::cos(Xi(i)*0.2)*std::cos(math::pi<double>()*2.*Vk(k)/20);
      f0[k][i] = std::exp( -SQ(Xi(i)-1.)/2. )*std::exp( -Vk(k)*Vk(k)/1.);
      //f[k][i] = ( std::exp(-0.5*SQ(Vk(k)))*np/std::sqrt(2.*math::pi<double>()) )*(1.+0.04*std::cos(0.3*Xi(i)));
    }
  }
  f0.write("init.dat");


  fft::spectrum_ hxf(Nx);
  fft::spectrum_ hvf(Nv);


  for ( int nb_iter=10 ; nb_iter<100 ; nb_iter+=10 ){
    double dt = Tf/nb_iter;

    int i_t=0;
    field<double,1> f=f0;
    while ( i_t*dt < Tf ) {
      //std::cout<<" \r"<<i_t<<" "<<std::flush;

      field<double,1> f1 = f;
      field<double,1> f2 = f;
  	  
  	  for ( std::size_t k=0 ; k<Nv ; ++k ) {
  	  	hxf.fft(&(f[k][0]));
  	  	for ( std::size_t i=0 ; i<Nx ; ++i ) { hxf[i] = std::exp( -I*v(k)*kx[i]*0.5*dt )*hxf[i]; }
  	  	hxf.ifft(&(f1[k][0]));
  	  }

      for ( std::size_t i=0 ; i<Nx ; ++i ) {
        std::valarray<double> fi(Nv); for ( std::size_t k=0 ; k<Nv ; ++k ) { fi[k] = f1[k][i]; }
        hvf.fft(&(fi[0]));
        for ( std::size_t k=0 ; k<Nv ; ++k ) { hvf[k] = std::exp( -I*E(i)*kv[k]*dt )*hvf[k]; }
        hvf.ifft(&(fi[0]));
        for ( std::size_t k=0 ; k<Nv ; ++k ) { f2[k][i] = fi[k]; }
      }

  	  for ( auto k=0 ; k<f.size(0) ; ++k ) {
  	  	hxf.fft(&(f2[k][0]));
  	  	for ( auto i=0 ; i<Nx ; ++i ) { hxf[i] = std::exp( -I*v(k)*kx[i]*0.5*dt )*hxf[i]; }
  	  	hxf.ifft(&(f[k][0]));
  	  }

      ++i_t;
  	}
    //std::cout << " \r" << nb_iter << " : " << i_t << std::endl;
    field<double,1> diff=f;
    for ( std::size_t k=0 ; k<Nv ; ++k ) {
      for ( std::size_t i=0 ; i<Nx ; ++i ) {
        diff[k][i] -= f0[k][i];
      }
    }
    double e_1  = std::accumulate( diff.origin() , diff.origin()+diff.num_elements() , 0. , [&](double a,double b){return a+std::abs(b)*f0.step.dx*f0.step.dv;} );
    double e_oo = std::abs(*std::max_element( diff.origin() , diff.origin()+diff.num_elements() , [](double a,double b){return (std::abs(a) < std::abs(b));} ));
    std::cout << dt << " " << e_1 << " " << e_oo << std::endl;
  }

	return 0;
}
