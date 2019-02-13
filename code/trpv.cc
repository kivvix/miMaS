#include <iostream>
#include <algorithm>
#include <iterator>
#include <cmath>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "miMaS/field.h"
#include "miMaS/weno.h"
#include "miMaS/fft.h"
#include "miMaS/array_view.h"
#include "miMaS/poisson.h"

namespace math = boost::math::constants;

int main(int,char**)
{
	std::size_t Nx = 125, Nv = 64;
	field<double,1> f(boost::extents[Nv][Nx]);

	f.range.v_min = 0.; f.range.v_max = 2.*math::pi<double>();
	f.step.dv = (f.range.v_max-f.range.v_min)/Nv;
	f.range.x_min = 0.; f.range.x_max = 1.;
	f.step.dx = (f.range.x_max-f.range.x_min)/Nx;
	double dt = 0.5*f.step.dv;

  field<double,1> f_sol = f;
	
	ublas::vector<double> E (Nx);
  for ( std::size_t i=0 ; i<Nx ; ++i ) { E[i] = 1.; }

#define X(i) (i*f.step.dx+f.range.x_min)
#define V(k) (k*f.step.dv+f.range.v_min)
  for (field<double,2>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,2>::size_type i=0 ; i<f.size(1) ; ++i ) {
      f[k][i] = cos(V(k));
      f_sol[k][i] = cos(V(k)-dt*20);
    }
  }
#undef X
  f.write("init.dat");

  for ( auto t=0 ; t<20 ; ++t ) {
  	if (t%32==0) { std::cout<<"\r"<<t<<" "<<std::flush ; }
	  field<double,1> Edvf = weno::trp_v(f,E);

    for ( auto k=0 ; k<f.size(0) ; ++k ) {
      for ( auto i=0 ; i<f.size(1) ; ++i ) {
        f[k][i] = f[k][i] - dt*Edvf[k][i];
      }
    }
  }

  f_sol.write("sol.dat");
  f.write("vp.dat");

  return 0;
}

