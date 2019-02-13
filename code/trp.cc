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

int main(int,char**)
{
  const int Nx = 9;
  double l=2.*math::pi<double>();
  double dx = l/Nx, dt = 0.5*dx;
  ublas::vector<double> v(Nx),kx(Nx);
  double vk=-0.5;

	for ( auto i=0 ; i<Nx/2+1 ; ++i )   { kx[i] = 2.*math::pi<double>()*i/l; }
	for ( auto i=0 ; i<((Nx/2)) ; ++i ) { kx[i+Nx/2+1] = -kx[Nx/2-i]; }

  for ( auto i=0 ; i<Nx ; ++i ) {
    v[i] = std::cos((double)i*dx);
  }

  std::ofstream f0("init.dat");
  for ( auto i=0 ; i<Nx ; ++i ) { f0 << i*dx << " " << v[i] << "\n"; }
  f0 << std::endl;
  f0.close();

  fft::spectrum s(Nx);

std::cout << std::endl;
  int Nb_iter = 100;
  for ( auto t=0 ; t<Nb_iter ; ++t ) {
    std::cout << t*dt << "----\n";
    s.fft(&v[0]);
    for ( auto i=0 ; i<Nx ; ++i ) {
      std::cout << kx[i] << " " << s[i][fft::re] << " " << s[i][fft::im] <<  "\n";
    }
    std::cout << "\n----\n" ; 

    for ( auto i=0 ; i<Nx ; ++i ) {
      double re = s[i][fft::re], im = s[i][fft::im];
      s[i][fft::re] = std::cos(kx[i]*dt*vk)*re + std::sin(kx[i]*dt*vk)*im;
      s[i][fft::im] = std::cos(kx[i]*dt*vk)*im - std::sin(kx[i]*dt*vk)*re;
    }
    s.ifft(&v[0]);
  }
std::cout << std::endl;

  std::ofstream f1("vp.dat");
  for ( auto i=0 ; i<Nx ; ++i ) {
    f1 << i*dx << " " << v[i] << " " << std::cos(i*dx-vk*(Nb_iter)*dt) << "\n";
  }
  f1 << std::endl;
  f1.close();

  return 0;
}








