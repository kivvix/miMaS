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
#define Xi(i) (i*f.step.dx+f.range.x_min)
#define Vk(k) (k*f.step.dv+f.range.v_min)

int main(int,char**)
{
  const std::size_t NumDimV = 1;
	const int Nx = 64, Nv = 128, Nb_iter=100;

	field<double,NumDimV> f( boost::extents[Nv][Nx] );

	f.range.v_min = -10.; f.range.v_max = 10;
	f.step.dv = (f.range.v_max-f.range.v_min)/Nv;
	f.range.x_min = 0.; f.range.x_max = 20.*math::pi<double>();
	f.step.dx = (f.range.x_max-f.range.x_min)/Nx;

	const double dt = 1.606*f.step.dv/0.6; //0.5*6.*math::pi<double>()/(Nv*f.range.v_max);

	ublas::vector<double> v (Nv,1.);
  ublas::vector<double> E (Nx,1.),rho(Nx);
  for ( std::size_t k=0 ; k<Nv ; ++k ) { v[k] = Vk(k); }
  //for ( std::size_t i=0 ; i<Nx ; ++i ) { E[i] = -Xi(i); }

	const double lx = f.range.x_max-f.range.x_min;
  const double lv = f.range.v_max-f.range.v_min;
	ublas::vector<double> kx(Nx),kv(Nv);
  for ( int i=0  ; i<Nx/2 ; ++i ) { kx[i]    = 2.*math::pi<double>()*i/lx; }
  for ( int i=-Nx/2 ; i<0 ; ++i ) { kx[Nx+i] = 2.*math::pi<double>()*i/lx; }

  for ( int k=0  ; k<Nv/2 ; ++k ) { kv[k]    = 2.*math::pi<double>()*k/lv; }
  for ( int k=-Nv/2 ; k<0 ; ++k ) { kv[Nv+k] = 2.*math::pi<double>()*k/lv; }
	
  double np = 0.9 , nb = 0.2 , ui = 4.5;
  for (field<double,NumDimV>::size_type k=0 ; k<f.size(0) ; ++k ) {
    for (field<double,NumDimV>::size_type i=0 ; i<f.size(1) ; ++i ) {
      f[k][i] = ( std::exp(-0.5*SQ(Vk(k)))*np/std::sqrt(2.*math::pi<double>()) + nb/std::sqrt(2.*math::pi<double>())*std::exp(-0.5*SQ(Vk(k)-ui)/0.25) )*(1.+0.04*std::cos(0.3*Xi(i)));
      //f[k][i] = std::cos(Xi(i)*0.2)*std::cos(math::pi<double>()*2.*Vk(k)/20);
      //f[k][i] = std::exp( -SQ(Xi(i)-1.)/2. )*std::exp( -Vk(k)*Vk(k)/1.);
      //f[k][i] = ( std::exp(-0.5*SQ(Vk(k)))*np/std::sqrt(2.*math::pi<double>()) )*(1.+0.04*std::cos(0.3*Xi(i)));
    }
  }
  f.write("init.dat");

  poisson<double> poisson_solver(Nx,lx);

  double Tf = 60.;//2*math::pi<double>();
  int i_t=0;

  std::cout << f.size(0) << "x" << f.size(1) << std::endl;
  std::cout << "dt " << dt << std::endl;
  std::cout << "dx " << f.step.dx << std::endl;
  std::cout << "dv " << f.step.dv << std::endl;
  std::cout << "Tf " << Tf << std::endl;

  ublas::vector<double> ee(int(Tf/dt)+1.);
  ublas::vector<double> Emax(int(Tf/dt)+1.);
  ublas::vector<double> H(int(Tf/dt)+1.);

  while ( i_t*dt < Tf ) {
    std::cout<<" \r"<<i_t<<" "<<std::flush;

    field<double,1> f1 = f;
    field<double,1> f2 = f;
	  
    fft::spectrum_ hxf(Nx);
	  for ( std::size_t k=0 ; k<Nv ; ++k ) {
	  	hxf.fft(&(f[k][0]));

	  	for ( std::size_t i=0 ; i<Nx ; ++i ) {
        hxf[i] = std::exp( -I*v(k)*kx[i]*0.5*dt )*hxf[i];
	  	}

	  	hxf.ifft(&(f1[k][0]));
	  }

    fft::spectrum_ hvf(Nv);
    rho = f1.density();
    E = poisson_solver(rho);
    for ( std::size_t i=0 ; i<Nx ; ++i ) {
      std::valarray<double> fi(Nv); for ( std::size_t k=0 ; k<Nv ; ++k ) { fi[k] = f1[k][i]; }
      hvf.fft(&(fi[0]));

      for ( std::size_t k=0 ; k<Nv ; ++k ) {
        hvf[k] = std::exp( -I*E(i)*kv[k]*dt )*hvf[k];
      }

      hvf.ifft(&(fi[0]));
      for ( std::size_t k=0 ; k<Nv ; ++k ) { f2[k][i] = fi[k]; }
    }

	  for ( auto k=0 ; k<f.size(0) ; ++k ) {
	  	hxf.fft(&(f2[k][0]));

	  	for ( auto i=0 ; i<Nx ; ++i ) {
        hxf[i] = std::exp( -I*v(k)*kx[i]*0.5*dt )*hxf[i];
	  	}

	  	hxf.ifft(&(f[k][0]));
	  }



    rho = f.density();
    E = poisson_solver(rho);
    ee(i_t) = 0.;
    for ( auto i=0 ; i<Nx ; ++i ) {
      ee(i_t) += SQ(E(i))*f.step.dx;
    }
    H(i_t) = energy(f,E);
    ++i_t;
	}

  std::cout << " \r" << i_t << std::endl;


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
  of.close();
  of.open("Emax.dat");
  std::transform( ee.begin() , ee.end() , std::ostream_iterator<std::string>(of,"\n") , dt_y );
  of.close();

	return 0;
}
