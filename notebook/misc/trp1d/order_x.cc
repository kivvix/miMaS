#include <iostream>
#include <valarray>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <numeric>

#include "bweno.h"
#include "weno.h"
#include "wenol.h"
#include "rk.h"

template < class Time_Method , typename Space_Method >
std::valarray<double>
trp_error ( std::valarray<double> const& u0 , std::valarray<double> const& uf , Time_Method rksn , Space_Method L , double v , double Tf , double dx , double dt )
{
  std::valarray<double> u = u0;
  for ( auto i_t=0 ; i_t*dt < Tf ; ++i_t ) { u = rksn(u,L,v,dx,dt); }
  
  // compute error
  std::valarray<double> err = uf-u;
  return err;
}

struct err {
  template<typename It>
  err(It first , It last, double dx)
  {
    e_1  = std::abs(*std::max_element( first , last , [](double a,double b){return (std::abs(a) < std::abs(b));} ));
    e_oo =   std::accumulate( first , last , 0. , [&](double a,double b){return a+std::abs(b)*dx; } );
  }

  double e_1,e_oo;
};

std::ostream&
operator << ( std::ostream & os , err const& e )
{
  os << e.e_1 << " " << e.e_oo; 
  return os;
}

int
main (int,char**)
{
  const std::size_t N=201;
  const double pi = 3.141592653589;
  const double dx = 2.*pi/N , dt = .1*dx/N;
  const double Tf = pi; // u_0 = cos(2x) donc u_f = u_0 avec Tf = pi

	std::size_t count;
  auto init_cos = [&,count=0]()mutable{ return std::cos(2.*(count++)*dx); };

  auto Lweno  = [](std::valarray<double> const& u,double v,double dx)->std::valarray<double>{ return  -weno::du(u,v,dx); };
  auto Lbweno = [](std::valarray<double> const& u,double v,double dx)->std::valarray<double>{ return -bweno::du(u,v,dx); };
  auto Lwenol = [](std::valarray<double> const& u,double v,double dx)->std::valarray<double>{ return -wenol::du(u,v,dx); };
  typedef decltype(Lweno)  weno_t;
  typedef decltype(Lbweno) bweno_t;
  typedef decltype(Lwenol) wenol_t;

  for ( int n=10;n<N;n+=10 ) {
    double _dx = 2.*pi/n;
    std::valarray<double> u0(n);
    std::generate(std::begin(u0),std::end(u0),[&,count=0]()mutable{ return std::cos(2.*(count++)*_dx); });
    std::valarray<double> err_weno = trp_error( u0 , u0 , rk::rk33<double,bweno_t> , Lbweno , 1. , Tf , _dx , dt );
    std::cout << u0.size() << " "  << _dx << " " << dt << " " << err(std::begin(err_weno),std::end(err_weno),_dx) << std::endl;
  }

  return 0;
}
