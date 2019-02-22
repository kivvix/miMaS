#include <iostream>
#include <valarray>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <sstream>

#include "bweno.h"
#include "weno.h"
#include "wenol.h"
#include "rk.h"

int
main (int,char**)
{
  const std::size_t N=100;
  const double pi = 3.141592653589;
  const double dx = 2.*pi/N , dt = 1.*dx;
  const double Tf = 2.*pi;

  std::valarray<double> u(N);
	std::size_t count;

  auto cren_0 = [&](std::size_t n)->std::valarray<double> {
    std::valarray<double> u0(n);
    int i=0;
    for ( i=0;i*dx<1.;++i ) { u0[i] = dx*i; }
    for (    ;i*dx<4.;++i ) { u0[i] = 1.;   }
    for (    ;i<n    ;++i ) { u0[i] = 0.;   }
    return u0;
  };
  auto cos_0 = [&](std::size_t n)->std::valarray<double> {
    std::valarray<double> u0(n);
    for ( auto i=0;i<n;++i ) { u0[i] = std::cos(2.*i*dx); }
    return u0;
  };
  auto cosHF_0 = [&](std::size_t n)->std::valarray<double> {
    std::valarray<double> u0(n);
    for ( auto i=0;i<n;++i ) { u0[i] = std::cos(10.*i*dx); }
    return u0;
  };

  auto chapi_chap0 = [&](std::size_t n)->std::valarray<double> {
    std::valarray<double> u0(n);
    int i=0;
    for ( i=0;i*dx<1.;++i ) { u0[i] = 0.;       }
    for (    ;i*dx<2.;++i ) { u0[i] =  dx*i-1.; }
    for (    ;i*dx<3.;++i ) { u0[i] = -dx*i+3.; }
    for (    ;i<n    ;++i ) { u0[i] = 0.;       }
    return u0;
  };
  auto cren0 = [&,count=0]()mutable->double {
    if ( (count++)*dx < 1. ) { return 2.; }
    return 1.;
  };

  auto u_0 = cren0;

  auto dx_y = [&,count=0]( auto const& x ) mutable { std::stringstream ss; ss<<dx*count++<<" "<<x; return ss.str(); };
	std::ofstream of;
    
  //u = u_0(N);
  std::generate(std::begin(u),std::end(u),u_0);
  of.open("init.dat");
	std::transform( std::begin(u) , std::end(u) , std::ostream_iterator<std::string>(of,"\n") , dx_y );
	of.close();

  auto Lweno  = [](std::valarray<double> const& u,double v,double dx)->std::valarray<double>{ return  -weno::du(u,v,dx); };
  auto Lbweno = [](std::valarray<double> const& u,double v,double dx)->std::valarray<double>{ return -bweno::du(u,v,dx); };
  auto Lwenol = [](std::valarray<double> const& u,double v,double dx)->std::valarray<double>{ return -wenol::du(u,v,dx); };

  std::cout << dx << " " << dt << std::endl;

  std::cout << "WENO" << std::endl;
  for ( auto i_t=0 ; i_t*dt < Tf ; ++i_t )
    { u = rk::rk53(u,Lweno,1.,dx,dt); }
  of.open("vp_weno.dat");
	std::transform( std::begin(u) , std::end(u) , std::ostream_iterator<std::string>(of,"\n") , dx_y );
	of.close();

  std::cout << "BWENO" << std::endl;
  //u = u_0(N);
  std::generate(std::begin(u),std::end(u),u_0);
  for ( auto i_t=0 ; i_t*dt < Tf ; ++i_t )
    { u = rk::rk33(u,Lbweno,1.,dx,dt); }
  of.open("vp_bweno.dat");
	std::transform( std::begin(u) , std::end(u) , std::ostream_iterator<std::string>(of,"\n") , dx_y );
	of.close();

  std::cout << "WENOl" << std::endl;
  //u = u_0(N);
  std::generate(std::begin(u),std::end(u),u_0);
  for ( auto i_t=0 ; i_t*dt < Tf ; ++i_t )
    { u = rk::rk33(u,Lwenol,1.,dx,dt); }
  of.open("vp_wenol.dat");
	std::transform( std::begin(u) , std::end(u) , std::ostream_iterator<std::string>(of,"\n") , dx_y );
	of.close();


  return 0;
}
