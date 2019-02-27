#include <iostream>
#include <valarray>
#include <algorithm>
#include <iterator>
#include <fstream>
#include <numeric>
#include <functional>
#include <vector>
#include <sstream>

#include "bweno.h"
#include "weno.h"
#include "wenol.h"
#include "rk.h"

template < typename _T , class Time_Method , typename Space_Method >
std::valarray<_T>
trp_error ( std::valarray<_T> const& u0 , std::valarray<_T> const& uf , Time_Method rksn , Space_Method L , _T v , _T Tf , _T dx , _T dt )
{
  std::valarray<_T> u = u0;
  for ( auto i_t=0 ; i_t*dt < Tf ; ++i_t ) { u = rksn(u,L,v,dx,dt); }
  
  // compute error
  std::valarray<_T> err = (u-uf)/uf;
  return err;
}

template < typename _T >
struct error {
  template<typename It>
  error(It first , It last, _T dx)
  {
    e_1  = std::accumulate( first , last , _T{0.} , [&](_T a,_T b){return a+std::abs(b)*dx; } );
    e_oo = std::abs(*std::max_element( first , last , [](_T a,_T b){return (std::abs(a) < std::abs(b));} ));
  }

  _T e_1,e_oo;
};

template < typename _T >
std::ostream&
operator << ( std::ostream & os , error<_T> const& err )
{
  os << err.e_1 << " " << err.e_oo; 
  return os;
}

auto Lweno  = [](std::valarray<long double> const& u,long double v,long double dx)->std::valarray<long double>{ return  -weno::du(u,v,dx); };
typedef decltype(Lweno) weno_t;

struct rksn {
  std::string name;
  std::size_t stage;
  std::function<std::valarray<long double>(std::valarray<long double>const&,weno_t,long double,long double,long double)> func;
  std::string tag;

  template < typename FUNC >
  rksn ( std::string name_ , std::size_t stage_ , FUNC func_ , std::string tag_ )
    : name(name_) , stage(stage_) , func(func_) , tag(tag_)
  {}
};

int
main (int,char**)
{
  const std::size_t N=10000;
  const long double pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628620899862803482534211706798214808651328230664709384460955058223172535940812848111745028410270l;
  const long double dx = 2.l*pi/N , dt = 10.l/N;
  const long double Tf = 1e0*dt;
  const std::size_t n_iter = 5;

	std::size_t count;
  auto cos0 = [&,count=0]()mutable{ return std::cos(7.l*(count++)*dx); };
  auto cosf = [&,count=0]()mutable{ return std::cos(7.l*((count++)*dx-Tf)); };

  auto chapi_chap0 = [&,count=0]()mutable->long double {
    if ( count*dx<1. ) { ++count; return 1.;       }
    if ( count*dx<2. ) { return dx*(count++);  }
    if ( count*dx<3. ) { return -dx*(count++)+4.; }
    ++count; return 1.;
  };
  auto cren0 = [&,count=0]()mutable->long double {
    if ( (count++)*dx < 1. ) { return 2.; }
    return 1.;
  };

  //auto init = cos0;

  std::valarray<long double> u0(N),uf(N);
  std::generate(std::begin(u0),std::end(u0),cos0);
  std::generate(std::begin(uf),std::end(uf),cosf);

  std::cout << Tf << std::endl;
  std::vector<rksn> time_methods;
  time_methods.push_back(rksn(    "RK (8,6)",8,rk::rk86<long double,weno_t>,"rk86"));
  time_methods.push_back(rksn(    "RK (7,6)",7,rk::rk76<long double,weno_t>,"rk76"));
  time_methods.push_back(rksn( "RK SSP(3,3)",3,rk::rk33<long double,weno_t>,"rk33"));
  time_methods.push_back(rksn( "RK SSP(4,3)",4,rk::rk43<long double,weno_t>,"rk43"));
  time_methods.push_back(rksn("RK NSSP(5,3)",5,rk::rk53<long double,weno_t>,"rk53"));
  time_methods.push_back(rksn( "RK SSP(4,4)",4,rk::rk44<long double,weno_t>,"rk44"));

  std::ofstream of;
 /* 
  auto dx_y = [&,count=0]( auto const& x ) mutable { std::stringstream ss; ss<<dx*count++<<" "<<x; return ss.str(); };
  of.open("init.dat");
  std::transform(std::begin(u0),std::end(u0),std::ostream_iterator<std::string>(of,"\n"),dx_y);
  of.close();
  of.open("sol.dat");
  std::transform(std::begin(uf),std::end(uf),std::ostream_iterator<std::string>(of,"\n"),dx_y);
  of.close();

  std::valarray<double> u=u0;
  for ( int i_t=0 ; i_t*dt*0.1<Tf ; ++i_t ) {
    u = rk::rk33(u,Lweno,1.,dx,dt*0.1);
  }
  of.open("vp.dat");
  std::transform(std::begin(u),std::end(u),std::ostream_iterator<std::string>(of,"\n"),dx_y);
  of.close();
*/
  for ( auto const& rk : time_methods ) {
    std::cout << rk.name << std::endl; of.open("rk_err/"+rk.tag+".dat");
    for ( int n=1;n<=n_iter;++n ) {
      long double _dt = dt/n;
      std::valarray<long double> err_weno = trp_error( u0,uf, rk.func , Lweno , 1.l , Tf , dx , _dt );
      std::cout << "\r " << n << std::flush;
      of << n << " " << dx << " " << _dt/rk.stage << " " << error<long double>(std::begin(err_weno),std::end(err_weno),dx) << std::endl;
    }
    std::cout << std::endl; of.close();
  }

  return 0;
}
