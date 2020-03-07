#include <ostream>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <functional>
#include <cmath>
#include <sstream>

#include <boost/numeric/ublas/matrix.hpp>

#include "miMaS/complex_field.h"

namespace iteration {
using namespace boost::numeric;

template < typename _T >
struct iteration {
  std::size_t iter;
  _T dt;
  _T current_time;
  _T Lhfh=0., LE=0.;
  bool success;

  template < typename Iterator >
  static _T
  error ( Iterator first1 , Iterator last1 , Iterator first2 , _T integrator_step )
  {
    return std::sqrt(std::inner_product(
      first1 , last1 , first2 , _T{0.} ,
      std::plus<_T>{} ,
      [&] ( const auto & a , const auto & b ) { return std::pow(std::abs( a - b ),2)*integrator_step; }
    ));
  }

  //template < std::size_t NumDimsV >
  _T
  hfh_error ( const complex_field<_T,1> & hfh1 , const complex_field<_T,1> & hfh2 , _T integrator_step )
  {
    Lhfh = error( hfh1.origin() , hfh1.origin()+hfh1.num_elements() , hfh2.origin() , integrator_step );
    success = false;
    return Lhfh;
  }

  _T
  E_error ( const ublas::vector<_T> & E1 , const ublas::vector<_T> & E2 , _T integrator_step )
  {
    LE = error( E1.begin() , E1.end() , E2.begin() , integrator_step );
    success = false;
    return LE;
  }

  const _T &
  increment ()
  {
    current_time += dt;
    return current_time;
  }

  const _T &
  error () const
  { return Lhfh; }

};

template < typename CharT , typename Traits = std::char_traits<CharT> , typename _T >
std::basic_ostream<CharT,Traits> &
operator << ( std::basic_ostream<CharT,Traits> & os , const iteration<_T> & iter )
{
  os << iter.iter << " " << iter.dt << " " << iter.current_time << " " << iter.Lhfh << " " << iter.LE << " " << std::noboolalpha << iter.success;
  return os;
}

template< typename _T >
struct __time_iteration
{
  std::size_t i;
  _T t;
  _T dt;

  __time_iteration ( std::size_t _i , _T _t , _T _dt )
    : i(_i) , t(_t) , dt(_dt)
  { ; }
};
template < typename CharT , typename Traits = std::char_traits<CharT> , typename _T >
std::basic_ostream<CharT,Traits> &
operator << ( std::basic_ostream<CharT,Traits> & os , const __time_iteration<_T> & ti )
{
  os << " [" << std::setw(6) << ti.i << "] " << std::setw(8) << ti.t << " (" << std::setw(9) << ti.dt << ")";
  return os;
}

template < typename _T >
__time_iteration<_T>
time ( const iteration<_T> & iter )
{
  return __time_iteration(iter.iter,iter.current_time,iter.dt);
}

template< typename _T >
struct __error_iteration
{
  bool success;
  _T Lhfh,LE;

  __error_iteration ( bool _success , _T _Lhfh , _T _LE )
    : success(_success) , Lhfh(_Lhfh) , LE(_LE)
  { ; }
};
template < typename CharT , typename Traits = std::char_traits<CharT> , typename _T >
std::basic_ostream<CharT,Traits> &
operator << ( std::basic_ostream<CharT,Traits> & os , const __error_iteration<_T> & ei )
{
  os << ((ei.success)?"\033[92m"s:"\033[31m"s) << std::setw(10) << ei.Lhfh << " " << std::setw(10) << ei.LE << "\033[0m";
  return os;
}

template < typename _T >
__error_iteration<_T>
error ( const iteration<_T> & iter )
{
  return __error_iteration(iter.success,iter.Lhfh,iter.LE);
}

} // namespace iteration
