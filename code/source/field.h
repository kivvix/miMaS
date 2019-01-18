#ifndef _FIELD_H_
#define _FIELD_H_

#include <iostream>
#include <algorithm>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "matrix_proxy.h"

/**
# `field<T>`

This class represent a field in phase space in 2 dimensions (1D $x$ and 1D $v$). To simplify implementation, `field<T>` class is a daughter of `boost::numeric::ublas::maxtrix<T>`. If necessary, precision can be choose by choosing `T` type as `double` or `float` (in fact this is just easier for me to write everything inside header file).

All storing stuff is done by `ublas::matrix<T>` class.

## Helper classes

* `step`: represent space and velocity step, it store $\Delta x$ and $\Delta v$ as `dx` and `dv` as public value (you can access to the space step by `f.step.dx` for example).
* `range`: represent domain in phase space just as theirs limits $\[x_{\text{min}},x_{\text{max}}\]\times\[v_{\text{min}},v_{\text{max}}\]$ as `x_min`, `x_max`, `v_min` and `v_max` as public value (you can access to a born with `f.range.v_max` for example).

## Constructors

Default value for step is $\Delta x = 1$ and $\Delta v = 1$. Default value for range domain is $x_{\text{min}}=0$, $x_{\text{max}}=1$ for $x$ domain and $v_{\text{min}}=-1$, $v_{\text{max}}=1$ for $v$ domain.

* `field()`: default constructor
* `field ( typename ublas::matrix<T>::size_type size1 , typename ublas::matrix<T>::size_type size2 )`: `size1` size of domain in $x$ direction, `size2` size of domain in $v$ direction
* `field ( const field & )`: copy constructor
* `field ( const ublas::matrix_expression<AE> & )`: constructor from *matrix expression*
* `field ( typename ublas::matrix<T>::size_type size1 , typename ublas::matrix<T>::size_type size2 , T dx , T dv )`: `size1` size of domain in $x$ direction, `size2` size of domain in $v$ direction, `dx` and `dv` is respectively $\Delta x$ and $\Delta v$.

## Methods

* `auto density() const`: return $\rho$ as a `ublas::vector<T>`
* `auto moments() const`: return $\sum_k m(v_k)f_{i,k}\Delta v$ as a `std::array<ublas::vector<T>,3>`
* `auto flux() const`: return $\sum_k v_km(v_k)f_{i,k}\Delta v$ as a `std::array<ublas::vector<T>,3>`
**/

using namespace boost::numeric;
template <typename _T>
struct field : public ublas::matrix<_T>
{
  struct step
  {
    _T dx=1.,dv=1.;
    step()            = default;
    step(const step&) = default;
    step(step&&)      = default;
    step(_T dx_,_T dv_)
      : dx(dx_) , dv(dv_)
    {}
    step &
    operator=(const step & rhs)
    {
      dx = rhs.dx; dv = rhs.dv;
      return *this;
    }
  } step;
  struct range
  {
    _T x_min=0. ,x_max=1.;
    _T v_min=-1.,v_max=1.;
    range()             = default;
    range(const range&) = default;
    range(range&&)      = default;
    range(_T x_min_,_T x_max_,_T v_min_,_T v_max_)
      : x_min(x_min_) , x_max(x_max_) ,
        v_min(v_min_) , v_max(v_max_)
    {}
    range &
    operator=(const range & rhs)
    {
      x_min = rhs.x_min; x_max = rhs.x_max;
      v_min = rhs.v_min; v_max = rhs.v_max;
      return *this;
    }
  } range;

  field ()
    : ublas::matrix<_T>()
  {}
  
  field ( typename ublas::matrix<_T>::size_type size1 , typename ublas::matrix<_T>::size_type size2 )
    : ublas::matrix<_T>(size1,size2)
  {}

  field ( const field &m )
    : ublas::matrix<_T>(m) , step(m.step) , range(m.range)
  {}

  template<class AE>
  field ( const ublas::matrix_expression<AE> &ae )
    : ublas::matrix<_T>(ae)
  {}

  field ( typename ublas::matrix<_T>::size_type size1 , typename ublas::matrix<_T>::size_type size2 , _T dx , _T dv )
    : ublas::matrix<_T>(size1,size2) , step(dx,dv) 
  {}

  ~field ()
  {}

  field &
  operator = ( const field & rhs )
  {
    if ( this != &rhs )
    {
      this->ublas::matrix<_T>::operator=(rhs);
      range = rhs.range;
      step = rhs.step;
    }
    return *this;
  }

  auto
  density () const
  {
    static const ublas::scalar_vector<_T> v(this->size2(),step.dv);
    return prod(*this,v);
  }

#define SQ(X) ((X)*(X))
#define CU(X) ((X)*(X)*(X))
  auto
  moments () const
  { 
    _T v;
    std::array<ublas::vector<_T>,3> U = { ublas::vector<_T>(this->size1(),0.) , ublas::vector<_T>(this->size1(),0.) ,
                                          ublas::vector<_T>(this->size1(),0.) };

    for ( auto i=0 ; i<this->size1() ; ++i ) {
      for ( auto k=0 ; k<this->size2() ; ++k )
      {
        v = k*step.dv+range.v_min;
        U[0][i] += this->operator()(i,k)*step.dv;
        U[1][i] += this->operator()(i,k)*v*step.dv;
        U[2][i] += 0.5*this->operator()(i,k)*SQ(v)*step.dv;
      }
    }
    return U;
  }

  auto
  flux () const
  {
    _T v;
    std::array<ublas::vector<_T>,3> U = { ublas::vector<_T>(this->size1(),0.) , ublas::vector<_T>(this->size1(),0.) ,
                                          ublas::vector<_T>(this->size1(),0.) };

    for ( auto i=0 ; i<this->size1() ; ++i ) {
      for ( auto k=0 ; k<this->size2() ; ++k ) {
        v = k*step.dv+range.v_min;
        U[0](i) += this->operator()(i,k)*v*step.dv;
        U[1](i) += this->operator()(i,k)*SQ(v)*step.dv;
        U[2](i) += 0.5*this->operator()(i,k)*CU(v)*step.dv;
      }
    }
    return U;
  }
#undef SQ 
#undef CU 

  template <direction::direction D>
  inline typename std::enable_if< D == direction::x ,
  typename ublas::matrix_vector_slice_periodic<const ublas::matrix<_T>,direction::x> >::type
  stencil ( typename ublas::matrix<_T>::size_type i , typename ublas::matrix<_T>::size_type k ) const
  {
    return ublas::matrix_vector_slice_periodic<const ublas::matrix<_T>,direction::x>(*this,ublas::slice((i-2+this->size1())%this->size1(),1,5),ublas::slice(k,0,5));
  }
  template <direction::direction D>
  inline typename std::enable_if< D == direction::x ,
  typename ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,direction::x> >::type
  stencil ( typename ublas::matrix<_T>::size_type i , typename ublas::matrix<_T>::size_type k )
  {
    return ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,direction::x>(*this,ublas::slice((i-2+this->size1())%this->size1(),1,5),ublas::slice(k,0,5));
  }

  template <direction::direction D>
  inline typename std::enable_if< D == direction::v ,
  typename ublas::matrix_vector_slice_periodic<const ublas::matrix<_T>,direction::v> >::type
  stencil ( typename ublas::matrix<_T>::size_type i , typename ublas::matrix<_T>::size_type k ) const
  {
    return ublas::matrix_vector_slice_periodic<const ublas::matrix<_T>,direction::v>(*this,ublas::slice(i,0,5),ublas::slice((k-2+this->size2())%this->size2(),1,5));
  }
  template <direction::direction D>
  inline typename std::enable_if< D == direction::v ,
  typename ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,direction::v> >::type
  stencil ( typename ublas::matrix<_T>::size_type i , typename ublas::matrix<_T>::size_type k )
  {
    return ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,direction::v>(*this,ublas::slice(i,0,5),ublas::slice((k-2+this->size2())%this->size2(),1,5));
  }

/*
  template <direction::direction D>
  ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,D>
  stencil ( typename ublas::matrix<_T>::size_type i , typename ublas::matrix<_T>::size_type k );
*/
};

/*
template < typename _T >
template <>
ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,D>
field<_T>::stencil<direction::x>( typename ublas::matrix<_T>::size_type i , typename ublas::matrix<_T>::size_type k )
{
  return ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,direction::x>(*this,ublas::slice((i-2+this->size1())%this->size1(),1,5),ublas::slice(k,0,5));
}

template < typename _T >
template <>
ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,D>
field<_T>::stencil<direction::x>( typename ublas::matrix<_T>::size_type i , typename ublas::matrix<_T>::size_type k )
{
  return ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,direction::v>(*this,ublas::slice(i,0,5),ublas::slice((k-2+this->size2())%this->size2(),1,5));
}
*/
  
/*
template<>
template<typename _T> 
typename ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,direction::x>
field<_T>::stencil<direction::x> ( typename ublas::matrix<_T>::size_type i , typename ublas::matrix<_T>::size_type k )
{
  return ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,direction::x>(*this,ublas::slice(i-2,1,5),ublas::slice(k,0,5));
}
template<>
template<typename _T>
typename ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,direction::v>
field<_T>::stencil<direction::v> ( typename ublas::matrix<_T>::size_type i , typename ublas::matrix<_T>::size_type k )
{
  return ublas::matrix_vector_slice_periodic<ublas::matrix<_T>,direction::v>(*this,ublas::slice(i,0,5),ublas::slice(k-2,1,5));
}
*/

template <typename _T>
struct variabili
{
  ublas::vector<_T> rho;
  ublas::vector<_T> u;
  ublas::vector<_T> T;

  variabili ( std::array<ublas::vector<_T>,3> && U )
    : rho(std::move(U[0])) , u(std::move(U[1])) , T(U[2].size())
  {
    for ( auto i=0 ; i<rho.size() ; ++i )
    {
      u(i) /= rho(i);
      T(i)  = (U[2](i)-U[1](i)*U[1](i)/(2*rho(i)))/(0.5*rho(i));
    }
  }
};

// faire une factory initialisée par un step et un range (précalcul du vecteur de vitesses)
// puis génération pour chaque triplet de variabili du champ de la maxwelienne associée
template <typename _T>
field<_T>
maxwellian ( const variabili<_T> &v , typename field<_T>::step step , typename field<_T>::range range )
{
  field<_T> M((range.x_max-range.x_min)/step.dx,(range.v_max-range.v_min)/step.dv);
  M.step = step; M.range = range;
  
  for ( auto i=0 ; i<M.size1() ; ++i )
  {
    auto tmp = v.rho(i)/sqrt(2.*boost::math::constants::pi<_T>()*v.T(i));
    for ( auto k=0 ; k<M.size2() ; ++k )
    {
      auto vrel = (M.step.dv*k+M.range.v_min) - v.u(i);
      auto alpha = -0.5*vrel*vrel/v.T(i);
      M(i,k) = tmp*exp(alpha);
    }
  }

  return M;
}


#endif

