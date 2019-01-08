#ifndef _WENO_H_
#define _WENO_H_

#include <iostream>
#include <algorithm>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "matrix_proxy.h"
#include "field.h"

namespace weno {
using namespace boost::numeric;
#define SQ(X) ((X)*(X))

template < typename _T >
_T
local_flux_p ( const typename ublas::matrix_vector_slice<ublas::matrix<_T>> & f_loc )
{
  static const _T epsi = 1e-6;
  _T w0,w1,w2;
  w0 = 13./12.*SQ( f_loc(0) - 2.*f_loc(1) + f_loc(2) ) + 0.25*SQ(    f(0) - 4.*f(1) + 3.*f(2) );
  w1 = 13./12.*SQ( f_loc(1) - 2.*f_loc(2) + f_loc(3) ) + 0.25*SQ( f(1) - f(3) );
  w2 = 13./12.*SQ( f_loc(2) - 2.*f_loc(3) + f_loc(4) ) + 0.25*SQ( 3.*f(2) - 4.*f(3) +    f(4) );

  w0 = 0.1/SQ(epsi+w0);
  w1 = 0.6/SQ(epsi+w1);
  w2 = 0.3/SQ(epsi+w2);

  _T sum_w = w0+w1+w2;
  w0 /= sum_w;
  w1 /= sum_w;
  w2 /= sum_w;

  return w0*( (2./6.)*f_loc(0) - (7./6.)*f(1) + (11./6.)*f(2) )
       + w1*(-(1./6.)*f_loc(1) + (5./6.)*f(2) +  (2./6.)*f(3) )
       + w2*( (2./6.)*f_loc(2) + (5./6.)*f(3) +  (1./6.)*f(4) );
}

template < typename _T >
_T
local_flux_m ( const typename ublas::matrix_vector_slice<ublas::matrix<_T>> & f_loc )
{
  static const _T epsi = 1e-6;
  _T w0,w1,w2;
  w0 = 13./12.*SQ( f_loc(2) - 2.*f_loc(3) + f_loc(4) ) + 0.25*SQ( 3.*f(2) - 4.*f(3) +    f(4) );
  w1 = 13./12.*SQ( f_loc(1) - 2.*f_loc(2) + f_loc(3) ) + 0.25*SQ( f(1) - f(3) );
  w2 = 13./12.*SQ( f_loc(0) - 2.*f_loc(1) + f_loc(2) ) + 0.25*SQ(    f(0) - 4.*f(1) + 3.*f(2) );

  w0 = 0.1/SQ(epsi+w0);
  w1 = 0.6/SQ(epsi+w1);
  w2 = 0.3/SQ(epsi+w2);

  _T sum_w = w0+w1+w2;
  w0 /= sum_w;
  w1 /= sum_w;
  w2 /= sum_w;

  return w2*(-(1./6.)*f_loc(0) + (5./6.)*f(1) +  (2./6.)*f(2) )
       + w1*( (2./6.)*f_loc(1) + (5./6.)*f(2) -  (1./6.)*f(3) )
       + w2*((11./6.)*f_loc(2) - (7./6.)*f(3) +  (2./6.)*f(4) );
}

#undef SQ

template < typename _T >
auto
flux_p ( const field<_T> & u , const ublas::vector<_T> & v )
{
  field<_T> fip12p(u.size1(),u.size2());
  fip12p.step = u.step; fip12p.range = u.range;

  for ( auto i=2 ; i<u.size1()-2 ; ++i ) {
    for ( auto k=0 ; k<u.size2() ; ++k )
      { fip12p(i,k) = local_flux_p( v(k)*u.stencil<direction::x>(i,k) ); }
  }

  for ( auto i=0 ; i<2 ; ++i ) {
    typename ublas::matrix<_T>::size_type im2 = (i-2+u.size1())%u.size1() , im1 = (i-1+u.size1())%u.size1();
    for ( auto k=0 ; k<u.size2() ; ++k )
     { fip12p(i,k) = local_flux_p( v(k)*{ u(im2,k) , u(im1,k) , u(i,k) , u(i+1,k) , u(i+2,k) } ); }
  }

  for ( auto i=u.size1()-2 ; i<u.size1() ; ++i ) {
    typename ublas::matrix<_T>::size_type ip1 = (i+1)%u.size1() , ip2 = (i+2)%u.size1();
    for ( auto k=0 ; k<u.size2() ; ++k )
     { fip12p(i,k) = local_flux_p( v(k)*{ u(i-2,k) , u(i-1,k) , u(i,k) , u(ip1,k) , u(ip2,k) } ); }
  }

  return fip12p;
}

template < typename _T >
auto
flux_m ( const field<_T> & u , const ublas::vector<_T> & v )
{
  field<_T> fip12p(u.size1(),u.size2());
  fip12p.step = u.step; fip12p.range = u.range;

  for ( auto i=1 ; i<u.size1()-3 ; ++i ) {
    for ( auto k=0 ; k<u.size2() ; ++k )
      { fip12m(i,k) = local_flux_p( v(k)*u.stencil<direction::x>(i+1,k) ); }
  }

  { auto i=0;
    typename ublas::matrix<_T>::size_type im1 = u.size1()-1;
    for ( auto k=0 ; k<u.size2() ; ++k )
     { fip12m(i,k) = local_flux_p( v(k)*{ u(im1,k) , u(i,k) , u(i+1,k) , u(i+2,k) , u(i+3,k) } ); }
  }

  for ( auto i=u.size1()-3 ; i<u.size1() ; ++i ) {
    typename ublas::matrix<_T>::size_type ip1 = (i+1)%u.size1() , ip2 = (i+2)%u.size1() , ip3 = (i+3)%u.size1();
    for ( auto k=0 ; k<u.size2() ; ++k )
     { fip12m(i,k) = local_flux_p( v(k)*{ u(i-1,k) , u(i,k) , u(ip1,k) , u(ip2,k) , u(ip3,k) } ); }
  }

  return fip12m;
}


}

#endif

