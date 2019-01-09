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

template < typename _T , direction::direction D >
_T
local_flux_p ( const typename ublas::matrix_vector_slice_periodic<const ublas::matrix<_T>,D> & f_loc )
{
  static const _T epsi = 1e-6;
  _T w0,w1,w2;
  w0 = 13./12.*SQ( f_loc(0) - 2.*f_loc(1) + f_loc(2) ) + 0.25*SQ(    f_loc(0) - 4.*f_loc(1) + 3.*f_loc(2) );
  w1 = 13./12.*SQ( f_loc(1) - 2.*f_loc(2) + f_loc(3) ) + 0.25*SQ( f_loc(1) - f_loc(3) );
  w2 = 13./12.*SQ( f_loc(2) - 2.*f_loc(3) + f_loc(4) ) + 0.25*SQ( 3.*f_loc(2) - 4.*f_loc(3) +    f_loc(4) );

  w0 = 0.1/SQ(epsi+w0);
  w1 = 0.6/SQ(epsi+w1);
  w2 = 0.3/SQ(epsi+w2);

  _T sum_w = w0+w1+w2;
  w0 /= sum_w;
  w1 /= sum_w;
  w2 /= sum_w;

  return w0*( (2./6.)*f_loc(0) - (7./6.)*f_loc(1) + (11./6.)*f_loc(2) )
       + w1*(-(1./6.)*f_loc(1) + (5./6.)*f_loc(2) +  (2./6.)*f_loc(3) )
       + w2*( (2./6.)*f_loc(2) + (5./6.)*f_loc(3) +  (1./6.)*f_loc(4) );
}

template < typename _T , direction::direction D >
_T
local_flux_m ( const typename ublas::matrix_vector_slice_periodic<const ublas::matrix<_T>,D> & f_loc )
{
  static const _T epsi = 1e-6;
  _T w0,w1,w2;
  w0 = 13./12.*SQ( f_loc(2) - 2.*f_loc(3) + f_loc(4) ) + 0.25*SQ( 3.*f_loc(2) - 4.*f_loc(3) +    f_loc(4) );
  w1 = 13./12.*SQ( f_loc(1) - 2.*f_loc(2) + f_loc(3) ) + 0.25*SQ( f_loc(1) - f_loc(3) );
  w2 = 13./12.*SQ( f_loc(0) - 2.*f_loc(1) + f_loc(2) ) + 0.25*SQ(    f_loc(0) - 4.*f_loc(1) + 3.*f_loc(2) );

  w0 = 0.1/SQ(epsi+w0);
  w1 = 0.6/SQ(epsi+w1);
  w2 = 0.3/SQ(epsi+w2);

  _T sum_w = w0+w1+w2;
  w0 /= sum_w;
  w1 /= sum_w;
  w2 /= sum_w;

  return w2*(-(1./6.)*f_loc(0) + (5./6.)*f_loc(1) +  (2./6.)*f_loc(2) )
       + w1*( (2./6.)*f_loc(1) + (5./6.)*f_loc(2) -  (1./6.)*f_loc(3) )
       + w2*((11./6.)*f_loc(2) - (7./6.)*f_loc(3) +  (2./6.)*f_loc(4) );
}

#undef SQ

template < direction::direction D >
int &
index_direction ( int & i , int & k );

template <>
int &
index_direction<direction::x> ( int & i , int & k )
{ return i; }
template <>
int &
index_direction<direction::v> ( int & i , int & k )
{ return k; }

template < typename _T , direction::direction D >
auto flux_p ( const field<_T> & u , const ublas::vector<_T> & v )
{
  field<_T> fip12p(u.size1(),u.size2());
  fip12p.step = u.step; fip12p.range = u.range;

  for ( auto i=0 ; i<u.size1() ; ++i ) {
    for ( auto k=0 ; k<u.size2() ; ++k )
      { fip12p(i,k) = v(index_direction<D>(i,k))*local_flux_p( u.field<_T>::template stencil<D>(i,k) ); }
  }

  return fip12p;
}
template < typename _T , direction::direction D >
auto flux_m ( const field<_T> & u , const ublas::vector<_T> & v )
{
  field<_T> fip12m(u.size1(),u.size2());
  fip12m.step = u.step; fip12m.range = u.range;

  for ( auto i=0 ; i<u.size1() ; ++i ) {
    for ( auto k=0 ; k<u.size2() ; ++k )
      { fip12m(i,k) = v(index_direction<D>(i,k))*local_flux_m( u.field<_T>::template stencil<D>(i+1,k) ); }
  }

  return fip12m;
}
/*
template < typename _T >
auto
flux_p<direction::x> ( const field<_T> & u , const ublas::vector<_T> & v )
{
  field<_T> fip12p(u.size1(),u.size2());
  fip12p.step = u.step; fip12p.range = u.range;

  for ( auto i=0 ; i<u.size1() ; ++i ) {
    for ( auto k=0 ; k<u.size2() ; ++k )
      { fip12p(i,k) = local_flux_p( v(k)*u.stencil<direction::x>(i,k) ); }
  }

  return fip12p;
}

template < typename _T >
auto
flux_m<direction::x> ( const field<_T> & u , const ublas::vector<_T> & v )
{
  field<_T> fip12m(u.size1(),u.size2());
  fip12m.step = u.step; fip12m.range = u.range;

  for ( auto i=0 ; i<u.size1() ; ++i ) {
    for ( auto k=0 ; k<u.size2() ; ++k )
      { fip12m(i,k) = local_flux_m( v(k)*u.stencil<direction::x>(i+1,k) ); }
  }

  return fip12m;
}

template < typename _T >
auto
flux_p<direction::v> ( const field<_T> & u , const ublas::vector<_T> & v )
{
  field<_T> fip12p(u.size1(),u.size2());
  fip12p.step = u.step; fip12p.range = u.range;

  for ( auto i=0 ; i<u.size1() ; ++i ) {
    for ( auto k=0 ; k<u.size2() ; ++k )
      { fip12p(i,k) = local_flux_p( v(k)*u.stencil<direction::v>(i,k) ); }
  }

  return fip12p;
}

template < typename _T >
auto
flux_m<direction::v> ( const field<_T> & u , const ublas::vector<_T> & v )
{
  field<_T> fip12m(u.size1(),u.size2());
  fip12m.step = u.step; fip12m.range = u.range;

  for ( auto i=0 ; i<u.size1() ; ++i ) {
    for ( auto k=0 ; k<u.size2() ; ++k )
      { fip12m(i,k) = local_flux_m( v(k)*u.stencil<direction::v>(i,k+1) ); }
  }

  return fip12m;
}
*/
template < typename _T >
auto
trp2D ( const field<_T> & u , const ublas::vector<_T> & v , const ublas::vector<_T> & E )
{
  field<_T> trp(u.size1(),u.size2());
  trp.step = u.step; trp.range = u.range;

  ublas::vector<_T> vm (v.size()) , vp(v.size());
  ublas::vector<_T> Em (E.size()) , Ep(E.size());

  for ( auto k=0 ; k<v.size() ; ++k ) {
    vp(k) = std::max(v(k),0.);
    vm(k) = std::min(v(k),0.);
  }
  for ( auto i=0 ; i<E.size() ; ++i ) {
    Ep(i) = std::max(E(i),0.);
    Em(i) = std::min(E(i),0.);
  }

  field<_T> fp_x = flux_p<_T,direction::x>(u,vp);
  field<_T> fm_x = flux_m<_T,direction::x>(u,vm);
  field<_T> fp_v = flux_p<_T,direction::v>(u,Ep);
  field<_T> fm_v = flux_m<_T,direction::v>(u,Em);

  for ( auto i=1 ; i<trp.size1() ; ++i )
  {
    for ( auto k=1 ; k<trp.size2() ; ++k )
    {
      trp(i,k) = 1./trp.step.dx*( fp_x(i,k)-fp_x(i-1,k) + fm_x(i,k)-fm_x(i-1,k) )
               + 1./trp.step.dv*( fp_v(i,k)-fp_v(i,k-1) + fm_v(i,k)-fm_v(i,k-1) );
    }
    { auto k=0, km1=trp.size2()-1;
      trp(i,k) = 1./trp.step.dx*( fp_x(i,k)-fp_x(i-1,k) + fm_x(i,k)-fm_x(i-1,k) )
               + 1./trp.step.dv*( fp_v(i,k)-fp_v(i,km1) + fm_v(i,k)-fm_v(i,km1) );
    }
  }

  { auto i=0, im1=trp.size1()-1;
    for ( auto k=1 ; k<trp.size2() ; ++k )
    {
      trp(i,k) = 1./trp.step.dx*( fp_x(i,k)-fp_x(im1,k) + fm_x(i,k)-fm_x(im1,k) )
               + 1./trp.step.dv*( fp_v(i,k)-fp_v(i,k-1) + fm_v(i,k)-fm_v(i,k-1) );
    }
    { auto k=0, km1=trp.size2()-1;
      trp(i,k) = 1./trp.step.dx*( fp_x(i,k)-fp_x(im1,k) + fm_x(i,k)-fm_x(im1,k) )
               + 1./trp.step.dv*( fp_v(i,k)-fp_v(i,km1) + fm_v(i,k)-fm_v(i,km1) );
    }
  }

  return trp;
}

} //namespace weno

#endif

