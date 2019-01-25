#ifndef _WENO_H_
#define _WENO_H_

#include <algorithm>
#include <iostream>
#include <iterator>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multi_array.hpp>

#include "direction.h"
#include "field.h"
#include "array_view.h"

namespace weno {

using namespace boost::numeric;
#define SQ(X) ((X)*(X))

template < class Collection >
auto
local_flux_p ( const Collection & f_loc )
{
  typedef typename std::remove_cv<typename Collection::value_type>::type _T;
  const static _T epsi = 1e-4;
  _T b0 = (13./12.)*SQ( f_loc[0] - 2.*f_loc[1] + f_loc[2] ) + 0.25*SQ(    f_loc[0] - 4.*f_loc[1] + 3.*f_loc[2] );
  _T b1 = (13./12.)*SQ( f_loc[1] - 2.*f_loc[2] + f_loc[3] ) + 0.25*SQ(    f_loc[1]               -    f_loc[3] );
  _T b2 = (13./12.)*SQ( f_loc[2] - 2.*f_loc[3] + f_loc[4] ) + 0.25*SQ( 3.*f_loc[2] - 4.*f_loc[3] +    f_loc[4] );

  _T a0 = 0.1/(SQ( epsi + b0 ));
  _T a1 = 0.6/(SQ( epsi + b1 ));
  _T a2 = 0.3/(SQ( epsi + b2 ));

  _T sum_w = a0+a1+a2;

  _T w0 = a0/sum_w;
  _T w1 = a1/sum_w;
  _T w2 = a2/sum_w;

  //std::copy( f_loc.begin() , f_loc.end() , std::ostream_iterator<_T>(std::cout," ") );
  //std::cout << "\t" << w0 << " " << w1 << " " << w2 << std::endl;
  //w0 = 0.1; w1 = 0.6; w2 = 0.3;

  return w0*( (2./6.)*f_loc[0] - (7./6.)*f_loc[1] + (11./6.)*f_loc[2] )
       + w1*(-(1./6.)*f_loc[1] + (5./6.)*f_loc[2] +  (2./6.)*f_loc[3] )
       + w2*( (2./6.)*f_loc[2] + (5./6.)*f_loc[3] -  (1./6.)*f_loc[4] );
}

template < class Collection >
auto
local_flux_m ( const Collection & f_loc )
{
  typedef typename std::remove_cv<typename Collection::value_type>::type _T;
  const static _T epsi = 1e-6;
  _T w0 = (13./12.)*SQ( f_loc[2] - 2.*f_loc[3] + f_loc[4] ) + 0.25*SQ( 3.*f_loc[2] - 4.*f_loc[3] +    f_loc[4] );
  _T w1 = (13./12.)*SQ( f_loc[1] - 2.*f_loc[2] + f_loc[3] ) + 0.25*SQ(    f_loc[1]               -    f_loc[3] );
  _T w2 = (13./12.)*SQ( f_loc[0] - 2.*f_loc[1] + f_loc[2] ) + 0.25*SQ(    f_loc[0] - 4.*f_loc[1] + 3.*f_loc[2] );

  w0 = 0.1/SQ(w0+epsi);
  w1 = 0.6/SQ(w1+epsi);
  w2 = 0.3/SQ(w2+epsi);

  _T sum_w = w0+w1+w2;
  w0 /= sum_w;
  w1 /= sum_w;
  w2 /= sum_w;

  //w0 = 0.1; w1 = 0.6; w2 = 0.3;

  return w2*(-(1./6.)*f_loc[0] + (5./6.)*f_loc[1] + (2./6.)*f_loc[2] )
       + w1*( (2./6.)*f_loc[1] + (5./6.)*f_loc[2] - (1./6.)*f_loc[3] )
       + w0*((11./6.)*f_loc[2] - (7./6.)*f_loc[3] + (2./6.)*f_loc[4] );
}

#undef SQ

template < direction::direction D> inline std::size_t index_velocity ( std::size_t i , std::size_t k );
template <> inline std::size_t index_velocity<direction::x> (std::size_t i , std::size_t k ) { return k; }
template <> inline std::size_t index_velocity<direction::v> (std::size_t i , std::size_t k ) { return i; }

template < direction::direction D > inline std::size_t index_i_m ( std::size_t i );
template <> inline std::size_t index_i_m<direction::x> ( std::size_t i ) { return i+1; }
template <> inline std::size_t index_i_m<direction::v> ( std::size_t i ) { return i;   }
template < direction::direction D > inline std::size_t index_k_m ( std::size_t k );
template <> inline std::size_t index_k_m<direction::x> ( std::size_t k ) { return k;   }
template <> inline std::size_t index_k_m<direction::v> ( std::size_t k ) { return k+1; }


template < typename _T , std::size_t NumDims , direction::direction D >
typename std::enable_if< D == direction::x ,
         boost::multi_array<_T,NumDims> >::type
flux_p ( field<_T,NumDims> const& u , ublas::vector<_T> const& velocity )
{
  boost::multi_array<_T,NumDims> fip12(tools::array_view<const std::size_t>(u.shape(),NumDims));

  for ( std::size_t i=0 ; i<2 ; ++i ) {
    for ( std::size_t k=0 ; k<u.shape()[1] ; ++k ) {
      fip12[i][k] = velocity(index_velocity<D>(i,k))*local_flux_p( u.field<_T,NumDims>::template stencil_border<D>(i,k) );
    }
  }
  for ( std::size_t i=2 ; i<u.shape()[0]-2 ; ++i ) {
    for ( std::size_t k=0 ; k<u.shape()[1] ; ++k ) {
      fip12[i][k] = velocity(index_velocity<D>(i,k))*local_flux_p( u.field<_T,NumDims>::template stencil<D>(i,k) );
    }
  }

  for ( std::size_t i=u.shape()[0]-2 ; i<u.shape()[0] ; ++i ) {
    for ( std::size_t k=0 ; k<u.shape()[1] ; ++k ) {
      fip12[i][k] = velocity(index_velocity<D>(i,k))*local_flux_p( u.field<_T,NumDims>::template stencil_border<D>(i,k) );
    }
  }

  return fip12;
}

template < typename _T , std::size_t NumDims , direction::direction D >
typename std::enable_if< D == direction::v ,
         boost::multi_array<_T,NumDims> >::type
flux_p ( field<_T,NumDims> const& u , ublas::vector<_T> const& velocity )
{
  boost::multi_array<_T,NumDims> fip12(tools::array_view<const std::size_t>(u.shape(),NumDims));

  for ( std::size_t i=0 ; i<u.shape()[0] ; ++i ) {
    for ( std::size_t k=0 ; k<2 ; ++k ) {
      fip12[i][k] = velocity(index_velocity<D>(i,k))*local_flux_p( u.field<_T,NumDims>::template stencil_border<D>(i,k) );
    }
    for ( std::size_t k=2 ; k<u.shape()[1]-2 ; ++k ) {
      fip12[i][k] = velocity(index_velocity<D>(i,k))*local_flux_p( u.field<_T,NumDims>::template stencil<D>(i,k) );
    }
    for ( std::size_t k=u.shape()[1]-2 ; k<u.shape()[1] ; ++k ) {
      fip12[i][k] = velocity(index_velocity<D>(i,k))*local_flux_p( u.field<_T,NumDims>::template stencil_border<D>(i,k) );
    }
  }

  return fip12;
}


template < typename _T , std::size_t NumDims , direction::direction D >
typename std::enable_if< D == direction::x ,
         boost::multi_array<_T,NumDims> >::type
flux_m ( field<_T,NumDims> const& u , ublas::vector<_T> const& velocity )
{
  boost::multi_array<_T,NumDims> fim12(tools::array_view<const std::size_t>(u.shape(),NumDims));

  for ( std::size_t i=0 ; i<1 ; ++i ) {
    for ( std::size_t k=0 ; k<u.shape()[1] ; ++k ) {
      fim12[i][k] = velocity(index_velocity<D>(i,k))*local_flux_m( u.field<_T,NumDims>::template stencil_border<D>(i+1,k) );
    }
  }
  for ( std::size_t i=1 ; i<u.shape()[0]-3 ; ++i ) {
    for ( std::size_t k=0 ; k<u.shape()[1] ; ++k ) {
      fim12[i][k] = velocity(index_velocity<D>(i,k))*local_flux_m( u.field<_T,NumDims>::template stencil<D>(i+1,k) );
    }
  }
  for ( std::size_t i=u.shape()[0]-3 ; i<u.shape()[0] ; ++i ) {
    for ( std::size_t k=0 ; k<u.shape()[1] ; ++k ) {
      fim12[i][k] = velocity(index_velocity<D>(i,k))*local_flux_m( u.field<_T,NumDims>::template stencil_border<D>(i+1,k) );
    }
  }

  return fim12;
}

template < typename _T , std::size_t NumDims , direction::direction D >
typename std::enable_if< D == direction::v ,
         boost::multi_array<_T,NumDims> >::type
flux_m ( field<_T,NumDims> const& u , ublas::vector<_T> const& velocity )
{
  boost::multi_array<_T,NumDims> fim12(tools::array_view<const std::size_t>(u.shape(),NumDims));

  for ( std::size_t i=0 ; i<u.shape()[0] ; ++i ) {
    for ( std::size_t k=0 ; k<1 ; ++k ) {
      fim12[i][k] = velocity(index_velocity<D>(i,k))*local_flux_m( u.field<_T,NumDims>::template stencil_border<D>(i,k+1) );
    }
    for ( std::size_t k=1 ; k<u.shape()[1]-3 ; ++k ) {
      fim12[i][k] = velocity(index_velocity<D>(i,k))*local_flux_m( u.field<_T,NumDims>::template stencil<D>(i,k+1) );
    }
    for ( std::size_t k=u.shape()[1]-3 ; k<u.shape()[1] ; ++k ) {
      fim12[i][k] = velocity(index_velocity<D>(i,k))*local_flux_m( u.field<_T,NumDims>::template stencil_border<D>(i,k+1) );
    }
  }

  return fim12;
}
/*
template < typename _T , std::size_t NumDims , direction::direction D >
auto
flux_m ( field<_T,NumDims> const& u , ublas::vector<_T> const& v )
{
  boost::multi_array<_T,NumDims> fim12(array_view<const std::size_t>(u.shape(),NumDims));

  for ( std::size_t i=0 ; i<u.shape()[0] ; ++i ) {
    for ( std::size_t k=0 ; k<u.shape()[1] ; ++k ) {
      fim12[i][k] = v(index_velocity<D>(i,k))*local_flux_m( u.field<_T,NumDims>::template stencil<D>(index_i_m<D>(i),index_k_m<D>(k)) );
    }
  }

  return fim12;
}
*/
template < typename _T , std::size_t NumDims >
auto
trp2D ( field<_T,NumDims> const& u , ublas::vector<_T> const& v , ublas::vector<_T> const& E /*, typename field<_T,NumDims>::step s*/ )
{
  field<_T,NumDims> trp(tools::array_view<const std::size_t>(u.shape(),NumDims));

  ublas::vector<_T> vm(v.size()),vp(v.size());
  ublas::vector<_T> Em(E.size()),Ep(E.size());

  for ( auto k=0 ; k<v.size() ; ++k) {
    vp(k) = std::max(v(k),0.);
    vm(k) = std::min(v(k),0.);
  }
  for ( auto i=0 ; i<E.size() ; ++i) {
    Ep(i) = std::max(E(i),0.);
    Em(i) = std::min(E(i),0.);
  }

  auto fp_x = flux_p<_T,NumDims,direction::x>(u,vp);
  auto fm_x = flux_m<_T,NumDims,direction::x>(u,vm);
  auto fp_v = flux_p<_T,NumDims,direction::v>(u,Ep);
  auto fm_v = flux_m<_T,NumDims,direction::v>(u,Em);


  for ( auto i=1 ; i<trp.size(0) ; ++i ) {
    for ( auto k=1 ; k<trp.size(1) ; ++k ) {
      trp[i][k] = (1./u.step.dx)*( fp_x[i][k]-fp_x[i-1][k] + fm_x[i][k]-fm_x[i-1][k] )
                + (1./u.step.dv)*( fp_v[i][k]-fp_v[i][k-1] + fm_v[i][k]-fm_v[i][k-1] );
    } { auto k=0,km1=trp.size(1)-1;
      trp[i][k] = (1./u.step.dx)*( fp_x[i][k]-fp_x[i-1][k] + fm_x[i][k]-fm_x[i-1][k] )
                + (1./u.step.dv)*( fp_v[i][k]-fp_v[i][km1] + fm_v[i][k]-fm_v[i][km1] );
    }
  }

  { auto i=0, im1=trp.size(0)-1;
    for ( auto k=1 ; k<trp.size(1) ; ++k ) {
      trp[i][k] = (1./u.step.dx)*( fp_x[i][k]-fp_x[im1][k] + fm_x[i][k]-fm_x[im1][k] )
                + (1./u.step.dv)*( fp_v[i][k]-fp_v[i][k-1] + fm_v[i][k]-fm_v[i][k-1] );
    } { auto k=0,km1=trp.size(1)-1;
      trp[i][k] = (1./u.step.dx)*( fp_x[i][k]-fp_x[im1][k] + fm_x[i][k]-fm_x[im1][k] )
                + (1./u.step.dv)*( fp_v[i][k]-fp_v[i][km1] + fm_v[i][k]-fm_v[i][km1] );
    }
  }

/*
  for ( auto i=0 ; i<trp.size(0)-1 ; ++i ) {
    for ( auto k=0 ; k<trp.size(1)-1 ; ++k ) {
      trp[i][k] = (1./u.step.dx)*( fp_x[i+1][k]-fp_x[i][k] + fm_x[i+1][k]-fm_x[i][k] )
               + (1./u.step.dv)*( fp_v[i][k+1]-fp_v[i][k] + fm_v[i][k+1]-fm_v[i][k] );
    } { auto k=trp.size(1)-1 , kp1=0;
      trp[i][k] = (1./u.step.dx)*( fp_x[i+1][k]-fp_x[i][k] + fm_x[i+1][k]-fm_x[i][k] )
               + (1./u.step.dv)*( fp_v[i][kp1]-fp_v[i][k-1] + fm_v[i][kp1]-fm_v[i][k-1] );
    }
  }

  { auto i=trp.size(0)-1, ip1=0;
    for ( auto k=0 ; k<trp.size(1)-1 ; ++k ) {
      trp[i][k] = (1./u.step.dx)*( fp_x[ip1][k]-fp_x[i][k] + fm_x[ip1][k]-fm_x[i][k] )
               + (1./u.step.dv)*( fp_v[i][k+1]-fp_v[i][k] + fm_v[i][k+1]-fm_v[i][k] );
    } { auto k=trp.size(1)-1 , kp1=0;
      trp[i][k] = (1./u.step.dx)*( fp_x[ip1][k]-fp_x[i][k] + fm_x[ip1][k]-fm_x[i][k] )
               + (1./u.step.dv)*( fp_v[i][kp1]-fp_v[i][k] + fm_v[i][kp1]-fm_v[i][k] );
    }
  }
*/

  return trp;
}
/*
template < typename _T , std::size_t NumDims >
auto
trp2D ( field<_T,NumDims> const& u , ublas::vector<_T> const& v , ublas::vector<_T> const& E )
{ return trp2D(u,v,E,u.step); }
*/
} // namespace weno

#endif
