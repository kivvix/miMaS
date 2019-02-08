#ifndef _WENO_H_
#define _WENO_H_

#include <algorithm>
#include <iostream>
#include <iterator>
#include <utility>
#include <type_traits>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multi_array.hpp>

#include <boost/iterator/zip_iterator.hpp>
#include "boost/tuple/tuple.hpp"

#include "field.h"
#include "array_view.h"

#include <boost/type_index.hpp>


namespace weno {

using namespace boost::numeric;
#define SQ(X) ((X)*(X))

/**
  @fn template < class IteratorTuple > auto local_flux ( boost::zip_iterator<IteratorTuple> & f_it )
  @brief Compute local flux $\left( f_{i,k+\frac{1}{2}}^2 , f_{i,k+\frac{1}{2}}^-\right)$
  @param `boost::zip_iterator<IteratorTuple> f_it` a `boost::zip_iterator` of length 6 on $(f_{i,k-2},\dots,f_{i,k+3})$ values
  @return a pair of flux plus and minus
**/
template < typename IteratorTuple >
auto
local_flux ( boost::zip_iterator<IteratorTuple> const& f_it )
{
  typedef typename std::remove_cv<typename std::remove_reference<decltype(f_it->template get<0>())>::type>::type _T; // récupération du type stocké dans le premier itérateur du zip_iterator (sans const ni reference)

  const static _T epsi = 1e-6;

  // NB pour optimisation minuscule : 2 termes sont calculés 2 fois :
  // * `(13./12.)*SQ( f_it->template get<1>() - 2.*f_it->template get<2>() + f_it->template get<3>() )` dans `w0p` et `w2m`
  // * `(13./12.)*SQ( f_it->template get<2>() - 2.*f_it->template get<3>() + f_it->template get<4>() )` dans `w2p` et `w1m`

  // $f_{i,k+1/2}^+$
  _T w0p = (13./12.)*SQ( f_it->template get<0>() - 2.*f_it->template get<1>() + f_it->template get<2>() ) + 0.25*SQ(    f_it->template get<0>() - 4.*f_it->template get<1>() + 3.*f_it->template get<2>() );
  _T w1p = (13./12.)*SQ( f_it->template get<1>() - 2.*f_it->template get<2>() + f_it->template get<3>() ) + 0.25*SQ(    f_it->template get<1>()                              -    f_it->template get<3>() );
  _T w2p = (13./12.)*SQ( f_it->template get<2>() - 2.*f_it->template get<3>() + f_it->template get<4>() ) + 0.25*SQ( 3.*f_it->template get<2>() - 4.*f_it->template get<3>() +    f_it->template get<4>() );

  w0p = 0.1/(SQ( epsi + w0p )); w1p = 0.6/(SQ( epsi + w1p )); w2p = 0.3/(SQ( epsi + w2p ));

  _T sum_w = w0p+w1p+w2p;
  w0p /= sum_w; w1p /= sum_w; w2p /= sum_w;

  _T fikp12p = w0p*( (2./6.)*f_it->template get<0>() - (7./6.)*f_it->template get<1>() + (11./6.)*f_it->template get<2>() )
             + w1p*(-(1./6.)*f_it->template get<1>() + (5./6.)*f_it->template get<2>() +  (2./6.)*f_it->template get<3>() )
             + w2p*( (2./6.)*f_it->template get<2>() + (5./6.)*f_it->template get<3>() -  (1./6.)*f_it->template get<4>() );

  // $f_{i,k+1/2}^-$
  _T w0m = (13./12.)*SQ( f_it->template get<3>() - 2.*f_it->template get<4>() + f_it->template get<5>() ) + 0.25*SQ( 3.*f_it->template get<3>() - 4.*f_it->template get<4>() +    f_it->template get<5>() );
  _T w1m = (13./12.)*SQ( f_it->template get<2>() - 2.*f_it->template get<3>() + f_it->template get<4>() ) + 0.25*SQ(    f_it->template get<2>()                              -    f_it->template get<4>() );
  _T w2m = (13./12.)*SQ( f_it->template get<1>() - 2.*f_it->template get<2>() + f_it->template get<3>() ) + 0.25*SQ(    f_it->template get<1>() - 4.*f_it->template get<2>() + 3.*f_it->template get<3>() );

  w0m = 0.1/SQ(w0m+epsi); w1m = 0.6/SQ(w1m+epsi); w2m = 0.3/SQ(w2m+epsi);

  sum_w = w0m+w1m+w2m;
  w0m /= sum_w; w1m /= sum_w; w2m /= sum_w;

  _T fikp12m = w2m*(-(1./6.)*f_it->template get<1>() + (5./6.)*f_it->template get<2>() + (2./6.)*f_it->template get<3>() )
             + w1m*( (2./6.)*f_it->template get<2>() + (5./6.)*f_it->template get<3>() - (1./6.)*f_it->template get<4>() )
             + w0m*((11./6.)*f_it->template get<3>() - (7./6.)*f_it->template get<4>() + (2./6.)*f_it->template get<5>() );
  
  // return a pair of fluxes, one loop on the data structure for two fluxes
  return std::make_pair(fikp12p,fikp12m);
}

#undef SQ


/**
  @fn template < typename _T , std::size_t NumDimsV > boost::multi_array<std::pair<_T,_T>,NumDimsV+1> flux ( flied<_T,NumDimsV> const& u )
  @brief Compute flux array $\left( f_{i,k+\frac{1}{2}}^2 , f_{i,k+\frac{1}{2}}^-\right)_{i,k}$
  @param `flied<_T,NumDimsV> const& u` the field to transport
  @return a `multi_array<std::pair<_T,_T>,NumDimsV+1>` of all flux (plus and minus)
**/
template < typename _T , std::size_t NumDimsV >
boost::multi_array<std::pair<_T,_T>,NumDimsV+1>
flux ( field<_T,NumDimsV> const& u )
{
  boost::multi_array<std::pair<_T,_T>,NumDimsV+1> fikp12(tools::array_view<const std::size_t>(u.shape(),NumDimsV));

  for ( std::size_t k=0 ; k<2 ; ++k ) {
    std::size_t i=0;
    for ( auto it = u.begin_border_stencil(k) ; it != u.end_border_stencil(k) ; ++it , ++i ) {
      fikp12[k][i] = local_flux(it);
    }
  }
  for ( std::size_t k=2 ; k<u.size(0)-3 ; ++k ) {
    std::size_t i=0;
    for ( auto it = u.begin_stencil(k) ; it != u.end_stencil(k) ; ++it , ++i ) {
      fikp12[k][i] = local_flux(it);
    }
  }
  for ( std::size_t k=u.size(0)-3 ; k<u.size(0) ; ++k ) {
    std::size_t i=0;
    for ( auto it = u.begin_border_stencil(k) ; it != u.end_border_stencil(k) ; ++it , ++i ) {
      fikp12[k][i] = local_flux(it);
    }
  }
  return fikp12;
}

/**
  @fn template < typename _T , std::size_t NumDimsV > auto trp_v ( field<_T,NumDimsV> const& u , ublas::vector<_T> const& E )
  @brief Transport a field `u` at velocity `E`
  @param `flied<_T,NumDimsV> const& u` the field to transport
  @param `ublas::vector<_T> const& E` the velocity in $v$-direction
  @return a `field<_T,NumDimsV>` a field with transported data at velocity `E` in $v$-direction
**/
template < typename _T , std::size_t NumDimsV >
auto
trp_v ( field<_T,NumDimsV> const& u , ublas::vector<_T> const& E )
{
  field<_T,NumDimsV> trp(tools::array_view<const std::size_t>(u.shape(),NumDimsV));

  ublas::vector<_T> Em(E.size()) , Ep(E.size());

  for ( auto i=0 ; i<E.size() ; ++i ) {
    Ep(i) = std::max(E(i),0.);
    Em(i) = std::min(E(i),0.);
  }

  auto fikp12 = flux(u);

  { auto k=0,km1=trp.size(0)-1;
    for ( auto i=0 ; i<trp.size(1) ; ++i ) {
      trp[k][i] = ( Ep(i)*(fikp12[k][i].first  - fikp12[km1][i].first)
                  + Em(i)*(fikp12[k][i].second - fikp12[km1][i].second) )/u.step.dv;
    }
  }
  for ( auto k=1 ; k<trp.size(0) ; ++k ) {
    for ( auto i=0 ; i<trp.size(1) ; ++i ) {
      trp[k][i] = ( Ep(i)*(fikp12[k][i].first  - fikp12[k-1][i].first)
                  + Em(i)*(fikp12[k][i].second - fikp12[k-1][i].second) )/u.step.dv;
    }
  }

  return trp;
}

} // namespace weno

#endif
