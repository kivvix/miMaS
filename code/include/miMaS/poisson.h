#ifndef _POISSON_H_
#define _POISSON_H_

#include <algorithm>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include <fftw3.h>

#include "fft.h"

namespace math = boost::math::constants;

template < typename _T >
struct poisson
{
  fft::spectrum     hE;
  ublas::vector<_T> kx;

  poisson  ( std::size_t nx , _T l )
    : hE(nx) , kx(nx)
  {
    /* TODO: trouver un meilleur truc qu'avoir les modes stockés dans la classe de poisson */
    //for ( auto i=0 ; i<nx/2+1 ; ++i )   { kx[i] = 2.*math::pi<double>()*i/l; }
    //for ( auto i=0 ; i<((nx/2)) ; ++i ) { kx[i+nx/2+1] = -kx[nx/2-i]; }
    for ( auto i=0 ; i<nx/2 ; ++i ) { kx[i] = 2.*math::pi<double>()*i/l; }
    for ( int i=-nx/2 ; i<0 ; ++i ) { kx[nx+i] = 2.*math::pi<double>()*i/l; }
  }

  ~poisson ()
  {}

  ublas::vector<_T>
  operator () ( ublas::vector<_T> const& rho )
  {
    ublas::vector<_T> E(rho.size());
    std::transform(rho.begin(),rho.end(),E.begin(),[]( _T const& r ){ return r-1.; } );

    hE.fft(&E[0]);
    hE[0][fft::re] = 0.; hE[0][fft::im] = 0.;
    for ( std::size_t i=1 ; i<hE.size() ; ++i ) {
      _T re = hE[i][fft::re]; _T im = hE[i][fft::im];
      hE[i][fft::re] =  im/kx[i];
      hE[i][fft::im] = -re/kx[i];
    }
    hE.ifft(&E[0]);
    return E;
  }
};

#endif

