#include <algorithm>
#include <iterator>
#include <cmath>
#include <valarray>
#include <sstream>
#include <complex>
#include <tuple>
#include <numeric>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "miMaS/field.h"
#include "miMaS/complex_field.h"
#include "miMaS/lagrange5.h"
#include "miMaS/fft.h"
#include "miMaS/array_view.h"
#include "miMaS/poisson.h"

using namespace boost::numeric;

/*
template < typename _T , std::size_t NumDimsV >
struct U_type {
  std::tuple< ublas::vector<_T> , ublas::vector<_T> , complex_field<_T,NumDimsV> > data;

  //U_type ( std::tuple< ublas::vector<_T> , ublas::vector<_T> , complex_field<_T,NumDimsV> > & t )
  //  : data(t)
  //{ ; }

  const auto &
  uc () const
  { return std::get<0>(data); }
  auto &
  uc ()
  { return std::get<0>(data); }

  const auto &
  E () const
  { return std::get<1>(data); }
  auto &
  E ()
  { return std::get<1>(data); }

  const auto &
  hfh () const
  { return std::get<2>(data); }
  auto &
  hfh ()
  { return std::get<2>(data); }
};
*/

template < typename _T , std::size_t NumDimsV >
struct splitting
{
  //typedef std::tuple< ublas::vector<_T>,ublas::vector<_T>,field<_T,NumDimsV> > U_type;
  typedef void U_type;

  _T _dx,_dv;
  _T _v_min , _x_max;
  _T _rho_c;
  field<_T,NumDimsV> _tmpf0;
  field<_T,NumDimsV> _tmpf1;
  std::size_t _Nx,_Nv;
  ublas::vector<_T> kx;

  splitting ( field<_T,NumDimsV> const & fh , _T l , _T rho_c) :
    _dx(fh.step.dx) , _dv(fh.step.dv) ,
    _v_min(fh.range.v_min) , _x_max(fh.range.x_max) ,
    _rho_c(rho_c) ,
    _tmpf0(tools::array_view<const std::size_t>(fh.shape(),2)) ,
    _tmpf1(tools::array_view<const std::size_t>(fh.shape(),2)) ,
    _Nx(fh.shape()[1]) ,
    _Nv(fh.shape()[0]) ,
    kx(fh.shape()[1],0.)
  {
    _tmpf0.step = fh.step; _tmpf0.range = fh.range;
    _tmpf1.step = fh.step; _tmpf1.range = fh.range;
    //kx[0] = 1.;
    /*
    for ( auto i=1 ; i<_Nx/2 ; ++i ) { kx[i] = 2.*math::pi<double>()*i/l; }
    for ( int i=-_Nx/2 ; i<0 ; ++i ) { kx[_Nx+i] = 2.*math::pi<double>()*i/l; }
    */
    for ( auto i=0 ; i<_Nx/2 ; ++i ) { kx[i]     = 2.*math::pi<double>()*i/l; }
    for ( int i=-_Nx/2 ; i<0 ; ++i ) { kx[_Nx+i] = 2.*math::pi<double>()*i/l; }
  }

  U_type
  phi_a ( _T dt , ublas::vector<_T> & uc , ublas::vector<_T> & E , complex_field<_T,NumDimsV> & hfh )
  {
    // equivalent of hf of Nicolas
    const std::complex<double> & I = std::complex<double>(0.,1.);

    // compute hdiffrho and update hfh
    ublas::vector<std::complex<_T>> hdiffrho(_Nx,0.);
    for ( auto k=0 ; k<_Nv ; ++k ) {
      double vk = k*_dv + _v_min;
      for ( auto i=0 ; i<_Nx ; ++i ) {
        // hrho_n
        hdiffrho[i] += hfh[k][i]*_dv;
        // update hfh
        hfh[k][i] *= std::exp(-I*vk*kx[i]*dt);
        // hrho_{n+1}
        hdiffrho[i] -= hfh[k][i]*_dv;
      }
    }

    // update E
    fft::spectrum_ hE(_Nx); hE.fft(E.begin());
    hE[0] = 0.;
    for ( auto i=1 ; i<_Nx ; ++i ) {
      hE[i] += I/kx[i]*hdiffrho[i];
    }
    hE.ifft(E.begin());
    

    /*
    for ( auto k=0 ; k < _Nv ; ++k )
      { fft::ifft( hfh[k].begin() , hfh[k].end() , _tmpf0[k].begin() ); }
    ublas::vector<_T> diffrho = _tmpf0.density();

    // update hfh
    for ( auto k=0 ; k<_Nv ; ++k ) {
      double vk = k*_dv + _v_min;
      for ( auto i=0 ; i<_Nx ; ++i ) {
        hfh[k][i] *= std::exp(-I*vk*kx[i]*dt);
      }
    }

    for ( auto k=0 ; k < _Nv ; ++k )
      { fft::ifft( hfh[k].begin() , hfh[k].end() , _tmpf0[k].begin() ); }
    diffrho -= _tmpf0.density();

    fft::spectrum_ hdiffrho(_Nx); hdiffrho.fft(diffrho.begin());

    // update E
    fft::spectrum_ hE(_Nx); hE.fft(E.begin());
    hE[0] = 0.;
    for ( auto i=1 ; i<_Nx ; ++i ) {
      hE[i] += I/kx[i]*hdiffrho[i];
    }
    hE.ifft(E.begin());
    */
  }

  U_type
  phi_b ( _T dt , ublas::vector<_T> & uc , ublas::vector<_T> const & E , complex_field<_T,NumDimsV> & hfh )
  {
    // equivalent of HE of Nicolas

    for ( auto k=0 ; k < _Nv ; ++k )
      { fft::ifft( &(hfh[k][0]) , &(hfh[k][0])+_Nx , &(_tmpf0[k][0]) ); }

    // faire des trucs sur `_tmpf0` : $f(x,v) = f(x,v-dt*E)$
    for ( auto k=0 ; k<_Nv ; ++k )
    {
      for ( auto i=0 ; i<_Nx ; ++i ) {
        _T  vstar = (k*_dv + _v_min) - dt*E[i]; // vstar = v - dt*E
        int kstar = std::ceil((vstar - _v_min)/_dv);
        //std::cout << i <<","<<k << " > " << kstar << " : " << (kstar-3+_Nv)%_Nv << " , " << (kstar-2+_Nv)%_Nv << " , " << (kstar-1+_Nv)%_Nv << " , " << (kstar+_Nv)%_Nv << " , " << (kstar+1+_Nv)%_Nv << " , " << (kstar+2+_Nv)%_Nv << "\n";

        auto N = lagrange5::generator(
                    _tmpf0[(kstar-3+_Nv)%_Nv][i],_tmpf0[(kstar-2+_Nv)%_Nv][i],_tmpf0[(kstar-1+_Nv)%_Nv][i],_tmpf0[(kstar+_Nv)%_Nv][i],_tmpf0[(kstar+1+_Nv)%_Nv][i],_tmpf0[(kstar+2+_Nv)%_Nv][i] ,
                    _dv , kstar*_dv + _v_min
                  );
        _tmpf1[k][i] = N(vstar);
      }
    }
    // update of hfh
    for ( auto k=0 ; k < _Nv ; ++k )
      { fft::fft( &(_tmpf1[k][0]) , &(_tmpf1[k][0])+_Nx , &(hfh[k][0]) ); }

    // update uc
    uc += dt*E;

    //return {uc,E,hfh};
  }

  U_type
  phi_c ( _T dt , ublas::vector<_T> & uc , ublas::vector<_T> & E , complex_field<_T,NumDimsV> & hfh )
  {
    // equivalent of Hu of Nicolas
    _T rho_c=_rho_c;
    _T curtot = std::accumulate(
        uc.begin() , uc.end() ,
        0. ,
        [rho_c]( _T s , _T ui ) { return s + rho_c*ui; }
      ) * _dx / _x_max;

    // update E
    std::transform(
        E.begin()  , E.end() ,
        uc.begin() ,
        E.begin()  ,
        [&]( _T ei , _T ui ) { return ei - dt*(_rho_c*ui - curtot); }
      );

    //return {uc,E,hfh};
  }
};

