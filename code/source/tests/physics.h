#ifndef _PHYSICS_H_
#define _PHYSICS_H_

#include <algorithm>
#include <iostream>
#include <iterator>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>

#include "field.h"

using namespace boost::numeric;
namespace math = boost::math::constants;

template <typename _T>
struct variabili
{
  ublas::vector<_T> rho;
  ublas::vector<_T> u;
  ublas::vector<_T> T;

  variabili ( std::array<ublas::vector<_T>,3> && U )
    : rho(std::move(U[0])) , u(std::move(U[1])) , T(U[2].size())
  {
    for ( std::size_t i=0 ; i<rho.size() ; ++i ) {
      u(i) /= rho(i);
      T(i)  = (U[2](i) - U[1](i)*U[1](i)/(2.*rho(i)))/(0.5*rho(i));
    }
  }
};

template <typename _T,std::size_t NumDims>
struct maxwellian
{
  maxwellian ( typename field<_T,NumDims>::step s , typename field<_T,NumDims>::range r )
    : step(s) , range(r) , v((r.v_max-r.v_min)/s.dv ) ,
      data(boost::extents[(r.x_max-r.x_min)/s.dx][(r.v_max-r.v_min)/s.dv])
  {
    for ( auto i=0 ; i<v.size() ; ++i )
      { v(i) = i*step.dv + range.v_min; }
  }

  void
  update ( variabili<_T> const& var )
  { update(var.rho,var.u,var.T); }

  void
  update ( ublas::vector<_T> const& rho , ublas::vector<_T> const& u , ublas::vector<_T> const& T )
  {
    for ( std::size_t i=0 ; i<data.shape()[0] ; ++i ) {
      for ( std::size_t k=0 ; k<data.shape()[1] ; ++k ) {
        data[i][k] = rho(i)/std::sqrt(2.*math::pi<_T>*T(i))*std::exp(-0.5*(v(k)-u(i))*(v(k)-u(i))/T(i));
      }
    }
  }
  
  _T const&
  operator () ( std::size_t i , std::size_t k ) const
  { return data[i][k]; }
  _T &
  operator () ( std::size_t i , std::size_t k )
  { return data[i][k]; }

  typename field<_T,NumDims>::step step;
  typename field<_T,NumDims>::range range;
  ublas::vector<_T> v;
  boost::multi_array<_T,NumDims> data;
};

#endif
