#ifndef _FTTRP_H_
#define _FTTRP_H_

#include <algorithm>
#include <iostream>
#include <iterator>
#include <vector>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/multi_array.hpp>

#include <fftw3.h>

#include "direction.h"
#include "field.h"
#include "array_view.h"

#define REAL 0
#define IMAG 1

using namespace boost::numeric;
namespace math = boost::math::constants;

namespace fft {

enum complex_part { re=0, im=1 };

struct spectrum
  : public tools::array_view<fftw_complex>
{
  spectrum ( std::size_t n )
    : tools::array_view<fftw_complex>( (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*n) , n )
  {}

  ~spectrum ()
  { fftw_free((fftw_complex*)this->data()); }

  template < typename Iterator >
  void
  fft ( Iterator it )
  {
    fftw_plan p = fftw_plan_dft_r2c_1d(this->size(),it,(fftw_complex*)this->front(),FFTW_ESTIMATE);
    fftw_execute(p); fftw_destroy_plan(p);
  }

  template < typename Iterator >
  void
  ifft ( Iterator it )
  {
    fftw_plan pI = fftw_plan_dft_c2r_1d(this->size(),(fftw_complex*)this->front(),it,FFTW_PRESERVE_INPUT);
    fftw_execute(pI); fftw_destroy_plan(pI);
    for ( std::size_t i=0 ; i<this->size() ; ++i,++it )
      { *it /= this->size(); }
  }
};


template < typename _T , std::size_t NumDimsV >
struct fft_x
{
  boost::multi_array<fftw_complex * , NumDimsV> spectrum;
  ublas::vector<_T> kx;
  _T l;
  std::size_t nx;

  template < typename ExtentList >
  fft_x ( ExtentList const& sizes , std::size_t const Nx , _T L )
    : spectrum(sizes) , kx(Nx) , l(L) , nx(Nx)
  {
    for ( auto k=0 ; k<spectrum.shape()[0] ; ++k ) {
      spectrum[k] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(Nx));
    }

    //std::generate(kx.begin(),kx.end(),[](){return 42.;});
    for ( auto i=0 ; i<Nx/2+1 ; ++i )   { kx[i] = i/l; }
    for ( auto i=0 ; i<((Nx/2)) ; ++i ) { kx[i+Nx/2+1] = -kx[Nx/2-i]; }
    /*kx[Nx/2] = 0.;
    for ( auto i=0 ; i<((Nx/2)-1) ; ++i )
      { kx[i+1+Nx/2] = -kx[Nx/2-i-1]; }*/
    //for ( auto i=0 ; i<Nx ; ++i ) { std::cout << kx[i] << " "; } std::cout << std::endl;
  }

  ~fft_x ()
  {
    for ( auto k=0 ; k<spectrum.shape()[0] ; ++k )
      { fftw_free(spectrum[k]); }
  }
/*
  auto
  fft ()
  {
    boost::multi_array<fftw_complex * , NumDimsV> s;
    for ( auto k=0 ; k<spectrum.shape()[0] ; ++k ) {
      s[k] = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(nx));
    }
    fftw_plan p = fftw_plan_dft_r2c_1d(kx.size(),&uk[0],s[k],FFTW_ESTIMATE);
    fftw_execute(p);
    return s;
  }

  auto
  ifft ()
  {}
*/
  auto
  operator () ( field<_T,NumDimsV> const& u , ublas::vector<_T> const& v , _T dt )
  {
    field<_T,NumDimsV> trp(tools::array_view<const std::size_t>(u.shape(),NumDimsV),u.step,u.range);

    for ( auto k=0 ; k<u.size(0) ; ++k ) {
      ublas::vector<_T> uk(kx.size()); std::copy(u[k].begin(),u[k].end(),uk.begin());
      fftw_plan p = fftw_plan_dft_r2c_1d(kx.size(),&uk[0],spectrum[k],FFTW_ESTIMATE);
      fftw_execute(p);

      // do something with spectrum
      
      for ( auto i=0 ; i<kx.size() ; ++i ) {
        auto re = spectrum[k][i][REAL], im = spectrum[k][i][IMAG];
        spectrum[k][i][REAL] = std::cos(kx[i]*v[k]*dt)*re + std::sin(kx[i]*v[k]*dt)*im; // Re(spectrum[k][i])
        spectrum[k][i][IMAG] = std::cos(kx[i]*v[k]*dt)*im - std::sin(kx[i]*v[k]*dt)*re; // Im(spectrum[k][i])
      }

      fftw_plan pI = fftw_plan_dft_c2r_1d(kx.size(),spectrum[k],&trp[k][0],FFTW_PRESERVE_INPUT);
      fftw_execute(pI);

      for ( auto i=0 ; i<trp[k].size() ; ++i ) { // t=trp[k].begin() ; it!=trp[k].end() ; ++it ) {
        //*it /= l;
        //*it /= trp[k].size();
        trp[k][i] = trp[k][i]/trp[k].size();
      }
      //std::cout << trp[k].size() << " " << kx.size() << " " << nx << " " << u.size(NumDimsV) << " " << spectrum.shape()[1] << "\n";

      fftw_destroy_plan(p); fftw_destroy_plan(pI);
    }

    return trp;
  }
};

} // namespace fft

#endif
