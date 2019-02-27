#ifndef _FTTRP_H_
#define _FTTRP_H_

#include <fftw3.h>

#include "array_view.h"

namespace fft {

enum complex_part { re=0, im=1 };

/*
  class: spectrum
  brief: la classe `spectrum` est un *wrapper* pour un signal réel 1D depuis `fftw`

  La classe contient seulement le tableau des coefficients de la FFT.
  On initialise le spectre en fonction de la taille des données en
  entrée, puis on appelle la méthode `spectrum::fft()` sur le signal
  dont on veut la transformée. On peut utiliser la structure `spectrum`
  comme un tableau (héritage de la structure `tools::array_view<T>`
  avec le type `T=fftw_complex` qui représente un complexe dans `fftw`.
  L'accès à ce complexe s'effectue par `spectre[i]` on accède à la
  partie réelle ou imaginaire à l'aide de l'énumération `complex_part`,
  `spectre[i][fft::re]` ou `spectre[i][fft:im]` (avec `fft::re=0` et
  `fft::im=1`).
  La gestion de la mémoire (allocation et libération) est transparente
  pour l'utilisateur, tout s'effectue dans le constructeur et
  destructeur.

  TODO: il semblerait en voyant le code de Lukas que std::complex<double> et fftw_complex se castent facilement ensemble, cela peut simplifier l'écriture des manipulations des coefficients de Fourier (pour Poisson et Lawson)
*/
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
    fftw_plan pI = fftw_plan_dft_c2r_1d(this->size(),(fftw_complex*)this->front(),it,FFTW_ESTIMATE);
    fftw_execute(pI); fftw_destroy_plan(pI);
    for ( std::size_t i=0 ; i<this->size() ; ++i,++it )
      { *it /= this->size(); }
  }
};

} // namespace fft

#endif
