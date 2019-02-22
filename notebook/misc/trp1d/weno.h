#ifndef MISC_TRP1D_WENO_H_
#define MISC_TRP1D_WENO_H_

#include <utility>

namespace weno
{

#define SQ(X) ((X)*(X))

template < typename _T >
std::pair<_T,_T>
local_flux ( _T uim2 , _T uim1 , _T ui , _T uip1 , _T uip2 , _T uip3 )
{
  _T w0p = (_T{13.}/_T{12.})*SQ(uim2 - _T{2.}*uim1 + ui  ) + _T{0.25}*SQ(       uim2 - _T{4.}*uim1 + _T{3.}*ui   );
  _T w1p = (_T{13.}/_T{12.})*SQ(uim1 - _T{2.}*ui   + uip1) + _T{0.25}*SQ(       uim1               -        uip1 );
  _T w2p = (_T{13.}/_T{12.})*SQ(ui   - _T{2.}*uip1 + uip2) + _T{0.25}*SQ(_T{3.}*ui   - _T{4.}*uip1 +        uip2 );

  w0p = _T{0.1}/SQ(_T{1e-6}+w0p); w1p = _T{0.6}/SQ(_T{1e-6}+w1p); w2p = _T{0.3}/SQ(_T{1e-6}+w2p);
  _T sum_wp = w0p+w1p+w2p;
  w0p /= sum_wp; w1p /= sum_wp; w2p /= sum_wp;

  _T fip12p = w0p*( (_T{2.}/_T{6.})*uim2 - (_T{7.}/_T{6.})*uim1 + (_T{11.}/_T{6.})*ui   )
            + w1p*(-(_T{1.}/_T{6.})*uim1 + (_T{5.}/_T{6.})*ui   +  (_T{2.}/_T{6.})*uip1 )
            + w2p*( (_T{2.}/_T{6.})*ui   + (_T{5.}/_T{6.})*uip1 -  (_T{1.}/_T{6.})*uip2 );

  _T w0m = (_T{13.}/_T{12.})*SQ(uip1 - _T{2.}*uip2 + uip3) + _T{0.25}*SQ(_T{3.}*uip1 - _T{4.}*uip2 +        uip3 );
  _T w1m = (_T{13.}/_T{12.})*SQ(ui   - _T{2.}*uip1 + uip2) + _T{0.25}*SQ(       ui                 -        uip2 );
  _T w2m = (_T{13.}/_T{12.})*SQ(uim1 - _T{2.}*ui   + uip1) + _T{0.25}*SQ(       uim1 - _T{4.}*ui   + _T{3.}*uip1 );

  w0m = _T{0.1}/SQ(_T{1e-6}+w0m); w1m = _T{0.6}/SQ(_T{1e-6}+w1m); w2m = _T{0.3}/SQ(_T{1e-6}+w2m);
  _T sum_wm = w0m+w1m+w2m;
  w0m /= sum_wm; w1m /= sum_wm; w2m /= sum_wm;

  _T fip12m = w2m*(-(_T{1.}/_T{6.})*uim1 + (_T{5.}/_T{6.})*ui   + (_T{2.}/_T{6.})*uip1 )
            + w1m*( (_T{2.}/_T{6.})*ui   + (_T{5.}/_T{6.})*uip1 - (_T{1.}/_T{6.})*uip2 )
            + w0m*((_T{11.}/_T{6.})*uip1 - (_T{7.}/_T{6.})*uip2 + (_T{2.}/_T{6.})*uip3 );

  return std::make_pair(fip12p,fip12m);
}

template <typename Container,typename _T>
Container
du ( Container const& u , _T const& v , _T dx )
{
  Container du(u.size());
  std::size_t N = u.size();

  std::pair<_T,_T> fim12,fip12;
  _T vp = std::max(v,_T{0.});
  _T vm = std::min(v,_T{0.});

  // i=0;
  fim12 = local_flux( u[(-3+N)%N],u[(-2+N)%N],u[(-1+N)%N],u[0],u[1],u[2] );
  for ( std::size_t i=0 ; i<3 ; ++i ) {
    fip12 = local_flux( u[(i-2+N)%N],u[(i-1+N)%N],u[i],u[i+1],u[i+2],u[i+3] );
    du[i] = ( vp*(fip12.first-fim12.first) + vm*(fip12.second-fim12.second) )/dx;
    fim12 = fip12;
  }
  for ( std::size_t i=3 ; i<u.size()-3 ; ++i ) {
    fip12 = local_flux( u[i-2],u[i-1],u[i],u[i+1],u[i+2],u[i+3] );
    du[i] = ( vp*(fip12.first-fim12.first) + vm*(fip12.second-fim12.second) )/dx;
    fim12 = fip12;
  }
  for ( std::size_t i=u.size()-3 ; i<u.size() ; ++i ) {
    fip12 = local_flux( u[i-2],u[i-1],u[i],u[(i+1)%N],u[(i+2)%N],u[(i+3)%N] );
    du[i] = ( vp*(fip12.first-fim12.first) + vm*(fip12.second-fim12.second) )/dx;
    fim12 = fip12;
  }

  return du;
}

} // namespace weno

#endif

