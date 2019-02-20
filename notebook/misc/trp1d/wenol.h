#ifndef MISC_TRP1D_WENOL_H_
#define MISC_TRP1D_WENOL_H_

#include <utility>

namespace wenol
{

#define SQ(X) ((X)*(X))

template < typename _T >
std::pair<_T,_T>
local_flux ( _T uim2 , _T uim1 , _T ui , _T uip1 , _T uip2 , _T uip3 )
{
  _T w0p = 0.1, w1p = 0.6, w2p = 0.3;

  _T fip12p = w0p*( (2./6.)*uim2 - (7./6.)*uim1 + (11./6.)*ui   )
            + w1p*(-(1./6.)*uim1 + (5./6.)*ui   +  (2./6.)*uip1 )
            + w2p*( (2./6.)*ui   + (5./6.)*uip1 -  (1./6.)*uip2 );

  _T w0m = 0.1, w1m = 0.6, w2m = 0.3;

  _T fip12m = w2m*(-(1./6.)*uim1 + (5./6.)*ui   + (2./6.)*uip1 )
            + w1m*( (2./6.)*ui   + (5./6.)*uip1 - (1./6.)*uip2 )
            + w0m*((11./6.)*uip1 - (7./6.)*uip2 + (2./6.)*uip3 );

  return std::make_pair(fip12p,fip12m);
}

template <typename Container,typename _T>
Container
du ( Container const& u , _T const& v , _T dx )
{
  Container du(u.size());
  std::size_t N = u.size();

  std::pair<_T,_T> fim12,fip12;
  _T vp = std::max(v,0.);
  _T vm = std::min(v,0.);

  //i=0;
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

} // namespace wenol

#endif

