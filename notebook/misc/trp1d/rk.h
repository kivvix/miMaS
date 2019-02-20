#ifndef MISC_TRP1D_RK_H_
#define MISC_TRP1D_RK_H_

#include <valarray>

namespace rk
{

template < typename _T , typename Func >
std::valarray<_T>
euler ( std::valarray<_T> const& un , Func L , _T v , _T dx , _T dt )
{
  return un + dt*L(un,v,dx);
}

template < typename _T , typename Func >
std::valarray<_T>
rk33 ( std::valarray<_T> const& un , Func L , _T v , _T dx , _T dt )
{
  std::valarray<_T> u1 = un + dt*L(un,v,dx);
  std::valarray<_T> u2 = 0.75*un + 0.25*u1 + 0.25*dt*L(u1,v,dx);
  return (1./3.)*un + (2./3.)*u2 + (2./3.)*dt*L(u2,v,dx);
}

template < typename _T , typename Func >
std::valarray<_T>
rk43 ( std::valarray<_T> const& un , Func L , _T v , _T dx , _T dt )
{
  std::valarray<_T> u1 = un + 0.5*dt*L(un,v,dx);
  std::valarray<_T> u2 = u1 + 0.5*dt*L(u1,v,dx);
  std::valarray<_T> u3 = (2./3.)*un + (1./3.)*u2 + (1./6.)*dt*L(u2,v,dx);
  return u3 + 0.5*dt*L(u3,v,dx);
}

} // namespace rk

#endif

