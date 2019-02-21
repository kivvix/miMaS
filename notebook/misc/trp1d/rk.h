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

template < typename _T , typename Func >
std::valarray<_T>
rk76 ( std::valarray<_T> const& un , Func L , _T v , _T dx , _T dt )
{
  _T s21 = std::sqrt(21.) , nu=0.5;
  std::valarray<_T> k1 = dt*L(un,v,dx);
  std::valarray<_T> k2 = dt*L(un+nu*k1,v,dx);
  std::valarray<_T> k3 = dt*L(un+((4.*nu-1.)*k1+k2)/(8.*nu),v,dx);
  std::valarray<_T> k4 = dt*L(un+((10.*nu-2.)*k1+2.*k2+8.*nu*k3)/(27.*nu),v,dx);
  std::valarray<_T> k5 = dt*L(un+(-((77.*nu-56.)+(17.*nu-8.)*s21)*k1 - 8.*(7.+s21)*k2 + 48.*(7.+s21)*nu*k3 -3.*(21.+s21)*nu*k4)/(392.*nu),v,dx);
  std::valarray<_T> k6 = dt*L(un+(-5.*((287*nu-56.)-(59.*nu-8.)*s21)*k1 - 40.*(7.-s21)*k2 + 320.*s21*nu*k3 + 3.*(21.-121.*s21)*nu*k4 + 392.*(6.-s21)*nu*k5)/(1960.*nu),v,dx);
  std::valarray<_T> k7 = dt*L(un+(15.*((30.*nu-8.)-7.*nu*s21)*k1 + 120.*k2 - 40.*(5.+7.*s21)*nu*k3 + 63.*(2.+3.*s21)*nu*k4 - 14.*(49.-9.*s21)*nu*k5 + 70.*(7.+s21)*nu*k6)/(180.*nu),v,dx);
  return un + (9.*k1+64.*k3+49.*k5+49.*k6+9.*k7)/180.;
}

template < typename _T , typename Func >
std::valarray<_T>
rk44 ( std::valarray<_T> const& un , Func L , _T v , _T dx , _T dt )
{
  std::valarray<_T> u1 = un + 0.5*dt*L(un,v,dx);
  std::valarray<_T> u2 = un + 0.5*dt*L(u1,v,dx);
  std::valarray<_T> u3 = un + dt*L(u2,v,dx);
  return un + dt*((1./6.)*L(un,v,dx) + (1./3.)*L(u1,v,dx) + (1./3.)*L(u2,v,dx) + (1./6.)*L(u3,v,dx));
}

template < typename _T , typename Func >
std::valarray<_T>
rk53 ( std::valarray<_T> const& un , Func L , _T v , _T dx , _T dt )
{
  std::valarray<_T> u1 = un + (1./7.)*dt*L(un,v,dx);
  std::valarray<_T> u2 = un + (3./16.)*dt*L(u1,v,dx);
  std::valarray<_T> u3 = un + (1./3.)*dt*L(u2,v,dx);
  std::valarray<_T> u4 = un + (2./3.)*dt*L(u3,v,dx);
  return un + 0.25*dt*L(un,v,dx) + 0.75*dt*L(u4,v,dx);
}

template < typename _T , typename Func >
std::valarray<_T>
rk86 ( std::valarray<_T> const& un , Func L , _T v , _T dx , _T dt )
{//(./.)*dt*L(u,v,dx)
  std::valarray<_T> u1 = un +             (1./9.)*dt*L(un,v,dx);
  std::valarray<_T> u2 = un +            (1./24.)*dt*L(un,v,dx) +          (1./8.)*dt*L(u1,v,dx);
  std::valarray<_T> u3 = un +             (1./6.)*dt*L(un,v,dx) -          (1./2.)*dt*L(u1,v,dx) +           (2./3.)*dt*L(u2,v,dx);
  std::valarray<_T> u4 = un +        (935./2536.)*dt*L(un,v,dx) -    (2781./2536.)*dt*L(u1,v,dx) +       (309./317.)*dt*L(u2,v,dx) +       (321./1268.)*dt*L(u3,v,dx);
  std::valarray<_T> u5 = un -       (12710./951.)*dt*L(un,v,dx) +     (8287./317.)*dt*L(u1,v,dx) -        (40./317.)*dt*L(u2,v,dx) -       (6335./317.)*dt*L(u3,v,dx) +            8.*dt*L(u4,v,dx);
  std::valarray<_T> u6 = un + (5840285./3104064.)*dt*L(un,v,dx) -    (7019./2536.)*dt*L(u1,v,dx) -   (52213./86224.)*dt*L(u2,v,dx) + (1278709./517344.)*dt*L(u3,v,dx) -  (433./2448.)*dt*L(u4,v,dx) +  (33./1088.)*dt*L(u5,v,dx);
  std::valarray<_T> u7 = un - (5101675./1767592.)*dt*L(un,v,dx) + (112077./25994.)*dt*L(u1,v,dx) + (334875./441898.)*dt*L(u2,v,dx) -  (973617./883796.)*dt*L(u3,v,dx) - (1421./1394.)*dt*L(u4,v,dx) + (333./5576.)*dt*L(u5,v,dx) + (36./41.)*dt*L(u6,v,dx);
  return un + (41./840.)*dt*L(un,v,dx)  + (9./35.)*dt*L(u2,v,dx) + (9./280.)*dt*L(u3,v,dx) + (34./105.)*dt*L(u4,v,dx) + (9./280.)*dt*L(u5,v,dx) + (9./35.)*dt*L(u6,v,dx) + (41./840.)*dt*L(u7,v,dx) ;
}


} // namespace rk

#endif

