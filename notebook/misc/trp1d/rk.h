#ifndef MISC_TRP1D_RK_H_
#define MISC_TRP1D_RK_H_

#include <valarray>

/*
NOTE: all cast by `_T{}` have done to work also with `_T = float` or `_T = long double` (so I'm sorry for the readability)
*/

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
  std::valarray<_T> u2 = _T{0.75}*un + _T{0.25}*u1 + _T{0.25}*dt*L(u1,v,dx);
  return (_T{1.}/_T{3.})*un + (_T{2.}/_T{3.})*u2 + (_T{2.}/_T{3.})*dt*L(u2,v,dx);
}

template < typename _T , typename Func >
std::valarray<_T>
rk43 ( std::valarray<_T> const& un , Func L , _T v , _T dx , _T dt )
{
  std::valarray<_T> u1 = un + _T{0.5}*dt*L(un,v,dx);
  std::valarray<_T> u2 = u1 + _T{0.5}*dt*L(u1,v,dx);
  std::valarray<_T> u3 = (_T{2.}/_T{3.})*un + (_T{1.}/_T{3.})*u2 + (_T{1.}/_T{6.})*dt*L(u2,v,dx);
  return u3 + _T{0.5}*dt*L(u3,v,dx);
}

template < typename _T , typename Func >
std::valarray<_T>
rk76 ( std::valarray<_T> const& un , Func L , _T v , _T dx , _T dt )
{
  _T s21 = std::sqrt(_T{21.}) , nu=_T{0.5};
  std::valarray<_T> k1 = dt*L(un,v,dx);
  std::valarray<_T> k2 = dt*L(un+nu*k1,v,dx);
  std::valarray<_T> k3 = dt*L(un+((_T{4.}*nu-_T{1.})*k1+k2)/(_T{8.}*nu),v,dx);
  std::valarray<_T> k4 = dt*L(un+((_T{10.}*nu-_T{2.})*k1+_T{2.}*k2+_T{8.}*nu*k3)/(_T{27.}*nu),v,dx);
  std::valarray<_T> k5 = dt*L(un+(-((_T{77.}*nu-_T{56.})+(_T{17.}*nu-_T{8.})*s21)*k1 - _T{8.}*(_T{7.}+s21)*k2 + _T{48.}*(_T{7.}+s21)*nu*k3 -_T{3.}*(_T{21.}+s21)*nu*k4)/(_T{392.}*nu),v,dx);
  std::valarray<_T> k6 = dt*L(un+(-_T{5.}*((_T{287.}*nu-_T{56.})-(_T{59.}*nu-_T{8.})*s21)*k1 - _T{40.}*(_T{7.}-s21)*k2 + _T{320.}*s21*nu*k3 + _T{3.}*(_T{21.}-_T{121.}*s21)*nu*k4 + _T{392.}*(_T{6.}-s21)*nu*k5)/(_T{1960.}*nu),v,dx);
  std::valarray<_T> k7 = dt*L(un+(_T{15.}*((_T{30.}*nu-_T{8.})-_T{7.}*nu*s21)*k1 + _T{120.}*k2 - _T{40.}*(_T{5.}+_T{7.}*s21)*nu*k3 + _T{63.}*(_T{2.}+_T{3.}*s21)*nu*k4 - _T{14.}*(_T{49.}-_T{9.}*s21)*nu*k5 + _T{70.}*(_T{7.}+s21)*nu*k6)/(_T{180.}*nu),v,dx);
  return un + (_T{9.}*k1+_T{64.}*k3+_T{49.}*k5+_T{49.}*k6+_T{9.}*k7)/_T{180.};
}

template < typename _T , typename Func >
std::valarray<_T>
rk44 ( std::valarray<_T> const& un , Func L , _T v , _T dx , _T dt )
{
  std::valarray<_T> u1 = un + _T{0.5}*dt*L(un,v,dx);
  std::valarray<_T> u2 = un + _T{0.5}*dt*L(u1,v,dx);
  std::valarray<_T> u3 = un + dt*L(u2,v,dx);
  return un + dt*((_T{1.}/_T{6.})*L(un,v,dx) + (_T{1.}/_T{3.})*L(u1,v,dx) + (_T{1.}/_T{3.})*L(u2,v,dx) + (_T{1.}/_T{6.})*L(u3,v,dx));
}

template < typename _T , typename Func >
std::valarray<_T>
rk53 ( std::valarray<_T> const& un , Func L , _T v , _T dx , _T dt )
{
  std::valarray<_T> u1 = un +  (_T{1.}/_T{7.})*dt*L(un,v,dx);
  std::valarray<_T> u2 = un + (_T{3.}/_T{16.})*dt*L(u1,v,dx);
  std::valarray<_T> u3 = un +  (_T{1.}/_T{3.})*dt*L(u2,v,dx);
  std::valarray<_T> u4 = un +  (_T{2.}/_T{3.})*dt*L(u3,v,dx);
  return un + _T{0.25}*dt*L(un,v,dx) + _T{0.75}*dt*L(u4,v,dx);
}

template < typename _T , typename Func >
std::valarray<_T>
rk86 ( std::valarray<_T> const& un , Func L , _T v , _T dx , _T dt )
{
  std::valarray<_T> u1 = un +             (_T{1.}/_T{9.})*dt*L(un,v,dx);
  std::valarray<_T> u2 = un +            (_T{1.}/_T{24.})*dt*L(un,v,dx) +          (_T{1.}/_T{8.})*dt*L(u1,v,dx);
  std::valarray<_T> u3 = un +             (_T{1.}/_T{6.})*dt*L(un,v,dx) -          (_T{1.}/_T{2.})*dt*L(u1,v,dx) +           (_T{2.}/_T{3.})*dt*L(u2,v,dx);
  std::valarray<_T> u4 = un +        (_T{935.}/_T{2536.})*dt*L(un,v,dx) -    (_T{2781.}/_T{2536.})*dt*L(u1,v,dx) +       (_T{309.}/_T{317.})*dt*L(u2,v,dx) +       (_T{321.}/_T{1268.})*dt*L(u3,v,dx);
  std::valarray<_T> u5 = un -       (_T{12710.}/_T{951.})*dt*L(un,v,dx) +     (_T{8287.}/_T{317.})*dt*L(u1,v,dx) -        (_T{40.}/_T{317.})*dt*L(u2,v,dx) -       (_T{6335.}/_T{317.})*dt*L(u3,v,dx) +                _T{8.}*dt*L(u4,v,dx);
  std::valarray<_T> u6 = un + (_T{5840285.}/_T{3104064.})*dt*L(un,v,dx) -    (_T{7019.}/_T{2536.})*dt*L(u1,v,dx) -   (_T{52213.}/_T{86224.})*dt*L(u2,v,dx) + (_T{1278709.}/_T{517344.})*dt*L(u3,v,dx) -  (_T{433.}/_T{2448.})*dt*L(u4,v,dx) +  (_T{33.}/_T{1088.})*dt*L(u5,v,dx);
  std::valarray<_T> u7 = un - (_T{5101675.}/_T{1767592.})*dt*L(un,v,dx) + (_T{112077.}/_T{25994.})*dt*L(u1,v,dx) + (_T{334875.}/_T{441898.})*dt*L(u2,v,dx) -  (_T{973617.}/_T{883796.})*dt*L(u3,v,dx) - (_T{1421.}/_T{1394.})*dt*L(u4,v,dx) + (_T{333.}/_T{5576.})*dt*L(u5,v,dx) + (_T{36.}/_T{41.})*dt*L(u6,v,dx);
  return un + (_T{41.}/_T{840.})*dt*L(un,v,dx)  + (_T{9.}/_T{35.})*dt*L(u2,v,dx) + (_T{9.}/_T{280.})*dt*L(u3,v,dx) + (_T{34.}/_T{105.})*dt*L(u4,v,dx) + (_T{9.}/_T{280.})*dt*L(u5,v,dx) + (_T{9.}/_T{35.})*dt*L(u6,v,dx) + (_T{41.}/_T{840.})*dt*L(u7,v,dx) ;
}


} // namespace rk

#endif

