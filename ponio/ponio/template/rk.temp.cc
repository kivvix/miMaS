#ifndef MISC_TRP1D_RK_H_
#define MISC_TRP1D_RK_H_

#include <valarray>

namespace rk
{

/**---
  type: `LType<_T>`
  brief: type of a valid function `L` in $\dot{u} = L(u,t)$
  details: This type correspond to `L(u,t)` with `u` (`std::valarray<_T>`) input function and `t` (_T) current time.
---**/
template < typename _T >
using LType = std::function< std::valarray<_T> ( std::valarray<_T> const& , _T ) >;

/**---
  type: `RKType<_T>`
  brief: type of a valid function to represent a $RK(s,n)$ method to solve $\dot{u} = L(u,t)$
---**/
template < typename _T >
using RKType = std::function< std::valarray<_T> ( std::valarray<_T> const& , LType<_T> , _T , _T ) >;

/**
All flowing function of namepsace `rk` are generated automatically. All of this function have the same signature : `RKType<_T>`. This is easy to make an array of this type (I don't know why but it seems to work only when I encapsulate this RK functions in a lambda function) and make some test on every RK methods.
**/

{% for rk in rk_list %}
{# /* for this template need `rk={ label=str , stages=[str] , un1=str }` each `str` is a correct C++ expression */ #}
/**---
  function: `rk::{{ label }}`
  brief: compute one iteration of this RK(s,n) method
  arguments:
    - `u_n` (`valarray<_T>`): value of $u^n$ at time $t^n$
    - `L` (`LType<_T>`): function $L$ in $\dot{u} = L(u,t)$
    - `tn` (_T): time $t^n$ of iteration
    - `dt` (_T): time step $\Delta t$
---**/
template < typename _T >
std::valarray<_T>
{{ rk.label }} ( std::valarray<_T> const& u_n , LType<_T> L , _T tn , _T dt )
{
  {% for stage in rk.stages %}
  std::valarray<_T> u_{{ loop.index }} = {{ stage }};
  {% endfor %}
  return {{ rk.un1 }};
}
{% endfor %}

} // namespace rk

#endif

