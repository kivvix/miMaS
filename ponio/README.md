# `ponio` description

`ponio` is a *Python* library for study of numerical integrators for solve linear transport equation :

$$
  u_t + u_x = 0
$$

In time we consider only Runge-Kutta methods, of different order and number of stages. In space we consider first WENO methods, in particular WENO5, but tools develop here can be use in many linear space schemes. We can automatically compute CFL number of some couple of RK(s,n) - WENO5, and various other type of couple.

