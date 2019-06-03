#!/usr/bin/env gnuplot

set terminal epslatex color font "Helvetica,12" standalone
xs = 0.3
xe = 0.7

set yrange [-0.1:1.1]
set samples 1000
set label "$x_s$" at (xs-0.02),-0.05
set label "$x_e$" at (xe-0.02),-0.05
h(x) = ( xs < x && x < xe ) ? 1 : 0

set output "h_gate.tex"
plot [0:1] h(x)

pas = 0.1
xss = xs+pas
xee = xe-pas
set yrange [-0.1:1.1]
set samples 1000
set label "$x_s$" at (xs-0.02),-0.05
set label "$x_e$" at (xe-0.02),-0.05
set label "$x_s^*$" at (xss-0.02),-0.05
set label "$x_e^*$" at (xee-0.02),-0.05
h(x) = ( xss < x && x < xee )? 1 : (xs < x && x<xss )? 1/(pas)*x-(xs)/pas : (xee <x && x<xe) ? -1/pas*x+(xe)/pas : 0
set arrow from xs,0 to xss,0 heads
set label "$\\delta x$" at xs+pas/2,0.05

set output "h_trap.tex"
plot [0:1] h(x)

set output "c_epsidt.tex"
unset arrow
unset label
set samples 1000
dx = 0.1
v = 18
h(x) = dx/v*x/(x-0.5*dx/v)
set arrow from 0.5*dx/v,-1 to 0.5*dx/v,1 nohead dt 2 lt 8
set arrow from -0.1,dx/v to 1,dx/v nohead dt 2 lt 8
set label "$\\frac{\\Delta x}{v_{max}}$" at 0.17,dx/v-0.00025
set label "$\\frac{1}{2}\\frac{\\Delta x}{v_{max}}$" at 0.5*dx/v+0.005,0.0015
set xlabel "$\\varepsilon$"
set ylabel "$\\Delta t$"
set yrange [0.0001:0.01]
plot [-0.01:0.2] h(x) notitle lt 6

