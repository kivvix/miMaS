#!/usr/bin/env gnuplot

set terminal png size 800,600 enhanced

# plot just Euler+cinetic (upwind)
i=2
set output "neu_1_rho.png"
plot "diag_euler.dat" u 1:i with lines lt 4 title "Euler" , "diag_cine_upwind.dat" u 1:i with lines lt 6 title "Cinetique (upwind)"

i=3
set output "neu_1_u.png"
plot "diag_euler.dat" u 1:i with lines lt 4 title "Euler" , "diag_cine_upwind.dat" u 1:i with lines lt 6 title "Cinetique (upwind)"

i=4
set output "neu_1_T.png"
plot "diag_euler.dat" u 1:i with lines lt 4 title "Euler" , "diag_cine_upwind.dat" u 1:i with lines lt 6 title "Cinetique (upwind)"

# compare scheme on cinetic (upwind|compact|WENO)
i=2
set output "neu_2_rho.png"
plot "diag_cine_upwind.dat" u 1:i with lines lt 4 title "Cinetique (upwind)", "diag_cine_compact.dat" u 1:i with lines lt 6 title "Cinetique (compact)", "diag_cine_weno.dat" u 1:i with lines lt 2 title "Cinetique (WENO)"

i=3
set output "neu_2_u.png"
plot "diag_cine_upwind.dat" u 1:i with lines lt 4 title "Cinetique (upwind)", "diag_cine_compact.dat" u 1:i with lines lt 6 title "Cinetique (compact)", "diag_cine_weno.dat" u 1:i with lines lt 2 title "Cinetique (WENO)"

i=4
set output "neu_2_T.png"
plot "diag_cine_upwind.dat" u 1:i with lines lt 4 title "Cinetique (upwind)", "diag_cine_compact.dat" u 1:i with lines lt 6 title "Cinetique (compact)", "diag_cine_weno.dat" u 1:i with lines lt 2 title "Cinetique (WENO)"

# compare cinetic miMaS
i=2
set output "neu_3_rho.png"
plot "diag_cine_weno.dat" u 1:i with lines lt 4 title "Cinetique (WENO)" , "diag_mimas.dat" u 1:i with lines lt 6 title "micro-macro (WENO)"

i=3
set output "neu_3_u.png"
plot "diag_cine_weno.dat" u 1:i with lines lt 4 title "Cinetique (WENO)" , "diag_mimas.dat" u 1:i with lines lt 6 title "micro-macro (WENO)"

i=4
set output "neu_3_T.png"
plot "diag_cine_weno.dat" u 1:i with lines lt 4 title "Cinetique (WENO)" , "diag_mimas.dat" u 1:i with lines lt 6 title "micro-macro (WENO)"



