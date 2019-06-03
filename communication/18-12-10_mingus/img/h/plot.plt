#!/usr/bin/env gnuplot

set terminal png size 800,600 enhanced

# plot v3g_3 avec hv3g_3
set output "neuh_1_g.png"
plot "diag_g.dat" u 1:4 with lines lt 4 title "Flux numérique de référence (WENO)" , "diag_g_h0307.dat" u 1:4 with lines lt 11 title "Flux numérique approximé (WENO)"

set output "neuh_2_g.png"
plot "diag_g.dat" u 1:4 with lines lt 4 title "Flux numérique de référence (WENO)" , "diag_g_h035065.dat" u 1:4 with lines lt 11 title "Flux numérique approximé (WENO)" , "diag_g_compact_h035065.dat" u 1:4 with lines lt 6 title "Flux numérique approximé (compact)"

set output "neuh_3_g.png"
plot "diag_g.dat" u 1:4 with lines lt 4 title "Flux numérique de référence (WENO)" , "diag_g_htrap.dat" u 1:4 with lines lt 11 title "Flux numérique approximé (WENO)" , "diag_g_compact_htrap.dat" u 1:4 with lines lt 6 title "Flux numérique approximé (compact)"

set output "neuh_4_g.png"
plot "diag_g.dat" u 1:4 with lines lt 4 title "Flux numérique de référence (WENO)" , "diag_g_ht.dat" u 1:4 every 20 with points lt 11 lw 0 title "Flux numérique approximé (WENO)"

