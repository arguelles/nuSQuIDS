#!/usr/bin/env gnuplot
if ( GPVAL_VERSION >= 4.4 && strstrt(GPVAL_TERMINALS, 'wxt') > 0 ) set terminal wxt persist
if ( GPVAL_VERSION >= 4.4 && strstrt(GPVAL_TERMINALS, 'qt') > 0 ) set terminal qt persist
if ( GPVAL_VERSION >= 4.4 && strstrt(GPVAL_TERMINALS, 'wxt') == 0 && strstrt(GPVAL_TERMINALS, 'qt') == 0 ) print "wxt and qt terminals not available, proceeding with default"
if ( GPVAL_VERSION < 4.4 ) print "gnuplot is too old to check for available terminals" ; print "attempting to use wxt terminal and hoping for the best" ; set terminal wxt persist
set key box
set key opaque
set yrange [-0.1:1.1]
set xlabel "log_{10}(E/GeV)"
set ylabel "Muon Flux Ratio"
plot "DMpmns.txt" u 1:($4) w l title "Q1"

replot "DMpmns.txt" u 1:($6) w l title "Q2"

replot "DMpmns.txt" u 1:($8) w l title "Q3"

replot "DMpmns.txt" u 1:($10) w l title "S0"

replot "DMpmns.txt" u 1:($12) w l title "E"

replot "DMpmns.txt" u 1:($14) w l title "M"

replot "DMpmns.txt" u 1:($16) w l title "T"

replot "DMpmns.txt" u 1:($4+$6+$8+$10+$12+$14+$16) w l title "All"

set terminal postscript enhance eps color
set output "plot.eps"
replot
