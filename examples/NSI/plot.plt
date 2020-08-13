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
#plot "fluxes_flavor.txt" u 1:($3) w l title "Electron NSI"
#replot "fluxes_flavor.txt" u 1:($4) w l title "Electron noNSI"

set style line 1 lt 6 lc rgb "red" lw 1
set style line 2 lt 1 lc rgb "red" lw 3

set style line 3 lt 6 lc rgb "blue" lw 1
set style line 4 lt 1 lc rgb "blue" lw 3

plot "fluxes_flavor.txt" u 1:($5) w l ls 1 title "Muon NSI"
replot "fluxes_flavor.txt" u 1:($6) w l ls 2 title "Muon noNSI"
replot "fluxes_flavor.txt" u 1:($7) w l ls 3 title "Tau NSI"
replot "fluxes_flavor.txt" u 1:($8) w l ls 4 title "Tau noNSI"



set terminal postscript enhance eps color
set output "plot.eps"
replot


