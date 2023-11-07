#!/usr/bin/env gnuplot
if ( GPVAL_VERSION >= 4.4 && strstrt(GPVAL_TERMINALS, 'wxt') > 0 ) set terminal wxt persist
if ( GPVAL_VERSION >= 4.4 && strstrt(GPVAL_TERMINALS, 'qt') > 0 ) set terminal qt persist
if ( GPVAL_VERSION >= 4.4 && strstrt(GPVAL_TERMINALS, 'wxt') == 0 && strstrt(GPVAL_TERMINALS, 'qt') == 0 ) print "wxt and qt terminals not available, proceeding with default"
if ( GPVAL_VERSION < 4.4 ) print "gnuplot is too old to check for available terminals" ; print "attempting to use wxt terminal and hoping for the best" ; set terminal wxt persist

#set cbrange [0:1]
set logscale cb
set xlabel "log_{10}(E/GeV)"
set ylabel "Cos(zenith)"
set pm3d map
splot "fluxes_flavor.txt" u 1:2:5

set xrange [2.:6.]
set terminal png size 800,600 enhanced font 'Verdana,10'
set output "plot.png"
#set terminal postscript enhance eps color
#set output "plot.eps"
#set terminal svg size 350,262 fname 'Verdana' fsize 10
#set output "plot.svg"
replot


