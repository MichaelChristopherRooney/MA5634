#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color font 'Helvetica,10'
#####################
# first plot cos(x)
#####################
set output 'graphs/cosx-history.eps'
set xrange [0:1000]
set yrange [-10:10]
set xlabel 'Time'
plot 'data/cosx-1.txt' using 1:2 with lines lt rgb 'black' title 'Delta = 1.5',\
'data/cosx-2.txt' using 1:2 with lines lt rgb 'red' title 'Delta = 15',\
'data/cosx-3.txt' using 1:2 with lines lt rgb 'green' title 'Delta = 150'
#####################
# now plot x*x
#####################
set output 'graphs/xsqr-history.eps'
set xrange [0:1000]
set yrange [-10:10]
set xlabel 'Time'
plot 'data/xsqr-1.txt' using 1:2 with lines lt rgb 'black' title 'Delta = 1.5',\
'data/xsqr-2.txt' using 1:2 with lines lt rgb 'red' title 'Delta = 15',\
'data/xsqr-3.txt' using 1:2 with lines lt rgb 'green' title 'Delta = 150'
exit

