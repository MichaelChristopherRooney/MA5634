#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'graphs/bin_size_vs_variance.eps'

set xrange [5:100000]
set xlabel 'Bin size'
set ylabel 'Bin variance'
plot 'data/bin_size_vs_bin_variance.txt' using 1:2 with linespoints lt rgb 'black'

exit

