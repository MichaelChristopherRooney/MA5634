#!/usr/bin/gnuplot -persist
set terminal postscript eps enhanced color font 'Helvetica,10'
set output 'graphs/delta_vs_acceptance.eps'
set xlabel 'Delta'
set ylabel 'Acceptance rate (%)'
plot 'data/delta_vs_acceptance_rate.txt' using 1:2 with linespoints lt rgb 'black'
exit

