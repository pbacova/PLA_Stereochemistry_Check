set term postscript eps enhanced color font 'Helvetica,28'
set output "inner_dist.eps"
set rmargin 2.5
set tmargin 3.1
set bmargin 2.5
set lmargin 8
set mxtics 10
set mytics 10
set logscale x
set key right bottom
set xrange [1:]
set xlabel "distance n between the backbone atoms" offset 0,0.7
set ylabel "<R_n^2>/n [nm^2]" offset 2.1,-0.9

p "in_d_all" u 1:2 w l lw 6 title "PLA"
