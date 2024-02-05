#!/usr/bin/gnuplot

set term pngcairo size 1560,860 enhanced font 'Verdana, 15'
set output 'MapaCalor.png'
set border linewidth 1.5

set pm3d map interpolate 0,0
set xlabel 
set ylabel 
set cblabel 'W(x,y)'
set palette rgb 34,35,36

splot 'matriz_solucao.dat' matrix