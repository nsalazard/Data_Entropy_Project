set term pdf
set out 'plot1.pdf'
set xlabel 'Blocksize'
set ylabel 'Tiempo CPU'
set logscale
set title 'Tiempo de CPU vs Blocksize'
plot "bloc_16.txt" u 1:2 w lp lw 2 t "N=16","bloc_32.txt" u 1:2 w lp lw 2 t "N=32","bloc_64.txt" u 1:2 w lp lw 2 t "N=64"