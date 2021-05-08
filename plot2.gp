set term pdf
set out 'plot2.pdf'
set xlabel 'Blocksize'
set ylabel 'Tiempo CPU'
set logscale
set title 'Tiempo de CPU vs Blocksize'
plot "bloc_128.txt" u 1:2 w lp lw 2 t "N=128","bloc_256.txt" u 1:2 w lp lw 2 t "N=256"