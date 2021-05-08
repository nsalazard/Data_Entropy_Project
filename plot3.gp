set term pdf
set out 'plot3.pdf'
set xlabel 'Blocksize'
set ylabel 'Tiempo CPU'
set logscale
set title 'Tiempo de CPU vs Blocksize'
plot "bloc_512.txt" u 1:2 w lp lw 2 t "N=512","bloc_1024.txt" u 1:2 w lp lw 2 t "N=1024"