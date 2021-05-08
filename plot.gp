set term pdf
set out 'plot.pdf'
set xlabel 'Blocksize'
set ylabel 'Tiempo CPU'
set logscale
set title 'Tiempo de CPU vs Blocksize'
plot "bloc_2.txt" u 1:2 w lp lw 2 t "N=2", "bloc_4.txt" u 1:2 w lp lw 2 t "N=4","bloc_8.txt" u 1:2 w lp lw 2 t "N=8", "bloc_16.txt" u 1:2 w lp lw 2 t "N=16","bloc_32.txt" u 1:2 w lp lw 2 t "N=32","bloc_64.txt" u 1:2 w lp lw 2 t "N=64", "bloc_128.txt" u 1:2 w lp lw 2 t "N=128","bloc_256.txt" u 1:2 w lp lw 2 t "N=256", "bloc_512.txt" u 1:2 w lp lw 2 t "N=512","bloc_1024.txt" u 1:2 w lp lw 2 t "N=1024"

 