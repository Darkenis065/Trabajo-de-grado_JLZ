set term pngcairo size 1000,700 font 'Arial,10'
set output 'osciladorarmonico.png'
set title 'Harmonic oscilator' font 'Arial,12'
set xlabel 'Time (t)'
set ylabel 'X(t)'
set grid lc rgb '#d3d3d3' lt 1 lw 0.5
set style line 1 lc rgb '#0066cc' lt 1 lw 2
set style line 2 lc rgb '#cc0000' lt 1 lw 2
set style line 3 lc rgb '#404040' lt 0 lw 1.5
set style line 4 lc rgb '#808080' lt 0 lw 1.5
plot 'datos_analiticos.dat' using 1:2 with lines linestyle 3 title 'x(t) Analytical', \
     'datos_analiticos.dat' using 1:3 with lines linestyle 4 title 'v(t) Analytical', \
     'datos_rk4.dat' using 1:2 with lines linestyle 1 title 'x(t) Runge-Kutta 4', \
     'datos_rk4.dat' using 1:3 with lines linestyle 2 title 'v(t) Runge-Kutta 4'
