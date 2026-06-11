set term png
set output 'Error.png'
set xlabel 'pasos'
set ylabel 'error'
set xrange [0:10] 
set grid 
plot 'Datos_runge_kutta4.dat' u 1:6 with l title 'Error velocidad con Rungekutta', 'Datos_runge_kutta4.dat' u 1:7 with l title 'Error posición con Rungekutta'