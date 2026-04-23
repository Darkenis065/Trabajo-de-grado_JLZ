# ==============================================================================
# SCRIPT DE ANÁLISIS NUMÉRICO - ECUACIÓN DE DIRAC 1D
# ==============================================================================
# Configuración del motor de renderizado y salida gráfica
set terminal pngcairo size 1200,800 enhanced font 'Arial,11'
set output 'analisis_dirac.png'

# Configuración de márgenes globales para la matriz 2x2
set multiplot layout 2,2 title "Auditoría Dinámica y de Conservación - Simulación de Dirac 1D" font ",14" margins 0.1,0.95,0.1,0.9 spacing 0.12,0.15

# Archivo de datos fuente
DATA = 'dirac_observables.dat'

# ------------------------------------------------------------------------------
# PANEL 1: Verificación de Invariantes Físicos (Norma y Energía)
# Propósito: Auditar la estabilidad incondicional de Crank-Nicolson y hermiticidad.
# ------------------------------------------------------------------------------
set title "Conservación de Invariantes (Teorema de Noether Numérico)" font ",12"
set xlabel "Tiempo (t)"
set ylabel "Norma Total" textcolor rgb "blue"
set y2label "Energía Esperada <E>" textcolor rgb "red"
set y2tics
set ytics nomirror
set grid
set format y "%.7f"   # Alta precisión para detectar inestabilidades
set format y2 "%.5f"
# Forzamos rangos estrechos alrededor del valor esperado para visibilizar el ruido numérico
set yrange [0.999999:1.000001]
set y2range [5.08:5.11] # E = sqrt(5^2 + 1) = 5.099

plot DATA using 1:2 axes x1y1 with lines lw 2 lc rgb "blue" title "Norma ||Psi||^2", \
     DATA using 1:6 axes x1y2 with lines lw 2 lc rgb "red" title "<E>"

# ------------------------------------------------------------------------------
# PANEL 2: Cinemática del Paquete de Ondas (<x> y <p>)
# Propósito: Verificar el teorema de Ehrenfest para una partícula relativista libre.
# ------------------------------------------------------------------------------
unset y2tics
set ytics mirror
set format y "%g"
set autoscale y
set title "Teorema de Ehrenfest: Dinámica Esperada" font ",12"
set xlabel "Tiempo (t)"
set ylabel "Valor Esperado"
set key top left

plot DATA using 1:3 with lines lw 2 lc rgb "dark-green" title "Posición <x>", \
     DATA using 1:5 with lines lw 2 lc rgb "purple" title "Momento <p>"

# ------------------------------------------------------------------------------
# PANEL 3: Dispersión Temporal en la Malla Principal (Componente Mayor)
# Propósito: Evaluar el tránsito de la componente 1 sobre el baricentro geométrico.
# ------------------------------------------------------------------------------
set title "Módulo de la Componente Mayor |psi_1(0, t)|^2" font ",12"
set xlabel "Tiempo (t)"
set ylabel "Amplitud (Densidad Local)"
set key top right
set autoscale y

plot DATA using 1:9 with lines lw 2 lc rgb "dark-orange" title "|psi_1(x=0)|^2"

# ------------------------------------------------------------------------------
# PANEL 4: Dispersión Temporal en la Sub-malla Desplazada (Componente Menor)
# Propósito: Evaluar el tránsito de la componente 2 interpolada en el baricentro.
# ------------------------------------------------------------------------------
set title "Módulo de la Componente Menor |psi_2(0, t)|^2" font ",12"
set xlabel "Tiempo (t)"
set ylabel "Amplitud (Densidad Local)"
set autoscale y

plot DATA using 1:12 with lines lw 2 lc rgb "dark-cyan" title "|psi_2(x=0)|^2"

# Cierre ordenado del motor gráfico
unset multiplot
