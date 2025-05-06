set title "Porówanie metod numerycznych z analitycznym"
set xlabel "x"
set ylabel "wartość"
set grid
set term wxt
plot "errors.txt" using 1:2 with lines title "błąd thomas", \
 "errors.txt" using 1:3 with lines title "błąd metod strzałów", \