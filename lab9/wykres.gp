set title "Porówanie metod numerycznych z analitycznym"
set xlabel "x"
set ylabel "wartość"
set grid
set term wxt
plot "wyniki.txt" using 1:2 with lines title "Wynik analityczny", \
 "wyniki.txt" using 1:3 with lines title "Wynik trzypunktowej dyskretyzacji konwencjonalnej", \
 "wyniki.txt" using 1:4 with lines title "Wynik metody strzałów", \



