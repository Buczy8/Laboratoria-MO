set title "Porówanie metod numerycznych z analitycznym"
set xlabel "t"
set ylabel "wartość"
set grid
set term wxt
plot  "laasonen.txt" using 1:4 with lines title "Wynik analityczny", \
     "laasonen.txt" using 1:3 with points title "Wynik laasonen z algorytmem thomasa", \
     "laasonen.txt" using 1:2 with points title "Wynik laasonen z dekompozycja LU", \

