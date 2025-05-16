set title "Porówanie metody numerycznych z analitycznym"
set xlabel "t"
set ylabel "wartość"
set grid
set term wxt
plot  "KMB.txt" using 1:3 with lines title "Wynik analityczny", \
     "KMB.txt" using 1:2 with points title "Wynik KMB", \

