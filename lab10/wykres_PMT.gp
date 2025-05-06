set title "Porówanie metod numerycznych z analitycznym"
set xlabel "t"
set ylabel "wartość"
set grid
set term wxt
plot  "PMT.txt" using 1:3 with points title "Metoda pośrednia trapezów (PMT)", \
    "PMT.txt" using 1:2 with lines title "Wynik analityczny", \