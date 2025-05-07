set title "Porówanie metod numerycznych z analitycznym"
set xlabel "t"
set ylabel "wartość"
set grid
set term wxt
plot  "BME.txt" using 1:2 with lines title "Wynik analityczny", \
    "BME.txt" using 1:3 with points title "Metoda BME stablina", \
    "BME.txt" using 1:4 with points title "Metoda BME rozniestablina", \
    "BME.txt" using 1:5 with points title "Metoda BME niestablina", \