set title "Porówanie błedów metod"
set xlabel "log10l(dt)"
set ylabel "log10(err)"
set grid
set term wxt
plot  "errors.txt" using 1:2 with lines title "BME", \
    "errors.txt" using 1:3 with lines title "PME", \
    "errors.txt" using 1:4 with lines title "PMT", \