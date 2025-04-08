set title "Blad wzgledny"
set xlabel "log10(x)"
set ylabel "log10(error)"
set grid
set term wxt
plot "wyniki.txt" using 1:2 with lines title "standardowe sinh"

