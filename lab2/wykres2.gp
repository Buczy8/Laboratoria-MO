set title "Blad wzgledny"
set xlabel "log10(x)"
set ylabel "log10(error)"
set grid
set term wxt
plot "wyniki2.txt" with lines title "alternatywne sinh"
