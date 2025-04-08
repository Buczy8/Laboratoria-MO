set title "Blad wzgledny"
set xlabel "log10(x)"
set ylabel "log10(error)"
set grid
set term wxt
plot "wyniki3.txt" using 1:2 with lines title "standardowe sinh", \
      "wyniki3.txt" using 1:3 with lines title "alternatywne sinh"

