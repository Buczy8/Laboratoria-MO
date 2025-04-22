set title "Błędy przybliżonych wartości pierwszych pochodnych funkcji f(x) = cos(x) dla zmiennych rzeczywistych typu long double"
set xlabel "log10(h)"
set ylabel "log10(error)"
set grid
set term wxt
plot "long_double_errors.txt" using 1:2 with lines title "różnica dwupunktowa punktu początkowego", \
     "long_double_errors.txt" using 1:3 with lines title "różnica trzypunktowa punktu początkowego", \
     "long_double_errors.txt" using 1:4 with lines title "różnica progresywna punktu centralnego", \
     "long_double_errors.txt" using 1:5 with lines title "różnica wsteczna punktu centralnego", \
     "long_double_errors.txt" using 1:6 with lines title "różnica centralna punktu centralnego", \
     "long_double_errors.txt" using 1:7 with lines title "różnica dwupunktowa punktu końcowego", \
     "long_double_errors.txt" using 1:8 with lines title "różnica trzypunktowa punktu końcowego", \

