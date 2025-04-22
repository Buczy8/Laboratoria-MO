set title "Błędy przybliżonych wartości pierwszych pochodnych funkcji f(x) = cos(x) dla zmiennych rzeczywistych typu double"
set xlabel "log10(h)"
set ylabel "log10(error)"
set grid
set term wxt
plot "double_errors.txt" using 1:2 with lines title "różnica dwupunktowa punktu początkowego", \
     "double_errors.txt" using 1:3 with lines title "różnica trzypunktowa punktu początkowego", \
     "double_errors.txt" using 1:4 with lines title "różnica progresywna punktu centralnego", \
     "double_errors.txt" using 1:5 with lines title "różnica wsteczna punktu centralnego", \
     "double_errors.txt" using 1:6 with lines title "różnica centralna punktu centralnego", \
     "double_errors.txt" using 1:7 with lines title "różnica dwupunktowa punktu końcowego", \
     "double_errors.txt" using 1:8 with lines title "różnica trzypunktowa punktu końcowego", \

