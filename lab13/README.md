
# Laboratorium 13

Napisz program w języku „C/C++”, obliczający numerycznie wartości funkcji:

$\mathrm{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-y^2} \, dy ,\, dla\, x = 1.0, 2.0, 3.0.$ 

Zastosuj złożone kwadratury:

(a) prostokątów (wariant z węzłem interpolacji po lewej stronie przedziału)

(b) prostokątów (wariant z węzłem interpolacji po prawej stronie przedziału)

(c) prostokątów (wariant z węzłem interpolacji w środku przedziału)

(d) trapezów

(e) parabol

na sieci o stałym kroku h. Oblicz błąd względny wyniku dla x = 3.0 w funkcji kroku h i
pokaż, że rzędy dokładności zastosowanych kwadratur są zgodne z przewidywaniami teoretycznymi.
W tym celu wykonaj (na jednym rysunku) wykresy zależności log10|błędu| od log10 h. Na podstawie
wykresów wyznacz doświadczalnie rzędy dokładności kwadratur. Do obliczenia ścisłych wartości
funkcji erf(x) (z dokładnością zbliżoną do maszynowej) zastosuj pakiet CALERF udostępniony
przez prowadzącego zajęcia.