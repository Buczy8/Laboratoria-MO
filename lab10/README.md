
# Laboratorium 10


Napisz program w języku „C/C++”, rozwiązujący równanie różniczkowe y'(t) + ((100 t + 10)/(t + 1)) (y(t) - 1) = 0, określone dla zmiennej t  0,
z warunkiem początkowym y(0) = 2, za pomocą
metod:

(a) bezpośredniej Eulera

(b) pośredniej Eulera

(c) metody trapezów.

Dla metod (b) i (c) wykonaj oddzielne rysunki przedstawiające po dwa wykresy: wykres
przykładowego rozwiązania numerycznego oraz (dla porównania) wykres rozwiązania
analitycznego: y(t) = 1 + (1 + t)
90 exp(−100t). Oba wykresy winny przedstawiać zależność y od
zmiennej niezależnej t. Rozwiązania analityczne zaznacz linią ciągłą, a numeryczne punktami. W
przypadku metody (a) wykonaj trzy takie rysunki: jeden uzyskany w warunkach numerycznej
stabilności metody (dla dowolnego t), drugi w warunkach numerycznej niestabilności (też dla
dowolnego t), a trzeci w pozostałym przypadku. Wyjaśnij różnice pomiędzy uzyskanymi wykresami.

Pokaż, że rząd dokładności uzyskanych stabilnych rozwiązań numerycznych jest zgodny z
przewidywaniami teoretycznymi. W tym celu wykonaj (na jednym rysunku) wykresy
przedstawiające zależności maksymalnych błędów bezwzględnych rozwiązań uzyskanych trzema
metodami, od kroku sieci czasowej t, posługując się skalą logarytmiczną (tzn. wykresy zależności
log10|błędu| od log10 t ). Na podstawie wykresów wyznacz doświadczalnie rzędy dokładności
rozwiązań uzyskanych za pomocą różnych metod i porównaj je z rzędami teoretycznymi. O ile to
możliwe, zidentyfikuj też wartości kroku sieci poniżej których pojawia się wpływ błędów
maszynowych.