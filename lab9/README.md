
# Laboratorium 9


Napisz program w języku „C/C++”, rozwiązujący równanie różniczkowe zwyczajne drugiego rzędu:
U’’(x) + 2 U’(x) – 4 U(x) + x
3
/2 = 0, określone na przedziale 0 <= x <= 1,
z warunkami brzegowymi U(0) = 2, U(1) = -2. Zastosuj typ double oraz trzypunktową dyskretyzację
konwencjonalną oraz metodę strzałów, na sieci jednorodnej. Do rozwiązania układu liniowych
równań algebraicznych zastosuj algorytm Thomasa (patrz zajęcia nr 6). Wykonaj rysunek
przedstawiający porównanie uzyskanych wyników numerycznych z rozwiązaniem analitycznym

U(x) = -((9 – 95 exp((-1 - sqrt(5))(-1 + x)) + 55 exp((-1 + sqrt(5))x) +
95 exp(1 + sqrt(5) + (-1 + sqrt(5))x) - 55 exp(2 sqrt(5) - (1 + sqrt(5))x) +
2x (6 + x(3 + 2x)) - exp(2 sqrt(5)) (9 + 2x (6 + x (3 + 2x))))(-1 + coth(sqrt(5))))/ 64.

Pokaż, że rząd dokładności rozwiązań numerycznych jest zgodny z przewidywaniami
teoretycznymi. W tym celu wykonaj (na jednym rysunku) wykresy przedstawiające zależności
maksymalnego błędu bezwzględnego rozwiązań od kroku sieci h, posługując się skalą logarytmiczną
(tzn. wykresy zależności log10|błędu| od log10 h ). Na podstawie wykresów wyznacz doświadczalnie
rzędy dokładności rozwiązań uzyskanych za pomocą obu metod, i porównaj je z rzędami
teoretycznymi. Ponadto zidentyfikuj wartości kroku sieci poniżej których pojawia się wpływ błędów
maszynowych.