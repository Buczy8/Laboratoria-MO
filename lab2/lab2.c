#include <stdio.h>
#include <math.h>
#include <float.h>


long double f(double x) {
    return (x * x * x) / (6.L * (sinh(x) - x));
}

long double blad(long double f_approx, long double f_exact) {
    return fabsl((f_approx - f_exact) / f_exact);
}

long double silnia(int n) {
    long double result = 1;
    for (int i = 1; i <= n; i++) {
        result *= i;
    }
    return result;
}
//rozwiniecie funkcji sinh(x)-x w szereg
long double sin_h(double x) {
    long double result = 0.L;
    for (int i = 1; i < 800; i++) {
        result += powl(x, 2 * i + 1) / silnia(2 * i + 1);
    }
    return result;
}

long double g(double x) {
    return (x * x * x) / (6.L * (sin_h(x)));
}
int main() {
    FILE *file = fopen("dane_do_laboratorium_2.txt", "r");  // Otwórz plik do odczytu
    if (file == NULL) {
        printf("Błąd: Nie można otworzyć pliku.\n");
        return 1;
    }
    FILE *file2 = fopen("wyniki.txt", "w");  // Otwórz plik do zapisu
    if (file2 == NULL) {
        printf("Błąd: Nie można otworzyć pliku.\n");
        return 1;
    }
    FILE *file3 = fopen("wyniki2.txt", "w");  // Otwórz plik do zapisu
    if (file3 == NULL) {
        printf("Błąd: Nie można otworzyć pliku.\n");
        return 1;
    }
    FILE *file4 = fopen("wyniki3.txt", "w");  // Otwórz plik do zapisu
    if (file4 == NULL) {
        printf("Błąd: Nie można otworzyć pliku.\n");
        return 1;
    }
    float log10_x;
    double x;
    long double f_x;
    //printf("%Le\n", silnia(200));
    // Wczytanie danych z pliku i wyświetlenie ich
    int i = 0;
    while (fscanf(file, "%f %le %Le", &log10_x, &x, &f_x) == 3) {
        //printf("%10.5f %25.20e %1.5Le %1.5Le %1.5Le\n", log10_x ,x, f(x), f_x, g(x));
        fprintf(file2, "%.5e %.5Le \n", log10_x, log10l(blad(f(x), f_x))); // ostatnie to bład reprezentacni
        fprintf(file3, "%.5e %.5Le \n", log10_x, log10l(blad(g(x), f_x)));
        fprintf(file4, "%.5le %.20Le %.20Le \n", log10(x), log10l(blad(f(x), f_x)), log10l(blad(g(x), f_x)));
    }

    fclose(file);  // Zamknięcie pliku
    fclose(file2);  // Zamknięcie pliku
    fclose(file3);  // Zamknięcie pliku
    fclose(file4);  // Zamknięcie pliku

    return 0;
}

