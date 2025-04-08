#include <stdio.h>

#include <math.h>

long double a(long double x) {
    return tanhl(x) + 2 * (x - 1);
}

long double b(long double x) {
    return sinhl(x) + (x / 4) - 1;
}

long double a_phi(long double x) {
    return -0.5 * tanhl(x) + 1;
}

long double b_phi(long double x) {
    return asinhl(1 - x / 4);
}

long double a_prime(long double x) {
    return 2 + 1 / (coshl(x) * coshl(x));
}

long double b_prime(long double x) {
    return coshl(x) + 0.25;
}

long double
picard(long double tolx, long double tolf, int n_max, long double (*phi)(long double), long double (*f)(long double)) {
    long double x = 0;
    int n = 0;

    while (n < n_max) {
        long double x_n = phi(x);
        printf("Iteracja %d: x = %.21Lf  f(x) = %.5Le \n", n, x_n, fabsl(f(x_n)));

        if (fabsl(x_n - x) < tolx && fabsl(f(x_n)) < tolf) {
            return x_n;
        }

        x = x_n;
        n++;
    }

    printf("Nie osiagnieto wymaganej tolerancji w %d iteracjach.\n", n_max);
    return x;
}

long double
bisekcja(long double a, long double b, long double tolx, long double tolf, int n_max, long double (*f)(long double)) {
    long double fa = f(a);
    long double fb = f(b);

    if((fa>0 && fb>0) || (fa<0 && fb<0)){
        printf("Funkcja nie spelnia zalozen metody bisekcji.\n");
        return 0;
    }

    int n = 0;
    while (n < n_max) {
        long double x_n = (a + b) / 2;
        long double fx_n = f(x_n);
        printf("Iteracja %d: x = %.21Lf f(x) = %.5Le\n", n, x_n, fabsl(fx_n));

        if (fabsl(fx_n) < tolf && fabsl((b - a) / 2) < tolx) {
            return x_n;
        }

        if (fa * fx_n < 0) {
            b = x_n;
            fb = fx_n;
        } else {
            a = x_n;
            fa = fx_n;
        }

        n++;
    }

    printf("Nie osiagnieto wymaganej tolerancji w %d iteracjach.\n", n_max);
    return 0;
}

long double newton(long double tolx, long double tolf, int n_max, long double (*f)(long double),
                   long double (*f_prime)(long double)) {
    int n = 0;
    long double x = 0;

    while (n < n_max) {
        long double fx = f(x);
        long double fx_prime = f_prime(x);
        long double x_n = x - fx / fx_prime;
        printf("Iteracja %d: x = %.21Lf f(x) = %.5Le\n", n, x_n, fabsl(f(x_n)));

        if (fabsl(f(x_n)) < tolf && fabsl(x_n - x) < tolx) {
            return x_n;
        }

        x = x_n;
        n++;
    }
    printf("Nie osiagnieto wymaganej tolerancji w %d iteracjach.\n", n_max);
    return 0;
}

long double sieczne(long double tolx, long double tolf, int n_max, long double (*f)(long double)) {
    int n = 0;
    long double x0 = 0;
    long double x1 = 1;

    while (n < n_max) {
        long double fx0 = f(x0);
        long double fx1 = f(x1);
        long double x2 = x1 - fx1 * ((x1 - x0) / (fx1 - fx0));
        printf("Iteracja %d: x = %.21Lf f(x) = %.5Le\n", n, x2, fabsl(f(x2)));

        if (fabsl(f(x2)) < tolf && fabsl(x2 - x1) < tolx) {
            return x2;
        }

        x0 = x1;
        x1 = x2;
        n++;
    }

    printf("Nie osiagnieto wymaganej tolerancji w %d iteracjach.\n", n_max);
    return 0;
}

int main() {

    long double tolx = 1e-12;
    long double tolf = 1e-12;
    int n_max = 1000;

    printf("=======================================\n");
    printf("Metoda Picarda dla rownania a\n");
    picard(tolx, tolf, n_max, a_phi, a);
    printf("=======================================\n");
    printf("Metoda Picarda dla rownania b\n");
    picard(tolx, tolf, n_max, b_phi, b);
    printf("=======================================\n");

    printf("Metoda bisekcji dla rownania a\n");
    bisekcja(0, 1, tolx, tolf, n_max, a);
    printf("=======================================\n");
    printf("Metoda bisekcji dla rownania b\n");
    bisekcja(0, 1, tolx, tolf, n_max, b);
    printf("=======================================\n");

    printf("Metoda newtona dla rownania a\n");
    newton(tolx, tolf, n_max, a, a_prime);
    printf("=======================================\n");
    printf("Metoda newtona dla rownania b\n");
    newton(tolx, tolf, n_max, b, b_prime);
    printf("=======================================\n");

    printf("Metoda siecznych dla rownania a\n");
    sieczne(tolx, tolf, n_max, a);
    printf("=======================================\n");
    printf("Metoda siecznych dla rownania b\n");
    sieczne(tolx, tolf, n_max, b);
    printf("=======================================\n");
    return 0;
}
