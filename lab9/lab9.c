#include <stdio.h>
#include <stdlib.h>
#include <math.h>
int n = 1000; // Zmniejszona liczba punktów dla testów
double a = 0.0, b = 1.0;
double alfa = 0.0, beta = 1.0, gamma = -2.0;
double phi = 0.0, psi = 1.0, teta = 2.0;
double U1 = 2.0, U2 = -2.0;


void thomasAlgorithm(double L[], double D[], double U[], double b[], int n)
{
    for (int i = 1; i < n; i++)
    {
        double m = L[i] / D[i - 1];
        D[i] = D[i] - m * U[i - 1];
        b[i] = b[i] - m * b[i - 1];
    }

    b[n - 1] = b[n - 1] / D[n - 1];
    for (int i = n - 2; i >= 0; i--)
    {
        b[i] = (b[i] - U[i] * b[i + 1]) / D[i];
    }
}

double s(double x)
{
    return pow(x, 3) * 0.5;
}

double analytical_solution(double x)
{
    const double sqrt5 = sqrt(5.0);
    const double exp2sqrt5 = exp(2.0 * sqrt5);
    const double coth_sqrt5 = 1.0 / tanh(sqrt5);

    double term1 = -95.0 * exp((-1.0 - sqrt5) * (-1.0 + x));
    double term2 = 55.0 * exp((-1.0 + sqrt5) * x);
    double term3 = 95.0 * exp(1.0 + sqrt5 + (-1.0 + sqrt5) * x);
    double term4 = -55.0 * exp(2.0 * sqrt5 - (1.0 + sqrt5) * x);
    double term5 = 2.0 * x * (6.0 + x * (3.0 + 2.0 * x));
    double term6 = -exp2sqrt5 * (9.0 + 2.0 * x * (6.0 + x * (3.0 + 2.0 * x)));

    double numerator = 9.0 + term1 + term2 + term3 + term4 + term5 + term6;
    double denominator = -1.0 + coth_sqrt5;

    return -numerator * denominator / 64.0;
}

void setup_tridiagonal_system(double L[], double D[], double U[], double b[], double x[], int n, double h,
                              double alfa, double beta, double gamma,
                              double phi, double psi, double teta)
{
    // Wewnętrzne węzły (i = 1..n-1)
    for (int i = 1; i < n; i++)
    {
        L[i] = (1.0 / (h * h)) - (1.0 / h);
        D[i] = -4.0 - (2.0 / (h * h));
        U[i] = (1.0 / (h * h)) + (1.0 / h);
        b[i] = -s(x[i]);
    }

    // Warunki brzegowe
    D[0] = beta - (alfa / h);
    U[0] = alfa / h;
    b[0] = -gamma;

    D[n] = phi / h + psi;
    L[n] = -phi / h;
    b[n] = -teta;
}

double check_shot(double p, double a, double U1, double Ub, int n, double h, double* y)
{
    double xi = a;

    y[0] = U1;
    y[1] = h * p + U1;

    for (int i = 1; i < n; i++)
    {
        xi += h;
        y[i + 1] = (y[i] * (2 + 4 * h * h) + y[i - 1] * (h - 1) - (xi * xi * xi * h * h) / 2) / (h + 1);
    }

    return y[n] - Ub;
}

double* RR_shoot(double shot_l, double shot_h, int n, double h, const double tol, int max_iter)
{
    int i = 0;
    double shot_m, check_m, check_l, check_h, tmp_m;

    // Alokacja pamięci dla tablicy y
    double* y = (double*)malloc((n + 1) * sizeof(double));
    if (y == NULL)
    {
        printf("Błąd alokacji pamięci dla tablicy y.\n");
        exit(1);
    }

    check_l = check_shot(shot_l, a, U1, U2, n, h, y);
    check_h = check_shot(shot_h, a, U1, U2, n, h, y);

    if ((check_l > 0 && check_h > 0) || (check_l < 0 && check_h < 0))
    {
        printf("Brak zmiany znaku: check_l = %f, check_h = %f\n", check_l, check_h);
        free(y);
        return NULL;
    }

    do
    {
        i++;

        tmp_m = shot_m;
        shot_m = (shot_h + shot_l) / 2;
        check_m = check_shot(shot_m, a, U1, U2, n, h, y);

        if (fabs(check_m) < tol)
        {
            break;
        }

        if ((check_l > 0 && check_m < 0) || (check_l < 0 && check_m > 0))
        {
            shot_h = shot_m;
            check_h = check_m;
        }
        else
        {
            shot_l = shot_m;
            check_l = check_m;
        }
    }
    while (fabs(check_m) > tol && i < max_iter);

    if (i >= max_iter)
    {
        printf("Osiągnięto maksymalną liczbę iteracji.\n");
    }

    return y;
}

double maxAbsoluteError(double* analytical, double* numerical, int size)
{
    double maxError = 0.0;
    for (int i = 0; i < size; i++)
    {
        double error = fabs(analytical[i] - numerical[i]);
        if (error > maxError)
        {
            maxError = error;
        }
    }
    return maxError;
}

int main()
{
    FILE* file = fopen("wyniki.txt", "w");
    if (file == NULL)
    {
        printf("Błąd: Nie można otworzyć pliku.\n");
        return 1;
    }

    int n = 1000;
    double h = (b - a) / n;

    double* x = malloc((n + 1) * sizeof(double));
    double* analytical = malloc((n + 1) * sizeof(double));
    double* lower = malloc((n + 1) * sizeof(double));
    double* diag = malloc((n + 1) * sizeof(double));
    double* upper = malloc((n + 1) * sizeof(double));
    double* rhs = malloc((n + 1) * sizeof(double));

    for (int i = 0; i <= n; i++)
    {
        x[i] = a + i * h;
        analytical[i] = analytical_solution(x[i]);
    }

    setup_tridiagonal_system(lower, diag, upper, rhs, x, n, h, alfa, beta, gamma, phi, psi, teta);
    thomasAlgorithm(lower, diag, upper, rhs, n + 1);

    double* u = RR_shoot(-10.0, 30.0, n, h, 1e-12, 100);
    if (u == NULL)
    {
        printf("Błąd: Nie można znaleźć rozwiązania.\n");
        goto cleanup;
    }

    for (int i = 0; i <= n; i++)
    {
        fprintf(file, "%.6f\t%.6f\t%.6f\t%.6f\n", x[i], analytical[i], rhs[i], u[i]);
    }
    fclose(file);

    FILE* errorFile = fopen("errors.txt", "w");
    if (errorFile == NULL)
    {
        printf("Błąd: Nie można otworzyć pliku z błędami.\n");
        goto cleanup;
    }

    for (n = 100; n <= 100000; n += 100)
    {
        h = (b - a) / n;

        double* x_new = malloc((n + 1) * sizeof(double));
        double* analytical_new = malloc((n + 1) * sizeof(double));
        double* lower_new = malloc((n + 1) * sizeof(double));
        double* diag_new = malloc((n + 1) * sizeof(double));
        double* upper_new = malloc((n + 1) * sizeof(double));
        double* rhs_new = malloc((n + 1) * sizeof(double));

        for (int i = 0; i <= n; i++)
        {
            x_new[i] = a + i * h;
            analytical_new[i] = analytical_solution(x_new[i]);
        }

        setup_tridiagonal_system(lower_new, diag_new, upper_new, rhs_new, x_new, n, h, alfa, beta, gamma, phi, psi,
                                 teta);
        thomasAlgorithm(lower_new, diag_new, upper_new, rhs_new, n + 1);
        double* u_new = RR_shoot(-20.0, 20.0, n, h, 1e-8, 1000);

        // for (int i = 0; i <= n; i++) {
        //     printf("%.6f\t%.6f\t%.6f\t%.6f\n", x_new[i], analytical_new[i], rhs_new[i], u_new[i]);
        // }
        // printf("\n\n\n\n\n");

        if (u_new == NULL)
        {
            printf("Błąd: Nie można znaleźć rozwiązania dla n = %d.\n", n);
            free(x_new);
            free(analytical_new);
            free(lower_new);
            free(diag_new);
            free(upper_new);
            free(rhs_new);
            continue;
        }

        double thomasError = maxAbsoluteError(analytical_new, rhs_new, n + 1);
        double shootingError = maxAbsoluteError(analytical_new, u_new, n + 1);
        fprintf(errorFile, "%f\t%.6e\t%.6e\n", log10(h), log10(thomasError), log10(shootingError));

        free(x_new);
        free(analytical_new);
        free(lower_new);
        free(diag_new);
        free(upper_new);
        free(rhs_new);
        free(u_new);
    }

    fclose(errorFile);

cleanup:
    free(x);
    free(analytical);
    free(lower);
    free(diag);
    free(upper);
    free(rhs);
    free(u);
    return 0;
}
