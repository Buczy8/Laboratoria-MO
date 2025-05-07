#include <stdio.h>
#include <stdlib.h>
#include <math.h>

long double t0 = 0.0L, y_0 = 2.0L;

long double f(long double t, long double y);
long double a(long double t);
long double b(long double t);
long double analytic_solution(long double t);
long double maxAbsoluteError(long double* analytical, long double* numerical, int size);

int main()
{
    FILE* file = fopen("BME.txt", "w");
    FILE* file2 = fopen("PME.txt", "w");
    FILE* file3 = fopen("PMT.txt", "w");
    FILE* file4 = fopen("errors.txt", "w");

    if (file == NULL || file2 == NULL || file3 == NULL || file4 == NULL)
    {
        printf("Błąd: Nie można otworzyć pliku.\n");
        return 1;
    }

    long double dt = 1e-4L; // krok czasowy
    int n = 1000;

    long double* t = malloc((n + 1) * sizeof(long double));
    long double* y_bme_stabilne = malloc((n + 1) * sizeof(long double));
    long double* y_bme_rozniestabilne = malloc((n + 1) * sizeof(long double));
    long double* y_bme_niestabilne = malloc((n + 1) * sizeof(long double));
    long double* y_pme = malloc((n + 1) * sizeof(long double));
    long double* y_pmt = malloc((n + 1) * sizeof(long double));
    long double* y_exact = malloc((n + 1) * sizeof(long double));

    if (!t || !y_bme_stabilne || !y_bme_niestabilne || !y_pme || !y_pmt || !y_exact || !y_bme_rozniestabilne)
    {
        printf("Błąd alokacji pamięci.\n");
        return 1;
    }

    // Obliczenia dla ustalonego dt = 1e-4
    t[0] = t0;
    y_bme_stabilne[0] = y_0;
    y_bme_rozniestabilne[0] = y_0;
    y_bme_niestabilne[0] = y_0;
    y_pme[0] = y_0;
    y_pmt[0] = y_0;
    y_exact[0] = analytic_solution(t0);

    for (int i = 1; i <= n; i++)
    {
        t[i] = t[i - 1] + dt;
        y_bme_stabilne[i] = y_bme_stabilne[i - 1] + dt * f(t[i - 1], y_bme_stabilne[i - 1]);
        y_bme_rozniestabilne[i] = y_bme_rozniestabilne[i - 1] + 0.1L * f(t[i - 1], y_bme_rozniestabilne[i - 1]);
        y_bme_niestabilne[i] = y_bme_niestabilne[i - 1] + 0.5L * f(t[i - 1], y_bme_niestabilne[i - 1]);
        y_pme[i] = (y_pme[i - 1] + dt * b(t[i])) / (1.0L - dt * a(t[i]));
        y_pmt[i] = (y_pmt[i - 1] + dt / 2.0L * (a(t[i - 1]) * y_pmt[i - 1] + b(t[i - 1]) + b(t[i]))) / (1.0L - dt / 2.0L
            * a(t[i]));
        y_exact[i] = analytic_solution(t[i]);

        fprintf(file, "%.6Le\t%.6Le\t%.6Le\t%.6Le\t%.6Le\n", t[i], y_exact[i], y_bme_stabilne[i],y_bme_rozniestabilne[i], y_bme_niestabilne[i]);
        fprintf(file2, "%.6Le\t%.6Le\t%.6Le\n", t[i], y_exact[i], y_pme[i]);
        fprintf(file3, "%.6Le\t%.6Le\t%.6Le\n", t[i], y_exact[i], y_pmt[i]);
    }


    int num_steps = 10000;
    long double dt_min = 10e-14L;
    long double dt_max = 10e-3L;
    for (int j = 0; j < num_steps; j++)
    {
        long double exponent = j * log10l(dt_max / dt_min) / (num_steps - 1);
        long double dt2 = dt_min * powl(10.0L, exponent);

        long double* t_local = malloc((n + 1) * sizeof(long double));
        long double* y_bme_stabilne_l = malloc((n + 1) * sizeof(long double));
        long double* y_pme_l = malloc((n + 1) * sizeof(long double));
        long double* y_pmt_l = malloc((n + 1) * sizeof(long double));
        long double* y_exact_l = malloc((n + 1) * sizeof(long double));

        if (!t_local || !y_bme_stabilne_l || !y_pme_l || !y_pmt_l || !y_exact_l)
        {
            printf("Błąd alokacji pamięci dla dt = %.6Le.\n", dt2);
            return 1;
        }

        t_local[0] = t0;
        y_bme_stabilne_l[0] = y_0;
        y_pme_l[0] = y_0;
        y_pmt_l[0] = y_0;
        y_exact_l[0] = analytic_solution(t0);

        for (int i = 1; i <= n; i++)
        {
            t_local[i] = t_local[i - 1] + dt2;
            y_bme_stabilne_l[i] = y_bme_stabilne_l[i - 1] + dt2 * f(t_local[i - 1], y_bme_stabilne_l[i - 1]);
            y_pme_l[i] = (y_pme_l[i - 1] + dt2 * b(t_local[i])) / (1.0L - dt2 * a(t_local[i]));
            y_pmt_l[i] = (y_pmt_l[i - 1] + dt2 / 2.0L * (a(t_local[i - 1]) * y_pmt_l[i - 1] + b(t_local[i - 1]) +
                b(t_local[i]))) / (1.0L - dt2 / 2.0L * a(t_local[i]));
            y_exact_l[i] = analytic_solution(t_local[i]);
        }

        long double maxErrorBMEStabilne = maxAbsoluteError(y_exact_l, y_bme_stabilne_l, n + 1);
        long double maxErrorPME = maxAbsoluteError(y_exact_l, y_pme_l, n + 1);
        long double maxErrorPMT = maxAbsoluteError(y_exact_l, y_pmt_l, n + 1);

        fprintf(file4, "%.6Le\t%.6Le\t%.6Le\t%.6Le\n",
                log10l(dt2), log10l(maxErrorBMEStabilne), log10l(maxErrorPME), log10l(maxErrorPMT));

        free(t_local);
        free(y_bme_stabilne_l);
        free(y_pme_l);
        free(y_pmt_l);
        free(y_exact_l);
    }

    fclose(file);
    fclose(file2);
    fclose(file3);
    fclose(file4);

    free(t);
    free(y_bme_stabilne);
    free(y_bme_niestabilne);
    free(y_pme);
    free(y_pmt);
    free(y_exact);

    return 0;
}

long double f(long double t, long double y)
{
    return -((100.0L * t + 10.0L) / (t + 1.0L)) * (y - 1.0L);
}

long double a(long double t)
{
    return -((100.0L * t + 10.0L) / (t + 1.0L));
}

long double b(long double t)
{
    return ((100.0L * t + 10.0L) / (t + 1.0L));
}

long double analytic_solution(long double t)
{
    return 1.0L + powl(1.0L + t, 90.0L) * expl(-100.0L * t);
}

long double maxAbsoluteError(long double* analytical, long double* numerical, int size)
{
    long double maxError = 0.0L;
    for (int i = 0; i < size; i++)
    {
        long double error = fabsl(analytical[i] - numerical[i]);
        if (error > maxError)
        {
            maxError = error;
        }
    }
    return maxError;
}
