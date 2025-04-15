#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define OMEGA 0.5L

void printMatrix(long double** A, int size);
long double** allocateMatrix(int size);
void freeMatrix(long double** A, int size);
void printVector(long double* b, int size);
void jacobi(long double** A, long double* b, const long double* x0, int N, int max_iter, long double tolx,
            long double tolf);
void Gauss_Seidel(long double** A, long double* b, const long double* x0, int N, int max_iter, long double tolx,
                  long double tolf);
void sor(long double** A, long double* b, const long double* x0, int N, int max_iter, long double tolx,
         long double tolf);
long double maxVectorNorm(long double* x, int N);
long double residualNorom(long double** A, long double* b, long double* x, int N);

int main()
{
    long double tolx = 1e-5;
    long double tolf = 1e-5;
    int n_max = 100;

    int N;
    FILE* file = fopen("dane.txt", "r");
    if (file == NULL)
    {
        printf("Nie można otworzyć pliku dane.txt\n");
        return 1;
    }

    fscanf(file, "%d", &N);

    long double** A = allocateMatrix(N);

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            fscanf(file, "%Lf", &A[i][j]);
        }
    }

    long double* b = (long double*)malloc(N * sizeof(long double));
    long double* x0 = (long double*)malloc(N * sizeof(long double));

    for (int i = 0; i < N; i++)
    {
        fscanf(file, "%Lf", &b[i]);
    }

    for (int i = 0; i < N; i++)
    {
        fscanf(file, "%Lf", &x0[i]);
    }

    fclose(file);
    jacobi(A, b, x0, N, n_max, tolx, tolf);
    Gauss_Seidel(A, b, x0, N, n_max, tolx, tolf);
    sor(A, b, x0, N, n_max, tolx, tolf);

    freeMatrix(A, N);
    free(b);
    free(x0);

    return 0;
}

void printMatrix(long double** A, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            printf("%.2Lf ", A[i][j]);
        }
        printf("\n");
    }
}

long double** allocateMatrix(int size)
{
    long double** A = (long double**)malloc(size * sizeof(long double*));
    for (int i = 0; i < size; i++)
    {
        A[i] = (long double*)malloc(size * sizeof(long double));
    }
    return A;
}

void freeMatrix(long double** A, int size)
{
    for (int i = 0; i < size; i++)
    {
        free(A[i]);
    }
    free(A);
}

void printVector(long double* b, int size)
{
    printf("[ ");
    for (int i = 0; i < size; i++)
    {
        printf("%.2Lf ", b[i]);
    }
    printf("] ");
}

long double maxVectorNorm(long double* x, int N)
{
    long double max = fabsl(x[0]);
    for (int i = 1; i < N; i++)
    {
        long double abs_val = fabsl(x[i]);
        if (abs_val > max)
        {
            max = abs_val;
        }
    }
    return max;
}

long double residualNorom(long double** A, long double* b, long double* x, int N)
{
    long double* r = (long double*)malloc(N * sizeof(long double));
    for (int i = 0; i < N; i++)
    {
        r[i] = b[i];
        for (int j = 0; j < N; j++)
        {
            r[i] -= A[i][j] * x[j];
        }
    }
    long double norm = fabsl(maxVectorNorm(r, N));
    free(r);
    return norm;
}

void jacobi(long double** A, long double* b, const long double* x0, int N, int max_iter, long double tolx,
            long double tolf)
{
    printf("Metoda Jacobiego\n");
    long double** M = allocateMatrix(N);
    long double* c = (long double*)malloc(N * sizeof(long double));

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            if (i == j)
            {
                M[i][j] = 0.0;
            }
            else
            {
                M[i][j] = A[i][j] / (-A[i][i]);
            }
        }
        c[i] = b[i] / (A[i][i]);
    }

    long double* x = (long double*)malloc(N * sizeof(long double));
    long double* x_prev = (long double*)malloc(N * sizeof(long double));
    for (int i = 0; i < N; i++)
    {
        x_prev[i] = x0[i];
    }
    for (int i = 0; i < max_iter; i++)
    {
        for (int j = 0; j < N; j++)
        {
            x[j] = c[j];
            for (int k = 0; k < N; k++)
            {
                x[j] += M[j][k] * x_prev[k];
            }
        }
        printf("Iteracja %d: ", i);
        printVector(x, N);
        long double b_x = fabsl(maxVectorNorm(x, N) - maxVectorNorm(x_prev, N));
        long double b_f = residualNorom(A, b, x, N);
        printf("Estymator bledu:  %Le ", b_x);
        printf("Reziduum: %Le\n", b_f);
        if (b_x < tolx && b_f < tolf)
        {
            break;
        }
        for (int j = 0; j < N; j++)
        {
            x_prev[j] = x[j];
        }
    }
    freeMatrix(M, N);
    free(c);
    free(x);
    free(x_prev);
}

void Gauss_Seidel(long double** A, long double* b, const long double* x0, int N, int max_iter, long double tolx,
                  long double tolf)
{
    long double* x = (long double*)malloc(N * sizeof(long double));
    long double* x_prev = (long double*)malloc(N * sizeof(long double));
    long double* rhs = (long double*)malloc(N * sizeof(long double));

    for (int i = 0; i < N; i++)
    {
        x[i] = x0[i];
    }

    printf("Metoda Gaussa-Seidela\n");

    for (int k = 0; k < max_iter; k++)
    {
        for (int i = 0; i < N; i++)
        {
            x_prev[i] = x[i];
        }
        //prawea strona: rhs = -Ux_prev + b
        for (int i = 0; i < N; i++)
        {
            rhs[i] = b[i];
            for (int j = i + 1; j < N; j++)
            {
                rhs[i] -= A[i][j] * x_prev[j];
            }
        }
        // rozwiązanie układu z macierzą trójkątną dolną
        for (int i = 0; i < N; i++)
        {
            long double sum = 0.0;
            for (int j = 0; j < i; j++)
            {
                sum += A[i][j] * x[j];
            }
            x[i] = (rhs[i] - sum) / A[i][i];
        }

        printf("Iteracja %d: ", k);
        printVector(x, N);
        long double b_x = fabsl(maxVectorNorm(x, N) - maxVectorNorm(x_prev, N));
        long double b_f = residualNorom(A, b, x, N);
        printf("Estymator bledu:  %Le ", b_x);
        printf("Reziduum: %Le\n", b_f);
        if (b_x < tolx && b_f < tolf)
        {
            break;
        }
    }

    free(x);
    free(x_prev);
    free(rhs);
}

void sor(long double** A, long double* b, const long double* x0, int N, int max_iter, long double tolx,
         long double tolf)
{
    long double* x = (long double*)malloc(N * sizeof(long double));
    long double* x_prev = (long double*)malloc(N * sizeof(long double));
    long double* rhs = (long double*)malloc(N * sizeof(long double));
    for (int i = 0; i < N; i++)
    {
        x[i] = x0[i];
    }
    printf("Metoda SOR\n");
    for (int k = 0; k < max_iter; k++)
    {
        for (int i = 0; i < N; i++)
        {
            x_prev[i] = x[i];
        }
        //prawa strona: rhs = -[(1-1/OMEGA)D + U]x_prev + b
        for (int i = 0; i < N; i++)
        {
            rhs[i] = b[i];
            for (int j = i; j < N; j++)
            {
                if (i == j)
                {
                    rhs[i] -= (1 - (1 / OMEGA)) * A[i][j] * x_prev[j];
                }
                else
                {
                    rhs[i] -= A[i][j] * x_prev[j];
                }
            }
        }
        // rozwiązanie układu z macierzą trójkątną dolną
        for (int i = 0; i < N; i++)
        {
            long double sum = 0.0;
            for (int j = 0; j < i; j++)
            {
                sum += A[i][j] * x[j];
            }
            x[i] = (rhs[i] - sum) / (A[i][i] / OMEGA);
        }
        printf("Iteracja %d: ", k);
        printVector(x, N);
        long double b_x = fabsl(maxVectorNorm(x, N) - maxVectorNorm(x_prev, N));
        long double b_f = residualNorom(A, b, x, N);
        printf("Estymator bledu:  %Le ", b_x);
        printf("Reziduum: %Le\n", b_f);
        if (b_x < tolx && b_f < tolf)
        {
            break;
        }
    }
    free(x);
    free(x_prev);
    free(rhs);
}
