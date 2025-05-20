#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "CALERF.h"

using namespace std;

long double tmax = 2.0L;
long double D = 1.0L;
long double b = 1.0L;
long double a = ceill(b + 6.0L * sqrtl(D * tmax));
long double dt = 0.001L;

struct Result
{
    vector<long double> x;
    vector<vector<long double>> numerical; // Wszystkie kroki czasowe - rozwiązanie numeryczne
    vector<vector<long double>> analytical; // Wszystkie kroki czasowe - rozwiązanie analityczne
};

void thomas1(vector<long double>& U, vector<long double>& D, vector<long double>& L, int N);
void thomas2(vector<long double>& U, vector<long double>& D, vector<long double>& L, vector<long double>& b, int N);
std::vector<long double> generate_log_steps(long double h_max, long double h_min, int num_steps);
long double maxAbsoluteError(vector<long double>& analytical, vector<long double>& numerical);
void LUDecomposition(vector<vector<long double>>& A, vector<int>& index);
void solve(const vector<vector<long double>>& A, vector<long double>& b, const vector<int>& index);

Result computeExplicit(long double h)
{
    Result result;
    long double lambda = D * dt / (h * h);
    cout << "lambda KMB = " << lambda << endl;

    int nx = static_cast<int>((2.0L * a) / h) + 1;
    int nt = static_cast<int>(tmax / dt) + 1;

    result.x.resize(nx);
    result.numerical.resize(nx, vector<long double>(nt, 0.0L));
    result.analytical.resize(nx, vector<long double>(nt, 0.0L));

    // Inicjalizacja siatki i warunku początkowego
    for (int i = 0; i < nx; ++i)
    {
        result.x[i] = -a + i * h;
        if (fabsl(result.x[i]) < b)
            result.numerical[i][0] = 1.0L;
    }

    // Obliczenia czasowe dla metody explicit
    for (int k = 1; k < nt; ++k)
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            result.numerical[i][k] = lambda * result.numerical[i - 1][k - 1]
                + (1.0L - 2.0L * lambda) * result.numerical[i][k - 1]
                + lambda * result.numerical[i + 1][k - 1];
        }
        result.numerical[0][k] = 0.0L;
        result.numerical[nx - 1][k] = 0.0L;
    }

    // Obliczenie rozwiązania analitycznego dla wszystkich czasów
    for (int k = 0; k < nt; ++k)
    {
        long double t = k * dt;
        for (int i = 0; i < nx; ++i)
        {
            if (t > 0.0L)
            {
                result.analytical[i][k] = 0.5L * calerfpack::erf_LD((result.x[i] + b) / (2.0L * sqrtl(D * t)))
                    - 0.5L * calerfpack::erf_LD((result.x[i] - b) / (2.0L * sqrtl(D * t)));
            }
            else
            {
                result.analytical[i][k] = 0.0L;
            }
        }
    }

    return result;
}

Result computeLaasonenT(long double h)
{
    Result result;
    long double lambda = D * dt / (h * h);
    cout << "lambda laasonen thomas = " << lambda << endl;
    int nx = static_cast<int>((2.0L * a) / h) + 1;
    int nt = static_cast<int>(tmax / dt) + 1;

    result.x.resize(nx);
    result.numerical.resize(nx, vector<long double>(nt, 0.0L));
    result.analytical.resize(nx, vector<long double>(nt, 0.0L));

    vector<long double> L(nx, -lambda);
    vector<long double> Di(nx, 1.0L + 2.0L * lambda);
    vector<long double> U(nx, -lambda);
    vector<long double> sol(nx);

    // Warunki brzegowe w macierzy
    Di[0] = 1.0L;
    U[0] = 0.0L;
    L[0] = 0.0L;

    Di[nx - 1] = 1.0L;
    U[nx - 1] = 0.0L;
    L[nx - 1] = 0.0L;

    // Inicjalizacja warunku początkowego
    for (int i = 0; i < nx; ++i)
    {
        result.x[i] = -a + i * h;
        if (fabsl(result.x[i]) < b)
            result.numerical[i][0] = 1.0L;
    }

    thomas1(U, Di, L, nx);

    // Obliczenia czasowe
    for (int k = 1; k < nt; ++k)
    {
        for (int i = 0; i < nx; ++i)
            sol[i] = result.numerical[i][k - 1];

        sol[0] = 0.0L;
        sol[nx - 1] = 0.0L;

        thomas2(U, Di, L, sol, nx);

        for (int i = 0; i < nx; ++i)
            result.numerical[i][k] = sol[i];
    }

    // Obliczenie rozwiązania analitycznego dla wszystkich czasów
    for (int k = 0; k < nt; ++k)
    {
        long double t = k * dt;
        for (int i = 0; i < nx; ++i)
        {
            result.analytical[i][k] = 0.5L * calerfpack::erf_LD((result.x[i] + b) / (2.0L * sqrtl(D * t)))
                - 0.5L * calerfpack::erf_LD((result.x[i] - b) / (2.0L * sqrtl(D * t)));
        }
    }

    return result;
}

Result computeLaasonenLU(long double h)
{
    Result result;
    long double lambda = D * dt / (h * h);
    cout << "lambda  lu = " << lambda << endl;
    int nx = static_cast<int>((2.0L * a) / h) + 1;
    int nt = static_cast<int>(tmax / dt) + 1;

    result.x.resize(nx);
    result.numerical.resize(nx, vector<long double>(nt, 0.0L));
    result.analytical.resize(nx, vector<long double>(nt, 0.0L));

    vector<vector<long double>> A(nx, vector<long double>(nx, 0.0L));
    vector<int> index(nx);
    vector<long double> sol(nx);
    for (int i = 0; i < nx; ++i)
    {
        index[i] = i;
    }

    // Budowa macierzy A
    for (int i = 0; i < nx; ++i)
    {
        if (i > 0)
            A[i][i - 1] = -lambda;
        A[i][i] = 1.0L + 2.0L * lambda;
        if (i < nx - 1)
            A[i][i + 1] = -lambda;
    }

    // Warunki brzegowe
    A[0][0] = 1.0L;
    A[0][1] = 0.0L;
    A[nx - 1][nx - 1] = 1.0L;
    A[nx - 1][nx - 2] = 0.0L;

    // Inicjalizacja warunku początkowego
    for (int i = 0; i < nx; ++i)
    {
        result.x[i] = -a + i * h;
        if (fabsl(result.x[i]) < b)
            result.numerical[i][0] = 1.0L;
    }

    // Dekompozycja LU
    LUDecomposition(A, index);

    // Obliczenia czasowe
    for (int k = 1; k < nt; ++k)
    {
        for (int i = 0; i < nx; ++i)
            sol[i] = result.numerical[i][k - 1];

        sol[0] = 0.0L;
        sol[nx - 1] = 0.0L;

        solve(A, sol, index);

        for (int i = 0; i < nx; ++i)
            result.numerical[i][k] = sol[i];
    }

    // Obliczenie rozwiązania analitycznego
    for (int k = 0; k < nt; ++k)
    {
        long double t = k * dt;
        for (int i = 0; i < nx; ++i)
        {
            if (t > 0.0L)
            {
                result.analytical[i][k] = 0.5L * erf((result.x[i] + b) / (2.0L * sqrtl(D * t)))
                    - 0.5L * erf((result.x[i] - b) / (2.0L * sqrtl(D * t)));
            }
            else
            {
                result.analytical[i][k] = 0.0L;
            }
        }
    }

    return result;
}

int main()
{
    long double h_explicit = 0.05L;
    long double h_laasonen = sqrt(10.0L) / 100.0L;
    int nt = static_cast<int>(tmax / dt) + 1;

    Result explicitResult = computeExplicit(h_explicit);
    Result laasonenResultT = computeLaasonenT(h_laasonen);
    Result laasonenResultLU = computeLaasonenLU(h_laasonen);

    // Zapis wyników explicit
    ofstream fp2("KMB.txt");
    fp2 << scientific << setprecision(15);

    for (int i = 0; i < explicitResult.x.size(); ++i)
    {
        fp2 << explicitResult.x[i] << " "
            << explicitResult.numerical[i][nt - 1] << " "
            << explicitResult.analytical[i][nt - 1] << "\n";
    }
    fp2.close();

    // Zapis wyników Laasonen
    ofstream fp1("laasonen.txt");
    fp1 << scientific << setprecision(15);

    for (int i = 0; i < laasonenResultLU.x.size(); ++i)
    {
        fp1 << laasonenResultLU.x[i] << " "
            << laasonenResultLU.numerical[i][nt - 1] << " "
            << laasonenResultT.numerical[i][nt - 1] << " "
            << laasonenResultLU.analytical[i][nt - 1] << "\n";
    }
    fp1.close();

    return 0;
}

void thomas1(vector<long double>& U, vector<long double>& D, vector<long double>& L, int N)
{
    for (int i = 1; i < N; i++)
    {
        D[i] = D[i] - (L[i - 1] * (1.0L / D[i - 1]) * U[i - 1]);
    }
}

void thomas2(vector<long double>& U, vector<long double>& D, vector<long double>& L, vector<long double>& b, int N)
{
    for (int i = 1; i < N; i++)
    {
        b[i] = b[i] - (L[i - 1] * (1.0L / D[i - 1]) * b[i - 1]);
    }
    b[N - 1] = b[N - 1] / D[N - 1];

    for (int i = N - 2; i >= 0; i--)
    {
        b[i] = (b[i] - U[i] * b[i + 1]) / D[i];
    }
}

std::vector<long double> generate_log_steps(long double h_max, long double h_min, int num_steps)
{
    std::vector<long double> steps;
    steps.reserve(num_steps);

    long double ratio = pow(h_min / h_max, 1.0 / (num_steps - 1));
    long double current_h = h_max;

    for (int i = 0; i < num_steps; ++i)
    {
        steps.push_back(current_h);
        current_h *= ratio;
    }

    return steps;
}

long double maxAbsoluteError(vector<long double>& analytical, vector<long double>& numerical)
{
    long double maxError = 0.0L;
    size_t size = analytical.size();

    for (size_t i = 0; i < size; ++i)
    {
        long double error = fabsl(analytical[i] - numerical[i]);
        if (error > maxError)
        {
            maxError = error;
        }
    }

    return maxError;
}

void LUDecomposition(vector<vector<long double>>& A, vector<int>& index)
{
    int N = (int)A.size();
    index.resize(N);
    for (int i = 0; i < N; ++i) index[i] = i;

    int i, j, k = 0;
    for (k = 0; k < N - 1; k++)
    {
        for (i = k + 1; i < N; i++)
        {
            if (fabsl(A[index[k]][k]) <= 10e-10)
            {
                int maxIndex = k;
                for (j = k + 1; j < N; j++)
                {
                    if (fabsl(A[index[j]][k]) > fabsl(A[index[maxIndex]][k]))
                    {
                        maxIndex = j;
                    }
                }
                if (maxIndex != k)
                {
                    int tempIndex = index[k];
                    index[k] = index[maxIndex];
                    index[maxIndex] = tempIndex;
                }
            }
            long double factor = A[index[i]][k] / A[index[k]][k];
            for (j = k; j < N; j++)
            {
                A[index[i]][j] -= factor * A[index[k]][j];
            }
            A[index[i]][k] = factor;
        }
    }
}


void solve(const vector<vector<long double>>& A, vector<long double>& b, const vector<int>& index)
{
    int N = (int)A.size();

    //      y-1 = b_1
    //      l_21*y_1 + y_2 = b_2
    //      l_n1*y_1 + l_n2*y_2 + ... + y_n = b_n

    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            b[index[i]] -= A[index[i]][j] * b[index[j]];
        }
    }

    //    u_11*x_1 + u_12*x_2 + ... + x_n = y_1
    //    u_22*x_2 + ... + x_n = y_2
    //    u_nn*x_n = y_n
    for (int i = N - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
        {
            b[index[i]] -= A[index[i]][j] * b[index[j]];
        }
        b[index[i]] /= A[index[i]][i];
    }
}
