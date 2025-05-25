#include <iostream> // bibliotek do wejścia/wyjścia
#include <fstream> // bibliotek do obsługi plików
#include <vector> // bibliotek do przechowywania wektorów
#include <cmath> // bibliotek do funkcji matematycznych
#include <iomanip> // bibliotek do formatowania wyjścia
#include "CALERF.h" // bibliotek do funkcji błędu

using namespace std;

// Parametry podane w zadaniu
long double tmax = 2.0L; // maksymalny czas
long double D = 1.0L; // współczynnik dyfuzji
long double b = 1.0L; // grubość warstwy 
long double a = ceill(b + 6.0L * sqrtl(D * tmax)); // długość obszaru, w którym rozważamy dyfuzję
long double dt = 0.001L; // krok czasowy

// Struktura do przechowywania wyników
struct Result
{
    vector<long double> x; // Współrzędne przestrzenne
    vector<vector<long double>> numerical; // Wszystkie kroki czasowe - rozwiązanie numeryczne
    vector<vector<long double>> analytical; // Wszystkie kroki czasowe - rozwiązanie analityczne
};

// Deklaracje pomocniczych funkcji
void thomas1(vector<long double>& U, vector<long double>& D, vector<long double>& L, int N);
void thomas2(vector<long double>& U, vector<long double>& D, vector<long double>& L, vector<long double>& b, int N);
void LUDecomposition(vector<vector<long double>>& A);
void solve(const vector<vector<long double>>& A, vector<long double>& b);

// Funkcja do obliczeń numerycznych dla metody KMB
Result computeKMB(long double h)
{
    Result result;
    long double lambda = D * dt / (h * h); // Parametr lambda potrzebny dla metody KMB

    int nx = static_cast<int>((2.0L * a) / h) + 1; // Liczba punktów w przestrzeni
    int nt = static_cast<int>(tmax / dt) + 1; // Liczba kroków czasowych

    result.x.resize(nx);
    result.numerical.resize(nx, vector<long double>(nt, 0.0L));
    result.analytical.resize(nx, vector<long double>(nt, 0.0L));

    // Inicjalizacja siatki i warunku początkowego
    for (int i = 0; i < nx; ++i)
    {
        result.x[i] = -a + i * h; // Obliczenie współrzędnych przestrzennych

        // Ustawienie warunku początkowego
        if (fabsl(result.x[i]) < b)
            result.numerical[i][0] = 1.0L;
    }

    // Obliczenia czasowe dla metody KMB
    for (int k = 1; k < nt; ++k)
    {
        for (int i = 1; i < nx - 1; ++i)
        {
            // Obliczenie wartości w punkcie i dla kroku czasowego k według metody KMB
            result.numerical[i][k] = lambda * result.numerical[i - 1][k - 1]
                + (1.0L - 2.0L * lambda) * result.numerical[i][k - 1]
                + lambda * result.numerical[i + 1][k - 1];
        }
        // Warunki brzegowe
        result.numerical[0][k] = 0.0L;
        result.numerical[nx - 1][k] = 0.0L;
    }

    // Obliczenie rozwiązania analitycznego dla wszystkich czasów
    for (int k = 0; k < nt; ++k)
    {
        long double t = k * dt; // Obliczenie czasu dla kroku k
        for (int i = 0; i < nx; ++i)
        {
            result.analytical[i][k] = 0.5L * calerfpack::erf_LD((result.x[i] + b) / (2.0L * sqrtl(D * t)))
                - 0.5L * calerfpack::erf_LD((result.x[i] - b) / (2.0L * sqrtl(D * t)));
            // Obliczenie rozwiązania analitycznego za pommocą pakietu CALERF dołączonego przez prowadząćego
        }
    }

    return result;
}

// Funkcja do obliczeń numerycznych dla metody Laasonen z wykorzystaniem metody Thomasa
Result computeLaasonenT(long double h)
{
    Result result;
    long double lambda = D * dt / (h * h); // Parametr lambda potrzebny dla metody Laasonen
    int nx = static_cast<int>((2.0L * a) / h) + 1; // Liczba punktów w przestrzeni
    int nt = static_cast<int>(tmax / dt) + 1; // Liczba kroków czasowych

    result.x.resize(nx);
    result.numerical.resize(nx, vector<long double>(nt, 0.0L));
    result.analytical.resize(nx, vector<long double>(nt, 0.0L));

    vector<long double> L(nx, -lambda); // Dolna przekątna macierzy
    vector<long double> Di(nx, 1.0L + 2.0L * lambda); // Główna przekątna macierzy
    vector<long double> U(nx, -lambda); // Górna przekątna macierzy
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
        result.x[i] = -a + i * h; // Obliczenie współrzędnych przestrzennych

        // Ustawienie warunku początkowego
        if (fabsl(result.x[i]) < b)
            result.numerical[i][0] = 1.0L;
    }

    thomas1(U, Di, L, nx); // Dekompozycja macierzy dla metody Thomasa

    // Obliczenia czasowe
    for (int k = 1; k < nt; ++k)
    {
        // Przygotowanie wektora sol do obliczeń
        for (int i = 0; i < nx; ++i)
            sol[i] = result.numerical[i][k - 1];

        thomas2(U, Di, L, sol, nx); // Rozwiązanie układu równań za pomocą metody Thomasa

        // Ustawienie warunków brzegowych
        sol[0] = 0.0L;
        sol[nx - 1] = 0.0L;

        // Zapisanie rozwiązania numerycznego dla określonego czasu
        for (int i = 0; i < nx; ++i)
            result.numerical[i][k] = sol[i];
    }

    // Obliczenie rozwiązania analitycznego dla wszystkich czasów
    for (int k = 0; k < nt; ++k)
    {
        long double t = k * dt; // Obliczenie czasu dla kroku k
        for (int i = 0; i < nx; ++i)
        {
            result.analytical[i][k] = 0.5L * calerfpack::erf_LD((result.x[i] + b) / (2.0L * sqrtl(D * t)))
                - 0.5L * calerfpack::erf_LD((result.x[i] - b) / (2.0L * sqrtl(D * t)));
            // Obliczenie rozwiązania analitycznego za pommocą pakietu CALERF dołączonego przez prowadząćego
        }
    }

    return result;
}

// Funkcja do obliczeń numerycznych dla metody Laasonen z wykorzystaniem metody dekompozycji LU macierzy pełnej
Result computeLaasonenLU(long double h)
{
    Result result;
    long double lambda = D * dt / (h * h); // Parametr lambda potrzebny dla metody Laasonen
    int nx = static_cast<int>((2.0L * a) / h) + 1; // Liczba punktów w przestrzeni
    int nt = static_cast<int>(tmax / dt) + 1; // Liczba kroków czasowych

    result.x.resize(nx);
    result.numerical.resize(nx, vector<long double>(nt, 0.0L));
    result.analytical.resize(nx, vector<long double>(nt, 0.0L));

    vector<vector<long double>> A(nx, vector<long double>(nx, 0.0L)); // Macierz A dla dekompozycji LU
    vector<long double> sol(nx);

    // Budowa macierzy A
    for (int i = 0; i < nx; ++i)
    {
        // Wypełnienie macierzy A zgodnie z równaniem różnicowym Laasonen
        if (i > 0)
            A[i][i - 1] = -lambda; // element poniżej głównej przekątnej
        A[i][i] = 1.0L + 2.0L * lambda; // element na głównej przekątnej
        if (i < nx - 1)
            A[i][i + 1] = -lambda; // element powyżej głównej przekątnej
    }

    // Warunki brzegowe
    A[0][0] = 1.0L;
    A[0][1] = 0.0L;
    A[nx - 1][nx - 1] = 1.0L;
    A[nx - 1][nx - 2] = 0.0L;

    // Inicjalizacja warunku początkowego
    for (int i = 0; i < nx; ++i)
    {
        result.x[i] = -a + i * h; // Obliczenie współrzędnych przestrzennych

        if (fabsl(result.x[i]) < b)
            result.numerical[i][0] = 1.0L;
    }

    // Dekompozycja LU
    LUDecomposition(A);

    // Obliczenia czasowe
    for (int k = 1; k < nt; ++k)
    {
        // Przygotowanie wektora sol z poprzedniego kroku czasowego
        for (int i = 0; i < nx; ++i)
            sol[i] = result.numerical[i][k - 1];

        // Ustawienie warunków brzegowych
        sol[0] = 0.0L;
        sol[nx - 1] = 0.0L;

        // Rozwiązanie układu równań Ax = b
        solve(A, sol);

        for (int i = 0; i < nx; ++i)
            result.numerical[i][k] = sol[i];
    }

    // Obliczenie rozwiązania analitycznego
    for (int k = 0; k < nt; ++k)
    {
        long double t = k * dt; // Obliczenie czasu dla kroku k

        // Obliczenie rozwiązania analitycznego dla wszystkich punktów przestrzennych
        for (int i = 0; i < nx; ++i)
        {
            result.analytical[i][k] = 0.5L * erf((result.x[i] + b) / (2.0L * sqrtl(D * t)))
                - 0.5L * erf((result.x[i] - b) / (2.0L * sqrtl(D * t)));
        }
    }

    return result;
}

int main()
{
    long double h_KMB = 0.05L;
    // Krok przestrzenny dla metody KMB dobrany tak żeby spełniał warunek stabilności (lambda = 0.4)
    long double h_laasonen = sqrt(10.0L) / 100.0L; // Krok przestrzenny dla metody Laasonen  dobrany tak aby lambda = 1
    int nt = static_cast<int>(tmax / dt) + 1; // Liczba kroków czasowych

    Result KMBResult = computeKMB(h_KMB);
    Result laasonenResultT = computeLaasonenT(h_laasonen);
    Result laasonenResultLU = computeLaasonenLU(h_laasonen);

    // Zapis wyników KMB
    ofstream fp2("KMB.txt");
    fp2 << scientific << setprecision(15); // Ustawienie precyzji zapisu do pliku

    // Zapis współrzędnych przestrzennych i wartości numerycznych oraz analitycznych dla metody KMB
    for (int i = 0; i < KMBResult.x.size(); ++i)
    {
        fp2 << KMBResult.x[i] << " "
            << KMBResult.numerical[i][nt - 1] << " "
            << KMBResult.analytical[i][nt - 1] << "\n";
    }
    fp2.close();

    // Zapis wyników Laasonen
    ofstream fp1("laasonen.txt");
    fp1 << scientific << setprecision(15);

    // Zapis współrzędnych przestrzennych i wartości numerycznych oraz analitycznych dla metody Laasonen
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

// eliminacja Gaussa która przekształca tylko jeden wiersz w każdym kroku redukcji
void thomas1(vector<long double>& U, vector<long double>& D, vector<long double>& L, int N)
{
    for (int i = 1; i < N; i++)
    {
        D[i] = D[i] - (L[i - 1] * (1.0L / D[i - 1]) * U[i - 1]);
    }
}

void thomas2(vector<long double>& U, vector<long double>& D, vector<long double>& L, vector<long double>& b, int N)
{
    // Eliminacja w przód – modyfikacja prawej strony
    for (int i = 1; i < N; i++)
    {
        b[i] = b[i] - (L[i - 1] * (1.0L / D[i - 1]) * b[i - 1]);
    }

    // Rozwiązanie ostatniego równania (najniższego w macierzy)
    b[N - 1] = b[N - 1] / D[N - 1];

    // Podstawianie wstecz – obliczanie pozostałych niewiadomych
    for (int i = N - 2; i >= 0; i--)
    {
        b[i] = (b[i] - U[i] * b[i + 1]) / D[i];
    }
}

// Faktoryzacja LU macierzy A
void LUDecomposition(vector<vector<long double>>& A)
{
    int N = (int)A.size();

    int i, j, k = 0;
    for (k = 0; k < N - 1; k++)
    {
        // Eliminacja Gaussa macierzy A
        for (i = k + 1; i < N; i++)
        {
            long double factor = A[i][k] / A[k][k]; // współczynnik eliminacji
            for (j = k; j < N; j++)
            {
                A[i][j] -= factor * A[k][j]; // modyfikacja macierzy powyżej głównej przekątnej
            }
            A[i][k] = factor; // zapis współczynnika do macierzy poniżej głównej przekątnej
        }
    }
}

// Rozwiązywanie układu równań Ax = b przy użyciu wcześniej wykonanej faktoryzacji LU
void solve(const vector<vector<long double>>& A, vector<long double>& b)
{
    int N = (int)A.size();

    // Rozwiązanie układu Ly = b (eliminacja w przód)
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            b[i] -= A[i][j] * b[j];
        }
    }

    // Rozwiązanie układu Ux = y (podstawianie wstecz)
    for (int i = N - 1; i >= 0; i--)
    {
        for (int j = i + 1; j < N; j++)
        {
            b[i] -= A[i][j] * b[j];
        }
        b[i] /= A[i][i];
    }
}
