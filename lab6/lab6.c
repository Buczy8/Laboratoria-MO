#include <stdio.h>

#define N 5

void printArray(const long double arr[], int size) {
    for (int i = 0; i < size; i++) {
        printf("x[%d] = %.2Lf\n", i, arr[i]);
    }
    printf("\n");
}

int readArray(FILE* inputFile, long double arr[], int size) {
    for (int i = 0; i < size; i++) {
        if (fscanf(inputFile, "%Lf", &arr[i]) != 1) {
            fprintf(stderr, "Błąd odczytu danych!\n");
            return 0;
        }
    }
    return 1;
}

void thomas1(long double U[], long double D[], long double L[]) {
    for (int i = 1; i < N; i++) {
        D[i] = D[i] - (L[i - 1] * (1 / D[i - 1]) * U[i - 1]);
    }
}

void thomas2(long double U[], long double D[], long double L[], long double b[]) {
    for (int i = 1; i < N; i++) {
        b[i] = b[i] - (L[i - 1] * (1 / D[i - 1]) * b[i - 1]);
    }
    b[N - 1] = b[N - 1] / D[N - 1];

    for (int i = N - 2; i >= 0; i--) {
        b[i] = (b[i] - U[i] * b[i + 1]) / D[i];
    }
}

int main() {
    long double U[N - 1], D[N], L[N - 1], b[N];

    FILE* inputFile = fopen("dane.txt", "r");
    if (!inputFile) {
        fprintf(stderr, "Błąd otwarcia pliku!\n");
        return 1;
    }

    if (!readArray(inputFile, U, N - 1) ||
        !readArray(inputFile, D, N) ||
        !readArray(inputFile, L, N - 1) ||
        !readArray(inputFile, b, N)) {
        fclose(inputFile);
        return 1;
    }

    fclose(inputFile);

    thomas1(U, D, L);
    thomas2(U, D, L, b);

    printArray(b, N);

    return 0;
}
