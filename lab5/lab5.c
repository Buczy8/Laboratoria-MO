#include <stdio.h>
#include <math.h>

#define N 5
void printMatrix(long double A[N][N], int index[N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%.2Lf ", A[index[i]][j]);
        }
        printf("\n");
    }
}

void LU(long double A[N][N], int index[N]) {
    int i, j, k = 0;
    for (k = 0; k < N-1; k++) {
        for (i = k + 1; i < N; i++) {
            if (fabsl(A[index[k]][k]) <= 10e-10 ) {
                int maxIndex = k;
                for (j = k + 1; j < N; j++) {
                    if (fabsl(A[index[j]][k]) > fabsl(A[index[maxIndex]][k])) {
                        maxIndex = j;
                    }
                }
                if (maxIndex != k) {
                    int tempIndex = index[k];
                    index[k] = index[maxIndex];
                    index[maxIndex] = tempIndex;
                }
            }
            long double factor = A[index[i]][k] / A[index[k]][k];
            for (j = k; j < N; j++) {
                A[index[i]][j] -= factor * A[index[k]][j];
            }
            A[index[i]][k] = factor;
        }
    }
}

void rowanie(long double A[N][N], long double b[N], int index[N]) {

//      y-1 = b_1
//      l_21*y_1 + y_2 = b_2
//      l_n1*y_1 + l_n2*y_2 + ... + y_n = b_n

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < i; j++) {
            b[index[i]] -= A[index[i]][j] * b[index[j]];
        }
    }

//    u_11*x_1 + u_12*x_2 + ... + x_n = y_1
//    u_22*x_2 + ... + x_n = y_2
//    u_nn*x_n = y_n
    for (int i = N - 1; i >= 0; i--) {
        for (int j = i + 1; j < N; j++) {
            b[index[i]] -= A[index[i]][j] * b[index[j]];
        }
        b[index[i]] /= A[index[i]][i];
    }
}



int main() {
    int index[N] = {0};
    for (int i = 0; i < N; i++) {
        index[i] = i;
    }
    long double A[N][N] = {
            {5,  4, 3,  2, 1},
            {10, 8, 7,  6, 5},
            {-1, 2, -3, 4, -5},
            {6,  5, -4, 3, -2},
            {1,  2, 3,  4, 5}
    };
    long double b[N] = {37,99,-9,12,53};


    LU(A, index);
    rowanie(A, b, index);
    for (int i = 0; i < N; i++) {
        printf("x[%d] = %.2Lf\n", i, b[index[i]]);
    }
    return 0;
}
