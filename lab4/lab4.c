#include <stdio.h>
#include <math.h>


void f(const long double x[3], long double fx[3]){
    fx[0] = x[0] * x[0] + x[1] * x[1] + x[2] * x[2] - 4;
    fx[1] = x[0] * x[0] + (x[1] * x[1])/2 - 1;
    fx[2] = x[0] * x[1] - 0.5;
}
void jacobi(const long double x[3], long double J[3][3]){
    J[0][0] = 2 * x[0];
    J[0][1] = 2 * x[1];
    J[0][2] = 2 * x[2];

    J[1][0] = 2 * x[0];
    J[1][1] = x[1];
    J[1][2] = 0;

    J[2][0] = x[1];
    J[2][1] = x[0];
    J[2][2] = 0;
}
void rowanie_liniowe(long double J[3][3], long double fx[3], long double delta[3]) {
    int i, j, k, maxIndex;
    long double ratio, maxElement, temp;

    long double A[3][4];
    for (i = 0; i < 3; i++) {
        for (j = 0; j < 3; j++) {
            A[i][j] = J[i][j];
        }
        A[i][3] = fx[i];
    }
    //czesciowa zamiana elementu podstawowego
    for (i = 0; i < 3; i++) {
        maxElement = fabsl(A[i][i]);
        maxIndex = i;
        for (j = i + 1; j < 3; j++) {
            if (fabsl(A[j][i]) > maxElement) {
                maxElement = fabsl(A[j][i]);
                maxIndex = j;
            }
        }
        if (maxIndex != i) {
            for (k = 0; k < 4; k++) {
                temp = A[i][k];
                A[i][k] = A[maxIndex][k];
                A[maxIndex][k] = temp;
            }
        }
    }
    //gauss
    for (i = 0; i < 3; i++) {
        for (j = i + 1; j < 3; j++) {
            ratio = A[j][i] / A[i][i];
            for (k = 0; k < 4; k++) {
                A[j][k] -= ratio * A[i][k];
            }
        }
    }
//liczenie delta_n+1
    delta[2] = A[2][3] / A[2][2];
    delta[1] = (A[1][3] - A[1][2] * delta[2]) / A[1][1];
    delta[0] = (A[0][3] - A[0][1] * delta[1] - A[0][2] * delta[2]) / A[0][0];


}

void newton(long double x0[3], long double x[3], long double tolx, long double tolf, int max_n) {
    int n = 0;
    long double fx[3], J[3][3], delta[3];
    f(x0, fx);
    x[0] = x0[0];
    x[1] = x0[1];
    x[2] = x0[2];
    while (n<max_n) {

        jacobi(x, J);
        rowanie_liniowe(J, fx, delta);
        x[0] = x[0] - delta[0];
        x[1] = x[1] - delta[1];
        x[2] = x[2] - delta[2];
        f(x, fx);

        printf("Iteracja: %d zbierzonsc: %.21Lf %.21Lf %.21Lf\n", n,x[0], x[1], x[2]);
        if(fabsl(fx[0]) < tolf && fabsl(fx[1]) < tolf && fabsl(fx[2]) < tolf) {
            if(fabsl(delta[0]) < tolx && fabsl(delta[1]) < tolx && fabsl(delta[2]) < tolx) break;
        }

        n++;

    }
}

int main(void) {
    long double tolf = 1e-15, tolx = 1e-15;
    int max_n = 100;
    //long double x0[3] = {-10, -10, -10}, x[3];
    long double x0[3] = {1, 1, 1}, x[3];

    newton(x0, x, tolx, tolf, max_n);
    return 0;
}
