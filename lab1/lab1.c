#include <stdio.h>
#include <float.h>


int main(void) {
    printf("%d\n", FLT_DIG);
    float epsylon_float = 1.f;
    int t_float = 0;
    int l_float = 0;
    while (1.f + epsylon_float / 2.f > 1.f){
        epsylon_float = epsylon_float / 2.f;
        t_float++;
    }
    printf("float: \n");
    printf("Liczba bitow mantysy bez najbardziej znaczacego bitu = %d\n", t_float);
    printf("epsylon maszynowy = %.*f\n", t_float,epsylon_float);
    while (epsylon_float * 10 < 1.f){
        epsylon_float = epsylon_float * 10;
        l_float++;
    }
    printf("Liczba cyfr znaczacych epsylona maszynowego = %d\n", t_float-l_float);

    double epsylon_double = 1.0;
    int t_double = 0;
    int l_double = 0;
    while(1.0 + epsylon_double / 2.0 > 1.0) {
        epsylon_double = epsylon_double / 2.0;
        t_double++;
    }
    printf("\ndouble: \n");
    printf("Liczba bitow mantysy bez najbardziej znaczacego bitu = %d\n", t_double);
    printf("epsylon maszynowy = %.*f\n", t_double,epsylon_double);
    while (epsylon_double * 10 < 1.0){
        epsylon_double = epsylon_double * 10;
        l_double++;
    }
    printf("Liczba cyfr znaczacych epsylona maszynowego = %d\n", t_double-l_double);

    long double epsylon_long_double = 1.L;
    int t_long_double = 0;
    int l_long_double = 0;
    while(1.l + epsylon_long_double / 2.0 > 1.l) {
        epsylon_long_double = epsylon_long_double / 2.l;
        t_long_double++;
    }
    printf("\nlong double: \n");
    printf("Liczba bitow mantysy bez najbardziej znaczacego bitu = %d\n", t_long_double);
    printf("epsylon maszynowy = %.*Lf\n", t_long_double,epsylon_long_double);
    while (epsylon_long_double * 10 < 1.l){
        epsylon_long_double = epsylon_long_double * 10;
        l_long_double++;
    }
    printf("Liczba cyfr znaczacych epsylona maszynowego = %d\n", t_long_double-l_long_double);

    return 0;
}
//epsylon =2*ni=2*2^-(t+1)= 2^(-t)