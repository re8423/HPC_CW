#include <stdio.h>				// needed for printing
#include <math.h>				// needed for tanh, used in init function
#include "params.h"				// model & simulation parameters
#include <omp.h> //openmp header file

void funcA( int a[N][N], int b, int c[N][N] ) {
    #pragma omp for schedule( static )
    for (int ii = 0; ii < b; ii++) {
        for (int jj = 0; jj < b; jj++) {
          a = a + 1;
          c = c + 1;
        }
    }
}

void funcB( int a[N][N], int b, int c[N][N] ) {
    #pragma omp for schedule( static )
    for (int ii = 0; ii < b; ii++) {
        for (int jj = 0; jj < b; jj++) {
          a = a + 2;
          c = c + 2;
        }
    }
}


double funcC (int a[N][N], int b, int c[N][N]){
    double k = 0;
    #pragma omp parallel for shared(a,b,c) reduction(+:k)
    for (int ii = 0; ii < b; ii++){
        for (int jj = 0; jj < b; jj++){
          // alter values of a and c
            k += a[ii][jj]*c[ii][jj];
    }
    }
    return k;
}


int main(int argc, char** argv){

int a[N][N], c[N][N];
int b = N;
double ans = 0;

#pragma omp parallel shared( a, b, c )
for (int i = 0; i < 4000; i++){
    funcA(a,b,c);
    funcB(a,b,c);
    ans = funcC(a,b,c);
}
printf("%d",ans);
return 0;
}
