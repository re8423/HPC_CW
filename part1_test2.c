#include <stdio.h>				// needed for printing
#include <math.h>				// needed for tanh, used in init function
#include "params.h"				// model & simulation parameters
#include <omp.h> //openmp header file

void funcA( int a[N][N], int b, int c[N][N] ) {
    #pragma omp for schedule( static )
    for (int ii = 0; ii < b; ii++) {
        for (int jj = 0; jj < b; jj++) {
          a[ii][jj] += 1;
          c[ii][jj] += 1;
        }
    }
}

void funcB( int a[N][N], int b, int c[N][N] ) {
    #pragma omp for schedule( static )
    for (int ii = 0; ii < b; ii++) {
        for (int jj = 0; jj < b; jj++) {
          a[ii][jj] += 1;
          c[ii][jj] += 1;
        }
    }
}

// void funcC( int a[N][N], int b, int c[N][N] ) {
//     #pragma omp for schedule( static )
//     for (int ii = 0; ii < b; ii++) {
//         for (int jj = 0; jj < b; jj++) {
//           a[ii][jj] += 1;
//           c[ii][jj] += 1;
//         }
//     }
// }
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

int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}

int main(int argc, char** argv){
double ans = 0;
// omp_set_num_threads(4);
int a[N][N], c[N][N];
int b=N;
printf("This program uses %d threads.\n", omp_thread_count());

#pragma omp parallel shared( a, b, c )

for (int i = 0; i < M; i++){
    funcA(a,b,c);
    funcB(a,b,c);
    ans = funcC(a,b,c);
    if (i%m == 0){
        funcC(a,b,c);
        printf("%d\t%f\n",i, ans);
    }
}
// printf("%d", ans);


return 0;
}
