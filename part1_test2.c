#include <stdio.h>				// needed for printing
#include <math.h>				// needed for tanh, used in init function
#include "params.h"				// model & simulation parameters
#include <omp.h> //openmp header file

void funcA( double u[N][N], int b, double v[N][N] ) {
    #pragma omp for schedule( static )
    for (int ii = 0; ii < b; ii++) {
        for (int jj = 0; jj < b; jj++) {
          u[ii][jj] += 1;
          v[ii][jj] += 1;
        }
    }
}

void funcB( double u[N][N], int b, double v[N][N] ) {
    #pragma omp for schedule( static )
    for (int ii = 0; ii < b; ii++) {
        for (int jj = 0; jj < b; jj++) {
          u[ii][jj] += 1;
          v[ii][jj] += 1;
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
double funcC (double x[N][N]){
    double k = 0;
    #pragma omp parallel for reduction(+:k)
    for (int ii = 0; ii < N; ii++){
        for (int jj = 0; jj < N; jj++){
          // alter values of a and c
            k += x[ii][jj]*x[ii][jj];
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
double u[N][N], v[N][N];
int b=N;
printf("This program uses %d threads.\n", omp_thread_count());

#pragma omp parallel shared( u, b, v )

for (int i = 0; i < M; i++){
    funcA(u,b,v);
    funcB(u,b,v);
    ans = funcC(u,b,v);
    if (i%m == 0){
        funcC(u,b,v);
        printf("%d\t%f\n",i, ans);
    }
}
// printf("%d", ans);


return 0;
}
