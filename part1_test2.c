#include <stdio.h>				// needed for printing
#include <math.h>				// needed for tanh, used in init function
#include "params.h"				// model & simulation parameters
#include <omp.h> //openmp header file

void init(double u[N][N], double v[N][N]){
	double uhi, ulo, vhi, vlo;
	uhi = 0.5; ulo = -0.5; vhi = 0.1; vlo = -0.1; // these are shared vars

	//collapse wouldnt make difference if using 128 threads
	#pragma omp parrallel for schedule( static )
	for (int i=0; i < N; i++){ //vars declared in loop are private
		for (int j=0; j < N; j++){
			// u[i][j] = ulo + (uhi-ulo)*0.5*(1.0 + tanh((i-N/2)/16.0)); 
			// v[i][j] = vlo + (vhi-vlo)*0.5*(1.0 + tanh((j-N/2)/16.0));
			u[i][j] = ulo + (uhi-ulo)*0.5*(1.0 + tanh((i-N/2)*0.0625)); 
			v[i][j] = vlo + (vhi-vlo)*0.5*(1.0 + tanh((j-N/2)*0.0625));
		}
	}
}
void funcA( double u[N][N], int b, double v[N][N] ) {
    double lapu, lapv;
	int up, down, left, right;

    #pragma omp for schedule( static )
    for (int ii = 0; ii < b; ii++) {
        for (int jj = 0; jj < b; jj++) {
            if (i == 0){
                down = i;
            }
            else{
                down = i-1;
            }
            if (i == N-1){
                up = i;
            }
            else{
                up = i+1;
            }

            if (j == 0){
                left = j;
            }
            else{
                left = j-1;
            }

            if (j == N-1){
                right = j;
            }
            else{
                right = j+1;
            }
        lapu = u[up][j] + u[down][j] + u[i][left] + u[i][right] + -4.0*u[i][j];
        lapv = v[up][j] + v[down][j] + v[i][left] + v[i][right] + -4.0*v[i][j];
        du[i][j] = DD*lapu + u[i][j]*(1.0 - u[i][j])*(u[i][j]-b) - v[i][j];
        dv[i][j] = d*DD*lapv + c*(a*u[i][j] - v[i][j]);	
        }
    }
}

void funcB( double u[N][N], int b, double v[N][N] ) {
    #pragma omp for schedule( static )
    for (int i = 0; i < b; i++) {
        for (int j = 0; j < b; j++) {
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
double u[N][N], v[N][N], du[N][N], dv[N][N];
int b=N;
printf("This program uses %d threads.\n", omp_thread_count());


init(u, v);
#pragma omp parallel shared( u, b, v )

for (int i = 0; i < M; i++){
    funcA(u,b,v);
    funcB(u,b,v);
    if (i%m == 0){
        ans = funcC(u);
        ans = funcC(v);
        printf("%d\t%f\n",i, ans);
    }
}
// printf("%d", ans);


return 0;
}
