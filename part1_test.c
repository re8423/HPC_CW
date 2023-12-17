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


int main(int argc, char** argv){

// double t = 0.0, nrmu, nrmv;
// double u[N][N], v[N][N], du[N][N], dv[N][N];

int u[N][N], v[N][N];
int b = N;
double ans = 0;
init(u, v);


#pragma omp parallel shared( u, b, v )
for (int k = 0; k < M; k++){
    funcA(u,b,v);
    funcB(u,b,v);
    
    if (i%m == 0){
        ans = funcC(u,b,v);
        ans = funcC(u,b,v);

        // printf("t = %2.1d\tv-norm = %2.5f\n", i, ans);
    }
}
printf("%d",ans);
return 0;
}
