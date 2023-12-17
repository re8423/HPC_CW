#include <stdio.h>				// needed for printing
#include <math.h>				// needed for tanh, used in init function
#include "params.h"				// model & simulation parameters
#include <omp.h> //openmp header file

void funcA( double a[N][N], double c[N][N] ) {
    double lapu, lapv;
	int up, down, left, right;
    #pragma omp for schedule( static )
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
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

void funcB( double a[N][N], double c[N][N] ) {
    #pragma omp for schedule( static )
    for (int ii = 0; ii < N; ii++) {
        for (int jj = 0; jj < N; jj++) {
          a = a + 2;
          c = c + 2;
        }
    }
}


double funcC (double a[N][N], double c[N][N]){
    double nrxw = 0;
    #pragma omp parallel for shared(a,c) reduction(+:nrxw)
    for (int ii = 0; ii < N; ii++){
        for (int jj = 0; jj < N; jj++){
          // alter values of a and c
            nrxw += a[ii][jj]*c[ii][jj];
    }
    }
    return nrxw;
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

double t = 0.0, nrmu, nrmv;
double u[N][N], v[N][N], du[N][N], dv[N][N];

// double u[N][N], v[N][N];
// int b = N;
// double ans = 0;
init(u, v);


#pragma omp parallel shared( u, v )
for (int k = 0; k < M; k++){
    funcA(u,v);
    funcB(u,v);
    
    if (k%m == 0){
        ans = funcC(u,v);
        ans = funcC(u,v);

        // printf("t = %2.1d\tv-norm = %2.5f\n", i, ans);
    }
}
printf("%d",ans);
return 0;
}
