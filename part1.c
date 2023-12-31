#include <stdio.h>				// needed for printing
#include <math.h>				// needed for tanh, used in init function
#include "params.h"				// model & simulation parameters
#include <omp.h> //openmp header file

void init(double u[N][N], double v[N][N]){
	double uhi, ulo, vhi, vlo;
	uhi = 0.5; ulo = -0.5; vhi = 0.1; vlo = -0.1; // these are shared vars
	//collapse wouldnt make difference if using 128 threads
	#pragma omp parrallel //128 is grain size shared(u, v) default(shared) schedule(static, 128)
	{
		#pragma omp for schedule(dynamic) 
			for (int i=0; i < N; i++){ //vars declared in loop are private
				for (int j=0; j < N; j++){
					u[i][j] = ulo + (uhi-ulo)*0.5*(1.0 + tanh((i-N/2)/16.0)); 
					v[i][j] = vlo + (vhi-vlo)*0.5*(1.0 + tanh((j-N/2)/16.0));
				}
			}
	}
}

void dxdt(double du[N][N], double dv[N][N], double u[N][N], double v[N][N]){ // u,v are not being changed (no need for reduction)
	
	// #pragma omp parrallel for schedule(static, 128)
	double lapu, lapv;
	int up, down, left, right;
	#pragma omp parrallel default(none) private(lapu, lapv, up, down, left, right) shared(du, dv, u, v, N, a, b, c, d)
	{ //Each value in grid is diff, hence computing du, dv for each i,v would be different, hence dynamic (i.e. work stealing is better)
		#pragma omp for schedule(dynamic, 64)  //64 + simple if would be faster than overhead of spawning a task (i.e. faster for same thread to do it)
			for (int i = 0; i < N; i++){
				for (int j = 0; j < N; j++){
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
}

void step(double du[N][N], double dv[N][N], double u[N][N], double v[N][N]){
	#pragma omp parrallel
	{
		#pragma omp for  //64 seems to work best, why? schedule(static, 64)
			for (int i = 0; i < N; i++){
				for (int j = 0; j < N; j++){
					u[i][j] += dt*du[i][j];
					v[i][j] += dt*dv[i][j];
				}
			}
	}

	// for (int i = 0; i < N; i++){
	// 	#pragma omp parallel  for
	// 	for (int j = 0; j < N; j++){
	// 		u[i][j] += dt*du[i][j];
	// 		v[i][j] += dt*dv[i][j];
	// 	}
	// }
}

double norm(double x[N][N]){
	double nrmx = 0.0;
	// #pragma omp parrallel for collapse(2)

	// #pragma omp parallel 
	// {
	// 	int partialsum = 0; #partial sum is slow because you would need 128 different partial sums
	// 	#pragma omp for
	// 	for (int i = 0; i < N; i++){
	// 		for (int j = 0; j < N; j++){
	// 			partialsum += x[i][j]*x[i][j];
	// 		}
	// 	}
	// 	#pragma omp atomic
	// 	nrmx += partialsum;
	// }
	#pragma omp parrallel for reduction(+:nrmx) schedule(static)
	for (int i = 0; i < N; i++){
		for (int j = 0; j < N; j++){
			nrmx += x[i][j]*x[i][j];
		}
		// printf("This program uses %d threads.\n", omp_get_num_threads());
	}
	// #pragma omp parrallel for reduction(+:nrmx) schedule(static)
	// for (int i = 0; i < N; i++){
	// 	for (int j = 0; j < N; j++){
	// 		nrmx += x[i][j]*x[i][j];
	// 	}
	// 	// printf("This program uses %d threads.\n", omp_get_num_threads());
	// }

	return nrmx;
}

int omp_thread_count() {
    int n = 0;
    #pragma omp parallel reduction(+:n)
    n += 1;
    return n;
}

int main(int argc, char** argv){
	
	double t = 0.0, nrmu, nrmv;
	double u[N][N], v[N][N], du[N][N], dv[N][N];
	
	FILE *fptr = fopen("nrms.txt", "w");
	fprintf(fptr, "#t\t\tnrmu\t\tnrmv\n");
	
	// initialize the state
	init(u, v);
	
	// time-loop
	for (int k=0; k < M; k++){
		// printf("This program uses %d threads.\n", omp_thread_count());
		// track the time
		t = dt*k;
		// evaluate the PDE
		dxdt(du, dv, u, v);
		// update the state variables u,v
		step(du, dv, u, v);
		if (k%m == 0){
			// calculate the norms
			nrmu = norm(u);
			nrmv = norm(v);
			printf("t = %2.1f\tu-norm = %2.5f\tv-norm = %2.5f\n", t, nrmu, nrmv);
			fprintf(fptr, "%f\t%f\t%f\n", t, nrmu, nrmv);
		}
	}
	
	fclose(fptr);
	return 0;
}