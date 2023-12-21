#include <stdio.h>				// needed for printing
#include <math.h>				// needed for tanh, used in init function
#include "params.h"				// model & simulation parameters
#include "mpi.h"

void init(double u[N][(N/4)+2], double v[N][(N/4)+2]){
	int rank, size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Status status;
	if (size != 4){	// Hardcoding a four-process decomposition
    	MPI_Abort( MPI_COMM_WORLD, 1 );
	}

	double uhi, ulo, vhi, vlo;
	uhi = 0.5; ulo = -0.5; vhi = 0.1; vlo = -0.1;

	int j_first, j_last;

	j_first = 1;
	j_last = N/4;
	if(rank==0){
		j_first++;
	}
	if(rank==size-1){
		j_last--;
	}
	
	for (int i=0; i < N; i++){
		for (int j=j_first; j <= j_last; j++){
			u[i][j] = ulo + (uhi-ulo)*0.5*(1.0 + tanh((i-N/2)/16.0));
			v[i][j] = vlo + (vhi-vlo)*0.5*(1.0 + tanh((j-N/2)/16.0));
		}
	}
	for(int i =0; i<1000; i++){
		printf("DJFIOPDJFIOPDSJFIOS");
	}

}

void dxdt(double du[N][(N/4)+2], double dv[N][(N/4)+2], double u[N][(N/4)+2], double v[N][(N/4)+2]){
	double lapu, lapv;
	int up, down, left, right;
	int rank, size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Status status;


	double du_edge[N+2], dv_edge[N+2];
	int j_first, j_last;
	j_first = 1;
	j_last = N/4;
	if(rank==0){
		j_first++;
	}
	if(rank==size-1){
		j_last--;
	}

	//Send right
	if(rank<size-1){
		for(int i=0; i<N; i++){
			du_edge[i] = du[i][j_last];
			dv_edge[i] = dv[i][j_last];
		}
		printf("Sending column %d from rank %d to rank %d\n", j_last, rank, rank+1);
		MPI_Send(du_edge, N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
		MPI_Send(dv_edge, N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
	}

	//Rec from left
	if(rank>0){
		printf("Recieving column %d from rank %d to rank %d\n", j_first-1, rank-1, rank);
		MPI_Recv(du_edge, N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD,&status);
		MPI_Recv(dv_edge, N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD,&status);
		for(int i=0; i<N; i++){
			du[i][j_first-1] = du_edge[i];
			dv[i][j_first-1] = dv_edge[i];
		}
	}
	//Send left
	if (rank>0){
		for(int i=0; i<N; i++){
			du_edge[i] = du[i][j_first];
			dv_edge[i] = dv[i][j_first];
		}
		printf("Sending column %d from rank %d to rank %d\n", j_first, rank, rank-1);
		MPI_Send(du_edge, N, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
		MPI_Send(dv_edge, N, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
	}

	//Rec from right
	if(rank<size-1){
		printf("Recieving column %d from rank %d to rank %d\n", j_last+1, rank+1, rank);
		MPI_Recv( du_edge, N, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status );
		MPI_Recv( dv_edge, N, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status );
		for(int i=0; i<N; i++){
			du[i][j_last+1] = du_edge[i];
			dv[i][j_last+1] = dv_edge[i];
		}
	}

	for (int i = 0; i < N; i++){
		for (int j = j_first; j <=j_last; j++){
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

void step(double du[N][(N/4)+2], double dv[N][(N/4)+2], double u[N][(N/4)+2], double v[N][(N/4)+2]){
	int rank, size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Status status;

	double u_edge[N+2], v_edge[N+2];
	int j_first, j_last;
	j_first = 1;
	j_last = N/4;
	if(rank==0){
		j_first++;
	}
	if(rank==size-1){
		j_last--;
	}

	//Send right
	if(rank<size-1){
		for(int i=0; i<N; i++){
			u_edge[i] = u[i][j_last];
			v_edge[i] = v[i][j_last];
		}
		printf("Sending column %d from rank %d to rank %d\n", j_last, rank, rank+1);
		MPI_Send(u_edge, N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
		MPI_Send(v_edge, N, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
	}

	//Rec from left
	if(rank>0){
		printf("Recieving column %d from rank %d to rank %d\n", j_first-1, rank-1, rank);
		MPI_Recv(u_edge, N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD,&status);
		MPI_Recv(v_edge, N, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD,&status);
		for(int i=0; i<N; i++){
			u[i][j_first-1] = u_edge[i];
			v[i][j_first-1] = v_edge[i];
		}
	}
	//Send left
	if (rank>0){
		for(int i=0; i<N; i++){
			u_edge[i] = u[i][j_first];
			v_edge[i] = v[i][j_first];
		}
		printf("Sending column %d from rank %d to rank %d\n", j_first, rank, rank-1);
		MPI_Send(u_edge, N, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
		MPI_Send(v_edge, N, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD);
	}

	//Rec from right
	if(rank<size-1){
		printf("Recieving column %d from rank %d to rank %d\n", j_last+1, rank+1, rank);
		MPI_Recv( u_edge, N, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status );
		MPI_Recv( v_edge, N, MPI_DOUBLE, rank + 1, 1, MPI_COMM_WORLD, &status );
		for(int i=0; i<N; i++){
			u[i][j_last+1] = u_edge[i];
			v[i][j_last+1] = v_edge[i];
		}
	}
	


	for (int i = 0; i < N; i++){
		for (int j = j_first; j <=j_last; j++){
			u[i][j] += dt*du[i][j];
			v[i][j] += dt*dv[i][j];
		}
	}
}

double norm(double x[N][(N/4)+2]){
	int rank, size;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );
	MPI_Status status;

	double nrmx = 0.0;
	int j_first, j_last;
	j_first = 1;
	j_last = N/4;
	if(rank==0){
		j_first++;
	}
	if(rank==size-1){
		j_last--;
	}

	for (int i = 0; i < N; i++){
		for (int j = j_first; j <=j_last; j++){
			nrmx += x[i][j]*x[i][j];
		}
	}
	return nrmx;
}

int main(int argc, char** argv){
	
	double t = 0.0, nrmu, nrmv, gnrmu, gnrmv;
	double u[N][(N/4)+2], v[N][(N/4)+2], du[N][(N/4)+2], dv[N][(N/4)+2];
	
	FILE *fptr = fopen("nrms.txt", "w");
	fprintf(fptr, "#t\t\tnrmu\t\tnrmv\n");
	MPI_Init( &argc, &argv );
	
	
	
	// initialize the state
	init(u, v);
	// for(int i =0; i<1000; i++){
	// 	printf("DJFIOPDJFIOPDSJFIOS");
	// }
	// time-loop
	for (int k=0; k < M; k++){
		// track the time
		t = dt*k;
		// evaluate the PDE
		dxdt(du, dv, u, v);
		// update the state variables u,v
		step(du, dv, u, v);
		if (k%m == 0){
			// calculate the norms
			nrmu = norm(u);
			MPI_Allreduce(&nrmu, &gnrmu, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			nrmv = norm(v);
			MPI_Allreduce(&nrmv, &gnrmv, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
			
			printf("t = %2.1f\tu-norm = %2.5f\tv-norm = %2.5f\n", t, gnrmu, gnrmv);
			fprintf(fptr, "%f\t%f\t%f\n", t, gnrmu, gnrmv);
		}
	}
	
	fclose(fptr);
	return 0;
}