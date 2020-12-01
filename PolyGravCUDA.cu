#include <iostream>
#include <cmath>
#include <fstream>
#include "Polygon.hpp"
#include "Matrix.hpp"
#include "Vect.hpp"
#include "book.h"

const int blocksPerGrid = 32;
const int threadsPerBlock = 256;


__global__ void FaceTerm(int* NOF, double* ListTri, double* F, double* X, double* Y, double* Z, double* Xsc, double* acc)
{

__shared__ double ax[threadsPerBlock];
__shared__ double ay[threadsPerBlock];
__shared__ double az[threadsPerBlock];

int i = threadIdx.x + blockIdx.x*blockDim.x;
int cacheIndex = threadIdx.x;

int v1, v2, v3;

double r1x, r1y, r1z;
double r2x, r2y, r2z;
double r3x, r3y, r3z;

double R1, R2, R3;

double dot_num, dot_23, dot_31, dot_12;

double wf;

double tempx = 0.0;
double tempy = 0.0;
double tempz = 0.0;

while(i < *NOF)
{

	v1 = (int) ListTri[3*i];
	v2 = (int) ListTri[3*i + 1];
	v3 = (int) ListTri[3*i + 2];

	//Define First Vertex
	r1x = X[v1-1] - Xsc[0];
	r1y = Y[v1-1] - Xsc[1];
	r1z = Z[v1-1] - Xsc[2];
	R1  = sqrt(r1x*r1x + r1y*r1y + r1z*r1z);

	//Define Second Vertex
	r2x = X[v2-1] - Xsc[0];
	r2y = Y[v2-1] - Xsc[1];
	r2z = Z[v2-1] - Xsc[2];
	R2  = sqrt(r2x*r2x + r2y*r2y + r2z*r2z);

	//Define Third Vertex
	r3x = X[v3-1] - Xsc[0];
	r3y = Y[v3-1] - Xsc[1];
	r3z = Z[v3-1] - Xsc[2];
	R3  = sqrt(r3x*r3x + r3y*r3y + r3z*r3z);


	dot_num = r1x*(r2y*r3z - r2z*r3y) + r1y*(r2z*r3x - r2x*r3z) + r1z*(r2x*r3y - r2y*r3x);

	dot_23 = r2x*r3x + r2y*r3y + r2z*r3z;
	dot_31 = r3x*r1x + r3y*r1y + r3z*r1z;
	dot_12 = r1x*r2x + r1y*r2y + r1z*r2z;

	//Define w_f
	wf = 2*atan2(dot_num, R1*R2*R3 + R1*dot_23 + R2*dot_31 + R3*dot_12);

	//Store acceleration
	tempx += wf*(F[9*i + 0]*r1x + F[9*i + 1]*r1y + F[9*i + 2]*r1z);
	tempy += wf*(F[9*i + 3]*r1x + F[9*i + 4]*r1y + F[9*i + 5]*r1z);
	tempz += wf*(F[9*i + 6]*r1x + F[9*i + 7]*r1y + F[9*i + 8]*r1z);

	i += blockDim.x * gridDim.x;
}

ax[cacheIndex] = tempx;
ay[cacheIndex] = tempy;
az[cacheIndex] = tempz;

// synchronize threads in the block
__syncthreads();


int j = blockDim.x/2;

while(j != 0){

	if(cacheIndex < j)
	{
		ax[cacheIndex] += ax[cacheIndex + j];
		ay[cacheIndex] += ay[cacheIndex + j];
		az[cacheIndex] += az[cacheIndex + j];
	}
	__syncthreads();

	j /= 2;
}

if(cacheIndex == 0)
{
	acc[3*blockIdx.x + 0] = ax[0];
	acc[3*blockIdx.x + 1] = ay[0];
	acc[3*blockIdx.x + 2] = az[0];
}


}



__global__ void EdgeTerm(int* NOE, double* ListE, double* E, double* X, double* Y, double* Z, double* Xsc, double* acc)
{

__shared__ double ax[threadsPerBlock];
__shared__ double ay[threadsPerBlock];
__shared__ double az[threadsPerBlock];

	int i = threadIdx.x + blockIdx.x*blockDim.x;
	int cacheIdx = threadIdx.x;

	int e1, e2;

	double r1x, r1y, r1z, R1;
	double r2x, r2y, r2z, R2;

	double dex, dey, dez, de;

	double Le;

	double tempx = 0.0;
	double tempy = 0.0;
	double tempz = 0.0;

	while(i < *NOE)
	{

		e1 = (int) ListE[2*i];
		e2 = (int) ListE[2*i + 1];

		//Define First Vertex
		r1x = X[e1-1] - Xsc[0];
		r1y = Y[e1-1] - Xsc[1];
		r1z = Z[e1-1] - Xsc[2];
		R1  = sqrt(r1x*r1x + r1y*r1y + r1z*r1z);

		//Define Second Vertex
		r2x = X[e2-1] - Xsc[0];
		r2y = Y[e2-1] - Xsc[1];
		r2z = Z[e2-1] - Xsc[2];
		R2  = sqrt(r2x*r2x + r2y*r2y + r2z*r2z);

		dex = r2x - r1x;
		dey = r2y - r1y;
		dez = r2z - r1z;
		de = sqrt(dex*dex + dey*dey + dez*dez);

		Le  = log((R1 + R2 + de)/(R1 + R2 - de));

		tempx += Le*(E[9*i + 0]*r1x + E[9*i + 1]*r1y + E[9*i + 2]*r1z);
		tempy += Le*(E[9*i + 3]*r1x + E[9*i + 4]*r1y + E[9*i + 5]*r1z);
		tempz += Le*(E[9*i + 6]*r1x + E[9*i + 7]*r1y + E[9*i + 8]*r1z);

		i += blockDim.x * gridDim.x;

	}

	ax[cacheIdx] = tempx;
	ay[cacheIdx] = tempy;
	az[cacheIdx] = tempz;

	__syncthreads();


	int j = blockDim.x/2;

	while(j!=0)
	{
	if(cacheIdx < j)
	{
	 ax[cacheIdx] += ax[cacheIdx + j];
	 ay[cacheIdx] += ay[cacheIdx + j];
	 az[cacheIdx] += az[cacheIdx + j];
	}

	__syncthreads();
	j /= 2;

	}

	if(cacheIdx == 0){
	 acc[3*blockIdx.x + 0] = ax[0];
	 acc[3*blockIdx.x + 1] = ay[0];
	 acc[3*blockIdx.x + 2] = az[0];
	} 

}



int main()
{

std::ifstream GRAV ("GravityFile.txt");

Polygon CG67P(GRAV);

int NOV = CG67P.GetNOV();

// Vertex Coordinates
Vect X = CG67P.GetX();
Vect Y = CG67P.GetY();
Vect Z = CG67P.GetZ();

// Load Vertex Coordinates on GPU
double* d_X;
double* d_Y;
double* d_Z;

HANDLE_ERROR(cudaMalloc((void**)&d_X, NOV*sizeof(double)));
HANDLE_ERROR(cudaMalloc((void**)&d_Y, NOV*sizeof(double)));
HANDLE_ERROR(cudaMalloc((void**)&d_Z, NOV*sizeof(double)));

HANDLE_ERROR(cudaMemcpy(d_X, &X[0], NOV*sizeof(double), cudaMemcpyHostToDevice));
HANDLE_ERROR(cudaMemcpy(d_Y, &Y[0], NOV*sizeof(double), cudaMemcpyHostToDevice));
HANDLE_ERROR(cudaMemcpy(d_Z, &Z[0], NOV*sizeof(double), cudaMemcpyHostToDevice));


// Facet Term
int NOF = CG67P.GetNOF();


Vect ListTri = CG67P.GetListTri();
//Matrix ListN   = KW4A.GetListN();
Vect F       = CG67P.GetF();


//Copy Variables on GPU
int*    d_NOF;
double* d_ListTri;
//double* d_ListN;
double* d_F;


HANDLE_ERROR(cudaMalloc((void**) &d_NOF, sizeof(int)));
HANDLE_ERROR(cudaMalloc((void**) &d_ListTri, NOF*3*sizeof(double)));
//cudaMalloc((void**) &d_ListN, NOF*3*sizeof(double));
HANDLE_ERROR(cudaMalloc((void**) &d_F, NOF*9*sizeof(double)));

HANDLE_ERROR(cudaMemcpy(d_NOF, &NOF, sizeof(int), cudaMemcpyHostToDevice));
HANDLE_ERROR(cudaMemcpy(d_ListTri, &ListTri[0], NOF*3*sizeof(double), cudaMemcpyHostToDevice));
//cudaMemcpy(d_ListN, &ListN[0], NOF*3*sizeof(double), cudaMemcpyHostToDevice);
HANDLE_ERROR(cudaMemcpy(d_F, &F[0], NOF*9*sizeof(double), cudaMemcpyHostToDevice));


//Define Spacecraft position vector
Vect Xsc(3);

Xsc[0] = 1000.0;
Xsc[1] = 0.;
Xsc[2] = 0.;


double* d_Xsc;

HANDLE_ERROR(cudaMalloc((void**) &d_Xsc, 3*sizeof(double)));

HANDLE_ERROR(cudaMemcpy(d_Xsc, &Xsc[0], 3*sizeof(double), cudaMemcpyHostToDevice));

double* d_aF;

HANDLE_ERROR(cudaMalloc((void**) &d_aF, blocksPerGrid*3*sizeof(double)));

FaceTerm<<<blocksPerGrid,threadsPerBlock>>>(d_NOF, d_ListTri, d_F, d_X, d_Y, d_Z, d_Xsc, d_aF);

Vect accF(blocksPerGrid*3);

HANDLE_ERROR(cudaMemcpy(&accF[0], d_aF, blocksPerGrid*3*sizeof(double), cudaMemcpyDeviceToHost));



// Edge Term
int  NOE    = CG67P.GetNOE();
Vect ListE  = CG67P.GetListE();
Vect E      = CG67P.GetE();

int* d_NOE;
double* d_ListE;
double* d_E;

HANDLE_ERROR(cudaMalloc((void**) &d_NOE, sizeof(int)));
HANDLE_ERROR(cudaMalloc((void**) &d_ListE, NOE*2*sizeof(double)));
HANDLE_ERROR(cudaMalloc((void**) &d_E, NOE*9*sizeof(double)));

HANDLE_ERROR(cudaMemcpy(d_NOE, &NOE, sizeof(int), cudaMemcpyHostToDevice));
HANDLE_ERROR(cudaMemcpy(d_ListE, &ListE[0], NOE*2*sizeof(double), cudaMemcpyHostToDevice));
HANDLE_ERROR(cudaMemcpy(d_E, &E[0], NOE*9*sizeof(double), cudaMemcpyHostToDevice));


double* d_aE;

HANDLE_ERROR(cudaMalloc((void**) &d_aE, blocksPerGrid*3*sizeof(double)));

EdgeTerm<<<blocksPerGrid,threadsPerBlock>>>(d_NOE, d_ListE, d_E, d_X, d_Y, d_Z, d_Xsc, d_aE);

Vect accE(blocksPerGrid*3);

HANDLE_ERROR(cudaMemcpy(&accE[0], d_aE, blocksPerGrid*3*sizeof(double),cudaMemcpyDeviceToHost));

double ax = 0.0;
double ay = 0.0;
double az = 0.0;

for(int i = 0; i < blocksPerGrid; i++)
{
ax += accF[3*i]   - accE[3*i];
ay += accF[3*i+1] - accE[3*i+1];
az += accF[3*i+2] - accE[3*i+2];
}


double Gs = CG67P.GetGs();

//std::cout << "Gs: " << Gs << std::endl;

Vect acc(3);
acc[0] = Gs*ax;
acc[1] = Gs*ay;
acc[2] = Gs*az;


std::cout << "Final Acceleration:" << std::endl; 
disp(acc);



// Free CUDA MEMORY
//Spacecraft Coordinates
cudaFree(d_Xsc);

//Vertex Coordinates
cudaFree(d_X);
cudaFree(d_Y);
cudaFree(d_Z);

//Facet Term
cudaFree(d_NOF);
cudaFree(d_ListTri);
//cudaFree(d_ListN);
cudaFree(d_F);
cudaFree(d_aF);

//Edge Term
cudaFree(d_NOE);
cudaFree(d_ListE);
cudaFree(d_E);
cudaFree(d_aE);

return 0;

}
