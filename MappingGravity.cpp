#include <iostream>
#include <cassert>
#include <cmath>
#include "Polygon.hpp"
#include "Matrix.hpp"
#include "Vect.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <iomanip> 

#include <stdio.h>
// For MPI
#include <mpi.h>

// RETURNS GRAVITY FOR A FILE OF POINTS//
// USES MPI TO EVALUATE SEVERAL POINTS SIMULTANEOUSLY//
int main(int argc, char** argv)
{  
	int p_id, n_procs, flag = 1;

	MPI_Init(NULL, NULL); // initialize MPI environment

	// READ GRAVITY FILE
	std::ifstream GravityFile("GravityFile.txt");
    Polygon polygon(GravityFile);
	Vect grav_output(12); // Acceleration + gravity matrix + position


	// READ FILE WITH MAPPING POINTS
	std::ifstream MapGravity("GravityMapping/Points.txt");
	assert(MapGravity.is_open());
	int nPoints;
	MapGravity >> nPoints;

	// MPI VARIABLES
	int world_size; // number of processes
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);
	int world_rank; // the rank of the process
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);   

	char processor_name[MPI_MAX_PROCESSOR_NAME]; // gets the name of the processor
	int name_len;
	MPI_Get_processor_name(processor_name, &name_len);
	int elperproc = nPoints/world_size;

	// Read S/C POSITIONS FROM THE .TXT FILE //
	double* Xsc = NULL;
	double* MapGrav = NULL;
	if (world_rank == 0) {
		Xsc = new double[3*nPoints];
		
		for (int l = 0; l < 3*nPoints; l++)
		{
			MapGravity>>Xsc[l];
		}

		MapGrav = new double[12*nPoints];
	}
	
	double* MapGravS;
	MapGravS = new double[12*elperproc];
	
	double* XscS;
	XscS = new double[3*elperproc];

	// SPLIT UP POINTS WITH PROCESSORS // 
	MPI_Scatter(Xsc, 3*elperproc, MPI_DOUBLE, XscS, 3*elperproc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	
	Vect pos(3);
	for (int ii = 0; ii < elperproc; ii++)
	{
		pos[0] = XscS[ii*3];
		pos[1] = XscS[ii*3+1];
		pos[2] = XscS[ii*3+2];
		grav_output = PolyGrav(pos,polygon);

		// Acceleration
		MapGravS[ii*12]   = grav_output[0];
		MapGravS[ii*12+1] = grav_output[1];
		MapGravS[ii*12+2] = grav_output[2];

		// Gradient
		MapGravS[ii*12+3]  = grav_output[3];
		MapGravS[ii*12+4]  = grav_output[4];
		MapGravS[ii*12+5]  = grav_output[5];
		MapGravS[ii*12+6]  = grav_output[6];
		MapGravS[ii*12+7]  = grav_output[7];
		MapGravS[ii*12+8]  = grav_output[8];
		MapGravS[ii*12+9]  = grav_output[9];
		MapGravS[ii*12+10] = grav_output[10];
		MapGravS[ii*12+11] = grav_output[11];

		// //Position
		// MapGravS[ii*15+12] = pos[0];
		// MapGravS[ii*15+13] = pos[1];
		// MapGravS[ii*15+14] = pos[2];
	}
	
	MPI_Gather(MapGravS, 12*elperproc, MPI_DOUBLE, MapGrav, 12*elperproc, MPI_DOUBLE, 0, MPI_COMM_WORLD);

	delete XscS;
	delete MapGravS;

	if (world_rank == 0)
	{	
		std::ofstream GravityFile ("GravityMapping/GravityMap.txt");
		assert(GravityFile.is_open());
		for (size_t i = 0; i < nPoints; i++)
		{
			GravityFile<<std::fixed<<std::setprecision(16)<<MapGrav[i*12]<<'\t'<<MapGrav[i*12+1]<<'\t'<<MapGrav[i*12+2]<<'\n';	
		}
		GravityFile.close();

		std::ofstream GravityGradientFile ("GravityMapping/GravityGradient.txt");
		assert(GravityGradientFile.is_open());
		for (size_t k = 0; k < nPoints; k++)
		{
			for (size_t i = 0; i < 3; i++)
			{
				GravityGradientFile<<MapGrav[3+3*i+12*k]<<'\t'<<MapGrav[4+3*i+12*k]<<'\t'<<MapGrav[5+3*i+12*k]<<'\n';	
			}
		
		}
		GravityGradientFile.close();

		// std::ofstream GravityPositionFile ("GravityMapping/GravityPosition.txt");
		// assert(GravityPositionFile.is_open());
		// for (size_t i = 0; i < nPoints; i++)
		// {
		// 	GravityPositionFile<<Xsc[i*3]<<'\t'<<Xsc[i*3+1]<<'\t'<<Xsc[i*3+2]<<'\n';	
		// }
		// GravityPositionFile.close();

		delete Xsc;
		delete MapGrav;
	}

	MPI_Finalize();

	return 0;
}



