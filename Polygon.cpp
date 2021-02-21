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

Polygon::Polygon(){}

Polygon::Polygon(std::ifstream& GravityFile)
{

  assert(GravityFile.is_open());

  double Gs;

  GravityFile >> Gs;

  mGs = Gs;

  int NOV, NOF, NOE;

  GravityFile >> NOV >> NOF >> NOE;

  mNOV = NOV;
  mNOF = NOF;
  mNOE = NOE;

  assert(mNOF == 2*mNOV - 4);
  assert(mNOE == 3*(mNOV - 2));

  mX = new double [mNOV];
  mY = new double [mNOV];
  mZ = new double [mNOV];

  for(int i = 0; i < mNOV; i++)
    {
      GravityFile >> mX[i] >> mY[i] >> mZ[i];   
    } 

  
  mListTri = new int* [mNOF];
  mListN = new double* [mNOF];
  mF = new double* [mNOF];

  for(int i = 0; i < mNOF; i++)
    {
      mListTri[i] = new int [3];
      mListN[i]   = new double [3];
      
      mF[i] = new double [9];

      GravityFile >> mListTri[i][0] >> mListTri[i][1] >> mListTri[i][2];
      GravityFile >> mListN[i][0] >> mListN[i][1] >> mListN[i][2];
      GravityFile >> mF[i][0] >> mF[i][1] >> mF[i][2] >> mF[i][3] >> mF[i][4] >> mF[i][5] >> mF[i][6] >> mF[i][7] >> mF[i][8];
    }

  mListE = new int* [mNOE];
  
  mE = new double* [mNOE];

  for(int i = 0; i < mNOE; i++)
    {
      mListE[i] = new int [2];

      mE[i] = new double [9];

      GravityFile >> mListE[i][0] >> mListE[i][1];
      GravityFile >> mE[i][0] >> mE[i][1] >> mE[i][2] >> mE[i][3] >> mE[i][4] >> mE[i][5] >> mE[i][6] >> mE[i][7] >> mE[i][8];
    }
}


Polygon::~Polygon()
{

  delete[] mX;
  delete[] mY;
  delete[] mZ;

  for(int i = 0; i < mNOF; i++)
    {
      delete[] mListTri[i];
      delete[] mListN[i];
      delete[] mF[i];
    }

  delete[] mListTri;
  delete[] mListN;
  delete[] mF;

  for(int i = 0; i < mNOE; i++)
    {
      delete[] mListE[i];
      delete[] mE[i];
    }

  delete[] mListE;
  delete[] mE;

}
   

// Get G*sigma value
double Polygon::GetGs() const
{
  return mGs;
}

// Get No. of Vertices
int Polygon::GetNOV() const
{
  return mNOV;
}

// Get No. of Facets
int Polygon::GetNOF() const
{
  return mNOF;
}

// Get No. of Edges
int Polygon::GetNOE() const
{
  return mNOE;
}


// Get Vertex X coordinates
Vect Polygon::GetX() const
{
  Vect x(mNOV);
  
  for(int i = 0; i < mNOV; i++)
    {
      x[i] = mX[i];
    }

  return x;
}

// Get Vertex Y coordinates
Vect Polygon::GetY() const
{
  Vect y(mNOV);
  
  for(int i = 0; i < mNOV; i++)
    {
      y[i] = mY[i];
    }

  return y;
}

// Get Vertex Z coordinates
Vect Polygon::GetZ() const
{
  Vect z(mNOV);
  
  for(int i = 0; i < mNOV; i++)
    {
      z[i] = mZ[i];
    }

  return z;
}


Vect Polygon::GetListTri() const
{

	Vect ListTri(mNOF*3);

	for(int i = 0; i < mNOF; i++)
	{
		for(int j = 0; j < 3; j++)
		{
		  ListTri[3*i+j] = (double) mListTri[i][j];
		}
	}

	return ListTri;
}


Vect Polygon::GetListN() const
{

	Vect ListN(mNOF*3);

	for(int i = 0; i < mNOF; i++)
	{
		for(int j = 0; j < 3; j++)
		{
			ListN[3*i+j] = mListN[i][j];
		}
	}

	return ListN;
}


Vect Polygon::GetF() const
{

	Vect F(mNOF*9);

	for(int i = 0; i < mNOF; i++)
	{
		for(int j = 0; j < 9; j++)
		{
			F[9*i+j] = mF[i][j];
		}
	}

	return F;
}


Vect Polygon::GetListE() const
{

	Vect ListE(mNOE*2);

	for(int i = 0; i < mNOE; i++)
	{
		for(int j = 0; j < 2; j++)
		{
			ListE[2*i+j] = mListE[i][j];
		}
	}

	return ListE;
}


Vect Polygon::GetE() const
{

	Vect E(mNOE*9);

	for(int i = 0; i < mNOE; i++)
	{
		for(int j = 0; j < 9; j++)
		{
			E[9*i+j] = mE[i][j];
		}
	}

	return E;
}


Vect PolyGrav(Vect& Xsc, Polygon& Body)
{

  Vect a_grav(3); // Gravity
  Matrix g_grav(3,3); //Gravity Gradient Matrix

  int v1, v2, v3;
  Vect r1(3), r2(3), r3(3);
  double R1, R2, R3;
  Vect  nf(3);
  Matrix F(3,3);
  int l;
  double wf;

  //Sum over Polyhedron Facets
  for(int i = 0; i < Body.mNOF; i++)
    {
      // Vertex no. 1
      v1 = Body.mListTri[i][0];
      r1[0] = Body.mX[v1-1];
      r1[1] = Body.mY[v1-1];
      r1[2] = Body.mZ[v1-1];

      // Vertex no. 2
      v2 = Body.mListTri[i][1];
      r2[0] = Body.mX[v2-1];
      r2[1] = Body.mY[v2-1];
      r2[2] = Body.mZ[v2-1];

      // Vertex no. 3
      v3 = Body.mListTri[i][2];
      r3[0] = Body.mX[v3-1];
      r3[1] = Body.mY[v3-1];
      r3[2] = Body.mZ[v3-1];

      // Normal Unit Vector
      nf[0] = Body.mListN[i][0];
      nf[1] = Body.mListN[i][1];
      nf[2] = Body.mListN[i][2];

      // Face Dyad
      l = 0;
      for(int j=0; j < 3; j++)
        {
          for(int k=0; k < 3; k++)
            {
              F(j,k) = Body.mF[i][l];
              l++;
            }
        }

      //Pos vector wrt S/C
      r1 = r1 - Xsc;
      r2 = r2 - Xsc;
      r3 = r3 - Xsc;

      R1 = norm(r1);
      R2 = norm(r2);
      R3 = norm(r3);

      wf = 2*atan2(dot(r1,cross(r2,r3)), R1*R2*R3 + R1*dot(r2,r3) + R2*dot(r3,r1) + R3*dot(r1,r2));

      a_grav = a_grav + F*r1*wf;
      g_grav = g_grav - F*wf;
    }

  // Sum over the polyhedron edges
  int e1, e2;
  double Re, Le;
  Matrix E(3,3);

  for(int i = 0; i < Body.mNOE; i++)
    {
     
      e1 = Body.mListE[i][0];
      e2 = Body.mListE[i][1];

      //Vertex 1 pos vector from spacecraft
      r1[0] = Body.mX[e1-1];
      r1[1] = Body.mY[e1-1];
      r1[2] = Body.mZ[e1-1];
      r1 = r1 - Xsc;

      R1 = norm(r1);

      //Vertex 2 pos vector from s/c
      r2[0] = Body.mX[e2-1];
      r2[1] = Body.mY[e2-1];
      r2[2] = Body.mZ[e2-1];
      r2 = r2 - Xsc;

      R2 = norm(r2);

      Re = norm(r2-r1);
      
      Le = log((R1 + R2 + Re)/(R1 + R2 - Re));

      // Edge dyad
      l = 0;
      for(int j = 0; j < 3; j++)
	      {
	        for(int k = 0; k < 3; k++)
	          {
	            E(j,k) = Body.mE[i][l];
	            l++;
	          }
	      }

      // Gravitational Acceleration
      a_grav = a_grav - E*r1*Le;
      g_grav = g_grav + E*Le;
  
    }
  g_grav = Body.mGs*g_grav;
  a_grav = Body.mGs*a_grav;

  Vect Output(12);
  for (int i = 0; i < 3; i++)
  {
    Output[i] = a_grav[i];

    for (int j = 0; j < 3; j++)
      {
          Output[3+3*i+j] = g_grav(i,j);

      }
  } 

  return Output;

}





