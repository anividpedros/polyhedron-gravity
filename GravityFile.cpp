// GIVEN A DENSITY AND SHAPE MODEL, CREATES GRAVITY FILE//
#include <iostream>
#include <cassert>
#include <cmath>
#include "Polygon.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

void Polygon::CreateGravFile(std::string filename, double density)
{

  double Gs = 6.67384E-20;

  mGs = density * Gs;
  
  // Get faces and vertices
  ReadPolygon(filename);

  mNOV = vertices.size();
  mNOF = faces.size();

  // Verify Polyhedron property
  assert(mNOF == 2*mNOV - 4);

  mNOE = 3*(mNOV - 2);

  // Vertex coordinates
  mX = new double [mNOV];
  mY = new double [mNOV];
  mZ = new double [mNOV];

  for (int i = 0; i < mNOV; i++)
  {
    mX[i] = vertices[i].x;
    mY[i] = vertices[i].y;
    mZ[i] = vertices[i].z;
  } 

  mListTri = new int* [mNOF];
  mListN = new double* [mNOF];
  mF = new double* [mNOF];

  // To determine edges (including repeated ones)
  int** edgeIndexAll;
  // Declare New Size w nF 
  edgeIndexAll = new int*[3*mNOF];

  for (int i = 0; i < 3*mNOF; i++) {
      edgeIndexAll[i] = new int[2];
      }

  // vertex identifier
  int a, b, c;
  double normN;

  for (int i = 0; i < mNOF; i++)
  {
    mListTri[i] = new int [3];
    mListN[i]   = new double [3];
    mF[i] = new double [9];

    a = faces[i].i;
    b = faces[i].j;
    c = faces[i].k;

    mListTri[i][0] = a;
    mListTri[i][1] = b;
    mListTri[i][2] = c;


    // Vertex2 - Vertex1
    double d1[3] = {vertices[b-1].x-vertices[a-1].x,vertices[b-1].y-vertices[a-1].y,vertices[b-1].z-vertices[a-1].z};
    // Vertex3 - Vertex1
    double d2[3] = {vertices[c-1].x-vertices[a-1].x,vertices[c-1].y-vertices[a-1].y,vertices[c-1].z-vertices[a-1].z};
    
    // Face Normal (unitary)
    mListN[i][0] = d1[1]*d2[2] - d1[2]*d2[1];
    mListN[i][1] = d1[2]*d2[0] - d1[0]*d2[2];
    mListN[i][2] = d1[0]*d2[1] - d1[1]*d2[0];

    normN = sqrt(mListN[i][0]*mListN[i][0] + mListN[i][1]*mListN[i][1] + mListN[i][2]*mListN[i][2]);
    mListN[i][0] = mListN[i][0]/normN; mListN[i][1] = mListN[i][1]/normN; mListN[i][2] = mListN[i][2]/normN;  
    
    // Face Dyads
    mF[i][0] = mListN[i][0]*mListN[i][0];
    mF[i][1] = mListN[i][0]*mListN[i][1];
    mF[i][2] = mListN[i][0]*mListN[i][2];
    mF[i][3] = mListN[i][1]*mListN[i][0];
    mF[i][4] = mListN[i][1]*mListN[i][1];
    mF[i][5] = mListN[i][1]*mListN[i][2];
    mF[i][6] = mListN[i][2]*mListN[i][0];
    mF[i][7] = mListN[i][2]*mListN[i][1];
    mF[i][8] = mListN[i][2]*mListN[i][2];


    // Edges
    int edgeAll[3][2] = {{a,b},{b,c},{c,a}};
    // Three edges per facet, and two vertices per edge
    // if (a>b){
    //     edgeAll[0][0] = b;
    //     edgeAll[0][1] = a;
        
    //     }
    // if (b>c){
    //     edgeAll[1][0] = c;
    //     edgeAll[1][1] = b;
        
    //     }
    // if (c<a){
    //     edgeAll[2][0] = a;
    //     edgeAll[2][1] = c;
        
    //     } 
    // Connectivity matrix
    for (int jj = 0; jj < 3; jj++) {
      for (int kk = 0; kk < 2; kk++) {
        edgeIndexAll[(3*i)+jj][kk] = edgeAll[jj][kk];
      } 
    }

  }

  mListE = new int* [mNOE];
  mE = new double* [mNOE];

  // Delete Repeated Edges
  bool repeated = 0; //set to false by default
  int numEdges = 0; // total number of non-repeated edges
  std::vector<int> index; // indices for non-repeated edges
  for (int nn = 0; nn < 3*mNOF; nn++) {
      repeated = 0;
      for (int mm = 0; mm < nn; mm++) {
          if ((edgeIndexAll[mm][0] == edgeIndexAll[nn][0])&&(edgeIndexAll[mm][1] == edgeIndexAll[nn][1])&&(nn!=mm)){
              repeated = 1;
          }
          else if((edgeIndexAll[mm][0] == edgeIndexAll[nn][1])&&(edgeIndexAll[mm][1] == edgeIndexAll[nn][0])&&(nn!=mm)){
              repeated = 1;
          }
          
      }
      if (!repeated){
          index.push_back(nn);
          mListE[numEdges] = new int [2];
          mListE[numEdges][0] = edgeIndexAll[nn][0];
          mListE[numEdges][1] = edgeIndexAll[nn][1];

          numEdges++;
      }

    }

  assert (numEdges == mNOE);

// Edge direction
// Set positive or negative orientation for the edge facets
  int edgeDirection, indexVertex1, indexVertex2; // either -1 or 1
  bool edgeAssig; // when edge exists
  int compt = 0;
  int** edgeFacet;

  edgeFacet = new int*[mNOE]; // Face1, Face2 per edge

  for (int i = 0; i < mNOE; i++) {
      edgeFacet[i] = new int[2];

  }

  Vect nhatA(3), nhatB(3), nhatA12(3), nhatB21(3), edgeLineA(3), edgeLineB(3);

  for (int jj = 0; jj < mNOE; jj++){
      // Here get the index of vertices per edge
      indexVertex1 = mListE[jj][0];
      indexVertex2 = mListE[jj][1];
      compt = 0;
      for (int ii = 0; ii < mNOF; ii++) {
          edgeAssig = 0; //set to false initially
          if ((faces[ii].i == indexVertex1) && (faces[ii].j == indexVertex2))
          {
              edgeAssig = 1;
              edgeDirection = 1;
              compt++;
          }else if ((faces[ii].j== indexVertex1) && (faces[ii].i == indexVertex2))
          { 
              //inverted -> -1
              edgeAssig = 1;
              edgeDirection = -1;
              compt++;
          }
          else if ((faces[ii].i== indexVertex1) && (faces[ii].k == indexVertex2))
          {
              edgeAssig = 1;
              edgeDirection = 1;
              compt++;
          }else if ((faces[ii].k== indexVertex1) && (faces[ii].i == indexVertex2))
          {
              //inverted -> -1
              edgeAssig = 1;
              edgeDirection = -1;
              compt++;
          }else if ((faces[ii].j== indexVertex1) && (faces[ii].k == indexVertex2))
          {
              edgeAssig = 1;
              edgeDirection = 1;
              compt++;
          }else if ((faces[ii].k== indexVertex1) && (faces[ii].j == indexVertex2))
          {
              //inverted -> -1
              edgeAssig = 1;
              edgeDirection = -1;
              compt++;
          }
          
          //If edge assigned, compute direction
          if (edgeAssig){
              if (compt == 1){
                  edgeFacet[jj][0] = ii;

                  edgeLineA[0] = vertices[indexVertex2-1].x - vertices[indexVertex1-1].x;
                  edgeLineA[1] = vertices[indexVertex2-1].y - vertices[indexVertex1-1].y;
                  edgeLineA[2] = vertices[indexVertex2-1].z - vertices[indexVertex1-1].z;

                  // if (edgeDirection == -1){
                  //     edgeLineA[0] = - 1 * edgeLineA[0];
                  //     edgeLineA[1] = - 1 * edgeLineA[1];
                  //     edgeLineA[2] = - 1 * edgeLineA[2];
                  // }

              } 
              else if (compt == 2){
                  edgeFacet[jj][1] = ii;

                  edgeLineB[0] = vertices[indexVertex2-1].x - vertices[indexVertex1-1].x;
                  edgeLineB[1] = vertices[indexVertex2-1].y - vertices[indexVertex1-1].y;
                  edgeLineB[2] = vertices[indexVertex2-1].z - vertices[indexVertex1-1].z;

                  // if (edgeDirection == -1){
                  //     edgeLineB[0] = - 1 * edgeLineB[0];
                  //     edgeLineB[1] = - 1 * edgeLineB[1];
                  //     edgeLineB[2] = - 1 * edgeLineB[2];
                  // }

              }
          }
      }
  

      // Edge Dyads
      mE[jj] = new double [9];
      for (int i = 0; i < 3; i++)
      {
        nhatA[i] = mListN[edgeFacet[jj][0]][i];
        nhatB[i] = mListN[edgeFacet[jj][1]][i];
      }

      nhatA12 = cross(edgeLineA,nhatA);
      nhatA12 = nhatA12/norm(nhatA12);
      nhatB21 = cross(nhatB,edgeLineB);
      nhatB21 = nhatB21/norm(nhatB21);
      compt = 0;
      for (int nn = 0; nn < 3; nn++) {
          for (int mm = 0; mm < 3; mm++) {
                mE[jj][compt] = nhatA[nn]*nhatA12[mm] + nhatB[nn]*nhatB21[mm];
                compt++;  
            } 
      }

  }



// Memory handling
  for (int i = 0; i < 3*mNOF; i++) {
      delete edgeIndexAll[i];
  }
  delete edgeIndexAll;   
    
  for (int i = 0; i < mNOF; i++) {
      delete edgeFacet[i];
  }
  delete edgeFacet;
  
  // Create Gravity File
  std::ofstream GravityFile("GravityFile.txt");
  GravityFile << mGs << std::endl;
  GravityFile << mNOV << std::endl;
  GravityFile << mNOF << std::endl;
  GravityFile << mNOE << std::endl;

  for(int i = 0; i < mNOV; i++)
    {
      GravityFile << mX[i] << "\t"<< mY[i] << "\t"<< mZ[i]<< std::endl;
 
    } 

  for(int i = 0; i < mNOF; i++)
    {
      GravityFile<< mListTri[i][0] <<"\t"<< mListTri[i][1] <<"\t"<< mListTri[i][2]<< std::endl;
      GravityFile<< mListN[i][0] << "\t"<<mListN[i][1] <<"\t"<< mListN[i][2]<< std::endl;  
      GravityFile<< mF[i][0] <<"\t"<< mF[i][1] <<"\t"<< mF[i][2] <<"\t"<< mF[i][3] << "\t"<<mF[i][4] << "\t"<<mF[i][5] << "\t"<<mF[i][6] << "\t"<<mF[i][7] << "\t"<<mF[i][8]<< std::endl;  
       
    }

  for(int i = 0; i < mNOE; i++)
    {
      GravityFile <<  mListE[i][0] << "\t"<< mListE[i][1]<< std::endl;
      GravityFile <<  mE[i][0] <<"\t"<<  mE[i][1] << "\t"<< mE[i][2] <<"\t"<<  mE[i][3] << "\t"<< mE[i][4] <<"\t"<<  mE[i][5] << "\t"<<mE[i][6] << "\t"<< mE[i][7] <<"\t"<<  mE[i][8]<< std::endl;
      

    }

  GravityFile.close();

}


// METHODS 
// Read Shape Model
void Polygon::ReadPolygon(std::string filename)
        {
            faces.clear();
            vertices.clear();                      

            Point Vertex;
            Triangle Face;

            std::string line;
            char c;
            int i, j, k;
            double x, y, z;
            std::string si, sj, sk;

            std::ifstream in(filename);
            while ( getline( in, line ) )                           
            {
                if ( line.find_first_of( "vVfF" ) == std::string::npos ) continue;     // skip pointless lines

                std::istringstream ss( line );                            // put line into a stream for input
                ss >> c;                                             // get first character
                switch ( c )
                {
                    case 'v':                                         // vertices
                    case 'V':                                         // (either case)
                        ss >> x >> y >> z;                             // read from internal stream
                        Vertex.x  = x;
                        Vertex.y  = y;
                        Vertex.z  = z;                     // read from internal stream
                        vertices.push_back(Vertex);             // add to vector
                        break;

                    case 'f':                                         // faces
                    case 'F':                                         // (either case)
                        ss >> si >> sj >> sk;                          // FIRST, read into individual strings
                        i = stoi(si);  j = stoi(sj);  k = stoi(sk);  // Get the FIRST integer from each string
                        Face.i = i;
                        Face.j = j;
                        Face.k = k;
                        faces.push_back(Face);
                        break;
                }
            }
            in.close();
        }      


int main()
{  

    // std::string filename = "GravityMapping/Bennu50kfixed.obj";
    std::string filename = "GravityMapping/Bennu50kfixed.obj";
    double density = 1.26e12; //kg/km3

    // Polygon polygon(filename, density);
    Polygon p;
    p.CreateGravFile(filename,density);

    return 0;

}
