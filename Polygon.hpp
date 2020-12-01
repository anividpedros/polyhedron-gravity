#ifndef POLYGONHEADERDEF
#define POLYGONHEADERDEF

#include <fstream>
#include "Matrix.hpp"
#include "Vect.hpp"
#include <string>
#include <vector>

class Polygon{
    private:
    struct Point{ double x, y, z; };       // Spatial coordinates
    struct Triangle{ int i, j, k; };       // Vertex numbers for a triangle

    double mGs;   // G*dens product

    int mNOV; //number of vertex
    int mNOF; // number of facets
    int mNOE; // mumber of edges

    // Vertex coordinates
    double* mX;
    double* mY;
    double* mZ;


    // Triangle Vertex List
    int** mListTri;

    // Face Normal Unit Vectors
    double** mListN;

    // Face Dyads
    double** mF;

    // Edge Vertex List
    int** mListE;

    // Edge Dyads
    double** mE;

    // Faces
    std::vector<Triangle> faces; // Indices for vertices per face

    // Vertices 
    std::vector<Point> vertices;

    public:
    
    Polygon();
    Polygon(std::ifstream& filename);
    ~Polygon();

    void CreateGravFile(std::string filename, double density);
    void ReadPolygon (std::string filename);

    double GetGs() const;

    int GetNOV() const;
    int GetNOF() const;
    int GetNOE() const;

    Vect GetX() const;
    Vect GetY() const;
    Vect GetZ() const;

    Vect GetListTri() const;
    Vect GetListN() const;
    Vect GetF() const;

    Vect GetListE() const;
    Vect GetE() const;

    friend Vect PolyGrav(Vect& Xsc, Polygon& Body);
    
};

#endif