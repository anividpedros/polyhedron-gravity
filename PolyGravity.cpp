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

// RETURNS GRAVITY FOR A SINGLE POINT //
int main()
{  

    std::ifstream GravityFile("GravityFile.txt");
    Polygon polygon(GravityFile);
    
    Vect Xsc(3);
    Xsc [0] = 200;
    Xsc [1] = 0;
    Xsc [2] = 0;    

    Vect gravity(12);
    gravity = PolyGrav(Xsc, polygon);
    disp(gravity);
    
    return 0;
}
