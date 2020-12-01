Matrix.o: Matrix.cpp Matrix.hpp 
	g++ -g -c Matrix.cpp
 
Vect.o: Vect.cpp Vect.hpp
	g++ -g -c Vect.cpp

Polygon.o: Vect.cpp Vect.hpp Matrix.cpp Matrix.hpp Polygon.cpp Polygon.hpp
	g++ -g -c Polygon.cpp

GravityFile.o: Vect.cpp Vect.hpp Matrix.cpp Matrix.hpp GravityFile.cpp Polygon.hpp
	g++ -g -c GravityFile.cpp

PolyGravity.o: Vect.cpp Vect.hpp Matrix.cpp Matrix.hpp PolyGravity.cpp Polygon.hpp
	g++ -g -c PolyGravity.cpp

MappingGravity.o: Vect.cpp Vect.hpp Matrix.cpp Matrix.hpp MappingGravity.cpp Polygon.hpp
	mpic++ -c MappingGravity.cpp

CreateGravFile: Vect.o Matrix.o Polygon.o GravityFile.o
	g++ -o CreateGravFile Vect.o Matrix.o Polygon.o GravityFile.o

PolyGrav: Vect.o Matrix.o Polygon.o PolyGravity.o
	g++ -g -o PolyGrav PolyGravity.o Polygon.o Vect.o Matrix.o 

MapGravityMPI: Vect.o Matrix.o MappingGravity.o Polygon.o
	mpic++ -o MapGravityMPI MappingGravity.o Vect.o Matrix.o Polygon.o

PolyGravCUDA: Vect.o Matrix.o Polygon.o
	nvcc PolyGravCUDA.cu -o PolyGravCUDA Polygon.o Vect.o Matrix.o

clean:
	rm *.o
