all: bsconnect-test

GEOM=../libgeom
EIGEN=/usr/include/eigen3

INCLUDES=-I$(GEOM) -I$(EIGEN)
LIBS=-L$(GEOM)/debug -lgeom

CXXFLAGS=-Wall -std=c++17 -pedantic -g -DDEBUG $(INCLUDES)

bsconnect-test: bsconnect-test.o bsconnect.o
	$(CXX) -o $@ $^ $(LIBS)
