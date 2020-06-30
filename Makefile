all: bsconnect-test

INCLUDES=-I../libgeom
LIBS=-L../libgeom/debug -lgeom

CXXFLAGS=-Wall -std=c++17 -pedantic -g -DDEBUG $(INCLUDES)

bsconnect-test: bsconnect-test.o bsconnect.o
	$(CXX) -o $@ $^ $(LIBS)
