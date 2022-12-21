all: nielson-test

# https://github.com/salvipeter/libgeom/
LIBGEOM=../libgeom

INCLUDES=-I$(LIBGEOM)
LIBS=-L$(LIBGEOM)/release -lgeom

CXXFLAGS=-std=c++17 -Wall -pedantic -O3 $(INCLUDES)

nielson-test: nielson.o nielson-test.o
	$(CXX) -o $@ $^ $(LIBS)
