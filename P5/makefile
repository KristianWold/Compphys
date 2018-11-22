CC = gcc
CXX = g++
CPPflags = -g -Wall -std=c++11 -O3

UNAME := $(shell uname)

ifeq ($(UNAME), Darwin)
INCL := -I /usr/local/Cellar/armadillo/9.100.5/include
LIB := -DARMA_DONT_USE_WRAPPER -framework Accelerate
endif

ifeq ($(UNAME), Linux)
INCL :=
LIB := -DARMA_DONT_USE_WRAPPER -lblas -llapack
endif

qmc.x: qmc.o
	${CXX} ${CPPflags} -o $@ $^

qmc.o: qmc.cpp
	${CXX} ${CPPflags} -c $^

clean:
	rm -f *.o *.txt *.x
