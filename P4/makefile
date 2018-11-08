CC = gcc
CXX = mpic++
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

main.x: main.o spins.o montecarlo.o
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

main.o: main.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

spins.o: classes/spins.cpp
	g++ ${CPPflags} -c $^ ${INCL} ${LIB}

montecarlo.o: classes/montecarlo.cpp
	g++ ${CPPflags} -c $^ ${INCL} ${LIB}

test.x: test.o
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

test.o: test.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

clean:
	rm -f *.o *.txt *.x
