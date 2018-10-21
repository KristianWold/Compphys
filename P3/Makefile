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

OBJS = planet.o solver.o

TARGETS = main.x allPlanets.x test.x

all: $(TARGETS)

main.x: main.o $(OBJS)
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

main.o: main.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}


allPlanets.x: allPlanets.o $(OBJS)
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

allPlanets.o: allPlanets.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}


escapeVel.x: escapeVel.o $(OBJS)
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

escapeVel.o: escapeVel.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}


earthAndJupiter.x: earthAndJupiter.o $(OBJS)
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

earthAndJupiter.o: earthAndJupiter.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}


relativistic.x: relativistic.o $(OBJS)
	${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}

relativistic.o: relativistic.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}


solver.o: Classes/solver.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

planet.o: Classes/planet.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

test.x: test.o $(OBJS)
		${CXX} ${CPPflags} -o $@ $^ ${INCL} ${LIB}
		./test.x

test.o: test.cpp
	${CXX} ${CPPflags} -c $^ ${INCL} ${LIB}

clean:
	rm -f *.o *.txt *.x
