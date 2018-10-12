CPPflags = g++ -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

twoBody: twoBody.o
	${CPPflags} -o twoBody.x twoBody.cpp ${LIB}

main: main.o planet.o solver.o
	${CPPflags} -o main.x main.cpp planet.cpp solver.cpp ${LIB}

test: test.o planet.o solver.o
	${CPPflags} -o test.x test.cpp planet.cpp solver.cpp ${LIB}
	./test.x



clean:
	rm -f *.o *.txt *.x
