CPPflags = g++ -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

twoBody: twoBody.o
	${CPPflags} -o twoBody.x twoBody.cpp ${LIB}

Planet: Planet.o
	${CPPflags} -o Planet.x Planet.cpp ${LIB}

clean:
	rm -f *.o *.txt *.x
