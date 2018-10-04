CPPflags = g++ -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

twoBody: twoBody.o
	${CPPflags} -o twoBody.x twoBody.cpp ${LIB}

clean:
	rm -f *.o *.txt *.x
