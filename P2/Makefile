CPPflags = g++ -O3

LIB = -DARMA_DONT_USE_WRAPPER -lblas -llapack

quantumDot: quantumDot.o jacobi.o
	${CPPflags} -o quantumDot.x quantumDot.cpp jacobi.cpp ${LIB}

eigenvectors: eigenvectors.o jacobi.o
	${CPPflags} -o eigenvectors.x eigenvectors.cpp jacobi.cpp ${LIB}

test: test.o jacobi.o
	${CPPflags} -o test.x test.cpp jacobi.cpp ${LIB}
	./test.x

benchmark: benchmark.o jacobi.o
	${CPPflags} -o benchmark.x benchmark.cpp jacobi.cpp ${LIB}
	./benchmark.x
	python plot.py benchmark

clean:
	rm -f *.o *.txt *.x
