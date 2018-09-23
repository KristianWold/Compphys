CCflags = g++ -O3

LIB = -larmadillo -DARMA_DONT_USE_WRAPPER -llapack -lblas

PROG = testcode

${PROG} :	main.o func.o
			${CPPflags} main.o func.o ${LIB} -o ${PROG}

main.o :	main.cpp
			${CPPflags} -c main.cpp

func.o :	func.cpp
			${CPPflags} -c func.cpp
