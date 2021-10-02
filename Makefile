#CHANGE PATHS AS NEEDED:
INCLUDE_PATH = ${CONDA_PREFIX}/include
LIBRARY_PATH = ${CONDA_PREFIX}/lib

#compiler: gcc for C programs, g++ for C++ programs
XX = g++
CC = gcc

#compiler flags
CFLAGS = -g
GSLFLAGS = -lgsl -lgslcblas

#subset of files
NMF = modules/initialization.cpp modules/nmf.cpp modules/utils.cpp modules/io.cpp

all: clean run_mvnmf

matf:
	$(CC) -c -o modules/random_svd/matrix_funcs.o modules/random_svd/matrix_vector_functions_gsl.c -I${INCLUDE_PATH}

rsvd:
	$(CC) -c -o modules/random_svd/rsvd.o modules/random_svd/low_rank_svd_algorithms_gsl.c -I${INCLUDE_PATH}

run_mvnmf:
	$(XX) run_mvnmf.cpp modules/*.cpp modules/random_svd/*.o -o run_mvnmf $(CFLAGS) -L${LIBRARY_PATH} ${GSLFLAGS} -I${INCLUDE_PATH}

run_nmf:
	$(XX) run_nmf.cpp modules/*.cpp modules/random_svd/*.o -o run_nmf $(CFLAGS) -L${LIBRARY_PATH} ${GSLFLAGS} -I${INCLUDE_PATH}

clean:
	rm run_mvnmf
