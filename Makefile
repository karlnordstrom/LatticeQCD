INCDIR=include
SRCDIR=src
LIBDIR=lib
CC=gcc
CFLAGS=-std=c99 -Ofast -march=native -fPIC
INCS=${HOME}/local/include
LIBS=${HOME}/local/lib

all: test.exe

${LIBDIR}/metropolis.o:	${SRCDIR}/metropolis.c ${INCDIR}/metropolis.h
	${CC} ${CFLAGS} -o $@ -c $< -I${INCDIR} -lm -I${INCS} -L${LIBS} -lmatrix

test.exe: ${SRCDIR}/test.c ${LIBDIR}/metropolis.o
	${CC} ${CFLAGS} -o $@ $^ -I${INCDIR} -lm -I${INCS} -L${LIBS} -lmatrix

clean:
	@rm ${LIBDIR}/metropolis.o
	@rm test.exe

