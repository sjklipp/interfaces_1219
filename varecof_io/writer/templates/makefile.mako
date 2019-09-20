# ******************* FORTRAN Compiler ********************
FC = gfortran
FFLAGS = -O3 -fPIC

# ****************** C Compiler *********************
CC = gcc
CFLAGS = -O3 -pedantic

FOR_SOURCES = pot_aux.f dummy_corr.f ${mol_corr_potentials} 

OBJECTS = $(FOR_SOURCES:.f=.o)

.PHONY: clean tags 

libcorrpot.so : ${OBJECTS} 
	${FC} -shared -o $@ ${OBJECTS}  

clean:
	rm -f *.o
	rm -f *~
	rm -f *.d
	rm -f a.out

tags:
	etags *.cc *.hh *.f

