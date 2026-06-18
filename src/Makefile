# Makefile

FC = gfortran
FFLAGS = -O3 -Wall -Wextra -std=f2018 -pedantic

OBJ = \
	  analysis.o \
	  conductivity_function.o \
	  constants.o \
	  fileio.o \
	  finite_difference.o \
	  geometry.o \
	  input.o \
	  kind.o \
	  linalg.o \
	  main.o \
	  output.o \
	  source_function.o

MOD = $(OBJ:o=mod)

LIBS =

EXEC = thermos.x

.PHONY : default all clean

default : all

all : $(EXEC)

$(EXEC) : $(OBJ)
	$(FC) $(FFLAGS) $^ -o $@ $(LIBS)

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	$(RM) $(OBJ)
	$(RM) $(MOD)
	$(RM) $(EXEC)

main.o : kind.o input.o geometry.o output.o finite_difference.o source_function.o conductivity_function.o analysis.o
kind.o : 
input.o : kind.o fileio.o output.o finite_difference.o
fileio.o :
geometry.o : kind.o output.o
output.o : fileio.o kind.o
finite_difference.o : kind.o linalg.o output.o source_function.o conductivity_function.o
linalg.o : kind.o
constants.o : kind.o
source_function.o : kind.o constants.o output.o
conductivity_function.o : kind.o output.o
analysis.o : kind.o output.o linalg.o
