# Makefile

FC = gfortran
FFLAGS = -O3 -Wall -Wextra -std=f2018 -pedantic

OBJ = \
	  fileio.o \
	  geometry.o \
	  input.o \
	  kind.o \
	  main.o \
	  output.o

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

main.o : kind.o input.o geometry.o output.o
kind.o : 
input.o : kind.o fileio.o output.o
fileio.o :
geometry.o : kind.o output.o
output.o : fileio.o
