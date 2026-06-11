# Makefile

FC = gfortran
FFLAGS = -O3 -Wall -Wextra -std=f2018 -pedantic

OBJ = \
	  fileio.o \
	  input.o \
	  kind.o \
	  main.o

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

main.o : kind.o input.o
kind.o : 
input.o : kind.o fileio.o
