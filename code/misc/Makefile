#.SUFFIXES: .o .f90 .F90

SRC_dir = source
OBJ_dir = obj

MOD = numeric weno time
SRC = $(MOD:%=%.f90)
OBJ = $(MOD:%=$(OBJ_dir)/%.o)
.PRECIOUS: $(OBJ_dir)/%.o

FC = gfortran-8 -c -cpp
LD = gfortran-8

#FFLAGS = -O0 -Wall -Wunderflow -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow -g -p -fbounds-check -fbacktrace -fdump-core
FFLAGS = -O0 -Wall -Wunderflow -fcheck=all -ffpe-trap=invalid -fbounds-check -fbacktrace -fdump-core

Fmain = $(wildcard *.f90)
Fbin  = $(Fmain:%.f90=%.out)

all: | $(OBJ_dir) $(OBJ)

$(OBJ_dir)/numeric.o : numeric.f90
	@ echo "\033[36m$@\033[0m"
	$(FC) $(FFLAGS) -I$(OBJ_dir) -J$(OBJ_dir) -c $< -o $@

$(OBJ_dir)/weno.o : weno.f90 $(OBJ_dir)/numeric.o
	@ echo "\033[37m$@\033[0m"
	$(FC) $(FFLAGS) -I$(OBJ_dir) -J$(OBJ_dir) -c $< -o $@

$(OBJ_dir)/time.o : time.f90 $(OBJ_dir)/numeric.o
	@ echo "\033[37m$@\033[0m"
	$(FC) $(FFLAGS) -I$(OBJ_dir) -J$(OBJ_dir) -c $< -o $@

$(OBJ_dir)/%.o : %.f90 $(OBJ)
	@ echo "\033[37m$@\033[0m"
	$(FC) $(FFLAGS) -I$(OBJ_dir) -J$(OBJ_dir) -c $< -o $@

%.out : $(OBJ_dir)/%.o $(OBJ)
	@ echo "\033[33m$@\033[0m"
	$(LD) -L/usr/local/lib -ldfftpack $(FFLAGS) $^ -o $@



$(OBJ_dir) :
	@ echo "\033[95m$@\033[0m"
	mkdir -p $@

opt: clean
	$(MAKE) $(OBJ) $(filter-out $@,$(MAKECMDGOALS)) FFLAGS="-O3"

clean:
	rm -f $(OBJ_dir)/*.o $(OBJ_dir)/*.mod

mrproper: clean
	rm -f *.out
	rm -Rf *.out.dSYM


