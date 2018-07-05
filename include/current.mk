#File containing compiler options

#F90C = ifort
#OPTS = -O2 -ip -xAVX -assume byterecl -convert big_endian #-traceback -check all
#F90FLAGS = -cpp  $(OPTS)

F90C = gfortran
OPTS = -O -Wall -fcheck=all -g -fbacktrace -ffree-line-length-0
#OPTS = -g -Wall -Wextra -Wconversion -fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all -ffpe-trap=zero,overflow,underflow -finit-real=nan
F90FLAGS = -cpp  $(OPTS)
LINK = -L../lib/ -lsrc
INCLUDE =

%.o: %.f90
	$(F90C) $(F90FLAGS) $(INCLUDE) -o $@ -c $<
%.x: %.f90
	$(F90C) $(F90FLAGS) $(INCLUDE) -cpp -o $@ $< $(OBJS) $(LINK)
