
include ../include/current.mk

# Objects, automatically includes all .f90 files in this dir
EXES = $(patsubst ./%.f90,./%.x,$(wildcard ./*.f90))

INCLUDE = -I../src/
LINK = -L../lib/ -lsrc

all:
	make -C ../src/
	make $(EXES)

clean:
	rm -f *.o *.mod *.x

include depends.mk

