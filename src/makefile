include ../include/current.mk

COMDIR = .

# Objects list
OBJS =   \
	kinds.o \
	universalconstants.o \
	parameters.o \
	integratephispherical.o \
	iteration.o \
	lambdas.o

all:
	make ILsimulationssrclib.o
	make libsrc.a

ILsimulationssrclib.o:  $(OBJS)

libsrc.a:  $(OBJS)
	ar rvs ../lib/libsrc.a $(OBJS) ILsimulationssrclib.o

# make clean for new make
clean:
	rm -f *.o *.mod *~ *.x

include depends.mk
