include ../include/current.mk

COMDIR = .

# Objects list
OBJS =   \
	kinds.o \
	universalconstants.o \
	parameters.o \
	integratephispherical.o \
	integratezcylindrical.o \
	iteration.o \
	lambdas.o \
	helpers.o \
	io.o \
	normalisation.o \
	surfaceforces.o \
	functionalderivatives.o

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
