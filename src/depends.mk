kinds.o:
universalconstants.o: kinds.o
parameters.o: kinds.o
integratephispherical.o: kinds.o
iteration.o: kinds.o parameters.o
ILsimulationssrclib.o: kinds.o universalconstants.o parameters.o integratephispherical.o iteration.o
