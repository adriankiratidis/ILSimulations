kinds.o:
universalconstants.o: kinds.o
parameters.o: kinds.o
integratephispherical.o: kinds.o parameters.o universalconstants.o
integratezcylindrical.o: kinds.o parameters.o universalconstants.o
iteration.o: kinds.o parameters.o
lambdas.o: kinds.o parameters.o universalconstants.o
ILsimulationssrclib.o: kinds.o universalconstants.o parameters.o integratephispherical.o integratezcylindrical.o iteration.o lambdas.o
