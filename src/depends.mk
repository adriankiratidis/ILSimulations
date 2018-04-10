kinds.o:
helpers.o: parameters.o
universalconstants.o: kinds.o
parameters.o: kinds.o
integratephispherical.o: kinds.o parameters.o universalconstants.o helpers.o
integratezcylindrical.o: kinds.o parameters.o universalconstants.o helpers.o
iteration.o: kinds.o parameters.o
lambdas.o: kinds.o parameters.o universalconstants.o helpers.o
ILsimulationssrclib.o: kinds.o universalconstants.o parameters.o integratephispherical.o integratezcylindrical.o iteration.o lambdas.o helpers.o
