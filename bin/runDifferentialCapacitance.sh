#### Shell script to run runSingleSphere.csh
#!/bin/bash

#declare -a charges=(-0.008 -0.0075 -0.007 -0.0065 -0.006 -0.0055 -0.005 -0.0045 -0.004 -0.0035 -0.003 -0.0025 -0.002 -0.0015 -0.001 -0.0005 0.0 0.0005 0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005 0.0055 0.006 0.0065 0.007 0.0075 0.008)

declare -a charges=(-0.0065 -0.006 -0.0055 -0.005 -0.0045 -0.004 -0.0035 -0.003 -0.0025 -0.002 -0.0015 -0.001 -0.0005 0.0 0.0005 0.001 0.0015 0.002 0.0025 0.003 0.0035 0.004 0.0045 0.005 0.0055 0.006 0.0065)

for i in "${charges[@]}"
do
surface_charge_density_left_wall="$i"

bash bin/runSimulation14.sh $surface_charge_density_left_wall


mv /home/adriank/postdoc/2018/branch/ILSimulations/run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-53.27_C4MIM_BF4--0.005-0-${surface_charge_density_left_wall}-hs1.0TESTINGDIFFERENTIALCAPACITANCE3/testing14-a2-electric_potential_and_charge.txt ~/postdoc/2018/branch/ILSimulations

done
