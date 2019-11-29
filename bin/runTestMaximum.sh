#### Shell script to run runSingleSphere.csh
#!/bin/bash


#previous_success=120
#previous_failure=140

#for i in "${charges[@]}"
#do

#bash bin/runSimulation15.sh $separation





size_of_empty_string=0

is_file_there=$(ls /home/adriank/postdoc/2018/branch/ILSimulations/run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/35.51-88.78_PositiveNeutralDoubleDimerMinusDimer-0.005-9.5238--0.003125 | grep "testing15-a2-potential-per-unit-area.txt")
size_of_string=${#is_file_there}

#echo $size_of_string

if [ "$size_of_string" -eq "$size_of_empty_string" ]; then
    echo "found nothing"
else
    echo "Finished correctly"
fi


#done
