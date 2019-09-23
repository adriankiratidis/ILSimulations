#### Shell script to run runSingleSphere.csh
#!/bin/bash

declare -a epsilons=(0)

for i in "${epsilons[@]}"
do
Epsilon_LJ_particle_wall="$i"

test_file_stub="testing8-a2"

#ionic_liquid_name="NeutralDimers"
#ionic_liquid_name="SingleNeutralSpheres"
ionic_liquid_name="C4MIM_BF4-"
#ionic_liquid_name="PositiveMinusSpheres"
#ionic_liquid_name="PositiveNeutralMinusSpheres"
#ionic_liquid_name="PositiveNeutralDimerMinusSpheres"
#ionic_liquid_name="C6MIM+_TFSI-_model1"
#ionic_liquid_name="C4MIM+_TFSI-_model1"
#ionic_liquid_name="PositiveNeutralDoubleDimerMinusDimer"
#ionic_liquid_name="Heptamers"
#ionic_liquid_name="Heptamer_SingleSphere"
#ionic_liquid_name="Hexamer_SingleSphere"
chi_parameter=0.71
Epsilon_r=14.0
Epsilon_LJ_particle_particle=35.51 # ( = epsilon_LJ * k_{B})
epsilon_eighth_power_const=140
type_of_iteration="r4"
#Epsilon_LJ_particle_wall=100
#mica_density=0.00505293988 #Mica density in particles/angstrom^3
mica_density=0.10611173748 #Mica density in particles/angstrom^3
#surface_charge_density_left_wall=-312500000000000000
#surface_charge_density_right_wall=-312500000000000000
#hs_diameter=0.00000000024
surface_charge_density_left_wall=-0.003125
surface_charge_density_right_wall=-0.003125
hs_diameter=2.4
hs_fraction=1.0
a_term_index=2
#bulk_density=0.04475			# ( =n_{b} * [hs_diameter**3] )
#bulk_density=0.005			# ( =n_{b} * [hs_diameter**3] )
bulk_density=0.01467			# ( =n_{b} * [hs_diameter**3] )
#bulk_density=0.000005			# ( =n_{b} * [hs_diameter**3] )
temp=294.0			        #Temperature in Kelvin
alpha_mixing_for_update=0.01
slope_for_initial_guess=0.000
#slope_for_initial_guess=0.0000
n_charge_iterations=1
positive_bead_charge=0.2
negative_bead_charge=-0.2
Donan_potential_intial_guess=0.0
string_length=1.2			#
n_points_per_hs_diameter=50			#number of discretised points
max_iteration_limit=10000
iterative_tolerance=0.0000000000001

starting_plate_separation=3
number_of_separations=350  #The number of plate separations we'll calculate.
number_of_plate_separations_in_hs_diameter=50

params_file=$test_file_stub.params

cat <<EOF > $params_file
${ionic_liquid_name}
${chi_parameter}
${Epsilon_r}
${Epsilon_LJ_particle_particle}
${Epsilon_LJ_particle_wall}
${epsilon_eighth_power_const}
${type_of_iteration}
${mica_density}
${surface_charge_density_left_wall}
${surface_charge_density_right_wall}
${hs_diameter}
${a_term_index}
${bulk_density}
${temp}
${alpha_mixing_for_update}
${slope_for_initial_guess}
${n_charge_iterations}
${positive_bead_charge}
${negative_bead_charge}
${Donan_potential_intial_guess}
${string_length}
${n_points_per_hs_diameter}
${max_iteration_limit}
${iterative_tolerance}
${number_of_separations}
EOF

d=$(echo print "1.0/${number_of_plate_separations_in_hs_diameter}" | python2.7)

counter=0
plate=${starting_plate_separation}

while [ $counter -le $number_of_separations ]
do
    
cat <<EOF >> ${params_file}
$plate
EOF

plate=$(echo "$plate + $d" | bc)
((counter++))
done

#./bin/runSingleSphere.x <<EOF

./bin/runC4MIMBF4.x <<EOF
$test_file_stub
EOF

write_dir="./run_results/compare_epsilonLJ/CentreCentrePotential/minuswalls/${Epsilon_LJ_particle_particle}-${Epsilon_LJ_particle_wall}_${ionic_liquid_name}-${bulk_density}-${epsilon_eighth_power_const}-${type_of_iteration}-${surface_charge_density_left_wall}-hs${hs_fraction}CHARGE"
mkdir -p ${write_dir}
mv ${test_file_stub}* ${write_dir}

done
