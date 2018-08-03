#### Shell script to run runSingleSphere.csh

test_file_stub="testing"

#ionic_liquid_name="NeutralDimers"
#ionic_liquid_name="SingleNeutralSpheres"
#ionic_liquid_name="C4MIM_BF4-"
ionic_liquid_name="PositiveMinusSpheres"
#ionic_liquid_name="PositiveNeutralDimerMinusSpheres"
#ionic_liquid_name="PositiveNeutralDoubleDimerMinusDimer"
chi_parameter=1.0
Epsilon_r=2.3
Epsilon_LJ=100                         # ( = epsilon_LJ * k_{B})
surface_charge_density_left_wall=0.003125
surface_charge_density_right_wall=-0.003125
hs_diameter=2.4
a_term_index=1
bulk_density=0.2			# ( =n_{b} * [hs_diameter**3] )
temp=294.0			        #Temperature in Kelvin
alpha_mixing_for_update=1.0
slope_for_initial_guess=0.00021
#slope_for_initial_guess=0.0000
positive_bead_charge=1.0
negative_bead_charge=-1.0
string_length=1.2			#
n_points_per_hs_diameter=50			#number of discretised points
max_iteration_limit=800
iterative_tolerance=0.0000000001

starting_plate_separation=5
number_of_separations=50   #The number of plate separations we'll calculate.
number_of_plate_separations_in_hs_diameter=50

params_file=$test_file_stub.params

cat <<EOF > $params_file
${ionic_liquid_name}
${chi_parameter}
${Epsilon_r}
${Epsilon_LJ}
${surface_charge_density_left_wall}
${surface_charge_density_right_wall}
${hs_diameter}
${a_term_index}
${bulk_density}
${temp}
${alpha_mixing_for_update}
${slope_for_initial_guess}
${positive_bead_charge}
${negative_bead_charge}
${string_length}
${n_points_per_hs_diameter}
${max_iteration_limit}
${iterative_tolerance}
${number_of_separations}
EOF

d=$(echo print "1.0/${number_of_plate_separations_in_hs_diameter}" | python)

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
