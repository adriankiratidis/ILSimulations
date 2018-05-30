#### Shell script to run runSingleSphere.csh
#!/bin/bash

test_file_stub="testing"

ionic_liquid_name="NeutralDimers"
#ionic_liquid_name="SingleNeutralSpheres"
chi_parameter=0.71
Epsilon_r=2.3
Epsilon_LJ=0.4                         # ( = epsilon_LJ * k_{B})
surface_charge_density=0.003125
hs_diameter=2.4
a_term_index=2
bulk_density=0.04			# ( =n_{b} * [hs_diameter**3] )
temp=294.0			        #Temperature in Kelvin
bead_charge=0.2
string_length=1.2			#
n_points_per_hs_diameter=100			#number of discretised points
max_iteration_limit=100
iterative_tolerance=0.00000000001

starting_plate_separation=5
number_of_separations=400   #The number of plate separations we'll calculate.
number_of_plate_separations_in_hs_diameter=10

params_file=$test_file_stub.params

cat <<EOF > $params_file
${ionic_liquid_name}
${chi_parameter}
${Epsilon_r}
${Epsilon_LJ}
${surface_charge_density}
${hs_diameter}
${a_term_index}
${bulk_density}
${temp}
${bead_charge}
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

./bin/runSingleSphere.x <<EOF
$test_file_stub
EOF
