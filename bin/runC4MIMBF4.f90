!Program to run the hard sphere simulation for the ionic liquid [C4MIM+][BF4-]
program runSingleSphere
  use ILsimulationssrclib
  implicit none

  character(len=256) :: file_stub
  integer            :: iteration, ith_separation

  real(dp), dimension(:), allocatable :: n_plus, n_neutral, n_minus! bead densities
  real(dp), dimension(:), allocatable :: n_plus_previous, n_neutral_previous, n_minus_previous! bead densities
  real(dp), dimension(:), allocatable :: n_plus_updated, n_neutral_updated, n_minus_updated ! bead densities

  !Here we define lambda_{i} = e^{l^{i}_{b} - l^{i}(r) where l^{i}(r) = dF/dn_{i} for example
  !and l^{i}_{b} is the value in the bulk.
  real(dp), dimension(:), allocatable :: lambda_plus, lambda_neutral, lambda_minus

  real(dp), dimension(:), allocatable :: lambda_cation_centre, lambda_anion_centre
  real(dp), dimension(:), allocatable :: n_cation_centre, n_anion_centre

  real(dp), dimension(:), allocatable :: grand_potential_per_unit_area, grand_potential_per_unit_area_in_bulk
  real(dp), dimension(:), allocatable :: normal_pressure_left_wall, normal_pressure_right_wall
  real(dp), dimension(:), allocatable :: negative_deriv_of_potential

  real(dp), dimension(:), allocatable :: dispersion_particle_particle_adjust_to_contact_thm
  real(dp), dimension(:), allocatable :: Centre_Separation_Metric
  
  integer :: ij
  integer :: icharge
  real(dp) :: intermediate, ratio
  integer :: end_size
  logical :: abort_now

  abort_now = .false.
  ! We use the standard notation of cj/aj to denote the contribution from bead j to the cation/anion
  ! as described in J. Phys. Chem C 2017, 121, 1742-1751. DOI: 10.1021/acs.jpcc.6b11491
  ! Here we use the notation c8c1 for example to denote the fact that due to symmetry c8 and c1
  ! are the same and consequently we don't need to calculate the same thing multiple times.
  ! Similarly, the contributions a1, a2, a3 and a4 are all identical due to symmetry.
  ! Note that we adopt the primed convention used in the document writeup associated with this code,
  ! (which largely follows that of the aforementioned paper), and label these contributions with a
  ! trailing 'p'.
  ! The cation beads are numbered as (with the symbol in brackets denoting the sign of the charge.)
  !
  !                 9(+)
  !                 |
  ! 1(0) -- 2(+) -- 3(+) -- 4(+) -- 5(0) --- 6(0) -- 7(0) -- 8(0)
  !                 |
  !                10(+)
  !
  !while the anions are numbered as
  !
  !        2(-)
  !        |
  !1(-) -- 5(-) -- (3)(-)
  !        |
  !        4(-)
  !
  !***********************************************
  !***********************************************
  !***************BEGIN EXECUTION*****************
  !***********************************************
  !***********************************************

  print *, "Please Enter the data file prefix"
  read(*,*) file_stub

  print *, "Reading and initialising model parameters"
  print *, "This includes discretisation params and simulation params"
  call InitialiseModelParameters(trim(file_stub))

  print *, "Initialisiong grand potential and variables for contact theorem check."
  call InitialisePotentialAndContactTheoremVariables(grand_potential_per_unit_area, grand_potential_per_unit_area_in_bulk, &
       normal_pressure_left_wall, normal_pressure_right_wall, negative_deriv_of_potential, dispersion_particle_particle_adjust_to_contact_thm, Centre_Separation_Metric)

  do ith_separation = 1, size(plate_separations)

     print *, ""
     print *, "****************************************"
     print *, "Starting calcaluation of ith_separation = ", ith_separation
     print *, "****************************************"
     print *, ""
     print *, "Initialising/ReInitialising Discretistion and setting integration ansatz."
     print *, "Doing this for the densities"
     call InitialiseDensityDiscretisationAndSetIntegrationAnsatz(ith_separation, n_plus, n_neutral, n_minus, n_cation_centre, n_anion_centre)

     print *, "Initialise/ReInitialise Discretisation for all the temperary variables we need."
     call InitialiseVariableDiscretisation(ith_separation, n_plus_updated, lambda_plus, &
          n_neutral_updated, lambda_neutral, n_minus_updated, lambda_minus, n_plus_previous, n_neutral_previous, n_minus_previous, lambda_cation_centre, lambda_anion_centre)
     call SetToZero(n_plus_updated, lambda_plus, n_neutral_updated, lambda_neutral, n_minus_updated, lambda_minus, lambda_cation_centre, lambda_anion_centre)


     if(ith_separation == 1) then
        !call ImposeChargeNeutrality(n_plus, n_neutral, n_minus, Donnan_potential, abort_now)
     end if

     call InitialiseChargeIncrement()
     do icharge = 1, n_charge_iterations
        if(icharge > 1) call UpdateChargeIncrement()


        print *, "Starting iteration.  Searching for convergence of density profiles."
        iteration = 0
        do while (iteration < MAX_ITERATION_LIMIT)
           iteration = iteration + 1

           !n_neutral = 0.0_dp
           
           !print *, "n_plus = ", n_plus
           !print *, ""
           !print *, "lambda_plus 1 = ", lambda_plus
           !call sleep(1)
           !print *, "n_neutral = ", n_neutral
           !print *, "n_minus = ", n_minus
           !call abort()

           !print *, "*************************************"
           !print *, "iteration = ", iteration
           !print *, "*************************************"
           !print *, "n_plus integral = ", integrate_z_cylindrical(positive_bead_charge*n_plus, "all_z")
           !print *, "n_neutral integral = ", n_neutral_updated
           !print *, "n_minus integral = ", integrate_z_cylindrical(negative_bead_charge*n_minus, "all_z")
           !print *, "ITERATION**********************************************", iteration

           call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, &
                lambda_cation_centre, n_cation_centre, lambda_anion_centre, n_anion_centre, ith_separation)

           !print *, "lambda_plus 2 = ", lambda_plus
           !call abort()

           !lambda_plus = 0.0_dp
           !lambda_neutral = 0.0_dp
           !lambda_minus = 0.0_dp

           ! print *, "lambda_plus = ", lambda_plus
           ! print *, ""
           ! print *, "lambda_neutral = ", lambda_neutral
           ! print *, ""
           ! print *, "lambda_minus = ", lambda_minus
           ! print *, ""
           ! print *, "" 
           ! call abort()


           call UpdateDensities(n_plus, n_neutral, n_minus, lambda_plus, n_plus_updated, lambda_neutral, n_neutral_updated, lambda_minus, n_minus_updated, &
                lambda_cation_centre, n_cation_centre, lambda_anion_centre, n_anion_centre, Donnan_potential, iteration, abort_now)


           !call GetLambdaContributionsByTerm(ith_separation, n_plus, n_neutral, n_minus, n_cation_centre, n_anion_centre, &
           !     lambda_plus, lambda_neutral, lambda_minus, lambda_cation_centre, lambda_anion_centre)
           
           !Found some problem, but still want to print what we've calculated so far.
           if(abort_now) exit


           do ij = 1, size(n_plus_updated)
              if(isnan((n_plus_updated(ij)))) then
                 print *, "n_plus = ", n_plus
                 print *, ""
                 print *, "n_plus_updated = ", n_plus_updated
                 print *, "iteration = ", iteration
                 print *, "DENSITY HAS A NAN......ABORTING..."
                 call abort()
              end if
           end do



           ! Now test convergence
           if(converged(n_plus_updated, n_plus, n_neutral_updated, n_neutral, n_minus_updated, n_minus)) then

              print *, ""
              print *, "************************************************************"
              print *, "runSingleSphere.x: Density calculations successfully converged."
              print *, "took ", iteration, " iterations."
              print *, "writing out density values to file"
              print *, "************************************************************"
              print *, ""
              print *, "charge = ", positive_bead_charge, negative_bead_charge


              call GetLambdaContributionsByTerm(ith_separation, n_plus, n_neutral, n_minus, n_cation_centre, n_anion_centre, &
                   lambda_plus, lambda_neutral, lambda_minus, lambda_cation_centre, lambda_anion_centre)

              !print *, "ns_bulk = ", bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads
              !print *, "ns = ", n_plus_updated + n_neutral_updated + n_minus_updated
              
              !Perform this update if we get the solution in one iteration.
              !Possible in principle because we rescale the solution at the previous separation.
              n_plus = n_plus_updated
              n_neutral = n_neutral_updated
              n_minus = n_minus_updated
              
              !call CalculateDonnanPotential(n_plus, n_minus, Donnan_potential)

              call WriteOutputFormattedAsFunctionOfPosition(n_plus_updated, trim(file_stub), &
                   "n_plus_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge)))
              call WriteOutputFormattedAsFunctionOfPosition(n_neutral_updated, trim(file_stub), &
                   "n_neutral_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge)))
              call WriteOutputFormattedAsFunctionOfPosition(n_minus_updated, trim(file_stub), &
                   "n_minus_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge)))
              call WriteOutputFormattedAsFunctionOfPosition(n_cation_centre, trim(file_stub), &
                   "n_cation_centre"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge)))
              call WriteOutputFormattedAsFunctionOfPosition(n_anion_centre, trim(file_stub), &
                   "n_anion_centre"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge)))

              

              
              ! !Also print out the sum, n_s
              ! call WriteOutputFormattedAsFunctionOfPosition(n_plus_updated + n_neutral_updated + n_minus_updated, trim(file_stub), &
              !      "n_s_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge)))

              ! !Now print out n_s_bar for checking
              ! call WriteOutputFormattedAsFunctionOfPosition(calculate_n_sbar(n_plus_updated + n_neutral_updated + n_minus_updated), trim(file_stub), &
              !      "n_sbar_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge)))


              exit

           else if(iteration == MAX_ITERATION_LIMIT) then

              print *, "runSingleSphere.x: iteration == MAX_ITERATION_LIMIT"
              print *, "Hit the iteration limit without converging"
              print *, "Increase the iteration limit"
              abort_now = .true.
              !call abort()

           else if(iteration > MAX_ITERATION_LIMIT) then

              print *, "runSingleSphere.x: iteration > MAX_ITERATION_LIMIT"
              print *, "This should never happen"
              print *, "Coding error...aborting..."
              call abort()

           else !Update and proceed to the next iteration

              ! call WriteOutputFormattedAsFunctionOfPosition(n_plus_updated, trim(file_stub), &
              !      "n_plus_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge))//"iteration"//trim(str(iteration)))
              ! call WriteOutputFormattedAsFunctionOfPosition(n_neutral_updated, trim(file_stub), &
              !      "n_neutral_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge))//"iteration"//trim(str(iteration)))
              ! call WriteOutputFormattedAsFunctionOfPosition(n_minus_updated, trim(file_stub), &
              !      "n_minus_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge))//"iteration"//trim(str(iteration)))

              !call WriteOutputFormattedAsFunctionOfPosition(n_plus_updated + n_neutral_updated + n_minus_updated, trim(file_stub), &
              !     "n_s_separation"//trim(str(plate_separations(ith_separation)))//"iteration"//trim(str(iteration)))


              n_plus_previous = n_plus
              n_neutral_previous = n_neutral
              n_minus_previous = n_minus

              if(iteration < 50000) then
                 do ij = 1, size(n_plus)
                    if(n_plus_updated(ij) .gt. n_plus(ij)) then
                       intermediate = 2*n_plus(ij) - (n_plus(ij)**2)/n_plus_updated(ij)
                       ratio = n_cation_centre(ij)/n_plus_updated(ij)
                       
                       if(intermediate .lt. n_plus_updated(ij)) then
                          n_plus_updated(ij) = intermediate
                          n_cation_centre(ij) = ratio * n_plus_updated(ij)
                       end if
                    end if
                 end do
                 do ij = 1, size(n_neutral)
                    if(n_neutral_updated(ij) .gt. n_neutral(ij)) then
                       intermediate = 2*n_neutral(ij) - (n_neutral(ij)**2)/n_neutral_updated(ij)
                       if(intermediate .lt. n_neutral_updated(ij)) then
                          n_neutral_updated(ij) = intermediate
                       end if
                    end if
                 end do
                 do ij = 1, size(n_minus)
                    if(n_minus_updated(ij) .gt. n_minus(ij)) then
                       intermediate = 2*n_minus(ij) - (n_minus(ij)**2)/n_minus_updated(ij)
                       ratio = n_anion_centre(ij)/n_minus_updated(ij)
                       
                       if(intermediate .lt. n_minus_updated(ij)) then
                          n_minus_updated(ij) = intermediate
                          n_anion_centre(ij) = ratio * n_minus_updated(ij)
                       end if
                    end if
                 end do
              end if


              if(iteration >= 0) then
                 n_plus_updated = (alpha_mixing_for_update * n_plus_updated) + (1.0_dp - alpha_mixing_for_update) * n_plus_previous
                 n_neutral_updated = (alpha_mixing_for_update * n_neutral_updated) + (1.0_dp - alpha_mixing_for_update) * n_neutral_previous
                 n_minus_updated = (alpha_mixing_for_update * n_minus_updated) + (1.0_dp - alpha_mixing_for_update) * n_minus_previous
              end if


              n_plus = n_plus_updated
              n_neutral = n_neutral_updated
              n_minus = n_minus_updated


           end if

        end do !end iteration loop
        !Found some problem, but still want to print what we've calculated so far.
        if(abort_now) exit
     end do !end charge increment loop

     !call CalculateDonnanPotential(n_plus, n_minus, Donnan_potential)

     if(abort_now) exit


     print *, "Calculating cation-anion centre length difference"
     call CalculateCentreLengthDifference(Centre_Separation_Metric(ith_separation), n_cation_centre, n_anion_centre)
     
     print *, "Calculating grand potential per unit area value."
     call CalculateGrandPotentialValuePerUnitArea(ith_separation, grand_potential_per_unit_area(ith_separation), &
          size(n_neutral_updated), n_plus_updated, n_neutral_updated, n_minus_updated, n_cation_centre, n_anion_centre, Donnan_potential)

     print *, "Calculating normal pressure from the contact theorem"
     call CalculateNormalPressureFromContactTheorem(n_plus_updated, n_neutral_updated, n_minus_updated, &
          normal_pressure_left_wall(ith_separation), normal_pressure_right_wall(ith_separation), &
          dispersion_particle_particle_adjust_to_contact_thm(ith_separation))

     print *, "integral_plus = ", integrate_z_cylindrical(n_plus_updated * (hs_diameter**2), unity_function)
     print *, "integral_neutral = ", integrate_z_cylindrical(n_neutral_updated * (hs_diameter**2), unity_function)
     print *, "integral_minus = ", integrate_z_cylindrical(n_minus_updated * (hs_diameter**2), unity_function)

     print *, "Calculating Potential for Differential Capacitance"
     call CalculateAndWriteElectricPotential(trim(file_stub), plate_separations(ith_separation), n_plus_updated, n_minus_updated)

  end do !end loop over plate separation

  print *, "Separations...."
  print *, Centre_Separation_Metric
  
  !call MakeContactTheoremAdjustmentFromParticleParticleDispersion(normal_pressure_left_wall, normal_pressure_right_wall, dispersion_particle_particle_adjust_to_contact_thm)  

  if(abort_now) then
     end_size = ith_separation - 1
  else
     end_size = size(plate_separations)
  end if
  print *, "end_size = ", end_size

  negative_deriv_of_potential(1:end_size) = CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation(grand_potential_per_unit_area(1:end_size))

  do ith_separation = 1, end_size
     grand_potential_per_unit_area(ith_separation) = grand_potential_per_unit_area(ith_separation) + &
          (negative_deriv_of_potential(end_size) * (plate_separations(ith_separation)) * &
          (hs_diameter))
  end do

  do ith_separation = 1, end_size
     grand_potential_per_unit_area(ith_separation) = grand_potential_per_unit_area(ith_separation) !- &
     !( grand_potential_per_unit_area(size(negative_deriv_of_potential)))
  end do

  negative_deriv_of_potential(1:end_size) = CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation(grand_potential_per_unit_area(1:end_size))
  print *, "grand potential = ", grand_potential_per_unit_area(1:end_size)

  call WriteOutputFormattedAsFunctionOfPlateSeparation(grand_potential_per_unit_area(1:end_size), trim(file_stub), "potential-per-unit-area")
  call WriteOutputFormattedAsFunctionOfPlateSeparation(normal_pressure_left_wall(1:end_size), trim(file_stub), "normal-pressure-left-wall")
  call WriteOutputFormattedAsFunctionOfPlateSeparation(normal_pressure_right_wall(1:end_size), trim(file_stub), "normal-pressure-right-wall")
  call WriteOutputFormattedAsFunctionOfPlateSeparation(negative_deriv_of_potential(1:end_size), trim(file_stub), "negative_deriv_of_potential")
  
  call DeAllocateModelParams()
  call DeAllocateLocalVariables()

  print *, ""
  print *, "************************************************"
  print *, "runC4MIMBF4.x completed running succesfully."
  print *, "Completed ", size(plate_separations), " different plate separations."
  print *, "Now check output files to verify results via the contact theorem."
  print *, "runC4MIMBF4.x completed running succesfully."
  print *, "************************************************"
  print *, ""

contains

  subroutine DeAllocateLocalVariables()

    if(allocated(n_plus)) deallocate(n_plus)
    if(allocated(n_neutral)) deallocate(n_neutral)
    if(allocated(n_minus)) deallocate(n_minus)
    if(allocated(n_plus_updated)) deallocate(n_plus_updated)
    if(allocated(n_neutral_updated)) deallocate(n_neutral_updated)
    if(allocated(n_minus_updated)) deallocate(n_minus_updated)

    if(allocated(lambda_plus)) deallocate(lambda_plus)
    if(allocated(lambda_neutral)) deallocate(lambda_neutral)
    if(allocated(lambda_minus)) deallocate(lambda_minus)
    if(allocated(lambda_cation_centre)) deallocate(lambda_cation_centre)
    if(allocated(lambda_anion_centre)) deallocate(lambda_anion_centre)
    if(allocated(n_cation_centre)) deallocate(n_cation_centre)
    if(allocated(n_anion_centre)) deallocate(n_anion_centre)
    if(allocated(grand_potential_per_unit_area)) deallocate(grand_potential_per_unit_area)
    if(allocated(grand_potential_per_unit_area_in_bulk)) deallocate(grand_potential_per_unit_area_in_bulk)
    if(allocated(normal_pressure_left_wall)) deallocate(normal_pressure_left_wall)
    if(allocated(normal_pressure_right_wall)) deallocate(normal_pressure_right_wall)
    if(allocated(negative_deriv_of_potential)) deallocate(negative_deriv_of_potential)
    if(allocated(dispersion_particle_particle_adjust_to_contact_thm)) deallocate(dispersion_particle_particle_adjust_to_contact_thm)
    if(allocated(Centre_Separation_Metric)) deallocate(Centre_Separation_Metric)
    if(allocated(n_plus_previous)) deallocate(n_plus_previous)
    if(allocated(n_neutral_previous)) deallocate(n_neutral_previous)
    if(allocated(n_minus_previous)) deallocate(n_minus_previous)

  end subroutine DeAllocateLocalVariables

end program runSingleSphere
