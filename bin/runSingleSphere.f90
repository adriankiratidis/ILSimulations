!Program to run the hard sphere simulation for the ionic liquid [C4MIM+][BF4-]
program run_C4MIM_BF4
  use ILsimulationssrclib
  implicit none

  character(len=256) :: file_stub
  integer            :: iteration, ith_separation

  real(dp), dimension(:), allocatable :: n_neutral! bead densities
  real(dp), dimension(:), allocatable :: n_neutral_updated ! bead densities

  !Need zero array as some contributions depend on n_s = n_+ + n_0 + n_-
  !and n_+ and n_- and their associated lambdas are zero.
  real(dp), dimension(:), allocatable :: zero_array 

  !Here we define lambda_{i} = e^{l^{i}_{b} - l^{i}(r) where l^{i}(r) = dF/dn_{i} for example
  !and l^{i}_{b} is the value in the bulk.
  real(dp), dimension(:), allocatable :: lambda_neutral

  real(dp), dimension(:), allocatable :: grand_potential_per_unit_area
  real(dp), dimension(:), allocatable :: normal_pressure_left_wall, normal_pressure_right_wall
  real(dp), dimension(:), allocatable :: negative_deriv_of_potential

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

  integer :: start_z_index
  integer :: end_z_index

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
  call InitialisePotentialAndContactTheoremVariables(grand_potential_per_unit_area, &
       normal_pressure_left_wall, normal_pressure_right_wall, negative_deriv_of_potential)

  do ith_separation = 1, size(plate_separations)

     print *, ""
     print *, "****************************************"
     print *, "Starting calcaluation of ith_separation = ", ith_separation
     print *, "****************************************"
     print *, ""
     print *, "Initialising/ReInitialising Discretistion and setting integration ansatz."
     print *, "Doing this for the densities"
     call InitialiseDensityDiscretisationAndSetIntegrationAnsatz(ith_separation, n_neutral)

     print *, "Initialise/ReInitialise Discretisation for all the temperary variables we need."
     call InitialiseVariableDiscretisation(ith_separation, n_neutral_updated, lambda_neutral, zero_array)

     print *, "Starting iteration.  Searching for convergence of density profiles."
     iteration = 0
     do while (iteration < MAX_ITERATION_LIMIT)
        iteration = iteration + 1

        !Calculate the lambdas from the densities.
        zero_array = 0.0_dp
        call CalculateLambdas(zero_array, zero_array, lambda_neutral, n_neutral, zero_array, zero_array, ith_separation)

        call UpdateDensities(lambda_neutral, n_neutral_updated)

        ! Now test convergence
        if(converged(n_neutral_updated, n_neutral)) then

           print *, ""
           print *, "************************************************************"
           print *, "runC4MIMBF4.x: Density calculations successfully converged."
           print *, "took ", iteration, " iterations."
           print *, "writing out density values to file"
           print *, "************************************************************"
           print *, ""
           call WriteOutputFormattedAsFunctionOfPosition(n_neutral_updated, trim(file_stub), "n_neutral")
           exit

        else if(iteration == MAX_ITERATION_LIMIT) then

           print *, "runC4MIMBf4.x: iteration == MAX_ITERATION_LIMIT"
           print *, "Hit the iteration limit without converging"
           print *, "Increase the iteration limit"
           call abort

        else if(iteration > MAX_ITERATION_LIMIT) then

           print *, "runC4MIMBf4.x: iteration > MAX_ITERATION_LIMIT"
           print *, "This should never happen"
           print *, "Coding error...aborting..."
           call abort

        else !Update and proceed to the next iteration

           n_neutral = n_neutral_updated

        end if

     end do !end iteration loop

     print *, "Calculating grand potential per unit area value."
     zero_array = 0.0_dp
     call CalculateGrandPotentialValuePerUnitArea(zero_array, n_neutral_updated, zero_array, &
          ith_separation, grand_potential_per_unit_area(ith_separation))

     print *, "Calculating normal from from the contact theorem"
     zero_array = 0.0_dp
     call CalculateNormalPressureFromContactTheorem(zero_array, n_neutral_updated, zero_array, &
          normal_pressure_left_wall(ith_separation), normal_pressure_right_wall(ith_separation))

  end do !end loop over plate separation

  call CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation(grand_potential_per_unit_area, negative_deriv_of_potential)

  call WriteOutputFormattedAsFunctionOfPlateSeparation(grand_potential_per_unit_area, trim(file_stub), "potential-per-unit-area")
  call WriteOutputFormattedAsFunctionOfPlateSeparation(normal_pressure_left_wall, trim(file_stub), "normal-pressure-left-wall")
  call WriteOutputFormattedAsFunctionOfPlateSeparation(normal_pressure_right_wall, trim(file_stub), "normal-pressure-right-wall")
  call WriteOutputFormattedAsFunctionOfPlateSeparation(negative_deriv_of_potential, trim(file_stub), "negative_deriv_of_potential")

  call DeAllocateModelParams()
  call DeAllocateLocalVariables()

  print *, "runSingleSphere.x completed running succesfully."
  print *, ""

contains

  subroutine DeAllocateLocalVariables()

    if(allocated(n_neutral)) deallocate(n_neutral)
    if(allocated(n_neutral_updated)) deallocate(n_neutral_updated)
    if(allocated(lambda_neutral)) deallocate(lambda_neutral)
    if(allocated(grand_potential_per_unit_area)) deallocate(grand_potential_per_unit_area)
    if(allocated(normal_pressure_left_wall)) deallocate(normal_pressure_left_wall)
    if(allocated(normal_pressure_right_wall)) deallocate(normal_pressure_right_wall)
    if(allocated(negative_deriv_of_potential)) deallocate(negative_deriv_of_potential)
    
  end subroutine DeAllocateLocalVariables

end program run_C4MIM_BF4
