!Program to run the hard sphere simulation for the ionic liquid [C4MIM+][BF4-]
program runSingleSphere
  use ILsimulationssrclib
  implicit none

  character(len=256) :: file_stub
  integer            :: iteration, ith_separation

  real(dp), dimension(:), allocatable :: n_plus, n_neutral, n_minus! bead densities
  real(dp), dimension(:), allocatable :: n_plus_updated, n_neutral_updated, n_minus_updated ! bead densities

  !Here we define lambda_{i} = e^{l^{i}_{b} - l^{i}(r) where l^{i}(r) = dF/dn_{i} for example
  !and l^{i}_{b} is the value in the bulk.
  real(dp), dimension(:), allocatable :: lambda_plus, lambda_neutral, lambda_minus

  real(dp), dimension(:), allocatable :: grand_potential_per_unit_area, grand_potential_per_unit_area_in_bulk
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
       normal_pressure_left_wall, normal_pressure_right_wall, negative_deriv_of_potential)

  do ith_separation = 1, size(plate_separations)

     print *, ""
     print *, "****************************************"
     print *, "Starting calcaluation of ith_separation = ", ith_separation
     print *, "****************************************"
     print *, ""
     print *, "Initialising/ReInitialising Discretistion and setting integration ansatz."
     print *, "Doing this for the densities"
     call InitialiseDensityDiscretisationAndSetIntegrationAnsatz(ith_separation, n_plus, n_neutral, n_minus)

     print *, "Initialise/ReInitialise Discretisation for all the temperary variables we need."
     call InitialiseVariableDiscretisation(ith_separation, n_plus_updated, lambda_plus, &
          n_neutral_updated, lambda_neutral, n_minus_updated, lambda_minus)
     call SetToZero(n_plus_updated, lambda_plus, n_neutral_updated, lambda_neutral, n_minus_updated, lambda_minus)

     print *, "Starting iteration.  Searching for convergence of density profiles."
     iteration = 0
     do while (iteration < MAX_ITERATION_LIMIT)
        iteration = iteration + 1

        !Calculate the lambdas from the densities.
        !call SetToZero(zero_array1, zero_array2)

        call CalculateLambdas(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_separation, iteration)

        !print *, "lambda_neutral(1:50) = ",  lambda_neutral(1:50) 

        call UpdateDensities(lambda_plus, n_plus_updated, lambda_neutral, n_neutral_updated, lambda_minus, n_minus_updated)


        ! if(iteration == 10) then
        !print *, "n_neutral_updated = ", n_neutral_updated
        !print *, "n_neutral = ", n_neutral
        !    call abort()
        ! end if



        ! Now test convergence
        if(converged(n_neutral_updated, n_neutral)) then

           print *, ""
           print *, "************************************************************"
           print *, "runSingleSphere.x: Density calculations successfully converged."
           print *, "took ", iteration, " iterations."
           print *, "writing out density values to file"
           print *, "************************************************************"
           print *, ""

           !Perform this update if we get the solution in one iteration.
           !Possible in principle because we rescale the solution at the previous separation.
           n_plus = n_plus_updated
           n_neutral = n_neutral_updated
           n_minus = n_minus_updated

           call WriteOutputFormattedAsFunctionOfPosition(n_plus_updated, trim(file_stub), &
                "n_plus_separation"//str(plate_separations(ith_separation)))
           call WriteOutputFormattedAsFunctionOfPosition(n_neutral_updated, trim(file_stub), &
                "n_neutral_separation"//str(plate_separations(ith_separation)))
           call WriteOutputFormattedAsFunctionOfPosition(n_minus_updated, trim(file_stub), &
                "n_minus_separation"//str(plate_separations(ith_separation)))
           exit

        else if(iteration == MAX_ITERATION_LIMIT) then

           print *, "runSingleSphere.x: iteration == MAX_ITERATION_LIMIT"
           print *, "Hit the iteration limit without converging"
           print *, "Increase the iteration limit"
           call abort()

        else if(iteration > MAX_ITERATION_LIMIT) then

           print *, "runSingleSphere.x: iteration > MAX_ITERATION_LIMIT"
           print *, "This should never happen"
           print *, "Coding error...aborting..."
           call abort()

        else !Update and proceed to the next iteration

           call WriteOutputFormattedAsFunctionOfPosition(n_plus_updated, trim(file_stub), &
                "n_plus_separation"//trim(str(plate_separations(ith_separation)))//"iteration"//trim(str(iteration)))
           call WriteOutputFormattedAsFunctionOfPosition(n_neutral_updated, trim(file_stub), &
                "n_neutral_separation"//trim(str(plate_separations(ith_separation)))//"iteration"//trim(str(iteration)))
           call WriteOutputFormattedAsFunctionOfPosition(n_minus_updated, trim(file_stub), &
                "n_minus_separation"//trim(str(plate_separations(ith_separation)))//"iteration"//trim(str(iteration)))

           n_plus = n_plus_updated
           n_neutral = n_neutral_updated
           n_minus = n_minus_updated

        end if

     end do !end iteration loop

     print *, "Calculating grand potential per unit area value."
     call CalculateGrandPotentialValuePerUnitArea(ith_separation, grand_potential_per_unit_area(ith_separation), &
          size(n_neutral_updated), n_plus_updated, n_neutral_updated, n_minus_updated)

     print *, "Calculating normal pressure from the contact theorem"
     call CalculateNormalPressureFromContactTheorem(n_plus_updated, n_neutral_updated, n_minus_updated, &
          normal_pressure_left_wall(ith_separation), normal_pressure_right_wall(ith_separation))

     !print *, "n_neutral_updated = ", n_neutral_updated
     !print *, "normal_pressure_left_wall = ", normal_pressure_left_wall
     !call abort()
  end do !end loop over plate separation

  call CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation(grand_potential_per_unit_area, negative_deriv_of_potential)

  call WriteOutputFormattedAsFunctionOfPlateSeparation(grand_potential_per_unit_area, trim(file_stub), "potential-per-unit-area")
  call WriteOutputFormattedAsFunctionOfPlateSeparation(normal_pressure_left_wall, trim(file_stub), "normal-pressure-left-wall")
  call WriteOutputFormattedAsFunctionOfPlateSeparation(normal_pressure_right_wall, trim(file_stub), "normal-pressure-right-wall")
  call WriteOutputFormattedAsFunctionOfPlateSeparation(negative_deriv_of_potential, trim(file_stub), "negative_deriv_of_potential")

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
    if(allocated(grand_potential_per_unit_area)) deallocate(grand_potential_per_unit_area)
    if(allocated(grand_potential_per_unit_area_in_bulk)) deallocate(grand_potential_per_unit_area_in_bulk)
    if(allocated(normal_pressure_left_wall)) deallocate(normal_pressure_left_wall)
    if(allocated(normal_pressure_right_wall)) deallocate(normal_pressure_right_wall)
    if(allocated(negative_deriv_of_potential)) deallocate(negative_deriv_of_potential)

  end subroutine DeAllocateLocalVariables

end program runSingleSphere
