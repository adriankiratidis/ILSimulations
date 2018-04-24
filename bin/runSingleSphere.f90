!Program to run the hard sphere simulation for the ionic liquid [C4MIM+][BF4-]
program run_C4MIM_BF4
  use ILsimulationssrclib
  implicit none

  character(len=256) :: file_stub
  integer            :: iteration, ith_separation

  real(dp), dimension(:), allocatable :: n_plus, n_minus, n_neutral! bead densities
  real(dp), dimension(:), allocatable :: n_plus_updated, n_minus_updated, n_neutral_updated ! bead densities

  !Here we define lambda_{i} = e^{l^{i}_{b} - l^{i}(r) where l^{i}(r) = dF/dn_{i} for example
  !and l^{i}_{b} is the value in the bulk.
  real(dp), dimension(:), allocatable :: lambda_plus, lambda_minus, lambda_neutral

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

  do ith_separation = 1, size(plate_separations)

     print *, "Initialising/ReInitialising Discretistion and setting integration ansatz."
     print *, "Doing this for the densities"
     call InitialiseDensityDiscretisationAndSetIntegrationAnsatz(ith_separation, n_plus, n_minus, n_neutral)

     print *, "Initialise/ReInitialise Discretisation for all the temperary variables we need."
     call InitialiseVariableDiscretisation(ith_separation, n_plus_updated, n_minus_updated, n_neutral_updated, &
          lambda_plus, lambda_minus, lambda_neutral)

     print *, "We don't calculate with hs_diameter/2 of the wall.  Therefore set it zero."
     print *, "This aids with plotting ease."
     call setNonCalculatedRegionToZero(n_plus, n_minus, n_neutral)

     print *, "Setting Bead Densities from the Bulk Ion Density."
     call SetSingleSphereBeadDensityFromBulkIonDensity()

     iteration = 0
     do while (iteration < MAX_ITERATION_LIMIT)
        iteration = iteration + 1

        !Calculate the lambdas from the densities.
        call CalculateLambdas(lambda_plus, lambda_neutral, lambda_minus, n_plus, n_neutral, n_minus, ith_separation)

        call UpdateSingleSphereArrangment(lambda_plus, lambda_neutral, lambda_minus, n_plus_updated, n_neutral_updated, n_minus_updated)

        call ReNormaliseToBulkDensity(n_plus_updated, n_neutral_updated, n_minus_updated)        

        ! Now test convergence
        if(converged(n_plus_updated, n_neutral_updated, n_minus_updated, n_plus, n_neutral, n_minus)) then

           print *, ""
           print *, "************************************************************"
           print *, "runC4MIMBF4.x: Density calculations successfully converged."
           print *, "took ", iteration, " iterations."
           print *, "writing out density values to file"
           print *, "************************************************************"
           print *, ""
           call WriteDensityOutputFormatted(n_plus_updated, trim(file_stub), "n_plus")
           call WriteDensityOutputFormatted(n_neutral_updated, trim(file_stub), "n_neutral")
           call WriteDensityOutputFormatted(n_minus_updated, trim(file_stub), "n_minus")
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

           n_plus = n_plus_updated
           n_neutral = n_neutral_updated
           n_minus = n_minus_updated

        end if

     end do !end iteration loop

     !call CalculateGrandPotentialValue(n_plus_updated, n_neutral_updated, n_minus_updated, grand_potential(ith_separation))

  end do !end loop over plate separation

  call DeAllocateModelParams()
  call DeAllocateLocalVariables()

contains

  subroutine DeAllocateLocalVariables()

    if(allocated(n_plus)) deallocate(n_plus)
    if(allocated(n_minus)) deallocate(n_minus)
    if(allocated(n_neutral)) deallocate(n_neutral)
    if(allocated(n_plus_updated)) deallocate(n_plus_updated)
    if(allocated(n_minus_updated)) deallocate(n_minus_updated)
    if(allocated(n_neutral_updated)) deallocate(n_neutral_updated)
    if(allocated(lambda_plus)) deallocate(lambda_plus)
    if(allocated(lambda_minus)) deallocate(lambda_minus)
    if(allocated(lambda_neutral)) deallocate(lambda_neutral)

  end subroutine DeAllocateLocalVariables

end program run_C4MIM_BF4
