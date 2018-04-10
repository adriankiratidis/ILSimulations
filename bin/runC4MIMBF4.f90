!Program to run the hard sphere simulation for the ionic liquid [C4MIM+][BF4-]
program run_C4MIM_BF4
  use ILsimulationssrclib
  implicit none

  character(len=256) :: file_stub
  integer            :: iteration, id

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
  real(dp), dimension(:), allocatable :: c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp ! Need for +ve beads
  real(dp), dimension(:), allocatable :: c2p, c4p, c5p, c6p, c7p ! neutral beads need these as well.
  real(dp), dimension(:), allocatable :: a1a2a3a4, a5p ! contributions needed for the -ve beads.

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

  do id = 1, size(plate_separations)

     print *, "Initialising/ReInitialising Discretistion and setting integration ansatz."
     print *, "Doing this for the densities"
     call InitialiseDensityDiscretisationAndSetIntegrationAnsatz(id, n_plus, n_minus, n_neutral)

     print *, "Initialise/ReInitialise Discretisation for all the temperary variables we need."
     call InitialiseVariableDiscretisation(id, n_plus_updated, n_minus_updated, n_neutral_updated, &
          lambda_plus, lambda_minus, lambda_neutral, c8c1, c9c10, c7, c6, c5, c4, c2, c3p, &
          c3pp, c3ppp, c2p, c4p, c5p, c6p, c7p, a1a2a3a4, a5p)

     iteration = 0
     do while (iteration < MAX_ITERATION_LIMIT)
        iteration = iteration + 1

        !Calculate the lambdas from the densities.
        call CalculateLambdas(lambda_plus, lambda_neutral, lambda_minus, n_plus, n_neutral, n_minus, id)

        ! print *, "lambda_plus = ", lambda_plus(1:10)
        ! print *, "lambda_neutral = ", lambda_neutral(1:10)
        ! print *, "lambda_minus = ", lambda_minus(1:10)
        ! print *, "n_plus(1) = ", n_plus(1)
        ! print *, "n_neutral(1) = ", n_neutral(1)
        ! print *, "n_minus(1) = ", n_minus(1)
        ! call abort()

        ! First we calculate the contributions from the required sites.

        ! Start with the calculation required for the positive bead site.
        c9c10 = integrate_phi_spherical(lambda_plus)

        c8c1 = integrate_phi_spherical(lambda_neutral)

        c7 = integrate_phi_spherical(lambda_neutral * c8c1)

        c6 = integrate_phi_spherical(lambda_neutral * c7)

        c5 = integrate_phi_spherical(lambda_neutral * c6)

        c4 = integrate_phi_spherical(lambda_plus * c5)

        c3p = integrate_phi_spherical(lambda_plus * c4 * c9c10 * c9c10)

        c2 = integrate_phi_spherical(lambda_plus * c8c1)

        c3pp = integrate_phi_spherical(lambda_plus * c4 * c2 * c9c10)

        c3ppp = integrate_phi_spherical(lambda_plus * c2 * c9c10 * c9c10)

        !Calculate the resulting positive bead densities.
        !n_plus_updated = nc2 + nc3 + nc4 + nc9 + nc10 =  nc2 + nc3 + nc4 + 2*nc9
        n_plus_updated = bulk_density * ( (lambda_plus * c8c1 * c3p) + (lambda_plus * c2 * c9c10 * c9c10 * c4) + &
             (lambda_plus * c3ppp * c5) + (2 * (lambda_plus * c3pp)) ) 

        !Now we calculate the neutral beads
        c2p = integrate_phi_spherical(lambda_plus * c3p)

        c4p = integrate_phi_spherical(lambda_plus * c3ppp)

        c5p = integrate_phi_spherical(lambda_neutral * c4p)

        c6p = integrate_phi_spherical(lambda_neutral * c5p)

        c7p = integrate_phi_spherical(lambda_neutral * c6p)

        !Calculate the resulting neutral bead densities.
        !n_zero_updated = nc1 + nc5 + nc6 + nc7 + nc8
        n_neutral_updated = bulk_density * ( (lambda_neutral * c2p) + (lambda_neutral * c4p * c6) + (lambda_neutral * c5p * c7) &
             + (lambda_neutral * c6p * c8c1) + (lambda_neutral * c7p) )

        !Calculate the required contributions for the anion
        a1a2a3a4 = integrate_phi_spherical(lambda_minus)

        a5p = integrate_phi_spherical(lambda_minus * (a1a2a3a4 ** 3.0_dp))

        !Calculate the resulting negative bead densities.
        !n_minus_updated = na1 + na2 + na3 + na4 + na5 = 4*na1 + na5
        n_minus_updated = bulk_density * ( 4.0_dp*(lambda_minus * a5p) + (lambda_minus * (a1a2a3a4**4.0_dp)) )

        print *, "n_plus(1) = ", n_plus(1:10)
        print *, "n_neutral(1) = ", n_neutral(1:10)
        print *, "n_minus(1) = ", n_minus(1:10)
        print *, "n_plus_updated(1) = ", n_plus_updated(1:10)
        print *, "n_netural_updated(1) = ", n_neutral_updated(1:10)
        print *, "n_minus_updated(1) = ", n_minus_updated(1:10)
        !call abort()
        
        ! Now test convergence
        if(converged(n_plus_updated, n_neutral_updated, n_minus_updated, n_plus, n_neutral, n_minus)) then

           print *, ""
           print *, "************************************************************"
           print *, "runC4MIMBF4.x: Density calculations successfully converged."
           print *, "writing out density values to file"
           print *, "************************************************************"
           print *, ""
           !call WriteDensityOutputFormatted(n_plus_updated, trim(file_stub), "n_plus")
           !call WriteDensityOutputFormatted(n_neutral_updated, trim(file_stub), "n_neutral")
           !call WriteDensityOutputFormatted(n_minus_updated, trim(file_stub), "n_minus")
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
    if(allocated(c8c1)) deallocate(c8c1)
    if(allocated(c9c10)) deallocate(c9c10)
    if(allocated(c7)) deallocate(c7)
    if(allocated(c6)) deallocate(c6)
    if(allocated(c5)) deallocate(c5)
    if(allocated(c4)) deallocate(c4)
    if(allocated(c2)) deallocate(c2)
    if(allocated(c3p)) deallocate(c3p)
    if(allocated(c3pp)) deallocate(c3pp)
    if(allocated(c3ppp)) deallocate(c3ppp)
    if(allocated(c2p)) deallocate(c2p)
    if(allocated(c4p)) deallocate(c4p)
    if(allocated(c5p)) deallocate(c5p)
    if(allocated(c6p)) deallocate(c6p)
    if(allocated(c7p)) deallocate(c7p)
    if(allocated(a1a2a3a4)) deallocate(a1a2a3a4)
    if(allocated(a5p)) deallocate(a5p)

  end subroutine DeAllocateLocalVariables

end program run_C4MIM_BF4
