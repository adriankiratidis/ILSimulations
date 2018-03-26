!Program to run the hard sphere simulation for the ionic liquid [C4MIM+][BF4-]
program run_C4MIM_BF4
  use ILsimulationssrclib
  implicit none

  character(len=256) :: input_file_stub
  
  real(dp), dimension(:), allocatable :: n_plus, n_minus, n_neutral, ns_bar ! bead densities
  real(dp), dimension(:), allocatable :: lambda_plus, lambda_minus, lambda_neutral !dF/dn_{i}

  ! We use the standard notation of cj/aj to denote the contribution from bead j to the cation/anion
  ! as described in J. Phys. Chem C 2017, 121, 1742-1751. DOI: 10.1021/acs.jpcc.6b11491
  ! Here we use the notation c8c1 for example to denote the fact that due to symmetry c8 and c1
  ! are the same and consequently we don't need to calculate the same thing multiple times.
  ! Similarly, the contributions a1, a2, a3 and a4 are all identical due to symmetry.
  real(dp), dimension(:), allocatable :: c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp
  real(dp), dimension(:), allocatable :: a1a2a3a4, a5
  
  print *, "Please Enter the data file prefix"
  read(*,*) input_file_stub

  print *, "Reading and initialising model parameters"
  print *, "This includes discretisation params and simulation params"
  call InitialiseModelParameters(trim(input_file_stub))

  print *, "Initialising Discretistion"
  call InitialiseDiscretisation()
  
  print *, "Calculating Lambda_{+}, Lambda_{-} and Lambda_{0}"
  call CalculateLambdas()

  iteration = 0
  do while (iteration < MAX_ITERATION_LIMIT)
     iteration = iteration + 1
     
     c8c1 = integrate_phi_spherical(lambda_zero)

     c7 = integrate_phi_spherical(lambda_zero * c8c1)

     c6 = integrate_phi_spherical(lambda_zero * c7)

     c5 = integrate_phi_spherical(lambda_zero * c6)

     c4 = integrate_phi_spherical()
     
     call TestConvergence()
     

  end do


end program run_C4MIM_BF4



  
