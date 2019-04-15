!Program to run the hard sphere simulation for the ionic liquid [C4MIM+][BF4-]
program runSimulation
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

  real(dp), dimension(:), allocatable :: grand_potential_per_unit_area, grand_potential_per_unit_area_in_bulk
  real(dp), dimension(:), allocatable :: normal_pressure_left_wall, normal_pressure_right_wall
  real(dp), dimension(:), allocatable :: negative_deriv_of_potential

  real(dp), dimension(:), allocatable :: dispersion_particle_particle_adjust_to_contact_thm

  integer :: ij
  integer :: icharge

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
       normal_pressure_left_wall, normal_pressure_right_wall, negative_deriv_of_potential, dispersion_particle_particle_adjust_to_contact_thm)

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
          n_neutral_updated, lambda_neutral, n_minus_updated, lambda_minus, n_plus_previous, n_neutral_previous, n_minus_previous)
     call SetToZero(n_plus_updated, lambda_plus, n_neutral_updated, lambda_neutral, n_minus_updated, lambda_minus)

     call InitialiseChargeIncrement()
     do icharge = 1, n_charge_iterations
        if(icharge > 1) call UpdateChargeIncrement()


        print *, "Starting iteration.  Searching for convergence of density profiles."
        iteration = 0
        do while (iteration < MAX_ITERATION_LIMIT)
           iteration = iteration + 1

           if(iteration > 1) then
              n_plus = (alpha_mixing_for_update * n_plus) + (1.0_dp - alpha_mixing_for_update) * n_plus_previous
              n_neutral = (alpha_mixing_for_update * n_neutral) + (1.0_dp - alpha_mixing_for_update) * n_neutral_previous
              n_minus = (alpha_mixing_for_update * n_minus) + (1.0_dp - alpha_mixing_for_update) * n_minus_previous
           end if

           print *, "*************************************"
           print *, "iteration = ", iteration
           print *, "*************************************"
           call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_separation)
   
           call UpdateDensities(lambda_plus, n_plus_updated, lambda_neutral, n_neutral_updated, lambda_minus, n_minus_updated)

           do ij = 1, size(n_plus_updated)
              if(isnan((n_plus_updated(ij)))) then
                 print *, "n_plus_updated = ", n_plus_updated
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

              !Perform this update if we get the solution in one iteration.
              !Possible in principle because we rescale the solution at the previous separation.
              n_plus = n_plus_updated
              n_neutral = n_neutral_updated
              n_minus = n_minus_updated

              call WriteOutputFormattedAsFunctionOfPosition(n_plus_updated, trim(file_stub), &
                   "n_plus_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge)))
              call WriteOutputFormattedAsFunctionOfPosition(n_neutral_updated, trim(file_stub), &
                   "n_neutral_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge)))
              call WriteOutputFormattedAsFunctionOfPosition(n_minus_updated, trim(file_stub), &
                   "n_minus_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge)))

              !Also print out the sum, n_s
              call WriteOutputFormattedAsFunctionOfPosition(n_plus_updated + n_neutral_updated + n_minus_updated, trim(file_stub), &
                   "n_s_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge)))

              !Now print out n_s_bar for checking
              call WriteOutputFormattedAsFunctionOfPosition(calculate_n_sbar(n_plus_updated + n_neutral_updated + n_minus_updated), trim(file_stub), &
                   "n_sbar_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge)))


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
                   "n_plus_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge))//"iteration"//trim(str(iteration)))
              call WriteOutputFormattedAsFunctionOfPosition(n_neutral_updated, trim(file_stub), &
                   "n_neutral_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge))//"iteration"//trim(str(iteration)))
              call WriteOutputFormattedAsFunctionOfPosition(n_minus_updated, trim(file_stub), &
                   "n_minus_separation"//trim(str(plate_separations(ith_separation)))//"charge"//trim(str(icharge))//"iteration"//trim(str(iteration)))


              n_plus_previous = n_plus
              n_neutral_previous = n_neutral
              n_minus_previous = n_minus

              n_plus = n_plus_updated
              n_neutral = n_neutral_updated
              n_minus = n_minus_updated


           end if

        end do !end iteration loop

     end do !end charge increment loop


     print *, "Calculating grand potential per unit area value."
     call CalculateGrandPotentialValuePerUnitArea(ith_separation, grand_potential_per_unit_area(ith_separation), &
          size(n_neutral_updated), n_plus_updated, n_neutral_updated, n_minus_updated)

     print *, "Calculating normal pressure from the contact theorem"
     call CalculateNormalPressureFromContactTheorem(n_plus_updated, n_neutral_updated, n_minus_updated, &
          normal_pressure_left_wall(ith_separation), normal_pressure_right_wall(ith_separation), &
          dispersion_particle_particle_adjust_to_contact_thm(ith_separation))

     print *, "integral_plus = ", integrate_z_cylindrical(n_plus_updated * (hs_diameter**2), unity_function)
     print *, "integral_neutral = ", integrate_z_cylindrical(n_neutral_updated * (hs_diameter**2), unity_function)
     print *, "integral_minus = ", integrate_z_cylindrical(n_minus_updated * (hs_diameter**2), unity_function)


  end do !end loop over plate separation

  call MakeContactTheoremAdjustmentFromParticleParticleDispersion(normal_pressure_left_wall, normal_pressure_right_wall, dispersion_particle_particle_adjust_to_contact_thm)  

  negative_deriv_of_potential = CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation(grand_potential_per_unit_area)

  do ith_separation = 1, size(plate_separations)
     grand_potential_per_unit_area(ith_separation) = grand_potential_per_unit_area(ith_separation) + &
          (negative_deriv_of_potential(size(negative_deriv_of_potential)) * (plate_separations(ith_separation)) * &
          (hs_diameter))
  end do

  do ith_separation = 1, size(plate_separations)
     grand_potential_per_unit_area(ith_separation) = grand_potential_per_unit_area(ith_separation) - &
          ( grand_potential_per_unit_area(size(negative_deriv_of_potential)))
  end do

  negative_deriv_of_potential = CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation(grand_potential_per_unit_area)


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
    if(allocated(dispersion_particle_particle_adjust_to_contact_thm)) deallocate(dispersion_particle_particle_adjust_to_contact_thm)
    if(allocated(n_plus_previous)) deallocate(n_plus_previous)
    if(allocated(n_neutral_previous)) deallocate(n_neutral_previous)
    if(allocated(n_minus_previous)) deallocate(n_minus_previous)

  end subroutine DeAllocateLocalVariables

end program runSimulation

module charge
  use kinds
  use parameters
  use helpers
  use integratezcylindrical
  implicit none
  private

  public :: ImposeChargeNeutrality
  
  public :: InitialiseChargeIncrement
  public :: UpdateChargeIncrement
  
  real(dp) :: positive_bead_increment
  real(dp) :: negative_bead_increment

  real(dp) :: left_plate_increment
  real(dp) :: right_plate_increment
  
contains
  
  subroutine InitialiseChargeIncrement()

    positive_bead_increment = positive_bead_charge / real(n_charge_iterations,dp)
    negative_bead_increment = negative_bead_charge / real(n_charge_iterations,dp)

    positive_bead_charge = positive_bead_increment 
    negative_bead_charge = negative_bead_increment

    !left_plate_increment = surface_charge_density_left_wall / real(n_charge_iterations,dp)
    !right_plate_increment = surface_charge_density_right_wall / real(n_charge_iterations,dp)

    !surface_charge_density_left_wall = left_plate_increment
    !surface_charge_density_right_wall = right_plate_increment
    
  end subroutine InitialiseChargeIncrement

  
  subroutine UpdateChargeIncrement()
    
    positive_bead_charge = positive_bead_charge + positive_bead_increment
    negative_bead_charge = negative_bead_charge + negative_bead_increment

    !surface_charge_density_left_wall = surface_charge_density_left_wall + left_plate_increment
    !surface_charge_density_right_wall = surface_charge_density_right_wall + right_plate_increment
    
  end subroutine UpdateChargeIncrement

  subroutine ImposeChargeNeutrality(n_plus, n_minus)
    real(dp), dimension(:) :: n_plus
    real(dp), dimension(:) :: n_minus

    real(dp) :: a, b, c, d
    real(dp) :: y1, y2 !two solutions for two roots
    real(dp) :: Donnan_potential
    
    !In order to solve for the Donnan potential and impose charge neutrality we must
    !solve the equation
    !
    !
    if( ((surface_charge_density_left_wall + surface_charge_density_right_wall)/electric_charge) < 1.0E-10_dp) then
       !The walls have equal and opposite charge.  By symmetry arguments the density profiles are already electroneutral.
       !=> Donnan potential is zero.  Therefore, there is nothing to do.
    else if(((positive_bead_charge + negative_bead_charge)/electric_charge) < 1.0E-10_dp) then
       !We have the special case of the positive and negative bead charge being the same.
       !We can solve this analytically, as q_{-} = -q_{+}.

       a = integrate_z_cylindrical(n_plus*positive_bead_charge, "all_z")
       b = integrate_z_cylindrical(n_minus*negative_bead_charge, "all_z")
       c = surface_charge_density_left_wall + surface_charge_density_right_wall

       !Finding the Donnan potential to impose electroneutrality is now equivalent to solving
       !a*exp(x) + b*exp(-x) + c = 0; where x = exp(beta*positive_bead_charge*Psi_{Donnan}).

       !Letting y = exp(x), and using the quadratic fromula we get

       y1 = (-1.0_dp*c + sqrt(c**2 - 4.0_dp*a*b))/(2.0_dp*a)
       y2 = (-1.0_dp*c - sqrt(c**2 - 4.0_dp*a*b))/(2.0_dp*a)

       if(y1 < 0 .and. y2 < 0) then
          print *, "charge.f90: ImposeChargeNeutrality: "
          print *, "Both solutions for y1 and y2 are negative."
          print *, "Since this is supposed to be exp(x), this is not possible...aborting..."
          call abort()
       else if(y1 > 0 .and. y2 > 0) then
          print *, "charge.f90: ImposeChargeNeutrality: "
          print *, "Both solutions for y1 and y2 are positive."
          print *, "Don't know which one to pick...aborting..."
          call abort()
       else
          y1 = max(y1, y2)
       end if

       Donnan_potential = log(y1)/(beta*positive_bead_charge)
       
       n_plus = n_plus*exp(beta*positive_bead_charge*Donnan_potential)
       n_minus = n_minus*exp(beta*negative_bead_charge*Donnan_potential)
       
    else
       !Use Newton's method to solve. 
       call SolveWithNewtonsMethod(a,b,c, Donnan_potential)
    end if

  end subroutine ImposeChargeNeutrality
end module charge


!Contains routines to contstruct the various types of oligomers
!that we may wish to calculate.  These routines are present in the module
!as they need to be called in two different places in the code.
!Once during iteration in order to calculate the updated density,
!and then again post density calculation in order to calculate the
!ideal chain contribution to the free energy.
module constructoligomers
  use kinds
  use parameters
  use helpers
  use integratephispherical
  use integratezcylindrical
  use normalisation
  use lambdas
  use charge
  implicit none
  private

  public :: UpdateDensities

  public :: calculate_single_neutral_sphere_ideal_chain_term
  public :: calculate_neutral_dimers_ideal_chain_term
  public :: calculate_C4MIMBF4_ideal_chain_term
  public :: calculate_PositiveMinusSpheres_ideal_chain_term
  public :: calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term
  public :: calculate_PositiveNeutralDoubleDimerMinusDimer_ideal_chain_term

  public :: calculate_chem_potential_term_neutral_spheres
  public :: calculate_chem_potential_term_neutral_dimers
  public :: calculate_chem_potential_C4MIMBF4
  public :: calculate_chem_potential_PositiveMinusSpheres
  public :: calculate_chem_potential_PositiveNeutralDimerMinusSpheres
  public :: calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer
  !Private subroutines
  !UpdateSinglePositiveSphereDensity
  !UpdateSingleNeutralSphereDensity
  !UpdateSinglePositiveSphereDensity
  !UpdateSinglePositiveNeutralMinusSphereDensities
  !UpdateC4MIMPositiveBeadDensities
  !UpdateC4MINNeutralBeadDensities
  !UpdateC4MINNegativeBeadDensities

contains

  subroutine UpdateDensities(lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated)
    real(dp), dimension(:), intent(in) :: lambda1
    real(dp), dimension(:), intent(out) :: n1_updated

    real(dp), dimension(:), intent(in), optional :: lambda2
    real(dp), dimension(:), intent(out), optional :: n2_updated

    real(dp), dimension(:), intent(in), optional :: lambda3
    real(dp), dimension(:), intent(out), optional :: n3_updated

    integer :: ij

    if(trim(ionic_liquid_name) == "SingleNeutralSpheres") then

       if(present(lambda2) .or. present(n2_updated) .or. present(lambda3) .or. present(n3_updated)) then
          print *, "When doing single spheres should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       end if

       call UpdateSingleNeutralSphereDensity(lambda1, n1_updated)

    else if(trim(ionic_liquid_name) == "SinglePositiveSpheres") then

       if(present(lambda2) .or. present(n2_updated) .or. present(lambda3) .or. present(n3_updated)) then
          print *, "When doing single spheres should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       end if

       call UpdateSinglePositiveSphereDensity(lambda1, n1_updated)

    else if(trim(ionic_liquid_name) == "SingleNegativeSpheres") then

       if(present(lambda2) .or. present(n2_updated) .or. present(lambda3) .or. present(n3_updated)) then
          print *, "When doing single spheres should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       end if

       call UpdateSingleNegativeSphereDensity(lambda1, n1_updated)

    else if(trim(ionic_liquid_name) == "SinglePositiveNeutralMinusSpheres") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateSinglePositiveNeutralMinusSphereDensities(lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated)
       end if


    else if(trim(ionic_liquid_name) == "PositiveMinusSpheres") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdatePositiveMinusSphereDensities(lambda1, n1_updated, lambda3, n3_updated)
       end if

    else if(trim(ionic_liquid_name) == "PositiveNeutralDimerMinusSpheres") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdatePositiveNeutralDimerMinusSphereDensities(lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated)
       end if

    else if(trim(ionic_liquid_name) == "PositiveNeutralDoubleDimerMinusDimer") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdatePositiveNeutralDoubleDimerMinusDimerDensities(lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated)
       end if

    else if(trim(ionic_liquid_name) == "NeutralDimers") then

       if(present(lambda2) .or. present(n2_updated) .or. present(lambda3) .or. present(n3_updated)) then
          print *, "When doing Neutral Dimers should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       else
          call UpdateNeutralDimerDensity(lambda1, n1_updated)
       end if

    else if(trim(ionic_liquid_name) == "C4MIM_BF4-") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC4MIMPositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC4MIMNeutralBeadDensities(lambda1, lambda2, n2_updated)
          call UpdateBF4NegativeBeadDensities(lambda3, n3_updated)
       end if

    else

       print *, "constructoligomers.f90: UpdateDensities: "
       print *, "Unsupported ionic_liquid_name value of ", trim(ionic_liquid_name)
       print *, "...aborting..."
       call abort()
    end if

    !Imposes charge neutrality, by solving for the Donnan potential.
    if(present(n3_updated)) then !we have positive, neutral and negative spheres all present
       !Note: if only doing a single charge species then this conditional call needs updating.
       call ImposeChargeNeutrality(n1_updated, n3_updated)
    end if

    !Clutch at some straws - introduce an artificial cutoff.
    ! do ij = 1, size(n1_updated)
    !    if(n1_updated(ij)*(hs_diameter**3) > 0.7) n1_updated(ij) = 0.7/(hs_diameter**3)
    ! end do
    ! if(present(n2_updated)) then
    !    do ij = 1, size(n2_updated)
    !       if(n2_updated(ij)*(hs_diameter**3) > 0.7) n2_updated(ij) = 0.7/(hs_diameter**3)
    !    end do
    ! end if
    ! if(present(n3_updated)) then
    !    do ij = 1, size(n3_updated)
    !       if(n3_updated(ij)*(hs_diameter**3) > 0.7) n3_updated(ij) = 0.7/(hs_diameter**3)
    !    end do
    ! end if


  end subroutine UpdateDensities

  subroutine UpdateSinglePositiveSphereDensity(lambda_plus, n_plus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated

    if(size(lambda_plus) /= size(n_plus_updated)) then
       print *, "constructoligomers.f90: ConstructSingleSphereArrangment: "
       print *, "Size mismatch.  size(lambda_plus) /= size(n_plus_updated)."
       print *, "can't update...aborting..."
       call abort()
    else
       n_plus_updated = bulk_density_positive_beads * exp(lambda_plus)
    end if

    call setNonCalculatedRegionToZero(n_plus_updated)

  end subroutine UpdateSinglePositiveSphereDensity

  subroutine UpdateSingleNeutralSphereDensity(lambda_neutral, n_neutral_updated)
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    if(size(lambda_neutral) /= size(n_neutral_updated)) then
       print *, "constructoligomers.f90: UpdateSingleNeutralSphereDensity: "
       print *, "Size mismatch.  size(lambda_neutral) /= size(n_plus_neutral)."
       print *, "can't update...aborting..."
       call abort()
    else
       n_neutral_updated = bulk_density_neutral_beads * exp(lambda_neutral)
    end if

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateSingleNeutralSphereDensity

  subroutine UpdateSingleNegativeSphereDensity(lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    if(size(lambda_minus) /= size(n_minus_updated)) then
       print *, "constructoligomers.f90: UpdateSingleNegativeSphereDensity: "
       print *, "Size mismatch.  size(lambda_plus) /= size(n_plus_updated)."
       print *, "can't update...aborting..."
       call abort()
    else
       n_minus_updated = bulk_density_negative_beads * exp(lambda_minus)
    end if

    call setNonCalculatedRegionToZero(n_minus_updated)

  end subroutine UpdateSingleNegativeSphereDensity

  subroutine UpdateNeutralDimerDensity(lambda_neutral, n_neutral_updated)
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    real(dp), dimension(:), allocatable :: c1

    ! First check the input variables are the same size
    if( (size(lambda_neutral) == size(n_neutral_updated)) ) then
       allocate(c1(size(lambda_neutral)))
    else
       print *, "constructoligomers.f90: UpdateNeutralDimerDensity:"
       print *, "Size mismatch. size(lambda_neutral) /= size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    ! Note that there is a factor of 1.0_dp/(4.0_dp * pi * hs_diameter**2)
    ! from the delta function bond, a factor of 2 * pi from the theta integral
    ! and a factor of hs_diameter**2 from the jacobian.  Combining these factors
    ! gives the required 0.5_dp factor that we are multiplying by.
    !c1 = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral))
    c1 = integrate_phi_spherical(exp(lambda_neutral))
    !c1 = 0.5_dp * (exp(lambda_neutral))

    !Note the factor of 2 is present as this formulation calculates the density of
    !an individual bead.  By symmetry we need both the beads are the same and to
    !get the total bead density we must add them.  Hence the factor of 2.
    n_neutral_updated = 2.0_dp * bulk_density * exp(lambda_neutral) * c1

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(c1)

    call setNonCalculatedRegionToZero(n_neutral_updated)

  end subroutine UpdateNeutralDimerDensity

  subroutine UpdatePositiveMinusSphereDensities(lambda_plus, n_plus_updated, lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    integer :: ij

    n_plus_updated = bulk_density * exp(lambda_plus)
    n_minus_updated = bulk_density * exp(lambda_minus)

    do ij = 1, (size(n_plus_updated) - 1)/2
       n_plus_updated(ij) = n_plus_updated(size(n_plus_updated) - ij + 1)
       n_minus_updated(ij) = n_minus_updated(size(n_minus_updated) - ij + 1)
    end do
    
    call setNonCalculatedRegionToZero(n_plus_updated)
    call setNonCalculatedRegionToZero(n_minus_updated)

  end subroutine UpdatePositiveMinusSphereDensities

  subroutine UpdatePositiveNeutralDimerMinusSphereDensities(lambda_plus, n_plus_updated, lambda_neutral, n_neutral_updated, lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    real(dp), dimension(:), allocatable :: c1, c2

    ! First check the input variables are the same size
    if( (size(lambda_plus) == size(n_plus_updated)) .and. (size(lambda_neutral) == size(n_neutral_updated)) .and. (size(lambda_minus) == size(n_minus_updated)) ) then
       allocate(c1(size(lambda_plus)))
       allocate(c2(size(lambda_plus)))
    else
       print *, "constructoligomers.f90: UpdatePositiveNeutralDimerMinusSphereDensities:"
       print *, "Size mismatch. The following expression is false."
       print *, "(size(lambda_plus) == size(n_plus_updated)) .and. (size(lambda_neutral) == size(n_neutral_updated)) .and. (size(lambda_minus) == size(n_minus_updated))"
       print *, "...aborting..."
       call abort()
    end if

    c1 = integrate_phi_spherical(exp(lambda_plus))
    c2 = integrate_phi_spherical(exp(lambda_neutral))

    n_plus_updated = bulk_density * exp(lambda_plus) * c2
    n_neutral_updated = bulk_density * exp(lambda_neutral) * c1

    n_minus_updated = bulk_density * exp(lambda_minus)

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(c1, c2)

    call setNonCalculatedRegionToZero(n_plus_updated)
    call setNonCalculatedRegionToZero(n_neutral_updated)
    call setNonCalculatedRegionToZero(n_minus_updated)

  end subroutine UpdatePositiveNeutralDimerMinusSphereDensities

  subroutine UpdatePositiveNeutralDoubleDimerMinusDimerDensities(lambda_plus, n_plus_updated, lambda_neutral, n_neutral_updated, lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    real(dp), dimension(:), allocatable :: c1, c2, c3, c4

    ! First check the input variables are the same size
    if( (size(lambda_plus) == size(n_plus_updated)) .and. (size(lambda_neutral) == size(n_neutral_updated)) .and. (size(lambda_minus) == size(n_minus_updated)) ) then
       allocate(c1(size(lambda_plus)))
       allocate(c2(size(lambda_plus)))
       allocate(c3(size(lambda_plus)))
       allocate(c4(size(lambda_plus)))       
    else
       print *, "constructoligomers.f90: UpdatePositiveNeutralDimerMinusSphereDensities:"
       print *, "Size mismatch. The following expression is false."
       print *, "(size(lambda_plus) == size(n_plus_updated)) .and. (size(lambda_neutral) == size(n_neutral_updated)) .and. (size(lambda_minus) == size(n_minus_updated))"
       print *, "...aborting..."
       call abort()
    end if

    c1 = integrate_phi_spherical(exp(lambda_minus))
    n_minus_updated = 2.0_dp * bulk_density * exp(lambda_minus) * c1

    c4 = integrate_phi_spherical(exp(lambda_neutral))
    c3 = integrate_phi_spherical(exp(lambda_neutral) * c4)
    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c1 = integrate_phi_spherical(exp(lambda_plus))

    n_plus_updated = bulk_density * ((exp(lambda_plus) * c2) + (exp(lambda_plus) * c3 * c1))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c1)
    c3 = integrate_phi_spherical(exp(lambda_neutral) * c2)

    n_neutral_updated = bulk_density * ( (exp(lambda_neutral) * c3) + (exp(lambda_neutral) * c2 * c4) )

    deallocate(c1, c2, c3, c4)

    call setNonCalculatedRegionToZero(n_plus_updated)
    call setNonCalculatedRegionToZero(n_neutral_updated)
    call setNonCalculatedRegionToZero(n_minus_updated)

  end subroutine UpdatePositiveNeutralDoubleDimerMinusDimerDensities


  function calculate_chem_potential_term_neutral_spheres(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_term_neutral_spheres

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp) :: lambda

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_term_neutral_spheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    calculate_chem_potential_term_neutral_spheres = (1.0_dp/beta) * (log(bulk_density_neutral_beads) + lambda) * &
         integrate_z_cylindrical(n_neutral, unity_function)

  end function calculate_chem_potential_term_neutral_spheres


  function calculate_single_neutral_sphere_ideal_chain_term(n_neutral)
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp) :: calculate_single_neutral_sphere_ideal_chain_term

    real(dp), dimension(size(n_neutral)) :: integrand

    integer :: start_z_index
    integer :: end_z_index

    integrand(:) = 0.0_dp
    call get_allowed_z_values(start_z_index, end_z_index, size(n_neutral))

    integrand(start_z_index:end_z_index) = n_neutral(start_z_index:end_z_index) * (log(n_neutral(start_z_index:end_z_index)) - 1.0_dp)

    calculate_single_neutral_sphere_ideal_chain_term = integrate_z_cylindrical(integrand, unity_function) / beta

  end function calculate_single_neutral_sphere_ideal_chain_term

  function calculate_chem_potential_term_neutral_dimers(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_term_neutral_dimers

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp) :: lambda

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_term_neutral_spheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)

    calculate_chem_potential_term_neutral_dimers = (1.0_dp/beta) * (log(bulk_density) + (2.0_dp * lambda)) * &
         integrate_z_cylindrical(0.5_dp * n_neutral, unity_function)

  end function calculate_chem_potential_term_neutral_dimers

  function calculate_neutral_dimers_ideal_chain_term(lambda_neutral)
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp) :: calculate_neutral_dimers_ideal_chain_term

    real(dp), dimension(size(lambda_neutral)) :: integrand
    real(dp), dimension(size(lambda_neutral)) :: integrand_with_lambda

    integer :: start_z_index
    integer :: end_z_index

    integrand(:) = 0.0_dp
    integrand_with_lambda(:) = 0.0_dp

    integrand(:) = integrate_phi_spherical(exp(lambda_neutral))
    integrand_with_lambda(:) = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    calculate_neutral_dimers_ideal_chain_term = bulk_density * integrate_z_cylindrical(&
         (integrand_with_lambda(:) + (integrand(:) * (lambda_neutral + log(bulk_density) - 1.0_dp))) * exp(lambda_neutral), unity_function ) / beta

  end function calculate_neutral_dimers_ideal_chain_term


  function calculate_chem_potential_PositiveMinusSpheres(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_PositiveMinusSpheres

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp) :: lambda_plus_bulk, lambda_minus_bulk

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_plus(start_z_index:end_z_index) - lambda_plus(start_z_index) < 0.000001_dp)) then
       lambda_plus_bulk = lambda_plus(start_z_index)
    else
       print *, "lambda_plus = ", lambda_plus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_minus(start_z_index:end_z_index) - lambda_minus(start_z_index) < 0.000001_dp)) then
       lambda_minus_bulk = lambda_minus(start_z_index)
    else
       print *, "lambda_minus = ", lambda_minus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    calculate_chem_potential_PositiveMinusSpheres = &
         ((1.0_dp/beta) * (log(bulk_density) + lambda_plus_bulk) * &
         integrate_z_cylindrical(n_plus, unity_function)) + &
         ((1.0_dp/beta) * (log(bulk_density) + lambda_minus_bulk) * &
         integrate_z_cylindrical(n_minus, unity_function))

  end function calculate_chem_potential_PositiveMinusSpheres

  function calculate_PositiveMinusSpheres_ideal_chain_term(n_plus, n_minus)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_minus

    real(dp) :: calculate_PositiveMinusSpheres_ideal_chain_term

    real(dp), dimension(size(n_plus)) :: integrand_plus, integrand_minus

    integer :: start_z_index
    integer :: end_z_index

    integrand_plus(:) = 0.0_dp
    integrand_minus(:) = 0.0_dp
    call get_allowed_z_values(start_z_index, end_z_index, size(n_plus))

    integrand_plus(start_z_index:end_z_index) = n_plus(start_z_index:end_z_index) * (log(n_plus(start_z_index:end_z_index)) - 1.0_dp)
    integrand_minus(start_z_index:end_z_index) = n_minus(start_z_index:end_z_index) * (log(n_minus(start_z_index:end_z_index)) - 1.0_dp)

    calculate_PositiveMinusSpheres_ideal_chain_term = integrate_z_cylindrical(integrand_plus, unity_function) / beta +&
         integrate_z_cylindrical(integrand_minus, unity_function) / beta

  end function calculate_PositiveMinusSpheres_ideal_chain_term


  function calculate_chem_potential_PositiveNeutralDimerMinusSpheres(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_PositiveNeutralDimerMinusSpheres

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp) :: lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_plus(start_z_index:end_z_index) - lambda_plus(start_z_index) < 0.000001_dp)) then
       lambda_plus_bulk = lambda_plus(start_z_index)
    else
       print *, "lambda_plus = ", lambda_plus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda_neutral_bulk = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if


    if(all(lambda_minus(start_z_index:end_z_index) - lambda_minus(start_z_index) < 0.000001_dp)) then
       lambda_minus_bulk = lambda_minus(start_z_index)
    else
       print *, "lambda_minus = ", lambda_minus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    calculate_chem_potential_PositiveNeutralDimerMinusSpheres = &
         ((1.0_dp/beta) * (log(bulk_density) + lambda_plus_bulk + lambda_neutral_bulk) * &
         integrate_z_cylindrical((n_plus + n_neutral)/2.0_dp, unity_function)) + &
         ((1.0_dp/beta) * (log(bulk_density) + lambda_minus_bulk) * &
         integrate_z_cylindrical(n_minus, unity_function))

  end function calculate_chem_potential_PositiveNeutralDimerMinusSpheres


  function calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term(lambda_plus, lambda_neutral, n_minus)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    real(dp) :: calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term

    real(dp), dimension(size(lambda_neutral)) :: integrand
    real(dp), dimension(size(lambda_neutral)) :: integrand_with_lambda

    integer :: start_z_index
    integer :: end_z_index

    integrand(:) = 0.0_dp
    call get_allowed_z_values(start_z_index, end_z_index, size(n_minus))

    integrand(start_z_index:end_z_index) = n_minus(start_z_index:end_z_index) * (log(n_minus(start_z_index:end_z_index)) - 1.0_dp)

    calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term = integrate_z_cylindrical(integrand, unity_function) / beta

    integrand(:) = 0.0_dp
    integrand_with_lambda(:) = 0.0_dp

    integrand(:) = integrate_phi_spherical(exp(lambda_neutral))
    integrand_with_lambda(:) = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term = calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term + (bulk_density * integrate_z_cylindrical(&
         (integrand_with_lambda(:) + (integrand(:) * (lambda_plus + log(bulk_density) - 1.0_dp))) * exp(lambda_plus), unity_function ) / beta)


  end function calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term







  function calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer

    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp) :: lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_plus(start_z_index:end_z_index) - lambda_plus(start_z_index) < 0.000001_dp)) then
       lambda_plus_bulk = lambda_plus(start_z_index)
    else
       print *, "lambda_plus = ", lambda_plus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda_neutral_bulk = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_minus(start_z_index:end_z_index) - lambda_minus(start_z_index) < 0.000001_dp)) then
       lambda_minus_bulk = lambda_minus(start_z_index)
    else
       print *, "lambda_minus = ", lambda_minus
       print *, "constructoligomers.f90: calculate_chem_potential_PositiveMinusSpheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)

    calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer = (1.0_dp/beta) * (&
         (log(bulk_density) + ((2.0_dp * lambda_plus_bulk) + (2.0_dp * lambda_neutral_bulk))) * &
         (integrate_z_cylindrical((n_plus + n_neutral)/4.0_dp, unity_function)) + &
         (integrate_z_cylindrical(n_minus/2.0_dp, unity_function) * (log(bulk_density) + 2.0_dp * lambda_minus_bulk)))


  end function calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer

  function calculate_PositiveNeutralDoubleDimerMinusDimer_ideal_chain_term(lambda_plus, lambda_neutral, lambda_minus)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp) :: calculate_PositiveNeutralDoubleDimerMinusDimer_ideal_chain_term

    real(dp), dimension(size(lambda_neutral)) :: c1, c2, c3, c4, a1, a2
    real(dp), dimension(size(lambda_neutral)) :: c1_lambda, c2_lambda, c3_lambda, c4_lambda, a1_lambda, a2_lambda
    real(dp), dimension(size(lambda_neutral)) :: cation_integrand, anion_integrand

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    cation_integrand  = 0.0_dp
    anion_integrand = 0.0_dp
    
    c4 = integrate_phi_spherical(exp(lambda_neutral))
    c4_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c3 = integrate_phi_spherical(exp(lambda_neutral) * c4)
    c3_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c4_lambda + c4*lambda_neutral))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp) + c2*lambda_plus
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_plus) * cation_integrand, unity_function)
    
    
    a1 = integrate_phi_spherical(exp(lambda_minus))
    a1_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)
    
    anion_integrand = a1_lambda + a1*(log(bulk_density) - 1.0_dp) + a1*lambda_minus
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)
    
    calculate_PositiveNeutralDoubleDimerMinusDimer_ideal_chain_term = anion_contribution + cation_contribution

  end function calculate_PositiveNeutralDoubleDimerMinusDimer_ideal_chain_term


  function calculate_chem_potential_C4MIMBF4(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_C4MIMBF4


    real(dp), dimension(size(n_plus)) :: lambda_plus
    real(dp), dimension(size(n_neutral)) :: lambda_neutral
    real(dp), dimension(size(n_minus)) :: lambda_minus

    real(dp) :: lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)

    !Check that lambda_bulk is the same everywhere.
    if(all(lambda_plus(start_z_index:end_z_index) - lambda_plus(start_z_index) < 0.000001_dp)) then
       lambda_plus_bulk = lambda_plus(start_z_index)
    else
       print *, "lambda_plus = ", lambda_plus
       print *, "constructoligomers.f90: calculate_chem_potential_C4MIMBF4: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_neutral(start_z_index:end_z_index) - lambda_neutral(start_z_index) < 0.000001_dp)) then
       lambda_neutral_bulk = lambda_neutral(start_z_index)
    else
       print *, "lambda_neutral = ", lambda_neutral
       print *, "constructoligomers.f90: calculate_chem_potential_C4MIMBF4: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    if(all(lambda_minus(start_z_index:end_z_index) - lambda_minus(start_z_index) < 0.000001_dp)) then
       lambda_minus_bulk = lambda_minus(start_z_index)
    else
       print *, "lambda_minus = ", lambda_minus
       print *, "constructoligomers.f90: calculate_chem_potential_C4MIMBF4: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)


    calculate_chem_potential_C4MIMBF4 = (1.0_dp/beta) * (&
         (log(bulk_density) + ((5.0_dp * lambda_plus_bulk) + (5.0_dp * lambda_neutral_bulk))) * &
         (integrate_z_cylindrical((n_plus + n_neutral)/10.0_dp, unity_function)) + &
         (integrate_z_cylindrical(n_minus/5.0_dp, unity_function) * (log(bulk_density) + 5.0_dp * lambda_minus_bulk)))


  end function calculate_chem_potential_C4MIMBF4

  function calculate_C4MIMBF4_ideal_chain_term(lambda_plus, lambda_neutral, lambda_minus)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus

    real(dp) :: calculate_C4MIMBF4_ideal_chain_term

    real(dp) :: cation_contribution
    real(dp) :: anion_contribution

    integer :: array_size

    real(dp), dimension(:), allocatable :: c8, c7, c6, c5, c4, c3, c2
    real(dp), dimension(:), allocatable :: c910, c4p
    real(dp), dimension(:), allocatable :: c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, c3_lambda, c2_lambda
    real(dp), dimension(:), allocatable :: c910_lambda, c4p_lambda

    real(dp), dimension(:), allocatable :: a1234
    real(dp), dimension(:), allocatable :: a1234_lambda

    real(dp), dimension(:), allocatable :: anion_integrand, cation_integrand

    array_size = size(lambda_plus)

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) == size(lambda_neutral)) .or. (size(lambda_plus) == size(lambda_minus))) then
       allocate(c8(array_size), c7(array_size), c6(array_size), c5(array_size), c4(array_size), c3(array_size), c2(array_size), &
            c910(array_size), c4p(array_size), c8_lambda(array_size), c7_lambda(array_size), c6_lambda(array_size), c5_lambda(array_size), &
            c4_lambda(array_size), c3_lambda(array_size), c2_lambda(array_size), c910_lambda(array_size), c4p_lambda(array_size), &
            a1234(array_size), a1234_lambda(array_size), anion_integrand(array_size), cation_integrand(array_size))
    else
       print *, "constructoligomers.f90: calculate_C4MIMBF4_ideal_chain_term:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    a1234 = integrate_phi_spherical(exp(lambda_minus))
    a1234_lambda = integrate_phi_spherical(exp(lambda_minus) * lambda_minus)

    anion_integrand = ( 4.0_dp*(a1234**3) * a1234_lambda ) + (a1234**4)*(log(bulk_density) - 1.0_dp) + ((a1234**4)*lambda_minus)
    anion_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_minus) * anion_integrand, unity_function)

    c8 = integrate_phi_spherical(exp(lambda_neutral))
    c8_lambda = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8)
    c7_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c8_lambda + c8*lambda_neutral))

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)
    c6_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c7_lambda + c7*lambda_neutral))

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)
    c5_lambda = integrate_phi_spherical(exp(lambda_neutral) * (c6_lambda + c6*lambda_neutral))

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)
    c4_lambda = integrate_phi_spherical(exp(lambda_plus) * (c5_lambda + c5*lambda_plus))

    c910 = integrate_phi_spherical(exp(lambda_plus))
    c910_lambda = integrate_phi_spherical(exp(lambda_plus) * lambda_plus)

    c4p = c4*(c910**2)
    c4p_lambda = c4_lambda*(c910**2) + 2.0_dp * (c4*c910_lambda*c910)

    c3 = integrate_phi_spherical(exp(lambda_plus) * c4p)
    c3_lambda = integrate_phi_spherical(exp(lambda_plus) * (c4p_lambda + c4p*lambda_plus))

    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c2_lambda = integrate_phi_spherical(exp(lambda_plus) * (c3_lambda + c3*lambda_plus))

    cation_integrand = c2_lambda + c2*(log(bulk_density) - 1.0_dp) + c2*lambda_neutral
    cation_contribution = (bulk_density/beta) * integrate_z_cylindrical(exp(lambda_neutral) * cation_integrand, unity_function)

    calculate_C4MIMBF4_ideal_chain_term =  cation_contribution + anion_contribution

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(c8, c7, c6, c5, c4, c3, c2, c910, c4p, c8_lambda, c7_lambda, c6_lambda, c5_lambda, c4_lambda, &
         c3_lambda, c2_lambda, c910_lambda, c4p_lambda, a1234, a1234_lambda, anion_integrand, cation_integrand)

  end function calculate_C4MIMBF4_ideal_chain_term


  subroutine UpdateSinglePositiveNeutralMinusSphereDensities(lambda_plus, n_plus_updated, &
       lambda_neutral, n_neutral_updated, lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated

    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: UpdateSinglePositiveNeutralMinusSphereDensities:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    call UpdateSinglePositiveSphereDensity(lambda_plus, n_plus_updated)
    call UpdateSingleNeutralSphereDensity(lambda_neutral, n_neutral_updated)
    call UpdateSingleNegativeSphereDensity(lambda_minus, n_minus_updated)

  end subroutine UpdateSinglePositiveNeutralMinusSphereDensities

  subroutine UpdateC4MIMPositiveBeadDensities(lambda_plus, lambda_neutral, n_plus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_plus_updated

    real(dp), dimension(:), allocatable :: c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp
    integer :: array_size

    ! First check the input variables are the same size
    if( (size(lambda_plus) == size(lambda_neutral)) .and. (size(lambda_plus) == size(n_plus_updated))) then

       array_size = size(lambda_plus)
       allocate(c8c1(array_size), c9c10(array_size), c7(array_size), c6(array_size), &
            c5(array_size), c4(array_size), c2(array_size), c3p(array_size), c3pp(array_size), c3ppp(array_size))
    else
       print *, "constructoligomers.f90: UpdateC4MIMPositiveBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_plus_updated)"
       print *, "...aborting..."
       call abort()
    end if

    ! Note that there is a factor of 1.0_dp/(4.0_dp * pi * hs_diameter**2)
    ! from the delta function bond, a factor of 2 * pi from the theta integral
    ! and a factor of hs_diameter**2 from the jacobian.  Combining these factors
    ! gives the required 0.5_dp factor.  This is taken care of inside the integrate_phi_spherical routine.
    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c8c1 = integrate_phi_spherical(exp(lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8c1)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c8c1)

    c3pp = integrate_phi_spherical(exp(lambda_plus) * c4 * c2 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    !Calculate the resulting positive bead densities.
    !n_plus_updated = nc2 + nc3 + nc4 + nc9 + nc10 =  nc2 + nc3 + nc4 + 2*nc9
    n_plus_updated = bulk_density * ( (exp(lambda_plus) * c8c1 * c3p) + (exp(lambda_plus) * c2 * c9c10 * c9c10 * c4) + &
         (exp(lambda_plus) * c3ppp * c5) + (2.0_dp * (exp(lambda_plus) * c3pp)) )

    call setNonCalculatedRegionToZero(n_plus_updated)

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp)

  end subroutine UpdateC4MIMPositiveBeadDensities

  subroutine UpdateC4MIMNeutralBeadDensities(lambda_plus, lambda_neutral, n_neutral_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    real(dp), dimension(:), allocatable :: c3p, c3ppp, c9c10, c8c1, c7, c6, c5, c4, c2 !contributions also used in +ve beads.
    real(dp), dimension(:), allocatable :: c2p, c4p, c5p, c6p, c7p !extra contributions for -ve beads.

    integer :: array_size

    ! First check the input variables are the same size
    if( (size(lambda_plus) == size(lambda_neutral)) .and. (size(lambda_plus) == size(n_neutral_updated))) then

       array_size = size(lambda_plus)
       allocate(c3p(array_size), c3ppp(array_size), c8c1(array_size), c9c10(array_size), c2p(array_size), &
            c7(array_size), c6(array_size), c5(array_size), c4(array_size), c2(array_size), &
            c4p(array_size), c5p(array_size), c6p(array_size), c7p(array_size))
    else
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    c9c10 = integrate_phi_spherical(exp(lambda_plus))

    c8c1 = integrate_phi_spherical(exp(lambda_neutral))

    c7 = integrate_phi_spherical(exp(lambda_neutral) * c8c1)

    c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = integrate_phi_spherical(exp(lambda_plus) * c8c1)

    c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

    !Calculate the resulting neutral bead densities.
    !n_zero_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8c1) + (exp(lambda_neutral) * c7p) )

    call setNonCalculatedRegionToZero(n_neutral_updated)

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(c3p, c3ppp, c8c1)
    deallocate(c2p, c4p, c5p, c6p, c7p)

  end subroutine UpdateC4MIMNeutralBeadDensities

  subroutine UpdateBF4NegativeBeadDensities(lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    real(dp), dimension(:), allocatable :: a1a2a3a4, a5p
    integer :: array_size

    ! First check the input variables are the same size
    if(size(lambda_minus) == size(n_minus_updated)) then
       array_size = size(lambda_minus)
       allocate(a1a2a3a4(array_size))
       allocate(a5p(array_size))

       a1a2a3a4(:) = 0.0_dp
       a5p(:) = 0.0_dp
    else
       print *, "constructoligomers.f90: UpdateBF4NegativeBeads:"
       print *, "Size mismatch. size(lambda_minus) /= size(n_minus_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !Calculate the required contributions for the anion
    a1a2a3a4 = integrate_phi_spherical(exp(lambda_minus))

    a5p = integrate_phi_spherical(exp(lambda_minus) * (a1a2a3a4 ** 3.0_dp))

    !Calculate the resulting negative bead densities.
    !n_minus_updated = na1 + na2 + na3 + na4 + na5 = 4*na1 + na5
    n_minus_updated = bulk_density * ( 4.0_dp*(exp(lambda_minus) * a5p) + (exp(lambda_minus) * (a1a2a3a4**4.0_dp)) )

    call setNonCalculatedRegionToZero(n_minus_updated)

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(a1a2a3a4, a5p)

  end subroutine UpdateBF4NegativeBeadDensities

end module constructoligomers


!Contains routines neccesary to compare results with the contact theorem
module contacttheorem
  use kinds
  use helpers
  use parameters
  use discretederivatives
  use integratezcylindrical
  use functionalderivatives
  implicit none
  private

  public :: CalculateNormalPressureFromContactTheorem
  public :: InitialisePotentialAndContactTheoremVariables
  public :: CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation
  public :: MakeContactTheoremAdjustmentFromParticleParticleDispersion

contains

  subroutine CalculateNormalPressureFromContactTheorem(n_plus, n_neutral, n_minus, &
       normal_pressure_left_wall, normal_pressure_right_wall, dispersion_particle_particle_adjust_to_contact_thm)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    real(dp), intent(out) :: normal_pressure_left_wall
    real(dp), intent(out) :: normal_pressure_right_wall
    real(dp), intent(out) :: dispersion_particle_particle_adjust_to_contact_thm

    real(dp), dimension(size(n_plus)) :: n_s

    real(dp), dimension(size(n_plus)) :: left_wall_dispersion_integrand
    real(dp), dimension(size(n_plus)) :: right_wall_dispersion_integrand

    real(dp) :: maxwell_stress_term_left_wall
    real(dp) :: maxwell_stress_term_right_wall

    integer :: start_z_index
    integer :: end_z_index

    !First check the sizes are the same
    if( (size(n_plus) /= size(n_neutral)) .or. (size(n_plus) /= size(n_minus)) ) then
       print *, "contacttheorem.f90:CalculateNormalPressureFromContactTheorem: "
       print *, "(size(n_plus) /= size(n_neutral)) .or. (size(n_plus) /= size(n_minus))"
       print *, "Density array size mismatch...aborting..."
       call abort()
    end if

    left_wall_dispersion_integrand = 0.0_dp
    right_wall_dispersion_integrand = 0.0_dp

    maxwell_stress_term_left_wall = 0.0_dp
    maxwell_stress_term_right_wall = 0.0_dp

    n_s = n_plus + n_neutral + n_minus
    call get_allowed_z_values(start_z_index, end_z_index, size(n_s))

    call CalculateDerivOfWallTerm(left_wall_dispersion_integrand, right_wall_dispersion_integrand)

    call CalculateDispersionAdjustment(dispersion_particle_particle_adjust_to_contact_thm, n_s)

    call CalculateMaxwellStressTerm(maxwell_stress_term_left_wall, maxwell_stress_term_right_wall)

    !left_wall_dispersion_integrand = 0.0_dp
    !right_wall_dispersion_integrand = 0.0_dp

    !maxwell_stress_term_left_wall = 0.0_dp
    !maxwell_stress_term_right_wall = 0.0_dp


    normal_pressure_left_wall =  (n_s(start_z_index) + integrate_z_cylindrical(n_s * left_wall_dispersion_integrand, unity_function) - maxwell_stress_term_left_wall) / beta
    normal_pressure_right_wall =  (n_s(end_z_index) + integrate_z_cylindrical(n_s * right_wall_dispersion_integrand, unity_function) - maxwell_stress_term_right_wall) / beta

  end subroutine CalculateNormalPressureFromContactTheorem

  subroutine InitialisePotentialAndContactTheoremVariables(grand_potential_per_unit_area, grand_potential_per_unit_area_in_bulk, &
       normal_pressure_left_wall, normal_pressure_right_wall, negative_deriv_of_potential, dispersion_particle_particle_adjust_to_contact_thm)
    real(dp), dimension(:), allocatable :: grand_potential_per_unit_area
    real(dp), dimension(:), allocatable :: grand_potential_per_unit_area_in_bulk
    real(dp), dimension(:), allocatable :: normal_pressure_left_wall
    real(dp), dimension(:), allocatable :: normal_pressure_right_wall
    real(dp), dimension(:), allocatable :: negative_deriv_of_potential
    real(dp), dimension(:), allocatable :: dispersion_particle_particle_adjust_to_contact_thm

    allocate(grand_potential_per_unit_area(size(plate_separations)))
    allocate(grand_potential_per_unit_area_in_bulk(size(plate_separations)))
    allocate(normal_pressure_left_wall(size(plate_separations)))
    allocate(normal_pressure_right_wall(size(plate_separations)))
    allocate(negative_deriv_of_potential(size(plate_separations)))
    allocate(dispersion_particle_particle_adjust_to_contact_thm(size(plate_separations)))

    grand_potential_per_unit_area(:) = 0.0_dp
    grand_potential_per_unit_area_in_bulk(:) = 0.0_dp
    normal_pressure_left_wall(:) = 0.0_dp
    normal_pressure_right_wall(:) = 0.0_dp
    negative_deriv_of_potential(:) = 0.0_dp
    dispersion_particle_particle_adjust_to_contact_thm(:) = 0.0_dp

  end subroutine InitialisePotentialAndContactTheoremVariables

  function CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation(grand_potential_per_unit_area)
    real(dp), dimension(:), intent(in) :: grand_potential_per_unit_area
    real(dp), dimension(size(grand_potential_per_unit_area)) :: CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation

    CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation = -1.0_dp * calculate_central_difference(grand_potential_per_unit_area)

  end function CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation

  subroutine CalculateMaxwellStressTerm(maxwell_stress_term_left_wall, maxwell_stress_term_right_wall)
    real(dp), intent(out) :: maxwell_stress_term_left_wall
    real(dp), intent(out) :: maxwell_stress_term_right_wall

    maxwell_stress_term_left_wall = 2.0_dp * pi * (surface_charge_density_left_wall**2)/ (epsilonr * epsilon0) / beta
    maxwell_stress_term_right_wall = 2.0_dp * pi * (surface_charge_density_right_wall**2)/ (epsilonr * epsilon0) / beta

  end subroutine CalculateMaxwellStressTerm

  subroutine CalculateDerivOfWallTerm(left_wall_term, right_wall_term)
    real(dp), dimension(:) :: left_wall_term
    real(dp), dimension(:) :: right_wall_term

    integer :: ij
    real(dp) :: distance_from_left_wall
    real(dp) :: distance_from_right_wall

    do ij = 2, size(left_wall_term) - 1
       
       distance_from_left_wall = (ij - 1)*hs_diameter/real(n_discretised_points_z,dp)
       distance_from_right_wall = (size(left_wall_term) - ij)*hs_diameter/real(n_discretised_points_z,dp)
       
       left_wall_term(ij) = (beta) * (2.0_dp * pi * epsilon_LJ * ((0.4_dp * (hs_diameter/distance_from_left_wall)**9) - (hs_diameter/distance_from_left_wall)**3))/distance_from_left_wall
       right_wall_term(ij) = (beta) * (2.0_dp * pi * epsilon_LJ * ((0.4_dp * (hs_diameter/distance_from_right_wall)**9) - (hs_diameter/distance_from_right_wall)**3))/distance_from_right_wall
    end do


  end subroutine CalculateDerivOfWallTerm

  subroutine CalculateDispersionAdjustment(dispersion_particle_particle_adjust_to_contact_thm, n_s)
    real(dp), intent(out) :: dispersion_particle_particle_adjust_to_contact_thm
    real(dp), dimension(:), intent(in) :: n_s

    dispersion_particle_particle_adjust_to_contact_thm = integrate_z_cylindrical(0.5_dp * n_s * calculate_vanderWaals_functional_deriv(n_s), unity_function)
    
  end subroutine CalculateDispersionAdjustment

  subroutine MakeContactTheoremAdjustmentFromParticleParticleDispersion(normal_pressure_left_wall, normal_pressure_right_wall, dispersion_particle_particle_adjust_to_contact_thm)
    real(dp), dimension(:) :: normal_pressure_left_wall
    real(dp), dimension(:) :: normal_pressure_right_wall
    real(dp), dimension(:) :: dispersion_particle_particle_adjust_to_contact_thm

    integer :: ij

    dispersion_particle_particle_adjust_to_contact_thm(:) = CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation(dispersion_particle_particle_adjust_to_contact_thm)
    
    do ij = 1, size(normal_pressure_left_wall)
       normal_pressure_left_wall(ij) =  normal_pressure_left_wall(ij) + dispersion_particle_particle_adjust_to_contact_thm(ij)
       normal_pressure_right_wall(ij) =  normal_pressure_right_wall(ij) + dispersion_particle_particle_adjust_to_contact_thm(ij)
    end do

  end subroutine MakeContactTheoremAdjustmentFromParticleParticleDispersion

end module contacttheorem


!Contains routines that perform the discrete derivative
module discretederivatives
  use kinds
  use parameters
  implicit none
  private

  public :: calculate_central_difference
  public :: calculate_backward_difference
  public :: calculate_forward_difference

contains

  !Calculate the central difference, unless we are at the beginning
  !or end of the walls in which case we take the appropriate forward
  !or backward difference.
  function calculate_central_difference(input_array)
    real(dp), dimension(:), intent(in) :: input_array
    real(dp), dimension(size(input_array)) :: calculate_central_difference

    integer :: ij

    calculate_central_difference = 0.0_dp

    !Check input_array is the correct size
    if(size(input_array) /= size(plate_separations)) then
       print *, "discretederivatives.f90: calculate_central_difference: "
       print *, "size(input_array) /= size(plate_separations)"
       print *, "size mismatch....likely coding error...aborting..."
       call abort()
    end if

    !If the input_array size is one we can't calculate any derivatives
    !but we may still want the density profile. So print a warning and set the derivative to 0.
    if(size(input_array) == 1) then
       print *, "discretederivatives.f90: calculate_central_difference: "
       print *, "trying to calculate a discrete derivative with only one point"
       print *, "not possible.  Setting the derivative to 0 and returning as we"
       print *, "may want the density profile."
       calculate_central_difference = 0.0_dp
       return
    end if

    !Calculate the central difference whenever we can
    do ij = 2, size(input_array) - 1
       calculate_central_difference(ij) = ( input_array(ij+1) - input_array(ij-1) ) / &
            ( abs(plate_separations(ij+1) - plate_separations(ij-1)) * hs_diameter)
    end do

    !If we are at the first array point, we can't do a central difference,
    !therefore, we calculate the forward difference.
    calculate_central_difference(1) = ( input_array(2) - input_array(1) ) / &
         ( abs(plate_separations(2) - plate_separations(1)) * hs_diameter)

    !If we are at the last array point, we can't do a central difference,
    !therefore, we calculate the backward difference.
    calculate_central_difference(size(input_array)) = &
         ( input_array(size(input_array)) - input_array(size(input_array) - 1) ) / &
         ( abs(plate_separations(size(input_array)) - plate_separations(size(input_array) - 1)) * hs_diameter)

  end function calculate_central_difference

  !Calculate the backward difference unless we are at the start point in
  !which case we calculate the forward difference.
  function calculate_backward_difference(input_array)
    real(dp), dimension(:), intent(in) :: input_array
    real(dp), dimension(size(input_array)) :: calculate_backward_difference

    integer :: ij

    !Check input_array is the correct size
    if(size(input_array) /= size(plate_separations)) then
       print *, "discretederivatives.f90: calculate_backward_difference: "
       print *, "size(input_array) /= size(plate_separations)"
       print *, "size mismatch....likely coding error...aborting..."
       call abort()
    end if

    !If the input_array size is one we can't calculate any derivatives
    !but we may still want the density profile. So print a warning and set the derivative to 0.
    if(size(input_array) == 1) then
       print *, "discretederivatives.f90: calculate_central_difference: "
       print *, "trying to calculate a discrete derivative with only one point"
       print *, "not possible.  Setting the derivative to 0 and returning as we"
       print *, "may want the density profile."
       calculate_backward_difference = 0.0_dp
       return
    end if

    calculate_backward_difference = 0.0_dp

    !Calculat the backward difference whenever we can.
    do ij = 2, size(input_array)
       calculate_backward_difference(ij) = ( input_array(ij) - input_array(ij-1) ) / &
            ( abs(plate_separations(ij) - plate_separations(ij-1)) * hs_diameter)
    end do

    !On the first point, we can't calculate the backward difference so calculate the forward difference.
    calculate_backward_difference(1) = ( input_array(2) - input_array(1) ) / &
         ( abs(plate_separations(2) - plate_separations(1)) * hs_diameter)

  end function calculate_backward_difference

  !Calculate the forward difference unless we are at the endpoint
  !in which case we calculate the backward difference.
  function calculate_forward_difference(input_array)
    real(dp), dimension(:), intent(in) :: input_array
    real(dp), dimension(size(input_array)) :: calculate_forward_difference

    integer :: ij

    !Check input_array is the correct size
    if(size(input_array) /= size(plate_separations)) then
       print *, "discretederivatives.f90: calculate_forward_difference: "
       print *, "size(input_array) /= size(plate_separations)"
       print *, "size mismatch....likely coding error...aborting..."
       call abort()
    end if

    !If the input_array size is one we can't calculate any derivatives
    !but we may still want the density profile. So print a warning and set the derivative to 0.
    if(size(input_array) == 1) then
       print *, "discretederivatives.f90: calculate_central_difference: "
       print *, "trying to calculate a discrete derivative with only one point"
       print *, "not possible.  Setting the derivative to 0 and returning as we"
       print *, "may want the density profile."
       calculate_forward_difference = 0.0_dp
       return
    end if

    calculate_forward_difference = 0.0_dp

    !Calculate the forward difference whenever we can
    do ij = 1, size(input_array) - 1
       calculate_forward_difference(ij) = ( input_array(ij + 1) - input_array(ij) ) / &
            ( abs(plate_separations(ij+1) - plate_separations(ij)) * hs_diameter)
    end do

    !On the end point, we can't do the forward difference, therefore do the backward difference.
    calculate_forward_difference(size(input_array)) = &
         ( input_array(size(input_array)) - input_array(size(input_array) - 1) ) / &
         ( abs(plate_separations(size(input_array)) - plate_separations(size(input_array) - 1)) * hs_diameter)

  end function calculate_forward_difference

end module discretederivatives

!The notation used in this module follows that of the paper
!"Evaluating the accuracy of a density functional theory of polymer
!solutions with additive hard sphere diameters", by Forsman and Woodward
!Journal of Chemical Physics, volume 120, number 1, Jan 2004.
module excessenergyfunctionalparameters
  use kinds
  use universalconstants
  use parameters
  implicit none
  private

  public :: GetAEx
  public :: GetAExDerivIntegrand

  real(dp) :: sigma_monomer
  real(dp) :: sigma_solvent

contains

  subroutine GetAExTermWeightings(n_mbar, n_sbar, term_index, Y_special, Psi, X, W, Z)
    real(dp), dimension(:), intent(in)              :: n_mbar
    real(dp), dimension(:), intent(in)              :: n_sbar
    integer, intent(in)                             :: term_index !0,1,2 as discussed in J.Chem.Phys. vol 120, #1, 1994.
    real(dp), dimension(size(n_mbar)), intent(out)  :: Y_special
    real(dp), dimension(size(n_mbar)), intent(out)  :: Psi
    real(dp), dimension(size(n_mbar)), intent(out)  :: W
    real(dp), dimension(size(n_mbar)), intent(out)  :: X
    real(dp), dimension(size(n_mbar)), intent(out)  :: Z

    real(dp) :: q

    real(dp) :: r_0, r_1, r_2
    real(dp) :: s_0, s_1, s_2
    real(dp) :: b_0, b_1, b_2

    real(dp), dimension(size(n_mbar)) :: x_d_hat
    real(dp), dimension(size(n_mbar)) :: x_m_hat

    real(dp), dimension(size(n_mbar)) :: rhat_2, rhat_1, rhat_0
    real(dp), dimension(size(n_mbar)) :: shat_2, shat_1, shat_0 
    real(dp), dimension(size(n_mbar)) :: bhat_2, bhat_1, bhat_0

    real(dp), dimension(size(n_mbar)) :: c_0, c_1, c_2

    real(dp), dimension(size(n_mbar)) :: omega_0, omega_1, omega_2
    real(dp), dimension(size(n_mbar)) :: phi_0, phi_1, phi_2
    real(dp), dimension(size(n_mbar)) :: omega_prime_0, omega_prime_1, omega_prime_2 
    real(dp), dimension(size(n_mbar)) :: W_0, W_1, W_2
    real(dp), dimension(size(n_mbar)) :: X_0, X_1, X_2
    real(dp), dimension(size(n_mbar)) :: Y_0, Y_1, Y_2
    real(dp), dimension(size(n_mbar)) :: Z_0, Z_1, Z_2
    real(dp), dimension(size(n_mbar)) :: Y_special_0, Y_special_1, Y_special_2
    real(dp), dimension(size(n_mbar)) :: Psi_0, Psi_1, Psi_2

    real(dp), dimension(size(n_mbar)) :: eta

    if(size(n_mbar) /= size(n_sbar)) then
       print *, "excessenergyfunctionalparameters.f90: GetAExTermWeightings"
       print *, "size(n_mbar) /= size(n_sbar)"
       print *, "input variable size mismatch...aborting..."
       call abort()
    end if

    call InitialiseHardSphereDiameters(sigma_monomer, sigma_solvent)
    
    q = sigma_solvent / sigma_monomer

    eta(:) = pi * (n_mbar(:) + n_sbar(:)) * (hs_diameter**3) / 6.0_dp

    r_2 = 0.75_dp * sigma_monomer
    r_1 = 0.5_dp * sigma_monomer
    r_0 = 0.5_dp * sigma_monomer * q

    s_2 = 2.0_dp * pi * (sigma_monomer**2)
    s_1 = pi * (sigma_monomer**2)
    s_0 = pi * (sigma_monomer**2) * (q**2)

    b_2 = (pi/3.0_dp) * (sigma_monomer**3)
    b_1 = (pi/6.0_dp) * (sigma_monomer**3)
    b_0 = (pi/6.0_dp) * (sigma_monomer**3) * (q**3)

    x_d_hat(:) = (0.5_dp * n_mbar(:)) / ((0.5_dp * n_mbar(:)) + n_sbar(:))
    x_m_hat(:) = n_mbar(:) / (n_mbar(:) + n_sbar(:))

    rhat_2(:) = (x_d_hat(:) * r_2) + ((1.0_dp - x_d_hat(:)) * r_0)
    rhat_1(:) = (x_m_hat(:) * r_1) + ((1.0_dp - x_m_hat(:)) * r_0)
    rhat_0(:) = rhat_2(:)

    shat_2(:) = (x_d_hat(:) * s_2) + ((1.0_dp - x_d_hat(:)) * s_0)
    shat_1(:) = (x_m_hat(:) * s_1) + ((1.0_dp - x_m_hat(:)) * s_0)
    shat_0(:) = shat_2(:)

    bhat_2(:) = (x_d_hat(:) * b_2) + ((1.0_dp - x_d_hat(:)) * b_0)
    bhat_1(:) = (x_m_hat(:) * b_1) + ((1.0_dp - x_m_hat(:)) * b_0)
    bhat_0(:) = bhat_2(:)

    c_2(:) = (x_d_hat(:) * (r_2**2)) + ((1.0_dp - x_d_hat(:)) * (r_0**2))
    c_1(:) = (x_m_hat(:) * (r_1**2)) + ((1.0_dp - x_m_hat(:)) * (r_0**2))
    c_0(:) = c_2(:)

    omega_2(:) = (rhat_2(:) * shat_2(:)) / bhat_2(:)
    omega_1(:) = (rhat_1(:) * shat_1(:)) / bhat_1(:)
    omega_0(:) = (rhat_0(:) * shat_0(:)) / bhat_0(:)

    phi_2(:) = (c_2(:) * (shat_2(:)**2)) / (9.0_dp * (bhat_2(:)**2))
    phi_1(:) = (c_1(:) * (shat_1(:)**2)) / (9.0_dp * (bhat_1(:)**2))
    phi_0(:) = (c_0(:) * (shat_0(:)**2)) / (9.0_dp * (bhat_0(:)**2))

    omega_prime_2(:) = ((r_2 * shat_2(:)) + (rhat_2(:) * shat_2(:)) - (b_2 * (1.0_dp / bhat_2(:)) * rhat_2(:) * shat_2(:))) &
         / (bhat_2(:))
    omega_prime_1(:) = ((r_1 * shat_1(:)) + (rhat_1(:) * shat_1(:)) - (b_1 * (1.0_dp / bhat_1(:)) * rhat_1(:) * shat_1(:))) &
         / (bhat_1(:))
    omega_prime_0(:) = ((r_0 * shat_0(:)) + (rhat_0(:) * shat_0(:)) - (b_0 * (1.0_dp / bhat_0(:)) * rhat_0(:) * shat_0(:))) &
         / (bhat_0(:))

    W_2(:) = ((2.0_dp * b_2 * (1.0_dp / bhat_2(:)) * c_2(:) * (shat_2(:)**2)) - (2.0_dp * c_2 * shat_2(:) * s_2) - &
         ((rhat_2(:)**2) * (shat_2(:)**2))) / (9.0_dp * (bhat_2(:)**2))
    W_1(:) = ((2.0_dp * b_1 * (1.0_dp / bhat_1(:)) * c_1(:) * (shat_1(:)**2)) - (2.0_dp * c_1 * shat_1(:) * s_1) - &
         ((rhat_1(:)**2) * (shat_1(:)**2))) / (9.0_dp * (bhat_1(:)**2))
    W_0(:) = ((2.0_dp * b_0 * (1.0_dp / bhat_0(:)) * c_0(:) * (shat_0(:)**2)) - (2.0_dp * c_0 * shat_0(:) * s_0) - &
         ((rhat_0(:)**2) * (shat_0(:)**2))) / (9.0_dp * (bhat_0(:)**2))

    X_2(:) = W_2(:) - omega_prime_2(:) - ( (omega_2(:) + 1.0_dp) * (b_2) / (bhat_2(:)) )
    X_1(:) = W_1(:) - omega_prime_1(:) - ( (omega_1(:) + 1.0_dp) * (b_1) / (bhat_1(:)) )
    X_0(:) = W_0(:) - omega_prime_0(:) - ( (omega_0(:) + 1.0_dp) * (b_0) / (bhat_0(:)) )

    Y_2(:) = (-1.0_dp * W_2(:)) + (2.0_dp * omega_prime_2(:)) - (3.0_dp * (phi_2(:) - omega_2(:) - 2.0_dp) * b_2 / bhat_2(:))
    Y_1(:) = (-1.0_dp * W_1(:)) + (2.0_dp * omega_prime_1(:)) - (3.0_dp * (phi_1(:) - omega_1(:) - 2.0_dp) * b_1 / bhat_1(:))
    Y_0(:) = (-1.0_dp * W_0(:)) + (2.0_dp * omega_prime_0(:)) - (3.0_dp * (phi_0(:) - omega_0(:) - 2.0_dp) * b_0 / bhat_0(:))

    Z_2(:) = ( (phi_2(:) - 1.0_dp) * b_2 / bhat_2(:) ) - omega_prime_2(:)
    Z_1(:) = ( (phi_1(:) - 1.0_dp) * b_1 / bhat_1(:) ) - omega_prime_1(:)
    Z_0(:) = ( (phi_0(:) - 1.0_dp) * b_0 / bhat_0(:) ) - omega_prime_0(:)

    Y_special_2(:) = W_2(:) + Y_2(:) + ( 3.0_dp * Z_2(:) )
    Y_special_1(:) = W_1(:) + Y_1(:) + ( 3.0_dp * Z_1(:) )
    Y_special_0(:) = W_0(:) + Y_0(:) + ( 3.0_dp * Z_0(:) )

    Psi_2(:) = ( 0.5_dp * X_2(:) ) + ( 2.5_dp * Y_2(:) ) + ( 5.5_dp * Z_2(:) ) + ( 3.0_dp * W_2(:) )
    Psi_1(:) = ( 0.5_dp * X_1(:) ) + ( 2.5_dp * Y_1(:) ) + ( 5.5_dp * Z_1(:) ) + ( 3.0_dp * W_1(:) )
    Psi_0(:) = ( 0.5_dp * X_0(:) ) + ( 2.5_dp * Y_0(:) ) + ( 5.5_dp * Z_0(:) ) + ( 3.0_dp * W_0(:) )  

    !For ease of readability and comparison with the paper we first calculate all the values of the
    !term weightings and then return only the ones that we have asked for.  If calculation speed becomes
    !a problem we can change this.
    if(term_index == 0) then

       Y_special(:) = Y_special_0(:)
       Psi(:) = Psi_0(:)
       Z(:) = Z_0(:)
       X(:) = X_0(:)
       W(:) = W_0(:)

    else if(term_index == 1) then

       Y_special(:) = Y_special_1(:)
       Psi(:) = Psi_1(:)
       Z(:) = Z_1(:)
       X(:) = X_1(:)
       W(:) = W_1(:)

    else if(term_index == 2) then

       Y_special(:) = Y_special_2(:)
       Psi(:) = Psi_2(:)
       Z(:) = Z_2(:)
       X(:) = X_2(:)
       W(:) = W_2(:)

    else
       print *, "excessenergyfunctionalparameters.f90: GetAExTermWeightings"
       print *, "Value of term index = ", term_index
       print *, "Only supports values of 0, 1 or 2 as per the paper"
       print *, "Coding error...aborting..."
    end if

  end subroutine GetAExTermWeightings

  function GetAEx(n_mbar, n_sbar, term_index)
    real(dp), dimension(:), intent(in)  :: n_mbar
    real(dp), dimension(:), intent(in)  :: n_sbar
    integer, intent(in)                 :: term_index
    real(dp), dimension(size(n_mbar)) :: GetAEx

    real(dp), dimension(size(n_mbar)) :: Y_special, Psi, X, W, Z
    real(dp), dimension(size(n_mbar)) :: eta_bar

    if(size(n_mbar) /= size(n_sbar)) then
       print *, "excessenergyfunctionalparameters.f90: GetAEx"
       print *, "size(n_mbar) /= size(n_sbar)"
       print *, "input variable size mismatch...aborting..."
       call abort()
    end if

    call GetAExTermWeightings(n_mbar, n_sbar, term_index, Y_special, Psi, X, W, Z)

    eta_bar(:) = (pi * (n_mbar(:) + n_sbar(:)) * hs_diameter**3) / 6.0_dp

    GetAEx(:) = (( (1.0_dp - eta_bar(:)) / eta_bar(:) ) * log(1.0_dp - eta_bar(:)) ) + 1.0_dp + &
         ( Y_special(:) * ( (log(1.0_dp - eta_bar(:)) * ( ((1.0_dp - eta_bar(:))/(eta_bar(:))) + 1.0_dp )) - (1.0_dp/(1.0_dp - eta_bar(:))) -&
         (1.0_dp/(2.0_dp * ((1.0_dp - eta_bar(:))**2))) + (5.0_dp/2.0_dp) ) ) + &
         ( (Psi(:) + (2.0_dp * Z(:)) - X(:)) * (((1.0_dp)/(2.0_dp * ((1.0_dp - eta_bar(:))**2))) - 0.5_dp) ) + &
         (W(:) * ((1.0_dp/(2.0_dp * ((1.0_dp - eta_bar(:))**2))) - (2.0_dp/(1.0_dp - eta_bar(:))) - log(1.0_dp - eta_bar(:)) + 1.5_dp )) +&
         (Psi(:) * ((1.0_dp / (1.0_dp - eta_bar(:))) - (1.0_dp / (2.0_dp * ((1.0_dp - eta_bar(:))**2))) - 0.5_dp) )

  end function GetAEx

  function GetAExDerivIntegrand(n_mbar, n_sbar, term_index)
    real(dp), dimension(:), intent(in) :: n_mbar
    real(dp), dimension(:), intent(in) :: n_sbar
    integer, intent(in)                :: term_index
    real(dp), dimension(size(n_mbar)) :: GetAExDerivIntegrand

    real(dp), dimension(size(n_mbar)) :: eta_bar

    real(dp), dimension(size(n_mbar)) :: Y_special, Psi, X, W, Z

    if(size(n_mbar) /= size(n_sbar)) then
       print *, "excessenergyfunctionalparameters.f90: GetAEx"
       print *, "size(n_mbar) /= size(n_sbar)"
       print *, "input variable size mismatch...aborting..."
       call abort()
    end if

    call GetAExTermWeightings(n_mbar, n_sbar, term_index, Y_special, Psi, X, W, Z)

    eta_bar(:) = (pi * (n_mbar(:) + n_sbar(:)) * hs_diameter**3) / 6.0_dp

    GetAExDerivIntegrand(:) = ( (-1.0_dp / 8.0_dp) * ((log(1.0_dp - eta_bar(:))/(eta_bar(:)**2)) + (1.0_dp/eta_bar(:))) ) + & !"T1"
         ( ((-Y_special(:))/8.0_dp) * ((1.0_dp/eta_bar(:)) + (1.0_dp/(1.0_dp - eta_bar(:))) + (1.0_dp/((1.0_dp - eta_bar(:))**2)) +&
         (1.0_dp/((1.0_dp - eta_bar(:))**3)) + (log(1.0_dp - eta_bar(:))/eta_bar(:)**2)) ) + & !"T2"
         (((Psi(:) + (2.0_dp * Z(:)) - X(:))/8.0_dp) * (1.0_dp / ((1.0_dp - eta_bar(:))**3))) +& !"T3"
         ((W(:)/8.0_dp) * ((1.0_dp/((1.0_dp - eta_bar(:))**3)) - (2.0_dp/((1.0_dp - eta_bar(:))**2)) + &
         (1.0_dp/(1.0_dp - eta_bar(:))))) +& !"T4"
         ((Psi(:)/8.0_dp) * ((1.0_dp/((1.0_dp - eta_bar(:))**2)) - (1.0_dp/((1.0_dp - eta_bar(:))**3)))) !"T5"

  end function GetAExDerivIntegrand

  function GetYMix(n_mbar, n_sbar, alpha)
    real(dp), dimension(:), intent(in) :: n_mbar
    real(dp), dimension(:), intent(in) :: n_sbar
    character(len=1), intent(in) :: alpha
    real(dp), dimension(size(n_mbar)) :: GetYMix

    real(dp), dimension(size(n_mbar)) :: phi_m

    real(dp) :: q
    real(dp) :: sigma_alpha

    call InitialiseHardSphereDiameters(sigma_monomer, sigma_solvent)

    ! if(trim(alpha) == 'm') then
    !    sigma_alpha = sigma_monomer
    ! else if(trim(alpha) == 's') then
    !    sigma_alpha = sigma_solvent
    ! else
    !    print *, "excessenergyfunctionalparameters.f90: GetYMix:"
    !    print *, "alpha input variable must be either 'm' or 's'"
    !    print *, "instead, trim(alpha) = ", trim(alpha), "...aborting..."
    !    call abort()
    ! end if

    ! q = sigma_solvent / sigma_monomer

    ! phi_m(:) = n_mbar / (n_mbar + n_sbar*(q**3))

    ! v_1_alpha = pi * ((1.0_dp + (sigma_alpha/sigma_monomer))**3) / (6.0_dp * (sigma_monomer**3))

    ! v_2_alpha = (pi * (2.0_dp + (6.0_dp*(sigma_alpha/sigma_monomer)) + (4.5_dp*((sigma_alpha/sigma_monomer)**2)) + ((sigma_alpha/sigma_monomer)**3))) &
    !      / (6.0_dp * (sigma_monomer**3))

    ! v_3_alpha = (1.57_dp + (4.75_dp*(sigma_alpha/sigma_monomer)) + 2.99_dp*((sigma_alpha/sigma_monomer)**2) + 0.52_dp*((sigma_alpha/sigma_monomer)**3)  ) / (sigma_monomer**3)

    ! v_r_alpha = v_3_alpha + (r - 3.0_dp)*(v_3_alpha - v_2_alpha) - &
    !      0.04915_dp*((r - 3.0_dp)**1.09_dp)*((sigma_alpha/sigma_monomer)**2.71_dp)*(sigma_monomer**3)

    ! GetYMix(:) = (v_r_alpha - v_2_alpha) / (v_2_alpha - v_1_alpha)

  end function GetYMix

  subroutine InitialiseHardSphereDiameters(sigma_monomer, sigma_solvent)
    real(dp), intent(out) :: sigma_monomer
    real(dp), intent(out) :: sigma_solvent
    
    sigma_monomer = hs_diameter
    sigma_solvent = hs_diameter
    
  end subroutine InitialiseHardSphereDiameters
  
end module excessenergyfunctionalparameters

!Contains routines to calculate the functional derivative of each
!external (non-ideal) term in the grand potential.
module functionalderivatives
  use kinds
  use universalconstants
  use parameters
  use helpers
  use integratezcylindrical
  use io
  use excessenergyfunctionalparameters
  implicit none
  private

  public :: calculate_vanderWaals_functional_deriv
  public :: calculate_surface_dispersion_functional_deriv
  public :: calculate_hardsphere_functional_deriv
  public :: calculate_surface_electrostatic_functional_deriv

  public :: calculate_electrostatic_like_term_functional_deriv
  public :: calculate_electrostatic_unlike_term_functional_deriv

  public :: calculate_n_sbar

  real(dp) :: CURRENT_BULK_BEAD_DENSITY

contains

  function calculate_hardsphere_functional_deriv(n_s, calculate_bulk)
    real(dp), dimension(:), intent(in) :: n_s
    logical, intent(in)                :: calculate_bulk
    real(dp), dimension(size(n_s)) :: calculate_hardsphere_functional_deriv

    real(dp), dimension(size(n_s)) :: n_mbar, n_sbar
    real(dp), dimension(size(n_s)) :: integrand, integral

    integer :: start_z_index
    integer :: end_z_index

    integrand(:) = 0.0_dp
    integral(:) = 0.0_dp

    n_mbar(:) = 0.0_dp
    n_sbar(:) = 0.0_dp
    calculate_hardsphere_functional_deriv(:) = 0.0_dp

    n_mbar(:) = calculate_n_sbar(n_s(:))

    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    call get_allowed_z_values(start_z_index, end_z_index, size(n_s))

    if(calculate_bulk) then !n_s input parameter is the bulk value so

       calculate_hardsphere_functional_deriv(start_z_index:end_z_index) = (0.5_dp /beta) * (&
            GetAEx(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), a_term_index) + &
            (((4.0_dp * pi * (hs_diameter**3))/3.0_dp) * ((n_s(start_z_index:end_z_index)) * &
            GetAExDerivIntegrand(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), a_term_index))))

    else

       integrand(start_z_index:end_z_index) = n_s(start_z_index:end_z_index) * &
            GetAExDerivIntegrand(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), a_term_index)

       integral(:) = calculate_n_sbar(integrand(:))

       calculate_hardsphere_functional_deriv(start_z_index:end_z_index) = (0.5_dp / beta) * (&
            GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), a_term_index) + &
            (((4.0_dp * pi * (hs_diameter**3))/3.0_dp) * integral(start_z_index:end_z_index)))



    end if
    !calculate_hardsphere_functional_deriv = 0.0_dp

  end function calculate_hardsphere_functional_deriv

  function calculate_surface_dispersion_functional_deriv(ith_plate_separation, array_size)
    integer, intent(in) :: ith_plate_separation
    integer, intent(in) :: array_size
    real(dp), dimension(array_size) :: calculate_surface_dispersion_functional_deriv

    integer :: start_z_index
    integer :: end_z_index
    integer :: iz

    real(dp) :: hs_d_divide_z
    real(dp) :: hs_d_divide_h_minus_z

    calculate_surface_dispersion_functional_deriv(:) = 0.0_dp

    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    call get_allowed_z_values(start_z_index, end_z_index, array_size)

    !Note that we exclude the points that are calculate on the wall (at r_{z} = 0 or r_{z} = h)
    !as this leads to a singularity.
    do iz = start_z_index, end_z_index

       hs_d_divide_z = (real(n_discretised_points_z, dp) / real((iz - 1), dp))
       hs_d_divide_h_minus_z = (1.0_dp/(real(plate_separations(ith_plate_separation),dp) - &
            ( real((iz - 1),dp) / real((n_discretised_points_z),dp) )))

       calculate_surface_dispersion_functional_deriv(iz) = 2.0_dp * pi * epsilon_LJ * (1.0_dp) *(&
            ( (2.0_dp/45.0_dp)* (hs_d_divide_z**9.0_dp) ) - &
            ( (1.0_dp/3.0_dp) * (hs_d_divide_z**3.0_dp) ) + &
            ( (2.0_dp/45.0_dp)* (hs_d_divide_h_minus_z**9.0_dp) ) - &
            ( (1.0_dp/3.0_dp) * (hs_d_divide_h_minus_z**3.0_dp) ))
    end do

    call setNonCalculatedRegionToZero(calculate_surface_dispersion_functional_deriv)
   
    !calculate_surface_dispersion_functional_deriv = 0.0_dp

  end function calculate_surface_dispersion_functional_deriv

  function calculate_vanderWaals_functional_deriv(n_s)
    real(dp), dimension(:), intent(in) :: n_s
    real(dp), dimension(size(n_s)) :: calculate_vanderWaals_functional_deriv

    calculate_vanderWaals_functional_deriv = -4.0_dp * epsilon_LJ * (hs_diameter**6.0_dp) * (1.0_dp/beta) * &
         2.0_dp * pi * integrate_z_cylindrical(n_s, van_der_waals_density_indept_integrand, "all_z")

    call setNonCalculatedRegionToZero(calculate_vanderWaals_functional_deriv)
    calculate_vanderWaals_functional_deriv = 0.0_dp

  end function calculate_vanderWaals_functional_deriv

  function calculate_surface_electrostatic_functional_deriv(size_array, charge)
    integer, intent(in) :: size_array
    real(dp), intent(in) :: charge

    real(dp), dimension(size_array) :: calculate_surface_electrostatic_functional_deriv

    integer :: ij, iz
    integer :: start_z_index
    integer :: end_z_index

    real(dp) :: d_to_left_wall
    real(dp) :: d_to_right_wall

    call get_allowed_z_values(start_z_index, end_z_index, size_array)

    do ij = start_z_index, end_z_index

       iz = ij !- start_z_index + 1
       d_to_left_wall = real((iz-1) * hs_diameter, dp) / real(n_discretised_points_z, dp)
       d_to_right_wall = real((size_array-iz) * hs_diameter, dp) / real(n_discretised_points_z, dp)


       calculate_surface_electrostatic_functional_deriv(ij) = (-1.0_dp / ( 2.0_dp * epsilon0 * epsilonr)) * (&
            charge * surface_charge_density_left_wall * d_to_left_wall + &
            charge * surface_charge_density_right_wall * d_to_right_wall)
    end do

    call setNonCalculatedRegionToZero(calculate_surface_electrostatic_functional_deriv)
    !calculate_surface_electrostatic_functional_deriv = 0.0_dp

  end function calculate_surface_electrostatic_functional_deriv

  function calculate_electrostatic_like_term_functional_deriv(n, charge, calculate_bulk)
    real(dp), dimension(:), intent(in) :: n
    real(dp), intent(in) :: charge
    logical, intent(in) :: calculate_bulk
    real(dp), dimension(size(n)) :: calculate_electrostatic_like_term_functional_deriv

    real(dp) :: lambda
    real(dp), dimension(size(n)) :: unit_array
    
    ! integer :: start_z_index
    ! integer :: end_z_index
    ! call get_allowed_z_values(start_z_index, end_z_index, size(n))

    if(charge > 0.0_dp) then
       CURRENT_BULK_BEAD_DENSITY = bulk_density_positive_beads
    else if(charge < 0.0_dp) then
       CURRENT_BULK_BEAD_DENSITY = bulk_density_negative_beads
    else if(charge == 0.0_dp) then
       !Note this is a float comparison so should never happen.
       !In any case if charge = 0.0_dp then we get an answer of zero anyway.
       !But we want to set it to something, to protect against zero * {some unset variable},
       !which I don't know what may happen.
       CURRENT_BULK_BEAD_DENSITY = bulk_density_neutral_beads
    end if

    unit_array(:) = 1.0_dp
    
    if(calculate_bulk) then
       calculate_electrostatic_like_term_functional_deriv(:) = 0.0_dp

       !calculate_electrostatic_like_term_functional_deriv(:) = (-1.0_dp / (2.0_dp * epsilon0 * epsilonr)) * (charge**2) * &
       !     CURRENT_BULK_BEAD_DENSITY * integrate_z_cylindrical(unit_array, electrostatic_like_integrand, "all_z")


    else
       calculate_electrostatic_like_term_functional_deriv(:) = (-1.0_dp / (2.0_dp * epsilon0 * epsilonr)) * (charge**2) * &
            integrate_z_cylindrical(n, electrostatic_unlike_integrand, "all_z")


    end if

    call setNonCalculatedRegionToZero(calculate_electrostatic_like_term_functional_deriv)

    !print *, "calculate_electrostatic_like_term_functional_deriv = ", n(200), calculate_electrostatic_like_term_functional_deriv(200)
    !calculate_electrostatic_like_term_functional_deriv = 0.0_dp

  end function calculate_electrostatic_like_term_functional_deriv

  function calculate_electrostatic_unlike_term_functional_deriv(n, charge1, charge2, calculate_bulk)
    real(dp), dimension(:), intent(in) :: n
    real(dp), intent(in) :: charge1
    real(dp), intent(in) :: charge2
    logical, intent(in) :: calculate_bulk
    real(dp), dimension(size(n)) :: calculate_electrostatic_unlike_term_functional_deriv

    real(dp) :: d, h

    !d = chi_parameter*hs_diameter
    !h = (size(n) - hs_diameter/n_discretised_points_z) * (hs_diameter/n_discretised_points_z)

    if(calculate_bulk) then
       calculate_electrostatic_unlike_term_functional_deriv(:) = 0.0_dp

       !calculate_electrostatic_unlike_term_functional_deriv(:) = (-1.0_dp / (2.0_dp * epsilon0 * epsilonr)) * charge1 * charge2 * &
       !     CURRENT_BULK_BEAD_DENSITY * ( d**2 )

    else

       calculate_electrostatic_unlike_term_functional_deriv(:) = (-1.0_dp / (2.0_dp * epsilon0 * epsilonr)) * charge1 * charge2 * &       
            integrate_z_cylindrical(n, electrostatic_unlike_integrand, "all_z")

    end if


    call setNonCalculatedRegionToZero(calculate_electrostatic_unlike_term_functional_deriv)
    !calculate_electrostatic_unlike_term_functional_deriv = 0.0_dp

  end function calculate_electrostatic_unlike_term_functional_deriv

  function electrostatic_unlike_integrand(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp) :: electrostatic_unlike_integrand

    integer :: xi_int
    real(dp) :: xi_real

    xi_int = xi - z ! centre xi on z
    xi_real = real(xi_int,dp) * hs_diameter / real(n_discretised_points_z,dp)

    if(abs(xi_int) >= (chi_parameter * n_discretised_points_z)) then
       electrostatic_unlike_integrand = abs(xi_real)
    else
       electrostatic_unlike_integrand = chi_parameter * hs_diameter
    end if

  end function electrostatic_unlike_integrand


  function electrostatic_like_integrand(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp) :: electrostatic_like_integrand

    real(dp) :: lambda
    real(dp) :: abs_xi

    lambda = get_lambda()
    abs_xi = (abs(z - xi)*hs_diameter/real(n_discretised_points_z,dp))

    electrostatic_like_integrand = (abs_xi) + ((exp(-(abs_xi*lambda)))/lambda)

  end function electrostatic_like_integrand

  function get_lambda()
    real(dp) :: get_lambda
    real(dp) :: s
    
    s = (3.0_dp/(4.0_dp * pi * CURRENT_BULK_BEAD_DENSITY))**(1.0_dp/3.0_dp)
    get_lambda = sqrt(2.0_dp)/s

  end function get_lambda

  function van_der_waals_density_indept_integrand(z, xi_in)
    integer, intent(in) :: z
    integer, intent(in) :: xi_in
    real(dp)             :: van_der_waals_density_indept_integrand

    integer :: xi_int
    real(dp) :: xi_real

    xi_int = xi_in - z ! centre xi on z
    xi_real = real(xi_int,dp) * hs_diameter / real(n_discretised_points_z,dp)

    if(abs(xi_int) >= n_discretised_points_z) then
       van_der_waals_density_indept_integrand = 1.0_dp / (4.0_dp * (real(xi_real,dp)**4.0_dp))
    else
       van_der_waals_density_indept_integrand = 1.0_dp / &
            (4.0_dp * ( (real(hs_diameter,dp) * cos(asin(real(xi_real,dp)/real(hs_diameter,dp))))**2.0_dp &
            + real(xi_real,dp)**2.0_dp )**2.0_dp)
    end if

    return
  end function van_der_waals_density_indept_integrand

  function calculate_n_sbar(array_to_integrate)
    real(dp), dimension(:), intent(in) :: array_to_integrate ! typically n_s, but sometime n_s / n_sbar
    real(dp), dimension(size(array_to_integrate)) :: calculate_n_sbar

    calculate_n_sbar(:) = 0.0_dp

    calculate_n_sbar(:) = (3.0_dp * ( integrate_z_cylindrical(array_to_integrate, n_sbar_integrand, "z_lteq_hs_diameter") ))&
         /(4.0_dp * pi * (hs_diameter**3.0_dp))

  end function calculate_n_sbar

  pure function n_sbar_integrand(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp)            :: n_sbar_integrand

    n_sbar_integrand = 2.0_dp * pi * (hs_diameter**2.0_dp - ((z - xi)*hs_diameter/real(n_discretised_points_z,dp))**2.0_dp) / 2.0_dp

  end function n_sbar_integrand

end module functionalderivatives

! Module containing so-called helper functions.
! Used to perform neccesary logic on data structures.
module helpers
  use kinds
  use parameters
  implicit none
  private

  public :: get_allowed_z_values
  public :: setNonCalculatedRegionToZero
  public :: unity_function
  public :: get_bulk_density
  public :: str
  public :: SetToZero

  interface str
     module procedure str_int
     module procedure str_real_dp
     module procedure str_real_sp
  end interface str

contains

  subroutine get_allowed_z_values(start_z_index, end_z_index, total_points_z)
    integer, intent(out) :: start_z_index
    integer, intent(out) :: end_z_index
    integer, intent(in)  :: total_points_z

    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    !Note we use integer division here
    if(is_even(n_discretised_points_z)) then
       start_z_index = (n_discretised_points_z/2)  + 1
       end_z_index = total_points_z - (n_discretised_points_z/2)
    else
       start_z_index = (n_discretised_points_z/2)  + 2
       end_z_index = total_points_z - (n_discretised_points_z/2) - 1
    end if

    !start_z_index = 1
    !end_z_index = total_points_z
    
  end subroutine get_allowed_z_values

  subroutine setNonCalculatedRegionToZero(array1, array2, array3)
    real(dp), dimension(:), intent(inout) :: array1
    real(dp), dimension(:), intent(inout), optional :: array2
    real(dp), dimension(:), intent(inout), optional :: array3

    integer :: start_z_index
    integer :: end_z_index

    !First ensure all input arguments are the same size.
    if(present(array2)) then
       if(size(array1) /= size(array2)) then
          print *, "helpers.f90: setNonCalculatedRegionToZero: "
          print *, "size mismatch.  Size of all input arguments must match."
          print *, "size(array1) /= size(array2)...aborting..."
          call abort()
       end if
    end if

    if(present(array3)) then
       if(size(array1) /= size(array3)) then
          print *, "helpers.f90: setNonCalculatedRegionToZero: "
          print *, "size mismatch.  Size of all input arguments must match."
          print *, "size(array1) /= size(array3)...aborting..."
          call abort()
       end if
    end if

    call get_allowed_z_values(start_z_index, end_z_index, size(array1))

    array1(1:start_z_index-1) = 0.0_dp
    array1(end_z_index+1:size(array1)) = 0.0_dp

    if(present(array2)) then
       array2(1:start_z_index-1) = 0.0_dp
       array2(end_z_index+1:size(array2)) = 0.0_dp
    end if

    if(present(array3)) then
       array3(1:start_z_index-1) = 0.0_dp
       array3(end_z_index+1:size(array3)) = 0.0_dp
    end if

  end subroutine setNonCalculatedRegionToZero

  pure function is_even(an_int)
    integer, intent(in) :: an_int
    logical             :: is_even

    if(modulo(an_int, 2) == 0) then
       is_even = .true.
    else
       is_even = .false.
    end if

    return
  end function is_even

  pure function unity_function(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp) :: unity_function

    unity_function = 1.0_dp
  end function unity_function

  function get_bulk_density(lambda) result(reslt)
    real(dp), dimension(:), intent(in)  :: lambda
    real(dp) :: reslt

    integer :: start_z_integrate
    integer :: end_z_integrate

    !Find the maximum range over which we are going to integrate
    call  get_allowed_z_values(start_z_integrate, end_z_integrate, size(lambda))

    reslt = sum(lambda(start_z_integrate:end_z_integrate))/real(size(lambda(start_z_integrate:end_z_integrate)), dp)

    return
  end function get_bulk_density

  !Casts an int to a string
  function str_int(x)
    integer, intent(in) :: x
    character(len=64) :: str_int

    write(str_int, *) x

    str_int = adjustl(str_int)
  end function str_int

  !Casts a double precision real to a string
  function str_real_dp(x)
    real(dp), intent(in) :: x
    character(len=64) :: str_real_dp

    write(str_real_dp, '(f10.5)') x

    str_real_dp = adjustl(str_real_dp)
  end function str_real_dp

  !Casts a single precision real to a string
  function str_real_sp(x)
    real(sp), intent(in) :: x
    character(len=64) :: str_real_sp

    write(str_real_sp, '(f10.5)') x

    str_real_sp = adjustl(str_real_sp)
  end function str_real_sp

  subroutine SetToZero(array1, array2, array3, array4, array5, array6)
    real(dp), dimension(:) :: array1
    real(dp), dimension(:), optional :: array2
    real(dp), dimension(:), optional :: array3
    real(dp), dimension(:), optional :: array4
    real(dp), dimension(:), optional :: array5
    real(dp), dimension(:), optional :: array6

    array1(:) = 0.0_dp
    if(present(array2)) array2(:) = 0.0_dp
    if(present(array3)) array3(:) = 0.0_dp
    if(present(array4)) array4(:) = 0.0_dp
    if(present(array5)) array5(:) = 0.0_dp
    if(present(array6)) array6(:) = 0.0_dp

  end subroutine SetToZero

end module helpers

!This module contains the prescription of how to perform an integral
!over the spherical phi coordinate given that the z coordinate in
!cylindrical coordinates in the one that is discretised.
module integratephispherical
  use kinds
  use universalconstants
  use parameters
  use helpers
  implicit none
  private

  public :: integrate_phi_spherical

contains

  !A function that integrates over all possible allowed values of phi
  !in spherical coordinates.  That is, from 0 to pi, with the restriction
  !that we can't integrate through walls. Noting that the integrand is a 1-D disretised
  !function of only z, we follow the algorithm...
  !
  !1. Loop over all values of z_{i} in the integrand array.
  !
  !2. For each z_{i} find the minimum and maximum values of z (z_min, z_max) corresponding
  !to a sphere of radius sigma around the initial z point.
  !
  !3. Form an array subsection of allowed values and apply the trapezoidal rule.  We calculate
  !the jacobian of an associated angle from the difference between z_{i} and the value of z being
  !summed over.
  !
  function integrate_phi_spherical(integrand_array) result(reslt)
    real(dp), dimension(:), intent(in) :: integrand_array
    real(dp), dimension(size(integrand_array)) :: reslt

    integer :: lower_z_limit
    integer :: upper_z_limit

    integer :: relative_z_index
    integer :: ires

    integer :: start_z_index
    integer :: end_z_index

    real(dp) :: surface_area_fraction
    
    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    call get_allowed_z_values(start_z_index, end_z_index, size(reslt))
    
    do ires = start_z_index, end_z_index

       call get_integrand_array_section_limits_and_surface_area_fraction(ires, size(reslt), lower_z_limit, upper_z_limit, &
            relative_z_index, surface_area_fraction)
       
       reslt(ires) = surface_area_fraction * apply_trapezoidal_rule(integrand_array(lower_z_limit:upper_z_limit), relative_z_index)

    end do

    reslt(1:start_z_index-1) = 0.0_dp
    reslt(end_z_index+1:size(reslt)) = 0.0_dp
    
  end function integrate_phi_spherical

  !Returns the Jacobian's phi dependence in spherical coordinates.
  pure function Jacobian_phi_dependence(theta)
    real(dp), intent(in) :: theta
    real(dp) :: Jacobian_phi_dependence

    Jacobian_phi_dependence = sin(theta)

  end function Jacobian_phi_dependence

  !Given an array section 'integrand_array' corresponding to the array section
  !we wish to integrate over and a 'z_index' corresponding to the fixed index that the
  !integral corresponds to, apply the trapezoidal rule.
  function apply_trapezoidal_rule(integrand_array, z_index) result(reslt)
    real(dp), dimension(:), intent(in) :: integrand_array
    integer, intent(in)                :: z_index
    real(dp)                           :: reslt

    real(dp) :: theta_i
    real(dp) :: theta_ip1
    integer  :: i
    reslt = 0.0_dp

    do i = 1, size(integrand_array) - 1

       theta_i =  get_angle_from_z_separation(i, z_index)
       theta_ip1 = get_angle_from_z_separation(i+1, z_index)

       reslt = reslt + ( (integrand_array(i)*Jacobian_phi_dependence(theta_i) + &
            integrand_array(i+1)*Jacobian_phi_dependence(theta_ip1) ) * abs(theta_i - theta_ip1)) &
            /2.0_dp
    end do

  end function apply_trapezoidal_rule

  !Given the fixed value of z from which the functional dependence originates and the value of
  !z being integrated over, we calculate the assocated angle phi (in spherical coordinates)
  !corresponding to their separation.
  function get_angle_from_z_separation(z_integrated, z_fixed) result(angle)
    integer, intent(in) :: z_integrated
    integer, intent(in) :: z_fixed

    real(dp)            :: angle

    if(z_integrated <= z_fixed) then
       angle = pi - acos(real((z_fixed - z_integrated),dp)/real(n_discretised_points_z, dp))
    else
       angle = acos(real((z_integrated - z_fixed),dp)/real(n_discretised_points_z, dp))
    end if

  end function get_angle_from_z_separation
  
  !Given a 'z_index' indicating the value of z that introduces the functional dependence, 'h' = number of valid z
  !discretised values between the plates, this routine calculates the upper and lower possible z limits for the
  !spherical integral and calculates the value of 'z_index' relative to this array section, storing the result in
  !'relative_z_index'.
  subroutine get_integrand_array_section_limits_and_surface_area_fraction(z_index, h, lower_z_limit, upper_z_limit, relative_z_index, surface_area_fraction)
    integer, intent(in)   :: z_index
    integer, intent(in)   :: h
    integer, intent(out)  :: lower_z_limit
    integer, intent(out)  :: upper_z_limit
    integer, intent(out)  :: relative_z_index
    real(dp), intent(out) :: surface_area_fraction

    integer :: lowest_z_calculated
    integer :: highest_z_calculated

    call get_allowed_z_values(lowest_z_calculated, highest_z_calculated, h)

    if(h < 2*n_discretised_points_z + 1) then
       print *, "integratephispherical.f90: get_integrand_array_section_limits:"
       print *, "distance between plates less than twice the hard sphere diameter"
       print *, "Do you really want that small a plate separation?"
       print *, "These short plate separations are currently not supported...aborting..."
       call abort()
    end if

    if(z_index < 3*(n_discretised_points_z)/2 + 1) then

       lower_z_limit = lowest_z_calculated
       upper_z_limit = z_index + n_discretised_points_z
       relative_z_index = z_index - lowest_z_calculated + 1

       !surface_area_fraction = 1.0_dp / (1.0_dp - cos(pi - acos((z_index - lowest_z_calculated)/real(n_discretised_points_z, dp))))
       surface_area_fraction = 0.5_dp
    else if((z_index >= (3*n_discretised_points_z/2) + 1) .and.  ((h - z_index) >= 3*n_discretised_points_z/2)) then

       lower_z_limit = z_index - n_discretised_points_z
       upper_z_limit = z_index + n_discretised_points_z
       relative_z_index = n_discretised_points_z + 1

       surface_area_fraction = 0.5_dp
       
    else if((h - z_index) < 3*(n_discretised_points_z)/2) then

       lower_z_limit = z_index - n_discretised_points_z
       upper_z_limit = highest_z_calculated
       relative_z_index = n_discretised_points_z + 1

       !surface_area_fraction = 1.0_dp / (1.0_dp + ((highest_z_calculated - z_index)/real(n_discretised_points_z, dp)))
       surface_area_fraction = 0.5_dp

    else
       print *, "integratephispherical.f90:get_integrand_array_section_limits:"
       print *, "Invalid values of z_index/h."
       print *, "z_index must be <= h."
       print *, "Almost certainly a coding error as opposed to an input error...aborting..."
       call abort()
    end if

  end subroutine get_integrand_array_section_limits_and_surface_area_fraction

end module integratephispherical

!A module containing the routines to integrate over all z
!in cylindrical coordinates.  Note that it is this variable
!that is discretised.  Also note that the Jacobian
!from  cylindrical coordinates is not present as the 'r'
!direction has already been integrated out.
module integratezcylindrical
  use kinds
  use parameters
  use universalconstants
  use helpers
  implicit none
  private

  public :: integrate_z_cylindrical

  interface integrate_z_cylindrical
     module procedure integrate_z_cylindrical_with_range_array
     module procedure integrate_z_cylindrical_with_range_real
     module procedure integrate_z_cylindrical_without_range
  end interface integrate_z_cylindrical

contains
  
  !When calculating the lambdas for example the integrand array is density dependent, while the integrand function
  !is density independent.
  function integrate_z_cylindrical_with_range_array(integrand_array, integrand_function, integration_range) result(reslt)
    real(dp), dimension(:), intent(in) :: integrand_array
    real(dp), external                 :: integrand_function
    character(len=*)                   :: integration_range

    real(dp), dimension(size(integrand_array)) :: reslt

    integer :: ires
    integer :: lower_z_limit
    integer :: upper_z_limit
    integer :: relative_z_index

    integer :: start_z_index
    integer :: end_z_index

    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    call get_allowed_z_values(start_z_index, end_z_index, size(reslt))

    do ires = start_z_index, end_z_index

       call get_integrand_array_section_limits(trim(integration_range), ires, size(reslt), &
            lower_z_limit, upper_z_limit, relative_z_index)   

       reslt(ires) = apply_trapezoidal_rule(integrand_array(lower_z_limit:upper_z_limit), &
            integrand_function, relative_z_index)

    end do

    reslt(1:start_z_index-1) = 0.0_dp
    reslt(end_z_index+1:size(reslt)) = 0.0_dp
    return

  end function integrate_z_cylindrical_with_range_array

  function integrate_z_cylindrical_with_range_real(integrand_array, integration_range) result(reslt)
    real(dp), dimension(:), intent(in) :: integrand_array
    character(len=*)                   :: integration_range

    real(dp) :: reslt
    real(dp), dimension(size(integrand_array)) :: reslt_array

    reslt_array = integrate_z_cylindrical_with_range_array(integrand_array, unity_function, integration_range)
    reslt = reslt_array(size(reslt_array)/2)
    
  end function integrate_z_cylindrical_with_range_real


  !When calculating the lambdas for example the integrand array is density dependent, while the integrand function
  !is density independent.
  function integrate_z_cylindrical_without_range(integrand_array, integrand_function) result(reslt)
    real(dp), dimension(:), intent(in) :: integrand_array
    real(dp), external                 :: integrand_function

    real(dp) :: reslt

    integer :: start_z_index
    integer :: end_z_index

    integer :: dummy_z_index

    !Doesn't matter what this index is. It's passed to the trapezoidal rule
    !but all calls to this routine are integrating an array over all z and
    !therefore want a number, and are independent of the so-called 'integrand_function'
    !which is this case is just the 'unity_function'.
    dummy_z_index = 1

    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    call get_allowed_z_values(start_z_index, end_z_index, size(integrand_array))

    reslt = apply_trapezoidal_rule(integrand_array(start_z_index:end_z_index), &
         integrand_function, dummy_z_index)

    return

  end function integrate_z_cylindrical_without_range

  !Given an array section 'z_dep_integrand' corresponding to the array section
  !we wish to integrate over and a 'z_index' corresponding to the fixed index that the
  !integral corresponds to, apply the trapezoidal rule.
  function apply_trapezoidal_rule(integrand_array, integrand_function, z_index) result(reslt)
    real(dp), dimension(:), intent(in) :: integrand_array
    real(dp), external                 :: integrand_function
    integer, intent(in)                :: z_index
    real(dp)                           :: reslt

    integer :: ixi
    reslt = 0.0_dp

    do ixi = 1, size(integrand_array) - 1
       reslt = reslt + ( (integrand_array(ixi)*integrand_function(z_index, ixi) +  &
            integrand_array(ixi + 1)*integrand_function(z_index, ixi + 1)) &
            * (hs_diameter/real(n_discretised_points_z,dp)))/2.0_dp
    end do

    return
  end function apply_trapezoidal_rule


  !Given an 'integration_range' that tells us whether we are integrating over everything or values >= or <= the
  !hs_diameter, a 'z_index' indicating the value of z that introduces the functional dependence and 'h' = number of valid z
  !discretised values between the plates, this routine calculates the upper and lower possible z limits for the
  !integral and calculates the value of 'z_index' relative to this array section, storing the result in
  !'relative_z_index'.
  subroutine get_integrand_array_section_limits(integration_range, z_index, h, lower_z_limit, upper_z_limit, relative_z_index)
    character(len=*), intent(in) :: integration_range
    integer, intent(in)          :: z_index
    integer, intent(in)          :: h
    integer, intent(out)         :: lower_z_limit
    integer, intent(out)         :: upper_z_limit
    integer, intent(out)         :: relative_z_index

    integer :: lowest_z_calculated
    integer :: highest_z_calculated

    call get_allowed_z_values(lowest_z_calculated, highest_z_calculated, h)

    if(h < 2*n_discretised_points_z + 1) then
       print *, "integratezcylindrical.f90: get_integrand_array_section_limits:"
       print *, "distance between plates less than twice the hard sphere diameter"
       print *, "Do you really want that small a plate separation?"
       print *, "h = ", h, "2*n_discretised_points_z + 1 = ", 2*n_discretised_points_z + 1
       print *, "These short plate separations are currently not supported...aborting..."
       call abort()
    end if

    if(trim(integration_range) == "all_z" .or. &
         trim(integration_range) == "z_gteq_hs_diameter") then

       lower_z_limit = lowest_z_calculated
       upper_z_limit = highest_z_calculated
       relative_z_index = z_index - lowest_z_calculated + 1

    else if(trim(integration_range) == "z_lteq_hs_diameter") then

       if(z_index < lowest_z_calculated + n_discretised_points_z) then !0 <= z < hs_diameter

          lower_z_limit = lowest_z_calculated
          upper_z_limit = z_index + n_discretised_points_z
          relative_z_index = z_index - lowest_z_calculated + 1

       else if((z_index >= lowest_z_calculated + n_discretised_points_z) .and.  &
            ((highest_z_calculated - z_index) >= n_discretised_points_z)) then !hs_diameter <= z <= h-hs_diameter

          lower_z_limit = z_index - n_discretised_points_z
          upper_z_limit = z_index + n_discretised_points_z
          relative_z_index = n_discretised_points_z + 1

       else if((highest_z_calculated - z_index) >= 0) then !h-hs_diameter < z <= h

          lower_z_limit = z_index - n_discretised_points_z
          upper_z_limit = highest_z_calculated
          relative_z_index = n_discretised_points_z + 1

       else ! z > h which is unphysical
          print *, "integratezcylindrical.f90:get_integrand_array_section_limits:"
          print *, "Invalid values of z_index/h."
          print *, "z_index must be <= h."
          print *, "Almost certainly a coding error as opposed to an input error...aborting..."
          call abort()
       end if

    else
       print *, "integratezcylindrical.f90:get_integrand_array_section_limits: "
       print *, "unsupported value of 'integration_range'"
       print *, "integration range = ", trim(integration_range)
       print *, "Almost certainly a coding error as opposed to an input error"
       print *, "Aborting..."
       call Abort()
    end if
  end subroutine get_integrand_array_section_limits

end module integratezcylindrical

!Contains all io routines
module io
  use kinds
  use parameters
  implicit none
  private

  public :: WriteOutputFormattedAsFunctionOfPosition
  public :: WriteOutputFormattedAsFunctionOfPlateSeparation
  
contains

  !Writes array out to txt file
  subroutine WriteOutputFormattedAsFunctionOfPosition(output_array, file_stub, file_suffix)
    real(dp), dimension(:), intent(in) :: output_array
    character(len=*), intent(in)       :: file_stub
    character(len=*), intent(in)       :: file_suffix

    integer :: iz
    integer :: file_unit
    file_unit = 171

    open(file_unit, file=trim(file_stub)//"-"//trim(file_suffix)//".txt", action='write')
    do iz = 1, size(output_array)
       write(file_unit, *) real((iz-1),dp)/real(n_discretised_points_z,dp), output_array(iz)*(hs_diameter**3)
    end do
    close(file_unit)

  end subroutine WriteOutputFormattedAsFunctionOfPosition

  subroutine WriteOutputFormattedAsFunctionOfPlateSeparation(output, file_stub, file_suffix)
    real(dp), dimension(:), intent(in) :: output
    character(len=*), intent(in)       :: file_stub
    character(len=*), intent(in)       :: file_suffix

    integer :: id
    integer :: file_unit
    file_unit = 173

    !First check the size is correct
    if(size(output) /= size(plate_separations)) then
       print *, "io.f90: WriteGrandPotentialOutputFormatted:"
       print *, "size(output) /= size(plate_separations)"
       print *, "should only have one potential value per plate separation"
       print *, "size mismatch...aborting..."
       call abort()
    else
       open(file_unit, file=trim(file_stub)//"-"//trim(file_suffix)//".txt", action='write')
       do id = 1, size(output)
          write(file_unit, *) plate_separations(id), output(id)
       end do
       close(file_unit)
    end if

  end subroutine WriteOutputFormattedAsFunctionOfPlateSeparation

end module io

!Contains relevant routines for our iterative procedure, such as convergence check.
module iteration
  use kinds
  use parameters
  use helpers
  implicit none
  private

  public :: converged
  public :: InitialiseDensityDiscretisationAndSetIntegrationAnsatz
  public :: InitialiseVariableDiscretisation
  
contains

  !Routine that returns true iff the average squared difference for all beads is
  !less than the iterative_tolerance input by the user.
  function converged(n1_updated, n1, n2_updated, n2, n3_updated, n3)
    real(dp), dimension(:), intent(in) :: n1_updated
    real(dp), dimension(:), intent(in) :: n1

    real(dp), dimension(:), optional, intent(in) :: n2_updated
    real(dp), dimension(:), optional, intent(in) :: n2

    real(dp), dimension(:), optional, intent(in) :: n3_updated
    real(dp), dimension(:), optional, intent(in) :: n3

    logical                            :: converged

    logical :: n1_converged
    logical :: n2_converged
    logical :: n3_converged

    real(dp) :: av_sq_diff
    real(dp) :: bulk_value

    integer :: start_z_index
    integer :: end_z_index

    converged = .false.
    n1_converged = .false.
    n2_converged = .false.
    n3_converged = .false.

    !Get the range of non-zero elements
    call get_allowed_z_values(start_z_index, end_z_index, size(n1))

    !Ensure that we recieve arrays of the correct size
    if(size(n1_updated) /= size(n1)) then
       print *, "iteration.f90: converged:"
       print *, "size(n1_updated) /= size(n1)"
       print *, "array size mismatch...aborting..."
       call abort()
    else

       av_sq_diff = sum((n1_updated(start_z_index:end_z_index)*(hs_diameter**3) - n1(start_z_index:end_z_index)*(hs_diameter**3))**2)&
            /real(size(n1(start_z_index:end_z_index)),dp)
       bulk_value = sum(n1(start_z_index:end_z_index)*(hs_diameter**3))/real(size(n1(start_z_index:end_z_index)),dp)

       ! if(av_sq_diff == 0.0_dp) then
       !    print *, "iteration.f90: converged: "
       !    print *, "average squared difference between dneisty and updated_density == 0.0_dp"
       !    print *, "would appear to be a coding bug...aborting..."
       !    call abort()
       ! end if
       !print *, "n1 = ", n1
       !print *, "n1 = ", n1_updated

       if(av_sq_diff <= iterative_tolerance*bulk_value) then
          n1_converged = .true.
       else
          n1_converged = .false.
       end if
    end if

    !Now do n2
    if(present(n2_updated)) then

       if(.not. present(n2))then
          print *, "iteration.f90: converged:"
          print *, "optional argument error.  n2_updated present while n2 isn't present"
          print *, "...aborting..."
          call abort()
       else

          if(size(n2_updated) /= size(n2)) then
             print *, "iteration.f90: converged:"
             print *, "size(n2_updated) /= size(n2)"
             print *, "array size mismatch...aborting..."
             call abort()
          else

             av_sq_diff = sum((n2_updated(start_z_index:end_z_index)*(hs_diameter**3) - n2(start_z_index:end_z_index)*(hs_diameter**3))**2)/&
                  size(n2(start_z_index:end_z_index))
             bulk_value = sum(n2(start_z_index:end_z_index)*(hs_diameter**3))/size(n2(start_z_index:end_z_index))

             ! if(av_sq_diff == 0.0_dp) then
             !    print *, "iteration.f90: converged: "
             !    print *, "average squared difference between dneisty and updated_density == 0.0_dp"
             !    print *, "would appear to be a coding bug...aborting..."
             !    !call abort()
             ! end if


             !print *, "n1 = ", n2
             !print *, "n1 = ", n2_updated

             if(av_sq_diff <= iterative_tolerance*bulk_value) then
                n2_converged = .true.
             else
                n2_converged = .false.
             end if
          end if
       end if

    else !if n2 not present, it shouldn't stop convergence
       n2_converged = .true.
    end if

    !Now do n3
    if(present(n3_updated)) then

       if(.not. present(n3)) then
          print *, "iteration.f90: converged:"
          print *, "size(n3_updated) /= size(n3)"
          print *, "array size mismatch...aborting..."
          call abort() 
       else

          if(size(n3_updated) /= size(n3)) then
             print *, "iteration.f90: converged:"
             print *, "size(n3_updated) /= size(n3)"
             print *, "array size mismatch...aborting..."
             call abort()
          else

             av_sq_diff = sum((n3_updated(start_z_index:end_z_index)*(hs_diameter**3) - n3(start_z_index:end_z_index)*(hs_diameter**3))**2)/ &
                  size(n3(start_z_index:end_z_index))
             bulk_value = sum(n3(start_z_index:end_z_index)*(hs_diameter**3))/size(n3(start_z_index:end_z_index))

             if(av_sq_diff <= iterative_tolerance*bulk_value) then
                n3_converged = .true.
             else
                n3_converged = .false.
             end if
          end if
       end if

    else
       n3_converged = .true.
    end if

    !Check for convergence
    if(n1_converged .and. n2_converged .and. n3_converged)then

       print *, "iteration.f90: converged: Iterative scheme successfully converged."
       converged = .true.
    else
       print *, "convergence = ", n1_converged, n2_converged, n3_converged
       converged = .false.
    end if

  end function converged

  !Intialises the bead density.  The variables are allocated/reallocated to the appropriate size based on which
  !plate separation 'ith_plate-separation' corresponds to.  In the case of it being the first run through
  !(i.e. ith_plate_separation = 1) they are initialised to be constant.  In the case of subsequent run throughs
  !(i.e. ith_plate_separation > 1) they are rescaled based on the previously converged value by the routine
  !'ReScaleArray'.
  subroutine InitialiseDensityDiscretisationAndSetIntegrationAnsatz(ith_plate_separation, n1, n2, n3)
    integer, intent(in) :: ith_plate_separation
    real(dp), dimension(:), allocatable, intent(inout) :: n1
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n2
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n3

    integer :: new_array_size
    
    !Note we add on one to include values at both endpoints/the walls.
    new_array_size = nint(plate_separations(ith_plate_separation) *  n_discretised_points_z) + 1

    call UpdateArraySize(n1, new_array_size, rescale=allocated(n1), charge='+')
    if(present(n2)) call UpdateArraySize(n2, new_array_size, rescale=allocated(n2), charge='0')
    if(present(n3)) call UpdateArraySize(n3, new_array_size, rescale=allocated(n3), charge='-')

    !We don't calculate with hs_diameter/2 of the wall.  Therefore set it zero.
    !This aids with plotting ease.
    call setNonCalculatedRegionToZero(n1)
    if(present(n2)) call setNonCalculatedRegionToZero(n2)
    if(present(n3)) call setNonCalculatedRegionToZero(n3)
    
  end subroutine InitialiseDensityDiscretisationAndSetIntegrationAnsatz

  !Intialises the all variables (other than the bead densities) that are functions of z.
  !The variables are allocated/reallocated to the appropriate size based on which
  !plate separation 'ith_plate-separation' corresponds to.  They are all subsequently initialised
  !to be a constant.  Note that this constant value is (at the time of writing) always overwritten.
  !
  !Note:
  !In order to support studies with a requirement for a large number of variables for clarity sakes,
  !we need to support an input list of a large number of variables.  As variables argument lists are
  !not allowed in Fortran, as far as I can tell there are two options.
  !
  !1. Explicitly include a long list of optional parameters.
  !2. Put them in an array and pass the single array in.
  !
  !The problem with 2 is that is needs to be an array of an explicit type, which in our case is an
  !array of real numbers.  We therefore need to define a derived type to store this information and hence
  !the main program would look substantially messier due to the need to perform array%val for example.
  !We have therefore chosen to aviod this complication in the main program and choose the unweidly option 1.
  !Perhaps we may switch to option 2 if there ever becomes a requirement for more arguments than we currently
  !have listed.
  subroutine InitialiseVariableDiscretisation(ith_plate_separation,a1,a2,a3,a4,a5,a6,a7, &
       a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28, &
       a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48,a49, &
       a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63)
    integer, intent(in) :: ith_plate_separation
    real(dp), dimension(:), allocatable, intent(inout) :: a1
    real(dp), dimension(:), allocatable, intent(inout), optional :: a2
    real(dp), dimension(:), allocatable, intent(inout), optional :: a3
    real(dp), dimension(:), allocatable, intent(inout), optional :: a4
    real(dp), dimension(:), allocatable, intent(inout), optional :: a5
    real(dp), dimension(:), allocatable, intent(inout), optional :: a6
    real(dp), dimension(:), allocatable, intent(inout), optional :: a7
    real(dp), dimension(:), allocatable, intent(inout), optional :: a8
    real(dp), dimension(:), allocatable, intent(inout), optional :: a9
    real(dp), dimension(:), allocatable, intent(inout), optional :: a10
    real(dp), dimension(:), allocatable, intent(inout), optional :: a11
    real(dp), dimension(:), allocatable, intent(inout), optional :: a12
    real(dp), dimension(:), allocatable, intent(inout), optional :: a13
    real(dp), dimension(:), allocatable, intent(inout), optional :: a14
    real(dp), dimension(:), allocatable, intent(inout), optional :: a15
    real(dp), dimension(:), allocatable, intent(inout), optional :: a16
    real(dp), dimension(:), allocatable, intent(inout), optional :: a17
    real(dp), dimension(:), allocatable, intent(inout), optional :: a18
    real(dp), dimension(:), allocatable, intent(inout), optional :: a19
    real(dp), dimension(:), allocatable, intent(inout), optional :: a20
    real(dp), dimension(:), allocatable, intent(inout), optional :: a21
    real(dp), dimension(:), allocatable, intent(inout), optional :: a22
    real(dp), dimension(:), allocatable, intent(inout), optional :: a23
    real(dp), dimension(:), allocatable, intent(inout), optional :: a24
    real(dp), dimension(:), allocatable, intent(inout), optional :: a25
    real(dp), dimension(:), allocatable, intent(inout), optional :: a26
    real(dp), dimension(:), allocatable, intent(inout), optional :: a27
    real(dp), dimension(:), allocatable, intent(inout), optional :: a28
    real(dp), dimension(:), allocatable, intent(inout), optional :: a29
    real(dp), dimension(:), allocatable, intent(inout), optional :: a30
    real(dp), dimension(:), allocatable, intent(inout), optional :: a31
    real(dp), dimension(:), allocatable, intent(inout), optional :: a32
    real(dp), dimension(:), allocatable, intent(inout), optional :: a33
    real(dp), dimension(:), allocatable, intent(inout), optional :: a34
    real(dp), dimension(:), allocatable, intent(inout), optional :: a35
    real(dp), dimension(:), allocatable, intent(inout), optional :: a36
    real(dp), dimension(:), allocatable, intent(inout), optional :: a37
    real(dp), dimension(:), allocatable, intent(inout), optional :: a38
    real(dp), dimension(:), allocatable, intent(inout), optional :: a39
    real(dp), dimension(:), allocatable, intent(inout), optional :: a40
    real(dp), dimension(:), allocatable, intent(inout), optional :: a41
    real(dp), dimension(:), allocatable, intent(inout), optional :: a42
    real(dp), dimension(:), allocatable, intent(inout), optional :: a43
    real(dp), dimension(:), allocatable, intent(inout), optional :: a44
    real(dp), dimension(:), allocatable, intent(inout), optional :: a45
    real(dp), dimension(:), allocatable, intent(inout), optional :: a46
    real(dp), dimension(:), allocatable, intent(inout), optional :: a47
    real(dp), dimension(:), allocatable, intent(inout), optional :: a48
    real(dp), dimension(:), allocatable, intent(inout), optional :: a49
    real(dp), dimension(:), allocatable, intent(inout), optional :: a50
    real(dp), dimension(:), allocatable, intent(inout), optional :: a51
    real(dp), dimension(:), allocatable, intent(inout), optional :: a52
    real(dp), dimension(:), allocatable, intent(inout), optional :: a53
    real(dp), dimension(:), allocatable, intent(inout), optional :: a54
    real(dp), dimension(:), allocatable, intent(inout), optional :: a55
    real(dp), dimension(:), allocatable, intent(inout), optional :: a56
    real(dp), dimension(:), allocatable, intent(inout), optional :: a57
    real(dp), dimension(:), allocatable, intent(inout), optional :: a58
    real(dp), dimension(:), allocatable, intent(inout), optional :: a59
    real(dp), dimension(:), allocatable, intent(inout), optional :: a60
    real(dp), dimension(:), allocatable, intent(inout), optional :: a61
    real(dp), dimension(:), allocatable, intent(inout), optional :: a62
    real(dp), dimension(:), allocatable, intent(inout), optional :: a63

    integer :: new_array_size

    !Note we add on one to include values at both endpoints/the walls.
    new_array_size = nint(plate_separations(ith_plate_separation) *  n_discretised_points_z) + 1

    call UpdateArraySize(a1, new_array_size)
    if(present(a2)) call UpdateArraySize(a2, new_array_size)
    if(present(a3)) call UpdateArraySize(a3, new_array_size)
    if(present(a4)) call UpdateArraySize(a4, new_array_size)
    if(present(a5)) call UpdateArraySize(a5, new_array_size)
    if(present(a6)) call UpdateArraySize(a6, new_array_size)
    if(present(a7)) call UpdateArraySize(a7, new_array_size)
    if(present(a8)) call UpdateArraySize(a8, new_array_size)
    if(present(a9)) call UpdateArraySize(a9, new_array_size)
    if(present(a10)) call UpdateArraySize(a10, new_array_size)
    if(present(a11)) call UpdateArraySize(a11, new_array_size)
    if(present(a12)) call UpdateArraySize(a12, new_array_size)
    if(present(a13)) call UpdateArraySize(a13, new_array_size)
    if(present(a14)) call UpdateArraySize(a14, new_array_size)
    if(present(a15)) call UpdateArraySize(a15, new_array_size)
    if(present(a16)) call UpdateArraySize(a16, new_array_size)
    if(present(a17)) call UpdateArraySize(a17, new_array_size)
    if(present(a18)) call UpdateArraySize(a18, new_array_size)
    if(present(a19)) call UpdateArraySize(a19, new_array_size)
    if(present(a20)) call UpdateArraySize(a20, new_array_size)
    if(present(a21)) call UpdateArraySize(a21, new_array_size)
    if(present(a22)) call UpdateArraySize(a22, new_array_size)
    if(present(a23)) call UpdateArraySize(a23, new_array_size)
    if(present(a24)) call UpdateArraySize(a24, new_array_size)
    if(present(a25)) call UpdateArraySize(a25, new_array_size)
    if(present(a26)) call UpdateArraySize(a26, new_array_size)
    if(present(a27)) call UpdateArraySize(a27, new_array_size)
    if(present(a28)) call UpdateArraySize(a28, new_array_size)
    if(present(a29)) call UpdateArraySize(a29, new_array_size)
    if(present(a30)) call UpdateArraySize(a30, new_array_size)
    if(present(a31)) call UpdateArraySize(a31, new_array_size)
    if(present(a32)) call UpdateArraySize(a32, new_array_size)
    if(present(a33)) call UpdateArraySize(a33, new_array_size)
    if(present(a34)) call UpdateArraySize(a34, new_array_size)
    if(present(a35)) call UpdateArraySize(a35, new_array_size)
    if(present(a36)) call UpdateArraySize(a36, new_array_size)
    if(present(a37)) call UpdateArraySize(a37, new_array_size)
    if(present(a38)) call UpdateArraySize(a38, new_array_size)
    if(present(a39)) call UpdateArraySize(a39, new_array_size)
    if(present(a40)) call UpdateArraySize(a40, new_array_size)
    if(present(a41)) call UpdateArraySize(a41, new_array_size)
    if(present(a42)) call UpdateArraySize(a42, new_array_size)
    if(present(a43)) call UpdateArraySize(a43, new_array_size)
    if(present(a44)) call UpdateArraySize(a44, new_array_size)
    if(present(a45)) call UpdateArraySize(a45, new_array_size)
    if(present(a46)) call UpdateArraySize(a46, new_array_size)
    if(present(a47)) call UpdateArraySize(a47, new_array_size)
    if(present(a48)) call UpdateArraySize(a48, new_array_size)
    if(present(a49)) call UpdateArraySize(a49, new_array_size)
    if(present(a50)) call UpdateArraySize(a50, new_array_size)
    if(present(a51)) call UpdateArraySize(a51, new_array_size)
    if(present(a52)) call UpdateArraySize(a52, new_array_size)
    if(present(a53)) call UpdateArraySize(a53, new_array_size)
    if(present(a54)) call UpdateArraySize(a54, new_array_size)
    if(present(a55)) call UpdateArraySize(a55, new_array_size)
    if(present(a56)) call UpdateArraySize(a56, new_array_size)
    if(present(a57)) call UpdateArraySize(a57, new_array_size)
    if(present(a58)) call UpdateArraySize(a58, new_array_size)
    if(present(a59)) call UpdateArraySize(a59, new_array_size)
    if(present(a60)) call UpdateArraySize(a60, new_array_size)
    if(present(a61)) call UpdateArraySize(a62, new_array_size)
    if(present(a62)) call UpdateArraySize(a62, new_array_size)
    if(present(a63)) call UpdateArraySize(a63, new_array_size)

  end subroutine InitialiseVariableDiscretisation

  !Routine that either resizes or rescales the array, based on the presence and value
  !of the optional parameter 'rescale'.  Also, initialises the array values.
  subroutine UpdateArraySize(array, new_size, rescale, charge)
    real(dp), dimension(:), allocatable, intent(inout) :: array
    integer                                            :: new_size
    logical, intent(in), optional                      :: rescale
    character(len=1), intent(in), optional             :: charge

    if(allocated(array)) then

       if(size(array) == new_size) then
          print *, "iteration.f90:UpdateArraySize:"
          print *, "Attempting to repeat the calculation at the same plate separation."
          print *, "There is no reason to do this."
          print *, "Please ensure the list of plate separations do not contain duplicates"
          call abort()   
       else

          if(present(rescale)) then

             if(rescale) then
                call ReScaleArray(array, new_size)
             else
                call ReSizeArray(array, new_size)
             end if

          else
             !call ReSizeArray(array, new_size)
             call ReScaleArray(array, new_size)
          end if

       end if

    else
       if(present(rescale)) then
          if(.not. present(charge)) then
             print *, "iteration.f90: UpdateArraySize:"
             print *, "If present(rescale) then charge must also be present...aborting..."
             call abort()
          end if
          if(rescale) then
             print *, "iteration.f90: UpdateArraySize:"
             print *, "Can't rescale an array that has not yet been allocated"
             call abort()
          end if
       end if

       allocate(array(new_size)) !size = M x N
       if(.not. present(charge)) then
          array(:) = 0.0_dp
       else
          call InitialiseIntegrationAnsatz(array, charge)
       end if
       
    end if

  end subroutine UpdateArraySize

  !Routine that resizes 'array' to 'new_size' and initialises all elements to 0.
  subroutine ReSizeArray(array, new_size)
    real(dp), dimension(:), allocatable, intent(inout) :: array
    integer, intent(in)                                :: new_size

    if(size(array) == new_size) then
       print *, "iteration.f90:ReSizeArray:"
       print *, "The size of the array is equal to the proposed new size."
       print *, "No need to call resize."
       print *, "Probable coding error...aborting..."
       call abort()
    end if

    if(allocated(array)) deallocate(array)
    allocate(array(new_size))
    array(:) = 0.0_dp

  end subroutine ReSizeArray

  !Rescales 'array' to an array of size 'new_size'.  The ith value of the new array
  !(the array value at the end of the routine) is given by the value of the old array
  !(the array value at the start of the routine) that is the closest to the same fraction
  !of the total size as the value in the new array that it is setting.
  subroutine ReScaleArray(array, new_size)
    real(dp), dimension(:), allocatable, intent(inout) :: array
    integer, intent(in)                                :: new_size

    real(dp), dimension(size(array)) :: old_array_values
    integer :: ith_component
    integer :: old_value_index

    integer :: start_z_index_new
    integer :: end_z_index_new

    integer :: start_z_index_old
    integer :: end_z_index_old
    
    if(.not. allocated(array)) then
       print *, "iteration.f90: Unable to rescale array that isn't allocated"
       print *, "This is almost certainly a coding error"
       print *, "as opposed to an error with the input parms"
       call abort()
    end if

    !Store the values of the array to be deallocated
    old_array_values = array

    !Reallocate the array to be the new size
    if(allocated(array)) deallocate(array)
    allocate(array(new_size))
    array = 0.0_dp

    !Find the range over which the values are non zero
    call get_allowed_z_values(start_z_index_new, end_z_index_new, size(array))
    call get_allowed_z_values(start_z_index_old, end_z_index_old, size(old_array_values))
    
    !Set the elements of the new array by rescaling the old values
    do ith_component = start_z_index_new, end_z_index_new
       old_value_index = start_z_index_old - 1 +  &
            ceiling((end_z_index_old - start_z_index_old + 1) * ( real(ith_component - start_z_index_new + 1, dp) &
            / real(end_z_index_new - start_z_index_new + 1, dp) ) )
       array(ith_component) = old_array_values(old_value_index)
    end do
  end subroutine ReScaleArray

  !Routine to initialise our ansatz for our integrative scheme to an arbitrary constant.
  subroutine InitialiseIntegrationAnsatz(array, charge)
    real(dp), dimension(:) :: array
    character(len=1), intent(in) :: charge

    integer :: midpoint
    integer :: ij
    
    
    if(charge == '+') then

       midpoint = size(array)/2

       do ij = 1, size(array)
          array(ij) = bulk_density_positive_beads + ((ij - midpoint)*hs_diameter/real(n_discretised_points_z,dp))*slope_for_initial_guess
       end do

    else if(charge == '0') then
       array(:) = bulk_density_neutral_beads

    else if(charge == '-') then

       midpoint = size(array)/2

       do ij = 1, size(array)
          array(ij) = bulk_density_negative_beads - ((ij - midpoint)*hs_diameter/real(n_discretised_points_z, dp))*slope_for_initial_guess
       end do
       
    else
       print *, "iteration.f90: InitialiseIntegrationAnsatz"
       print *, "charge has an illegal value of ", charge, "...aborting..."
       call abort()
    end if
    
    !array(:) = bulk_density_positive_beads
  end subroutine InitialiseIntegrationAnsatz

end module iteration

!A module that stores kind type parameters corresponding to the various precisions we may need.
module Kinds
  implicit none
  public

  integer, parameter :: i32  = selected_int_kind(9)  ! 32 bit integer.
  integer, parameter :: i64  = selected_int_kind(18) ! 64 bit integer.

  integer, parameter :: sp = kind(1.0) ! Single precision real scalars.
  integer, parameter :: dp = kind(1.0D0) ! Double precision real scalars.
  !kind(1.0D0) ! Double precision real scalars.
  !selected_real_kind(33, 4931) !Quadruple precision real scalars

  integer, parameter :: sc = kind((1.0,1.0)) ! Single precision complex scalars.
  integer, parameter :: dc = kind((1.0D0,1.0D0)) ! Double precision complex scalars.
  
end module Kinds

!This is a module that calculates the lambdas of our model.
!Recall that lambda is defined to be the paritial derivative of
!all non-ideal contributions to the grand potential functional
!w.r.t the bead densities.
module lambdas
  use kinds
  use helpers
  use parameters
  use functionalderivatives
  implicit none
  private

  public :: CalculateLambdasDifference
  public :: CalculateLambdasBulk
  public :: CalculateLambdasDifference_old

contains
  
  subroutine CalculateLambdas(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus
    real(dp), dimension(:), intent(in)  :: n_plus
    
    real(dp), dimension(:), intent(out) :: lambda_neutral
    real(dp), dimension(:), intent(in)  :: n_neutral
    
    real(dp), dimension(:), intent(out) :: lambda_minus
    real(dp), dimension(:), intent(in)  :: n_minus
    integer, intent(in)                 :: ith_plate_separation
    
    real(dp), dimension(size(lambda_plus)) :: lambda_common_terms
    
    integer :: input_array_size

    ! First ensure that the sizes of all the input variables are the same.
    input_array_size = size(lambda_plus)
    if((size(lambda_neutral) /= input_array_size) .or. &
         (size(lambda_minus) /= input_array_size) .or. &
         (size(n_plus) /= input_array_size) .or. &
         (size(n_neutral) /= input_array_size) .or. &
         (size(n_minus) /= input_array_size)) then

       print *, "lambda.f90: CalculateLambdas:"
       print *, "Input arrays must all be the same size"
       print *, "Almost certainly coding error as opposed to input error...aborting..."
       call abort()
    end if

    lambda_plus(:) = 0.0_dp
    lambda_neutral(:) = 0.0_dp
    lambda_minus(:) = 0.0_dp

    
    ! Now calculate the terms in common to all the lambdas
    lambda_common_terms = CalculateLambdaCommonTerms(n_plus, n_neutral, n_minus, ith_plate_separation, .false.)
    
    ! Now calculate our lambdas(r)
    lambda_plus = beta * (lambda_common_terms + CalculateLambdaPlusSpecificTerms(n_plus, n_minus, .false.))
    lambda_neutral = beta * (lambda_common_terms + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .false.))
    lambda_minus = beta * (lambda_common_terms + CalculateLambdaMinusSpecificTerms(n_plus, n_minus, .false.))
    
  end subroutine CalculateLambdas

  subroutine CalculateLambdasBulk(lambda_plus_bulk, n_plus, lambda_neutral_bulk, n_neutral, lambda_minus_bulk, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus_bulk
    real(dp), dimension(:), intent(in)  :: n_plus

    real(dp), dimension(:), intent(out) :: lambda_neutral_bulk
    real(dp), dimension(:), intent(in)  :: n_neutral

    real(dp), dimension(:), intent(out) :: lambda_minus_bulk
    real(dp), dimension(:), intent(in)  :: n_minus
    integer, intent(in)                 :: ith_plate_separation

    real(dp), dimension(size(lambda_plus_bulk)) :: lambda_common_terms_bulk
    real(dp), dimension(size(n_neutral)) :: n_neutral_array
    
    integer :: input_array_size

    ! First ensure that the sizes of all the input variables are the same.
    input_array_size = size(lambda_plus_bulk)
    if((size(lambda_neutral_bulk) /= input_array_size) .or. &
         (size(lambda_minus_bulk) /= input_array_size) .or. &
         (size(n_plus) /= input_array_size) .or. &
         (size(n_neutral) /= input_array_size) .or. &
         (size(n_minus) /= input_array_size)) then

       print *, "lambda.f90: CalculateLambdasBulk:"
       print *, "Input arrays must all be the same size"
       print *, "Almost certainly coding error as opposed to input error...aborting..."
       call abort()
    end if

    lambda_plus_bulk(:) = 0.0_dp
    lambda_neutral_bulk(:) = 0.0_dp
    lambda_minus_bulk(:) = 0.0_dp
    
    lambda_common_terms_bulk = CalculateLambdaCommonTerms(n_plus, n_neutral, n_minus, ith_plate_separation, .true.)

    lambda_plus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaPlusSpecificTerms(n_plus, n_minus, .true.))
    lambda_neutral_bulk = beta * (lambda_common_terms_bulk + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .true.))
    lambda_minus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaMinusSpecificTerms(n_plus, n_minus, .true.))

  end subroutine CalculateLambdasBulk

  subroutine CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus
    real(dp), dimension(:), intent(in)  :: n_plus

    real(dp), dimension(:), intent(out) :: lambda_neutral
    real(dp), dimension(:), intent(in)  :: n_neutral

    real(dp), dimension(:), intent(out) :: lambda_minus
    real(dp), dimension(:), intent(in)  :: n_minus
    integer, intent(in)                 :: ith_plate_separation

    real(dp), dimension(size(lambda_plus)) :: lambda_plus_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_neutral_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_minus_bulk

    integer :: input_array_size

    ! First ensure that the sizes of all the input variables are the same.
    input_array_size = size(lambda_plus)
    if((size(lambda_neutral) /= input_array_size) .or. &
         (size(lambda_minus) /= input_array_size) .or. &
         (size(n_plus) /= input_array_size) .or. &
         (size(n_neutral) /= input_array_size) .or. &
         (size(n_minus) /= input_array_size)) then

       print *, "lambda.f90: CalculateLambdasBulk:"
       print *, "Input arrays must all be the same size"
       print *, "Almost certainly coding error as opposed to input error...aborting..."
       call abort()
    end if

    call CalculateLambdasBulk(lambda_plus_bulk, n_plus, lambda_neutral_bulk, n_neutral, lambda_minus_bulk, n_minus, ith_plate_separation)
    call CalculateLambdas(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)

    lambda_plus(:) = lambda_plus_bulk(:) - lambda_plus(:)
    lambda_neutral(:) = lambda_neutral_bulk(:) - lambda_neutral(:)
    lambda_minus(:) = lambda_minus_bulk(:) - lambda_minus(:)

  end subroutine CalculateLambdasDifference



  subroutine CalculateLambdasDifference_old(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus
    real(dp), dimension(:), intent(in)  :: n_plus

    real(dp), dimension(:), intent(out) :: lambda_neutral
    real(dp), dimension(:), intent(in)  :: n_neutral

    real(dp), dimension(:), intent(out) :: lambda_minus
    real(dp), dimension(:), intent(in)  :: n_minus
    integer, intent(in)                 :: ith_plate_separation


    real(dp), dimension(size(lambda_plus)) :: lambda_common_terms
    real(dp), dimension(size(lambda_plus)) :: lambda_common_terms_bulk

    real(dp), dimension(size(lambda_plus)) :: lambda_plus_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_neutral_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_minus_bulk

    integer :: input_array_size

    ! First ensure that the sizes of all the input variables are the same.
    input_array_size = size(lambda_plus)
    if((size(lambda_neutral) /= input_array_size) .or. &
         (size(lambda_minus) /= input_array_size) .or. &
         (size(n_plus) /= input_array_size) .or. &
         (size(n_neutral) /= input_array_size) .or. &
         (size(n_minus) /= input_array_size)) then

       print *, "lambda.f90: CalculateLambdas:"
       print *, "Input arrays must all be the same size"
       print *, "Almost certainly coding error as opposed to input error...aborting..."
       call abort()
    end if

    lambda_plus(:) = 0.0_dp
    lambda_neutral(:) = 0.0_dp
    lambda_minus(:) = 0.0_dp

    lambda_plus_bulk(:) = 0.0_dp
    lambda_neutral_bulk(:) = 0.0_dp
    lambda_minus_bulk(:) = 0.0_dp

    ! Now calculate the terms in common to all the lambdas
    lambda_common_terms = CalculateLambdaCommonTerms(n_plus, n_neutral, n_minus, ith_plate_separation, .false.)

    lambda_common_terms_bulk = CalculateLambdaCommonTerms(n_plus, n_neutral, n_minus, ith_plate_separation, .true.)

    ! Now calculate our lambdas(r)
    lambda_plus = beta * (lambda_common_terms + CalculateLambdaPlusSpecificTerms(n_plus, n_minus, .false.))
    lambda_neutral = beta * (lambda_common_terms + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .false.))
    lambda_minus = beta * (lambda_common_terms + CalculateLambdaMinusSpecificTerms(n_plus, n_minus, .false.))

    lambda_plus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaPlusSpecificTerms(n_plus, n_minus, .true.))
    lambda_neutral_bulk = beta * (lambda_common_terms_bulk + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .true.))
    lambda_minus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaMinusSpecificTerms(n_plus, n_minus, .true.))


    !Now calculate the bulk value, lambda^{bulk} and in order to return e^{lambda^{bulk} - lambda(r)}
    lambda_plus = lambda_plus_bulk - lambda_plus
    lambda_neutral = lambda_neutral_bulk - lambda_neutral
    lambda_minus = lambda_minus_bulk - lambda_minus


  end subroutine CalculateLambdasDifference_old

  function CalculateLambdaCommonTerms(n_plus, n_neutral, n_minus, ith_plate_separation, calculate_bulk, iteration)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
    real(dp), dimension(:), intent(in)  :: n_minus
    integer, intent(in)                 :: ith_plate_separation
    logical, intent(in)                 :: calculate_bulk
    integer, optional :: iteration

    real(dp), dimension(size(n_plus)) :: CalculateLambdaCommonTerms

    real(dp), dimension(size(n_plus)) :: hs_term
    real(dp), dimension(size(n_plus)) :: van_der_waals_term
    real(dp), dimension(size(n_plus)) :: surface_fluid_dispersion_term

    real(dp), dimension(size(n_plus)) :: n_s

    hs_term = 0.0_dp
    van_der_waals_term = 0.0_dp
    surface_fluid_dispersion_term = 0.0_dp

    if(calculate_bulk) then
       n_s(:) = bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads

       hs_term = 1.0_dp * calculate_hardsphere_functional_deriv(n_s, .true.)
       van_der_waals_term = calculate_vanderWaals_functional_deriv(n_s)
       
    else
       n_s = n_plus + n_neutral + n_minus
       hs_term = 1.0_dp * calculate_hardsphere_functional_deriv(n_s, .false.)

       surface_fluid_dispersion_term = calculate_surface_dispersion_functional_deriv(&
            ith_plate_separation, size(surface_fluid_dispersion_term))

       van_der_waals_term = calculate_vanderWaals_functional_deriv(n_s)
    end if


    CalculateLambdaCommonTerms = hs_term + surface_fluid_dispersion_term + van_der_waals_term



   end function CalculateLambdaCommonTerms


  function CalculateLambdaPlusSpecificTerms(n_plus, n_minus, calculate_bulk)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_minus
    logical, intent(in) :: calculate_bulk

    real(dp), dimension(size(n_plus)) :: CalculateLambdaPlusSpecificTerms

    real(dp), dimension(size(n_plus)) :: surface_electrostatic_term
    real(dp), dimension(size(n_plus)) :: like_electrostatic_term
    real(dp), dimension(size(n_plus)) :: unlike_electrostatic_term

    real(dp), dimension(size(n_plus))  :: n_plus_in
    real(dp), dimension(size(n_minus))  :: n_minus_in

    surface_electrostatic_term = 0.0_dp
    like_electrostatic_term = 0.0_dp
    unlike_electrostatic_term = 0.0_dp

    if(calculate_bulk) then

       n_plus_in(:) = bulk_density_positive_beads
       n_minus_in(:) = bulk_density_negative_beads

       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_plus_in, positive_bead_charge, .true.)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_minus_in, positive_bead_charge, negative_bead_charge, .true.)

    else
       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_plus, positive_bead_charge, .false.)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_minus, positive_bead_charge, negative_bead_charge, .false.)

       surface_electrostatic_term = calculate_surface_electrostatic_functional_deriv(size(n_plus), positive_bead_charge)


       print *, "surface_electrostatic_term = ",  surface_electrostatic_term(51), surface_electrostatic_term(size(surface_electrostatic_term) - 50), &
            surface_electrostatic_term(size(surface_electrostatic_term) - 50) - surface_electrostatic_term(51)
       print *, "like_electrostatic_term plus = ",  like_electrostatic_term(51), like_electrostatic_term(size(like_electrostatic_term) - 50), &
            like_electrostatic_term(size(like_electrostatic_term) - 50) - like_electrostatic_term(51)
       print *, "unlike_electrostatic_term plus = ", unlike_electrostatic_term(51), unlike_electrostatic_term(size(unlike_electrostatic_term) - 50), &
            unlike_electrostatic_term(size(unlike_electrostatic_term) - 50) - unlike_electrostatic_term(51)

    end if



    CalculateLambdaPlusSpecificTerms = surface_electrostatic_term + like_electrostatic_term + unlike_electrostatic_term


  end function CalculateLambdaPlusSpecificTerms


  function CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, calculate_bulk)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
    
    real(dp), dimension(:), intent(in)  :: n_minus
    logical, intent(in) :: calculate_bulk

    real(dp), dimension(size(n_plus)) :: CalculateLambaNeutralSpecificTerms

    CalculateLambaNeutralSpecificTerms = 0.0_dp

  end function CalculateLambaNeutralSpecificTerms


  function CalculateLambdaMinusSpecificTerms(n_plus, n_minus, calculate_bulk)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_minus
    logical, intent(in) :: calculate_bulk

    real(dp), dimension(size(n_plus)) :: CalculateLambdaMinusSpecificTerms

    real(dp), dimension(size(n_plus)) :: surface_electrostatic_term
    real(dp), dimension(size(n_plus)) :: like_electrostatic_term
    real(dp), dimension(size(n_plus)) :: unlike_electrostatic_term

    real(dp), dimension(size(n_plus))  :: n_plus_in
    real(dp), dimension(size(n_minus))  :: n_minus_in

    surface_electrostatic_term = 0.0_dp
    like_electrostatic_term = 0.0_dp
    unlike_electrostatic_term = 0.0_dp

    if(calculate_bulk) then

       n_plus_in(:) = bulk_density_positive_beads
       n_minus_in(:) = bulk_density_negative_beads

       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_minus_in, negative_bead_charge, .true.)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_plus_in, positive_bead_charge, negative_bead_charge, .true.)

    else

       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_minus, negative_bead_charge, .false.)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_plus, positive_bead_charge, negative_bead_charge, .false.)

       surface_electrostatic_term = calculate_surface_electrostatic_functional_deriv(size(n_minus), negative_bead_charge)


    end if

    CalculateLambdaMinusSpecificTerms = surface_electrostatic_term + like_electrostatic_term + unlike_electrostatic_term

  end function CalculateLambdaMinusSpecificTerms

end module lambdas

!Module storing routines that read and initialise all required parameters.
module parameters
  use kinds
  use universalconstants
  implicit none
  public

  public :: InitialiseModelParameters

  public :: ionic_liquid_name 

  ! 'd' = chi_parameter * hs_diameter provides the only adjusted
  ! parameter in the model. It is fitted to give the correct bulk
  ! ion density for the particlar RTIL at a given pressure and temp.
  ! Only present in the unlike electric correlation term.
  public :: chi_parameter

  ! Number of discretised points in the z direction
  ! over range of hs_diameter. (i.e. >= x and < x + hs_diameter)
  public :: n_discretised_points_z

  public :: epsilonr
  public :: epsilon_LJ
  public :: surface_charge_density_left_wall
  public :: surface_charge_density_right_wall
  public :: hs_diameter
  public :: a_term_index !Index on a_ex in hard sphere term
  public :: plate_separations ! in multiples of hs_diameter
  public :: bulk_density
  public :: temperature
  public :: positive_bead_charge
  public :: negative_bead_charge
  public :: string_length
  public :: iterative_tolerance
  public :: max_iteration_limit
  public :: beta
  public :: alpha_mixing_for_update

  
  public :: bulk_density_positive_beads
  public :: bulk_density_neutral_beads
  public :: bulk_density_negative_beads

  public :: slope_for_initial_guess
  public :: n_charge_iterations
  
  character(len=256) :: ionic_liquid_name 
  real(dp) :: chi_parameter
  integer  :: n_discretised_points_z
  real(dp) :: epsilonr
  real(dp) :: epsilon_LJ
  real(dp) :: surface_charge_density_left_wall
  real(dp) :: surface_charge_density_right_wall
  real(dp) :: hs_diameter
  integer  :: a_term_index

  !array of plater separations in multiples of hs_diameter
  real(dp), dimension(:), allocatable  :: plate_separations

  real(dp) :: bulk_density
  real(dp) :: temperature
  real(dp) :: positive_bead_charge
  real(dp) :: negative_bead_charge
  real(dp) :: string_length
  real(dp) :: iterative_tolerance
  integer  :: max_iteration_limit
  real(dp) :: beta
  real(dp) :: alpha_mixing_for_update
  
  real(dp) :: bulk_density_positive_beads
  real(dp) :: bulk_density_neutral_beads
  real(dp) :: bulk_density_negative_beads

  real(dp) :: slope_for_initial_guess
  integer :: n_charge_iterations
  
contains

  !Subroutine that reads in all required params
  subroutine InitialiseModelParameters(file_stub)
    character(len=*), intent(in) :: file_stub

    integer :: iseparation
    integer :: n_plate_separations
    integer :: file_unit
    file_unit = 171

    open(file_unit, file=trim(file_stub)//".params", action='read')

    read(file_unit, *) ionic_liquid_name
    read(file_unit, *) chi_parameter
    read(file_unit, *) epsilonr
    read(file_unit, *) epsilon_LJ !
    read(file_unit, *) surface_charge_density_left_wall
    read(file_unit, *) surface_charge_density_right_wall
    read(file_unit, *) hs_diameter
    read(file_unit, *) a_term_index
    read(file_unit, *) bulk_density !
    read(file_unit, *) temperature
    read(file_unit, *) alpha_mixing_for_update
    read(file_unit, *) slope_for_initial_guess
    read(file_unit, *) n_charge_iterations
    read(file_unit, *) positive_bead_charge
    read(file_unit, *) negative_bead_charge
    read(file_unit, *) string_length
    read(file_unit, *) n_discretised_points_z
    read(file_unit, *) max_iteration_limit
    read(file_unit, *) iterative_tolerance
    read(file_unit, *) n_plate_separations

    allocate(plate_separations(n_plate_separations))
    do iseparation = 1, n_plate_separations
       read(file_unit, *) plate_separations(iseparation)
    end do

    close(file_unit)

    ! Set derived parameters
    beta = 1.0_dp / (k_B * temperature)
    
    ! Apply unit transformations
    epsilon_LJ = epsilon_LJ * k_B
    bulk_density = bulk_density / (hs_diameter**3.0_dp)
    
    positive_bead_charge = positive_bead_charge * electric_charge
    negative_bead_charge = negative_bead_charge * electric_charge

    surface_charge_density_left_wall = surface_charge_density_left_wall * electric_charge
    surface_charge_density_right_wall = surface_charge_density_right_wall * electric_charge
    

    print *,  "Succesfully set the following values"
    print *,  "ionic_liquid_name = ", ionic_liquid_name
    print *,  "chi_parameter = ", chi_parameter
    print *,  "epsilonr = ", epsilonr
    print *,  "epsilon_LJ = ", epsilon_LJ
    print *,  "surface_charge_density_left_wall = ", surface_charge_density_left_wall
    print *,  "surface_charge_density_right_wall = ", surface_charge_density_right_wall
    print *,  "hs_diameter = ", hs_diameter
    print *,  "a_term_index = ", a_term_index
    print *,  "bulk_density = ", bulk_density
    print *,  "temperature = ", temperature
    print *,  "alpha_mixing_for_update = ", alpha_mixing_for_update
    print *,  "slope_for_initial_guess = ", slope_for_initial_guess
    print *,  "n_charge_iterations = ", n_charge_iterations
    print *,  "positive_bead_charge = ", positive_bead_charge
    print *,  "negative_bead_charge = ", negative_bead_charge
    print *,  "string_length = ", string_length
    print *,  "n_discretised_points_z = ", n_discretised_points_z
    print *,  "max_iteration_limit = ", max_iteration_limit
    print *,  "iterative tolerance = ", iterative_tolerance
    print *,  "n_plate_separations = ", n_plate_separations
    print *,  "plate separations are: ", plate_separations
    print *,  "thermodynamic beta = ", beta

    print *, "Setting Bead Densities from the Bulk Ion Density."
    call SetBeadDensityFromBulkIonDensity(ionic_liquid_name)

    print *, "Checking plate separations are valid given discretisation scheme."
    call CheckValidityOfPlateSeparations()

  end subroutine InitialiseModelParameters

  subroutine DeAllocateModelParams()

    if(allocated(plate_separations)) deallocate(plate_separations)

  end subroutine DeAllocateModelParams

  subroutine SetBeadDensityFromBulkIonDensity(ionic_liquid_name)
    character(len=*), intent(in) :: ionic_liquid_name

    if(trim(ionic_liquid_name) == "SingleNeutralSpheres") then
       call SetSingleNeutralSphereBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "SinglePositiveSpheres") then
       call SetSinglePositiveSphereBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "SingleNegativeSpheres") then
       call SetSingleNegativeSphereBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "NeutralDimers") then
       call SetNeutralDimerDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "C4MIM_BF4-") then
       call SetC4MIN_BF4BeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "PositiveMinusSpheres") then
       call SetPositiveMinusBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "PositiveNeutralDimerMinusSpheres") then
       call SetPlusNeutralDimerMinusSpheresBeadDensityFromBulkIonDensity()

    else if(trim(ionic_liquid_name) == "PositiveNeutralDoubleDimerMinusDimer") then
       call SetDimerDoubleDimerBeadDensityFromBulkIonDensity()


    else
       print *, "parameters.f90: SetBeadDensityFromBulkIonDensity:"
       print *, "Unsupported 'ionic_liquid_name' value of ", trim(ionic_liquid_name)
       print *, "...aborting..."
       call abort()
    end if

  end subroutine SetBeadDensityFromBulkIonDensity

  subroutine SetSingleNeutralSphereBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 0.0_dp
    bulk_density_neutral_beads = bulk_density
    bulk_density_negative_beads = 0.0_dp
  end subroutine SetSingleNeutralSphereBeadDensityFromBulkIonDensity

  subroutine SetSinglePositiveSphereBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = bulk_density
    bulk_density_neutral_beads = 0.0_dp
    bulk_density_negative_beads = 0.0_dp
  end subroutine SetSinglePositiveSphereBeadDensityFromBulkIonDensity

  subroutine SetSingleNegativeSphereBeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 0.0_dp
    bulk_density_neutral_beads = 0.0_dp
    bulk_density_negative_beads = bulk_density
  end subroutine SetSingleNegativeSphereBeadDensityFromBulkIonDensity

  subroutine SetNeutralDimerDensityFromBulkIonDensity()

    bulk_density_positive_beads = 0.0_dp
    bulk_density_neutral_beads = 2.0_dp * bulk_density
    bulk_density_negative_beads = 0.0_dp
  end subroutine SetNeutralDimerDensityFromBulkIonDensity

  subroutine SetC4MIN_BF4BeadDensityFromBulkIonDensity()

    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 5.0_dp * bulk_density
    bulk_density_negative_beads = 5.0_dp * bulk_density
  end subroutine SetC4MIN_BF4BeadDensityFromBulkIonDensity

  subroutine SetPositiveMinusBeadDensityFromBulkIonDensity()
    
    bulk_density_positive_beads = bulk_density
    bulk_density_neutral_beads = 0.0_dp
    bulk_density_negative_beads = bulk_density
  end subroutine SetPositiveMinusBeadDensityFromBulkIonDensity

  subroutine SetPlusNeutralDimerMinusSpheresBeadDensityFromBulkIonDensity()
    
    bulk_density_positive_beads = bulk_density
    bulk_density_neutral_beads = bulk_density
    bulk_density_negative_beads = bulk_density

  end subroutine SetPlusNeutralDimerMinusSpheresBeadDensityFromBulkIonDensity

  subroutine SetDimerDoubleDimerBeadDensityFromBulkIonDensity()
    
    bulk_density_positive_beads = 2.0_dp * bulk_density
    bulk_density_neutral_beads = 2.0_dp * bulk_density
    bulk_density_negative_beads = 2.0_dp * bulk_density

  end subroutine SetDimerDoubleDimerBeadDensityFromBulkIonDensity

  subroutine CheckValidityOfPlateSeparations()

    integer :: ith_separation

    real(dp) :: distance_between_z_values
    real(dp) :: distance_beyond_hs_diameter_multiple
    real(dp) :: n_extra_discretisation_points

    distance_between_z_values = hs_diameter / n_discretised_points_z

    !Check that all the separations are a multiple of the distance between allowed points
    do ith_separation = 1, size(plate_separations)

       distance_beyond_hs_diameter_multiple = hs_diameter * &
            (plate_separations(ith_separation) - floor(plate_separations(ith_separation)))

       n_extra_discretisation_points = distance_beyond_hs_diameter_multiple / distance_between_z_values

       if(abs(n_extra_discretisation_points - nint(n_extra_discretisation_points)) > 1.0E-6) then
          print *, "parameter.f90: CheckValidityOfPlateSeparations: "
          print *, "abs(n_extra_discretisation_points - nint(n_extra_discretisation_points)) > 1.0E-6"
          print *, "abs(n_extra_discretisation_points - nint(n_extra_discretisation_points)) = ", &
               abs(n_extra_discretisation_points - nint(n_extra_discretisation_points))
          print *, "distance_beyond_hs_diameter_multiple = ", distance_beyond_hs_diameter_multiple
          print *, "n_extra_discretisation_points = ", n_extra_discretisation_points
          print *, "Total distance between plates is not an integer multiple of the distance between"
          print *, "discretisations.  Input error.  Change plate separations...aborting..."
          call abort()
       end if

    end do

  end subroutine CheckValidityOfPlateSeparations

end module parameters

!Contains routines neccesary to calculate surface forces.
module surfaceforces
  use kinds
  use parameters
  use integratezcylindrical
  use lambdas
  use helpers
  use functionalderivatives
  use constructoligomers
  use excessenergyfunctionalparameters
  implicit none
  private

  public :: CalculateGrandPotentialValuePerUnitArea

contains

  !Subroutine that calculates the 
  subroutine CalculateGrandPotentialValuePerUnitArea(ith_plate_separation, grand_potential_per_unit_area, &
       size_of_ns_array, n_plus, n_neutral, n_minus)
    integer, intent(in)   :: ith_plate_separation
    real(dp), intent(out) :: grand_potential_per_unit_area
    integer, intent(in) :: size_of_ns_array
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus

    real(dp), dimension(size_of_ns_array) :: n_s
    real(dp), dimension(size_of_ns_array) :: n_sbar

    real(dp), dimension(size_of_ns_array) :: n_plus_input
    real(dp), dimension(size_of_ns_array) :: n_neutral_input
    real(dp), dimension(size_of_ns_array) :: n_minus_input

    real(dp) :: chemical_potential_term

    real(dp) :: potential_per_unit_area_not_in_bulk
    real(dp) :: potential_per_unit_area_in_bulk

    real(dp) :: F_ideal_chain
    real(dp) :: F_van_der_waals

    real(dp) :: F_hard_sphere
    real(dp) :: F_surface_disp
    real(dp) :: F_surface_electro
    real(dp) :: F_electric_like
    real(dp) :: F_electric_unlike

    F_ideal_chain = 0.0_dp
    F_van_der_waals = 0.0_dp
    F_hard_sphere = 0.0_dp
    F_surface_disp = 0.0_dp
    F_surface_electro = 0.0_dp
    F_electric_like = 0.0_dp
    F_electric_unlike = 0.0_dp


    n_s(:) = n_plus(:) + n_neutral(:) + n_minus(:)
    n_sbar(:) = calculate_n_sbar(n_s(:))

    n_plus_input(:) = n_plus(:)
    n_neutral_input(:) = n_neutral(:)
    n_minus_input(:) = n_minus(:)

    F_surface_disp = integrate_z_cylindrical(n_s * &
         calculate_surface_dispersion_functional_deriv(ith_plate_separation, size(n_s)), unity_function)

    F_ideal_chain = calculate_ideal_chain_term_per_unit_area(size(n_plus_input), n_plus_input, n_neutral_input, n_minus_input, ith_plate_separation)

    F_hard_sphere = calculate_hardsphere_term_per_unit_area(n_s, n_sbar)

    F_van_der_waals = integrate_z_cylindrical(0.5_dp * n_s * calculate_vanderWaals_functional_deriv(n_s), unity_function)

    F_surface_electro = integrate_z_cylindrical(n_plus_input(:) * calculate_surface_electrostatic_functional_deriv(size(n_plus_input), positive_bead_charge), unity_function) +&
         integrate_z_cylindrical(n_minus_input * calculate_surface_electrostatic_functional_deriv(size(n_minus_input), negative_bead_charge), unity_function)


    F_electric_like = 0.5_dp * (&
         integrate_z_cylindrical(n_plus * calculate_electrostatic_like_term_functional_deriv(n_plus, positive_bead_charge, .false.), unity_function) +&
         integrate_z_cylindrical(n_minus * calculate_electrostatic_like_term_functional_deriv(n_minus, negative_bead_charge, .false.), unity_function))

    F_electric_unlike = 0.5_dp * (&
         integrate_z_cylindrical(n_plus * calculate_electrostatic_unlike_term_functional_deriv(n_minus, positive_bead_charge, negative_bead_charge, .false.), unity_function) +&
         integrate_z_cylindrical(n_minus * calculate_electrostatic_unlike_term_functional_deriv(n_plus, positive_bead_charge, negative_bead_charge, .false.), unity_function))

    
    potential_per_unit_area_not_in_bulk = (F_ideal_chain + F_van_der_waals + F_surface_disp + F_hard_sphere + &
         F_surface_electro + F_electric_like + F_electric_unlike)


    chemical_potential_term = calculate_chem_potential_term(n_plus_input, n_neutral_input, n_minus_input, ith_plate_separation)


    grand_potential_per_unit_area = potential_per_unit_area_not_in_bulk - chemical_potential_term

  end subroutine CalculateGrandPotentialValuePerUnitArea


  function calculate_chem_potential_term(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_term

    if(trim(ionic_liquid_name) == "SingleNeutralSpheres") then
       calculate_chem_potential_term = calculate_chem_potential_term_neutral_spheres(n_plus, n_neutral, n_minus, ith_plate_separation)
    else if(trim(ionic_liquid_name) == "NeutralDimers") then
       calculate_chem_potential_term = calculate_chem_potential_term_neutral_dimers(n_plus, n_neutral, n_minus, ith_plate_separation)
    else if(trim(ionic_liquid_name) == "C4MIM_BF4-") then
       calculate_chem_potential_term = calculate_chem_potential_C4MIMBF4(n_plus, n_neutral, n_minus, ith_plate_separation)
    else if(trim(ionic_liquid_name) == "PositiveMinusSpheres") then
       calculate_chem_potential_term = calculate_chem_potential_PositiveMinusSpheres(n_plus, n_neutral, n_minus, ith_plate_separation)
    else if(trim(ionic_liquid_name) == "PositiveNeutralDimerMinusSpheres") then
       calculate_chem_potential_term = calculate_chem_potential_PositiveNeutralDimerMinusSpheres(n_plus, n_neutral, n_minus, ith_plate_separation)
    else if(trim(ionic_liquid_name) == "PositiveNeutralDoubleDimerMinusDimer") then
       calculate_chem_potential_term = calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer(n_plus, n_neutral, n_minus, ith_plate_separation)
    else
       print *, "surfaceforces.f90: CalculateChemicalPotentialTerm:"
       print *, "Unsupported 'ionic_liquid_name' of ", trim(ionic_liquid_name)
    end if

  end function Calculate_Chem_Potential_Term


  !Calculates the value of the ideal chain term per unit area.  Note that all the input arguments are optional.
  !In the case that none are present then calculate the bulk value. (by setting lambda_b - lambda = 0).
  !In the case that they're all present calulculate the value based on the input densities.
  !In any other case, abort.
  function calculate_ideal_chain_term_per_unit_area(n_points, n_plus, n_neutral, n_minus, ith_plate_separation)
    integer, intent(in) :: n_points !Number of points to calculate
    real(dp), dimension(:), intent(in), optional :: n_plus
    real(dp), dimension(:), intent(in), optional :: n_neutral
    real(dp), dimension(:), intent(in), optional :: n_minus
    integer, intent(in), optional :: ith_plate_separation
    real(dp) :: calculate_ideal_chain_term_per_unit_area

    real(dp), dimension(n_points) :: lambda_plus, lambda_neutral, lambda_minus
    real(dp), dimension(n_points) :: n_plus_input, n_neutral_input, n_minus_input

    calculate_ideal_chain_term_per_unit_area = 0.0_dp

    lambda_plus(:) = 0.0_dp
    lambda_neutral(:) = 0.0_dp
    lambda_minus(:) = 0.0_dp

    if(present(n_plus) .and. present(n_neutral) .and. present(n_minus) .and. present(ith_plate_separation)) then

       if((size(n_plus) /= size(n_neutral)) .or. (size(n_plus) /= size(n_minus))) then
          print *, "surfaceforces.f90.f90: calculate_ideal_chain_term_per_unit_area"
          print *, "Size mismatch..."
          print *, "(size(n_plus) /= size(n_neutral)) .or. (size(n_plus) /= size(n_minus))"
          print *, "coding error...aborting..."
          call abort()
       end if

       call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)
       n_plus_input(:) = n_plus(:)
       n_neutral_input(:) = n_neutral(:)
       n_minus_input(:) = n_minus(:)
    else if(present(n_plus) .or. present(n_neutral) .or. present(n_minus) .or. present(ith_plate_separation)) then
       print *, "surfaceforces.f90: calculate_ideal_chain_term_per_unit_area"
       print *, "routine has optional argument list.  They must be all present or all not present."
       print *, "This rule has been violated...coding error...aborting..."
       call abort()
    else !else we want to calculate the value in the bulk so we leave lambda_b - lambda = 0
       n_plus_input(:) = bulk_density_positive_beads
       n_neutral_input(:) = bulk_density_neutral_beads
       n_minus_input(:) = bulk_density_negative_beads
    end if

    if(trim(ionic_liquid_name) == "SingleNeutralSpheres") then
       calculate_ideal_chain_term_per_unit_area = calculate_single_neutral_sphere_ideal_chain_term(n_neutral_input)
    else if(trim(ionic_liquid_name) == "NeutralDimers") then
       calculate_ideal_chain_term_per_unit_area = calculate_neutral_dimers_ideal_chain_term(lambda_neutral)
    else if(trim(ionic_liquid_name) == "C4MIM_BF4-") then
       calculate_ideal_chain_term_per_unit_area = calculate_C4MIMBF4_ideal_chain_term(lambda_plus, lambda_neutral, lambda_minus)
    else if(trim(ionic_liquid_name) == "PositiveMinusSpheres") then
       calculate_ideal_chain_term_per_unit_area = calculate_PositiveMinusSpheres_ideal_chain_term(n_plus_input, n_minus_input)
    else if(trim(ionic_liquid_name) == "PositiveNeutralDimerMinusSpheres") then
       calculate_ideal_chain_term_per_unit_area = calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term(lambda_plus, lambda_neutral, n_minus_input)
    else if(trim(ionic_liquid_name) == "PositiveNeutralDoubleDimerMinusDimer") then
       calculate_ideal_chain_term_per_unit_area = calculate_PositiveNeutralDoubleDimerMinusDimer_ideal_chain_term(lambda_plus, lambda_neutral, lambda_minus)
    else
       print *, "surfaceforces.f90: calculate_ideal_chain_term_per_unit_area:"
       print *, "Unsupported 'ionic_liquid_name' of ", trim(ionic_liquid_name)
    end if

  end function calculate_ideal_chain_term_per_unit_area

  !Calculates the hard sphere contribution to the grand potential.
  !Note that the integration over the angle and radial directions in cylindrical
  !co-ordinates cancel, as we are interested in the value of the term
  !PER UNIT AREA.
  function calculate_hardsphere_term_per_unit_area(n_s, n_sbar)
    real(dp), dimension(:), intent(in) :: n_s
    real(dp), dimension(:), intent(in) :: n_sbar
    real(dp) :: calculate_hardsphere_term_per_unit_area

    integer :: start_z_index
    integer :: end_z_index

    real(dp), dimension(size(n_s)) :: hs_integrand
    real(dp), dimension(size(n_s)) :: n_solvent

    hs_integrand(:) = 0.0_dp
    n_solvent(:) = 0.0_dp

    call get_allowed_z_values(start_z_index, end_z_index, size(n_s))

    hs_integrand(start_z_index:end_z_index) = (0.5_dp / beta) * n_s(start_z_index:end_z_index) * &
         GetAEx(n_sbar(start_z_index:end_z_index), n_solvent(start_z_index:end_z_index), a_term_index)

    calculate_hardsphere_term_per_unit_area = integrate_z_cylindrical(hs_integrand, unity_function)

  end function calculate_hardsphere_term_per_unit_area

end module surfaceforces

!A module containing the values of universal constants.
module universalconstants
  use kinds
  implicit none
  public

  real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp ! pi in double precision.
  real(dp), parameter :: epsilon0=8.854187817620389850536563031710750260608E-22_dp ! epsilon_0 in double precision [Farads/Angstrom] ([Farad] = [C^2]/[J])
  real(dp), parameter :: k_B=1.38064852E-23_dp !([J]/[K])
  real(dp), parameter :: electric_charge=1.6021766E-19_dp !the electric charge ([C])

  complex(sc), parameter :: i_sc = ( 0.0_sp, 1.0_sp ) ! The square root of -1, in single precision.
  complex(dc), parameter :: i_dc = ( 0.0_dp, 1.0_dp ) ! The square root of -1, in double precision.

end module universalconstants
