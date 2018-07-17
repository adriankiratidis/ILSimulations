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

  subroutine UpdateDensities(lambda1, n1_updated_total, n1_updated_end, lambda2, n2_updated_total, n2_updated_end, lambda3, n3_updated_total, n3_updated_end)
    real(dp), dimension(:), intent(in) :: lambda1
    real(dp), dimension(:), intent(out) :: n1_updated_total
    real(dp), dimension(:), intent(out) :: n1_updated_end

    real(dp), dimension(:), intent(in), optional :: lambda2
    real(dp), dimension(:), intent(out), optional :: n2_updated_total
    real(dp), dimension(:), intent(out), optional :: n2_updated_end

    real(dp), dimension(:), intent(in), optional :: lambda3
    real(dp), dimension(:), intent(out), optional :: n3_updated_total
    real(dp), dimension(:), intent(out), optional :: n3_updated_end

    real(dp), dimension(size(n1_updated_end)) :: n1_updated_non_end
    real(dp), dimension(size(n1_updated_end)) :: n2_updated_non_end
    real(dp), dimension(size(n1_updated_end)) :: n3_updated_non_end

    if(trim(ionic_liquid_name) == "SingleNeutralSpheres") then

       if(present(lambda2) .or. present(n2_updated_end) .or. present(lambda3) .or. present(n3_updated_end) .or. &
            present(n2_updated_total) .or. present(n3_updated_total)) then
          print *, "When doing single spheres should not have more than one lambda/density pair present"
          print *, "...aborting..."
          !call abort()
       end if

       call UpdateSingleNeutralSphereDensity(lambda1, n1_updated_end)
       n1_updated_total = n1_updated_end

    else if(trim(ionic_liquid_name) == "SinglePositiveSpheres") then

       if(present(lambda2) .or. present(n2_updated_end) .or. present(lambda3) .or. present(n3_updated_end) .or. &
            present(n2_updated_total) .or. present(n3_updated_total)) then
          print *, "When doing single spheres should not have more than one lambda/density pair present"
          print *, "...aborting..."
          !call abort()
       end if

       call UpdateSinglePositiveSphereDensity(lambda1, n1_updated_end)
       n1_updated_total = n1_updated_end

    else if(trim(ionic_liquid_name) == "SingleNegativeSpheres") then

       if(present(lambda2) .or. present(n2_updated_end) .or. present(lambda3) .or. present(n3_updated_end) .or. &
            present(n2_updated_total) .or. present(n3_updated_total)) then
          print *, "When doing single spheres should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       end if

       call UpdateSingleNegativeSphereDensity(lambda1, n1_updated_end)
       n1_updated_total = n1_updated_end

    else if(trim(ionic_liquid_name) == "SinglePositiveNeutralMinusSpheres") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated_end)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated_end)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateSinglePositiveNeutralMinusSphereDensities(lambda1, n1_updated_end, lambda2, n2_updated_end, lambda3, n3_updated_end)
       end if
       n1_updated_total = n1_updated_end
       n2_updated_total = n2_updated_end
       n3_updated_total = n3_updated_end

    else if(trim(ionic_liquid_name) == "PositiveMinusSpheres") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated_end)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated_end)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdatePositiveMinusSphereDensities(lambda1, n1_updated_end, lambda3, n3_updated_end)
       end if

       n1_updated_total = n1_updated_end
       n3_updated_total = n3_updated_end

    else if(trim(ionic_liquid_name) == "PositiveNeutralDimerMinusSpheres") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated_end)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated_end)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdatePositiveNeutralDimerMinusSphereDensities(lambda1, n1_updated_end, lambda2, n2_updated_end, lambda3, n3_updated_end)
       end if

       n1_updated_total = n1_updated_end 
       n2_updated_total = n2_updated_end 
       n3_updated_total = n3_updated_end 

    else if(trim(ionic_liquid_name) == "PositiveNeutralDoubleDimerMinusDimer") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated_end)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated_end)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdatePositiveNeutralDoubleDimerMinusDimerDensities(lambda1, n1_updated_end, n1_updated_non_end, lambda2, &
               n2_updated_end, n2_updated_non_end, lambda3, n3_updated_end, n3_updated_non_end)
       end if

       n1_updated_total = n1_updated_end + n1_updated_non_end 
       n2_updated_total = n2_updated_end + n2_updated_non_end 
       n3_updated_total = n3_updated_end + n3_updated_non_end 

    else if(trim(ionic_liquid_name) == "NeutralDimers") then

       if(present(lambda2) .or. present(n2_updated_end) .or. present(lambda3) .or. present(n3_updated_end) .or. &
            present(n2_updated_total) .or. present(n3_updated_total)) then
          print *, "When doing Neutral Dimers should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       else
          !call UpdateSingleNeutralSphereDensity(lambda1, n1_updated)
          call UpdateNeutralDimerDensity(lambda1, n1_updated_end)
       end if

       n1_updated_total = n1_updated_end

    else if(trim(ionic_liquid_name) == "C4MIM_BF4-") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated_end)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated_end)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC4MIMPositiveBeadDensities(lambda1, lambda2, n1_updated_end, n1_updated_non_end)
          call UpdateC4MIMNeutralBeadDensities(lambda1, lambda2, n2_updated_end, n2_updated_non_end)
          call UpdateBF4NegativeBeadDensities(lambda3, n3_updated_end, n3_updated_non_end)
       end if

       n1_updated_total = n1_updated_end + n1_updated_non_end
       n2_updated_total = n2_updated_end + n2_updated_non_end
       n3_updated_total = n3_updated_end + n3_updated_non_end

    else

       print *, "constructoligomers.f90: UpdateDensities: "
       print *, "Unsupported ionic_liquid_name value of ", trim(ionic_liquid_name)
       print *, "...aborting..."
       call abort()
    end if

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
    !call RenormaliseToBulkDensity(n_plus_updated, "n+")

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
    !call RenormaliseToBulkDensity(n_neutral_updated, "n0")

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
    !call RenormaliseToBulkDensity(n_minus_updated, "n-")

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

    n_plus_updated = bulk_density * exp(lambda_plus)
    n_minus_updated = bulk_density * exp(lambda_minus)

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

  subroutine UpdatePositiveNeutralDoubleDimerMinusDimerDensities(lambda_plus, n_plus_updated_end, n_plus_updated_non_end, &
       lambda_neutral, n_neutral_updated_end, n_neutral_updated_non_end, lambda_minus, n_minus_updated_end, n_minus_updated_non_end)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated_end
    real(dp), dimension(:), intent(out) :: n_plus_updated_non_end
    
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated_end
    real(dp), dimension(:), intent(out) :: n_neutral_updated_non_end
    
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated_end
    real(dp), dimension(:), intent(out) :: n_minus_updated_non_end

    real(dp), dimension(:), allocatable :: c1, c2, c3, c4

    ! First check the input variables are the same size
    if( (size(lambda_plus) == size(n_plus_updated_end)) .and. (size(lambda_neutral) == size(n_neutral_updated_end)) .and. &
         (size(lambda_minus) == size(n_minus_updated_end)) .and. (size(lambda_plus) == size(n_plus_updated_non_end)) .and. &
         (size(lambda_neutral) == size(n_neutral_updated_non_end)) .and. (size(lambda_minus) == size(n_minus_updated_non_end))) then
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
    n_minus_updated_end = 2.0_dp * bulk_density * exp(lambda_minus) * c1
    n_minus_updated_non_end = 0.0_dp
    
    c4 = integrate_phi_spherical(exp(lambda_neutral))
    c3 = integrate_phi_spherical(exp(lambda_neutral) * c4)
    c2 = integrate_phi_spherical(exp(lambda_plus) * c3)
    c1 = integrate_phi_spherical(exp(lambda_plus))

    !n_plus_updated = bulk_density * ((exp(lambda_plus) * c2) + (exp(lambda_plus) * c3 * c1))
    n_plus_updated_end = bulk_density * (exp(lambda_plus) * c2)
    n_plus_updated_non_end = bulk_density * (exp(lambda_plus) * c3 * c1)
    
    c2 = integrate_phi_spherical(exp(lambda_plus) * c1)
    c3 = integrate_phi_spherical(exp(lambda_neutral) * c2)
    
    !n_neutral_updated = bulk_density * ( (exp(lambda_neutral) * c3) + (exp(lambda_neutral) * c2 * c4) )
    n_neutral_updated_end = bulk_density * (exp(lambda_neutral) * c3)
    n_neutral_updated_non_end = (exp(lambda_neutral) * c2 * c4)

    deallocate(c1, c2, c3, c4)

    call setNonCalculatedRegionToZero(n_plus_updated_end)
    call setNonCalculatedRegionToZero(n_plus_updated_non_end)
    
    call setNonCalculatedRegionToZero(n_neutral_updated_end)
    call setNonCalculatedRegionToZero(n_neutral_updated_non_end)
    
    call setNonCalculatedRegionToZero(n_minus_updated_end)
    call setNonCalculatedRegionToZero(n_minus_updated_non_end)

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

  function calculate_chem_potential_term_neutral_dimers(n_plus_total, n_neutral_total, n_minus_total, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus_total
    real(dp), dimension(:), intent(in) :: n_neutral_total
    real(dp), dimension(:), intent(in) :: n_minus_total
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_term_neutral_dimers

    real(dp), dimension(size(n_plus_total)) :: lambda_plus
    real(dp), dimension(size(n_neutral_total)) :: lambda_neutral
    real(dp), dimension(size(n_minus_total)) :: lambda_minus

    real(dp) :: lambda

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus_total, lambda_neutral, n_neutral_total, lambda_minus, n_minus_total, ith_plate_separation)

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

    calculate_chem_potential_term_neutral_dimers = (1.0_dp/beta) * (log(bulk_density) + (2.0_dp * lambda)) * &
         integrate_z_cylindrical(0.5_dp * n_neutral_total, unity_function)

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

    ! integrand(:) = 4.0_dp * pi * (hs_diameter**2) * integrate_phi_spherical(exp(lambda_neutral))
    ! integrand_with_lambda(:) = 4.0_dp * pi * (hs_diameter**2) * integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

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

    ! integrand(:) = 4.0_dp * pi * (hs_diameter**2) * integrate_phi_spherical(exp(lambda_neutral))
    ! integrand_with_lambda(:) = 4.0_dp * pi * (hs_diameter**2) * integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    integrand(:) = integrate_phi_spherical(exp(lambda_neutral))
    integrand_with_lambda(:) = integrate_phi_spherical(exp(lambda_neutral) * lambda_neutral)

    calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term = calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term + (bulk_density * integrate_z_cylindrical(&
         (integrand_with_lambda(:) + (integrand(:) * (lambda_plus + log(bulk_density) - 1.0_dp))) * exp(lambda_plus), unity_function ) / beta)


  end function calculate_PositiveNeutralDimerMinusSpheres_ideal_chain_term







  function calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer(n_plus_total, n_plus_end, n_neutral_total, n_neutral_end, n_minus_total, n_minus_end, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus_total
    real(dp), dimension(:), intent(in) :: n_plus_end
    real(dp), dimension(:), intent(in) :: n_neutral_total
    real(dp), dimension(:), intent(in) :: n_neutral_end
    real(dp), dimension(:), intent(in) :: n_minus_total
    real(dp), dimension(:), intent(in) :: n_minus_end
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer

    real(dp), dimension(size(n_plus_total)) :: lambda_plus
    real(dp), dimension(size(n_neutral_total)) :: lambda_neutral
    real(dp), dimension(size(n_minus_total)) :: lambda_minus

    real(dp) :: lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus_total, lambda_neutral, n_neutral_total, lambda_minus, n_minus_total, ith_plate_separation)

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

    call CalculateLambdasDifference(lambda_plus, n_plus_total, n_plus_end, lambda_neutral, n_neutral_total, n_neutral_end, lambda_minus, n_minus_total, n_minus_end, ith_plate_separation)

    calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer = (1.0_dp/beta) * (&
         (log(bulk_density) + ((2.0_dp * lambda_plus_bulk) + (2.0_dp * lambda_neutral_bulk))) * &
         (integrate_z_cylindrical((n_plus_total + n_neutral_total)/4.0_dp, unity_function)) + &
         (integrate_z_cylindrical(n_minus_total/2.0_dp, unity_function) * (log(bulk_density) + 2.0_dp * lambda_minus_bulk)))

    ! calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer = (1.0_dp/beta) * (&
    !      (log(bulk_density) + ((2.0_dp * lambda_plus_bulk) + (2.0_dp * lambda_neutral_bulk))) * &
    !      (integrate_z_cylindrical((n_plus_total + n_neutral_total)/4.0_dp, unity_function)))

    ! calculate_chem_potential_PositiveNeutralDoubleDimerMinusDimer = (1.0_dp/beta) * &
    !      (log(bulk_density) + ((2.0_dp * lambda_minus_bulk))) * &
    !      (integrate_z_cylindrical(0.5_dp * n_minus_total, unity_function))


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






  function calculate_chem_potential_C4MIMBF4(n_plus_total, n_plus_end, n_neutral_total, n_neutral_end, n_minus_total, n_minus_end, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus_total
    real(dp), dimension(:), intent(in) :: n_plus_end
    real(dp), dimension(:), intent(in) :: n_neutral_total
    real(dp), dimension(:), intent(in) :: n_neutral_end
    real(dp), dimension(:), intent(in) :: n_minus_total
    real(dp), dimension(:), intent(in) :: n_minus_end
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chem_potential_C4MIMBF4


    real(dp), dimension(size(n_plus_total)) :: lambda_plus
    real(dp), dimension(size(n_neutral_total)) :: lambda_neutral
    real(dp), dimension(size(n_minus_total)) :: lambda_minus

    real(dp) :: lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk

    integer :: start_z_index, end_z_index

    call get_allowed_z_values(start_z_index, end_z_index, size(lambda_neutral))

    call CalculateLambdasBulk(lambda_plus, n_plus_total, lambda_neutral, n_neutral_total, lambda_minus, n_minus_total, ith_plate_separation)

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

    call CalculateLambdasDifference(lambda_plus, n_plus_total, n_plus_end, lambda_neutral, n_neutral_total, n_neutral_end, lambda_minus, n_minus_total, n_minus_end, ith_plate_separation)

    !calculate_chem_potential_C4MIMBF4 = 0.0_dp

    calculate_chem_potential_C4MIMBF4 = (1.0_dp/beta) * (&
         (log(bulk_density) + ((5.0_dp * lambda_plus_bulk) + (5.0_dp * lambda_neutral_bulk))) * &
         (integrate_z_cylindrical((n_plus_total + n_neutral_total)/10.0_dp, unity_function)) + &
         (integrate_z_cylindrical(n_minus_total/5.0_dp, unity_function) * (log(bulk_density) + 5.0_dp * lambda_minus_bulk)))

    ! calculate_chem_potential_C4MIMBF4 = (1.0_dp/beta) * (1.0_dp*log(bulk_density) + &
    !      ((5.0_dp * lambda_plus_bulk) + (5.0_dp * lambda_neutral_bulk))) * &
    !      (integrate_z_cylindrical((n_plus_total + n_neutral_total)/10.0_dp, unity_function))


    ! calculate_chem_potential_C4MIMBF4 = (1.0_dp/beta) * (1.0_dp*log(bulk_density) + &
    !      (5.0_dp * lambda_minus_bulk)) * &
    !      (integrate_z_cylindrical(n_minus_total/5.0_dp, unity_function))


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
    !print *, "calculate_C4MIMBF4_ideal_chain_term = ", calculate_C4MIMBF4_ideal_chain_term
    !call abort()

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

  subroutine UpdateC4MIMPositiveBeadDensities(lambda_plus, lambda_neutral, n_plus_updated_end, n_plus_updated_non_end)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_plus_updated_end
    real(dp), dimension(:), intent(out) :: n_plus_updated_non_end

    real(dp), dimension(:), allocatable :: c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp
    integer :: array_size

    ! First check the input variables are the same size
    if( (size(lambda_plus) == size(lambda_neutral)) .and. (size(lambda_plus) == size(n_plus_updated_end)) .and. &
        (size(lambda_plus) == size(n_plus_updated_non_end))) then

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
    ! gives the required 0.5_dp factor that we are multiplying by.
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
    
    !n_plus_updated = bulk_density * ( (exp(lambda_plus) * c8c1 * c3p) + (exp(lambda_plus) * c2 * c9c10 * c9c10 * c4) + &
    !     (exp(lambda_plus) * c3ppp * c5) + (2.0_dp * (exp(lambda_plus) * c3pp)) )

    !call setNonCalculatedRegionToZero(n_plus_updated)

    n_plus_updated_non_end = bulk_density * ((exp(lambda_plus) * c8c1 * c3p) + (exp(lambda_plus) * c2 * c9c10 * c9c10 * c4) + (exp(lambda_plus) * c3ppp * c5))
    n_plus_updated_end = bulk_density * (2.0_dp * (exp(lambda_plus) * c3pp))

    call setNonCalculatedRegionToZero(n_plus_updated_end)
    call setNonCalculatedRegionToZero(n_plus_updated_non_end)
    
    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp)

  end subroutine UpdateC4MIMPositiveBeadDensities

  subroutine UpdateC4MIMNeutralBeadDensities(lambda_plus, lambda_neutral, n_neutral_updated_end, n_neutral_updated_non_end)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated_end
    real(dp), dimension(:), intent(out) :: n_neutral_updated_non_end

    real(dp), dimension(:), allocatable :: c3p, c3ppp, c9c10, c8c1, c7, c6, c5, c4, c2 !contributions also used in +ve beads.
    real(dp), dimension(:), allocatable :: c2p, c4p, c5p, c6p, c7p !extra contributions for -ve beads.

    integer :: array_size

    ! First check the input variables are the same size
    if( (size(lambda_plus) == size(lambda_neutral)) .and. (size(lambda_plus) == size(n_neutral_updated_end)) .and. &
         (size(lambda_plus) == size(n_neutral_updated_non_end))) then

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
    
    !n_neutral_updated = bulk_density * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
    !     (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8c1) + (exp(lambda_neutral) * c7p) )
    !call setNonCalculatedRegionToZero(n_neutral_updated)

    n_neutral_updated_end = bulk_density * ((exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c7p))
    n_neutral_updated_non_end = bulk_density * ((exp(lambda_neutral) * c4p * c6) + (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8c1))

    call setNonCalculatedRegionToZero(n_neutral_updated_end)
    call setNonCalculatedRegionToZero(n_neutral_updated_non_end)
    
    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(c3p, c3ppp, c8c1)
    deallocate(c2p, c4p, c5p, c6p, c7p)

  end subroutine UpdateC4MIMNeutralBeadDensities

  subroutine UpdateBF4NegativeBeadDensities(lambda_minus, n_minus_updated_end, n_minus_updated_non_end)
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated_end
    real(dp), dimension(:), intent(out) :: n_minus_updated_non_end

    real(dp), dimension(:), allocatable :: a1a2a3a4, a5p
    integer :: array_size

    ! First check the input variables are the same size
    if((size(lambda_minus) == size(n_minus_updated_end)) .and. (size(lambda_minus) == size(n_minus_updated_non_end))) then
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

    !print *, "lambda_minus = ", lambda_minus
    !print *, "exp(lambda_minus) = ", exp(lambda_minus)

    !Calculate the required contributions for the anion
    a1a2a3a4 = integrate_phi_spherical(exp(lambda_minus))

    !print *, "a1a2a3a4 = ", a1a2a3a4 

    a5p = integrate_phi_spherical(exp(lambda_minus) * (a1a2a3a4 ** 3.0_dp))

    !Calculate the resulting negative bead densities.
    !n_minus_updated = na1 + na2 + na3 + na4 + na5 = 4*na1 + na5

    n_minus_updated_end = bulk_density * ( 4.0_dp*(exp(lambda_minus) * a5p))
    n_minus_updated_non_end = bulk_density * (exp(lambda_minus) * (a1a2a3a4**4.0_dp))

    !n_minus_updated = bulk_density * ( 4.0_dp*(exp(lambda_minus) * a5p) + (exp(lambda_minus) * (a1a2a3a4**4.0_dp)) )
    !call setNonCalculatedRegionToZero(n_minus_updated)

    call setNonCalculatedRegionToZero(n_minus_updated_end)
    call setNonCalculatedRegionToZero(n_minus_updated_non_end)
    !n_minus_updated = 0.0_dp

    !print *, "n_minus_updated = ", n_minus_updated

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(a1a2a3a4, a5p)

  end subroutine UpdateBF4NegativeBeadDensities

  !   subroutine UpdateC4MIMPositiveBeadDensities(lambda_plus, lambda_neutral, n_plus_updated)
  !   real(dp), dimension(:), intent(in) :: lambda_plus
  !   real(dp), dimension(:), intent(in) :: lambda_neutral
  !   real(dp), dimension(:), intent(out) :: n_plus_updated

  !   real(dp), dimension(:), allocatable :: c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp
  !   integer :: array_size

  !   ! First check the input variables are the same size
  !   if( (size(lambda_plus) == size(lambda_neutral)) .and. (size(lambda_plus) == size(n_plus_updated))) then

  !      array_size = size(lambda_plus)
  !      allocate(c8c1(array_size), c9c10(array_size), c7(array_size), c6(array_size), &
  !           c5(array_size), c4(array_size), c2(array_size), c3p(array_size), c3pp(array_size), c3ppp(array_size))
  !   else
  !      print *, "constructoligomers.f90: UpdateC4MIMPositiveBeads:"
  !      print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
  !      print *, "or size(lambda_plus) == size(n_plus_updated)"
  !      print *, "...aborting..."
  !      call abort()
  !   end if

  !   ! Note that there is a factor of 1.0_dp/(4.0_dp * pi * hs_diameter**2)
  !   ! from the delta function bond, a factor of 2 * pi from the theta integral
  !   ! and a factor of hs_diameter**2 from the jacobian.  Combining these factors
  !   ! gives the required 0.5_dp factor that we are multiplying by.
  !   c9c10 = integrate_phi_spherical(exp(lambda_plus))

  !   c8c1 = integrate_phi_spherical(exp(lambda_neutral))

  !   c7 = integrate_phi_spherical(exp(lambda_neutral) * c8c1)

  !   c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

  !   c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

  !   c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

  !   c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

  !   c2 = integrate_phi_spherical(exp(lambda_plus) * c8c1)

  !   c3pp = integrate_phi_spherical(exp(lambda_plus) * c4 * c2 * c9c10)

  !   c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

  !   !Calculate the resulting positive bead densities.
  !   !n_plus_updated = nc2 + nc3 + nc4 + nc9 + nc10 =  nc2 + nc3 + nc4 + 2*nc9
  !   n_plus_updated = bulk_density * ( (exp(lambda_plus) * c8c1 * c3p) + (exp(lambda_plus) * c2 * c9c10 * c9c10 * c4) + &
  !        (exp(lambda_plus) * c3ppp * c5) + (2.0_dp * (exp(lambda_plus) * c3pp)) )

  !   !n_plus_updated = 0.0_dp

  !   !call RenormaliseToBulkDensity(n_plus_updated, "n+")
  !   call setNonCalculatedRegionToZero(n_plus_updated)

  !   !Don't check if the variables have been allocated.
  !   !We want an error thrown if they haven't been (which should never happen anyway).
  !   deallocate(c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp)

  ! end subroutine UpdateC4MIMPositiveBeadDensities

  ! subroutine UpdateC4MIMNeutralBeadDensities(lambda_plus, lambda_neutral, n_neutral_updated)
  !   real(dp), dimension(:), intent(in) :: lambda_plus
  !   real(dp), dimension(:), intent(in) :: lambda_neutral
  !   real(dp), dimension(:), intent(out) :: n_neutral_updated

  !   real(dp), dimension(:), allocatable :: c3p, c3ppp, c9c10, c8c1, c7, c6, c5, c4, c2 !contributions also used in +ve beads.
  !   real(dp), dimension(:), allocatable :: c2p, c4p, c5p, c6p, c7p !extra contributions for -ve beads.

  !   integer :: array_size

  !   ! First check the input variables are the same size
  !   if( (size(lambda_plus) == size(lambda_neutral)) .and. (size(lambda_plus) == size(n_neutral_updated))) then

  !      array_size = size(lambda_plus)
  !      allocate(c3p(array_size), c3ppp(array_size), c8c1(array_size), c9c10(array_size), c2p(array_size), &
  !           c7(array_size), c6(array_size), c5(array_size), c4(array_size), c2(array_size), &
  !           c4p(array_size), c5p(array_size), c6p(array_size), c7p(array_size))
  !   else
  !      print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
  !      print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
  !      print *, "or size(lambda_plus) == size(n_neutral_updated)"
  !      print *, "...aborting..."
  !      call abort()
  !   end if

  !   c9c10 = integrate_phi_spherical(exp(lambda_plus))

  !   c8c1 = integrate_phi_spherical(exp(lambda_neutral))

  !   c7 = integrate_phi_spherical(exp(lambda_neutral) * c8c1)

  !   c6 = integrate_phi_spherical(exp(lambda_neutral) * c7)

  !   c5 = integrate_phi_spherical(exp(lambda_neutral) * c6)

  !   c4 = integrate_phi_spherical(exp(lambda_plus) * c5)

  !   c2 = integrate_phi_spherical(exp(lambda_plus) * c8c1)

  !   c3p = integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

  !   c3ppp = integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

  !   c2p = integrate_phi_spherical(exp(lambda_plus) * c3p)

  !   c4p = integrate_phi_spherical(exp(lambda_plus) * c3ppp)

  !   c5p = integrate_phi_spherical(exp(lambda_neutral) * c4p)

  !   c6p = integrate_phi_spherical(exp(lambda_neutral) * c5p)

  !   c7p = integrate_phi_spherical(exp(lambda_neutral) * c6p)

  !   !Calculate the resulting neutral bead densities.
  !   !n_zero_updated = nc1 + nc5 + nc6 + nc7 + nc8
  !   n_neutral_updated = bulk_density * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
  !        (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8c1) + (exp(lambda_neutral) * c7p) )

  !   !n_neutral_updated = 0.0_dp

  !   !call RenormaliseToBulkDensity(n_neutral_updated, "n0")
  !   call setNonCalculatedRegionToZero(n_neutral_updated)

  !   !Don't check if the variables have been allocated.
  !   !We want an error thrown if they haven't been (which should never happen anyway).
  !   deallocate(c3p, c3ppp, c8c1)
  !   deallocate(c2p, c4p, c5p, c6p, c7p)

  ! end subroutine UpdateC4MIMNeutralBeadDensities

  ! subroutine UpdateBF4NegativeBeadDensities(lambda_minus, n_minus_updated)
  !   real(dp), dimension(:), intent(in) :: lambda_minus
  !   real(dp), dimension(:), intent(out) :: n_minus_updated

  !   real(dp), dimension(:), allocatable :: a1a2a3a4, a5p
  !   integer :: array_size

  !   ! First check the input variables are the same size
  !   if(size(lambda_minus) == size(n_minus_updated)) then
  !      array_size = size(lambda_minus)
  !      allocate(a1a2a3a4(array_size))
  !      allocate(a5p(array_size))

  !      a1a2a3a4(:) = 0.0_dp
  !      a5p(:) = 0.0_dp
  !   else
  !      print *, "constructoligomers.f90: UpdateBF4NegativeBeads:"
  !      print *, "Size mismatch. size(lambda_minus) /= size(n_minus_updated)"
  !      print *, "...aborting..."
  !      call abort()
  !   end if

  !   !print *, "lambda_minus = ", lambda_minus
  !   !print *, "exp(lambda_minus) = ", exp(lambda_minus)

  !   !Calculate the required contributions for the anion
  !   a1a2a3a4 = integrate_phi_spherical(exp(lambda_minus))

  !   !print *, "a1a2a3a4 = ", a1a2a3a4 

  !   a5p = integrate_phi_spherical(exp(lambda_minus) * (a1a2a3a4 ** 3.0_dp))

  !   !print *, "a5p = ", a5p

  !   !Calculate the resulting negative bead densities.
  !   !n_minus_updated = na1 + na2 + na3 + na4 + na5 = 4*na1 + na5
  !   n_minus_updated = bulk_density * ( 4.0_dp*(exp(lambda_minus) * a5p) + (exp(lambda_minus) * (a1a2a3a4**4.0_dp)) )

  !   !call RenormaliseToBulkDensity(n_minus_updated, "n-")
  !   call setNonCalculatedRegionToZero(n_minus_updated)

  !   !n_minus_updated = 0.0_dp

  !   !print *, "n_minus_updated = ", n_minus_updated

  !   !Don't check if the variables have been allocated.
  !   !We want an error thrown if they haven't been (which should never happen anyway).
  !   deallocate(a1a2a3a4, a5p)

  ! end subroutine UpdateBF4NegativeBeadDensities
  
end module constructoligomers
