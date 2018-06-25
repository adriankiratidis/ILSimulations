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
  public :: calculate_chemical_potential_term_neutral_spheres
  public :: calculate_chemical_potential_term_neutral_dimers
  
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

    else if(trim(ionic_liquid_name) == "NeutralDimers") then

       if(present(lambda2) .or. present(n2_updated) .or. present(lambda3) .or. present(n3_updated)) then
          print *, "When doing Neutral Dimers should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       else
          !call UpdateSingleNeutralSphereDensity(lambda1, n1_updated)
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

  
  function calculate_chemical_potential_term_neutral_spheres(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chemical_potential_term_neutral_spheres

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
       print *, "constructoligomers.f90: calculate_chemical_potential_term_neutral_spheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    calculate_chemical_potential_term_neutral_spheres = (1.0_dp/beta) * (log(bulk_density_neutral_beads) + lambda) * &
         integrate_z_cylindrical(n_neutral, unity_function)

  end function calculate_chemical_potential_term_neutral_spheres

  
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

    !integer, intent(in) :: ith_plate_separation
    ! real(dp), dimension(size(n_neutral)) :: zero_array1, zero_array2
    ! real(dp), dimension(size(n_neutral)) :: dummy_array1, dummy_array2
    ! real(dp), dimension(size(n_neutral)) :: lambda_neutral

    ! call SetToZero(zero_array1, zero_array2)
    ! call CalculateLambdas(dummy_array1, zero_array1, lambda_neutral, n_neutral, dummy_array2, zero_array2, ith_plate_separation)

    ! calculate_single_neutral_sphere_ideal_chain_term = integrate_z_cylindrical(bulk_density_neutral_beads * exp(lambda_neutral) * &
    !      ((log(bulk_density_neutral_beads) + lambda_neutral) - 1.0_dp), unity_function) / beta

  end function calculate_single_neutral_sphere_ideal_chain_term

  function calculate_chemical_potential_term_neutral_dimers(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chemical_potential_term_neutral_dimers

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
       print *, "constructoligomers.f90: calculate_chemical_potential_term_neutral_spheres: "
       print *, "When calculating lambda bulk all the values of lambda should be the same"
       print *, "but they aren't, they are...(printed above)...aborting"
       call abort()
    end if

    call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)

    calculate_chemical_potential_term_neutral_dimers = (1.0_dp/beta) * (log(bulk_density) + (2.0_dp * lambda)) * bulk_density * &
        integrate_z_cylindrical(integrate_phi_spherical(exp(lambda_neutral)) * exp(lambda_neutral), unity_function)

    ! calculate_chemical_potential_term_neutral_dimers = (1.0_dp/beta) * (log(bulk_density) + (2.0_dp * lambda)) * &
    !      integrate_z_cylindrical(0.5_dp * n_neutral, unity_function)

  end function calculate_chemical_potential_term_neutral_dimers

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
    ! gives the required 0.5_dp factor that we are multiplying by.
    c9c10 = 0.5_dp * integrate_phi_spherical(exp(lambda_plus))

    c8c1 = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral))

    c7 = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral) * c8c1)

    c6 = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = 0.5_dp * integrate_phi_spherical(exp(lambda_plus) * c5)

    c3p = 0.5_dp * integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c2 = 0.5_dp * integrate_phi_spherical(exp(lambda_plus) * c8c1)

    c3pp = 0.5_dp * integrate_phi_spherical(exp(lambda_plus) * c4 * c2 * c9c10)

    c3ppp = 0.5_dp * integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    !Calculate the resulting positive bead densities.
    !n_plus_updated = nc2 + nc3 + nc4 + nc9 + nc10 =  nc2 + nc3 + nc4 + 2*nc9
    n_plus_updated = bulk_density * ( (exp(lambda_plus) * c8c1 * c3p) + (exp(lambda_plus) * c2 * c9c10 * c9c10 * c4) + &
         (exp(lambda_plus) * c3ppp * c5) + (2.0_dp * (exp(lambda_plus) * c3pp)) )

    !call RenormaliseToBulkDensity(n_plus_updated, "n+")
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

    c9c10 = 0.5_dp * integrate_phi_spherical(exp(lambda_plus))

    c8c1 = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral))

    c7 = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral) * c8c1)

    c6 = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral) * c7)

    c5 = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral) * c6)

    c4 = 0.5_dp * integrate_phi_spherical(exp(lambda_plus) * c5)

    c2 = 0.5_dp * integrate_phi_spherical(exp(lambda_plus) * c8c1)

    c3p = 0.5_dp * integrate_phi_spherical(exp(lambda_plus) * c4 * c9c10 * c9c10)

    c3ppp = 0.5_dp * integrate_phi_spherical(exp(lambda_plus) * c2 * c9c10 * c9c10)

    c2p = 0.5_dp * integrate_phi_spherical(exp(lambda_plus) * c3p)

    c4p = 0.5_dp * integrate_phi_spherical(exp(lambda_plus) * c3ppp)

    c5p = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral) * c4p)

    c6p = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral) * c5p)

    c7p = 0.5_dp * integrate_phi_spherical(exp(lambda_neutral) * c6p)

    !Calculate the resulting neutral bead densities.
    !n_zero_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * ( (exp(lambda_neutral) * c2p) + (exp(lambda_neutral) * c4p * c6) + &
         (exp(lambda_neutral) * c5p * c7) + (exp(lambda_neutral) * c6p * c8c1) + (exp(lambda_neutral) * c7p) )

    !call RenormaliseToBulkDensity(n_neutral_updated, "n0")
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

    print *, "lambda_minus = ", lambda_minus
    print *, "exp(lambda_minus) = ", exp(lambda_minus)

    !Calculate the required contributions for the anion
    a1a2a3a4 = 0.5_dp * integrate_phi_spherical(exp(lambda_minus))

    print *, "a1a2a3a4 = ", a1a2a3a4 

    a5p = 0.5_dp * integrate_phi_spherical(exp(lambda_minus) * (a1a2a3a4 ** 3.0_dp))

    print *, "a5p = ", a5p

    !Calculate the resulting negative bead densities.
    !n_minus_updated = na1 + na2 + na3 + na4 + na5 = 4*na1 + na5
    n_minus_updated = bulk_density * ( 4.0_dp*(exp(lambda_minus) * a5p) + (exp(lambda_minus) * (a1a2a3a4**4.0_dp)) )

    !call RenormaliseToBulkDensity(n_minus_updated, "n-")
    call setNonCalculatedRegionToZero(n_minus_updated)

    print *, "n_minus_updated = ", n_minus_updated

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(a1a2a3a4, a5p)

  end subroutine UpdateBF4NegativeBeadDensities

end module constructoligomers
