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
  
  subroutine CalculateLambdas(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, lambda_hs_end_cation, n_hs_end_cation, lambda_hs_nonend_cation, n_hs_nonend_cation, lambda_hs_end_anion, n_hs_end_anion, lambda_hs_nonend_anion, n_hs_nonend_anion, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus
    real(dp), dimension(:), intent(in)  :: n_plus

    real(dp), dimension(:), intent(out) :: lambda_neutral
    real(dp), dimension(:), intent(in)  :: n_neutral

    real(dp), dimension(:), intent(out) :: lambda_minus
    real(dp), dimension(:), intent(in)  :: n_minus

    !Variables for hs calculation.  End and non-end monomer density needed.
    real(dp), dimension(:), intent(out) :: lambda_hs_end_cation
    real(dp), dimension(:), intent(in)  :: n_hs_end_cation
    real(dp), dimension(:), intent(out) :: lambda_hs_nonend_cation
    real(dp), dimension(:), intent(in)  :: n_hs_nonend_cation

    real(dp), dimension(:), intent(out) :: lambda_hs_end_anion
    real(dp), dimension(:), intent(in)  :: n_hs_end_anion
    real(dp), dimension(:), intent(out) :: lambda_hs_nonend_anion
    real(dp), dimension(:), intent(in)  :: n_hs_nonend_anion

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


    ! Now calculate the terms in common to all the lambdas.
    lambda_common_terms = CalculateLambdaCommonTerms(n_plus, n_neutral, n_minus, ith_plate_separation, .false.)

    ! Now calculate our lambdas(r)
    lambda_plus = beta * (lambda_common_terms + CalculateLambdaPlusSpecificTerms(n_plus, n_minus, .false.))
    lambda_neutral = beta * (lambda_common_terms + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .false.))
    lambda_minus = beta * (lambda_common_terms + CalculateLambdaMinusSpecificTerms(n_plus, n_minus, .false.))

    ! lambda_hs_end and lambda_hs_nonend terms correspond to the contributions from the hard sphere term (only) from 
    ! end and nonend monomers respectively.
    call CalculateHSEndAndNonEndHSFunctionalDeriv(n_plus + n_neutral + n_minus, lambda_hs_end_cation, n_hs_end_cation, lambda_hs_nonend_cation, n_hs_nonend_cation, &
         lambda_hs_end_anion, n_hs_end_anion, lambda_hs_nonend_anion, n_hs_nonend_anion, .false.)

    lambda_hs_end_cation = beta * lambda_hs_end_cation
    lambda_hs_nonend_cation = beta * lambda_hs_nonend_cation
    lambda_hs_end_anion = beta * lambda_hs_end_anion
    lambda_hs_nonend_anion = beta * lambda_hs_nonend_anion

    !print *, "TESTING non bulk", lambda_hs_end_cation
    !print *, "n_plus", n_plus
    !print *, "n_neutral", n_neutral
    !print *, "n_minus", n_minus

  end subroutine CalculateLambdas

  subroutine CalculateLambdasBulk(lambda_plus_bulk, n_plus, lambda_neutral_bulk, n_neutral, lambda_minus_bulk, n_minus, lambda_hs_end_cation_bulk, lambda_hs_nonend_cation_bulk, &
       lambda_hs_end_anion_bulk, lambda_hs_nonend_anion_bulk, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus_bulk
    real(dp), dimension(:), intent(in)  :: n_plus

    real(dp), dimension(:), intent(out) :: lambda_neutral_bulk
    real(dp), dimension(:), intent(in)  :: n_neutral

    real(dp), dimension(:), intent(out) :: lambda_minus_bulk
    real(dp), dimension(:), intent(in)  :: n_minus

    !Variables for hs calculation.  End and non-end monomer density needed.
    real(dp), dimension(:), intent(out) :: lambda_hs_end_cation_bulk
    real(dp), dimension(:), intent(out) :: lambda_hs_nonend_cation_bulk
    real(dp), dimension(:), intent(out) :: lambda_hs_end_anion_bulk
    real(dp), dimension(:), intent(out) :: lambda_hs_nonend_anion_bulk

    integer, intent(in)                 :: ith_plate_separation

    real(dp), dimension(size(n_plus)) :: n_hs_end_cation_bulk
    real(dp), dimension(size(n_plus)) :: n_hs_nonend_cation_bulk
    real(dp), dimension(size(n_plus)) :: n_hs_end_anion_bulk
    real(dp), dimension(size(n_plus)) :: n_hs_nonend_anion_bulk

    real(dp), dimension(size(lambda_plus_bulk)) :: lambda_common_terms_bulk
    real(dp), dimension(size(n_neutral)) :: n_neutral_array

    real(dp), dimension(size(n_plus)) :: n_s_bulk

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

    ! lambda_hs_end and lambda_hs_nonend terms correspond to the contributions from the hard sphere term (only) from 
    ! end and nonend monomers respectively.
    n_s_bulk(:) = bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads
    n_hs_end_cation_bulk = n_end_cation_bulk
    n_hs_nonend_cation_bulk = n_nonend_cation_bulk
    n_hs_end_anion_bulk = n_end_anion_bulk
    n_hs_nonend_anion_bulk = n_nonend_anion_bulk

    call CalculateHSEndAndNonEndHSFunctionalDeriv(n_s_bulk, lambda_hs_end_cation_bulk, n_hs_end_cation_bulk, lambda_hs_nonend_cation_bulk, n_hs_nonend_cation_bulk, &
         lambda_hs_end_anion_bulk, n_hs_end_anion_bulk, lambda_hs_nonend_anion_bulk, n_hs_nonend_anion_bulk,.true.)

    lambda_hs_end_cation_bulk = beta * lambda_hs_end_cation_bulk
    lambda_hs_nonend_cation_bulk = beta * lambda_hs_nonend_cation_bulk
    lambda_hs_end_anion_bulk = beta * lambda_hs_end_anion_bulk
    lambda_hs_nonend_anion_bulk = beta * lambda_hs_nonend_anion_bulk

    !print *, "TESTING bulk", lambda_hs_end_cation_bulk, bulk_density_positive_beads, bulk_density_neutral_beads, bulk_density_negative_beads
    

  end subroutine CalculateLambdasBulk

  subroutine CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, lambda_hs_end_cation, n_hs_end_cation, &
       lambda_hs_nonend_cation, n_hs_nonend_cation, lambda_hs_end_anion, n_hs_end_anion, lambda_hs_nonend_anion, n_hs_nonend_anion, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus
    real(dp), dimension(:), intent(in)  :: n_plus

    real(dp), dimension(:), intent(out) :: lambda_neutral
    real(dp), dimension(:), intent(in)  :: n_neutral

    real(dp), dimension(:), intent(out) :: lambda_minus
    real(dp), dimension(:), intent(in)  :: n_minus

    !Variables for hs calculation.  End and non-end monomer density needed.
    real(dp), dimension(:), intent(out) :: lambda_hs_end_cation
    real(dp), dimension(:), intent(in)  :: n_hs_end_cation
    real(dp), dimension(:), intent(out) :: lambda_hs_nonend_cation
    real(dp), dimension(:), intent(in)  :: n_hs_nonend_cation

    real(dp), dimension(:), intent(out) :: lambda_hs_end_anion
    real(dp), dimension(:), intent(in)  :: n_hs_end_anion
    real(dp), dimension(:), intent(out) :: lambda_hs_nonend_anion
    real(dp), dimension(:), intent(in)  :: n_hs_nonend_anion

    integer, intent(in)                 :: ith_plate_separation

    real(dp), dimension(size(lambda_plus)) :: lambda_plus_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_neutral_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_minus_bulk

    real(dp), dimension(size(lambda_plus)) :: lambda_hs_end_cation_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_hs_nonend_cation_bulk

    real(dp), dimension(size(lambda_plus)) :: lambda_hs_end_anion_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_hs_nonend_anion_bulk

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

    call CalculateLambdasBulk(lambda_plus_bulk, n_plus, lambda_neutral_bulk, n_neutral, lambda_minus_bulk, n_minus, lambda_hs_end_cation_bulk, lambda_hs_nonend_cation_bulk, &
         lambda_hs_end_anion_bulk, lambda_hs_nonend_anion_bulk, ith_plate_separation)
    call CalculateLambdas(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, lambda_hs_end_cation, n_hs_end_cation, lambda_hs_nonend_cation, n_hs_nonend_cation, &
         lambda_hs_end_anion, n_hs_end_anion, lambda_hs_nonend_anion, n_hs_nonend_anion, ith_plate_separation)

    lambda_plus(:) = lambda_plus_bulk(:) - lambda_plus(:)
    lambda_neutral(:) = lambda_neutral_bulk(:) - lambda_neutral(:)
    lambda_minus(:) = lambda_minus_bulk(:) - lambda_minus(:)

    lambda_hs_end_cation(:) = lambda_hs_end_cation_bulk(:) - lambda_hs_end_cation(:)
    lambda_hs_nonend_cation(:) = lambda_hs_nonend_cation_bulk(:) - lambda_hs_nonend_cation(:)
    
    lambda_hs_end_anion(:) = lambda_hs_end_anion_bulk(:) - lambda_hs_end_anion(:)
    lambda_hs_nonend_anion(:) = lambda_hs_nonend_anion_bulk(:) - lambda_hs_nonend_anion(:)

    ! lambda_hs_end_cation(:) = 0.0_dp
    ! lambda_hs_nonend_cation(:) = 0.0_dp
    ! lambda_hs_end_anion(:) = 0.0_dp
    ! lambda_hs_nonend_anion(:) = 0.0_dp
    
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

    !print *, "lambda_common_terms = ", lambda_common_terms
    !print *, "lambda_common_terms_bulk = ", lambda_common_terms_bulk
    !call abort()

    ! Now calculate our lambdas(r)
    lambda_plus = beta * (lambda_common_terms + CalculateLambdaPlusSpecificTerms(n_plus, n_minus, .false.))
    lambda_neutral = beta * (lambda_common_terms + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .false.))
    lambda_minus = beta * (lambda_common_terms + CalculateLambdaMinusSpecificTerms(n_plus, n_minus, .false.))

    lambda_plus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaPlusSpecificTerms(n_plus, n_minus, .true.))
    lambda_neutral_bulk = beta * (lambda_common_terms_bulk + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .true.))
    lambda_minus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaMinusSpecificTerms(n_plus, n_minus, .true.))

    !print *, "lambda_neutral_bulk = ", lambda_neutral_bulk
    !print *, "lambda_neutral = ", lambda_neutral

    !Now calculate the bulk value, lambda^{bulk} and in order to return e^{lambda^{bulk} - lambda(r)}
    lambda_plus = lambda_plus_bulk - lambda_plus
    lambda_neutral = lambda_neutral_bulk - lambda_neutral
    lambda_minus = lambda_minus_bulk - lambda_minus

    !print *, "lambda_neutral diff = ", lambda_neutral

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

    real(dp) :: allowed_distance_between_plates

    hs_term = 0.0_dp
    van_der_waals_term = 0.0_dp
    surface_fluid_dispersion_term = 0.0_dp

    if(calculate_bulk) then
       n_s(:) = bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads

       allowed_distance_between_plates = (((size(n_s) - 1)*hs_diameter/n_discretised_points_z) + ((hs_diameter)*1.0_dp))/2.0_dp
       
       !hs_term = 1.0_dp * calculate_hardsphere_functional_deriv(n_s, .true.)
       !print *, "he_term_in_bulk", beta * calculate_hardsphere_functional_deriv(n_s, .true.)
       
       van_der_waals_term = (-8.0_dp * epsilon_LJ_particle_particle * (hs_diameter**6) * pi * n_s(:)) * ( &
            (2.0_dp/(3.0_dp*(hs_diameter**3))))


       call setNonCalculatedRegionToZero(n_s)

    else
       n_s = n_plus + n_neutral + n_minus
       !hs_term = 1.0_dp * calculate_hardsphere_functional_deriv(n_s, .false.)
       !print *, "he_term_not_in_bulk", beta * calculate_hardsphere_functional_deriv(n_s, .false.)

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

       !print *, "n_plus_in = ", n_plus_in
       !print *, "n_minus_in = ", n_minus_in

       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_plus_in, positive_bead_charge, .true.)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_minus_in, positive_bead_charge, negative_bead_charge, .true.)

       !print *, "like_electrostatic_term = ", like_electrostatic_term
       !print *, "unlike_electrostatic_term = ", unlike_electrostatic_term
       !call abort()
    else
       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_plus, positive_bead_charge, .false.)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_minus, positive_bead_charge, negative_bead_charge, .false.)

       surface_electrostatic_term = calculate_surface_electrostatic_functional_deriv(size(n_plus), positive_bead_charge)


       !print *, "surface_electrostatic_term = ",  surface_electrostatic_term(51), surface_electrostatic_term(size(surface_electrostatic_term) - 50), &
       !     surface_electrostatic_term(size(surface_electrostatic_term) - 50) - surface_electrostatic_term(51)
       !print *, "like_electrostatic_term plus = ",  like_electrostatic_term(51), like_electrostatic_term(size(like_electrostatic_term) - 50), &
       !     like_electrostatic_term(size(like_electrostatic_term) - 50) - like_electrostatic_term(51)
       !print *, "unlike_electrostatic_term plus = ", unlike_electrostatic_term(51), unlike_electrostatic_term(size(unlike_electrostatic_term) - 50), &
       !     unlike_electrostatic_term(size(unlike_electrostatic_term) - 50) - unlike_electrostatic_term(51)


       !print *, ""
       !print *, "sum =",  (surface_electrostatic_term(26) + like_electrostatic_term(26) + &
       !     unlike_electrostatic_term(26)) * beta
       !print *, "surface_electrostatic_term = ", surface_electrostatic_term
       !call abort()

    end if



    CalculateLambdaPlusSpecificTerms = surface_electrostatic_term + like_electrostatic_term + unlike_electrostatic_term
    !if(calculate_bulk .eqv. .false.) then
    !   print *, "CalculateLambdaPlusSpecificTerms = ", CalculateLambdaPlusSpecificTerms
    !   call abort()
    !end if

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

       !print *, "surface_electrostatic_term = ",  surface_electrostatic_term(26), surface_electrostatic_term(size(surface_electrostatic_term) - 25), &
       !     surface_electrostatic_term(size(surface_electrostatic_term) - 25) - surface_electrostatic_term(26)
       !print *, "like_electrostatic_term minus = ",  like_electrostatic_term(26), like_electrostatic_term(size(like_electrostatic_term) - 25), &
       !     like_electrostatic_term(size(like_electrostatic_term) - 25) - like_electrostatic_term(26)
       !print *, "unlike_electrostatic_term minus = ", unlike_electrostatic_term(26), unlike_electrostatic_term(size(unlike_electrostatic_term) - 25), &
       !     unlike_electrostatic_term(size(unlike_electrostatic_term) - 25) - unlike_electrostatic_term(26)
       !print *, "sum = ", (surface_electrostatic_term(26) + like_electrostatic_term(26) + unlike_electrostatic_term(26))*beta
       !print *, ""
       
       !call abort()


    end if



    ! print *, "minus terms"
    ! print *, "calculate_bulk = ", calculate_bulk
    ! print *, "surface_electrostatic_term(100) = ", surface_electrostatic_term(100)
    ! print *, "like_electrostatic_term(100) = ", like_electrostatic_term(100)
    ! print *, "unlike_electrostatic_term(100) = ", unlike_electrostatic_term(100)
    !call abort()

    CalculateLambdaMinusSpecificTerms = surface_electrostatic_term + like_electrostatic_term + unlike_electrostatic_term

  end function CalculateLambdaMinusSpecificTerms

end module lambdas
