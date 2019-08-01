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
  !public :: CalculateLambdasDifference_old

contains

  subroutine CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, &
       lambda_cation_centre, n_cation_centre, lambda_anion_centre, n_anion_centre, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus
    real(dp), dimension(:), intent(in)  :: n_plus

    real(dp), dimension(:), intent(out) :: lambda_neutral
    real(dp), dimension(:), intent(in)  :: n_neutral

    real(dp), dimension(:), intent(out) :: lambda_minus
    real(dp), dimension(:), intent(in)  :: n_minus

    real(dp), dimension(:), intent(out) :: lambda_cation_centre
    real(dp), dimension(:), intent(in)  :: n_cation_centre

    real(dp), dimension(:), intent(out) :: lambda_anion_centre
    real(dp), dimension(:), intent(in)  :: n_anion_centre

    integer, intent(in)                 :: ith_plate_separation

    real(dp), dimension(size(lambda_plus)) :: lambda_plus_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_neutral_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_minus_bulk

    real(dp), dimension(size(lambda_plus)) :: lambda_cation_centre_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_anion_centre_bulk

    integer :: input_array_size

    lambda_cation_centre = 0.0_dp
    lambda_anion_centre = 0.0_dp

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

    call CalculateLambdasBulk(lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk, lambda_cation_centre_bulk, lambda_anion_centre_bulk, ith_plate_separation)
    call CalculateLambdas(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, &
         lambda_cation_centre, n_cation_centre, lambda_anion_centre, n_anion_centre, ith_plate_separation)

    !print *, "lambda_plus_bulk = ", lambda_plus_bulk
    !print *, "lambda_minus_bulk = ", lambda_minus_bulk(100)
    !print *, "lambda_plus = ", lambda_plus
    !print *, "lambda_minus = ", lambda_minus()
    !call abort()

    lambda_plus(:) = lambda_plus_bulk(:) - lambda_plus(:)
    lambda_neutral(:) = lambda_neutral_bulk(:) - lambda_neutral(:)
    lambda_minus(:) = lambda_minus_bulk(:) - lambda_minus(:)

    lambda_cation_centre(:) = lambda_cation_centre_bulk(:) - lambda_cation_centre(:)
    lambda_anion_centre(:) = lambda_anion_centre_bulk(:) - lambda_anion_centre(:)

    !print *, "lambda_plus after= ", lambda_plus
    !print *, "lambda_minus after = ", lambda_minus(100)
    !call abort()

    !print *, "lambda_plus = ", lambda_plus
    !print *, "lambda_neutral = ", lambda_neutral
    !print *, "lambda_minus = ", lambda_minus
    !print *, "lambda_cation_centre = ", lambda_cation_centre
    !print *, "lambda_anion_centre = ", lambda_anion_centre
    
  end subroutine CalculateLambdasDifference


  subroutine CalculateLambdas(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, &
       lambda_cation_centre, n_cation_centre, lambda_anion_centre, n_anion_centre, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus
    real(dp), dimension(:), intent(in)  :: n_plus

    real(dp), dimension(:), intent(out) :: lambda_neutral
    real(dp), dimension(:), intent(in)  :: n_neutral

    real(dp), dimension(:), intent(out) :: lambda_minus
    real(dp), dimension(:), intent(in)  :: n_minus

    real(dp), dimension(:), intent(out) :: lambda_cation_centre
    real(dp), dimension(:), intent(in) :: n_cation_centre

    real(dp), dimension(:), intent(out) :: lambda_anion_centre
    real(dp), dimension(:), intent(in) :: n_anion_centre

    integer, intent(in)                 :: ith_plate_separation

    real(dp), dimension(size(lambda_plus)) :: lambda_common_terms
    real(dp), dimension(size(n_plus)) :: n_s
    
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

    n_s(:) = n_plus(:) + n_neutral(:) + n_minus(:)

    ! Now calculate the terms in common to all the lambdas
    lambda_common_terms = CalculateLambdaCommonTerms(.false., size(n_plus), n_s, ith_plate_separation)

    ! Now calculate our lambdas(r)
    lambda_plus = beta * (lambda_common_terms + CalculateLambdaPlusSpecificTerms(.false., size(lambda_plus), n_plus, n_minus))
    lambda_neutral = beta * (lambda_common_terms + CalculateLambaNeutralSpecificTerms(.false., size(lambda_neutral), n_plus, n_neutral, n_minus))
    lambda_minus = beta * (lambda_common_terms + CalculateLambdaMinusSpecificTerms(.false., size(lambda_minus), n_plus, n_minus))

    lambda_cation_centre = beta * CalculateLambdaCentreToCentreCorrection(.false., size(lambda_cation_centre), 'c', n_cation_centre, n_anion_centre)
    lambda_anion_centre = beta * CalculateLambdaCentreToCentreCorrection(.false., size(lambda_anion_centre), 'a', n_cation_centre, n_anion_centre)

    !print *, "lambda_plus = ", lambda_plus
    !print *, "lambda_neutral = ", lambda_neutral
    !print *, "lambda_minus = ", lambda_minus
    !print *, "lambda_cation_centre = ", lambda_cation_centre
    !print *, "lambda_anion_centre = ", lambda_anion_centre
    !call abort()

  end subroutine CalculateLambdas

  subroutine CalculateLambdasBulk(lambda_plus_bulk, lambda_neutral_bulk, lambda_minus_bulk, lambda_cation_centre_bulk, lambda_anion_centre_bulk, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus_bulk
    real(dp), dimension(:), intent(out) :: lambda_neutral_bulk
    real(dp), dimension(:), intent(out) :: lambda_minus_bulk
    real(dp), dimension(:), intent(out) :: lambda_cation_centre_bulk
    real(dp), dimension(:), intent(out) :: lambda_anion_centre_bulk
    integer, intent(in)                 :: ith_plate_separation

    real(dp), dimension(size(lambda_plus_bulk)) :: lambda_common_terms_bulk
    real(dp), dimension(size(lambda_plus_bulk)) :: n_neutral_array
    
    integer :: input_array_size

    ! First ensure that the sizes of all the input variables are the same.
    input_array_size = size(lambda_plus_bulk)
    if((size(lambda_neutral_bulk) /= input_array_size) .or. &
         (size(lambda_minus_bulk) /= input_array_size) .or. &
         (size(lambda_cation_centre_bulk) /= input_array_size) .or. &
         (size(lambda_anion_centre_bulk) /= input_array_size)) then

       print *, "lambda.f90: CalculateLambdasBulk:"
       print *, "Input arrays must all be the same size"
       print *, "Almost certainly coding error as opposed to input error...aborting..."
       call abort()
    end if

    lambda_plus_bulk(:) = 0.0_dp
    lambda_neutral_bulk(:) = 0.0_dp
    lambda_minus_bulk(:) = 0.0_dp
    lambda_cation_centre_bulk(:) = 0.0_dp
    lambda_anion_centre_bulk(:) = 0.0_dp
    
    lambda_common_terms_bulk = CalculateLambdaCommonTerms(.true., size(lambda_plus_bulk))

    lambda_plus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaPlusSpecificTerms(.true., size(lambda_plus_bulk)))
    lambda_neutral_bulk = beta * (lambda_common_terms_bulk + CalculateLambaNeutralSpecificTerms(.true., size(lambda_neutral_bulk)))
    lambda_minus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaMinusSpecificTerms(.true., size(lambda_minus_bulk)))

    lambda_cation_centre_bulk = beta * CalculateLambdaCentreToCentreCorrection(.true., size(lambda_cation_centre_bulk), 'c')
    lambda_anion_centre_bulk = beta * CalculateLambdaCentreToCentreCorrection(.true., size(lambda_anion_centre_bulk), 'a')
    
  end subroutine CalculateLambdasBulk

  !This routine is the correction for the centre-atom to centre-atom potential.
  !It is written such that the lambda can be added to lambda_plus (for example) for the
  !positive bead on the cation and it will add on the different potential and subtract off the bit that
  !isn't applicable for that bead, but is for all other beads.  This improves code readability.
  function CalculateLambdaCentreToCentreCorrection(calculate_bulk, size_of_terms, cation_or_anion, n_cation_centre, n_anion_centre)
    logical, intent(in) :: calculate_bulk
    integer, intent(in) :: size_of_terms
    character(len=*), intent(in) :: cation_or_anion
    real(dp), dimension(:), intent(in), optional :: n_cation_centre
    real(dp), dimension(:), intent(in), optional :: n_anion_centre

    real(dp), dimension(size_of_terms) :: CalculateLambdaCentreToCentreCorrection

    real(dp), dimension(size_of_terms) :: hs_term
    real(dp), dimension(size_of_terms) :: van_der_waals_term
    
    if(calculate_bulk) then
       if(present(n_cation_centre) .or. present(n_cation_centre)) then
          print *, "lambdas.f90:CalculateLambdaCationCentreCorrection:"
          print *, "if we're calculating the value in the bulk then NO need to included densities."
          call abort()
       end if

       CalculateLambdaCentreToCentreCorrection = &!-1.0_dp * CalculateLambdaCommonTerms(calculate_bulk, size_of_terms) + &
            calculate_centre_to_centre_functional_deriv(calculate_bulk, size_of_terms, trim(cation_or_anion), type_of_interaction_r4r8)

    else
       if((.not. present(n_cation_centre)) .or. (.not. present(n_cation_centre))) then
          print *, "lambdas.f90:CalculateLambdaCationCentreCorrection:"
          print *, "if we're calculating the value NOT in the bulk then need to included densities."
          call abort()
       end if

       !hs_term = 1.0_dp * calculate_hardsphere_functional_deriv(n_s, .false.)
       !van_der_waals_term = calculate_vanderWaals_functional_deriv(n_s)

       CalculateLambdaCentreToCentreCorrection = &!-1.0_dp * (hs_term + van_der_waals_term) + &
            calculate_centre_to_centre_functional_deriv(calculate_bulk, size_of_terms, trim(cation_or_anion), type_of_interaction_r4r8, n_cation_centre, n_anion_centre)

    end if

    !CalculateLambdaCentreToCentreCorrection = 0.0_dp
  end function CalculateLambdaCentreToCentreCorrection
  
  ! subroutine CalculateLambdasDifference_old(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)
  !   real(dp), dimension(:), intent(out) :: lambda_plus
  !   real(dp), dimension(:), intent(in)  :: n_plus

  !   real(dp), dimension(:), intent(out) :: lambda_neutral
  !   real(dp), dimension(:), intent(in)  :: n_neutral

  !   real(dp), dimension(:), intent(out) :: lambda_minus
  !   real(dp), dimension(:), intent(in)  :: n_minus
  !   integer, intent(in)                 :: ith_plate_separation

  !   real(dp), dimension(size(lambda_plus)) :: lambda_plus_bulk
  !   real(dp), dimension(size(lambda_plus)) :: lambda_neutral_bulk
  !   real(dp), dimension(size(lambda_plus)) :: lambda_minus_bulk

  !   real(dp), dimension(size(lambda_plus)) :: lambda_common_terms
  !   real(dp), dimension(size(lambda_plus)) :: lambda_common_terms_bulk

  !   real(dp), dimension(size(lambda_plus)) :: lambda_plus_bulk
  !   real(dp), dimension(size(lambda_plus)) :: lambda_neutral_bulk
  !   real(dp), dimension(size(lambda_plus)) :: lambda_minus_bulk

  !   integer :: input_array_size

  !   ! First ensure that the sizes of all the input variables are the same.
  !   input_array_size = size(lambda_plus)
  !   if((size(lambda_neutral) /= input_array_size) .or. &
  !        (size(lambda_minus) /= input_array_size) .or. &
  !        (size(n_plus) /= input_array_size) .or. &
  !        (size(n_neutral) /= input_array_size) .or. &
  !        (size(n_minus) /= input_array_size)) then

  !      print *, "lambda.f90: CalculateLambdas:"
  !      print *, "Input arrays must all be the same size"
  !      print *, "Almost certainly coding error as opposed to input error...aborting..."
  !      call abort()
  !   end if

  !   lambda_plus(:) = 0.0_dp
  !   lambda_neutral(:) = 0.0_dp
  !   lambda_minus(:) = 0.0_dp

  !   lambda_plus_bulk(:) = 0.0_dp
  !   lambda_neutral_bulk(:) = 0.0_dp
  !   lambda_minus_bulk(:) = 0.0_dp

  !   ! Now calculate the terms in common to all the lambdas
  !   lambda_common_terms = CalculateLambdaCommonTerms(.false., n_plus, n_neutral, n_minus, ith_plate_separation)

  !   lambda_common_terms_bulk = CalculateLambdaCommonTerms(.true.)

  !   !print *, "lambda_common_terms = ", lambda_common_terms
  !   !print *, "lambda_common_terms_bulk = ", lambda_common_terms_bulk
  !   !call abort()

  !   ! Now calculate our lambdas(r)
  !   lambda_plus = beta * (lambda_common_terms + CalculateLambdaPlusSpecificTerms(n_plus, n_minus, .false.))
  !   lambda_neutral = beta * (lambda_common_terms + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .false.))
  !   lambda_minus = beta * (lambda_common_terms + CalculateLambdaMinusSpecificTerms(n_plus, n_minus, .false.))

  !   lambda_plus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaPlusSpecificTerms(n_plus, n_minus, .true.))
  !   lambda_neutral_bulk = beta * (lambda_common_terms_bulk + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .true.))
  !   lambda_minus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaMinusSpecificTerms(n_plus, n_minus, .true.))

  !   !print *, "lambda_neutral_bulk = ", lambda_neutral_bulk
  !   !print *, "lambda_neutral = ", lambda_neutral

  !   !Now calculate the bulk value, lambda^{bulk} and in order to return e^{lambda^{bulk} - lambda(r)}
  !   lambda_plus = lambda_plus_bulk - lambda_plus
  !   lambda_neutral = lambda_neutral_bulk - lambda_neutral
  !   lambda_minus = lambda_minus_bulk - lambda_minus

  !   !print *, "lambda_neutral diff = ", lambda_neutral

  ! end subroutine CalculateLambdasDifference_old

  function CalculateLambdaCommonTerms(calculate_bulk, size_of_terms, n_s, ith_plate_separation, iteration)
    logical, intent(in) :: calculate_bulk
    integer, intent(in) :: size_of_terms
    real(dp), dimension(:), intent(in), optional :: n_s
    integer, intent(in), optional :: ith_plate_separation
    integer, optional :: iteration

    real(dp), dimension(size_of_terms) :: CalculateLambdaCommonTerms

    real(dp), dimension(size_of_terms) :: hs_term
    real(dp), dimension(size_of_terms) :: van_der_waals_term
    real(dp), dimension(size_of_terms) :: surface_fluid_dispersion_term
    real(dp), dimension(size_of_terms) :: n_s_bulk

    real(dp) :: allowed_distance_between_plates

    hs_term = 0.0_dp
    van_der_waals_term = 0.0_dp
    surface_fluid_dispersion_term = 0.0_dp

    if(calculate_bulk) then
       if(present(n_s) .or. present(ith_plate_separation)) then
          print *, "lambdas.f90: CalculateLambdaCommonTerms:"
          print *, "When calculate_bulk is present, don't pass in densities as they're not needed.  Coding bug."
          call abort()
       end if
    else
       if((.not. present(n_s)) .or. (.not. present(ith_plate_separation))) then
          print *, "lambdas.f90: CalculateLambdaCommonTerms:"
          print *, "When calculate_bulk is NOT present, densities are required to be passed in.  Coding bug."
          call abort()
       end if
    end if


    if(calculate_bulk) then

       n_s_bulk(:) = bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads
       allowed_distance_between_plates = (((size(n_s) - 1)*hs_diameter/n_discretised_points_z) + ((hs_diameter)*1.0_dp))/2.0_dp

       hs_term = 1.0_dp * calculate_hardsphere_functional_deriv(n_s_bulk, .true.)
       !van_der_waals_term = calculate_vanderWaals_functional_deriv(n_s)

       van_der_waals_term = (-8.0_dp * epsilon_LJ_particle_particle * (hs_diameter**6) * pi * n_s_bulk(:)) * ( &
            (2.0_dp/(3.0_dp*(hs_diameter**3))))

       !van_der_waals_term = (-2.0_dp * epsilon_LJ_particle_particle * (hs_diameter**6) * pi * n_s_bulk(:)) * ( &
       !     (-1.0_dp / (3.0_dp * ((allowed_distance_between_plates)**3))) + (2.0_dp/(hs_diameter**3))) !+ &


       !call setNonCalculatedRegionToZero(n_s_bulk)

    else

       hs_term = 1.0_dp * calculate_hardsphere_functional_deriv(n_s, .false.)

       surface_fluid_dispersion_term = calculate_surface_dispersion_functional_deriv(&
            ith_plate_separation, size(surface_fluid_dispersion_term))

       van_der_waals_term = calculate_vanderWaals_functional_deriv(n_s)
    end if

    CalculateLambdaCommonTerms = hs_term + surface_fluid_dispersion_term + van_der_waals_term

  end function CalculateLambdaCommonTerms


  function CalculateLambdaPlusSpecificTerms(calculate_bulk, size_of_terms, n_plus, n_minus)
    logical, intent(in) :: calculate_bulk
    integer, intent(in) :: size_of_terms
    real(dp), dimension(:), intent(in), optional  :: n_plus
    real(dp), dimension(:), intent(in), optional  :: n_minus


    real(dp), dimension(size_of_terms) :: CalculateLambdaPlusSpecificTerms

    real(dp), dimension(size_of_terms) :: surface_electrostatic_term
    real(dp), dimension(size_of_terms) :: like_electrostatic_term
    real(dp), dimension(size_of_terms) :: unlike_electrostatic_term

    real(dp), dimension(size_of_terms)  :: n_plus_in
    real(dp), dimension(size_of_terms)  :: n_minus_in

    real(dp), dimension(size_of_terms) :: n_s

    surface_electrostatic_term = 0.0_dp
    like_electrostatic_term = 0.0_dp
    unlike_electrostatic_term = 0.0_dp

    if(calculate_bulk) then
       if( present(n_plus) .or. present(n_minus) ) then
          print *, "lambdas.f90: CalculateLambdaCommonTerms:"
          print *, "When calculate_bulk is present, don't pass in densities as they're not needed.  Coding bug."
          call abort()
       end if
    else
       if( (.not. present(n_plus)) .or. (.not. present(n_minus)) ) then
          print *, "lambdas.f90: CalculateLambdaCommonTerms:"
          print *, "When calculate_bulk is NOT present, densities are required to be passed in.  Coding bug."
          call abort()
       end if
    end if

    if(calculate_bulk) then

       n_plus_in(:) = bulk_density_positive_beads
       n_minus_in(:) = bulk_density_negative_beads

       !print *, "n_plus_in = ", n_plus_in
       !print *, "n_minus_in = ", n_minus_in

       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_plus_in, positive_bead_charge, .true.)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_minus_in, positive_bead_charge, negative_bead_charge, .true.)

       !n_s(:) = bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads
       !hs_term(:) = calculate_hard_sphere_functional_deriv(n_s, n_plus, n_cation_end, n_cation_nonend, n_anion_end, n_anion_nonend, .true.)
       
       !print *, "like_electrostatic_term = ", like_electrostatic_term
       !print *, "unlike_electrostatic_term = ", unlike_electrostatic_term
       !call abort()
    else
       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_plus, positive_bead_charge, .false.)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_minus, positive_bead_charge, negative_bead_charge, .false.)

       surface_electrostatic_term = calculate_surface_electrostatic_functional_deriv(size_of_terms, positive_bead_charge)

       !n_s(:) = n_plus + n_neutral + n_minus
       !print *, "calling n_plus"
       !hs_term(:) = calculate_hard_sphere_functional_deriv(n_s, n_plus, n_cation_end, n_cation_nonend, n_anion_end, n_anion_nonend, .false.)


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

    !print *, "calculate_bulk", calculate_bulk
    !print *, "surface_electrostatic_term", surface_electrostatic_term
    !print *, "like_electrostatic_term", like_electrostatic_term
    !print *, "unlike_electrostatic_term", unlike_electrostatic_term
    !print *, "hs_term", hs_term
    !call abort()

    CalculateLambdaPlusSpecificTerms = surface_electrostatic_term + like_electrostatic_term + unlike_electrostatic_term! + hs_term
    !if(calculate_bulk .eqv. .false.) then
    !   print *, "CalculateLambdaPlusSpecificTerms = ", CalculateLambdaPlusSpecificTerms
    !   call abort()
    !end if

  end function CalculateLambdaPlusSpecificTerms


  function CalculateLambaNeutralSpecificTerms(calculate_bulk, size_of_terms, n_plus, n_neutral, n_minus)
    logical, intent(in) :: calculate_bulk
    integer, intent(in) :: size_of_terms
    real(dp), dimension(:), intent(in), optional  :: n_plus
    real(dp), dimension(:), intent(in), optional  :: n_neutral
    real(dp), dimension(:), intent(in), optional  :: n_minus

    real(dp), dimension(size_of_terms) :: CalculateLambaNeutralSpecificTerms

    if(calculate_bulk) then
       if( present(n_plus) .or. present(n_neutral) .or. present(n_minus)) then
          print *, "lambdas.f90: CalculateLambdaNeutralSpecificTerms:"
          print *, "When calculate_bulk is present, don't pass in densities as they're not needed.  Coding bug."
          call abort()
       end if
    else
       if( (.not. present(n_plus)) .or. (.not. present(n_neutral)) .or. (.not. present(n_minus)) ) then
          print *, "lambdas.f90: CalculateLambdaNeutralSpecificTerms:"
          print *, "When calculate_bulk is NOT present, densities are required to be passed in.  Coding bug."
          call abort()
       end if
    end if


    CalculateLambaNeutralSpecificTerms = 0.0_dp

  end function CalculateLambaNeutralSpecificTerms


  function CalculateLambdaMinusSpecificTerms(calculate_bulk, size_of_terms, n_plus, n_minus)
    logical, intent(in) :: calculate_bulk
    integer, intent(in) :: size_of_terms
    real(dp), dimension(:), intent(in), optional  :: n_plus
    real(dp), dimension(:), intent(in), optional  :: n_minus

    real(dp), dimension(size_of_terms) :: CalculateLambdaMinusSpecificTerms

    real(dp), dimension(size_of_terms) :: surface_electrostatic_term
    real(dp), dimension(size_of_terms) :: like_electrostatic_term
    real(dp), dimension(size_of_terms) :: unlike_electrostatic_term
    
    real(dp), dimension(size_of_terms)  :: n_plus_in
    real(dp), dimension(size_of_terms)  :: n_minus_in

    real(dp), dimension(size_of_terms) :: n_s

    surface_electrostatic_term = 0.0_dp
    like_electrostatic_term = 0.0_dp
    unlike_electrostatic_term = 0.0_dp
    

    if(calculate_bulk) then
       if( present(n_plus) .or. present(n_minus) ) then
          print *, "lambdas.f90: CalculateLambdaMinusSpecificTerms:"
          print *, "When calculate_bulk is present, don't pass in densities as they're not needed.  Coding bug."
          call abort()
       end if
    else
       if( (.not. present(n_plus)) .or. (.not. present(n_minus)) ) then
          print *, "lambdas.f90: CalculateLambdaMinusSpecificTerms:"
          print *, "When calculate_bulk is NOT present, densities are required to be passed in.  Coding bug."
          call abort()
       end if
    end if

    if(calculate_bulk) then

       n_plus_in(:) = bulk_density_positive_beads
       n_minus_in(:) = bulk_density_negative_beads

       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_minus_in, negative_bead_charge, .true.)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_plus_in, positive_bead_charge, negative_bead_charge, .true.)

       !n_s(:) = bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads
       !hs_term(:) = calculate_hard_sphere_functional_deriv(n_s, n_minus, n_cation_end, n_cation_nonend, n_anion_end, n_anion_nonend, .true.)

    else

       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_minus, negative_bead_charge, .false.)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_plus, positive_bead_charge, negative_bead_charge, .false.)

       surface_electrostatic_term = calculate_surface_electrostatic_functional_deriv(size_of_terms, negative_bead_charge)

       !n_s(:) = n_plus + n_neutral + n_minus
       !print *, "calling n_minus"
       !hs_term(:) = calculate_hard_sphere_functional_deriv(n_s, n_minus, n_cation_end, n_cation_nonend, n_anion_end, n_anion_nonend, .false.)

       !print *, "surface_electrostatic_term = ",  surface_electrostatic_term(26), surface_electrostatic_term(size(surface_electrostatic_term) - 25), &
       !     surface_electrostatic_term(size(surface_electrostatic_term) - 25) - surface_electrostatic_term(26)
       !print *, ""
       !print *, "like_electrostatic_term minus = ",  like_electrostatic_term(26), like_electrostatic_term(size(like_electrostatic_term) - 25), &
       !     like_electrostatic_term(size(like_electrostatic_term) - 25) - like_electrostatic_term(26)
       !print *, ""
       !print *, "unlike_electrostatic_term minus = ", unlike_electrostatic_term(26), unlike_electrostatic_term(size(unlike_electrostatic_term) - 25), &
       !     unlike_electrostatic_term(size(unlike_electrostatic_term) - 25) - unlike_electrostatic_term(26)
       !print *, ""

       !call abort()


    end if



    ! print *, "minus terms"
    ! print *, "calculate_bulk = ", calculate_bulk
    ! print *, "surface_electrostatic_term = ", surface_electrostatic_term
    ! print *, ""
    ! print *, "like_electrostatic_term = ", like_electrostatic_term
    ! print *, ""
    ! print *, "unlike_electrostatic_term = ", unlike_electrostatic_term
    ! print *, ""
    ! print *, "hs_term = ", hs_term
    ! print *, ""
    ! print *, "**********"
    ! print *, ""
    !call abort()

    CalculateLambdaMinusSpecificTerms = surface_electrostatic_term + like_electrostatic_term + unlike_electrostatic_term ! + hs_term

  end function CalculateLambdaMinusSpecificTerms

end module lambdas
