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
    lambda_plus = beta * (lambda_common_terms + CalculateLambdaPlusSpecificTerms(n_plus, n_neutral, n_minus, .false.))
    lambda_neutral = beta * (lambda_common_terms + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .false.))
    lambda_minus = beta * (lambda_common_terms + CalculateLambdaMinusSpecificTerms(n_plus, n_neutral, n_minus, .false.))
    
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

    lambda_plus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaPlusSpecificTerms(n_plus, n_neutral, n_minus, .true.))
    lambda_neutral_bulk = beta * (lambda_common_terms_bulk + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .true.))
    lambda_minus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaMinusSpecificTerms(n_plus, n_neutral, n_minus, .true.))

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

    !print *, "lambda_common_terms = ", lambda_common_terms
    !print *, "lambda_common_terms_bulk = ", lambda_common_terms_bulk
    !call abort()

    ! Now calculate our lambdas(r)
    lambda_plus = beta * (lambda_common_terms + CalculateLambdaPlusSpecificTerms(n_plus, n_neutral, n_minus, .false.))
    lambda_neutral = beta * (lambda_common_terms + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .false.))
    lambda_minus = beta * (lambda_common_terms + CalculateLambdaMinusSpecificTerms(n_plus, n_neutral, n_minus, .false.))

    lambda_plus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaPlusSpecificTerms(n_plus, n_neutral, n_minus, .true.))
    lambda_neutral_bulk = beta * (lambda_common_terms_bulk + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, .true.))
    lambda_minus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaMinusSpecificTerms(n_plus, n_neutral, n_minus, .true.))

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

    hs_term = 0.0_dp
    van_der_waals_term = 0.0_dp
    surface_fluid_dispersion_term = 0.0_dp

    if(calculate_bulk) then
       n_s(:) = bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads

       hs_term = calculate_hardsphere_functional_deriv(n_s, .true.)

       !print *, "hs_term bulk = ", beta * hs_term
    else
       n_s = n_plus + n_neutral + n_minus
       hs_term = calculate_hardsphere_functional_deriv(n_s, .false.)

       surface_fluid_dispersion_term = calculate_surface_dispersion_functional_deriv(&
            ith_plate_separation, size(surface_fluid_dispersion_term))


       !print *, "hs_term non bulk = ", beta * hs_term
       !print *, "surface_fluid_dispersion_term non bulk = ", beta * surface_fluid_dispersion_term

       !van_der_waals_term = calculate_vanderWaals_functional_deriv(n_s)
    end if



    CalculateLambdaCommonTerms = hs_term + surface_fluid_dispersion_term + van_der_waals_term

    !print *, "common terms = ", CalculateLambdaCommonTerms

  end function CalculateLambdaCommonTerms


  function CalculateLambdaPlusSpecificTerms(n_plus, n_neutral, n_minus, calculate_bulk)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
    real(dp), dimension(:), intent(in)  :: n_minus
    logical, intent(in) :: calculate_bulk

    real(dp), dimension(size(n_plus)) :: CalculateLambdaPlusSpecificTerms

    CalculateLambdaPlusSpecificTerms = 0.0_dp
    !CalculateLambdaPlusSpecificTerms = bead_charge * /(4.0_dp * pi * epsilonr * epsilon0)

  end function CalculateLambdaPlusSpecificTerms


  function CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus, calculate_bulk)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
    real(dp), dimension(:), intent(in)  :: n_minus
    logical, intent(in) :: calculate_bulk

    real(dp), dimension(size(n_plus)) :: CalculateLambaNeutralSpecificTerms

    CalculateLambaNeutralSpecificTerms = 0.0_dp

  end function CalculateLambaNeutralSpecificTerms


  function CalculateLambdaMinusSpecificTerms(n_plus, n_neutral, n_minus, calculate_bulk)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
    real(dp), dimension(:), intent(in)  :: n_minus
    logical, intent(in) :: calculate_bulk

    real(dp), dimension(size(n_plus)) :: CalculateLambdaMinusSpecificTerms

    CalculateLambdaMinusSpecificTerms = 0.0_dp

  end function CalculateLambdaMinusSpecificTerms

end module lambdas
