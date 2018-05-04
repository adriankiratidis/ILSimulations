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

  public :: CalculateLambdas

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

    ! Now calculate the terms in common to all the lambdas
    lambda_common_terms = CalculateLambdaCommonTerms(n_plus, n_neutral, n_minus, ith_plate_separation)

    !print *, "n_plus = ", n_plus
    !print *, "n_neutral = ", n_neutral
    !print *, "n_minus = ", n_minus
    !print *, "ith_plate_separation = ", ith_plate_separation  
    !print *, "lambda_common_terms = ", lambda_common_terms

    ! Now calculate our lambdas(r)
    lambda_plus = beta * (lambda_common_terms + CalculateLambdaPlusSpecificTerms(n_plus, n_neutral, n_minus))
    lambda_neutral = beta * (lambda_common_terms + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus))
    lambda_minus = beta * (lambda_common_terms + CalculateLambdaMinusSpecificTerms(n_plus, n_neutral, n_minus))

    !Now calculate the bulk value, lambda^{bulk} and in order to return e^{lambda^{bulk} - lambda(r)}
    lambda_plus = get_bulk_density(lambda_plus) - lambda_plus
    lambda_neutral = get_bulk_density(lambda_neutral) - lambda_neutral
    lambda_minus = get_bulk_density(lambda_minus) - lambda_minus

  end subroutine CalculateLambdas

  function CalculateLambdaCommonTerms(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
    real(dp), dimension(:), intent(in)  :: n_minus
    integer, intent(in)                 :: ith_plate_separation

    real(dp), dimension(size(n_plus)) :: CalculateLambdaCommonTerms

    real(dp), dimension(size(n_plus)) :: hs_term
    real(dp), dimension(size(n_plus)) :: van_der_waals_term
    real(dp), dimension(size(n_plus)) :: surface_fluid_dispersion_term

    real(dp), dimension(size(n_plus)) :: n_s

    n_s = n_plus + n_neutral + n_minus

    !van_der_waals_term = calculate_vanderWaals_functional_deriv(n_s)
                         
    !surface_fluid_dispersion_term = calculate_surface_dispersion_functional_deriv(&
    !     ith_plate_separation, size(surface_fluid_dispersion_term))

    hs_term = calculate_hardsphere_functional_deriv(n_s)

    !print *, "hs_term= ", hs_term
    !print *, "van_der_waals_term = ", van_der_waals_term
    !print *, "surface_fluid_dispersion_term = ",&
    !     surface_fluid_dispersion_term

    CalculateLambdaCommonTerms = hs_term !+ van_der_waals_term + surface_fluid_dispersion_term

  end function CalculateLambdaCommonTerms


  function CalculateLambdaPlusSpecificTerms(n_plus, n_neutral, n_minus)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
    real(dp), dimension(:), intent(in)  :: n_minus

    real(dp), dimension(size(n_plus)) :: CalculateLambdaPlusSpecificTerms

    CalculateLambdaPlusSpecificTerms = 0.0_dp
    !CalculateLambdaPlusSpecificTerms = bead_charge * /(4.0_dp * pi * epsilonr * epsilon0)

  end function CalculateLambdaPlusSpecificTerms


  function CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
    real(dp), dimension(:), intent(in)  :: n_minus

    real(dp), dimension(size(n_plus)) :: CalculateLambaNeutralSpecificTerms

    CalculateLambaNeutralSpecificTerms = 0.0_dp

  end function CalculateLambaNeutralSpecificTerms


  function CalculateLambdaMinusSpecificTerms(n_plus, n_neutral, n_minus)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
    real(dp), dimension(:), intent(in)  :: n_minus

    real(dp), dimension(size(n_plus)) :: CalculateLambdaMinusSpecificTerms

    CalculateLambdaMinusSpecificTerms = 0.0_dp

  end function CalculateLambdaMinusSpecificTerms

end module lambdas
