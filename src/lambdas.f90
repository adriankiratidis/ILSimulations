!This is a module that calculates the lambdas of our model.
!Recall that lambda is defined to be the paritial derivative of
!all non-ideal contributions to the grand potential functional
!w.r.t the bead densities.
module lambdas
  use kinds
  use parameters
  implicit none

contains

  subroutine CalculateLambdas(lambda_plus, lambda_neutral, lambda_minus, n_plus, n_neutral, n_minus)
    real(dp), dimension(:), intent(out) :: lambda_plus
    real(dp), dimension(:), intent(out) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: lambda_minus
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    
  end subroutine CalculateLambdas


  subroutine CalculateCommonTerms(n_plus, n_neutral, n_minus)
  end subroutine CalculateCommonTerms

  subroutine CalculatePlusSpecificTerms(n_plus, n_neutral, n_minus)
  end subroutine CalculatePlusSpecificTerms

  subroutine CalculateMinusSpecificTerms(n_plus, n_neutral, n_minus)
  end subroutine CalculateMinusSpecificTerms

  subroutine CalculateNeutralSpecificTerms(n_plus, n_neutral, n_minus)
  end subroutine CalculateNeutralSpecificTerms

end module lambdas
