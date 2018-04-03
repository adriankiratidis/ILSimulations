!This is a module that calculates the lambdas of our model.
!Recall that lambda is defined to be the paritial derivative of
!all non-ideal contributions to the grand potential functional
!w.r.t the bead densities.
module lambdas
  use kinds
  use parameters
  use universalconstants
  implicit none
  private

  public :: CalculateLambdas
  
contains

  subroutine CalculateLambdas(lambda_plus, lambda_neutral, lambda_minus, n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus
    real(dp), dimension(:), intent(out) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: lambda_minus
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
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

    ! Now calculate our lambdas(r)
    lambda_plus = lambda_common_terms + CalculateLambdaPlusSpecificTerms(n_plus, n_neutral, n_minus)
    lambda_neutral = lambda_common_terms + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus)
    lambda_minus = lambda_common_terms + CalculateLambdaMinusSpecificTerms(n_plus, n_neutral, n_minus)


    !Now calculate the bulk value, lambda^{bulk} and in order to return e^{lambda^{bulk} - lambda(r)}
    lambda_plus = exp(get_bulk_density(lambda_plus) - lambda_plus)
    lambda_neutral = exp(get_bulk_density(lambda_neutral) - lambda_neutral)
    lambda_minus = exp(get_bulk_density(lambda_minus) - lambda_minus)
    
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
    
    real(dp), dimension(size(n_plus) :: n_s
    
    integer :: iz

    surface_fluid_dispersion_term = 0.0_dp
    !Note that we exclude the points that are calculate on the wall (at r_{z} = 0 or r_{z} = h)
    !as this leads to a singularity.
    do iz = 2, size(surface_fluid_dispersion_term) - 1

       surface_fluid_dispersion_term(iz) = 2.0_dp * pi * epsilon_LJ * (&
            ( ( (2.0_dp/45.0_dp)* ( (n_discretised_points_z/(iz - 1) )**9.0_dp) ) - &
            ( (1.0_dp/3.0_dp)*( (n_discretised_points_z/(iz - 1) )**3.0_dp)) )  + &
            ( ( (2.0_dp/45.0_dp)* ( (1.0_dp/(plate_separations(ith_plate_separation) - &
            ((iz - 1)/(n_discretised_points_z))) ))**9.0_dp) ) - &
            ( (1.0_dp/3.0_dp)*( (1.0_dp/(plate_separations(ith_plate_separation) - &
            ((iz - 1)/(n_discretised_points_z))) ))**3.0_dp) )
    end do

    n_s = n_plus + n_neutral + n_minus
    
    van_der_waals_term = -4.0_dp * epsilon_LJ * (hs_diameter**6.0_dp) * integrate_z_cylindrical(n_s, gt_sigma = .true.)

    CalculateLambdaCommonTerms = hs_term + van_der_waals_term + surface_fluid_dispersion_term

  end function CalculateLambdaCommonTerms

  function CalculateLambdaPlusSpecificTerms(n_plus, n_neutral, n_minus)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
    real(dp), dimension(:), intent(in)  :: n_minus

    real(dp), dimension(size(n_plus)) :: CalculateLambdaPlusSpecificTerms

  end function CalculateLambdaPlusSpecificTerms

  
  function CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
    real(dp), dimension(:), intent(in)  :: n_minus

    real(dp), dimension(size(n_plus)) :: CalculateLambaNeutralSpecificTerms

  end function CalculateLambaNeutralSpecificTerms


  function CalculateLambdaMinusSpecificTerms(n_plus, n_neutral, n_minus)
    real(dp), dimension(:), intent(in)  :: n_plus
    real(dp), dimension(:), intent(in)  :: n_neutral
    real(dp), dimension(:), intent(in)  :: n_minus

    real(dp), dimension(size(n_plus)) :: CalculateLambdaMinusSpecificTerms

  end function CalculateLambdaMinusSpecificTerms

  function get_bulk_density(lambda)
    real(dp), dimension(:), intent(in)  :: lambda

    real(dp), dimension(size(lambda)) :: get_bulk_density

    get_bulk_density = integrate_z_cylindrical(lambda, "all_z") / real(size(lambda), dp)

  end function get_bulk_density

end module lambdas
