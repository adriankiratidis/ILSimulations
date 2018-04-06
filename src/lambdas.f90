!This is a module that calculates the lambdas of our model.
!Recall that lambda is defined to be the paritial derivative of
!all non-ideal contributions to the grand potential functional
!w.r.t the bead densities.
module lambdas
  use kinds
  use parameters
  use universalconstants
  use integratezcylindrical
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

    real(dp), dimension(size(n_plus)) :: n_s
    real(dp), dimension(size(n_plus)) :: n_sbar

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

    van_der_waals_term = -4.0_dp * epsilon_LJ * (hs_diameter**6.0_dp) * &
         integrate_z_cylindrical(n_s, van_der_waals_density_indept_integrand, "all_z")

    n_sbar = (3.0_dp * ( integrate_z_cylindrical(n_s, n_sbar_integrand, "z_lteq_hs_diameter") ))/(2.0_dp * (hs_diameter**3.0_dp))
    
    hs_term = (3.0_dp * (log( (1.0_dp - (hs_diameter**3)*n_sbar)/(n_sbar) ) - &
         (1.0_dp)/(1.0_dp - (hs_diameter**3.0_dp) * n_sbar)) ) / &
         (4.0_dp * pi * (hs_diameter**3.0_dp))

    
    CalculateLambdaCommonTerms = hs_term + van_der_waals_term + surface_fluid_dispersion_term

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

  
  function get_bulk_density(lambda)
    real(dp), dimension(:), intent(in)  :: lambda

    real(dp), dimension(size(lambda)) :: get_bulk_density

    !Setting the bulk density to be int(lambda)/plate_separation
    !We could of course calculate the total plate separation by hs_diameter * plate_separations(ith_separation)
    !but we choose the current version so we can calculate it without passing in an extra parameter.
    get_bulk_density = integrate_z_cylindrical(lambda, unity_function, "all_z") / &
         ( (real(size(lambda) - 1, dp) * hs_diameter)/n_discretised_points_z )

  end function get_bulk_density

  pure function van_der_waals_density_indept_integrand(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp)             :: van_der_waals_density_indept_integrand

    if(abs(z-xi) >= n_discretised_points_z) then
       van_der_waals_density_indept_integrand = 1.0_dp / (4.0_dp * (xi**4.0_dp))
    else
       van_der_waals_density_indept_integrand = 1.0_dp / &
            (4.0_dp * ( (hs_diameter * cos(asin(xi/hs_diameter)))**2.0_dp + xi**2.0_dp ))
    end if

  end function van_der_waals_density_indept_integrand

  pure function n_sbar_integrand(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp)            :: n_sbar_integrand
    
    n_sbar_integrand = (hs_diameter**2.0_dp + (z - xi)**2.0_dp) / 2.0_dp

  end function n_sbar_integrand
  
  pure function unity_function(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp) :: unity_function

    unity_function = 1.0_dp
  end function unity_function
  
end module lambdas
