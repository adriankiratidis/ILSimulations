!This is a module that calculates the lambdas of our model.
!Recall that lambda is defined to be the paritial derivative of
!all non-ideal contributions to the grand potential functional
!w.r.t the bead densities.
module lambdas
  use kinds
  use parameters
  use universalconstants
  use integratezcylindrical
  use helpers
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

    !print *, "n_plus = ", n_plus
    !print *, "n_neutral = ", n_neutral
    !print *, "n_minus = ", n_minus
    !print *, "ith_plate_separation = ", ith_plate_separation  
    !print *, "lambda_common_terms = ", lambda_common_terms

    ! Now calculate our lambdas(r)
    lambda_plus = beta * (lambda_common_terms + CalculateLambdaPlusSpecificTerms(n_plus, n_neutral, n_minus))
    lambda_neutral = beta * (lambda_common_terms + CalculateLambaNeutralSpecificTerms(n_plus, n_neutral, n_minus))
    lambda_minus = beta * (lambda_common_terms + CalculateLambdaMinusSpecificTerms(n_plus, n_neutral, n_minus))

    ! print *, "error = ", get_bulk_density(lambda_plus) - lambda_plus
    ! call abort()
    print *, "lambda_plus = ", lambda_plus
    print *, "lambda_plus bulk density = ", get_bulk_density(lambda_plus)
    print *, "lambda_plus_diff = ", get_bulk_density(lambda_plus) - lambda_plus
    
    !Now calculate the bulk value, lambda^{bulk} and in order to return e^{lambda^{bulk} - lambda(r)}
    lambda_plus = exp(get_bulk_density(lambda_plus) - lambda_plus)
    lambda_neutral = exp(get_bulk_density(lambda_neutral) - lambda_neutral)
    lambda_minus = exp(get_bulk_density(lambda_minus) - lambda_minus)

    print *, "lambda_plus = ", lambda_plus
    ! !print *, "size(lambda_plus) = ", size(lambda_plus)
    print *, "lambda_neutral = ", lambda_neutral
    print *, "lambda_minus = ", lambda_minus
    !call abort()


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
    real(dp) :: hs_d_divide_z
    real(dp) :: hs_d_divide_h_minus_z

    integer :: start_z_index
    integer :: end_z_index

    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    call get_allowed_z_values(start_z_index, end_z_index, size(surface_fluid_dispersion_term))

    surface_fluid_dispersion_term = 0.0_dp
    !Note that we exclude the points that are calculate on the wall (at r_{z} = 0 or r_{z} = h)
    !as this leads to a singularity.
    do iz = start_z_index, end_z_index

       hs_d_divide_z = (real(n_discretised_points_z, dp) / real((iz - 1), dp))
       hs_d_divide_h_minus_z = (1.0_dp/(real(plate_separations(ith_plate_separation),dp) - &
            ( real((iz - 1),dp) / real((n_discretised_points_z),dp) )))

       surface_fluid_dispersion_term(iz) = 2.0_dp * pi * epsilon_LJ * (&
            ( (2.0_dp/45.0_dp)* (hs_d_divide_z**9.0_dp) ) - &
            ( (1.0_dp/3.0_dp) * (hs_d_divide_z**3.0_dp) ) + &
            ( (2.0_dp/45.0_dp)* (hs_d_divide_h_minus_z**9.0_dp) ) - &
            ( (1.0_dp/3.0_dp) * (hs_d_divide_h_minus_z**3.0_dp) ))
    end do

    n_s = n_plus + n_neutral + n_minus

    van_der_waals_term = -4.0_dp * epsilon_LJ * (hs_diameter**6.0_dp) * &
         integrate_z_cylindrical(n_s, van_der_waals_density_indept_integrand, "all_z")

    !print *, "n_s = ", n_s
    
    n_sbar = (3.0_dp * ( integrate_z_cylindrical(n_s, n_sbar_integrand, "z_lteq_hs_diameter") ))&
         /(4.0_dp * pi * (hs_diameter**3.0_dp))

    !print *, "n_sbar = ", n_sbar
    !call abort
    
    hs_term(start_z_index:end_z_index) = (-1.0_dp / beta) * &
         (3.0_dp * (log( (1.0_dp - (hs_diameter**3)*n_sbar(start_z_index:end_z_index))/(n_sbar(start_z_index:end_z_index)) ) - &
         (1.0_dp)/(1.0_dp - (hs_diameter**3.0_dp) * n_sbar(start_z_index:end_z_index))) ) / &
         (4.0_dp * pi * (hs_diameter**3.0_dp))

    print *, "hs_term= ",hs_term
    print *, "van_der_waals_term = ", van_der_waals_term
    print *, "surface_fluid_dispersion_term = ",&
         surface_fluid_dispersion_term
    !call abort()

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

  
  function get_bulk_density(lambda) result(reslt)
    real(dp), dimension(:), intent(in)  :: lambda

    real(dp), dimension(size(lambda)) :: reslt
    
    !Setting the bulk density to be int(lambda)/plate_separation
    !We could of course calculate the total plate separation by hs_diameter * plate_separations(ith_separation)
    !but we choose the current version so we can calculate it without passing in an extra parameter.
    reslt = integrate_z_cylindrical(lambda, unity_function, "all_z") / &
         ( (real(size(lambda) - 1, dp) * hs_diameter)/real(n_discretised_points_z,dp) )
    return
  end function get_bulk_density

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

    ! print *, "van_der_waals_density_indept_integrand = ", z, xi_in, xi_int, z - xi_int, van_der_waals_density_indept_integrand
    ! print *, ""
    ! if(isnan(van_der_waals_density_indept_integrand)) then
    !    print *, xi_real, hs_diameter
    !    !print *, asin(3.0_dp/2.4_dp)
    !    call abort
    ! end if
    return
  end function van_der_waals_density_indept_integrand

  pure function n_sbar_integrand(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp)            :: n_sbar_integrand
    
    n_sbar_integrand = 2.0_dp * pi * (hs_diameter**2.0_dp - ((z - xi)*hs_diameter/real(n_discretised_points_z,dp))**2.0_dp) / 2.0_dp

  end function n_sbar_integrand
  
  pure function unity_function(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp) :: unity_function

    unity_function = 1.0_dp
  end function unity_function
  
end module lambdas
