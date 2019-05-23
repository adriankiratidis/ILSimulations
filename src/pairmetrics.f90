module pairmetrics
  use kinds
  use parameters
  use integratezcylindrical
  use helpers
  implicit none
  private

  !public :: pair_density_evaluation_integrand
  public :: CalculateCentreLengthDifference
  

contains

  subroutine CalculateCentreLengthDifference(Centre_to_Centre_distance_metric, n_cation_centre, n_anion_centre)
    real(dp), intent(out) :: Centre_to_Centre_distance_metric
    real(dp), dimension(:), intent(in) :: n_cation_centre
    real(dp), dimension(:), intent(in) :: n_anion_centre

    Centre_to_Centre_distance_metric = integrate_z_cylindrical(n_cation_centre * &
         integrate_z_cylindrical(n_anion_centre, pair_density_evaluation_integrand, "all_z"), unity_function)
    
  end subroutine CalculateCentreLengthDifference

  function pair_density_evaluation_integrand(z, xi_in)
    integer, intent(in) :: z
    integer, intent(in) :: xi_in
    real(dp)             :: pair_density_evaluation_integrand

    integer :: xi_int
    real(dp) :: xi_real

    xi_int = xi_in - z ! centre xi on z
    xi_real = real(xi_int,dp) * hs_diameter / real(n_discretised_points_z,dp)

    if(abs(xi_int) >= n_discretised_points_z) then
       pair_density_evaluation_integrand = -2.0_dp / (3.0_dp * (real(xi_real,dp)**3.0_dp))
    else
       ! pair_density_evaluation_integrand = 1.0_dp / &
       !      (4.0_dp * ( (real(hs_diameter,dp) * cos(asin(real(xi_real,dp)/real(hs_diameter,dp))))**2.0_dp &
       !      + real(xi_real,dp)**2.0_dp )**2.0_dp)
       pair_density_evaluation_integrand = -2.0_dp / (3.0_dp * (hs_diameter**3.0_dp))
    end if

  end function pair_density_evaluation_integrand

end module pairmetrics
