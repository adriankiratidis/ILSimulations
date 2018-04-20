!Contains routines to calculate the functional derivative of each
!external (non-ideal) term in the grand potential.
module functionalderivatives
  use kinds
  use universalconstants
  use parameters
  use helpers
  use integratezcylindrical
  implicit none
  private

  public :: calculate_vanderWaals_functional_deriv
  public :: calculate_surface_dispersion_functional_deriv
  public :: calculate_hardsphere_functional_deriv

  !public :: calculate_surface_electrostatic_functional_deriv
  !public :: calculate_electrostatic_like_term_functional_deriv
  !public :: calculate_electrostatic_unlike_term_functional_deriv
  public :: calculate_n_sbar

contains

  function calculate_hardsphere_functional_deriv(n_s)
    real(dp), dimension(:) :: n_s
    real(dp), dimension(size(n_s)) :: calculate_hardsphere_functional_deriv

    real(dp), dimension(size(n_s)) :: n_sbar

    integer :: start_z_index
    integer :: end_z_index

    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    call get_allowed_z_values(start_z_index, end_z_index, size(n_s))

    n_sbar = calculate_n_sbar(n_s)

    calculate_hardsphere_functional_deriv(start_z_index:end_z_index) = (-1.0_dp / beta) * &
         (3.0_dp * (log( (1.0_dp - (hs_diameter**3)*n_sbar(start_z_index:end_z_index))/(n_sbar(start_z_index:end_z_index)) ) - &
         (1.0_dp)/(1.0_dp - (hs_diameter**3.0_dp) * n_sbar(start_z_index:end_z_index))) ) / &
         (4.0_dp * pi * (hs_diameter**3.0_dp))

    call setNonCalculatedRegionToZero(calculate_hardsphere_functional_deriv)

  end function calculate_hardsphere_functional_deriv

  function calculate_surface_dispersion_functional_deriv(ith_plate_separation, array_size)
    integer, intent(in) :: ith_plate_separation
    integer, intent(in) :: array_size
    real(dp), dimension(array_size) :: calculate_surface_dispersion_functional_deriv

    integer :: start_z_index
    integer :: end_z_index
    integer :: iz

    real(dp) :: hs_d_divide_z
    real(dp) :: hs_d_divide_h_minus_z

    calculate_surface_dispersion_functional_deriv = 0.0_dp

    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    call get_allowed_z_values(start_z_index, end_z_index, array_size)

    !Note that we exclude the points that are calculate on the wall (at r_{z} = 0 or r_{z} = h)
    !as this leads to a singularity.
    do iz = start_z_index, end_z_index

       hs_d_divide_z = (real(n_discretised_points_z, dp) / real((iz - 1), dp))
       hs_d_divide_h_minus_z = (1.0_dp/(real(plate_separations(ith_plate_separation),dp) - &
            ( real((iz - 1),dp) / real((n_discretised_points_z),dp) )))

       calculate_surface_dispersion_functional_deriv(iz) = 2.0_dp * pi * epsilon_LJ * (&
            ( (2.0_dp/45.0_dp)* (hs_d_divide_z**9.0_dp) ) - &
            ( (1.0_dp/3.0_dp) * (hs_d_divide_z**3.0_dp) ) + &
            ( (2.0_dp/45.0_dp)* (hs_d_divide_h_minus_z**9.0_dp) ) - &
            ( (1.0_dp/3.0_dp) * (hs_d_divide_h_minus_z**3.0_dp) ))
    end do

    call setNonCalculatedRegionToZero(calculate_surface_dispersion_functional_deriv)

  end function calculate_surface_dispersion_functional_deriv
           
  function calculate_vanderWaals_functional_deriv(n_s)
    real(dp), dimension(:), intent(in) :: n_s
    real(dp), dimension(size(n_s)) :: calculate_vanderWaals_functional_deriv

    calculate_vanderWaals_functional_deriv = -4.0_dp * epsilon_LJ * (hs_diameter**6.0_dp) * &
         2.0_dp * pi * integrate_z_cylindrical(n_s, van_der_waals_density_indept_integrand, "all_z")

    call setNonCalculatedRegionToZero(calculate_vanderWaals_functional_deriv)

  end function calculate_vanderWaals_functional_deriv

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

    return
  end function van_der_waals_density_indept_integrand

  function calculate_n_sbar(n_s)
    real(dp), dimension(:), intent(in) :: n_s
    real(dp), dimension(size(n_s)) :: calculate_n_sbar

    calculate_n_sbar = (3.0_dp * ( integrate_z_cylindrical(n_s, n_sbar_integrand, "z_lteq_hs_diameter") ))&
         /(4.0_dp * pi * (hs_diameter**3.0_dp))

  end function calculate_n_sbar

  pure function n_sbar_integrand(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp)            :: n_sbar_integrand

    n_sbar_integrand = 2.0_dp * pi * (hs_diameter**2.0_dp - ((z - xi)*hs_diameter/real(n_discretised_points_z,dp))**2.0_dp) / 2.0_dp

  end function n_sbar_integrand

end module functionalderivatives