!Contains routines to calculate the functional derivative of each
!external (non-ideal) term in the grand potential.
module functionalderivatives
  use kinds
  use universalconstants
  use parameters
  use helpers
  use integratezcylindrical
  use io
  use excessenergyfunctionalparameters
  implicit none
  private

  public :: calculate_vanderWaals_functional_deriv
  public :: calculate_surface_dispersion_functional_deriv
  public :: calculate_hardsphere_functional_deriv
  public :: calculate_surface_electrostatic_functional_deriv

  public :: calculate_electrostatic_like_term_functional_deriv
  public :: calculate_electrostatic_unlike_term_functional_deriv

  public :: calculate_n_sbar

  real(dp) :: CURRENT_BULK_BEAD_DENSITY

contains

  function calculate_hardsphere_functional_deriv(n_s, calculate_bulk)
    real(dp), dimension(:), intent(in) :: n_s
    logical, intent(in)                :: calculate_bulk
    real(dp), dimension(size(n_s)) :: calculate_hardsphere_functional_deriv

    real(dp), dimension(size(n_s)) :: n_mbar, n_sbar
    real(dp), dimension(size(n_s)) :: integrand, integral

    integer :: start_z_index
    integer :: end_z_index

    integrand(:) = 0.0_dp
    integral(:) = 0.0_dp

    n_mbar(:) = 0.0_dp
    n_sbar(:) = 0.0_dp
    calculate_hardsphere_functional_deriv(:) = 0.0_dp

    n_mbar(:) = calculate_n_sbar(n_s(:))

    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    call get_allowed_z_values(start_z_index, end_z_index, size(n_s))

    if(calculate_bulk) then !n_s input parameter is the bulk value so


       calculate_hardsphere_functional_deriv(start_z_index:end_z_index) = (0.5_dp / beta) * (&
            GetAEx(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), a_term_index) + &
            (((4.0_dp * pi * (hs_diameter**3))/3.0_dp) * ((n_s(start_z_index:end_z_index)) * &
            GetAExDerivIntegrand(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), a_term_index))))




    else

       integrand(start_z_index:end_z_index) = n_s(start_z_index:end_z_index) * &
            GetAExDerivIntegrand(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), a_term_index)

       integral(:) = calculate_n_sbar(integrand(:))

       calculate_hardsphere_functional_deriv(start_z_index:end_z_index) = (0.5_dp / beta) * (&
            GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), a_term_index) + &
            (((4.0_dp * pi * (hs_diameter**3))/3.0_dp) * integral(start_z_index:end_z_index)))



    end if


    ! if(calculate_bulk) then !n_s input parameter is the bulk value so, 

    !    ! calculate_hardsphere_functional_deriv(start_z_index:end_z_index) = (1.0_dp / beta) * (&
    !    !      log(n_s(start_z_index:end_z_index) / (1.0_dp - ((hs_diameter**3)*n_s(start_z_index:end_z_index)))) + &
    !    !      ((n_s(start_z_index:end_z_index) * (hs_diameter**3))/(1.0_dp - ((hs_diameter**3)*n_s(start_z_index:end_z_index)))) + &
    !    !      1.0_dp )

    !    extra_integral_contribution2(:) = 0.0_dp
    !    extra_integral_contribution2(start_z_index:end_z_index) = (hs_diameter**3)/(1.0_dp - (hs_diameter**3)*n_s(start_z_index:end_z_index))
    !    extra_integral_contribution2 = calculate_n_sbar(extra_integral_contribution2)

    !    calculate_hardsphere_functional_deriv(start_z_index:end_z_index) = (1.0_dp / beta) * (&
    !         log(n_s(start_z_index:end_z_index) / (1.0_dp - ((hs_diameter**3)*n_s(start_z_index:end_z_index)))) + &
    !         extra_integral_contribution2(start_z_index:end_z_index) + &
    !         1.0_dp )


    ! else
    !    n_sbar = calculate_n_sbar(n_s)

    !    extra_integral_contribution(:) = 0.0_dp
    !    extra_integral_contribution(start_z_index:end_z_index) = n_s(start_z_index:end_z_index)/n_sbar(start_z_index:end_z_index)
    !    extra_integral_contribution = calculate_n_sbar(extra_integral_contribution)

    !    extra_integral_contribution2(:) = 0.0_dp
    !    extra_integral_contribution2(start_z_index:end_z_index) = (hs_diameter**3)/(1.0_dp - (hs_diameter**3)*n_sbar(start_z_index:end_z_index))
    !    extra_integral_contribution2 = calculate_n_sbar(extra_integral_contribution2)

    !    calculate_hardsphere_functional_deriv(start_z_index:end_z_index) = (1.0_dp / beta) * (&
    !         log(n_sbar(start_z_index:end_z_index) / (1.0_dp - ((hs_diameter**3)*n_sbar(start_z_index:end_z_index)))) + &
    !         extra_integral_contribution2(start_z_index:end_z_index) + &
    !         extra_integral_contribution(start_z_index:end_z_index) )



    !    ! end if
    !    !call abort()

    ! end if
    !calculate_hardsphere_functional_deriv = 0.0_dp
    !calculate_hardsphere_functional_deriv = 0.5_dp * calculate_hardsphere_functional_deriv
    !print *, "bulk density n_s = ", get_bulk_density(n_s)
    !print *, "n_sbar = ", n_sbar
    !print *, "bulk density n_sbar = ", get_bulk_density(n_sbar)
    !print *, "f_hs = ", calculate_hardsphere_functional_deriv

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

    calculate_surface_dispersion_functional_deriv(:) = 0.0_dp

    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    call get_allowed_z_values(start_z_index, end_z_index, array_size)

    !Note that we exclude the points that are calculate on the wall (at r_{z} = 0 or r_{z} = h)
    !as this leads to a singularity.
    do iz = start_z_index, end_z_index

       ! hs_d_divide_z = (real(n_discretised_points_z, dp) / real((iz - start_z_index), dp))
       ! hs_d_divide_h_minus_z = (1.0_dp/(real(plate_separations(ith_plate_separation),dp) - &
       !      ( real((iz - start_z_index),dp) / real((n_discretised_points_z),dp) )))

       hs_d_divide_z = (real(n_discretised_points_z, dp) / real((iz - 1), dp))
       hs_d_divide_h_minus_z = (1.0_dp/(real(plate_separations(ith_plate_separation),dp) - &
            ( real((iz - 1),dp) / real((n_discretised_points_z),dp) )))

       calculate_surface_dispersion_functional_deriv(iz) = 2.0_dp * pi * epsilon_LJ * (1.0_dp / beta) *(&
            ( (2.0_dp/45.0_dp)* (hs_d_divide_z**9.0_dp) ) - &
            ( (1.0_dp/3.0_dp) * (hs_d_divide_z**3.0_dp) ) + &
            ( (2.0_dp/45.0_dp)* (hs_d_divide_h_minus_z**9.0_dp) ) - &
            ( (1.0_dp/3.0_dp) * (hs_d_divide_h_minus_z**3.0_dp) ))
    end do

    call setNonCalculatedRegionToZero(calculate_surface_dispersion_functional_deriv)

    ! calculate_surface_dispersion_functional_deriv = 0.0_dp

  end function calculate_surface_dispersion_functional_deriv

  function calculate_vanderWaals_functional_deriv(n_s)
    real(dp), dimension(:), intent(in) :: n_s
    real(dp), dimension(size(n_s)) :: calculate_vanderWaals_functional_deriv

    calculate_vanderWaals_functional_deriv = -4.0_dp * epsilon_LJ * (hs_diameter**6.0_dp) * (1.0_dp/beta) * &
         2.0_dp * pi * integrate_z_cylindrical(n_s, van_der_waals_density_indept_integrand, "all_z")

    call setNonCalculatedRegionToZero(calculate_vanderWaals_functional_deriv)
    !calculate_vanderWaals_functional_deriv = 0.0_dp

  end function calculate_vanderWaals_functional_deriv

  function calculate_surface_electrostatic_functional_deriv(size_array, charge)
    integer, intent(in) :: size_array
    real(dp), intent(in) :: charge

    real(dp), dimension(size_array) :: calculate_surface_electrostatic_functional_deriv

    integer :: ij

    real(dp) :: d_to_left_wall
    real(dp) :: d_to_right_wall

    do ij = 1, size_array
       d_to_left_wall = (ij-1) * hs_diameter / n_discretised_points_z
       d_to_right_wall = (size_array-ij) * hs_diameter / n_discretised_points_z

       calculate_surface_electrostatic_functional_deriv(ij) = (1.0_dp / (4.0_dp * pi * epsilon0 * epsilonr)) * (&
            -2.0_dp * charge * pi * surface_charge_density_left_wall * d_to_left_wall &
            -2.0_dp * charge * pi * surface_charge_density_right_wall * d_to_right_wall)
    end do

  end function calculate_surface_electrostatic_functional_deriv

  function calculate_electrostatic_like_term_functional_deriv(n, charge)
    real(dp), dimension(:), intent(in) :: n
    real(dp), intent(in) :: charge
    real(dp), dimension(size(n)) :: calculate_electrostatic_like_term_functional_deriv

    if(charge > 0.0_dp) then
       CURRENT_BULK_BEAD_DENSITY = bulk_density_positive_beads
    else if(charge < 0.0_dp) then
       CURRENT_BULK_BEAD_DENSITY = bulk_density_negative_beads
    else if(charge < 0.0_dp) then
       !Note this is a float comparison so should never happen.
       !In any case if charge = 0.0_dp then we get an answer of zero anyway.
       !But we want to set it to something, to protect against zero * {some unset variable},
       !which I don't know what may happen.
       CURRENT_BULK_BEAD_DENSITY = bulk_density_neutral_beads
    end if

    calculate_electrostatic_like_term_functional_deriv(:) = (-1.0_dp / (4.0_dp * pi * epsilon0 * epsilonr)) * pi * (charge**2) * &
         integrate_z_cylindrical(n, electrostatic_like_integrand, "all_z")

    calculate_electrostatic_like_term_functional_deriv = 0.0_dp

  end function calculate_electrostatic_like_term_functional_deriv

  function calculate_electrostatic_unlike_term_functional_deriv(n, charge1, charge2)
    real(dp), dimension(:), intent(in) :: n
    real(dp), intent(in) :: charge1
    real(dp), intent(in) :: charge2
    real(dp), dimension(size(n)) :: calculate_electrostatic_unlike_term_functional_deriv

    calculate_electrostatic_unlike_term_functional_deriv(:) = (-1.0_dp / (4.0_dp * pi * epsilon0 * epsilonr)) * pi * charge1 * charge2 * &
         integrate_z_cylindrical(n, electrostatic_unlike_integrand, "z_gteq_hs_diameter")

    calculate_electrostatic_unlike_term_functional_deriv = 0.0_dp

  end function calculate_electrostatic_unlike_term_functional_deriv

  function electrostatic_unlike_integrand(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp) :: electrostatic_unlike_integrand

    integer :: xi_int
    real(dp) :: xi_real

    xi_int = xi - z ! centre xi on z
    xi_real = real(xi_int,dp) * hs_diameter / real(n_discretised_points_z,dp)

    if(abs(xi_int) >= n_discretised_points_z) then
       electrostatic_unlike_integrand = abs(xi_real)
    else
       electrostatic_unlike_integrand = hs_diameter
    end if

  end function electrostatic_unlike_integrand


  function electrostatic_like_integrand(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp) :: electrostatic_like_integrand

    real(dp) :: lambda
    real(dp) :: s

    s = (3.0_dp/(4.0_dp * pi * CURRENT_BULK_BEAD_DENSITY))**(1.0_dp/3.0_dp)
    lambda = sqrt(2.0_dp)/s

    electrostatic_like_integrand = ((abs(z - xi)*hs_diameter/real(n_discretised_points_z,dp)) + (1.0_dp/lambda))

  end function electrostatic_like_integrand


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

  function calculate_n_sbar(array_to_integrate)
    real(dp), dimension(:), intent(in) :: array_to_integrate ! typically n_s, but sometime n_s / n_sbar
    real(dp), dimension(size(array_to_integrate)) :: calculate_n_sbar

    calculate_n_sbar(:) = 0.0_dp

    calculate_n_sbar(:) = (3.0_dp * ( integrate_z_cylindrical(array_to_integrate, n_sbar_integrand, "z_lteq_hs_diameter") ))&
         /(4.0_dp * pi * (hs_diameter**3.0_dp))

  end function calculate_n_sbar

  pure function n_sbar_integrand(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp)            :: n_sbar_integrand

    n_sbar_integrand = 2.0_dp * pi * (hs_diameter**2.0_dp - ((z - xi)*hs_diameter/real(n_discretised_points_z,dp))**2.0_dp) / 2.0_dp

  end function n_sbar_integrand

end module functionalderivatives
