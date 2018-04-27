!Contains routines neccesary to calculate surface forces.
module surfaceforces
  use kinds
  use parameters
  use integratezcylindrical
  use lambdas
  use helpers
  use functionalderivatives
  use constructoligomers
  implicit none
  private

  public :: CalculateGrandPotentialValuePerUnitArea

contains

  !Subroutine that calculates the 
  subroutine CalculateGrandPotentialValuePerUnitArea(n_plus, n_neutral, n_minus, &
       ith_plate_separation, grand_potential_per_unit_area)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in)   :: ith_plate_separation
    real(dp), intent(out) :: grand_potential_per_unit_area

    real(dp), dimension(size(n_plus)) :: n_s
    real(dp), dimension(size(n_plus)) :: n_sbar

    real(dp) :: F_ideal_chain
    real(dp) :: F_van_der_waals
    real(dp) :: F_hard_sphere
    real(dp) :: F_surface_disp
    real(dp) :: F_surface_electro
    real(dp) :: F_electric_like
    real(dp) :: F_electric_unlike

    !If sizes of densities aren't the same then abort
    if((size(n_plus) /= size(n_neutral)) .or. (size(n_plus) /= size(n_minus))) then
       print *, "surfaceforces.f90:CalculateGrandPotentialValue:"
       print *, "sizess of input density arrays must match"
       print *, "coding error...aborting..."
       call abort()
    end if

    n_s = n_plus + n_neutral + n_minus

    n_sbar = calculate_n_sbar(n_s)

    F_ideal_chain = calculate_ideal_chain_term_per_unit_area(n_plus, n_neutral, n_minus, ith_plate_separation, ionic_liquid_name)

    !Note in the following we don't need to integrate over the theta and rho (r) directions
    !as we are only interested in the value PER UNIT AREA.
    F_hard_sphere = calculate_hardsphere_term_per_unit_area(n_s, n_sbar)

    F_van_der_waals = integrate_z_cylindrical(0.5_dp * n_s * calculate_vanderWaals_functional_deriv(n_s), unity_function)

    F_surface_disp = integrate_z_cylindrical(n_s * &
         calculate_surface_dispersion_functional_deriv(ith_plate_separation, size(n_s)), unity_function)

    F_surface_electro = 0.0_dp; F_electric_like = 0.0_dp; F_electric_unlike = 0.0_dp;
    !F_surface_electro = calculate_surface_electrostatic_functional_deriv()

    !F_electric_like = calculate_electrostatic_like_term_functional_deriv()

    !F_electric_unlike = calculate_electrostatic_unlike_term_functional_deriv()

    grand_potential_per_unit_area = F_ideal_chain + F_van_der_waals + F_van_der_waals + F_hard_sphere + &
         F_surface_electro + F_electric_like + F_electric_unlike

  end subroutine CalculateGrandPotentialValuePerUnitArea

  function calculate_ideal_chain_term_per_unit_area(n_plus, n_neutral, n_minus, ith_plate_separation, ionic_liquid_name)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation
    character(len=*), intent(in) :: ionic_liquid_name
    real(dp) :: calculate_ideal_chain_term_per_unit_area

    calculate_ideal_chain_term_per_unit_area = 0.0_dp

    if(trim(ionic_liquid_name) == "SingleNeutralSpheres") then
       calculate_ideal_chain_term_per_unit_area = calculate_single_neutral_sphere_ideal_chain_term(n_neutral, ith_plate_separation)
    else
       print *, "surfaceforces.f90: calculate_ideal_chain_term_per_unit_area:"
       print *, "Unsupported 'ionic_liquid_name' of ", trim(ionic_liquid_name)
    end if

  end function calculate_ideal_chain_term_per_unit_area

  !Calculates the hard sphere contribution to the grand potential.
  !Note that the integration over the angle and radial directions in cylindrical
  !co-ordinates cancel, as we are interested in the value of the term
  !PER UNIT AREA.
  function calculate_hardsphere_term_per_unit_area(n_s, n_sbar)
    real(dp), dimension(:), intent(in) :: n_s
    real(dp), dimension(:), intent(in) :: n_sbar
    real(dp) :: calculate_hardsphere_term_per_unit_area

    integer :: start_z_index
    integer :: end_z_index

    real(dp), dimension(size(n_s)) :: hs_integrand
    hs_integrand = 0.0_dp
    
    call get_allowed_z_values(start_z_index, end_z_index, size(n_s))

    hs_integrand(start_z_index:end_z_index) = (-1.0_dp / beta) * n_sbar(start_z_index:end_z_index) * &
         log( (1.0_dp - (hs_diameter**3)*n_sbar(start_z_index:end_z_index))/(n_sbar(start_z_index:end_z_index)) )

    call setNonCalculatedRegionToZero(hs_integrand)

    calculate_hardsphere_term_per_unit_area = integrate_z_cylindrical(n_sbar * hs_integrand, unity_function)
    
  end function calculate_hardsphere_term_per_unit_area
  
  
end module surfaceforces
