!Contains routines neccesary to calculate surface forces.
module surfaceforces
  use kinds
  use parameters
  use integratezcylindrical
  use lambdas
  use helpers
  use functionalderivatives
  use constructoligomers
  use excessenergyfunctionalparameters
  implicit none
  private

  public :: CalculateGrandPotentialValuePerUnitArea

contains

  !Subroutine that calculates the 
  subroutine CalculateGrandPotentialValuePerUnitArea(ith_plate_separation, grand_potential_per_unit_area, &
       size_of_ns_array, n_plus, n_neutral, n_minus)
    integer, intent(in)   :: ith_plate_separation
    real(dp), intent(out) :: grand_potential_per_unit_area
    integer, intent(in) :: size_of_ns_array
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus

    real(dp), dimension(size_of_ns_array) :: n_s
    real(dp), dimension(size_of_ns_array) :: n_sbar
    
    real(dp), dimension(size_of_ns_array) :: n_plus_input
    real(dp), dimension(size_of_ns_array) :: n_neutral_input
    real(dp), dimension(size_of_ns_array) :: n_minus_input

    real(dp) :: chemical_potential_term

    real(dp) :: potential_per_unit_area_not_in_bulk
    real(dp) :: potential_per_unit_area_in_bulk

    real(dp) :: F_ideal_chain
    real(dp) :: F_van_der_waals

    real(dp) :: F_hard_sphere
    real(dp) :: F_surface_disp
    real(dp) :: F_surface_electro
    real(dp) :: F_electric_like
    real(dp) :: F_electric_unlike

    F_ideal_chain = 0.0_dp
    F_van_der_waals = 0.0_dp
    F_hard_sphere = 0.0_dp
    F_surface_disp = 0.0_dp
    F_surface_electro = 0.0_dp
    F_electric_like = 0.0_dp
    F_electric_unlike = 0.0_dp

    
    n_s(:) = n_plus(:) + n_neutral(:) + n_minus(:)
    n_sbar(:) = calculate_n_sbar(n_s(:))

    n_plus_input(:) = n_plus(:)
    n_neutral_input(:) = n_neutral(:)
    n_minus_input(:) = n_minus(:)

    F_surface_disp = integrate_z_cylindrical(n_s * &
         calculate_surface_dispersion_functional_deriv(ith_plate_separation, size(n_s)), unity_function)

    F_ideal_chain = calculate_ideal_chain_term_per_unit_area(size(n_plus_input), n_plus_input, n_neutral_input, n_minus_input, ith_plate_separation)

    F_hard_sphere = calculate_hardsphere_term_per_unit_area(n_s, n_sbar)

    potential_per_unit_area_not_in_bulk = (F_ideal_chain + F_van_der_waals + F_surface_disp + F_hard_sphere + &
         F_surface_electro + F_electric_like + F_electric_unlike)

    !potential_per_unit_area_not_in_bulk = 0.0_dp

    ! !Now calculate the value of the potential in the bulk
    ! F_ideal_chain = 0.0_dp
    ! F_van_der_waals = 0.0_dp
    ! F_hard_sphere = 0.0_dp
    ! F_surface_disp = 0.0_dp
    ! F_surface_electro = 0.0_dp
    ! F_electric_like = 0.0_dp
    ! F_electric_unlike = 0.0_dp

    ! n_s(:) = bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads
    ! n_sbar(:) = n_s(:)

    ! !Note: Only passing onw input parameters calculates the value in the bulk.
    ! F_ideal_chain = calculate_ideal_chain_term_per_unit_area(size(n_plus_input))

    ! !Note in the following we don't need to integrate over the theta and rho (r) directions
    ! !as we are only interested in the value PER UNIT AREA.
    ! F_hard_sphere = calculate_hardsphere_term_per_unit_area(n_s, n_sbar)

    ! potential_per_unit_area_in_bulk = (F_ideal_chain + F_van_der_waals + F_surface_disp + F_hard_sphere + &
    !      F_surface_electro + F_electric_like + F_electric_unlike)



    chemical_potential_term = calculate_chemical_potential_term(n_plus_input, n_neutral_input, n_minus_input, ith_plate_separation)

    !F_van_der_waals = integrate_z_cylindrical(0.5_dp * n_s * calculate_vanderWaals_functional_deriv(n_s), unity_function)
    !F_surface_electro = calculate_surface_electrostatic_functional_deriv()
    !F_electric_like = calculate_electrostatic_like_term_functional_deriv()
    !F_electric_unlike = calculate_electrostatic_unlike_term_functional_deriv()

    !grand_potential_per_unit_area = potential_per_unit_area_not_in_bulk - potential_per_unit_area_in_bulk
    grand_potential_per_unit_area = potential_per_unit_area_not_in_bulk - chemical_potential_term

    
    !grand_potential_per_unit_area = potential_per_unit_area_in_bulk 
    !print *, "grand_potential_per_unit_area = ", grand_potential_per_unit_area
    !call abort()

  end subroutine CalculateGrandPotentialValuePerUnitArea


  function calculate_chemical_potential_term(n_plus, n_neutral, n_minus, ith_plate_separation)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    integer, intent(in) :: ith_plate_separation

    real(dp) :: calculate_chemical_potential_term

    if(trim(ionic_liquid_name) == "SingleNeutralSpheres") then
       calculate_chemical_potential_term = calculate_chemical_potential_term_neutral_spheres(n_plus, n_neutral, n_minus, ith_plate_separation)
    else if(trim(ionic_liquid_name) == "NeutralDimers") then
       !calculate_chemical_potential_term = calculate_chemical_potential_term_neutral_dimers()
    else
       print *, "surfaceforces.f90: CalculateChemicalPotentialTerm:"
       print *, "Unsupported 'ionic_liquid_name' of ", trim(ionic_liquid_name)
    end if


  end function Calculate_Chemical_Potential_Term



  !Calculates the value of the ideal chain term per unit area.  Note that all the input arguments are optional.
  !In the case that none are present then calculate the bulk value. (by setting lambda_b - lambda = 0).
  !In the case that they're all present calulculate the value based on the input densities.
  !In any other case, abort.
  function calculate_ideal_chain_term_per_unit_area(n_points, n_plus, n_neutral, n_minus, ith_plate_separation)
    integer, intent(in) :: n_points !Number of points to calculate
    real(dp), dimension(:), intent(in), optional :: n_plus
    real(dp), dimension(:), intent(in), optional :: n_neutral
    real(dp), dimension(:), intent(in), optional :: n_minus
    integer, intent(in), optional :: ith_plate_separation
    real(dp) :: calculate_ideal_chain_term_per_unit_area

    real(dp), dimension(n_points) :: lambda_plus, lambda_neutral, lambda_minus
    real(dp), dimension(n_points) :: n_plus_input, n_neutral_input, n_minus_input

    calculate_ideal_chain_term_per_unit_area = 0.0_dp

    lambda_plus(:) = 0.0_dp
    lambda_neutral(:) = 0.0_dp
    lambda_minus(:) = 0.0_dp

    if(present(n_plus) .and. present(n_neutral) .and. present(n_minus) .and. present(ith_plate_separation)) then

       if((size(n_plus) /= size(n_neutral)) .or. (size(n_plus) /= size(n_minus))) then
          print *, "surfaceforces.f90.f90: calculate_ideal_chain_term_per_unit_area"
          print *, "Size mismatch..."
          print *, "(size(n_plus) /= size(n_neutral)) .or. (size(n_plus) /= size(n_minus))"
          print *, "coding error...aborting..."
          call abort()
       end if

       call CalculateLambdasDifference(lambda_plus, n_plus, lambda_neutral, n_neutral, lambda_minus, n_minus, ith_plate_separation)
       n_plus_input(:) = n_plus(:)
       n_neutral_input(:) = n_neutral(:)
       n_minus_input(:) = n_minus(:)
    else if(present(n_plus) .or. present(n_neutral) .or. present(n_minus) .or. present(ith_plate_separation)) then
       print *, "surfaceforces.f90: calculate_ideal_chain_term_per_unit_area"
       print *, "routine has optional argument list.  They must be all present or all not present."
       print *, "This rule has been violated...coding error...aborting..."
       call abort()
    else !else we want to calculate the value in the bulk so we leave lambda_b - lambda = 0
       n_plus_input(:) = bulk_density_positive_beads
       n_neutral_input(:) = bulk_density_neutral_beads
       n_minus_input(:) = bulk_density_negative_beads
    end if

    if(trim(ionic_liquid_name) == "SingleNeutralSpheres") then
       calculate_ideal_chain_term_per_unit_area = calculate_single_neutral_sphere_ideal_chain_term(n_neutral_input)
    else if(trim(ionic_liquid_name) == "NeutralDimers") then
       calculate_ideal_chain_term_per_unit_area = calculate_neutral_dimers_ideal_chain_term(lambda_neutral)       
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
    real(dp), dimension(size(n_s)) :: n_solvent

    hs_integrand(:) = 0.0_dp
    n_solvent(:) = 0.0_dp

    call get_allowed_z_values(start_z_index, end_z_index, size(n_s))

    hs_integrand(start_z_index:end_z_index) = (0.5_dp / beta) * n_s(start_z_index:end_z_index) * &
         GetAEx(n_sbar(start_z_index:end_z_index), n_solvent(start_z_index:end_z_index), a_term_index)

    !call setNonCalculatedRegionToZero(hs_integrand)

    calculate_hardsphere_term_per_unit_area = integrate_z_cylindrical(hs_integrand, unity_function)

  end function calculate_hardsphere_term_per_unit_area

end module surfaceforces
