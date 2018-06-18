!Contains routines neccesary to compare results with the contact theorem
module contacttheorem
  use kinds
  use helpers
  use parameters
  use discretederivatives
  use integratezcylindrical
  implicit none
  private

  public :: CalculateNormalPressureFromContactTheorem
  public :: InitialisePotentialAndContactTheoremVariables
  public :: CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation

contains

  subroutine CalculateNormalPressureFromContactTheorem(n_plus, n_neutral, n_minus, &
       normal_pressure_left_wall, normal_pressure_right_wall)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    real(dp), intent(out) :: normal_pressure_left_wall
    real(dp), intent(out) :: normal_pressure_right_wall

    real(dp), dimension(size(n_plus)) :: n_s

    real(dp), dimension(size(n_plus)) :: left_wall_dispersion_integrand
    real(dp), dimension(size(n_plus)) :: right_wall_dispersion_integrand

    integer :: start_z_index
    integer :: end_z_index

    !First check the sizes are the same
    if( (size(n_plus) /= size(n_neutral)) .or. (size(n_plus) /= size(n_minus)) ) then
       print *, "contacttheorem.f90:CalculateNormalPressureFromContactTheorem: "
       print *, "(size(n_plus) /= size(n_neutral)) .or. (size(n_plus) /= size(n_minus))"
       print *, "Density array size mismatch...aborting..."
       call abort()
    end if

    left_wall_dispersion_integrand = 0.0_dp
    right_wall_dispersion_integrand = 0.0_dp
    n_s = n_plus + n_neutral + n_minus
    call get_allowed_z_values(start_z_index, end_z_index, size(n_s))

    call CalculateDerivOfWallTerm(left_wall_dispersion_integrand, right_wall_dispersion_integrand)

    !left_wall_dispersion_integrand = 0.0_dp
    !right_wall_dispersion_integrand = 0.0_dp


    !right_wall_dispersion_integrand = 0.0_dp
    !left_wall_dispersion_integrand = 0.0_dp

    ! print *, "n_s = ", n_s
    ! print *, "left_wall_dispersion_integrand = ", left_wall_dispersion_integrand
    ! print *, "right_wall_dispersion_integrand = ", right_wall_dispersion_integrand

    ! print *, "n_s left = ", n_s * left_wall_dispersion_integrand
    ! print *, "n_s right = ", n_s * right_wall_dispersion_integrand
    ! call abort()

    normal_pressure_left_wall =  (n_s(start_z_index) + integrate_z_cylindrical(n_s * left_wall_dispersion_integrand, unity_function)) / beta
    normal_pressure_right_wall =  (n_s(end_z_index) + integrate_z_cylindrical(n_s * right_wall_dispersion_integrand, unity_function)) / beta

  end subroutine CalculateNormalPressureFromContactTheorem

  subroutine InitialisePotentialAndContactTheoremVariables(grand_potential_per_unit_area, grand_potential_per_unit_area_in_bulk, &
       normal_pressure_left_wall, normal_pressure_right_wall, negative_deriv_of_potential)
    real(dp), dimension(:), allocatable :: grand_potential_per_unit_area
    real(dp), dimension(:), allocatable :: grand_potential_per_unit_area_in_bulk
    real(dp), dimension(:), allocatable :: normal_pressure_left_wall
    real(dp), dimension(:), allocatable :: normal_pressure_right_wall
    real(dp), dimension(:), allocatable :: negative_deriv_of_potential

    allocate(grand_potential_per_unit_area(size(plate_separations)))
    allocate(grand_potential_per_unit_area_in_bulk(size(plate_separations)))
    allocate(normal_pressure_left_wall(size(plate_separations)))
    allocate(normal_pressure_right_wall(size(plate_separations)))
    allocate(negative_deriv_of_potential(size(plate_separations)))

    grand_potential_per_unit_area(:) = 0.0_dp
    grand_potential_per_unit_area_in_bulk(:) = 0.0_dp
    normal_pressure_left_wall(:) = 0.0_dp
    normal_pressure_right_wall(:) = 0.0_dp
    negative_deriv_of_potential(:) = 0.0_dp

  end subroutine InitialisePotentialAndContactTheoremVariables

  subroutine CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation(grand_potential_per_unit_area, &
       negative_deriv_of_potential)
    real(dp), dimension(:), intent(in) :: grand_potential_per_unit_area
    real(dp), dimension(:), intent(out) :: negative_deriv_of_potential

    negative_deriv_of_potential = -1.0_dp * calculate_central_difference(grand_potential_per_unit_area)

  end subroutine CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation

  subroutine CalculateMaxwellStressTerm(maxwell_stress_term_left_wall, maxwell_stress_term_right_wall)
    real(dp), intent(out) :: maxwell_stress_term_left_wall
    real(dp), intent(out) :: maxwell_stress_term_right_wall

    maxwell_stress_term_left_wall = 2.0_dp * pi * (surface_charge_density_left_wall**2)/ (epsilonr * epsilon0)
    maxwell_stress_term_right_wall = 2.0_dp * pi * (surface_charge_density_right_wall**2)/ (epsilonr * epsilon0)

  end subroutine CalculateMaxwellStressTerm

  subroutine CalculateDerivOfWallTerm(left_wall_term, right_wall_term)
    real(dp), dimension(:) :: left_wall_term
    real(dp), dimension(:) :: right_wall_term

    integer :: ij
    real(dp) :: distance_from_left_wall
    real(dp) :: distance_from_right_wall

    do ij = 2, size(left_wall_term) - 1
       
       distance_from_left_wall = (ij - 1)*hs_diameter/real(n_discretised_points_z,dp)
       distance_from_right_wall = (size(left_wall_term) - ij)*hs_diameter/real(n_discretised_points_z,dp)

       !print *, "number = ", distance_from_left_wall, distance_from_right_wall, distance_from_right_wall + distance_from_left_wall
       
       left_wall_term(ij) = (2.0_dp * pi * epsilon_LJ * ((0.4_dp * (hs_diameter/distance_from_left_wall)**9) - (hs_diameter/distance_from_left_wall)**3))/distance_from_left_wall
       right_wall_term(ij) = (2.0_dp * pi * epsilon_LJ * ((0.4_dp * (hs_diameter/distance_from_right_wall)**9) - (hs_diameter/distance_from_right_wall)**3))/distance_from_right_wall
    end do


  end subroutine CalculateDerivOfWallTerm



end module contacttheorem
