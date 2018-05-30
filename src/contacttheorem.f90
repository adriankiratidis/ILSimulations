!Contains routines neccesary to compare results with the contact theorem
module contacttheorem
  use kinds
  use helpers
  use parameters
  use discretederivatives
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

    integer :: start_z_index
    integer :: end_z_index
    
    !First check the sizes are the same
    if( (size(n_plus) /= size(n_neutral)) .or. (size(n_plus) /= size(n_minus)) ) then
       print *, "contacttheorem.f90:CalculateNormalPressureFromContactTheorem: "
       print *, "(size(n_plus) /= size(n_neutral)) .or. (size(n_plus) /= size(n_minus))"
       print *, "Density array size mismatch...aborting..."
       call abort()
    end if

    n_s = n_plus + n_neutral + n_minus
    call get_allowed_z_values(start_z_index, end_z_index, size(n_s))
    
    normal_pressure_left_wall =  n_s(start_z_index) / beta
    normal_pressure_right_wall =  n_s(end_z_index) / beta
    
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


end module contacttheorem
