!Contains routines neccesary to compare results with the contact theorem
module contacttheorem
  use kinds
  use helpers
  use parameters
  use discretederivatives
  use integratezcylindrical
  use functionalderivatives
  implicit none
  private

  public :: CalculateNormalPressureFromContactTheorem
  public :: InitialisePotentialAndContactTheoremVariables
  public :: CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation
  public :: MakeContactTheoremAdjustmentFromParticleParticleDispersion

contains

  subroutine CalculateNormalPressureFromContactTheorem(n_plus, n_neutral, n_minus, &
       normal_pressure_left_wall, normal_pressure_right_wall, dispersion_particle_particle_adjust_to_contact_thm)
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    real(dp), intent(out) :: normal_pressure_left_wall
    real(dp), intent(out) :: normal_pressure_right_wall
    real(dp), intent(out) :: dispersion_particle_particle_adjust_to_contact_thm

    real(dp), dimension(size(n_plus)) :: n_s

    real(dp), dimension(size(n_plus)) :: left_wall_dispersion_integrand
    real(dp), dimension(size(n_plus)) :: right_wall_dispersion_integrand

    real(dp) :: maxwell_stress_term_left_wall
    real(dp) :: maxwell_stress_term_right_wall

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

    maxwell_stress_term_left_wall = 0.0_dp
    maxwell_stress_term_right_wall = 0.0_dp

    n_s = n_plus + n_neutral + n_minus
    call get_allowed_z_values(start_z_index, end_z_index, size(n_s))

    call CalculateDerivOfWallTerm(left_wall_dispersion_integrand, right_wall_dispersion_integrand)

    call CalculateDispersionAdjustment(dispersion_particle_particle_adjust_to_contact_thm, n_s)

    call CalculateMaxwellStressTerm(maxwell_stress_term_left_wall, maxwell_stress_term_right_wall)

    !left_wall_dispersion_integrand = 0.0_dp
    !right_wall_dispersion_integrand = 0.0_dp

    !maxwell_stress_term_left_wall = 0.0_dp
    !maxwell_stress_term_right_wall = 0.0_dp

    ! print *, "n_s = ", n_s
    ! print *, "left_wall_dispersion_integrand = ", left_wall_dispersion_integrand
    ! print *, "right_wall_dispersion_integrand = ", right_wall_dispersion_integrand

    ! print *, "n_s left = ", n_s * left_wall_dispersion_integrand
    ! print *, "n_s right = ", n_s * right_wall_dispersion_integrand
    ! call abort()

    !print *, "1", n_s(start_z_index)/beta
    !print *, integrate_z_cylindrical(n_s * left_wall_dispersion_integrand, unity_function)
    !print *, maxwell_stress_term_left_wall
    !call abort()
    
    normal_pressure_left_wall =  (n_s(start_z_index)/beta + integrate_z_cylindrical(n_s * left_wall_dispersion_integrand, unity_function) - maxwell_stress_term_left_wall)
    normal_pressure_right_wall =  (n_s(end_z_index)/beta + integrate_z_cylindrical(n_s * right_wall_dispersion_integrand, unity_function) - maxwell_stress_term_right_wall)

  end subroutine CalculateNormalPressureFromContactTheorem

  subroutine InitialisePotentialAndContactTheoremVariables(grand_potential_per_unit_area, grand_potential_per_unit_area_in_bulk, &
       normal_pressure_left_wall, normal_pressure_right_wall, negative_deriv_of_potential, dispersion_particle_particle_adjust_to_contact_thm, Centre_Separation_Metric)
    real(dp), dimension(:), allocatable :: grand_potential_per_unit_area
    real(dp), dimension(:), allocatable :: grand_potential_per_unit_area_in_bulk
    real(dp), dimension(:), allocatable :: normal_pressure_left_wall
    real(dp), dimension(:), allocatable :: normal_pressure_right_wall
    real(dp), dimension(:), allocatable :: negative_deriv_of_potential
    real(dp), dimension(:), allocatable :: dispersion_particle_particle_adjust_to_contact_thm
    real(dp), dimension(:), allocatable :: Centre_Separation_Metric

    allocate(grand_potential_per_unit_area(size(plate_separations)))
    allocate(grand_potential_per_unit_area_in_bulk(size(plate_separations)))
    allocate(normal_pressure_left_wall(size(plate_separations)))
    allocate(normal_pressure_right_wall(size(plate_separations)))
    allocate(negative_deriv_of_potential(size(plate_separations)))
    allocate(dispersion_particle_particle_adjust_to_contact_thm(size(plate_separations)))
    allocate(Centre_Separation_Metric(size(plate_separations)))
    
    grand_potential_per_unit_area(:) = 0.0_dp
    grand_potential_per_unit_area_in_bulk(:) = 0.0_dp
    normal_pressure_left_wall(:) = 0.0_dp
    normal_pressure_right_wall(:) = 0.0_dp
    negative_deriv_of_potential(:) = 0.0_dp
    dispersion_particle_particle_adjust_to_contact_thm(:) = 0.0_dp
    Centre_Separation_Metric(:) = 0.0_dp
    
  end subroutine InitialisePotentialAndContactTheoremVariables

  function CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation(grand_potential_per_unit_area)
    real(dp), dimension(:), intent(in) :: grand_potential_per_unit_area
    real(dp), dimension(size(grand_potential_per_unit_area)) :: CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation

    CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation = -1.0_dp * calculate_central_difference(grand_potential_per_unit_area)

  end function CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation

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
       
       left_wall_term(ij) = (2.0_dp * pi * epsilon_LJ_particle_wall * ((0.4_dp * (hs_diameter/distance_from_left_wall)**9) - (hs_diameter/distance_from_left_wall)**3))/distance_from_left_wall
       right_wall_term(ij) = (2.0_dp * pi * epsilon_LJ_particle_wall * ((0.4_dp * (hs_diameter/distance_from_right_wall)**9) - (hs_diameter/distance_from_right_wall)**3))/distance_from_right_wall
    end do


  end subroutine CalculateDerivOfWallTerm

  subroutine CalculateDispersionAdjustment(dispersion_particle_particle_adjust_to_contact_thm, n_s)
    real(dp), intent(out) :: dispersion_particle_particle_adjust_to_contact_thm
    real(dp), dimension(:), intent(in) :: n_s

    dispersion_particle_particle_adjust_to_contact_thm = integrate_z_cylindrical(0.5_dp * n_s * calculate_vanderWaals_functional_deriv(n_s), unity_function)
    
  end subroutine CalculateDispersionAdjustment

  subroutine MakeContactTheoremAdjustmentFromParticleParticleDispersion(normal_pressure_left_wall, normal_pressure_right_wall, dispersion_particle_particle_adjust_to_contact_thm)
    real(dp), dimension(:) :: normal_pressure_left_wall
    real(dp), dimension(:) :: normal_pressure_right_wall
    real(dp), dimension(:) :: dispersion_particle_particle_adjust_to_contact_thm

    integer :: ij

    dispersion_particle_particle_adjust_to_contact_thm(:) = CalculateNegativeDerivOfPotentialPerUnitAreaWRTSeparation(dispersion_particle_particle_adjust_to_contact_thm)
    
    do ij = 1, size(normal_pressure_left_wall)
       normal_pressure_left_wall(ij) =  normal_pressure_left_wall(ij) - dispersion_particle_particle_adjust_to_contact_thm(ij)
       normal_pressure_right_wall(ij) =  normal_pressure_right_wall(ij) - dispersion_particle_particle_adjust_to_contact_thm(ij)
    end do

  end subroutine MakeContactTheoremAdjustmentFromParticleParticleDispersion

end module contacttheorem
