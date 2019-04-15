module diagnostics
  use kinds
  use parameters
  use functionalderivatives
  implicit none
  private

  public :: GetLambdaContributionsByTerm
            
contains

  subroutine GetLambdaContributionsByTerm(ith_plate_separation, n_plus, n_neutral, n_minus, lambda_plus, lambda_neutral, lambda_minus, lambda_cation_centre, lambda_anion_centre)
    integer, intent(in) :: ith_plate_separation
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(in) :: lambda_cation_centre
    real(dp), dimension(:), intent(in) :: lambda_anion_centre

    real(dp), dimension(size(n_plus)) :: hs_term
    real(dp), dimension(size(n_plus)) :: electrostatic_like_term
    real(dp), dimension(size(n_plus)) :: electrostatic_unlike_term
    real(dp), dimension(size(n_plus)) :: particle_particle_vanDerWaals_term
    real(dp), dimension(size(n_plus)) :: particle_surface_dispersion_term
    real(dp), dimension(size(n_plus)) :: particle_surface_electrostatic_term

    real(dp), dimension(size(n_plus)) :: n_s
    real(dp), dimension(size(n_plus)) :: n_s_bulk

    real(dp), dimension(size(n_plus)) :: n_plus_in, n_minus_in
    real(dp), dimension(size(n_plus)) :: TOTAL
    
    n_s(:) = n_plus + n_neutral + n_minus
    n_s_bulk(:) = bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads
    n_plus_in(:) = bulk_density_positive_beads
    n_minus_in(:) = bulk_density_negative_beads

    hs_term(:) = calculate_hardsphere_functional_deriv(n_s_bulk, .true.) - calculate_hardsphere_functional_deriv(n_s, .false.)

    particle_particle_vanDerWaals_term(:) = &
         (-8.0_dp * epsilon_LJ_particle_particle * (hs_diameter**6) * pi * n_s_bulk(:)) * ((2.0_dp/(3.0_dp*(hs_diameter**3)))) &
         - calculate_vanderWaals_functional_deriv(n_s)

    particle_surface_dispersion_term(:) = -calculate_surface_dispersion_functional_deriv(ith_plate_separation, size(n_plus))

    particle_surface_electrostatic_term(:) = & !-calculate_surface_electrostatic_functional_deriv(size(n_plus), negative_bead_charge) &
         - calculate_surface_electrostatic_functional_deriv(size(n_plus), positive_bead_charge)

    electrostatic_like_term(:) = calculate_electrostatic_like_term_functional_deriv(n_plus_in, positive_bead_charge, .true.) &
         - calculate_electrostatic_like_term_functional_deriv(n_plus, positive_bead_charge, .false.)
    !+ calculate_electrostatic_like_term_functional_deriv(n_minus_in, negative_bead_charge, .true.) &
    !- calculate_electrostatic_like_term_functional_deriv(n_minus, negative_bead_charge, .false.)

    electrostatic_unlike_term(:) = &!calculate_electrostatic_unlike_term_functional_deriv(n_minus_in, positive_bead_charge, negative_bead_charge, .true.) &
                                !- calculate_electrostatic_unlike_term_functional_deriv(n_minus, positive_bead_charge, negative_bead_charge, .false.) &
         + calculate_electrostatic_unlike_term_functional_deriv(n_plus_in, positive_bead_charge, negative_bead_charge, .true.) &
         - calculate_electrostatic_unlike_term_functional_deriv(n_plus, positive_bead_charge, negative_bead_charge, .false.)

    TOTAL(:) = (hs_term(:) + particle_particle_vanDerWaals_term(:) + particle_surface_dispersion_term(:) + particle_surface_electrostatic_term(:) + &
         electrostatic_like_term(:) + electrostatic_unlike_term(:)) 

    
    print *, "******************************************************"
    print *, "******************************************************"
    print *, "******************************************************"
    print *, "******************************************************"
    print *, "Contributions to lambda plus"
    print *, "hs_term = ", hs_term(:)/TOTAL(:) * 100
    print *, ""
    print *, "particle_particle_vanDerWaals_term = ", particle_particle_vanDerWaals_term(:)/TOTAL(:) * 100
    print *, ""
    print *, "particle_surface_dispersion_term = ", particle_surface_dispersion_term(:)/TOTAL(:) * 100
    print *, ""
    print *, "particle_surface_electrostatic_term = ", particle_surface_electrostatic_term(:)/TOTAL(:) * 100
    print *, ""
    print *, "electrostatic particle_particle term = ", (electrostatic_like_term(:) + electrostatic_unlike_term(:))/TOTAL(:) * 100
    
    
    !print *, ""
    !print *, "electrostatic_like_term = ", electrostatic_like_term(:)/TOTAL(:) * 100
    !print *, ""
    !print *, "electrostatic_unlike_term = ", electrostatic_unlike_term(:)/TOTAL(:) * 100

    print *, "centre potential adjustment =  ", lambda_cation_centre(:) / lambda_plus(:) * 100
    call abort()



  end subroutine GetLambdaContributionsByTerm

end module diagnostics
