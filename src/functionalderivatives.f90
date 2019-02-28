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
  use discretederivatives
  implicit none
  private

  public :: calculate_vanderWaals_functional_deriv
  public :: calculate_surface_dispersion_functional_deriv
  public :: calculate_hardsphere_functional_deriv
  public :: calculate_hard_sphere_functional_deriv !This one is the latest one
  public :: calculate_surface_electrostatic_functional_deriv

  public :: calculate_electrostatic_like_term_functional_deriv
  public :: calculate_electrostatic_unlike_term_functional_deriv

  !public :: CalculateHSEndAndNonEndHSFunctionalDeriv
  
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

  end function calculate_hardsphere_functional_deriv

  ! subroutine CalculateHSEndAndNonEndHSFunctionalDeriv(n_s, lambda_hs_end_cation, n_hs_end_cation, lambda_hs_nonend_cation, n_hs_nonend_cation, &
  !      lambda_hs_end_anion, n_hs_end_anion, lambda_hs_nonend_anion, n_hs_nonend_anion, calculate_bulk)
  !   real(dp), dimension(:), intent(in) :: n_s

  !   real(dp), dimension(:), intent(out) :: lambda_hs_end_cation
  !   real(dp), dimension(:), intent(in)  :: n_hs_end_cation
  !   real(dp), dimension(:), intent(out) :: lambda_hs_nonend_cation
  !   real(dp), dimension(:), intent(in)  :: n_hs_nonend_cation

  !   real(dp), dimension(:), intent(out) :: lambda_hs_end_anion
  !   real(dp), dimension(:), intent(in)  :: n_hs_end_anion
  !   real(dp), dimension(:), intent(out) :: lambda_hs_nonend_anion
  !   real(dp), dimension(:), intent(in)  :: n_hs_nonend_anion

  !   logical, intent(in) :: calculate_bulk

  !   real(dp), dimension(size(n_s)) :: n_mbar, n_sbar
  !   real(dp), dimension(size(n_s)) :: integrand1, integrand2
  !   real(dp), dimension(size(n_s)) :: integral1_end, integral1_nonend, integral2_end, integral2_nonend

  !   integer :: start_z_index
  !   integer :: end_z_index

  !   integrand1(:) = 0.0_dp
  !   integrand2(:) = 0.0_dp

  !   integral1_end(:) = 0.0_dp
  !   integral1_nonend(:) = 0.0_dp

  !   integral2_end(:) = 0.0_dp
  !   integral2_nonend(:) = 0.0_dp

  !   n_mbar(:) = 0.0_dp
  !   n_sbar(:) = 0.0_dp

  !   n_mbar(:) = calculate_n_sbar(n_s(:))

  !   !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
  !   call get_allowed_z_values(start_z_index, end_z_index, size(n_s))

  !   !print *, "n_mbar", n_mbar, calculate_bulk

  !   !integral1(:) =  integrand1!calculate_n_sbar(integrand1(:))
  !   !integral2(:) =  integrand2!calculate_n_sbar(integrand2(:))
  !   !print *, integrand1
  !   !call abort()

  !   !Replace
  !   !GetYMix(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index)
  !   !with
  !   !GetYMix(n_hs_end_cation(start_z_index:end_z_index) + n_hs_nonend_cation(start_z_index:end_z_index), n_hs_end_anion(start_z_index:end_z_index) + n_hs_nonend_anion(start_z_index:end_z_index)

  !   if(calculate_bulk) then            


  !      integrand1(start_z_index:end_z_index) = ((4.0_dp * pi * (hs_diameter**3))/3.0_dp) * &
  !           GetAExDerivIntegrand(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 1)
  !      integrand2(start_z_index:end_z_index) = ((4.0_dp * pi * (hs_diameter**3))/3.0_dp) * &
  !           GetAExDerivIntegrand(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2)

  !      !First we do the cation
  !      lambda_hs_end_cation(start_z_index:end_z_index) = (1.0_dp / (beta * num_end_monomers_cation)) * (&
  !           GetAEx(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2) + &
  !           ((n_hs_end_cation(start_z_index:end_z_index) * integrand2(start_z_index:end_z_index)))) + &
  !           (n_hs_nonend_cation(start_z_index:end_z_index) * (GetYMix(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_cation)/(r_cation - num_end_monomers_cation)) * &
  !           ((integrand2(start_z_index:end_z_index) * n_hs_nonend_cation(start_z_index:end_z_index)) - &
  !           (integrand1(start_z_index:end_z_index) * n_hs_nonend_cation(start_z_index:end_z_index))))/beta

  !      lambda_hs_nonend_cation(start_z_index:end_z_index) = (1.0_dp / (beta * num_end_monomers_cation)) * (&
  !           ((n_hs_end_cation(start_z_index:end_z_index) * integrand2(start_z_index:end_z_index)))) + &
  !           (n_hs_nonend_cation(start_z_index:end_z_index) * (GetYMix(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_cation)/(r_cation - num_end_monomers_cation)) * &
  !           ((integrand2(start_z_index:end_z_index) * n_hs_nonend_cation(start_z_index:end_z_index)) - &
  !           (integrand1(start_z_index:end_z_index) * n_hs_nonend_cation(start_z_index:end_z_index))))/beta + &
  !           ((GetYMix(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_cation)/(r_cation - num_end_monomers_cation)) * GetAEx(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2) - &
  !           GetAEx(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 1))/beta

  !      !Now we do the anion
  !      lambda_hs_end_anion(start_z_index:end_z_index) = (1.0_dp / (beta * num_end_monomers_anion)) * (&
  !           GetAEx(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2) + &
  !           ((n_hs_end_anion(start_z_index:end_z_index) * integrand2(start_z_index:end_z_index)))) + &
  !           (n_hs_nonend_anion(start_z_index:end_z_index) * (GetYMix(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_anion)/(r_anion - num_end_monomers_anion)) * &
  !           ((integrand2(start_z_index:end_z_index) * n_hs_nonend_anion(start_z_index:end_z_index)) - &
  !           (integrand1(start_z_index:end_z_index) * n_hs_nonend_anion(start_z_index:end_z_index))))/beta

  !      lambda_hs_nonend_anion(start_z_index:end_z_index) = (1.0_dp / (beta * num_end_monomers_anion)) * (&
  !           ((n_hs_end_anion(start_z_index:end_z_index) * integrand2(start_z_index:end_z_index)))) + &
  !           (n_hs_nonend_anion(start_z_index:end_z_index) * (GetYMix(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_anion)/(r_anion - num_end_monomers_anion)) * &
  !           ((integrand2(start_z_index:end_z_index) * n_hs_nonend_anion(start_z_index:end_z_index)) - &
  !           (integrand1(start_z_index:end_z_index) * n_hs_nonend_anion(start_z_index:end_z_index))))/beta + &
  !           ((GetYMix(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_anion)/(r_anion - num_end_monomers_anion)) * GetAEx(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2) - &
  !           GetAEx(n_s(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 1))/beta       


  !   else

  !      integrand1(start_z_index:end_z_index) = ((4.0_dp * pi * (hs_diameter**3))/3.0_dp) * &
  !           GetAExDerivIntegrand(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 1)
  !      integrand2(start_z_index:end_z_index) = ((4.0_dp * pi * (hs_diameter**3))/3.0_dp) * &
  !           GetAExDerivIntegrand(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2)

  !      integral2_end(:) = calculate_n_sbar(n_hs_end_cation(:) * integrand2(:))
  !      integral2_nonend(:) = calculate_n_sbar(n_hs_nonend_cation(:) * integrand2(:))

  !      integral1_end(:) = calculate_n_sbar(n_hs_end_cation(:) * integrand1(:))
  !      integral1_nonend(:) = calculate_n_sbar(n_hs_nonend_cation(:) * integrand1(:))

  !      !First we do the cation
  !      lambda_hs_end_cation(start_z_index:end_z_index) = (1.0_dp / (beta * num_end_monomers_cation)) * (&
  !           GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2) + &
  !           (integral2_end(start_z_index:end_z_index))) + &
  !           ((GetYMix(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_cation)/(r_cation - num_end_monomers_cation)) * &
  !           (integral2_nonend(start_z_index:end_z_index) - &
  !           integral1_nonend(start_z_index:end_z_index)))/beta

  !      lambda_hs_nonend_cation(start_z_index:end_z_index) = (1.0_dp / (beta * num_end_monomers_cation)) * (&
  !           (integral2_end(start_z_index:end_z_index))) + &
  !           ((GetYMix(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_cation)/(r_cation - num_end_monomers_cation)) * &
  !           (integral2_nonend(start_z_index:end_z_index) - &
  !           integral1_nonend(start_z_index:end_z_index)))/beta + &
  !           ((GetYMix(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_cation)/(r_cation - num_end_monomers_cation)) * GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2) - &
  !           GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 1))/beta


  !      integral2_end(:) = calculate_n_sbar(n_hs_end_anion(:) * integrand2(:))
  !      integral2_nonend(:) = calculate_n_sbar(n_hs_nonend_anion(:) * integrand2(:))

  !      integral1_end(:) = calculate_n_sbar(n_hs_end_anion(:) * integrand1(:))
  !      integral1_nonend(:) = calculate_n_sbar(n_hs_nonend_anion(:) * integrand1(:))

  !      !Now we do the anion
  !      lambda_hs_end_anion(start_z_index:end_z_index) = (1.0_dp / (beta * num_end_monomers_anion)) * (&
  !           GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2) + &
  !           (integral2_end(start_z_index:end_z_index))) + &
  !           ((GetYMix(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_anion)/(r_anion - num_end_monomers_anion)) * &
  !           (integral2_nonend(start_z_index:end_z_index) - &
  !           integral1_nonend(start_z_index:end_z_index)))/beta

  !      lambda_hs_nonend_anion(start_z_index:end_z_index) = (1.0_dp / (beta * num_end_monomers_anion)) * (&
  !           (integral2_end(start_z_index:end_z_index))) + &
  !           ((GetYMix(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_anion)/(r_anion - num_end_monomers_anion)) * &
  !           (integral2_nonend(start_z_index:end_z_index) - &
  !           integral1_nonend(start_z_index:end_z_index)))/beta + &
  !           ((GetYMix(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_anion)/(r_anion - num_end_monomers_anion)) * GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2) - &
  !           GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 1))/beta

  !   end if

  !   ! lambda_hs_end_cation = 0.0_dp
  !   ! lambda_hs_nonend_cation = 0.0_dp
  !   ! lambda_hs_end_anion = 0.0_dp
  !   ! lambda_hs_nonend_anion = 0.0_dp

  ! end subroutine CalculateHSEndAndNonEndHSFunctionalDeriv

  function calculate_hard_sphere_functional_deriv(n_s, n_deriv_wrt, n_cation_end, n_cation_nonend, n_anion_end, n_anion_nonend, calculate_bulk)
    real(dp), dimension(:), intent(in) :: n_s
    real(dp), dimension(:), intent(in) :: n_deriv_wrt
    real(dp), dimension(:), intent(in) :: n_cation_end
    real(dp), dimension(:), intent(in) :: n_cation_nonend
    real(dp), dimension(:), intent(in) :: n_anion_end
    real(dp), dimension(:), intent(in) :: n_anion_nonend
    logical, intent(in)                :: calculate_bulk
    real(dp), dimension(size(n_deriv_wrt)) :: calculate_hard_sphere_functional_deriv

    real(dp), dimension(size(n_s)) :: n_mbar, n_sbar, n_mbar_bulk
    real(dp), dimension(size(n_s)) :: integrand1, integrand2

    real(dp), dimension(size(n_cation_end)) :: n_cation_end_bulk
    real(dp), dimension(size(n_cation_nonend)) :: n_cation_nonend_bulk
    real(dp), dimension(size(n_anion_end)) :: n_anion_end_bulk
    real(dp), dimension(size(n_anion_nonend)) :: n_anion_nonend_bulk

    real(dp), dimension(size(n_s)) :: ratio_cation_end, ratio_cation_nonend
    real(dp), dimension(size(n_s)) :: ratio_anion_end, ratio_anion_nonend
    real(dp), dimension(size(n_s)) :: deriv_wrt

    integer :: start_z_index
    integer :: end_z_index
    integer :: ij

    ratio_cation_end(:) = 0.0_dp
    ratio_cation_nonend(:) = 0.0_dp
    ratio_anion_end(:) = 0.0_dp
    ratio_anion_nonend(:) = 0.0_dp

    n_cation_end_bulk(:) = 0.0_dp
    n_cation_nonend_bulk(:) = 0.0_dp
    n_anion_end_bulk(:) = 0.0_dp
    n_anion_nonend_bulk(:) = 0.0_dp

    deriv_wrt(:) = 0.0_dp

    n_mbar(:) = 0.0_dp
    n_sbar(:) = 0.0_dp

    n_mbar(:) = calculate_n_sbar(n_s(:))       

    call get_allowed_z_values(start_z_index, end_z_index, size(n_s))

    calculate_hard_sphere_functional_deriv(:) = 0.0_dp

    if(calculate_bulk) then

       !First calculate the bulk properties
       n_cation_end_bulk(start_z_index:end_z_index) = n_plus_cation_end_bulk + n_neutral_cation_end_bulk + n_minus_cation_end_bulk
       n_cation_nonend_bulk(start_z_index:end_z_index) = n_plus_cation_nonend_bulk + n_neutral_cation_nonend_bulk + n_minus_cation_nonend_bulk
       n_anion_end_bulk(start_z_index:end_z_index) = n_plus_anion_end_bulk + n_neutral_anion_end_bulk + n_minus_anion_end_bulk
       
       n_anion_nonend_bulk(start_z_index:end_z_index) = bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads
       n_anion_nonend_bulk(start_z_index:end_z_index) = n_anion_nonend_bulk(start_z_index:end_z_index) - n_cation_end_bulk(start_z_index:end_z_index) - &
            n_cation_nonend_bulk(start_z_index:end_z_index) - n_anion_end_bulk(start_z_index:end_z_index)

       n_mbar_bulk(start_z_index:end_z_index) = bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads

       !Now calculate the functional derivative
       calculate_hard_sphere_functional_deriv(start_z_index:end_z_index) = &
            (n_cation_nonend_bulk(start_z_index:end_z_index) * &
            GetYMix(n_mbar_bulk(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_cation)/(r_cation - num_end_monomers_cation) * &
            (((4.0_dp * pi * (hs_diameter**3)/3.0_dp) * GetAExDerivIntegrand(n_mbar_bulk(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2)) - &
            ((4.0_dp * pi * (hs_diameter**3)/3.0_dp) * GetAExDerivIntegrand(n_mbar_bulk(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 1))))/beta + &
                                !
            (n_cation_end_bulk(start_z_index:end_z_index)/num_end_monomers_cation * &
            ((4.0_dp * pi * (hs_diameter**3)/3.0_dp) * GetAExDerivIntegrand(n_mbar_bulk(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2)))/beta + &
                                !
                                !Now we've done the cation, we can do the anion.
                                !
            (n_anion_nonend_bulk(start_z_index:end_z_index) * &
            GetYMix(n_mbar_bulk(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_anion)/(r_anion - num_end_monomers_anion) * &
            (((4.0_dp * pi * (hs_diameter**3)/3.0_dp) * GetAExDerivIntegrand(n_mbar_bulk(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2)) - &
            ((4.0_dp * pi * (hs_diameter**3)/3.0_dp) * GetAExDerivIntegrand(n_mbar_bulk(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 1))))/beta + &
                                !
            (n_anion_end_bulk(start_z_index:end_z_index)/num_end_monomers_anion * &
            ((4.0_dp * pi * (hs_diameter**3)/3.0_dp) * GetAExDerivIntegrand(n_mbar_bulk(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2)))/beta


    else

       integrand1(:) = calculate_n_sbar(GetAExDerivIntegrand(n_mbar(:), n_sbar(:), 1))
       integrand2(:) = calculate_n_sbar(GetAExDerivIntegrand(n_mbar(:), n_sbar(:), 2))

       !print *, "integrand1", integrand1
       !print *, "integrand2", integrand2
       !print *, "3", calculate_central_difference(n_cation_nonend(start_z_index:end_z_index), hs_diameter/n_discretised_points_z)
       !call abort()

       ratio_cation_end(start_z_index:end_z_index) = calculate_central_difference(n_cation_end(start_z_index:end_z_index), hs_diameter/n_discretised_points_z)
       ratio_cation_nonend(start_z_index:end_z_index) = calculate_central_difference(n_cation_nonend(start_z_index:end_z_index), hs_diameter/n_discretised_points_z)
       ratio_anion_end(start_z_index:end_z_index) = calculate_central_difference(n_anion_end(start_z_index:end_z_index), hs_diameter/n_discretised_points_z)
       ratio_anion_nonend(start_z_index:end_z_index) = calculate_central_difference(n_anion_nonend(start_z_index:end_z_index), hs_diameter/n_discretised_points_z)
       deriv_wrt(start_z_index:end_z_index) = calculate_central_difference(n_deriv_wrt(start_z_index:end_z_index), hs_diameter/n_discretised_points_z)


       !if(all(ratio_cation_end(start_z_index:end_z_index) == 0.0_dp)) ratio_cation_end(start_z_index:end_z_index) = 1.0_dp
       !if(all(ratio_cation_nonend(start_z_index:end_z_index) == 0.0_dp)) ratio_cation_nonend(start_z_index:end_z_index) = 1.0_dp
       !if(all(ratio_anion_end(start_z_index:end_z_index) == 0.0_dp)) ratio_anion_end(start_z_index:end_z_index) = 1.0_dp
       !if(all(ratio_anion_nonend(start_z_index:end_z_index) == 0.0_dp)) ratio_anion_nonend(start_z_index:end_z_index) = 1.0_dp

       !This should only happen on the first iteration
       if(all(deriv_wrt(start_z_index:end_z_index) == 0.0_dp)) deriv_wrt(start_z_index:end_z_index) = 1.0_dp

       print *, "STARTING..."
       print *, "n_deriv_wrt = ", n_deriv_wrt(start_z_index:end_z_index)
       print *, ""
       print *, "deriv_wrt = ", deriv_wrt
       print *, ""
       print *, "n_cation_end = ", n_cation_end
       print *, ""
       print *, "n_cation_nonend = ", n_cation_nonend
       print *, ""
       print *, "ratio_cation_end = ", ratio_cation_end
       print *, ""
       print *, "ratio_cation_nonend = ", ratio_cation_nonend
       print *, ""

       do ij = start_z_index, end_z_index
          if(deriv_wrt(ij) == 0.0_dp) then
             if(ij < size(deriv_wrt)/2) then
                deriv_wrt(ij) = deriv_wrt(ij+1)
             else
                deriv_wrt(ij) = deriv_wrt(ij-1)
             end if
             !print *, "HELLLLLLLLLLLLLLLLLLLOOOOOOOOOOOO", ij

          end if
          ! if(ratio_cation_end(ij) == 0.0_dp) then
          !    ratio_cation_end(ij) = ratio_cation_end(ij-1)
          ! end if
          ! if(ratio_cation_nonend(ij) == 0.0_dp) then
          !    ratio_cation_nonend(ij) = ratio_cation_nonend(ij-1)
          ! end if
          ! if(ratio_anion_end(ij) == 0.0_dp) then
          !    ratio_anion_end(ij) = ratio_anion_end(ij-1)
          ! end if
          ! if(ratio_anion_nonend(ij) == 0.0_dp) then
          !    ratio_anion_nonend(ij) = ratio_anion_nonend(ij-1)
          ! end if
       end do

       !print *, "deriv_wrt2 = ", deriv_wrt

       !if(any(deriv_wrt(start_z_index:end_z_index) == 0.0_dp)) then

       !print *, "***************************************"
       !print *, "***************************************"
       !print *, "**********WARNING - DENOMINATOR ZERO***"
       !print *, "***************************************"
       !print *, "***************************************"
       !print *, "deriv_wrt(start_z_index:end_z_index) = ", deriv_wrt(start_z_index:end_z_index)
       !print *, ""
       !print *, ""
       !print *, ""
       !   ratio_cation_end(start_z_index:end_z_index) = 1.0_dp
       !   ratio_cation_nonend(start_z_index:end_z_index) = 1.0_dp
       !   ratio_anion_end(start_z_index:end_z_index) = 1.0_dp
       !   ratio_anion_nonend(start_z_index:end_z_index) = 1.0_dp
       !else
       ratio_cation_end(start_z_index:end_z_index) = ratio_cation_end(start_z_index:end_z_index) / deriv_wrt(start_z_index:end_z_index)
       ratio_cation_nonend(start_z_index:end_z_index) = ratio_cation_nonend(start_z_index:end_z_index) / deriv_wrt(start_z_index:end_z_index)
       ratio_anion_end(start_z_index:end_z_index) = ratio_anion_end(start_z_index:end_z_index) / deriv_wrt(start_z_index:end_z_index)
       ratio_anion_nonend(start_z_index:end_z_index) = ratio_anion_nonend(start_z_index:end_z_index) / deriv_wrt(start_z_index:end_z_index)
       !end if


       calculate_hard_sphere_functional_deriv(start_z_index:end_z_index) = &
            (ratio_cation_nonend(start_z_index:end_z_index) * &
            GetYMix(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_cation)/(r_cation - num_end_monomers_cation) * &
            (GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2) - GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 1)))/beta + &
                                !
            (ratio_cation_end(start_z_index:end_z_index) * &
            (1.0_dp / num_end_monomers_cation) * &
            (GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2)))/beta + &
                                !
            (n_cation_nonend(start_z_index:end_z_index) * &
            GetYMix(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_cation)/(r_cation - num_end_monomers_cation) * &
            (((4.0_dp * pi * (hs_diameter**3)/3.0_dp) * integrand2(start_z_index:end_z_index)) - &
            ((4.0_dp * pi * (hs_diameter**3)/3.0_dp) * integrand1(start_z_index:end_z_index))))/beta + &
                                !
            (n_cation_end(start_z_index:end_z_index)/num_end_monomers_cation * &
            ((4.0_dp * pi * (hs_diameter**3)/3.0_dp) * integrand2(start_z_index:end_z_index)))/beta + &
                                !
                                !Now we've done the cation, we can do the anion.
                                !
            (ratio_anion_nonend(start_z_index:end_z_index) * &
            GetYMix(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_anion)/(r_anion - num_end_monomers_anion) * &
            (GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2) - GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 1)))/beta + &
                                !
            (ratio_anion_end(start_z_index:end_z_index) * &
            (1.0_dp / num_end_monomers_anion) * &
            (GetAEx(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 2)))/beta + &
                                !
            (n_anion_nonend(start_z_index:end_z_index) * &
            GetYMix(n_mbar(start_z_index:end_z_index), n_sbar(start_z_index:end_z_index), 'm', r_anion)/(r_anion - num_end_monomers_anion) * &
            (((4.0_dp * pi * (hs_diameter**3)/3.0_dp) * integrand2(start_z_index:end_z_index)) - &
            ((4.0_dp * pi * (hs_diameter**3)/3.0_dp) * integrand1(start_z_index:end_z_index))))/beta + &
                                !
            (n_anion_end(start_z_index:end_z_index)/num_end_monomers_anion * &
            ((4.0_dp * pi * (hs_diameter**3)/3.0_dp) * integrand2(start_z_index:end_z_index)))/beta


    end if

    !print *, "calculate_hard_sphere_functional_deriv = ", calculate_hard_sphere_functional_deriv
    !calculate_hard_sphere_functional_deriv = 0.0_dp

  end function calculate_hard_sphere_functional_deriv


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

       calculate_surface_dispersion_functional_deriv(iz) = 2.0_dp * pi * epsilon_LJ_particle_wall * (1.0_dp) *(&
            ( (2.0_dp/45.0_dp)* (hs_d_divide_z**9.0_dp) ) - &
            ( (1.0_dp/3.0_dp) * (hs_d_divide_z**3.0_dp) ) + &
            ( (2.0_dp/45.0_dp)* (hs_d_divide_h_minus_z**9.0_dp) ) - &
            ( (1.0_dp/3.0_dp) * (hs_d_divide_h_minus_z**3.0_dp) ))
    end do

    call setNonCalculatedRegionToZero(calculate_surface_dispersion_functional_deriv)

    !print *, "surface_dispersion = ", calculate_surface_dispersion_functional_deriv
    !call abort()
    !calculate_surface_dispersion_functional_deriv = 0.0_dp

  end function calculate_surface_dispersion_functional_deriv

  function calculate_vanderWaals_functional_deriv(n_s)
    real(dp), dimension(:), intent(in) :: n_s
    real(dp), dimension(size(n_s)) :: calculate_vanderWaals_functional_deriv

    calculate_vanderWaals_functional_deriv = -4.0_dp * epsilon_LJ_particle_particle * (hs_diameter**6.0_dp) * (1.0_dp) * &
         2.0_dp * pi * integrate_z_cylindrical(n_s, van_der_waals_density_indept_integrand, "all_z")

    call setNonCalculatedRegionToZero(calculate_vanderWaals_functional_deriv)

    !print *, "calculate_vanderWaals_functional_deriv = ", calculate_vanderWaals_functional_deriv
    !call abort()
    !calculate_vanderWaals_functional_deriv = 0.0_dp

  end function calculate_vanderWaals_functional_deriv

  function calculate_surface_electrostatic_functional_deriv(size_array, charge)
    integer, intent(in) :: size_array
    real(dp), intent(in) :: charge

    real(dp), dimension(size_array) :: calculate_surface_electrostatic_functional_deriv

    integer :: ij, iz
    integer :: start_z_index
    integer :: end_z_index

    real(dp) :: d_to_left_wall
    real(dp) :: d_to_right_wall

    call get_allowed_z_values(start_z_index, end_z_index, size_array)

    do ij = start_z_index, end_z_index
       iz = ij !- start_z_index + 1
       d_to_left_wall = real((iz-1) * hs_diameter, dp) / real(n_discretised_points_z, dp)
       d_to_right_wall = real((size_array-iz) * hs_diameter, dp) / real(n_discretised_points_z, dp)

       !print *, "charge = ", charge
       !print *, "surface_charge_density_left_wall = ", surface_charge_density_left_wall
       !print *, "surface_charge_density_rightt_wall = ", surface_charge_density_right_wall

       !call abort()
       !print *,  "d_to_left_wall, d_to_right_wall= ", d_to_left_wall, d_to_right_wall
       !call abort()
       calculate_surface_electrostatic_functional_deriv(ij) = (-1.0_dp / ( 2.0_dp * epsilon0 * epsilonr)) * (&
            charge * surface_charge_density_left_wall * d_to_left_wall + &
            charge * surface_charge_density_right_wall * d_to_right_wall)
    end do

    call setNonCalculatedRegionToZero(calculate_surface_electrostatic_functional_deriv)
    !print *, surface_charge_density_left_wall
    !print *, surface_charge_density_right_wall
    !print *, charge
    !print *, epsilon0, epsilonr
    !print *, "calculate_surface_electrostatic_functional_deriv = ", calculate_surface_electrostatic_functional_deriv
    !call abort()
    !calculate_surface_electrostatic_functional_deriv = 0.0_dp

  end function calculate_surface_electrostatic_functional_deriv

  function calculate_electrostatic_like_term_functional_deriv(n, charge, calculate_bulk)
    real(dp), dimension(:), intent(in) :: n
    real(dp), intent(in) :: charge
    logical, intent(in) :: calculate_bulk
    real(dp), dimension(size(n)) :: calculate_electrostatic_like_term_functional_deriv

    real(dp) :: lambda
    real(dp), dimension(size(n)) :: unit_array

    ! integer :: start_z_index
    ! integer :: end_z_index
    ! call get_allowed_z_values(start_z_index, end_z_index, size(n))

    if(charge > 0.0_dp) then
       CURRENT_BULK_BEAD_DENSITY = bulk_density_positive_beads
    else if(charge < 0.0_dp) then
       CURRENT_BULK_BEAD_DENSITY = bulk_density_negative_beads
    else if(charge == 0.0_dp) then
       !Note this is a float comparison so should never happen.
       !In any case if charge = 0.0_dp then we get an answer of zero anyway.
       !But we want to set it to something, to protect against zero * {some unset variable},
       !which I don't know what may happen.
       CURRENT_BULK_BEAD_DENSITY = bulk_density_neutral_beads
    end if

    ! if(calculate_bulk) then
    !    calculate_electrostatic_like_term_functional_deriv(:) = (-2.0_dp / (4.0_dp * pi * epsilon0 * epsilonr)) * pi * (charge**2) * &
    !         n(size(n)/2) * ((((end_z_index - start_z_index)*hs_diameter/n_discretised_points_z)/get_lambda()) - (hs_diameter**2)) / beta
    ! else
    !    calculate_electrostatic_like_term_functional_deriv(:) = (-2.0_dp / (4.0_dp * pi * epsilon0 * epsilonr)) * pi * (charge**2) * &
    !         integrate_z_cylindrical(n, electrostatic_like_integrand, "all_z") / beta
    ! end if

    !lambda = get_lambda()
    unit_array(:) = 1.0_dp

    if(calculate_bulk) then
       !calculate_electrostatic_like_term_functional_deriv(:) = 0.0_dp

       !calculate_electrostatic_like_term_functional_deriv(:) = (-1.0_dp / (2.0_dp * epsilon0 * epsilonr)) * (charge**2) * &
       !     CURRENT_BULK_BEAD_DENSITY * integrate_z_cylindrical(unit_array, electrostatic_like_integrand, "all_z")

       !print *, calculate_electrostatic_like_term_functional_deriv
       !call abort()

       calculate_electrostatic_like_term_functional_deriv(:) = 0.0_dp
       
    else
       calculate_electrostatic_like_term_functional_deriv(:) = (-1.0_dp / (2.0_dp * epsilon0 * epsilonr)) * (charge**2) * &
            integrate_z_cylindrical(n, electrostatic_like_integrand, "all_z")

       !print *, "n = ", n
       !print *, "integrate_z_cylindrical(n, electrostatic_unlike_integrand, 'all_z') = ", integrate_z_cylindrical(n, electrostatic_unlike_integrand, "all_z")
       !call abort()

       !print *, "integrate_z_cylindrical = ", integrate_z_cylindrical(n, electrostatic_unlike_integrand, "all_z")
       !print *, "calculate_electrostatic_like_term_functional_deriv = ", calculate_electrostatic_like_term_functional_deriv
       !call abort()


    end if

    call setNonCalculatedRegionToZero(calculate_electrostatic_like_term_functional_deriv)

    !print *, "calculate_electrostatic_like_term_functional_deriv = ", n(200), calculate_electrostatic_like_term_functional_deriv(200)
    !calculate_electrostatic_like_term_functional_deriv = 0.0_dp

  end function calculate_electrostatic_like_term_functional_deriv

  function calculate_electrostatic_unlike_term_functional_deriv(n, charge1, charge2, calculate_bulk)
    real(dp), dimension(:), intent(in) :: n
    real(dp), intent(in) :: charge1
    real(dp), intent(in) :: charge2
    logical, intent(in) :: calculate_bulk
    real(dp), dimension(size(n)) :: calculate_electrostatic_unlike_term_functional_deriv

    real(dp) :: d, h

    d = chi_parameter*hs_diameter
    !h = (size(n) - hs_diameter/n_discretised_points_z) * (hs_diameter/n_discretised_points_z)

    if(charge1 > 0.0_dp) then
       CURRENT_BULK_BEAD_DENSITY = bulk_density_positive_beads
    else if(charge1 < 0.0_dp) then
       CURRENT_BULK_BEAD_DENSITY = bulk_density_negative_beads
    else if(charge1 == 0.0_dp) then
       !Note this is a float comparison so should never happen.
       !In any case if charge = 0.0_dp then we get an answer of zero anyway.
       !But we want to set it to something, to protect against zero * {some unset variable},
       !which I don't know what may happen.
       CURRENT_BULK_BEAD_DENSITY = bulk_density_neutral_beads
    end if




    if(calculate_bulk) then
       !calculate_electrostatic_unlike_term_functional_deriv(:) = 0.0_dp

       calculate_electrostatic_unlike_term_functional_deriv(:) = (-1.0_dp / (2.0_dp * epsilon0 * epsilonr)) * charge1 * charge2 * &
            CURRENT_BULK_BEAD_DENSITY * ( d**2 )

    else

       !print *, "n = ", charge1, charge2, n
       calculate_electrostatic_unlike_term_functional_deriv(:) = (-1.0_dp / (2.0_dp * epsilon0 * epsilonr)) * charge1 * charge2 * &       
            integrate_z_cylindrical(n, electrostatic_unlike_integrand, "all_z")


       !print *, "integrate_z_cylindrical = ", integrate_z_cylindrical(n, electrostatic_unlike_integrand, "all_z")
       !call abort()
    end if


    call setNonCalculatedRegionToZero(calculate_electrostatic_unlike_term_functional_deriv)
    !calculate_electrostatic_unlike_term_functional_deriv = 0.0_dp

  end function calculate_electrostatic_unlike_term_functional_deriv

  function electrostatic_unlike_integrand(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp) :: electrostatic_unlike_integrand

    integer :: xi_int
    real(dp) :: xi_real

    xi_int = xi - z ! centre xi on z
    xi_real = real(xi_int,dp) * hs_diameter / real(n_discretised_points_z,dp)

    if(abs(xi_int) >= (chi_parameter * n_discretised_points_z)) then
       electrostatic_unlike_integrand = abs(xi_real)
    else
       electrostatic_unlike_integrand = chi_parameter * hs_diameter
    end if

  end function electrostatic_unlike_integrand


  function electrostatic_like_integrand(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp) :: electrostatic_like_integrand

    real(dp) :: lambda
    real(dp) :: abs_xi

    lambda = get_lambda()
    abs_xi = (abs(z - xi)*hs_diameter/real(n_discretised_points_z,dp))

    electrostatic_like_integrand = (abs_xi) + ((exp(-(abs_xi*lambda)))/lambda)

  end function electrostatic_like_integrand

  function get_lambda()
    real(dp) :: get_lambda
    real(dp) :: s

    s = (3.0_dp/(4.0_dp * pi * CURRENT_BULK_BEAD_DENSITY))**(1.0_dp/3.0_dp)
    get_lambda = sqrt(2.0_dp)/s

  end function get_lambda

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
       ! van_der_waals_density_indept_integrand = 1.0_dp / &
       !      (4.0_dp * ( (real(hs_diameter,dp) * cos(asin(real(xi_real,dp)/real(hs_diameter,dp))))**2.0_dp &
       !      + real(xi_real,dp)**2.0_dp )**2.0_dp)
       van_der_waals_density_indept_integrand = 1.0_dp / (4.0_dp * (hs_diameter**4.0_dp))
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
