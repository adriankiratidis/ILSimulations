!This is a module that calculates the lambdas of our model.
!Recall that lambda is defined to be the paritial derivative of
!all non-ideal contributions to the grand potential functional
!w.r.t the bead densities.
module lambdas
  use kinds
  use helpers
  use parameters
  use functionalderivatives
  implicit none
  private

  public :: CalculateLambdasDifference
  public :: CalculateLambdasBulk
  !public :: CalculateLambdasDifference_old

contains
  
  subroutine CalculateLambdas(lambda_plus, n_plus_total, n_plus_end, lambda_neutral, n_neutral_total, n_neutral_end, lambda_minus, n_minus_total, n_minus_end, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus
    real(dp), dimension(:), intent(in)  :: n_plus_total
    real(dp), dimension(:), intent(in)  :: n_plus_end

    real(dp), dimension(:), intent(out) :: lambda_neutral
    real(dp), dimension(:), intent(in)  :: n_neutral_total
    real(dp), dimension(:), intent(in)  :: n_neutral_end

    real(dp), dimension(:), intent(out) :: lambda_minus
    real(dp), dimension(:), intent(in)  :: n_minus_total
    real(dp), dimension(:), intent(in)  :: n_minus_end
    integer, intent(in)                 :: ith_plate_separation

    real(dp), dimension(size(lambda_plus)) :: lambda_common_terms

    integer :: input_array_size

    ! First ensure that the sizes of all the input variables are the same.
    input_array_size = size(lambda_plus)
    if((size(lambda_neutral) /= input_array_size) .or. &
         (size(lambda_minus) /= input_array_size) .or. &
         (size(n_plus_total) /= input_array_size) .or. &
         (size(n_neutral_total) /= input_array_size) .or. &
         (size(n_minus_total) /= input_array_size)) then

       print *, "lambda.f90: CalculateLambdas:"
       print *, "Input arrays must all be the same size"
       print *, "Almost certainly coding error as opposed to input error...aborting..."
       call abort()
    end if

    lambda_plus(:) = 0.0_dp
    lambda_neutral(:) = 0.0_dp
    lambda_minus(:) = 0.0_dp

    ! Now calculate the terms in common to all the lambdas
    lambda_common_terms = CalculateLambdaCommonTerms(n_plus_total, n_plus_end, n_neutral_total, n_neutral_end, n_minus_total, n_minus_end, ith_plate_separation, .false.)

    ! Now calculate our lambdas(r)
    lambda_plus = beta * (lambda_common_terms + &
         CalculateLambdaPlusSpecificTerms(n_plus_total, n_plus_end, n_minus_total, n_plus_total + n_neutral_total + n_minus_total, .false.))
    lambda_neutral = beta * (lambda_common_terms + &
         CalculateLambaNeutralSpecificTerms(n_neutral_total, n_neutral_end, n_plus_total + n_neutral_total + n_minus_total, .false.))
    lambda_minus = beta * (lambda_common_terms + &
         CalculateLambdaMinusSpecificTerms(n_minus_total, n_minus_end, n_plus_total, n_plus_total + n_neutral_total + n_minus_total, .false.))

  end subroutine CalculateLambdas

  subroutine CalculateLambdasBulk(lambda_plus_bulk, n_plus_total, lambda_neutral_bulk, n_neutral_total, lambda_minus_bulk, &
       n_minus_total, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus_bulk
    real(dp), dimension(:), intent(in)  :: n_plus_total

    real(dp), dimension(:), intent(out) :: lambda_neutral_bulk
    real(dp), dimension(:), intent(in)  :: n_neutral_total

    real(dp), dimension(:), intent(out) :: lambda_minus_bulk
    real(dp), dimension(:), intent(in)  :: n_minus_total
    integer, intent(in)                 :: ith_plate_separation

    real(dp), dimension(size(lambda_plus_bulk)) :: lambda_common_terms_bulk
    real(dp), dimension(size(n_neutral_total)) :: n_neutral_total_array

    integer :: input_array_size

    real(dp), dimension(size(n_plus_total)) :: n_plus_end
    real(dp), dimension(size(n_neutral_total)) :: n_neutral_end
    real(dp), dimension(size(n_minus_total)) :: n_minus_end

    ! First ensure that the sizes of all the input variables are the same.
    input_array_size = size(lambda_plus_bulk)
    if((size(lambda_neutral_bulk) /= input_array_size) .or. &
         (size(lambda_minus_bulk) /= input_array_size) .or. &
         (size(n_plus_total) /= input_array_size) .or. &
         (size(n_neutral_total) /= input_array_size) .or. &
         (size(n_minus_total) /= input_array_size)) then

       print *, "lambda.f90: CalculateLambdasBulk:"
       print *, "Input arrays must all be the same size"
       print *, "Almost certainly coding error as opposed to input error...aborting..."
       call abort()
    end if

    lambda_plus_bulk(:) = 0.0_dp
    lambda_neutral_bulk(:) = 0.0_dp
    lambda_minus_bulk(:) = 0.0_dp

    lambda_common_terms_bulk = CalculateLambdaCommonTerms(n_plus_total, n_plus_end, n_neutral_total, n_neutral_end, &
         n_minus_total, n_minus_end, ith_plate_separation, .true.)

    call SetBulkDensityOfEndMonomers(n_plus_total, n_plus_end, n_neutral_total, n_neutral_end, n_minus_total, n_minus_end)
    
    lambda_plus_bulk = beta * (lambda_common_terms_bulk + &
         CalculateLambdaPlusSpecificTerms(n_plus_total, n_plus_end, n_minus_total, n_plus_total + n_neutral_total + n_minus_total, .true.))
    lambda_neutral_bulk = beta * (lambda_common_terms_bulk + &
         CalculateLambaNeutralSpecificTerms(n_neutral_total, n_neutral_end, n_plus_total + n_neutral_total + n_minus_total, .true.))
    lambda_minus_bulk = beta * (lambda_common_terms_bulk + &
         CalculateLambdaMinusSpecificTerms(n_minus_total, n_minus_end, n_plus_total, n_plus_total + n_neutral_total + n_minus_total, .true.))

  end subroutine CalculateLambdasBulk

  subroutine CalculateLambdasDifference(lambda_plus, n_plus_total, n_plus_end, lambda_neutral, n_neutral_total, n_neutral_end, lambda_minus, n_minus_total, n_minus_end, ith_plate_separation)
    real(dp), dimension(:), intent(out) :: lambda_plus
    real(dp), dimension(:), intent(in)  :: n_plus_total
    real(dp), dimension(:), intent(in)  :: n_plus_end

    real(dp), dimension(:), intent(out) :: lambda_neutral
    real(dp), dimension(:), intent(in)  :: n_neutral_total
    real(dp), dimension(:), intent(in)  :: n_neutral_end

    real(dp), dimension(:), intent(out) :: lambda_minus
    real(dp), dimension(:), intent(in)  :: n_minus_total
    real(dp), dimension(:), intent(in)  :: n_minus_end

    integer, intent(in)                 :: ith_plate_separation

    real(dp), dimension(size(lambda_plus)) :: lambda_plus_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_neutral_bulk
    real(dp), dimension(size(lambda_plus)) :: lambda_minus_bulk

    integer :: input_array_size

    ! First ensure that the sizes of all the input variables are the same.
    input_array_size = size(lambda_plus)
    if((size(lambda_neutral) /= input_array_size) .or. &
         (size(lambda_minus) /= input_array_size) .or. &
         (size(n_plus_total) /= input_array_size) .or. &
         (size(n_neutral_total) /= input_array_size) .or. &
         (size(n_minus_total) /= input_array_size)) then

       print *, "lambda.f90: CalculateLambdasBulk:"
       print *, "Input arrays must all be the same size"
       print *, "Almost certainly coding error as opposed to input error...aborting..."
       call abort()
    end if

    call CalculateLambdasBulk(lambda_plus_bulk, n_plus_total, lambda_neutral_bulk, n_neutral_total, lambda_minus_bulk, &
         n_minus_total, ith_plate_separation)
    call CalculateLambdas(lambda_plus, n_plus_total, n_plus_end, lambda_neutral, n_neutral_total, n_neutral_end, lambda_minus, n_minus_total, &
         n_minus_end, ith_plate_separation)

    lambda_plus(:) = lambda_plus_bulk(:) - lambda_plus(:)
    lambda_neutral(:) = lambda_neutral_bulk(:) - lambda_neutral(:)
    lambda_minus(:) = lambda_minus_bulk(:) - lambda_minus(:)

  end subroutine CalculateLambdasDifference


  ! subroutine CalculateLambdasDifference_old(lambda_plus, n_plus_total, lambda_neutral, n_neutral_total, lambda_minus, n_minus_total, ith_plate_separation)
  !   real(dp), dimension(:), intent(out) :: lambda_plus
  !   real(dp), dimension(:), intent(in)  :: n_plus_total

  !   real(dp), dimension(:), intent(out) :: lambda_neutral
  !   real(dp), dimension(:), intent(in)  :: n_neutral_total

  !   real(dp), dimension(:), intent(out) :: lambda_minus
  !   real(dp), dimension(:), intent(in)  :: n_minus_total
  !   integer, intent(in)                 :: ith_plate_separation


  !   real(dp), dimension(size(lambda_plus)) :: lambda_common_terms
  !   real(dp), dimension(size(lambda_plus)) :: lambda_common_terms_bulk

  !   real(dp), dimension(size(lambda_plus)) :: lambda_plus_bulk
  !   real(dp), dimension(size(lambda_plus)) :: lambda_neutral_bulk
  !   real(dp), dimension(size(lambda_plus)) :: lambda_minus_bulk

  !   integer :: input_array_size

  !   ! First ensure that the sizes of all the input variables are the same.
  !   input_array_size = size(lambda_plus)
  !   if((size(lambda_neutral) /= input_array_size) .or. &
  !        (size(lambda_minus) /= input_array_size) .or. &
  !        (size(n_plus_total) /= input_array_size) .or. &
  !        (size(n_neutral_total) /= input_array_size) .or. &
  !        (size(n_minus_total) /= input_array_size)) then

  !      print *, "lambda.f90: CalculateLambdas:"
  !      print *, "Input arrays must all be the same size"
  !      print *, "Almost certainly coding error as opposed to input error...aborting..."
  !      call abort()
  !   end if

  !   lambda_plus(:) = 0.0_dp
  !   lambda_neutral(:) = 0.0_dp
  !   lambda_minus(:) = 0.0_dp

  !   lambda_plus_bulk(:) = 0.0_dp
  !   lambda_neutral_bulk(:) = 0.0_dp
  !   lambda_minus_bulk(:) = 0.0_dp

  !   ! Now calculate the terms in common to all the lambdas
  !   lambda_common_terms = CalculateLambdaCommonTerms(n_plus_total, n_neutral_total, n_minus_total, ith_plate_separation, .false.)

  !   lambda_common_terms_bulk = CalculateLambdaCommonTerms(n_plus_total, n_neutral_total, n_minus_total, ith_plate_separation, .true.)

  !   ! Now calculate our lambdas(r)
  !   lambda_plus = beta * (lambda_common_terms + CalculateLambdaPlusSpecificTerms(n_plus_total, n_minus_total, .false.))
  !   lambda_neutral = beta * (lambda_common_terms + CalculateLambaNeutralSpecificTerms(n_plus_total, n_neutral_total, n_minus_total, .false.))
  !   lambda_minus = beta * (lambda_common_terms + CalculateLambdaMinusSpecificTerms(n_plus_total, n_minus_total, .false.))

  !   lambda_plus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaPlusSpecificTerms(n_plus_total, n_minus_total, .true.))
  !   lambda_neutral_bulk = beta * (lambda_common_terms_bulk + CalculateLambaNeutralSpecificTerms(n_plus_total, n_neutral_total, n_minus_total, .true.))
  !   lambda_minus_bulk = beta * (lambda_common_terms_bulk + CalculateLambdaMinusSpecificTerms(n_plus_total, n_minus_total, .true.))

  !   !Now calculate the bulk value, lambda^{bulk} and in order to return e^{lambda^{bulk} - lambda(r)}
  !   lambda_plus = lambda_plus_bulk - lambda_plus
  !   lambda_neutral = lambda_neutral_bulk - lambda_neutral
  !   lambda_minus = lambda_minus_bulk - lambda_minus

  ! end subroutine CalculateLambdasDifference_old

  
  function CalculateLambdaCommonTerms(n_plus_total, n_plus_end, n_neutral_total, n_neutral_end, &
       n_minus_total, n_minus_end, ith_plate_separation, calculate_bulk, iteration)
    real(dp), dimension(:), intent(in)  :: n_plus_total
    real(dp), dimension(:), intent(in)  :: n_plus_end
    real(dp), dimension(:), intent(in)  :: n_neutral_total
    real(dp), dimension(:), intent(in)  :: n_neutral_end
    real(dp), dimension(:), intent(in)  :: n_minus_total
    real(dp), dimension(:), intent(in)  :: n_minus_end
    integer, intent(in)                 :: ith_plate_separation
    logical, intent(in)                 :: calculate_bulk
    integer, optional :: iteration

    real(dp), dimension(size(n_plus_total)) :: CalculateLambdaCommonTerms

    real(dp), dimension(size(n_plus_total)) :: hs_term
    real(dp), dimension(size(n_plus_total)) :: van_der_waals_term
    real(dp), dimension(size(n_plus_total)) :: surface_fluid_dispersion_term

    real(dp), dimension(size(n_plus_total)) :: n_s
    real(dp), dimension(size(n_plus_total)) :: n_s_end

    hs_term = 0.0_dp
    van_der_waals_term = 0.0_dp
    surface_fluid_dispersion_term = 0.0_dp

    if(calculate_bulk) then
       n_s(:) = bulk_density_positive_beads + bulk_density_neutral_beads + bulk_density_negative_beads

       !hs_term = calculate_hardsphere_functional_deriv(n_s, .true.)
       !van_der_waals_term = calculate_vanderWaals_functional_deriv(n_s)

    else
       n_s = n_plus_total + n_neutral_total + n_minus_total
       !n_s_end = n_plus_end + n_neutral_end + n_minus_end
       !hs_term = calculate_hardsphere_functional_deriv(n_s, n_s_end, .false.)

       surface_fluid_dispersion_term = calculate_surface_dispersion_functional_deriv(&
            ith_plate_separation, size(surface_fluid_dispersion_term))

       !van_der_waals_term = calculate_vanderWaals_functional_deriv(n_s)
    end if

    !print *, "hs_term = ", hs_term
    !print *, "surface_fluid_dispersion_term = ", surface_fluid_dispersion_term
    !print *, "van_der_waals_term = ", van_der_waals_term
    !call abort()

    CalculateLambdaCommonTerms = hs_term + surface_fluid_dispersion_term + van_der_waals_term

    !print *, "common terms = ", CalculateLambdaCommonTerms

  end function CalculateLambdaCommonTerms


  function CalculateLambdaPlusSpecificTerms(n_plus_total, n_plus_end, n_minus_total, n_s, calculate_bulk)
    real(dp), dimension(:), intent(in)  :: n_plus_total
    real(dp), dimension(:), intent(in) :: n_plus_end
    real(dp), dimension(:), intent(in)  :: n_minus_total
    real(dp), dimension(:), intent(in) :: n_s
    logical, intent(in) :: calculate_bulk

    real(dp), dimension(size(n_plus_total)) :: CalculateLambdaPlusSpecificTerms

    real(dp), dimension(size(n_plus_total)) :: surface_electrostatic_term
    real(dp), dimension(size(n_plus_total)) :: like_electrostatic_term
    real(dp), dimension(size(n_plus_total)) :: unlike_electrostatic_term

    real(dp), dimension(size(n_plus_total)) :: hard_sphere_term
    
    real(dp), dimension(size(n_plus_total))  :: n_plus_total_in
    real(dp), dimension(size(n_minus_total))  :: n_minus_total_in

    surface_electrostatic_term = 0.0_dp
    like_electrostatic_term = 0.0_dp
    unlike_electrostatic_term = 0.0_dp
    hard_sphere_term = 0.0_dp
    
    if(calculate_bulk) then

       n_plus_total_in(:) = bulk_density_positive_beads
       n_minus_total_in(:) = bulk_density_negative_beads
       
       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_plus_total_in, positive_bead_charge)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_minus_total_in, positive_bead_charge, negative_bead_charge)

       hard_sphere_term = calculate_hardsphere_functional_deriv(n_plus_total, n_plus_end, n_s, calculate_bulk, r_plus)
       
    else
       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_plus_total, positive_bead_charge)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_minus_total, positive_bead_charge, negative_bead_charge)

       surface_electrostatic_term = calculate_surface_electrostatic_functional_deriv(size(n_plus_total), positive_bead_charge)

       hard_sphere_term = calculate_hardsphere_functional_deriv(n_plus_total, n_plus_end, n_s, calculate_bulk, r_plus)
       
    end if

    CalculateLambdaPlusSpecificTerms = hard_sphere_term !surface_electrostatic_term + like_electrostatic_term + unlike_electrostatic_term + hard_sphere_term

  end function CalculateLambdaPlusSpecificTerms


  function CalculateLambaNeutralSpecificTerms(n_neutral_total, n_neutral_end, n_s, calculate_bulk)
    real(dp), dimension(:), intent(in)  :: n_neutral_total
    real(dp), dimension(:), intent(in)  :: n_neutral_end
    real(dp), dimension(:), intent(in)  :: n_s
    logical, intent(in) :: calculate_bulk

    real(dp), dimension(size(n_neutral_total)) :: CalculateLambaNeutralSpecificTerms
    real(dp), dimension(size(n_neutral_total)) :: hard_sphere_term

    hard_sphere_term(:) = 0.0_dp
    
    !print *, "calculate_bulk = ", calculate_bulk
    !print *, "n_neutral_total = ", n_neutral_total
    !print *, "n_neutral_end = ", n_neutral_end
    
    hard_sphere_term = calculate_hardsphere_functional_deriv(n_neutral_total, n_neutral_end, n_s, calculate_bulk, r_neutral)

    !print *, hard_sphere_term
    
    !call abort()
    
    CalculateLambaNeutralSpecificTerms = hard_sphere_term !surface_electrostatic_term + like_electrostatic_term + unlike_electrostatic_term

  end function CalculateLambaNeutralSpecificTerms


  function CalculateLambdaMinusSpecificTerms(n_minus_total, n_minus_end, n_plus_total, n_s, calculate_bulk)
    real(dp), dimension(:), intent(in)  :: n_minus_total
    real(dp), dimension(:), intent(in)  :: n_minus_end
    real(dp), dimension(:), intent(in)  :: n_plus_total
    real(dp), dimension(:), intent(in)  :: n_s
    logical, intent(in) :: calculate_bulk

    real(dp), dimension(size(n_plus_total)) :: CalculateLambdaMinusSpecificTerms

    real(dp), dimension(size(n_plus_total)) :: surface_electrostatic_term
    real(dp), dimension(size(n_plus_total)) :: like_electrostatic_term
    real(dp), dimension(size(n_plus_total)) :: unlike_electrostatic_term
    
    real(dp), dimension(size(n_plus_total)) :: hard_sphere_term
    
    real(dp), dimension(size(n_plus_total))  :: n_plus_total_in
    real(dp), dimension(size(n_minus_total))  :: n_minus_total_in

    surface_electrostatic_term = 0.0_dp
    like_electrostatic_term = 0.0_dp
    unlike_electrostatic_term = 0.0_dp
    hard_sphere_term = 0.0_dp
    
    if(calculate_bulk) then

       n_plus_total_in(:) = bulk_density_positive_beads
       n_minus_total_in(:) = bulk_density_negative_beads

       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_minus_total_in, positive_bead_charge)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_plus_total_in, positive_bead_charge, negative_bead_charge)

       hard_sphere_term = calculate_hardsphere_functional_deriv(n_minus_total, n_minus_end, n_s, calculate_bulk, r_plus)
       
    else
       like_electrostatic_term = calculate_electrostatic_like_term_functional_deriv(n_minus_total, positive_bead_charge)
       unlike_electrostatic_term = calculate_electrostatic_unlike_term_functional_deriv(n_plus_total, positive_bead_charge, negative_bead_charge)

       surface_electrostatic_term = calculate_surface_electrostatic_functional_deriv(size(n_minus_total), negative_bead_charge)

       hard_sphere_term = calculate_hardsphere_functional_deriv(n_minus_total, n_minus_end, n_s, calculate_bulk, r_plus)
       
    end if

    CalculateLambdaMinusSpecificTerms = hard_sphere_term !surface_electrostatic_term + like_electrostatic_term + unlike_electrostatic_term

  end function CalculateLambdaMinusSpecificTerms

end module lambdas
