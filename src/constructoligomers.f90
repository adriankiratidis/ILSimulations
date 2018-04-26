!Contains routines to contstruct the various types of oligomers
!that we may wish to calculate.  These routines are present in the module
!as they need to be called in two different places in the code.
!Once during iteration in order to calculate the updated density,
!and then again post density calculation in order to calculate the
!ideal chain contribution to the free energy.
module constructoligomers
  use kinds
  use parameters
  use integratephispherical
  use normalisation
  implicit none
  private

  public :: UpdateDensities

  !Private subroutines
  !UpdateSinglePositiveSphereDensity
  !UpdateSingleNeutralSphereDensity
  !UpdateSinglePositiveSphereDensity
  !UpdateSinglePositiveNeutralMinusSphereDensities
  !UpdateC4MIMPositiveBeadDensities
  !UpdateC4MINNeutralBeadDensities
  !UpdateC4MINNegativeBeadDensities
  
contains

  subroutine UpdateDensities(lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated)
    real(dp), dimension(:), intent(in) :: lambda1
    real(dp), dimension(:), intent(out) :: n1_updated

    real(dp), dimension(:), intent(in), optional :: lambda2
    real(dp), dimension(:), intent(out), optional :: n2_updated

    real(dp), dimension(:), intent(in), optional :: lambda3
    real(dp), dimension(:), intent(out), optional :: n3_updated


    if(trim(ionic_liquid_name) == "SingleNeutralSpheres") then

       if(present(lambda2) .or. present(n2_updated) .or. present(lambda3) .or. present(n3_updated)) then
          print *, "When doing single spheres should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       end if

       call UpdateSingleNeutralSphereDensity(lambda1, n1_updated)

    else if(trim(ionic_liquid_name) == "SinglePositiveSpheres") then

       if(present(lambda2) .or. present(n2_updated) .or. present(lambda3) .or. present(n3_updated)) then
          print *, "When doing single spheres should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       end if

       call UpdateSinglePositiveSphereDensity(lambda1, n1_updated)

    else if(trim(ionic_liquid_name) == "SingleNegativeSpheres") then

       if(present(lambda2) .or. present(n2_updated) .or. present(lambda3) .or. present(n3_updated)) then
          print *, "When doing single spheres should not have more than one lambda/density pair present"
          print *, "...aborting..."
          call abort()
       end if

       call UpdateSingleNegativeSphereDensity(lambda1, n1_updated)

    else if(trim(ionic_liquid_name) == "SinglePositiveNeutralMinusSpheres") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateSinglePositiveNeutralMinusSphereDensities(lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated)
       end if

    else if(trim(ionic_liquid_name) == "C4MIM_BF4") then

       if( (.not. present(lambda2)) .or. (.not. present(n2_updated)) .or. &
            (.not. present(lambda3)) .or. (.not. present(n3_updated)) ) then
          print *, "constructoligomers.90: UpdateDensities:"
          print *, "Trying to Update positive, neutral and negative spheres"
          print *, "but haven't been passed the appropriate number of arguments"
          print *, "coding error...aborting..."
          call abort()
       else
          call UpdateC4MIMPositiveBeadDensities(lambda1, lambda2, n1_updated)
          call UpdateC4MIMNeutralBeadDensities(lambda1, lambda2, n2_updated)
          call UpdateBF4NegativeBeadDensities(lambda3, n3_updated)
       end if

    end if

  end subroutine UpdateDensities

  subroutine UpdateSinglePositiveSphereDensity(lambda_plus, n_plus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated

    if(size(lambda_plus) /= size(n_plus_updated)) then
       print *, "constructoligomers.f90: ConstructSingleSphereArrangment: "
       print *, "Size mismatch.  size(lambda_plus) /= size(n_plus_updated)."
       print *, "can't update...aborting..."
       call abort()
    else
       n_plus_updated = lambda_plus
    end if

    call RenormaliseToBulkDensity(n_plus_updated, "n+")

  end subroutine UpdateSinglePositiveSphereDensity

  subroutine UpdateSingleNeutralSphereDensity(lambda_neutral, n_neutral_updated)
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    if(size(lambda_neutral) /= size(n_neutral_updated)) then
       print *, "constructoligomers.f90: ConstructSingleSphereArrangment: "
       print *, "Size mismatch.  size(lambda_neutral) /= size(n_plus_neutral)."
       print *, "can't update...aborting..."
       call abort()
    else
       n_neutral_updated = lambda_neutral
    end if

    call RenormaliseToBulkDensity(n_neutral_updated, "n0")

  end subroutine UpdateSingleNeutralSphereDensity

  subroutine UpdateSingleNegativeSphereDensity(lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    if(size(lambda_minus) /= size(n_minus_updated)) then
       print *, "constructoligomers.f90: ConstructSingleSphereArrangment: "
       print *, "Size mismatch.  size(lambda_plus) /= size(n_plus_updated)."
       print *, "can't update...aborting..."
       call abort()
    else
       n_minus_updated = lambda_minus
    end if

    call RenormaliseToBulkDensity(n_minus_updated, "n-")

  end subroutine UpdateSingleNegativeSphereDensity

  subroutine UpdateSinglePositiveNeutralMinusSphereDensities(lambda_plus, n_plus_updated, &
       lambda_neutral, n_neutral_updated, lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(out) :: n_plus_updated

    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    !First ensure all the lambda sizes are the same
    if((size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))) then
       print *, "constructoligomers.f90: UpdateSinglePositiveNeutralMinusSphereDensities:"
       print *, "(size(lambda_plus) /= size(lambda_neutral)) .or. (size(lambda_plus) /= size(lambda_minus))"
       print *, "Input variable size mismatch...aborting..."
       call abort()
    end if

    call UpdateSinglePositiveSphereDensity(lambda_plus, n_plus_updated)
    call UpdateSingleNeutralSphereDensity(lambda_neutral, n_neutral_updated)
    call UpdateSingleNegativeSphereDensity(lambda_minus, n_minus_updated)

  end subroutine UpdateSinglePositiveNeutralMinusSphereDensities

  subroutine UpdateC4MIMPositiveBeadDensities(lambda_plus, lambda_neutral, n_plus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_plus_updated

    real(dp), dimension(:), allocatable :: c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp
    integer :: array_size

    ! First check the input variables are the same size
    if( (size(lambda_plus) == size(lambda_neutral)) .and. (size(lambda_plus) == size(n_plus_updated))) then

       array_size = size(lambda_plus)
       allocate(c8c1(array_size), c9c10(array_size), c7(array_size), c6(array_size), &
            c5(array_size), c4(array_size), c2(array_size), c3p(array_size), c3pp(array_size), c3ppp(array_size))
    else
       print *, "constructoligomers.f90: UpdateC4MIMPositiveBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_plus_updated)"
       print *, "...aborting..."
       call abort()
    end if

    ! Note that there is a factor of 1.0_dp/(4.0_dp * pi * hs_diameter**2)
    ! from the delta function bond, a factor of 2 * pi from the theta integral
    ! and a factor of hs_diameter**2 from the jacobian.  Combining these factors
    ! gives the required 0.5_dp factor that we are multiplying by.
    c9c10 = 0.5_dp * integrate_phi_spherical(lambda_plus)

    c8c1 = 0.5_dp * integrate_phi_spherical(lambda_neutral)

    c7 = 0.5_dp * integrate_phi_spherical(lambda_neutral * c8c1)

    c6 = 0.5_dp * integrate_phi_spherical(lambda_neutral * c7)

    c5 = 0.5_dp * integrate_phi_spherical(lambda_neutral * c6)

    c4 = 0.5_dp * integrate_phi_spherical(lambda_plus * c5)

    c3p = 0.5_dp * integrate_phi_spherical(lambda_plus * c4 * c9c10 * c9c10)

    c2 = 0.5_dp * integrate_phi_spherical(lambda_plus * c8c1)

    c3pp = 0.5_dp * integrate_phi_spherical(lambda_plus * c4 * c2 * c9c10)

    c3ppp = 0.5_dp * integrate_phi_spherical(lambda_plus * c2 * c9c10 * c9c10)

    !Calculate the resulting positive bead densities.
    !n_plus_updated = nc2 + nc3 + nc4 + nc9 + nc10 =  nc2 + nc3 + nc4 + 2*nc9
    n_plus_updated = bulk_density * ( (lambda_plus * c8c1 * c3p) + (lambda_plus * c2 * c9c10 * c9c10 * c4) + &
         (lambda_plus * c3ppp * c5) + (2.0_dp * (lambda_plus * c3pp)) )

    call RenormaliseToBulkDensity(n_plus_updated, "n+")

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp)

  end subroutine UpdateC4MIMPositiveBeadDensities

  subroutine UpdateC4MIMNeutralBeadDensities(lambda_plus, lambda_neutral, n_neutral_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    real(dp), dimension(:), allocatable :: c3p, c3ppp, c9c10, c8c1, c7, c6, c5, c4, c2 !contributions also used in +ve beads.
    real(dp), dimension(:), allocatable :: c2p, c4p, c5p, c6p, c7p !extra contributions for -ve beads.

    integer :: array_size

    ! First check the input variables are the same size
    if( (size(lambda_plus) == size(lambda_neutral)) .and. (size(lambda_plus) == size(n_neutral_updated))) then

       array_size = size(lambda_plus)
       allocate(c3p(array_size), c3ppp(array_size), c8c1(array_size), c9c10(array_size), c2p(array_size), &
            c7(array_size), c6(array_size), c5(array_size), c4(array_size), c2(array_size), &
            c4p(array_size), c5p(array_size), c6p(array_size), c7p(array_size))
    else
       print *, "constructoligomers.f90: UpdateC4MINNeutralBeads:"
       print *, "Size mismatch. size(lambda_plus) /= size(lambda_neutral)"
       print *, "or size(lambda_plus) == size(n_neutral_updated)"
       print *, "...aborting..."
       call abort()
    end if

    c9c10 = 0.5_dp * integrate_phi_spherical(lambda_plus)

    c8c1 = 0.5_dp * integrate_phi_spherical(lambda_neutral)

    c7 = 0.5_dp * integrate_phi_spherical(lambda_neutral * c8c1)

    c6 = 0.5_dp * integrate_phi_spherical(lambda_neutral * c7)

    c5 = 0.5_dp * integrate_phi_spherical(lambda_neutral * c6)

    c4 = 0.5_dp * integrate_phi_spherical(lambda_plus * c5)

    c2 = 0.5_dp * integrate_phi_spherical(lambda_plus * c8c1)

    c3p = 0.5_dp * integrate_phi_spherical(lambda_plus * c4 * c9c10 * c9c10)

    c3ppp = 0.5_dp * integrate_phi_spherical(lambda_plus * c2 * c9c10 * c9c10)

    c2p = 0.5_dp * integrate_phi_spherical(lambda_plus * c3p)

    c4p = 0.5_dp * integrate_phi_spherical(lambda_plus * c3ppp)

    c5p = 0.5_dp * integrate_phi_spherical(lambda_neutral * c4p)

    c6p = 0.5_dp * integrate_phi_spherical(lambda_neutral * c5p)

    c7p = 0.5_dp * integrate_phi_spherical(lambda_neutral * c6p)

    !Calculate the resulting neutral bead densities.
    !n_zero_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * ( (lambda_neutral * c2p) + (lambda_neutral * c4p * c6) + (lambda_neutral * c5p * c7) &
         + (lambda_neutral * c6p * c8c1) + (lambda_neutral * c7p) )

    call RenormaliseToBulkDensity(n_neutral_updated, "n0")

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(c3p, c3ppp, c8c1)
    deallocate(c2p, c4p, c5p, c6p, c7p)

  end subroutine UpdateC4MIMNeutralBeadDensities

  subroutine UpdateBF4NegativeBeadDensities(lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated

    real(dp), dimension(:), allocatable :: a1a2a3a4, a5p
    integer :: array_size

    ! First check the input variables are the same size
    if(size(lambda_minus) == size(n_minus_updated)) then
       array_size = size(lambda_minus)
       allocate(a1a2a3a4(array_size))
       allocate(a5p(array_size))
    else
       print *, "constructoligomers.f90: UpdateBF4NegativeBeads:"
       print *, "Size mismatch. size(lambda_minus) /= size(n_minus_updated)"
       print *, "...aborting..."
       call abort()
    end if

    !Calculate the required contributions for the anion
    a1a2a3a4 = 0.5_dp * integrate_phi_spherical(lambda_minus)

    a5p = 0.5_dp * integrate_phi_spherical(lambda_minus * (a1a2a3a4 ** 3.0_dp))

    !Calculate the resulting negative bead densities.
    !n_minus_updated = na1 + na2 + na3 + na4 + na5 = 4*na1 + na5
    n_minus_updated = bulk_density * ( 4.0_dp*(lambda_minus * a5p) + (lambda_minus * (a1a2a3a4**4.0_dp)) )

    call RenormaliseToBulkDensity(n_minus_updated, "n-")

    !Don't check if the variables have been allocated.
    !We want an error thrown if they haven't been (which should never happen anyway).
    deallocate(a1a2a3a4, a5p)

  end subroutine UpdateBF4NegativeBeadDensities

end module constructoligomers
