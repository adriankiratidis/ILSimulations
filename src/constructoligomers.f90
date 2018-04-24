!Contains routines to contstruct the various types of oligomers
!that we may wish to calculate.  These routines are present in the module
!as they need to be called in two different places in the code.
!Once during iteration in order to calculate the updated density,
!and then again post density calculation in order to calculate the
!ideal chain contribution to the free energy.
module constructoligomers
  use kinds
  implicit none
  private

  public :: UpdateSingleSphereArrangment
  public :: UpdateC4MIMPositiveBeads
  public :: UpdateC4MINNeutralBeads
  public :: UpdateBF4NegativeBeads

contains

  subroutine UpdateSingleSphereArrangment(lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated)
    real(dp), dimension(:), intent(in) :: lambda1
    real(dp), dimension(:), intent(out) :: n1_updated

    real(dp), dimension(:), intent(in), optional :: lambda2
    real(dp), dimension(:), intent(out), optional :: n2_updated

    real(dp), dimension(:), intent(in), optional :: lambda3
    real(dp), dimension(:), intent(out), optional :: n3_updated

    !First do the n1 calculation.
    if(size(lambda1) /= size(n1_updated)) then
       print *, "constructoligomers.f90: ConstructSingleSphereArrangment: "
       print *, "Size mismatch.  size(lambda1) /= size(n1_updated)."
       print *, "can't update...aborting..."
       call abort()
    else
       n1_updated = lambda1
    end if

    !Now do the n2 calculation
    if(present(lambda2)) then
       if(.not. present(n2_updated)) then
          print *, "constructoligomers.f90: ConstructSingleSphereArrangment: "
          print *, "Must have n2_updated present if lambda2 is present"
          print *, "coding error...aborting..."
          call abort()
       else
          if(size(lambda2) /= size(n2_updated)) then
             print *, "constructoligomers.f90: ConstructSingleSphereArrangment: "
             print *, "Size mismatch.  size(lambda1) /= size(n1_updated)."
             print *, "can't update...aborting..."
             call abort()
          else
             n2_updated = lambda2
          end if
       end if
    end if

    !Now do the n3 calculation
     if(present(lambda3)) then
       if(.not. present(n3_updated)) then
          print *, "constructoligomers.f90: ConstructSingleSphereArrangment: "
          print *, "Must have n3_updated present if lambda3 is present"
          print *, "coding error...aborting..."
          call abort()
       else
          if(size(lambda3) /= size(n3_updated)) then
             print *, "constructoligomers.f90: ConstructSingleSphereArrangment: "
             print *, "Size mismatch.  size(lambda3) /= size(n3_updated)."
             print *, "can't update...aborting..."
             call abort()
          else
             n3_updated = lambda3
          end if
       end if
    end if

  end subroutine ConstructSingleSphereArrangment
  
  subroutine UpdateC4MIMPositiveBeads(lambda_plus, lambda_neutral, n_plus_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_plus_updated

    real(dp), dimension(:), allocatable :: c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp

    call InitialiseVariableDiscretisation(c8c1, c9c10, c7, c6, c5, c4, c2, c3p, c3pp, c3ppp)

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

  end subroutine UpdateC4MIMPositiveBeads

  subroutine UpdateC4MINNeutralBeads(lambda_plus, lambda_neutral, n_neutral_updated)
    real(dp), dimension(:), intent(in) :: lambda_plus
    real(dp), dimension(:), intent(in) :: lambda_neutral
    real(dp), dimension(:), intent(out) :: n_neutral_updated

    real(dp), dimension(:), allocatable :: c2p, c4p, c5p, c6p, c7p

    call InitialiseVariableDiscretisation(c2p, c4p, c5p, c6p, c7p)
    

    c2p = 0.5_dp * integrate_phi_spherical(lambda_plus * c3p)

    c4p = 0.5_dp * integrate_phi_spherical(lambda_plus * c3ppp)

    c5p = 0.5_dp * integrate_phi_spherical(lambda_neutral * c4p)

    c6p = 0.5_dp * integrate_phi_spherical(lambda_neutral * c5p)

    c7p = 0.5_dp * integrate_phi_spherical(lambda_neutral * c6p)

    !Calculate the resulting neutral bead densities.
    !n_zero_updated = nc1 + nc5 + nc6 + nc7 + nc8
    n_neutral_updated = bulk_density * ( (lambda_neutral * c2p) + (lambda_neutral * c4p * c6) + (lambda_neutral * c5p * c7) &
         + (lambda_neutral * c6p * c8c1) + (lambda_neutral * c7p) )

  end subroutine UpdateC4MINNeutralBeads

  subroutine UpdateBF4NegativeBeads(lambda_minus, n_minus_updated)
    real(dp), dimension(:), intent(in) :: lambda_minus
    real(dp), dimension(:), intent(out) :: n_minus_updated
    
    real(dp), dimension(:), allocatable :: a1a2a3a4, a5p


    call InitialiseVariableDiscretisation(a1a2a3a4, a5p)
    
    !Calculate the required contributions for the anion
    a1a2a3a4 = 0.5_dp * integrate_phi_spherical(lambda_minus)

    a5p = 0.5_dp * integrate_phi_spherical(lambda_minus * (a1a2a3a4 ** 3.0_dp))

    !Calculate the resulting negative bead densities.
    !n_minus_updated = na1 + na2 + na3 + na4 + na5 = 4*na1 + na5
    n_minus_updated = bulk_density * ( 4.0_dp*(lambda_minus * a5p) + (lambda_minus * (a1a2a3a4**4.0_dp)) )

  end subroutine UpdateBF4NegativeBeads



end module constructoligomers
