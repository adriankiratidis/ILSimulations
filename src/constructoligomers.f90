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

  public :: ConstructSingleSphereArrangment
  public :: 

contains

  subroutine ConstructSingleSphereArrangment(lambda1, n1_updated, lambda2, n2_updated, lambda3, n3_updated)
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


end module constructoligomers
