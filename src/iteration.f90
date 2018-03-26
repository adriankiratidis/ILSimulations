!Contains relevant routines for our iterative procedure.  I.e. Allocation of
!Variables, intialisation of ansatz and convergence check.
module iteration
  use kinds
  use parameters
  implicit none
  private

  public :: converged
  public :: InitialiseIntegrationAnsatz
  
contains

  !Routine that returns true iff the average squared difference for all beads is
  !less than the iterative_tolerance input by the user.
  function converged(n_plus_updated, n_neutral_updated, n_minus_updated, n_plus, n_neutral, n_minus)
    real(dp), dimension(:), intent(in) :: n_plus_updated
    real(dp), dimension(:), intent(in) :: n_neutral_updated
    real(dp), dimension(:), intent(in) :: n_minus_updated
    real(dp), dimension(:), intent(in) :: n_plus
    real(dp), dimension(:), intent(in) :: n_neutral
    real(dp), dimension(:), intent(in) :: n_minus
    logical                            :: converged

    real(dp) :: av_sq_diff_plus
    real(dp) :: av_sq_diff_minus
    real(dp) :: av_sq_diff_neutral
    
    converged = .false.

    !Ensure that we recieve arrays of the correct size
    if(size(n_plus_updated) /= size(n_plus) .or. &
         size(n_neutral_updated) /= size(n_neutral) .or. &
         size(n_minus_updated) /= size(n_minus)) then
       print *, "iteration.f90: converged needs to compare arrays of same size for convergence."
       print *, "array size mismatch...aborting..."
       call abort
    end if

    !Find the average difference of least squares for the 3 bead types
    av_sq_diff_plus = sum((n_plus_updated - n_plus)**2)/size(n_plus)
    av_sq_diff_neutral = sum((n_neutral_updated - n_neutral)**2)/size(n_plus)
    av_sq_diff_minus = sum((n_minus_updated - n_minus)**2)/size(n_plus)

    !Check for convergence
    if((av_sq_diff_plus <= iterative_tolerance) .and. &
       (av_sq_diff_minus <= iterative_tolerance) .and. &
       (av_sq_diff_neutral <= iterative_tolerance))then

       print *, "Iterative scheme successfully converged."
       converged = .true.
    end if

  end function converged

  !Routine to initialise our ansatz for our integrative scheme
  subroutine InitialiseIntegrationAnsatz(n_plus, n_neutral, n_minus)
    real(dp), dimension(:) :: n_plus
    real(dp), dimension(:) :: n_neutral
    real(dp), dimension(:) :: n_minus

    n_plus(:) = 1.0_dp
    n_neutral(:) = 1.0_dp
    n_minus(:) = 1.0_dp
    
  end subroutine InitialiseIntegrationAnsatz

end module iteration
