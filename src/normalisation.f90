!Contains the routines relevant to the normalisation of arrays
!after each iteration.
module normalisation
  use kinds
  use parameters
  use lambdas
  use helpers
  implicit none
  private

  public :: ReNormaliseToBulkDensity
  public :: SetSingleSphereBeadDensityFromBulkIonDensity
  public :: SetC4MIN_BF4BeadDensityFromBulkIonDensity

  real(dp) :: bulk_density_positive_beads
  real(dp) :: bulk_density_neutral_beads
  real(dp) :: bulk_density_negative_beads
  
contains

  !The second two arguments are optional soas to allow the call of the routine with only
  !one parameter, if for example we are studying the single sphere  case for example.
  subroutine ReNormaliseToBulkDensity(n_plus, n_neutral, n_minus)
    real(dp), dimension(:) :: n_plus
    real(dp), dimension(:), optional :: n_neutral
    real(dp), dimension(:), optional :: n_minus

    integer :: start_z_index
    integer :: end_z_index

    real(dp), dimension(size(n_plus)) :: bulk_plus
    real(dp), dimension(size(n_plus)) :: bulk_neutral
    real(dp), dimension(size(n_plus)) :: bulk_minus

    call get_allowed_z_values(start_z_index, end_z_index, size(n_plus))

    bulk_plus = get_bulk_density(n_plus)
    n_plus(start_z_index:end_z_index) = bulk_density_positive_beads * n_plus(start_z_index:end_z_index) &
         / bulk_plus(start_z_index:end_z_index)

    if(present(n_neutral)) then
       if(size(n_neutral) /= size(n_plus)) then
          print *, "normalisation.f90: ReNormaliseToBulkDensity:"
          print *, "Size mismatch. size(n_neutral) /= size(n_plus)"
          print *, "...aborting..."
          call abort()
       else
          bulk_neutral = get_bulk_density(n_neutral)
          n_neutral(start_z_index:end_z_index) = bulk_density_neutral_beads * n_neutral(start_z_index:end_z_index) &
               / bulk_neutral(start_z_index:end_z_index)
       end if
    end if

    if(present(n_minus)) then
       if(size(n_minus) /= size(n_plus)) then
          print *, "normalisation.f90: ReNormaliseToBulkDensity:"
          print *, "Size mismatch. size(n_minus) /= size(n_plus)"
          print *, "...aborting..."
          call abort()
       else
          bulk_minus = get_bulk_density(n_minus)
          n_minus(start_z_index:end_z_index) = bulk_density_minus_beads * n_minus(start_z_index:end_z_index) &
               / bulk_minus(start_z_index:end_z_index)
       end if
    end if

    call setNonCalculatedRegionToZero(n_plus, n_minus, n_neutral)

  end subroutine ReNormaliseToBulkDensity

  subroutine SetSingleSphereBeadDensityFromBulkIonDensity()
    
    bulk_density_positive_beads = bulk_density
    bulk_density_neutral_beads = bulk_density
    bulk_density_negative_beads = bulk_density
  end subroutine SetSingleSphereBeadDensityFromBulkIonDensity

  subroutine SetC4MIN_BF4BeadDensityFromBulkIonDensity()
    
    bulk_density_positive_beads = 5.0_dp * bulk_density
    bulk_density_neutral_beads = 5.0_dp * bulk_density
    bulk_density_negative_beads = 5.0_dp * bulk_density
  end subroutine SetC4MIN_BF4BeadDensityFromBulkIonDensity

end module normalisation
