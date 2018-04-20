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

contains

  subroutine ReNormaliseToBulkDensity(n_plus, n_neutral, n_minus)
    real(dp), dimension(:) :: n_plus
    real(dp), dimension(:) :: n_neutral
    real(dp), dimension(:) :: n_minus

    integer :: start_z_index
    integer :: end_z_index

    real(dp), dimension(size(n_plus)) :: bulk_plus
    real(dp), dimension(size(n_neutral)) :: bulk_neutral
    real(dp), dimension(size(n_minus)) :: bulk_minus

    !Calculate the value in the bulk
    bulk_plus = get_bulk_density(n_plus)
    bulk_neutral = get_bulk_density(n_neutral)
    bulk_minus = get_bulk_density(n_minus)

    call get_allowed_z_values(start_z_index, end_z_index, size(n_plus))

    n_plus(start_z_index:end_z_index) = bulk_density * n_plus(start_z_index:end_z_index) &
         / bulk_plus(start_z_index:end_z_index)

    n_neutral(start_z_index:end_z_index) = bulk_density * n_neutral(start_z_index:end_z_index) &
         / bulk_neutral(start_z_index:end_z_index)

    n_minus(start_z_index:end_z_index) = bulk_density * n_minus(start_z_index:end_z_index) &
         / bulk_minus(start_z_index:end_z_index)

  end subroutine ReNormaliseToBulkDensity

end module normalisation
