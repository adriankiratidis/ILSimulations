!Contains the routines relevant to the normalisation of arrays
!after each iteration.
module normalisation
  use kinds
  use lambdas
  use helpers
  use parameters
  implicit none
  private

  public :: ReNormaliseToBulkDensity
  
contains

  subroutine ReNormaliseToBulkDensity(n, bead_charge)
    real(dp), dimension(:) :: n !density
    character(len=*), intent(in) :: bead_charge

    integer :: start_z_index
    integer :: end_z_index

    real(dp), dimension(size(n)) :: n_bulk
    real(dp) :: bead_bulk_density

    call get_allowed_z_values(start_z_index, end_z_index, size(n))

    n_bulk = get_bulk_density(n)
    if(trim(bead_charge) == "n+") then
       
       bead_bulk_density = bulk_density_positive_beads
       
    else if(trim(bead_charge) == "n0") then
       
       bead_bulk_density = bulk_density_neutral_beads
       
    else if(trim(bead_charge) == "n-") then
       
       bead_bulk_density = bulk_density_negative_beads
       
    else
       print *, "normalisation.f90: ReNormaliseToBulkDensity: "
       print *, "Unsupported bead_charge type of ", trim(bead_charge)
       print *, "...aborting..."
       call abort()
    end if

    n(start_z_index:end_z_index) = bead_bulk_density * n(start_z_index:end_z_index) &
         / n_bulk(start_z_index:end_z_index)

    call setNonCalculatedRegionToZero(n)

  end subroutine ReNormaliseToBulkDensity

end module normalisation
