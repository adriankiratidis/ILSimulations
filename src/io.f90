!Contains all io routines
module io
  use kinds
  use parameters
  implicit none
  private

  public :: WriteDensityOutputFormatted
  
contains

  !Writes array out to txt file
  subroutine WriteDensityOutputFormatted(output_array, file_stub, file_suffix)
    real(dp), dimension(:), intent(in) :: output_array
    character(len=*), intent(in)       :: file_stub
    character(len=*), intent(in)       :: file_suffix

    integer :: iz
    integer :: file_unit
    file_unit = 171

    open(file_unit, file=trim(file_stub)//"-"//trim(file_suffix)//".txt", action='write')
    do iz = 1, size(output_array)
       write(file_unit, *) real((iz-1),dp)/real(n_discretised_points_z,dp), output_array(iz)*(hs_diameter**3)
    end do
    close(file_unit)


  end subroutine WriteDensityOutputFormatted

end module io
