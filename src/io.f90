!Contains all io routines
module io
  use kinds
  use parameters
  implicit none
  private

  public :: WriteOutputFormattedAsFunctionOfPosition
  public :: WriteOutputFormattedAsFunctionOfPlateSeparation
  
contains

  !Writes array out to txt file
  subroutine WriteOutputFormattedAsFunctionOfPosition(output_array, file_stub, file_suffix)
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

  end subroutine WriteOutputFormattedAsFunctionOfPosition

  subroutine WriteOutputFormattedAsFunctionOfPlateSeparation(output, file_stub, file_suffix)
    real(dp), dimension(:), intent(in) :: output
    character(len=*), intent(in)       :: file_stub
    character(len=*), intent(in)       :: file_suffix

    integer :: id
    integer :: file_unit
    file_unit = 173

    !First check the size is correct
    if(size(output) /= size(plate_separations)) then
       print *, "io.f90: WriteGrandPotentialOutputFormatted:"
       print *, "size(output) /= size(plate_separations)"
       print *, "should only have one potential value per plate separation"
       print *, "size mismatch...aborting..."
       call abort()
    else
       open(file_unit, file=trim(file_stub)//"-"//trim(file_suffix)//".txt", action='write')
       do id = 1, size(output)
          write(file_unit, *) plate_separations(id), output
       end do
       close(file_unit)
    end if

  end subroutine WriteOutputFormattedAsFunctionOfPlateSeparation

end module io
