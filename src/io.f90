!Contains all io routines
module io
  use kinds
  implicit none
  private

  public :: WriteDensityOutput
  
contains

  !Writes array out to txt file
  subroutine WriteDensityOutputFormatted(ouput_array, file_stub, file_suffix)
    real(dp), dimension(:), intent(in) :: output_array
    character(len=*), intent(in)       :: file_stub
    character(len=*), intent(in)       :: file_suffix
    
  end subroutine WriteDensityOutputFormatted

end module io
