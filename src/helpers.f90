! Module containing so-called helper functions.
! Used to perform neccesary logic on data structures.
module helpers
  use parameters
  implicit none
  private

  public :: get_allowed_z_values

contains

  subroutine get_allowed_z_values(start_z_index, end_z_index, total_points_z)
    integer, intent(out) :: start_z_index
    integer, intent(out) :: end_z_index
    integer, intent(in)  :: total_points_z
    
    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    !Note we use integer division here
    if(is_even(n_discretised_points_z)) then
       start_z_index = (n_discretised_points_z/2)  + 1
       end_z_index = total_points_z - (n_discretised_points_z/2)
    else
       start_z_index = (n_discretised_points_z/2)  + 2
       end_z_index = total_points_z - (n_discretised_points_z/2) - 1
    end if

  end subroutine get_allowed_z_values

  pure function is_even(an_int)
    integer, intent(in) :: an_int
    logical             :: is_even

    if(modulo(an_int, 2) == 0) then
       is_even = .true.
    else
       is_even = .false.
    end if

    return
  end function is_even

end module helpers
