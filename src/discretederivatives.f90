!Contains routines that perform the discrete derivative
module discretederivatives
  use kinds
  use parameters
  implicit none
  private

  public :: calculate_central_difference
  public :: calculate_backward_difference
  public :: calculate_forward_difference


  interface calculate_central_difference
     module procedure calculate_central_difference_plate_separation
     module procedure calculate_central_difference_against_constant
  end interface
  
contains

  function calculate_central_difference_against_constant(input_array, distance_between_points)
    real(dp), dimension(:), intent(in) :: input_array
    real(dp), intent(in)               :: distance_between_points
    real(dp), dimension(size(input_array)) :: calculate_central_difference_against_constant

    integer :: ij
    
    !If the input_array size is one we can't calculate any derivatives
    !but we may still want the density profile. So print a warning and set the derivative to 0.
    if(size(input_array) == 1) then
       print *, "discretederivatives.f90: calculate_central_difference: "
       print *, "trying to calculate a discrete derivative with only one point"
       print *, "not possible.  Setting the derivative to 0 and returning as we"
       print *, "may want the density profile."
       calculate_central_difference_against_constant = 0.0_dp
       return
    end if

    !Calculate the central difference whenever we can
    do ij = 2, size(input_array) - 1
       calculate_central_difference_against_constant(ij) = ( input_array(ij+1) - input_array(ij-1) ) / &
            (2.0_dp * distance_between_points)
    end do

    !If we are at the first array point, we can't do a central difference,
    !therefore, we calculate the forward difference.
    calculate_central_difference_against_constant(1) = ( input_array(2) - input_array(1) ) / &
         (distance_between_points)

    !If we are at the last array point, we can't do a central difference,
    !therefore, we calculate the backward difference.
    calculate_central_difference_against_constant(size(input_array)) = &
         ( input_array(size(input_array)) - input_array(size(input_array) - 1) ) / &
         (distance_between_points)

    ! do ij =1, size(calculate_central_difference_against_constant)
    !    if(calculate_central_difference_against_constant(ij) <= 1E-15) then
    !       calculate_central_difference_against_constant(ij) = 1E-15
    !    end if
    ! end do
    
  end function calculate_central_difference_against_constant


  !Calculate the central difference, unless we are at the beginning
  !or end of the walls in which case we take the appropriate forward
  !or backward difference.
  function calculate_central_difference_plate_separation(input_array)
    real(dp), dimension(:), intent(in) :: input_array
    real(dp), dimension(size(input_array)) :: calculate_central_difference_plate_separation

    integer :: ij

    calculate_central_difference_plate_separation = 0.0_dp

    !Check input_array is the correct size
    if(size(input_array) /= size(plate_separations)) then
       print *, "discretederivatives.f90: calculate_central_difference: "
       print *, "size(input_array) /= size(plate_separations)"
       !print *, "size mismatch....likely coding error...aborting..."
       !call abort()
    end if

    !If the input_array size is one we can't calculate any derivatives
    !but we may still want the density profile. So print a warning and set the derivative to 0.
    if(size(input_array) == 1) then
       print *, "discretederivatives.f90: calculate_central_difference: "
       print *, "trying to calculate a discrete derivative with only one point"
       print *, "not possible.  Setting the derivative to 0 and returning as we"
       print *, "may want the density profile."
       calculate_central_difference_plate_separation = 0.0_dp
       return
    end if

    !Calculate the central difference whenever we can
    do ij = 2, size(input_array) - 1
       calculate_central_difference_plate_separation(ij) = ( input_array(ij+1) - input_array(ij-1) ) / &
            ( abs(plate_separations(ij+1) - plate_separations(ij-1)) * hs_diameter)
    end do

    !If we are at the first array point, we can't do a central difference,
    !therefore, we calculate the forward difference.
    calculate_central_difference_plate_separation(1) = ( input_array(2) - input_array(1) ) / &
         ( abs(plate_separations(2) - plate_separations(1)) * hs_diameter)

    !If we are at the last array point, we can't do a central difference,
    !therefore, we calculate the backward difference.
    calculate_central_difference_plate_separation(size(input_array)) = &
         ( input_array(size(input_array)) - input_array(size(input_array) - 1) ) / &
         ( abs(plate_separations(size(input_array)) - plate_separations(size(input_array) - 1)) * hs_diameter)

  end function calculate_central_difference_plate_separation

  !Calculate the backward difference unless we are at the start point in
  !which case we calculate the forward difference.
  function calculate_backward_difference(input_array)
    real(dp), dimension(:), intent(in) :: input_array
    real(dp), dimension(size(input_array)) :: calculate_backward_difference

    integer :: ij

    !Check input_array is the correct size
    if(size(input_array) /= size(plate_separations)) then
       print *, "discretederivatives.f90: calculate_backward_difference: "
       print *, "size(input_array) /= size(plate_separations)"
       !print *, "size mismatch....likely coding error...aborting..."
       !call abort()
    end if

    !If the input_array size is one we can't calculate any derivatives
    !but we may still want the density profile. So print a warning and set the derivative to 0.
    if(size(input_array) == 1) then
       print *, "discretederivatives.f90: calculate_central_difference: "
       print *, "trying to calculate a discrete derivative with only one point"
       print *, "not possible.  Setting the derivative to 0 and returning as we"
       print *, "may want the density profile."
       calculate_backward_difference = 0.0_dp
       return
    end if

    calculate_backward_difference = 0.0_dp

    !Calculat the backward difference whenever we can.
    do ij = 2, size(input_array)
       calculate_backward_difference(ij) = ( input_array(ij) - input_array(ij-1) ) / &
            ( abs(plate_separations(ij) - plate_separations(ij-1)) * hs_diameter)
    end do

    !On the first point, we can't calculate the backward difference so calculate the forward difference.
    calculate_backward_difference(1) = ( input_array(2) - input_array(1) ) / &
         ( abs(plate_separations(2) - plate_separations(1)) * hs_diameter)

  end function calculate_backward_difference

  !Calculate the forward difference unless we are at the endpoint
  !in which case we calculate the backward difference.
  function calculate_forward_difference(input_array)
    real(dp), dimension(:), intent(in) :: input_array
    real(dp), dimension(size(input_array)) :: calculate_forward_difference

    integer :: ij

    !Check input_array is the correct size
    if(size(input_array) /= size(plate_separations)) then
       print *, "discretederivatives.f90: calculate_forward_difference: "
       print *, "size(input_array) /= size(plate_separations)"
       !print *, "size mismatch....likely coding error...aborting..."
       !call abort()
    end if

    !If the input_array size is one we can't calculate any derivatives
    !but we may still want the density profile. So print a warning and set the derivative to 0.
    if(size(input_array) == 1) then
       print *, "discretederivatives.f90: calculate_central_difference: "
       print *, "trying to calculate a discrete derivative with only one point"
       print *, "not possible.  Setting the derivative to 0 and returning as we"
       print *, "may want the density profile."
       calculate_forward_difference = 0.0_dp
       return
    end if

    calculate_forward_difference = 0.0_dp

    !Calculate the forward difference whenever we can
    do ij = 1, size(input_array) - 1
       calculate_forward_difference(ij) = ( input_array(ij + 1) - input_array(ij) ) / &
            ( abs(plate_separations(ij+1) - plate_separations(ij)) * hs_diameter)
    end do

    !On the end point, we can't do the forward difference, therefore do the backward difference.
    calculate_forward_difference(size(input_array)) = &
         ( input_array(size(input_array)) - input_array(size(input_array) - 1) ) / &
         ( abs(plate_separations(size(input_array)) - plate_separations(size(input_array) - 1)) * hs_diameter)

  end function calculate_forward_difference

end module discretederivatives
