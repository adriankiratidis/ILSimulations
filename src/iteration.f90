
!Contains relevant routines for our iterative procedure, such as convergence check.
module iteration
  use kinds
  use parameters
  use helpers
  implicit none
  private

  public :: converged
  public :: InitialiseDensityDiscretisationAndSetIntegrationAnsatz
  public :: InitialiseVariableDiscretisation
  
contains

  !Routine that returns true iff the average squared difference for all beads is
  !less than the iterative_tolerance input by the user.
  function converged(n1_updated, n1, n2_updated, n2, n3_updated, n3)
    real(dp), dimension(:), intent(in) :: n1_updated
    real(dp), dimension(:), intent(in) :: n1

    real(dp), dimension(:), optional, intent(in) :: n2_updated
    real(dp), dimension(:), optional, intent(in) :: n2

    real(dp), dimension(:), optional, intent(in) :: n3_updated
    real(dp), dimension(:), optional, intent(in) :: n3

    logical                            :: converged

    logical :: n1_converged
    logical :: n2_converged
    logical :: n3_converged

    real(dp) :: av_sq_diff
    real(dp) :: bulk_value

    integer :: start_z_index
    integer :: end_z_index

    converged = .false.
    n1_converged = .false.
    n2_converged = .false.
    n3_converged = .false.

    !Get the range of non-zero elements
    call get_allowed_z_values(start_z_index, end_z_index, size(n1))

    !Ensure that we recieve arrays of the correct size
    if(size(n1_updated) /= size(n1)) then
       print *, "iteration.f90: converged:"
       print *, "size(n1_updated) /= size(n1)"
       print *, "array size mismatch...aborting..."
       call abort()
    else

       av_sq_diff = sum((n1_updated(start_z_index:end_z_index)*(hs_diameter**3) - n1(start_z_index:end_z_index)*(hs_diameter**3))**2)&
            /real(size(n1(start_z_index:end_z_index)),dp)
       bulk_value = sum(n1(start_z_index:end_z_index)*(hs_diameter**3))/real(size(n1(start_z_index:end_z_index)),dp)

       ! if(av_sq_diff == 0.0_dp) then
       !    print *, "iteration.f90: converged: "
       !    print *, "average squared difference between dneisty and updated_density == 0.0_dp"
       !    print *, "would appear to be a coding bug...aborting..."
       !    call abort()
       ! end if
       !print *, "n1 = ", n1
       !print *, "n1 = ", n1_updated

       if(av_sq_diff <= iterative_tolerance*bulk_value) then
          n1_converged = .true.
       else
          n1_converged = .false.
       end if
    end if

    !Now do n2
    if(present(n2_updated)) then

       if(.not. present(n2))then
          print *, "iteration.f90: converged:"
          print *, "optional argument error.  n2_updated present while n2 isn't present"
          print *, "...aborting..."
          call abort()
       else

          if(size(n2_updated) /= size(n2)) then
             print *, "iteration.f90: converged:"
             print *, "size(n2_updated) /= size(n2)"
             print *, "array size mismatch...aborting..."
             call abort()
          else

             av_sq_diff = sum((n2_updated(start_z_index:end_z_index)*(hs_diameter**3) - n2(start_z_index:end_z_index)*(hs_diameter**3))**2)/&
                  size(n2(start_z_index:end_z_index))
             bulk_value = sum(n2(start_z_index:end_z_index)*(hs_diameter**3))/size(n2(start_z_index:end_z_index))

             ! if(av_sq_diff == 0.0_dp) then
             !    print *, "iteration.f90: converged: "
             !    print *, "average squared difference between dneisty and updated_density == 0.0_dp"
             !    print *, "would appear to be a coding bug...aborting..."
             !    !call abort()
             ! end if


             !print *, "n1 = ", n2
             !print *, "n1 = ", n2_updated

             if(av_sq_diff <= iterative_tolerance*bulk_value) then
                n2_converged = .true.
             else
                n2_converged = .false.
             end if
          end if
       end if

    else !if n2 not present, it shouldn't stop convergence
       n2_converged = .true.
    end if

    !Now do n3
    if(present(n3_updated)) then

       if(.not. present(n3)) then
          print *, "iteration.f90: converged:"
          print *, "size(n3_updated) /= size(n3)"
          print *, "array size mismatch...aborting..."
          call abort() 
       else

          if(size(n3_updated) /= size(n3)) then
             print *, "iteration.f90: converged:"
             print *, "size(n3_updated) /= size(n3)"
             print *, "array size mismatch...aborting..."
             call abort()
          else

             av_sq_diff = sum((n3_updated(start_z_index:end_z_index)*(hs_diameter**3) - n3(start_z_index:end_z_index)*(hs_diameter**3))**2)/ &
                  size(n3(start_z_index:end_z_index))
             bulk_value = sum(n3(start_z_index:end_z_index)*(hs_diameter**3))/size(n3(start_z_index:end_z_index))

             if(av_sq_diff <= iterative_tolerance*bulk_value) then
                n3_converged = .true.
             else
                n3_converged = .false.
             end if
          end if
       end if

    else
       n3_converged = .true.
    end if

    !Check for convergence
    if(n1_converged .and. n2_converged .and. n3_converged)then

       print *, "iteration.f90: converged: Iterative scheme successfully converged."
       converged = .true.
    else
       !print *, "convergence = ", n1_converged, n2_converged, n3_converged
       converged = .false.
    end if

  end function converged

  !Intialises the bead density.  The variables are allocated/reallocated to the appropriate size based on which
  !plate separation 'ith_plate-separation' corresponds to.  In the case of it being the first run through
  !(i.e. ith_plate_separation = 1) they are initialised to be constant.  In the case of subsequent run throughs
  !(i.e. ith_plate_separation > 1) they are rescaled based on the previously converged value by the routine
  !'ReScaleArray'.
  subroutine InitialiseDensityDiscretisationAndSetIntegrationAnsatz(ith_plate_separation, n_plus, n_neutral, n_minus, n_plus_cation_end, n_neutral_cation_end, n_minus_cation_end, &
       n_plus_cation_nonend, n_neutral_cation_nonend, n_minus_cation_nonend, n_plus_anion_end, n_neutral_anion_end, n_minus_anion_end)
    integer, intent(in) :: ith_plate_separation
    real(dp), dimension(:), allocatable, intent(inout) :: n_plus
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n_neutral
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n_minus
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n_plus_cation_end
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n_neutral_cation_end
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n_minus_cation_end
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n_plus_cation_nonend
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n_neutral_cation_nonend
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n_minus_cation_nonend
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n_plus_anion_end
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n_neutral_anion_end
    real(dp), dimension(:), allocatable, optional, intent(inout) :: n_minus_anion_end

    integer :: new_array_size

    !Note we add on one to include values at both endpoints/the walls.
    new_array_size = nint(plate_separations(ith_plate_separation) *  n_discretised_points_z) + 1

    call UpdateArraySize(n_plus, new_array_size, rescale=allocated(n_plus), bead_type='+')
    if(present(n_neutral)) call UpdateArraySize(n_neutral, new_array_size, rescale=allocated(n_neutral), bead_type='0')
    if(present(n_minus)) call UpdateArraySize(n_minus, new_array_size, rescale=allocated(n_minus), bead_type='-')

    if(present(n_plus_cation_end)) call UpdateArraySize(n_plus_cation_end, new_array_size, rescale=allocated(n_plus_cation_end), bead_type='pce')
    if(present(n_neutral_cation_end)) call UpdateArraySize(n_neutral_cation_end, new_array_size, rescale=allocated(n_neutral_cation_end), bead_type='nce')
    if(present(n_minus_cation_end)) call UpdateArraySize(n_minus_cation_end, new_array_size, rescale=allocated(n_minus_cation_end), bead_type='mce')
    if(present(n_plus_cation_nonend)) call UpdateArraySize(n_plus_cation_nonend, new_array_size, rescale=allocated(n_plus_cation_nonend), bead_type='pcne')
    if(present(n_neutral_cation_nonend)) call UpdateArraySize(n_neutral_cation_nonend, new_array_size, rescale=allocated(n_neutral_cation_nonend), bead_type='ncne')
    if(present(n_minus_cation_nonend)) call UpdateArraySize(n_minus_cation_nonend, new_array_size, rescale=allocated(n_minus_cation_nonend), bead_type='mcne')
    if(present(n_plus_anion_end)) call UpdateArraySize(n_plus_anion_end, new_array_size, rescale=allocated(n_plus_anion_end), bead_type='pae')
    if(present(n_neutral_anion_end)) call UpdateArraySize(n_neutral_anion_end, new_array_size, rescale=allocated(n_neutral_anion_end), bead_type='nae')
    if(present(n_minus_anion_end)) call UpdateArraySize(n_minus_anion_end, new_array_size, rescale=allocated(n_minus_anion_end), bead_type='mae')

    !We don't calculate with hs_diameter/2 of the wall.  Therefore set it zero.
    !This aids with plotting ease.
    call setNonCalculatedRegionToZero(n_plus)
    if(present(n_neutral)) call setNonCalculatedRegionToZero(n_neutral)
    if(present(n_minus)) call setNonCalculatedRegionToZero(n_minus)
    if(present(n_plus_cation_end)) call setNonCalculatedRegionToZero(n_plus_cation_end)
    if(present(n_neutral_cation_end)) call setNonCalculatedRegionToZero(n_neutral_cation_end)
    if(present(n_minus_cation_end)) call setNonCalculatedRegionToZero(n_minus_cation_end)
    if(present(n_plus_cation_nonend)) call setNonCalculatedRegionToZero(n_plus_cation_nonend)
    if(present(n_neutral_cation_nonend)) call setNonCalculatedRegionToZero(n_neutral_cation_nonend)
    if(present(n_minus_cation_nonend)) call setNonCalculatedRegionToZero(n_minus_cation_nonend)
    if(present(n_plus_anion_end)) call setNonCalculatedRegionToZero(n_plus_anion_end)
    if(present(n_neutral_anion_end)) call setNonCalculatedRegionToZero(n_neutral_anion_end)
    if(present(n_minus_anion_end)) call setNonCalculatedRegionToZero(n_minus_anion_end)
    
  end subroutine InitialiseDensityDiscretisationAndSetIntegrationAnsatz

  !Intialises the all variables (other than the bead densities) that are functions of z.
  !The variables are allocated/reallocated to the appropriate size based on which
  !plate separation 'ith_plate-separation' corresponds to.  They are all subsequently initialised
  !to be a constant.  Note that this constant value is (at the time of writing) always overwritten.
  !
  !Note:
  !In order to support studies with a requirement for a large number of variables for clarity sakes,
  !we need to support an input list of a large number of variables.  As variables argument lists are
  !not allowed in Fortran, as far as I can tell there are two options.
  !
  !1. Explicitly include a long list of optional parameters.
  !2. Put them in an array and pass the single array in.
  !
  !The problem with 2 is that is needs to be an array of an explicit type, which in our case is an
  !array of real numbers.  We therefore need to define a derived type to store this information and hence
  !the main program would look substantially messier due to the need to perform array%val for example.
  !We have therefore chosen to aviod this complication in the main program and choose the unweidly option 1.
  !Perhaps we may switch to option 2 if there ever becomes a requirement for more arguments than we currently
  !have listed.
  subroutine InitialiseVariableDiscretisation(ith_plate_separation,a1,a2,a3,a4,a5,a6,a7, &
       a8,a9,a10,a11,a12,a13,a14,a15,a16,a17,a18,a19,a20,a21,a22,a23,a24,a25,a26,a27,a28, &
       a29,a30,a31,a32,a33,a34,a35,a36,a37,a38,a39,a40,a41,a42,a43,a44,a45,a46,a47,a48,a49, &
       a50,a51,a52,a53,a54,a55,a56,a57,a58,a59,a60,a61,a62,a63)
    integer, intent(in) :: ith_plate_separation
    real(dp), dimension(:), allocatable, intent(inout) :: a1
    real(dp), dimension(:), allocatable, intent(inout), optional :: a2
    real(dp), dimension(:), allocatable, intent(inout), optional :: a3
    real(dp), dimension(:), allocatable, intent(inout), optional :: a4
    real(dp), dimension(:), allocatable, intent(inout), optional :: a5
    real(dp), dimension(:), allocatable, intent(inout), optional :: a6
    real(dp), dimension(:), allocatable, intent(inout), optional :: a7
    real(dp), dimension(:), allocatable, intent(inout), optional :: a8
    real(dp), dimension(:), allocatable, intent(inout), optional :: a9
    real(dp), dimension(:), allocatable, intent(inout), optional :: a10
    real(dp), dimension(:), allocatable, intent(inout), optional :: a11
    real(dp), dimension(:), allocatable, intent(inout), optional :: a12
    real(dp), dimension(:), allocatable, intent(inout), optional :: a13
    real(dp), dimension(:), allocatable, intent(inout), optional :: a14
    real(dp), dimension(:), allocatable, intent(inout), optional :: a15
    real(dp), dimension(:), allocatable, intent(inout), optional :: a16
    real(dp), dimension(:), allocatable, intent(inout), optional :: a17
    real(dp), dimension(:), allocatable, intent(inout), optional :: a18
    real(dp), dimension(:), allocatable, intent(inout), optional :: a19
    real(dp), dimension(:), allocatable, intent(inout), optional :: a20
    real(dp), dimension(:), allocatable, intent(inout), optional :: a21
    real(dp), dimension(:), allocatable, intent(inout), optional :: a22
    real(dp), dimension(:), allocatable, intent(inout), optional :: a23
    real(dp), dimension(:), allocatable, intent(inout), optional :: a24
    real(dp), dimension(:), allocatable, intent(inout), optional :: a25
    real(dp), dimension(:), allocatable, intent(inout), optional :: a26
    real(dp), dimension(:), allocatable, intent(inout), optional :: a27
    real(dp), dimension(:), allocatable, intent(inout), optional :: a28
    real(dp), dimension(:), allocatable, intent(inout), optional :: a29
    real(dp), dimension(:), allocatable, intent(inout), optional :: a30
    real(dp), dimension(:), allocatable, intent(inout), optional :: a31
    real(dp), dimension(:), allocatable, intent(inout), optional :: a32
    real(dp), dimension(:), allocatable, intent(inout), optional :: a33
    real(dp), dimension(:), allocatable, intent(inout), optional :: a34
    real(dp), dimension(:), allocatable, intent(inout), optional :: a35
    real(dp), dimension(:), allocatable, intent(inout), optional :: a36
    real(dp), dimension(:), allocatable, intent(inout), optional :: a37
    real(dp), dimension(:), allocatable, intent(inout), optional :: a38
    real(dp), dimension(:), allocatable, intent(inout), optional :: a39
    real(dp), dimension(:), allocatable, intent(inout), optional :: a40
    real(dp), dimension(:), allocatable, intent(inout), optional :: a41
    real(dp), dimension(:), allocatable, intent(inout), optional :: a42
    real(dp), dimension(:), allocatable, intent(inout), optional :: a43
    real(dp), dimension(:), allocatable, intent(inout), optional :: a44
    real(dp), dimension(:), allocatable, intent(inout), optional :: a45
    real(dp), dimension(:), allocatable, intent(inout), optional :: a46
    real(dp), dimension(:), allocatable, intent(inout), optional :: a47
    real(dp), dimension(:), allocatable, intent(inout), optional :: a48
    real(dp), dimension(:), allocatable, intent(inout), optional :: a49
    real(dp), dimension(:), allocatable, intent(inout), optional :: a50
    real(dp), dimension(:), allocatable, intent(inout), optional :: a51
    real(dp), dimension(:), allocatable, intent(inout), optional :: a52
    real(dp), dimension(:), allocatable, intent(inout), optional :: a53
    real(dp), dimension(:), allocatable, intent(inout), optional :: a54
    real(dp), dimension(:), allocatable, intent(inout), optional :: a55
    real(dp), dimension(:), allocatable, intent(inout), optional :: a56
    real(dp), dimension(:), allocatable, intent(inout), optional :: a57
    real(dp), dimension(:), allocatable, intent(inout), optional :: a58
    real(dp), dimension(:), allocatable, intent(inout), optional :: a59
    real(dp), dimension(:), allocatable, intent(inout), optional :: a60
    real(dp), dimension(:), allocatable, intent(inout), optional :: a61
    real(dp), dimension(:), allocatable, intent(inout), optional :: a62
    real(dp), dimension(:), allocatable, intent(inout), optional :: a63

    integer :: new_array_size

    !Note we add on one to include values at both endpoints/the walls.
    new_array_size = nint(plate_separations(ith_plate_separation) *  n_discretised_points_z) + 1

    call UpdateArraySize(a1, new_array_size)
    if(present(a2)) call UpdateArraySize(a2, new_array_size)
    if(present(a3)) call UpdateArraySize(a3, new_array_size)
    if(present(a4)) call UpdateArraySize(a4, new_array_size)
    if(present(a5)) call UpdateArraySize(a5, new_array_size)
    if(present(a6)) call UpdateArraySize(a6, new_array_size)
    if(present(a7)) call UpdateArraySize(a7, new_array_size)
    if(present(a8)) call UpdateArraySize(a8, new_array_size)
    if(present(a9)) call UpdateArraySize(a9, new_array_size)
    if(present(a10)) call UpdateArraySize(a10, new_array_size)
    if(present(a11)) call UpdateArraySize(a11, new_array_size)
    if(present(a12)) call UpdateArraySize(a12, new_array_size)
    if(present(a13)) call UpdateArraySize(a13, new_array_size)
    if(present(a14)) call UpdateArraySize(a14, new_array_size)
    if(present(a15)) call UpdateArraySize(a15, new_array_size)
    if(present(a16)) call UpdateArraySize(a16, new_array_size)
    if(present(a17)) call UpdateArraySize(a17, new_array_size)
    if(present(a18)) call UpdateArraySize(a18, new_array_size)
    if(present(a19)) call UpdateArraySize(a19, new_array_size)
    if(present(a20)) call UpdateArraySize(a20, new_array_size)
    if(present(a21)) call UpdateArraySize(a21, new_array_size)
    if(present(a22)) call UpdateArraySize(a22, new_array_size)
    if(present(a23)) call UpdateArraySize(a23, new_array_size)
    if(present(a24)) call UpdateArraySize(a24, new_array_size)
    if(present(a25)) call UpdateArraySize(a25, new_array_size)
    if(present(a26)) call UpdateArraySize(a26, new_array_size)
    if(present(a27)) call UpdateArraySize(a27, new_array_size)
    if(present(a28)) call UpdateArraySize(a28, new_array_size)
    if(present(a29)) call UpdateArraySize(a29, new_array_size)
    if(present(a30)) call UpdateArraySize(a30, new_array_size)
    if(present(a31)) call UpdateArraySize(a31, new_array_size)
    if(present(a32)) call UpdateArraySize(a32, new_array_size)
    if(present(a33)) call UpdateArraySize(a33, new_array_size)
    if(present(a34)) call UpdateArraySize(a34, new_array_size)
    if(present(a35)) call UpdateArraySize(a35, new_array_size)
    if(present(a36)) call UpdateArraySize(a36, new_array_size)
    if(present(a37)) call UpdateArraySize(a37, new_array_size)
    if(present(a38)) call UpdateArraySize(a38, new_array_size)
    if(present(a39)) call UpdateArraySize(a39, new_array_size)
    if(present(a40)) call UpdateArraySize(a40, new_array_size)
    if(present(a41)) call UpdateArraySize(a41, new_array_size)
    if(present(a42)) call UpdateArraySize(a42, new_array_size)
    if(present(a43)) call UpdateArraySize(a43, new_array_size)
    if(present(a44)) call UpdateArraySize(a44, new_array_size)
    if(present(a45)) call UpdateArraySize(a45, new_array_size)
    if(present(a46)) call UpdateArraySize(a46, new_array_size)
    if(present(a47)) call UpdateArraySize(a47, new_array_size)
    if(present(a48)) call UpdateArraySize(a48, new_array_size)
    if(present(a49)) call UpdateArraySize(a49, new_array_size)
    if(present(a50)) call UpdateArraySize(a50, new_array_size)
    if(present(a51)) call UpdateArraySize(a51, new_array_size)
    if(present(a52)) call UpdateArraySize(a52, new_array_size)
    if(present(a53)) call UpdateArraySize(a53, new_array_size)
    if(present(a54)) call UpdateArraySize(a54, new_array_size)
    if(present(a55)) call UpdateArraySize(a55, new_array_size)
    if(present(a56)) call UpdateArraySize(a56, new_array_size)
    if(present(a57)) call UpdateArraySize(a57, new_array_size)
    if(present(a58)) call UpdateArraySize(a58, new_array_size)
    if(present(a59)) call UpdateArraySize(a59, new_array_size)
    if(present(a60)) call UpdateArraySize(a60, new_array_size)
    if(present(a61)) call UpdateArraySize(a62, new_array_size)
    if(present(a62)) call UpdateArraySize(a62, new_array_size)
    if(present(a63)) call UpdateArraySize(a63, new_array_size)

  end subroutine InitialiseVariableDiscretisation

  !Routine that either resizes or rescales the array, based on the presence and value
  !of the optional parameter 'rescale'.  Also, initialises the array values.
  subroutine UpdateArraySize(array, new_size, rescale, bead_type)
    real(dp), dimension(:), allocatable, intent(inout) :: array
    integer                                            :: new_size
    logical, intent(in), optional                      :: rescale
    character(len=*), intent(in), optional             :: bead_type

    if(allocated(array)) then

       if(size(array) == new_size) then
          print *, "iteration.f90:UpdateArraySize:"
          print *, "Attempting to repeat the calculation at the same plate separation."
          print *, "There is no reason to do this."
          print *, "Please ensure the list of plate separations do not contain duplicates"
          call abort()   
       else

          if(present(rescale)) then

             if(rescale) then
                call ReScaleArray(array, new_size)
             else
                call ReSizeArray(array, new_size)
             end if

          else
             !call ReSizeArray(array, new_size)
             call ReScaleArray(array, new_size)
          end if

       end if

    else
       if(present(rescale)) then
          if(.not. present(bead_type)) then
             print *, "iteration.f90: UpdateArraySize:"
             print *, "If present(rescale) then bead_type must also be present...aborting..."
             call abort()
          end if
          if(rescale) then
             print *, "iteration.f90: UpdateArraySize:"
             print *, "Can't rescale an array that has not yet been allocated"
             call abort()
          end if
       end if

       allocate(array(new_size)) !size = M x N
       if(.not. present(bead_type)) then
          array(:) = 0.0_dp
       else
          call InitialiseIntegrationAnsatz(array, trim(bead_type))
       end if
       
    end if

  end subroutine UpdateArraySize

  !Routine that resizes 'array' to 'new_size' and initialises all elements to 0.
  subroutine ReSizeArray(array, new_size)
    real(dp), dimension(:), allocatable, intent(inout) :: array
    integer, intent(in)                                :: new_size

    if(size(array) == new_size) then
       print *, "iteration.f90:ReSizeArray:"
       print *, "The size of the array is equal to the proposed new size."
       print *, "No need to call resize."
       print *, "Probable coding error...aborting..."
       call abort()
    end if

    if(allocated(array)) deallocate(array)
    allocate(array(new_size))
    array(:) = 0.0_dp

  end subroutine ReSizeArray

  !Rescales 'array' to an array of size 'new_size'.  The ith value of the new array
  !(the array value at the end of the routine) is given by the value of the old array
  !(the array value at the start of the routine) that is the closest to the same fraction
  !of the total size as the value in the new array that it is setting.
  subroutine ReScaleArray(array, new_size)
    real(dp), dimension(:), allocatable, intent(inout) :: array
    integer, intent(in)                                :: new_size

    real(dp), dimension(size(array)) :: old_array_values
    integer :: ith_component
    integer :: old_value_index

    integer :: start_z_index_new
    integer :: end_z_index_new

    integer :: start_z_index_old
    integer :: end_z_index_old
    
    if(.not. allocated(array)) then
       print *, "iteration.f90: Unable to rescale array that isn't allocated"
       print *, "This is almost certainly a coding error"
       print *, "as opposed to an error with the input parms"
       call abort()
    end if

    !Store the values of the array to be deallocated
    old_array_values = array

    !Reallocate the array to be the new size
    if(allocated(array)) deallocate(array)
    allocate(array(new_size))
    array = 0.0_dp

    !Find the range over which the values are non zero
    call get_allowed_z_values(start_z_index_new, end_z_index_new, size(array))
    call get_allowed_z_values(start_z_index_old, end_z_index_old, size(old_array_values))
    
    !Set the elements of the new array by rescaling the old values
    do ith_component = start_z_index_new, end_z_index_new
       old_value_index = start_z_index_old - 1 +  &
            ceiling((end_z_index_old - start_z_index_old + 1) * ( real(ith_component - start_z_index_new + 1, dp) &
            / real(end_z_index_new - start_z_index_new + 1, dp) ) )
       array(ith_component) = old_array_values(old_value_index)
    end do
    !array(:) = bulk_density_positive_beads
  end subroutine ReScaleArray

  !Routine to initialise our ansatz for our integrative scheme to an arbitrary constant.
  subroutine InitialiseIntegrationAnsatz(array, bead_type)
    real(dp), dimension(:) :: array
    character(len=*), intent(in) :: bead_type

    integer :: midpoint
    integer :: ij

    !print *, "bead_type = ", bead_type
    !print *, "size(array) = ", size(array), size(array)/2


    if(trim(bead_type) == '+') then

       midpoint = size(array)/2

       do ij = 1, size(array)
          array(ij) = bulk_density_positive_beads + ((ij - midpoint)*hs_diameter/real(n_discretised_points_z,dp))*slope_for_initial_guess
       end do

    else if(trim(bead_type) == '0') then
       array(:) = bulk_density_neutral_beads

    else if(trim(bead_type) == '-') then

       midpoint = size(array)/2

       do ij = 1, size(array)
          array(ij) = bulk_density_negative_beads - ((ij - midpoint)*hs_diameter/real(n_discretised_points_z, dp))*slope_for_initial_guess
       end do

    else if(trim(bead_type) == 'pce') then !bead_type = positive bead in cation end density
       array(:) = n_plus_cation_end_bulk

    else if(trim(bead_type) == 'nce') then !bead_type = neutral bead in cation end density
       array(:) = n_neutral_cation_end_bulk

    else if(trim(bead_type) == 'mce') then !bead_type = minus bead in cation end density
       array(:) = n_minus_cation_end_bulk

    else if(trim(bead_type) == 'pcne') then !bead_type = positive bead in cation nonend density
       array(:) = n_plus_cation_nonend_bulk

    else if(trim(bead_type) == 'ncne') then !bead_type = neutral bead in cation nonend density
       array(:) = n_neutral_cation_nonend_bulk

    else if(trim(bead_type) == 'mcne') then !bead_type = minus bead in cation nonend density
       array(:) = n_minus_cation_nonend_bulk

    else if(trim(bead_type) == 'pae') then !bead_type = positive bead in anion end density
       array(:) = n_plus_anion_end_bulk

    else if(trim(bead_type) == 'nae') then !bead_type = neutral bead in anion end density
       array(:) = n_neutral_anion_end_bulk

    else if(trim(bead_type) == 'mae') then !bead_type = minus bead in anion end density
       array(:) = n_minus_anion_end_bulk

    else
       print *, "iteration.f90: InitialiseIntegrationAnsatz"
       print *, "bead_type has an illegal value of ", trim(bead_type), "...aborting..."
       call abort()
    end if
    !print *, "printing bead_type and array" 
    !print *, bead_type, array(midpoint)*(hs_diameter**3), array*(hs_diameter**3)

    !array(:) = bulk_density_positive_beads
  end subroutine InitialiseIntegrationAnsatz

end module iteration
