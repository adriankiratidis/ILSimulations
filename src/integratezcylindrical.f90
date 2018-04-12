!A module containing the routines to integrate over all z
!in cylindrical coordinates.  Note that it is this variable
!that is discretised.  Also note that the Jacobian
!from  cylindrical coordinates is not present as the 'r'
!direction has already been integrated out.
module integratezcylindrical
  use kinds
  use parameters
  use universalconstants
  use helpers
  implicit none
  private

  public :: integrate_z_cylindrical
  
contains

  !When calculating the lambdas for example the integrand array is density dependent, while the integrand function
  !is density independent.
  function integrate_z_cylindrical(integrand_array, integrand_function, integration_range) result(reslt)
    real(dp), dimension(:), intent(in) :: integrand_array
    real(dp), external                 :: integrand_function
    character(len=*)                   :: integration_range

    real(dp), dimension(size(integrand_array)) :: reslt

    integer :: ires
    integer :: lower_z_limit
    integer :: upper_z_limit
    integer :: relative_z_index

    integer :: start_z_index
    integer :: end_z_index

    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    call get_allowed_z_values(start_z_index, end_z_index, size(reslt))

    do ires = start_z_index, end_z_index

       call get_integrand_array_section_limits(trim(integration_range), ires, size(reslt), &
            lower_z_limit, upper_z_limit, relative_z_index)

       !print *,"helloooooooooooooooo = ", ires, lower_z_limit, upper_z_limit, relative_z_index

       reslt(ires) = apply_trapezoidal_rule(integrand_array(lower_z_limit:upper_z_limit), &
            integrand_function, relative_z_index)

       !print *, reslt(ires)

    end do

    reslt(1:start_z_index-1) = 0.0_dp
    reslt(end_z_index+1:size(reslt)) = 0.0_dp
    return
    !print *, "reslt = ", reslt

  end function integrate_z_cylindrical


  !Given an array section 'z_dep_integrand' corresponding to the array section
  !we wish to integrate over and a 'z_index' corresponding to the fixed index that the
  !integral corresponds to, apply the trapezoidal rule.
  function apply_trapezoidal_rule(integrand_array, integrand_function, z_index) result(reslt)
    real(dp), dimension(:), intent(in) :: integrand_array
    real(dp), external                 :: integrand_function
    integer, intent(in)                :: z_index
    real(dp)                           :: reslt

    integer :: ixi
    reslt = 0.0_dp

    do ixi = 1, size(integrand_array) - 1
       reslt = reslt + ( (integrand_array(ixi)*integrand_function(z_index, ixi) +  &
            integrand_array(ixi + 1)*integrand_function(z_index, ixi + 1)) &
            * (hs_diameter/real(n_discretised_points_z,dp)))/2.0_dp
    end do

    return
  end function apply_trapezoidal_rule


  !Given an 'integration_range' that tells us whether we are integrating over everything or values >= or <= the
  !hs_diameter, a 'z_index' indicating the value of z that introduces the functional dependence and 'h' = number of valid z
  !discretised values between the plates, this routine calculates the upper and lower possible z limits for the
  !integral and calculates the value of 'z_index' relative to this array section, storing the result in
  !'relative_z_index'.
  subroutine get_integrand_array_section_limits(integration_range, z_index, h, lower_z_limit, upper_z_limit, relative_z_index)
    character(len=*), intent(in) :: integration_range
    integer, intent(in)          :: z_index
    integer, intent(in)          :: h
    integer, intent(out)         :: lower_z_limit
    integer, intent(out)         :: upper_z_limit
    integer, intent(out)         :: relative_z_index

    integer :: lowest_z_calculated
    integer :: highest_z_calculated

    call get_allowed_z_values(lowest_z_calculated, highest_z_calculated, h)
    
    if(h < 2*n_discretised_points_z + 1) then
       print *, "integratephispherical.f90: get_integrand_array_section_limits:"
       print *, "distance between plates less than twice the hard sphere diameter"
       print *, "Do you really want that small a plate separation?"
       print *, "These short plate separations are currently not supported...aborting..."
       call abort()
    end if

    if(trim(integration_range) == "all_z" .or. &
         trim(integration_range) == "z_gteq_hs_diameter") then

       lower_z_limit = lowest_z_calculated
       upper_z_limit = highest_z_calculated
       relative_z_index = z_index - lowest_z_calculated + 1

    else if(trim(integration_range) == "z_lteq_hs_diameter") then

       if(z_index < n_discretised_points_z + 1) then !0 <= z < hs_diameter

          lower_z_limit = lowest_z_calculated
          upper_z_limit = z_index + n_discretised_points_z
          relative_z_index = z_index - lowest_z_calculated + 1

       else if((z_index >= n_discretised_points_z + 1) .and.  ((h - z_index) >= n_discretised_points_z)) then !hs_diameter <= z <= h-hs_diameter

          lower_z_limit = z_index - n_discretised_points_z
          upper_z_limit = z_index + n_discretised_points_z
          relative_z_index = n_discretised_points_z + 1

       else if((h - z_index) < n_discretised_points_z) then !h-hs_diameter < z <= h

          lower_z_limit = z_index - n_discretised_points_z
          upper_z_limit = highest_z_calculated
          relative_z_index = n_discretised_points_z + 1

       else ! z > h which is unphysical
          print *, "integratephispherical.f90:get_integrand_array_section_limits:"
          print *, "Invalid values of z_index/h."
          print *, "z_index must be <= h."
          print *, "Almost certainly a coding error as opposed to an input error...aborting..."
          call abort()
       end if

    else
       print *, "integratezcylindrical.f90:get_integrand_array_section_limits: "
       print *, "unsupported value of 'integration_range'"
       print *, "integration range = ", trim(integration_range)
       print *, "Almost certainly a coding error as opposed to an input error"
       print *, "Aborting..."
       call Abort()
    end if
  end subroutine get_integrand_array_section_limits

end module integratezcylindrical
