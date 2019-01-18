!This module contains the prescription of how to perform an integral
!over the spherical phi coordinate given that the z coordinate in
!cylindrical coordinates in the one that is discretised.
module integratephispherical
  use kinds
  use universalconstants
  use parameters
  use helpers
  implicit none
  private

  public :: integrate_phi_spherical

contains

  !A function that integrates over all possible allowed values of phi
  !in spherical coordinates.  That is, from 0 to pi, with the restriction
  !that we can't integrate through walls. Noting that the integrand is a 1-D disretised
  !function of only z, we follow the algorithm...
  !
  !1. Loop over all values of z_{i} in the integrand array.
  !
  !2. For each z_{i} find the minimum and maximum values of z (z_min, z_max) corresponding
  !to a sphere of radius sigma around the initial z point.
  !
  !3. Form an array subsection of allowed values and apply the trapezoidal rule.  We calculate
  !the jacobian of an associated angle from the difference between z_{i} and the value of z being
  !summed over.
  !
  function integrate_phi_spherical(integrand_array) result(reslt)
    real(dp), dimension(:), intent(in) :: integrand_array
    real(dp), dimension(size(integrand_array)) :: reslt

    integer :: lower_z_limit
    integer :: upper_z_limit

    integer :: relative_z_index
    integer :: ires

    integer :: start_z_index
    integer :: end_z_index

    real(dp) :: surface_area_fraction
    
    !Ensure that we only integrate from hs_diameter/2 up to h - hs_diameter/2
    call get_allowed_z_values(start_z_index, end_z_index, size(reslt))
    
    do ires = start_z_index, end_z_index

       call get_integrand_array_section_limits_and_surface_area_fraction(ires, size(reslt), lower_z_limit, upper_z_limit, &
            relative_z_index, surface_area_fraction)
       
       reslt(ires) = surface_area_fraction * apply_trapezoidal_rule(integrand_array(lower_z_limit:upper_z_limit), relative_z_index)

    end do

    reslt(1:start_z_index-1) = 0.0_dp
    reslt(end_z_index+1:size(reslt)) = 0.0_dp
    
  end function integrate_phi_spherical

  !Returns the Jacobian's phi dependence in spherical coordinates.
  pure function Jacobian_phi_dependence(theta)
    real(dp), intent(in) :: theta
    real(dp) :: Jacobian_phi_dependence

    Jacobian_phi_dependence = sin(theta)

  end function Jacobian_phi_dependence

  !Given an array section 'integrand_array' corresponding to the array section
  !we wish to integrate over and a 'z_index' corresponding to the fixed index that the
  !integral corresponds to, apply the trapezoidal rule.
  function apply_trapezoidal_rule(integrand_array, z_index) result(reslt)
    real(dp), dimension(:), intent(in) :: integrand_array
    integer, intent(in)                :: z_index
    real(dp)                           :: reslt

    real(dp) :: theta_i
    real(dp) :: theta_ip1
    integer  :: i
    reslt = 0.0_dp

    do i = 1, size(integrand_array) - 1

       theta_i =  get_angle_from_z_separation(i, z_index)
       theta_ip1 = get_angle_from_z_separation(i+1, z_index)

       reslt = reslt + ( (integrand_array(i)*Jacobian_phi_dependence(theta_i) + &
            integrand_array(i+1)*Jacobian_phi_dependence(theta_ip1) ) * abs(theta_i - theta_ip1)) &
            /2.0_dp
    end do

  end function apply_trapezoidal_rule

  !Given the fixed value of z from which the functional dependence originates and the value of
  !z being integrated over, we calculate the assocated angle phi (in spherical coordinates)
  !corresponding to their separation.
  function get_angle_from_z_separation(z_integrated, z_fixed) result(angle)
    integer, intent(in) :: z_integrated
    integer, intent(in) :: z_fixed

    real(dp)            :: angle

    if(z_integrated <= z_fixed) then
       angle = pi - acos(real((z_fixed - z_integrated),dp)/real(n_discretised_points_z, dp))
    else
       angle = acos(real((z_integrated - z_fixed),dp)/real(n_discretised_points_z, dp))
    end if

  end function get_angle_from_z_separation
  
  !Given a 'z_index' indicating the value of z that introduces the functional dependence, 'h' = number of valid z
  !discretised values between the plates, this routine calculates the upper and lower possible z limits for the
  !spherical integral and calculates the value of 'z_index' relative to this array section, storing the result in
  !'relative_z_index'.
  subroutine get_integrand_array_section_limits_and_surface_area_fraction(z_index, h, lower_z_limit, upper_z_limit, relative_z_index, surface_area_fraction)
    integer, intent(in)   :: z_index
    integer, intent(in)   :: h
    integer, intent(out)  :: lower_z_limit
    integer, intent(out)  :: upper_z_limit
    integer, intent(out)  :: relative_z_index
    real(dp), intent(out) :: surface_area_fraction

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

    if(z_index < 3*(n_discretised_points_z)/2 + 1) then

       lower_z_limit = lowest_z_calculated
       upper_z_limit = min(z_index + n_discretised_points_z, highest_z_calculated)
       relative_z_index = z_index - lowest_z_calculated + 1

       ! lower_z_limit = 1
       ! upper_z_limit = z_index + n_discretised_points_z
       ! relative_z_index = z_index 

       !surface_area_fraction = 1.0_dp / (1.0_dp - cos(pi - acos((z_index - lowest_z_calculated)/real(n_discretised_points_z, dp))))
       surface_area_fraction = 0.5_dp
    else if((z_index >= (3*n_discretised_points_z/2) + 1) .and.  ((h - z_index) >= 3*n_discretised_points_z/2)) then

       lower_z_limit = z_index - n_discretised_points_z
       upper_z_limit = z_index + n_discretised_points_z
       relative_z_index = n_discretised_points_z + 1

       surface_area_fraction = 0.5_dp
       
    else if((h - z_index) < 3*(n_discretised_points_z)/2) then

       lower_z_limit = max(z_index - n_discretised_points_z, lowest_z_calculated)
       upper_z_limit = highest_z_calculated
       relative_z_index = n_discretised_points_z + 1

       ! lower_z_limit = z_index - n_discretised_points_z
       ! upper_z_limit = h
       ! relative_z_index = n_discretised_points_z + 1

       !surface_area_fraction = 1.0_dp / (1.0_dp + ((highest_z_calculated - z_index)/real(n_discretised_points_z, dp)))
       surface_area_fraction = 0.5_dp

    else
       print *, "integratephispherical.f90:get_integrand_array_section_limits:"
       print *, "Invalid values of z_index/h."
       print *, "z_index must be <= h."
       print *, "Almost certainly a coding error as opposed to an input error...aborting..."
       call abort()
    end if

  end subroutine get_integrand_array_section_limits_and_surface_area_fraction

  
  ! subroutine get_integral_limits(z_index, h, lower_phi_limit, upper_phi_limit)
  !   integer, intent(in)   :: z_index
  !   integer, intent(in)   :: h
  !   real(dp), intent(out) :: lower_phi_limit
  !   real(dp), intent(out) :: upper_phi_limit

  !   if(h < 2*n_discretised_points_z + 1) then
  !      print *, "integratephispherical.f90: get_integral_limits:"
  !      print *, "distance between plates less than twice the hard sphere diameter"
  !      print *, "Do you really want that small a plate separation?"
  !      print *, "These short plate separations are currently not supported...aborting..."
  !      call abort()
  !   end if

  !   if(z_index < n_discretised_points_z + 1) then

  !      lower_phi_limit = 0.0_dp
  !      !upper_phi_limit = pi - acos( ((z_index - 1)*(sigma/N))/sigma )
  !      !where sigma = hs_diameter and N = n_discretised_points_z
  !      upper_phi_limit = pi - acos( real((z_index - 1.0_dp),dp)/ real(n_discretised_points_z,dp) )

  !   else if((z_index >= n_discretised_points_z + 1) .and.  ((h - z_index) >= n_discretised_points_z)) then

  !      lower_phi_limit = 0.0_dp
  !      upper_phi_limit = pi

  !   else if((h - z_index) < n_discretised_points_z) then

  !      !lower_phi_limit = pi - acos( ((h - z_index)*(sigma/N))/sigma )
  !      !where sigma = hs_diameter and N = n_discretised_points_z

  !      lower_phi_limit = acos(real((h - z_index), dp)/real(n_discretised_points_z,dp))
  !      upper_phi_limit = pi
  !   else
  !      print *, "integratephispherical.f90:get_integral_limits:"
  !      print *, "Invalid values of z_index/h."
  !      print *, "z_index must be <= h."
  !      print *, "Almost certainly a coding error as opposed to an input error...aborting..."
  !      call abort()
  !   end if

  ! end subroutine get_integral_limits

end module integratephispherical

