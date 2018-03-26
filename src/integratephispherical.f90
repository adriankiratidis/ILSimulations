!This module contains the prescription of how to perform an integral
!over the spherical phi coordinate given that the z coordinate in
!cylindrical coordinates in the one that is discretised.
module integratephispherical
  use kinds
  implicit none
  private

  public :: integrate_phi_spherical

contains

  function integrate_phi_spherical(z_dep_integrand) result(res)
    real(dp), dimension(:), intent(in) :: z_dep_integrand
    real(dp), dimension(size(z_dep_integrand)) :: res
    
    res = 0
    
  end function integrate_phi_spherical

  function Jacobian(r, theta)
    real(dp), intent(in) :: r
    real(dp), intent(in) :: theta
    real(dp) :: Jacobian

    Jacobian = (r**2) * sin(theta)
    
  end function Jacobian

  ! function trapezoidal_rule()

  ! end function trapezoidal_rule

  ! subroutine get_integral_limits(z)
    
  ! end subroutine get_integral_limits
  
end module integratephispherical
  
