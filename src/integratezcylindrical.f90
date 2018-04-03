!A module containing the routines to integrate over all z
!in cylindrical coordinates.  Note that it is this variable
!that is discretised.  Also note that the Jacobian
!from  cylindrical coordinates is not present as the 'r'
!direction has already been integrated out.
module integratezcylindrical
  use kinds
  use parameters
  use universalconstants
  implicit none
  private

  public :: integrate_z_cylindrical
  
contains
  
  function integrate_z_cylindrical(z_dep_integrand, integration_range) result(reslt)
    real(dp), dimension(:), intent(in) :: z_dep_integrand
    character(len=*), intent(in) :: integration_range
    real(dp), dimension(size(z_dep_integrand)) :: reslt

    if(trim(integration_range) == "all_z") then !integrate over all possible values of z
       
       
    else if(trim(integration_range) == "z_lteq_hs_diameter") then !integrate over all z <= hs_diameter
       
    else if(trim(integration_range) == "z_gteq_hs_diameter") then !integrate over all z >= hs_diameter

    end if


  end function integrate_z_cylindrical


  !Given an array section 'z_dep_integrand' corresponding to the array section
  !we wish to integrate over and a 'z_index' corresponding to the fixed index that the
  !integral corresponds to, apply the trapezoidal rule.
  function apply_trapezoidal_rule(z_dep_integrand) result(reslt)
    real(dp), dimension(:), intent(in) :: z_dep_integrand
    real(dp)                           :: reslt

    integer :: iz
    reslt = 0.0_dp

    do iz = 1, size(z_dep_integrand) - 1
       reslt = ( (z_dep_integrand(iz) +  z_dep_integrand(iz + 1)) * (hs_diameter/n_discretised_points_z))/2.0_dp
    end do

  end module integratezcylindrical
