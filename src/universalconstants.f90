!A module containing the values of universal constants.
module universalconstants
  use kinds
  implicit none
  public

  real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp ! pi in double precision.
  real(dp), parameter :: epsilon0=8.854187817620389850536563031710750260608_dp ! epsilon_0 in double precision
  real(dp), parameter :: k_B=1.38064852E-23_dp

  complex(sc), parameter :: i_sc = ( 0.0_sp, 1.0_sp ) ! The square root of -1, in single precision.
  complex(dc), parameter :: i_dc = ( 0.0_dp, 1.0_dp ) ! The square root of -1, in double precision.

end module universalconstants
  
