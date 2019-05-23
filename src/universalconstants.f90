!A module containing the values of universal constants.
module universalconstants
  use kinds
  implicit none
  public

  real(dp), parameter :: pi=3.141592653589793238462643383279502884197_dp ! pi in double precision.
  real(dp), parameter :: epsilon0=8.854187817620389850536563031710750260608E-22_dp ! epsilon_0 in double precision [Farads/Angstrom] ([Farad] = [C^2]/[J])
  real(dp), parameter :: k_B=1.38064852E-23_dp !([J]/[K])
  real(dp), parameter :: electric_charge=1.6021766E-19_dp !the electric charge ([C])

  complex(sc), parameter :: i_sc = ( 0.0_sp, 1.0_sp ) ! The square root of -1, in single precision.
  complex(dc), parameter :: i_dc = ( 0.0_dp, 1.0_dp ) ! The square root of -1, in double precision.

end module universalconstants
  
