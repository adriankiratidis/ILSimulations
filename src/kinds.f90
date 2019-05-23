!A module that stores kind type parameters corresponding to the various precisions we may need.
module Kinds
  implicit none
  public

  integer, parameter :: i32  = selected_int_kind(9)  ! 32 bit integer.
  integer, parameter :: i64  = selected_int_kind(18) ! 64 bit integer.

  integer, parameter :: sp = kind(1.0) ! Single precision real scalars.
  integer, parameter :: dp = kind(1.0D0) ! Double precision real scalars.
  !kind(1.0D0) ! Double precision real scalars.
  !selected_real_kind(33, 4931) !Quadruple precision real scalars

  integer, parameter :: sc = kind((1.0,1.0)) ! Single precision complex scalars.
  integer, parameter :: dc = kind((1.0D0,1.0D0)) ! Double precision complex scalars.
  
end module Kinds
