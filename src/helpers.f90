! Module containing so-called helper functions.
! Used to perform neccesary logic on data structures.
module helpers
  use kinds
  use parameters
  implicit none
  private

  public :: get_allowed_z_values
  public :: setNonCalculatedRegionToZero
  public :: unity_function
  public :: get_bulk_density
  public :: str
  public :: SetToZero
  public :: SolveWithNewtonsMethod

  interface str
     module procedure str_int
     module procedure str_real_dp
     module procedure str_real_sp
  end interface str

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

    !start_z_index = 1
    !end_z_index = total_points_z

  end subroutine get_allowed_z_values

  subroutine setNonCalculatedRegionToZero(array1, array2, array3)
    real(dp), dimension(:), intent(inout) :: array1
    real(dp), dimension(:), intent(inout), optional :: array2
    real(dp), dimension(:), intent(inout), optional :: array3

    integer :: start_z_index
    integer :: end_z_index

    !First ensure all input arguments are the same size.
    if(present(array2)) then
       if(size(array1) /= size(array2)) then
          print *, "helpers.f90: setNonCalculatedRegionToZero: "
          print *, "size mismatch.  Size of all input arguments must match."
          print *, "size(array1) /= size(array2)...aborting..."
          call abort()
       end if
    end if

    if(present(array3)) then
       if(size(array1) /= size(array3)) then
          print *, "helpers.f90: setNonCalculatedRegionToZero: "
          print *, "size mismatch.  Size of all input arguments must match."
          print *, "size(array1) /= size(array3)...aborting..."
          call abort()
       end if
    end if

    call get_allowed_z_values(start_z_index, end_z_index, size(array1))

    array1(1:start_z_index-1) = 0.0_dp
    array1(end_z_index+1:size(array1)) = 0.0_dp

    if(present(array2)) then
       array2(1:start_z_index-1) = 0.0_dp
       array2(end_z_index+1:size(array2)) = 0.0_dp
    end if

    if(present(array3)) then
       array3(1:start_z_index-1) = 0.0_dp
       array3(end_z_index+1:size(array3)) = 0.0_dp
    end if

  end subroutine setNonCalculatedRegionToZero

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

  pure function unity_function(z, xi)
    integer, intent(in) :: z
    integer, intent(in) :: xi
    real(dp) :: unity_function

    unity_function = 1.0_dp
  end function unity_function

  function get_bulk_density(lambda) result(reslt)
    real(dp), dimension(:), intent(in)  :: lambda
    real(dp) :: reslt

    integer :: start_z_integrate
    integer :: end_z_integrate

    !Find the maximum range over which we are going to integrate
    call  get_allowed_z_values(start_z_integrate, end_z_integrate, size(lambda))

    reslt = sum(lambda(start_z_integrate:end_z_integrate))/real(size(lambda(start_z_integrate:end_z_integrate)), dp)

    !Setting the bulk density to be int(lambda)/plate_separation
    !We could of course calculate the total plate separation by hs_diameter * plate_separations(ith_separation)
    !but we choose the current version so we can calculate it without passing in an extra parameter.
    ! result_array = integrate_z_cylindrical(lambda, unity_function, "all_z") / &
    !      ( (real(size(lambda(start_z_integrate:end_z_integrate)) - 1, dp) * hs_diameter)/real(n_discretised_points_z,dp) )

    ! if(count(result_array(start_z_integrate:end_z_integrate) /= result_array(start_z_integrate)) /= 0) then
    !    print *, "lambdas.f90: get_bulk_density:"
    !    print *, "more than one different value when calculating bulk density"
    !    print *, "almost certainly coding error...aborting..."
    !    call abort
    ! else
    !    reslt = result_array(start_z_integrate)
    ! end if

    return
  end function get_bulk_density

  !Casts an int to a string
  function str_int(x)
    integer, intent(in) :: x
    character(len=64) :: str_int

    write(str_int, *) x

    str_int = adjustl(str_int)
  end function str_int

  !Casts a double precision real to a string
  function str_real_dp(x)
    real(dp), intent(in) :: x
    character(len=64) :: str_real_dp

    write(str_real_dp, '(f10.5)') x

    str_real_dp = adjustl(str_real_dp)
  end function str_real_dp

  !Casts a single precision real to a string
  function str_real_sp(x)
    real(sp), intent(in) :: x
    character(len=64) :: str_real_sp

    write(str_real_sp, '(f10.5)') x

    str_real_sp = adjustl(str_real_sp)
  end function str_real_sp

  subroutine SetToZero(array1, array2, array3, array4, array5, array6)
    real(dp), dimension(:) :: array1
    real(dp), dimension(:), optional :: array2
    real(dp), dimension(:), optional :: array3
    real(dp), dimension(:), optional :: array4
    real(dp), dimension(:), optional :: array5
    real(dp), dimension(:), optional :: array6

    array1(:) = 0.0_dp
    if(present(array2)) array2(:) = 0.0_dp
    if(present(array3)) array3(:) = 0.0_dp
    if(present(array4)) array4(:) = 0.0_dp
    if(present(array5)) array5(:) = 0.0_dp
    if(present(array6)) array6(:) = 0.0_dp

  end subroutine SetToZero

  !Given equation to solve of the form
  !
  !a*exp(\beta*{q_+}*Donnan_potential) + b*exp(\beta*{q_-}*Donnan_potential) + c = 0
  !
  subroutine SolveWithNewtonsMethod(a,b,c,Donnan_potential, abort_now)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp), intent(in) :: c
    real(dp), intent(inout) :: Donnan_potential
    logical :: abort_now

    real(dp) :: tolerance
    integer :: ij

    real(dp) :: xi, xip1
    logical :: converged

    converged = .false.
    tolerance = abs(0.0000001_dp * electric_charge)
    xi = Donnan_potential !Start with initial guess = previous Donnan potential.

    do ij = 1, MAX_ITERATION_LIMIT

       xip1 = xi - Donnan_func_to_solve(a,b,c,xi)/Donnan_func_to_solve_deriv(a,b,xi)

       !print *, Donnan_func_to_solve(a,b,c,xi), Donnan_func_to_solve_deriv(a,b,xi), Donnan_func_to_solve(a,b,c,xi)/Donnan_func_to_solve_deriv(a,b,xi)
       !print *, xip1, xi, abs(xip1 - xi)
       !print *, a, b, c, abs(a*exp(beta*xip1*positive_oligomer_charge) + b*exp(beta*xip1*negative_oligomer_charge) + c), tolerance
       !print *, a*exp(beta*xip1*positive_oligomer_charge)/positive_bead_charge, b*exp(beta*xip1*negative_oligomer_charge)/negative_bead_charge
       !print *, ""

       !print *, beta*positive_oligomer_charge, beta*negative_oligomer_charge
       !call abort()

       if(isnan(xip1)) then
          print *, "found a nan"
          print *, xip1, xi, abs(xip1 - xi)
          print *, a, b, c, abs(a*exp(beta*xip1*positive_oligomer_charge) + b*exp(beta*xip1*negative_oligomer_charge) + c), tolerance
          print *, a*exp(beta*xip1*positive_oligomer_charge)/positive_bead_charge, 5*b*exp(beta*xip1*negative_oligomer_charge)/negative_bead_charge
          print *, ""
          print *, "found a nan"
          !Found a problem, but still want to print out what we've calculated so far.
          abort_now = .true.
          !call abort()
       end if

       !if((abs(xip1 - xi) < tolerance)) then

       if(abs(Donnan_func_to_solve(a,b,c,xip1)) < tolerance) then
          !print *, Donnan_func_to_solve(a,b,c,xip1), "sadkjfh"
          if(Donnan_func_to_solve(a,b,c,xip1) > 1.0E-4_dp) then
             print *, "helpers.f90:SolveWithNewtonsMethod: Converged to wrong solution"
             print *, "in calculation of the Donnan potential."
             print *, "Function value (which should be zero) = ", Donnan_func_to_solve(a,b,c,xip1)
             !Found a problem, but still want to print out what we've calculated so far.
             abort_now = .true.
             !call abort()
          end if

          converged = .true.
          exit
       else
          xi = xip1
       end if
    end do

    if(converged) then
       Donnan_potential = xip1
       !print *, "value of Donnan function = ",
       !print *, xip1, xi, abs(xip1 - xi)
       !print *, a, b, c, abs(a*exp(beta*xip1*positive_oligomer_charge) + b*exp(beta*xip1*negative_oligomer_charge) + c), tolerance
       !print *, a*exp(beta*xip1*positive_oligomer_charge)/positive_bead_charge, 5*b*exp(beta*xip1*negative_oligomer_charge)/negative_bead_charge
       !print *, ""

    else
       print *, "helpers.f90:SolveWithNewtonsMethod: Hit the iteration limit before convergence"
       print *, "in calculation of Donnan potential"
       print *, xip1, xi, abs(xip1 - xi)
       print *, a, b, c, abs(a*exp(beta*xip1*positive_oligomer_charge) + b*exp(beta*xip1*negative_oligomer_charge) + c), tolerance
       print *, a*exp(beta*xip1*positive_oligomer_charge)/positive_bead_charge, 5*b*exp(beta*xip1*negative_oligomer_charge)/negative_bead_charge
       print *, ""
       !Found a problem, but still want to print out what we've calculated so far.
       abort_now = .true.
       !call abort()
    end if

  end subroutine SolveWithNewtonsMethod

  function Donnan_func_to_solve(a,b,c,Donnan_potential)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp), intent(in) :: c
    real(dp), intent(in) :: Donnan_potential
    real(dp) :: Donnan_func_to_solve

    Donnan_func_to_solve = a*exp(beta*positive_oligomer_charge*Donnan_potential) + b*exp(beta*negative_oligomer_charge*Donnan_potential) + c
    
    ! Donnan_func_to_solve = a*exp(beta*Donnan_potential*(positive_bead_charge - negative_bead_charge)) &
    !      + c*exp(-1.0_dp*beta*Donnan_potential*negative_bead_charge) + b

  end function Donnan_func_to_solve

  function Donnan_func_to_solve_deriv(a,b,Donnan_potential)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: b
    real(dp), intent(in) :: Donnan_potential
    real(dp) :: Donnan_func_to_solve_deriv

    Donnan_func_to_solve_deriv = a*beta*positive_oligomer_charge*exp(beta*positive_oligomer_charge*Donnan_potential) &
         + b*beta*negative_oligomer_charge*exp(beta*negative_oligomer_charge*Donnan_potential) 
    
    
    ! Donnan_func_to_solve_deriv = a*exp(beta*Donnan_potential*(positive_bead_charge - negative_bead_charge))*beta*(positive_bead_charge - negative_bead_charge) &
    !      - c*exp(-1.0_dp*beta*Donnan_potential*negative_bead_charge)*beta*negative_bead_charge

  end function Donnan_func_to_solve_deriv

end module helpers
