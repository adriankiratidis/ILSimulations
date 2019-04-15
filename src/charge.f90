module charge
  use kinds
  use parameters
  use helpers
  use integratezcylindrical
  implicit none
  private

  public :: ImposeChargeNeutrality
  public :: CalculateDonnanPotential
  
  public :: InitialiseChargeIncrement
  public :: UpdateChargeIncrement
  
  real(dp) :: positive_bead_increment
  real(dp) :: negative_bead_increment

  real(dp) :: left_plate_increment
  real(dp) :: right_plate_increment
  
contains
  
  subroutine InitialiseChargeIncrement()

    positive_bead_increment = positive_bead_charge / real(n_charge_iterations,dp)
    negative_bead_increment = negative_bead_charge / real(n_charge_iterations,dp)

    positive_bead_charge = positive_bead_increment 
    negative_bead_charge = negative_bead_increment

    !left_plate_increment = surface_charge_density_left_wall / real(n_charge_iterations,dp)
    !right_plate_increment = surface_charge_density_right_wall / real(n_charge_iterations,dp)

    !surface_charge_density_left_wall = left_plate_increment
    !surface_charge_density_right_wall = right_plate_increment
    
  end subroutine InitialiseChargeIncrement

  
  subroutine UpdateChargeIncrement()
    
    positive_bead_charge = positive_bead_charge + positive_bead_increment
    negative_bead_charge = negative_bead_charge + negative_bead_increment

    !surface_charge_density_left_wall = surface_charge_density_left_wall + left_plate_increment
    !surface_charge_density_right_wall = surface_charge_density_right_wall + right_plate_increment
    
  end subroutine UpdateChargeIncrement

  subroutine ImposeChargeNeutrality(n_plus, n_neutral, n_minus, Donnan_potential, abort_now)
    real(dp), dimension(:) :: n_plus
    real(dp), dimension(:) :: n_neutral
    real(dp), dimension(:) :: n_minus
    real(dp) :: Donnan_potential
    logical :: abort_now

    real(dp) :: a, b, c, d
    real(dp) :: y1, y2 !two solutions for two roots
    !real(dp) :: Donnan_potential

    !print *, "n_plus = ", n_plus
    !print *, "n_minus = ", n_minus

    !print *, "positive integral 1= ", integrate_z_cylindrical(n_plus, "all_z")
    !print *, "neutral integral 1= ", integrate_z_cylindrical(n_neutral, "all_z")

    !In order to solve for the Donnan potential and impose charge neutrality we must
    !solve the equation
    !
    !

    !if( ((surface_charge_density_left_wall + surface_charge_density_right_wall)/electric_charge) < 1.0E-10_dp) then
    !The walls have equal and opposite charge.  By symmetry arguments the density profiles are already electroneutral.
    !=> Donnan potential is zero.  Therefore, there is nothing to do.
    !else
    if(abs((positive_bead_charge + negative_bead_charge)/electric_charge) < 1.0E-10_dp) then
       !We have the special case of the positive and negative bead charge being the same.
       !We can solve this analytically, as q_{-} = -q_{+}.

       a = integrate_z_cylindrical(n_plus*positive_bead_charge, "all_z")
       b = integrate_z_cylindrical(n_minus*negative_bead_charge, "all_z")
       c = surface_charge_density_left_wall + surface_charge_density_right_wall

       !print *, "a,b,c", a, b, c

       !Finding the Donnan potential to impose electroneutrality is now equivalent to solving
       !a*exp(x) + b*exp(-x) + c = 0; where x = exp(beta*positive_bead_charge*Psi_{Donnan}).

       !Letting y = exp(x), and using the quadratic fromula we get

       y1 = (-1.0_dp*c + sqrt(c**2 - 4.0_dp*a*b))/(2.0_dp*a)
       y2 = (-1.0_dp*c - sqrt(c**2 - 4.0_dp*a*b))/(2.0_dp*a)

       if(y1 < 0 .and. y2 < 0) then
          print *, "charge.f90: ImposeChargeNeutrality: "
          print *, "Both solutions for y1 and y2 are negative."
          print *, "Since this is supposed to be exp(x), this is not possible...aborting..."
          call abort()
       else if(y1 > 0 .and. y2 > 0) then
          print *, "charge.f90: ImposeChargeNeutrality: "
          print *, "Both solutions for y1 and y2 are positive."
          print *, "Don't know which one to pick...aborting..."
          call abort()
       else
          y1 = max(y1, y2)
       end if

       Donnan_potential = log(y1)/(beta*positive_oligomer_charge)
       !print *, y1, y2
       !print *, "Donnan = ", Donnan_potential
       !print *, "Donnan potential = ", exp(beta*positive_oligomer_charge*Donnan_potential)
       !call abort()
       n_plus = n_plus*exp(beta*positive_oligomer_charge*Donnan_potential)
       n_neutral = n_neutral*exp(beta*positive_oligomer_charge*Donnan_potential)
       n_minus = n_minus*exp(beta*negative_oligomer_charge*Donnan_potential)

       !a = integrate_z_cylindrical(n_plus*positive_bead_charge, "all_z")
       !b = integrate_z_cylindrical(n_minus*negative_bead_charge, "all_z")
       !c = surface_charge_density_left_wall + surface_charge_density_right_wall

       !print *, "**************************************"
       !print *, a, b, c, abs(a + b + c)
       !print *, "**************************************"
       
       ! call SolveWithNewtonsMethod(a,b,c,Donnan_potential)

       ! n_plus = n_plus*exp(beta*positive_oligomer_charge*Donnan_potential)
       ! n_neutral = n_neutral*(exp(beta*positive_oligomer_charge*Donnan_potential))
       ! n_minus = n_minus*exp(beta*negative_oligomer_charge*Donnan_potential)


    else
       !Use Newton's method to solve.
       a = integrate_z_cylindrical(n_plus*positive_bead_charge, "all_z")
       b = integrate_z_cylindrical(n_minus*negative_bead_charge, "all_z")
       c = surface_charge_density_left_wall + surface_charge_density_right_wall

       call SolveWithNewtonsMethod(a,b,c,Donnan_potential, abort_now)

       !print *, "should never get here"
       !call abort()

       n_plus = n_plus*exp(beta*positive_oligomer_charge*Donnan_potential)
       n_neutral = n_neutral*(exp(beta*positive_oligomer_charge*Donnan_potential)*(5.0_dp/19.0_dp) + exp(beta*negative_oligomer_charge*Donnan_potential)*(14.0_dp/19.0_dp))
       n_minus = n_minus*exp(beta*negative_oligomer_charge*Donnan_potential)
       
       
    end if

    n_plus = 5.0_dp * n_minus
    
    !(surface_charge_density_left_wall + surface_charge_density_right_wall)

    !print *, "positive integral 2= ", integrate_z_cylindrical(n_plus, "all_z")
    !print *, "neutral integral 2= ", integrate_z_cylindrical(n_neutral, "all_z")
    !print *, "minus integral 2= ", integrate_z_cylindrical(n_minus, "all_z")

  end subroutine ImposeChargeNeutrality
  
  subroutine CalculateDonnanPotential(n_plus, n_minus, Donnan_potential, abort_now)
    real(dp), dimension(:) :: n_plus
    real(dp), dimension(:) :: n_minus
    real(dp) :: Donnan_potential
    logical :: abort_now

    real(dp) :: a, b, c, d
    real(dp) :: y1, y2 !two solutions for two roots

    !call abort()
    
    !if( ((surface_charge_density_left_wall + surface_charge_density_right_wall)/electric_charge) < 1.0E-10_dp) then
       !The walls have equal and opposite charge.  By symmetry arguments the density profiles are already electroneutral.
       !=> Donnan potential is zero.  Therefore, there is nothing to do.
    if(abs((positive_bead_charge + negative_bead_charge)/electric_charge) < 1.0E-10_dp) then
       !We have the special case of the positive and negative bead charge being the same.
       !We can solve this analytically, as q_{-} = -q_{+}.

       a = integrate_z_cylindrical(n_plus*positive_bead_charge, "all_z")
       b = integrate_z_cylindrical(n_minus*negative_bead_charge, "all_z")
       c = surface_charge_density_left_wall + surface_charge_density_right_wall

       !print *, "a,b,c", a, b, c

       !Finding the Donnan potential to impose electroneutrality is now equivalent to solving
       !a*exp(x) + b*exp(-x) + c = 0; where x = exp(beta*positive_bead_charge*Psi_{Donnan}).

       !Letting y = exp(x), and using the quadratic fromula we get

       y1 = (-1.0_dp*c + sqrt(c**2 - 4.0_dp*a*b))/(2.0_dp*a)
       y2 = (-1.0_dp*c - sqrt(c**2 - 4.0_dp*a*b))/(2.0_dp*a)

       if(y1 < 0 .and. y2 < 0) then
          print *, "charge.f90: ImposeChargeNeutrality: "
          print *, "Both solutions for y1 and y2 are negative."
          print *, "Since this is supposed to be exp(x), this is not possible...aborting..."
          call abort()
       else if(y1 > 0 .and. y2 > 0) then
          print *, "charge.f90: ImposeChargeNeutrality: "
          print *, "Both solutions for y1 and y2 are positive."
          print *, "Don't know which one to pick...aborting..."
          call abort()
       else
          y1 = max(y1, y2)
       end if

       Donnan_potential = log(y1)/(beta*positive_oligomer_charge)
       !print *, y1, y2
       !print *, "CALCULATED DONANN = ", Donnan_potential
       !print *, "Donnan potential = ", exp(beta*positive_oligomer_charge*Donnan_potential)
       !call abort()
       !n_plus = n_plus*exp(beta*positive_oligomer_charge*Donnan_potential)
       !n_neutral = n_neutral*exp(beta*positive_oligomer_charge*Donnan_potential)
       !n_minus = n_minus*exp(beta*negative_oligomer_charge*Donnan_potential)

    else
       !Use Newton's method to solve.
       a = integrate_z_cylindrical(n_plus*positive_bead_charge, "all_z")
       b = integrate_z_cylindrical(n_minus*negative_bead_charge, "all_z")
       c = surface_charge_density_left_wall + surface_charge_density_right_wall

       !print *, "c = ", c
       call SolveWithNewtonsMethod(a,b,c,Donnan_potential, abort_now)
       !print *, Donnan_potential
       !call abort()

       !n_plus = n_plus*exp(beta*positive_oligomer_charge*Donnan_potential)
       !n_neutral = n_neutral*(exp(beta*positive_oligomer_charge*Donnan_potential)*(5.0_dp/19.0_dp) + exp(beta*negative_oligomer_charge*Donnan_potential)*(14.0_dp/19.0_dp))
       !n_minus = n_minus*exp(beta*negative_oligomer_charge*Donnan_potential)

    end if


    ! if( ((surface_charge_density_left_wall + surface_charge_density_right_wall)/electric_charge) < 1.0E-10_dp) then
    !    !The walls have equal and opposite charge.  By symmetry arguments the density profiles are already electroneutral.
    !    !=> Donnan potential is zero.  Therefore, there is nothing to do.
    ! else if(((positive_bead_charge + negative_bead_charge)/electric_charge) < 1.0E-10_dp) then
    !    !We have the special case of the positive and negative bead charge being the same.
    !    !We can solve this analytically, as q_{-} = -q_{+}.

    !    a = integrate_z_cylindrical(n_plus*positive_bead_charge, "all_z")
    !    b = integrate_z_cylindrical(n_minus*negative_bead_charge, "all_z")
    !    c = surface_charge_density_left_wall + surface_charge_density_right_wall

    !    !Finding the Donnan potential to impose electroneutrality is now equivalent to solving
    !    !a*exp(x) + b*exp(-x) + c = 0; where x = exp(beta*positive_bead_charge*Psi_{Donnan}).

    !    !Letting y = exp(x), and using the quadratic fromula we get

    !    y1 = (-1.0_dp*c + sqrt(c**2 - 4.0_dp*a*b))/(2.0_dp*a)
    !    y2 = (-1.0_dp*c - sqrt(c**2 - 4.0_dp*a*b))/(2.0_dp*a)

    !    if(y1 < 0 .and. y2 < 0) then
    !       print *, "charge.f90: ImposeChargeNeutrality: "
    !       print *, "Both solutions for y1 and y2 are negative."
    !       print *, "Since this is supposed to be exp(x), this is not possible...aborting..."
    !       call abort()
    !    else if(y1 > 0 .and. y2 > 0) then
    !       print *, "charge.f90: ImposeChargeNeutrality: "
    !       print *, "Both solutions for y1 and y2 are positive."
    !       print *, "Don't know which one to pick...aborting..."
    !       call abort()
    !    else
    !       y1 = max(y1, y2)
    !    end if

    !    Donnan_potential = log(y1)/(beta*positive_oligomer_charge)

    !    print *, "CALCULATED_DONNAN = ", Donnan_potential
    !    !print *, "Donnan potential = ", exp(beta*positive_oligomer_charge*Donnan_potential)

    !    !n_plus = n_plus*exp(beta*positive_oligomer_charge*Donnan_potential)
    !    !n_neutral = n_neutral*exp(beta*positive_oligomer_charge*Donnan_potential)
    !    !n_minus = n_minus*exp(beta*negative_oligomer_charge*Donnan_potential)

    ! else
    !    !Use Newton's method to solve.
    !    a = integrate_z_cylindrical(n_plus*positive_bead_charge, "all_z")
    !    b = integrate_z_cylindrical(n_minus*negative_bead_charge, "all_z")
    !    c = surface_charge_density_left_wall + surface_charge_density_right_wall

    !    call SolveWithNewtonsMethod(a,b,c,Donnan_potential)

    !    !print *, "Newton's Method not yet implemented"
    !    !call abort()
    ! end if

  end subroutine CalculateDonnanPotential

end module charge
