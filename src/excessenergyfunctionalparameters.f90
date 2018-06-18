!The notation used in this module follows that of the paper
!"Evaluating the accuracy of a density functional theory of polymer
!solutions with additive hard sphere diameters", by Forsman and Woodward
!Journal of Chemical Physics, volume 120, number 1, Jan 2004.
module excessenergyfunctionalparameters
  use kinds
  use universalconstants
  use parameters
  implicit none
  private

  public :: GetAEx
  public :: GetAExDerivIntegrand

  real(dp) :: sigma_monomer
  real(dp) :: sigma_solvent

contains

  subroutine GetAExTermWeightings(n_mbar, n_sbar, term_index, Y_special, Psi, X, W, Z)
    real(dp), dimension(:), intent(in)              :: n_mbar
    real(dp), dimension(:), intent(in)              :: n_sbar
    integer, intent(in)                             :: term_index !0,1,2 as discussed in J.Chem.Phys. vol 120, #1, 1994.
    real(dp), dimension(size(n_mbar)), intent(out)  :: Y_special
    real(dp), dimension(size(n_mbar)), intent(out)  :: Psi
    real(dp), dimension(size(n_mbar)), intent(out)  :: W
    real(dp), dimension(size(n_mbar)), intent(out)  :: X
    real(dp), dimension(size(n_mbar)), intent(out)  :: Z

    real(dp) :: q

    real(dp) :: r_0, r_1, r_2
    real(dp) :: s_0, s_1, s_2
    real(dp) :: b_0, b_1, b_2

    real(dp), dimension(size(n_mbar)) :: x_d_hat
    real(dp), dimension(size(n_mbar)) :: x_m_hat

    real(dp), dimension(size(n_mbar)) :: rhat_2, rhat_1, rhat_0
    real(dp), dimension(size(n_mbar)) :: shat_2, shat_1, shat_0 
    real(dp), dimension(size(n_mbar)) :: bhat_2, bhat_1, bhat_0

    real(dp), dimension(size(n_mbar)) :: c_0, c_1, c_2

    real(dp), dimension(size(n_mbar)) :: omega_0, omega_1, omega_2
    real(dp), dimension(size(n_mbar)) :: phi_0, phi_1, phi_2
    real(dp), dimension(size(n_mbar)) :: omega_prime_0, omega_prime_1, omega_prime_2 
    real(dp), dimension(size(n_mbar)) :: W_0, W_1, W_2
    real(dp), dimension(size(n_mbar)) :: X_0, X_1, X_2
    real(dp), dimension(size(n_mbar)) :: Y_0, Y_1, Y_2
    real(dp), dimension(size(n_mbar)) :: Z_0, Z_1, Z_2
    real(dp), dimension(size(n_mbar)) :: Y_special_0, Y_special_1, Y_special_2
    real(dp), dimension(size(n_mbar)) :: Psi_0, Psi_1, Psi_2

    real(dp), dimension(size(n_mbar)) :: eta

    if(size(n_mbar) /= size(n_sbar)) then
       print *, "excessenergyfunctionalparameters.f90: GetAExTermWeightings"
       print *, "size(n_mbar) /= size(n_sbar)"
       print *, "input variable size mismatch...aborting..."
       call abort()
    end if

    call InitialiseHardSphereDiameters(sigma_monomer, sigma_solvent)
    
    q = sigma_solvent / sigma_monomer

    eta(:) = pi * (n_mbar(:) + n_sbar(:)) * (hs_diameter**3) / 6.0_dp

    r_2 = 0.75_dp * sigma_monomer
    r_1 = 0.5_dp * sigma_monomer
    r_0 = 0.5_dp * sigma_monomer * q

    s_2 = 2.0_dp * pi * (sigma_monomer**2)
    s_1 = pi * (sigma_monomer**2)
    s_0 = pi * (sigma_monomer**2) * (q**2)

    b_2 = (pi/3.0_dp) * (sigma_monomer**3)
    b_1 = (pi/6.0_dp) * (sigma_monomer**3)
    b_0 = (pi/6.0_dp) * (sigma_monomer**3) * (q**3)

    x_d_hat(:) = (0.5_dp * n_mbar(:)) / ((0.5_dp * n_mbar(:)) + n_sbar(:))
    x_m_hat(:) = n_mbar(:) / (n_mbar(:) + n_sbar(:))

    rhat_2(:) = (x_d_hat(:) * r_2) + ((1.0_dp - x_d_hat(:)) * r_0)
    rhat_1(:) = (x_m_hat(:) * r_1) + ((1.0_dp - x_m_hat(:)) * r_0)
    rhat_0(:) = rhat_2(:)

    shat_2(:) = (x_d_hat(:) * s_2) + ((1.0_dp - x_d_hat(:)) * s_0)
    shat_1(:) = (x_m_hat(:) * s_1) + ((1.0_dp - x_m_hat(:)) * s_0)
    shat_0(:) = shat_2(:)

    bhat_2(:) = (x_d_hat(:) * b_2) + ((1.0_dp - x_d_hat(:)) * b_0)
    bhat_1(:) = (x_m_hat(:) * b_1) + ((1.0_dp - x_m_hat(:)) * b_0)
    bhat_0(:) = bhat_2(:)

    c_2(:) = (x_d_hat(:) * (r_2**2)) + ((1.0_dp - x_d_hat(:)) * (r_0**2))
    c_1(:) = (x_m_hat(:) * (r_1**2)) + ((1.0_dp - x_m_hat(:)) * (r_0**2))
    c_0(:) = c_2(:)

    omega_2(:) = (rhat_2(:) * shat_2(:)) / bhat_2(:)
    omega_1(:) = (rhat_1(:) * shat_1(:)) / bhat_1(:)
    omega_0(:) = (rhat_0(:) * shat_0(:)) / bhat_0(:)

    phi_2(:) = (c_2(:) * (shat_2(:)**2)) / (9.0_dp * (bhat_2(:)**2))
    phi_1(:) = (c_1(:) * (shat_1(:)**2)) / (9.0_dp * (bhat_1(:)**2))
    phi_0(:) = (c_0(:) * (shat_0(:)**2)) / (9.0_dp * (bhat_0(:)**2))

    omega_prime_2(:) = ((r_2 * shat_2(:)) + (rhat_2(:) * shat_2(:)) - (b_2 * (1.0_dp / bhat_2(:)) * rhat_2(:) * shat_2(:))) &
         / (bhat_2(:))
    omega_prime_1(:) = ((r_1 * shat_1(:)) + (rhat_1(:) * shat_1(:)) - (b_1 * (1.0_dp / bhat_1(:)) * rhat_1(:) * shat_1(:))) &
         / (bhat_1(:))
    omega_prime_0(:) = ((r_0 * shat_0(:)) + (rhat_0(:) * shat_0(:)) - (b_0 * (1.0_dp / bhat_0(:)) * rhat_0(:) * shat_0(:))) &
         / (bhat_0(:))

    W_2(:) = ((2.0_dp * b_2 * (1.0_dp / bhat_2(:)) * c_2(:) * (shat_2(:)**2)) - (2.0_dp * c_2 * shat_2(:) * s_2) - &
         ((rhat_2(:)**2) * (shat_2(:)**2))) / (9.0_dp * (bhat_2(:)**2))
    W_1(:) = ((2.0_dp * b_1 * (1.0_dp / bhat_1(:)) * c_1(:) * (shat_1(:)**2)) - (2.0_dp * c_1 * shat_1(:) * s_1) - &
         ((rhat_1(:)**2) * (shat_1(:)**2))) / (9.0_dp * (bhat_1(:)**2))
    W_0(:) = ((2.0_dp * b_0 * (1.0_dp / bhat_0(:)) * c_0(:) * (shat_0(:)**2)) - (2.0_dp * c_0 * shat_0(:) * s_0) - &
         ((rhat_0(:)**2) * (shat_0(:)**2))) / (9.0_dp * (bhat_0(:)**2))

    X_2(:) = W_2(:) - omega_prime_2(:) - ( (omega_2(:) + 1.0_dp) * (b_2) / (bhat_2(:)) )
    X_1(:) = W_1(:) - omega_prime_1(:) - ( (omega_1(:) + 1.0_dp) * (b_1) / (bhat_1(:)) )
    X_0(:) = W_0(:) - omega_prime_0(:) - ( (omega_0(:) + 1.0_dp) * (b_0) / (bhat_0(:)) )

    Y_2(:) = (-1.0_dp * W_2(:)) + (2.0_dp * omega_prime_2(:)) - (3.0_dp * (phi_2(:) - omega_2(:) - 2.0_dp) * b_2 / bhat_2(:))
    Y_1(:) = (-1.0_dp * W_1(:)) + (2.0_dp * omega_prime_1(:)) - (3.0_dp * (phi_1(:) - omega_1(:) - 2.0_dp) * b_1 / bhat_1(:))
    Y_0(:) = (-1.0_dp * W_0(:)) + (2.0_dp * omega_prime_0(:)) - (3.0_dp * (phi_0(:) - omega_0(:) - 2.0_dp) * b_0 / bhat_0(:))

    Z_2(:) = ( (phi_2(:) - 1.0_dp) * b_2 / bhat_2(:) ) - omega_prime_2(:)
    Z_1(:) = ( (phi_1(:) - 1.0_dp) * b_1 / bhat_1(:) ) - omega_prime_1(:)
    Z_0(:) = ( (phi_0(:) - 1.0_dp) * b_0 / bhat_0(:) ) - omega_prime_0(:)

    Y_special_2(:) = W_2(:) + Y_2(:) + ( 3.0_dp * Z_2(:) )
    Y_special_1(:) = W_1(:) + Y_1(:) + ( 3.0_dp * Z_1(:) )
    Y_special_0(:) = W_0(:) + Y_0(:) + ( 3.0_dp * Z_0(:) )

    Psi_2(:) = ( 0.5_dp * X_2(:) ) + ( 2.5_dp * Y_2(:) ) + ( 5.5_dp * Z_2(:) ) + ( 3.0_dp * W_2(:) )
    Psi_1(:) = ( 0.5_dp * X_1(:) ) + ( 2.5_dp * Y_1(:) ) + ( 5.5_dp * Z_1(:) ) + ( 3.0_dp * W_1(:) )
    Psi_0(:) = ( 0.5_dp * X_0(:) ) + ( 2.5_dp * Y_0(:) ) + ( 5.5_dp * Z_0(:) ) + ( 3.0_dp * W_0(:) )  

    !For ease of readability and comparison with the paper we first calculate all the values of the
    !term weightings and then return only the ones that we have asked for.  If calculation speed becomes
    !a problem we can change this.
    if(term_index == 0) then

       Y_special(:) = Y_special_0(:)
       Psi(:) = Psi_0(:)
       Z(:) = Z_0(:)
       X(:) = X_0(:)
       W(:) = W_0(:)

    else if(term_index == 1) then

       Y_special(:) = Y_special_1(:)
       Psi(:) = Psi_1(:)
       Z(:) = Z_1(:)
       X(:) = X_1(:)
       W(:) = W_1(:)

    else if(term_index == 2) then

       Y_special(:) = Y_special_2(:)
       Psi(:) = Psi_2(:)
       Z(:) = Z_2(:)
       X(:) = X_2(:)
       W(:) = W_2(:)

    else
       print *, "excessenergyfunctionalparameters.f90: GetAExTermWeightings"
       print *, "Value of term index = ", term_index
       print *, "Only supports values of 0, 1 or 2 as per the paper"
       print *, "Coding error...aborting..."
    end if

    ! print *, "Y_special = ", Y_special
    ! print *, "Psi = ", Psi
    ! print *, "Z = ", Z
    ! print *, "X = ", X
    ! print *, "W = ", W
    ! call abort()

  end subroutine GetAExTermWeightings

  function GetAEx(n_mbar, n_sbar, term_index)
    real(dp), dimension(:), intent(in)  :: n_mbar
    real(dp), dimension(:), intent(in)  :: n_sbar
    integer, intent(in)                 :: term_index
    real(dp), dimension(size(n_mbar)) :: GetAEx

    real(dp), dimension(size(n_mbar)) :: Y_special, Psi, X, W, Z
    real(dp), dimension(size(n_mbar)) :: eta_bar

    if(size(n_mbar) /= size(n_sbar)) then
       print *, "excessenergyfunctionalparameters.f90: GetAEx"
       print *, "size(n_mbar) /= size(n_sbar)"
       print *, "input variable size mismatch...aborting..."
       call abort()
    end if

    call GetAExTermWeightings(n_mbar, n_sbar, term_index, Y_special, Psi, X, W, Z)

    eta_bar(:) = (pi * (n_mbar(:) + n_sbar(:)) * hs_diameter**3) / 6.0_dp

    GetAEx(:) = (( (1.0_dp - eta_bar(:)) / eta_bar(:) ) * log(1.0_dp - eta_bar(:)) ) + 1.0_dp + &
         ( Y_special(:) * ( (log(1.0_dp - eta_bar(:)) * ( ((1.0_dp - eta_bar(:))/(eta_bar(:))) + 1.0_dp )) - (1.0_dp/(1.0_dp - eta_bar(:))) -&
         (1.0_dp/(2.0_dp * ((1.0_dp - eta_bar(:))**2))) + (5.0_dp/2.0_dp) ) ) + &
         ( (Psi(:) + (2.0_dp * Z(:)) - X(:)) * (((1.0_dp)/(2.0_dp * ((1.0_dp - eta_bar(:))**2))) - 0.5_dp) ) + &
         (W(:) * ((1.0_dp/(2.0_dp * ((1.0_dp - eta_bar(:))**2))) - (2.0_dp/(1.0_dp - eta_bar(:))) - log(1.0_dp - eta_bar(:)) + 1.5_dp )) +&
         (Psi(:) * ((1.0_dp / (1.0_dp - eta_bar(:))) - (1.0_dp / (2.0_dp * ((1.0_dp - eta_bar(:))**2))) - 0.5_dp) )

  end function GetAEx

  function GetAExDerivIntegrand(n_mbar, n_sbar, term_index)
    real(dp), dimension(:), intent(in) :: n_mbar
    real(dp), dimension(:), intent(in) :: n_sbar
    integer, intent(in)                :: term_index
    real(dp), dimension(size(n_mbar)) :: GetAExDerivIntegrand

    real(dp), dimension(size(n_mbar)) :: eta_bar

    real(dp), dimension(size(n_mbar)) :: Y_special, Psi, X, W, Z

    if(size(n_mbar) /= size(n_sbar)) then
       print *, "excessenergyfunctionalparameters.f90: GetAEx"
       print *, "size(n_mbar) /= size(n_sbar)"
       print *, "input variable size mismatch...aborting..."
       call abort()
    end if

    call GetAExTermWeightings(n_mbar, n_sbar, term_index, Y_special, Psi, X, W, Z)

    eta_bar(:) = (pi * (n_mbar(:) + n_sbar(:)) * hs_diameter**3) / 6.0_dp

    GetAExDerivIntegrand(:) = ( (-1.0_dp / 8.0_dp) * ((log(1.0_dp - eta_bar(:))/(eta_bar(:)**2)) + (1.0_dp/eta_bar(:))) ) + & !"T1"
         ( ((-Y_special(:))/8.0_dp) * ((1.0_dp/eta_bar(:)) + (1.0_dp/(1.0_dp - eta_bar(:))) + (1.0_dp/((1.0_dp - eta_bar(:))**2)) +&
         (1.0_dp/((1.0_dp - eta_bar(:))**3)) + (log(1.0_dp - eta_bar(:))/eta_bar(:)**2)) ) + & !"T2"
         (((Psi(:) + (2.0_dp * Z(:)) - X(:))/8.0_dp) * (1.0_dp / ((1.0_dp - eta_bar(:))**3))) +& !"T3"
         ((W(:)/8.0_dp) * ((1.0_dp/((1.0_dp - eta_bar(:))**3)) - (2.0_dp/((1.0_dp - eta_bar(:))**2)) + &
         (1.0_dp/(1.0_dp - eta_bar(:))))) +& !"T4"
         ((Psi(:)/8.0_dp) * ((1.0_dp/((1.0_dp - eta_bar(:))**2)) - (1.0_dp/((1.0_dp - eta_bar(:))**3)))) !"T5"

  end function GetAExDerivIntegrand

  ! function GetYMix(n_mbar, n_sbar)
  !   real(dp), dimension(:), intent(in) :: n_mbar
  !   real(dp), dimension(:), intent(in) :: n_sbar
  !   real(dp), dimension(size(n_mbar)) :: GetYMix

  !   real(dp), dimension(size(n_mbar)) :: phi_m

  !   real(dp) :: q
  !   call InitialiseHardSphereDiameters(sigma_monomer, sigma_solvent)
  !   q = sigma_solvent / sigma_monomer

  !   phi_m(:) = n_mbar / (n_mbar + n_sbar*q)

    
  ! end function GetYMix

  subroutine InitialiseHardSphereDiameters(sigma_monomer, sigma_solvent)
    real(dp), intent(out) :: sigma_monomer
    real(dp), intent(out) :: sigma_solvent
    
    sigma_monomer = hs_diameter
    sigma_solvent = hs_diameter
    
  end subroutine InitialiseHardSphereDiameters
  
end module excessenergyfunctionalparameters
