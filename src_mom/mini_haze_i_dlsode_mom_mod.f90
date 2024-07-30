module mini_haze_i_dlsode_mom_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomended if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128

  !! Constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: R_gas = 8.31446261815324e7_dp
  real(dp), parameter :: amu = 1.66053906892e-24_dp

  !! Global variables
  real(dp) :: T, nd_atm, rho, P_cgs, grav_cgs

  !! Variables needed to be sent in from outside
  character(len=20) :: prod_scheme
  real(dp) :: Prod_in ! Mass mixing ratio production rate of precurser molecules
  real(dp) :: r_mon ! Haze particle monomer size
  real(dp) :: rho_h ! bulk density of haze [g cm-3]
  real(dp) :: p_deep ! Deep removal pressure level
  real(dp) :: tau_loss ! Deep removal loss timescale
  real(dp) :: tau_decay, tau_act, tau_form

  real(dp) :: V_mon ! Haze particle monomer volume

  real(dp), parameter :: d_H2 = 2.827e-8_dp, LJ_H2 = 59.7_dp * kb, molg_H2 = 2.01588_dp
  real(dp), parameter :: d_He = 2.511e-8_dp, LJ_He = 10.22_dp * kb, molg_He = 4.002602_dp
  real(dp) :: mfp, eta

  public :: mini_haze_i_dlsode_mom, RHS_mom, jac_dum

  contains

  subroutine mini_haze_i_dlsode_mom(n_eq, T_in, P_in, mu_in, grav_in, t_end, q)
    implicit none

    ! Input variables
    integer, intent(in) :: n_eq
    real(dp), intent(in) :: T_in, P_in, mu_in, grav_in, t_end

    ! Input/Output tracer values
    real(dp), dimension(n_eq), intent(inout) :: q

    integer :: ncall

    ! Time controls
    real(dp) :: t_now

    ! DLSODE variables
    real(dp), allocatable, dimension(:) :: y
    real(dp), allocatable, dimension(:) :: rwork, rtol, atol
    integer, allocatable, dimension(:) :: iwork
    integer :: itol, itask, istate, iopt, mf
    integer :: rworkdim, iworkdim

    !! Find the number density of the atmosphere
    T = T_in
    P_cgs = P_in * 10.0_dp   ! Convert pascal to dyne cm-2

    !! Change gravity to cgs
    grav_cgs = grav_in * 100.0_dp

    !! Number density [cm-3] of layer
    nd_atm = P_cgs/(kb*T_in)  

    !! Mass density of layer
    rho = (P_cgs*mu_in*amu)/(kb * T_in) ! Mass density [g cm-3]

    !! Calculate dynamical viscosity for this layer
    eta = (5.0_dp/16.0_dp) * (sqrt(pi*(molg_H2*amu)*kb*T_in)/(pi*d_H2**2)) &
      & * ((((kb*T_in)/LJ_H2)**(0.16_dp))/1.22_dp)

    !! Calculate mean free path for this layer
    mfp = (2.0_dp*eta/rho) * sqrt((pi * mu_in)/(8.0_dp*R_gas*T_in))

    !! Haze particle monomer volume
    V_mon = 4.0_dp/3.0_dp * pi * r_mon**3


    ! -----------------------------------------
    ! ***  parameters for the DLSODE solver  ***
    ! -----------------------------------------

    itask = 1
    istate = 1
    iopt = 1

    ! Problem is stiff (usual)
    ! mf = 21 - full jacobian matrix with jacobian save
    ! mf = 22 - internal calculated jacobian
    mf = 22
    rworkdim = 22 + 9*n_eq + n_eq**2
    iworkdim = 20 + n_eq
    allocate(rtol(n_eq), atol(n_eq), rwork(rworkdim), iwork(iworkdim))

    itol = 4
    rtol(:) = 1.0e-3_dp           ! Relative tolerances for each scalar
    atol(:) = 1.0e-99_dp               ! Absolute tolerance for each scalar (floor value)

    rwork(:) = 0.0_dp
    iwork(:) = 0

    rwork(1) = 0.0_dp               ! Critical T value (don't integrate past time here)
    rwork(5) = 0.0_dp              ! Initial starting timestep (start low, will adapt in DVODE)
    rwork(6) = 0.0_dp       ! Maximum timestep

    iwork(5) = 0               ! Max order required
    iwork(6) = 100000               ! Max number of internal steps
    iwork(7) = 1                ! Number of error messages


    allocate(y(n_eq))

    !! Give tracer values to y - convert to number density
    y(1) = q(1)*nd_atm ! Convert to real number density
    y(2) = q(2)*rho   ! Convery to real mass
    y(3) = q(3)
    y(4) = q(4)

    t_now = 0.0_dp

    ! Set the printing flag
    ! 0 = no printing, 1 = printing
    call xsetf(1)

    ncall = 0

    do while (t_now < t_end)

      call DLSODE (RHS_mom, n_eq, y, t_now, t_end, itol, rtol, atol, itask, &
        & istate, iopt, rwork, rworkdim, iwork, iworkdim, jac_dum, mf)

      ncall = ncall + 1

      if (mod(ncall,50) == 0) then
        istate = 1
      else  if (istate == -1) then
        istate = 2
      else if (istate < -1) then
        print*, 'dlsode: ', istate
        exit
      end if

    end do

    !print*, t_now, y(:), ((3.0_dp*y(2)/y(1))/(4.0_dp*pi*rho_h))**(1.0_dp/3.0_dp) * 1e4_dp

    !! Give y values to tracers - convert to number mixing ratio + mass mixing ratio
    q(1) = y(1)/nd_atm
    q(2) = y(2)/rho
    q(3) = y(3)
    q(4) = y(4)

    deallocate(y, rtol, atol, rwork, iwork)

  end subroutine mini_haze_i_dlsode_mom

  subroutine RHS_mom(n_eq, time, y, f)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), intent(inout) :: time
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq), intent(inout) :: f

    real(dp) :: f_coal, f_coag
    real(dp) :: f_act, f_decay_pre, f_decay_act, f_form
    real(dp), dimension(n_eq) :: f_loss, f_prod

    !! In this routine, you calculate the new fluxes (f) for each moment
    !! The values of each bin (y) are typically kept constant
    !! Basically, you solve for the RHS of the ODE for each moment

    f(:) = 0.0_dp

    !! Find the RHS of the ODE for each particle bin size

    !! Calculate the coalesence loss rate for the zeroth moment
    call calc_coal(n_eq, y, f_coal)

    !! Calculate the coagulation loss rate for the zeroth moment
    call calc_coag(n_eq, y, f_coag)

    !! Add thermal decomposition loss term if above given pressure level (pa)
    if (P_cgs >= p_deep) then
      f_loss(:) = y(:)/tau_loss
    else
      f_loss(:) = 0.0_dp
    end if

    !! Calculate precursor and activated molecules rates
    f_act = y(3)/tau_act
    f_decay_pre = y(3)/tau_decay
    f_decay_act = y(4)/tau_decay
    f_form = y(4)/tau_form

    !! Add a production rate source term for each moment and precursor activation
    call find_production_rate(n_eq, f_prod, f_form)

    f(1) = f_prod(1) - f_coal - f_coag - f_loss(1)
    f(2) = f_prod(2) - f_loss(2)
    f(3) = f_prod(3) - f_act - f_decay_pre
    f(4) = f_prod(4) + f_act - f_form - f_decay_act

    !print*, 'y', time, y(:), ((3.0_dp*y(2)/y(1))/(4.0_dp*pi*rho_h))**(1.0_dp/3.0_dp) * 1e4_dp
    !print*, 'f', f(:), f_prod(:)
    !print*, 'f2', f_coal, f_coag, f_loss, f_act, f_decay_pre, f_decay_act, f_form
      

  end subroutine RHS_mom


  !! Routine called by the mini-haze integrator - calculates production rates from mass mixing ratio rate
  subroutine find_production_rate(n_eq, f_prod, f_form)
    implicit none

    integer, intent(in) :: n_eq

    real(dp), intent(in) :: f_form

    real(dp), dimension(n_eq), intent(out) :: f_prod

    f_prod(1) = f_form * (rho/(V_mon*rho_h))
    f_prod(2) = f_form * rho
    f_prod(3) = Prod_in
    f_prod(4) = 0.0_dp

  end subroutine find_production_rate

  subroutine calc_coag(n_eq, y, f_coag)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y 

    real(dp), intent(inout) :: f_coag

    real(dp) :: m_c, r, Kn, beta, regime1, regime2

    if (y(1) <= 1.0e-20_dp) then
      f_coag = 0.0_dp
      return
    end if

    !! Condensed mass
    m_c = y(2)/y(1)

    !! Mean radius
    r = ((3.0_dp*m_c)/(4.0_dp*pi*rho_h))**(1.0_dp/3.0_dp)

    regime1 = 8.0_dp * sqrt((pi*kb*T)/(m_c)) * r**2 * y(1)**2

    !! Knudsen number
    Kn = mfp/r

    !! Cunningham slip factor
    beta = 1.0_dp + Kn*(1.257_dp + 0.4_dp * exp(-1.1_dp/Kn))

    regime2 = (4.0_dp * kb * T * beta)/(3.0_dp * eta) * y(1)**2

    f_coag = min(regime1, regime2)

  end subroutine calc_coag

  subroutine calc_coal(n_eq, y, f_coal)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y 

    real(dp), intent(inout) :: f_coal

    real(dp) :: beta, Kn, r, m_c, d_vf, vf, StK, E
    real(dp), parameter :: eps = 0.5_dp

    if (y(1) <= 1.0e-20_dp) then
      f_coal = 0.0_dp
      return
    end if

    !! Condensed mass
    m_c = y(2)/y(1)

    !! Mean radius
    r = ((3.0_dp*m_c)/(4.0_dp*pi*rho_h))**(1.0_dp/3.0_dp)

    !! Knudsen number
    Kn = mfp/r

    !! Cunningham slip factor
    beta = 1.0_dp + Kn*(1.257_dp + 0.4_dp * exp(-1.1_dp/Kn))

    !! Settling velocity
    vf = (2.0_dp * beta * grav_cgs * r**2 * rho_h)/(9.0_dp * eta) & 
      & * (1.0_dp + & 
      & ((0.45_dp*grav_cgs*r**3*rho*rho_h)/(54.0_dp*eta**2))**(2.0_dp/5.0_dp))**(-5.0_dp/4.0_dp)

    !! Estimate differential velocity
    d_vf = eps * vf

    !! Calculate E
    if (Kn > 1.0_dp) then
      !! E = 1 when Kn > 1
      E = 1.0_dp
    else
      !! Calculate Stokes number
      StK = (vf * d_vf)/(grav_cgs * r)
      E = max(0.0_dp,1.0_dp - 0.42_dp*Stk**(-0.75_dp))
    end if

    !! Finally calculate the loss flux term
    f_coal = 2.0_dp*pi*r**2*y(1)**2*d_vf*E

  end subroutine calc_coal

  subroutine jac_dum (NEQ, X, Y, ML, MU, PD, NROWPD)
    integer, intent(in) :: NEQ, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(NEQ), intent(in) :: Y
    real(dp), dimension(NROWPD, NEQ), intent(inout) :: PD
  end subroutine jac_dum

end module mini_haze_i_dlsode_mom_mod