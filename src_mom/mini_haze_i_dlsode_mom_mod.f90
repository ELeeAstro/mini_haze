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
  real(dp) :: Prod_in ! Mass mixing ratio production rate of precurser molecules
  real(dp) :: r_mon ! Haze particle monomer size
  real(dp) :: rho_h ! bulk density of haze [g cm-3]
  real(dp) :: p_deep ! Deep removal pressure level
  real(dp) :: tau_loss ! Deep removal loss timescale
  real(dp) :: tau_decay, tau_act, tau_form ! Precurser timescales

  real(dp) :: V_mon ! Haze particle monomer volume

  real(dp) :: mfp, eta

  !! Diameter, LJ potential and molecular weight for background gases
  real(dp), parameter :: d_OH = 3.06e-8_dp, LJ_OH = 100.0_dp * kb, molg_OH = 17.00734_dp  ! estimate
  real(dp), parameter :: d_H2 = 2.827e-8_dp, LJ_H2 = 59.7_dp * kb, molg_H2 = 2.01588_dp
  real(dp), parameter :: d_H2O = 2.641e-8_dp, LJ_H2O = 809.1_dp * kb, molg_H2O = 18.01528_dp
  real(dp), parameter :: d_H = 2.5e-8_dp, LJ_H = 30.0_dp * kb, molg_H = 1.00794_dp
  real(dp), parameter :: d_CO = 3.690e-8_dp, LJ_CO = 91.7_dp * kb, molg_CO = 28.0101_dp
  real(dp), parameter :: d_CO2 = 3.941e-8_dp, LJ_CO2 = 195.2_dp * kb, molg_CO2 = 44.0095_dp
  real(dp), parameter :: d_O = 2.66e-8_dp, LJ_O = 70.0_dp * kb, molg_O = 15.99940_dp
  real(dp), parameter :: d_CH4 = 3.758e-8_dp, LJ_CH4 = 148.6_dp * kb, molg_CH4 = 16.0425_dp
  real(dp), parameter :: d_C2H2 = 4.033e-8_dp, LJ_C2H2 = 231.8_dp * kb, molg_C2H2 = 26.0373_dp
  real(dp), parameter :: d_NH3 = 2.900e-8_dp, LJ_NH3 = 558.3_dp * kb, molg_NH3 = 17.03052_dp
  real(dp), parameter :: d_N2 = 3.798e-8_dp, LJ_N2 = 71.4_dp * kb, molg_N2 = 14.0067_dp
  real(dp), parameter :: d_HCN = 3.630e-8_dp, LJ_HCN = 569.1_dp * kb, molg_HCN = 27.0253_dp
  real(dp), parameter :: d_He = 2.511e-8_dp, LJ_He = 10.22_dp * kb, molg_He = 4.002602_dp

  !! Constuct required arrays for calculating gas mixtures
  real(dp), dimension(3) :: d_g = (/d_H2, d_He, d_H/)
  real(dp), dimension(3) :: d_LJ = (/LJ_H2, LJ_He, LJ_H/)
  real(dp), dimension(3) :: d_molg = (/molg_H2, molg_He, molg_H/)

  public :: mini_haze_i_dlsode_mom, RHS_mom, jac_dum

  contains

  subroutine mini_haze_i_dlsode_mom(n_eq, T_in, P_in, mu_in, grav_in, t_end, q, n_gas, g_VMR)
    implicit none

    ! Input variables
    integer, intent(in) :: n_eq, n_gas
    real(dp), intent(in) :: T_in, P_in, mu_in, grav_in, t_end

    ! Input/Output tracer values
    real(dp), dimension(n_eq), intent(inout) :: q
    real(dp), dimension(n_gas), intent(in) :: g_VMR

    integer :: ncall, n

    ! Time controls
    real(dp) :: t_now

    ! DLSODE variables
    real(dp), allocatable, dimension(:) :: y
    real(dp), allocatable, dimension(:) :: rwork, rtol, atol
    integer, allocatable, dimension(:) :: iwork
    integer :: itol, itask, istate, iopt, mf
    integer :: rworkdim, iworkdim

    !! Work variables
    real(dp), dimension(n_gas) :: g_eta
    real(dp) :: bot, top

    !! Find the number density of the atmosphere
    T = T_in
    P_cgs = P_in * 10.0_dp   ! Convert pascal to dyne cm-2

    !! Change gravity to cgs
    grav_cgs = grav_in * 100.0_dp

    !! Number density [cm-3] of layer
    nd_atm = P_cgs/(kb*T_in)  

    !! Mass density of layer
    rho = (P_cgs*mu_in*amu)/(kb * T_in) ! Mass density [g cm-3]

    !! Calculate dynamical viscosity for this layer - do square root mixing law from Rosner 2012
    do n = 1, n_gas
      g_eta(n) = (5.0_dp/16.0_dp) * (sqrt(pi*(d_molg(n)*amu)*kb*T_in)/(pi*d_g(n)**2)) &
        & * ((((kb*T_in)/d_LJ(n))**(0.16_dp))/1.22_dp)
    end do

    !! Mass square root mixing law
    top = 0.0_dp
    bot = 0.0_dp
    do n = 1, n_gas
      top = top + sqrt(d_molg(n)*amu)*g_VMR(n)*g_eta(n)
      bot = bot + sqrt(d_molg(n)*amu)*g_VMR(n)
    end do

    !! Mixture dynamical viscosity
    eta = top/bot

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

    !! Give tracer values to y
    y(:) = q(:)

    t_now = 0.0_dp

    ! Set the printing flag
    ! 0 = no printing, 1 = printing
    call xsetf(1)

    ncall = 0

    do while ((t_now < t_end) .and. (ncall < 100))

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

    !! Give y values to tracers
    q(:) = y(:)

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

    !! Convert y to real numbers to calculate f
    y(1) = y(1)*nd_atm ! Convert to real number density
    y(2) = y(2)*rho   ! Convert to real mass density

    !! Find the RHS of the ODE for each particle bin size
    f(:) = 0.0_dp

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

    !! Calculate final net flux rate for each tracer
    f(1) = f_prod(1) - f_coal - f_coag - f_loss(1)
    f(2) = f_prod(2) - f_loss(2)
    f(3) = f_prod(3) - f_act - f_decay_pre
    f(4) = f_prod(4) + f_act - f_form - f_decay_act

    !! Convert f to ratios
    f(1) = f(1)/nd_atm
    f(2) = f(2)/rho
      
    !! Convert y back to ratios
    y(1) = y(1)/nd_atm
    y(2) = y(2)/rho 

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

    !! Mean mass of particle
    m_c = y(2)/y(1)

    !! Mass weighted mean radius of particle
    r = ((3.0_dp*m_c)/(4.0_dp*pi*rho_h))**(1.0_dp/3.0_dp)

    regime1 = 8.0_dp * sqrt((pi*kb*T)/m_c) * r**2 * y(1)**2

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

    !! Mean mass of particle
    m_c = y(2)/y(1)

    !! Mass weighted mean radius of particle
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
    if (Kn >= 1.0_dp) then
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