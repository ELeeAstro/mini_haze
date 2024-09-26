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
  real(dp), parameter :: third = 1.0_dp/3.0_dp

  !! Global variables
  real(dp) :: T, mu, nd_atm, rho, p, grav

  !! Variables needed to be sent in from outside
  real(dp) :: Prod_in ! Mass mixing ratio production rate of precurser molecules
  real(dp) :: r_mon ! Haze particle monomer size
  real(dp) :: rho_d ! bulk density of haze [g cm-3]
  real(dp) :: p_deep ! Deep removal pressure level
  real(dp) :: tau_loss ! Deep removal loss timescale
  real(dp) :: tau_decay, tau_act, tau_form ! Precurser timescales

  real(dp) :: V_mon ! Haze particle monomer volume
  real(dp) :: m_mon ! Haze particle monomer mass

  real(dp) :: mfp, eta, nu, eps_d

  !! Kolmogorov scale
  real(dp), parameter :: l_k = 1.0_dp

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
  real(dp), allocatable, dimension(:) :: d_g, LJ_g, molg_g

  public :: mini_haze_i_dlsode_mom, RHS_mom, jac_dum
  private :: calc_coal, calc_coag, calc_turb_acc, calc_turb_shear, calc_turb_couple, &
    & find_production_rate, eta_construct

  contains

  subroutine mini_haze_i_dlsode_mom(n_eq, T_in, P_in, mu_in, grav_in, bg_VMR_in, t_end, sp_bg, q)
    implicit none

    ! Input variables
    integer, intent(in) :: n_eq
    character(len=20), dimension(:), intent(in) :: sp_bg
    real(dp), intent(in) :: T_in, P_in, mu_in, grav_in, t_end
    real(dp), dimension(:), intent(in) :: bg_VMR_in

    ! Input/Output tracer values
    real(dp), dimension(n_eq), intent(inout) :: q

    integer :: ncall, n

    ! Time controls
    real(dp) :: t_now

    ! DLSODE variables
    real(dp), allocatable, dimension(:) :: y
    real(dp), allocatable, dimension(:) :: rwork
    integer, allocatable, dimension(:) :: iwork
    integer :: itol, itask, istate, iopt, mf
    integer :: rworkdim, iworkdim
    real(dp) :: rtol, atol

    !! Work variables
    integer :: n_gas
    real(dp) :: eta_g, bot, top
    real(dp), allocatable, dimension(:) :: VMR_g

    !! Alter input values to mini-haze units
    !! (note, some are obvious not not changed in case specific models need different conversion factors)

    !! Find the number density of the atmosphere
    T = T_in             ! Convert temperature to K
    p = P_in * 10.0_dp   ! Convert pascal to dyne cm-2

    !! Background gas array
    n_gas = size(bg_VMR_in)
    allocate(VMR_g(n_gas))
    VMR_g(:) = bg_VMR_in(:)

    !! Change mu_in to mu
    mu = mu_in ! Convert mean molecular weight to mu [g mol-1]

    !! Change gravity to cgs [cm s-2]
    grav = grav_in * 100.0_dp

    !! Number density [cm-3] of layer
    nd_atm = p/(kb*T)  

    !! Mass density of layer
    rho = (p*mu*amu)/(kb * T) ! Mass density [g cm-3]

    !! Calculate dynamical viscosity for this layer - do square root mixing law from Rosner 2012
    call eta_construct(sp_bg)
    top = 0.0_dp
    bot = 0.0_dp
    do n = 1, n_gas
      eta_g = (5.0_dp/16.0_dp) * (sqrt(pi*(molg_g(n)*amu)*kb*T)/(pi*d_g(n)**2)) &
        & * ((((kb*T)/LJ_g(n))**(0.16_dp))/1.22_dp)
      top = top + sqrt(molg_g(n))*VMR_g(n)*eta_g
      bot = bot + sqrt(molg_g(n))*VMR_g(n)
    end do

    !! Mixture dynamical viscosity
    eta = top/bot

    !! Mixture kinematic viscosity
    nu = eta/rho

    !! dissipation of turbulent kinetic energy
    eps_d = nu**3/l_k**4

    !! Calculate mean free path for this layer
    mfp = (2.0_dp*eta/rho) * sqrt((pi * mu)/(8.0_dp*R_gas*T))

    !! Haze particle monomer volume
    V_mon = 4.0_dp/3.0_dp * pi * r_mon**3

    !! Haze particle momomer mass
    m_mon = rho_d * V_Mon


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
    allocate(rwork(rworkdim), iwork(iworkdim))

    itol = 1
    rtol = 1.0e-3_dp           ! Relative tolerances for each scalar
    atol = 1.0e-30_dp               ! Absolute tolerance for each scalar (floor value)

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
    y(:) = max(q(:),1e-30_dp)

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

    !! Give y values to tracers
    q(:) = max(y(:),1e-30_dp)

    deallocate(y, rwork, iwork, d_g, LJ_g, molg_g)

  end subroutine mini_haze_i_dlsode_mom

  subroutine RHS_mom(n_eq, time, y, f)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), intent(inout) :: time
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq), intent(inout) :: f

    real(dp) :: w_sh2, w_acc2, w_co2
    real(dp) :: f_coal, f_coag, f_turb
    real(dp) :: f_act, f_decay_pre, f_decay_act, f_form
    real(dp), dimension(2) :: f_loss
    real(dp), dimension(n_eq) ::  f_prod
    real(dp) :: m_h, r_h, Kn, beta, vf

    !! In this routine, you calculate the new fluxes (f) for each moment
    !! The values of each moment (y) are typically kept constant
    !! Basically, you solve for the RHS of the ODE for each moment

    !! Convert y to real numbers to calculate f
    y(1) = y(1)*nd_atm ! Convert to real number density
    y(2) = y(2)*rho   ! Convert to real mass density

    !! Find the RHS of the ODE for each moment

    !! Calculate coagulation and coalesence fluxes

    !! Mean mass of particle
    m_h = max(y(2)/y(1), m_mon)

    !! Mass weighted mean radius of particle
    r_h = max(((3.0_dp*m_h)/(4.0_dp*pi*rho_d))**(third), r_mon)

    !! Knudsen number
    Kn = mfp/r_h

    !! Cunningham slip factor
    beta = 1.0_dp + Kn*(1.257_dp + 0.4_dp * exp(-1.1_dp/Kn))

    !! Settling velocity
    vf = (2.0_dp * beta * grav * r_h**2 * rho_d)/(9.0_dp * eta) & 
      & * (1.0_dp &
      & + ((0.45_dp*grav*r_h**3*rho*rho_d)/(54.0_dp*eta**2))**(0.4_dp))**(-1.25_dp)

    !! Calculate the coagulation loss rate for the zeroth moment
    call calc_coag(n_eq, y, m_h, r_h, beta, f_coag)
    f_coag = max(1e-30_dp,f_coag)

    !! Calculate the coalesence loss rate for the zeroth moment
    call calc_coal(n_eq, y, r_h, Kn, vf, f_coal)
    f_coal = max(1e-30_dp,f_coal)

    !! Calculate the turbulent shear collision velocity
    call calc_turb_shear(n_eq, y, r_h, w_sh2)

    !! Calculate the turbulent acceleration collision velocity
    call calc_turb_acc(n_eq, y, vf, w_acc2)

    !! Calculate the turbulent fluid coupling collision velocity
    call calc_turb_couple(n_eq, y, r_h, vf, w_co2)

    !! Combine turbulent velocities into collision rate using total kernel
    f_turb = 4.0_dp * pi * r_h**2 * sqrt(2.0_dp/pi) * sqrt(w_sh2 + w_acc2 + w_co2) * y(1)**2
    f_turb = max(1e-30_dp,f_turb)

    !! Add thermal decomposition loss term if above given pressure level (pa)
    if (p > p_deep) then
      f_loss(1) = y(1)/tau_loss
      f_loss(2) = y(2)/tau_loss
    else
      f_loss(1) = 0.0_dp
      f_loss(2) = 0.0_dp
    end if
    f_loss = max(1e-30_dp,f_loss)

    !! Calculate precursor and activated molecules rates
    f_act = y(3)/tau_act
    f_decay_pre = y(3)/tau_decay
    f_decay_act = y(4)/tau_decay
    f_form = y(4)/tau_form

    !! Add a production rate source term for each moment and precursor activation
    call find_production_rate(n_eq, f_prod, f_form)

    !! Calculate final net flux rate for each tracer
    f(1) = f_prod(1) - f_coal - f_coag - f_loss(1) - f_turb
    f(2) = f_prod(2) - f_loss(2)
    f(3) = f_prod(3) - f_act - f_decay_pre
    f(4) = f_prod(4) + f_act - f_form - f_decay_act

    !! Convert f to ratios
    f(1) = f(1)/nd_atm
    f(2) = f(2)/rho
      
    !! Convert y back to ratios
    y(1) = y(1)/nd_atm
    y(2) = y(2)/rho 

    !print*, 'y', time, y(:), ((3.0_dp*y(2)/y(1))/(4.0_dp*pi*rho_d))**(1.0_dp/3.0_dp) * 1e4_dp
    !print*, 'f', f(:), f_prod(:)
    !print*, 'f2', f_coal, f_coag, f_loss, f_act, f_decay_pre, f_decay_act, f_form

  end subroutine RHS_mom

  !! Particle-particle turbulent shear collisions
  subroutine calc_turb_shear(n_eq, y, r_h, w_sh2)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y 
    real(dp), intent(in) :: r_h

    real(dp), intent(out) :: w_sh2

    !! Shear relative velocity (Wang et al. 1998)
    w_sh2 = 1.0_dp/15.0_dp * (2.0_dp*r_h)**2 * eps_d/nu

  end subroutine calc_turb_shear

  !! Particle-particle turbulent acceleration collisions
  subroutine calc_turb_acc(n_eq, y, vf, w_acc2)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y 
    real(dp), intent(in) :: vf

    real(dp), intent(out) :: w_acc2

    real(dp), parameter :: gam = 10.0_dp
    real(dp) :: u2, tau_i, b, T_l, th_i, c1, c2

    !! Follows Wang et al. (2001) 
    !! - see also Park et el. (2002) and Kruis & Kusters (1997)

    u2 = (gam * sqrt(eps_d*nu))/0.183_dp
    b = (3.0_dp*rho)/(2.0_dp*rho_d + rho)
    T_L = (0.4_dp*u2)/eps_d

    tau_i = vf/grav
    th_i = tau_i/T_L

    c1 = sqrt((1.0_dp + th_i + th_i)/((1.0_dp + th_i)*(1.0_dp + th_i)))
    c2 = (1.0_dp/(((1.0_dp + th_i)*(1.0_dp + th_i))) - 1.0_dp/(((1.0_dp + gam*th_i)*(1.0_dp + gam*th_i))))

    w_acc2 = 3.0_dp*(1.0_dp-b)**2*u2*(gam/(gam-1.0_dp)) & 
      & * (((th_i + th_i)**2 - 4.0_dp*th_i*th_i*c1)/(th_i + th_i)) * c2

  end subroutine calc_turb_acc

  !! Particle-particle turbulent fluid coupling collisions
  subroutine calc_turb_couple(n_eq, y, r_h, vf, w_co2)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y 
    real(dp), intent(in) :: r_h, vf

    real(dp), intent(out) :: w_co2

    real(dp), parameter :: lam_d = 10.0_dp
    real(dp) :: dudt2, tau

    !! Fluid coupling turbulent term (Wang et al. 1998)
    tau = vf/grav
    dudt2 = 1.16_dp * eps_d**(1.5_dp)/sqrt(nu)
    w_co2 = 2.0_dp*(1.0_dp - rho/rho_d)**2 * tau**2 * dudt2 * ((2.0_dp*r_h)**2/lam_d**2)

  end subroutine calc_turb_couple

  !! Particle-particle Brownian coagulation
  subroutine calc_coag(n_eq, y, m_h, r_h, beta, f_coag)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y 
    real(dp), intent(in) :: m_h, r_h, beta

    real(dp), intent(inout) :: f_coag

    real(dp) :: phi, del_r, D_r, V_r, lam_r

    !! Particle diffusion rate
    D_r = (kb*T*beta)/(6.0_dp*pi*eta*r_h)

    !! Thermal velocity limit rate
    V_r = sqrt((8.0_dp*kb*T)/(pi*m_h))

    !! Ratio fraction
    lam_r = (8.0_dp*D_r)/(pi*V_r)

    !! Interpolation function
    del_r = ((2.0_dp*r_h + lam_r)**3 - (4.0_dp*r_h**2 + lam_r**2)**(1.5_dp))/(6.0_dp*r_h*lam_r) &
      & - 2.0_dp*r_h

    !! Correction factor
    phi = 2.0_dp*r_h/(2.0_dp*r_h + sqrt(2.0_dp)*del_r) + (4.0_dp*D_r)/(r_h*sqrt(2.0_dp)*V_r) 

    !! Coagulation flux (Zeroth moment) [cm-3 s-1]
    f_coag = (4.0_dp * kb * T * beta)/(3.0_dp * eta * phi) * y(1)**2

  end subroutine calc_coag

  !! Particle-particle gravitational coalesence
  subroutine calc_coal(n_eq, y, r_h, Kn, vf, f_coal)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), dimension(n_eq), intent(in) :: y 
    real(dp), intent(in) :: r_h, Kn, vf

    real(dp), intent(out) :: f_coal

    real(dp) :: d_vf, Stk, E
    real(dp), parameter :: eps = 0.5_dp

    !! Estimate differential velocity
    d_vf = eps * vf

    !! Calculate E
    if (Kn >= 1.0_dp) then
      !! E = 1 when Kn > 1
      E = 1.0_dp
    else
      !! Calculate Stokes number
      Stk = (vf * d_vf)/(grav * r_h)
      E = max(1.0e-6_dp,1.0_dp - 0.42_dp*Stk**(-0.75_dp))
    end if

    !! Finally calculate the loss flux term
    f_coal = 2.0_dp*pi*r_h**2*y(1)**2*d_vf*E

  end subroutine calc_coal

  !! Monomer production rates from mass mixing ratio rate passed into module
  subroutine find_production_rate(n_eq, f_prod, f_form)
    implicit none

    integer, intent(in) :: n_eq

    real(dp), intent(in) :: f_form

    real(dp), dimension(n_eq), intent(out) :: f_prod

    f_prod(1) = f_form * (rho/(V_mon*rho_d))
    f_prod(2) = f_form * rho
    f_prod(3) = Prod_in
    f_prod(4) = 0.0_dp

  end subroutine find_production_rate

  !! Dummy jacobian subroutine required for dlsode
  subroutine jac_dum (NEQ, X, Y, ML, MU, PD, NROWPD)
    integer, intent(in) :: NEQ, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(NEQ), intent(in) :: Y
    real(dp), dimension(NROWPD, NEQ), intent(inout) :: PD
  end subroutine jac_dum

  !! eta for background gas
  subroutine eta_construct(sp_bg)
    implicit none

    character(len=20), dimension(:), intent(in) :: sp_bg
    
    integer :: n_bg, i

    n_bg = size(sp_bg)

    allocate(d_g(n_bg), LJ_g(n_bg), molg_g(n_bg))

    do i = 1, n_bg
      select case(sp_bg(i))

      case('OH')
        d_g(i) = d_OH
        LJ_g(i) = LJ_OH
        molg_g(i) = molg_OH
      case('H2')
        d_g(i) = d_H2
        LJ_g(i) = LJ_H2
        molg_g(i) = molg_H2
      case('H2O')
        d_g(i) = d_H2O
        LJ_g(i) = LJ_H2O
        molg_g(i) = molg_H2O
      case('H')
        d_g(i) = d_H
        LJ_g(i) = LJ_H
        molg_g(i) = molg_H
      case('CO')
        d_g(i) = d_CO
        LJ_g(i) = LJ_CO
        molg_g(i) = molg_CO
      case('CO2')
        d_g(i) = d_CO2
        LJ_g(i) = LJ_CO2
        molg_g(i) = molg_CO2
      case('O')
        d_g(i) = d_O
        LJ_g(i) = LJ_O
        molg_g(i) = molg_O
      case('CH4')
        d_g(i) = d_CH4
        LJ_g(i) = LJ_CH4
        molg_g(i) = molg_CH4
      case('C2H2')
        d_g(i) = d_C2H2
        LJ_g(i) = LJ_C2H2
        molg_g(i) = molg_C2H2
      case('NH3')
        d_g(i) = d_NH3
        LJ_g(i) = LJ_NH3
        molg_g(i) = molg_NH3
      case('N2')
        d_g(i) = d_N2
        LJ_g(i) = LJ_N2
        molg_g(i) = molg_N2 
      case('HCN')
        d_g(i) = d_HCN
        LJ_g(i) = LJ_HCN
        molg_g(i) = molg_HCN
      case('He')
        d_g(i) = d_He
        LJ_g(i) = LJ_He
        molg_g(i) = molg_He
      case default
        print*, 'Background gas species data not found: ', trim(sp_bg(i)), 'STOP'
        stop
      end select

    end do
    
  end subroutine eta_construct

end module mini_haze_i_dlsode_mom_mod
