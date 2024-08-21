module mini_haze_i_dlsode_bin_mod
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
  real(dp) :: T, mu, nd_atm, rho, p, grav, mfp, eta
  logical :: first_call = .True.

  !! Variables needed to be sent in from outside
  real(dp) :: Prod_in ! Mass mixing ratio production rate of precurser molecules
  real(dp) :: r_mon ! Haze particle monomer size
  real(dp) :: rho_d ! bulk density of haze [g cm-3]
  real(dp) :: p_deep ! Deep removal pressure level
  real(dp) :: tau_loss ! Deep removal loss timescale
  real(dp) :: tau_decay, tau_act, tau_form ! Precurser timescales

  real(dp) :: V_mon ! Haze particle monomer volume
  real(dp) :: m_mon ! Haze particle monomer mass

  !! Bin calculation variables
  integer :: n_bin
  real(dp) :: r_min, r_max, m_min, m_max ! Max and min bin radius and mas
  real(dp), allocatable, dimension(:) :: re, r, me, m ! radius and mass bin edges and center
  real(dp), allocatable, dimension(:,:) :: Kr




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
  real(dp), dimension(3) :: LJ_g = (/LJ_H2, LJ_He, LJ_H/)
  real(dp), dimension(3) :: molg_g = (/molg_H2, molg_He, molg_H/)

  public :: mini_haze_i_dlsode_bin, RHS_bin, RHS_analy, jac_dum
  private :: find_production_rate

  contains

  subroutine mini_haze_i_dlsode_bin(n_eq, T_in, P_in, mu_in, grav_in, t_end, q, vf, n_gas, VMR_g)
    implicit none

    ! Input variables
    integer, intent(in) :: n_eq, n_gas
    real(dp), intent(in) :: T_in, P_in, mu_in, grav_in, t_end
    real(dp), dimension(n_gas), intent(in) :: VMR_g

    ! Input/Output tracer values
    real(dp), dimension(n_eq), intent(inout) :: q
    real(dp), dimension(n_eq-2), intent(out) :: vf

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
    integer :: n_bin
    real(dp) :: eta_g, bot, top
    real(dp), allocatable, dimension(:) ::  Kn, beta

    !! To initialisations for first call
    if (first_call .eqv. .True.) then
      n_bin = n_eq - 2
      call init_mini_haze_bin(n_eq)
      call calc_eps()
      first_call = .False.
    end if

    !! Alter input values to mini-haze units
    !! (note, some are obvious not not changed in case specific models need different conversion factors)

    !! Find the number density of the atmosphere
    T = T_in             ! Convert temperature to K
    p = P_in * 10.0_dp   ! Convert pascal to dyne cm-2

    !! Change mu_in to mu
    mu = mu_in ! Convert mean molecular weight to mu [g mol-1]

    !! Change gravity to cgs [cm s-2]
    grav = grav_in * 100.0_dp

    !! Number density [cm-3] of layer
    nd_atm = p/(kb*T)  

    !! Mass density of layer
    rho = (p*mu*amu)/(kb * T) ! Mass density [g cm-3]

    !! Calculate dynamical viscosity for this layer - do square root mixing law from Rosner 2012
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

    !! Calculate mean free path for this layer
    mfp = (2.0_dp*eta/rho) * sqrt((pi * mu)/(8.0_dp*R_gas*T))

    !! Haze particle monomer volume
    V_mon = 4.0_dp/3.0_dp * pi * r_mon**3

    !! Haze particle momomer mass
    m_mon = rho_d * V_Mon

    !! Since we know the bin centers, we can pre-calculate some values for each bin
    !! Knudsen number
    allocate(Kn(n_bin))
    Kn(:) = mfp/r(:)

    !! Cunningham slip factor
    allocate(beta(n_bin))
    beta(:) = 1.0_dp + Kn(:)*(1.257_dp + 0.4_dp * exp(-1.1_dp/Kn(:)))

    !! Settling velocity
    vf(:) = (2.0_dp * beta * grav * r(:)**2 * rho_d)/(9.0_dp * eta) & 
      & * (1.0_dp + & 
      & ((0.45_dp*grav*r(:)**3*rho*rho_d)/(54.0_dp*eta**2))**(0.4_dp))**(-1.25_dp)


    !! Due to the equation set, we can pre-calculate the collison Kernel before integration
    !! Analytical collisional kernel
    allocate(Kr(n_bin,n_bin))
    Kr(:,:) = 1.0_dp

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
    atol(:) = 1.0e-30_dp               ! Absolute tolerance for each scalar (floor value)

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

      call DLSODE (RHS_analy, n_eq, y, t_now, t_end, itol, rtol, atol, itask, &
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
    q(:) = max(y(:), 1e-30_dp)

    deallocate(y, Kn, Kr, beta, rtol, atol, rwork, iwork)

    end subroutine mini_haze_i_dlsode_bin

  subroutine RHS_bin(n_eq, time, y, f)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), intent(inout) :: time
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq), intent(inout) :: f

    real(dp) :: f_coal, f_coag
    real(dp) :: f_act, f_decay_pre, f_decay_act, f_form
    real(dp), dimension(n_eq-2) :: f_loss
    real(dp), dimension(2) ::  f_prod

    !! In this routine, you calculate the new fluxes (f) for each bin
    !! The values of each bin (y) are typically kept constant
    !! Basically, you solve for the RHS of the ODE for each moment

    !! Add thermal decomposition loss term if above given pressure level (pa)
    if (p > p_deep) then
      f_loss(:) = y(1:n_bin)/tau_loss
    else
      f_loss(:) = 0.0_dp
    end if

    !! Calculate precursor and activated molecules rates
    f_act = y(n_eq-1)/tau_act
    f_decay_pre = y(n_eq-1)/tau_decay
    f_decay_act = y(n_eq)/tau_decay
    f_form = y(n_eq)/tau_form

    !! Add a production rate source term for each moment and precursor activation
    call find_production_rate(f_prod, f_form)

    !! Calculate final net flux rate for each bin
    !! First bin contains flux of monomers from production rate
    f(1) = f_prod(1)

    !! Rest contain net flux in-out of bin


    !! Calculate final net flux rate for precusor molecules and activation
    f(n_eq-1) = f_prod(2) - f_act - f_decay_pre
    f(n_eq) = f_act - f_form - f_decay_act


    !print*, 'y', time, y(:), ((3.0_dp*y(2)/y(1))/(4.0_dp*pi*rho_d))**(1.0_dp/3.0_dp) * 1e4_dp
    !print*, 'f', f(:), f_prod(:)
    !print*, 'f2', f_coal, f_coag, f_loss, f_act, f_decay_pre, f_decay_act, f_form

  end subroutine RHS_bin

  !! Analytical testing RHS scheme
  subroutine RHS_analy(n_eq, time, y, f)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), intent(inout) :: time
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq), intent(inout) :: f




  end subroutine RHS_analy

  !! Routine called by the mini-haze integrator - calculates production rates from mass mixing ratio rate
  subroutine find_production_rate(f_prod, f_form)
    implicit none

    real(dp), intent(in) :: f_form

    real(dp), dimension(2), intent(out) :: f_prod

    f_prod(1) = 0.0_dp
    f_prod(2) = Prod_in

  end subroutine find_production_rate

  subroutine init_mini_haze_bin(n_eq)
    implicit none

     integer, intent(in) :: n_eq

    integer :: u_nml, i
    real(dp) :: lrmin, lrmax, rmin, rmax
    real(dp) :: lmmin, lmmax, mmin, mmax

    ! Read the namelist to get the mini-haze parameters

    !! Allocate other arrays
    allocate(re(n_bin+1), r(n_bin), me(n_bin+1), m(n_bin))

    !! For testing do other way round
    lmmin = log10(mmin)
    lmmax = log10(mmax)
    do i = 1, n_bin+1
      me(i) = 10.0_dp**((lmmax-lmmin) * real(i-1,dp) / real(n_bin+1-1,dp) + lmmin)
    end do

    !! Bin centers are the central value of each bin edge
    m(:) = (me(1:n_bin) + me(2:n_bin+1))/2.0_dp

    !! Calculate radii of bin edges and center
    re(:) = ((3.0*me(:))/(4.0_dp*pi*rho_d))**(1.0_dp/3.0_dp)
    r(:) = ((3.0*m(:))/(4.0_dp*pi*rho_d))**(1.0_dp/3.0_dp) 

    print*,re(:)*1e4_dp
    print*,r(:)*1e4_dp

    print*, me(:)
    print*, m(:)

    return

    !! Calculate log spaced values between rmin and rmax
    !! for the bin edges - convert to cm here
    rmin = rmin * 1e-4_dp
    rmax = rmax * 1e-4_dp
    lrmin = log10(rmin)
    lrmax = log10(rmax)
    do i = 1, n_bin+1
      re(i) = 10.0_dp**((lrmax-lrmin) * real(i-1,dp) / real(n_bin+1-1,dp) + lrmin)
    end do

    !! Bin centers are the central value of each bin edge
    r(:) = (re(1:n_bin) + re(2:n_bin+1))/2.0_dp

    !! Calculate masses of bin edges and center
    me(:) = 4.0_dp/3.0_dp * pi * re(:)**3 * rho_h
    m(:) = 4.0_dp/3.0_dp * pi * r(:)**3 * rho_h

    print*, rmin*1e4_dp, rmax*1e4_dp
    print*,re(:)*1e4_dp
    print*,r(:)*1e4_dp

    print*, me(:)
    print*, m(:)

  end subroutine init_mini_haze_bin

  subroutine jac_dum (NEQ, X, Y, ML, MU, PD, NROWPD)
    integer, intent(in) :: NEQ, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(NEQ), intent(in) :: Y
    real(dp), dimension(NROWPD, NEQ), intent(inout) :: PD
  end subroutine jac_dum

end module mini_haze_i_dlsode_bin_mod