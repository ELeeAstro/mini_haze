module mini_haze_integration
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomended if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128

  ! Constants
  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: R_gas = 8.31446261815324e7_dp
  real(dp), parameter :: amu = 1.66053906892e-24_dp

  ! Global variables
  real(dp) :: nd_atm, rho, P_cgs, Rd
  logical :: first_call = .True.

  integer :: n_bin
  character(len=20) :: prod_scheme
  real(dp), allocatable, dimension(:) :: re, r ! Radii of bin edge and bin center
  real(dp), allocatable, dimension(:,:) :: K ! Collisional kernel
  real(dp) :: rho_h ! bulk density of haze [g cm-3]

  public :: RHS_bin, jac_dum
  private :: dp, init_mini_haze_bin

  contains

  subroutine mini_haze_i_dlsode(n_eq, T_in, P_in, mu_in, t_end, q)
    implicit none

    ! Input variables
    integer, intent(in) :: n_eq
    real(dp), intent(in) :: T_in, P_in, mu_in, t_end

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

    !! To initialisations for first call
    if (first_call .eqv. .True.) then
      call init_mini_haze_bin(n_eq)
      first_call = .False.
    end if

    !! Find the number density of the atmosphere
    P_cgs = P_in * 10.0_dp   ! Convert pascal to dyne cm-2
    nd_atm = P_cgs/(kb*T_in)  ! Find initial number density [cm-3] of atmosphere

    !! Find mass density of atmosphere [g cm-3]
    Rd = R_gas/mu_in ! Specific gas constant [erg g-1 K-1]
    rho = P_cgs/(Rd * T_in) ! Mass density [g cm-3]

    call calc_kernel(T_in, rho, mu_in)

    stop

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
    rworkdim = 22 +  9*n_eq + n_eq**2
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


    allocate(y(n_bin))

    !! Give tracer values to y - convert to number density
    y(:) = q(:)*nd_atm

    t_now = 0.0_dp

    ! Set the printing flag
    ! 0 = no printing, 1 = printing
    call xsetf(1)

    ncall = 0

    do while (t_now < t_end)

      call DLSODE (RHS_bin, n_eq, y, t_now, t_end, itol, rtol, atol, itask, &
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

    !! Give y values to tracers - convert to number mixing ratio
    q(:) = y(:)/nd_atm

  end subroutine mini_haze_i_dlsode

  subroutine RHS_bin(n_eq, time, y, f)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), intent(inout) :: time
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq), intent(inout) :: f

    integer :: i

    !! In this routine, you calculate the new fluxes (f) into each bin
    !! The values of each bin (y) are typically kept constant
    !! Basically, you solve for the RHS of the ODE for each size bin

    f(:) = 0.0_dp


    !! Find the RHS of the ODE for each particle bin size
    do i = 1, n_bin

      if (i == 1) then
        !! For the first bin, add a production rate source term
        call find_production_rate(prod_scheme, f(i))
      end if


      !! The rate depends on the kernel * i number density * j number density

      end do
    end do



  end subroutine RHS_bin

  subroutine jac_dum (NEQ, X, Y, ML, MU, PD, NROWPD)
    integer, intent(in) :: NEQ, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(NEQ), intent(in) :: Y
    real(dp), dimension(NROWPD, NEQ), intent(inout) :: PD
  end subroutine jac_dum

  subroutine calc_kernel(T, rho, mu)
    implicit none

    real(dp), intent(in) :: T, rho, mu

    integer :: i, j
    real(dp), dimension(n_bin) :: D, V, lam, del, m, beta, Kn
    real(dp) :: eta, mfp

    real(dp), parameter :: d_H2 = 2.827e-8_dp, LJ_H2 = 59.7_dp * kb, molg_H2 = 2.01588_dp

    !! Here we calculate the kernel function between each particle size going up in size

    !! First we calculate the diffusion coefficents, thermal velocities and interpolation factor
    !! for each bin size

    !! Calculate dynamical viscosity for this layer
    eta = (5.0_dp/16.0_dp) * (sqrt(pi*(molg_H2*amu)*kb*T)/(pi*d_H2**2)) &
      & * (((kb*T)/LJ_H2)**(0.16_dp)/1.22_dp)

    print*, eta

    !! Calculate mean free path for this layer
    mfp = (2.0_dp*eta/rho) * sqrt((pi * mu)/(8.0_dp*R_gas*T))

    print*, mfp

    do i = 1, n_bin

      !! Knudsen number
      Kn(i) = mfp/r(i)

      print*, Kn(i)

      !! Cunningham slip factor
      beta(i) = 1.0_dp + 1.246_dp*Kn(i) + 0.42_dp*Kn(i) * exp(-0.87_dp/Kn(i))

      print*,beta(i)

      D(i) = (kb*T*beta(i))/(6.0_dp*pi*eta*r(i))

      print*, D(i)

      !! calculate mass of the particle

      m(i) = 4.0_dp/3.0_dp * pi * r(i)**3 * rho_h

      V(i) = sqrt((8.0_dp*kb*T)/(pi*m(i)))

      print*, V(i)

      lam(i) = (8.0_dp*D(i))/(pi * V(i))

      print*, lam(i)

      del(i) = ((2.0_dp*r(i) + lam(i))**3 - (4.0_dp*r(i)**2 + lam(i)**2)**(3.0_dp/2.0_dp)) &
        & / (6.0_dp * r(i) * lam(i)) - 2.0_dp * r(i) 

      print*, i, eta, mfp, Kn(i), beta(i), D(i), m(i), V(i), lam(i), del(i)

    end do

    !! Calculate the kernel matrix
    do i = 1, n_bin
      do j = 1, n_bin

        if (j < i) then
          !! We don't calculate the kernel for particle sizes less the current bin
          cycle
        end if

        K(i,j) = (4.0_dp*pi*(D(i) + D(j))*(r(i) + r(j))) & 
          & / ((r(i) + r(j))/(r(i) + r(j) + sqrt(del(i)**2 + del(j)**2)) &
          & + (4.0_dp*(D(i) + D(j)))/((r(i) + r(j)) * sqrt(V(i)**2 + V(j)**2)))

        print*, i, j, K(i,j)

      end do
    end do

  end subroutine calc_kernel

  subroutine init_mini_haze_bin(n_eq)
    implicit none

     integer, intent(in) :: n_eq

    integer :: u_nml, i
    real(dp) :: lrmin, lrmax, rmin, rmax

    namelist /mini_haze_nml/ prod_scheme, rho_h, rmin, rmax

    n_bin = n_eq

    open(newunit=u_nml, file='mini_haze.nml', status='old', action='read')
    read(u_nml, nml=mini_haze_nml)
    close(u_nml)

    ! Read the namelist to get the mini-haze parameters

    !! Allocate the collisional kernel array
    allocate(K(n_bin,n_bin))
    
    !! Allocate other arrays
    allocate(re(n_bin+1), r(n_bin))

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

    print*, rmin*1e4_dp, rmax*1e4_dp
    print*,re(:)*1e4_dp
    print*,r(:)*1e4_dp

  end subroutine init_mini_haze_bin

end module mini_haze_integration