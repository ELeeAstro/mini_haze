module mini_haze_i_dlsode_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_haze_production_mod, only : find_production_rate
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
  real(dp) :: nd_atm, rho, P_cgs, Rd, grav_cgs
  logical :: first_call = .True.

  integer :: n_bin
  character(len=20) :: prod_scheme
  real(dp), allocatable, dimension(:) :: re, r ! Radii of bin edge and bin center (cm)
  real(dp), allocatable, dimension(:) :: me, m ! Mass of bin edges and bin center (g)
  real(dp), allocatable, dimension(:,:) :: Kr ! Collisional kernel
  real(dp), allocatable, dimension(:,:,:) :: C ! C value (Brauer et al. 2008)
  real(dp) :: rho_h ! bulk density of haze [g cm-3]

  public :: mini_haze_i_dlsode, RHS_bin, jac_dum
  private :: dp, init_mini_haze_bin, calc_kernel, calc_eps

  contains

  subroutine mini_haze_i_dlsode(n_eq, T_in, P_in, mu_in, grav_in, t_end, vf, q)
    implicit none

    ! Input variables
    integer, intent(in) :: n_eq
    real(dp), intent(in) :: T_in, P_in, mu_in, grav_in, t_end

    ! Input/Output tracer values and settling velocity
    real(dp), dimension(n_eq), intent(inout) :: q
    real(dp), dimension(n_eq), intent(out) :: vf

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
      call calc_eps()
      first_call = .False.
    end if

    !! Find the number density of the atmosphere
    P_cgs = P_in * 10.0_dp   ! Convert pascal to dyne cm-2
    nd_atm = P_cgs/(kb*T_in)  ! Find initial number density [cm-3] of atmosphere

    !! Change gravity to cgs
    grav_cgs = grav_in * 100.0_dp

    !! Find mass density of atmosphere [g cm-3]
    Rd = R_gas/mu_in ! Specific gas constant [erg g-1 K-1]
    rho = P_cgs/(Rd * T_in) ! Mass density [g cm-3]

    call calc_kernel(T_in, rho, mu_in, grav_cgs, vf)
    !call calc_kernel_analy()

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
    atol(:) = 1.0e-30_dp               ! Absolute tolerance for each scalar (floor value)

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

    print*, t_end * 4.0_dp/3.0_dp * pi * r(1)**3 * rho_h

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

    print*, sum(y(:) * 4.0_dp/3.0_dp * pi * r(:)**3 * rho_h)

    stop

  end subroutine mini_haze_i_dlsode

  subroutine RHS_bin(n_eq, time, y, f)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), intent(inout) :: time
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq), intent(inout) :: f

    integer :: i, j, k

    !! In this routine, you calculate the new fluxes (f) into each bin
    !! The values of each bin (y) are typically kept constant
    !! Basically, you solve for the RHS of the ODE for each size bin

    f(:) = 0.0_dp

    !! Find the RHS of the ODE for each particle bin size

    !! For the first bin, add a production rate source term
    call find_production_rate(prod_scheme, f(1))

    !! Calculate the third term in Brauer et al. 2008 (Equation A.7)
    !! This is loss from a bin from collisions
    !! The loss rate depends on the kernel * i number density * j number density
    do k = 1, n_bin
      do j = 1, n_bin
        f(k) = f(k) - y(j) * y(k) * Kr(j,k)
      end do
    end do
    
    !! Calculate the first term in Brauer et al. 2008 (Equation A.7)
    !! This is gain into a bin for self size collisions
    do k = 2, n_bin
      do i = 1, k
        f(k) = f(k) + 0.5_dp * y(i) * y(i) * Kr(i,i) * C(i,i,k)
      end do
    end do

    !! Calculate the second term in Brauer et al. 2008 (Equation A.7)
    !! This is gain into a bin for different sized collisions
    do k = 2, n_bin
      do i = 1, k
        do j = 1, i-1
          f(k) = f(k) + y(i) * y(j) * Kr(i,j) * C(i,j,k)
        end do
      end do
    end do

    !print*,time,f(:)

  end subroutine RHS_bin

  subroutine jac_dum (NEQ, X, Y, ML, MU, PD, NROWPD)
    integer, intent(in) :: NEQ, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(NEQ), intent(in) :: Y
    real(dp), dimension(NROWPD, NEQ), intent(inout) :: PD
  end subroutine jac_dum

  subroutine calc_kernel(T, rho, mu, grav, vf)
    implicit none

    real(dp), intent(in) :: T, rho, mu, grav

    real(dp), dimension(n_bin), intent(out) :: vf

    integer :: i, j
    real(dp), dimension(n_bin,n_bin) :: K_coag, K_coal
    real(dp), dimension(n_bin) :: D, V, lam, del, beta, Kn
    real(dp) :: eta, mfp

    real(dp) :: d_vf, Stk, E

    real(dp), parameter :: d_H2 = 2.827e-8_dp, LJ_H2 = 59.7_dp * kb, molg_H2 = 2.01588_dp

    !! Here we calculate the kernel function between each particle size going up in size

    !! First we calculate the diffusion coefficents, thermal velocities and interpolation factor
    !! for each bin size

    !! Calculate dynamical viscosity for this layer
    eta = (5.0_dp/16.0_dp) * (sqrt(pi*(molg_H2*amu)*kb*T)/(pi*d_H2**2)) &
      & * (((kb*T)/LJ_H2)**(0.16_dp)/1.22_dp)

    !print*, eta

    !! Calculate mean free path for this layer
    mfp = (2.0_dp*eta/rho) * sqrt((pi * mu)/(8.0_dp*R_gas*T))

    !print*, mfp

    do i = 1, n_bin

      !! Knudsen number
      Kn(i) = mfp/r(i)

      !print*, Kn(i)

      !! Cunningham slip factor
      beta(i) = 1.0_dp + 1.246_dp*Kn(i) + 0.42_dp*Kn(i) * exp(-0.87_dp/Kn(i))

      !print*,beta(i)

      D(i) = (kb*T*beta(i))/(6.0_dp*pi*eta*r(i))

      !print*, D(i)

      V(i) = sqrt((8.0_dp*kb*T)/(pi*m(i)))

      !print*, V(i)

      lam(i) = (8.0_dp*D(i))/(pi * V(i))

      !print*, lam(i)

      del(i) = ((2.0_dp*r(i) + lam(i))**3 - (4.0_dp*r(i)**2 + lam(i)**2)**(3.0_dp/2.0_dp)) &
        & / (6.0_dp * r(i) * lam(i)) - 2.0_dp * r(i) 

      !print*, i, eta, mfp, Kn(i), beta(i), D(i), m(i), V(i), lam(i), del(i)

    end do

    K_coag(:,:) = 0.0_dp

    !! Calculate the kernel matrix for Brownian motion 
    !! (coagulation)
    do i = 1, n_bin
      do j = 1, n_bin

        ! if (j < i) then
        !   !! We don't calculate the kernel for particle sizes less the current bin
        !   cycle
        ! end if

        K_coag(i,j) = (4.0_dp*pi*(D(i) + D(j))*(r(i) + r(j))) & 
          & / ((r(i) + r(j))/(r(i) + r(j) + sqrt(del(i)**2 + del(j)**2)) &
          & + (4.0_dp*(D(i) + D(j)))/((r(i) + r(j)) * sqrt(V(i)**2 + V(j)**2)))

        !print*, i, j, Kr(i,j)

      end do
    end do

    !! Calculate settling velocity for each size bin
    vf(:) = 0.0_dp
    do i = 1, n_bin

      beta(i) = 1.0_dp + Kn(i)*(1.257_dp + 0.4_dp * exp(-1.1_dp/Kn(i)))

      vf(i) = (2.0_dp * beta(i) * grav * r(i)**2 * rho_h)/(9.0_dp * eta) & 
        & * (1.0_dp + & 
        & ((0.45_dp*grav*r(i)**3*rho*rho_h)/(54.0_dp*eta**2))**(2.0_dp/5.0_dp))**(-5.0_dp/4.0_dp)

    end do

    K_coal(:,:) = 0.0_dp
 
    !! Calculate the kernel matrix for settling velocity diferentiation
    !! (coalescence)
    do i = 1, n_bin
      do j = 1, n_bin

        if (j >= i) then
          !! Particle size is greater or equal than current bin
          cycle
        end if

        d_vf = abs(vf(i) - vf(j))

        Stk = (vf(j)*d_vf)/(r(i)*grav)
        if (Kn(j) >= 1.0_dp) then
          E = 1.0_dp
        else 
          E = max(0.0_dp,1.0_dp - 0.42_dp*Stk**(-0.75_dp))
        end if

        K_coal(i,j) = 0.0_dp !pi * (r(i) + r(j))**2 * d_vf * E

      end do
    end do

    !! Total kernel is coagulation + coalescence
    Kr(:,:) = K_coag(:,:) + K_coal(:,:)

  end subroutine calc_kernel

  subroutine calc_kernel_analy()
    implicit none

    integer :: i, j

    do i = 1, n_bin
      do j = 1, n_bin
        Kr(i,j) = 1.0_dp!m(i) + m(j)
      end do
    end do

  end subroutine calc_kernel_analy

  subroutine calc_eps()
    implicit none

    integer :: i, j, k

    real(dp) :: eps, msum

    !! Allocate eps array
    allocate(C(n_bin,n_bin,n_bin))
    !! Calculate the epsilon (C_ijk) matrix from Brauer et al. 2008

    C(:,:,:) = 0.0_dp

    do i = 1, n_bin
      do j = 1, n_bin
        msum = m(i) + m(j)
        do k = 2, n_bin
          !! Find nearest neighbours
          if ((msum > m(k-1)) .and. (msum <= m(k))) then
            !! Sum of collision mass lies between k-1 and k mass bin
            eps = (m(k) - msum)/(m(k) - m(k-1))
            C(i,j,k) = 1.0_dp - eps
            C(i,j,k-1) = eps
            
          end if
        end do
      end do
    end do

  end subroutine calc_eps

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
    allocate(Kr(n_bin,n_bin))
    
    !! Allocate other arrays
    allocate(re(n_bin+1), r(n_bin), me(n_bin+1), m(n_bin))

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

end module mini_haze_i_dlsode_mod