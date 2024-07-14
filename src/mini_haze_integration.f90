module mini_haze_integration
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomended if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128

  ! Constants
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: R = 8.31446261815324e7_dp
  real(dp), parameter :: P0 = 1.0e6_dp

  ! Global variables
  real(dp) :: nd_atm, P_cgs
  logical :: first_call = .True.

  public :: RHS_bin, jac_dum
  private :: dp, init_mini_haze_bin

  contains

  subroutine mini_haze_i_dlsode(n_eq, T_in, P_in, t_end, q)
    implicit none

    ! Input variables
    integer, intent(in) :: n_eq
    real(dp), intent(in) :: T_in, P_in, t_end

    ! Input/Output tracer values
    real(dp), dimension(n_eq), intent(inout) :: q

    integer :: ncall

    ! Time controls
    real(dp) :: t_now

    ! DLSODE variables
    real(dp), dimension(n_eq) :: y
    real(dp), allocatable, dimension(:) :: rwork, rtol, atol
    integer, allocatable, dimension(:) :: iwork
    integer :: itol, itask, istate, iopt, mf
    integer :: rworkdim, iworkdim

    !! To initialisations for first call
    if (first_call .eqv. .True.) then
      call init_mini_haze_bin()
      first_call = .False.
    end if

    !! Find the number density of the atmosphere
    P_cgs = P_in * 10.0_dp   ! Convert pascal to dyne cm-2
    nd_atm = P_cgs/(kb*T_in)  ! Find initial number density [cm-3] of atmosphere

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

  end subroutine mini_haze_i_dlsode



  subroutine RHS_bin(n_eq, time, y, f)
    implicit none

    integer, intent(in) :: n_eq
    real(dp), intent(inout) :: time
    real(dp), dimension(n_eq), intent(inout) :: y
    real(dp), dimension(n_eq), intent(inout) :: f

  end subroutine RHS_bin

  subroutine jac_dum (NEQ, X, Y, ML, MU, PD, NROWPD)
    integer, intent(in) :: NEQ, ML, MU, NROWPD
    real(dp), intent(in) :: X
    real(dp), dimension(NEQ), intent(in) :: Y
    real(dp), dimension(NROWPD, NEQ), intent(inout) :: PD
  end subroutine jac_dum

  subroutine init_mini_haze_bin()
    implicit none
  end subroutine init_mini_haze_bin

end module mini_haze_integration