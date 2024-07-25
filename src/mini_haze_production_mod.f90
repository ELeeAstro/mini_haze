module mini_haze_production_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomonded if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)

  real(dp), allocatable, dimension(:) :: f_prod

  logical :: first_call = .True.

  public :: mini_haze_prod, find_production_rate
  private :: dp

contains

  !! Landing routine for the production rate schemes
  !! Called before main mini-haze integration to save production rates 
  subroutine mini_haze_prod(n_eq, f_prod_in)
    implicit none

    integer, intent(in) :: n_eq

    real(dp), dimension(n_eq), intent(in) :: f_prod_in

    if (first_call .eqv. .True.) then
      allocate(f_prod(n_eq))
      first_call = .False.
    end if

    f_prod(:) = f_prod_in(:)

  end subroutine mini_haze_prod

  !! Routine called by the mini-haze integrator
  subroutine find_production_rate(n_eq, f_prod_int)
    implicit none

    integer, intent(in) :: n_eq

    real(dp), dimension(n_eq), intent(out) :: f_prod_int

    f_prod_int(:) = f_prod(:)

  end subroutine find_production_rate

end module mini_haze_production_mod