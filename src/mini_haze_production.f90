module mini_haze_production
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomonded if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128


  public :: find_production_rate
  private :: dp, constant_rate

contains

  !! Landing routine for the production rate schemes
  subroutine find_production_rate(prod_scheme, f1)
    implicit none

    character(len=20), intent(in) :: prod_scheme

    real(dp), intent(out) :: f1


    select case (prod_scheme)
    case('constant')
      call constant_rate()
    end select

  end subroutine find_production_rate

  subroutine constant_rate()
    implicit none


  end subroutine constant_rate


end module mini_haze_production