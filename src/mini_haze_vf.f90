module mini_haze_vf
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomonded if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128


  public :: calc_haze_vf
  private :: dp

  contains

  subroutine calc_haze_vf()
    implicit none

  end subroutine calc_haze_vf


end module mini_haze_vf