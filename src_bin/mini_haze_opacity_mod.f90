module mini_haze_opacity_mod
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomonded if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128

  public :: mini_haze_bin_opacity
  private :: dp

  contains

  subroutine mini_haze_bin_opacity()
    implicit none

  end subroutine mini_haze_bin_opacity

end module mini_haze_opacity_mod