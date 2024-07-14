module mini_haze_production
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomonded if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128


  public :: constant_rate
  private :: dp

contains

  subroutine constant_rate
    implicit none


  end subroutine constant_rate


end module mini_haze_production