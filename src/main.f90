program mini_haze_main
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_haze_integration, only : mini_haze_i_dlsode
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomonded if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128


  integer :: n_bin
  real(dp) :: temp, p, mu, t_step
  
  real(dp), allocatable, dimension(:) :: q

  temp = 1000.0_dp
  p = 1e5_dp
  n_bin = 5
  mu = 2.33

  allocate(q(n_bin))
  q(:) = 1e-30_dp

  call mini_haze_i_dlsode(n_bin, temp, p, mu, t_step, q)


end program mini_haze_main