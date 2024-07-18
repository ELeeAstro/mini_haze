program mini_haze_main
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_haze_i_dlsode_mod, only : mini_haze_i_dlsode
  use mini_haze_production_mod
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomonded if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128

  real(dp), parameter :: kb = 1.380649e-16_dp


  integer :: n_bin, n_steps, n
  real(dp) :: temp, p, mu, g, t_step, nd_atm
  
  real(dp), allocatable, dimension(:) :: q, vf

  temp = 500.0_dp
  p = 1e-1_dp
  n_bin = 31
  mu = 2.33_dp
  n_steps = 1
  g = 10.0_dp

  t_step = 30.0_dp

  allocate(q(n_bin), vf(n_bin))
  q(:) = 0.0_dp
  vf(:) = 0.0_dp


  !nd_atm = (p*10.0_dp)/(kb*temp)

  !q(1) = 1.0_dp/nd_atm

  do n = 1, n_steps

    ! Set constant rate -> stored in production module
    f_const = 1.0_dp

    call mini_haze_i_dlsode(n_bin, temp, p, mu, g, t_step, vf, q)

    print*, n*t_step, n, q(:), sum(q(:))

  end do


end program mini_haze_main