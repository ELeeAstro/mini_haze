program mini_haze_main
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_haze_i_dlsode_mod, only : mini_haze_i_dlsode
  use mini_haze_production_mod, only : f_const
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomonded if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-16_dp


  integer :: n_bin, n_steps, n
  real(dp) :: temp, p, mu, g, t_step, nd_atm
  real(dp) :: F0, mu_z, sig, m0, Prod, VMR
  
  real(dp), allocatable, dimension(:) :: q, vf


  !! Mock atmosphere conditions
  temp = 500.0_dp
  p = 1e-1_dp
  n_bin = 5
  mu = 2.33_dp
  n_steps = 10
  g = 10.0_dp

  !! Timestep (of GCM usually)
  t_step = 30.0_dp

  !! Allocate tracers and fall velocities
  !! q here is the volume mixing ratio
  allocate(q(n_bin), vf(n_bin))
  q(:) = 0.0_dp
  vf(:) = 0.0_dp


  !! Number density of atmosphere (cgs)
  nd_atm = (p*10.0_dp)/(kb*temp)

  !! Set up a mock production rate using Steinrueck et al. (2023)
  mu_z = 1.0_dp
  F0 = 1.0e-11_dp
  sig = 0.25_dp*log(10.0_dp)
  m0 = 2.0e-6_dp * 1e5_dp

  !! Production rate is dimensionless and is a mass mixing ratio
  Prod = F0 * g * mu_z * (1.0_dp/(sqrt(2.0_dp*pi) * p * sig)) &
    & * exp(-(log(p/m0)**2)/(2.0_dp*sig**2))

  !! The resulting production rate is the mass mixing ratio, 
  !! now turn this into a volume mixing ratio
  !! Assume a molecular weight
  VMR = (Prod  * mu)/ (200.0_dp)

  !! f_const is now the constant production rate in monomers cm-3 s-1 
  f_const = VMR * nd_atm

  print*, Prod, VMR, f_const

  do n = 1, n_steps

    call mini_haze_i_dlsode(n_bin, temp, p, mu, g, t_step, vf, q)

    print*, n*t_step, n, q(:), sum(q(:))

  end do


end program mini_haze_main