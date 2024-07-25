program mini_haze_main
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  !use mini_haze_i_dlsode_mod, only : mini_haze_i_dlsode, m
  use mini_haze_i_dlsode_mom_mod, only : mini_haze_i_dlsode_mom, rho_h, p_deep, tau_loss
  use mini_haze_production_mod, only : mini_haze_prod
  implicit none

  ! Fortran 2008 intrinsic precisions - reccomonded if possible
  !integer, parameter :: sp = REAL32
  integer, parameter :: dp = REAL64
  !integer, parameter :: qp = REAL128

  real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
  real(dp), parameter :: kb = 1.380649e-16_dp
  real(dp), parameter :: amu = 1.66054e-24_dp


  integer :: n_bin, n_steps, n, u
  real(dp) :: temp, p, mu, g, t_step, nd_atm, rho
  real(dp) :: F0, mu_z, sig, m0, Prod

  real(dp) :: bulk_den, a_seed, r
  
  real(dp), allocatable, dimension(:) :: q, vf, f_prod

  !! Mock atmosphere conditions
  temp = 700.0_dp ! Temperature [K]
  p = 2e-1_dp ! Pressure [pa]
  mu = 2.33_dp ! mean molecular weight [g mol-1]
  g = 10.0_dp ! Gravity [m s-2]

  n_bin = 2 ! Number of size bins = 2 for moment method
  n_steps = 1000 ! Number of time steps
 
  rho = (p*10.0_dp*mu*amu)/(kb * temp) ! Mass density [g cm-3]
  nd_atm = (p*10.0_dp)/(kb*temp) ! Number density [cm-3]

  print*, 'nd: ', nd_atm, 'rho: ', rho

  !! Timestep (of GCM usually)
  t_step = 100.0_dp

  !! Allocate tracers and fall velocities
  !! q here is the volume mixing ratio
  allocate(q(n_bin), vf(n_bin))
  q(:) = 1.0e-30_dp
  vf(:) = 0.0_dp

  !! Give initial values analytical solution

  !! Set up a mock production rate using Steinrueck et al. (2023)
  mu_z = 1.0_dp
  F0 = 1.0e-11_dp
  sig = 0.25_dp*log(10.0_dp)
  m0 = 2.0e-6_dp * 1e5_dp

  !! Production rate is dimensionless and is a mass mixing ratio
  Prod = F0 * g * mu_z * (1.0_dp/(sqrt(2.0_dp*pi) * p * sig)) &
    & * exp(-(log(p/m0)**2)/(2.0_dp*sig**2))

  !! Haze parameters
  bulk_den = 2.0_dp
  a_seed = 1e-7_dp
  
  allocate(f_prod(n_bin))
  !! Find the production rate for each moments/bin size 
  f_prod(1) = (3.0_dp * Prod)/(4.0_dp*pi*a_seed**3*bulk_den) * rho
  f_prod(2) = Prod * rho

  print*, f_prod(:), a_seed * 1e4_dp

  !! Send some parameters to the integration module
  rho_h = 2.0_dp
  p_deep = 1e5_dp
  tau_loss = 1e3_dp

  open(newunit=u,file='test.txt',action='readwrite')

  write(u,*) q(:)

  do n = 1, n_steps

    call mini_haze_prod(n_bin, f_prod)
    call mini_haze_i_dlsode_mom(n_bin, temp, p, mu, g, t_step, vf, q)

    r = ((3.0_dp*(q(2)*rho)/(q(1)*nd_atm))/(4.0_dp*pi*bulk_den))**(1.0_dp/3.0_dp)

    print*, n*t_step, n, q(:), r * 1e4_dp

    write(u,*) q(:)

  end do

end program mini_haze_main