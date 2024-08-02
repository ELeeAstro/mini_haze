program mini_haze_main
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_haze_i_dlsode_mom_mod, only : mini_haze_i_dlsode_mom,  &
    & Prod_in, r_mon, rho_h, p_deep, tau_loss, tau_act, tau_decay, tau_form 
  implicit none


  double precision, parameter :: pi = 4.0 * atan(1.0)
  double precision, parameter :: kb = 1.380649e-16
  double precision, parameter :: amu = 1.66054e-24


  integer :: n_eq, n_bin, n_steps, n, u
  double precision :: temp, p, mu, g, t_step, nd_atm, rho
  double precision :: F0, mu_z, sig, m0

  double precision :: bulk_den, r
  
  double precision, allocatable, dimension(:) :: q

  !! Mock atmosphere conditions
  temp = 700.0! Temperature [K]
  p = 2e-1 ! Pressure [pa]
  mu = 2.33 ! mean molecular weight [g mol-1]
  g = 10.0 ! Gravity [m s-2]

  n_bin = 2 ! Number of size bins = 2 for moment method
  n_eq = n_bin + 2 ! Number of size bins + precursor tracers
  n_steps = int(1e5) ! Number of time steps
 
  rho = (p*10.0*mu*amu)/(kb * temp) ! Mass density [g cm-3]
  nd_atm = (p*10.0)/(kb*temp) ! Number density [cm-3]

  print*, 'nd: ', nd_atm, 'rho: ', rho

  !! Timestep (of GCM usually)
  t_step = 100.0

  !! Allocate tracers and fall velocities
  !! q here is the volume mixing ratio
  allocate(q(n_eq))
  q(:) = 1.0e-30

  !! Give initial values analytical solution

  !! Set up a mock production rate using Steinrueck et al. (2023)
  mu_z = 1.0
  F0 = 1.0e-11
  sig = 0.25*log(10.0)
  m0 = 2.0e-6 * 1e5

  !! Production rate is dimensionless and is a mass mixing ratio
  Prod_in = F0 * g * mu_z * (1.0/(sqrt(2.0*pi) * p * sig)) &
    & * exp(-(log(p/m0)**2)/(2.0*sig**2))

  print*, Prod_in

  !! Send some parameters to the integration module
  r_mon = 1e-7
  rho_h = 2.0
  p_deep = 1e5
  tau_loss = 1e3
  tau_act = 1.0
  tau_decay = 10000.0
  tau_form = 1.0

  bulk_den = rho_h

  open(newunit=u,file='mom_test.txt',action='readwrite')

  write(u,*) q(:)

  do n = 1, n_steps

    call mini_haze_i_dlsode_mom(n_eq, temp, p, mu, g, t_step, q)

    r = ((3.0*(q(2)*rho)/(q(1)*nd_atm))/(4.0*pi*bulk_den))**(1.0/3.0)

    print*, n*t_step, n, q(:), r * 1e4, Prod_in

    write(u,*) q(:)

  end do

end program mini_haze_main