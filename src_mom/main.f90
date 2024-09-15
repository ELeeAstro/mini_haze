program mini_haze_main
  use, intrinsic :: iso_fortran_env ! Requires fortran 2008
  use mini_haze_i_dlsode_mom_mod, only : mini_haze_i_dlsode_mom,  &
    & Prod_in, r_mon, rho_d, p_deep, tau_loss, tau_act, tau_decay, tau_form
  use mini_haze_vf_mod, only : mini_haze_vf
  use mini_haze_opac_mie_mod, only : opac_mie
  implicit none


  double precision, parameter :: pi = 4.0 * atan(1.0)
  double precision, parameter :: kb = 1.380649e-16
  double precision, parameter :: amu = 1.66054e-24


  integer :: n_eq, n_bin, n_steps, n, u
  double precision :: temp, p, mu, grav, t_step, nd_atm, rho
  double precision :: F0, mu_z, sig, m0

  double precision :: r_h, m_h, v_f
  
  character(len=20) :: sp
  character(len=20) :: sp_bg(3)
  double precision ::  VMR_in(3)
  double precision, allocatable, dimension(:) :: q

  integer :: n_wl
  double precision, allocatable, dimension(:) :: wl_e, wl, k_ext, ssa, g

  !! Mock atmosphere conditions
  temp = 950.0! Temperature [K]
  p = 2e-1 ! Pressure [pa]
  mu = 2.33 ! mean molecular weight [g mol-1]
  grav = 10.0 ! Gravity [m s-2]

  n_bin = 2 ! Number of size bins = 2 for moment method
  n_eq = n_bin + 2 ! Number of size bins + precursor tracers
  n_steps = int(1e5) ! Number of time steps
 
  rho = (p*10.0*mu*amu)/(kb * temp) ! Mass density [g cm-3]
  nd_atm = (p*10.0)/(kb*temp) ! Number density [cm-3]

  print*, 'nd: ', nd_atm, 'rho: ', rho

  !! Timestep (of GCM usually)
  t_step = 60.0

  !! Allocate tracers and fall velocities
  !! q(1) = zeroth moment volume mixing ratio
  !! q(2) = first moment mass mixing ratio
  !! q(3) = Non-activated precursor mass mixing ratio
  !! q(4) = Activated precursor mass mixing ratio
  allocate(q(n_eq))
  q(:) = 1.0e-30

  !! Assume constant H2, He and H background VMR @ approx solar
  sp_bg = (/'H2','He','H '/)
  VMR_in(1) = 0.85
  VMR_in(2) = 0.15
  VMR_in(3) = 1e-6

  !! Setup wavelength grid for opacity calculations
  n_wl = 11
  allocate(wl_e(n_wl+1), wl(n_wl), k_ext(n_wl), ssa(n_wl), g(n_wl))

  !! Wavelengths to calculate opacity
  wl_e = (/0.260, 0.420, 0.610, 0.850, 1.320, 2.020,2.500,3.500,4.400,8.70,20.00,324.68 /)
  wl(:) = (wl_e(2:n_wl+1) +  wl_e(1:n_wl))/ 2.0

  !! Set up a mock production rate using Steinrueck et al. (2023)
  mu_z = 1.0
  F0 = 1.0e-11
  sig = 0.25*log(10.0)
  m0 = 2.0e-6 * 1e5

  !! Production rate is dimensionless and is a mass mixing ratio
  Prod_in = F0 * grav * mu_z * (1.0/(sqrt(2.0*pi) * p * sig)) &
    & * exp(-(log(p/m0)**2)/(2.0*sig**2))

  print*, Prod_in

  !! Species of haze (for opacity calculation)
  sp = 'Tholin'

  !! Send some parameters to the integration module
  r_mon = 1e-7
  rho_d = 1.0
  p_deep = 1e5
  tau_loss = 1e3
  tau_act = 100.0
  tau_decay = 1.0
  tau_form = 100.0

  open(newunit=u,file='mom_test.txt',action='readwrite')

  write(u,*) q(:)

  do n = 1, n_steps

    !! Integrate moments
    call mini_haze_i_dlsode_mom(n_eq, temp, p, mu, grav, VMR_in, t_step, sp_bg, q)

    !! Calculate settling velocity for this layer (v_f [cm s-1])
    call mini_haze_vf(temp, p, grav, mu, VMR_in, rho_d, sp_bg, q(1), q(2), v_f)

    !! Calculate the opacity at the wavelength grid
    call opac_mie(1, sp, temp, mu, p, q(1), q(2), rho_d, n_wl, wl, k_ext, ssa, g)

    !! Mean mass [g]
    m_h = (q(2)*rho)/(q(1)*nd_atm)

    !! Mean radius [cm]
    r_h = ((3.0*m_h)/(4.0*pi*rho_d))**(1.0/3.0)

    ! Shut off production after 1 timestep
    !Prod_in = 0.0

    print*, n*t_step/86400.0
    print*, 'q', n, q(:), v_f
    print*, 'r', n, m_h, r_h * 1e4
    print*, 'o', n, k_ext(1), ssa(1), g(1), k_ext(n_wl), ssa(n_wl), g(n_wl)

    !print*, n*t_step/86400, n, q(:), Prod_in, tau_decay/tau_act + tau_decay/tau_form !* Prod_in

    write(u,*) q(:)

  end do

end program mini_haze_main