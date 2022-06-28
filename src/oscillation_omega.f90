!
!---------------------------------------------!
! Computation of the expected oscillation due !
! to coupling between the core and the mantle !
!---------------------------------------------!
!
! Written by Sevan Adourian and Harriet Lau
! January 2021
! email: sevan.adourian@berkeley.edu
!
!-----------
program main
  !-----------
  use gravitational_potential
  use volumetric_integral
  use moment_of_inertia
  !
  implicit none
  !
  ! argv, argc
  character(len=100) :: arg4, arg5

  character(len=100) :: file_radius
  character(len=100) :: file_model_1d
  character(len=100) :: file_del_rho
  character(len=100) :: file_del_topo

  integer :: nlayer, l, m, i
  integer :: ndisc, NR, N_ICB, N_CMB

  integer, allocatable  :: disc(:)
  real*8,  allocatable  :: rdisc(:)

  real*8, allocatable :: rho_1d(:), rad_norm(:), rho_prime(:)
  real*8, allocatable :: s1(:), s2(:), s3(:)
  real*8              :: phi

  ! Kernel
  complex*16, allocatable :: delta_rho(:), delta_d(:)
  integer :: which_chunk
  integer :: layer_perturb_start, layer_perturb_end, disc_perturb, iostat
  real*8  :: value_perturb, epsilon_CMB, epsilon_pert
  logical :: kern_rho, kern_topo
  ! DEBUG
  
  ! Integration
  real*8 :: moment, MOI_FLUID, MOI_NO_FLUID
  real*8 :: total_potential, delta_v_dim, delta_v_norm
  real*8 :: epsilon_surf, prefactor, rho_eps, sig, sig_sq
  real*8 :: ACC_NORM
  real*8, parameter :: PI = 3.1415927, R_EARTH = 6371.d3, RHO_AV = 5510.d0, GRAV_CST = 6.67408d-11
  integer :: l_chos = 2, m_chos = 2


  ! if (iargc) etc etc
  call get_command_argument(1, file_model_1d) ! 1-d model to consider
  call get_command_argument(2, file_del_rho) ! file for delta_rho
  call get_command_argument(3, file_del_topo) ! file for delta_topo
  call get_command_argument(4, arg4) ! compute kernel density flag
  call get_command_argument(5, arg5) ! compute kernel topography flag
  
  ! Convert arg4 and arg5 to logical
  read(arg4,*,iostat=iostat)  kern_rho
  if (iostat .ne. 0) then
     print *,'Argument conversion for ',arg4,' failed!'
     stop
  endif
  read(arg5,*,iostat=iostat)  kern_topo
  if (iostat .ne. 0) then
     print *,'Argument conversion for ',arg5,' failed!'
     stop
  endif
  
  if (kern_rho) then
     open(100, file='./kernel_rho.txt', form='formatted', status='unknown',&
          position='append')
  endif

  if (kern_topo) then
     open(100, file='./kernel_topo.txt', form='formatted', status='unknown',&
          position='append')
  endif

  ! Initialize 1D Earth model and find dicontinuities
  call find_disc(file_model_1d, disc, rdisc, ndisc, NR)
  print*, disc, ndisc
  ! Get 1D Earth profile of density, and get its first spatial derivative
  call get_1d_rho(file_model_1d, rad_norm, rho_1d, NR, N_ICB, N_CMB, rdisc, ndisc)

  ! Get 1D Earth profile of density, and get its first spatial derivative
  call get_delta_rho(delta_rho, file_del_rho, layer_perturb_start, layer_perturb_end,&
       value_perturb, NR, kern_rho)

  ! Get 1D Earth profile of density, and get its first spatial derivative
  ! call get_delta_topography(delta_d, file_del_topo, disc_perturb, value_perturb, ndisc, kern_comp)
  call get_delta_topography(delta_d, file_del_topo, disc_perturb, value_perturb, ndisc, kern_topo)

  ! Compute delta V from eq. (15) from Buffett, 1996.
  prefactor = 4*PI*GRAV_CST/5
  ! print*, "Entering compute_rho_epsilon_average"
  call compute_rho_epsilon_average(rho_eps, rad_norm, rho_1d, delta_rho, delta_d, &
       l_chos, m_chos, ndisc, disc, n_icb, n_cmb, NR, kern_rho, kern_topo)
  ! print*, "Passed compute_rho_epsilon_average", rho_eps*RHO_AV

  call compute_epsilon(epsilon_surf, l_chos, m_chos, rad_norm, rho_1d, delta_rho, delta_d, &
       N_ICB, ndisc, disc, n_cmb, NR, kern_rho, kern_topo)
  ! print*, "Passed compute_epsilon for ICB", epsilon_surf*rho_1d(N_ICB)*RHO_AV, epsilon_surf

  call compute_equatorial_moment_of_inertia_IC(MOI_no_fluid, rho_1d, rad_norm, delta_rho, &
       delta_d, disc, ndisc, n_cmb, n_icb, NR, l_chos, m_chos, .false., kern_rho, kern_topo)
  ! print*, "Passed compute_equatorial_moment_of_inertia_IC", MOI_no_fluid*R_EARTH**5*RHO_AV

  call compute_equatorial_moment_of_inertia_IC(MOI_fluid, rho_1d, rad_norm, delta_rho, delta_d, &
       disc, ndisc, n_cmb, n_icb, NR, l_chos, m_chos, .true., kern_rho, kern_topo)
  ! print*, "Passed compute_equatorial_moment_of_inertia_IC for fluid", MOI_fluid*R_EARTH**5*RHO_AV
 
  delta_V_norm = prefactor * (abs(rho_eps+epsilon_surf)) * (abs(MOI_no_fluid-MOI_fluid))
  
  ! Re-dimensionalization
  ! ACC_NORM = (PI * GRAV_CST * RHO_AV * R_EARTH)
  ACC_NORM = R_EARTH**5 * RHO_AV**2
  delta_V_dim = ACC_NORM * delta_V_norm ! DIMENSION ?!?
  !
  sig_sq = (4*delta_V_dim/5.87d34) * (1.d0+5.87d34/7.12d37)
  sig = ((2*PI)/abs(sig_sq)**(0.5))/(86400*365)

  if (kern_rho) then
     call compute_epsilon(epsilon_CMB, l_chos, m_chos, rad_norm, rho_1d, delta_rho, delta_d, &
          N_CMB, ndisc, disc, n_cmb, NR, kern_rho, kern_topo)
     ! call compute_epsilon(epsilon_pert, l_chos, m_chos, rad_norm, rho_1d, delta_rho, delta_d, &
     !      layer_perturb, ndisc, disc, n_cmb, NR, kern_rho, kern_topo)
     epsilon_CMB = epsilon_surf*rho_1d(N_CMB)*RHO_AV
     ! epsilon_pert = epsilon_surf*rho_1d(layer_perturb)*RHO_AV
     ! write(100,'(I8,F8.2)') floor(rad_norm(layer_perturb)*R_EARTH/1000), sig
     write(100,'(I8,3E12.4)') floor(rad_norm(layer_perturb_start)*R_EARTH/1000), delta_V_dim, &
          epsilon_CMB
     write(100,'(I8,3E12.4)') floor(rad_norm(layer_perturb_end)*R_EARTH/1000), delta_V_dim, &
          epsilon_CMB

     close(100)
  elseif (kern_topo) then
     write(100,'(I8,3E12.4)') floor(rdisc(layer_perturb_start)/1000), delta_V_dim, &
          epsilon_CMB
     close(100)
  end if
  
  print*, "The period is: ", sig, " years"
  ! write(10,*) rad_norm(layer_perturb)*R_EARTH/1000, l_perturb, m_perturb, &
  !      value_perturb, sig
  !----------------------------------------
  deallocate(disc, rdisc)
  !
  stop  
  !
end program main
  
