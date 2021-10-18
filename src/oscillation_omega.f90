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
  character :: file_radius*100
  ! character :: file_rho_ylm_im*100
  ! character :: file_rho_ylm_re*100
  character :: file_model_1d*100

  integer :: nlayer, l, m, i
  integer :: ndisc, NR

  real*8 :: lat_test, lon_test
  integer, allocatable :: disc(:)
  real*8,  allocatable  :: rdisc(:)

  real*8, allocatable :: rho_1d(:), rad_norm(:), rho_prime(:)
  real*8, allocatable :: s1(:), s2(:), s3(:)
  real*8, allocatable :: phi_zero(:), g_zero(:)
  real*8              :: phi

  ! DEBUG
  real*8, allocatable :: phi_cmb(:,:)
  
  ! Integration
  real*8 :: total_potential, moment
  real*8 :: ACC_NORM
  real*8, parameter :: PI = 3.1415927, R_EARTH = 6371.d3, RHO_AV = 5510.d0, GRAV_CST = 6.67408d-11
  integer :: lmax = 16
  integer :: lmax_model = 20

  ! Kernel
  real*8, allocatable :: kernel_grav(:,:)
  real*8 :: kernel_grav_r
  file_radius = '../../../../data/make_s20rts/rho_ulm/rad.dat'
  file_model_1d = '../../../../data/make_s20rts/isoprem808.md'

  open(unit=10, file='kernel.txt', ACTION="write", STATUS="replace")

  ! Initialize 1D Earth model and find dicontinuities
  call find_disc(file_model_1d, disc, rdisc, ndisc, NR)

  ! Get 1D Earth profile of density, and get its first spatial derivative
  call get_1d_rho(file_model_1d, rad_norm, rho_1d, NR, rdisc, ndisc)
  call deriv(rho_1d,rho_prime,NR,rad_norm,ndisc,disc,s1,s2,s3)

  ! Compute gravitational potential in unperturbed Earth
  call compute_phi_zero(phi_zero, g_zero, rho_1d, rad_norm, disc, ndisc, NR)

  ! Compute moment of inertia
  ! call compute_moment_of_inertia(moment, rho_1d, &
  !      rad_norm, 1, NR-26, lmax, lmax_model, disc, ndisc, NR)

  
  ! allocate(kernel_grav(NR,6))
  ! kernel_grav(1,:) = 0.d0
  ! do ii = 2,NR
  !    do l = 1,6
  !       call compute_kernel_grav(kernel_grav_r, rho_1d, rad_norm, ii, l, NR, ndisc, disc)
  !       kernel_grav(ii,l) = kernel_grav_r
  !    end do
  !    write(10,"(6E10.2)")kernel_grav(ii,:) 
  ! end do
  ! deallocate(kernel_grav)
  !
  ! DEBUG CHECK THE GRAVITATIONAL POTENTIAL
  ! allocate(phi_cmb(nx,ny))
  ! do ii = 1,nx
  !    lat_test = ((ii*dlat)-91)*pi/180
  !    do jj = 1,ny
  !       lon_test = (jj*dlon - 1)*pi/180
  !       ! call compute_gravitational_potential(phi, rho_1d, rho_prime, rad_norm, &
  !       !      nlayer, lat_test, lon_test, disc, ndisc, NR)
  !       call compute_delta_phi(phi, phi_zero, rho_1d, rad_norm, rho_prime, 20, nlayer, (pi/2)-lat_test, lon_test, NR)
  !       phi_cmb(ii,jj) = phi
  !       write(10,*)lat_test, lon_test, phi
  !       print*, lat_test, lon_test
  !    end do
  ! end do
  
  ! Compute the integral of the gravitational potential
  call compute_grav_pot_volumetric_integral(total_potential, rho_1d, rad_norm,&
       lmax, lmax_model, disc, ndisc, NR)
  ! TEMPORARY
  ! ACC_NORM = (PI * GRAV_CST * RHO_AV * R_EARTH)
  ! print*, R_EARTH**4 * RHO_AV * ACC_NORM*total_potential
  ! total_potential =  R_EARTH**4 * RHO_AV * ACC_NORM * &
  !      (4*total_potential/5.87d34) * (1.d0+5.87d34/7.12d37)
  total_potential =  R_EARTH**5 * RHO_AV**2 * &
       (4*total_potential/5.87d34) * (1.d0+5.87d34/7.12d37)

  print*, total_potential
  ! total_potential =  R_EARTH**5 * RHO_AV * moment

  ! print*, total_potential
  !----------------------------------------
  ! PRINTING SPACE FOR DEBUGGING THE MAIN

  ! do i = 1,NR
  !    print*,rad_norm(i)*6371, phi_zero(i)! , phi(i)
  ! enddo
  !----------------------------------------
  ! deallocate(disc, rdisc)
  ! deallocate(phi_cmb)
  close(10)
  stop  
  
end program main
  
