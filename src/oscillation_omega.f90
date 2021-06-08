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
  real*8, allocatable :: phi_zero(:)
  real*8              :: phi

  ! Integration
  real*8 :: total_potential
  real*8, parameter :: PI = 3.1415927, R_EARTH = 6371.d3
    
  file_radius = '../../../data/make_s20rts/rho_ulm/rad.dat'
  file_model_1d = '../../../data/make_s20rts/isoprem808.md'

  ! Initialize 1D Earth model and find dicontinuities
  call find_disc(file_model_1d, disc, rdisc, ndisc, NR)

  ! Get 1D Earth profile of density, and get its first spatial derivative
  call get_1d_rho(file_model_1d, rad_norm, rho_1d, NR, rdisc, ndisc)
  call deriv(rho_1d,rho_prime,NR,rad_norm,ndisc,disc,s1,s2,s3)

  ! Compute gravitational potential in unperturbed Earth
  call compute_phi_zero(phi_zero, rho_1d, rad_norm, disc, ndisc, NR)
  
  ! Compute the integral of the gravitational potential
  call compute_volumetric_integral(total_potential, phi_zero, rad_norm, disc, ndisc, NR)  
  total_potential = total_potential * R_EARTH**3
  
  !----------------------------------------
  ! PRINTING SPACE FOR DEBUGGING THE MAIN
  print*, 'total potential', total_potential

  ! do i = 1,NR
  !    print*,rad_norm(i)*6371, phi_zero(i)! , phi(i)
  ! enddo
  !----------------------------------------
  ! deallocate(disc, rdisc)
  stop  
  
end program main
  
