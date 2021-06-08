module gravitational_potential
  !
  use ylm_procedures
  use read_model
  use discontinuities
  !
  implicit none
  !
contains

  ! subroutine compute_gravitational_potential(phi, rho, r, &
  !      nlayer, lat, lon, disc, ndisc, NR)
  !   !
  !   implicit none
  !   !
  !   real*8, allocatable, intent(in)  :: rho(:), r(:)
  !   integer, allocatable, intent(in) :: disc(:)
  !   integer, intent(in)              :: NR, ndisc
  !   !
  !   real*8, intent(out) :: phi


  subroutine compute_phi_zero(phi_zero, rho, r, disc, ndisc, NR)

    implicit none
    !
    real*8, allocatable, intent(in)  :: rho(:), r(:)
    integer, intent(in) :: NR, ndisc
    integer, allocatable, intent(in) :: disc(:)
    !
    real*8, allocatable, intent(out) :: phi_zero(:)
    !
    real*8, parameter :: R_EARTH = 6371.d3
    real*8, parameter :: GRAV_CST = 6.67408d-11
    real*8, parameter :: PI = 3.1415927
    real*8, allocatable :: integrd_1(:), integrd_2(:)
    real*8 :: integrl_1, integrl_2, cste
    integer :: i

    allocate(phi_zero(NR))
    allocate(integrd_1(NR), integrd_2(NR))
    !
    integrd_1(:) = rho(:) * r(:)*r(:)
    integrd_2(:) = rho(:) * r(:)

    phi_zero(:) = 0.d0

    do i = 1,NR
       call intgrl_disc(integrl_1,NR,r(1:NR), disc,ndisc,1,i,integrd_1(1:NR))
       call intgrl_disc(integrl_2,NR,r(1:NR),disc,ndisc,i,NR,integrd_2(1:NR))
       cste = (-4 *  PI * GRAV_CST) 
       phi_zero(i) = cste * ( (integrl_1 / (r(i))) + integrl_2) * R_EARTH * R_EARTH 
    end do

    deallocate(integrd_1)
    deallocate(integrd_2)

    return
  end subroutine compute_phi_zero

  subroutine compute_delta_phi(phi, phi_zero, rad_norm, rho_prime, nlayer, lat, lon, NR)
    !
    implicit none
    !
    real*8, allocatable, intent(in) :: rho_prime(:), rad_norm(:)
    real*8, allocatable, intent(in) :: phi_zero(:)
    real*8, intent(in) :: lat, lon
    integer, intent(in) :: nlayer, NR
    ! 
    real*8, intent(out) :: phi
    ! 
    integer :: i, l, m
    complex*8 :: delta_rho, chi_lm, chi_lm_ylm_sum, chi_lm_ylm, chi_lm_ylm_m

    ! TEMPORARY, BEFORE WE HAVE A SETUP FILE WHERE TO PUT THE ROOT OF THOSE NAMES
    character :: file_rho_ylm_im*100
    character :: file_rho_ylm_re*100
    ! allocate(phi(NR))

    ! form the name of files to scan for density perturbations
    write(file_rho_ylm_im,"(a,I3.3,a)")"../../data/make_s20rts/rho_ulm/rho_ulm_im_lay",&
         nlayer,".dat"
    write(file_rho_ylm_re,"(a,I3.3,a)")"../../data/make_s20rts/rho_ulm/rho_ulm_re_lay",&
         nlayer,".dat"

    ! Start loop on l and m
    chi_lm_ylm_sum = 0
    do l = 1,6
       chi_lm_ylm_m = 0
       do m = 1,l
          ! Get the density perturbation
          call get_delta_rho(file_rho_ylm_re, file_rho_ylm_im, &
               l, m, delta_rho)
          ! Compute the
          call compute_chi(chi_lm, rho_prime(nlayer), rad_norm(nlayer), delta_rho, NR)
          ! Expand in spatial domain from spherical harmonics
          call ylm(lat, lon, l, m, chi_lm_ylm)
          ! sum chi lm over m
          ! 2* real part except m = 0
          if (m == 0) then
             chi_lm_ylm_m = chi_lm_ylm_m + realpart(chi_lm_ylm * chi_lm)
          else
             chi_lm_ylm_m = chi_lm_ylm_m + 2*(realpart(chi_lm_ylm * chi_lm))
          end if
       end do
       chi_lm_ylm_sum = chi_lm_ylm_sum + chi_lm_ylm_m
    end do
    phi = phi_zero(nlayer) + chi_lm_ylm_sum 

    return

  end subroutine compute_delta_phi

  subroutine compute_chi(chi_lm, rho_prime_of_ind, r_of_ind, del_rho_l_m, NR)

    implicit none
    !
    real*8, intent(in)  :: rho_prime_of_ind, r_of_ind
    complex*8, intent(in)  :: del_rho_l_m
    integer, intent(in) :: NR
    !
    complex, intent(out) :: chi_lm
    !
    double precision, parameter :: R_EARTH = 6371.d3

    chi_lm = del_rho_l_m/(rho_prime_of_ind*(r_of_ind*R_EARTH))

    return
  end subroutine compute_chi

end module gravitational_potential
