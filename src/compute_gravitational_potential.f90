module gravitational_potential
  !
  use ylm_procedures
  use read_model
  use discontinuities
  !
  implicit none
  !
contains

  subroutine compute_gravitational_potential(del_phi, rho, rho_prime, rad_norm,&
       g_zero,nlayer, lat, lon, disc, ndisc, NR)
    !
    implicit none
    !
    real*8, allocatable, intent(in)  :: rho(:), rho_prime(:), rad_norm(:)
    real*8, intent(in)               :: lat, lon
    integer, allocatable, intent(in) :: disc(:)
    integer, intent(in)              :: NR, ndisc, nlayer
    !
    real*8, intent(out) :: del_phi
    !
    real*8, allocatable :: phi_zero(:), g_zero(:)
    
    integer :: lmax = 20
    ! call compute_phi_zero(phi_zero, g_zero, rho, rad_norm, disc, ndisc, NR)
    ! print*,"[compute_gravitational_potential] phi_zero", rho(nlayer),rad_norm(nlayer), phi_zero(nlayer)
    ! Compute the gravitational potential at lat/lon
    call compute_delta_phi(del_phi, g_zero, rho, rad_norm, rho_prime, lmax, nlayer, lat, lon, NR)
    ! print*, phi
  end subroutine compute_gravitational_potential

  subroutine compute_phi_zero(phi_zero, g_zero, rho, r, disc, ndisc, NR)

    implicit none
    !
    real*8, allocatable, intent(in)  :: rho(:), r(:)
    integer, intent(in) :: NR, ndisc
    integer, allocatable, intent(in) :: disc(:)
    !
    real*8, allocatable, intent(out) :: phi_zero(:), g_zero(:)
    !
    real*8, parameter :: R_EARTH = 6371.d3, RHO_AV = 5510.d0
    real*8, parameter :: GRAV_CST = 6.67408d-11
    real*8, parameter :: PI = 3.1415927
    real*8, allocatable :: integrd_1(:), integrd_2(:)
    real*8 :: integrl_1, integrl_2, cste
    integer :: i

    allocate(phi_zero(NR))
    allocate(g_zero(NR))
    allocate(integrd_1(NR), integrd_2(NR))
    !
    integrd_1(:) = rho(:) * r(:)*r(:)
    integrd_2(:) = rho(:) * r(:)

    phi_zero(:) = 0.d0

    open(11,file="g_zero.txt",form="formatted")
    do i = 1,NR
       call intgrl_disc(integrl_1,NR,r(1:NR), disc,ndisc,1,i,integrd_1(1:NR))
       call intgrl_disc(integrl_2,NR,r(1:NR),disc,ndisc,i,NR,integrd_2(1:NR))
       cste = (-4 *  PI * GRAV_CST) 
       ! phi_zero(i) = cste * ( (integrl_1 / (r(i))) + integrl_2) * R_EARTH * R_EARTH
       phi_zero(i) = cste * ( (integrl_1 / (r(i))) + integrl_2) !* R_EARTH * R_EARTH
       g_zero(i) = -cste*integrl_1 / (r(i)*r(i)) / (PI * GRAV_CST)
       write(11,*)g_zero(i) * (PI * GRAV_CST * RHO_AV * R_EARTH)
    end do

    deallocate(integrd_1)
    deallocate(integrd_2)
    close(11)
    return
  end subroutine compute_phi_zero

  subroutine compute_delta_phi(phi, rho_norm, rad_norm, rho_prime, g_zero, &
       lmax, nlayer, lat, lon, NR)
    !
    implicit none
    !
    real*8, allocatable, intent(in) :: rho_prime(:), rad_norm(:), rho_norm(:), g_zero(:)
    ! real*8, allocatable, intent(in) :: phi_zero(:)
    real*8, intent(in) :: lat, lon
    integer, intent(in) :: nlayer, lmax, NR
    ! 
    real*8, intent(out) :: phi
    ! 
    integer :: i, l, m
    complex*8 :: chi_lm, chi_lm_ylm_sum, chi_lm_ylm, chi_lm_ylm_m
    complex*8 :: delta_rho(0:lmax,0:lmax)

    ! TEMPORARY, BEFORE WE HAVE A SETUP FILE WHERE TO PUT THE ROOT OF THOSE NAMES
    character :: file_rho_ylm_im*100
    character :: file_rho_ylm_re*100
    ! allocate(phi(NR))

    ! form the name of files to scan for density perturbations
    write(file_rho_ylm_im,"(a,I3.3,a)")"../../../../data/make_s20rts/rho_ulm/rho_ulm_im_lay",&
         nlayer,".dat"
    write(file_rho_ylm_re,"(a,I3.3,a)")"../../../../data/make_s20rts/rho_ulm/rho_ulm_re_lay",&
         nlayer,".dat"

    ! Get the density perturbation
    call get_delta_rho(delta_rho, file_rho_ylm_re, file_rho_ylm_im, lmax)

    ! Start loop on l and m
    chi_lm_ylm_sum = 0
    do l = 1,lmax
       chi_lm_ylm_m = 0
       do m = 0,l
          ! Compute the
          call compute_chi(chi_lm, rho_norm(nlayer), rho_prime(nlayer), rad_norm(nlayer), &
               g_zero(NR), delta_rho(l,m), NR)
          !
          ! MAYBE MOVE THIS TO A SEPARATE SUBROUTINE IN ORDER TO USE IT TO EXPAND RHO AS WELL
          ! Expand in spatial domain from spherical harmonics
          call ylm(lat, lon, l, m, chi_lm_ylm)
          ! sum chi lm over m
          if (m == 0) then
             chi_lm_ylm_m = chi_lm_ylm_m + realpart(chi_lm_ylm * chi_lm)
          else
             chi_lm_ylm_m = chi_lm_ylm_m + 2*(realpart(chi_lm_ylm * chi_lm))
          end if
       end do
       chi_lm_ylm_sum = chi_lm_ylm_sum + chi_lm_ylm_m
    end do
    ! phi = phi_zero(nlayer) + chi_lm_ylm_sum
    phi = chi_lm_ylm_sum
    ! <SA> DEBUG
    ! print*, "[compute_delta_phi]", phi
    return

  end subroutine compute_delta_phi

  subroutine compute_chi(chi_lm, rho_of_ind, rho_prime_of_ind, r_of_ind, g_zero_of_ind,&
       del_rho_l_m, NR)

    implicit none
    !
    real*8, intent(in)  :: rho_prime_of_ind, r_of_ind, rho_of_ind, g_zero_of_ind
    complex*8, intent(in)  :: del_rho_l_m
    integer, intent(in) :: NR
    !
    complex, intent(out) :: chi_lm
    !
    real*8 :: absolute_rho_ND
    !
    ! Make it absolute perturbation, add rho as an input argument of the subroutine
    absolute_rho_ND = (del_rho_l_m * rho_of_ind)
    ! chi_lm = (absolute_rho_ND/(rho_prime_of_ind)*(r_of_ind))
    chi_lm = (absolute_rho_ND/(rho_prime_of_ind))*g_zero_of_ind

    return
  end subroutine compute_chi

end module gravitational_potential
