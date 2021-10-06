module gravitational_potential
  !
  use ylm_procedures
  use read_model
  use discontinuities
  !
  implicit none
  !
contains

  ! subroutine compute_gravitational_potential(del_phi, rho, rho_prime, rad_norm,&
  !      g_zero,nlayer, lat, lon, disc, ndisc, NR)
  !   !
  !   implicit none
  !   !
  !   real*8, allocatable, intent(in)  :: rho(:), rho_prime(:), rad_norm(:)
  !   real*8, intent(in)               :: lat, lon
  !   integer, allocatable, intent(in) :: disc(:)
  !   integer, intent(in)              :: NR, ndisc, nlayer
  !   !
  !   real*8, intent(out) :: del_phi
  !   !
  !   real*8, allocatable :: phi_zero(:), g_zero(:)
    
  !   ! integer :: lmax = 20
  !   ! call compute_phi_zero(phi_zero, g_zero, rho, rad_norm, disc, ndisc, NR)
  !   ! print*,"[compute_gravitational_potential] phi_zero", rho(nlayer),rad_norm(nlayer), phi_zero(nlayer)
  !   ! Compute the gravitational potential at lat/lon
  !   ! call compute_delta_phi(del_phi, g_zero, rho, rad_norm, rho_prime, lmax, nlayer, lat, lon, NR)
  !   ! print*, phi
  ! end subroutine compute_gravitational_potential

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

    ! open(11,file="g_zero.txt",form="formatted")
    do i = 1,NR
       call intgrl_disc(integrl_1,NR,r(1:NR), disc,ndisc,1,i,integrd_1(1:NR))
       call intgrl_disc(integrl_2,NR,r(1:NR),disc,ndisc,i,NR,integrd_2(1:NR))
       cste = (-4 *  PI * GRAV_CST) 
       ! phi_zero(i) = cste * ( (integrl_1 / (r(i))) + integrl_2) * R_EARTH * R_EARTH
       phi_zero(i) = cste * ( (integrl_1 / (r(i))) + integrl_2) !* R_EARTH * R_EARTH
       g_zero(i) = cste*integrl_1 / (r(i)*r(i)) / (PI * GRAV_CST)
       ! write(11,*)g_zero(i) * (PI * GRAV_CST * RHO_AV * R_EARTH)
    end do

    deallocate(integrd_1)
    deallocate(integrd_2)
    ! close(11)
    return
  end subroutine compute_phi_zero

  subroutine compute_delta_phi(del_phi_sum, rho_abs_lm, rho_norm, rad_norm, &
       delta_rho, kern_grav_nlayer, lmax, nlayer, lat, lon, NR)
    !
    implicit none
    !
    real*8, allocatable, intent(in) :: rad_norm(:), rho_norm(:)
    real*8, allocatable, intent(in) :: kern_grav_nlayer(:)
    real*8, intent(in) :: lat, lon
    complex*8, intent(in) :: delta_rho(:,:)
    integer, intent(in) :: nlayer, lmax, NR
    ! 
    real*8, intent(out) :: del_phi_sum
    real*8, intent(out) :: rho_abs_lm
    ! 
    integer :: i, l, m
    ! integer :: nlayer_for_file
    ! complex*8 :: chi_lm, chi_lm_ylm_sum, chi_lm_ylm, chi_lm_ylm_m
    complex*8 :: del_phi_lm, del_phi_m
    complex*8 :: y_lm

    ! TEMPORARY, BEFORE WE HAVE A SETUP FILE WHERE TO PUT THE ROOT OF THOSE NAMES
    ! character :: file_rho_ylm_im*100
    ! character :: file_rho_ylm_re*100
    real*8, parameter :: GRAV_CST = 6.67408d-11

    ! allocate(phi(NR))

    ! <SA> Moved this to the get_delta_rho subroutine to avoid many useless computations.
    ! form the name of files to scan for density perturbations
    ! nlayer_for_file = nlayer - 330
    ! write(file_rho_ylm_im,"(a,I3.3,a)")"../../../../data/make_s20rts/rho_ulm/rho_ulm_im_lay",&
    !      nlayer_for_file,".dat"
    ! write(file_rho_ylm_re,"(a,I3.3,a)")"../../../../data/make_s20rts/rho_ulm/rho_ulm_re_lay",&
    !      nlayer_for_file,".dat"

    ! Get the density perturbation
    ! <SA> 10/2021 Moved this to the integral computation, in order to read it only once
    ! per layer.

    ! Start loop on l and m
    del_phi_sum = 0
    do l = 0,lmax
       del_phi_m = 0
       ! No need to do it here, call it directly from the volumetric integration routine
       ! so we don't do the same computations for each layers.
       ! call compute_kernel_grav()
       do m = 0,l
          ! MAYBE MOVE THIS TO A SEPARATE SUBROUTINE IN ORDER TO USE IT TO EXPAND RHO AS WELL
          ! Expand in spatial domain from spherical harmonics
          call ylm(lat, lon, l, m, y_lm)
          ! Computing here absolute value for rho for a given l and m, used to compute
          ! the integral.
          rho_abs_lm = delta_rho(l,m) * rho_norm(nlayer) + rho_norm(nlayer)
          del_phi_lm = delta_rho(l,m) * rho_norm(nlayer) * kern_grav_nlayer(l)
          ! sum chi lm over m
          if (m == 0) then
             del_phi_m = del_phi_m + realpart(del_phi_lm * y_lm)
          else
             del_phi_m = del_phi_m + 2 *(realpart(del_phi_lm * y_lm))
          end if
       end do
       del_phi_sum = del_phi_sum + del_phi_m
    end do
    del_phi_sum = -1 * GRAV_CST * del_phi_sum 
    ! phi = phi_zero(nlayer) + chi_lm_ylm_sum
    ! phi = chi_lm_ylm_sum
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

  subroutine compute_kernel_grav(k_g_lm, rho, r, lmax, NR, ndisc, disc)
    !
    implicit none
    !
    real*8, allocatable, intent(in)  :: rho(:), r(:)
    integer, intent(in) :: NR, ndisc, lmax
    integer, allocatable, intent(in) :: disc(:)
    !
    real*8, allocatable, intent(out) :: k_g_lm(:,:)
    !
    real*8, parameter :: R_EARTH = 6371.d3, RHO_AV = 5510.d0
    real*8, parameter :: GRAV_CST = 6.67408d-11
    real*8, parameter :: PI = 3.1415927
    real*8, allocatable :: integrd_1(:), R_rl(:)
    real*8 :: integrl_1, k_g_lm_ind
    integer :: ii, l
    !
    allocate(R_rl(NR))
    allocate(integrd_1(NR))
    !
    allocate(k_g_lm(1:NR,0:lmax))

    open(11,file="kernel.txt", form="formatted",status="replace")
    do ii = 2,NR
       do l = 0,lmax
          call compute_one_over_r_sph_harm(R_rl, r, ii, l, NR)
          integrd_1(:) = rho(:) * R_rl(:)
          integrd_1(1) = 0

          ! Integration over the radius
          call intgrl_disc(k_g_lm_ind, NR,r(1:NR), disc,ndisc,1,ii,integrd_1(1:NR))
          k_g_lm(ii,l) = k_g_lm_ind
          
          ! open(11,file="k_g_lm.txt",form="formatted")

          ! do i= 1,NR
          ! end do
          ! do i = 1,NR
          ! write(11,*)k_g_lm ! * (PI * GRAV_CST * RHO_AV * R_EARTH)
          ! end do
       end do
       write(11,"(21F10.2)") k_g_lm(ii,:)
    end do
    
    deallocate(R_rl)
    deallocate(integrd_1)
    
    ! close(11)
    return
  end subroutine compute_kernel_grav

  subroutine compute_one_over_r_sph_harm(one_over_r, r, ind_r, l, NR)
    !
    implicit none
    !
    real*8, allocatable, intent(in) :: r(:)
    integer, intent(in) :: ind_r, l, NR
    !
    real*8, allocatable, intent(out) :: one_over_r(:)
    !
    real*8 :: ratio, const
    integer :: i
    real*8, parameter :: PI = 3.1415927
    !
    allocate(one_over_r(NR))

    const = (4*PI)/(2*l+1)
    do i=1,NR
       if (i < ind_r) then
          ratio = (r(i)/r(ind_r))**l
          one_over_r(i) = (1.d0/r(i)) * const * ratio
       else
          ratio = (r(ind_r)/r(i))**l
          one_over_r(i) = (1.d0/r(ind_r)) * const * ratio
       end if
    end do

    return
  end subroutine compute_one_over_r_sph_harm
  
end module gravitational_potential
