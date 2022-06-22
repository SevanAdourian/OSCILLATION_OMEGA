module gravitational_potential
  !
  use ylm_procedures
  use read_model
  use discontinuities
  !
  implicit none
  !
contains

  subroutine compute_epsilon(epsilon_l_m, l, m, rad, rho, del_rho, del_d, &
       ind_r, ndisc, disc, N_CMB, NR, kern_rho, kern_topo)
    !
    implicit none
    !
    complex*16, allocatable, intent(in) :: del_rho(:), del_d(:)
    real*8, allocatable, intent(in)  :: rad(:), rho(:)
    integer, allocatable, intent(in) :: disc(:)
    integer, intent(in) :: l,m,ind_r,ndisc, N_CMB, NR
    logical, intent(in) :: kern_rho, kern_topo
    !
    real*8, intent(out) :: epsilon_l_m
    !
    complex*16, allocatable :: delta_rho(:,:)
    real*8, allocatable :: integrd1_re(:), integrd1_im(:), integrd2_re(:), integrd2_im(:)
    real*8 :: diff_rho
    real*8 :: integrl1_re, integrl1_im, integrl2_re, integrl2_im
    real*8 :: r_nat, rho_layer
    real*8 :: sum_disc, sum_disc_ind, prefactor
    integer :: i, j 
    !
    allocate(integrd1_re(NR),integrd2_re(NR))
    allocate(integrd1_im(NR),integrd2_im(NR))
    integrd1_re(:) = 0.d0
    integrd2_re(:) = 0.d0
    integrd1_im(:) = 0.d0
    integrd2_im(:) = 0.d0
    sum_disc = 0.d0
    !
    r_nat = rad(ind_r)
    do i = 2,NR
       rho_layer = rho(i)
       if (i >= N_CMB) then
          ! if (i <= ind_r) then ! <= or < ?
             integrd1_re(i) = realpart(del_rho(i)) * rho_layer *&
                  ((rad(i)**(l+2))/(rad(ind_r)**(l+1)))
             integrd1_im(i) = imagpart(del_rho(i)) * rho_layer *&
                  ((rad(i)**(l+2))/(rad(ind_r)**(l+1)))
          ! else
             integrd2_re(i) = realpart(del_rho(i)) * rho_layer *&
                  ((rad(ind_r)**(l))/(rad(i)**(l-1)))
             integrd2_im(i) = imagpart(del_rho(i)) * rho_layer *&
                  ((rad(ind_r)**(l))/(rad(i)**(l-1)))
          ! end if
       end if
    end do
    !
    call intgrl_disc(integrl1_re, NR, rad, disc, ndisc, &
         1, ind_r, integrd1_re)
    call intgrl_disc(integrl2_re, NR, rad, disc, ndisc, &
         ind_r+1, NR, integrd2_re)
    call intgrl_disc(integrl1_im, NR, rad, disc, ndisc, &
         1, ind_r, integrd1_im)
    call intgrl_disc(integrl2_im, NR, rad, disc, ndisc, &
         ind_r+1, NR, integrd2_im)

    ! Summing discontinuities contributions
    do j = 2,ndisc-1
       diff_rho = rho(disc(j+1))-rho(disc(j))
       if (r_nat .lt. rad(disc(j))) then
          sum_disc_ind = diff_rho * (rad(disc(j))**l)/(r_nat**(l+1)) * del_d(j)
       else
          sum_disc_ind = diff_rho * (r_nat**l)/(rad(disc(i))**(l+1)) * del_d(j)
       end if
       sum_disc = sum_disc + sum_disc_ind
       print*, sum_disc_ind
    end do
    !
    ! Summing all contributions
    ! sum_disc = 0.d0
    prefactor = 1 / (rho(ind_r)*r_nat**l)
    if (kern_rho .and. .not. kern_topo) then
       epsilon_l_m = prefactor * (integrl1_re + integrl1_im + integrl2_re + integrl2_im)
    elseif (kern_topo .and. .not. kern_rho) then
       epsilon_l_m = prefactor * sum_disc
    elseif ( .not. kern_topo .and. .not. kern_rho) then
       epsilon_l_m = prefactor * (integrl1_re + integrl1_im + integrl2_re + integrl2_im + sum_disc)
    end if

    print*, epsilon_l_m
    !
    ! write(111,*) ind_r, epsilon_l_m, integrl1_re, integrl1_im, integrl2_re, integrl2_im, rho(ind_r), r_nat
    ! if (ind_r .eq. 400) close(112)
    deallocate(integrd1_re,integrd2_re,integrd1_im,integrd2_im)
    !
    return
    !
  end subroutine compute_epsilon

  subroutine compute_rho_epsilon_average(rho_epsilon_average, rad, rho, del_rho, del_d, &
       l, m, ndisc, disc, N_ICB, N_CMB, NR, kern_rho, kern_topo)
    !
    implicit none
    !
    complex*16, allocatable, intent(in) :: del_rho(:), del_d(:)
    real*8, allocatable, intent(in) :: rad(:), rho(:)
    integer, allocatable, intent(in) :: disc(:)
    integer, intent(in) :: ndisc, NR, l, m, N_CMB, N_ICB
    logical, intent(in) :: kern_rho, kern_topo
    !
    real*8, intent(out) :: rho_epsilon_average
    !
    real*8,  allocatable :: epsilon_a_all(:), epsilon_a_prime(:), integrd(:)
    real*8, allocatable :: s1(:), s2(:), s3(:)
    real*8 :: epsilon_a
    integer :: i
    !
    allocate(epsilon_a_all(NR), epsilon_a_prime(NR), integrd(NR))
    epsilon_a_all(:) = 0.d0
    integrd(:) = 0.d0

    ! Equation (14) - Buffett, 1996
    ! Gather epsilon from all radii
    ! Starting at 2 to avoid NaN at the center 

    ! do i = N_ICB,NR
    do i = 2,NR
       call compute_epsilon(epsilon_a, l, m, rad, rho, del_rho, del_d, i, &
            ndisc, disc, N_CMB, NR, kern_rho, kern_topo)
       epsilon_a_all(i) = epsilon_a
    end do
    !
    call deriv(epsilon_a_all, epsilon_a_prime, NR, rad, ndisc, disc, s1, s2, s3)
    !
    ! Gathering of the integrand in eq. (14)
    ! do i = N_ICB,NR
    do i = 2,NR
       integrd(i) = epsilon_a_prime(i)*rho(i)
    end do
    !
    ! Integration over the whole radius
    call intgrl_disc(rho_epsilon_average, NR, rad, disc, ndisc, &
         N_ICB, NR, integrd)
    !
    deallocate(epsilon_a_all, epsilon_a_prime, integrd)
    return
    !
  end subroutine compute_rho_epsilon_average
    
end module gravitational_potential
