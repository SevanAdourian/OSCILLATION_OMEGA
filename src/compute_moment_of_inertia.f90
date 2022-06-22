module moment_of_inertia
  !
  use discontinuities
  use gravitational_potential
  !
  implicit none
  !
contains
  !
  ! eq. (12), Buffett, 1996
  subroutine compute_equatorial_moment_of_inertia_IC(moment_of_inertia, rho_norm, rad, &
       delta_rho, delta_d, disc, ndisc, N_CMB, N_ICB, NR, l, m, flag_fluid, kern_rho, kern_topo)
    !
    implicit none
    !
    complex*16, allocatable, intent(in) :: delta_rho(:), delta_d(:)
    real*8, allocatable, intent(in) :: rho_norm(:), rad(:)
    integer, allocatable, intent(in) :: disc(:)
    integer, intent(in) :: ndisc, NR, N_ICB, N_CMB, l, m
    logical, intent(in) :: flag_fluid, kern_rho, kern_topo
    !
    real*8, intent(out) :: moment_of_inertia
    !
    real*8, allocatable :: epsilon_a_all(:), epsilon_a_prime(:), integrd(:)
    real*8, allocatable :: s1(:), s2(:), s3(:)
    real*8 ::epsilon_a, rho_fluid_core, prefactor
    integer :: i

    real*8, parameter :: PI = 3.1415926535
    !
    allocate(integrd(NR), epsilon_a_all(NR), epsilon_a_prime(NR))
    integrd(:)         = 0.d0
    epsilon_a_all(:)   = 0.d0
    epsilon_a_prime(:) = 0.d0
    
    ! DEBUG
    ! open(1,file='epsilon_a.txt')
    rho_fluid_core = rho_norm(N_ICB) ! N_ICB? N_ICB+1?
    do i = 2,NR
       call compute_epsilon(epsilon_a, l, m, rad, rho_norm, delta_rho, delta_d, &
            i, ndisc, disc, N_CMB, NR, kern_rho, kern_topo)
       epsilon_a_all(i) = epsilon_a * rad(i)**5
    end do
    !
    call deriv(epsilon_a_all, epsilon_a_prime, NR, rad, ndisc, disc, s1, s2, s3)
    !
    ! Element by element multiplication, depending on fluid or not
    if (flag_fluid .eqv. .true.) then
       integrd(:) = epsilon_a_prime(:) * rho_fluid_core
    else
       integrd(:) = epsilon_a_prime(:) * rho_norm(:)
    end if
    !
    ! DEBUG
    ! do i = 1,N_ICB-1
    !    write(1,*) i, N_ICB, epsilon_a_all(i), epsilon_a_prime(i), rho_norm(i), integrd(i)
    ! enddo
    !
    ! close(1)
    !
    call intgrl_disc(moment_of_inertia, NR, rad, disc, ndisc, &
         1, N_ICB-1, integrd)
    !
    prefactor = 8.d0*PI/15.d0
    moment_of_inertia = prefactor * moment_of_inertia
    return
    !
  end subroutine compute_equatorial_moment_of_inertia_IC
  
end module moment_of_inertia
