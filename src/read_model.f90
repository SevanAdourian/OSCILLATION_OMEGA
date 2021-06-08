module read_model
  !
  ! subroutine find_disc(model_file, disc, rdisc, ndisc, NR)
  ! subroutine get_1d_rho(model_1d, rad_norm, rho_1d, NR, rdisc, ndisc)
  ! subroutine get_delta_rho(fname_rho_ylm_re, fname_rho_ylm_im, )
  implicit none
  !
contains
  
  subroutine get_1d_rho(model_1d, rad_norm, rho_1d, NR, rdisc, ndisc)
    !-----------------------------------------
    !
    ! reads the rho profile from the 1D file
    !
    implicit none
    !
    character, intent(in) :: model_1d*100
    integer, intent(in) :: NR, ndisc
    real*8, allocatable, intent(in) :: rdisc(:)
    double precision, allocatable, intent(out) :: rho_1d(:), rad_norm(:)
    !
    integer :: iomod=11, i, j
    integer:: nb_layers
    character :: title*100
    double precision :: dummy(3), radius(9), rho
    double precision, parameter :: R_EARTH = 6371.d3, RHO_AV = 5510.d0
    !
    open(iomod,file=model_1d,status='old',form='formatted',action='read')
    read(iomod,*) title
    read(iomod,*) dummy
    read(iomod,*) dummy

    nb_layers = int(dummy(1))
    print*, '[get_1d_rho]', nb_layers
    j = 1
    allocate(rho_1d(nb_layers))
    allocate(rad_norm(nb_layers))

    do i = 1, NR
       read(iomod,*) radius(:) !HL
       rad_norm(i-j+1) = radius(1)/R_EARTH
       rho_1d(i-j+1) = radius(2)
    enddo

    close(iomod)
    return

  end subroutine get_1d_rho



  subroutine get_delta_rho(fname_rho_ylm_re, fname_rho_ylm_im, &
       l, m, delta_rho)
    !---------------------------------------
    !
    ! reads the rho perturbation due to 3D
    ! lateral heterogeneties
    !
    implicit none
    !
    character, intent(in) :: fname_rho_ylm_im*100
    character, intent(in) :: fname_rho_ylm_re*100
    integer, intent(in) :: l, m

    complex, intent(out) :: delta_rho

    integer :: ii
    real*8  :: mat_im(0:19), mat_re(0:19)

    open(1,file=fname_rho_ylm_re,status='old',form='formatted',action='read')
    open(2,file=fname_rho_ylm_im,status='old',form='formatted',action='read')

    do ii = 0,l-1
       read(1,*)
       read(2,*)
    enddo

    read(1,*) mat_re
    read(2,*) mat_im

    delta_rho = cmplx(mat_re(m), mat_im(m))

    close(1)
    close(2)

    return
  end subroutine get_delta_rho

  ! subroutine find_disc(model_file, disc, rdisc, ndisc, NR)

  !   implicit none

  !   character, intent(in) :: model_file*100
  !   integer, intent(out) :: ndisc, NR
  !   integer, allocatable, intent(out) :: disc(:)
  !   real*8, allocatable, intent(out) :: rdisc(:)
  !   !
  !   real*8 :: radius_last,radius,eps
  !   integer :: jdisc,iomod,i
  !   character*80 :: junk ! HL
  !   double precision :: dummy(3)

  !   ! HL changed model format
  !   ! no need for '.dk' files, just '.md' files

  !   radius_last = -1.0
  !   eps = 1.0d-8
  !   jdisc = 0
  !   iomod = 11
  !   open(iomod,file=model_file,status='old')
  !   ! First pass to determine the number of discontinuities
  !   read(iomod,*) junk !HL
  !   read(iomod,*) dummy !HL, SA
  !   read(iomod,*) dummy !HL, SA

  !   NR = int(dummy(1))

  !   do i = 1, NR
  !      read(iomod,*) radius !HL
  !      if (abs(radius-radius_last) .lt. eps) then
  !         jdisc = jdisc + 1
  !      end if
  !      radius_last = radius
  !   enddo
  !   ndisc = jdisc

  !   allocate(disc(ndisc))
  !   allocate(rdisc(ndisc))

  !   ! 2nd pass to fill in the radii of the discontinuities
  !   jdisc = 0
  !   radius_last = -1.0

  !   rewind(iomod)
  !   read(iomod,*) junk !HL
  !   read(iomod,*) junk !HL
  !   read(iomod,*) junk !HL
  !   do i = 1, NR
  !      read(iomod,*) radius !HL
  !      if (abs(radius-radius_last) .lt. eps) then
  !         jdisc = jdisc + 1
  !         rdisc(jdisc) = radius
  !         disc(jdisc) = i - 1
  !      endif
  !      radius_last = radius
  !   enddo
  !   close(iomod)

  ! end subroutine find_disc
end module read_model
