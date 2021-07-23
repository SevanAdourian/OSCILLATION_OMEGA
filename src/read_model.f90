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
       ! NON-DIMENSIONALIZE RHO AS WELL
       ! rho_1d(i-j+1) = radius(2)
       rho_1d(i-j+1) = radius(2)/RHO_AV
    enddo

    close(iomod)
    return

  end subroutine get_1d_rho



  subroutine get_delta_rho(delta_rho_mat, fname_rho_ylm_re, fname_rho_ylm_im, lmax)
    !---------------------------------------
    !
    ! reads the rho perturbation due to 3D
    ! lateral heterogeneties
    !
    ! <SA> 06/2021 Returns a matrix for the entire layer instead of a specific line,
    ! to avoid I/O overload.
    ! </SA>
    implicit none
    !
    character, intent(in) :: fname_rho_ylm_im*100
    character, intent(in) :: fname_rho_ylm_re*100
    ! integer, intent(in) :: l, m
    integer, intent(in) :: lmax

    complex, intent(out) :: delta_rho_mat(0:lmax,0:lmax)

    integer :: ii
    real*8  :: mat_im(0:lmax,0:lmax), mat_re(0:lmax,0:lmax)

    open(1,file=fname_rho_ylm_re,status='old',form='formatted',action='read')
    open(2,file=fname_rho_ylm_im,status='old',form='formatted',action='read')

    do ii = 0,lmax
       read(1,*) mat_re(ii,:)
       read(2,*) mat_im(ii,:)
    end do

    delta_rho_mat = cmplx(mat_re,mat_im)
    
    ! do ii = 0,l-1
    !    read(1,*)
    !    read(2,*)
    ! enddo

    ! read(1,*) mat_re
    ! read(2,*) mat_im

    ! delta_rho = cmplx(mat_re(m), mat_im(m))

    close(1)
    close(2)

    return
  end subroutine get_delta_rho

  subroutine get_interpolate_rho(rho_value, rho_1d, rad_norm, lat, lon, depth)
    !
    ! Interpolates all the density for any lat lon depth in order to compute
    ! an integral on consecutive cylinders without depending on the spherical
    ! geometry of the model and the problem.
    !
    implicit none
    !
    ! Normalize depth
    r_norm = depth/R_EARTH
    if (r_norm > 1) then
       exit
    end if
    !
    ! Find the index layer

    ! 
    
  end subroutine get_interpolate_rho
    
end module read_model
