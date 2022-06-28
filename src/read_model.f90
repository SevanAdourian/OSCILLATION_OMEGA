module read_model
  !
  ! subroutine find_disc(model_file, disc, rdisc, ndisc, NR)
  ! subroutine get_1d_rho(model_1d, rad_norm, rho_1d, NR, rdisc, ndisc)
  ! subroutine get_delta_rho(fname_rho_ylm_re, fname_rho_ylm_im, )
  implicit none
  !
contains
  
  subroutine get_1d_rho(model_1d, rad_norm, rho_1d, NR, n_icb, n_cmb, rdisc, ndisc)
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
    integer, intent(out) :: n_icb, n_cmb
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
    n_icb = int(dummy(2))+1
    n_cmb = int(dummy(3))+1

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
  !
  subroutine get_delta_rho_s20rts(delta_rho_mat, nlayer, lmax_model)
    !---------------------------------------
    !
    ! reads the rho perturbation due to 3D
    ! lateral heterogeneties
    !
    ! <SA> 06/2021 Returns a matrix for the entire layer instead of a specific line,
    ! to avoid I/O overload.
    ! </SA>
    ! <SA> 10/2021 Moved the naming of the files here for better readability and ensure
    ! that we access the files only once.
    ! </SA>
    implicit none
    !
    ! integer, intent(in) :: l, m
    integer, intent(in) :: lmax_model, nlayer

    ! complex*16, intent(out) :: delta_rho_mat(0:lmax_model,0:lmax_model)
    complex*16, allocatable, intent(out) :: delta_rho_mat(:,:)
    ! 
    character(len=100) :: fname_rho_ylm_im
    character(len=100) :: fname_rho_ylm_re
    ! 
    integer :: ii, nlayer_for_file
    real*8  :: mat_im(0:lmax_model,0:lmax_model), mat_re(0:lmax_model,0:lmax_model)

    nlayer_for_file = nlayer - 330
    allocate(delta_rho_mat(0:lmax_model,0:lmax_model))

    write(fname_rho_ylm_im,"(a,I3.3,a)")"../data/make_s20rts/rho_ulm/rho_ulm_im_lay",&
         nlayer_for_file,".dat"
    write(fname_rho_ylm_re,"(a,I3.3,a)")"../data/make_s20rts/rho_ulm/rho_ulm_re_lay",&
         nlayer_for_file,".dat"

    open(1,file=fname_rho_ylm_re,status='old',form='formatted',action='read')
    open(2,file=fname_rho_ylm_im,status='old',form='formatted',action='read')

    do ii = 0,lmax_model
       read(1,*) mat_re(ii,:)
       read(2,*) mat_im(ii,:)
    end do

    ! Here we correct for the fact that delta_rho(0,0) in the files is supposed to be 0
    ! but it is not because of some artifacts of S20RTS. Should be ok with any other model
    ! as physically, delta_rho(0,0) should always be 0 as it would change the mass of the
    ! Earth otherwise.
    mat_re(0,:) = 0
    mat_im(0,:) = 0
    delta_rho_mat = cmplx(mat_re,mat_im)

    close(1)
    close(2)

    return
  end subroutine get_delta_rho_s20rts

  subroutine get_delta_rho(delta_rho_mat, fname, layer_pert_start, layer_pert_end, &
       value_pert, NR, kernel_or_not)
    !---------------------------------------
    !
    ! reads the rho perturbation due to 3D
    ! lateral heterogeneties
    !
    ! <SA> 06/2021 Returns a matrix for the entire layer instead of a specific line,
    ! to avoid I/O overload.
    ! </SA>
    ! <SA> 10/2021 Moved the naming of the files here for better readability and ensure
    ! that we access the files only once.
    ! </SA>
    implicit none
    !
    character(len=100), intent(in) :: fname
    integer, intent(in) :: NR, layer_pert_start, layer_pert_end
    real*8, intent(in) :: value_pert
    logical, intent(in) :: kernel_or_not
    !
    complex*16, allocatable, intent(out) :: delta_rho_mat(:)
    ! 
    ! 
    integer :: ii, layer_mid
    real*8, allocatable  :: mat_im(:), mat_re(:)
    !
    allocate(delta_rho_mat(NR), mat_re(NR), mat_im(NR))

    ! Here we correct for the fact that delta_rho(0,0) in the files is supposed to be 0
    ! but it is not because of some artifacts of S20RTS. Should be ok with any other model
    ! as physically, delta_rho(0,0) should always be 0 as it would change the mass of the
    ! Earth otherwise.
    !
        ! if (kernel_or_not) then
    !    delta_rho_mat(:) = cmplx(0,0)
    !    ! delta_rho_mat(layer_pert_start:layer_pert_end) = cmplx(value_pert, 0)
    !    layer_mid = (layer_pert_start+layer_pert_end)/2
    !    print*, layer_mid
    !    ! delta_rho_mat(layer_mid-3:layer_mid+3) = cmplx(value_pert, 0)
    ! else
    !    open(1,file=fname,status='old',form='formatted',action='read')
    !    do ii = 1,NR
    !       read(1,*) mat_re(ii), mat_im(ii)
    !    end do
    ! endif

    open(1,file=fname,status='old',form='formatted',action='read')
    do ii = 1,NR
       read(1,*) mat_re(ii), mat_im(ii)
    end do
    !
    close(1)
    !
    delta_rho_mat = cmplx(mat_re,mat_im)
    !
    ! print*, layer_pert_start, layer_pert_end
    ! print*, delta_rho_mat(layer_pert_start:layer_pert_end)
    deallocate(mat_re)
    deallocate(mat_im)
    !
    return
  end subroutine get_delta_rho


  subroutine get_delta_topography(delta_d_mat, fname, disc_pert, value_pert, &
       ndisc, kernel_or_not)
    !---------------------------------------
    !
    ! reads the topography perturbation at the discontinuities specified in a 1D model
    !
    implicit none
    !
    character(len=100), intent(in) :: fname
    integer, intent(in) :: ndisc, disc_pert
    real*8, intent(in) :: value_pert
    logical, intent(in) :: kernel_or_not
    !
    complex*16, allocatable, intent(out) :: delta_d_mat(:)
    ! 
    integer :: ii
    real*8, allocatable  :: mat_im(:), mat_re(:)
    !
    allocate(delta_d_mat(ndisc), mat_re(ndisc), mat_im(ndisc))

    ! if (kernel_or_not) then
    !    delta_d_mat(:) = cmplx(0,0)
    !    delta_d_mat(disc_pert) = cmplx(value_pert, 0)
    ! else
    !    open(2,file=fname,status='old',form='formatted',action='read')
    !    do ii = 1,ndisc
    !       read(2,*) mat_re(ii), mat_im(ii)
    !    end do
    !    delta_d_mat = cmplx(mat_re, mat_im)
    ! endif
    
    open(2,file=fname,status='old',form='formatted',action='read')
    do ii = 1,ndisc
       read(2,*) mat_re(ii), mat_im(ii)
    end do
    delta_d_mat = cmplx(mat_re, mat_im)
    
    close(2)
    !
    deallocate(mat_re)
    deallocate(mat_im)
    !
    return
  end subroutine get_delta_topography

subroutine construct_rho_map(delta_rho_all_r, r_min, r_max, lmax_model)
    !
    implicit none
    !
    complex*16, intent(out), allocatable :: delta_rho_all_r(:,:,:)
    !
    integer, intent(in) :: lmax_model, r_min, r_max
    !
    ! complex*16 :: delta_rho(0:lmax_model,0:lmax_model)
    complex*16, allocatable :: delta_rho(:,:)
    integer :: nlayer
    !
    ! Allocating and initialize the output matrix
    allocate(delta_rho_all_r(r_min:r_max,0:lmax_model,0:lmax_model))
    delta_rho_all_r(:,:,:) = 0.d0
    !
    do nlayer = r_min, r_max
       call get_delta_rho_s20rts(delta_rho, nlayer, lmax_model)
       delta_rho_all_r(nlayer,:,:) = delta_rho
    end do 
    
    return
    !
  end subroutine construct_rho_map
  
end module read_model
