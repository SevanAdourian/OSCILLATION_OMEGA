module volumetric_integral
  !
  use legendre_quadrature
  use discontinuities
  use gravitational_potential
  !
  implicit none
  !
contains

  subroutine compute_volumetric_integral_gravitational_energy(grav_potential, rho_norm, &
       rad, lmax, lmax_model,disc, ndisc, NR)

    use legendre_quadrature
    implicit none
    !
    integer, intent(in) :: NR, ndisc, lmax, lmax_model
    integer, allocatable, intent(in) :: disc(:)
    real*8, allocatable, intent(in) :: rad(:), rho_norm(:)
    !
    ! real*8, intent(out) :: integral_value
    real*8, intent(out) :: grav_potential
    !
    integer :: kind, order
    real*8, allocatable :: abs_int(:), w_int(:)
    real*8 :: a, b, alpha, beta
    complex*16 :: delta_rho(0:lmax_model,0:lmax_model)
    complex*16, allocatable :: delta_rho_all_r(:,:,:) !!!
    ! 
    integer :: i, j, nlayer, nlayer_r_prime, l, m
    real*8, allocatable :: integral_lon(:)
    real*8 :: rho_abs_lm
    real*8 :: colat, lon
    real*8 :: integral_lat
    complex*16 :: rho_abs, y_lm
    real*8, allocatable :: r_l(:)
    complex*16 :: integral_r, integral_r_r_prime_lm
    complex*16 :: del_rho_lm
    complex*16 :: integral_r_r_prime_l
    complex*16, allocatable :: integrd_r(:), integrd_lm_r_prime(:)
    real*8, allocatable :: integrd_r_re(:), integrd_r_im(:)
    real*8 :: integral_r_re, integral_r_im
    real*8 :: integral_r_r_prime_lm_re, integral_r_r_prime_lm_im
    
    real*8, parameter :: PI = 3.1415927, GRAV_CST = 6.67408d-11

    ! Parameters for the interation rule (move to arg of subroutine eventually)
    ! order = 3
    ! kind = 1
    ! alpha = 0.d0
    ! beta = 0.d0
    ! a = -1
    ! b = 1

    ! Allocation
    ! allocate(integral_lon(NR))
    ! allocate(abs_int(order), w_int(order))
    allocate(r_l(NR))
    allocate(integrd_r(NR), integrd_lm_r_prime(NR))
    allocate(integrd_r_re(NR), integrd_r_im(NR))

    r_l(:) = 0.d0
    integrd_r(:) = 0.d0
    integrd_lm_r_prime(:) = 0.d0
    integrd_r_re(:) = 0.d0
    integrd_r_im(:) = 0.d0
    grav_potential = cmplx(0.d0,0.d0)

    !
    call construct_rho_map(delta_rho_all_r, 331, 782, lmax_model)
    ! Major loop in l and m, we moved the summation outside the integral computation
    integral_r_r_prime_l = cmplx(0.d0,0.d0)
    do l = 2,2
       print*, "Computing l value ", l,"/",lmax
       integral_r_r_prime_lm = cmplx(0.d0,0.d0)
       ! do m = -l,l
       do m = -2,2,4
       print*, "m = ", m
          ! Second and final radial integration
          !do nlayer_r_prime = 331,332
          do nlayer_r_prime = 2,782
             ! Initialization for a given layer
             ! First volume integration, in r' and omega
             ! Starting at 2 to avoid the singularity at the center, ending right before the crust
             ! do nlayer = 331, 332
             do nlayer = 2, 780
                ! get the density perturbation, make a map of rho
                if (nlayer.ge. 331 .and. nlayer .lt. 782)then
                   ! call get_delta_rho(delta_rho, nlayer, lmax_model)
                   delta_rho(:,:) = delta_rho_all_r(nlayer,:,:)
                else
                   delta_rho(:,:) = cmplx(0.d0,0.d0)
                end if
                ! print*, "After delta rho"
                   
                ! No integration in lat/lon
                call compute_one_over_r_sph_harm(r_l, rad, nlayer, l, NR)
                ! print*, "After one over r"               
                if (l .eq. 0) del_rho_lm = cmplx(rho_norm(nlayer),0)
                if (m < 0) then
                   del_rho_lm = (-1.d0)**m * conjg(delta_rho(l,abs(m)))
                else
                   del_rho_lm = delta_rho(l,abs(m))
                end if
                
                integrd_r(nlayer) = r_l(nlayer) * conjg(del_rho_lm) * rad(nlayer)**2
                
             end do ! nlayer, r'         
             ! Radial integration
             integrd_r_re = realpart(integrd_r)
             integrd_r_im = dimag(integrd_r)
             call intgrl_disc(integral_r_re, NR,rad(1:NR), disc,ndisc,&
                  1,NR,integrd_r_re)
             call intgrl_disc(integral_r_im, NR,rad(1:NR), disc,ndisc,&
                  1,NR,integrd_r_im)
             integral_r = cmplx(integral_r_re, integral_r_im)

             if (l .eq. 0) del_rho_lm = cmplx(0,0)
             
             integrd_lm_r_prime(nlayer_r_prime) = integral_r * del_rho_lm *&
                  rho_norm(nlayer_r_prime) * rad(nlayer_r_prime)**2
          end do
          
          integrd_r_re = realpart(integrd_lm_r_prime)
          integrd_r_im = dimag(integrd_lm_r_prime)
          call intgrl_disc(integral_r_r_prime_lm_re, NR,rad(1:NR),disc,ndisc, &
               1,NR,integrd_r_re)
          call intgrl_disc(integral_r_r_prime_lm_im, NR,rad(1:NR),disc,ndisc, &
               1,NR,integrd_r_im)
          integral_r_r_prime_lm = integral_r_r_prime_lm +&
               cmplx(integral_r_r_prime_lm_re,integral_r_r_prime_lm_im)
       end do
       integral_r_r_prime_l = integral_r_r_prime_l + integral_r_r_prime_lm
    end do
    grav_potential = realpart(integral_r_r_prime_l)
    ! 
    return    
    !
  end subroutine compute_volumetric_integral_gravitational_energy

  subroutine compute_volumetric_integral_gravitational_energy_for_kernel(grav_potential, &
       rho_norm, rad, layer_perturb, l_perturb, m_perturb, value_perturb, &
       lmax, lmax_model,disc, ndisc, NR)

    use legendre_quadrature
    implicit none
    !
    integer, intent(in) :: NR, ndisc, lmax, lmax_model
    integer, allocatable, intent(in) :: disc(:)
    real*8, allocatable, intent(in) :: rad(:), rho_norm(:)
    ! debug
    integer, intent(in) :: layer_perturb, l_perturb, m_perturb
    real*8, intent(in) :: value_perturb
    !
    ! real*8, intent(out) :: integral_value
    real*8, intent(out) :: grav_potential
    !
    integer :: kind, order
    real*8, allocatable :: abs_int(:), w_int(:)
    real*8 :: a, b, alpha, beta
    complex*16 :: delta_rho(0:lmax_model,0:lmax_model)
    complex*16, allocatable :: delta_rho_all_r(:,:,:) !!!
    ! 
    integer :: i, j, nlayer, nlayer_r_prime, l, m
    real*8, allocatable :: integral_lon(:)
    real*8 :: rho_abs_lm
    real*8 :: colat, lon
    real*8 :: integral_lat
    complex*16 :: rho_abs, y_lm
    real*8, allocatable :: r_l(:)
    complex*16 :: integral_r, integral_r_r_prime_lm
    complex*16 :: del_rho_lm
    complex*16 :: integral_r_r_prime_l
    complex*16, allocatable :: integrd_r(:), integrd_lm_r_prime(:)
    real*8, allocatable :: integrd_r_re(:), integrd_r_im(:)
    real*8 :: integral_r_re, integral_r_im
    real*8 :: integral_r_r_prime_lm_re, integral_r_r_prime_lm_im
    
    real*8, parameter :: PI = 3.1415927, GRAV_CST = 6.67408d-11

    ! Parameters for the interation rule (move to arg of subroutine eventually)
    ! order = 3
    ! kind = 1
    ! alpha = 0.d0
    ! beta = 0.d0
    ! a = -1
    ! b = 1

    ! Allocation
    ! allocate(integral_lon(NR))
    ! allocate(abs_int(order), w_int(order))
    allocate(r_l(NR))
    allocate(integrd_r(NR), integrd_lm_r_prime(NR))
    allocate(integrd_r_re(NR), integrd_r_im(NR))

    r_l(:) = 0.d0
    integrd_r(:) = 0.d0
    integrd_lm_r_prime(:) = 0.d0
    integrd_r_re(:) = 0.d0
    integrd_r_im(:) = 0.d0
    grav_potential = cmplx(0.d0,0.d0)

    !
    ! call construct_rho_map(delta_rho_all_r, 331, 782, lmax_model)
    ! Major loop in l and m, we moved the summation outside the integral computation
    integral_r_r_prime_l = cmplx(0.d0,0.d0)
    do l = l_perturb,l_perturb
    ! l = l_perturb
    print*, "Computing l value ", l,"/",lmax
    integral_r_r_prime_lm = cmplx(0.d0,0.d0)
    ! do m = -l,l
    ! do m = -m_perturb,m_perturb,2*m_perturb
    do m = 0,0
       print*, "m = ", -m_perturb,m_perturb,2*m_perturb
          ! Second and final radial integration
          !do nlayer_r_prime = 331,332
          do nlayer_r_prime = 2,782
             ! Initialization for a given layer
             ! First volume integration, in r' and omega
             ! Starting at 2 to avoid the singularity at the center, ending right before the crust
             ! do nlayer = 331, 332
             do nlayer = layer_perturb,layer_perturb
                ! get the density perturbation, make a map of rho
                ! if (nlayer.ge.  .and. nlayer .lt. 782)then
                   ! call get_delta_rho(delta_rho, nlayer, lmax_model)
                ! delta_rho(:,:) = delta_rho_all_r(nlayer,:,:)
                delta_rho(:,:) = cmplx(value_perturb,value_perturb)
                ! else
                !    delta_rho(:,:) = cmplx(0.d0,0.d0)
                ! end if
                ! print*, "After delta rho"
                   
                ! No integration in lat/lon
                call compute_one_over_r_sph_harm(r_l, rad, nlayer, l, NR)
                ! print*, "After one over r"               
                if (l .eq. 0) del_rho_lm = cmplx(rho_norm(nlayer),0)
                if (m < 0) then
                   del_rho_lm = (-1.d0)**m * conjg(delta_rho(l,abs(m)))
                else
                   del_rho_lm = delta_rho(l,abs(m))
                end if
                
                integrd_r(nlayer) = r_l(nlayer) * conjg(del_rho_lm) * rad(nlayer)**2
                
             end do ! nlayer, r'         
             ! Radial integration
             integrd_r_re = realpart(integrd_r)
             integrd_r_im = dimag(integrd_r)
             call intgrl_disc(integral_r_re, NR,rad(1:NR), disc,ndisc,&
                  1,NR,integrd_r_re)
             call intgrl_disc(integral_r_im, NR,rad(1:NR), disc,ndisc,&
                  1,NR,integrd_r_im)
             integral_r = cmplx(integral_r_re, integral_r_im)

             if (l .eq. 0) del_rho_lm = cmplx(0,0)
             
             integrd_lm_r_prime(nlayer_r_prime) = integral_r * del_rho_lm *&
                  rho_norm(nlayer_r_prime) * rad(nlayer_r_prime)**2
          end do
          
          integrd_r_re = realpart(integrd_lm_r_prime)
          integrd_r_im = dimag(integrd_lm_r_prime)
          call intgrl_disc(integral_r_r_prime_lm_re, NR,rad(1:NR),disc,ndisc, &
               1,NR,integrd_r_re)
          call intgrl_disc(integral_r_r_prime_lm_im, NR,rad(1:NR),disc,ndisc, &
               1,NR,integrd_r_im)
          integral_r_r_prime_lm = integral_r_r_prime_lm +&
               cmplx(integral_r_r_prime_lm_re,integral_r_r_prime_lm_im)
       end do
       integral_r_r_prime_l = integral_r_r_prime_l + integral_r_r_prime_lm
    end do
    grav_potential = realpart(integral_r_r_prime_l)
    ! 
    return    
    !
  end subroutine compute_volumetric_integral_gravitational_energy_for_kernel


  ! subroutine compute_volumetric_integral_gravitational_energy_for_kernel(grav_potential, &
  !      rho_norm, rad, layer_perturb, l_perturb, m_perturb, value_perturb, &
  !      lmax, lmax_model,disc, ndisc, NR)

  !   use legendre_quadrature
  !   implicit none
  !   !
  !   integer, intent(in) :: NR, ndisc, lmax, lmax_model
  !   integer, intent(in) :: layer_perturb, l_perturb, m_perturb
  !   real*8, intent(in)  :: value_perturb
  !   integer, allocatable, intent(in) :: disc(:)
  !   real*8, allocatable, intent(in) :: rad(:), rho_norm(:)
  !   !
  !   ! real*8, intent(out) :: integral_value
  !   real*8, intent(out) :: grav_potential
  !   !
  !   integer :: kind, order
  !   real*8, allocatable :: abs_int(:), w_int(:)
  !   real*8 :: a, b, alpha, beta
  !   complex*16 :: delta_rho(0:lmax_model,0:lmax_model)
  !   complex*16, allocatable :: delta_rho_all_r(:,:,:) !!!
  !   ! 
  !   integer :: i, j, nlayer, nlayer_r_prime, l, m
  !   real*8, allocatable :: integral_lon(:)
  !   real*8 :: rho_abs_lm
  !   real*8 :: colat, lon
  !   real*8 :: integral_lat
  !   complex*16 :: rho_abs, y_lm
  !   real*8, allocatable :: r_l(:)
  !   complex*16 :: integral_r, integral_r_r_prime_lm
  !   complex*16 :: del_rho_lm
  !   complex*16 :: integral_r_r_prime_l
  !   complex*16, allocatable :: integrd_r(:), integrd_lm_r_prime(:)
  !   real*8, allocatable :: integrd_r_re(:), integrd_r_im(:)
  !   real*8 :: integral_r_re, integral_r_im
  !   real*8 :: integral_r_r_prime_lm_re, integral_r_r_prime_lm_im
    
  !   real*8, parameter :: PI = 3.1415927, GRAV_CST = 6.67408d-11

  !   ! Parameters for the interation rule (move to arg of subroutine eventually)
  !   ! order = 3
  !   ! kind = 1
  !   ! alpha = 0.d0
  !   ! beta = 0.d0
  !   ! a = -1
  !   ! b = 1

  !   ! Allocation
  !   ! allocate(integral_lon(NR))
  !   ! allocate(abs_int(order), w_int(order))
  !   allocate(r_l(NR))
  !   allocate(integrd_r(NR), integrd_lm_r_prime(NR))
  !   allocate(integrd_r_re(NR), integrd_r_im(NR))

  !   r_l(:) = 0.d0
  !   integrd_r(:) = 0.d0
  !   integrd_lm_r_prime(:) = 0.d0
  !   integrd_r_re(:) = 0.d0
  !   integrd_r_im(:) = 0.d0
  !   grav_potential = cmplx(0.d0,0.d0)

  !   !
  !   call construct_rho_map(delta_rho_all_r, 331, 782, lmax_model)
  !   ! Major loop in l and m, we moved the summation outside the integral computation
  !   integral_r_r_prime_l = cmplx(0.d0,0.d0)
  !   ! do l = l_perturb,l_perturb+1
  !   do l = 0,lmax
  !      print*, "Computing l value ", l,"/",lmax
  !      integral_r_r_prime_lm = cmplx(0.d0,0.d0)
  !      do m = -l,l
  !      ! do m = -m_perturb,m_perturb,2*m_perturb
  !      ! do m = m_perturb,m_perturb
  !         print*, "m = ", m
  !         ! Second and final radial integration
  !         !do nlayer_r_prime = 331,332
  !         do nlayer_r_prime = 2,782
  !            ! Initialization for a given layer
  !            ! First volume integration, in r' and omega
  !            ! Starting at 2 to avoid the singularity at the center, ending right before the crust
  !            ! do nlayer = 331, 332
  !            do nlayer = 2, 782
  !               ! get the density perturbation, make a map of rho
  !               ! if (nlayer .eq. layer_perturb)then
  !               if (nlayer .ge. layer_perturb .and. nlayer .lt. layer_perturb+1)then
  !                  ! call get_delta_rho(delta_rho, nlayer, lmax_model)
  !                  delta_rho(:,:) = delta_rho_all_r(nlayer,:,:)
  !                  ! CHECK THAT IT IS THE CORRECT VALUE FOR THE MATRIX
  !                  ! delta_rho(l_perturb,m_perturb) = cmplx(value_perturb, value_perturb)
  !                  ! print*, layer_perturb, delta_rho(l_perturb, m_perturb)
  !               else
  !                  delta_rho(:,:) = cmplx(0.d0,0.d0)
  !               end if
  !               ! print*, "After delta rho"
                   
  !               ! No integration in lat/lon
  !               call compute_one_over_r_sph_harm(r_l, rad, nlayer, l, NR)
  !               ! print*, "After one over r"               
  !               if (l .eq. 0) del_rho_lm = cmplx(rho_norm(nlayer),0)
  !               if (m < 0) then
  !                  del_rho_lm = (-1.d0)**m * conjg(delta_rho(l,abs(m)))
  !               else
  !                  del_rho_lm = delta_rho(l,abs(m))
  !               end if
                
  !               integrd_r(nlayer) = r_l(nlayer) * conjg(del_rho_lm) * rad(nlayer)**2
                
  !            end do ! nlayer, r'         
  !            ! Radial integration
  !            integrd_r_re = realpart(integrd_r)
  !            integrd_r_im = dimag(integrd_r)
  !            call intgrl_disc(integral_r_re, NR,rad(1:NR), disc,ndisc,&
  !                 1,NR,integrd_r_re)
  !            call intgrl_disc(integral_r_im, NR,rad(1:NR), disc,ndisc,&
  !                 1,NR,integrd_r_im)
  !            integral_r = cmplx(integral_r_re, integral_r_im)

  !            if (l .eq. 0) del_rho_lm = cmplx(0,0)
             
  !            integrd_lm_r_prime(nlayer_r_prime) = integral_r * del_rho_lm *&
  !                 rho_norm(nlayer_r_prime) * rad(nlayer_r_prime)**2
  !         end do
          
  !         integrd_r_re = realpart(integrd_lm_r_prime)
  !         integrd_r_im = dimag(integrd_lm_r_prime)
  !         call intgrl_disc(integral_r_r_prime_lm_re, NR,rad(1:NR),disc,ndisc, &
  !              1,NR,integrd_r_re)
  !         call intgrl_disc(integral_r_r_prime_lm_im, NR,rad(1:NR),disc,ndisc, &
  !              1,NR,integrd_r_im)
  !         integral_r_r_prime_lm = integral_r_r_prime_lm +&
  !              cmplx(integral_r_r_prime_lm_re,integral_r_r_prime_lm_im)
  !      end do
  !      integral_r_r_prime_l = integral_r_r_prime_l + integral_r_r_prime_lm
  !      ! DEBUG
  !      print*, integral_r_r_prime_l
  !   end do
  !   grav_potential = realpart(integral_r_r_prime_l)
  !   ! 
  !   return    
  !   !
  ! end subroutine compute_volumetric_integral_gravitational_energy_for_kernel


  
  ! <SA> THIS SUBROUTINE COMPUTES THE VOLUMIC INTEGRAL OF THE GRAVITATIONAL POTENTIAL ONLY
  ! Because the gravitational potential is a 3D value, we need to call a function to calculate
  ! it at each lat/lon.
  
  subroutine compute_grav_pot_volumetric_integral(integral_value, rho_norm, &
       rad, lmax, lmax_model,disc, ndisc, NR)

    use legendre_quadrature
    implicit none
    !
    integer, intent(in) :: NR, ndisc, lmax, lmax_model
    integer, allocatable, intent(in) :: disc(:)
    real*8, allocatable, intent(in) :: rad(:), rho_norm(:)
    !
    real*8, intent(out) :: integral_value
    !
    integer :: kind, order
    real*8, allocatable :: abs_int(:), w_int(:)
    real*8, allocatable :: kernel_grav_all_l_nlayer(:)
    real*8, allocatable :: kernel_grav_all_l(:,:)
    real*8 :: rho_lm
    complex*16 :: delta_rho(0:lmax_model,0:lmax_model)
    ! 
    integer :: i, j, nlayer, l, m,ii
    real*8 :: a, b, alpha, beta
    real*8, allocatable :: integral_lon(:)
    real*8 :: lat, lon
    real*8 :: integral_lat, del_phi
    real*8, parameter :: PI = 3.1415927

    ! Parameters for the integration rule (move to arg of subroutine eventually)
    order = 200
    ! order = 30
    kind = 1
    alpha = 0.d0
    beta = 0.d0
    a = -1
    b = 1

    ! Compute knots and weights for the Gauss-Legendre quadrature
    allocate(abs_int(order), w_int(order))
    call cgqf( order, kind, alpha, beta, a, b, abs_int, w_int)

    ! Initialization
    allocate(integral_lon(NR))
    allocate(kernel_grav_all_l(1:NR,0:lmax))
    allocate(kernel_grav_all_l_nlayer(0:lmax))
    integral_value = 0.d0
    integral_lon(:) = 0.d0

    ! Computation of the kernels for all l.
    call compute_kernel_grav(kernel_grav_all_l, rho_norm, rad, lmax, NR, ndisc, disc)
    ! <SA> DEBUG LATLON
    open(41,file="lat_lon.txt",form="formatted",status="replace")

    ! Integration on the volume per se
    ! do nlayer = 331, NR-26
    do nlayer = 600, 600
       ! Starting above the core mantle boundary as the core is 1D,
       ! hence making delta_phi = 0
       ! Get the kernel for all l for this specfic layer
       kernel_grav_all_l_nlayer = kernel_grav_all_l(nlayer,:)
       ! Get the density perturbation
       call get_delta_rho(delta_rho, nlayer, lmax_model)
       ! <SA> DEBUG
       ! if (nlayer .eq. 600) then
       !    print*, 'coucou', nlayer
       !    do l = 0,lmax
       !       do m = 0,lmax
       !          print*, l, m, delta_rho(l,m)
       !       end do
       !    end do          
       ! end if
       ! </SA>

       print*, '[compute_volume_integral] layer', nlayer,'/',NR
       ! Loop on longitude between -1 and 1 (transformation)
       do j = 1, 2*order
          lon = ((j-1/2) * PI / order) ! Make things cleaner here
          integral_lat = 0.d0
          do i = 1, order
             ! transform into geographical latitude (radians)
             lat = acos(abs_int(i))
             ! lat = abs_int(i)
             if (nlayer .eq. 331)then
                write(41,"(3F10.5)") abs_int(i),lat,lon
             endif
             ! Compute delta phi, summed over l and m.
             ! call compute_delta_phi(del_phi, rho_lm, rho_norm, rad, delta_rho, &
             !      kernel_grav_all_l_nlayer, lmax, lmax_model, nlayer, lat, lon, NR)
             ! integral_lat = integral_lat + w_int(i) * rho_lm * del_phi
             ! integral_lat = integral_lat + w_int(i) * (del_phi * (1 - abs_int(i)**2)**(0.5))
             del_phi = cexp(cmplx(0,lat))!  * cmplx((1-abs_int(i)**2)**(0.5))
             ! del_phi = del_phi * ((1-abs_int(i)**2)**(0.5))
             ! integral_lat = integral_lat + w_int(i) * abs(del_phi)
             integral_lat = integral_lat + w_int(i) * del_phi
          end do ! longitude
          ! Sum over all latitudes of integration
          integral_lon(nlayer) = integral_lon(nlayer) + integral_lat * (PI/order) 
       end do ! latitude
       ! Normalize the surface integral and multiply by dr**2 to prepare the radial integration
       ! integral_lon(nlayer) = integral_lon(nlayer) * rad(nlayer)**2
       print*, "[compute_volumetric_integral.f90]", integral_lon(nlayer)
    end do ! nlayer

    ! <SA> DEBUG
    close(41)
    ! Radial integration
    ! call intgrl_disc(integral_value,NR,rad(1:NR), disc,ndisc,1,NR,integral_lon(1:NR))
    ! integral_value = integral_value/2
    ! print*, "[compute_volumetric_integral.f90]",integral_value

    deallocate(integral_lon)
    deallocate(kernel_grav_all_l)
    deallocate(kernel_grav_all_l_nlayer)
    deallocate(abs_int, w_int)

    return

  end subroutine compute_grav_pot_volumetric_integral

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module volumetric_integral
