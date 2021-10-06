module volumetric_integral
  !
  use legendre_quadrature
  use discontinuities
  use gravitational_potential
  !
  implicit none
  !
contains

  ! subroutine compute_volumetric_integral(integral_value, func_to_integrate, rad, disc, ndisc, NR)

  !   use legendre_quadrature
  !   implicit none
  !   !
  !   integer, intent(in) :: NR, ndisc
  !   integer, allocatable, intent(in) :: disc(:)
  !   real*8, allocatable, intent(in) :: rad(:), func_to_integrate(:)
  !   !
  !   real*8, intent(out) :: integral_value
  !   !
  !   integer :: kind, order
  !   real*8, allocatable :: abs_int(:), w_int(:)
  !   real*8 :: a, b, alpha, beta
  !   ! 
  !   integer :: i, j, nlayer
  !   real*8, allocatable :: integral_lon(:)
  !   real*8 :: lat, lon
  !   real*8 :: integral_lat
  !   real*8, parameter :: PI = 3.1415927

  !   ! Parameters for the integration rule (move to arg of subroutine eventually)
  !   order = 8
  !   kind = 1
  !   alpha = 0.d0
  !   beta = 0.d0
  !   a = -1
  !   b = 1

  !   ! Compute knots and weights for the Gauss-Legendre quadrature
  !   allocate(abs_int(order), w_int(order))
  !   call cgqf( order, kind, alpha, beta, a, b, abs_int, w_int)

  !   ! Initialization
  !   allocate(integral_lon(NR))
  !   integral_value = 0.d0
  !   integral_lon(:) = 0.d0

  !   ! Integration on the volume per se
  !   do nlayer = 2, NR ! Starting at 2 to avoid the singularity at the center
  !      ! Loop on longitude between -1 and 1 (transformation)
  !      print*, '[compute_volume_integral] layer', nlayer,'/',NR
  !      do j = 1, 2*order
  !         lon = (j * PI / order) ! Make things cleaner here
  !         integral_lat = 0.d0
  !         do i = 1, order
  !            ! transform into geographical latitude (radians)
  !            lat = acos(abs_int(i))
  !            ! call compute_phi(phi, phi_zero, nlayer, lat, lon, NR)
  !            integral_lat = integral_lat + w_int(i) * func_to_integrate(nlayer)
  !            ! integral_lat = integral_lat + w_int(i) * phi 
  !         end do ! longitude
  !         ! Sum over all latitudes of integration
  !         integral_lon(nlayer) = integral_lon(nlayer) + integral_lat 
  !      end do ! latitude
  !      ! Normalize the surface integral and multiply by dr**2 to prepare the radial integration
  !      integral_lon(nlayer) = (PI/order) * integral_lon(nlayer) * rad(nlayer)**2

  !   end do ! nlayer

  !   ! Radial integration
  !   call intgrl_disc(integral_value,NR,rad(1:NR), disc,ndisc,1,NR,integral_lon(1:NR))

  !   return

  ! end subroutine compute_volumetric_integral

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
    complex*8 :: delta_rho(0:lmax_model,0:lmax_model)
    ! 
    integer :: i, j, nlayer, l, m,ii
    real*8 :: a, b, alpha, beta
    real*8, allocatable :: integral_lon(:)
    real*8 :: lat, lon
    real*8 :: integral_lat, del_phi
    real*8, parameter :: PI = 3.1415927

    ! Parameters for the integration rule (move to arg of subroutine eventually)
    order = (lmax*2+1)
    kind = 1
    alpha = 0.d0
    beta = 0.d0
    a = -1
    b = 1

    ! Compute knots and weights for the Gauss-Legendre quadrature
    allocate(abs_int(order), w_int(order))
    call cgqf( order, kind, alpha, beta, a, b, abs_int, w_int)

    ! <SA> DEBUG INTEGRATION
    open(21,file="int_nodes.txt", form="formatted",status="replace")
 
    do ii = 1,order
       write(21,"(2F10.2)")abs_int(ii), w_int(ii)
    end do
    close(21)
    ! </SA>
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
    do nlayer = 331, NR-26
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
          lon = (j * PI / order) ! Make things cleaner here
          integral_lat = 0.d0
          do i = 1, order
             ! transform into geographical latitude (radians)
             lat = acos(abs_int(i))
             if (nlayer .eq. 331)then
                write(41,"(3F10.5)") abs_int(i),lat,lon
             endif
             ! Compute delta phi, summed over l and m.
             call compute_delta_phi(del_phi, rho_lm, rho_norm, rad, delta_rho, &
                  kernel_grav_all_l_nlayer, lmax, nlayer, lat, lon, NR)
             ! integral_lat = integral_lat + w_int(i) * rho_lm * del_phi
             integral_lat = integral_lat + w_int(i) * del_phi
          end do ! longitude
          ! Sum over all latitudes of integration
          integral_lon(nlayer) = integral_lon(nlayer) + integral_lat 
       end do ! latitude
       ! Normalize the surface integral and multiply by dr**2 to prepare the radial integration
       integral_lon(nlayer) = (PI/order) * integral_lon(nlayer) * rad(nlayer)**2
       print*, "[compute_volumetric_integral.f90]", integral_lon(nlayer)
    end do ! nlayer

    ! <SA> DEBUG
    close(41)
    ! Radial integration
    call intgrl_disc(integral_value,NR,rad(1:NR), disc,ndisc,1,NR,integral_lon(1:NR))
    integral_value = integral_value/2
    print*, "[compute_volumetric_integral.f90]",integral_value

    deallocate(integral_lon)
    deallocate(kernel_grav_all_l)
    deallocate(kernel_grav_all_l_nlayer)
    deallocate(abs_int, w_int)

    return

  end subroutine compute_grav_pot_volumetric_integral

end module volumetric_integral
