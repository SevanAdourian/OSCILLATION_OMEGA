module moment_of_inertia
  !
  use legendre_quadrature
  use discontinuities
  use gravitational_potential
  !
  implicit none
  !
contains
  !
  subroutine compute_moment_of_inertia(moment_of_inertia, rho_norm, &
       rad, layer_min, layer_max, lmax, lmax_model, disc, ndisc, NR)

    use legendre_quadrature
    implicit none
    !
    integer, intent(in) :: NR, ndisc, lmax, lmax_model, layer_min, layer_max
    integer, allocatable, intent(in) :: disc(:)
    real*8, allocatable, intent(in) :: rad(:), rho_norm(:)
    !
    real*8, intent(out) :: moment_of_inertia
    !
    complex*8 :: delta_rho(0:lmax_model,0:lmax_model)
    !
    integer :: kind, order
    real*8, allocatable :: abs_int(:), w_int(:)
    !
    integer :: i, j, nlayer, l, m, number_of_layers
    real*8 :: a, b, alpha, beta
    real*8, allocatable :: integral_lon(:), radial_integrd(:)
    real*8 :: lat, lon
    real*8 :: radial_integrl_1, integral_lat, integrd
    real*8 :: moment_of_inertia_hetero, moment_of_inertia_radial
    ! 
    real*8, parameter :: PI = 3.1415927

    ! Parameters for the integration rule (move to arg of subroutine eventually - no need
    ! it just has to depend on l <SA>)
    order = (lmax*2+1)
    kind = 1
    alpha = 0.d0
    beta = 0.d0
    a = -1
    b = 1

    ! Compute knots and weights for the Gauss-Legendre quadrature
    number_of_layers = layer_max - layer_min + 1
    allocate(abs_int(order), w_int(order))
    call cgqf( order, kind, alpha, beta, a, b, abs_int, w_int)

    ! Initialization
    allocate(radial_integrd(number_of_layers))
    allocate(integral_lon(NR))
    integral_lon(:) = 0.d0
    radial_integrd(:) = 0.d0
    moment_of_inertia = 0.d0

    ! <SA> DEBUG delta_rho = 0.d0
    delta_rho(:,:) = 0.d0
    ! </SA>
    ! Computation of the kernels for all l.
    ! call compute_kernel_moi(kernel_moi_all_l, rho_norm, rad, lmax, NR, ndisc, disc)
    ! <SA> DEBUG LATLON
    ! open(41,file="lat_lon.txt",form="formatted",status="replace")

    ! Compute 1D part of moment of inertia
    radial_integrd(:) = rho_norm(layer_min:layer_max) * rad(layer_min:layer_max)**4
    call intgrl_disc(radial_integrl_1,NR,rad(1:NR), disc,ndisc,layer_min,layer_max, &
         radial_integrd(1:NR))
    moment_of_inertia_radial = radial_integrl_1 * 8.d0 * PI / 3.d0
    
    ! Integration on the volume per se
    if (layer_max > 331) then
       do nlayer = 331, layer_max
          ! Starting above the core mantle boundary as the core is 1D,
          ! hence making delta_rho = 0
          
          ! Get the density perturbation
          call get_delta_rho(delta_rho, nlayer, lmax_model)
          
          print*, '[compute_volume_integral] layer', nlayer,'/',NR
          ! Loop on longitude between -1 and 1 (transformation)
          do j = 1, 2*order
             lon = (j * PI / order) ! Make things cleaner here
             integral_lat = 0.d0
             do i = 1, order
                ! transform into geographical latitude (radians)
                lat = acos(abs_int(i))
                ! Compute integrand of moment of inertia, summed over l and m.
                call compute_moi_hetero(integrd, rho_norm, rad, delta_rho, &
                     lmax, nlayer, lat, lon, NR)
                integral_lat = integral_lat + w_int(i) * integrd
             end do ! latitude
             ! Sum over all latitudes of integration
             integral_lon(nlayer) = integral_lon(nlayer) + integral_lat * (PI/order) 
          end do ! longitude
          ! Normalize the surface integral and multiply by dr**2 to prepare the radial integration
          ! integral_lon(nlayer) = integral_lon(nlayer)
          print*, "[compute_volumetric_integral.f90]", integral_lon(nlayer)
       end do ! nlayer
       ! <SA> DEBUG
       close(41)
       ! Radial integration
       call intgrl_disc(moment_of_inertia_hetero,NR,rad(1:NR), disc,ndisc,&
            layer_min,layer_max,integral_lon(1:NR))
    else
       moment_of_inertia_hetero = 0.d0
    end if
    moment_of_inertia = moment_of_inertia_radial + moment_of_inertia_hetero
    
    print*, "[compute_volumetric_integral.f90]",moment_of_inertia, moment_of_inertia_radial, &
         moment_of_inertia_hetero

    deallocate(integral_lon)
    deallocate(radial_integrd)
    deallocate(abs_int, w_int)

    return

  end subroutine compute_moment_of_inertia

  subroutine compute_moi_hetero(integrd_sum, rho, r, delta_rho, &
       lmax, nlayer, lat, lon, NR)
    !
    implicit none
    !
    integer, intent(in) :: lmax, nlayer, NR
    real*8, intent(in)  :: lat, lon
    real*8, allocatable, intent(in) :: rho(:), r(:)
    complex*8, intent(in) :: delta_rho(:,:)
    !
    real*8, intent(out) :: integrd_sum
    ! 
    complex*8 :: integrd_lm
    complex*8 :: y_lm
    real*8 :: integrd_m
    integer :: l, m
    
    
    integrd_sum = 0
    do l = 0,lmax
       integrd_m = 0
       ! No need to do it here, call it directly from the volumetric integration routine
       ! so we don't do the same computations for each layers.
       ! call compute_kernel_grav()
       do m = 0,l
          ! MAYBE MOVE THIS TO A SEPARATE SUBROUTINE IN ORDER TO USE IT TO EXPAND RHO AS WELL
          ! Expand in spatial domain from spherical harmonics
          call ylm(lat, lon, l, m, y_lm)
          ! Computing here absolute value for rho for a given l and m, used to compute
          ! the integral.
          integrd_lm = delta_rho(l,m) * rho(nlayer) * r(nlayer)**4 * sin(lat)**3
          ! sum chi lm over m
          if (m == 0) then
             integrd_m = integrd_m + realpart(integrd_lm * y_lm)
          else
             integrd_m = integrd_m + 2 *(realpart(integrd_lm * y_lm))
          end if
       end do
       integrd_sum = integrd_sum + integrd_m
    end do
    !
  end subroutine compute_moi_hetero
  !
end module moment_of_inertia
