module volumetric_integral
  !
  use legendre_quadrature
  use discontinuities
  !
  implicit none
  !
contains

  subroutine compute_volumetric_integral(integral_value, func_to_integrate, rad, disc, ndisc, NR)

    use legendre_quadrature
    implicit none
    !
    integer, intent(in) :: NR, ndisc
    integer, allocatable, intent(in) :: disc(:)
    real*8, allocatable, intent(in) :: rad(:), func_to_integrate(:)
    !
    real*8, intent(out) :: integral_value
    !
    integer :: kind, order
    real*8, allocatable :: abs_int(:), w_int(:)
    real*8 :: a, b, alpha, beta
    ! 
    integer :: i, j, nlayer
    real*8, allocatable :: integral_lon(:)
    real*8 :: lat, lon
    real*8 :: integral_lat
    real*8, parameter :: PI = 3.1415927

    ! Parameters for the integration rule (move to arg of subroutine eventually)
    order = 8
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
    integral_value = 0.d0
    integral_lon(:) = 0.d0

    ! Integration on the volume per se
    do nlayer = 2, NR ! Starting at 2 to avoid the singularity at the center
       ! Loop on longitude between -1 and 1 (transformation)
       print*, '[compute_volume_integral] layer', nlayer,'/',NR
       do j = 1, 2*order
          lon = (j * PI / order) ! Make things cleaner here
          integral_lat = 0.d0
          do i = 1, order
             ! transform into geographical latitude (radians)
             lat = acos(abs_int(i))
             ! call compute_phi(phi, phi_zero, nlayer, lat, lon, NR)
             integral_lat = integral_lat + w_int(i) * func_to_integrate(nlayer)
             ! integral_lat = integral_lat + w_int(i) * phi 
          end do ! longitude
          ! Sum over all latitudes of integration
          integral_lon(nlayer) = integral_lon(nlayer) + integral_lat 
       end do ! latitude
       ! Normalize the surface integral and multiply by dr**2 to prepare the radial integration
       integral_lon(nlayer) = (PI/order) * integral_lon(nlayer) * rad(nlayer)**2

    end do ! nlayer

    ! Radial integration
    call intgrl_disc(integral_value,NR,rad(1:NR), disc,ndisc,1,NR,integral_lon(1:NR))

    return

  end subroutine compute_volumetric_integral
end module volumetric_integral
