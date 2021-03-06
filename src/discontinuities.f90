module discontinuities
  !
contains
  !
  subroutine find_disc(model_file, disc, rdisc, ndisc, NR)

    implicit none

    character, intent(in) :: model_file*100
    integer, intent(out) :: ndisc, NR
    integer, allocatable, intent(out) :: disc(:)
    real*8, allocatable, intent(out)  :: rdisc(:)
    !
    real*8 :: radius_last,radius,eps
    integer :: jdisc,iomod,i
    integer :: ierr
    character*80 :: junk ! HL
    double precision :: dummy(3)

    ! HL changed model format
    ! no need for '.dk' files, just '.md' files

    radius_last = -1.0
    eps = 1.0d-8
    jdisc = 0
    iomod = 11
    open(iomod,file=model_file,status='old')
    ! First pass to determine the number of discontinuities
    read(iomod,*) junk !HL
    read(iomod,*) dummy !HL, SA
    read(iomod,*) dummy !HL, SA

    NR = int(dummy(1))

    do i = 1, NR
       read(iomod,*) radius !HL
       if (abs(radius-radius_last) .lt. eps) then
          jdisc = jdisc + 1
       end if
       radius_last = radius
    enddo
    ndisc = jdisc
    allocate(disc(ndisc))
    allocate(rdisc(ndisc))

    ! 2nd pass to fill in the radii of the discontinuities
    jdisc = 0
    radius_last = -1.0

    rewind(iomod)
    read(iomod,*) junk !HL
    read(iomod,*) junk !HL
    read(iomod,*) junk !HL
    do i = 1, NR
       read(iomod,*) radius !HL
       if (abs(radius-radius_last) .lt. eps) then
          jdisc = jdisc + 1
          rdisc(jdisc) = radius
          disc(jdisc) = i - 1
       endif
       radius_last = radius
    enddo
    close(iomod)

  end subroutine find_disc

  subroutine deriv(y,yprime,n,r,ndis,kdis,s1,s2,s3)
    !     slightly altered dp version of subroutine rspln.

    implicit none
    !
    integer, intent(in) :: n, ndis
    integer, intent(in) :: kdis(:)
    real*8, intent(in) :: y(:), r(:)
    real*8, allocatable, intent(out) :: s1(:), s2(:), s3(:), yprime(:)
    ! 
    integer :: nd, ndp, i, j, j1, j2, k
    real*8 :: f(3,1000), yy(3), y0
    real*8 :: a0, h, ha, h2, h2a, h2b, h3a, b0, b1, s32, s21, s13

    equivalence (yy(1),y0)
    data yy/3*0.0/

    allocate(s1(n), s2(n), s3(n),yprime(n))
    ndp=ndis+1
    do nd=1,ndp
       if(nd.eq.1) go to 4
       if(nd.eq.ndp) go to 5
       j1=kdis(nd-1)+1
       j2=kdis(nd)-2
       go to 6
4      j1=1
       j2=kdis(1)-2
       go to 6
5      j1=kdis(ndis)+1
       j2=n-2
6      if((j2+1-j1).gt.0) go to 11
       j2=j2+2
       y0=(y(j2)-y(j1))/(r(j2)-r(j1))
       s1(j1)=yy(1)
       s1(j2)=yy(1)
       s2(j1)=yy(2)
       s2(j2)=yy(2)
       s3(j1)=yy(3)
       s3(j2)=yy(3)
       go to 3
11     a0=0.0
       if(j1.eq.1) go to 7
       h=r(j1+1)-r(j1)
       h2=r(j1+2)-r(j1)
       y0=h*h2*(h2-h)
       h=h*h
       h2=h2*h2
       b0=(y(j1)*(h-h2)+y(j1+1)*h2-y(j1+2)*h)/y0
       go to 8
7      b0=0.0
8      b1=b0
       if(j2 .gt. 1000)write(0,'("error:deriv:j2= ",i5)')j2
       do i=j1,j2
          h=r(i+1)-r(i)
          y0=y(i+1)-y(i)
          h2=h*h
          ha=h-a0
          h2a=h-2.0*a0
          h3a=2.0*h-3.0*a0
          h2b=h2*b0
          s1(i)=h2/ha
          s2(i)=-ha/(h2a*h2)
          s3(i)=-h*h2a/h3a
          f(1,i)=(y0-h*b0)/(h*ha)
          f(2,i)=(h2b-y0*(2.0*h-a0))/(h*h2*h2a)
          f(3,i)=-(h2b-3.0*y0*ha)/(h*h3a)
          a0=s3(i)
          b0=f(3,i)
       enddo
       i=j2+1
       h=r(i+1)-r(i)
       y0=y(i+1)-y(i)
       h2=h*h
       ha=h-a0
       h2a=h*ha
       h2b=h2*b0-y0*(2.*h-a0)
       s1(i)=h2/ha
       f(1,i)=(y0-h*b0)/h2a
       ha=r(j2)-r(i+1)
       y0=-h*ha*(ha+h)
       ha=ha*ha
       y0=(y(i+1)*(h2-ha)+y(i)*ha-y(j2)*h2)/y0
       s3(i)=(y0*h2a+h2b)/(h*h2*(h-2.0*a0))
       !the following statements were expanded to prevent register overflow on unix
       s13=s1(i)*s3(i)
       s2(i)=f(1,i)-s13
       do  j=j1,j2
          k=i-1
          s32=s3(k)*s2(i)
          s1(i)=f(3,k)-s32
          s21=s2(k)*s1(i)
          s3(k)=f(2,k)-s21
          s13=s1(k)*s3(k)
          s2(k)=f(1,k)-s13
          i=k
       enddo
       s1(i)=b1
       j2=j2+2
       s1(j2)=yy(1)
       s2(j2)=yy(2)
       s3(j2)=yy(3)
3   enddo
    do  i=1,n
       yprime(i)=s1(i)
    enddo

    return

  end subroutine deriv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine intgrl_disc(sum,n,r,idisc,ndis,nir,ner,f)

    !-----------------------------------------------------------
    !  Usage: Computes the integral of f[i] from i=nir to i=ner for
    !  Input : n -- number of total layers
    !          r (1:n) -- normalized radius
    !          idisc(1:ndis) -- discontinuity indices
    !          ndis -- number of discontinuities
    !          nir,ner -- start and end indices
    !          f -- function to be integrated.
    !  Output: sum -- the integral
    !--------------------------------------------------------------
    !
    implicit none
    !
    real*8,  intent(in) :: f(:), r(:)
    integer, intent(in) :: n, nir, ner, ndis
    integer, intent(in) :: idisc(:)
    !
    real*8, intent(out) :: sum
    !
    real*8, allocatable :: yprime(:), s1(:), s2(:), s3(:)
    integer :: kdis(ndis)
    integer :: nir1, i, j
    real*8  :: third, fifth, sixth, rji

    kdis = 0
    kdis(1:ndis) = idisc(1:ndis)

    third = 1.0/3.0
    fifth = 1.0/5.0
    sixth = 1.0/6.0

    call deriv(f,yprime,n,r,ndis,kdis,s1,s2,s3)
    nir1 = nir + 1
    sum = 0.0

    do i=nir1,ner
       j = i-1
       rji = r(i) - r(j)
       sum=sum+rji*(f(j)+rji*(.50*s1(j)+rji*(third*s2(j)+rji* &
            .250*s3(j))))+2.0*r(j)*rji*rji*(.50*f(j)+rji*(third*s1(j)+rji* &
            (.250*s2(j)+rji*fifth*s3(j))))+rji*rji*rji*(third*f(j)+rji*  &
            (.250*s1(j)+rji*(fifth*s2(j)+rji*sixth*s3(j))))
    enddo

  end subroutine intgrl_disc


end module discontinuities
