module ylm_procedures
  implicit none

contains
  
  subroutine ylm(clat,phi,l,msign,yval)

    !-----------------------------------------------------------
    ! Usage: computes the spherical harmonic (normalization of Edmonds, 1960)
    ! Input:  clat -> colat (rad)
    !         phi -> longitude (rad, 0-2pi)
    !         l,msign -> degree and order
    ! Output: yval -> ylm value at the given site.

    implicit none

    ! inputs
    integer, intent(in) :: l,msign
    real*8, intent (in) :: clat
    real*8, intent (in) :: phi
    ! outputs
    complex*8, intent(out) :: yval
    ! others:
    integer :: m,l0
    real*8 :: cost,sint,mm
    real*8 :: x(l+1),dx(l+1)
    complex*16 :: imphi

    m = abs(msign)
    l0 = l
    cost = cos(clat)
    sint = sin(clat)
    call lgndr(l0,cost,sint,x,dx)
    mm=real(m)
    imphi = cmplx(0.0,mm*phi)
    yval=exp(imphi)*x(m+1)
    if (msign.lt.0) yval=conjg(yval)*(-1.**msign)

    return

  end subroutine ylm



  subroutine lgndr(l,c,s,x,dx)
    !
    !    computes legendre function x(l,m,theta)
    !    theta=colatitude,c=cos(theta),s=sin(theta),l=angular order,
    !    sin(theta) restricted so that sin(theta).ge.1.e-7
    !    x(1) contains m=0, x(2) contains m=1, x(k+1) contains m=k
    !    m=azimuthal(longitudenal) order 0.le.m.le.l
    !    dx=dx/dtheta
    !    SUBROUTINE originally came from Physics Dept. Princeton through
    !    Peter Davis
    !    modified to run stably on the Perkin-Elmer by Jeffrey Park 11/23/85
    !
    !----------------------------------------------------
    !    calculate X_lm(theta)= 
    !             (-1)^m * sqrt((2l+1)/4*pi) * sqrt[(l-m)!/(l+m)!] * P_lm(cos_theta)
    !                                            (B.58)
    !                 where P_lm is the "associated lengendre function"
    !
    !              dX_lm/d_theta
    !----------------------------------------------------
    IMPLICIT NONE
    INTEGER,INTENT(INOUT) :: l
    REAL*8,   INTENT(IN   ) :: c
    REAL*8,   INTENT(INOUT) :: s
    REAL*8,   INTENT(  OUT) :: x(:),dx(:)  !x(l+1),dx(l+1)

    REAL*8 :: tol, rfpi, root3, boeing

    REAL*8 :: c1,c2,sos,cot,ct,ss,g3,g2,g1,f3,f2,f1,w,v,y,z,d,t,&
         stom,fac
    INTEGER :: lsave,lp1,i,mp1,m,lpsafex,lpsafe,mmm,maxsin

    tol=1.e-5
    rfpi=0.282094791773880  !sqrt(1/4(pi))
    root3=1.73205080756890
    boeing=0.707106781186550
    x(:)=0.0 ! set all to zero
    dx(:)=0.0 ! set all to zero

    if(s.ge.1.0-tol) s=1.0-tol
    lsave=l
    if(l.lt.0) THEN
       l=-1-l
       WRITE(*,*) 'are you sure l is negtive in SUB lgndr?? STOP'
       STOP
    END IF
    if(l.gt.0) go to 1
    x(1)=rfpi
    dx(1)=0.0
    l=lsave
    return
1   if(l.ne.1) go to 2
    c1=root3*rfpi
    c2=boeing*c1
    x(1)=c1*c
    x(2)=-c2*s
    dx(1)=-c1*s
    dx(2)=-c2*c
    l=lsave
    return
2   sos=s
    if(s.lt.tol) s=tol
    cot=c/s
    ct=2.0*c
    ss=s*s
    lp1=l+1
    g3=0.0
    g2=1.0
    f3=0.0
    !  evaluate m=l value, sans (sin(theta))**l
    do 100 i=1,l
100    g2=g2*(1.0-1.0/(2.0*i))
       g2=rfpi*sqrt((2*l+1)*g2)
       f2=l*cot*g2
       x(lp1)=g2
       dx(lp1)=f2
       w=0.0
       v=1.0
       y=2.0*l
       z=y+1.0
       d=sqrt(v*y)
       t=0.0
       mp1=l
       m=l-1
       !  these recursions are similar to ordinary m-recursions, but since we
       !  have taken the s**m factor out of the xlm's, the recursion has the powers
       !  of sin(theta) instead
3      g1=-(ct*mp1*g2+ss*t*g3)/d
       f1=(mp1*(2.0*s*g2-ct*f2)-t*ss*(f3+cot*g3))/d-cot*g1
       x(mp1)=g1
       dx(mp1)=f1
       if(m.eq.0) go to 4
       mp1=m
       m=m-1
       v=v+1.0
       y=y-1.0
       t=d
       d=sqrt(v*y)
       g3=g2
       g2=g1
       f3=f2
       f2=f1
       go to 3
4      maxsin=-72.0/log10(s)
       !  maxsin is the max exponent of sin(theta) without underflow
       lpsafe=min0(lp1,maxsin)
       stom=1.0
       fac=sign(1.0,(l/2)*2-l+.50)
       !  multiply xlm by sin**m
       do 5 m=1,lpsafe
          x(m)=fac*x(m)*stom
          dx(m)=fac*dx(m)*stom
          stom=stom*s
5         continue
          !  set any remaining xlm to zero
          if(maxsin.le.l) then
             mmm=maxsin+1
             do 200 m=mmm,lp1
                x(m)=0.0
200             dx(m)=0.0
             endif
             s=sos
             l=lsave
             return
           END SUBROUTINE lgndr


end module
