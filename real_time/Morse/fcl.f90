        subroutine derivs(x,y,Ntraj,dvdx,dvdy,pe)
        implicit none
        integer*8,intent(IN) :: Ntraj
        real*8,intent(IN) :: x(Ntraj),y(Ntraj)
        integer*4 :: i
        real*8 :: De,a,k1,k2,x0,y0,b,l,xi(Ntraj),pi
        real*8,intent(OUT) :: dvdx(Ntraj),dvdy(Ntraj),pe(Ntraj)
	common/group/x0,y0
!       harmonic potential for x coordinate, Eckart potential for y
!       direction  
        De=7d0
        !a=1.04435d0
        k1=2d0
        k2=2d0
	A=1d0
	B=-2d0
	l=1d0
	pi=3.1415926d0
        do i=1,Ntraj
	!dvdx(i)=2*a*De*(1-exp(-a*x(i)))*exp(-a*x(i))
        dvdx(i)=k1*x(i)
        !pe(i)=k1*x(i)**2/2d0+k2*y(i)**2/2d0
        !dvdy(i)=-43.5968*sinh(1.3624*y(i))/
     !& (cosh(1.3624*y(i)))**3
        !pe(i)=De*(1-exp(-a*x(i)))**2+16d0/(cosh(1.3624*y(i)))**2
	xi(i)=-exp(2d0*pi*y(i)/l)
	dvdy(i)=2d0*pi*xi(i)*(A-A*xi(i)+B+
     &  B*xi(i))/((1-xi(i))**3*l)
	pe(i)=k1*x(i)**2/2d0+A*xi(i)/(1d0-xi(i))+B*xi(i)/(1d0-xi(i))**2
	
        end do
        return
        end subroutine
