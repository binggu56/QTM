        program main 
        implicit real*8(a-h,o-z)
      	real*8 :: t,k1,k2,m1,m2,xav,yav,ki
	integer*8 :: Ntraj
     	real*8,allocatable :: ke(:),pe(:),x(:),y(:),px(:),py(:),w(:)
        real*8,allocatable :: qfx(:),qfy(:),qp(:),fx(:),fy(:),s(:)
     	real*8 :: dv,dt,dt2,f(3),qx0,qy0,prob
	integer :: k,kmax,i,j,ipot
	complex*16 :: im

      real*8 ::ax,ay,px0,py0,xx,yy
	real :: gasdev
      integer :: idum1,idum2
      COMMON/PES/ipot
        open(100,file='energy.dat')
        open(101,file='xoutput')
        open(102,file='xav')
	open(103,file='prob')
	open(104,file='youtput') 
        open(105,file='ini')
        open(5,file='IN')
	read(5,*) Ntraj
        read(5,*) kmax,dt
        read(5,*) m1,m2
	read(5,*) idum1,idum2
        read(5,*) ax,ay
        read(5,*) qx0,qy0
        read(5,*) px0,py0
        read(5,*) ipot 
        close(5)
        write(*,*) 'Quantum Trajectories'

	allocate(ke(Ntraj),pe(Ntraj),x(Ntraj),
     & y(Ntraj),px(Ntraj),py(Ntraj),s(Ntraj),
     & w(Ntraj),qfx(Ntraj),qfy(Ntraj),
     & qp(Ntraj),fx(Ntraj),fy(Ntraj))

        dt2=dt/2d0
        t=0d0
	im=(0d0,1d0)
!C        m1 = m1*2d0*m2/(m1+2d0*m2)
!C        m2 = m2/2d0
!C        m1 = m1*1836d0
!C        m2 = m2*1836d0
        !open(103,file='xyoutput')
	do i=1,Ntraj
	x(i)=gasdev(idum1)
	y(i)=gasdev(idum2)
        x(i)=x(i)/sqrt(4d0*ax)+qx0
        y(i)=y(i)/sqrt(4d0*ay)+qy0
	write (105,*) x(i),y(i)
	end do 
       
        px = px0
        py = py0
        
! 	center for the potential energy
!	x0=0d0
!	y0=0d0

        write(*,*) 'Initial Conditions'
        write(*,1001) 'ax =',ax
        write(*,1001) 'ay =',ay
        write(*,1002) 'Ntraj =', Ntraj
        write(*,1002) 'Time steps =',kmax
        write(*,1001) 'Time step =',dt
        write(*,1001) 'Mass =',m1,m2
        print *,'Initial Momentum'
        write(*,1001) 'px=',px0
        write(*,1001) 'py=',py0
1001    format(A20,20(f10.6))
1002    format(A20,20(I6))
       
        w=1d0/Ntraj
c	call derivs(x,y,Ntraj,fx,fy,pe)
c        call quantum_potential(m1,m2,x,y,w,Ntraj,qp,qfx,qfy)

c 	initial value for action function s
c	do i=1,Ntraj
c	s(i)=px(i)*(x(i)-x0)+py(i)*(y(i)-y0)
c        enddo
! begin the time step loop
        do 10 k=1,kmax

	  t=t+dt2
	  call derivs(kmax,dt,t,x,y,Ntraj,fx,fy,pe)
          call quantum_potential(m1,m2,x,y,w,Ntraj,qp,qfx,qfy)
! begin the trajectory loop
          do 11 i=1,Ntraj
! half-step increments of moment
            px(i)=px(i)+(fx(i)+qfx(i))*dt2
            py(i)=py(i)+(fy(i)+qfy(i))*dt2

! full step increment of positions
            x(i)=x(i)+px(i)*dt/m1
            y(i)=y(i)+py(i)*dt/m2
11        enddo
        
          t = t+dt2
          call derivs(kmax,dt,t,x,y,Ntraj,fx,fy,pe)
          call quantum_potential(m1,m2,x,y,w,Ntraj,qp,qfx,qfy)

          do 12 i=1,Ntraj
! half-step increments of momenta
  	    px(i)=px(i)+(fx(i)+qfx(i))*dt2
  	    py(i)=py(i)+(fy(i)+qfy(i))*dt2
  	
  ! update potential, kinetic, and total energy for each trajectory	
  	    ke(i)=px(i)**2/(2d0*m1)+py(i)**2/(2d0*m2)
12	  enddo

! calculate action function
	
	do i=1,Ntraj
	  s(i)=s(i)+(ke(i)-pe(i)-qp(i))*dt
	end do

! print out a random trajectory	
        
	write(101,10000) t,(x(i),i=1,40)
	write(104,10000) t,(y(i),i=1,40)

!       calculate the expectation value of x,y 
        call aver(Ntraj,x,w,xav)
        call aver(Ntraj,y,w,yav)
        xx=0d0
        yy=0d0
!        prob=0d0
!       do i=1,Ntraj
!        if (x(i) .lt. -5d0) then
!       prob=prob+w(i)
!        endif
!        xx=x(i)**2*w(i)+xx
!        yy=y(i)**2*w(i)+yy
!        enddo
!      write(102,10000) t,xav,yav,xx,yy
!	write(103,10000) t,prob


! calculate the total energy, the sum of all the trajectories
        call aver(Ntraj,pe,w,po)
        call aver(Ntraj,ke,w,ki)
        call aver(Ntraj,qp,w,qu)
        en=po+ki+qu

        write(100,10000) t,po,ki,qu,en
! end of time loop
10	end do

10000   format(2000(e13.6,1x))   
	end program main
!----------------------------------------------
! average 
!---------------------------------------------
        subroutine aver(Ntraj,x,w,xav)
        implicit real*8(a-h,o-z)
        integer*4,intent(in)::Ntraj
        real*8,intent(in)::x(Ntraj),w(Ntraj)
        real*8,intent(out)::xav 
        
        xav = 0d0
        do i=1,Ntraj
          xav = xav+w(i)*x(i)
        enddo
        
        return
        end subroutine
       

