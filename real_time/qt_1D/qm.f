        program main
C******************************
C Friction 1-D
C****************************** 
        implicit real*8(a-h,o-z)
	real*8 :: ki,po,qu,tot,x0,y0
      	real*8 :: t,k1,k2,m1,m2,xav,yav
     	real*8,allocatable :: ke(:),pe(:),x(:),y(:),px(:),py(:),w(:)
        real*8,allocatable ::  qfx(:),qfy(:),qp(:),fx(:),
     & fy(:),s(:)
     	real*8 :: dv,dt,dt2,f(3),qx0,qy0,prob
	complex*16 :: im
	real :: gasdev

        open(100,file='energy.dat')
        open(101,file='xoutput')
        open(102,file='xav')  
	open(103,file='prob')

        open(5,file='IN')
	read(5,*) Ntraj
        read(5,*) kmax,dt
        read(5,*) m1
	read(5,*) idum1
        read(5,*) ax
        read(5,*) qx0
        read(5,*) px0
        read(5,*) gamma 
        close(5)
        write(*,*) 'Quantum Trajectories'

	allocate(ke(Ntraj),pe(Ntraj),x(Ntraj),
     & y(Ntraj),px(Ntraj),py(Ntraj),s(Ntraj),
     & w(Ntraj),qfx(Ntraj),qfy(Ntraj),
     & qp(Ntraj),fx(Ntraj),fy(Ntraj))

        dt2=dt/2d0
        t=0d0
	im=(0d0,1d0)
C	m1 = m1*2d0*m2/(m1+2d0*m2)
C	m2 = m2/2d0
C	m1 = m1*1836d0
C	m2 = m2*1836d0
        !open(103,file='xyoutput')
	do i=1,Ntraj
          x(i)=gasdev(idum1)
          x(i)=x(i)/sqrt(4d0*ax)+qx0
	  !write (103,*) x(i),y(i)
	enddo 
       
        px = px0
        
!       center for the potential energy
	x0=0d0

c        print out the initial conditions        
        write(*,*) 'Initial Conditions'
        print *,'ax =',ax
        write(*,*) 'Number of trajectories =', Ntraj
        write(*,*) 'Time steps =',kmax,dt
        write(*,*) 'Mass =',m1,m2
        print *,'Initial Momentum'
        print *,'px=',px0
	print *,'gamma=',gamma
       
        w=1d0/Ntraj
!	call derivs(x,y,Ntraj,fx,fy,pe)
!        call quantum_potential(m1,m2,x,y,w,Ntraj,qp,qfx,qfy)

! 	initial value for action function s
        aver_s = 0d0    
	do i=1,Ntraj
	  s(i)=px(i)*(x(i)-qx0)
          aver_s=aver_s+w(i)*s(i)
        enddo

! begin the time step loop
        do 10 k=1,kmax

! increase t by dt
	t=t+dt2
! UPDATE AVERAGE ACTION FUNCTION
!	aver_s = 0d0	
!	do i=1,Ntraj
!	  aver_s=aver_s+w(i)*s(i)
!	enddo
        
	call derivs(kmax,dt,t,x,Ntraj,fx,pe)
        call quantum_potential(m1,x,w,Ntraj,qp,qfx)
! begin the trajectory loop
        do 11 i=1,Ntraj
! half-step increments of moment
          px(i)=px(i)+(fx(i)+qfx(i)-gamma*px(i))*dt2
! full step increment of positions
          x(i)=x(i)+px(i)*dt/m1
11      enddo
        
        t = t+dt2
        call derivs(kmax,dt,t,x,Ntraj,fx,pe)
        call quantum_potential(m1,x,w,Ntraj,qp,qfx)

        do 12 i=1,Ntraj
! half-step increments of momenta
	px(i)=px(i)+(fx(i)+qfx(i))*dt2
! update potential, kinetic, and total energy for each trajectory	
	ke(i)=px(i)**2/(2d0*m1)
12	end do

! calculate action function
	
	do i = 1,Ntraj
	s(i)=s(i)+(ke(i)-pe(i)-qp(i)-gamma*(s(i)-aver_s))*dt
	end do

        aver_s = 0d0    
        do i=1,Ntraj
          aver_s=aver_s+w(i)*s(i)
        enddo

! print out a random trajectory	
        
	write(101,10000) t,(x(i),i=1,40)

c       calculate the expectation value of x,y 
        xav=0d0
        xx=0d0
        prob=0d0
        do i=1,Ntraj
        xav=w(i)*x(i)+xav
        xx=x(i)**2*w(i)+xx
        enddo
        write(102,10000) t,xav,xx
	write(103,10000) t,aver_s


! calculate the total energy, the sum of all the trajectories	
	po=0d0
	ki=0d0
        qu=0d0
	do i=1,Ntraj
        po=po+pe(i)*w(i)
        ki=ki+ke(i)*w(i)
        qu=qu+qp(i)*w(i)
        end do
        tot=po+ki+qu
        write(100,10000) t,po,ki,qu,tot
10000   format(2000(e13.6,2x))   

! end of time loop
10	end do
		
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

