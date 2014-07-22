        program main 
        implicit none
	real*8 :: ki,po,qu,tot,x0,y0
      	real*8 :: t,k1,k2,m1,m2,xav,yav,lambda
	integer*8 :: Ntraj
     	real*8,allocatable :: ke(:),pe(:),x(:),y(:),px(:),py(:),w(:)
        real*8,allocatable ::  qfx(:),qfy(:),qp(:),fx(:),
     & fy(:),s(:)
     	real*8 :: dv,dt,dt2,f(3),qx0,qy0,prob,hbar
      real*8 cor_d,an0,pi
	integer :: k,kmax,i,j,ipot
	complex*16 :: im,auc
      real*8 ::ax,ay,px0,py0,xx,yy
	real :: gasdev
      integer :: idum1,idum2
      COMMON/PES/ipot
        open(100,file='output')
        open(101,file='xoutput')
        open(102,file='xav')
	open(103,file='prob')
	open(104,file='youtput') 
      open(105,file='ini')
      open(106,file='cor')
      open(107,file='cor2')
        open(5,file='IN')
	read(5,*) Ntraj
      read(5,*) kmax,dt
        read(5,*) m1,m2
	read(5,*) idum1,idum2
        read(5,*) ax,ay
        read(5,*) qx0,qy0
        read(5,*) px0,py0
        read(5,*) lambda
      read(5,*) ipot 
        close(5)
        write(*,*) 'Quantum Trajectories'

	allocate(ke(Ntraj),pe(Ntraj),x(Ntraj),
     & y(Ntraj),px(Ntraj),py(Ntraj),s(Ntraj),
     & w(Ntraj),qfx(Ntraj),qfy(Ntraj),
     & qp(Ntraj),fx(Ntraj),fy(Ntraj))

        dt2=dt/2d0
        t=0d0
      hbar=1d0
      PI=4.0*atan(1d0)
      an0=pi/2d0/dsqrt(ax*ay)
      an0=1d0/dsqrt(an0)
	im=(0d0,1d0)
C        m1 = m1*2d0*m2/(m1+2d0*m2)
C        m2 = m2/2d0
C        m1 = m1*1836d0
C        m2 = m2*1836d0
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
        
c 	center for the potential energy
	x0=0d0
	y0=0d0

c        print out the initial conditions        
        write(*,*) 'Initial Conditions'
        print *,'ax =',ax
        print *,'ay =',ay
        write(*,*) 'Number of trajectories =', Ntraj
        write(*,*) 'Time steps =',kmax,dt
        write(*,*) 'Mass =',m1,m2
        print *,'Initial Momentum'
        print *,'px=',px0
        print *,'py=',py0
      print *,'qx0=',qx0
      print *,'qy0=',qy0
	print *,'lambda=',lambda
       
        w=1d0/Ntraj
c	call derivs(x,y,Ntraj,fx,fy,pe)
c        call quantum_potential(m1,m2,x,y,w,Ntraj,qp,qfx,qfy)

c 	initial value for action function s
	do i=1,Ntraj
	s(i)=px(i)*(x(i)-qx0)+py(i)*(y(i)-qy0)
      enddo
! begin the time step loop
        do 10 k=1,kmax

! increase t by dt
	t=t+dt	


	call derivs(lambda,x,y,Ntraj,fx,fy,pe)
        call quantum_potential(m1,m2,x,y,w,Ntraj,qp,qfx,qfy)
! begin the trajectory loop
        do 11 i=1,Ntraj
! half-step increments of moment
        px(i)=px(i)+(fx(i)+qfx(i))*dt2
        py(i)=py(i)+(fy(i)+qfy(i))*dt2

! full step increment of positions
        x(i)=x(i)+px(i)*dt/m1
        y(i)=y(i)+py(i)*dt/m2
11      end do
        

        call derivs(lambda,x,y,Ntraj,fx,fy,pe)
        call quantum_potential(m1,m2,x,y,w,Ntraj,qp,qfx,qfy)

        do 12 i=1,Ntraj
! half-step increments of momenta
	px(i)=px(i)+(fx(i)+qfx(i))*dt2
	py(i)=py(i)+(fy(i)+qfy(i))*dt2
	

! update potential, kinetic, and total energy for each trajectory	
        

	ke(i)=px(i)**2/(2d0*m1)+py(i)**2/(2d0*m2)

! end of trajectory loop
12	end do
! calculate action function
	
	do i=1,Ntraj
	s(i)=s(i)+(ke(i)-pe(i)-qp(i))*dt
	end do

! print out a random trajectory	
        
	write(101,10000) t,(x(i),i=1,40)
	write(104,10000) t,(y(i),i=1,40)

c       calculate the expectation value of x,y 
        xav=0d0
        yav=0d0
        xx=0d0
        yy=0d0
        prob=0d0
      cor_d=0d0
      auc=(0d0,0d0)
        do i=1,Ntraj
        if (x(i) .lt. -5d0) then
        prob=prob+w(i)
        endif
        xav=w(i)*x(i)+xav
        yav=w(i)*y(i)+yav
        xx=x(i)**2*w(i)+xx
        yy=y(i)**2*w(i)+yy
      auc=auc+w(i)*exp(2d0*im*s(i))
      cor_d=cor_d+w(i)*abs(an0*exp(-ax*(x(i)-qx0)**2+
     &      im*px(i)*(x(i)-qx0)/hbar
     &      -ay*(y(i)-qy0)**2+im*py(i)*(y(i)-qy0)/hbar))**2
        enddo
      write(106,10000) t,cor_d
      write(107,10000) 2d0*t,auc
      write(102,10000) t,xav,yav,xx,yy
	write(103,10000) t,prob


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

        

