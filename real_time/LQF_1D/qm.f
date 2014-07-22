      program main
      use cdat
C******************************
C     LQF 1D
C****************************** 
      implicit real*8(a-h,o-z)
     	real*8,allocatable :: ke(:),v(:),x(:),y(:),p(:),w(:)
      real*8,allocatable :: qfx(:),qfy(:),qp(:),fx(:),
     &                      ddv(:),s(:)
      complex*16,allocatable :: h(:,:),mat(:,:)
     	real*8 :: f(3),ki
      complex*16,allocatable :: c(:),c0(:),hc(:)
      real :: gasdev
      complex*16 :: psi
      common/wave/al,q0,p0

      open(100,file='energy.dat')
      open(101,file='x.dat')
      open(102,file='aver.dat')  
      open(103,file='wf.dat')
      open(104,file='wf0.dat')

      open(5,file='IN')
      read(5,*) Ntraj
      read(5,*) kmax,dt
      read(5,*) am
      read(5,*) al
      read(5,*) q0
      read(5,*) p0
      read(5,*) idum1
      close(5)

      write(*,1002) 
1002  format(/'1D QTM code with LQF, version 1.0'//)

      allocate(ke(Ntraj),v(Ntraj),x(Ntraj),
     & p(Ntraj),s(Ntraj),w(Ntraj),qfx(Ntraj),
     & qp(Ntraj),fx(Ntraj),ddv(Ntraj))

      allocate(mat(Ntraj,Ntraj))
      allocate(c(Ntraj),c0(Ntraj),h(Ntraj,Ntraj),hc(Ntraj))     

      dt2=dt/2d0
      t=0d0
      Nb = Ntraj

!     grid size
      xmin = -4d0
      xmax = 4d0
      Np = 1000
      dx = (xmax-xmin)/(np-1)

	do i=1,Ntraj
          x(i) = gasdev(idum1)
          x(i)=x(i)/sqrt(4d0*al)+q0
	enddo 
      p = p0   

!     print out the initial conditions        
      write(*,1001) al,Ntraj,kmax,dt,am,p0,q0
1001  format('Initial Conditions'// ,
     +       'ax    = ', f10.6/, 
     +       'Ntraj = ', i6/ , 
     +       'Kmax  = ', i6/ ,
     +       'dt    = ', f10.6/ ,
     +       'Mass  = ', f10.6/ , 
     +       'p0    = ', f10.6/ ,
     +       'q0    = ', f10.6/)
       
      w=1d0/Ntraj


!     begin the time step loop
      do 10 k=1,kmax

        t = t+dt2
!       propogate trajectories with LQF
!       -------------------------------
        call derivs(x,Ntraj,fx,v,ddv)
        call quantum_potential(am,x,w,Ntraj,qp,qfx)

!       half-step increments of moment, full step increment of positions
        do 11 i=1,Ntraj
          p(i)=p(i)+(fx(i)+qfx(i))*dt2
          x(i)=x(i)+p(i)*dt/am
11      enddo
        
        t = t+dt2

        call derivs(x,Ntraj,fx,v,ddv)
        call quantum_potential(am,x,w,Ntraj,qp,qfx)

        do 12 i=1,Ntraj
!         half-step increments of momenta
          p(i)=p(i)+(fx(i)+qfx(i))*dt2
          ke(i)=p(i)**2/(2d0*am)
12      enddo


!     calculate action function
      
!	 do i = 1,Ntraj
!	  s(i)=s(i)+(ke(i)-pe(i)-qp(i))*dt
!	 end do

!       print out a random trajectory	
        write(101,10000) t,(x(i),i=1,20)


!       calculate the expectation values
!       call aver(Ntraj,x,w,xav)
        x2 = 0d0
        do i=1,Ntraj
          x2 = x2+x(i)**2*w(i)
        enddo
          
        write(102,10000) t,x2


!     end of time loop
10    end do

10000 format(2000(e14.7,1x))
      end program main

        
!     ----------------------------------------------
!     average
!     ---------------------------------------------
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
