        program main
!**********************************************************************
!	N dimensional quantum trajectory code with linear quantum force
!	test with coupled oscilator for each DOF (3 dim)
!*********************************************************************
        implicit real*8(a-h,o-z)
	real*8 :: ki,po,qu,tot
      	real*8 :: t,k1,k2,px,py,ax,ay
	integer*4 :: Ntraj,counter
	real*8 :: dt,dt2,s
       	real*8,dimension(:,:),allocatable :: dvdx,
     &  x,p,qf
	integer,dimension(:),allocatable :: idum
        real*8,dimension(:),allocatable :: am,p0,w,alpha,
     &  x0,pe,ke,qp,cf
	real :: gasdev

! the initial conditions
        open(100,file='output')
	open(101,file='xoutput')
	open(102,file='traj4')
	
        open(10,file='IN')
	read(10,*) Ntraj
	read(10,*) Ndim
        read(10,*) kmax,dt

        close(10)

	allocate(dvdx(Ndim,Ntraj),pe(Ntraj),ke(Ntraj),
     &           x(Ndim,Ntraj),p(Ndim,Ntraj),
     &           qp(Ntraj),qf(Ndim,Ntraj),w(Ntraj))
	allocate(p0(Ndim),idum(Ndim),alpha(Ndim),am(Ndim),
     &           x0(Ndim),cf(Ndim)) 
!	read(10,*) (p0(i),i=1,Ndim)
!	read(10,*) (idum(i),i=1,Ndim)
!	read(10,*) (alpha(i),i=1,Ndim)
!	read(10,*) (am(i),i=1,Ndim)
!	read(10,*) (x0(i),i=1,Ndim)


        dt2=dt/2d0
        t=0d0
	p0 = 0d0
	x0 = 0d0
	cf = 0.0d0
	am = 1d0
	alpha = 0.5d0

	call seed(idum,Ndim)


!       write(*,*) 'Quantum Trajectories'
!	print *,x0

        write(*,*) 'Initial Conditions'
        write(*,*) 'Ntraj=',Ntraj
	print *,'Ndim=',Ndim
	print *,'kmax=',kmax
	print *,'dt=',dt

! initial grid points 
        do i=1,Ntraj
	  do j=1,Ndim
            x(j,i)=gasdev(idum(j))
            x(j,i)=x(j,i)/sqrt(4d0*alpha(j))+x0(j)
	  enddo
	enddo
!--------------------------------------------
!	do i=1,Ntraj
!	  write (*,10000) (x(j,i),j=1,Ndim)
!	enddo

! initial momentum 
	do i=1,Ntraj
	  do j=1,Ndim
            p(j,i)=p0(j)
          enddo
	enddo
!-----------------------------------------	
        w=1d0/Ntraj


! time propagation	
        do 10 counter=1,kmax
! increase t by dt
	  t=t+dt	
	  call derivs(x,x0,dvdx,Ndim,Ntraj,pe)
	  call qpot(am,x,w,Ndim,Ntraj,qp,qf)

! half-step increments of momenta & full step increment of positions
	  do i=1,Ntraj
	    do j=1,Ndim
              p(j,i)=p(j,i)+(-dvdx(j,i)+qf(j,i)-cf(j)*p(j,i))*dt2
              x(j,i)=x(j,i)+p(j,i)*dt/am(j)
     	    enddo   
          enddo 
    
          call derivs(x,x0,dvdx,Ndim,Ntraj,pe)
	  call qpot(am,x,w,Ndim,Ntraj,qp,qf)

! half-step increments of momenta
          do i=1,Ntraj
	    do j=1,Ndim
	      p(j,i)=p(j,i)+(-dvdx(j,i)+qf(j,i)-cf(j)*p(j,i))*dt2
	    enddo
          enddo
!------------------------------------------------------------
! update potential, kinetic, and total energy
          ke = 0d0
          do i=1, Ntraj
            do j=1,Ndim	
	      ke(i)=p(j,i)**2/(2d0*am(j))+ke(i)
            enddo
	  enddo
	
          write(101,10000) t,(x(1,i),i=1,40),(x(2,i),i=1,40)
	  write(102,10000) t,(x(j,4),j=1,Ndim),(x(j,4),j=1,Ndim)

	
	  call aver(Ntraj,w,pe,po)
	  call aver(Ntraj,w,ke,ki)
	  call aver(Ntraj,w,qp,qu)

	  tot=po+ki+qu
	  s=0d0
	  write(100,10000) t,po,ki,qu,tot

10      enddo

10000   format(100(e14.7,1x))
	
	print *,'Total Energy =', tot	
	end program 

!---------------------------------------------
! average over trajectories y = sum(w(i)*x(i)) 
!---------------------------------------------
	subroutine aver(Ntraj,w,x,y)
	implicit real*8(a-h,o-z)
	real*8 :: x(Ntraj),w(Ntraj)

	y = 0d0

	do i=1,Ntraj
	  y = y+x(i)*w(i) 
	enddo
	
	return
	end subroutine
!----------------------------------------------
! assign idum[Ndim]
	subroutine seed(idum,Ndim)
	implicit real*8(a-h,o-z)
	integer*4, intent(IN) :: Ndim
	integer*4, intent(OUT) :: idum(Ndim)
	do i=1,Ndim
	  idum(i) = 5 + i
	enddo

	return
	end subroutine
!----------------------------------------------

