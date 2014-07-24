        program main
!**********************************************************************
!	N dimensional quantum trajectory code with imaginary time
!	test with coupled oscilator for each DOF
!	Eulerian Frame
!*********************************************************************
        implicit real*8(a-h,o-z)
	real*8 :: ki,po,qu,tot
      	real*8 :: t,k1,k2,px,py,ax,ay
	integer*4 :: Ntraj,counter
       	real*8,dimension(:,:),allocatable :: dvdx,
     &  x,p,qf,dp
	integer,dimension(:),allocatable :: idum
        real*8,dimension(:),allocatable :: am,p0,w,alpha,
     &  x0,pe,ke,qp,s,s0
	real :: gasdev

! the initial conditions
        open(100,file='ener.dat')
	open(101,file='px.dat')
	open(102,file='traj4')
	open(103,file='pot.dat')
        open(104,file='weight.dat')
        open(105,file='phase.dat')	
        open(10,file='IN')
	read(10,*) Ntraj
	read(10,*) Ndim
        read(10,*) kmax,dt
	read(10,*) am0
	read(10,*) x00
        read(10,*) a0
	read(10,*) ipot
        read(10,*) iqp

        close(10)

	allocate(dvdx(Ndim,Ntraj),pe(Ntraj),ke(Ntraj),
     &           x(Ndim,Ntraj),p(Ndim,Ntraj),
     &           qp(Ntraj),qf(Ndim,Ntraj),w(Ntraj),dp(Ndim,Ntraj))
	allocate(p0(Ndim),idum(Ndim),alpha(Ndim),am(Ndim),
     &           x0(Ndim))
        allocate(s(Ntraj),s0(Ntraj))
!	read(10,*) (p0(i),i=1,Ndim)
!	read(10,*) (idum(i),i=1,Ndim)
!	read(10,*) (alpha(i),i=1,Ndim)
!	read(10,*) (am(i),i=1,Ndim)
!	read(10,*) (x0(i),i=1,Ndim)


!        dt2=dt/2d0
        t=0d0
	x0 = x00
	am = am0
	alpha = a0

!------------------------------------
! initial grid points 
!------------------------------------
!	do i=1,Ndim
!	  idum(i) = 5 + i
!	enddo

!        do i=1,Ntraj
!	  do j=1,Ndim
!            x(j,i)=gasdev(idum(j))
!            x(j,i)=x(j,i)/sqrt(4d0*alpha(j))+x0(j)
!	  enddo
!	enddo
!--------------------------------------
        write(*,*) 'Initial Conditions'
        write(*,*) 'Ntraj=',Ntraj
	print *,'Ndim=',Ndim
	print *,'kmax=',kmax
	print *,'dt=',dt
	write(*,1001) 'Mass = ',(am(i),i=1,Ndim)
	write(*,1001) 'Alpha = ',(alpha(i),i=1,Ndim)
1001    format(A10,20(f10.6,1x))
	call sample(Ndim,Ntraj,alpha,x,x0)

!--------------------------------------------
!	do i=1,Ntraj
!	  write (*,10000) (x(j,i),j=1,Ndim)
!	enddo

! initial momentum
	s0 = 0d0
	do i=1,Ntraj
	  do j=1,Ndim
            p(j,i)=2d0*alpha(j)*(x(j,i)-x0(j))
	    s0(i) = alpha(j)*(x(j,i)-x0(j))**2 + s0(i)
          enddo
	enddo

	s = 0d0
      w=1d0/Ntraj
	if(iqp .eq. 1) then
	  Nb = Ndim+1
	elseif(iqp .eq. 2) then
	  Nb = (Ndim+1)*(Ndim+2)/2
	else
	  write(*,*) 'ILLEGAL IQP!!!'
	  stop
	endif


	if(ipot .eq. 3) then
	  call derivs3(x,dvdx,Ndim,Ntraj,pe)
	elseif(ipot .eq. 2) then
	  call derivs2(x,dvdx,Ndim,Ntraj,pe)
	elseif(ipot .eq. 1) then
	  call derivs(x,dvdx,Ndim,Ntraj,pe)
	endif

	do i=1,Ntraj
	  write(103,10000) x(1,i),pe(i)
	enddo

	write(*,1001) 'V = ', (pe(i),i=1,10)
! time propagation	
        do 10 counter=1,kmax
! increase t by dt
	  t=t+dt

	  call qpot(am,x,p,w,Ndim,Ntraj,qp,qf,dp,Nb)

!----one-step increment of action function, S(Ntraj)----
    	  do i=1, Ntraj
	    do j=1,Ndim
!----one-step increment of momenta, p(Ndim,Ntraj)----
              p(j,i)=p(j,i)+(dvdx(j,i)-qf(j,i))*dt
     &	             -dp(j,i)*dt
     	    enddo
	  enddo

	  do i=1,Ntraj
	    do j=1,Ndim
	      s(i)=s(i)+(-p(j,i)**2/2d0/am(j))*dt
     	    enddo
	      s(i) = s(i)+(pe(i)+qp(i))*dt
	      w(i) = exp(-2d0*(s(i)))*1d0/Ntraj
	  enddo

!------------------------------------------------------------
! update potential, kinetic, and total energy
          ke = 0d0
          do i=1, Ntraj
            do j=1,Ndim	
	      ke(i)=p(j,i)**2/(2d0*am(j))+ke(i)
            enddo
	  enddo
	
          write(101,10000) t,(p(1,i),i=1,40)
!	  write(102,10000) t,(x(j,4),j=1,Ndim),(x(j,4),j=1,Ndim)

	
	  call aver(Ntraj,w,pe,po)
	  call aver(Ntraj,w,ke,ki)
	  call aver(Ntraj,w,qp,qu)


	  ww = 0d0
	  do i=1,Ntraj
            ww = ww+w(i)
          enddo
          tot=(-ki+po+qu)/ww
          
          write(100,10000) t,po/ww,ki/ww,qu/ww,tot
          write(104,10000) t,(w(i),i=1,10)

10      enddo
          do i=1,200
            write(105,10000) x(1,i),p(1,i),s(i)+s0(i),w(i)/ww
          enddo
          

10000   format(100(e14.7,1x))
	
	print *,'Total Energy =', tot	
	end program 

!----------------------------------------------
! average over trajectories y = sum(w(i)**x(i))
!----------------------------------------------
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
! guassian smapling
!---------------------------------------------
	subroutine sample(Ndim,Ntraj,a0,x,x0)
	implicit real*8(a-h,o-z)
	real*8, intent(in) :: x0(Ndim),a0(Ndim)
     	integer :: idum(Ndim)
	real*8, intent(OUT) :: x(Ndim,Ntraj)
        real gasdev
        
        xcut = 0.2d0

	do i=1,Ndim
	  idum(i) = 6 + i
	enddo

! initial grid points 
        do i=1,Ntraj
	  do j=1,Ndim
100         x(j,i)=gasdev(idum(j))
            x(j,i)=x(j,i)/sqrt(4d0*a0(j))+x0(j)
            if(abs(x(j,i)) > xcut) goto 100
	  enddo
	enddo

	return 
	end subroutine
