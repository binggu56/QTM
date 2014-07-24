        program main
!**********************************************************************
!	N dimensional quantum trajectory code with imaginary time
!	test with coupled oscilator for each DOF
!	Eulerian Frame
!*********************************************************************
        implicit real*8(a-h,o-z)
	real*8 :: ki,po,qu,tot
      	real*8 :: t,k1,k2
	integer*4 :: Ntraj,counter
       	real*8,dimension(:,:),allocatable :: dvdx,
     &  x,p,qf,dp
	integer,dimension(:),allocatable :: idum
        real*8,dimension(:),allocatable :: am,p0,w,alpha,
     &  x0,pe,ke,qp,s,s0,w0
        character*4,dimension(:),allocatable :: atom
! the initial conditions
        open(100,file='energy.dat')
	open(101,file='px.dat')
	open(102,file='traj4')
	open(103,file='pot.dat')
        open(104,file='samp.dat')
        open(105,file='weight.dat')

        open(10,file='IN')
	read(10,*) Ntraj
	read(10,*) Ndim
        read(10,*) kmax,dt
!	read(10,*) am0
!	read(10,*) x00
!        read(10,*) a0
!	read(10,*) ipot
        read(10,*) iqp

        close(10)
        
        np = ndim/3

	allocate(dvdx(Ndim,Ntraj),pe(Ntraj),ke(Ntraj),
     &           x(Ndim,Ntraj),p(Ndim,Ntraj),
     &           qp(Ntraj),qf(Ndim,Ntraj),w(Ntraj),dp(Ndim,Ntraj))
	allocate(p0(Ndim),idum(Ndim),alpha(Ndim),am(Ndim),
     &           x0(Ndim))
        allocate(s(Ntraj),s0(Ntraj),w0(Ntraj),atom(Np))
!	read(10,*) (p0(i),i=1,Ndim)
!	read(10,*) (idum(i),i=1,Ndim)
!	read(10,*) (alpha(i),i=1,Ndim)
!	read(10,*) (am(i),i=1,Ndim)
!	read(10,*) (x0(i),i=1,Ndim)


!        dt2=dt/2d0
        t=0d0
	am = 1836.18d0
	alpha = 8d0

        open(11,file='h2.xyz')
        do i=1,np
          read(11,*) atom(i),x0(3*i-2),x0(3*i-1),x0(3*i)
          write(*,*) atom(i),x0(3*i-2),x0(3*i-1),x0(3*i)
        enddo
        close(11)

        

!------------------------------------
! initial grid points 
!------------------------------------
!	do i=1,Ndim
!	  idum(i) = 5 + i
!	enddo
!
!        do i=1,Ntraj
!	  do j=1,Ndim
!            x(j,i)=gasdev(idum(j))
!            x(j,i)=x(j,i)/sqrt(4d0*alpha(j))+x0(j)
!	  enddo
!	enddo
        call gsamp(Ntraj,Ndim,alpha,x0,x,w)
!        call usamp(Ntraj,Ndim,alpha,x0,x,w)
        
        w0 = w
!--------------------------------------
        write(*,*) 'Initial Conditions'
        write(*,1002) 'Ntraj=',Ntraj
	write(*,1002) 'Ndim=',Ndim
	write(*,1002) 'kmax=',kmax
	write(*,'(A10,f10.6)') 'dt=',dt
	write(*,*) 'Mass = ',(am(i),i=1,Ndim)
	write(*,1001) 'Alpha = ',(alpha(i),i=1,Ndim)
1001    format(A10,20(f10.6,1x))
1002    format(A10,I10)

!--------------------------------------------
	do i=1,Ntraj
	  write (104,10000) (x(j,i),j=1,Ndim)
	enddo

! initial momentum
	s0 = 0d0
	do i=1,Ntraj
	  do j=1,Ndim
            p(j,i)=2d0*alpha(j)*(x(j,i)-x0(j))
	    s0(i) = alpha(j)*(x(j,i)-x0(j))**2 + s0(i)
          enddo
	enddo

	s = 0d0
!        w=1d0/Ntraj
	if(iqp .eq. 1) then
	  Nb = Ndim+1
	elseif(iqp .eq. 2) then
	  Nb = (Ndim+1)*(Ndim+2)/2
	else
	  write(*,*) 'ILLEGAL IQP!!!'
	  stop
	endif


!	if(ipot .eq. 3) then
!	  call derivs3(x,dvdx,Ndim,Ntraj,pe)
!	elseif(ipot .eq. 2) then
!	  call derivs2(x,dvdx,Ndim,Ntraj,pe)
!	elseif(ipot .eq. 1) then
!	  call derivs(x,x0,dvdx,Ndim,Ntraj,pe)
!	endif


        call derivs(Ndim,Ntraj,x,pe,dvdx)

	do i=1,Ntraj
	  write(103,10000) dsqrt((x(4,i)-x(1,i))**2+(x(5,i)-x(2,i))**2+
     &                   (x(6,i)-x(3,i))**2),pe(i)
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
	      w(i) = exp(-2d0*(s(i)))*w0(i)
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
10      enddo
        
        do i=1,Ntraj
          write(105,10000) x(1,i),p(1,i),w(i)/ww
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
	subroutine gsamp(Ntraj,Ndim,alpha,x0,x,w)
	implicit real*8(a-h,o-z)
	real*8,intent(in)::alpha(Ndim),x0(Ndim)
     	integer*4 :: idum(Ndim)
	real*8,intent(OUT) :: x(Ndim,Ntraj),w(Ntraj)
        real :: gasdev

	do i=1,Ndim
	  idum(i) = 6 + i
	enddo

! initial grid points 
        do i=1,Ntraj
	  do j=1,Ndim
            x(j,i)=gasdev(idum(j))
            x(j,i)=x(j,i)/sqrt(4d0*alpha(j))+x0(j)
	  enddo
	enddo

        w = 1d0/Ntraj

	return 
	end subroutine
!-----------------------------------------------
! uniform smapling
!----------------------------------------------
        subroutine usamp(Ntraj,Ndim,alpha,x0,x,w)
        implicit real*8(a-h,o-z)
        integer*4,intent(in)::Ntraj,Ndim
        real*8,intent(in)::alpha(Ndim),x0(Ndim)
        real*8,intent(out)::x(Ndim,Ntraj),w(Ntraj)
        real*8::dx(Ndim),xmin(Ndim),xmax(Ndim),s(Ntraj)
        
!        cut = 1d-6
        pi = 4d0*atan(1.0d0)
        xmin = -1.5d0
        xmax = 1.5d0

	do i=1,Ndim
!	  idum(i) = 5 + i
!          dx(i) = sqrt(-log(cut*sqrt(alpha(i)/pi))
!     &            /2d0/alpha(i))
          dx(i) = (xmax(i)-xmin(i))/Ntraj
	enddo
        
        w = 1d0
        s = 0d0

        do i=1,Ntraj
          do j=1,Ndim
!            x(j,i) = ran2(idum(j))
!            x(j,i) = 2d0*x(j,i)-1d0
!            x(j,i) = x0(j)+x(j,i)*dx(j)
            x(j,i) = xmin(j)+(i-1)*dx(j)
            s(i) = s(i)+alpha(j)*(x(j,i)-x0(j))**2
!     &             dx(j)*dsqrt(2d0*alpha(j)/pi)
          enddo
        enddo
        
        wt = 0d0
        do i=1,Ntraj
          w(i) = exp(-2d0*s(i))
          wt = wt+w(i)
        enddo

        do i=1,Ntraj
          w(i) = w(i)/wt
        enddo
                
        return
        end subroutine
!----------------------------------------------
! uniform random sampling between [0..1]
!---------------------------------------------
      FUNCTION ran1(idum)
C-----------------------------------------------------------------
C    DOUBLE PRECISION VERSION OF THE NUMREC ROUTINE
C    FRANK GROSSMANN, 12.9.1994
C-----------------------------------------------------------------
      INTEGER*4 idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8 ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1.d0/IM,IQ=127773,IR=2836,
     *NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2d-14,RNMX=1.0d0-EPS)
      INTEGER*4 j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/
      if (idum.le.0.or.iy.eq.0) then
        idum=max(-idum,1)
        do 11 j=NTAB+8,1,-1
          k=idum/IQ
          idum=IA*(idum-k*IQ)-IR*k
          if (idum.lt.0) idum=idum+IM
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ
      idum=IA*(idum-k*IQ)-IR*k
      if (idum.lt.0) idum=idum+IM
      j=1+iy/NDIV
      iy=iv(j)
      iv(j)=idum
      ran1=min(AM*iy,RNMX)
      return
      END
!-----------------------------------------------
