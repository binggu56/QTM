!----linear basis f = (x(1),x(2),...,x(Ndim),1)----
	subroutine qpot(am,x,p,w,Ndim,Ntraj,qp,qf,dp,Nb)
        implicit real*8(a-h,o-z)
        integer*4,intent(IN) :: Ntraj,Ndim
        real*8,intent(IN) :: w(Ntraj),x(Ndim,Ntraj),am(Ndim)
        real*8, intent(OUT) :: qf(Ndim,Ntraj),qp(Ntraj),dp(Ndim,Ntraj)
        real*8 :: f(Ndim+1,Ntraj),r(Ndim,Ntraj),p(Ndim,Ntraj)
	real*8 :: M(Ndim+1,Ndim+1),c(Ndim+1,Ndim),s(Ntraj)
! define c matrix 
!	do j=1,Ndim
!	c(j,j)=-0.5d0
!	enddo

!-------choose which quantum potential to use----------
!	if(iqpot .eq. 1) then
!	Nf = Ndim+1
!	endif
!	

! basis vector f = (x(1),x(2),...,x(Ndim),1) for each trajectory
	f=0d0
        do i=1,Ntraj
	  do j=1,Ndim
            f(j,i)=x(j,i)
	  enddo
          f(Ndim+1,i)=1d0
        enddo

!---------c = p*f-------------------
	c=0d0
	do k1=1,Ndim+1
	  do k2=1,Ndim
	    do i=1,Ntraj
	      c(k1,k2)=c(k1,k2)+p(k2,i)*f(k1,i)*w(i)
	    enddo
	  enddo
	enddo

! Matrix M=f X f [Ndim+1,Ndim+1]
        M = 0d0
        do k1=1,Ndim+1
          do k2=1,Ndim+1
            do i=1,Ntraj
              M(k1,k2)=w(i)*f(k1,i)*f(k2,i)+M(k1,k2)
            enddo
          enddo
        enddo
!-------print out matrix equation---------
        write(*,*) 'M = '
        do i=1,Nb
          write(*,1002) (M(i,j),j=1,Nb)
        enddo
	write(*,*) 'C='
        do i=1,Nb
          write(*,1002) (C(i,j),j=1,Ndim)
        enddo
1002    format(20f10.6)

! calculate iatrix c(t)
        call DPOSV('U',Ndim+1,Ndim,M,Ndim+1,c,Ndim+1,INFO)
        if(INFO/=0) then
        print *, "info=",info
        print *, "matrix fails"
        stop
        end if

	! the momentum operator r=cf
	r=0d0
	qp=0d0
	qf=0d0

!        do i=1,Ntraj
!          do j=1,Ndim
!	    do n=1,Ndim+1
!              r(j,i)=c(n,j)*f(n,i)+r(j,i)
!        !r(2,i)=c(1,2)*x(1,i)+c(2,2)*x(2,i)+c(3,2)
!	    enddo
!! calculate quantum potential
!            qp(i)=-r(j,i)**2/(2d0*am(j))
!     &            -c(j,j)/(2d0*am(j))+qp(i)
!	  enddo
!	enddo

	do i=1,Ntraj
	  do j=1,Ndim
	      dp(j,i) = c(j,j)
	      qp(i) = qp(i) + dp(j,i)/2d0/am(j)
	  enddo
	enddo


! quantum force
!	do i=1,Ntraj
!          do j=1,Ndim
!	    do n=1,Ndim
!	      qf(j,i)=r(n,i)*c(j,n)/am(n)+qf(j,i)
!	    enddo
!          enddo
!	enddo
	qf = 0d0

        return
        end subroutine
