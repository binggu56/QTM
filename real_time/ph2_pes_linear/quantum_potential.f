        subroutine qpot(am,x,w,Ndim,Ntraj,qp,qf)
        implicit real*8(a-h,o-z)
        integer*4,intent(IN) :: Ntraj,Ndim
        real*8,intent(IN) :: w(Ntraj),x(Ndim,Ntraj),am(Ndim)
        real*8, intent(OUT) :: qf(Ndim,Ntraj),qp(Ntraj)
        real*8 :: f(Ndim+1,Ntraj),r(Ndim,Ntraj)
	real*8 :: s(Ndim+1,Ndim+1),c(Ndim+1,Ndim)
! define c matrix 
	c=0d0
	do j=1,Ndim
	c(j,j)=-0.5d0
	enddo
	
        !c(1,1)=-0.5d0
        !c(1,2)=0d0
        !c(2,1)=0d0
        !c(2,2)=-0.5d0
        !c(3,1)=0d0
        !c(3,2)=0d0

! basis vector f = (1,x(1),x(2),...,x(Ndim)) for each trajectory
        s=0d0
	f=0d0
        do i=1,Ntraj
	  do j=1,Ndim
            f(j,i)=x(j,i)
	  enddo
          f(Ndim+1,i)=1d0
        enddo
!	f(Ndim+1)=1d0

! Matrix S=f X f [Ndim+1,Ndim+1]
        s = 0d0
        do k1=1,Ndim+1
          do k2=1,Ndim+1
            do i=1,Ntraj
              s(k1,k2)=w(i)*f(k1,i)*f(k2,i)+s(k1,k2)
            enddo
          enddo
        enddo

!       	print *,f
!	print *,s 

! calculate matrix c(t)
        call DPOSV('U',Ndim+1,Ndim,s,Ndim+1,c,Ndim+1,INFO)
        if(INFO/=0) then
        print *, "info=",info
        print *, "matrix fails"
        stop
        end if

! the momentum operator r=cf
	r=0d0
	qp=0d0
	qf=0d0

        do i=1,Ntraj
          do j=1,Ndim
	    do n=1,Ndim+1
              r(j,i)=c(n,j)*f(n,i)+r(j,i)
        !r(2,i)=c(1,2)*x(1,i)+c(2,2)*x(2,i)+c(3,2)
	    enddo
! calculate quantum potential
            qp(i)=-r(j,i)**2/(2d0*am(j))
     &            -c(j,j)/(2d0*am(j))+qp(i)
	  enddo
	enddo

! quantum force
	do i=1,Ntraj
          do j=1,Ndim
	    do n=1,Ndim
        !qfx(i)=rx(i)*c(1,1)/m1+ry(i)*c(1,2)/m2
        !qfy(i)=rx(i)*c(2,1)/m1+ry(i)*c(2,2)/m2
	      qf(j,i)=r(n,i)*c(j,n)/am(n)+qf(j,i)
	!qf(2,i)=r(2,i)*c(2,1)/m(1)+r(2,i)*c(2,2)/m(2)
	    enddo
          enddo
	enddo

        return
        end subroutine
