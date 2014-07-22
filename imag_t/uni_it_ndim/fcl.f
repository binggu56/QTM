	subroutine derivs(x,x0,dvdx,Ndim,Ntraj,pe)
        implicit real*8(a-h,o-z)
        integer*4,intent(IN) :: Ntraj,Ndim
        real*8,intent(IN) :: x(Ndim,Ntraj),x0(Ndim)
        real*8,intent(OUT) :: dvdx(Ndim,Ntraj),pe(Ntraj)
	real*8 :: hess(Ndim,Ndim)

!	call hessian(hess,Ndim)

        pi=4.d0*atan(1.d0)
        pe=0d0
	dvdx = 0d0

        do i=1,Ntraj
              dvdx(1,i) = 4d0*x(1,i)**3
              pe(i) = x(1,i)**4
        enddo

        return
        end subroutine


!------------------------------------------------------
! Morse oscilator 
!------------------------------------------------------
	subroutine derivs2(x,dvdx,Ndim,Ntraj,pe)
        implicit real*8(a-h,o-z)
        integer*4,intent(IN) :: Ntraj,Ndim
        real*8,intent(IN) :: x(Ndim,Ntraj)
        real*8,intent(OUT) :: dvdx(Ndim,Ntraj),pe(Ntraj)
        real*8 :: xe(Ndim),hess(Ndim,Ndim)
        pi=4.d0*atan(1.d0)
        pe=0d0
        dvdx = 0d0
        xe = 0d0
        De = 0.176d0
        a0 = 1.02
        
	call hessian(hess,Ndim)

        do i=1,Ntraj
          do j=1,Ndim
 	    do k=1,Ndim
              dvdx(j,i) = dvdx(j,i)+hess(j,k)*2d0*De*a0*
     &                    exp(-a0*(x(k,i)-xe(k)))
     &                   *(1d0-exp(-a0*(x(k,i)-xe(k))))
              pe(i) = pe(i)+De*(1d0-exp(-a0*(x(k,i)-xe(k))))*hess(k,j)
     &               *(1d0-exp(-a0*(x(j,i)-xe(j))))
     	    enddo
          enddo
        enddo

        return
        end subroutine
!----------------------------------------------
	subroutine derivs3(x,dvdx,Ndim,Ntraj,pe)
	implicit real*8(a-h,o-z)
	real*8,intent(in) :: x(Ndim,Ntraj)
	real*8,intent(out):: dvdx(Ndim,Ntraj),pe(Ntraj)
	real*8 :: hess(Ndim,Ndim),xe(Ndim)
	
	xe = 0d0
	dvdx = 0d0
	pe = 0d0

	do i=1,Ntraj
	  do j=1,Ndim
	    dvdx(j,i) = 0.36622d0*x(j,i)-0.5603178240*x(j,i)**2
	    pe(i) = pe(i)+0.18311d0*x(j,i)**2
     &             -0.186773d0*x(j,i)**3
!	    dvdx(j,i) = 0.36622d0*x(j,i)
!	    pe(i) = pe(i)+0.18311d0*x(j,i)**2
	  enddo
	enddo


	return 
	end subroutine
!--------------------------------------
! hessian matrix
!--------------------------------------
	subroutine hessian(hess,Ndim)
	implicit real*8(a-h,o-z)
	real*8 :: hess(Ndim,Ndim)

	hess = 0d0
!	hess(1,1) = 1d0
!	hess(1,2) = 0.4d0
!	hess(2,2) = 1d0
!	hess(1,3) = 0.2d0
!	hess(2,3) = 0.4d0
!	hess(3,3) = 1d0
	do i=1,Ndim
	 hess(i,i) = 1d0
	enddo

	do i=1,Ndim
	  do j=1,Ndim
	    if( j .ne. i) then
	      hess(i,j) = 0.d0
	    endif
	  enddo
	enddo

! symmetric matrix
!	do i=1,Ndim
!	  do j=1,Ndim
!	    if(j < i) then
!	      hess(i,j) = hess(j,i)
!	    endif
!	  enddo
!	enddo

	return
	end subroutine
