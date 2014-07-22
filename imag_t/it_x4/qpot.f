!----quadratic basis f = (x(1),x(2),...,x(Ndim),x(1)**2,...,x(Nd)**2,1)----
	subroutine qpot(am,x,p,w,Ndim,Ntraj,qp,qf,dp,Nb)
        implicit real*8(a-h,o-z)
        integer*4,intent(IN) :: Ntraj,Ndim
        real*8,intent(IN) :: w(Ntraj),x(Ndim,Ntraj),am(Ndim)
        real*8, intent(OUT) :: qf(Ndim,Ntraj),qp(Ntraj),dp(Ndim,Ntraj)
        real*8 :: f(Nb,Ntraj),r(Ndim,Ntraj),p(Ndim,Ntraj)
	real*8 :: M(Nb,Nb),c(Nb,Ndim),s(Ntraj),df(Nb,Ndim),
     &            ddf(Nb,Ndim,Ndim)
  
        call basis(Nb,Ndim,Ntraj,x,f)

!---------c = p*f-------------------
	c=0d0
	do k1=1,Nb
	  do k2=1,Ndim
	    do i=1,Ntraj
	      c(k1,k2)=c(k1,k2)+f(k1,i)*p(k2,i)*w(i)
	    enddo
	  enddo
	enddo

! Matrix M=f x f [Nb,Nb]
        M = 0d0
        do k1=1,Nb
          do k2=1,Nb
            do i=1,Ntraj
              M(k1,k2)=w(i)*f(k1,i)*f(k2,i)+M(k1,k2)
            enddo
          enddo
        enddo

!	write(*,'(A10)') 'M = '
!	do i=1,Nb
!	  write(*,1002) (M(i,j),j=1,Nb)
!	enddo
!	write(*,'(A10)') 'C before = '
!	do i=1,Nb
!	  write(*,1002) (C(i,j),j=1,Ndim)
!	enddo
!1002	format(20f10.6)
	
! calculate matrix c(t)
        call DPOSV('U',Nb,Ndim,M,Nb,c,Nb,INFO)
        if(INFO/=0) then
        print *, "info=",info
        print *, "matrix fails"
        stop
        end if
	
!	write(*,'(A10)') 'C after = '
!	do i=1,Nb
!	  write(*,1002) (C(i,j),j=1,Ndim)
!	enddo
	
!-------qradratic basis---------------------
!        do i=1,Ntraj
!          do j=1,Ndim
!	    do n=1,Nb
!              r(j,i)=c(n,j)*f(n,i)+r(j,i)
!	    enddo
!	  enddo
!	enddo
	if(Nb .eq. (Ndim+1)*(Ndim+2)/2) then

!	do 11 i=1,Ntraj
!	  do 12 j=1,Ndim
! calculate quantum potential
!            qp(i)=-r(j,i)**2/(2d0*am(j))
!     &            -c(j,j)/(2d0*am(j))+qp(i)
!-------------------------------------------
!	    kk = 2*Ndim
!	    do k1=1,Ndim-1
!	      do k2=k1+1,Ndim
!	        kk = kk+1
!	        if(k1 .eq. j) then
!		  qp(i) = qp(i)+c(kk,j)*x(k2,i)/2d0/am(j)
!	          dp(j,i) = dp(j,i)+(c(j,j)+2d0*c(j+Ndim,j)*x(j,i)
!     &                   +c(kk,j)*x(k2,i)
!	        elseif(k2 .eq. j) then
!		  qp(i) = qp(i)+c(kk,j)*x(k1,i)/2d0/am(j)
!	          dp(j,i) = dp(j,i)+c(j,j)+2d0*c(j+Ndim,j)*x(j,i)
!     &                     +c(kk,j)*x(k1,i)
!     	        endif
!
!		qp(i) = qp(i)+c(j,j)/2d0/am(j)+c(j+Ndim)*2d0*x(j,i)
!	      enddo
!	    enddo
	
	call ggf(Ntraj,Ndim,Nb,ddf)

!	do k=1,Nb
!	  write(*,1003) 'DDF =',(ddf(k,i,1),i=1,Ndim)
!	enddo

	qp = 0d0
	do i=1,Ntraj
	  call gradf(x,i,Ntraj,Ndim,Nb,df)

	  do j=1,Nb
	    do k=1,Ndim
	      qp(i) = qp(i)+df(j,k)*c(j,k)/2d0/am(k)
!	      qf(j,i) = qf(j,i)+ddf(j,i)*c(k,j)/2d0/am(j)
	    enddo
	  enddo
	enddo

!	write(*,1003) 'x =',(x(j,Ntraj),j=1,Ndim)
!	do k=1,Nb
!	  write(*,1003) 'DF = ',(df(k,i),i=1,Ndim)
!	enddo
!
!	write(*,1003) 'QPOT = ', (qp(i),i=1,10)

	qf = 0d0	
	do i=1,Ntraj
	  do j=1,Ndim
	    do k=1,Ndim
	      do n=1,Nb  
	        qf(j,i) = qf(j,i)-ddf(n,k,j)*c(n,k)/2d0/am(k)
	      enddo
	    enddo
	  enddo
	enddo

!	write(*,1003) 'QFORCE = ',(qf(1,i),i=1,10)
	


	dp = 0d0
	do i=1,Ntraj
	  call gradf(x,i,Ntraj,Ndim,Nb,df)
	  do j=1,Ndim
	    do n=1,Ndim
	      do k=1,Nb
	        dp(j,i) = dp(j,i)+df(k,j)*c(k,n)*p(n,i)/am(n)
	      enddo
	    enddo
	  enddo
	enddo

!	write(*,1003) 'PdP = ', (dp(1,i),i=1,10)
!-----------------------------------------
!	    qf(j,i) = 2d0*c(j+Ndim,j)
!	    qp(i) = qp(i) + dp(j,i)/2d0/am(j)

!12	  enddo
!11	enddo

	elseif(Nb .eq. Ndim+1) then

! quantum force
!	do i=1,Ntraj
!          do j=1,Ndim
!	    do n=1,Ndim
!	      qf(j,i)=r(n,i)*c(j,n)/am(n)+qf(j,i)
!	    enddo
!          enddo
!	enddo

!-------linear basis----------------------------
	dp = 0d0
	qp = 0d0

	do i=1,Ntraj
	  do j=1,Ndim
	    do k=1,Ndim
	      dp(j,i) = dp(j,i)+p(k,i)*c(j,k)/am(k)
	    enddo
	      qp(i) = qp(i) + c(j,j)/2d0/am(j)
	  enddo
	enddo

	qf = 0d0
!----------------------------------------------
	

!	write(*,1003) 'QPOT = ', (qp(i),i=1,10)
!	write(*,1003) 'QFORCE =', (qf(1,j),j=1,10)
!	write(*,'(A10)') 'PdP = '
!	do i=1,10
!	  write(*,1002) (dp(j,i),j=1,Ndim)
!	enddo
	
	endif

1003 	format(A10,20f10.6) 
        return
        end subroutine

!----------------------------------------
! basis
!----------------------------------------
	subroutine basis(Nb,Ndim,Ntraj,x,f)
	implicit real*8(a-h,o-z)
	integer*4,intent(IN) :: Nb,Ndim,Ntraj
	real*8,intent(IN) :: x(Ndim,Ntraj)
	real*8,intent(OUT):: f(Nb,Ntraj)

	f = 0d0

	if(Nb .eq. Ndim+1) then
! basis vector f = ((x(i),i=1,Ndim),1) for each trajectory
          do i=1,Ntraj
  	    do j=1,Ndim
              f(j,i)=x(j,i)
  	    enddo
              f(Ndim+1,i)=1d0
          enddo
	elseif(Nb .eq. (Ndim+1)*(Ndim+2)/2) then
! basis vector f = (x(i),x(i)**2,x(i)*x(j),1)
	  do 10 i=1,Ntraj
	    do j=1,Ndim
	      f(j,i) = x(j,i)
	    enddo

	    do j=1,Ndim
	      f(j+Ndim,i) = x(j,i)**2
	    enddo

	    mm = 2*Ndim
	    do j=1,Ndim-1
	      do k=j+1,Ndim
	        mm = mm+1
		f(mm,i) = x(j,i)*x(k,i)
	      enddo
	    enddo
	    f(Nb,i) = 1d0

!	    mm = 2*Ndim
!	    do j=1,Ndim
!	      do k1=1,Ndim-1
!	        do k2=k1+1,Ndim
!		  if(k1 .eq. j) then
!		    df(mm,i) = x(k2,i)
!		  elseif(k2 .eq. j) then
!		    df(mm,i) = x(k1,i)
!		  endif
!		enddo
!	      enddo
!	    enddo
10	  enddo
	  
	else
	  write(*,*) "WARNING: ILLEGAL BASIS SIZE!!!"
	endif

	return
	end subroutine
!--------------------------------------------------
! derivative of basis of qradratic basis
!-------------------------------------------------
	subroutine gradf(x,ntr,Ntraj,Ndim,Nb,df)
	implicit real*8(a-h,o-z)
	integer*4,intent(in)::ntr,Ntraj,Ndim,Nb
	real*8,intent(in) :: x(Ndim,Ntraj)
	real*8 :: df(Nb,Ndim)
	
	df = 0d0
	ddf = 0d0
	do i=1,Ndim
 	  df(i,i) = 1d0
	  df(i+Ndim,i) = 2d0*x(i,ntr)
	enddo
	
	do i=1,Ndim
	mm = 2*Ndim
	do k1=1,Ndim-1
	  do k2=k1+1,Ndim
	    mm = mm+1
	    if(k1 .eq. i) then
	      df(mm,i) = x(k2,ntr)
	    elseif(k2 .eq. i) then
	      df(mm,i) = x(k1,ntr)
	    endif
	  enddo
	enddo
	enddo

	return
	end subroutine
!---------------------------------------
! second derivative of quadratic basis
!---------------------------------------
	subroutine ggf(Ntraj,Ndim,Nb,ddf)
	implicit real*8(a-h,o-z)
!	integer*4,intent(in) :: ntr,Ntraj,Ndim,Nb
!	real*8 :: x(Ndim,Ntraj)
	real*8 :: ddf(Nb,Ndim,Ndim)

	ddf = 0d0

	do i=1,Ndim
	  do j=1,Ndim
	    if(j .eq. i) then
	      ddf(i+Ndim,i,j) = 2d0
	    endif

  	    mm = 2*Ndim
	    if(j .ne. i) then
	      do k1=1,Ndim-1
	        do k2=k1+1,Ndim
	          mm = mm+1
	          if(k1.eq.i .and. k2.eq.j) then
	            ddf(mm,i,j) = 1d0
	          elseif(k1.eq.j .and. k2.eq.i) then
	            ddf(mm,i,j) = 1d0
	          endif
	        enddo
	      enddo
	    endif
	  enddo
	enddo
	     
	return
	end subroutine
!-----------------------------------------------
	      
	

	    

