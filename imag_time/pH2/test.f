        
	program main
	implicit real*8(a-h,o-z)
        parameter(Ndim=1,Ntraj=100)
	dimension idum(Ndim),x(Ndim,Ntraj)
        REAL ran2
!	integer*4, parameter :: Ntraj=10
!	real*8 :: x(Ntraj),y(Ntraj)
!        x(1)=0
!        y(1)=0
!        do i=1,Ntraj-1
	
!       x(i+1)=x(i)+1d0/Ntraj
!        y(i+1)=y(i)+1d0/Ntraj
       
	
!        real*8 :: gasdev,random(100) 
        open(1000,file='testoutput')
        do i=1,Ndim   
          idum=5+i
        enddo

        do i=1,Ntraj
          do j=1,Ndim
            x(j,i)=ran2(idum(j))
          enddo
            write(1000,*) i,(x(j,i),j=1,Ndim)        
        enddo

        return
	end program main
