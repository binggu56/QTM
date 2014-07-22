!        implicit real*8(a-h,o-z)
!        integer*4,intent(in)::m,n,Ndim,Ntraj
!        real*8,intent(out)::fx,v
!        real*8 :: y(Ndim,Ntraj)
!        call ph2(m,Ndim,Ntraj,x,v)
!        dx = 1d-3
!        y = x
!        y(n,m) = y(n,m)+dx
!        call ph2(m,Ndim,Ntraj,y,v1)
!        fx = -(v1-v)/dx
!        
!        return
!        end subroutine        
!-----------------------------------------
! potential for mth trajectory   
!-----------------------------------------    
        subroutine pot(d,m,Ntraj,Ndim,x,v)
        implicit real*8(a-h,o-z)
        integer*4,intent(in)::m,Ntraj,Ndim
        real*8,intent(in) :: x(Ndim,Ntraj),d
        real*8,intent(out) :: v
	real*8 :: f(3),y(27*Ndim)

        Np = Ndim/3
	Npt = 27*Np
        c6 = 12.14d0
        c8 = 215.2d0
        c10 = 4813.9d0
        c9 = 143.1d0
        al = 1.713
        be = 1.5671
        ga = 0.00993d0
        rm = 3.44d0/0.52918d0
        v = 0d0

	f = (/0d0,1d0,-1d0/)
	i = 0
        do k1=1,3
          do k2 = 1,3
            do k3 = 1,3
              i = i+1
              do j=1,Np
                z = x(3*j-2,m)+f(k1)*d
                if(abs(z) > 3d0*d/2d0) z = z-dsign(3d0*d,z)
                y((i-1)*Ndim+3*j-2) = z
                z = x(3*j-1,m)+f(k2)*d
                if(abs(z) > 3d0*d/2d0) z = z-dsign(3d0*d,z)
                y((i-1)*Ndim+3*j-1) = z
                z = x(3*j,m)+f(k3)*d
                if(abs(z) > 3d0*d/2d0) z = z-dsign(3d0*d,z)
                y((i-1)*ndim+3*j) = z
              enddo
            enddo
          enddo
        enddo
	
	do i=1,Npt
	  write(105,*) 'H',y(3*i-2),y(3*i-1),y(3*i)
	enddo

	do i=1,4
          do j=i+1,Npt
              r = dsqrt((y(3*j-2)-y(3*i-2))**2+
     &            (y(3*j-1)-y(3*i-1))**2+
     &            (y(3*j)-y(3*i))**2)
            if(r < 1.28*rm) then
              fc = exp(-(1.28d0*rm/r-1d0)**2)
            else
              fc = 1.d0
            endif
            v = v+exp(al-be*r-ga*r**2)-
     &            (c6/r**6+c8/r**8+c10/r**10)*fc+
     &            (c9/r**9)*fc
          enddo
        enddo
	
1000	format(20(e14.7,1x))
        return
        end subroutine
