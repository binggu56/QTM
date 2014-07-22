        subroutine derivs(Ndim,Ntraj,x,v,dv)
        implicit real*8(a-h,o-z)
        integer*4,intent(in)::Ndim,Ntraj
        real*8,intent(in)::x(Ndim,Ntraj)
        real*8,intent(out)::v(Ntraj),dv(Ndim,Ntraj)
        real*8 :: y(Ndim,Ntraj)
       
        v = 0d0
        dv = 0d0
        y = x
        dx = 5d-4
       

        do m=1,Ntraj
          call pot(m,Ntraj,Ndim,x,v0)
          v(m) = v0
          do n=1,Ndim
            y(n,m) = y(n,m)+dx
!            do i=1,Np
!              do j=i+1,Np
!                rx = dsqrt((x(3*j-2,m)-x(3*i-2,m))**2+
!     &                 (x(3*j-1,m)-x(3*i-1,m))**2+
!     &                 (x(3*j,m)-x(3*i,m))**2)
!                ry = dsqrt((y(3*j-2,m)-x(3*i-2,m))**2+
!                       (y(3*j-1,m)-y(3*i-1))**2+
!                       (y(3*j,m)-y(3*i,m))**2)
            call pot(m,Ntraj,Ndim,y,v1)  
            dv(n,m) = (v1-v0)/dx
            y(n,m) = x(n,m)
          enddo
        enddo  
    
        return
        end subroutine
        
!        subroutine cforce(m,n,Ndim,Ntraj,x,fx,v)
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
        subroutine pot(m,Ntraj,Ndim,x,v)
        implicit real*8(a-h,o-z)
        integer*4,intent(in)::m,Ntraj,Ndim
        real*8,intent(in) :: x(Ndim,Ntraj)
        real*8,intent(out) :: v

        Np = Ndim/3
        c6 = 12.14d0
        c8 = 215.2d0
        c10 = 4813.9d0
        c9 = 143.1d0
        al = 1.713
        be = 1.5671
        ga = 0.00993d0
        rm = 3.44d0/0.52918d0
        v = 0d0
        
        do i=1,Np
          do j=i+1,Np
            r = dsqrt((x(3*j-2,m)-x(3*i-2,m))**2+
     &        (x(3*j-1,m)-x(3*i-1,m))**2+
     &        (x(3*j,m)-x(3*i,m))**2)        
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
       
        return
        end subroutine
