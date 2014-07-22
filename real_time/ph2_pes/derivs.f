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

