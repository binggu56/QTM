        subroutine derivs(x,Ntraj,fx,v,ddv)
        implicit real*8(a-h,o-z)
        integer*4,intent(IN) :: Ntraj
        real*8,intent(IN) :: x(Ntraj)
        real*8,intent(OUT) :: fx(Ntraj),v(Ntraj),ddv(Ntraj)
        
        do i=1,Ntraj
C          fx(i) = -x(i)
C          v(i)  = x(i)**2/2d0
C          ddv(i) = 1d0
          fx(i) = 0d0
          v(i) = 0d0
          ddv(i) = 0d0
        enddo

        return
        end subroutine

