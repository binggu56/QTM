        subroutine derivs(Nt,dt,t,x,Ntraj,fx,pe)
        implicit real*8(a-h,o-z)
        integer*4,intent(IN) :: Ntraj,Nt
        real*8,intent(IN) :: x(Ntraj),dt,t
        real*8 :: De,a,k1,k2,x0
        real*8,intent(OUT) :: fx(Ntraj),pe(Ntraj)
c   morse potential for x coordinate, Eckart potential for y direction  
        
!        call swit(Nt,dt,t,ct)
        do i=1,Ntraj
!          fx(i) = -x(i)+ct*(0.15d0*x(i)**2)
!          pe(i) = x(i)**2/2d0+ct*(-0.05d0*x(i)**3)
           pe(i) = x(i)**2/2d0
           fx(i) = -x(i)
        enddo

        return
        end subroutine

!----------------------------------------------------------
! switching function
!----------------------------------------------------------
        subroutine swit(Nt,dt,t,ct)
        implicit real*8(a-h,o-z)
        integer*4,intent(in)::Nt
        real*8,intent(in)::dt
        ct = t/(Nt*dt)
        return
        end subroutine

