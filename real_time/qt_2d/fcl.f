        subroutine derivs(Nt,dt,t,x,y,Ntraj,fx,fy,pe)
        implicit real*8(a-h,o-z)
        integer*4,intent(IN) :: Ntraj,Nt
        real*8,intent(IN) :: x(Ntraj),y(Ntraj),t,dt
        real*8 :: k1,k2
        real*8,intent(OUT) :: fx(Ntraj),fy(Ntraj),pe(Ntraj)
        common/pes/ipot
c   morse potential for x coordinate, Eckart potential for y direction  
      De = 160d0
      xe=1.4d0
C
        xmin=0.3342534836
        if (ipot .eq. 1) then
        do i=1,Ntraj
        fx(i)=0.0908295*exp(-1.5d0*(x(i)-y(i)*lambda/2d0))/
     & (1d0+exp(-1.5d0*(x(i)-y(i)*lambda/2d0)))**2
     & -0.09246906195d0*sinh(0.75d0*(x(i)-y(i)*lambda/2d0))/
     & (cosh(0.75d0*(x(i)-y(i)*lambda/2d0)))**3

        fy(i)=-k2*y(i)-0.04541475d0*lambda*exp(-1.5d0*x(i)+
     & 0.75d0*y(i)*lambda)/(1d0+exp(-1.5d0*x(i)+0.75d0*y(i)*lambda))**2+
     & 0.04623453098d0*sinh(0.75d0*x(i)-0.375d0*y(i)*lambda)*lambda/
     & cosh(0.75d0*x(i)-0.375d0*y(i)*lambda)**3

       pe(i)=0.063740d0-0.060553d0/(1d0+exp(-1.5d0*(x(i)-
     & y(i)*lambda/2d0)))
     & - 0.06164604130d0/cosh(0.75d0*(x(i)-y(i)*lambda/2d0))**2+
     & k2*y(i)**2/2d0
        end do
        endif
C Harmonic 
        if(ipot .eq. 2) then
          call swit(Nt,dt,t,ct)
          do i=1,Ntraj
            fx(i) = -x(i)-0.5d0*ct*x(i)
            fy(i) = -y(i)-0.5d0*ct*y(i)
            pe(i) = x(i)**2/2d0+y(i)**2/2d0+
     &              ct*(0.5d0*x(i)*y(i))
          enddo
        endif
C Morse
      if(ipot .eq. 3) then
      do i=1,Ntraj
      fx(i)=-2d0*De*exp(-(x(i)-xe))*(1d0-exp(-(x(i)-xe)))
      fy(i)=-ky*y(i)
      pe(i)=De*(1d0-exp(-(x(i)-xe)))**2+ky*y(i)**2/2d0
      enddo
      endif
C Henon-Hailes
      IF(ipot .eq. 4) THEN
      k1=1d0
      k2=1d0
      do i=1,Ntraj
      fx(i)=-k1*x(i)-2d0*x(i)*y(i)
      fy(i)=-k2*y(i)-x(i)**2+y(i)**2
      pe(i)=k1*x(i)**2/2d0+k2*y(i)**2/2d0+
     &      (x(i)**2*y(i)-y(i)**3/3d0)
      ENDDO
      ENDIF

      RETURN
      END subroutine

!----------------------------------------------------------
! switching function
!----------------------------------------------------------
        subroutine swit(Nt,dt,t,ct)
        implicit real*8(a-h,o-z)
        integer*4,intent(in)::Nt
        real*8,intent(in)::dt,t
        ct = t/(Nt*dt)
        return
        end subroutine

