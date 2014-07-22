        subroutine derivs(lambda,x,y,Ntraj,fx,fy,pe)
        implicit none
        integer*4,intent(IN) :: Ntraj
        real*8,intent(IN) :: x(Ntraj),y(Ntraj),lambda
        integer*4 :: ipot,i
        real*8 :: De,a,x0,y0,b,c,d,e,v,f,w,ky,xmin
        real*8,intent(OUT) :: fx(Ntraj),fy(Ntraj),pe(Ntraj)
      real*8 xe,k1,k2
      common/pes/ipot
c        common/group/qx0,qy0
c   morse potential for x coordinate, Eckart potential for y direction  
      k2=0.3696726511
      De = 160d0
      ky=1d0
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
        do i=1,Ntraj
        fx(i)=-x(i)
        fy(i)=-y(i)
        pe(i)=x(i)**2/2d0+y(i)**2/2d0
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
