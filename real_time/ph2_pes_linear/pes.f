        program main 
        implicit real*8(a-h,o-z)
        parameter(Ntraj=1,Ndim=21)
        real*8 :: x(Ndim,Ntraj)
        real*8  :: v
        open(100,file='pes.dat')

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
        m = 1
        
        x = 0d0
        x(1,1) = -2d0
        x(4,1) = rm
        x(7,1) = -rm
        x(10,1) = 2d0*rm
        x(13,1) = -2d0*rm
        x(16,1) = 3d0*rm
        x(19,1) = -3d0*rm
        
        do n=1,100
          x(1,1) = x(1,1)+0.04d0
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
        write(100,1000) x(1,1),v
        enddo
        
1000    format(20(e14.7,1x))
        end program
