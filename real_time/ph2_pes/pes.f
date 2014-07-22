        program main 
        implicit real*8(a-h,o-z)
        real*8,dimension(:,:),allocatable :: x
        open(100,file='pes.dat')
	open(105,file='cris.xyz')

        Np = Ndim/3
	Ntraj = 1
	Ndim = 12
	am = 2d0*1836.18d0
        c6 = 12.14d0
        c8 = 215.2d0
        c10 = 4813.9d0
        c9 = 143.1d0
        al = 1.713
        be = 1.5671
        ga = 0.00993d0
        rm = 3.44d0/0.52918d0
        v = 0d0
	rho = 0.02509d0
        d = (4d0/rho)**(1d0/3d0)/0.52918d0*dsqrt(2d0)/2d0
        m = 1
	dx = 0.06d0
        
	allocate(x(Ndim,Ntraj))
        x = 0d0
	x(1,1) = d/2d0-3d0
	x(2,1) = 0d0
	x(4,1) = 0d0
	x(5,1) = d/2d0
	x(7,1) = 0d0
	x(9,1) = -d/2d0
	x(10,1) = d/2d0
	x(11,1) = d/2d0
	x(12,1) = -d/2d0
	
        write(*,*) 'd/2d0 =', d/2d0 
        do n=1,100
	  x(1,1) = x(1,1)+dx
	  call pot(d,m,Ntraj,Ndim,x,v)
          write(100,1000) x(1,1),v
        enddo

	 ! call pot(d,m,Ntraj,Ndim,x,v0)
         ! x(1,1) = x(1,1)+dx
	 ! call pot(d,m,Ntraj,Ndim,x,v1)
	 ! f1 = (v1-v0)/dx
	 ! x(1,1) = x(1,1)-2d0*dx
	 ! call pot(d,m,ntraj,Ndim,x,v2)
	 ! f2 = (v0-v2)/dx
	 ! al = (f1-f2)/dx
	 ! write(*,*) 'alpha = ',dsqrt(al*am)/2d0

        
1000    format(20(e14.7,1x))
        end program
