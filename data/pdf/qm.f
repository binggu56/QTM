      program main

! --- Compute pair distribution function 
! --- input : postion array x(ndim, ntraj), or read from data file 
! --- compute over bins
! --- g(r) = sum over trajectries r(ij)

      implicit real*8(a-h,o-z)

      include 'sizes.h'

      include 'qsats.h'

      integer*4, parameter :: ntraj = 19200, ndim = 540, nbins = 2000

      real*8 :: x(ndim, ntraj), gr(nbins),w(ntraj)

      real*8, dimension(NATOMS,NATOMS) :: r12

      open(10, file = 'gr.dat')

      if(ndim/3 /= NATOMS ) then
        write(*,*) 'ERROR: illegal number of DoFs'
      endif

      open(11, file = '/home/bing/data/1.4.2/dat1/temp.dat')
!      open(11, file= '/home/bing/data/darter/dat2/temp.dat')
      open(11, file = '/home/bing/data/darter/1.4.4/dat4/temp.dat')
      do i=1,ntraj
        do j=1,ndim
          read(11, *) x(j,i),a,b
        enddo
      enddo
      close(11)
! --- read crystal lattice points.

      w = 1d0/dble(Ntraj)

      den = 5.231d-3

      ltfile = 'ltfile'

      write (6, 6200) ltfile
6200  format ('READING crystal lattice from ', a16/)

      open (8, file = ltfile, status='old')

      read (8, *) nlpts

      if (nlpts.ne.NATOMS) then
         write (6, *) 'ERROR: number of atoms in lattice file = ', 
     +                 nlpts
         write (6, *) 'number of atoms in source code = ', NATOMS
         stop
      end if

! --- read the edge lengths of the supercell.
      
      read (8, *) xlen, ylen, zlen

! --- compute a distance scaling factor.

      den0=dble(NATOMS)/(xlen*ylen*zlen)

! --- scale is a distance scaling factor, computed from the atomic
!     number density specified by the user.

      scale=exp(dlog(den/den0)/3.0d0)

      write (6, 6300) scale
6300  format ('supercell scaling factor computed from density = ',
     +         f12.8/)

      xlen=xlen/scale
      ylen=ylen/scale
      zlen=zlen/scale

      write (6, 6310) xlen, ylen, zlen
6310  format ('supercell edge lengths [bohr]         = ', 3f10.5/)

      dxmax=half*xlen
      dymax=half*ylen
      dzmax=half*zlen

      do i=1, NATOMS

         read (8, *) xtal(i, 1), xtal(i, 2), xtal(i, 3)

         xtal(i, 1)=xtal(i, 1)/scale
         xtal(i, 2)=xtal(i, 2)/scale
         xtal(i, 3)=xtal(i, 3)/scale

      end do

      close (8)

      bin = 1d-2 
      gr = 0d0

!      allocate(r12(natoms,natoms))

      traj: do k=1,ntraj

        do i=2,NATOMS 
          do j=1,i-1 
!        do n=1, nvpair2

!          i = ivpair2(1,n)
!          j = ivpair2(2,n)

!      r12(j,i) =  (x(3*i-2,k)-x(3*j-2,k))**2+(x(3*i-1,k)-x(3*j-1,k))**2+(x(3*i,k)-x(3*j,k))**2
          dx=(x(3*j-2,k)+xtal(j,1)-(x(3*i-2,k)+xtal(i,1)))
          dy=(x(3*j-1,k)+xtal(j,2)-(x(3*i-1,k)+xtal(i,2)))
          dz=(x(3*j,k)  +xtal(j,3)-(x(3*i,k)+xtal(i,3)))

          r = dsqrt(dx*dx+dy*dy+dz*dz)

          ibin = INT(r/bin) 

          if(ibin < NBINS) then
            gr(ibin) = gr(ibin)+w(k)/r**2
          endif

          enddo

        enddo 

      end do traj

      do i=1,nbins
        write(10,1000) i*bin, gr(i)
      enddo

      close(10)

1000  format(20(e14.7,1x))
      end program main
      
