      module cdat

        implicit real*8(a-h,o-z)
      
        real*8, public, parameter :: pi = 4.d0*atan(1.d0), half = 0.5d0,  & 
                  one = 1.0d0 

        complex*16, public, parameter :: im=(0d0,1d0)

        integer*4 :: nb 
!        integer*4 :: ntraj,ndim,nb,kmax

!        real*8    :: dt,dt2,am0,x00,p00,cf0 

      save

      contains

        double precision function trace(n,A)

        real*8, dimension(n,n), intent(in) :: A

        trace = 0d0 

        do i=1,n
          trace = trace + A(i,i)
        enddo 

        end function trace  

! ----- average over trajectories y = sum(w(i)*x(i)) 
        double precision function aver_traj(ntraj,w,x)

        integer*4, intent(in) :: Ntraj
        
        real*8 :: x(Ntraj),w(Ntraj)

        aver_traj = 0d0

        do i=1,Ntraj
          aver_traj = aver_traj + x(i)*w(i) 
        enddo

        end function aver_traj
        
! ----- random number seeds
        subroutine seed(idum,Ndim)
        
        implicit real*8(a-h,o-z)
        
        integer*4, intent(in) :: Ndim
        
        integer*4, intent(out) :: idum(Ndim)

        do i=1,Ndim
          idum(i) = 5 + i
        enddo

        return
        end subroutine

      end module cdat
