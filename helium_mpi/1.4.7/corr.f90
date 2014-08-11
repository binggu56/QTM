! --- subroutine to compute correlation function C(AB)
      subroutine corr(am,dt,nb,ndim,ntraj,w,x,g,cor)
     
      implicit real*8 (a-h,o-z)

      real*8, intent(in) :: am(ndim),w(ntraj),x(ndim,ntraj)

      integer*4 ,intent(in) :: nb,ntraj,ndim 

      real*8 :: f(nb)

      complex*16,intent(in) :: g(nb)

      complex*16,intent(out) :: cor 
      
      complex*16 :: z0

      cor = (0d0,0d0) 

      do i=1,ntraj 
        do j=1,ndim 
          f(j) = x(j,i)
        enddo 
        f(ndim+1) = 1d0 
        
        z0 = dot_product(f,g)

        do j=1,ndim
            cor = cor+ conjg(z0)*x(j,i)*w(i)
        enddo 

      enddo
      
      return 
      end subroutine 

! --- get time derivative of coefficient vector c  
      subroutine prop_c(s,mat1,am,dt,ndim,ntraj,c)

      use cdat, only : im, nb  

      implicit real*8 (a-h, o-z) 

      integer*4, intent(in) :: ndim,ntraj
      
      real*8,    intent(in) :: dt,am(ndim),mat1(nb,nb)

      real*8 :: f(nb),df(ndim,nb),s(nb,nb),sinv(nb,nb),mat2(nb,nb), mats(nb,nb)

      complex*16 :: cor,mat(nb,nb),c(nb),dc(nb)

      am0 = am(1)

! --- get df(ndim,nb)
      df = 0d0 
      do k=1,ndim 
        df(k,k) = 1d0 
      enddo 

! --- matrix mat2 = <df|df> 

      mat2 = matmul(transpose(df),df)

! --- matrix matrix mat = <f|f> 
!      s = 0d0 
!      do i=1,ntraj 
!
!        do j=1,ndim 
!          f(j) = x(j,i)
!        enddo 
!        f(ndim+1) = 1d0 
!        
!        do k=1,nb
!          do j=1,k  
!            s(j,k) = s(j,k)+f(j)*f(k)*w(i)
!          enddo 
!        enddo 
!      enddo 
!
!      do k=1,nb 
!        do j=1,k
!          s(k,j) = s(j,k)
!        enddo 
!      enddo       


! --- matrix mat1 = <p*f|df>, simplified version for linear basis 
!     needs modification if using higher order basis 
!     computation of mat1 will be distributed  
!      mat1 = (0d0,0d0)
!
!      traj: do i=1,ntraj
!
!        do j=1,ndim 
!          f(j) = x(j,i) 
!        enddo 
!        f(ndim+1) = 1d0 
!
!        do j=1,ndim 
!          do k=1,nb 
!              mat1(k,j) = mat1(k,j)+p(j,i)*f(k)*df(j,j)*w(i)
!          enddo 
!        enddo 
!        
!      enddo traj 

      mat = 2d0*mat1 + im*mat2
      mat = -mat/2d0/am0

      dc = matmul(mat,c)

      mats = s 
      
      call inverse(mats,sinv,nb)
      
      dc = matmul(sinv,dc)
 
      do j=1,nb 
        c(j) = c(j)+dc(j)*dt 
      enddo 
           
      return 
      end subroutine 

