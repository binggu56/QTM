      SUBROUTINE rkdumb(vstart,nvar,x1,x2,nstep,derivs)
C------------------------------------------------------------
C    DOUBLE PRECISION VERSION, UP TO 200 VARIABLES
C------------------------------------------------------------
      INTEGER nstep,nvar,NMAX,NSTPMX
      PARAMETER (NMAX=200,NSTPMX=200)
      REAL*8 x1,x2,vstart(nvar),xx(NSTPMX),y(NMAX,NSTPMX)
      EXTERNAL derivs
      COMMON /path/ xx,y
CU    USES rk4
      INTEGER i,k
      REAL*8 h,x,dv(NMAX),v(NMAX)
      do 11 i=1,nvar
        v(i)=vstart(i)
        y(i,1)=v(i)
11    continue
      xx(1)=x1
      x=x1
      h=(x2-x1)/nstep
      do 13 k=1,nstep
        call derivs(x,v,dv)
        call rk4(v,dv,nvar,x,h,v,derivs)
        if(x+h.eq.x)pause 'stepsize not significant in rkdumb'
        x=x+h
        xx(k+1)=x
        do 12 i=1,nvar
          y(i,k+1)=v(i)
12      continue
13    continue
      return
      END
