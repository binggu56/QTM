      SUBROUTINE rkfr(vstart,nvar,x1,x2,nstep,derivs)
C------------------------------------------------------------
C    DOUBLE PRECISION VERSION, UP TO 200 VARIABLES
C    THE RESULT IS IN vstart! NO STORAGE OF INTERMEDIATE VALUES!
C------------------------------------------------------------
      INTEGER nstep,nvar,NMAX,NSTPMX
      PARAMETER (NMAX=200,NSTPMX=200)
      REAL*8 x1,x2,vstart(nvar)
      EXTERNAL derivs
CU    USES rk4
      INTEGER i,k
      REAL*8 h,x,dv(NMAX),v(NMAX)
C
      do 11 i=1,nvar
        v(i)=vstart(i)
11    continue
      x=x1
      h=(x2-x1)/nstep
      do 13 k=1,nstep
        call derivs(x,v,dv)
        call rk4(v,dv,nvar,x,h,v,derivs)
        if(x+h.eq.x)pause 'stepsize not significant in rkdumb'
        x=x+h
13    continue
      do 14, i=1,nvar
       vstart(i)=v(i)
14    continue
      return
      END
