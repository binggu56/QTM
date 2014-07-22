      SUBROUTINE QSIMP(FUNC,A,B,S)
C--------------------------------------------
C    DOPPELT GENAUE VERSION DER SIMPSONREGEL
C--------------------------------------------
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      EXTERNAL FUNC
C
      PARAMETER (EPS=1.D-6, JMAX=200)
      OST=-1.D30
      OS= -1.D30
      DO 11 J=1,JMAX
        CALL TRAPZD(FUNC,A,B,ST,J)
        S=(4.*ST-OST)/3.
        IF (ABS(S-OS).LT.EPS*ABS(OS)) RETURN
        OS=S
        OST=ST
11    CONTINUE
      PAUSE 'Too many steps.'
      END
