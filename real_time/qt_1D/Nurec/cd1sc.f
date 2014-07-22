         PROGRAM CD1SC
C------------------------------------------------------------
C       CALCULATES SEMICLASSICAL CORRELATION FUNCTION
C       FOR SCATTERING WITH ERIC HELLER'S CELLULAR DYNAMICS
C       FRANK GROSSMANN, 7. FEBRUAR 1994 (last modified:)
C------------------------------------------------------------
	 IMPLICIT NONE
	 INTEGER I,K,N,M,NMAX,MMAXL,MPS,IRK,NPER,NSCH,MSW,MM,MT
         INTEGER NU,NUM(-100:100),NUT(-100:100)
	 DOUBLE PRECISION YS(7),YM(7,-100:100),YT(7,-100:100)
         DOUBLE PRECISION YIM(2,-100:100),YIT(2,-100:100)
	 DOUBLE PRECISION AL,PI,BE,DBE
	 DOUBLE PRECISION QA,QB,P0,PT,DP,QM,PM,PL
	 DOUBLE PRECISION D,Z,TAU,T
	 COMPLEX*16 CORR,CEX,IM,ALS
C
         COMMON /PARAMS/ D,Z
	 COMMON /RECH1/ NPER,NSCH,IRK
         COMMON /RECH3/ AL,BE,QA,QB,P0
         COMMON /RECH4/ YIM,YIT
C
         WRITE(*,*)'D,Z'
         WRITE(*,*)'ECKART BARRIER: x,y'
         WRITE(*,*)'FREE PARTICLE: 0,0'
         READ(*,*)D,Z
         WRITE(*,*)'ANFANGS-ENDZUSTAND: QA,QB,P0,AL'
         READ(*,*)QA,QB,P0,AL
         WRITE(*,*)'GRIDPARAMETER: BE,DBE,MMAXL'
         READ(*,*)BE,DBE,MMAXL
         WRITE(*,*)'RUNGE-KUTTA-GEN.: IRK'
         READ(*,*)IRK
	 WRITE(*,*)'TAU,NPER,NSCH'
	 READ(*,*)TAU,NPER,NSCH
C
	 OPEN(UNIT=30,FILE='ce.redat')
	 OPEN(UNIT=31,FILE='ce.imdat')
	 OPEN(UNIT=32,FILE='tr.dat')
C
	 WRITE(30,*)'PROGRAM CESC'
         WRITE(30,*)D,Z
         WRITE(30,*)QA,QB,P0,AL
         WRITE(30,*)BE,DBE,MMAXL
         WRITE(30,*)IRK
	 WRITE(30,*)TAU,NPER,NSCH
C
         PI=4*ATAN(1.D0)
         IM=(0.D0,1.D0)
         T=2**NPER*TAU
         IF (BE.NE.0) DP=(DBE/BE)**.5
C
*****************************************************************
******************BACKWARD PROPAGATION***************************
*****************************************************************
C
          TAU=-TAU
          DO K=1,2**NPER
	   DO M=-MMAXL,MMAXL
            IF (K.EQ.1) THEN
             PM=P0+M*DP
             QM=QB+T*PM
             YIM(1,M)=QM
             YIM(2,M)=PM
             YS(1)=QM
             YS(2)=PM
             YS(3)=0
             YS(4)=1
             YS(5)=0
             YS(6)=0
             YS(7)=1
             NU=0
            ELSE
             DO I=1,7
              YS(I)=YM(I,M)
             ENDDO
             NU=NUM(M)
            ENDIF
C
            CALL KLASSIK(YS,NU,TAU)
C
            IF (K.EQ.2**NPER) WRITE(32,*)YS(1),YS(2)
            DO I=1,7
             YM(I,M)=YS(I)
            ENDDO
            NUM(M)=NU
           ENDDO
          ENDDO
C
*****************FORWARD PROPAGATION****************************
C
          TAU=-TAU
          DO K=1,2*2**NPER
	   DO M=-MMAXL,MMAXL
            IF (K.EQ.1) THEN
             PM=P0+M*DP
             QM=QA-T*PM
             YIT(1,M)=QM
             YIT(2,M)=PM
             YS(1)=QM
             YS(2)=PM
             YS(3)=0
             YS(4)=1
             YS(5)=0
             YS(6)=0
             YS(7)=1
             NU=0
            ELSE
             DO I=1,7
              YS(I)=YT(I,M)
             ENDDO
             NU=NUT(M)
            ENDIF
C
            CALL KLASSIK(YS,NU,TAU)
C
            IF (K/10.EQ.K/10.) WRITE(32,*)YS(1),YS(2)
            DO I=1,7
             YT(I,M)=YS(I)
            ENDDO
            NUT(M)=NU
           ENDDO
C
           CALL UEB(YM,YT,NUM,NUT,MMAXL,T,CORR)
           WRITE(30,*)(K-2**NPER)*TAU,REAL(CORR)*DBE/PI
           WRITE(31,*)(K-2**NPER)*TAU,IMAG(CORR)*DBE/PI
          ENDDO
****************************************************************
****************************************************************
C
	 CLOSE(30)
	 CLOSE(31)
         CLOSE(32)
C
         STOP
         END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
	 SUBROUTINE KLASSIK(YS,NU,TAU)
C---------------------------------------------------------------
C       CLASSICAL INFORMATION FOR ONE TRAJECTORY
C---------------------------------------------------------------
	 IMPLICIT NONE
	 INTEGER I,NPER,NSCH,IRK,KMAX,KOUNT,NVAR,NOK,NBAD
         INTEGER NU
	 DOUBLE PRECISION YS(7),Y(7),DYDX(7)
	 DOUBLE PRECISION XP(200),YP(200,200)
         DOUBLE PRECISION X,X1,X2,EPS,PI,DXSAV
	 DOUBLE PRECISION H1,HMIN
	 DOUBLE PRECISION ALTS11,ALTS21,ALTS12,ALTS22
	 DOUBLE PRECISION D,Z,TAU
	 EXTERNAL DERIVS,RKQC
C
         COMMON /PARAMS/ D,Z
	 COMMON /RECH1/ NPER,NSCH,IRK
	 COMMON /PATH/ KMAX,KOUNT,DXSAV,XP,YP
C
         PI=4*ATAN(1D0)
         EPS=10D0**(-IRK)
         NVAR=7
C
         DO I=1,2**NSCH
          ALTS11=YS(4)
	  ALTS12=YS(5)
	  ALTS21=YS(6)
          ALTS22=YS(7)
          X1=(I-1)*TAU/(2.**NSCH)
          X2=I*TAU/(2.**NSCH)
          H1=(X2-X1)/10.
          HMIN=(X2-X1)/1000000.
          CALL ODEINT(YS,NVAR,X1,X2,EPS,H1,HMIN,NOK,NBAD,DERIVS,RKQC)
C
	  IF (SIGN(1D0,YS(6)).NE.SIGN(1D0,ALTS21)) THEN
           IF (SIGN(1D0,YS(6)/YS(4)).EQ.1D0) THEN
	    NU=NU+1
           ELSE IF (SIGN(1D0,YS(6)/YS(4)).EQ.-1D0) THEN
            NU=NU-1
           ENDIF
	  ENDIF
         ENDDO
C
         RETURN
	 END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         SUBROUTINE DERIVS(X,Y,DYDX)
C----------------------------------------------------------
C       RIGHT HAND SIDE OF THE DIFFERENTIAL EQUATION
C----------------------------------------------------------
	 IMPLICIT NONE
         DOUBLE PRECISION Y(7),DYDX(7)
         DOUBLE PRECISION X,V,DVX,D2VXX
	 DOUBLE PRECISION D,Z
C
         COMMON /PARAMS/ D,Z
C
         V=D/((COSH(Z*Y(1)))**2)
         DVX=-2*Z*TANH(Z*Y(1))*V
         D2VXX=-2*Z**2*V*(1-3*TANH(Z*Y(1))**2)
C
         DYDX(1)=Y(2)
         DYDX(2)=-DVX
	 DYDX(3)=Y(2)**2/(2D0)-V
         DYDX(4)=-D2VXX*Y(6)
         DYDX(5)=-D2VXX*Y(7)
	 DYDX(6)=Y(4)
	 DYDX(7)=Y(5)
C
         RETURN
         END
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         SUBROUTINE UEB(YM,YT,NUM,NUT,MMAXL,T,CORR)
C--------------------------------------------------------------
C       CALCULATES THE CORRELATION FUNCTION 
C--------------------------------------------------------------
	 IMPLICIT NONE
         INTEGER MMAXL,M,N
         INTEGER NUM(-100:100),NUT(-100:100)
	 DOUBLE PRECISION YM(7,-100:100),YT(7,-100:100)
         DOUBLE PRECISION YIM(2,-100:100),YIT(2,-100:100)
	 DOUBLE PRECISION AL,BE,QA,QB,P0,T
	 DOUBLE PRECISION DQT,DQA,DQB,PI
         COMPLEX*16 CORR,DEN
	 COMPLEX*16 IM,EX,PRE,DET,A0,A11,A12,A13,A22,A33,B1,B2,B3
C
         COMMON /RECH3/ AL,BE,QA,QB,P0
         COMMON /RECH4/ YIM,YIT
C
         PI=4*ATAN(1.D0)
	 IM=(0.D0,1.D0)
         DEN=1-2*IM*AL*T
         CORR=(0.D0,0.D0)
C
         DO N=-MMAXL,MMAXL
          DQA=YIT(1,N)-QA
          DO M=-MMAXL,MMAXL
           DQB=YIM(1,M)-QB
           DQT=YM(1,M)-YT(1,N)
C
	   A0=IM*(YT(3,N)-YM(3,M))-IM*PI*(NUT(N)-NUM(M))/2D0
     +       +IM*YT(2,N)*DQT
     +       +(IM/2.*YT(4,N)/YT(6,N)-BE/YT(6,N)**2)*DQT**2
     +       -1/DEN*(AL*(DQA**2+DQB**2)-IM*P0*(DQA-DQB)-IM*T*P0**2)
C
           A11=-IM/2.*YT(4,N)/YT(6,N)+BE/YT(6,N)**2
     +         +IM/2.*YM(4,M)/YM(6,M)+BE/YM(6,M)**2
           A12=-IM/(2*YM(6,M))+2*BE*YM(7,M)/YM(6,M)**2
           A13=IM/(2*YT(6,N))+2*BE*YT(7,N)/YT(6,N)**2
           A22=AL/DEN+IM/2.*YM(7,M)/YM(6,M)+BE*(YM(7,M)/YM(6,M))**2
           A33=AL/DEN-IM/2.*YT(7,N)/YT(6,N)+BE*(YT(7,N)/YT(6,N))**2
C
           B1=IM*(YT(2,N)-YM(2,M))
     +       +(IM*YT(4,N)/YT(6,N)-2*BE/YT(6,N)**2)*DQT
           B2=IM*YIM(2,M)-1/DEN*(2*AL*DQB+IM*P0)
           B3=-IM*YIT(2,N)-1/DEN*(2*AL*DQA-IM*P0)
     +       -(IM/YT(6,N)+2*BE*YT(7,N)/YT(6,N)**2)*DQT
C
           DET=A11*A22*A33-A12**2*A33-A13**2*A22
           EX=A0+(B1**2*A22*A33-2*B1*B2*A12*A33-2*B1*B3*A13*A22
     +         +B2**2*(A11*A33-A13**2)+2*B2*B3*A13*A12
     +         +B3**2*(A11*A22-A12**2))/(4*DET)
           PRE=1/(2*DEN)*SQRT(2*AL/(DET*ABS(YM(6,M))*ABS(YT(6,N))))
C
           CORR=CORR+PRE*EXP(EX)
          ENDDO
         ENDDO
C
         RETURN
         END
