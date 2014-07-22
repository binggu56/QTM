      SUBROUTINE ADAMS6(Y,N,JC,X,HX,DERIVS,EPS,JARDT)
C
C  INTEGRATION PAR LA METHODE D'ADAMS-MOULTON : PREDICTEUR D'ORDRE 6 ET
C  CORRECTEUR D'ORDRE 7 (written by Leforestier, got from Marston)
C **********************************************************************
C *                                                                    *
C *                          SPECIFICATIONS                            *
C * Y : TABLEAU A N DIMENSIONS                                         *
C *      PREMIER APPEL : CONTIENT LES CONDITIONS INITALES              *
C *      APRES : CONTIENT LES VALEURS DES FONCTIONS INTEGREES          *
C * JC : INITIALISE A 1 PAR LE PROGAMME APPELANT LORS DU PREMIER APPEL *
C * N : DIMENSION DU SYSTEME D'EQUATIONS DIFFERENTIELLES               *
C * X : PREMIER APPEL VALEUR CORRESPONDANT AUX CONDITIONS INITIALES    *
C *      VALEUR EN RETOUR CORRESPONDANT AUX Y RETOURNES RETOURNES      *
C * HX : PAS UTILISE LORS DE L'INTEGRATION ADAMS-MOULTON               *
C * DERIVS :FONCTION DECLAREE EXTERNAL PAR LE PROGRAMME APPELANT       *
C *      ET QUI SPECIFIE LE SYSTEME DIFFERENTIEL                       *
C * EPS : ECART MAXIMAL RELATIF ENTRE DEUX VALEURS ITEREES ADAMS       *
C *      -MOULTON                                                      *
C * LM : NOMBRE D'ITERATIONS MAXIMUM PAR ADAMS-MOULTON                 *
C * $ : RETOUR NON STANDARD EN CAS DE NON-CONVERGENCE                  *
C * LK : NOMBRE D'ITERATIONS REELLEMENT EFFECTUEES                     *
C * M : DIMENSION MAXIMUM DU SYSTEME DIFFERENTIEL                      *
C * MM : RAPPORT ENTRE LE PAS SPECIFIE EN ARGUMENT ET LE PAS UTILISE   *
C * POUR LE CALCUL DE MM*6 POINTS ESPACES DE HX/MM                     *
C *                                                                    *
C **********************************************************************
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION XY(30,6),XZ(30,7),V(30),W(30),U(30,2),XX(30,6)
      DIMENSION C(6),D(7),B(3),BB(3),Y(1)
      EXTERNAL DERIVS
C
      DATA LM/8/,MM/5/
C
      DATA B,BB/.5E0,.5E0,1.E0,.5E0,0.E0,.5E0/
      DATA C/.3298611111111111    E0,.3486111111111111    E0,.375E0,.41
     L6666666666666    E0,.5E0,1.E0/,D/-.1426917989417989    E-1 ,-.187
     $5E-1,-.2638888888888888    E-1,-.4166666666666666    E-1,-.8
     $3333333333333333    E-1 ,-.5E0,1.E0/
   36 IC=JC
      JC=JC+1
      IF(IC.NE.1) GO TO 1
C INITIALISATION
      H=HX/MM
      INIT=0
      IF(MM.GT.1) INIT=1
      JM=5*MM+1
      CALL DERIVS(X,Y,XX)
      DO 30 K=1,N
   30 XY(K,1)=Y(K)
    1 IF(JC-7)100,3,12
C INTEGRATION PAR RUNGE-KUTTA D'ORDRE 4 POUR LES 5 PREMIERS POINTS POUR
C INITIALISER ADAMS-MOULTON
C CALCUL DE KV -> Y
  100 DO 5 K=1,N
    5 V(K)=XY(K,IC)
      J=0
    2 J=J+1
      CALL DERIVS(X,V,W)
      IF(J.GT.3) GO TO 6
      DO 7 K=1,N
      XZ(K,J)=H *W(K)
    7 V(K)=XY(K,IC)+B(J)*XZ(K,J)
      X=X+BB(J)*H
      GO TO 2
    6 DO 8 K=1,N
      XZ(K,4)=H *W(K)
    8 Y(K)=XY(K,IC)+(XZ(K,1)+XZ(K,2)*2.+XZ(K,3)*2.+XZ(K,4))/6.
C CALCUL DE F(N) -> XY(IC)
      XP=X-H
      CALL DERIVS(XP,XY(1,IC),V)
      DO 31 K=1,N
      XY(K,JC)=Y(K)
   31 XY(K,IC)=V(K)
      IF(INIT.EQ.1) GO TO 40
      RETURN
C CONSTRUCTION DES NABLA(I)F(N) INITIAUX
    3 DO 9 K=1,N
    9 V(K)=XY(K,IC)
      CALL DERIVS(X ,V,XY(1,IC))
      DO 10 K=1,N
      XY(K,1)=XY(K,6)-5*XY(K,5)+10*XY(K,4)-10*XY(K,3)+5*XY(K,2)-XY(K,1)
      XY(K,2)=XY(K,6)-4*XY(K,5)+6*XY(K,4)-4*XY(K,3)+XY(K,2)
      XY(K,3)=XY(K,6)-3*XY(K,5)+3*XY(K,4)-XY(K,3)
      XY(K,4)=XY(K,6)-2*XY(K,5)+XY(K,4)
   10 XY(K,5)=XY(K,6)-XY(K,5)
      GO TO 11
   41 INIT=0
      H=HX
C CALCUL DES SECONDS NABLA INITIAUX
      DO 42 K=1,N
      XY(K,1)=XX(K,6)-5*XX(K,5)+10*XX(K,4)-10*XX(K,3)+5*XX(K,2)-XX(K,1)
      XY(K,6)=XX(K,6)
      XY(K,2)=XX(K,6)-4*XX(K,5)+6*XX(K,4)-4*XX(K,3)+XX(K,2)
      XY(K,3)=XX(K,6)-3*XX(K,5)+3*XX(K,4)-XX(K,3)
      XY(K,4)=XX(K,6)-2*XX(K,5)+XX(K,4)
   42 XY(K,5)=XX(K,6)-XX(K,5)
      GO TO 11
C CALCUL DES NABLA(I)F(N) COURANTS
   12 IF(IC-JM)101,41,101
  101 DO 13 J=1,5
      DO 13 K=1,N
   13 XY(K,J)=XY(K,J+1)
      CALL DERIVS(X ,V,XY(1,6))
      DO 14 IX=1,5
      I=6-IX
      DO 14 K=1,N
   14 XY(K,I)=XY(K,I+1)-XY(K,I)
C CALCUL DE Y0(N+1) -> U(1)
   11 LL=0
      X=X+H
      DO 15 K=1,N
      A=0.D0
      DO 16 I=1,6
   16 A=A+C(I)*XY(K,I)
   15 U(K,1)=H*A+V(K)
C TEST DE CONVERGENCE
   23 IF(LL.EQ.0) GO TO 17
      DO 18 K=1,N
      IF(ABS(U(K,1)).LT.1E-10) GOTO 18
      IF((ABS((U(K,1)-U(K,2))/U(K,1))).GT.EPS) GOTO 17
   18 CONTINUE
      DO 19 K=1,N
      V(K)=U(K,IN)
   19 Y(K)=V(K)
      LK=LL
      IF(INIT.EQ.1) GO TO 40
      RETURN
C ITERATIONS
C CALCUL DES NABLA(I)FLL(N+1) -> XZ(7-I)
   17 IF(LL.GE.LM) JARDT=1
      IF(LL.GE.LM) RETURN
      IN=MOD(LL,2)+1
      CALL DERIVS(X,U(1,IN),XZ(1,7))
      DO 20 IX=1,6
      I=7-IX
      DO 20 K=1,N
   20 XZ(K,I)=XZ(K,I+1)-XY(K,I)
C CALCUL DE YLL(N+1) -> U
      LL=LL+1
      IN=MOD(LL,2)+1
      DO 21 K=1,N
      A=0.D0
      DO 22 I=1,7
   22 A=A+D(I)*XZ(K,I)
   21 U(K,IN)=H*A+V(K)
      GO TO 23
   40 JK=JC-1
      IF((MOD(JK,MM)).NE.0) GO TO 36
      JK=JK/MM+1
      CALL DERIVS(X,Y,XX(1,JK))
      RETURN
C  99 RETURN A1
   99 JARDT=1
      RETURN
      END
