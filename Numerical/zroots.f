      SUBROUTINE ZROOTS(A,M,ROOTS,POLISH)
      PARAMETER (EPS=1.E-6,MAXM=101)
      double complex A(*),ROOTS(M),AD(MAXM),X,B,C
      LOGICAL POLISH
      DO 11 J=1,M+1
        AD(J)=A(J)
11    CONTINUE
      DO 13 J=M,1,-1
        X=CMPLX(0.,0.)
        CALL LAGUER(AD,J,X,EPS,.FALSE.)
        IF(ABS(imagpart(X)).LE.2.*EPS**2*ABS(realpart(X))) 
     &          X=CMPLX(realpart(X),0.)
        ROOTS(J)=X
        B=AD(J+1)
        DO 12 JJ=J,1,-1
          C=AD(JJ)
          AD(JJ)=B
          B=X*B+C
12      CONTINUE
13    CONTINUE
      IF (POLISH) THEN
        DO 14 J=1,M
          CALL LAGUER(A,M,ROOTS(J),EPS,.TRUE.)
14      CONTINUE
      ENDIF
      DO 16 J=2,M
        X=ROOTS(J)
        DO 15 I=J-1,1,-1
          IF(realpart(ROOTS(I)).LE.realpart(X))GO TO 10
          ROOTS(I+1)=ROOTS(I)
15      CONTINUE
        I=0
10      ROOTS(I+1)=X
16    CONTINUE
      RETURN
      END
