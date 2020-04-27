      SUBROUTINE LAGUER(A,M,X,EPS,POLISH)
      double COMPLEX A(*),X,DX,X1,B,D,F,G,H,SQ,GP,GM,G2,ZERO,XX,WW
      LOGICAL POLISH
      PARAMETER (ZERO=(0.,0.),TINY=1.E-15,MAXIT=100)
      IF (POLISH) THEN
        DXOLD=ABS(X)
        NPOL=0
      ENDIF
      DO 12 ITER=1,MAXIT
        B=A(M+1)
        D=ZERO
        F=ZERO
        DO 11 J=M,1,-1
          F=X*F+D
          D=X*D+B
          B=X*B+A(J)
11      CONTINUE
        IF(ABS(B).LE.TINY) THEN
          DX=ZERO
        ELSE IF(ABS(D).LE.TINY.AND.ABS(F).LE.TINY)THEN
          DX=CMPLX(ABS(B/A(M+1))**(1./M),0.)
        ELSE
          G=D/B
          G2=G*G
          H=G2-2.*F/B
	  XX=(M-1)*(M*H-G2)
	  YY=ABS(realpart(XX))
	  ZZ=ABS(imagpart(XX))
	  IF(YY.LT.TINY.AND.ZZ.LT.TINY) THEN
	    SQ=ZERO
	  ELSE IF (YY.GE.ZZ) THEN
	    WW=(1.0/YY)*XX
	    SQ=SQRT(YY)*sqrt(WW)
	  ELSE
	    WW=(1.0/ZZ)*XX
	    SQ=SQRT(ZZ)*sqrt(WW)
	  ENDIF
          GP=G+SQ
          GM=G-SQ
          IF(ABS(GP).LT.ABS(GM)) GP=GM
          DX=M/GP
        ENDIF
        X1=X-DX
        IF(X.EQ.X1)RETURN
        X=X1
        IF (POLISH) THEN
          NPOL=NPOL+1
          CDX=ABS(DX)
          IF(NPOL.GT.9.AND.CDX.GE.DXOLD)RETURN
          DXOLD=CDX
        ELSE
          IF(ABS(DX).LE.EPS*ABS(X))RETURN
        ENDIF
12    CONTINUE
      PAUSE 'too many iterations'
      RETURN
      END
