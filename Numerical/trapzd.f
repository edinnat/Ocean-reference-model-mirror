      SUBROUTINE TRAPZD(FUNC,A,B,S,N,IT)

      external func
      double precision a, b, s, x, del, tnm, sum, func
      integer n, j
      integer IT
      
      IF (N.EQ.1) THEN
        S=0.5D0*(B-A)*(FUNC(A)+FUNC(B))
        IT=1
      ELSE
        TNM=IT
        DEL=(B-A)/TNM
        X=A+0.5D0*DEL
        SUM=0.D0
        DO  J=1,IT
          SUM=SUM+FUNC(X)
          X=X+DEL
        enddo
        S=0.5D0*(S+(B-A)*SUM/TNM)
        IT=2*IT
      ENDIF
      RETURN
      END
