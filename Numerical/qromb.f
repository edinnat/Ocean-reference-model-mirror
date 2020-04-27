        SUBROUTINE qromb(FUNC,A,B,SS, EPS)
        
        implicit none
        
        common /Ncount/ Ncount

        integer j, k, km, jmax, jmaxp, IT
        double precision eps, Ncount
        external func
        parameter (jmax=24,jmaxp=jmax+1,k=5,km=k-1)
        double precision h(jmaxp), s(jmaxp), a, b, ss, dss
        double precision func
        integer i

        i = 0
        H(1)=1.D0
        DO J=1,JMAX
        CALL trapzd(FUNC,A,B,S(J),J,IT)
        IF (J.GE.K) THEN
          CALL polint(H(J-KM),S(J-KM),K,0.D0,SS,DSS)
          IF (ABS(DSS).LT.EPS*ABS(SS)) then
                  i = i + 1
          else
                  i = 0
          endif
          IF (i.eq.2) RETURN
        ENDIF
        S(J+1)=S(J)
        H(J+1)=0.25D0*H(J)
        enddo
        PAUSE 'Too many steps.'
        END
