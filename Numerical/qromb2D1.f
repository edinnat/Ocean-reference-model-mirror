        SUBROUTINE qromb2D1(FUNC1,A1,B1,SS1, EPS1)
        
        implicit none
        
        common /Flag1/ Flag1

        integer j1, k1, km1, jmax1, jmaxp1, IT1
        double precision eps1, Ncount1
        external func1
        parameter (jmax1=15,jmaxp1=jmax1+1,k1=5,km1=k1-1)
        double precision h1(jmaxp1)
        double precision s1(jmaxp1,1:2)
        double precision sTemp1(1:2)
        double precision a1
        double precision b1
        double precision ss1(1:2)
        double precision dss1(1:2)
c        double precision func1
        integer Flag1(1:2)
        integer i1(1:2)
        integer indice1

        Ncount1 = 2
        do indice1 = 1, 2
                H1(1)          = 1.0D0
                i1(indice1)    = 0
                Flag1(indice1) = 1
        enddo
        DO J1 = 1, JMAX1
        CALL trapzd2D1(FUNC1, A1, B1, sTemp1, J1, IT1)
        do indice1 = 1, 2
        if (flag1(indice1).eq.1) then
                s1(J1, indice1) = sTemp1(indice1)
                IF (J1.GE.K1) THEN
                CALL polint(h1(J1-KM1), s1(J1-KM1, indice1), K1
     &                      ,0.D0, SS1(indice1), DSS1(indice1))
                IF (ABS(DSS1(indice1)).LT.EPS1*ABS(SS1(indice1)))then
                                i1(indice1) = i1(indice1) + 1
                        else
                                i1(indice1) = 0
                endif
                if (abs(SS1(indice1)).lt.1.0D-08) i1(indice1) = 2
                        IF (i1(indice1).eq.2) then
                                flag1(indice1) = 0
                                Ncount1 = Ncount1 - 1
                                if (Ncount1.eq.0) return 
                        endif
                endif
        ENDIF
        S1(J1+1, indice1) = S1(J1, indice1)
        H1(J1+1) = 0.25D0*H1(J1)
        enddo
        enddo
        do indice1 = 1,2
                if (Flag1(indice1).eq.1) then
                print *, '_____________  QROMB2D1 ________________'
                print *, 'l''intégrale #',indice1,' n''a pas convergé'
                print *, 'res = ',SS1(indice1),' Dres = ',DSS1(indice1)
                endif
        enddo
        
        END
