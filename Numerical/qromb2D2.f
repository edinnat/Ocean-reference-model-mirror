        SUBROUTINE qromb2D2(FUNC2, A2, B2, SS2, EPS2)
        
        implicit none
        
        common /Flag2/ Flag2

        integer j2, k2, km2, jmax2, jmaxp2, IT2
        double precision eps2, Ncount2
        external func2
        parameter (jmax2=15,jmaxp2=jmax2+1,k2=5,km2=k2-1)
        double precision h2(jmaxp2)
        double precision s2(jmaxp2, 1:2)
        double precision sTemp2(1:2)
        double precision a2
        double precision b2
        double precision ss2(1:2)
        double precision dss2(1:2)
c        double precision func2
        integer Flag2(1:2)
        integer i2(1:2)
        integer indice2

        Ncount2 = 2
        do indice2 = 1, 2
                H2(1)          = 1.0D0
                i2(indice2)    = 0
                Flag2(indice2) = 1
        enddo
        DO J2 = 1, JMAX2
        CALL trapzd2D2(FUNC2, A2, B2, sTemp2, J2, IT2)
        do indice2 = 1, 2
        if (flag2(indice2).eq.1) then
                s2(J2, indice2) = sTemp2(indice2)
                IF (J2.GE.K2) THEN
                CALL polint(h2(J2-KM2), s2(J2-KM2, indice2), K2
     &                      , 0.D0, SS2(indice2), DSS2(indice2))
                IF (ABS(DSS2(indice2)).LT.EPS2*ABS(SS2(indice2)))then
                                i2(indice2) = i2(indice2) + 1
                        else    
                                i2(indice2) = 0
                endif
                if (abs(SS2(indice2)).lt.1.0D-08) i2(indice2) = 2
                        IF (i2(indice2).eq.2) then
                                flag2(indice2) = 0
                                Ncount2 = Ncount2 - 1
                                if (Ncount2.eq.0) return 
                        endif
                endif
        ENDIF
        S2(J2+1,indice2) = S2(J2, indice2)
        H2(J2+1) = 0.25D0*H2(J2)
        enddo
        enddo
        do indice2 = 1, 2
                if (Flag2(indice2).eq.1) then
                print *, '____________   QROMB2D2 ______________'
                print *, 'l''intégrale #',indice2,' n''a pas convergé'
                print *, 'res = ',SS2(indice2),' Dres = ',DSS2(indice2)
                endif
        enddo
        
        END
