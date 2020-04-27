c______________________________________________________________________
c       Calcul de la constant c, paramètre du Spectre de
c       Durden & Vesecky. U12 est le vent à 12.5 mètres.
c______________________________________________________________________

        subroutine c_ (c, U12)

        implicit none
        external funcD2, funcD1
        double precision c, D1, D2, R, D, D11, D22, EPS, U12
        double precision funcD2, funcD1
        integer i
        parameter (EPS=1.D-05)

c  CALCUL DE D1
        i = 0
        D11 = 100.D10
        D1 = 0.D0
        do while (abs(D11).ge.EPS*D1)
        call qromb(funcD1,10.D0**i, 10.D0**(i+1), D11, 1.D-06)
        D1 = D1 + D11
        i=i-1
        enddo
        i = 1
        D11 = 100.D10
        do while (abs(D11).ge.EPS*D1)
        call qromb(funcD1,10.D0**i, 10.D0**(i+1), D11, 1.D-06)
        D1 = D1 + D11
        i=i+1
        enddo
c  CALCUL DE D2
        i = 0
        D22 = 100.D0
        D2 = 0.D0
        do while (abs(D22).ge.EPS*D2)
        call qromb(funcD2,10.D0**i, 10.D0**(i+1), D22, 1.D-06)
        D2 = D2 + D22
        i=i-1
        enddo
        i = 1
        D22 = 100.D0
        do while (abs(D22).ge.EPS*D2)
        call qromb(funcD2,10.D0**i, 10.D0**(i+1), D22, 1.D-06)
        D2 = D2 + D22
        i=i+1
        enddo
c  CALCUL DE D
        D = D1/D2
c  CALCUL DE R
        R = (3.D-03 + 1.92D-03*U12)/3.16D-03/U12
c  CALCUL DE c
        c = 2.D0*(1.D0-R)/(1.D0+R)/(1.D0-D)
        
        end

        function funcD1(x)
        implicit none
        external S_DV
        double precision x,S_DV, funcD1

        funcD1 = x*x*S_DV(x)*dexp(- 1.5D-04*x*x)
        return
        end

        function funcD2(x)
        implicit none
        external S_DV
        double precision x,S_DV, funcD2

        funcD2 = x*x*S_DV(x)
        return
        end

