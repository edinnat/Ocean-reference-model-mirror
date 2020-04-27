c______________________________________________________________________
c       Calcul de la variance des pentes dans la direction 
c       longitudinale à celle du vent.
c______________________________________________________________________

        subroutine Su_ (Su, kmin, kd, cspec)

        implicit none
       
        external funcSuEl
        external funcSuDV
        external funcSuPL
        external funcSuLem
        double precision Su
        double precision EPS
        double precision kmin
        double precision kd
        integer i
        integer cspec

        parameter (EPS=1.D-06)

        Su = 0.D0

c Integration du spectrum des pentes en upwind "funcSuX" de 1e-03 à kd
c X = KS pour le spectrum de Durden & Vesecky
c X = El pour le spectrum de Elfouhaily
        if (cspec.eq.1) then
                call qromb(funcSuDV,kmin, kd, Su, EPS)
        elseif (cspec.eq.2) then
                call qromb(funcSuEl,kmin, kd, Su, EPS)
        elseif (cspec.eq.3) then
                call qromb(funcSuPL,kmin, kd, Su, EPS)
        elseif (cspec.eq.4) then
                call qromb(funcSuLem,kmin, kd, Su, EPS)
        else
                print *, 'Su : cspec ', cspec
        endif
        Su = 5.D-01*sqrt(Su)

        end

        function funcSuDV (x)

        implicit none
        
        common /c/ c
        external S_DV


        double precision funcSuDV
        double precision x
        double precision S_DV
        double precision c

        funcSuDV = x*x*S_DV(x)*(2.D0 + c*(1.0D0 - dexp(-1.5D-04*x*x)))
        return

        end

        function funcSuEl (x)

        implicit none
        
        double precision funcSuEl
        double precision x
        double precision Wel
        double precision Bel
        double precision DeltaR

        call W_EL (x, 0.0D0, Wel, Bel, DeltaR)
        funcSuEl = Bel/x*(2.D0 + DeltaR)
        return

        end

        function funcSuPL (x)

        implicit none
        
        double precision funcSuPL
        double precision x
        double precision Wel
        double precision Bel
        double precision S_powerLaw
        double precision DeltaR

        call W_EL (x, 0.0D0, Wel, Bel, DeltaR)
        funcSuPL = x*x*S_powerLaw(x)*(2.D0 + DeltaR)
        return

        end

        function funcSuLem (x)

        implicit none
        
        double precision funcSuLem
        double precision x
        double precision WLem
        double precision SLem
        double precision DeltaR

        call SpectreWLemaire (x, 45, WLem, SLem, DeltaR)
        funcSuLem = SLem*x**2*(2.D0 + DeltaR) 
        return

        end

