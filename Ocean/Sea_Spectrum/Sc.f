c______________________________________________________________________
c       Calcul de la variance des pentes dans la direction 
c       orthogonale à celle du vent.
c______________________________________________________________________

        subroutine Sc_ (Sc, kmin, kd, cspec)

        implicit none
       
        external funcScEl
        external funcScDV
        external funcScPL
        external funcScLem
        double precision Sc 
        double precision EPS 
        double precision kmin
        double precision kd
        integer i
        integer cspec

        parameter (EPS=1.D-06)

        Sc = 0.D0
        
c Integration du spectrum des pentes en crosswind "funcSuX" de 1e-03 à kd
c X = DV pour le spectrum de Durden & Vesecky
c X = El pour le spectrum de Elfouhaily

        if (cSpec.eq.1) then
                call qromb(funcScDV, kmin, kd, Sc, EPS)
        elseif (cspec.eq.2) then
                call qromb(funcScEl, kmin, kd, Sc, EPS)
        elseif (cspec.eq.3) then
                call qromb(funcScPL, kmin, kd, Sc, EPS)
        elseif (cspec.eq.4) then
                call qromb(funcScLem, kmin, kd, Sc, EPS)                
        else
                print *, 'Sc : cspec = ',cspec
        endif
        Sc = 5.D-01*sqrt(Sc)

        end

        function funcScDV (x)

        implicit none
        
        common /c/ c
        external S_DV

        double precision funcScDV
        double precision x
        double precision S_DV
        double precision c
        
        funcScDV = x*x*S_DV(x)*(2.D0 - c*(1.D0 - dexp(-1.5D-04*x*x)))
        return

        end

        function funcScEl (x)

        implicit none
        
        double precision funcScEl
        double precision x
        double precision Wel
        double precision Bel
        double precision DeltaR
        
        call W_EL (x, 0.0D0, Wel, Bel, DeltaR)
        funcScEl = Bel/x*(2.0D0 - DeltaR)
        
        return

        end

        function funcScPL (x)

        implicit none
        
        double precision funcScPL
        double precision x
        double precision Wel
        double precision Bel
        double precision S_powerLaw
        double precision DeltaR
        
        call W_EL (x, 0.0D0, Wel, Bel, DeltaR)
        funcScPL = x*x*S_powerLaw(x)*(2.0D0 - DeltaR)
        return

        end

        function funcScLem (x)

        implicit none
        
        double precision funcScLem
        double precision x
        double precision WLem
        double precision SLem
        double precision DeltaR
        
        call SpectreWLemaire (x, 45, WLem, SLem, DeltaR)
        funcScLem = SLem*x**2*(2.0D0 - DeltaR)
        
        return

        end
