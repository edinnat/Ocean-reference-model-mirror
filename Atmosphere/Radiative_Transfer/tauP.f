c-----------------------------------------------------------------------
c       Calcul l'epaisseur optique entre la couche 
c       icouche et la couche limite couchlim  
c-----------------------------------------------------------------------

        function tauP(icouche, couchlim, dpath, KH2O, KO2)

        implicit none

        double precision tauP
        double precision KO2(*) ! coef absorption oxygene
        double precision KH2O(*)! coef absorption vapeur eau
        double precision theta  ! angle de visée
        double precision dpath(*)! epaisseur atmosphere elementaires tabulée
        integer icouche
        integer couchlim
        integer i
        
        tauP = 0.0D0
        if (icouche.ne.couchlim) then
        do i = min(icouche,couchlim), max(icouche-1,couchlim-1)
                tauP = tauP + 
     &          (KH2O(i)+KO2(i)+KH2O(i+1)+KO2(i+1))/2.0D0*dpath(i)
        enddo
        endif

        end
