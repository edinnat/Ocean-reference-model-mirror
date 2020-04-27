c calcul de la couverture d'ecume active et passive d'apres Monahan et Lu
c
c calcul de la viscosite cinematique fn de la T selon la formule donnee par
c E. Obligis (origine?) en accord avec la fig. 6 de Monahan et Muircheartaigh
c (1986)
c
      subroutine MonahanLu (U, T, Wa, Wb, W)
      
      implicit none

      double precision T ! temperature en degC
      double precision U ! vent a 10m en m/s
      double precision Fr ! fraction d'ecume
      double precision rnu ! viscosité cinématique      
      double precision Wa ! fraction d ecume active      
      double precision Wb ! fraction d ecume passive      
      double precision W ! fraction d ecume totale
      double precision Wb1 !       
c
        rnu=0.018D0-4.72D-4*T+5.24D-6*T*T
c
c Ecume passive:
        Wb1=5.21D-4*(U*((9.81D0*rnu*0.0001D0)**(-0.3333D0))-75.9D0)
        Wb=Wb1**3
        Wa=Wb/9.7D0
        W = Wa + Wb

      end
