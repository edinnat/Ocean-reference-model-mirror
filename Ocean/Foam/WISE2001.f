c calcul de la couverture d'ecume active et passive d'apres 
c les mesures d'A. Weill sur WISE 2001
c
      subroutine WISE2001 (U, W)
      
      implicit none

      double precision U ! vent a 10m en m/s
      double precision W ! fraction d ecume totale

        W =  10.0D0**(-5.2076)*U**2.8788

      end
