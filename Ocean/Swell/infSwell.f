c-----------------------------------------------------------------------
c  Calcul la variance des pentes de la houle
c
c     Var(1)    : variance des pentes dans la direction du vent (dir x)
c     Var(2)    : variance des pentes normales au vent (dir y)
c     N(1 ou 2) : nombre de points de tabulation du spectre de la houle
c              dans les direction x et y
c     h         : rms de la hauteur de la houle (m)
c     sigma(1 ou 2) : largeur à mi-hauteur du spectre en x ou y
c     KMax(1 ou 2)  : maximum de puissance du spectre en x ou y
c-----------------------------------------------------------------------

       subroutine infSwell(N, h, sigma, KMax, Var)

       implicit none

       double precision h ! rms des hauteurs de la houle
       double precision sigma(1:2) ! largeur a mi-hauteur
       double precision KMax(1:2) ! pics de la distribution de la houle
       double precision Swell(1:1000,1:1000) ! Spectre de puissance
       double precision Var(1:2) ! Variance des pentes la houle
       double precision pi

       integer N(1:2) ! nombre de points de tabulation de la houle

       pi = acos(-1.0D0)

       call TabSwell (h, sigma, KMax, N, Swell) ! Tabulation du spectre de puissance de la houle
       call intSwellX (sigma, KMax, N, Swell, Var(1)) ! integration de la projection du spectre des hauteur sur kx
       call intSwellY (sigma, KMax, N, Swell, Var(2)) ! integration de la projection du spectre des hauteur sur ky

       end
