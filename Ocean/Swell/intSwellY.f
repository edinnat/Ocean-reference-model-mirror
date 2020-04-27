c-----------------------------------------------------------------------
c  Integration du spectre des pentes le long de la direction du vent
c
c    sigma(1 ou 2) : largeur à mi-hauteur du spectre en x ou y
c    KMax(1 ou 2)  : maximum de puissance du spectre en x ou y
c    N(1 ou 2) : nombre de points de tabulation du spectre de la houle
c              dans les direction x et y
c    Swell : densite du spectre de puissance des hauteurs (m^4)
c    Var   : variance des pentes en y (dir normal au vent)
c-----------------------------------------------------------------------

       subroutine intSwellY (sigma, KMax, N, Swell, Var)

       implicit none
       
       double precision KMax(1:2)
       double precision sigma(1:2)
       double precision Swell (1:1000, 1:1000)
       double precision Var
       double precision kx, ky
       double precision dkx, dky
       double precision kxmin, kymin
       double precision kxmax, kymax

       integer N(1:2)
       integer ikx, iky

       kxmin = -5.0D0*sigma(1) + KMax(1)
       kxmax =  5.0D0*sigma(1) + KMax(1)
       kymin = -5.0D0*sigma(2) + KMax(2)
       kymax =  5.0D0*sigma(2) + KMax(2)
       dkx = (kxmax-kxmin)/(N(1) - 1)
       dky = (kymax-kymin)/(N(2) - 1)

       Var = 0.0D0
       do ikx = 1, N(1)
          kx = kxmin + (ikx-1)*dkx
          do iky = 1, N(2)
             ky = kymin + (iky-1)*dky
             Var = Var + Swell(ikx,iky)*ky**2*dkx*dky
             if ((kx**2+ky**2).gt.10) print*, '!!!'
          enddo
       enddo

       end
