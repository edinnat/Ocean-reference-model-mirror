c-----------------------------------------------------------------------
c  Tabulation du spectre de puissance de la houle
c
c     h         : rms de la hauteur de la houle (m)
c     sigma(1 ou 2) : largeur à mi-hauteur du spectre en x ou y
c     KMax(1 ou 2)  : maximum de puissance du spectre en x ou y
c     N(1 ou 2) : nombre de points de tabulation du spectre de la houle
c              dans les direction x et y
c     Swell : spectre de puissance des hauteur tabulé (m^4)
c-----------------------------------------------------------------------

       subroutine TabSwell (h, sigma, KMax, N, Swell)

       implicit none

       double precision h
       double precision sigma(1:2)
       double precision KMax(1:2)
       double precision Swell(1:1000,1:1000)
       double precision kx, ky
       double precision dkx, dky
       double precision kxmin, kymin
       double precision kxmax, kymax
       double precision pi

       integer N(1:2)
       integer ikx, iky

       kxmin = -5.0D0*sigma(1) + KMax(1) ! borne inf = -5x variance
       kxmax =  5.0D0*sigma(1) + KMax(1) ! borne sup = +5x variance
       kymin = -5.0D0*sigma(2) + KMax(2)
       kymax =  5.0D0*sigma(2) + KMax(2)
       dkx = (kxmax-kxmin)/(N(1) - 1) ! increment
       dky = (kymax-kymin)/(N(2) - 1)
       pi = dacos(-1.0D0)

       do ikx = 1, N(1)
          kx = kxmin + (ikx-1)*dkx
          do iky = 1, N(2)
             ky = kymin + (iky-1)*dky
             Swell(ikx, iky) = h*h/2.0D0/sigma(1)/sigma(2)/pi
     &         *dexp(-0.5D0*(((kx - KMax(1))/sigma(1))**2 +
     &                       ((ky - KMax(2))/sigma(2))**2 ))
          enddo
       enddo

       end
