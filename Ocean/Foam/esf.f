       subroutine esf(theta, nu, epsi_sf)

       implicit none

       double precision theta, nu, epsi_sf(2)
       double precision epsi0Tw, Fp(2)

       ! theta : incidence angle in degrees
       ! nu : frequency in GHz
       ! epsi_sf : foam brightness temperature in K for V (1) and H(2)
       !                pol

       epsi0Tw = 2.08D02 + 1.29D0*nu
c  Calcul du facteur de forme vertical   
       Fp(1) = 1.D0-9.946D-04*theta+3.218D-05*theta**2
     &         -1.187D-06*theta**3+7.D-20*theta**10 
c  Calcul du facteur de forme horizontal   
       Fp(2) = 1.-1.748D-03*theta-7.336D-05*
     &         theta**2+1.044D-07*theta**3

       epsi_sf(1) = epsi0Tw*Fp(1)
       epsi_sf(2) = epsi0Tw*Fp(2)

       end
