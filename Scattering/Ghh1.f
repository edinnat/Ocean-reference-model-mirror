c________________________________________________________________________
c
c       Calcul du coefficient bistatique hh d'ordre 1.
c       = coefficient from Yueh (1997)/(2 cos(theta_i))/(epsi-1)
c________________________________________________________________________

        subroutine Ghh1_(Ghh1)

        implicit none
        
        common /cosdphi/ cosdphi
        common /D1/      D1
        common /Di1/     Di1
        
        double precision cosdphi
        double complex   D1
        double complex   Di1
        double complex   Ghh1

        Ghh1 = cosdphi/D1/Di1

        return
        end
