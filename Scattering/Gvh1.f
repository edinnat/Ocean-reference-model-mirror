c________________________________________________________________________
c
c       Calcul du coefficient bistatique vh d'ordre 1.
c       = coefficient from Yueh (1997)/(2 cos(theta_i))
c________________________________________________________________________

        subroutine Gvh1_(Gvh1)

        implicit none
        
        common /sindphi/ sindphi
        common /C1/      C1
        common /D2/      D2
        common /Di1/     Di1
        
        double precision sindphi
        double complex   C1
        double complex   D2
        double complex   Di1
        double complex   Gvh1

        Gvh1 = C1*sindphi/D2/Di1
        
        return

        end
