c________________________________________________________________________
c
c       Calcul du coefficient bistatique vv d'ordre 1.
c       = coefficient from Yueh (1997)/(2 cos(theta_i))
c________________________________________________________________________

        subroutine Gvv1_(Gvv1)

        implicit none
        
        common /cosdphi/ cosdphi
        common /epsi/    epsi
        common /sin1/    sin1
        common /Z/       Z
        common /C1/      C1
        common /C3/      C3
        common /D2/      D2
        common /Di2/     Di2
        
        double precision cosdphi
        double complex   epsi
        double precision sin1
        double precision Z
        double complex   C1
        double complex   C3
        double complex   D2
        double complex   Di2
        double complex   Gvv1

        Gvv1 = (epsi*sin1*Z - C1*C3*cosdphi)/Di2/D2
        
        return

        end
