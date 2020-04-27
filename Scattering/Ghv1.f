c________________________________________________________________________
c
c       Calcul du coefficient bistatique hv d'ordre 1.
c________________________________________________________________________

        subroutine Ghv1_(Ghv1)

        implicit none
        
        common /sindphi/ sindphi
        common /C3/      C3
        common /D1/     D1
        common /Di2/     Di2
        
        double precision sindphi
        double complex   C3
        double complex   D1
        double complex   Di2
        double complex   Ghv1

        Ghv1 = C3*sindphi/D1/Di2
        
        return

        end
