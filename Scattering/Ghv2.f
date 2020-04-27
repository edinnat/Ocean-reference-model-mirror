        subroutine Ghv2_ (Ghv2)

        implicit none
                
        common /Z/       Z
        common /Z2/      Z2
        common /C1/      C1
        common /D1/      D1
        common /D2/      D2
        common /Dc1/     Dc1
        common /Dc2/     Dc2
        common /epsi/    epsi
        common /epsi_/   epsi_
        common /sin1/    sin1
        common /cos1/    cos1
        common /cosdphi/ cosdphi
        common /sindphi/ sindphi

        double precision cosdphi
        double precision sindphi
        double precision sin1
        double precision cos1
        double precision Z
        double precision Z2
        
        double complex C1
        double complex D1
        double complex D2
        double complex Dc1
        double complex Dc2
        double complex epsi
        double complex epsi_
        double complex Ghv2

        
        Ghv2 = -2.0D0*cos1*epsi_*sindphi/D1/D2/Dc1*(epsi*Z*sin1-
     &         epsi_*Z2*C1*cosdphi/Dc2)
        
        return
        
        end
