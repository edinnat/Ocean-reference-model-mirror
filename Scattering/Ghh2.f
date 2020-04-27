        subroutine Ghh2_(Ghh2) 

        implicit none

        common /cos1/     cos1
        common /cosdphi/  cosdphi
        common /epsi_/    epsi_
        common /Z2/       Z2
        common /C1/       C1
        common /C2/       C2
        common /C3/       C3
        common /D1/       D1
        common /Dc1/      Dc1
        common /Dc2/      Dc2

        double precision  cos1
        double precision  cosdphi
        double precision  Z2

        double complex    epsi_
        double complex    C1
        double complex    C2
        double complex    C3
        double complex    D1
        double complex    Dc1
        double complex    Dc2
        double complex    Ghh2

        
        Ghh2 = 2.0D0*cos1*epsi_/D1/D1*(C1-epsi_/Dc1/Dc2*(C2*C3+
     &          Z2*cosdphi*cosdphi))

        return
        
        end
