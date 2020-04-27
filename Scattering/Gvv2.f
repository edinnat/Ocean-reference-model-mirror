        subroutine Gvv2_(Gvv2) 

        implicit none

        common /epsi/     epsi
        common /epsi_/    epsi_
        common /cos1/     cos1
        common /sin1/     sin1
        common /cosdphi/  cosdphi
        common /Z/        Z
        common /Z2/       Z2
        common /C1/       C1
        common /D2/       D2
        common /Dc1/      Dc1
        common /Dc2/      Dc2

        double precision  cos1
        double precision  sin1
        double precision  cosdphi
        double precision  Z
        double precision  Z2
        
        double complex    epsi
        double complex    epsi_
        double complex    C1
        double complex    D2
        double complex    Dc1
        double complex    Dc2
        double complex    Gvv2

  
        Gvv2 = -2.0D0*cos1*epsi_*epsi/D2/D2*(epsi_*Z2*sin1*sin1/Dc1/Dc2+
     &          C1*(1.0D0-2.0D0*Z*sin1*cosdphi/Dc1)-(epsi-sin1*sin1)
     &          *epsi_/epsi/Dc2*(1.0D0-Z2*cosdphi*cosdphi/Dc1)) 

        return
        
        end
