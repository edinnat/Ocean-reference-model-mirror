        subroutine func2D22(x2, y2)

        implicit none

        external dcosd
        external dsind
        external dacosd
c Passage des variables
        common   /kd/           kd
        common   /c/            c
        common   /s/            s
        common   /Rvv0/         Rvv0
        common   /Rhh0/         Rhh0
        common   /C2/           C2
        common   /cos1/         cos1
        common   /epsi_/        epsi_
        common   /phi_1/        phi_1
        common   /krho1/        krho1
        common   /krho2/        krho2
        common   /flag2/        flag2
        common   /cosdphi/      cosdphi
        common   /sindphi/      sindphi
        common   /flagI/        flagI
        common   /cspec/        cspec
c Declaration des variables
        double precision kd
        double precision dcosd
        double precision dsind
        double precision dacosd
        double precision dphi
        double precision cosdphi
        double precision sindphi
        double precision c
        double precision s
        double precision S_DV
        double precision y2(1:2)
        double precision x2
        double precision cos1
        double precision krho1
        double precision krho2
        double precision phi_1
        double precision k
        double precision phi
        double precision Azimuth
        double precision W
        double precision Wel
        double precision modu5
        double precision Bel
        double precision SLem
        double precision DeltaR
        double precision SpecOut(1:6) ! output Spec Kudry
        
        double complex   Gvv1
        double complex   Ghh1
        double complex   Ghv1
        double complex   Gvh1
        double complex   Gvh2
        double complex   Ghv2
        double complex   Rvv0
        double complex   Rhh0
        double complex   C2
        double complex   epsi_
        double complex   conj1
        double complex   conj2
        
        integer          indice
        integer          flag2(1:2)
        integer          flagI
        integer          cspec
        integer          Nout

        do indice = 1, 2
                y2(indice) = 0.0D0
        enddo
        
        modu5   = abs(epsi_)
        modu5   = 4.0D0*modu5*modu5
        k       = sqrt(krho1*krho1 + krho2*krho2 - 
     &            2.0D0*krho1*krho2*dcosd(phi_1 - x2))
        phi     = dsign (dacosd((krho1*dcosd(phi_1)-krho2*dcosd(x2))/
     &            k*0.9999999D0),krho1*dsind(phi_1) - krho2*dsind(x2))
        dphi    = phi_1 - x2
        cosdphi = dcosd (dphi)
        sindphi = dsind (dphi)

c  CALCUL DES COEFFICEINTS COHERENTS 

        call Ghv2_ (Ghv2)
        Gvh2 = -Ghv2
        
        if (flagI.ne.0) then
                call Gvv1_ (Gvv1)
                call Ghh1_ (Ghh1)
                call Ghv1_ (Ghv1)
                call Gvh1_ (Gvh1)
        else
                Gvv1   = 0.0D0
                Ghh1   = 0.0D0
                Ghv1   = 0.0D0
                Gvh1   = 0.0D0
        endif

        conj1 = conjg(Ghh1)
        conj2 = conjg(Ghv1)

        if (cspec.eq.1) then
        Azimuth = 1.0D0 + c*(1.0D0-dexp(-1.5D-04*k*k))*dcosd(2.0D0*phi)
        W =  0.15915494312D0/k*S_DV(k)*Azimuth
        elseif (cspec.eq.2) then
                call W_El (k, phi, W, Bel, DeltaR)
        elseif (cspec.eq.4) then
                call spectreWLemaire(k, phi, W, SLem, DeltaR)
        elseif (cspec.eq.5) then
                call spectreWKud(k, phi, SpecOut, Nout)
                W = SpecOut(3)
        else
                print *,'func2D22 : cspec = ',cspec
                stop
        endif

        y2(1) = 2.0D0*realpart((Gvh2*conjg(Rhh0)+Rvv0*conjg(Ghv2))
     &          +(Gvh1*conj1+Gvv1*conj2)*cos1*C2*modu5
     &          )*1.7453293D-02*W
        y2(2) = 2.0D0*imagpart((Gvh2*conjg(Rhh0)+Rvv0*conjg(Ghv2))
     &          +(Gvh1*conj1+Gvv1*conj2)*cos1*C2*modu5
     &          )*1.7453293D-02*W
c        y2(1) = 0.0D0
c        y2(2) = 0.0D0
c        y2(1) = 2.0D0*realpart((conjg(Rvv0)-conjg(Rhh0))*Ghv2
c     &          +(Gvh1*conj1+Gvv1*conj2)*cos1*C2*modu5
c     &          )*1.7453293D-02*W
c        y2(2) = 2.0D0*imagpart((-conjg(Rhh0)-conjg(Rvv0))*Ghv2
c     &          +(Gvh1*conj1+Gvv1*conj2)*cos1*C2*modu5
c     &          )*1.7453293D-02*W
c        y2(1) = 2.0D0*realpart(!(Gvh2*conjg(Rhh0)+Rvv0*conjg(Ghv2))
c     &          (Gvh1*conj1+Gvv1*conj2)*cos1*C2*modu5
c     &          )*1.7453293D-02*W
c        y2(2) = 2.0D0*imagpart(!(Gvh2*conjg(Rhh0)+Rvv0*conjg(Ghv2))
c     &          (Gvh1*conj1+Gvv1*conj2)*cos1*C2*modu5
c     &          )*1.7453293D-02*W

        return
        
        end
