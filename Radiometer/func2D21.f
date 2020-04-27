        subroutine func2D21(x2, y2)

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
        double precision SLem
        double precision modu1
        double precision modu2
        double precision modu3
        double precision modu4
        double precision modu5
        double precision Bel
        double precision DeltaR
        double precision SpecOut(1:6) ! output spec Kudry
        
        double complex   Gvv1
        double complex   Ghh1
        double complex   Ghv1
        double complex   Gvh1
        double complex   Gvv2
        double complex   Ghh2
        double complex   Gvh2
        double complex   Ghv2
        double complex   Rvv0
        double complex   Rhh0
        double complex   C2
        double complex   epsi_
        
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

        if (flag2(1).ne.0) then
                call Gvv2_ (Gvv2)
                if (flagI.ne.0) then
                        call Gvv1_ (Gvv1)
                        call Gvh1_ (Gvh1)
                else
                        Gvv1   = 0.0D0
                        Gvh1   = 0.0D0
                endif
        else
                Gvv2 = 0.0D0
                Gvv1   = 0.0D0
                Gvh1   = 0.0D0
        endif
        if (flag2(2).ne.0) then
                call Ghh2_ (Ghh2)
                if (flagI.ne.0) then
                        call Ghh1_ (Ghh1)
                        call Ghv1_ (Ghv1)
                else
                        Ghh1   = 0.0D0
                        Ghv1   = 0.0D0
                endif
        else
                Ghh2 = 0.0D0
                Ghh1   = 0.0D0
                Ghv1   = 0.0D0
        endif
        
        modu1   = abs(Gvv1)
        modu2   = abs(Ghh1)
        modu3   = abs(Ghv1)
        modu4   = abs(Gvh1)

c Calcul du spectre, a partir du choix contenu dans cspec

        if (cspec.eq.1) then ! Spectre Durden et Vesecky
	    Azimuth = 1.0D0 + 
     &		c*(1.0D0-dexp(-1.5D-04*k*k))*dcosd(2.0D0*phi)
	    W =  0.15915494312D0/k*S_DV(k)*Azimuth 
        elseif (cspec.eq.2)then ! Spectre Elfouhaily
                call W_El (k, phi, W, Bel, DeltaR)
        elseif (cspec.eq.4)then ! Spectre Lemaire
               call spectreWLemaire (k, phi, W, SLem, DeltaR)
        elseif (cspec.eq.5)then ! Spectre Kudryavtsev
!!! Compute half of omnidirectional value !!!
!!! it is compute at phi = 0 and phi = 90 and then averaged
!!! so the average of 2 times 0.5 omni => omni
!!! This is because KM spectrum has 1 order component
!!! => omnidirectional component can't be derived from
!!! average of phi = 0 and phi = 45 !!!
!!! Omni for Kudry ~ 1/4(W@0 + 2*W@90 + W@180)
!!! because W ~ W0 + W1 cos phi + W2 cos 2phi + o(phi3)
               call spectreWKud (k, 0.D0, SpecOut, Nout)
               W = SpecOut(3)*0.25D0
               call spectreWKud (k, 90.0D0, SpecOut, Nout)
               W = W + SpecOut(3)*0.5D0
               call spectreWKud (k, 180.0D0, SpecOut, Nout)
               W = W + SpecOut(3)*0.25D0
        else 
                print *, 'func2D21 : cspec = ', cspec
                stop
        endif
c Spectre => W

        y2(1) = (2.0D0*realpart(conjg(Rvv0)*Gvv2)
     &          +(modu1*modu1+modu4*modu4)
     &          *cos1*C2*modu5)*1.7453293D-02*W
        y2(2) = (2.0D0*realpart(conjg(Rhh0)*Ghh2)
     &          +(modu2*modu2+modu3*modu3)
     &          *cos1*C2*modu5)*1.7453293D-02*W
         
c        print*, k, y2(1)
c	print*, krho1, phi_1, krho2, x2, Rvv0, Gvv2,modu1, modu4  




        return
        
        end
