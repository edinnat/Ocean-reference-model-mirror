c       Calcul l'epaisseur d'atmosphere traversee 'path' pour un point
c       visé à l'altitude 'z' dans la direction theta

        function path (R, theta, z)

        implicit none
        
        external dsind
        external dcosd
        
        double precision path  ! epaisseur atmosphere(km)
        double precision R      ! Rayon terrestre(km)
        double precision theta  ! angle d'incidence du satellite(°)
        double precision thetp  ! complementaire de theta à pi
        double precision z      ! altitude du point visé(km)
        double precision dsind
        double precision dcosd
        double precision phi
        double precision s

        double complex pathtab(1:4)! ensemble des solutions pour h
        double complex C(0:4)   ! Coefficients du polynome à annuler

        integer indice

        logical verif

c	if (theta.eq.90.0D0) then
c		path = sqrt(R*R + (R+z)**2)
c		return
c	endif

        path = 1.0D10
        thetp = 180.0D0 - theta
        !print*, z, R
        
        C(0) = 4.0D0*R*z*z*(z + R) + z**4
        C(1) = 0.0D0
        C(2) = -2.0D0*((R+z)**2 + R*R - 2.0D0*R*R*dsind(thetp)**2)
        C(3) = 0.0D0
        C(4) = 1.0D0

        call zroots(C, 4, pathtab, .false.)
        do indice = 1,4
                if (verif(R, z, thetp, pathtab(indice))) then
                        if (path.eq.1.0D10) then
                                path = realpart(pathtab(indice))
                        else
!                                print*, '------- ERREUR  -------------'
!                                print*, 'plusieurs solutions pour path'
!                                print*, path
                                if (realpart(pathtab(indice)).lt.path)
     &                                path = realpart(pathtab(indice))
!                                print*, path
                        endif
                endif
        enddo

        if (path.eq.1.0D10) then
                print*, '------ ERREUR  ---------------'
                print*, 'Pas de solution valable pour L'
                stop
        endif
        return 
        end

        function verif(R, z, thetp, path)

        implicit none

        external dcosd
        external dsind
        
        double precision R
        double precision z
        double precision L
        double precision thetp
        double precision phi
        double precision dcosd
        double precision dsind
        double precision acos1
        double precision asin1
        
        double complex path
        
        logical    verif

        verif = .false. 
c        if (imagpart(path).ne.0.0D0) return ! solution doit etre reelle
        L = realpart(path)      ! on prend la partie reelle de solution
        if (L.lt.0.0D0) return
c        if (abs((R*R + (R+z)**2 - L*L)/2/R/(R+z)).le.1.0D0) then
                acos1 = dacos((R*R + (R+z)**2 - L*L)/2/R/(R+z))
c        else
c                return
c        endif
c        if (abs(L*dsind(thetp)/(R+z)).le.1.0D0) then
                asin1 = dasin(L*dsind(thetp)/(R+z))
c        else
c                return
c        endif
c        if (abs(acos1-asin1).gt.1.0D-06) return ! phi independant de la methode
c        phi = asin1
c        if (abs(R*R - ((R+z)**2 + L*L - 2.0D0*(R+z)*L
c     &  *dcosd(180-thetp-phi))).gt.1.0D-06) return ! R doit etre coherent avec 3eme angle(=phi)
        verif = .true.
        
        return
        end
