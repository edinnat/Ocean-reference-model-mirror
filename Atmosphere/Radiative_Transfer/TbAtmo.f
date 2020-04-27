c Tabulation de la température de brillance de l'atmosphère (TbA) en fonction de l'angle d'incidence (Theta)

        subroutine TbAtmo(Param, nParam, TbAd, TbAu, tau)

        implicit none

        double precision Param(1:10)
        double precision TbAd(0:90) ! Tableau des températures de brillance de l'atmosphere vers le bas en fonction de l'angle d'incidence (K)
        double precision TbAu(0:90) ! Tableau des températures de brillance de l'atmosphere vers le haut en fonction de l'angle d'incidence (K)
        double precision P0 ! Pression au sol (mb)
        double precision T0 ! Température au sol (K)
        double precision R  ! Humidité relative (%)
        double precision Rt ! Rayon terrestre (km)
        double precision alt_z !
        double precision f !fréquence (Hz)
        double precision altmax ! altitude max de l'atmosphère (km)
        double precision KO2(1:10000) ! coef absorption Oxygene(Neper/km)
        double precision KH2O(1:10000)!coef absorption Vapeur d'eau(N/km)
        double precision P(1:10000) ! profil de pression (mb)
        double precision T(1:10000) ! profil de température (°C)
        double precision q(1:10000)      ! profil humidite absolue (g/m3)
        double precision Tbup ! Température rayonnée vers le haut (K)   
        double precision Tbdown ! Température rayonnée vers le bas (K)
        double precision tauP           ! tau' epaisseur optique "partielle" (de la couche i a la limite) (Néper)
        double precision tau(0:90)
        double precision dalt
        double precision dpath(1:10000)!epaisseur d'atmo elementaire (km)
        double precision path        ! epaisseur d'atmo en 1 point (km)
        double precision Theta
        double precision secTheta
        double precision qtot ! Humidite absolue integree (mm)
        double precision grad ! gradient de T dans la troposphere(K/km)

        integer ncouche
        integer i
        integer j
        integer nParam

        character*80 opt ! si opt = 'f' profil de temperature a partir un fichier
        character*80 fichier ! fichier contenant profil de temperature

        

        Rt = Param(1)
        P0 = Param(2)
        T0 = Param(3)
        R  = Param(4)
        f  = Param(5)
        altmax = Param(6)
        ncouche = int(Param(7))
        grad = Param(8)
        opt = 't'     ! si opt = 'f' profil de temperature a partir 
                      ! d'un fichier. Sinon profil lineaire de pente 
                      ! '-grad'
        fichier='temperatures.txt'
        write (*,*)
        print*, ' --> Tabulation des Tb atmosphériques'
        write (*,*)

c     calcul les profils de pression, température, coef. d'abosrption
c      de l'oxygene et de la vapeur d'eau et de l'humidite absolue.
        call param_atmo (f, T0, P0,R, altmax,
     &         opt, fichier, grad, ncouche, P, T, KO2, KH2O, q)

c >
        do theta = 0, 90
c construction de dpath
c     dpath = elt d'atmosphere traversee qd on fait varier l'altitude de
c             'altmax/ncouche' et que l'on regarde dans la direction 
c             theta. dpath prend en compte la courbure de la terre de 
c             rayon Rt. Il tend vers sec(theta) = 1/cos(theta) qd theta
c             tend vers 0.

c >>
        do j = 1, ncouche
               dpath(j) = path(Rt, theta,altmax/ncouche*(j+1))
     &                   - path(Rt, theta, altmax/ncouche*j)
        enddo ! j
c <<

c -------- integrale Tb down  ----------------------

        Tbdown= 0.0D0 ! initialisation

c >>
        do j = 2, ncouche
                Tbdown = Tbdown + (T(j)+273.15)*(KO2(j)+KH2O(j))
     &               *dexp(-tauP(j,1,dpath,KH2O,KO2))*dpath(j)
        enddo ! j
c <<
c     Contribution des extremité (sol et haut de l'atmosphere pondérées
c     par un facteur 1/2 (methode des trapèzes)
        Tbdown = Tbdown + (T(1)+273.15)*(KO2(1)+KH2O(1))
     &          /2.0D0*dpath(1)
        Tbdown = Tbdown + (T(ncouche+1)+273.15)*(KO2(ncouche+1)
     &          +KH2O(ncouche+1))
     &          *dexp(-tauP(1,ncouche+1,dpath,KH2O,KO2))/2.0D0
     &          +2.7*dexp(-tauP(1,ncouche+1,dpath,KH2O,KO2))

c -------- integrale Tb upward  ----------------------

        Tbup = 0.0D0 ! initialisation

c >>

        do j = 2, ncouche
                Tbup = Tbup + (T(j)+273.15)*(KO2(j)+KH2O(j))
     &        *dexp(-tauP(j,ncouche+1,dpath,KH2O,KO2))*dpath(j)
        enddo ! j
c <<
c     Contribution des extremité (sol et haut de l'atmosphere pondérées
c     par un facteur 1/2 (methode des trapèzes)
        Tbup = Tbup + (T(1)+273.15)*(KO2(1)+KH2O(1))
     &          *dexp(-tauP(1,ncouche+1,dpath,KH2O,KO2)*dpath(1))/2.0D0
     &          *dpath(1)
        Tbup = Tbup + (T(ncouche+1)+273.15)*(KO2(ncouche+1)
     &          +KH2O(ncouche+1))/2.0D0

c ----------- integrale humidite totale qtot  --------
c >>
        do j = 2, ncouche+1
                qtot = qtot + (q(j)+q(j-1))/2.0D0*dpath(j-1)
        enddo ! j
c <<

        TbAd(theta) = Tbdown
        TbAu(theta) = Tbup
        tau(theta) = tauP(1,ncouche+1,dpath,KH2O,KO2)
        enddo ! theta
	dalt = altmax/ncouche
c	do j = 1, ncouche+1
c		print*, dalt*(j-1), T(j), P(j), q(j)
c	enddo
c <
        print*, ' --> fin de la tabulation des Tb atmosphériques'

        end ! fin de programme
