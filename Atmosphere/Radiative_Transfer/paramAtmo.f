! Calcul les profil de pression, temperature et coef absorption

        subroutine param_atmo(f, T0, P0,R, altmax, 
     &                          opt, fichier, grad,
     &                          ncouche, P, T, KO2, KH2O, q)

        implicit none

        double precision f                                              ! frequence (Hz)
        double precision nu                                             ! frequence (GHz)
        double precision KO2(1:10000)                                         ! profil du coef d'absorption de O2 (neper/km)
        double precision KH2O(1:10000)                                        ! profil du coef d'absorption de H2O (neper/km)
        double precision T0                                             ! T0 temperature a la surface de la mer (°K)
        double precision P0                                             ! P0 pression a la surface de la mer (mbar)
        double precision R                                              ! Humidite relative (%) constante
        double precision altmax                                         ! altitude max de l'atmosphere (km)
        double precision altTab(1:10000)                                 ! tabulation altitudes des limites des couches
        double precision TTab (1:10000)                                  ! Temperature tabulees
        double precision altTTab(1:10000)                                ! altitudes des temperatures tabulees
        double precision T(1:10000)                                           ! profil de temperature (°C)
        double precision P(1:10000)                                           ! profil de pression (mbar) 
        double precision q(1:10000)                                           ! Profil d'humidite absolue
        double precision grad                                           ! gradient de T dans la troposphere(K/km)
        double precision dT
        double precision F0_O2(44)                                      ! Frequences raies O2 (GHz)
        double precision F0_H2O(35)                                     ! Frequences raies HO2 (GHz)
        double precision A_O2(6,44)                                     ! parametres raies O2
        double precision A_H2O(6,35)                                    ! parametres raies H2O
        integer nTab                                                    ! nombre de valeurs tabulee pour la temperature
        integer ncouche                                                 ! nb de couches pour les profils
                                                                        ! avec nombre altitudes =  ncouche + 1
        integer icouche                                                 ! numero de la couche atmospherique
        integer ii

        character*80 opt                                                ! si opt = 'f' profil ... 
                                                                        ! de temperature a partir un fichier
        character*80 fichier                                            ! fichier contenant profil de temperature
                                                                
        nu = f/1.0D09 

        do  icouche = 0, ncouche
                altTab (icouche) = altmax/ncouche*icouche
        enddo
        if ((opt.ne.'f').and.(opt.ne.'F')) then
                call profil_T(ncouche, altmax, T0, T, grad)             ! calcul le profil de temperature
        else    
                open (unit=1,file=fichier)
                read(1,*) nTab                                          ! lit le nombre de valeurs de T
                do ii = 1, nTab
                        read(1,*) altTTab(ii), TTab(ii)                 ! lit les altitudes et temperatures
                enddo
                close (1)
                do ii = 1, ncouche
                call polint(altTTab, TTab                               ! Calcul les temperature sur les ncouches
     &                          , nTab, altTab(ii), T(ii), dT)         
                T(ii) = T(ii) - 273.15D0
c                        write (*,*) altTab(ii), T(ii)
                enddo
        endif
        call profil_P(ncouche, altmax, P0, T0, P)                       ! calcul le profil de pression
        call lecture_tab(F0_O2,A_O2,F0_H2O,A_H2O)                       ! lit les tableaux de parametres raies
        call calcul_profil(nu, T, P, R, ncouche,q,                      ! calcul profil coef absoption
     &                  F0_O2,A_O2,F0_H2O,A_H2O,KO2,KH2O)    
        
        end


c ----------------   Profil de temperature  ------------------------
        subroutine profil_T(ncouche, altmax, T0, T, grad)

        implicit none
        
        double precision T(*)                                              
        double precision Tlim                                           ! temperature minimale (°K)
        double precision grad                                           ! taux de decroissance (°K/km)
        double precision altit                                          ! altitude (km)
        double precision T0
        double precision altmax
        integer ncouche
        integer ialt

        Tlim = 223.15D0                                                
        
        do ialt = 0, ncouche
                altit = altmax/ncouche*ialt
                T(ialt+1) = max(T0 - grad*altit, Tlim) - 273.15
        enddo
        
        end
c ----------------   Profil de pression  ------------------------
        subroutine profil_P(ncouche, altmax, P0, T0, P)

        implicit none
        
        double precision P(*)                                              
        double precision altit                                          ! altitude (km)
        double precision P0
        double precision T0
        double precision altmax
        double precision a
        double precision g
        double precision Ra
        double precision expo
        
        integer ncouche
        integer ialt

        a = 6.5D0 
        g = 9.80665                                                     ! acceleration de la pesanteur
        Ra=287
        expo = g/Ra/(0.001*a)                                           ! exposant pour la pression

        do ialt = 0, ncouche
                altit = altmax/ncouche*ialt
                P(ialt+1) = P0*(1-a*altit/T0)**expo;
c                P(ialt+1) = P0*dexp(-altit/7.7D0)
        enddo
        
        end
c ------------   Lecture des parametres des raies -------------
        subroutine lecture_tab(F0_O2,A_O2,F0_H2O,A_H2O)

        implicit none

        double precision F0_O2(*) 
        double precision F0_H2O(*) 
        double precision A_O2(6,*) 
        double precision A_H2O(6,*) 
		character*80 Fort
        integer i
        integer j

		call getenv('FORTRAN', Fort)
        open (unit=2,file=Fort(1:lnblnk(Fort))//'/Data/Atmosphere/WATER.
     &DAT')
        read(2,*) ! saute les 2 ligne de commentaires
        read(2,*)
        do i = 1, 35
                read (2,*) F0_H2O (i),(A_H2O(j,i),j=1,6)
        enddo
        close (2)
        open (unit=2,file=Fort(1:lnblnk(Fort))//'/Data/Atmosphere/OXYGEN
     &.DAT')
        read(2,*) ! saute les 2 ligne de commentaires
        read(2,*)
        do i = 1,44
                read (2,*) F0_O2(i),(A_O2(j,i),j=1,6)
        enddo
        close (2)

        end
c ------------  Calcul profils coef absorption  -------------
        subroutine calcul_profil(nu, T, P, R, ncouche, q
     &                          ,F0_O2,A_O2,F0_H2O,A_H2O,KO2,KH2O)

        implicit none

        double precision nu
        double precision T(*)
        double precision P(*)
        double precision R
        double precision F0_O2(*)
        double precision F0_H2O(*)
        double precision A_H2O(6,*)
        double precision A_O2(6,*)
        double precision KO2(*)
        double precision KH2O(*)
        double precision V
        double precision Sn
        double precision So
        double precision Gammao
        double precision Gamma
        double precision gamH
        double precision gamD2
        double precision delH
        double precision Delta
        double precision S
        double precision Pdry                                           ! Pression de l'air sec ?
        double precision E                                              ! Pression partielle de vapeur d'eau
        double precision Es
        double precision Y
        double precision q(*)                                           ! Profil d'humidite absolue
        
        double complex Zn
        double complex Zf
        double complex Zfo
        double complex Zfn
        
        integer ncouche
        integer ialt
        integer i

        
c   Iteration sur les couche de l'atmosphere
        do ialt = 0, ncouche
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c   Calcul la pression partielle de vapeur d'eau en mbar
        if (T(ialt+1).gt.-40.0D0) then
                Y = 373.16D0/(T(ialt+1)+273.16D0)
                Es = -7.90298D0*(Y - 1.0D0) + 5.02808D0*dlog10(Y)-
     &          1.3816D-7*(10**(11.344D0*(1.0D0 - (1.0D0/Y))) - 1.0D0)+
     &          8.1328D-3*(10**(-3.49149D0*(Y - 1.0D0)) - 1.0D0)
     &          +dlog10(1013.246D0)
        else
                Y = 273.16D0/(T(ialt+1) + 273.16D0)
                Es = -9.09718*(Y - 1.0D0) - 3.56654D0*dlog10(Y) +
     &          0.876793D0*(1.0D0 - (1.0D0/Y))+dlog10(6.1071D0)
        endif
        Es = 10**Es
c        Es = 2.408D11*V**5.0D0*dexp(-22.644*V)
        E = Es*R/100.0D0


        V = 300.0D0/(T(ialt+1) + 273.15D0)
        Pdry = P(ialt+1) - E
        if (Pdry.le.0) then
                E = P(ialt+1)
                Pdry = 0.0D0
        endif

        q(ialt+1) = 0.7223*E*v  ! humidite absolue de la couche ialt+1

        
c   ---         Oxygene          ------

        Zn = (0.0D0,0.0D0)
        do i = 1,44 ! somme contribution 44 raies O2
                Gamma = 0.0D0
                S = A_O2(1,i)*Pdry*V**3*exp(A_O2(2,i)*(1.0D0-V))*1.0D-6
                Gamma = A_O2(3,I)*(Pdry*V**(0.8D0 - A_O2(4,I))
     &                  + 1.1D0*E*V)*1.D-3
                Gamma = (Gamma**2 + (25.0D0*0.6D-4)**2)**0.5D0
                Delta = (A_O2(5,I) + 
     &                  A_O2(6,I)*V)*(Pdry + E)*(V**0.8D0)*1.0D-3
                Zf = nu/F0_O2(i)*(cmplx(1.0D0,-Delta)/
     &               cmplx(F0_O2(I)-nu,-Gamma)
     &               -cmplx(1.0D0, Delta)/cmplx(F0_O2(i) + nu, Gamma))
                Zn = Zn + S*Zf
        enddo
c  ---     raie absorption oxygen   ---
        KO2(ialt+1) = 0.182D0*nu*imagpart(Zn)
c   ne peut pas etre < 0
        if (KO2(ialt+1).lt.0.0D0) KO2(ialt+1) = 0.0D0

c   Continum air sec
        Zn = (0.0D0, 0.0D0)
        So = 6.14D-5*Pdry*V**2
        Gammao = 0.56D-3*(Pdry+E)*V**0.8D0
        Zfo = - nu/cmplx(nu,Gammao)
        Sn = 1.40D-12*Pdry**2*V**3.5D0
        Zfn = cmplx(0.D0, nu/(1.93D-5*nu**1.5D0 + 1.0D0))
        Zn = So*Zfo + Sn*Zfn
c   ajoute absoption non-resonnante air sec
        KO2(ialt+1) = KO2(ialt+1) + 0.182D0*nu*imagpart(Zn)
      
c   ---         Vapeur d'eau    -------   

        Zn = (0.0D0, 0.0D0)
        do i = 1,34 ! somme contribution 34 raies H2O
                gamH = 0.0D0
                S = A_H2O(1,i)*E*V**3.5D0*exp(A_H2O(2,i)*(1.0D0 - V))
c   approximation Doppler
                gamH = A_H2O(3,i)*(Pdry*V**A_H2O(5,i)
     &          +A_H2O(4,i)*E*V**A_H2O(6,i))*1.0D-3
                gamD2 = 1.0D-12/V*(1.46D0*F0_H2O(i))**2
                gamH = 0.535D0*gamH + (0.217D0*gamH**2 + gamD2)**0.5D0
                delH = 0.0D0
                Zf = nu/F0_H2O(i)*(cmplx(1.0D0,-delH)
     &               /cmplx(F0_H2O(i)-nu,-gamH)-
     &               cmplx(1.0D0,delH)/cmplx(F0_H2O(i) + nu, gamH))
                Zn = Zn + S*Zf
        enddo
C   Raie absorption vapeur d'eau
        KH2O(ialt+1) = 0.182D0*nu*imagpart(Zn)

        Zn = (0.0D0 , 0.0D0)
        do i = 35,35
                gamH = 0.0D0
                S = A_H2O(1,i)*E*V**3.5D0*exp(A_H2O(2,i)*(1.0D0 - V))
C   approximation Doppler
                gamH = A_H2O(3,i)*(Pdry*V**A_H2O(5,i)
     &                  +A_H2O(4,i)*E*V**A_H2O(6,i))*1.0D-3
                gamD2 = 1.0D-12/V*(1.46D0*F0_H2O(i))**2
                gamH = 0.535D0*gamH + (0.217D0*GAMH**2 + GAMD2)**0.5D0
                delH=0.0D0
                Zf = nu/F0_H2O(i)*(cmplx(1.0D0,-delH)
     &               /cmplx(F0_H2O(i)-nu,-GAMH)-
     &               cmplx(1.0D0,delH)/cmplx(F0_H2O(i) + nu,gamH))
                Zn = Zn + S*Zf
        enddo
        
c   ajoute 35eme raie absorption H2O
        KH2O(ialt+1) = KH2O(ialt+1) + 0.182D0*nu*imagpart(Zn)

c        conversion en neper/km
        KH2O(ialt+1) = KH2O(ialt+1) * 0.1D0*dlog(10.0D0)
        KO2(ialt+1) = KO2(ialt+1) * 0.1D0*dlog(10.0D0)

c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        enddo !fin boucle couche atmo


        end
