c-----------------------------------------------------------------------
c          Programme principal du modele de temperature de brillance
c          en fonction de l'etat de surface, de la SST, de la SSS,
c          de l'angle d'incidence, de l'azimuth par rapport a la
c          direction du vent.
c          Modele de Yueh, 1997: Modeling if Wind Direction Signals
c          in Polarimetric Sea Surface Brightness Temperatures, IEEE
c          Transactions on Geoscience and Remote Sensing.
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c   * Uz est la valeur du vent(m/s) connue a l'altitude z(metres).
c   * phi et theta sont l'azimuth/direction du vent et l'angle 
c   * d'incidence (degrés)
c   * SSS est la salinite en PSU
c   * SST est la temperature de surface en degres Centigrades
c   * ustar est le stress du vent
c-----------------------------------------------------------------------

       	implicit none

c Passage des parmaètres communs
        common /s/       s
       	common /c/       c
        common /g/       g
        common /kc/      kc
        common /ustar/   ustar
        common /ustar_c/   ustar_c
        common /b0/      b0
        common /B_/      B_
        common /k0/      k0
        common /kd/      kd
        common /kmax/    kmax
       	common /phi_1/   phi_1
        common /krho1/   krho1
       	common /theta_1/ theta_1
        common /cos1/    cos1
        common /sin1/    sin1
        common /h/       h
       	common /Uz/      Uz
        common /alt/     alt
        common /PI2/     PI2
       	common /epsi/    epsi
        common /epsi_/   epsi_
        common /C1/      C1
        common /D1/      D1
        common /D2/      D2
        common /Rvv0/    Rvv0
        common /Rhh0/    Rhh0
        common /flagI/   flagI
        common /U10/     U10
        common /cspec/   cspec
        common /Omega/   Omega
c Définition des fonction externes        
       	external funcUstar_y
       	external funcUstar_c
        external fkinf
        external fksup
     	external dcosd
     	external dsind
     	external dtand
     	external dacosd
     	external dasind
     	external epsilon_KS
     	external epsilon_El
     	external func2D1
c Déclaration des variables 
        character*80 fout1
        character*80 fout2
        character*1 cVar
        character*1  cepsi
        character*1  csigneeps
        character*1  cSpectre
        character*1  cType
        character*2  cCouvEcume
        character*1  cEmisEcume
        character*1  cCD
     	character*1  cAtmo
     	character*1  cSwell
        character*80 sautligne
        character*40 Modepsi
        character*40 ModVarPente
        character*40 ModSpectre
        character*40 TypeMod
        character*40 ModCouvEcume
        character*40 ModEmisEcume
        character*40 ModCD
        character*10 date
        character*9 time
	
        double precision kmin
        double precision kmax
        double precision phimin
        double precision phimax
        double precision dcosd
        double precision dsind
        double precision dtand
        double precision dacosd
        double precision dasind
        double precision Wind(1:100)
        double precision Windmin
        double precision Windmax
        double precision phi1
        double precision phi2
       	double precision theta(1:100)
        double precision thetamin
        double precision theta_max
        double precision SSS(1:100)
        double precision SSSmin
        double precision SSSmax
        double precision SST(1:100)
        double precision SSTmin
        double precision SSTmax
        double precision nu
        double precision Tw
        double precision kdtab(1:100)
        double precision kd
        double precision kdmin
        double precision kdmax
        double precision k0
        double precision k02
        double precision lambda
        double precision lambdad
        double precision phi
        double precision theta_1
        double precision theta_1min
        double precision thetaAtmo
        double precision cos1
        double precision sin1
       	double precision cophi1
        double precision siphi1
        double precision h
       	double precision phi_1
       	double precision freq
     	double precision Ir (0:2,1:4)
     	double precision dIr(1:4)
     	double precision Irh (0:2,1:4,0:180)
     	double precision dIrh(1:4)
        double precision PI2
     	double precision B_
        double precision g
        double precision s
        double precision R
        double precision D
        double precision kc
        double precision b0
        double precision c
        double precision k
        double precision kinf
        double precision ksup
        double precision kmid
        double precision ustar
        double precision ustar_y
        double precision ustar_c
        double precision Uz
        double precision alt
        double precision U19
        double precision U12
        double precision U10
        double precision krho
        double precision krho1
        double precision krhomin
        double precision krhomax
        double precision binf
        double precision THETAmax
        double precision TPHI_1 (0:2)
        double precision TPHI (1:4)
        double precision Su
        double precision Sc
        double precision sigma
        double precision Stab(1:100)
        double precision Stabmin
        double precision Stabmax
        double precision TbAd(0:90)
        double precision TbAu(0:90)
        double precision tau(0:90)
        double precision tau_
        double precision AtmoD
        double precision AtmoU
        double precision Param(1:10)
        double precision Fr
        double precision Fr_
        double precision C2
        double precision C7
        double precision Norm2
        double precision cothet
        double precision sithet
        double precision cophi
        double precision siphi
        double precision Sx
        double precision Sy
        double precision Sx_sup
        double precision Sy_sup
        double precision Sx_inf
        double precision Sy_inf
        double precision Sx_
        double precision del1
        double precision del2
        double precision limit
        double precision projec
        double precision thetal
        double precision phil
        double precision cosl
        double precision sinl
        double precision epsi_sf(1:2)
        double precision Is(1:4) 
        double precision Is_e(1:4) 
        double precision IsA(1:4)  ! Is avec atmosphere 
        double precision Is_eA(1:4) ! Is_e avec atmosphere 
        double precision Isl(1:4) 
        double precision Isl_e(1:2) 
        double precision IslA(1:4)  ! Isl avec atmosphere
        double precision Isl_eA(1:2) ! Isl_e avec atmosphere 
        double precision int1(1:4,0:40)
        double precision int1e(1:4,0:40)
        double precision int1A(1:4,0:40) ! int1 avec atmosphere
        double precision int1eA(1:4,0:40) ! int1e avec atmosphere
        double precision cosalpha
        double precision sinalpha 
        double precision P_Sx_Sy
        double precision Diff(1:4)
        double precision Poids
        double precision DTHETAtab
        double precision lambdamin
        double precision lambdamax
        double precision niVar
        double precision Tv0
        double precision Th0
        double precision Tv1
        double precision Th1
        double precision U1 
        double precision V1
        double precision Tv2
        double precision Th2
        double precision U2 
        double precision V2
        double precision Tvn
        double precision Thn
     	double precision RvEff
     	double precision RhEff
        double precision xVar
        double precision VarSwell(1:2) ! variances des pentes de la houle en x et y
        double precision hSwell ! rms des hauteurs de la houle
        double precision sigSwell(1:2) ! largeur à mi-puissance de la densite spectrale de la houle en x et y
        double precision KMaxSwell(1:2) ! pics de la densite spectrale de la houle en x et y
        double precision Omega

        
       	double complex epsi_KS
       	double complex epsi_El
       	double complex epsi
       	double complex epsi_
        double complex C1
     	double complex D1
        double complex D2
        double complex Rvv0
        double complex Rhh0
       	
        integer i 
        integer iWind
        integer nWind
        integer iphi
        integer j 
        integer jj 
        integer l 
        integer ll 
        integer nParam
        integer nphi 
        integer ntheta 
        integer nSST 
        integer nSSS
        integer Nstab
        integer nkd
        integer iSST 
       	integer iSSS 
        integer iStab
        integer itheta
        integer ikd
        integer flagB
        integer indice 
        integer flagI 
        integer compt 
        integer nTHETA_1
        integer nPHI_1
        integer N1
        integer N2
        integer ind
        integer cspec
        integer fDiff
        integer fGeo
        integer fVar
        integer fCouvEcume
        integer fEmisEcume
        integer fCD
        integer fSwell
        integer nsorties
        integer ufile
        integer sortiestab(1:2)
        integer NSwell (1:2) ! nb de point de tabulation de la houle
	
c--Paramètres----------------------------------------------------------
c       B_ : Amplitude du spectre de Durden Vesecky (0.004 à l'origine)
c       g  : Accélération de la gravitation
c       s  : Paramètre de la variation angulaire du spectre
c       N1 : le nombre de pts d'intégration sur les pentes en Sx
c       N2 : le nombre de pts d'intégration sur les pentes en Sy
c       PI2 = 2*Pi
c       THETAmax : Valeur max des theta pour la tabulation 
c                  de la diffusion
c       TPHI_1 : contient les 3 angle phi pour déterminer le fondamental
c               et la 2eme harmonique des 2 premiers parametres de 
c               Stokes (à 0 et 90°) et la 2eme harmonique des 3eme et
c               4eme parametres de Stokes (à 45°)
c       TPHI   : contient les 4 angle phi pour déterminer les 3 
c               coefficients de la decompostion en serie de fourier 
c               paire pour les 2 premiers parametres de Stokes (0, 90
c               et 180 °) et en serie impaire pour les 2 derniers 
c               parametres de Stokes(45 et 90°)

       	data g/9.81D0/
       	data s/1.5D-04/
        PI2=2.0D0*acos(-1.0D0)
        PI2=2.0D0*3.141592653589793238462643383279D0
        TPHI_1(0) = 0.0D0
        TPHI_1(1) = 90.0D0
        TPHI_1(2) = 45.0D0
        TPHI(1)   = 0.0D0 
        TPHI(2)   = 45.0D0 
        TPHI(3)   = 90.0D0 
        TPHI(4)   = 180.0D0 
        do i = 0, 2
                do j = 1, 180
                        do indice = 1, 4
                                Irh(i, j, indice) = 0.0D0
                        enddo
                enddo
        enddo

     	open (unit=30,file='Tb.p',status='unknown')
        open (unit=40, file='Diffusion.dat',status='unknown')
        open (unit=50, file='rapport.dat',status='unknown')
        
c----------------------------------------------------------------------

c------------------ ENTREE DES VARIABLES  -----------------------------
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,*) lambda
		write(*,*) 'ok' 
       	freq = 3.D08/lambda
       	nu = 3.D-01/lambda
        call readvec(kdtab, 30, nkd)
		write(*,*) 'ok' 
        call extr_val(kdtab, nkd, kdmin, kdmax)
		write(*,*) 'ok' 
        call readvec(theta, 30, ntheta)
		write(*,*) 'ok' 
        call extr_val(theta, ntheta, thetamin, theta_max)
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        call readvec(Wind, 30, nWind)
		write(*,*) 'ok' 
        call extr_val(Wind, nWind, Windmin, Windmax)
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
       	read (30,*) alt
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) cCD               ! cCD mod. coefficient trainee
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        call readvec(SST, 30, nSST)
		write(*,*) 'ok' 
        call extr_val(SST, nSST, SSTmin, SSTmax)
		write(*,*) 'ok' 
        call readvec(SSS, 30, nSSS)
		write(*,*) 'ok' 
        call extr_val(SSS, nSSS, SSSmin, SSSmax)
		write(*,*) 'ok' 
        call readvec(Stab, 30, nStab)
		write(*,*) 'ok' 
        call extr_val(Stab, nStab, Stabmin, Stabmax)
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') fout1
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') fout2
		write(*,*) 'ok' 
        nsorties = 0
        if (fout1.eq.'aucun') then
                print *, '  Pas de sortie fichier avec écume'
                if (fout2.eq.'aucun') then
                        print *, '  Pas de sortie fichier sans écume'
                        print *, '!!!!!!! ATTENTION  !!!!!!'
                        print *, 'Pas de sortie fichier prévue'
                else
        	   open (unit=20,file=fout2,status='unknown')
		   open (unit=25,file="Atm_"//fout2,status='unknown')
                        nsorties = 1
                        sortiestab(1) = 20
                endif
        else
                open (unit=10,file=fout1,status='unknown')
                open (unit=15,file="Atm_"//fout1,status='unknown')
                if (fout2.eq.'aucun') then
                        print *, '  Pas de sortie fichier sans écume'
                        nsorties = 1
                        sortiestab(1) = 10
                else
                   open (unit=20,file=fout2,status='unknown')
                   open (unit=25,file="Atm_"//fout2,status='unknown')
                        nsorties = 2
                        sortiestab(1) = 10
                        sortiestab(2) = 20
                endif
        endif
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,'(a)') cVar
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne     ! B_ Amplitude du spectre DV
		write(*,*) 'ok' 
        read (30,*) B_
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) cSpectre          ! cSpectre mod. de spectre
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) Omega          ! Inverse de l'age des vagues
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne     ! cCouvEcume mod. couverture ecume
		write(*,*) 'ok' 
        read (30,*) cCouvEcume
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) cEmisEcume        ! cEmisEcume mod. emissiv. ecume  
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) cepsi             ! cepsi mod. permitivite
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) csigneeps         ! csigneeps Signe Im(permitivite)
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) cSwell
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) NSwell(1)
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) NSwell(2)
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) hSwell
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) sigSwell(1)
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) sigSwell(2)
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) KMaxSwell(1)
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) KMaxSwell(2)
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) THETAmax, DTHETAtab
        if (DTHETAtab.lt.1.0D0) then
                write (*,*)
                print *, '!!!!!! ERREUR  !!!!!!!!!!!'
                print *, 'DTHETAtab doit etre au moins égal à 1'
                write (*,*)
                stop
        endif
     	read (30,'(a)') sautligne     ! lambdamin
		write(*,*) 'ok' 
        read (30,*) lambdamin
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne     ! lambdamax
		write(*,*) 'ok' 
        read (30,*) lambdamax
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne     ! niVar
		write(*,*) 'ok' 
        read (30,*) niVar
		write(*,*) 'ok' 
        read (30,'(a)') sautligne     ! cType Type de Modele 
		write(*,*) 'ok' 
        read (30,*) cType
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) N1                ! Nb integration sur Sx
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) N2                ! Nb integration sur Sy
		write(*,*) 'ok' 
     	read (30,'(a)') sautligne
		write(*,*) 'ok' 
        read (30,*) xVar              ! Limite integration sur pentes (xVar x Variance)
		write(*,*) 'ok' 
        
        call date_and_time(date, time)
        date = date(7:8)//'/'//date(5:6)//'/'//date(1:4)
        time = time(1:2)//':'//time(3:4)//':'//time(5:6)


c----------- VERIFICATION DES PARAMETRES D'ENTRÉE ----------------------

        if (int(THETAmax/DTHETAtab).ne.(THETAmax/DTHETAtab)) then
                write (*,*)
                print *, '!!!!!! ERREUR  !!!!!!!!!!!'
                print *, 'THETAmax n''est pas multiple de DTHETAtab'
                write (*,*) 'THETAmax = ', THETAmax
                write (*,*) 'DTHETAtab = ', DTHETAtab
                write (*,*)
                stop
        endif
        if ((theta_max.gt.THETAmax).and.(nTheta.ne.1)) then
                write (*,*)
                print *, '!!!!!!!! ERREUR  !!!!!!!!!!!!!!!!'
                print *, 'La valeur de theta est limitée à ', THETAmax
     &          ,' par THETAmax (valeur maximale tabulée)'
                write (*,*)
                stop
        endif
        if ((cSpectre.eq.'D').or.(cSpectre.eq.'d')) then
                cspec = 1
                ModSpectre = 'Durden & Vesecky (85)'
        elseif ((cSpectre.eq.'E').or.(cSpectre.eq.'e')) then
                cspec = 2
                ModSpectre = 'Elfouhaily (97)'
        else
                write (*,*)
                print *, '!!!!!!!! ERREUR  !!!!!!!!!!!!!!!!'
                print *, 'Le choix du spectre est ', cSpectre
                print *, 'Choix possibles : E pour Elfouhaily (97)'
                print *, '                : D pour Durden & Vesecky(85)'
                write (*,*)
                stop
        endif
        if (cType.eq.'2') then
                TypeMod = '2 Echelles'
                fGeo  = 1
                fDiff = 1 
        elseif ((cType.eq.'p').or.(cType.eq.'P')) then
                TypeMod = 'Petites Echelles (Diffusion)'
                fGeo  = 0
                fDiff = 1
        elseif ((cType.eq.'g').or.(cType.eq.'G')) then
                TypeMod = 'Grandes Echelles (Optique Geometrique)'
                fGeo  = 1
                fDiff = 0
        else
                write (*,*)
                print *, '!!!!!!!! ERREUR  !!!!!!!!!!!!!!!!'
                print *, 'Le choix du type de modele est ', cType
                print *, 'Choix possibles : 2 pour 2 echelles'
                print *, '                : p pour petites echelles'
                print *, '                : g pour grandes echelles'
                write (*,*)
                stop
        endif
        if ((cepsi.eq.'K').or.(cepsi.eq.'k')) then
                Modepsi = 'Klein & Swift (77)'
        elseif ((cepsi.eq.'e').or.(cepsi.eq.'E')) then
                Modepsi = 'Ellison (98)'
        else
                write(*,*)
                print *,'Mauvais choix pour constante dielectrique:'
                print *,'     K : Klein & Swift (77)'
                print *,'     E : Ellison (98)'
                print *,'Votre choix etait : ',cepsi
                write(*,*)
                stop
        endif
        if ((cVar.eq.'C').or.(cVar.eq.'c')) then
                ModVarPente = 'Cox & Munk (55)'
                fVar = 0
        elseif ((cVar.eq.'S').or.(cVar.eq.'s')) then
                ModVarPente = 'Spectre'
                fVar = 1
        else
                write(*,*)
                print *,'Mauvais choix pour le calcul des variances:'
                print *,'     C : Cox & Munk (55)'
                print *,'     S : calcul avec le spectre'
                print *,'Votre choix etait : ',cVar
                write(*,*)
                stop
        endif
        if ((cCouvEcume.eq.'A').or.(cCouvEcume.eq.'a')) then
                ModCouvEcume = 'Pas d''ecume'
                fCouvEcume = 0
        elseif ((cCouvEcume.eq.'M1').or.(cCouvEcume.eq.'m1')) then
                ModCouvEcume = 'Monahan (86, eq.5)'
                fCouvEcume = 1
        elseif ((cCouvEcume.eq.'M2').or.(cCouvEcume.eq.'m2')) then
                ModCouvEcume = 'Monahan moindres carres (86, eq. 3a)'
                fCouvEcume = 2
        else
                write(*,*)
                print *,'Mauvais choix pour la couverture d''ecume:'
                print *,'     A  : pas d''ecume'
                print *,'     M1 : Monahan (86, eq.5)'
                print *,'     M1 : Monahan moindre carres (86, eq. 3a)'
                print *,'Votre choix etait : ',cCouvEcume
                write(*,*)
                stop
        endif
        if ((cEmisEcume.eq.'S').or.(cEmisEcume.eq.'s')) then
                ModEmisEcume = 'Stogryn (72)'
                fEmisEcume = 1
        else
                write(*,*)
                print *,'Mauvais choix pour l''emissivite de l''ecume:'
                print *,'     S : Stogryn(72)'
                print *,'Votre choix etait : ',cEmisEcume
                write(*,*)
                stop
        endif
        if ((cCD.eq.'y').or.(cCD.eq.'Y')) then
                ModCD = 'Cardone (69)'
                fCD = 1
        elseif ((cCD.eq.'c').or.(cCD.eq.'C')) then
                ModCD = 'Charnock (55)'
                fCD = 2
        else
                write(*,*)
                print *,'Mauvais choix pour le coefficient de trainee:'
                print *,'     y : Cardone (69)'
                print *,'     c : Charnock (55)' 
                print *,'Votre choix etait : ',cCD
                write(*,*)
                stop
        endif
        if ((cSwell.eq.'O').or.(cCD.eq.'o')) then
                fSwell = 1
        elseif ((cSwell.eq.'N').or.(cSwell.eq.'n')) then
                fSwell = 0
        else
                write(*,*)
                print *,'Mauvais choix pour la presence de houle:'
                print *,'     O : oui'
                print *,'     N : non' 
                print *,'Votre choix etait : ',cSwell
                write(*,*)
                stop
        endif

c---------------- AFFICHAGE DES VARIABLES DANS LE FICHIER --------------
        do i = 1, nsorties
        ufile = sortiestab(i)
        write (ufile,*) 'date = '//date//' ; heure = '//time(1:8)//' ;'
       	k0 = 2.D0*3.141592653D0/lambda 
       	write (ufile,*) 'Lambda   = ',lambda,' ;'
       	if (kdmin.eq.0.)  then
                write (*,*) ' !!!! Erreur !!!!'
                write (*,*) 'kdmin = 0 !!!!!'
                stop
        endif
        write (ufile,*) 'Lambda_d min = ', PI2/kdmax,
     &               ' ; Lambda_d max = ', PI2/kdmin,
     &               ' ; Nlambda_d = ', nkd,' ;'


        write (ufile,*) 'Vent min.        = ', Windmin,
     &               ' ; Vent max.        = ', Windmax, 
     &               ' ; Nvent         = ', nWind, ' ;'
        write (ufile,*) 'Incidence min.   = ', thetamin, 
     &               ' ; Incidence max.   =' , theta_max,
     &               ' ; Nincidences   = ', ntheta, ' ;'
        write (ufile,*) 'Temperature min. = ', SSTmin, 
     &               ' ; Temperature max. =' ,SSTmax, 
     &               ' ; Ntemperatures = ', nSST, ' ;'
        write (ufile,*) 'Stabilite min.   = ', Stabmin,
     &               ' ; Stabilite max.   =' ,Stabmax,
     &               ' ; Nstabilites   = ', nStab, ' ;'
        write (ufile,*) 'Salinite min.    = ', SSSmin,
     &               ' ; Salinite max.    =' ,SSSmax,
     &               ' ; Nsalinites    = ', nSSS, ' ;'
        write (ufile,*) 'Modele constante dielectrique = ', Modepsi,' ;'
        write (ufile,*) 'Calcul de la Variance des pentes = ',
     &               ModVarPente,' ;'
        write (ufile,*) 'Type de modele                = ', TypeMod,' ;'
        write (ufile,*) 'Modele de spectre des vagues  = ',
     &                 ModSpectre,' ;'
        write (ufile,*) 'Nombre de points sur Gaussienne sur axe x = ',
     &                N1, ' ;' 
        write (ufile,*) 'Nombre de points sur Gaussienne sur axe y = ',
     &                N2, ' ;' 
        if (ufile.eq.10) then
                write (ufile,*) 'Modele de couverture d''ecume = ', 
     &                  ModCouvEcume, ' ;'
                write (ufile,*) 'Modele d''emissivite de l''ecume = ', 
     &                  ModEmisEcume, ' ;'
        else if (ufile.eq.20) then
                write (ufile,*) 'Modele de couverture d''ecume = ',
     &                 'Pas d''ecume;'
                write (ufile,*) 'Modele d''emissivite de l''ecume = ',
     &                 'Pas d''ecume;'
        else
                print *, 'Erreur, ne reconnait pas le fichier de sortie'
                stop
        endif
        write (ufile,*) 'Modele de coefficient de trainee = ',ModCD,' ;'
        write (ufile,*) '! freq     SST   SSS     U10   ustar theta phi'
     &  //' TvN    ThN     Tv0      Th0    Tv1    Th1    U1      V1    '
     &  //' Tv2   Th2    U2      V2    lam_d    Stab Re(epsi)I(epsi)'
     &  //' Var_u    Var_c      Foam;'
        write (ufile,*) '! constante dielectrique e = Re(e) + i Im(e) $'

        enddo
        
c---------------------------------------------------------------------
c------------------ PROGRAMME PRINCIPAL ------------------------------

c  Calcul des Variance de pentes de la houle

        write(50,*) 'Variances des pentes de la houle'
        if (fSwell.eq.1) then
                call infSwell(NSwell, hSwell, sigSwell , 
     &          KMaxSwell, VarSwell)
        else
                VarSwell(1) = 0.0D0
                VarSwell(2) = 0.0D0
        endif
        write(50,*) 'Variances des pentes de la houle terminée'
        write (50,*) 'SuSwell = ', sqrt(VarSwell(1))
        write (50,*) 'ScSwell = ', sqrt(VarSwell(2))

c Début Boucle sur SSS

        do iSSS=1, nSSS
       	   
c Début Boucle sur SST

       	do iSST = 1, nSST
       	   
       	   Tw = 2.7315D2 + SST(iSST)
           
c----------------------------------------------------------------------
c                       CONSTANTE DIELECTRIQUE : epsi
c
c epsi : constante dielectrique complexe relative à epsilon 0 
c        (ie permitivité du vide)
c       epsi_KS pour Klein & Swift 
c       epsi_El pour Ellison
c SST  : Temperature de surface de l'ocean en °K
c SSS  : Salinité de surface de l'ocean en ppm
c freq : Frequence electromagnetique du radiometre en Hz
        
        write (50,*) 'Calcul de constante dielectrique'
        call epsilon_KS (SST(iSST), SSS(iSSS), epsi_KS, freq)
        call epsilon_El (SST(iSST), SSS(iSSS), epsi_El, freq)
        if ((cepsi.eq.'K').or.(cepsi.eq.'k')) then
                epsi = epsi_KS
                Modepsi = 'Klein & Swift (77)'
        elseif ((cepsi.eq.'e').or.(cepsi.eq.'E')) then
                epsi = epsi_El
                Modepsi = 'Ellison (98)'
        else
                write(*,*)
                print *,'Mauvais choix pour constante dielectrique:'
                print *,'     K : Klein & Swift'
                print *,'     E : Ellison'
                print *,'Votre choix etait : ',cepsi
                write(*,*)
                stop
        endif
        if ((csigneeps.eq.'o').or.(csigneeps.eq.'O')) then
                continue
        elseif ((csigneeps.eq.'n').or.(csigneeps.eq.'N')) then
                epsi = conjg(epsi)
        else
                write(*,*)
                print *, 'Mauvais choix pour le signe de la '
                print *, 'constante dielectrique comme positif:'
                print *, '     O : oui'
                print *, '     N : non'
                print *, 'Votre choix etait : ',csigneeps
                write(*,*)
                stop
        endif

        !epsi = (1.0D0, -1.0D16)
        write (50,*) 'Calcul de constante dielectrique terminée'
        write (50,*) '  Klein & Swift = ', epsi_KS
        write (50,*) '  Ellison       = ', epsi_El
        write (50,*) '      =>       ',epsi

c Tabulation des Tb Atmosphérique en fonction de l'angle d'incidence 
c theta (theta de 0 à 90° par pas de 1°)
        
        nParam = 8
        Param(1) = 6370.0D0 ! rayon terrestre (km)
        Param(2) = 1013.0D0 ! pression atmosphérique au sol (mb)
        Param(3) = SST(iSST) + 273.15D0  ! Température au sol (K)
        Param(4) = 36.0D0  ! Humidité relative (%)
        Param(5) = freq  ! fréquence (Hz)
        Param(6) = 20.0D0  ! altitude maximum de l'atmosphère (km)
        Param(7) = 400  ! nb de couches d'atmosphere
        Param(8) = 6.5D0  ! module gradient de temperature (K/km)
        write (50,*) 'Tabulation de températures de brillance de l'''
     &   ,'atmosphère'
c        call TbAtmo(Param, nParam, TbAd, TbAu, tau) 
        ! TbAd Tb atmosphere vers le bas (downward)
        ! TbAu Tb atmosphere vers le haut (upward)
c	do i = 0, 90
c        write (*,*) i*1.0D0, TbAu(i), TbAd(i), tau(i)
c	enddo
        do i = 0, 90
                tau(i) = 0.0D0
                TbAd(i) = 0.0D0
                TbAu(i) = 0.0D0
        enddo
        write (50,*) 'Tabulation de températures de brillance de l'''
     &        ,'atmosphère terminée'


c----------------------------------------------------------------------
c Début Boucle sur Vent

       	do iWind= 1, nWind
                Uz = Wind(iWind)
 
c----------------------------------------------------------------------
c                            STRESS : ustar
c
c ustar est borne entre arg#1 et arg#2 et on impose une precision 
c de arg#3 sur le calcul de la racine. 
c
c Uz : Module du vent en m/s
c z  : Altitude pour le module du vent en metres
c ustar_c calculé avec Z0 Charnock 55
c ustar_y calculé avec Z0 Yueh 97

         write (50,*) 'Calcul du stress pour U(',alt,') = ', Uz
      	if (Uz.ne.0.) then
                call root (1.D-05, 1.5D1, 1.D-06,funcUstar_y, ustar_y)
                call root (1.D-05, 1.5D1, 1.D-06,funcUstar_c, ustar_c)
       	else
                ustar = 0.D0
       	endif
        if (fCD.eq.1) then
                ustar = ustar_y
        elseif (fCD.eq.2) then
                ustar = ustar_c
        else
                print *, 'Je n''ai pas compris quel parametrisation'
                print *, 'vous desirez pour le coefficient de trainee'
                print *, 'Votre choix etait : ', cCD
                print *, 'Choix Possibles : y pour Cardone (69)'
                print *, '                  c pour Charnock (55)'
                write (*,*)
                stop
        endif
        if (ustar.eq.1.0D04) print *, 'Changer les bornes passée à root'
     &          ,'pour le calcul de u*'
         write (50,*) 'Calcul du stress pour U(',alt,') = ', Uz
     &   , ' terminée'
         write (50,*) '   ustar = ', ustar

c       	write (10,*)
c        write (10,*) '________________________________________'
c        write (10,*) '       PARAMETRES VENT ET STRESS        '
c       	write (10,*)
c	write (10,*) " ---> U*(Yueh)     = ",  (ustar)
c	write (10,*) " ---> U*(Charnock) = ",  (ustar_c)

c----------------------------------------------------------------------
c               CONVERSION DU VENT A 10, 12.5 ET 19.5 METRES
c
c  Uz    : Module du vent en m/s
c  z     : Altitude pour le module du vent en metres
c  ustar : Stress du vent en m
c  U19   : Module du vent en m/s pour une hauteur de 19.5m
c  U12   : Module du vent en m/s pour une hauteur de 12.5m
c  U10   : Module du vent en m/s pour une hauteur de 10.0m

         write (50,*) 'Conversion du vent'
       	if (Uz.ne.0.) then
                call vent (alt,1.95D01,Uz,U19,ustar)
                call vent (alt,1.25D01,Uz,U12,ustar)
                call vent (alt,1.0D01 ,Uz,U10,ustar)
       	else 
                U19 = 0.
                U12 = 0.
       	endif
         
         write (50,*) 'Conversion du vent terminée'
         write (50,*) '   U10 = ', U10
         write (50,*) '   U12.5 = ', U12
         write (50,*) '   U19.5 = ', U19
c Début boucle sur kd

        do ikd = 1, nkd
        
                kd = kdtab(ikd)
                lambdad = PI2/kd
        
c----------------------------------------------------------------------
c       PARAMETRES DU MODELE DE SPECTRE DE DURDEN & VESECKY
c

       	kc   = g/U19/U19
        kmid = sqrt(1.48D0/3.0D0)*kc
       	b0   = B_*dexp(1.85D-01*kc*kc)
        write (50,*) 'Calcul de la borne inférieure du spectre'
        call root (1.0D-03, 2.0D0, 1.D-06,fkinf, kinf)
!        kinf = 1.0D-03 
        write (50,*) 'Calcul de la borne inférieure du spectre terminé'
        write (50,*) '   kinf = ', kinf
        kd = max(kinf, kd)
        write (50,*) ' => kd = ', kd
        write (50,*) 'Calcul de ''c'' du spectre de Durden & Vesecky'
        call c_ (c, U12)
        write (50,*) 'Calcul de ''c'' du spectre de Durden & Vesecky'
     &              ,' terminé'
        write (50,*) '   c = ',c
        kmax = sqrt(1.48D0/3.0D0)*kc ! k pour max du spectre
!        print*, kd
c        if (kinf.eq.1.0D04) print *, 'Changer les bornes passée à root '
c     &          ,'pour le calcul de kinf'
c        call root (2.0D0, 1.0D+05, 1.D-06,fksup, ksup)
c        if (ksup.eq.1.0D04) print *, 'Changer les bornes passée à root '
c     &          ,'pour le calcul de ksup'
c        print *, 'kinf = ', kinf, '  S_DV(kinf)/kinf', S_DV(kinf)/kinf
c        print *, 'ksup = ', ksup, '  S_DV(ksup)', S_DV(ksup)

c----------------------------------------------------------------------
c               VARIANCES DES PENTES DU SPECTRE
c
c  Su : Variances des pentes en direction "upwind"
c  Sc : Variances des pentes en direction "downwind"
c  kd : Nombre d'onde de coupure en rad/m
c  cspec : modele de spectre

        write (50,*) 'Calcul des variances des pentes'
        if (fVar.eq.1) then
                call Su_ (Su, 1.0D-03, kd, cspec)
                call Sc_ (Sc, 1.0D-03, kd, cspec)
        else if (fVar.eq.0) then
                Su = sqrt(3.16D-03*U12)
                Sc = sqrt(0.003D0 + 1.92D-03*U12)
        else
                write (*,*)
                print *, '!!!!!!!!!!  ERREUR !!!!!!!!!!!!!! '
                print *, 'Pas de Choix valide pour le modele'
                print *, 'des variances de pentes !!!!'
                print *, 'votre choix est ', fVar
                print *, 'Il doit etre C pour Cox & Munk'
                print *, '             S pour calcule par le spectre'
                write (*,*)
                stop
        endif
        Su = sqrt(Su*Su + VarSwell(1))
        Sc = sqrt(Sc*Sc + VarSwell(2))
        write (50,*) 'Calcul des variances des pentes terminé'
        write (50,*) '   Su = ', Su
        write (50,*) '   Sc = ', Sc

c----------------------------------------------------------------------
c               VARIANCES DES HAUTEURS DU SPECTRE DES PETITE VAGUES
c
c  sigma : Variances des hauteurs des petites échelles
c  kd : Nombre d'onde de coupure en rad/m
c  cspec : modele de spectre

        write (50,*) 'Calcul des variances des hauteurs'
        call sigma_ (sigma, kd, cspec)
        write (50,*) 'Calcul des variances des hauteurs terminé'
        write (50,*) '   sigma =', sigma

c----------------------------------------------------------------------
c                       FRACTION DE COUVERTURE D'ECUME
c

c  Couverture d'ecume avec dependance en Stab(stabilite atmos.)
        write (50,*) 'Calcul de la fraction d''écume'
        if (fCouvEcume.eq.1) then
                call foam (U10,Stab(1),Fr)
        elseif (fCouvEcume.eq.2) then
                Fr = 2.95D-6*U10**(3.52)
        elseif (fCouvEcume.eq.0) then
                Fr  = 0.0D0
        else
                        write (*,*)
                        print *, '!!!!!!!!!!  ERREUR !!!!!!!!!!!!!! '
                        print *, 'Je n''ai pas compris le modele'
                        print *, 'd''ecume que vous voulez!'
                        write (*,*)
                        stop
        endif
        Fr_ = 1.0D0 - Fr
        write (50,*) 'Calcul de la fraction d''écume terminé'
        write (50,*) '   Fraction d''écume = ', Fr

c  Couverture d'ecume par methode des moindres carres

c----------------------------------------------------------------------
c                              TABULATION DIFFUSION

        if (fDiff.eq.1) then
c   Calcul les constantes
        k02     = k0*k0
        epsi_   = (epsi - 1.0D0)
        
c    Tabulation de la diffusion en fonction de theta_local et phi_local
c       
c       La coefficient de reflexion dans le plan local Ir = Irc + Iri
c       On décompose Ir en fonction de l'azimuth local phi_1 comme suit:
c               Ir = Irh(0) + Irh(1)*cos (phi_1) + Irh(2)*cos(2*phi_1)
c       La première harmonique étant nulle, on a détermine Irh(0) 
c       et Irh(2) en phi_1 = 0° et phi_2 = 180°

        print *, 'Tabulation de la diffusion engagée ...'
c Echantillonnage en theta de 0 à THETAmax par pas de DTHETAtab

c Début BOUCLE THETA_1>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        do nTHETA_1 =  0, THETAmax/DTHETAtab
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        theta_1 = nTHETA_1*DTHETAtab
        cos1    = dcosd(theta_1)
        sin1    = dsind(theta_1)
        krho1   = k0*sin1
c   Calcul les constantes
        C1      = sqrt (epsi - sin1*sin1)
        D1      = cos1 + C1
        D2      = epsi*cos1 + C1
c   Calcul les coefficient de réflexion de Fresnel
        Rvv0    = (epsi*cos1-C1)/(epsi*cos1+C1)
        Rhh0    = (cos1 - C1)/(cos1+C1)

c Début BOUCLE PHI_1>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        do nPHI_1 = 0,1
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        phi_1   = TPHI_1(nPHI_1)
        cophi1  = dcosd(phi_1)
        siphi1  = dsind(phi_1)

c      initialisation
        do indice = 1,2
                Ir(nPHI_1, indice)  = 0.0D0
                dIr(indice) = 0.0D0
        enddo

        krhomin = max(kinf, (kd - krho1)*1.0000001D0)
c        krhomax = min(1.1D0*krhomin, k0*0.999999D0)
        krhomax = k0*0.999999D0
        
c  INTEGRALE SUR L'HEMISPHERE SUPERIEUR DE LA VAGUE
        flagI = 1
        if (krhomin.lt.krhomax) then
        flagB = 1
        dowhile(flagB.eq.1)
                call qromb2D1(func2D1, krhomin, krhomax, dIr, 1.0D-03)
                do indice = 1, 2
                   Ir(nPHI_1,indice) = Ir(nPHI_1,indice)+dIr(indice)*k02
                enddo
                if (abs(krhomax-(k0*0.999999D0)).lt.1.0D-06) flagB = 0
                krhomin   = krhomax
c               krhomax = min(1.1D0*krhomin, k0*0.999999D0)
                krhomax = k0*0.999999D0

        enddo
        endif
        
c  INTEGRALE DE LA DIFFUSION SUR LES FAIBLE LONGUEURES D'ONDE

        flagI = 0
        do i = 1,30
                krhomin = k0*1.000001D0*10.0D0**((i-1.0D0)/10.0D0)
                krhomax = k0*10.0D0**(i/10.0D0)
               call qromb2D1(func2D1, krhomin, krhomax, dIr, 1.0D-03)
                do indice = 1, 2
                   Ir(nPHI_1,indice) = Ir(nPHI_1,indice)+dIr(indice)*k02
                enddo
        enddo


c  Fin Boucle Phi_1<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        enddo
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c  Harmoniques 0 et 2 des Stokes 1 et 2
c
c       Irh(i,indice,nTHETA_1) est l'harmonique "i", 
c                              pour le Stokes "indice",
c                              en theta = nTHETA_1*DTHETAtab
c       Ir(0,indice) est la Tb pour le Stokes "indice" en phi=TPHI_1(0)
c       Ir(1,indice) est la Tb pour le Stokes "indice" en phi=TPHI_1(1)
        
        do indice = 1, 2
                Irh(0,indice,nTHETA_1)=0.5D0*(Ir(0,indice)+Ir(1,indice))
                Irh(1,indice,nTHETA_1)=0.0D0
                Irh(2,indice,nTHETA_1)=0.5D0*(Ir(0,indice)-Ir(1,indice))
        enddo
        write (40,*) nu, kd, SST(iSST), SSS(iSSS), U10, ustar, theta_1
     &                 , Irh(0,1,nTHETA_1)*Tw
     &                 , Irh(0,2,nTHETA_1)*Tw
     &                 , Irh(2,1,nTHETA_1)*Tw
     &                 , Irh(2,2,nTHETA_1)*Tw
c        write (*,*) theta_1, Ir(0,1)*Tw
c     &                 , Ir(1,1)*Tw

c ---------------  DEBUT CALCUL STOKES 3 ET 4  -----------------------
        
        nPHI_1  = 2
        phi_1   = TPHI_1(nPHI_1)
        cophi1  = dcosd(phi_1)
        siphi1  = dsind(phi_1)

c      initialisation
        do indice = 3,4
                Ir(nPHI_1, indice)  = 0.0D0
                dIr(indice) = 0.0D0
        enddo

        krhomin = max(kinf, (kd - krho1)*1.0000001D0)
        krhomax = min(1.1D0*krhomin, k0*0.999999D0)
        krhomax = k0*0.999999D0
        
c  INTEGRALE SUR L'HEMISPHERE SUPERIEUR DE LA VAGUE
        flagI = 1
        if (krhomin.lt.krhomax) then
        flagB = 1
        dowhile(flagB.eq.1)
                call qromb2D1(func2D1, krhomin, krhomax, dIr, 1.0D-03)
                do indice = 3, 4
                   Ir(nPHI_1,indice)=Ir(nPHI_1,indice)+dIr(indice-2)*k02
                enddo
                if (abs(krhomax-(k0*0.999999D0)).lt.1.0D-06) flagB = 0
                krhomin   = krhomax
c               krhomax = min(1.1D0*krhomin, k0*0.999999D0)
        krhomax = k0*0.999999D0

	enddo
	endif
        
c  INTEGRALE DE LA DIFFUSION SUR LES FAIBLE LONGUEURES D'ONDE

        flagI = 0
        do i = 1,30
                krhomin = k0*1.000001D0*10.0D0**((i-1.0D0)/10.0D0)
                krhomax = k0*10.0D0**(i/10.0D0)
                call qromb2D1(func2D1, krhomin, krhomax, dIr, 1.0D-03)
                do indice = 3, 4
                   Ir(nPHI_1,indice)=Ir(nPHI_1,indice)+dIr(indice-2)*k02
                enddo
        enddo
c -----------   FIN DU CALCUL STOKES 3 ET 4 -------------------------

c  sortie fichier de thetal, fondamental, seconde harmonique pour Vpol
c  et Hpol
        do indice = 3, 4
                Irh(0,indice,nTHETA_1) = 0.0D0
                Irh(1,indice,nTHETA_1) = 0.0D0
                Irh(2,indice,nTHETA_1) = Ir(2,indice)
        enddo
c        write (10,*) theta_1, Irh(2,3,nTHETA_1)*Tw,
c     &                   Irh(2,4,nTHETA_1)*Tw
        write (*,*) theta_1, Irh(0,1,nTHETA_1)*Tw
     &                 , Irh(0,2,nTHETA_1)*Tw
     &                 , Irh(2,1,nTHETA_1)*Tw
     &                 , Irh(2,2,nTHETA_1)*Tw
     &                 , Irh(2,3,nTHETA_1)*Tw
     &                 , Irh(2,4,nTHETA_1)*Tw

c  Fin Boucle Theta_1<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        enddo
        elseif (fDiff.eq.0) then
                write(*,*)
                print *, 'Pas de diffusion de Bragg!'
                write(*,*)
        else
                write(*,*)
                print *, 'Je n''ai pas compris si vous vouliez de la'
                print *, 'diffusion de Bragg'
                write (*,*)
        endif
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
c   Integration : PARTIE GEOMETRIQUE 
        
c Début boucle theta>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	do itheta =1, ntheta
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
     	C2 = dtand(theta(itheta))
     	if (C2.ne.0.D0) then	
		C7 = 1.D0/C2
     	else
		C7 = 1.D+99
     	endif
     	cothet = dcosd(theta(itheta))
     	sithet = dsind(theta(itheta))

c-----------------------------------------------------------------------
c                            COEFFICIENTS DE FRESNEL
c
c  Rvv0    : Coefficient de reflexion en polar V
c  Rhh0    : Coefficient de reflexion en polar H
c  epsi    : Constante diélectrique
c-----------------------------------------------------------------------

        Rvv0 = (epsi*cothet - sqrt(epsi-sithet*sithet))
     &        /(epsi*cothet + sqrt(epsi-sithet*sithet))
        Rhh0 = (cothet - sqrt(epsi-sithet*sithet))
     &        /(cothet + sqrt(epsi-sithet*sithet))
        Rvv0 = abs(Rvv0) 
        Rhh0 = abs(Rhh0) 
        Rvv0 = Rvv0*Rvv0
        Rhh0 = Rhh0*Rhh0

c-----------------------------------------------------------------------
c                       TEMPERATURE DE BRILLANCE SANS VENT
c
c  Tvn      : Temperature de brillance en polar verticale pour une 
c            incidence theta et une surface plane
c  Thn      : idem en polar horizontale 
c-----------------------------------------------------------------------

        Tvn  = Tw*(1.0D0 - Rvv0)
        Thn  = Tw*(1.0D0 - Rhh0)


c Début boucle phi>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     	do iphi = 0, 40
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

c        phi = TPHI(iphi)
        phi = 9.0D0*iphi
     	cophi = dcosd(phi)
     	siphi = dsind(phi)
       
c  Calcul des bornes de l'integrale: en +- 3*variance des pentes
c  avec pour condition Sx' < cotan(theta)

        Norm2 = 0.D0
        Sy_sup =  xVar*Sc
       	Sy_inf = -xVar*Sc
       	del1 = (Sy_sup-Sy_inf)/(N1-1)
       	compt = 0

c  INITIALISATION DU VECTEUR INTEGRATION
        do indice=1,4
                int1(indice,iphi) = 0.0D0
                int1e(indice,iphi) = 0.0D0
                int1A(indice,iphi) = 0.0D0
                int1eA(indice,iphi) = 0.0D0
        enddo
c  INITIALISATION DU VECTEUR DIFFUSION
        do indice=1,4
                Diff(indice) = 0.0D0
        enddo
        

c Début  boucle Sy
       	do i=0, N1-1
       		Sy = Sy_inf + i*del1
                if (cophi.eq.0.D0) then
                        Sx_ = Sy*siphi
                        if (Sx_.le.C7) then
                                Sx_sup =  xVar*Su
                                Sx_inf = -xVar*Su
                        else
                                Sx_sup =  66666.D0
                                Sx_inf =  66666.D0
                        endif
                else
                        limit = (C7-Sy*siphi)/cophi
                        if (cophi.gt.0.D0) then
                                Sx_sup = min(limit, xVar*Su)
c                                if (Sx_sup.eq.limit)
c     &                                  compt = compt + 1
                                Sx_sup = Sx_sup-abs(Sx_sup*1.0D-06)
                                Sx_inf = -xVar*Su
                        else
                                Sx_sup =  xVar*Su
                                Sx_inf = max(limit,-xVar*Su)*1.000001D0
                                Sx_inf = Sx_inf+abs(Sx_inf*1.0D-06)
                        endif
                endif
                del2 = (Sx_sup-Sx_inf)/(N2-1)
                if ((Sx_sup.ne.66666.D0).and.(Sx_inf.ne.66666.D0)) then
c Début boucle Sx
          	do j=0, N2-1
          		Sx = Sx_inf + j*del2
                        Sx_ = Sx*cophi + Sy*siphi
          		projec = 1.0D0-Sx_*C2

c-----------------------------------------------------------------------
c                   DENSITE DE PROBABILITE DES PENTES P_Sx_Sy
c
c  Sx      : Pente dans la direction du vent
c  Sy      : Pente dans la direction orthogonale à celle du vent
c  P_Sx_Sy : Densité de probabilité pour les pentes Sx et Sy
c  Su      : Variance des pentes dans la direction du vent
c  Sc      : Variance des pentes dans la direction orthogonale au vent
c  C7      : cotan(theta), limit de visibilité de la surface des vagues
c-----------------------------------------------------------------------

       			if (Sx_.gt.C7) then
				P_Sx_Sy = 0.0D0
                        else
                                call P (Sx, Sy, P_Sx_Sy, Su, Sc)
       			endif

c-----------------------------------------------------------------------
c      PASSAGE DES COORDONÉES LOCALES EN COORD GLOBALES 
c-----------------------------------------------------------------------

        call Glob2Loc(theta(itheta), phi, Sx, Sy, thetal, phil,
     &                cosAlpha, sinAlpha)
c        print*, thetaL, phiL, cosAlpha, sinAlpha
        sinl = dsind(thetal)
     	cosl = dcosd(thetal)
c        print*, thetal, phil

c-----------------------------------------------------------------------
c       ANGLE D'INCIDENCE DANS LE REPERE GLOBAL DE LA REFLEXION 
c       SPECULAIRE SUR LA VAGUE
c         le satellite est dans la direction theta(itheta) (rep global)
c                          dans la direction [thetal, phil] (rep local 
c                               de la vague)
c         le rayon incident réfléchi spéculairement dans le plan de 
c             la vague vient de [thetal,-phil] (rep local)
c             et il vient de l'angle thetaAtmo (dans le repère gloabal)
c-----------------------------------------------------------------------

c        thetaAtmo = dacosd(cothet + 
c     &                  2.0D0*sinl*siphil*sithetn*siphin*cosbeta)
c
c        AtmoD = TbAd(int(thetaAtmo))*(int(thetaAtmo)-thetaAtmo+1) 
c     &         + TbAd(int(thetaAtmo)+1)*(thetaAtmo-int(thetaAtmo))
        tau_ = tau(int(theta(itheta)))
     &                  *(int(theta(itheta))-theta(itheta)+1)
     &         + tau(int(theta(itheta))+1)
     &                  *(theta(itheta)-int(theta(itheta)))
	AtmoD = 0.0D0
	tau_ = 0.0D0

c-----------------------------------------------------------------------
c               COEFFICIENT DE MODULATION HYDRODYNAMIQUE h
c
c  Sx : Pente de la vague dans la direction du vent
c  Su : Variance de la distribution des vagues dans la direction du vent
c  h  : Coefficient de modulation hydrodynamique
c-----------------------------------------------------------------------

          		if (abs(Sx/Su).le.(1.25D0)) then
       				h = 1.D0- 4.D-01*Sx/Su
          		else
       				h = 1.D0 - 5.D-01*dsign(1.D0,Sx)
          		endif

c-----------------------------------------------------------------------
c                       EMISSIVITE DE L'ECUME epsi_sf
c
c  thetal  : Angle d'incidence dans le plan de la vague ( = local)
c  nu      : Fréquence électromagnétique en GHz
c  epsi_sf : Emissivité de l'écume en polars V (indice 1) et H(indice 2)
c-----------------------------------------------------------------------

       			call esf (thetal, nu, epsi_sf)

c-----------------------------------------------------------------------
c                      COEFFICIENTS DE FRESNEL LOCAUX
c
c  Rvv0    : Coefficient de reflexion en polar V
c  Rhh0    : Coefficient de reflexion en polar H
c  epsi    : Constante diélectrique
c  cosl    : Cosinus de l'incidence dans le plan de la vague ( = local)
c  sinl    : Sinus de l'incidence dans le plan de la vague ( = local)
c-----------------------------------------------------------------------

        Rvv0 = (epsi*cosl - sqrt(epsi-sinl*sinl))
     &        /(epsi*cosl + sqrt(epsi-sinl*sinl))
        Rhh0 = (cosl - sqrt(epsi-sinl*sinl))
     &        /(cosl + sqrt(epsi-sinl*sinl))
        Rvv0 = abs(Rvv0) 
        Rhh0 = abs(Rhh0) 
        Rvv0 = Rvv0*Rvv0
        Rhh0 = Rhh0*Rhh0

c-----------------------------------------------------------------------
c                       INTERPOLATION DE LA DIFFUSION
c
c-----------------------------------------------------------------------

c SI ON RESTE DANS LE DOMAINE TABULÉ
                        if (thetal.le.THETAmax) then
                                ind     = thetal/DTHETAtab
                                Poids   = (thetal/DTHETAtab-ind) 
c Fondamental
                        Diff(1) = Irh(0,1,ind)*(1.0D0-Poids)
     &                            + Irh(0,1,ind+1)*Poids
                        Diff(2) = Irh(0,2,ind)*(1.0D0-Poids)
     &                            + Irh(0,2,ind+1)*Poids
                        Diff(3) = Irh(0,3,ind)*(1.0D0-Poids)
     &                            + Irh(0,3,ind+1)*Poids
                        Diff(4) = Irh(0,4,ind)*(1.0D0-Poids)
     &                            + Irh(0,4,ind+1)*Poids
c 2eme Harmonique

                        Diff(1) = Diff(1) + (Irh(2,1,ind)*(1.0D0-Poids)
     &                                    + Irh(2,1,ind+1)*Poids
     &                                      )*dcosd(2.0D0*phiL)
                        Diff(2) = Diff(2) + (Irh(2,2,ind)*(1.0D0-Poids)
     &                                    + Irh(2,2,ind+1)*Poids
     &                                      )*dcosd(2.0D0*phiL)
                        Diff(3) = Diff(3) + (Irh(2,3,ind)*(1.0D0-Poids)
     &                                    + Irh(2,3,ind+1)*Poids
     &                                      )*dsind(2.0D0*phiL)
                        Diff(4) = Diff(4) + (Irh(2,4,ind)*(1.0D0-Poids)
     &                                    + Irh(2,4,ind+1)*Poids
     &                                      )*dsind(2.0D0*phiL)
c       print*, Diff(1),  Irh(0,1,ind), Irh(0,1,ind+1), Irh(2,1,ind)
c     &             , Irh(2,1,ind+1), dcosd(2.0D0*phil), phil
c SI ON SORT DU DOMAINE TABULÉ À LA PRECISION NUMÉRIQUE PRÈS
                        elseif ((P_Sx_Sy.eq.0.D0).or.
     &                          (thetal.le.THETAmax*1.0001D0)) then
                                do indice = 1, 4
                                        Diff(indice) = 0.0D0
                                enddo
                        else
c SI ON SORT DU DOMAINE, ERREUR => SORTIE DU PROGRAMME
                                print *, 'ERREUR !!!!!!!!!!!!!'
                                print *, 'thetal > THETAmax'
                                print *, thetal, ' > ', THETAmax
                                print *, 'Densité de Proba = ',P_Sx_Sy
                                stop
                        endif


c-----------------------------------------------------------------------c                       VECTEUR DE STOKES GLOBAL
c

c  CONTRIBUTION DE L'ECUME ET DE LA DIFFUSION DE BRAGG AU VECTEUR LOCAL

c--- Vecteur de stokes dans le repère local
c  		AVEC ATMOSPHERE
c                       avec ecume

          		Isl_eA(1) = (Tw*(1.0D0-(Rvv0+h*Diff(1)))*Fr_
     &                   + epsi_sf(1)*Fr + AtmoD*Rvv0)*dexp(-tau_)
          		Isl_eA(2) = (Tw*(1.0D0-(Rhh0+h*Diff(2)))*Fr_
     &                   + epsi_sf(2)*Fr + AtmoD*Rhh0)*dexp(-tau_)

c                       sans ecume

          		IslA(1) = (Tw*(1.0D0-(Rvv0+h*Diff(1)))
     &                                  + AtmoD*Rvv0)*dexp(-tau_) 
          		IslA(2) = (Tw*(1.0D0-(Rhh0+h*Diff(2)))
     &                                  + AtmoD*Rhh0)*dexp(-tau_) 
          		IslA(3) = (Tw*(0.D0-h*Diff(3)))*dexp(-tau_)
                        IslA(4) = (Tw*(0.D0-h*Diff(4)))*dexp(-tau_)
c  		SANS ATMOSPHERE
c                       avec ecume

          		Isl_e(1) = Tw*(1.0D0-(Rvv0+h*Diff(1)))*Fr_
     &                   + epsi_sf(1)*Fr
          		Isl_e(2) = Tw*(1.0D0-(Rhh0+h*Diff(2)))*Fr_
     &                   + epsi_sf(2)*Fr

c                       sans ecume

          		Isl(1) = Tw*(1.0D0-(Rvv0+h*Diff(1)))
          		Isl(2) = Tw*(1.0D0-(Rhh0+h*Diff(2)))
          		Isl(3) = Tw*(0.D0-h*Diff(3))
       			Isl(4) = Tw*(0.D0-h*Diff(4))

c--- Vecteur de stokes dans le repère global
c  		AVEC ATMOSPHERE
c                       avec ecume

                        Is_eA(1) = Isl_eA(1)*cosalpha**2
     &                           + Isl_eA(2)*sinalpha**2
     &                           + IslA(3)*cosalpha*sinalpha
          		Is_eA(2) = Isl_eA(1)*sinalpha**2
     &                           + Isl_eA(2)*cosalpha**2
     &          		 - IslA(3)*cosalpha*sinalpha
                        Is_eA(3) = IslA(3)*(cosalpha**2-sinalpha**2) -
     &                      (Isl_eA(1)-Isl_eA(2))*2.D0*sinalpha*cosalpha
                        Is_eA(4) = IslA(4)

c                       sans ecume

          		IsA(1) = IslA(1)*cosalpha**2 
     &                         + IslA(2)*sinalpha**2 
     &          	       + IslA(3)*cosalpha*sinalpha 
          		IsA(2) = IslA(1)*sinalpha**2 
     &                         + IslA(2)*cosalpha**2 
     &          	       - IslA(3)*cosalpha*sinalpha
                        IsA(3) = IslA(3)*(cosalpha**2-sinalpha**2) -
     &                      (IslA(1)-IslA(2))*2.D0*sinalpha*cosalpha
                        IsA(4) = IslA(4)

c  		SANS ATMOSPHERE
c                       avec ecume

          		Is_e(1) = Isl_e(1)*cosalpha**2 
     &                          + Isl_e(2)*sinalpha**2 
     &          		+ Isl(3)*cosalpha*sinalpha 
          		Is_e(2) = Isl_e(1)*sinalpha**2 
     &                          + Isl_e(2)*cosalpha**2 
     &          		- Isl(3)*cosalpha*sinalpha
                        Is_e(3) = Isl(3)*(cosalpha**2-sinalpha**2) -
     &                        (Isl_e(1)-Isl_e(2))*2.D0*sinalpha*cosalpha
                        Is_e(4) = Isl(4)

c                       sans ecume

          		Is(1) = Isl(1)*cosalpha**2 
     &                        + Isl(2)*sinalpha**2 
     &                        + Isl(3)*cosalpha*sinalpha 
          		Is(2) = Isl(1)*sinalpha**2 
     &                        + Isl(2)*cosalpha**2 
     &                        - Isl(3)*cosalpha*sinalpha
                        Is(3) = Isl(3)*(cosalpha**2-sinalpha**2) -
     &                        (Isl(1)-Isl(2))*2.D0*sinalpha*cosalpha
                        Is(4) = Isl(4)

c---Integrand Tb
c  		AVEC ATMOSPHERE
c                       avec ecume

          		int1eA(1,iphi) = int1eA(1,iphi) + 
     &                                 Is_eA(1)*projec*P_Sx_Sy*del1*del2
          		int1eA(2,iphi) = int1eA(2,iphi) + 
     &                                 Is_eA(2)*projec*P_Sx_Sy*del1*del2
          		int1eA(3,iphi) = int1eA(3,iphi) + 
     &                                 Is_eA(3)*projec*P_Sx_Sy*del1*del2
          		int1eA(4,iphi) = int1eA(4,iphi) + 
     &                                 Is_eA(4)*projec*P_Sx_Sy*del1*del2

c                       sans ecume

          		int1A(1,iphi) = int1A(1,iphi) + 
     &                                 IsA(1)*projec*P_Sx_Sy*del1*del2
          		int1A(2,iphi) = int1A(2,iphi) + 
     &                                 IsA(2)*projec*P_Sx_Sy*del1*del2
          		int1A(3,iphi) = int1A(3,iphi) + 
     &                                 IsA(3)*projec*P_Sx_Sy*del1*del2
          		int1A(4,iphi) = int1A(4,iphi) + 
     &                                 IsA(4)*projec*P_Sx_Sy*del1*del2

c  		SANS ATMOSPHERE
c                       avec ecume

          		int1e(1,iphi) = int1e(1,iphi) + 
     &                                 Is_e(1)*projec*P_Sx_Sy*del1*del2
          		int1e(2,iphi) = int1e(2,iphi) + 
     &                                 Is_e(2)*projec*P_Sx_Sy*del1*del2
          		int1e(3,iphi) = int1e(3,iphi) + 
     &                                 Is_e(3)*projec*P_Sx_Sy*del1*del2
          		int1e(4,iphi) = int1e(4,iphi) + 
     &                                 Is_e(4)*projec*P_Sx_Sy*del1*del2
c                       sans ecume

          		int1(1,iphi) = int1(1,iphi) + 
     &                                 Is(1)*projec*P_Sx_Sy*del1*del2
          		int1(2,iphi) = int1(2,iphi) + 
     &                                 Is(2)*projec*P_Sx_Sy*del1*del2
          		int1(3,iphi) = int1(3,iphi) + 
     &                                 Is(3)*projec*P_Sx_Sy*del1*del2
          		int1(4,iphi) = int1(4,iphi) + 
     &                                 Is(4)*projec*P_Sx_Sy*del1*del2


c  NORMALISATION DE LA DENSITE DE PROBABILITE : Norm2
                        Norm2 = Norm2 + projec*P_Sx_Sy*del1*del2
                enddo
                endif
                enddo
c  AVEC ATMO
c               avec ecume
                int1eA(1,iphi) = int1eA(1,iphi)/Norm2
                int1eA(2,iphi) = int1eA(2,iphi)/Norm2
                int1eA(3,iphi) = int1eA(3,iphi)/Norm2
                int1eA(4,iphi) = int1eA(4,iphi)/Norm2
c               sans ecume
                int1A(1,iphi) = int1A(1,iphi)/Norm2
                int1A(2,iphi) = int1A(2,iphi)/Norm2
                int1A(3,iphi) = int1A(3,iphi)/Norm2
                int1A(4,iphi) = int1A(4,iphi)/Norm2

c  SANS ATMO
c               avec ecume
                int1e(1,iphi) = int1e(1,iphi)/Norm2
                int1e(2,iphi) = int1e(2,iphi)/Norm2
                int1e(3,iphi) = int1e(3,iphi)/Norm2
                int1e(4,iphi) = int1e(4,iphi)/Norm2
c               sans ecume
                int1(1,iphi) = int1(1,iphi)/Norm2
                int1(2,iphi) = int1(2,iphi)/Norm2
                int1(3,iphi) = int1(3,iphi)/Norm2
                int1(4,iphi) = int1(4,iphi)/Norm2
        write (20,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
     &                  theta(itheta), phi, Tvn, Thn, int1(1,iphi)-Tvn,
     &               int1(2,iphi)-Thn ,int1(3,iphi),int1(4,iphi) ,0.0D0 
     &, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, lambdad
     &                  , Stab(1), realpart(epsi),imagpart(epsi), Su*Su,
     &                  Sc*Sc, Fr-Fr
        write (*,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
     &                  theta(itheta), phi, Tvn, Thn, int1(1,iphi)-Tvn,
     &               int1(2,iphi)-Thn,int1(3,iphi),int1(4,iphi) ,0.0D0 
     &, 0.0D0, 0.0D0, 0.0D0, 0.0D0, 0.0D0, lambdad
     &                  , Stab(1), realpart(epsi),imagpart(epsi), Su*Su,
     &                  Sc*Sc, Fr-Fr
c Fin boucle phi<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      		enddo
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c----------------------------------------------------------------------
c                       CALCUL DES HARMONIQUES AVEC ECUME
c
c int1(1,1) : Tb verticale en phi = TPHI(1) = 0°
c int1(1,2) : Tb verticale en phi = TPHI(2) = 45°
c int1(1,3) : Tb verticale en phi = TPHI(3) = 90°
c int1(1,4) : Tb verticale en phi = TPHI(4) = 180°
c Le premier indice de int1 defini l'indice du parametre de stokes et 
c varie de 1 à 4 respectivement pour Tv, Th, U et V.

        Tv0 = 0.25D0*(int1e(1,1) + 2.0D0*int1e(1,3) + int1e(1,4))
        Th0 = 0.25D0*(int1e(2,1) + 2.0D0*int1e(2,3) + int1e(2,4))
        Tv1 = 0.5D0 *(int1e(1,1) - int1e(1,4)) 
        Th1 = 0.5D0 *(int1e(2,1) - int1e(2,4))
        U1  = int1e(3,3) 
        V1  = int1e(4,3)
        Tv2 = 0.25D0*(int1e(1,1) - 2.0D0*int1e(1,3) + int1e(1,4))
        Th2 = 0.25D0*(int1e(2,1) - 2.0D0*int1e(2,3) + int1e(2,4))
        U2  = int1e(3,2) - int1e(3,3)
        V2  = int1e(4,2) - int1e(4,3)

c Ecriture fichier et ecran
c        write (10,1000) nu, SST(iSST), SSS(iSSS), U10, 
c     &                 ustar*100, theta(itheta), Tvn, Thn,
c     &                (Tv0-Tvn), (Th0-Thn), Tv1, Th1, U1, V1, Tv2, 
c     &                Th2, U2, V2, lambdad, Stab(1), realpart(epsi), 
c     &                  imagpart(epsi), Su*Su, Sc*Sc, Fr*100.0D0
c        write (*,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
c     &                 theta(itheta), Tvn, Thn,
c     &                (Tv0-Tvn), (Th0-Thn), Tv1, Th1, U1, V1, Tv2, 
c     &                Th2, U2, V2, lambdad, Stab(1), realpart(epsi), 
c     &                  imagpart(epsi), Su*Su, Sc*Sc, Fr*100.0D0

c----------------------------------------------------------------------
c                       CALCUL DES HARMONIQUES SANS ECUME
c
c int1(1,1) : Tb verticale en phi = TPHI(1) = 0°
c int1(1,2) : Tb verticale en phi = TPHI(2) = 45°
c int1(1,3) : Tb verticale en phi = TPHI(3) = 90°
c int1(1,4) : Tb verticale en phi = TPHI(4) = 180°
c Le premier indice de int1 defini l'indice du parametre de stokes et 
c varie de 1 à 4 respectivement pour Tv, Th, U et V.


        Tv0 = 0.25D0*(int1A(1,1)+2.0D0*int1A(1,3)+int1A(1,4))+AtmoU
        Th0 = 0.25D0*(int1A(2,1)+2.0D0*int1A(2,3)+int1A(2,4))+AtmoU
        Tv1 = 0.5D0 *(int1A(1,1) - int1A(1,4)) 
        Th1 = 0.5D0 *(int1A(2,1) - int1A(2,4))
        U1  = int1A(3,3) 
        V1  = int1A(4,3)
        Tv2 = 0.25D0*(int1A(1,1)-2.0D0*int1A(1,3)+int1A(1,4))
        Th2 = 0.25D0*(int1A(2,1)-2.0D0*int1A(2,3)+int1A(2,4))
        U2  = int1A(3,2) - int1A(3,3)
        V2  = int1A(4,2) - int1A(4,3)

c        Tv0 = 0.25D0*(int1e(1,1) + 2.0D0*int1e(1,3) + int1e(1,4))+AtmoU
c        Th0 = 0.25D0*(int1e(2,1) + 2.0D0*int1e(2,3) + int1e(2,4))+AtmoU
c        Tv1 = 0.5D0 *(int1e(1,1) - int1e(1,4)) 
c        Th1 = 0.5D0 *(int1e(2,1) - int1e(2,4))
c        U1  = int1e(3,3) 
c        V1  = int1e(4,3)
c        Tv2 = 0.25D0*(int1e(1,1) - 2.0D0*int1e(1,3) + int1e(1,4))
c        Th2 = 0.25D0*(int1e(2,1) - 2.0D0*int1e(2,3) + int1e(2,4))
c        U2  = int1e(3,2) - int1e(3,3)
c        V2  = int1e(4,2) - int1e(4,3)

c Ecriture fichier et ecran
c        write (20,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
c     &                  theta(itheta), Tvn, Thn, (Tv0-Tvn), (Th0-Thn),
c     &                  Tv1, Th1, U1, V1, Tv2, Th2, U2, V2, lambdad
c     &                  , Stab(1), realpart(epsi),imagpart(epsi), Su*Su,
c     &                  Sc*Sc, Fr-Fr
c        write (*,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
c     &                 theta(itheta), Tvn, Thn, (Tv0-Tvn), (Th0-Thn), 
c     &                 Tv1, Th1, U1, V1, Tv2, Th2, U2, V2, lambdad
c     &                 , Stab(1), realpart(epsi),imagpart(epsi), Su*Su,
c     &                 Sc*Sc, Fr-Fr

c Fin boucle Theta<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        enddo
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c  Fin Boucle kd
        enddo

c  Fin Boucle Vent
        enddo

c  Fin Boucle SST
        enddo


c  Fin Boucle SSS
        enddo


c----------------------------------------------------------------------
c                               FORMATS

  100     format (1x,f5.2,4(1x,f4.1),2(1x,f6.2),2(1x,f6.2)
     &          ,8(1x,f5.2))
 1000     format (1x,f9.5,1x,f5.2,1x,f6.3,2(1x,f6.2),1x,f4.1,4(1x,f7.3),
     &         8(1x,f6.3),1x,e9.3,1x,f5.2,2(1x,f6.2),2(1x,e9.3),1x,f5.1)

        close (10)
	close (15)
        close (20)
        close (25)
        close (30)
        close (40)

        end

c------------------------ Include  -------------------------------------

c        include "/home/dinnat/Code/Fortran/Modele2Echelles/readvec.f" 
c        include "/home/dinnat/Code/Fortran/Modele2Echelles/func2D1.f"
c        include "/home/dinnat/Code/Fortran/Numerical/extr_val.f"
c        include "/home/dinnat/Code/Fortran/Numerical/root.f"
c        include "/home/dinnat/Code/Fortran/Numerical/mylib.f"
c        include "/home/dinnat/Code/Fortran/Numerical/qromb2D1.f"
c        include "/home/dinnat/Code/Fortran/Permitivite/epsilon_KS.f"
c        include "/home/dinnat/Code/Fortran/Permitivite/epsilon_El.f"
c        include "/home/dinnat/Code/Fortran/Vent/vent.f"
c        include "/home/dinnat/Code/Fortran/Vent/funcUstar_y.f"
c        include "/home/dinnat/Code/Fortran/Vent/funcUstar_c.f"
c        include "/home/dinnat/Code/Fortran/SeaState/c.f"
c        include "/home/dinnat/Code/Fortran/SeaState/Su.f"
c        include "/home/dinnat/Code/Fortran/SeaState/Sc.f"
c        include "/home/dinnat/Code/Fortran/SeaState/P.f"
c        include "/home/dinnat/Code/Fortran/SeaState/Sigma.f"
c        include "/home/dinnat/Code/Fortran/Ecume/foam.f"
c        include "/home/dinnat/Code/Fortran/Ecume/esf.f"
c        include "/home/dinnat/Code/Fortran/Atmosphere/TbAtmo.f"
c        include "/home/dinnat/Code/Fortran/Geometrie/Glob2Loc.f"

c------------- Inlude des subroutines  ---------------------------------
c        include "/home/dinnat/Code/Fortran/SeaState/Swell/infSwell.f"
c dans func2D1.f
c        include "/home/dinnat/Code/Fortran/Modele2Echelles/func2D21.f"
c        include "/home/dinnat/Code/Fortran/Modele2Echelles/func2D22.f"
c        include "/home/dinnat/Code/Fortran/Numerical/qromb2D2.f"
c dans root.f
c        include "/home/dinnat/Code/Fortran/Numerical/dichoto.f"
c dans qromb2D1.f
c        include "/home/dinnat/Code/Fortran/Numerical/trapzd2D1.f"
c        include "/home/dinnat/Code/Fortran/Numerical/polint.f"
c dans c.f
c        include "/home/dinnat/Code/Fortran/Numerical/qromb.f"
c dans func2D21.f et func2D22.f
c        include "/home/dinnat/Code/Fortran/EM/Gvv2.f"
c        include "/home/dinnat/Code/Fortran/EM/Gvv1.f"
c        include "/home/dinnat/Code/Fortran/EM/Ghh2.f"
c        include "/home/dinnat/Code/Fortran/EM/Ghh1.f"
c        include "/home/dinnat/Code/Fortran/EM/Gvh1.f"
c        include "/home/dinnat/Code/Fortran/EM/Ghv1.f"
c        include "/home/dinnat/Code/Fortran/EM/Ghv2.f"
c        include "/home/dinnat/Code/Fortran/SeaState/spectre_DV.f"
c        include "/home/dinnat/Code/Fortran/SeaState/spectre_El.f"
c        include "/home/dinnat/Code/Fortran/SeaState/spectre_powerLaw.f"
c dans qromb.f
c        include "/home/dinnat/Code/Fortran/Numerical/trapzd.f"
c        include "/home/dinnat/Code/Fortran/SeaState/fksup.f"
c        include "/home/dinnat/Code/Fortran/SeaState/fkinf.f"
c dans qromb2D2.f
c        include "/home/dinnat/Code/Fortran/Numerical/trapzd2D2.f"
c dans TbAtmo.f
c        include "/home/dinnat/Code/Fortran/Atmosphere/paramAtmo.f" ! calcul parmetres atmospheriques
c        include "/home/dinnat/Code/Fortran/Atmosphere/tauP.f" ! epaisseur optique atmo (Neper)
c        include "/home/dinnat/Code/Fortran/Atmosphere/convDpath.f"   ! conversion (alt,theta) -> dpath
c dans convDpath.f
c        include "/home/dinnat/Code/Fortran/Numerical/zroots.f"
c dans zroots.f
c        include "/home/dinnat/Code/Fortran/Numerical/laguer.f"
