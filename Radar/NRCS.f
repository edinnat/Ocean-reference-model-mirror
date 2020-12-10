c Modif
c Norm du spectre = 1 (supprimer la normalisation, ne devrait pas exister)      

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
        common /paramLem/ paramLem
        common /paramKudryavtsev/ paramKudryavtsev
        common /filtersKudry/ filtersKudry
c Définition des fonction externes        
       	external funcUstar_y
       	external funcUstar_c
        external fkinf
        external fksup
        external datand
     	external dcosd
     	external dsind
     	external dtand
     	external dacosd
     	external dasind
     	external epsilon_KS
     	external epsilon_El
     	external func2D1
c Déclaration des variables 
        character*80 fout1 ! fichier de sortie avec ecume
        character*80 fout2 ! fichier de sortie sans ecume
        character*80 fin1  ! fichier d'entrée des parametres
        character*80 pathout ! chemin pour les fichiers de sortie
        character*80 KudryFilter1 ! filename of lookup table for Kudryavtsev filter function 1
        character*80 KudryFilter2 ! filename of lookup table for Kudryavtsev filter function 2
        character*80 KudryFilter3 ! filename of lookup table for Kudryavtsev filter function 3
        character*1 cVar
        character*1  cepsi
        character*1  csigneeps
        character*1  cSpectre
        character*1  cType
        character*16  cCouvEcume
        character*1  cEmisEcume
        character*1  cCD
     	character*1  cAtmo
     	character*1  cSwell
        character*80 jump_line
        character*40 Modepsi
        character*40 ModVarPente
        character*40 ModSpectre
        character*40 TypeMod
        character*40 ModCouvEcume
        character*40 ModEmisEcume
        character*40 ModCD
        character*10 date
        character*9 time
	
        double precision paramLem(10)
      double precision paramKudryavtsev(10)
      double precision filtersKudry(1:3,1:50003)
        double precision ss, Pi, beta
        double precision kmin
        double precision phimin
        double precision phimax
        double precision datand
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
     	double precision dIr(1:2)
     	double precision Irh (0:2,1:4,0:18000)
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
        double precision kinf
        double precision ksup
        double precision kmid
        double precision ustar
        double precision ustar_y
        double precision ustar_c
        double precision ustar_don
        double precision ustar_Kudry
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
        double precision sigmaQS
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
        double precision int1(1:4,1:4)
        double precision int1e(1:4,1:4)
        double precision int1A(1:4,1:4) ! int1 avec atmosphere
        double precision int1eA(1:4,1:4) ! int1e avec atmosphere
        double precision cosalpha
        double precision sinalpha 
        double precision P_Sx_Sy
        double precision Diff(1:4)
        double precision Poids
        double precision DTHETAtab
        double precision lambdamin
        double precision lambdamax
        double precision niVar
        double precision sigmaVV0
        double precision sigmaHH0
        double precision sigmaVV1
        double precision sigmaHH1
        double precision sigmaHV1 
        double precision sigmaVH1
        double precision sigmaVV2
        double precision sigmaHH2
        double precision sigmaHV2 
        double precision sigmaVH2
        double precision Tvn
        double precision Thn
     	double precision Rveff
     	double precision Rheff
     	double precision Reff
        double precision xVar
        double precision xSigma
        double precision VarSwell(1:2) ! variances des pentes de la houle en x et y
        double precision hSwell ! rms des hauteurs de la houle
        double precision sigSwell(1:2) ! largeur à mi-puissance de la densite spectrale de la houle en x et y
        double precision KMaxSwell(1:2) ! pics de la densite spectrale de la houle en x et y
        double precision temp
        double precision temp1
        double precision temp2
        double precision temp3
        double precision temp4
        double precision Omega
        
       	double complex epsi_KS
       	double complex epsi_El
       	double complex epsi_MW
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
c		sampling of the sea spectrum
		integer ik					! index of wavenumber sample
		integer nk					! numb. of samples in k (log spaced)
		double precision dlogk		! step in log(k)
		double precision k			! spectrum wavenumber 
		double precision S_DV		! function for spectrum
		double precision S_Yueh		! function for spectrum
		
        integer iIncid
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
        integer ikfilt ! index for Kudryavtsev filters
        
        logical fin1exist
	
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
c        PI2=2.0D0*acos(-1.0D0)
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
c----------- Process Time and Date of Program Execution -----------        
        call date_and_time(date, time)
        date = date(7:8)//'_'//date(5:6)//'_'//date(1:4)
        time = time(1:2)//'_'//time(3:4)//'_'//time(5:6)

c------------ Management of Input Output Files ----------------
c Lookup table for Kudryavtsev filters        
        KudryFilter1 = 'Data/Kudryavtsev_Spectrum/filterF.dat'
        KudryFilter2 = 'Data/Kudryavtsev_Spectrum/filterFres.dat' 
        KudryFilter3 = 'Data/Kudryavtsev_Spectrum/filterPhi.dat'  
c Identify file of input parameters
        call getarg (1, fin1)                   ! get 1st parameter
                                                ! passed to function

        inquire (file = fin1, exist = fin1exist)
        if (.not. fin1exist) then
	  call getenv ('inputFortran', fin1)
          if (fin1.eq.'') then
            print*, 'Pas de fichier de parametre specifié ou fichier non
     & valide. De plus la variable d''environnement ''inputFortran''
     & pour le fichier de sortie n''est pas définie.'
            print*, 'Fin du programme'
            stop
          else
            fin1 = fin1(1:lnblnk(fin1))//"/NRCS.p"
            inquire (file = fin1, exist = fin1exist)
            if (.not. fin1exist) then
            print*, 'Invalid input file and invalid default file: ',
     &  fin1
            print*, '-- ABORT -- '
            stop
            else    
            print*, 'Pas de nom de fichier input specifié ou valide.'
            print*, 'Utilisation de fichier par defaut ', fin1
            endif
          endif
	endif

     	open (unit=30,file=fin1,status='unknown',err=60, action='read')
        call getenv ('outputFortran', pathOut)
        if (pathOut.eq.'') then
            print*, 'La variable d''environnement ''outputFortran'' pour
     & le chemin des fichiers de sortie n''est pas définie.'
            print*, 'Fin du programme'
            stop
        endif
        pathout = pathout(1:lnblnk(pathout))//'/'
        open (unit=40, file=pathout(1:lnblnk(pathout))//"integPts.dat"
     &,status='unknown',err=70)
c        open (unit=40, file=pathout(1:lnblnk(pathout))//"Diffusion.dat"
c     &,status='unknown',err=70)
        open (unit=50, file=pathout(1:lnblnk(pathout))//"rapport.dat", 
     &status='unknown',err=80)

c----------------------------------------------------------------------

c------------------ ENTREE DES VARIABLES  -----------------------------
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,*) lambda
       	freq = 3.D08/lambda
       	nu = 3.D-01/lambda
        call readvec(kdtab, 30, nkd)
        call extr_val(kdtab, nkd, kdmin, kdmax)
        call readvec(theta, 30, ntheta)
        call extr_val(theta, ntheta, thetamin, theta_max)
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
        call readvec(Wind, 30, nWind)
        call extr_val(Wind, nWind, Windmin, Windmax)
     	read (30,'(a)') jump_line
       	read (30,*) alt
     	read (30,'(a)') jump_line
        read (30,*) cCD               ! cCD mod. coefficient trainee
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
        call readvec(SST, 30, nSST)
        call extr_val(SST, nSST, SSTmin, SSTmax)
        call readvec(SSS, 30, nSSS)
        call extr_val(SSS, nSSS, SSSmin, SSSmax)
        call readvec(Stab, 30, nStab)
        call extr_val(Stab, nStab, Stabmin, Stabmax)
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') fout1
     	read (30,'(a)') jump_line
     	read (30,'(a)') fout2
        nsorties = 0
        if (fout1.eq.'aucun') then
                print *, '  Pas de sortie fichier avec écume'
                if (fout2.eq.'aucun') then
                        print *, '  Pas de sortie fichier sans écume'
                        print *, '!!!!!!! ATTENTION  !!!!!!'
                        print *, 'Pas de sortie fichier prévue'
                else
                   open (unit=20,file=pathout(1:lnblnk(pathout))//fout2
     &,status='unknown')
c                   open (unit=25,file=pathout(1:lnblnk(pathout))//"Atm_"
c     &//fout2,status='unknown')
                        nsorties = 1
                        sortiestab(1) = 20
                endif
        else
                open (unit=10,file=pathout(1:lnblnk(pathout))//fout1,sta
     &tus='unknown')
                open (unit=15,file=pathout(1:lnblnk(pathout))//"Atm_"//
     &fout1,status='unknown')
                if (fout2.eq.'aucun') then
                        print *, '  Pas de sortie fichier sans écume'
                        nsorties = 1
                        sortiestab(1) = 10
                else
                   open (unit=20,file=pathout(1:lnblnk(pathout))//fout2,
     &status='unknown')
c                   open (unit=25,file=pathout(1:lnblnk(pathout))//"Atm_"
c     &//fout2,status='unknown')
                        nsorties = 2
                        sortiestab(1) = 10
                        sortiestab(2) = 20
                endif
        endif
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
        read (30,'(a)') cVar
     	read (30,'(a)') jump_line     ! B_ Amplitude du spectre DV
        read (30,*) B_
     	read (30,'(a)') jump_line
        read (30,*) cSpectre          ! cSpectre mod. de spectre
     	read (30,'(a)') jump_line
        read (30,*) Omega          ! Inverse de l'age des vagues
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line     ! cCouvEcume mod. couverture ecume
        read (30,*) cCouvEcume
     	read (30,'(a)') jump_line
        read (30,*) cEmisEcume        ! cEmisEcume mod. emissiv. ecume  
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
        read (30,*) cepsi             ! cepsi mod. permitivite
     	read (30,'(a)') jump_line
        read (30,*) csigneeps         ! csigneeps Signe Im(permitivite)
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
        read (30,*) cSwell
     	read (30,'(a)') jump_line
        read (30,*) NSwell(1)
     	read (30,'(a)') jump_line
        read (30,*) NSwell(2)
     	read (30,'(a)') jump_line
        read (30,*) hSwell
     	read (30,'(a)') jump_line
        read (30,*) sigSwell(1)
     	read (30,'(a)') jump_line
        read (30,*) sigSwell(2)
     	read (30,'(a)') jump_line
        read (30,*) KMaxSwell(1)
     	read (30,'(a)') jump_line
        read (30,*) KMaxSwell(2)
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
        read (30,*) THETAmax, DTHETAtab
        if (DTHETAtab.lt.0.010D0) then
                write (*,*)
                print *, '!!!!!! ERREUR  !!!!!!!!!!!'
                print *, 'DTHETAtab doit etre au moins égal à 1'
                write (*,*)
                stop
        endif
     	read (30,'(a)') jump_line     ! lambdamin
        read (30,*) lambdamin
     	read (30,'(a)') jump_line     ! lambdamax
        read (30,*) lambdamax
     	read (30,'(a)') jump_line     ! niVar
        read (30,*) niVar
        read (30,'(a)') jump_line     ! cType Type de Modele 
        read (30,*) cType
     	read (30,'(a)') jump_line
        read (30,*) N1                ! Nb integration sur Sx
     	read (30,'(a)') jump_line
        read (30,*) N2                ! Nb integration sur Sy
     	read (30,'(a)') jump_line
        read (30,*) xVar              ! Limite integration sur pentes (xVar x Variance)
        xSigma = sqrt(xVar) ! (limite = xSigma x Ecart type)
        
        call date_and_time(date, time)
        date = date(7:8)//'/'//date(5:6)//'/'//date(1:4)
        time = time(1:2)//':'//time(3:4)//':'//time(5:6)


c----------- VERIFICATION DES PARAMETRES D'ENTRÉE ----------------------

c        if (int(THETAmax/DTHETAtab).ne.(THETAmax/DTHETAtab)) then
c                write (*,*)
c                print *, '!!!!!! ERREUR  !!!!!!!!!!!'
c                print *, 'THETAmax n''est pas multiple de DTHETAtab'
c                write (*,*) 'THETAmax = ', THETAmax
c                write (*,*) 'DTHETAtab = ', DTHETAtab
c                write (*,*)
c                stop
c        endif
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
c        elseif ((cSpectre.eq.'P').or.(cSpectre.eq.'p')) then
c                cspec = 3
c                ModSpectre = 'Power law'
        elseif ((cSpectre.eq.'L').or.(cSpectre.eq.'l')) then
                cspec = 4
                ModSpectre = 'Lemaire et al. (1999)'                
        elseif ((cSpectre.eq.'K').or.(cSpectre.eq.'k')) then
                cspec = 5
                ModSpectre = 'Kudryavtsev et al. (2003)'
        else
                write (*,*)
                print *, '!!!!!!!! ERREUR  !!!!!!!!!!!!!!!!'
                print *, 'Le choix du spectre est ', cSpectre
                print *, 'Choix possibles : E pour Elfouhaily (97)'
                print *, '                : D pour Durden & Vesecky(85)'
                print *, '                : L pour Lemaire et al. (99)'
                print *, '                : K Kudryavtsev et al. (2003)'
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
        elseif ((cCouvEcume.eq.'M3').or.(cCouvEcume.eq.'m3')) then
                ModCouvEcume = 'Monahan et Lu (90, actif + passif)'
                fCouvEcume = 3
        elseif ((cCouvEcume.eq.'M4').or.(cCouvEcume.eq.'m4')) then
                ModCouvEcume = 'Monahan et Lu (90, actif)'
                fCouvEcume = 4
        elseif ((cCouvEcume.eq.'M5').or.(cCouvEcume.eq.'m5')) then
                ModCouvEcume = 'Données WISE 2001'
                fCouvEcume = 5
        else
                write(*,*)
                print *,'Mauvais choix pour la couverture d''ecume:'
                print *,'     A  : pas d''ecume'
                print *,'     M1 : Monahan (86, eq.5)'
                print *,'     M2 : Monahan moindre carres (86, eq. 3a)'
                print *,'     M3 : Monahan et Lu (90, actif + passif)'
                print *,'     M4 : Monahan et Lu (90, actif)'
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
        elseif ((cCD.eq.'d').or.(cCD.eq.'D')) then
                ModCD = 'Donelan (93)'
                fCD = 3
        else
                write(*,*)
                print *,'Mauvais choix pour le coefficient de trainee:'
                print *,'     y : Cardone (69)'
                print *,'     c : Charnock (55)' 
                print *,'     d : Donelan (93)' 
                print *,'Votre choix etait : ',cCD
                write(*,*)
                stop
        endif
        if ((cSwell.eq.'O').or.(cSwell.eq.'o')) then
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
        write(50,*) 'k0 = ', k0, ' ;'
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
        write (ufile,*) '! freq     SST   SSS     U10   ustar theta  '
     &  //' sigmaVV0   sigmaHH0  sigmaVV1  sigmaHH1 '
     &  //' sigmaHV1      sigmaVH1    '
     &  //' sigmaVV2   sigmaHH2    sigmaHV2      sigmaVH2'
     &  //' lam_d/lam_0 Reff/Rfr sigma_Bragg_VV sigma_Bragg_HH ;'
        write (ufile,*) '! constante dielectrique e = Re(e) + i Im(e) $'

        enddo
        
c---------------------------------------------------------------------
c------------------ PROGRAMME PRINCIPAL ------------------------------

c  Calcul des Variance de pentes de la houle

        write(50,*)
        write(50,*) 'c: ============== Report ========================='
        write(50,*) 'c: Swell Slope variances'
        if (fSwell.eq.1) then
                call infSwell(NSwell, hSwell, sigSwell , 
     &          KMaxSwell, VarSwell)
        else
                VarSwell(1) = 0.0D0
                VarSwell(2) = 0.0D0
        endif
        write(50,*) 'Ecart type des pentes de la houle terminée'
        write(50,*) 'sigma(Swell) Upwind    = ', sqrt(VarSwell(1)),
     &                ' => variance = ', VarSwell(1)
        write(50,*) 'sigma(Swell) Crosswind = ', sqrt(VarSwell(2)),
     &                ' => variance = ', VarSwell(2)

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
c       epsi_MW pour Meissner and Wentz 
c SST  : Temperature de surface de l'ocean en °C
c SSS  : Salinité de surface de l'ocean en ppm
c freq : Frequence electromagnetique du radiometre en Hz
        
        write(50,*) 'Dielectric constant'
        call epsilon_KS (SST(iSST), SSS(iSSS), epsi_KS, freq)
        call epsilon_El (SST(iSST), SSS(iSSS), epsi_El, freq)
        call Epsilon_MW (SST(iSST), SSS(iSSS), epsi_MW, freq)
        if ((cepsi.eq.'K').or.(cepsi.eq.'k')) then
                epsi = epsi_KS
                Modepsi = 'Klein & Swift (77)'
        elseif ((cepsi.eq.'e').or.(cepsi.eq.'E')) then
                epsi = epsi_El
                Modepsi = 'Ellison (98)'
        elseif ((cepsi.eq.'m').or.(cepsi.eq.'M')) then
                epsi = epsi_MW
                Modepsi = 'Meissner et al. (2004,2012,2014)'
        else
                write(*,*)
                print *,'Mauvais choix pour constante dielectrique:'
                print *,'     K : Klein & Swift'
                print *,'     E : Ellison'
                print *,'     M : Meissner et al. (2004,2012,2014)'
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
        write(50,*) 'Calcul de constante dielectrique terminée'
        write(50,*) '  Klein & Swift = ', epsi_KS
        write(50,*) '  Ellison       = ', epsi_El
        write (50,*) '  Meissner et al. = ', epsi_MW
        write(50,*) '      =>       ',epsi


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
c ustar_don calculé avec Z0 Donelan 93(inclue le Fetch)

         write(50,*) 'Calcul du stress pour U(',alt,') = ', Uz
      	if (Uz.ne.0.) then
                call root (1.D-05, 1.5D1, 1.D-06,funcUstar_y, ustar_y)
                call root (1.D-05, 1.5D1, 1.D-06,funcUstar_c, ustar_c)
       ustar_don = 0.4D0*Uz/dlog(alt/(3.7D-05*(Uz**2/g)
     &         *Omega**0.9D0));
       	else
                ustar = 0.D0
       	endif
        if (fCD.eq.1) then
                ustar = ustar_y
        elseif (fCD.eq.2) then
                ustar = ustar_c
        elseif (fCD.eq.3) then
                ustar = ustar_don
        else
                print *, 'Je n''ai pas compris quel parametrisation'
                print *, 'vous desirez pour le coefficient de trainee'
                print *, 'Votre choix etait : ', cCD
                print *, 'Choix Possibles : y pour Cardone (69)'
                print *, '                  c pour Charnock (55)'
                print *, '                  d pour Donelan (93)'
                write (*,*)
                stop
        endif
        if (ustar.eq.1.0D04) print *, 'Changer les bornes passée à root'
     &          ,'pour le calcul de u*'
         write(50,*) 'Calcul du stress pour U(',alt,') = ', Uz
     &   , ' terminée'
         write(50,*) '   ustar = ', ustar

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

         write(50,*) 'Conversion du vent'
       	if (Uz.ne.0.) then
                call vent (alt,1.95D01,Uz,U19,ustar)
                call vent (alt,1.25D01,Uz,U12,ustar)
                call vent (alt,1.0D01 ,Uz,U10,ustar)
       	else 
                U19 = 0.
                U12 = 0.
       	endif
         
         write(50,*) 'Conversion du vent terminée'
         write(50,*) '   U10 = ', U10
         write(50,*) '   U12.5 = ', U12
         write(50,*) '   U19.5 = ', U19
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
        write(50,*) 'Calcul de la borne inférieure du spectre'
        call root (1.0D-03, 2.0D0, 1.D-06,fkinf, kinf)
!        kinf = 1.0D-03 
		ksup = 1.0D+03
        write(50,*) 'Calcul de la borne inférieure du spectre terminé'
        write(50,*) '   kinf = ', kinf
        kd = max(kinf, kd)
        write(50,*) ' => kd = ', kd
        write(50,*) 'Calcul de ''c'' du spectre de Durden & Vesecky'
        call c_ (c, U12)
        write(50,*) 'Calcul de ''c'' du spectre de Durden & Vesecky'
     &              ,' terminé'
        write(50,*) '   c = ',c

c		Sampling of sea spectra

		nk = 100
		dlogk = (dlog10(ksup) - dlog10(kinf))/(nk-1)
       if (nk.lt.1) then
			print*, 'Error: nk must be greater or equal to 1'
			stop
       endif

c		write(50,*)
c		write(50, *) 'Spectra'
c		write(50, *) 'col 1: U m/s at 10 '
c		write(50, *) 'col 2: k (rad/m)'
c		write(50, *) 'col 3: Psi(k) = S(k)/2/pi/k'
c		write(50,*)


c /* Durden and Vesecky (1985) 
c		write(50, *) 'Durden and Vesecky (1985) spectrum'
c		write(50,*)

c		do ik = 1, nk
c			k = 10.0D0**(dlog10(kinf) + dlogk*(ik-1))
c			write(50, '(f5.2,2(1x,e12.5))') U10, k, S_DV(k)/k/PI2
c		enddo 

c /* Yueh (1997)

c		write(50, *) 'Yueh (1997) spectrum'
c		write(50,*)
c		do ik = 1, nk
c			k = 10.0D0**(dlog10(kinf) + dlogk*(ik-1))
c			write(50, '(f5.2,2(1x,e12.5))') U10, k, S_Yueh(k)/k/PI2
c		enddo 

!        print*, kd
c        if (kinf.eq.1.0D04) print *, 'Changer les bornes passée à root '
c     &          ,'pour le calcul de kinf'
c        call root (2.0D0, 1.0D+05, 1.D-06,fksup, ksup)
c        if (ksup.eq.1.0D04) print *, 'Changer les bornes passée à root '
c     &          ,'pour le calcul de ksup'
c        print *, 'kinf = ', kinf, '  S_DV(kinf)/kinf', S_DV(kinf)/kinf
c        print *, 'ksup = ', ksup, '  S_DV(ksup)', S_DV(ksup)

c----------------------------------------------------------------------
c       PARAMETRES DU MODELE DE LEMAIRE ET AL.
c

      ss = 0.007  
      Pi = dacos(-1.0d0)
      beta = (2*Pi*ss)**2*5.0d0 ! fully developed, Pierson
c params
      paramLem(1) = ss ! significant slope
      paramLem(2) = 5 ! shape factor m
      paramLem(3) = beta ! amplitude factor
      paramLem(4) = 5e6 ! dimensional fetch (m)
      paramLem(5) = U10
      paramLem(6) = U12
      paramLem(7) = ustar

c----------------------------------------------------------------------
c       PARAMETRES DU MODELE DE Kudryavtsev ET AL.
c

! compute ustar
      temp1 = 1d-3*(0.8d0 + 0.065*U10) ! cd10
      ustar_Kudry = dsqrt(temp1)*U10   !
      temp2 = 0d0 ! # iteration
      temp3 = 1d0 ! eps
      do while (temp3 .gt. 0.01)
        temp2 = temp2 +1d0 ! 
        temp4 = 0.1*14d-6/ustar_Kudry+0.018*ustar_Kudry**(2)/9.81d0 ! z0
        temp1 = (0.41d0/dlog(10.d0/temp4))**2
        temp = dsqrt(temp1) * U10 ! usn
        temp3 = dabs(ustar_Kudry - temp)/ustar_Kudry
        ustar_Kudry = temp
       enddo      
      paramKudryavtsev(1) = U10
      paramKudryavtsev(2) = ustar_Kudry 
! load lookup tables for Kudryavtsev filter functions
      open (unit=11,file=KudryFilter1,status='old',err=61)
      open (unit=12,file=KudryFilter2,status='old',err=62)
      open (unit=13,file=KudryFilter3,status='old',err=63)
      do ikfilt = 1, 50000
           read(11,*) temp, filtersKudry(1,ikfilt+3) ! filter F, temp is for k
           read(12,*) temp, filtersKudry(2,ikfilt+3) ! filter Fres
           read(13,*) temp, filtersKudry(3,ikfilt+3) ! filter Phi
      enddo
      filtersKudry(1,1) = 1d-3 ! kmin
      filtersKudry(1,2) = 5d4  ! kmax
      filtersKudry(1,3) = 5d4  ! nk
      close(11)
      close(12)
      close(13)



c----------------------------------------------------------------------
c               VARIANCES DES PENTES DU SPECTRE
c
c  Su : Variances des pentes en direction "upwind"
c  Sc : Variances des pentes en direction "downwind"
c  kd : Nombre d'onde de coupure en rad/m
c  cspec : modele de spectre

        write(50,*) 'Calcul des variances des pentes'
        if (fVar.eq.1) then
            if (cspec.eq.5) then
               call SuScKud (Su, Sc, 1d-3, kd)
            else
                call Su_ (Su, 1.0D-03, kd, cspec)
                call Sc_ (Sc, 1.0D-03, kd, cspec)
            endif
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
        write(50,*) 'Calcul des variances des pentes terminé'
        write(50,*) '   Su = ', Su, ' => ', 'Su2 = ', Su*Su
        write(50,*) '   Sc = ', Sc, ' => ', 'Sc2 = ', Sc*Sc

c----------------------------------------------------------------------
c               VARIANCES DES HAUTEURS DU SPECTRE DES PETITE VAGUES
c
c  sigma : Variances des hauteurs des petites échelles
c  kd : Nombre d'onde de coupure en rad/m
c  cspec : modele de spectre

        write(50,*) 'Calcul des variances des hauteurs'
        call sigma_ (sigma, kd, cspec)
        write(50,*) 'Calcul des variances des hauteurs terminé'
        write(50,*) '   sigma =', sigma
c----------------------------------------------------------------------
c   Calcul du coefficient de réflexion de Fresnel effectif au nadir
        Reff = abs((1-sqrt(epsi))/(1+sqrt(epsi)))**2
     &          *dexp(-4.0D0*(sigma*k0)**2)
        write(50,*) 'R(0) = ', abs((1-sqrt(epsi))/(1+sqrt(epsi)))**2 
        write(50,*) 'Reff = ', Reff
        write(50,*) 'Reff/R(0) = ', dexp(-4.0D0*(sigma*k0)**2)

c----------------------------------------------------------------------
c                       FRACTION DE COUVERTURE D'ECUME
c

c  Couverture d'ecume avec dependance en Stab(stabilite atmos.)
        write(50,*) 'Calcul de la fraction d''écume'
        if (fCouvEcume.eq.1) then
                call foam (U10,Stab(1),Fr)
        elseif (fCouvEcume.eq.2) then
                Fr = 2.95D-6*U10**(3.52)
        elseif (fCouvEcume.eq.3) then
                call MonahanLu(U10, SST(iSST), temp1, temp2, Fr)
        elseif (fCouvEcume.eq.4) then
                call MonahanLu(U10, SST(iSST), Fr, temp1, temp2)
        elseif (fCouvEcume.eq.5) then
                call WISE2001(U10, Fr)
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
        write(50,*) 'Calcul de la fraction d''écume terminé'
        write(50,*) '   Fraction d''écume = ', Fr

c  Couverture d'ecume par methode des moindres carres

c-----------------------------------------------------------------------
c                           Lecture de tabulations de la diffusion
c        print *, 'Lecture de la diffusion engagée ...'
c        do iIncid = 0,THETAmax/DTHETAtab
c          read (33,*) Irh(0, 1, iIncid)
c     &     , Irh(0, 2, iIncid)
c     &     , Irh(2, 1, iIncid)
c     &     , Irh(2, 2, iIncid)
c     &     , Irh(2, 3, iIncid)
c     &     , Irh(2, 4, iIncid)
c        Irh(0, 1, iIncid) = Irh(0, 1, iIncid)/288.15D0
c        Irh(0, 2, iIncid) = Irh(0, 2, iIncid)/288.15D0
c        Irh(2, 1, iIncid) = Irh(2, 1, iIncid)/288.15D0
c        Irh(2, 2, iIncid) = Irh(2, 2, iIncid)/288.15D0
c        Irh(2, 3, iIncid) = Irh(2, 3, iIncid)/288.15D0
c        Irh(2, 4, iIncid) = Irh(2, 4, iIncid)/288.15D0
        
c        enddo
c-----------------------------------------------------------------------
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

c        print *, 'Tabulation de la diffusion engagée ...'
c        write(50,*)
c        write(50,*) 'c: Lookup table for scattering coefficients'
c        write(50,*) 'c: Format:'
c        write(50,*) 'c: col 1: frequency (Hz)'
c        write(50,*) 'c: col 2: sea surface temperature (degC)'
c        write(50,*) 'c: col 3: sea surface salinity (psu)'
c        write(50,*) 'c: col 4: wind speed at 10 meters height (m/s)'
c        write(50,*) 'c: col 5: wind friction velocity (m/s)'
c        write(50,*) 'c: col 6: incidence angle (deg)'
c        write(50,*) 'c: col 7: scat. coef., VV pol., omnidirectional'
c        write(50,*) 'c: col 8: scat. coef., HH pol., omnidirectional'
c        write(50,*) 'c: col 9: scat. coef., VV pol., Up/Down Wind asym'
c        write(50,*) 'c: col 10: same,       HH pol., Up/Down Wind asym'
c        write(50,*) 'c: col 11: same,       HV pol., Up/Down Wind asym'
c        write(50,*) 'c: col 12: same,       VH pol., Up/Down Wind asym'
c        write(50,*) 'c: col 13: same,      VV pol., Up/Cross Wind asym'
c        write(50,*) 'c: col 14: same,      HH pol., Up/Cross Wind asym'
c        write(50,*) 'c: col 15: same,      HV pol., Up/Cross Wind asym'
c        write(50,*) 'c: col 16: same,      VH pol., Up/Cross Wind asym'
c
c Echantillonnage en theta de 0 à THETAmax par pas de DTHETAtab

c Début BOUCLE THETA_1>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        do nTHETA_1 =  0, int ( THETAmax/DTHETAtab )
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        theta_1 = nTHETA_1*DTHETAtab
        if ((theta_1.ne.0.0D0).and.(theta_1.ne.90.0D0)) then
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


c > START LOOP: PHI_1
c Computation of 'Samples' to derive 0th and 2nd 
c harmonic coefficients for 1st and 2nd Stokes parameters.

c  	Harmonics 0 and 2 of Stokes 1 and 2
c
c       Irh(i,indice,nTHETA_1) is the i-th harmonic, 
c                              for the indice-th Stokes parameter,
c                              at theta = nTHETA_1*DTHETAtab
c       Ir(0,indice) est la Tb pour le Stokes "indice" en phi=TPHI_1(0)
c       Ir(1,indice) est la Tb pour le Stokes "indice" en phi=TPHI_1(1)      do nPHI_1 = 0,1
      do nPHI_1 = 0,1
		phi_1   = TPHI_1(nPHI_1)
		cophi1  = dcosd(phi_1)
		siphi1  = dsind(phi_1)
c  		initialisation
		do indice = 1,2
			Ir(nPHI_1, indice)  = 0.0D0
		enddo
c		computation
		call func2D1(krho1, dIr)
      do indice = 1, 2
      	Ir(nPHI_1,indice) = dIr(indice)*k02
     &                      *4*dacos(-1.0D0)/dtand(theta_1)*k0*cos1
c        print*, 'integ 3 : ',  dIr(indice), k02, 
c     &                          dIr(indice)*k02,
c     &                          dacos(-1.0D0), dtand(theta_1), k0, cos1

      enddo
      enddo
c < END LOOP: Phi_1
c Computation of 'Samples' to derive 2nd 
c harmonic coefficient for 3rd and 4th Stokes parameters.

      nPHI_1  = 2
      phi_1   = TPHI_1(nPHI_1)
      cophi1  = dcosd(phi_1)
      siphi1  = dsind(phi_1)
c  		initialisation
		do indice = 3, 4
			Ir(nPHI_1, indice)  = 0.0D0
		enddo
c		computation
		call func2D1(krho1, dIr)
      do indice = 3, 3
      	Ir(nPHI_1,indice) = dIr(indice-2)*k02
     &                      *4*dacos(-1.0D0)/dtand(theta_1)*k0*cos1
      enddo

      
c Derive the harmonic coefficients from samples


      do indice = 1, 2
        Irh(0,indice,nTHETA_1)=0.5D0*(Ir(0,indice)+Ir(1,indice))
        Irh(1,indice,nTHETA_1)=0.0D0
        Irh(2,indice,nTHETA_1)=0.5D0*(Ir(0,indice)-Ir(1,indice))
      enddo
      do indice = 3, 4
        Irh(0,indice,nTHETA_1) = 0
        Irh(1,indice,nTHETA_1) = 0
        Irh(2,indice,nTHETA_1) = 0.0D0 !Ir(2,indice)
      enddo
      else ! theta_1 = 90 deg
          do indice = 1, 4
            Irh(0,indice,nTHETA_1) = 0.0D0
            Irh(1,indice,nTHETA_1) = 0.0D0
            Irh(2,indice,nTHETA_1) = 0.0D0
          enddo
      endif

c Print coefficients
c        write(50,2000) nu, SST(iSST), SSS(iSSS), U10, ustar, 
c     &  theta_1                                               
c     &                 , Irh(0,1,nTHETA_1)
c     &                 , Irh(0,2,nTHETA_1)
c     &                 , Irh(1,1,nTHETA_1)
c     &                 , Irh(1,2,nTHETA_1)
c     &                 , Irh(1,3,nTHETA_1)
c     &                 , Irh(1,4,nTHETA_1)
c     &                 , Irh(2,1,nTHETA_1)
c     &                 , Irh(2,2,nTHETA_1)
c     &                 , Irh(2,3,nTHETA_1)
c     &                 , Irh(2,4,nTHETA_1)

        write (*,2000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
     &  theta_1                                               
     &                 , Irh(0,1,nTHETA_1)
     &                 , Irh(0,2,nTHETA_1)
     &                 , Irh(1,1,nTHETA_1)
     &                 , Irh(1,2,nTHETA_1)
     &                 , Irh(1,3,nTHETA_1)
     &                 , Irh(1,4,nTHETA_1)
     &                 , Irh(2,1,nTHETA_1)
     &                 , Irh(2,2,nTHETA_1)
     &                 , Irh(2,3,nTHETA_1)
     &                 , Irh(2,4,nTHETA_1)

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

	if (fGeo.ne.0) then ! seulement si on veut les grandes vaugues
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
c   Integration : PARTIE GEOMETRIQUE 
        
c Début boucle theta>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	do itheta =1, ntheta
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

		ind = theta(itheta)/1 ! on prend la partie entière de theta_1
		Poids = (theta(itheta)/1 -ind) ! poids des diff Tb
		AtmoU = TbAu(ind)*(1 - Poids) + TbAu(ind+1)*Poids
        
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
     	do iphi = 1, 4
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

        phi = TPHI(iphi)
     	cophi = dcosd(phi)
     	siphi = dsind(phi)
       
c  Calcul des bornes de l'integrale: en +- xVar*variance des pentes
c  avec pour condition Sx' < cotan(theta)


        Norm2 = 0.D0
        Sy_sup =  xSigma*Sc
       	Sy_inf = -xSigma*Sc
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
        

c                write(50, *) 'Sy, Sx, theta, phi, thetal, phil, P, intVV
c     &, intHH'
c Début  boucle Sy
       	do i=0, N1-1
       		Sy = Sy_inf + i*del1
                if (cophi.eq.0.D0) then
                        Sx_ = Sy*siphi
                        if (Sx_.le.C7) then
                                Sx_sup =  xSigma*Su
                                Sx_inf = -xSigma*Su
                        else
                                Sx_sup =  66666.D0
                                Sx_inf =  66666.D0
                        endif
                else
                        limit = (C7-Sy*siphi)/cophi
                        if (cophi.gt.0.D0) then
                                Sx_sup = min(limit, xSigma*Su)
c                                if (Sx_sup.eq.limit)
c     &                                  compt = compt + 1
                                Sx_sup = Sx_sup-abs(Sx_sup*1.0D-06)
                                Sx_inf = -xSigma*Su
                        else
                                Sx_sup =  xSigma*Su
                                Sx_inf =max(limit,-xSigma*Su)*1.000001D0
                                Sx_inf = Sx_inf+abs(Sx_inf*1.0D-06)
                        endif
                endif
                del2 = (Sx_sup-Sx_inf)/(N2-1)
                if ((Sx_sup.ne.66666.D0).and.(Sx_inf.ne.66666.D0)) then
c Début boucle Sx
          	do j=0, N2-1
          		Sx = Sx_inf + j*del2
                        Sx_ = Sx*cophi + Sy*siphi
c          		projec = 1.0D0-Sx_*C2
                projec = 1.0d0/dcos(datan(Sx))/dcos(datan(Sy))
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                        
c                        projec = 1.0D0
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                        

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
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!
c Version Valenzuela thetal (change also computation of cosalpha,
c sinalpha in accordance)
c        cosl = dcosd(theta(itheta)+datand(Sx))*dcos(datan(Sy))
c       	sinl = dsqrt(1.0D0-cosl*cosl) 
c        thetal = dacosd(cosl)
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c        phil = 0
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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

        call Glob2Loc(theta(itheta), phi, Sx, Sy, thetaAtmo, temp1,
     &                temp2, temp3)


c-----------------------------------------------------------------------
c               COEFFICIENT DE MODULATION HYDRODYNAMIQUE h
c
c  Sx : Pente de la vague dans la direction du vent
c  Su : Variance de la distribution des vagues dans la direction du vent
c  h  : Coefficient de modulation hydrodynamique
c-----------------------------------------------------------------------


	   h = 1.D0 - 0.20D0*Sx/Su

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
c                            COEFFICIENTS DE FRESNEL "EFFECTIFS" LOCAUX
c
c  Rvv0    : Coefficient de reflexion en polar V
c  Rhh0    : Coefficient de reflexion en polar H
c  Rveff   : Coefficient de reflexion effectif en polar V
c  Rheff   : Coefficient de reflexion effectif en polar V

        Rveff = Rvv0*dexp(-4.0D0*(sigma*cosl*k0)**2)
        Rheff = Rhh0*dexp(-4.0D0*(sigma*cosl*k0)**2)

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
                        Diff(3) = Irh(2,3,ind)*(1.0D0-Poids)
     &                            + Irh(2,3,ind+1)*Poids
                        Diff(4) = Irh(0,4,ind)*(1.0D0-Poids)
     &                            + Irh(0,4,ind+1)*Poids
c 2eme Harmonique
                        Diff(1) = Diff(1) + (Irh(2,1,ind)*(1.0D0-Poids)
     &                                    + Irh(2,1,ind+1)*Poids
     &                                      )*dcosd(2.0D0*phiL)
                        Diff(2) = Diff(2) + (Irh(2,2,ind)*(1.0D0-Poids)
     &                                    + Irh(2,2,ind+1)*Poids
     &                                      )*dcosd(2.0D0*phiL)
                        Diff(3) = Diff(3) !+ (Irh(2,3,ind)*(1.0D0-Poids)
     &                                    !+ Irh(2,3,ind+1)*Poids
     &                                    !  )*dsind(2.0D0*phiL)
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


c-----------------------------------------------------------------------
c                       VECTEUR DE STOKES GLOBAL
c

c  CONTRIBUTION DE L'ECUME ET DE LA DIFFUSION DE BRAGG AU VECTEUR LOCAL

c--- Vecteur de stokes dans le repère local

c                       sans ecume

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c                        h = 1.0D0
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c h: coefficient hydrodynamique                        
          		Isl(1) = h*Diff(1)
          		Isl(2) = h*Diff(2)
          		Isl(3) = h*Diff(3)
       			Isl(4) = h*Diff(4)

c--- Vecteur de stokes dans le repère global
c                       sans ecume


c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c Version Valenzuela pour alpha (change also computation of thetal
c in accordance)
c                cosalpha = dsind(theta(itheta)
c     &                +datand(Sx))*dcosd(datand(Sy))/sinl
c                sinalpha = dsind(datand(Sy))/sinl
c                Isl(3) = 0.0D0
c                Isl(4) = 0.0D0
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c projection Valenzuela 1978               
                Is(1) = Isl(1)*cosalpha**4 
     &                        + Isl(2)*sinalpha**4 
c     &                        - Isl(3)*(cosalpha*sinalpha)**2
          		Is(2) = Isl(1)*sinalpha**4 
     &                        + Isl(2)*cosalpha**4 
c     &                        - Isl(3)*(cosalpha*sinalpha)**2
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                
                
c          		Is(1) = Isl(1)*cosalpha**2 
c     &                        + Isl(2)*sinalpha**2 
c     &                        + Isl(3)*cosalpha*sinalpha 
c          		Is(2) = Isl(1)*sinalpha**2 
c     &                        + Isl(2)*cosalpha**2 
c     &                        - Isl(3)*cosalpha*sinalpha
                        Is(3) = Isl(3)*(cosalpha**2-sinalpha**2) -
     &                        (Isl(1)-Isl(2))*2.D0*sinalpha*cosalpha
                        Is(4) = Isl(4)

c---Integrand Tb
c                       sans ecume

c!!!!!!!!!!!!!!!!!!!!!!!!!!                      
c       On force tout a 0 quand SPM pas valide                        
                        if ((2.0D0*k0*sinl).lt.kd) then
                            temp1 = 0.0D0
                        else
                            temp1 = 1.0D0
                        endif
                        do indice = 1,4
                         Is(indice) = Is(indice)*temp1
                        enddo
c!!!!!!!!!!!!!!!!!!!!!!!!!                         

          		int1(1,iphi) = int1(1,iphi) + 
     &                                 Is(1)*projec*P_Sx_Sy*del1*del2
          		int1(2,iphi) = int1(2,iphi) + 
     &                                 Is(2)*projec*P_Sx_Sy*del1*del2
          		int1(3,iphi) = int1(3,iphi) + 
     &                                 Is(3)*projec*P_Sx_Sy*del1*del2
          		int1(4,iphi) = int1(4,iphi) + 
     &                                 Is(4)*projec*P_Sx_Sy*del1*del2

c!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c                write(40, '(11(1x,e17.9))')  Sy, Sx, theta(itheta), phi,
c     &thetal, phil, P_Sx_Sy, Is(1)*P_Sx_Sy, Is(2)*P_Sx_Sy, cosalpha,
c     &sinalpha
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c  NORMALISATION DE LA DENSITE DE PROBABILITE : Norm2
c                        Norm2 = Norm2 + projec*P_Sx_Sy*del1*del2
                Norm2 = 1.0D0
                enddo
                endif
                enddo
c               sans ecume
                int1(1,iphi) = int1(1,iphi)/Norm2
                int1(2,iphi) = int1(2,iphi)/Norm2
                int1(3,iphi) = int1(3,iphi)/Norm2
                int1(4,iphi) = int1(4,iphi)/Norm2
c Fin boucle phi<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
      		enddo
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

c----------------------------------------------------------------------
c                       CALCUL DES HARMONIQUES SANS ECUME
c
c int1(1,1) : Tb verticale en phi = TPHI(1) = 0°
c int1(1,2) : Tb verticale en phi = TPHI(2) = 45°
c int1(1,3) : Tb verticale en phi = TPHI(3) = 90°
c int1(1,4) : Tb verticale en phi = TPHI(4) = 180°
c Le premier indice de int1 defini l'indice du parametre de stokes et 
c varie de 1 à 4 respectivement pour Tv, Th, U et V.


        sigmaVV0 = 0.25D0*(int1(1,1)+2.0D0*int1(1,3)+int1(1,4))
        sigmaHH0 = 0.25D0*(int1(2,1)+2.0D0*int1(2,3)+int1(2,4))
        sigmaVV1 = 0.5D0 *(int1(1,1) - int1(1,4)) 
        sigmaHH1 = 0.5D0 *(int1(2,1) - int1(2,4))
        sigmaHV1  = int1(3,3) 
        sigmaVH1  = int1(4,3)
        sigmaVV2 = 0.25D0*(int1(1,1)-2.0D0*int1(1,3)+int1(1,4))
        sigmaHH2 = 0.25D0*(int1(2,1)-2.0D0*int1(2,3)+int1(2,4))
        sigmaHV2  = int1(3,2) - int1(3,3)
        sigmaVH2  = int1(4,2) - int1(4,3)

        sigmaQS = Reff*dexp(-dtand(theta(itheta))**2/2/Su**2)/2/Su/Sc
     &              /dcosd(theta(itheta))**4        
        if (sigmaQS.le.1.0D-100) then
           sigmaQS = 0.0D0
        endif 
		write(*,*) 'Sans Atmo'

c   COMPUTE AND WRITE BRAGG SCATTERING
c!!!!        phiL = 0.d0

c SI ON RESTE DANS LE DOMAINE TABULÉ
                        if (theta(itheta).le.THETAmax) then
                                ind     = theta(itheta)/DTHETAtab
                                Poids   = (theta(itheta)/DTHETAtab-ind) 
c Fondamental
                        Diff(1) = Irh(0,1,ind)*(1.0D0-Poids)
     &                            + Irh(0,1,ind+1)*Poids
                        Diff(2) = Irh(0,2,ind)*(1.0D0-Poids)
     &                            + Irh(0,2,ind+1)*Poids
                        Diff(3) = Irh(2,3,ind)*(1.0D0-Poids)
     &                            + Irh(2,3,ind+1)*Poids
                        Diff(4) = Irh(0,4,ind)*(1.0D0-Poids)
     &                            + Irh(0,4,ind+1)*Poids
c 2eme Harmonique
                        Diff(1) = Diff(1) + (Irh(2,1,ind)*(1.0D0-Poids)
     &                                    + Irh(2,1,ind+1)*Poids
     &                                      )*dcosd(2.0D0*phiL)
                        Diff(2) = Diff(2) + (Irh(2,2,ind)*(1.0D0-Poids)
     &                                    + Irh(2,2,ind+1)*Poids
     &                                      )*dcosd(2.0D0*phiL)
                        Diff(3) = Diff(3) !+ (Irh(2,3,ind)*(1.0D0-Poids)
     &                                    !+ Irh(2,3,ind+1)*Poids
     &                                    !  )*dsind(2.0D0*phiL)
                        Diff(4) = Diff(4) + (Irh(2,4,ind)*(1.0D0-Poids)
     &                                    + Irh(2,4,ind+1)*Poids
     &                                      )*dsind(2.0D0*phiL)
c       print*, Diff(1),  Irh(0,1,ind), Irh(0,1,ind+1), Irh(2,1,ind)
c     &             , Irh(2,1,ind+1), dcosd(2.0D0*phil), phil
c SI ON SORT DU DOMAINE TABULÉ À LA PRECISION NUMÉRIQUE PRÈS
                        elseif ((P_Sx_Sy.eq.0.D0).or.
     &                     (theta(itheta).le.THETAmax*1.0001D0)) then
                                do indice = 1, 4
                                        Diff(indice) = 0.0D0
                                enddo
                        else
c SI ON SORT DU DOMAINE, ERREUR => SORTIE DU PROGRAMME
                                print *, 'ERREUR !!!!!!!!!!!!!'
                                print *, 'theta(itheta) > THETAmax'
                                print *, theta(itheta), ' > ', THETAmax
                                print *, 'Densité de Proba = ',P_Sx_Sy
                                stop
                        endif

c Compute NRCS for bragg scattering: interpolation of |Diff| at |theta(itheta)|
c Fundamental only
                        if (theta(itheta).le.THETAmax) then
                                ind     = theta(itheta)/DTHETAtab
                                Poids   = (theta(itheta)/DTHETAtab-ind) 
                             
                        temp1 = Irh(0,1,ind)*(1.0D0-Poids)
     &                            + Irh(0,1,ind+1)*Poids
                        temp2 = Irh(0,2,ind)*(1.0D0-Poids)
     &                            + Irh(0,2,ind+1)*Poids
                        else
                            temp1 = 999.99
                            temp2 = 999.99
                        endif

        
        write(20,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
     &                  theta(itheta), sigmaVV0, sigmaHH0,
     &                  sigmaVV1, sigmaHH1, sigmaHV1, sigmaVH1,
     &                  sigmaVV2, sigmaHH2, sigmaHV2, sigmaVH2, sigmaQS,
     &                   Diff(1), Diff(2), Diff(3), Diff(4),
     &                 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,  k0/kd,
     &                 dexp(-4.0D0*(sigma*k0)**2), realpart(epsi),
     &                 imagpart(epsi), temp1, temp2, Fr*100.0D0  
        write (*,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
     &                 theta(itheta), sigmaVV0, sigmaHH0, 
     &                 sigmaVV1, sigmaHH1, sigmaHV1, sigmaVH1,
     &                 sigmaVV2, sigmaHH2, sigmaHV2, sigmaVH2, sigmaQS,
     &                   Diff(1), Diff(2), Diff(3), Diff(4),
     &                 0.D0, 0.D0, 0.D0, 0.D0, 0.D0, 0.D0,  k0/kd,
     &                 dexp(-4.0D0*(sigma*k0)**2), realpart(epsi),
     &                 imagpart(epsi), temp1, temp2, Fr*100.0D0
c Fin boucle Theta<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        enddo
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	endif ! fin condition sur fGeo
c  Fin intégrale sur les grandes vagues

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

        close (10)
	close (15)
        close (20)
c        close (25)
        close (30)
        close (40)
        close (50)

  100     format (1x,f5.2,4(1x,f4.1),2(1x,f6.2),2(1x,f6.2)
     &          ,8(1x,f5.2))
 1000     format (1x,f9.5,1x,f5.2,1x,f6.3,2(1x,f6.2),1x,f4.1,
     &         28(1x,e12.5))
 2000     format (1x,f9.5,1x,f5.2,1x,f6.3,2(1x,f6.2),1x,f4.1,
     &         10(1x,e12.5))
       
        stop
   60      print*, 'Unable to read ', fin1
        stop
   61      print*, 'Unable to open filterF.dat'
        stop
   62      print*, 'Unable to open filterFres.dat'
        stop
   63      print*, 'Unable to open filterPhi.dat'
        stop        
   70      print*, 'Unable to create Output/integPts.dat'
        stop
c   70      print*, 'Unable to create Output/Diffusion.dat'
c        stop
   80      print*, 'Unable to create Output/rapport.dat'
        stop


        end

