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

c-------------------------------------------------------------------------------------
c   Modifications as of 10 Nov 2021 
c    by Magdalena Anguelova, NRL, Washington, DC, USA
C	* Frequency band (nuBand = L,C,X,K,Ka,W) is input from the config file Tb.p (Line 461)
c	* Introduce case 7 for foam emissivity tuned for foam parameters (line 1757)
c	* Call to subroutine for tuned foam emissivity (line 1759)
c	* Stability parameter for foam coverage M1 (Monahan, 1986, eq.5) with minus sign "-" (Line 1278)
c		as atm. stability from the config file Tb.p is Tair-Twater
c	* Added foam + atmopshere calculations (only foam without atmopshere in the original) (Lines 2082-2104)
c-------------------------------------------------------------------------------------


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
        character*200 logFile ! filename for program execution report
        character*1 cVar
        character*1  cepsi
        character*1  csigneeps
        character*1  cSpectre
        character*1  cType
        character*16 cCouvEcume
        character*16  cEmisEcume
        character*16 nuBand
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
        character*8 time
	
        double precision paramLem(10)
        double precision ss, Pi, beta
        double precision kmin
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
      double precision paramKudryavtsev(10)
      double precision filtersKudry(1:3,1:50003)
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
        double precision ustar_don
        double precision ustar_kudry
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
        double precision Is_GO(1:4) 
        double precision Is_e(1:4) 
        double precision IsA(1:4)  ! Is avec atmosphere 
        double precision Is_eA(1:4) ! Is_e avec atmosphere 
        double precision Isl(1:4) 
        double precision Isl_GO(1:4) 
        double precision Isl_e(1:2) 
        double precision IslA(1:4)  ! Isl avec atmosphere
        double precision Isl_eA(1:2) ! Isl_e avec atmosphere 
        double precision int1(1:4,1:4)
        double precision int1_GO(1:4,1:4)
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
        double precision Tv0_GO
        double precision Th0_GO
        double precision Tv1_GO
        double precision Th1_GO
        double precision T31_GO 
        double precision T41_GO
        double precision Tv2_GO
        double precision Th2_GO
        double precision T32_GO 
        double precision T42_GO
         double precision Tv0_SPM
        double precision Th0_SPM
        double precision T30_SPM
        double precision T40_SPM
        double precision Tv2_SPM
        double precision Th2_SPM
        double precision T32_SPM
        double precision T42_SPM
        double precision Tvn
        double precision Thn
     	double precision Rveff
     	double precision Rheff
        double precision xVar
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
       	double complex epsi_hifreq
       	double complex epsi
       	double complex epsi_
        double complex C1
     	double complex D1
        double complex D2
        double complex Rvv0
        double complex Rhh0
       	
        integer i 
        integer iWind
! TO BE REMOVED !!! (use for continuation of computations)
        integer iSSS0
        integer iSST0
        integer iWind0
        integer itheta0
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
        if (.not. fin1exist) then               ! if no parameter ...
          call getenv ('inputFortran', fin1)    ! get default path for
                                                ! input files. 
          if (fin1.eq.'') then                  ! if no default path ...
            print*, 'No valid input file for parameters'
            print*, 'Moreover, environment variable ''inputFortran''
     &               for default path not defined.'
            print*, 'Program aborted.'          !  => error, end of prog
            stop
          else
            fin1 = fin1(1:lnblnk(fin1))//"/Tb.p" ! if default path exist, 
                               ! then append it to default input filename
                               ! (check of existance of file performed
                               ! when opening file)
            print*, 'No valid filename for input file.'
            print*, 'Use (if exists) default file : '
            print*, fin1 
          endif
        endif
c Open file of input parameters        
        open (unit=30, file=fin1, status='old', err=60)
c Identify path of output files (data output file is read in the input
c file)       
        call getenv ('outputFortran', pathOut)
        if (pathOut.eq.'') then
            print*, 'Environnement  variable ''outputFortran'' for
     & path of output files not defined.'
            print*, 'Program aborted.'
            stop
        endif
        pathout = pathout(1:lnblnk(pathout))//'/'
c Open output file for small scale reflectivity
        open (unit=60, file=pathout(1:lnblnk(pathout))//"Diffusion.dat"
     &,status='unknown',err=70)
c Open output file for Tb Atmo
        open (unit=40, file=pathout(1:lnblnk(pathout))//"TbAtmo.dat"
     &,status='unknown',err=70)
c Open file for program execution report (i.e. log file !)
      logFile = pathout(1:lnblnk(pathout))//"Tb_"//date//"_"//time//
     &".log"
        open (unit=50, file=logFile, status='unknown', err=80)

c----------------------------------------------------------------------
c------------------ READ INPUT PARAMETERS -----------------------------
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
     	read (30,'(a)') jump_line
    	read (30,*) nuBand					! frequency Band
     	read (30,'(a)') jump_line
     	read (30,*) lambda
       	freq = 3.D08/lambda
       	nu = 3.D-01/lambda
       	k0 = 2.D0*3.141592653D0/lambda 
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
        if ((fout1.eq.'aucun').or.(fout1.eq.'none')) then ! test on 'aucun' kept for legacy
                print *, ' * No output file with foam selected.'
                if ((fout2.eq.'aucun').or.(fout2.eq.'none')) then
                      print *, '  No output file without foam selected.'
                      print *, '/!\ WARNING /!\'
                      print *, 'No output file selected. Results only'//
     &                ' reported on screen'
                else
                   open (unit=20,file=pathout(1:lnblnk(pathout))//fout2
     &,status='unknown')
! TO BE REMOVED !!! (use for continuation of computations)
!                   open (unit=20,file=pathout(1:lnblnk(pathout))//fout2
!     &,status='unknown', access='append')
                   open (unit=25,file=pathout(1:lnblnk(pathout))//"Atm_"
     &//fout2,status='unknown')
                        nsorties = 1
                        sortiestab(1) = 20
                   print *, ' * Results without foam are saved in file :
     & '//pathout(1:lnblnk(pathout))//fout2
                endif
        else
                open (unit=10,file=pathout(1:lnblnk(pathout))//fout1,sta
     &tus='unknown')
                open (unit=15,file=pathout(1:lnblnk(pathout))//"Atm_"//
     &fout1,status='unknown')
                if ((fout2.eq.'aucun').or.(fout2.eq.'none')) then
                   print *, ' * Results with foam are saved in file : '
     & //pathout(1:lnblnk(pathout))//fout1
                   print *, ' * No output file without foam selected.'
                        nsorties = 1
                        sortiestab(1) = 10
                else
                   open (unit=20,file=pathout(1:lnblnk(pathout))//fout2,
     &status='unknown')
                   open (unit=25,file=pathout(1:lnblnk(pathout))//"Atm_"
     &//fout2,status='unknown')
                        nsorties = 2
                        sortiestab(1) = 10
                        sortiestab(2) = 20
                   print *, ' * Results with foam are saved in file : '
     & //pathout(1:lnblnk(pathout))//fout1
                 print *, ' * Results without foam are saved in file : '
     & //pathout(1:lnblnk(pathout))//fout2
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
        if (DTHETAtab.lt.1.0D0) then
                write (*,*)
                print *, '!!!!!! ERROR  !!!!!!!!!!!'
                print *, 'DTHETAtab must be larger or equal to 1'
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
        


c----------- VERIFICATION DES PARAMETRES D'ENTRÉE ----------------------

        if (int(THETAmax/DTHETAtab).ne.(THETAmax/DTHETAtab)) then
                write (*,*)
                print *, '!!!!!! ERROR  !!!!!!!!!!!'
                print *, 'THETAmax is not a multiple of DTHETAtab'
                write (*,*) 'THETAmax = ', THETAmax
                write (*,*) 'DTHETAtab = ', DTHETAtab
                write (*,*)
                stop
        endif
        if ((theta_max.gt.THETAmax).and.(nTheta.ne.1)) then
                write (*,*)
                print *, '!!!!!!!! ERROR  !!!!!!!!!!!!!!!!'
                print *, 'theta must be less or equal to ', THETAmax
     &          ,' which is THETAmax (max value for lookup table)'
                write (*,*)
                stop
        endif
        if ((cSpectre.eq.'D').or.(cSpectre.eq.'d')) then
                cspec = 1
                ModSpectre = 'Durden and Vesecky (1985)'
        elseif ((cSpectre.eq.'E').or.(cSpectre.eq.'e')) then
                cspec = 2
                ModSpectre = 'Elfouhaily et al. (1997)'
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
                print *, '!!!!!!!! ERROR  !!!!!!!!!!!!!!!!'
                print *, 'Current choice for sea spectrum model is:', 
     &                  cSpectre
                print*,'Possible choices are:'
                print*,'   - E for Elfouhaily et al. (1997)'
                print*,'   - D for Durden and Vesecky (1985)'
                print*,'   - L for Lemaire et al. (1999)'
                print*,'   - K for Kudryavtsev et al. (2003)'
                write (*,*)
                stop
        endif
        if (cType.eq.'2') then
                TypeMod = 'Two-scale'
                fGeo  = 1
                fDiff = 1 
        elseif ((cType.eq.'p').or.(cType.eq.'P')) then
                TypeMod = 'Small scales (Scattering)'
                fGeo  = 0
                fDiff = 1
        elseif ((cType.eq.'g').or.(cType.eq.'G')) then
                TypeMod = 'Large scales (Geometric Optics)'
                fGeo  = 1
                fDiff = 0
        else
                write (*,*)
                print *, '!!!!!!!! ERROR  !!!!!!!!!!!!!!!!'
                print *, 'Current choice for the type of model is: '
     &                          , cType
                print *, 'Possible Choices are:'
                print*,  '  - 2 for  two-scale model'
                print *, '  - p for small scales scattering only'
                print *, '  - g for large scales only'
                write (*,*)
                stop
        endif
        if ((cepsi.eq.'K').or.(cepsi.eq.'k')) then
                Modepsi = 'Klein & Swift (77)'
        elseif ((cepsi.eq.'e').or.(cepsi.eq.'E')) then
                Modepsi = 'Ellison (98)'
        elseif ((cepsi.eq.'m').or.(cepsi.eq.'M')) then
                Modepsi = 'Meissner et al. (2004,2012,2014)'
        elseif ((cepsi.eq.'h').or.(cepsi.eq.'H')) then
                Modepsi = 'High frequency dataset'
        else
                write(*,*)
                print *,'Invalid choice for dielectric constant model:'
                print *,'     K : Klein & Swift (77)'
                print *,'     E : Ellison (98)'
                print *,'     M : Meissner et al. (2004,2012,2014)'
                print *,'     H : High frequency dataset'
                print *,'Current choice is: ',cepsi
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
                print *,'Invalid choice for the slope variance '//
     &                  'computation. Possible choices are:'
                print *,'   - C : Cox and Munk (1955)'
                print *,'   - S : Compute from sea spectrum model'
                print *,'Current choise is: ',cVar
                write(*,*)
                stop
        endif
        if ((cCouvEcume.eq.'A').or.(cCouvEcume.eq.'a')) then
                ModCouvEcume = 'No foam.'
                fCouvEcume = 0
        elseif ((cCouvEcume.eq.'M1').or.(cCouvEcume.eq.'m1')) then
                ModCouvEcume = 'Monahan (1986, eq.5)'
                fCouvEcume = 1
        elseif ((cCouvEcume.eq.'M2').or.(cCouvEcume.eq.'m2')) then
                ModCouvEcume = 'Monahan least square (1986, eq. 3a)'
                fCouvEcume = 2
        elseif ((cCouvEcume.eq.'M3').or.(cCouvEcume.eq.'m3')) then
                ModCouvEcume = 'Monahan and Lu (90, active + passive)'
                fCouvEcume = 3
        elseif ((cCouvEcume.eq.'M4').or.(cCouvEcume.eq.'m4')) then
                ModCouvEcume = 'Monahan and Lu (90, active)'
                fCouvEcume = 4
        elseif ((cCouvEcume.eq.'M5').or.(cCouvEcume.eq.'m5')) then
                ModCouvEcume = 'Data from WISE 2001 campaign'
                fCouvEcume = 5
        elseif ((cCouvEcume.eq.'M-Du-E')) then
                ModCouvEcume = 'Yin et al. 2016 - M-Du-E'
                fCouvEcume = 6
        elseif ((cCouvEcume.eq.'M-Du-E1')) then
                ModCouvEcume = 'Yin et al. 2016 - M-Du-E1'
                fCouvEcume = 7
        elseif ((cCouvEcume.eq.'M-Du-S')) then
                ModCouvEcume = 'Yin et al. 2016 - M-Du-S'
                fCouvEcume = 8
        elseif ((cCouvEcume.eq.'M-Ku-E')) then
                ModCouvEcume = 'Yin et al. 2016 - M-Ku-E'
                fCouvEcume = 9
        elseif ((cCouvEcume.eq.'M-Ku-S')) then
                ModCouvEcume = 'Yin et al. 2016 - M-Ku-S'
                fCouvEcume = 10
        else
                write(*,*)
                print *,'Invalid choice for foam fraction model.'//
     &                  ' Possible choices are:'
                print *,'  A  : No foam'
                print *,'  M1 : Monahan (1986, eq.5)'
                print *,'  M2 : Monahan least square (1986, eq. 3a)'
                print *,'  M3 : Monahan and Lu (1990, active + passive)'
                print *,'  M4 : Monahan and Lu (1990, active)'
                print *,'  M-Du-E : Yin et al. (2016) M-Du-E version'
              print *,'    M-Du-E1 : Yin et al. (2016) M-Du-E1 version'
              print *,'    M-Du-S : Yin et al. (2016) M-Du-S version'
              print *,'    M-Ku-E : Yin et al. (2016) M-Ku-E version'
              print *,'    M-Ku-S : Yin et al. (2016) M-Ku-S version'
                print *,'Current choice is: ',cCouvEcume
                write(*,*)
                stop
        endif
        if ((cEmisEcume.eq.'S').or.(cEmisEcume.eq.'s')) then
                ModEmisEcume = 'Stogryn (72)'
                fEmisEcume = 1
        elseif ((cEmisEcume.eq.'M-Du-E')) then
                ModEmisEcume = 'Yin et al. 2016 - M-Du-E'
                fEmisEcume = 2
        elseif ((cEmisEcume.eq.'M-Du-E1')) then
                ModEmisEcume = 'Yin et al. 2016 - M-Du-E1'
                fEmisEcume = 3
        elseif ((cEmisEcume.eq.'M-Du-S')) then
                ModEmisEcume = 'Yin et al. 2016 - M-Du-S'
                fEmisEcume = 4
        elseif ((cEmisEcume.eq.'M-Ku-E')) then
                ModEmisEcume = 'Yin et al. 2016 - M-Ku-E'
                fEmisEcume = 5
        elseif ((cEmisEcume.eq.'M-Ku-S')) then
                ModEmisEcume = 'Yin et al. 2016 - M-Ku-S'
                fEmisEcume = 6
	elseif ((cEmisEcume.eq.'M-Du-Tune')) then
                ModEmisEcume = 'Yin et al. 2016 - M-Du-Tune'
                fEmisEcume = 7	
        else
                write(*,*)
                print *,'Invalid choice for foam emissivity model.'//
     &                  ' Possible choices are:'
                print *,'     S : Stogryn (1972)'
                print *,'     M-Du-E : Yin et al. (2016) M-Du-E version'
              print *,'     M-Du-E1 : Yin et al. (2016) M-Du-E1 version'
                print *,'     M-Du-S : Yin et al. (2016) M-Du-S version'
                print *,'     M-Ku-E : Yin et al. (2016) M-Ku-E version'
                print *,'     M-Ku-S : Yin et al. (2016) M-Ku-S version'
                print *,'  M-Du-tune : Yin et al. (2016) tuned version'
                print *,'Current choice is: ',cEmisEcume
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
                print *,'Invalid choice for the drag coefficient'//
     &                  ' model. Possible choices are:'
                print *,'     y : Cardone (1969)'
                print *,'     c : Charnock (1955)' 
                print *,'     d : Donelan et al. (1993)' 
                print *,'Current choice is: ',cCD
                write(*,*)
                stop
        endif
        if ((cSwell.eq.'O').or.(cSwell.eq.'o')) then
                fSwell = 1
        elseif ((cSwell.eq.'N').or.(cSwell.eq.'n')) then
                fSwell = 0
        else
                write(*,*)
                print *,'Invalid choice for the inclusion of swell'//
     &                  'Possible choice are:'
                print *,'     O : swell is included'
                print *,'     N : No swell included' 
                print *,'Current choice is: ',cSwell
                write(*,*)
                stop
        endif

c---------------- AFFICHAGE DES VARIABLES DANS LE FICHIER --------------
        if (nsorties.ge.1) then

        do i = 1, nsorties
        ufile = sortiestab(i)
        write (ufile,*) 'Date = '//date//' ; Time = '//time(1:8)//' ;'
       	write (ufile,*) 'Lambda   = ',lambda,' ;'
       	if (kdmin.eq.0.)  then
                write (*,*) ' !!!! Error !!!!'
                write (*,*) 'kdmin = 0 !!!!!'
                stop
        endif
        write (ufile,*) 'Lambda_d min = ', PI2/kdmax,
     &               ' ; Lambda_d max = ', PI2/kdmin,
     &               ' ; Nlambda_d = ', nkd,' ;'


        write (ufile,*) 'Wind min.        = ', Windmin,
     &               ' ; Wind max.        = ', Windmax, 
     &               ' ; N_wind         = ', nWind, ' ;'
        write (ufile,*) 'Incidence angle min.   = ', thetamin, 
     &               ' ; Incidence angle max.   =' , theta_max,
     &               ' ; Nincidences   = ', ntheta, ' ;'
        write (ufile,*) 'Temperature min. = ', SSTmin, 
     &               ' ; Temperature max. =' ,SSTmax, 
     &               ' ; Ntemperatures = ', nSST, ' ;'
        write (ufile,*) 'Stability min.   = ', Stabmin,
     &               ' ; Stability max.   =' ,Stabmax,
     &               ' ; Nstabilities   = ', nStab, ' ;'
        write (ufile,*) 'Salinity min.    = ', SSSmin,
     &               ' ; Salinity max.    =' ,SSSmax,
     &               ' ; Nsalinities    = ', nSSS, ' ;'
        write (ufile,*) 'Model for dielectric constant of sea water = '
     &                  , Modepsi,' ;'
        write (ufile,*) 'Computation of slope variances = ',
     &               ModVarPente,' ;'
        write (ufile,*) 'Type of model  = ', TypeMod,' ;'
        write (ufile,*) 'Model for sea spectrum  = ',
     &                 ModSpectre,' ;'
        write (ufile,*) 'Number of integraton points on gaussian pdf'//
     &                  ' along x-axis = ', N1, ' ;' 
        write (ufile,*) 'Number of integraton points on gaussian pdf'//
     &                  ' along y-axis = ', N2, ' ;' 
        if (ufile.eq.10) then
                write (ufile,*) 'Model for foam fraction = ', 
     &                  ModCouvEcume, ' ;'
                write (ufile,*) 'Model for foam emissivity = ', 
     &                  ModEmisEcume, ' ;'
        else if (ufile.eq.20) then
                write (ufile,*) 'Model for foam fraction = ', 
     &                 'No foam ;'
                write (ufile,*) 'Model for foam emissivity = ', 
     &                 'No foam ;'
        else
                print *, 'Error, ambiguous output file.'
                stop
        endif
        write (ufile,*) 'Model for drag coefficient = ',ModCD,' ;'
        write (ufile,*) '! freq        SST   SSS     U10   ustar'
     &  //' theta  ' 
     &  //' TvN    ThN     Tv0      Th0    Tv1    Th1    U1      V1    '
     &  //' Tv2   Th2    U2      V2    lam_d    Stab Re(epsi)I(epsi)'
     &  //' Var_u    Var_c      Foam;'
        write (ufile,*) '! dielectric constant: e = Re(e) + i Im(e) $'

        enddo
        endif ! if (nsorties.ge.1)
        
c---------------------------------------------------------------------
c------------------ PROGRAMME PRINCIPAL ------------------------------

c  Calcul des Variance de pentes de la houle

        write(50,*) 'Computation of slope variances for swell ...'
        if (fSwell.eq.1) then
                call infSwell(NSwell, hSwell, sigSwell , 
     &          KMaxSwell, VarSwell)
        else
                VarSwell(1) = 0.0D0
                VarSwell(2) = 0.0D0
        endif
        write(50,*) ' ... Done.'
        write (50,*) 'SuSwell = ', sqrt(VarSwell(1))
        write (50,*) 'ScSwell = ', sqrt(VarSwell(2))

c Début Boucle sur SSS

! TO BE REMOVED !!! (use for continuation of computations)
!      iSSS0 = 4
!      iSST0 = 3
!      iWind0 = 1
!      iTheta0 = 13
      iSSS0 = 1
      iSST0 = 1
      iWind0 = 1
      iTheta0 = 1
        do iSSS=iSSS0, nSSS
       	   
c Début Boucle sur SST

       	do iSST = iSST0, nSST
       	   
       	   Tw = 2.7315D2 + SST(iSST)
           
c----------------------------------------------------------------------
c                       CONSTANTE DIELECTRIQUE : epsi
c
c epsi : constante dielectrique complexe relative à epsilon 0 
c        (ie permitivité du vide)
c       epsi_KS pour Klein & Swift 
c       epsi_El pour Ellison
c       epsi_MW pour Meissner and Wentz 
c       espi_hifreq [high frequency dataset]
c SST  : Temperature de surface de l'ocean en °C
c SSS  : Salinité de surface de l'ocean en ppm
c freq : Frequence electromagnetique du radiometre en Hz
        
        write (50,*) 'Computation of the sea water dielectric'//
     &               ' constant ...'
        call epsilon_KS (SST(iSST), SSS(iSSS), epsi_KS, freq)
        call epsilon_El (SST(iSST), SSS(iSSS), epsi_El, freq)
        call Epsilon_MW (SST(iSST), SSS(iSSS), epsi_MW, freq)
        call epsilon_hifreq (SST(iSST), SSS(iSSS), epsi_hifreq, freq)
        if ((cepsi.eq.'K').or.(cepsi.eq.'k')) then
                epsi = epsi_KS
                Modepsi = 'Klein & Swift (77)'
        elseif ((cepsi.eq.'e').or.(cepsi.eq.'E')) then
                epsi = epsi_El
                Modepsi = 'Ellison (98)'
        elseif ((cepsi.eq.'m').or.(cepsi.eq.'M')) then
                epsi = epsi_MW
                Modepsi = 'Meissner et al. (2004,2012,2014)'
        elseif ((cepsi.eq.'h').or.(cepsi.eq.'H')) then
                epsi = epsi_hifreq
                Modepsi = 'High frequency dataset'
        else
                write(*,*)
                print *,'Invalid choice for dielectric constant model.'
                print *,'Possible choices are :'
                print *,'     K : Klein and Swift (1977)'
                print *,'     E : Ellison et al. (1998)'
                print *,'     M : Meissner et al. (2004,2012,2014)'
                print *,'     H : High frequency dataset'
                print *,'Current choice is: ',cepsi
                write(*,*)
                stop
        endif
        if ((csigneeps.eq.'o').or.(csigneeps.eq.'O')) then
                continue
        elseif ((csigneeps.eq.'n').or.(csigneeps.eq.'N')) then
                epsi = conjg(epsi)
        else
                write(*,*)
                print *,'Invalid choice for the sign of the imaginary'//
     &                   ' part for the dielectric constant.'
                print *, 'Possible choices are:'
                print *, '     O : positive imaginary part'
                print *, '     N : negative imaginary part'
                print *, 'Current choice is: ',csigneeps
                write(*,*)
                stop
        endif

        !epsi = (1.0D0, -1.0D16)
        write (50,*) '   ... Done.'
        write (50,*) '  Klein & Swift = ', epsi_KS
        write (50,*) '  Ellison       = ', epsi_El
        write (50,*) '  Meissner et al. = ', epsi_MW
        write (50,*) '  Selected model    =>       ',epsi

c Tabulation des Tb Atmosphérique en fonction de l'angle d'incidence 
c theta (theta de 0 à 90° par pas de 1°)
        
        nParam = 8
        Param(1) = 6370.0D0 ! rayon terrestre (km)
        Param(2) = 1013.0D0 ! pression atmosphérique au sol (mb)
        Param(3) = SST(iSST) + 273.15D0  ! Température au sol (K)
        Param(4) = 70.0D0  ! Humidité relative (%)
        Param(5) = freq  ! fréquence (Hz)
        Param(6) = 20.0D0  ! altitude maximum de l'atmosphère (km)
        Param(7) = 400  ! nb de couches d'atmosphere
        Param(8) = 6.5D0  ! module gradient de temperature (K/km)
        write (50,*) 'Tabulation of atmospheric brightness tem'//
     &'peratures ...'
        call TbAtmo(Param, nParam, TbAd, TbAu, tau) 
        ! TbAd Tb atmosphere vers le bas (downward)
        ! TbAu Tb atmosphere vers le haut (upward)
c        do i = 0, 90
c                tau(i) = 0.0D0
c                TbAd(i) = 0.0D0
c                TbAu(i) = 0.0D0
c        enddo
c	Write Header for TbAtmo file
        write (40,*) 'P(0,mb), T(0,K) Hr(0,%), freq (Hz), grad(K/km), the
     &ta (deg), TbUp(k), TbDown(K), tau(neper)' 
        write (40,*) '$'   
        
        write(*,*) ''
        write(*,*) '----------   Atmo Tb Lookup Table  ----------------'
        write(*,*) '  incidence  Tb Up (K)  Tb down (K)    attenuation'
        write(*,*) '---------------------------------------------------'
	do i = 0, 90
        write (*,'(5x,f4.1,5x,f7.3,5x,f7.3,9x,f7.5)') i*1.0D0, TbAu(i),
     &        TbAd(i), tau(i)
        write (40,'(1x,f7.2,1x,f6.2,1x,f5.2,1x,e11.4,1x,f5.2,3(1x,f5.2),
     &e11.4)')
     &  Param(2), Param(3),  Param(4),Param(5),Param(8),
     &   i*1.0D0, TbAu(i), TbAd(i), tau(i)
	enddo
        write(*,*) '---------------------------------------------------'
        close (40) ! close Tb Atmo file
        write (50,*) ' ... Done.'


c----------------------------------------------------------------------
c Début Boucle sur Vent


       	do iWind= iWind0, nWind
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

         write (50,*) 'Computation of wind stress for wind speed U(',alt
     &       ,') = ', Uz
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
                print *, 'Invalid choice for drag coefficient model.'
                print *, 'Your choice: ', cCD
                print *, 'Possible choices are:'
                print*, '   y , Cardone (1969)'
                print *,'   c , Charnock (1955)'
                print *,'   d , Donelan et al. (1993)'
                write (*,*)
                stop
        endif
        if (ustar.eq.1.0D04) print *, '/!\ Change boundary values '//
     &   ' used in root function to compute wind stress u*, '//
     &   'computation failed. /!\'
         write (50,*) '   u* = ', ustar

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

         write (50,*) 'Conversion of wind speed at different'//
     &       '  reference heights ...'

       	if (Uz.ne.0.) then
                call vent (alt,1.95D01,Uz,U19,ustar)
                call vent (alt,1.25D01,Uz,U12,ustar)
                call vent (alt,1.0D01 ,Uz,U10,ustar)
       	else 
                U19 = 0.
                U12 = 0.
       	endif
         
         write (50,*) '   ... Done.'
         write (50,*) '   U10 = ', U10, ' m/s'
         write (50,*) '   U12.5 = ', U12, ' m/s'
         write (50,*) '   U19.5 = ', U19, ' m/s'
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
        write (50,*) 'Computation of lower boundary wavenumber for sea s
     &pectrum ...'
        call root (1.0D-04, 2.0D0, 1.D-06,fkinf, kinf)
!        kinf = 1.0D-03 
        write (50,*) '   ... Done.'
        write (50,*) '   kinf = ', kinf
        kd = max(kinf, kd)
        write (50,*) ' => kd = ', kd
        write (50,*) 'Computation of the coefficient ''c'' in Durden & V
     &esecky (1985) sea spectrum model ...'
        call c_ (c, U12)
        write (50,*) '   ... Done.'
        write (50,*) '   c = ',c
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

        write (50,*) 'Computation of slope variances ...'
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
                print *, '!!!!!!!!!!  ERROR  !!!!!!!!!!!!!! '
                print *, 'Invalid choice for the model of '
                print *, 'slope variances.'
                print *, 'Current choice is ', fVar
                print *, 'Possible choices are :'
                print *, '  - C for Cox and Munk (1955)'
                print *, '  - S for computing variances from the sea'//
     &                    ' spectrum model.'
                write (*,*)
                stop
        endif

        Su = sqrt(Su*Su + VarSwell(1))
        Sc = sqrt(Sc*Sc + VarSwell(2))

        write (50,*) '   ... Done.'
        write (50,*) '   Su = ', Su
        write (50,*) '   Sc = ', Sc

c----------------------------------------------------------------------
c               VARIANCES DES HAUTEURS DU SPECTRE DES PETITE VAGUES
c
c  sigma : Variances des hauteurs des petites échelles
c  kd : Nombre d'onde de coupure en rad/m
c  cspec : modele de spectre

        write (50,*) 'Computation of height variance ...'
        call sigma_ (sigma, kd, cspec)
        write (50,*) '   ... Done.'
        write (50,*) '   sigma = ', sigma

c----------------------------------------------------------------------
c                       FRACTION DE COUVERTURE D'ECUME
c

c  Couverture d'ecume avec dependance en Stab(stabilite atmos.)
        write (50,*) 'Computation of foam fraction ...'
        if (fCouvEcume.eq.1) then
                call foam (U10,-Stab(1),Fr)
        elseif (fCouvEcume.eq.2) then
                Fr = 2.95D-6*U10**(3.52)
        elseif (fCouvEcume.eq.3) then
                call MonahanLu(U10, SST(iSST), temp1, temp2, Fr)
        elseif (fCouvEcume.eq.4) then
                call MonahanLu(U10, SST(iSST), Fr, temp1, temp2)
        elseif (fCouvEcume.eq.5) then
                call WISE2001(U10, Fr)
        elseif ((fCouvEcume.ge.6).and.(fCouvEcume.le.10)) then
                call foam_fr_Yin16(U10, cCouvEcume, Fr)
        elseif (fCouvEcume.eq.0) then
                Fr  = 0.0D0
        else
                        write (*,*)
                        print *, '!!!!!!!!!!  ERROR !!!!!!!!!!!!!! '
                        print*,'Invalid choice for foam fraction model.'
                        write (*,*)
                        stop
        endif
        Fr_ = 1.0D0 - Fr
        write (50,*) ' ... Done.'
        write (50,*) '   Foam fraction = ', Fr
        close (50)

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

        print *, 'Tabulation of small scale scattering component started
     &...'
c Echantillonnage en theta de 0 à THETAmax par pas de DTHETAtab

c Début BOUCLE THETA_1>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        do nTHETA_1 =  0, int ( THETAmax/DTHETAtab )
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
c                print*, 'fin 12'
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
c                print*, 'fin 12'
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
c        write (40,*) nu, kd, SST(iSST), SSS(iSSS), U10, ustar, Omega, 
c     &  theta_1
c     &                 , Irh(0,1,nTHETA_1)*Tw
c     &                 , Irh(0,2,nTHETA_1)*Tw
c     &                 , Irh(2,1,nTHETA_1)*Tw
c     &                 , Irh(2,2,nTHETA_1)*Tw
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
c                print*, 'fin 34'
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
c                print*, 'fin 34'
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
   		write (60,*) theta_1, Irh(0,1,nTHETA_1)*Tw
     &                 , Irh(0,2,nTHETA_1)*Tw
     &                 , Irh(2,1,nTHETA_1)*Tw
     &                 , Irh(2,2,nTHETA_1)*Tw
     &                 , Irh(2,3,nTHETA_1)*Tw
     &                 , Irh(2,4,nTHETA_1)*Tw
	if (fGeo.eq.0) then
c Ecriture fichier et ecran
	Tvn = Tw*(1-abs(Rvv0)**2)
	Thn = Tw*(1-abs(Rhh0)**2)
	Tv0 = -Tw*Irh(0, 1, nTHETA_1)
	Th0 = -Tw*Irh(0, 2, nTHETA_1)
	Tv1 = -Tw*Irh(1, 1, nTHETA_1)
	Th1 = -Tw*Irh(1, 2, nTHETA_1)
	U1 = -Tw*Irh(1, 3, nTHETA_1)
	V1 = -Tw*Irh(1, 4, nTHETA_1)
	Tv2 = -Tw*Irh(2, 1, nTHETA_1)
	Th2 = -Tw*Irh(2, 2, nTHETA_1)
	U2 = -Tw*Irh(2, 3, nTHETA_1)
	V2 = -Tw*Irh(2, 4, nTHETA_1)

        write (20,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
     &                  theta_1, Tvn, Thn, Tv0, Th0,
     &                  Tv1, Th1, U1, V1, Tv2, Th2, U2, V2, lambdad
     &                  , Stab(1), realpart(epsi),imagpart(epsi), Su*Su,
     &                  Sc*Sc, Fr
        write (*,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
     &                 theta_1, Tvn, Thn, Tv0, Th0, 
     &                 Tv1, Th1, U1, V1, Tv2, Th2, U2, V2, lambdad
     &                 , Stab(1), realpart(epsi),imagpart(epsi), Su*Su,
     &                 Sc*Sc, Fr
	endif
c  Fin Boucle Theta_1<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        enddo
        elseif (fDiff.eq.0) then
                write(*,*)
                print *, '/!\ No Bragg scattering selected !'
                write(*,*)
        else
                write(*,*)
                print *, '/!\ Invalid choice for Bragg scattering '//
     &          'inclusion.'
                write (*,*)
        endif

	if (fGeo.ne.0) then ! seulement si on veut les grandes vaugues
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        
c   Integration : PARTIE GEOMETRIQUE 
        
c Début boucle theta>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
	do itheta = itheta0, ntheta
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

		ind = theta(itheta)/1 ! on prend la partie entière de theta_1
		Poids = (theta(itheta)/1 -ind) ! poids des diff Tb
		AtmoU = TbAu(ind)*(1 - Poids) + TbAu(ind+1)*Poids
		write(*,*) 'theta = ', theta(itheta)
		write(*,*) 'ind = ', ind
		write(*,*) 'Poids = ', Poids
		write(*,*) 'TbAu(ind), TbAu(ind+1)', TbAu(ind), TbAu(ind+1)
		write(*,*)  'AtmoU', AtmoU
        
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
                int1_GO(indice,iphi) = 0.0D0
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

        call Glob2Loc(theta(itheta), phi, Sx, Sy, thetaAtmo, temp1,
     &                temp2, temp3)

c TEST USE OF GLOBAL POLAR ANGLE OF INCIDENT ATMO RADIATION 
c INSTEAD OF LOCAL INCIDENCE ANGLE FOR TB ATMO
        call Loc2Glob(thetal, phil+180, Sx, Sy, thetaAtmo)
c        thetaAtmo = 180 - thetaAtmo
c		write(*,*) thetal, thetaAtmo
c ORIGINAL PHD VERSION (BUG??)
c NOTE : index out of bound (i.e theta > 90 deg) bring TbAd = 0
        AtmoD = TbAd(int(thetaAtmo))*(int(thetaAtmo)-thetaAtmo+1) 
     &         + TbAd(int(thetaAtmo)+1)*(thetaAtmo-int(thetaAtmo))
        tau_ = tau(int(theta(itheta)))
     &                  *(int(theta(itheta))-theta(itheta)+1)
     &         + tau(int(theta(itheta))+1)
     &                  *(theta(itheta)-int(theta(itheta)))

c	AtmoD = 0.0D0
c	tau_ = 0.0D0

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
c  epsi_sf : Brightness temperature of foam in polarisaion V (indice 1) et H(indice 2)
c-----------------------------------------------------------------------
      foamEmiss : select case (fEmisEcume)

         case (1)

                call esf (thetal, nu, epsi_sf)

        case (2:6)

                call foam_TbYin16_wrapper (fEmisEcume - 1, SSS(iSSS),
     &                          SST(iSST), nu, thetal, epsi_sf)

	case (7) 	! M-Du-Tune, Yin et al. (2016) tuned by freq and pol

                call foam_emiss(nuband, SSS(iSSS), SST(iSST),
     &                           nu, thetal, epsi_sf)

        end select foamEmiss

        !print *, '-->', thetal, epsi_sf(1), epsi_sf(2)

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
                                print *, 'ERROR !!!!!!!!!!!!!'
                                print *, 'thetal > THETAmax'
                                print *, thetal, ' > ', THETAmax
                                print *, 'Probability density of slopes 
     &= ',P_Sx_Sy
                                stop
                        endif


c-----------------------------------------------------------------------
c                       VECTEUR DE STOKES GLOBAL
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
c		For GO model only
          		Isl_GO(1) = Tw*(1.0D0-(Rvv0))
          		Isl_GO(2) = Tw*(1.0D0-(Rhh0))
          		Isl_GO(3) = 0.D0
       			Isl_GO(4) = 0.D0


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
c                     For GO model Only 

          		Is_GO(1) = Isl_GO(1)*cosalpha**2 
     &                        + Isl_GO(2)*sinalpha**2 
     &                        + Isl_GO(3)*cosalpha*sinalpha 
          		Is_GO(2) = Isl_GO(1)*sinalpha**2 
     &                        + Isl_GO(2)*cosalpha**2 
     &                        - Isl_GO(3)*cosalpha*sinalpha
                        Is_GO(3) = Isl_GO(3)*(cosalpha**2-sinalpha**2) -
     &                    (Isl_GO(1)-Isl_GO(2))*2.D0*sinalpha*cosalpha
                        Is_GO(4) = Isl_GO(4)


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
c                       For GO model only

          		int1_GO(1,iphi) = int1_GO(1,iphi) + 
     &                                 Is_GO(1)*projec*P_Sx_Sy*del1*del2
          		int1_GO(2,iphi) = int1_GO(2,iphi) + 
     &                                 Is_GO(2)*projec*P_Sx_Sy*del1*del2
          		int1_GO(3,iphi) = int1_GO(3,iphi) + 
     &                                 Is_GO(3)*projec*P_Sx_Sy*del1*del2
          		int1_GO(4,iphi) = int1_GO(4,iphi) + 
     &                                 Is_GO(4)*projec*P_Sx_Sy*del1*del2



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
c            for GO only
                int1_GO(1,iphi) = int1_GO(1,iphi)/Norm2
                int1_GO(2,iphi) = int1_GO(2,iphi)/Norm2
                int1_GO(3,iphi) = int1_GO(3,iphi)/Norm2
                int1_GO(4,iphi) = int1_GO(4,iphi)/Norm2
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

c SANS ATMO
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
        write (10,1000) nu, SST(iSST), SSS(iSSS), U10, 
     &                 ustar*100, theta(itheta), Tvn, Thn,
     &                (Tv0-Tvn), (Th0-Thn), Tv1, Th1, U1, V1, Tv2, 
     &                Th2, U2, V2, lambdad, Stab(1), realpart(epsi), 
     &                  imagpart(epsi), Su*Su, Sc*Sc, Fr*100.0D0
c        write (*,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
c     &                 theta(itheta), Tvn, Thn,
c     &                (Tv0-Tvn), (Th0-Thn), Tv1, Th1, U1, V1, Tv2, 
c     &                Th2, U2, V2, lambdad, Stab(1), realpart(epsi), 
c     &                  imagpart(epsi), Su*Su, Sc*Sc, Fr*100.0D0

c AVEC ATMO
        Tv0 = 0.25D0*(int1eA(1,1)+ 2.0D0*int1eA(1,3)+int1eA(1,4))+AtmoU
        Th0 = 0.25D0*(int1eA(2,1)+ 2.0D0*int1eA(2,3)+int1eA(2,4))+AtmoU
        Tv1 = 0.5D0 *(int1eA(1,1) - int1eA(1,4)) 
        Th1 = 0.5D0 *(int1eA(2,1) - int1eA(2,4))
        U1  = int1eA(3,3) 
        V1  = int1eA(4,3)
        Tv2 = 0.25D0*(int1eA(1,1) - 2.0D0*int1eA(1,3) + int1eA(1,4))
        Th2 = 0.25D0*(int1eA(2,1) - 2.0D0*int1eA(2,3) + int1eA(2,4))
        U2  = int1eA(3,2) - int1eA(3,3)
        V2  = int1eA(4,2) - int1eA(4,3)

c Ecriture fichier et ecran
        write (15,1000) nu, SST(iSST), SSS(iSSS), U10, 
     &                 ustar*100, theta(itheta), Tvn, Thn,
     &                (Tv0-Tvn), (Th0-Thn), Tv1, Th1, U1, V1, Tv2, 
     &                Th2, U2, V2, lambdad, Stab(1), realpart(epsi), 
     &                  imagpart(epsi), Su*Su, Sc*Sc, Fr*100.0D0
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

c AVEC ATMO
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
        write (25,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
     &                  theta(itheta), Tvn, Thn, (Tv0-Tvn), (Th0-Thn),
     &                  Tv1, Th1, U1, V1, Tv2, Th2, U2, V2, lambdad
     &                  , Stab(1), realpart(epsi),imagpart(epsi), Su*Su,
     &                  Sc*Sc, Fr-Fr
		write(*,*) 'With Atmosphere'
        write (*,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
     &                 theta(itheta), Tvn, Thn, (Tv0-Tvn), (Th0-Thn), 
     &                 Tv1, Th1, U1, V1, Tv2, Th2, U2, V2, lambdad
     &                 , Stab(1), realpart(epsi),imagpart(epsi), Su*Su,
     &                 Sc*Sc, Fr-Fr

c SANS ATMO
        Tv0 = 0.25D0*(int1(1,1)+2.0D0*int1(1,3)+int1(1,4))
        Th0 = 0.25D0*(int1(2,1)+2.0D0*int1(2,3)+int1(2,4))
        Tv1 = 0.5D0 *(int1(1,1) - int1(1,4)) 
        Th1 = 0.5D0 *(int1(2,1) - int1(2,4))
        U1  = int1(3,3) 
        V1  = int1(4,3)
        Tv2 = 0.25D0*(int1(1,1)-2.0D0*int1(1,3)+int1(1,4))
        Th2 = 0.25D0*(int1(2,1)-2.0D0*int1(2,3)+int1(2,4))
        U2  = int1(3,2) - int1(3,3)
        V2  = int1(4,2) - int1(4,3)
c For GO model only
        Tv0_GO = 0.25D0*(int1_GO(1,1)+2.0D0*int1_GO(1,3)+int1_GO(1,4))
        Th0_GO = 0.25D0*(int1_GO(2,1)+2.0D0*int1_GO(2,3)+int1_GO(2,4))
        Tv1_GO = 0.5D0 *(int1_GO(1,1) - int1_GO(1,4)) 
        Th1_GO = 0.5D0 *(int1_GO(2,1) - int1_GO(2,4))
        T31_GO  = int1_GO(3,3) 
        T41_GO  = int1_GO(4,3)
        Tv2_GO = 0.25D0*(int1_GO(1,1)-2.0D0*int1_GO(1,3)+int1_GO(1,4))
        Th2_GO = 0.25D0*(int1_GO(2,1)-2.0D0*int1_GO(2,3)+int1_GO(2,4))
        T32_GO  = int1_GO(3,2) - int1_GO(3,3)
        T42_GO  = int1_GO(4,2) - int1_GO(4,3)

        write(*,*) 'Without Atmosphere'

c   COMPUTE AND WRITE BRAGG SCATTERING

c SI ON RESTE DANS LE DOMAINE TABULÉ
                        if (theta(itheta).le.THETAmax) then
                                ind     = theta(itheta)/DTHETAtab
                                Poids   = (theta(itheta)/DTHETAtab-ind) 
c Fondamental
                        Tv0_SPM = Tw*(-(Irh(0,1,ind)*(1.0D0-Poids)
     &                            + Irh(0,1,ind+1)*Poids))
                        Th0_SPM = Tw*(-(Irh(0,2,ind)*(1.0D0-Poids)
     &                            + Irh(0,2,ind+1)*Poids))
                        T30_SPM = Tw*(-(Irh(0,3,ind)*(1.0D0-Poids)
     &                            + Irh(0,3,ind+1)*Poids))
                        T40_SPM = Tw*(-(Irh(0,4,ind)*(1.0D0-Poids)
     &                            + Irh(0,4,ind+1)*Poids))
c 2eme Harmonique
                        Tv2_SPM = Tw*(-(Irh(2,1,ind)*(1.0D0-Poids)
     &                                    + Irh(2,1,ind+1)*Poids))
                        Th2_SPM = Tw*(-(Irh(2,2,ind)*(1.0D0-Poids)
     &                                    + Irh(2,2,ind+1)*Poids))
                        T32_SPM = Tw*(-(Irh(2,3,ind)*(1.0D0-Poids)
     &                                   + Irh(2,3,ind+1)*Poids))
                        T42_SPM = Tw*(-(Irh(2,4,ind)*(1.0D0-Poids)
     &                                    + Irh(2,4,ind+1)*Poids))

c SI ON SORT DU DOMAINE TABULÉ À LA PRECISION NUMÉRIQUE PRÈS
                        elseif ((P_Sx_Sy.eq.0.D0).or.
     &                     (theta(itheta).le.THETAmax*1.0001D0)) then
                                do indice = 1, 4
                                        Diff(indice) = 0.0D0
                                enddo
                        else
c SI ON SORT DU DOMAINE, ERREUR => SORTIE DU PROGRAMME
                                print *, 'ERROR !!!!!!!!!!!!!'
                                print *, 'theta(itheta) > THETAmax'
                                print *, theta(itheta), ' > ', THETAmax
                                print *, 'Probability densoty of slopes 
     &= ',P_Sx_Sy
                                stop
                        endif



        write (20,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
     &                  theta(itheta), Tvn, Thn, (Tv0-Tvn), (Th0-Thn),
     &                  Tv1, Th1, U1, V1, Tv2, Th2, U2, V2, lambdad
     &                  , Stab(1), realpart(epsi),imagpart(epsi), Su*Su,
     &                  Sc*Sc, Fr-Fr, Tv0_SPM, Th0_SPM, T30_SPM
     &                  , T40_SPM, Tv2_SPM, Th2_SPM, T32_SPM,
     &                   T42_SPM, Tv0_GO, Th0_GO, Tv1_GO
     &                 , Th1_GO, T31_GO, T41_GO, Tv2_GO, Th2_GO, T32_GO,
     &                  T42_GO

        write (*,1000) nu, SST(iSST), SSS(iSSS), U10, ustar*100, 
     &                 theta(itheta), Tvn, Thn, (Tv0-Tvn), (Th0-Thn), 
     &                 Tv1, Th1, U1, V1, Tv2, Th2, U2, V2, lambdad
     &                 , Stab(1), realpart(epsi),imagpart(epsi), Su*Su,
     &                 Sc*Sc, Fr-Fr, Tv0_SPM, Th0_SPM, T30_SPM
     &                 , T40_SPM, Tv2_SPM, Th2_SPM, T32_SPM,
     &                  T42_SPM, Tv0_GO, Th0_GO, Tv1_GO
     &                 , Th1_GO, T31_GO, T41_GO, Tv2_GO, Th2_GO, T32_GO,
     &                  T42_GO


 
c Fin boucle Theta<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        enddo
! TO BE REMOVED !!! (use for continuation of computations)
      itheta0 = 1
c<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
	endif ! fin condition sur fGeo
c  Fin intégrale sur les grandes vagues

c  Fin Boucle kd
        enddo

c  Fin Boucle Vent
        enddo
! TO BE REMOVED !!! (use for continuation of computations)
      iWind0 = 1
c  Fin Boucle SST
        enddo
! TO BE REMOVED !!! (use for continuation of computations)
      iSST0 = 1

c  Fin Boucle SSS
        enddo
! TO BE REMOVED !!! (use for continuation of computations)
      iSSS0 = 1



        close (10)
        close (15)
        close (20)
        close (25)
        close (30)
        close (60)

c ---------- Formats --------------------------      
  100     format (1x,f5.2,4(1x,f4.1),2(1x,f6.2),2(1x,f6.2)
     &          ,8(1x,f5.2))
 1000     format (1x,f12.5,1x,f5.2,1x,f6.3,2(1x,f6.2),1x,f4.1,4(1x,f7.3),
     &        8(1x,f6.3),1x,e9.3,1x,f5.2,2(1x,f6.2),2(1x,e9.3),1x,f5.1,
     &       18(1x,f7.3))
        stop
c ---------- Error messages --------------------------      
   60      print*, 'Unable to open or read ', fin1
        stop
   61      print*, 'Unable to open ', KudryFilter1
        stop
   62      print*, 'Unable to open ', KudryFilter2
        stop
   63      print*, 'Unable to open ', KudryFilter3
        stop
        
   70      print*, 'Unable to create TbAtmo.dat in Output directory'
        stop
   80      print*, 'Unable to create ', logFile
        stop


        end

