subroutine foam_emiss(band, SalIn,SSTin,Frin,ThIn,ef_foam)

! -----------------------------------------------------------------------------------------------
! Created by Lise Kilic, France: Translated to Fortran from Yin et al. MATLAB code 
! Modified by Magdalena Anguelova, Remote Sensing Division, NRL, Washngton, DC, USA:
!	* foam parameters (thickness t and upper limit of the void fraction profile vaf) 
!		tuned for frequency and polarization
!	* t value for each freq 
!	* vaf values for V and H pol at each freq
! -----------------------------------------------------------------------------------------------

USE mod_foam_emiss
	
IMPLICIT NONE
     
! Inputs     
        character*16, intent(in) :: band
 	REAL (KIND=8) :: SalIn,SSTin,Frin,ThIn
  	REAL (KIND=8) :: sst, sss, freq, thetai 
! Internal
	integer :: i
  	REAL (KIND=8) :: evfoam, ehfoam, t, vafv, vafh	
  	DOUBLE COMPLEX :: perm

! Output				foam emiss in V and H pol for index 1 and 2, respectively
  	REAL (KIND=8) :: ef_foam(1:2)						 

! Input test
!	freq = 1.4D0
!  	thetai=55.0D0

!  	sss=34.0D0
!  	sst=20.0D0
!  	sst = sst + 2.7315D2

!  	t = 2.0
!  	vafv = 0.94D0
! 	vafh = 0.97D0 

	freq = FrIn
	thetai = ThIn
	
	sst = SSTin + 2.7315D2
	sss = SalIn

! Foam parameters by frequency and polarization
	if (band.eq.'L') then		! 1.4 GHz
		t = 2.0D0
		vafv = 0.95D0		! V pol	
		vafh = 0.95D0		! H pol
	elseif (band.eq.'C') then	! 6.9 GHz
		t = 0.6D0
		vafv = 0.95D0
		vafh = 0.96D0
	elseif (band.eq.'X') then	! 10.6 GHz
		t = 0.4D0
		vafv = 0.95D0
		vafh = 0.964D0
	elseif (band.eq.'K') then	! 18.7 GHz
		t = 0.2D0
		vafv = 0.95D0
		vafh = 0.968D0
	elseif (band.eq.'Ka') then	! 36.5 GHz
		t = 0.1D0
		vafv = 0.98D0
		vafh = 0.97D0
	elseif (band.eq.'W') then	! 89.0 GHz
		t = 0.1D0
		vafv = 0.97D0
		vafh = 0.98D0
        else
                write(*,*)
                print *,'Invalid choice for frequency band.'
     		print *,' Possible choices are:'
                print *,'  L : 1.4 GHz'
                print *,'  C : 6.9 GHz'
              	print *,'  X : 10.6 GHz'
                print *,'  K : 18.7 GHz'
                print *,'  Ka : 36.5-37 GHz'
                print *,'  W : 89.0 GHz'
                print *,'Current choice is: ',band
                write(*,*)
                stop
       endif
		!print *, t, vafv, vafh

! Computation of the water permittivity 
	call epsilon_MW (sst-2.7315D2, sss, perm, freq*1.D9)

! The original foam emiss subroutine uses one void fraction value to return 2 ef values, for V and H pol
! To use the same subroutine with minimal change, call it twice, 
! once with foam fraction upper limit for V pol, then for H pol
! For each call, only the ef value for the pol corresponding to the input is correct
 
	CALL esf_anguelova(thetai,freq,perm,t,vafv,evfoam,ehfoam)	! Only V pol is correct
		ef_foam(1) = evfoam
	CALL esf_anguelova(thetai,freq,perm,t,vafh,evfoam,ehfoam)	! Only H pol is correct
		ef_foam(2) = ehfoam

	ef_foam = ef_foam*sst
		!print *, thetai, freq, perm
		!print *, t, vafv, vafh, evfoam, ehfoam
		!print *, ef_foam
	
	RETURN

END subroutine
