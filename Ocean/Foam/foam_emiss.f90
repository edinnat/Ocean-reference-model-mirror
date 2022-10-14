subroutine foam_emiss(perm,SSTin,Frin,ThIn,ef_foam)

! -----------------------------------------------------------------------------------------------
! Created by Lise Kilic, France: Translated to Fortran from Yin et al. MATLAB code 
! Modified by Magdalena Anguelova, Remote Sensing Division, NRL, Washngton, DC, USA:
!	* foam parameters (thickness t and upper limit of the void fraction profile vaf) 
!		tuned for frequency and polarization
!	* t value for each freq 
!	* vaf values for V and H pol at each freq
!
! Modified by Emmanuel Dinnat: removed the additional freqency input ("band" variable) because it could be set 
! independently of the frequency used in the rest of the model ("Frin" variable here), it was limited to 6 
! discrete frequencies; it was used to define foam model parameters t (thickness) and foam fraction vafv & vafh.
! Instead, the routine now uses the same input freq as the rest of the model (Frin) as a unique frequency,
!  and it interpolates linearly the foam parameters t, vafv and vafh between the values set for the 6 reference
!  frequencies between 1.4 GHz and 89 GHz. If the freq is lower than the 1.4 GHz lower bound, values are set
!  the same as for 1.4 GHz, and if it is higher than 89 GHz, the 89 GHz values are used.
! -----------------------------------------------------------------------------------------------

USE mod_foam_emiss

IMPLICIT NONE
     
! Inputs     
 	REAL (KIND=8) :: SSTin,Frin,ThIn
  	REAL (KIND=8) :: sst, freq, thetai 
    DOUBLE COMPLEX :: perm
! Internal
    integer :: i, i_t
  	REAL (KIND=8) :: evfoam, ehfoam, t, vafv, vafh
  	REAL (KIND=8) :: vafv_t(1:6), vafh_t(1:6), t_t(1:6), f_t(1:6), wght

! Output				foam emiss in V and H pol for index 1 and 2, respectively
  	REAL (KIND=8) :: ef_foam(1:2)


    f_t = (/1.4D0, 6.9D0, 10.6D0, 18.7D0, 36.5D0, 89.0D0/)
    t_t = (/2.0D0, 0.6D0, 0.4D0, 0.2D0, 0.1D0, 0.1D0/)
    vafv_t = (/0.950D0, 0.950D0, 0.950D0, 0.950D0, 0.980D0, 0.970D0/)
    vafh_t = (/0.950D0, 0.960D0, 0.964D0, 0.968D0, 0.970D0, 0.980D0/)

    if (Frin.lt.f_t(1)) then

        t = t_t(1)
        vafv = vafv_t(1)
        vafh = vafh_t(1)
   
    endif

    if (Frin.gt.f_t(6)) then

        t = t_t(6)
        vafv = vafv_t(6)
        vafh = vafh_t(6)
  
    endif
    
    if ((Frin.ge.f_t(1)).and.(Frin.le.f_t(6))) then

    do i_t = 1, 6
        
        if ((Frin.ge.f_t(i_t)).and.(Frin.le.f_t(i_t+1))) then

            wght = (Frin - f_t(i_t))/(f_t(i_t+1) - f_t(i_t))
            t = t_t(i_t)*(1.D0 - wght) + t_t(i_t+1)*wght
            vafv = vafv_t(i_t)*(1.D0 - wght) + vafv_t(i_t+1)*wght
            vafh = vafh_t(i_t)*(1.D0 - wght) + vafh_t(i_t+1)*wght

        exit

        endif

    enddo

    endif


	freq = FrIn
	thetai = ThIn
	
	sst = SSTin + 2.7315D2


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
