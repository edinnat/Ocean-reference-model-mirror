Module mod_foam_emiss

! -----------------------------------------------------------------------------------------------
! Model of Anguelova & Gaiser (2013, RSE) with Yin et al. (2016) modifications for fast calculations 
! Fast calculations obtained by using semi-closed form of the incoherent approach (with only one
! integration over the foam layer) instead of the general form (with 2 integrations from each stratum
! upward and downward). 
!
! Created by Lise Kilic, France: Translated to Fortran from Yin et al. MATLAB code 
! Modified by Magdalena Anguelova, Remote Sensing Division, NRL, Washngton, DC, USA:
! 	* Semi-closed (Yin et al.) and general (nguelova & Gaiser) forms of the emissivity 
!	  incoherent approach verified 
!	* A few elements in the code corrected
!	* Differences between semi-closed and genral forms quantified 
!	* This comparison is documented in a white paper
! -----------------------------------------------------------------------------------------------

IMPLICIT NONE
DOUBLE PRECISION :: pi=3.141592653589793238462643383279D0
CONTAINS

	SUBROUTINE esf_anguelova(thetai,hz,seaDC,tin,vaf,evfoam,ehfoam)

	IMPLICIT NONE
! INPUT
	REAL (KIND=8) thetai 			! incidence angle (deg) from air onto foam 
	REAL (KIND=8) hz 			! frequency in gigahertz
	REAL (KIND=8) tin			! foam thickness in cm
	REAL (KIND=8) vaf			! upper limit of void fraction
	COMPLEX (KIND=8) seaDC			! seawater permittivity
! OUTPUT
	REAL (KIND=8) evfoam
	REAL (KIND=8) ehfoam
! INTERNAL
	REAL (KIND=8) aa
	REAL (KIND=8) vfw
	REAL (KIND=8) t
	REAL (KIND=8) m
	COMPLEX (KIND=8) sigfoam1
	COMPLEX (KIND=8) sigfoam2
	REAL (KIND=8) taof
	REAL (KIND=8) sigw_re
	REAL (KIND=8) sigw_im
! Variables for Simpson integration 
	REAL (KIND=8) a,b, aire 
	INTEGER n
	REAL (KIND=8) h,SommePaire, SommeImpaire
	INTEGER i

	aa=0
	
	vfw=0.01				! foam void fraction at the foam-seawater boundary set at 1%
	t=tin/100. 				! thickness input is in cm, convert here to m
	m=1					! shape of void fraction profile is nominal exponential function

! Mixing rule for foam dielectric constant using seawater permittivity and foam void fraction
	sigw_re=REAL(seaDC)
	sigw_im=AIMAG(seaDC)
	sigfoam1=(vaf+(1-vaf)*sqrt(Cmplx(sigw_re,sigw_im)))**2;
	sigfoam2=(vfw+(1-vfw)*sqrt(Cmplx(sigw_re,sigw_im)))**2;

! INTEGRAL of SIMPSON 
	a=0.D0
	b=t
	n=10

	aire = 0.D0
	h = (b-a)/(n*2)
	SommePaire = 0.D0
	SommeImpaire = 0.D0

! Calculation of the sum of the even indices
	DO i=1, n-1
	SommePaire = SommePaire + &
	       loss_factor((a+h*2.*i),t,thetai,vaf,vfw,m,hz,sigw_re,sigw_im)
    ENDDO

! Calculation of the sum of the odd indices
	DO i=1, n
	SommeImpaire = SommeImpaire + &		
	       loss_factor((a+h*(2.*i-1)),t,thetai,vaf,vfw,m,hz,sigw_re,sigw_im)
    ENDDO

! Final computation of the total area
	aire = h*(loss_factor(a,t,thetai,vaf,vfw,m,hz,sigw_re,sigw_im) + &
       loss_factor(b,t,thetai,vaf,vfw,m,hz,sigw_re,sigw_im) + &
       2.*SommePaire + 4.*SommeImpaire)/3.

	taof=aire		! total extinction in the foam layer

	CALL IncohEmissivity_TwoLayer(sigfoam1,sigfoam2,seaDC,thetai,aa,taof,evfoam, ehfoam)

	RETURN

	END 

!----------End of principal subroutine ------------------------------------


! %Description: Code computes incoherent emissivity of an inhomogeneous layer
! %separated by air on top and a homogeneous medium on the bottom, with
! %perfectly smooth parallel boundaries on both sides.

! %Input Variables: 
!    %eps2: dielectric constant of middle layer
!    %eps3: dielectric constant of bottom layer
!    %theta_i: incidence angle in air (deg)
!    %a: Single scattering albedo
!    %d: layer thickness (m)
!    %kappa_e: extinction rate through layer (Np/m)
    
! %Output Products:
!    %e_v_inc(theta_i): v-polarized emissivity
!    %e_h_inc(theta_i): h-polarized emissivity
    
! %Book Reference: Section 12-12.2

	SUBROUTINE IncohEmissivity_TwoLayer(eps2,eps3,epsw, theta_i,a, kappa_e,e_v_inc,e_h_inc) 

	IMPLICIT NONE
! ---INPUT---
	COMPLEX (KIND=8) eps2	
	COMPLEX (KIND=8) eps3
	COMPLEX (KIND=8) epsw
	REAL (KIND=8) theta_i
	REAL (KIND=8) a
	REAL (KIND=8) kappa_e
! ---OUTPUT---
	REAL (KIND=8) e_v_inc
	REAL (KIND=8) e_h_inc
! ---INTERNE---
	COMPLEX (KIND=8) eps1
	REAL (KIND=8) theta1
	REAL (KIND=8) theta2
	COMPLEX (KIND=8) n1
	COMPLEX (KIND=8) n2
	REAL (KIND=8) rhoh
	REAL (KIND=8) rhov
	REAL (KIND=8) gammah12
	REAL (KIND=8) gammav12
	REAL (KIND=8) gammah23
	REAL (KIND=8) gammav23    
	REAL (KIND=8) tauh
	REAL (KIND=8) tauv
	REAL (KIND=8) Th
	REAL (KIND=8) Tv
	REAL (KIND=8) trans

	eps1 = 1				! Air permittivity (layer 1)
	n1 = sqrt(eps1) 			! Index of refraction for air, no further use
	n2 = sqrt(eps2)				! Index of refraction for foam, no further use

	theta1 = theta_i*pi/180 				
	theta2 = asin(abs(n1/n2) * sin(theta1))*180/pi 	! Incidence angle in medium 2 (refracted)


! -- calculate reflectivies at the two interfaces

! 12 interface (air-foam in our case; use foam permittivity for the air-foam boundary) 
	CALL ReflTransm_PlanarBoundary(eps1,eps2,theta_i,rhoh,rhov,gammah12,gammav12,tauh,tauv,Th,Tv);

! 23 interface (foam-seawater in our case; use foam permittivity for the foam-seawater boundary)
	CALL ReflTransm_PlanarBoundary(eps3,epsw,theta2 ,rhoh,rhov,gammah23,gammav23,tauh,tauv,Th,Tv);

! extinction coefficient inside medium (foam layer in our case)
	trans = exp(-kappa_e);

	e_v_inc = ( (1.D0 - gammav12)/(1.D0 - gammav12*gammav23*trans*trans) ) * &
		  ( (1.D0 + gammav23*trans)*(1.D0 - a)*(1.D0 - trans) + &
		    (1.D0 - gammav23) * trans )

	e_h_inc = ( (1.D0 - gammah12)/(1.D0 - gammah12*gammah23*trans*trans) ) * &
     		  ( (1.D0 + gammah23*trans)*(1.D0 - a)*(1.D0 - trans) + &
		    (1.D0 - gammah23) * trans )
 
	RETURN 

	END 


! %Code 2.3: Oblique Reflection and Transmission @ Planar Boundry
! %Description: Code computes the reflection coefficients, transmission
! %coefficients, reflectivities and transmissivities for incidence in
! %medium (medium 1) upon the planar boundary of a lossless or lossy
! %medium (medium 2) at any incidence angle, for both h and v polarizations

! %Input Variables:
!    %eps1: eps1r -j*eps1i: relative dielectric constant of medium 1
!    %eps2 = eps2r-j*eps2i: relative dielectric constant of medium 2
!    %theta1d: incidence angle in medium 1 in degrees

!%Output Products:
!    %rhoh: reflection coefficient for h pol
!    %rhov: reflection coefficient for v pol
!    %gammah:reflectivity for h pol
!    %gammav: reflectivity for v pol
!    %tauh: transmission coefficient for h pol
!    %tauv: transmission coefficient for v pol
!    %Th: transmissivity for h pol
!    %Tv:transmissivity for v pol
!%Book Reference: Sections 2-7 & 2-8

!%Example call: [rhoh rhov gammah gammav tauh tauv Th Tv] = ReflTransm_PlanarBoundary(eps1, eps2, theta1d)
!%Computes 

	SUBROUTINE ReflTransm_PlanarBoundary(eps1, eps2, theta1d, rhoh, rhov, gammah, gammav, tauh, tauv, Th, Tv)

	IMPLICIT NONE 
! ---INPUT---
	COMPLEX (KIND=8) eps1
	COMPLEX (KIND=8) eps2
	REAL (KIND=8) theta1d
! ---OUTPUT---
	REAL (KIND=8) rhoh
	REAL (KIND=8) rhov
	REAL (KIND=8) gammah
	REAL (KIND=8) gammav
	REAL (KIND=8) tauh
	REAL (KIND=8) tauv
	REAL (KIND=8) Th
	REAL (KIND=8) Tv  
! ---INTERNE--- 
	REAL (KIND=8) theta1
	REAL (KIND=8) sin_theta2
	REAL (KIND=8) cos_theta2

    	theta1 = theta1d*pi/180
    
    	sin_theta2 = sqrt(eps1)/sqrt(eps2)*sin(theta1);
    	cos_theta2 = sqrt(1 - sin_theta2**2);
    
    	rhoh = (sqrt(eps1)*cos(theta1)-sqrt(eps2)*cos_theta2) / (sqrt(eps1)*cos(theta1) + sqrt(eps2)*cos_theta2);
    	rhov = (sqrt(eps1)*cos_theta2-sqrt(eps2)*cos(theta1)) / (sqrt(eps1)*cos_theta2 + sqrt(eps2)*cos(theta1));

    	tauh = 1 + rhoh;
    	tauv = (1 + rhov)*(cos(theta1)/cos_theta2);
       

    	gammah = abs(rhoh)**2;
    	gammav = abs(rhov)**2;
        
    	Th = 1-gammah;
    	Tv = 1-gammav;
    
	RETURN

  	END


! ---loss factor----------------------------------------------------------------

	REAL (Kind=8) FUNCTION loss_factor(z,t,thetai,vaf,vfw,m,hz,sigw_re,sigw_im)
! calul taof
	IMPLICIT NONE
! ---INPUT---
	REAL (Kind=8) z
	REAL (KIND=8)  t
	REAL (KIND=8)  thetai
	REAL (KIND=8)  vaf
	REAL (KIND=8)  vfw
	REAL (KIND=8)  m
	REAL (KIND=8)  hz
	REAL (KIND=8)  sigw_re
	REAL (KIND=8)  sigw_im
! ---OUTPUT---
!	REAL (KIND=8) :: loss_factor
! ---INTERNE---
	REAL (KIND=8) :: az
	REAL (KIND=8) :: tz

	az=alphaz(z,t,vaf,vfw,m,hz,sigw_re,sigw_im)
	tz=thetaz(z,t,thetai,vaf,vfw,m,hz,sigw_re,sigw_im)

	loss_factor=2*az/cos(tz)

	RETURN

	END FUNCTION loss_factor

! ---void fraction profile------------------------------------------------------

	COMPLEX (KIND=8) FUNCTION void_fraction_profile(z,t,vaf,vfw,m,hz,sigw_re,sigw_im)
! retourne la valeur de permitivité de l'écume sigfoam

	IMPLICIT NONE
! ---INPUT---
	REAL (KIND=8)  z
	REAL (KIND=8)  t
	REAL (KIND=8)  vaf
	REAL (KIND=8)  vfw
	REAL (KIND=8)  m
	REAL (KIND=8)  hz
	REAL (KIND=8)  sigw_re
	REAL (KIND=8)  sigw_im
! ---OUTPUT--- 
!	COMPLEX (KIND=8)  void_fraction_profile
! ---INTERNE---
	REAL (KIND=8)  av
	REAL (KIND=8)  bv
	REAL (KIND=8)  fa

	av=vaf+m
	bv=1./t*log((av-vfw)/m)
	fa=av-m*exp(bv*z)

! fa(fa<vfw)=vfw;
! Permittivity profile
! [sigw_re, sigw_im]=eps_klein_swift(T,S,hz);

	void_fraction_profile=(fa+(1-fa)*sqrt(Cmplx(sigw_re,sigw_im)))**2

	RETURN

	END FUNCTION void_fraction_profile


! ----------------thetaz--------------------------------------------------------
	REAL (KIND=8) FUNCTION thetaz(z,t,thetai,vaf,vfw,m,hz,sigw_re,sigw_im)
	IMPLICIT NONE
! ---INPUT---
	REAL (KIND=8)  z
	REAL (KIND=8)  t
	REAL (KIND=8)  thetai
	REAL (KIND=8)  vaf
	REAL (KIND=8)  vfw
	REAL (KIND=8)  m
	REAL (KIND=8)  hz
	REAL (KIND=8)  sigw_re
	REAL (KIND=8)  sigw_im
! ---OUTPUT---
!	REAL (KIND=8)  thetaz
! ---INTERNE---
	REAL (KIND=8)  k0
	REAL (KIND=8)  az
	REAL (KIND=8)  bz
	REAL (KIND=8)  pz
	REAL (KIND=8)  qz

	k0=2*pi*hz*1.D9/3.D8

	az=alphaz(z,t,vaf,vfw,m,hz,sigw_re,sigw_im)
	bz=betaz(z,t,vaf,vfw,m,hz,sigw_re,sigw_im)

	pz=2*az*bz
	!qz=bz**2-az**2-k0**2*sin(thetai*pi/180.)
	qz=bz**2-az**2-k0**2*(sin(thetai*pi/180.))**2

	thetaz=atan(sqrt(2.)*k0*sin(thetai*pi/180.)/sqrt(sqrt(pz**2+qz**2)+qz))

	RETURN

	END FUNCTION thetaz

! -------betaz------------------------------------------------------------------
	REAL (KIND=8) FUNCTION betaz(z,t,vaf,vfw,m,hz,sigw_re,sigw_im)

	IMPLICIT NONE
! ---INPUT---
	REAL (KIND=8)  z
	REAL (KIND=8)  t
	REAL (KIND=8)  vaf
	REAL (KIND=8)  vfw
	REAL (KIND=8)  m
	REAL (KIND=8)  hz
	REAL (KIND=8)  sigw_re
	REAL (KIND=8)  sigw_im
! ---OUTPUT---
!	REAL (KIND=8)  betaz
! ---INTERNE---
	REAL (KIND=8)  k0
	COMPLEX (KIND=8)  sigfoam

! The index of refraction, of air, is 1.000293. 
! The Speed of Light, in a vacuum, is defined to be 299792458 m/s 
! So the speed of light in air is 299792458/1.000293, or about 299704644.54 m/s


	k0=2*pi*hz*1.D9/3.D8;

	sigfoam=void_fraction_profile(z,t,vaf,vfw,m,hz,sigw_re,sigw_im)
	betaz=k0*real(sqrt(sigfoam))

	RETURN

	END FUNCTION betaz

! -----------------alphaz-------------------------------------------------------
	REAL (KIND=8) FUNCTION alphaz(z,t,vaf,vfw,m,hz,sigw_re,sigw_im)

	IMPLICIT NONE
! ---INPUT----
	REAL (KIND=8)  z
	REAL (KIND=8)  t
	REAL (KIND=8)  vaf
	REAL (KIND=8)  vfw
	REAL (KIND=8)  m
	REAL (KIND=8)  hz
	REAL (KIND=8)  sigw_re
	REAL (KIND=8)  sigw_im
! ---OUTPUT---
!	REAL (KIND=8)  alphaz
! ---INTERNE---
	REAL (KIND=8)  k0
	COMPLEX (KIND=8)  sigfoam

	k0=2*pi*hz*1.D9/3.D8

	sigfoam=void_fraction_profile(z,t,vaf,vfw,m,hz,sigw_re,sigw_im)
	alphaz=k0*abs(aimag(sqrt(sigfoam)))

	RETURN

	END FUNCTION alphaz

! -----END OF MODULE E SEA FOAM------------------------------------
END
