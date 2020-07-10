c
c                                             /--- alphaz
c                        /--- loss_factor ---/
c                       /                    \
c foam_emissivity -----/                      \--- thetaz --- betaz --- void_fraction_profile 
c                      \                                                            \   
c                       \                                                            \--- eps_klein_swift
c                        \                                                           /
c                         \                                                         /
c                          \--- IncohEmissivity_TwoLayer --------------------------/ 
c

c Description: Code computes incoherent emissivity of an inhomogeneous layer
c separated by air on top and a homogeneous medium on the bottom, with
c perfectly smooth parallel boundaries on both sides.

c Input Variables: 
c    %eps2: dielectric constant of middle layer
c    %eps3: dielectric constant of bottom layer
c    %theta_i: incidence angle in air (deg)
c    %a: Single scattering albedo
c    %d: layer thickness (m)
c    %kappa_e: extinction rate through layer (Np/m)
c    
c%Output Products:
c    %e_v_inc(theta_i): v-polarized emissivity
c    %e_h_inc(theta_i): h-polarized emissivity
c    


      subroutine IncohEmissivity_TwoLayer(eps2,eps3, thetai, a, kappa_e
     &                                    e_v, e_h)

      ! Inputs
      complex*16 :: eps2, eps3
      real*8 :: thetai, a, kappa_e
      ! Outputs
      real*8 :: ev, eh
      ! Other
      complex*16 :: n1, n2, rhoh, rhov, tauh, tauv
      complex*16 :: eps1 = ( 1.D0, 0.D0 )
      real*8 :: theta2, gammah12, gammav12, gammah23, gammav23, Th, Tv
      real*8 :: theta1, trans, PI = dacos(-1.D0)


      theta1 = thetai*PI/180.D0 ! transform to radians
      
      n1 = sqrt( eps2 ) ! calc index of refraction

      n2 = sqrt( eps3 )

      theta2 = asin(abs(n1/n2)*dsin ( theta1 ) ) ! incidence angle in medium 2 (refracted)

c-- calculate reflectivies at the two interfaces

c      !- 12 interface
      
      call ReflTransm_PlanarBoundary(eps1, eps2, thetai, rhoh, rhov,
     & gammah12, gammav12, tauh, tauv, Th, Tv)

c    !- 23 interface

      call ReflTransm_PlanarBoundary(eps2, eps3, theta2, rhoh, rhov,
     & gammah23, gammav23, tauh, tauv, Th, Tv)

c     ! extinction coefficient inside medium.
c     ! trans = exp(-kappa_e .* d ./cosd(theta2)); %% constant

      trans = exp( - kappa_e )

      e_v = ( (1.D0 - gammav12)/(1.D0 - gammav12*gammav23*trans*trans) )
     &  * ( (1.D0 + gammav23*trans)*(1.D0 - a)*(1.D0 - trans) 
     &  + ( 1.D0 - gammav23) * trans )

      e_h = ( (1.D0 - gammah12)/(1.D0 - gammah12*gammah23*trans*trans) )
     &  * ( (1.D0 + gammah23*trans)*(1.D0 - a)*(1.D0 - trans) 
     &  + ( 1.D0 - gammah23) * trans )

      end subroutine IncohEmissivity_TwoLayer


c %Code 2.3: Oblique Reflection and Transmission @ Planar Boundry
c %Description: Code computes the reflection coefficients, transmission
c %coefficients, reflectivities and transmissivities for incidence in
c %medium (medium 1) upon the planar boundary of a lossless or lossy
c %medium (medium 2) at any incidence angle, for both h and v polarizations

c %Input Variables:
c     %eps1: eps1r -j*eps1i: relative dielectric constant of medium 1
c     %eps2 = eps2r-j*eps2i: relative dielectric constant of medium 2
c     %theta1d: incidence angle in medium 1 in degrees

cOutput Products:
c   %rhoh: reflection coefficient for h pol
c   %rhov: reflection coefficient for v pol
c   %gammah:reflectivity for h pol
c   %gammav: reflectivity for v pol
c   %tauh: transmission coefficient for h pol
c   %tauv: transmission coefficient for v pol
c   %Th: transmissivity for h pol
c   %Tv:transmissivity for v pol
cBook Reference: Sections 2-7 & 2-8
c%Example call: [rhoh rhov gammah gammav tauh tauv Th Tv] = ReflTransm_PlanarBoundary(eps1, eps2, theta1d)
cComputes 
c%MATLAB Code

      subroutine ReflTransm_PlanarBoundary(eps1, eps2, theta1d, 
     &                   rhoh, rhov, gammah, gammav, tauh, tauv, Th, Tv)

      implicit none

      ! Inputs
      complex*16 :: eps1, eps2
      real*8 :: theta1d
      ! Outputs
      complex*16 :: rhoh, rhov, tauh, tauv
      real*8 :: gammah, gammav, Th, Tv
      ! Other
      real*8 :: theta1, sin_theta2, cos_theta2
      real*8 :: PI = dacos(-1.D0)


      theta1 = theta1d*PI/180.D0
    
      sin_theta2 = sqrt(eps1)/sqrt(eps2)*sin(theta1);
      cos_theta2 = sqrt(1 - sin_theta2*sin_theta2);
    
      rhoh = (sqrt(eps1)*cos(theta1)-sqrt(eps2)*cos_theta2) 
     &  / (sqrt(eps1)*cos(theta1) + sqrt(eps2)*cos_theta2)
      rhov = (sqrt(eps1)*cos_theta2-sqrt(eps2)*cos(theta1)) 
     &  / (sqrt(eps1)*cos_theta2 + sqrt(eps2)*cos(theta1))

      tauh =  1.D0 + rhoh
      tauv = (1.D0 + rhov)*(cos(theta1)/cos_theta2)
       
      gammah = abs(rhoh)**2
      gammav = abs(rhov)**2
        
      Th = 1.D0 - gammah
      Tv = 1.D0 - gammav
    
      end subroutine ReflTransm_PlanarBoundary


      
      function alphaz(z, thick, vaf, vfw, m, temp, sal, freq)

      implicit none

      real*8 :: z,thick,vaf,vfw,m,temp,sal,freq ! inputs
      real*8 :: alphaz ! output
      real*8 :: k0, sigfoam, void_fraction_profile
      real*8 :: PI = dacos(-1.D0)

      k0 = 2.D0*PI*freq*1.D9/3.D8

      sigfoam = void_fraction_profile(z, thick,vaf,vfw,m,temp,sal,freq)

      alphaz = k0*abs(imag(sqrt(complex(sigfoam,0.D0))))

      end function alphaz




      function thetaz(z, thick, thetai, vaf, vfw, m, temp, sal, freq)

      implicit none

      real*8 :: z, thick, thetai, vaf, vfw, m, temp, sal, freq ! inputs
      real*8 :: thetaz ! output
      real*8 :: k0, az, bz, pz, qz, alphaz, betaz
      real*8 :: PI = dacos(-1.D0)
      real*8 :: thetai_r

      thetai_r = thetai/180.D0*PI

      
      k0 = 2*PI*freq*1.D9/3.D8

      az = alphaz(z, thick, vaf, vfw, m, temp, sal, freq)
      
      bz = betaz (z, thick, vaf, vfw, m, temp, sal, freq)

      pz = 2*az*bz
     
      qz = bz*bz - az*az - k0*k0*dsin(thetai)

      thetaz = atan(sqrt(2.D0)*k0*dsin(thetai) 
     &         /sqrt ( sqrt ( pz*pz + qz*qz ) + qz ) )

      end function thetaz

      function betaz (z, thick, vaf, vfw, m, temp, sal, freq)

      implicit none

      real*8 :: z, thick, vaf, vfw, m, temp, sal, freq ! input
      real*8 :: betaz ! output
      real*8 :: k0, sigfoam, void_fraction_profile
      real*8 :: PI = dacos(-1.D0)

      k0 = 2*PI*freq*1e9/3e8

      sigfoam = void_fraction_profile (z, thick, vaf, vfw, m, temp, sal,
     &                                 freq)
      
      betaz = k0*real(sqrt(complex(sigfoam, 0.D0) ) ) 

      end function betaz

      function void_fraction_profile (z, thick, vaf, vfw, m, temp, sal,
     &                                 freq)

      implicit none

      real*8 :: z, thick, vaf, vfw, m, temp, sal, freq ! input
      real*8 :: void_fraction_profile ! output
      real*8 :: av, bv, fa
      complex*16 epsi_KS

      av = vaf + m
      bv = 1/thick*log( (av - vfw) / m )
      fa = av - m*exp ( bv*z )

c      % fa(fa<vfw)=vfw;

      ! Pemittivity profile
      
      call  epsilon_KS (temp, sal, epsi_KS, freq)

      void_fraction_profile = ( fa + (1 - fa)*sqrt( epsi_KS ) )**2

      end function void_fraction_profile


      function loss_factor(z, thick, thetai, vaf, vfw, m, temp, sal, 
     &                      freq) 

      implicit none

!      external alphaz
!      external thetaz

      real*8 :: z, thick, thetai, vaf, vfw, m, temp, sal, freq ! Inputs
      real*8 :: loss_factor ! output 
      real*8 :: az, tz, alphaz, thetaz

      az = alphaz(z,thick,vaf,vfw,m,temp,sal,freq)

      tz  = thetaz(z,thick,thetai,vaf,vfw,m,temp,sal,freq)

      loss_factor = 2*az/dcos(tz)

      end function loss_factor 



