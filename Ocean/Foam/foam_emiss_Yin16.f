c
c                                             /--- alphaz
c                        /--- loss_factor ---/
c                       /                    \
c foam_emissYin16 -----/                      \--- thetaz --- betaz --- void_fraction_profile 
c                      \                                                            \   
c                       \                                                            \--- eps_klein_swift
c                        \                                                           
c                         \                                                         
c                          \--- IncohEmissivity_TwoLayer ----------- ReflTransm_PlanarBoundary  
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

    
c   DEBUG MAIN ROUTINE

      module test

      implicit none

      real*8 :: h_fe,theta_i,Vaf,Vfw,m_fm,SST,SSS, freq

      contains

      function  loss_factor_int(z)

      implicit none

      real*8 :: z
      real*8 :: loss_factor_int, loss_factor


      loss_factor_int = loss_factor(z,h_fe,theta_i,Vaf,Vfw,m_fm,SST,SSS,
     &                                   freq/1e09)
c      print *, 'loss f --', z, h_fe, theta_i,Vaf,Vfw,m_fm,SST,SSS,freq
c      print *, 'loss f --', loss_factor_int

      end function loss_factor_int

      end module

c      program foam_test

c      implicit none

c      real*8 :: ev_foam, eh_foam
c      real*8 ::  test, loss_factor

c      test = loss_factor(0.61254272460936954d0,
c     & 1.77d0,40.D0,0.9d0,1.d-002,1.d0,15.d0,35.d0,1.3999999999999d0)

c      print *, 'test = ', test

c      stop

c      call foam_emissYin16(35.d0, 15.d0, 1.41309d9, 0.97d0, 0.01d0, 1.d0
c     &                          , 40.d0, 1.77d0, ev_foam, eh_foam)

c        print *,  ev_foam, eh_foam

c      end 

      subroutine foam_TbYin16_wrapper(modl, S, T, f, tht, e_foam)

        implicit none

        ! Inputs
        integer modl ! selector of model parameterization as row number from Table 1 in Yin et al. 2016
        real*8 :: S  ! SSS in psu
        real*8 :: T  ! SST in degC
        real*8 :: f  ! frequency in GHz 
        real*8 :: tht  ! incidence angle in degrees
        ! Outputs
        real*8 , INTENT(OUT) :: e_foam(1:2) ! foam Tb in v and h pol for index 1 and 2 respectively
        ! Intermediate
        real*8 :: ev_foam, eh_foam, Vaf0, Vfw0, m_fm0, h_fe0

        call foam_emissYin16_wrapper ( modl, S, T, f, tht, e_foam )

        e_foam(1) = e_foam(1) * ( T + 273.15D0 )
        e_foam(2) = e_foam(2) * ( T + 273.15D0 )

      end subroutine foam_TbYin16_wrapper

 
      subroutine foam_emissYin16_wrapper(modl, S, T, f, tht, e_foam)

        implicit none

        ! Inputs
        integer modl ! selector of model parameterization as row number from Table 1 in Yin et al. 2016
        real*8 :: S  ! SSS in psu
        real*8 :: T  ! SST in degC
        real*8 :: f  ! frequency in GHz 
        real*8 :: tht  ! incidence angle in degrees
        ! Outputs
        real*8 , INTENT(OUT) :: e_foam(1:2) ! foam emissivity in v and h pol for index 1 and 2 respectively
        ! Intermediate
        real*8 :: ev_foam, eh_foam, Vaf0, Vfw0, m_fm0, h_fe0

        Vfw0 = 0.01D0
        m_fm0 = 1.D0

        modelEmiss: select case (modl)
        case (1) ! model M-Du-E
               Vaf0 = 0.97D0
               h_fe0 = 1.85D0
        case (2) ! model M-Du-E1
               Vaf0 = 0.97D0
               h_fe0 = 1.77D0
        case (3) ! model M-Du-S
               Vaf0 = 0.98D0
               h_fe0 = 1.83D0
        case (4) ! model M-Ku-E
               Vaf0 = 0.94D0
               h_fe0 = 1.99D0
        case (5) ! model M-Ku-S
               Vaf0 = 0.96D0
               h_fe0 = 1.91D0
        case default
                print *, 'Invalid choice for foam emissivity model in'
                print *, 'Yin 2016'
                stop
        end select modelEmiss

        call foam_emissYin16 (S, T, f*1.D9, Vaf0, Vfw0, 
     &                           m_fm0, tht, h_fe0, ev_foam, eh_foam)

        e_foam(1) = ev_foam
        e_foam(2) = eh_foam


      end subroutine foam_emissYin16_wrapper

      subroutine foam_emissYin16 (SSS1, SST_dc, freq1, Vaf1, Vfw1, 
     &                           m_fm1, theta1, h_fe1, ev_foam, eh_foam)

      use test

      implicit none

      ! Inputs
      real*8 :: freq1   ! frequency in Hz
      real*8 :: SST_dc ! sea surface temperature in degreeC
      real*8 :: SSS1    ! sea surface salinity in psu
      real*8 :: Vaf1    ! void fraction at air/foam interface
      real*8 :: Vfw1    ! void fraction at foam/water interface
      real*8 :: m_fm1   ! profile shape factor 
      real*8 :: theta1  ! incidence angle in degrees
      real*8 :: h_fe1   ! foam layer effective thickness (cm)
      real*8 :: a       ! 
      ! Intermediate
      real*8 :: sigw_re, sigw_im, taof
      complex*16 :: sigfoam1, sigfoam2, epsks, i = (0.0, 1.0)
      ! Outputs
      real*8 :: ev_foam !foam emissivity at v-pol
      real*8 :: eh_foam !foam emissivity at h-pol

      freq = freq1   ! frequency in Hz
      SST = SST_dc ! sea surface temperature in degreeC
      SSS = SSS1    ! sea surface salinity in psu
      Vaf = Vaf1    ! void fraction at air/foam interface
      Vfw = Vfw1    ! void fraction at foam/water interface
      m_fm = m_fm1   ! profile shape factor 
      theta_i = theta1  ! incidence angle in degrees
      h_fe = h_fe1   ! foam layer effective thickness (cm)

      call epsilon_KS(SST_dc, SSS, epsks, freq)

      sigw_re = real(epsks) 
      sigw_im = imag(epsks)  
      
      sigfoam1 = (Vaf+(1-Vaf)*sqrt(sigw_re+sigw_im*i))**2
      sigfoam2 = (Vfw+(1-Vfw)*sqrt(sigw_re+sigw_im*i))**2       


      call qromb ( loss_factor_int, 0, h_fe, taof, 1.D-6)

      call IncohEmissivity_TwoLayer(sigfoam1, sigfoam2, theta1, a, taof
     &                              ,ev_foam, eh_foam)

      end subroutine foam_emissYin16

      subroutine IncohEmissivity_TwoLayer(eps2,eps3, thetai, a, kappa_e
     &                                    , e_v, e_h)

      implicit none

      ! Inputs
      complex*16 :: eps2, eps3
      real*8 :: thetai, a, kappa_e
      ! Outputs
      real*8 :: e_v, e_h
      ! Other
      complex*16 :: n1, n2, rhoh, rhov, tauh, tauv
      complex*16 :: eps1 = ( 1.D0, 0.D0 )
      real*8 :: theta2, gammah12, gammav12, gammah23, gammav23, Th, Tv
      real*8 :: theta1, trans, PI = dacos(-1.D0)


      theta1 = thetai*PI/180.D0 ! transform to radians
      
      n1 = sqrt( eps2 ) ! calc index of refraction

      n2 = sqrt( eps3 )

      theta2 = asin(abs(n1/n2)*dsin ( theta1 ) )/PI*180.D0 ! incidence angle in medium 2 (refracted)

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

c      print *, 'taof -- ', kappa_e
c      print *, 'gv12, gv23, gh12, gh23, trans -- ', gammav12, gammav23, 
c     &  gammah12, gammah23, trans

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
    

      sin_theta2 = abs(sqrt(eps1)/sqrt(eps2))*dsin(theta1);
      cos_theta2 = sqrt(1 - sin_theta2*sin_theta2);

c      print *, 'sin_theta2 -- ', sin_theta2
    
      rhoh = (sqrt(eps1)*dcos(theta1) - sqrt(eps2)*cos_theta2) 
     &     / (sqrt(eps1)*dcos(theta1) + sqrt(eps2)*cos_theta2)
      rhov = (sqrt(eps1)*cos_theta2   - sqrt(eps2)*dcos(theta1)) 
     &     / (sqrt(eps1)*cos_theta2   + sqrt(eps2)*dcos(theta1))

      tauh =  1.D0 + rhoh
      tauv = (1.D0 + rhov)*(dcos(theta1)/cos_theta2)
       
      gammah = abs(rhoh)**2
      gammav = abs(rhov)**2
        
      Th = 1.D0 - gammah
      Tv = 1.D0 - gammav

c      print *, 'ReflTransm_ -- ', eps1, eps2, dsin(theta1), rhoh, rhov,
c     &     (sqrt(eps1)*dcos(theta1) - sqrt(eps2)*cos_theta2),
c     &    (sqrt(eps1)*dcos(theta1) + sqrt(eps2)*cos_theta2)
    
      end subroutine ReflTransm_PlanarBoundary


      
      function alphaz(z, thick, vaf, vfw, m, temp, sal, freq)

      implicit none

      real*8 :: z,thick,vaf,vfw,m,temp,sal,freq ! inputs
      real*8 :: alphaz ! output
      real*8 :: k0
      complex*16 :: sigfoam, void_fraction_profile
      real*8 :: PI = dacos(-1.D0)

      k0 = 2.D0*PI*freq*1.D9/3.D8

      sigfoam = void_fraction_profile(z, thick,vaf,vfw,m,temp,sal,freq)

      alphaz = k0*abs(imag(sqrt(sigfoam)))

c      print *, 'alphaz -- ',  z, thick, vaf, vfw, m, temp, sal, freq,
c     &                       sigfoam, alphaz

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
     
      qz = bz*bz - az*az - k0*k0*dsin(thetai_r)

      thetaz = atan(sqrt(2.D0)*k0*dsin(thetai_r) 
     &         /sqrt ( sqrt ( pz*pz + qz*qz ) + qz ) )

c      print *, 'thetaz -- ', k0, az, bz, pz, qz, thetaz

      end function thetaz

      function betaz (z, thick, vaf, vfw, m, temp, sal, freq)

      implicit none

      real*8 :: z, thick, vaf, vfw, m, temp, sal, freq ! input
      real*8 :: betaz ! output
      real*8 :: k0
      complex*16 :: sigfoam, void_fraction_profile
      real*8 :: PI = dacos(-1.D0)

      k0 = 2*PI*freq*1e9/3e8

      sigfoam = void_fraction_profile (z, thick, vaf, vfw, m, temp, sal,
     &                                 freq)
      
      betaz = k0*real(sqrt(sigfoam ) ) 

      end function betaz

      function void_fraction_profile (z, thick, vaf, vfw, m, temp, sal,
     &                                 freq)

      implicit none

      ! Inputs
      ! z:
      ! thick
      ! vaf
      ! vfm
      ! m
      ! temp : temprature in degree C 
      ! sal  : salinity in psu
      ! freq : frequency in GHz
      real*8 :: z, thick, vaf, vfw, m, temp, sal, freq ! input
      real*8 :: av, bv, fa
      complex*16 :: epsi_KS
      complex*16 :: void_fraction_profile ! output

      av = vaf + m
      bv = 1.D0/thick*log( (av - vfw) / m )
      fa = av - m*exp ( bv*z )

c      % fa(fa<vfw)=vfw;

      ! Pemittivity profile
      
      call  epsilon_KS (temp, sal, epsi_KS, freq*1D09)

      void_fraction_profile = ( fa + (1 - fa)*sqrt( epsi_KS ) )**2

c      print *, 'void f -- ', freq, epsi_KS, void_fraction_profile

      end function void_fraction_profile


      function loss_factor(z, thick, thetai, vaf, vfw, m, temp, sal, 
     &                      freq) 

      implicit none

!      external alphaz
!      external thetaz

      real*8 :: z, thick, thetai, vaf, vfw, m, temp, sal, freq ! Inputs
      real*8 :: loss_factor ! output 
      real*8 :: az, tz, alphaz, thetaz

c      print *, z, thick, thetai, vaf, vfw, m, temp, sal, freq

      az = alphaz(z,thick,vaf,vfw,m,temp,sal,freq)

      tz  = thetaz(z,thick,thetai,vaf,vfw,m,temp,sal,freq)

      loss_factor = 2*az/dcos(tz)

c      print *, az , tz

      end function loss_factor 



