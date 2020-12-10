        subroutine epsilon_MW (SST, SSS, epsi_MW, f) 

! Computes the dielectric constant of sea water according to the model
! by Meissner and Wentz (2004) including corrections by Meissner and Wentz
! 2012 and Meissner et al. 2014.
!
!       Inputs:
!               SST : Sea surface temperature in degrees C
!               SSS  : Sea surface salinity in psu
!               f    : frequency in Hz
!
!       Outputs:
!               epsi_MW : complex dielectric constant of sea water
!
!
!       E. Dinnat

! Meissner, T., & Wentz, F. J. (2004). The complex dielectric constant of pure
! and sea water from microwave satellite observations. IEEE Transactions on
! Geoscience and Remote Sensing, 42(9), 18361849.
! https://doi.org/10.1109/TGRS.2004.831888
!
! Meissner, T., & Wentz, F. J. (2012). The Emissivity of the Ocean Surface
! Between 6 and 90 GHz Over a Large Range of Wind Speeds and Earth Incidence
! Angles. IEEE Transactions on Geoscience and Remote Sensing, 50(8), 30043026.
! https://doi.org/10.1109/TGRS.2011.2179662
!
! Meissner, T., Wentz, F. J., & Ricciardulli, L. (2014). The emission and
! scattering of L-band microwave radiation from rough ocean surfaces and wind
! speed measurements from the Aquarius sensor. Journal of Geophysical Research C:
! Oceans, 119(9), 64996522. https://doi.org/10.1002/2014JC009837


      implicit none

      ! Inputs
      REAL*8, INTENT(IN) :: SST, SSS, f
      ! Outputs
      COMPLEX*16, INTENT(OUT) :: epsi_MW

      ! Intermediate variables
!     ai coefficients from Table III in Meissner and Wentz 2004
      REAL*8, PARAMETER :: ai(11)=(/ 5.7230D+00, 2.2379D-02
     &      , -7.1237D-04, 5.0478D+00, -7.0315D-02, 6.0059D-04
     &      , 3.6143D+00, 2.8841D-02, 1.3652D-01,  1.4825D-03,
     &              2.4166D-04 /)
!     bi coefficients from Table VI in Meissner and Wentz 2004
      REAL*8, PARAMETER :: bi(13)=(/-3.56417D-03,  4.74868D-06
     &      , 1.15574D-05, 2.39357D-03, -3.13530D-05, 2.52477D-07
     &      , -6.28908D-03, 1.76032D-04, -9.22144D-05, -1.99723D-02,
     &            1.81176D-04, -2.04265D-03,  1.57883D-04/)
!    updated bi (i = 0 -> 2) coefficients from Table VI in Meissner and Wentz
!    2012
      REAL*8, PARAMETER :: bi_new1(3)=(/ -0.33330D-02,  4.74868D-06, 
     &                                  0.0D0/)
!    updated bi (i = 3 -> 5) coefficients from Table VII in Meissner and Wentz
!    2012 (now uses 5 coefficients instead of 3), includes change of
!    sign on coef#4 reported in Meissner et al. 2016
      REAL*8, PARAMETER :: bi_new2(5)=(/0.23232D-02, -0.79208D-04
     &      , 0.36764D-05, -0.35594D-06, 0.89795D-08/)
      REAL*8 :: fGHz, f0, SST2, SST3, SST4, SSS2
      REAL*8 :: es, e1, einf, nu1,nu2 ! Fresh water parameters
      REAL*8 :: sigma35, R15, alpha0, alpha1, RTR15, sigma
      REAL*8 :: es_s, e1_s, einf_s, nu1_s, nu2_s ! Salt water parameters
      REAL*8 :: c1, c2, c3, c4 ! conversion coef from fresh to salt water
      COMPLEX*16, PARAMETER :: j = (0.D0, 1.D0)

      if (SST.le.-40.D0) then

        print *, 'Error in epsilon_MW. SST is lower than 40 degree C.'
        error stop

      endif
        

      fGHz = f/1.D09 ! convert fequency from Hz to Ghz
      f0 = 17.97510D0 ! (1/(2*pi*epsi_0) term in conductivity term in GHz.m/S

      SST2 = SST * SST
      SST3 = SST2 * SST
      SST4 = SST3 * SST

      SSS2 = SSS * SSS
 
      ! Fresh water parameters
      ! es:  Static dielectric constant for pure water by Stogryn et al. 1995
      es = ( 3.70886D4 - 8.2168D1*SST ) / ( 4.21854D2 + SST )
      e1 = ai(1) + ai(2)*SST + ai(3)*SST2
      nu1 = ( 45.00D0 + SST ) / ( ai(4) + ai(5)*SST + ai(6)*SST2 )
      einf = ai(7) + ai(8)*SST
      nu2 = ( 45.00D0 + SST ) / ( ai(9) + ai(10)*SST + ai(11)*SST2 )
 
      ! Salt water
      ! Conductivity [s/m] by Stogryn et al. 1995
      sigma35 = 2.903602D0 + 8.60700D-2*SST + 4.738817D-4*SST2 
     &        - 2.9910D-6*SST3 + 4.3047D-9*SST4
      R15 = SSS*( 37.5109D0 + 5.45216D0*SSS + 1.4409D-2*SSS2 ) 
     &        / ( 1004.75D0 + 182.283D0*SSS + SSS2 )
      alpha0 = ( 6.9431D0 + 3.2841D0*SSS - 9.9486D-2*SSS2 )
     &        / (84.850D0 + 69.024D0*SSS + SSS2 )
      alpha1 = 49.843D0 - 0.2276D0*SSS + 0.198D-2*SSS2
      RTR15 = 1.0D0 + ( SST - 15.0D0 )*alpha0/ ( alpha1 + SST ) 
      sigma = sigma35 * R15 * RTR15

      ! Salt water parameters 
      es_s = es *exp ( bi_new1(1)*SSS + bi_new1(2)*SSS2 
     &                  + bi_new1(3)*SSS*SST ) 

      if ( SST .le. 30.D0 ) then

        c1 = 1.0D0 + SSS * ( bi_new2(1) + bi_new2(2)*SST 
     &          + bi_new2(3)*SST2 + bi_new2(4)*SST3 + bi_new2(5)*SST4 )

      else

        c1 = 1.0D0 + SSS * 
     &          ( 9.1873715D-04 + 1.5012396D-04*(SST - 30.D0 ) )

      endif

      nu1_s = nu1*c1

      c2  = exp ( bi(7)*SSS + bi(8)*SSS2 + bi(9)*SSS*SST ) 
      e1_s = e1*c2

!     c3 = 1.0D0 + SSS*(bi(10) + bi(11)*SST) ! Expression for c3  amended in Meissner et al. 2016
      c3 = 1.0D0 + SSS*(bi(10) + 0.5D0*bi(11) * (SST + 30.D0) )
      nu2_s = nu2*c3
      c4 = 1.0D0  + SSS * ( bi(12) + bi(13) * SST )
      einf_s = einf*c4  

      epsi_MW = ( es_s - e1_s  ) / ( 1.0D0 - j * ( fGHz/nu1_s ) ) 
     &        + ( e1_s - einf_s ) / ( 1.0D0 - j * ( fGHz/nu2_s ) )
     &        + einf_s + j*sigma*f0/fGHz

      end subroutine epsilon_MW  

