c---------------------------------------------------------------------
c Calcul de la constante dielectrique de l'eau de mer d'apres Klein et 
c Swift (1977) +++  25 Mars 1999  +++
c  * SST est la temperature de surface de la mer en degres celsius.
c  * SSS est la salinite de surface de la mer en PSU
c    (partie par millier).
c  * epsi est la constante dielectrique epsilon complexe, donnee par l'
c    equation de Debye.
c  * freq est la frequence electromagnetique du radiometre en Hertz.
c---------------------------------------------------------------------


       subroutine epsilon_KS (SST, SSS, epsi_KS, freq)

       implicit none
       
       double complex epsi_KS 
       double complex D 
       double complex im1 
       double complex im2
       double precision SST
       double precision sigma 
       double precision sigma_25 
       double precision SSS 
       double precision epsi_inf 
       double precision epsi_s 
       double precision ohmega 
       double precision freq
       double precision epsi_0 
       double precision a 
       double precision b 
       double precision delta 
       double precision beta
       double precision tau_0 
       double precision epsi_s_T
       double precision tau 
       double precision alpha

**********************************************************************
*                        INITIALISATION DES CONSTANTES               *
**********************************************************************
       data epsi_0/8.854D-12/
       data epsi_inf/4.9D0/
       data alpha/0.D0/
       ohmega = freq*(2.D0*3.141592653D0)

**********************************************************************
*                             CALCUL DE SIGMA                        *
**********************************************************************

       sigma_25 = SSS*(1.82521D-1 - 1.46192D-03*SSS + 2.09324D-05*SSS**2
     &            - 1.28205D-07*SSS**3)
       delta = 2.5D01 - SST
       beta = 2.033D-02 + 1.266D-04*delta + 2.464D-06*delta**2 -
     &        SSS*(1.849D-05 - 2.551D-07*delta + 2.551D-08*delta**2)
       sigma = sigma_25*dexp(-delta*beta)

**********************************************************************
*                           CALCUL DE EPSILON_S                      *
**********************************************************************

       epsi_s_T = 8.7134D1 - 1.949D-01*SST - 1.276D-02*SST**2 + 
     &            2.491D-04*SST**3
       a = 1.D0 + 1.613D-05*SSS*SST - 3.656D-03*SSS + 
     &     3.210D-05*SSS**2 - 4.232D-07*SSS**3
       epsi_s = epsi_s_T*a

**********************************************************************
*                             CALCUL DE TAU                          *
**********************************************************************
  
       tau_0 = 1.768D-11 - 6.086D-13*SST + 1.104D-14*SST**2 
     &         - 8.111D-17*SST**3
       b = 1.D00 + 2.282D-05*SSS*SST - 7.638D-04*SSS 
     &     - 7.760D-06*SSS**2 + 1.105D-08*SSS**3
       tau = tau_0*b

**********************************************************************
*                           CALCUL DE EPSILON                        *
**********************************************************************

       im1 = (0.D0,1.D0)*(ohmega * tau)
       im2 = (0.D0,1.D0)*(sigma/(ohmega*epsi_0))
       D = 1.D0 - im1**(1-alpha)
       epsi_KS = epsi_inf + (epsi_s - epsi_inf)/D + im2

       end




