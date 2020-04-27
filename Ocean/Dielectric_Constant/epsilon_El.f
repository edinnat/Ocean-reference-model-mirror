c---------------------------------------------------------------------
c Calcul de la constante dielectrique de l'eau de mer d'apres Ellison 
c (1998) +++ 10  Nov 1999  +++
c  * SST est la temperature de surface de la mer en degres celsius.
c  * SSS est la salinite de surface de la mer en PSU
c    (partie par millier).
c  * epsi est la constante dielectrique epsilon complexe, donnee par l'
c    equation de Debye.
c  * freq est la frequence electromagnetique du radiometre en Hertz.
c---------------------------------------------------------------------


	subroutine epsilon_El (SST, SSS, epsi_El, freq)

	implicit none
       
	double complex epsi_El
        double complex im
	double precision SST
        double precision SSS 
        double precision epsi_inf 
        double precision epsi_s
        double precision tau
        double precision freq
        double precision ohmega 
        double precision epsi_0
	double precision C1 
        double precision C2 
        double precision A1 
        double precision A2 
        double precision B1 
        double precision B2
	double precision sigma 

**********************************************************************
*                        INITIALISATION DES CONSTANTES               *
**********************************************************************
	data epsi_0/8.884D-12/
	ohmega = freq*(2.D0*3.141592653D0)

**********************************************************************
*                             CALCUL DE SIGMA                        *
**********************************************************************

	C1 = 8.6374D-02 + 3.0606D-02*SST - 4.121D-04*SST**2
	C2 = 7.7454D-02 + 1.687D-03*SST + 1.937D-05*SST**2
	sigma = C1 + C2*SSS

**********************************************************************
*                           CALCUL DE EPSILON_S                      *
**********************************************************************

	A1 = 8.182D01 - 6.0503D-02*SST - 3.1661D-02*SST**2 + 3.1097D-03
     &	*SST**3 - 1.1791D-04*SST**4 + 1.4838D-06*SST**5
	A2 = 1.2544D-01 + 9.4037D-03*SST - 9.5551D-04*SST**2 + 9.0888D-0
     &	5*SST**3 - 3.6011D-06*SST**4 + 4.7130D-08*SST**5
	epsi_s = A1 - A2*SSS

**********************************************************************
*                             CALCUL DE TAU                          *
**********************************************************************
  
	B1 = 1.7303D01 - 6.6651D-01*SST + 5.1482D-03*SST**2 + 1.2145D-03
     &	*SST**3 - 5.0325D-05*SST**4 + 5.8272D-07*SST**5
	B2 = -6.272D-03 + 2.357D-04*SST + 5.075D-04*SST**2 - 6.3983D-05
     &	*SST**3 + 2.463D-06*SST**4 - 3.0676D-08*SST**5
	tau = (B1 + B2*SSS)*1.D-12

**********************************************************************
*	Calcul de la permittivite aux hautes frequences	             *
**********************************************************************

	epsi_inf = 6.4587 - 0.04203*SST - 0.0065881*SST**2 + 6.4924D-04
     &	*SST**3 - 1.2328D-05*SST**4 + 5.0433D-08*SST**5

**********************************************************************
*                           CALCUL DE EPSILON                        *
**********************************************************************

	im = (0.,1.)
c  partie reelle d'epsilon
	epsi_El = epsi_inf+(epsi_s - epsi_inf)/(1.D0 + (ohmega*tau)**2) 
c  partie imaginaire  d'epsilon
	epsi_El = ((epsi_s-epsi_inf)*ohmega*tau/(1.D0+ohmega**2*tau**2)+
     &	sigma/ohmega/epsi_0)*im + epsi_El

	end




