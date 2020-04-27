      subroutine SpectreWLemaire (k, phi_ang, WLem, SLem, DeltaR)

c	Execution:
c		Compute the Lemaire et al. (1999) power spectra.
c		
c
c	Variables In:
c           k     : wave number (rad/m),
c           phi_ang : direction w/r wind vector (deg),
c					available models ('y' model recommended as used to make the spectrum model),
c           ss    : significant slope = Hs/lambda0 (typical 0.006-0.007),
c           m     : shape factor
c           beta  : amplitude factor      
c           x     : adimensionnal fetch (g*x'/wind**2, x' is dimensional fetch), optional, default value = 5e7.
c           U10   : wind speed (m/s) at 10 meters height      
c           U12   : wind speed (m/s) at 12.5 meters height      
c           us    : friction velocity (m/s)
c
c	Variables Out:
c			WLem  : sea waves power spectrum (m4).
c           DeltaR : delta ratio      
c
c	Author: E. P. Dinnat
c	Lastest update: 03/02/2004
c	Contact: http://www.lodyc.jussieu.fr/~dinnat

c Declarations      

      implicit none
      common /paramLem/ paramLem      
      
        double precision k, phi_ang, ss, m, beta, x, U10, U12, us, WLem
        double precision Root2, Pi, log2, X0, GAMMA, KM, K1, G, u0
        double precision lim_us
        double precision Xa, f0, w0, k0, kp, kcopy, kmax, k2km2, w, ww0
        double precision phi, dwdk, sg, k2, a0, b0, gs, scg, s1, s2sqr
        double precision sc, SLem, r, Fazi, DeltaR
        double precision dcosd, paramLem(10)
        external dcosd
c Constants
      Root2 = 1.41421356237309504880d0
      Pi = 3.14159265358979323846d0
      log2 = 0.69314718055994530942d0
      X0 = 2.2d4
      GAMMA = 7.4D-5
      KM = 363d0 ! high frequency peak, defined after eq. 2.12
      K1 = 100d0 ! cutoff wavenumber (rad/m !!) for gravity waves, eq 2.37
      G = 9.807d0 ! gravity cte
      u0 = 0.01d0
      lim_us = 0.08d0
c params
      ss = paramLem(1)
      m = paramLem(2)
      beta = paramLem(3)
      x = paramLem(4)
      U10 = paramLem(5)
      U12 = paramLem(6)
      us = paramLem(7)

c Peak Wavenumber
      Xa = x*G/U10/U10 ! fetch adim.
c      Xa = x ! Sinon pic ne correspond pas a Pierson !!!!!!!!!!!!!!!!!!!!
      Xa = (Xa/X0)**(0.4) !eq. 2.22
      f0 = 0.14*G/(U12*1.05)*dtanh(Xa)**(-0.7576) !eq 2.22, 0.14 = 0.877/2/pi 
      w0 = 2.0d0*Pi*f0
      k0 = w0*w0/G ! eq 2.23, max of frequency spectrum
c      k0 = 7.54/U12/U12/1.05/1.05
c      w0 = sqrt(G*k0)
      kp = dsqrt(m/m+1)*k0 ! max of wavenumber spectrum

c Gravity range
 
      kcopy = 0.0d0
      kmax = 10.0d0*kp ! limit for gravity range
      if (k.gt.kmax) then ! beyond gravity range
        kcopy = k
        k = kmax
      endif
 
      k2km2 = k*k/(KM*KM)
      w = dsqrt( G*k*( 1.0d0 + k2km2) )  ! eq. dispersion 2.12
      ww0 = w/w0
      phi = beta*G*G*w0**(-5) ! eq. 2.14 (note : extra 2pi w/r article removed)
      phi = phi*ww0**(-m) ! "
      phi = phi*dexp(-0.25d0*m*ww0**(-4.0d0)) ! "
      dwdk = 0.5d0*G*( 1.0d0 + 3.0d0*k2km2 )/w
      sg = phi*dwdk ! eq. 2.11

      if (kcopy.ne.0.0d0) then ! if beyond gravity range, trend towards k**-3 dependance
        k = kcopy
        sg = sg*(kmax/k)**(3.0)
      endif

c Capillary-gravity range
      k2 = k*k
      k2km2 = k*k/(KM*KM)
      if (us.lt.lim_us) then ! eq. 2.35
        a0 = 1.0d-6
        b0 = 5.2d0
      else
        a0 = 4.1d-3
        b0 = 1.2d0
      endif 
      gs = G*( 1.0d0 + k2km2 ) ! eq. 2.31
      scg = 0.5d0*a0*us**(b0)*u0**(1.0d0-b0)
      scg = scg*gs**(-1.5)*k**(-2.5)
      scg = scg*G*( 1.0d0 + 3.0d0*k2km2 )

c Stitch ranges together 1

      if (k.lt.kp) then
        SLem = sg
      else 
        sc = 1.473d-4*KM**(6)*us**(3)*k**(-9)*u0**(-3)
        s1 = sg*dexp(-k2/(K1*K1))
        s2sqr = 1.0d0/(1.0d0/(scg*scg)+1.0d0/(sc*sc))
        SLem = dsqrt( s1*s1 + s2sqr )
      endif

c /* stitch ranges together 2 */

c if (k .lt. kp) then
c   SLem = sg
c else 
c   s1 = sg*exp( -k2/(K1*K1) )
c   SLem = sqrt( s1*s1 + scg*scg )
c endif

c Stitch ranges together 3
c
c if ( k .le. kp ) then
c  *p = sg
c else 
c   fm = k/KM - 1.0
c   scg *= exp( -0.25*fm*fm )
c   s1 = sg*exp( -k2/(K1*K1) )
c   *p = sqrt( s1*s1 + scg*scg )
c   endif

c Azimuthal variation
      r = (3.0+1.92*U12)/(3.16*U12) ! Cox and Munk Ratio between Up/Crosswind Slope variances, eq 2.43
      DeltaR = 2.0d0*(1.0d0-r)/(1.0d0+r)
      Fazi = (1.0+DeltaR*dcosd(2.0*phi_ang))/2/Pi ! eq 2.39
c      Fazi = 1.0d0/2/Pi
      WLem = SLem*Fazi/k ! eq. 2.10


      return
      end
