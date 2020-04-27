c--------------------------------------------------------------------
c  Calcul de la valeur du spectrum de puissance de la mer, pour un 
c  nombre d onde k
c  modele : Yueh (1997).
c           description: based on Durden and Vesecky (1985) model 
c                   with following modifications:
c                   1) amplitude of the high frequency domain
c                      is not set but computed to insure continuity
c                      at k = 2 rad/m (i.e. limit between high and low
c                       frequency domains).
c                   2) amplitude of the low frequency spectrum is
c                      changed from 0.004 to 0.008
c                   3) mistypped equations are taken into account 
c                      (see spectrum_DV_Base.f)
c
c                       +++++  17 October 2003  +++++
c  * k est en m-1
c--------------------------------------------------------------------


        function S_Yueh (k)

		implicit none

		common /kc/ kc
		double precision k
        double precision a
        double precision b 
        double precision B_ 
        double precision S_Yueh 
        double precision b0
		double precision S_DV_Base
		double precision kc

        data a/0.225D0/
		data b/1.25D0/
		data B_/0.008/
        b0   = B_*dexp(1.85D-01*kc*kc)

		S_Yueh = S_DV_Base(k, a, b, B_, b0)
        return

        end
