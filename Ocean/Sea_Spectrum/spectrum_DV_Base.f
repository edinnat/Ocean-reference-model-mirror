c--------------------------------------------------------------------
c  Calcul de la valeur du spectrum de puissance de la mer, pour un 
c  nombre d onde k
c  modele : Durden et Vesecky (1985).
c           Description : modification with respect to the original
c           article introduced in Yueh (1997) (mistypes only !!!)
c              1) the aerodynamical roughness scale was mistyped
c              2) the angular speading coeeficients were mistyped
c                       +++++  17 Oct 2003  +++++
c  * k est en m-1
c--------------------------------------------------------------------


        function S_DV_Base (k, a, b, B_, b0)

		implicit none
        common /g/ g
        common /kc/ kc
        common /ustar/ ustar

		double precision k
        double precision a
        double precision b 
        double precision B_ 
        double precision kj 
        double precision gstar 
        double precision kc
        double precision g 
        double precision S_DV_Base
        double precision ustar 
        double precision gamma 
        double precision b0

c*********************************************************************
c*                   INITIALISATION DES CONSTANTE                    *
c*********************************************************************

		data gamma/7.25D-05/
		data kj/2.0D0/
        if (k.le.kj) then
                S_DV_Base = b0/k/k/k*dexp(-0.74D0*kc*kc/k/k)
				return
        else
	        gstar = g + gamma*k*k
            S_DV_Base = B_/k/k/k*(b*k*ustar*ustar/gstar)
     &          **(a*dlog10(k/kj))
			return
        endif

        end
