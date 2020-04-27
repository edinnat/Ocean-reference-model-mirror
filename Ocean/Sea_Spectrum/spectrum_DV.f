c--------------------------------------------------------------------
c  Calcul de la valeur du spectrum de puissance de la mer, pour un 
c  nombre d onde k
c  modele : Durden et Vesecky (1985).
c           description : see spectrum_DV_Base.f
c                       +++++  17 Octobre 2003  +++++
c  * k est en m-1
c--------------------------------------------------------------------


        function S_DV (k)
        
		implicit none

        common /B_/ B_
        common /b0/ b0

		double precision k
        double precision a
        double precision b 
        double precision S_DV 
        double precision B_ ! temp
        double precision b0
		double precision S_DV_Base

        data a/0.225D0/
		data b/1.25D0/
c		data b0/0.004/ ! temp
        

!temp		S_DV = S_DV_Base(k, a, b, b0, b0)
     	S_DV = S_DV_Base(k, a, b, B_, b0)
		return

        end
