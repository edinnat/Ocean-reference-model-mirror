c--------------------------------------------------------------------
c  Calcul de la valeur du spectrum de puissance de la mer, pour un 
c  nombre d'onde k 
c  model : loi de puissance en k-3
c                       +++++  25 Mars 1999  +++++
c  * k est en m-1
c !!!!!!!! U10 = 10 m/s
c--------------------------------------------------------------------


        function S_powerLaw (k)

	implicit none

	double precision k
	double precision S_powerLaw 

        if ((k.lt.1.0D0).or.(k.gt.125.0D0)) then
                S_powerLaw = 0.0D0
        else
	        S_powerLaw = 5.0D-03/k/k/k
        endif
        return

        end

