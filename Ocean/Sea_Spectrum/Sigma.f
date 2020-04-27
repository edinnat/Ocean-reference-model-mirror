c	Calcul de la variance des hauteurs du spectrum à partir de kd
c
c	sigma : variance des hauteurs
c	kd : nombre d'onde de coupure. le spectrum est intégré entre kd
c		et l'infini.
c	cspec : choix du spectrum
c		1 -> Durden & Vesecky
c		2 -> Elfouhaily

	subroutine sigma_ (sigma, kd, cspec)

	implicit none

	external funcSigDV
	external funcSigEl

	double precision sigma
	double precision kd
	double precision EPS
	integer cspec

        parameter (EPS=1.D-06)

c Integration du spectrum des hauteurs "funcSigX" de kd à 1000*kd
c X = KS pour le spectrum de Durden & Vesecky
c X = El pour le spectrum de Elfouhaily
        if (cspec.eq.1) then
                call qromb(funcSigDV,kd, 1000.0D0*kd, sigma, EPS)
        elseif (cspec.eq.2) then
                call qromb(funcSigEl,kd, 1000.0D0*kd, sigma, EPS)
        else
                print *, 'Su : cspec ', cspec
        endif

        sigma = sqrt(sigma)


	end

	function funcSigDV (x)

        implicit none
        
        common /c/ c
        external S_DV


        double precision funcSigDV
        double precision x
        double precision S_DV
        double precision c

        funcSigDV = S_DV(x)

        return

        end

	function funcSigEl (x)

        implicit none
        
        double precision funcSigEl
        double precision x
        double precision Wel
        double precision Bel
        double precision DeltaR

        call W_EL (x, 0.0D0, Wel, Bel, DeltaR)
        funcSigEl = Bel/x/x/x

        return

        end
