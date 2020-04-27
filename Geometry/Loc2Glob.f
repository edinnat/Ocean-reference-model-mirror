
      
c-----------------------------------------------------------------------
c	Passage des coordonnées du repere local de la vague vers le repere 
c       global (implemented for theta only TODO implement other angles)
c-----------------------------------------------------------------------

       	subroutine Loc2Glob(thetal, phil, Sx, Sy, thetaG)
	
       	implicit none

       	external dcosd
       	external dsind
       	external dacosd
       	external dasind
c Variables : entree
       	double precision thetal
       	double precision phil
       	double precision Sx
       	double precision Sy
c Variables : sortie
       	double precision thetaG
       	double precision phiG
c TODO       	double precision cosAlpha
c TODO      	double precision sinAlpha
c Variables : intermediaire
       	double precision sq
       	double precision Norm
       	double precision cosTheta
       	double precision sinTheta
       	double precision cosPhi
       	double precision sinPhi
       	double precision cosThetaN
       	double precision sinThetaN
       	double precision cosPhiN
       	double precision sinPhiN
       	double precision cosThetaG
       	double precision sinThetaG
       	double precision cosPhiG
       	double precision sinPhiG
       	double precision Cx
       	double precision Cy
       	double precision cosBeta
       	double precision sinBeta
c Fonction trigo
       	double precision dsind
       	double precision dcosd
       	double precision dacosd
       	double precision dasind

c       Input cartesian coordinates
       	cosTheta = dcosd(thetal)
       	sinTheta = dsind(thetal)
       	cosPhi   = dcosd(phil)
       	sinPhi   = dsind(phil)
       	Cx = sinTheta*cosPhi
       	Cy = sinTheta*sinPhi
c       Constants
       	sq = Sx*Sx + Sy*Sy
       	Norm = dsqrt(sq + 1.0D0)
c       
       	sinThetaN = dsqrt(sq)/Norm
       	cosThetaN = dsqrt(1.0D0 - sinThetaN*sinThetaN)
       	if (sq.ne.0.0D0) then
       		cosPhiN = (-Sx/dsqrt(sq))
       		sinPhiN = (-Sy/dsqrt(sq))
       	else
       		cosPhiN = cosPhi
       		sinPhiN = sinPhi
       	endif
       	cosBeta = 1.0D0/dsqrt(Sx*Sx + 1.0D0)  
       	sinBeta = -Sx/dsqrt(Sx*Sx + 1.0D0)

       	cosThetaG = (-sinBeta*Cx-sinThetaN*sinPhiN*cosBeta*Cy
     &  +cosThetaN*cosTheta)
       	sinThetaG = dsqrt(1.0D0-cosThetaG*cosThetaG) 
       thetaG = dsign(dacosd(cosThetaG),sinThetaG)	
c TODO       	if (sinThetaG.gt.1.0D-04) then
c TODO        		C = (Cy*sqrtcosBeta - cosTheta*sinBeta)/sinThetaG
c TODO        		cosPhiG = C
c TODO             	sinPhiG = (dsqrt(1.0D0-(sinThetaN*sinPhiN)**2)
c TODO      &			*Cy-sinThetaN*sinPhiN*(sinBeta*Cx+cosBeta*
c TODO      &			cosTheta))/sinThetaG
c TODO        	else	
c TODO                	cosPhiG = cosPhiN
c TODO 	      	sinPhiG = sinPhiN
c TODO         endif

       	end
