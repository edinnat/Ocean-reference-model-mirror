
      
c-----------------------------------------------------------------------
c	Passage des coordonnées dans le repere global vers le repere 
c       local de la vague
c-----------------------------------------------------------------------

c-----------------------------------------------------------------------
c-----------------------------------------------------------------------

       	subroutine Glob2Loc(thetaG, phiG, Sx, Sy, thetaL, phiL, 
     &			    cosAlpha, sinAlpha)
	
       	implicit none

       	external dcosd
       	external dsind
       	external dacosd
       	external dasind
c Variables : entree
       	double precision thetaG
       	double precision phiG
       	double precision Sx
       	double precision Sy
c Variables : sortie
       	double precision thetaL
       	double precision phiL
       	double precision cosAlpha
       	double precision sinAlpha
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
       	double precision cosThetaL
       	double precision sinThetaL
       	double precision cosPhiL
       	double precision sinPhiL
       	double precision Cx
       	double precision Cy
       	double precision cosBeta
       	double precision sinBeta
       	double precision C
       	double precision Hlx
       	double precision Hly
       	double precision Hx
       	double precision Hy
       	double precision Vx
       	double precision Vy
       	double precision Vz
c Fonction trigo
       	double precision dsind
       	double precision dcosd
       	double precision dacosd
       	double precision dasind

       	cosTheta = dcosd(thetaG)
       	sinTheta = dsind(thetaG)
       	cosPhi   = dcosd(phiG)
       	sinPhi   = dsind(phiG)
       	sq = Sx*Sx + Sy*Sy
       	Norm = dsqrt(sq + 1.0D0)
       	sinThetaN = dsqrt(sq)/Norm
       	cosThetaN = dsqrt(1.0D0 - sinThetaN*sinThetaN)
       	if (sq.ne.0.0D0) then
       		cosPhiN = (-Sx/dsqrt(sq))
       		sinPhiN = (-Sy/dsqrt(sq))
       	else
       		cosPhiN = cosPhi
       		sinPhiN = sinPhi
       	endif
       	Cx = sinTheta*cosPhi
       	Cy = sinTheta*sinPhi
       	cosThetaL = (-Sx*Cx-Sy*Cy+cosTheta)/Norm
       	sinThetaL = dsqrt(1.0D0-cosThetaL*cosThetaL) 
       	
       	cosBeta = 1.0D0/dsqrt(Sx*Sx + 1.0D0)  
       	sinBeta = -Sx/dsqrt(Sx*Sx + 1.0D0)
       	if (sinThetal.gt.1.0D-04) then
       		C = (Cx*cosBeta - cosTheta*sinBeta)/sinThetaL
       		cosPhiL = C
            	sinPhiL = (dsqrt(1.0D0-(sinThetaN*sinPhiN)**2)
     &			*Cy-sinThetaN*sinPhiN*(sinBeta*Cx+cosBeta*
     &			cosTheta))/sinThetaL
       	else	
               	cosPhiL = cosPhiN
	      	sinPhiL = sinPhiN
        endif
	
       	Hlx = sinPhiL 
       	Hly = -cosPhiL 
       	Hx = sinPhi
       	Hy = -cosPhi
       	Vx = -cosTheta*cosPhi
       	Vy = -cosTheta*sinPhi
       	Vz = sinTheta
       	cosAlpha = Hlx*(Hx*cosBeta)+Hly*(Hx*(-sinThetaN*sinPhiN*sinBeta)
     &			+Hy*dsqrt(1.0D0-sinThetaN**2*sinPhiN**2))
       	sinAlpha = Hlx*(Vx*cosBeta-Vz*sinBeta)+Hly*(Vx*(-sinThetaN
     &			*sinPhiN*sinBeta)+Vy*dsqrt(1.0D0-sinPhiN**2
     &			*sinThetaN**2)-Vz*sinThetaN*sinPhiN*cosBeta) 

       	thetaL = dsign(dacosd(cosThetaL),sinThetaL)
c       	phiL = dsign(dacosd(cosPhiL),sinPhiL)
        if (dabs(sinPhiL).lt.1.0D0) then
	       	if (cosPhiL.lt.0.0D0) then
        	     phiL = 180.0D0 - dasind(sinPhiL)
	        else
        	     phiL = dasind(sinPhiL)
	        endif
        elseif (dabs(sinPhiL).lt.1.000001D0) then
               phiL = dsign(90.0D0,sinPhiL)
        else
               print*, 'Erreur Glob2Loc.f'
               print*, 'Sin(phil) > 1.000001'
               stop
        endif
c        print*, thetal, phil, cosalpha, sinalpha, 
c     &          cosalpha**2 + sinalpha**2
c        write(*,'(7(1x,e13.5))') 
c     &      sinThetaL, thetaL,sinPhiN, thetaG, phiG, cosAlpha, sinAlpha 

       	end
