c               SPECTRE DIRECTIONNEL D'ELFOUHAILY
c Version avec Omega en parametre

        subroutine W_El (k,phid, PSI, B, delta_ratio)

        implicit none

        common /PI2/ PI2        ! PI2 = 2*pi
        common /U10/ U10        ! U10 = vent à 10m en m/s
        common /ustar/ ustar    ! ustar = stress en m/s
        common /Omega/ Omega    ! inverse de l'age des vagues

        double precision k      ! nombre d'onde en rad/m
        double precision phi    ! phi en radian
        double precision phid   ! phi en degres
        double precision U10    
        double precision ustar
        double precision PSI
        double precision Omega
        double precision B
        double precision delta_ratio
        double precision km
        double precision g
        double precision c
        double precision cm
        double precision cp
        double precision k0
        double precision kp
        double precision sigma
        double precision gamma
        double precision Jp
        double precision Bl
        double precision Bh
        double precision a0
        double precision ap
        double precision am
        double precision PI2 
        double precision alpha_p 
        double precision alpha_m 
        
c unified_spectrum   computes the directional unified spectrum for long
c and and short waves
c
c        input:

c            k:                wavenumber in rd/m
c            phi:               wave direction in rd
c            U10:              wind speed at 10m in m/s
c            ustar:             wind friction velocity in m/s
c            Omega=U10/Cp    the inverse wave-age
c

c         output:

c             PSI                direction wave spectrum
c             B:                 the curvature spectrum
c             delta__ratio:         the spreading factor
c
c By Tanos Elfouhaily March 1996.
        
c        Omega = 0.84D0
        phi  = phid*PI2/360.0D0

        km=363.0D0 ! gravity-capillary peak (rd/m)
        g=9.81D0   ! acceleration due to gravity (m/s/s)

        c=sqrt(g/k*(1.0D0+(k/km)**2)) ! phase speed of short and long waves
        cp=U10/Omega                  ! phase speed of dominant waves
        cm=sqrt(2.0D0*g/km)           ! phase speed at gravity-capillary peak
        k0=g/U10**2                   ! peak wavenumber
        kp=k0*Omega**2                ! wave-age peak wavenumber

c Peak enhancement according to JONSWAP
        sigma=0.08D0*(1.0D0+4.0D0/Omega**3)
        gamma=1.7D0
        if (Omega.gt.1.0D0) then
                gamma=1.7D0+6.0D0*dlog10(Omega)
        endif
        Jp=dexp(-(5.0D0/4.0D0)*(kp/k)**2)
        Jp=Jp*gamma**(dexp(-(sqrt(k/kp)-1.0D0)**2/2.0D0/sigma/sigma))

c Curvature spectrum of Long waves
        alpha_p=6.0D-3*sqrt(Omega)
        Bl=0.5D0*alpha_p*cp/c*dexp(-Omega/sqrt(10.0D0)
     &  *(sqrt(k/kp)-1.0D0))

c Curvature spectrum of Short waves

c        alpha_m=2.0D0*ustar/32.0D0          
c Phillips [1985]
c Jaehne and Riemer [1990] and Hara et al. [1994]
c alpha_m=2*max([1.4*dlog(ustar/cm)+1,4*dlog(ustar/cm)+1])*4*1e-3
        alpha_m=max(0.9D0*dlog(ustar/cm)+1.0D0,
     &          3.0D0*dlog(ustar/cm)+1.0D0)*1.0D-2
c        alpha_m=max(1.0D0*dlog(ustar/cm)+1.0D0,
c     &          3.0D0*dlog(ustar/cm)+1.0D0)*1.0D-2
        Bh=0.5D0*alpha_m*cm/c*dexp(-(k/km-1.0D0)**2/4.0D0)
        B=Jp*(Bl+Bh)               
c curvature spectrum

        a0=dlog(2.0D0)/4.0D0                   
c constant shift of spreading factor
        ap=4.0D0                         
c long wave factor
        am=0.13D0*ustar/cm             
c short wave factor
        delta_ratio=tanh(a0+ap*(c/cp)**2.5+am*(cm/c)**2.5D0)
        PSI=B/k**4/PI2*(1.0D0+delta_ratio*dcos(2.0D0*phi)) 
c directional spectrum

        return

        end
