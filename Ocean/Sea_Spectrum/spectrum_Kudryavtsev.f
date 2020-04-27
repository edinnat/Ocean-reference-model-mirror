        subroutine SpectreWKud (k, phi_ang, SpecOut, Nout)

        implicit none
        common /paramKudryavtsev/ paramKudryavtsev
        common /filtersKudry/ filtersKudry
      
        double precision filtersKudry(1:3,1:50003)
        double precision k, phi_ang, SpecOut(1:*),S,paramKudryavtsev(10)
        double precision rho, ck, g, Vw, Va, dep, Cst, const, nc, ng
        double precision a, fetch, kgam, w, c, cg, kres, wres, cres 
        double precision tau, ftch, age, agemax, Cp, kpd, kp, wp, hlow
        double precision fphi, al, si, ga, Cut, Gbig, blow1, Blow
        double precision hlowres, fphires, Cutres, Gres, Blow1res
        double precision Blowres 
        double precision fres, avCbeta, n0, lna, alpha
        double precision lnares, alphares, diss1, diss, diss1res
        double precision zgen, zgenres, angbeta, Cbeta, Cbetares, beta1
        double precision betaphi, beta1res, betaphires, grow ,growres
        double precision Bresrip, Bres, Brip, B, T, U10, ustar, z0
        double precision p_ang, cgres, n0res, dissres, pi, theta
        double precision F, FF, Phi, kmin, kmax, nkfilt
        integer Nout, ikfilt

        pi = 3.141592653589793d0
        rho = 1.D3/1.3d0 ! ratio water/air density
        ck = 0.41D0     ! Von Karman Constant
        g = 9.81D0
        T = 74D-6     ! gravity, surface tension/densite eau de mer 
        Vw = 1.0D-6
        Va = 14D-6    ! viscosity of water and air
        dep=10000     ! depth
        Cst=1.5       
        Cst = Cst/rho ! constant for Cbeta, parameter for wind growth rate Beta
        const = 0.018D0 ! Charnok constant
        nc =1 ! n exponent in (17) in capillarity domain 
        ng=5 ! n exponent in gravity domain
        p_ang=1d0
        a=2.5d-3
        fetch=1d7;

!%%%%%%%%%%%%%%GRID%%%%%%%%%%%%%%%%%%%
        theta = phi_ang*1.745329251994330d-02 ! conv in radian
        kgam = dsqrt(g/T)
        w = dsqrt(g*k + T*k**(3)) ! eq dispersion, T/g = kgam**2
        c = w / k ! phase speed
        cg=(g+3*T*k**(2))/(2*w) ! group speed 
        kres = g/T /k ! kg, gravity waves, eq (14)
        wres = dsqrt(g*kres + T*kres**(3)) ! eq dispersion, T/g = kgam**2
        cres = wres / kres
        cgres=(g+3*T*kres**(2))/(2*wres)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        U10  = paramKudryavtsev(1)
        ustar = paramKudryavtsev(2)
        ! filters
        kmin = filtersKudry(1, 1)
        kmax = filtersKudry(1, 2)
        nkfilt = filtersKudry(1, 3)
        ikfilt = (dlog(k)-dlog(kmin))/((dlog(kmax)-dlog(kmin))
     &/(nkfilt-1))+4 ! 4 is first index after 3 first elts for kmin, kmax, nk
         if (ikfilt .ge. (nkfilt+2)) then
              ikfilt = nkfilt - 2
          endif
          F = (filtersKudry(1,ikfilt)+filtersKudry(1,ikfilt+1))/2d0
          Fres = (filtersKudry(2,ikfilt)+filtersKudry(2,ikfilt+1))/2d0
          Phi = (filtersKudry(3,ikfilt)+filtersKudry(3,ikfilt+1))/2d0
        FF = Phi
        
        z0 = 0.1*Va/ustar+const*ustar**(2)/g !eq (21)
        tau=ustar**(2)
        ftch=9.8*fetch/u10**(2)
        age=2*pi*3.5*ftch**(-0.33)
        age=max(0.83,age)
        !age=0.84*(tanh((ftch/2.2E4)**0.4))**(-0.75);
        agemax=u10*dsqrt(pi/dep/g)
        age=max(age,agemax)


! Frequency Space
        Cp  = u10/age
        kpd = g/(1.2*u10)**2
        kp = g/Cp**2
        wp=g/Cp

                  !%LOW FREQUENCY SPECTRUM

        if (k/kp .lt. .31) then
            hlow = 1.24
        elseif (k/kp .lt. .90) then
            hlow = 2.61*(k/kp)**(.65)
        else
            hlow = 2.28*(k/kp)**(-.65)
        endif


!        fphi = 0.5*(ones(Nphi,1)*hlow)/(dcosh(theta*hlow)**(2))
        fphi = 0.5*hlow/(dcosh(theta*hlow)**(2))

!%DHH-spectrum Donelan, Hamilton, Hui
        al =6E-3 * age**(.55) ! al = rear phace or Phillips and Kitaigorodskii eq. range param. (see Elfouhaily et al. 1997)
        if (age .lt. 5) then
            si = .08 * (1 + 4 * age**(-3))
        else
            si=min(0.24,0.08+(age-5)*(0.24-0.08)/15)
        endif
        ga = 1.7 + 6*dlog10(max(1.d0,age)) ! gamma
        Cut=dexp(-age/dsqrt(10.d0)*max(0.d0,(dsqrt(k/kp)-1)))
        Gbig = dexp(-(dsqrt(k/kp) - 1)**(2) /(2*si**2)) ! Gamma in Elfouhaily et al.
!%Blow1 =0.5* al*dsqrt(k/kp) * dexp(-(kp/k)**(2)) * ga**(G)*Cut;
        Blow1 =0.5* al*dsqrt(k/kp) * dexp(-(kp/k)**(2))*ga**(Gbig)*Cut ! ga**G = Jp, exp(...) ~ LPM
!        Blow = ones(Nphi,1)*Blow1 *fphi
        Blow = Blow1*fphi

!%RESONANT LOW FREQ.SPECTRA:

        if (kres/kp .lt. .31) then
            hlowres = 1.24
        elseif (kres/kp .lt. .90) then
            hlowres = 2.61*(kres/kp)**(.65)
        else
            hlowres = 2.28*(kres/kp)**(-.65)
        endif
!        fphires = 0.5*(ones(Nphi,1)*hlowres)/(dcosh(theta*hlowres)**(2))
        fphires = 0.5*hlowres/(dcosh(theta*hlowres)**(2))

!DHH 
        Cutres=dexp(-age/dsqrt(10d0)*max(0.d0,(dsqrt(kres/kp)-1)))
        Gres = dexp(-(dsqrt(kres/kp) - 1)**(2)/(2*si**2))
        Blow1res=0.5* al*dsqrt(kres/kp) * dexp(-(kp/kres)**2) 
     &*ga**(Gres)*Cutres
!        Blowres = ones(Nphi,1)*Blow1res *fphires
        Blowres = Blow1res*fphires



        avCbeta=0.03 ! Average of Beta


!        n0=ones(Nphi,1)*(1/((1-1/ng)*f+1/ng)); ! eq (22)
        n0=(1/((1-1/ng)*F+1/ng)); ! eq (22)
        lna=dlog(a)-dlog(avCbeta)*1/n0 
        alpha=dexp(lna) ! % eq. (26), a = tuning constant

!        n0res=ones(Nphi,1)*(1/((1-1/ng)*fres+1/ng))
        n0res=(1/((1-1/ng)*fres+1/ng))
        lnares=dlog(a)-dlog(avCbeta)*1/n0res
        alphares=dexp(lnares) ! % eq. (26), a = tuning constant

        diss1 = 4*Vw*k**(2)/ w ! % rate of viscous dissipation
!        diss  = ones(Nphi,1)*diss1 !
        diss  = diss1 !
        diss1res = 4*Vw*kres**(2)/ wres
!        dissres  = ones(Nphi,1)*diss1res
        dissres  = diss1res

        zgen=pi/k
        zgenres=pi/kres
        angbeta=dexp(-theta**(2)/p_ang) ! % eq (20)

        Cbeta=max(1.d0,1/ck*dlog(1/z0*zgen))
        Cbeta=Cst*max(0.d0,Cbeta-c/ustar) ! % eq (19) Kuryatsev et al. 2003

        Cbetares=max(1.d0,1/ck*dlog(1/z0*zgenres))
        Cbetares=Cst*max(0.d0,Cbetares-cres/ustar)
        beta1 = Cbeta*(ustar/c)**(2)
        betaphi = angbeta* beta1 ! % eq (20)
        beta1res = Cbetares*(ustar/cres)**(2)
        betaphires  = angbeta*beta1res
!          grow = (betaphi-diss); 
        grow=angbeta*(beta1-diss1) ! % effective growth rate (just after (13)) BetaNu
! %         growres = (betaphires - dissres);
        growres=angbeta*(beta1res-diss1res)
        growres = max(0.d0, growres)
         
        Bresrip =alphares*(growres)**(1/n0res)
        Bres=Bresrip+Blowres
!  % eq (17), high f spectrum
        Brip = grow **(2)+4 /alpha *growres *Bres *FF; ! Ipc = growres*Bres*FF, BetaNu = grow, FF = Phi filter of parasit cap.
        Brip = ((grow + Brip**(0.5))/2)**(1/n0)
        Brip = alpha *Brip
!% Total Spectrum, eq (30)
!        B=Blow+Brip*(ones(Nphi,1)*(1-Cut))  
        B=Blow+Brip*(1-Cut)  
        SpecOut(1) = B/k**3 ! S
        SpecOut(2) = B      !
        SpecOut(3) = B/k**4 ! W
        SpecOut(4) = n0
        SpecOut(5) = growres
        SpecOut(6) = Bres
        Nout = 6
        end
