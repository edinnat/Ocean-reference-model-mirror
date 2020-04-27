        subroutine func2D1(x1, y1)

        implicit none

        external dacosd
        external func2D21
        external func2D22

        common /krho1/   krho1 
        common /phi_1/    phi_1 
        common /kd/      kd 
        common /epsi/    epsi
        common /k0/      k0
        common /Z/       Z
        common /Z2/      Z2
        common /C2/      C2
        common /C3/      C3
        common /Di1/     Di1
        common /Di2/     Di2
        common /Dc1/     Dc1
        common /Dc2/     Dc2
        common /krho2/   krho2
        common /kmin/   kmin
        
        double precision dacosd
        double precision y1(1:2)
        double precision x1
        double precision res1(1:2)
        double precision Z
        double precision Z2
        double precision lim
        double precision binf
        double precision bsup
        double precision delt
        double precision kd
        double precision k0
        double precision krho1
        double precision phi_1
        double precision krho2

        double complex   epsi
        double complex   C2
        double complex   C3
        double complex   Di1
        double complex   Di2
        double complex   Dc1
        double complex   Dc2
        
        integer indice
       
        double precision del, binff, bsupp, dres1(1:2)
        integer i, N
        double precision   kmin
        
        krho2 = x1
        do indice = 1, 2
                y1(indice)   = 0.0D0
                res1(indice) = 0.D0
        enddo
c   Calcul des variables ne dépendant pas de phi
c       Z est égal à sin (thetai)
        Z    = x1/k0
        Z2   = Z*Z
        C2   = 1.0D0 - Z2
        C2   = sqrt(C2)
        C3   = sqrt(epsi  - Z2)
c     Partie diffusion incohérente
        Di1  = C2 + C3
        Di2  = epsi*C2 + C3
c     Partie diffusion cohérente
        Dc1  = Z2 + C3*C2
        Dc2  = C3 + C2 
c   Calcul les bornes d'intégration en phi pour que k > kd
        if ((2.0D0*x1).lt.kd) then
	        do indice = 1, 2
                y1(indice) = 0.0D0
c				print*, 'k < kd'
c				print*, y1(indice)
   		     enddo
             return
        endif

c   Intègre la diffusion sur phi entre binf et bsup
                if ((phi_1.eq.0.0D0).or.(phi_1.eq.90.0D0)) then  ! pour les Stokes 1 et 2
                       call func2D21(phi_1 + 180.0D0, res1)
                elseif (phi_1.eq.45.0D0) then ! pour les Stokes 3 et 4
                        phi_1 = 90.0D0
                       call func2D22(phi_1 + 180.0D0, res1)
                else
                        print *, 'Peut pas déterminer quel param Stokes'
                        print *, 'phi_1 = ', phi_1, ' diff de 0, 45 90°'
                        stop
                endif
c   Multiplie le resultat par krho
         
        do indice = 1, 2
                y1(indice) = res1(indice)*krho2
c				print*, 'integ 2 : ', res1(indice), krho2, res1(indice)*krho2
        enddo

        return

        end
