      subroutine SuScKud (Su, Sc, kmin, kd)

      implicit none

      external dsind, dcosd
    
      double precision Su, Sc, kmin, kd, dsind, dcosd
      double precision Dphi, dlnk, k, Phi, SpecOut(1:6)
      integer nPhi, nK, iphi, ik, Nout

!                    double precision k, phi_ang, U10, us,SpecOut (10),z0, usn, it, ep
!      double precision param(10), Kudryavtsev, cd10, F, Fres, Phi
!      integer nK, nPhi

      
      nPhi = 36
      nK = 5000
      Dphi = 180d0/(nPhi - 1d0)
      dlnk = (dlog(kd) - dlog(kmin))/(nK - 1d0)
! initialisation
      Su = 0d0
      Sc = 0d0
      do ik = 1, (nK - 1)
        k = kmin*dexp(ik*dlnk)
        do iphi = 0, (nPhi - 1)
          Phi = Dphi*iphi     
          call SpectreWKud(k, Phi, SpecOut, Nout)
          Su = Su + dcosd(Phi)**2*SpecOut(2)
          Sc = Sc + dsind(Phi)**2*SpecOut(2)
        enddo
      enddo
      Su = dsqrt(2*dlnk*Dphi*1.745329251994330d-02*Su)
      Sc = dsqrt(2*dlnk*Dphi*1.745329251994330d-02*Sc)
      end
