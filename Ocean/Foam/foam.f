       subroutine foam (U_10, DeltaT,Fr)

       implicit none
       double precision U_10, DeltaT, Fr

       Fr = 1.95D-05*U_10**2.55D0*dexp(8.61D-02*DeltaT)

       end
