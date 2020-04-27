       subroutine P (Sx, Sy, P_Sx_Sy, Su, Sc)

       implicit none

       double precision Sx
       double precision Sy
       double precision P_Sx_Sy
       double precision Su
       double precision Sc

       P_Sx_Sy = 1.5915494D-01/Su/Sc*dexp(- 5.D-01*(Sx*Sx/Su/Su+
     &           Sy*Sy/Sc/Sc))

       end
