       subroutine vent (z1, z2, U_z1, U_z2, u_star)

       implicit none

       double precision z1,z2,U_z1,U_z2,u_star,Z_0, a, b, c, D

       data a/6.84D-05/
       data b/4.28D-03/
       data c/-4.43D-04/
      
       Z_0 = a/u_star + b * u_star**2 + c
       D = - dlog (Z_0)

       U_z2 = U_z1 * (dlog(z2)+D)/(dlog(z1)+D)

       end
