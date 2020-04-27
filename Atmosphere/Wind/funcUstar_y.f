c_____________________________________________________________
c  Calcul du stress du vent a partir d'un vent donne U_z, a 
c  une altitude donnee z, de la constante de Karman K, de
c  et de Z_0 d'apres Yueh (1997), expressions (59) et (60).
c  +++++  26 Mars 1999    +++++
c_____________________________________________________________

c      f est la fonction donnant le stress dont on cherche la
c        racine

 
       function funcUstar_y (u)

       implicit none
       
       common /Uz/ U_z
       common /alt/ alt

       double precision a
       double precision b
       double precision c
       double precision D
       double precision E
       double precision Z_0
       double precision alt
       double precision U_z
       double precision u
       double precision funcUstar_y

       data a/6.84D-05/
       data b/4.28D-03/
       data c/-4.43D-04/
       D = 4.D-01 * U_z
       E = -dlog(alt)
       Z_0 = a/u + b * u*u + c
       funcUstar_y = u*(dlog(Z_0)+E) + D 
       return

       end

