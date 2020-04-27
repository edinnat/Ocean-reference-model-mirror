c_____________________________________________________________
c  Calcul du stress du vent a partir d'un vent donne U_z, a 
c  une altitude donnee z, de la constante de Karman K, de
c  et de Z_0 d'apres Charnock (1955)
c  +++++  26 Mars 1999    +++++
c_____________________________________________________________

c      f est la fonction donnant le stress dont on cherche la
c        racine

 
       function funcUstar_c (u)

       implicit none
       
       common /Uz/ U_z
       common /alt/ alt

       double precision a
       double precision nu
       double precision zc
       double precision zs
       double precision D
       double precision E
       double precision Z_0
       double precision alt
       double precision U_z
       double precision u
       double precision funcUstar_c

       a   = 0.012D0
       nu  = 14.0D-06
       zc  = a*u*u/9.81 
       zs  = 0.11*nu/u
       Z_0 = zc + zs
       D   = 4.D-01 * U_z
       E   = -dlog(alt)
       funcUstar_c = u*(dlog(Z_0)+E) + D 
       return

       end

