       subroutine dichoto (x1, x2, eps, func, root)

       implicit none

       double precision b_inf, b_sup, eps, root, func,x1,x2

       b_inf = x1
       b_sup = x2
       if ((func(b_inf)*func(b_sup)).gt.0.) goto 2
    1  root = (b_inf+b_sup)/2.D0
       if (((b_sup-b_inf)/root).le.eps) goto 2
       if ((func(b_inf)*func(root)).le.0.) then
          b_sup = root
       else
          b_inf = root
       endif 
       goto 1
    2  end
