      Module interp_func

      contains

      subroutine interp(x, y, xi, yi)

      implicit none

      double precision, intent (in) :: x(:), y(:)
      double precision, intent (in) :: xi
      double precision, intent (out) :: yi
      double precision wght
      integer i_t, Ndat

      Ndat = size(x)

      if (xi.lt.x(1)) then

        yi = y(1)
   
      endif

      if (xi.gt.x(Ndat)) then

        yi = y(Ndat)
  
      endif
    
      if ((xi.ge.x(1)).and.(xi.le.x(Ndat))) then

      do i_t = 1, Ndat
        
        if ((xi.ge.x(i_t)).and.(xi.le.x(i_t+1))) then

            wght = (xi - x(i_t))/(x(i_t+1) - x(i_t))
            yi = y(i_t)*(1.D0 - wght) + y(i_t+1)*wght

        exit

      endif

      enddo

      endif 

      end subroutine interp

      end module interp_func
