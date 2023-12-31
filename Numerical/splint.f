c----   Subroutine splint  takes the data points, xa and ya vectors of length n, and the second
c       derivatives computed by spline, y2a, and it computes the y value for the specified x value.
c
        subroutine splint(xa,ya,y2a,n,x,y)
        dimension xa(n),ya(n),y2a(n)
        klo=1
        khi=n
  1     if (khi-klo.gt.1) then
           k=(khi+klo)/2
           if(xa(k).gt.x)then
              khi=k
           else
              klo=k
           endif
           goto 1
        endif
        h=xa(khi)-xa(klo)
        if (h.eq.0.) then
            write(*,*) 'Error from splint.f: bad xa input.'
            stop
        endif
        a=(xa(khi)-x)/h
        b=(x-xa(klo))/h
        y=a*ya(klo)+b*ya(khi)+
     &     ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
        return
        end
