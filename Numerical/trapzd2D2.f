        subroutine trapzd2D2(func2, a2,b2,s2,n2,it2)

        implicit none

c        external func2
        double precision a2, b2, s2(1:2), x2, del2, tnm2, sum2(1:2)
        double precision func1_2(1:2), func2_2(1:2)
        integer n2, indice_2, j2
        integer it2
      
        if (n2.eq.1) then
                call func2 (a2, func1_2)
                call func2 (b2, func2_2)
                do indice_2 = 1, 2
                        s2(indice_2) = 0.5D0*(b2-a2)*(func2_2(indice_2) 
     &                          + func2_2(indice_2))
                enddo
                it2 = 1
        else
                tnm2 = it2 
                del2 = (b2 - a2)/tnm2
                x2   = a2 + 0.5D0 * del2
                do indice_2 = 1, 2
                        sum2 (indice_2) = 0.D0
                enddo
                do  j2 = 1, it2
                        call func2 (x2, func1_2)
                        do indice_2 = 1, 2
                                sum2(indice_2) = sum2(indice_2)
     &                                          + func1_2(indice_2)
                        enddo
                        x2 = x2 + del2
                enddo
                do indice_2 = 1, 2
                        s2 (indice_2) = 0.5D0*(s2(indice_2)+
     &                          (b2-a2)*sum2(indice_2)/tnm2)
                enddo
                it2 = 2 * it2 
        endif
        return
        
        end
