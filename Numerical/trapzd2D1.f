        subroutine trapzd2D1(func1, a1,b1,s1,n1,it1)

        implicit none

c        external func1

        double precision a1, b1, s1(1:2), x1, del1, tnm1, sum1(1:2)
        double precision func1_1(1:2), func2_1(1:2)
        integer n1, indice_1, j1
        integer it1
      
        if (n1.eq.1) then
                call func1 (a1, func1_1)
                call func1 (b1, func2_1)
                do indice_1 = 1, 2
                        s1(indice_1) = 0.5D0*(b1-a1)*(func1_1(indice_1) 
     &                          + func2_1(indice_1))
                enddo
                it1 = 1
        else
                tnm1 = it1 
                del1 = (b1 - a1)/tnm1
                x1   = a1 + 0.5D0 * del1
                do indice_1 = 1, 2
                        sum1 (indice_1) = 0.D0
                enddo
                do  j1 = 1, it1
                        call func1 (x1, func1_1)
                        do indice_1 = 1, 2
                                sum1(indice_1) = sum1(indice_1)
     &                                          + func1_1(indice_1)
                        enddo
                        x1 = x1 + del1
                enddo
                do indice_1 = 1, 2
                        s1 (indice_1) = 0.5D0*(s1(indice_1)+
     &                          (b1-a1)*sum1(indice_1)/tnm1)
                enddo
                it1 = 2 * it1 
        endif
        return
        
        end
