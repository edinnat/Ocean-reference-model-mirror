c       Calcul les valeurs extreme contenues dans un vecteur

        subroutine extr_val (vector, size_vector, min_val, max_val)

        implicit none

        double precision vector(1:100)
        double precision min_val
        double precision max_val
        integer          size_vector
        integer          i

        min_val = vector(1)
        max_val = vector(1)
        
        do i=2, size_vector
                if (vector(i).lt.min_val) min_val = vector(i)
                if (vector(i).gt.max_val) max_val = vector(i)
        enddo
                
        end
