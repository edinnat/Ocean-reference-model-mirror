c   Lit des vecteurs de taille variable dans un fichier

        subroutine readvec (vector, unit_file, size_vector)

        implicit none

        double precision vector(1:100)
        double precision value
        integer          size_vector
        integer          i
        integer          unit_file
        character*80     sautligne

	read (unit_file,'(a)') sautligne
        do i = 1, 101
                read (unit_file,*) value
                if (value.ne.9999) then
                        if (i.ne.101) then
                                vector(i) = value
                        else
                                stop 'plus de 100 elts pour la variable'
                        endif
                else
                        size_vector = i-1
                        exit
                endif
        enddo

        end
