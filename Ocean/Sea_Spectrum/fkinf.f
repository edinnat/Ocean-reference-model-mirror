        function fkinf(x)

        implicit none

        
        double precision fkinf
        double precision x
        double precision S_DV

        fkinf = S_DV(x)/x-1.0D-20
        
        end
