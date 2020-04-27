        subroutine root (b_inf,b_sup,eps,func,root_)

        implicit none

        external func
        double precision b_inf, b_sup, eps, root_, func
 
        root_ = 1.0D04
        call dichoto (b_inf, b_sup, eps, func, root_)
        
        end
