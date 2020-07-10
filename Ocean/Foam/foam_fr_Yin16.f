c Compute foam fraction from Yin et al. 2016
c
c Input: 
c
c           U: Wind speed (m/s) WS
c           mod: character to define model to be used
c                Format is 'm_XX_YY'
c                XX is:
c                      Du for sea spectrum by Durden & Vesecky (1985)
c                      Ku for sea spectrum by Kudryavtsev et al. (2003)
c                YY is:
c                      E   collocated with ECMWF WS
c                      E1  collocated with ECMWF WS and model optimized for 20-50 deg incidence angles
c                      S   collocated with SSMIS WS
c                
c
c Output:
c           W : Foam fraction (percent)
      
      subroutine foam_fr_Yin16 (U, modl, W)

      implicit none

      double precision U, W, b, c
      character(*) modl

      select case ( modl )
      
      case ('M-Du-E')
       
       b = 1.12D-6 
       c = 3.15D0 

      case ('M-Du-E1')
        
       b = 1.29D-6 
       c = 3.10D0 

      case ('M-Du-S')
        
       b = 3.83D-6 ;
       c = 2.76D0 

      case ('M-Ku-E')
        
       b = 6.48D-7 
       c = 3.40D0 

      case ('M-Ku-S')
        
       b = 1.62D-6 
       c = 3.12D0

      case  default

          stop 'error'

      end select
      
      W = b*U**c 

      end
