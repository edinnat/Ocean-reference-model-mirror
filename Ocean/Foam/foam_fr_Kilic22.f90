! Computation of the foam coverage from Kilic et al., 2022
! Kilic, L., Prigent, C., Jimenez, C. (2022). Development of
! the SURface Fast Emissivity Model for Ocean (SURFEM-Ocean)
! (Report No. NWPSAF-EC_VS-060). EUMETSAT.
!
subroutine foam_fr_Kilic22 (U10, Fr)

implicit none

double precision U10 ! ocean wind speed at 10m in m/s
double precision Fr ! foam fraction

if (U10.le.20) then
    Fr = 6.25D-6*U10**(3)
else
    Fr = 6.25D-6*3.*20**2*(U10-20) + 6.25D-6*20**3
endif


end
