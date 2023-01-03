! Compute foam emissivity and brightness temperature using Kilic et al., 2022, parameterization
! based on the work of Yin et al., 2016, the foam emissivity model of Anguelova and Gaiser, 2013.

  subroutine foam_Tb_Kilic22 (perm, T, f, tht, tb_foam )

    implicit none

    ! Inputs
    REAL (KIND=8) :: T  ! SST in degC
    REAL (KIND=8) :: f  ! frequency in GHz
    REAL (KIND=8) :: tht  ! incidence angle in degrees
    COMPLEX (KIND=8) :: perm ! permittivity of sea water (usually from Meissner and Wentz, 2012 for this foam emissivity module)
    ! Outputs
    REAL (KIND=8) , INTENT(OUT) :: tb_foam(1:2) ! foam Tb in v and h pol for index 1 and 2 respectively
    ! Intermediate
    REAL (KIND=8) :: ev_foam, eh_foam, Vaf, h_fe
    REAL (KIND=8) :: e_foam(1:2)

    call foam_emiss_Kilic22 (perm, f, tht, e_foam )

    tb_foam(1) = e_foam(1) * ( T + 273.15D0 )
    tb_foam(2) = e_foam(2) * ( T + 273.15D0 )

  end subroutine foam_Tb_Kilic22

  subroutine foam_emiss_Kilic22 (perm, f, tht, e_foam )

  USE mod_foam_emiss

    implicit none

    ! Inputs
    REAL (KIND=8) :: f  ! frequency in GHz
    REAL (KIND=8) :: tht  ! incidence angle in degrees
    COMPLEX (KIND=8) :: perm
    ! Outputs
    REAL (KIND=8) , INTENT(OUT) :: e_foam(1:2) ! foam emissivity in v and h pol for index 1 and 2 respectively
    ! Intermediate
    REAL (KIND=8) :: ev_foam, eh_foam, Vaf, h_fe

    Vaf = 0.97D0
    h_fe = 0.2D0 !in centimeters


    call esf_anguelova (tht,f,perm,h_fe,Vaf,ev_foam,eh_foam)
    e_foam(1) = ev_foam
    e_foam(2) = eh_foam


  end subroutine foam_emiss_Kilic22
