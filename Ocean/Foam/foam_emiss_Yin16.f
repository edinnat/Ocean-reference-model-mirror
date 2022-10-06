c
c                                             /--- alphaz
c                        /--- loss_factor ---/
c                       /                    \
c foam_emissYin16 -----/                      \--- thetaz --- betaz --- void_fraction_profile 
c                      \                                                            \   
c                       \                                                            \--- eps_klein_swift
c                        \                                                           
c                         \                                                         
c                          \--- IncohEmissivity_TwoLayer ----------- ReflTransm_PlanarBoundary  
c

c Description: Code computes incoherent emissivity of an inhomogeneous layer
c separated by air on top and a homogeneous medium on the bottom, with
c perfectly smooth parallel boundaries on both sides.

c Input Variables: 
c    %eps2: dielectric constant of middle layer
c    %eps3: dielectric constant of bottom layer
c    %theta_i: incidence angle in air (deg)
c    %a: Single scattering albedo
c    %d: layer thickness (m)
c    %kappa_e: extinction rate through layer (Np/m)
c    
c%Output Products:
c    %e_v_inc(theta_i): v-polarized emissivity
c    %e_h_inc(theta_i): h-polarized emissivity

    

      subroutine foam_Tb_Yin16 ( modl, S, T, f, tht, e_foam )

        implicit none

        ! Inputs
        integer modl ! selector of model parameterization as row number from Table 1 in Yin et al. 2016
        real*8 :: S  ! SSS in psu
        real*8 :: T  ! SST in degC
        real*8 :: f  ! frequency in GHz 
        real*8 :: tht  ! incidence angle in degrees
        ! Outputs
        real*8 , INTENT(OUT) :: e_foam(1:2) ! foam Tb in v and h pol for index 1 and 2 respectively
        ! Intermediate
        real*8 :: ev_foam, eh_foam, Vaf0, Vfw0, m_fm0, h_fe0

        call foam_emiss_Yin16 ( modl, S, T, f, tht, e_foam )

        e_foam(1) = e_foam(1) * ( T + 273.15D0 )
        e_foam(2) = e_foam(2) * ( T + 273.15D0 )

      end subroutine foam_Tb_Yin16

      subroutine foam_emiss_Yin16( modl, S, T, f, tht, e_foam )

      USE mod_foam_emiss

        implicit none

        ! Inputs
        integer modl ! selector of model parameterization as row number from Table 1 in Yin et al. 2016
        real*8 :: S  ! SSS in psu
        real*8 :: T  ! SST in degC
        real*8 :: f  ! frequency in GHz 
        real*8 :: tht  ! incidence angle in degrees
        ! Outputs
        real*8 , INTENT(OUT) :: e_foam(1:2) ! foam emissivity in v and h pol for index 1 and 2 respectively
        ! Intermediate
        real*8 :: ev_foam, eh_foam, Vaf0, Vfw0, m_fm0, h_fe0
        complex*16 :: epsks

        Vfw0 = 0.01D0
        m_fm0 = 1.D0

        modelEmiss: select case (modl)
        case (1) ! model M-Du-E
               Vaf0 = 0.97D0
               h_fe0 = 1.85D0
        case (2) ! model M-Du-E1
               Vaf0 = 0.97D0
               h_fe0 = 1.77D0
        case (3) ! model M-Du-S
               Vaf0 = 0.98D0
               h_fe0 = 1.83D0
        case (4) ! model M-Ku-E
               Vaf0 = 0.94D0
               h_fe0 = 1.99D0
        case (5) ! model M-Ku-S
               Vaf0 = 0.96D0
               h_fe0 = 1.91D0
        case default
                print *, 'Invalid choice for foam emissivity model in'
                print *, 'Yin 2016'
                stop
        end select modelEmiss

        call epsilon_KS(S, T, epsks, f*1.D9)

        call esf_anguelova (tht,f*1.D9,epsks,h_fe0,Vaf0,ev_foam,eh_foam)

        e_foam(1) = ev_foam
        e_foam(2) = eh_foam


      end subroutine foam_emiss_Yin16
