      subroutine epsilon_hifreq (SST, SSS, epsi_hifreq, f)

!---------------------------------------------------------------------
! Computes the dielectric constant of seawater from tabulated data in
! the range 0.962-14999.6 inverse cm (frequency 28.8-449677 GHz)
!
!       Inputs:
!               SST : Sea surface temperature in degrees C
!               SSS : Sea surface salinity in psu
!               f   : frequency in Hz
!
!       Outputs:
!               epsi_hifreq : complex dielectric constant of sea water
!
!---------------------------------------------------------------------
! Tabulated data
!
! Data are read for the complex refractive index of seawater in the 
! range 0.962-14999.6 cm^-1 (28.8-449677 GHz). Pure water refractive 
! indices are based primarily on data at 273 K and 298 K from
! Rowe et al., J. Geophys. Res. 125 (2020) 10.1029/2020JD032624.
! Data from Pinkley and Williams, J. Opt. Soc. Am., 66, 554-558 (1976)
! are used to supply a correction to pure water values due to ocean
! salinity. The salinity correction is interpolated onto the
! Rowe et al. spectral grid - the correction is only considered
! reliable in the range 500-5000 cm^-1 and is elsewhere set to zero.
! Newman et al., Q. J. R. Meteorol. Soc., 131, 2539-2557 (2005) data
! (with temperature dependence) are used for the interval
! 770-1230 cm^-1. Newman et al. seawater data have been converted to
! pure water equivalents by subtracting the correction of
! Pinkley and Williams - these values are merged with the Rowe et al.
! data outside this range over an interval of 250 cm^-1.
! The refractive index is modelled to vary linearly with temperature.
! For water temperature T in degrees C,
!   pure water (real) n(T) = n_T0 + (T - 0.0)*coef_n
!   pure water (imag) k(T) = k_T0 + (T - 0.0)*coef_k
! where n_T0, k_T0 are the real and imaginary parts of the refractive
! index at 0.0 degrees C (273.15 K) and coef_n, coef_k are coefficients
! describing the linear change per 1.0 degree C change in temperature.
! The salinity correction of Pinkley and Williams is added to pure
! water values at temperature T:
!   salt water (real) n(T,SSS) = n(T) + n_salt_corr*(SSS/SSS_standard)
!   salt water (imag) k(T,SSS) = k(T) + k_salt_corr*(SSS/SSS_standard)
! where the sea surface salinity (SSS) is in psu (g/kg). Pinkley and
! Williams (1976) used standard sea water according to the prescription
! of Lyman and Fleming, J. Marine Res., 3, 134 (1940). Table VII in
! Lyman and Fleming gives a sum SSS_standard = 34.482 g/kg by which SSS
! is scaled in this routine.
! Newman et al. investigated how the difference in saltwater versus
! pure water emissivity varied with temperature and found no significant
! temperature dependence of the salinity correction in the infrared.
! Table columns, varying with wavenumber f_wn, are:
!  f_wn (cm^-1), n_T0, k_T0, coef_n, coef_k, n_salt_corr, k_salt_corr
!
!---------------------------------------------------------------------
! Computation
!
! For input SST and SSS, the temperature- and salinity-dependent
! complex refractive index is computed using the relations above.
! The complex refractive index is interpolated to input frequency f.
! For output, there is a conversion to complex dielectric constant.
!---------------------------------------------------------------------

      implicit none

! Inputs
      double precision, intent(in) :: SST, SSS, f
! Outputs
      double complex, intent(out) :: epsi_hifreq

      character*300       :: hiresfile, header_text, Fort
      
      integer            :: stat, i
      integer, parameter :: nwn = 15556

      double precision            :: table_wn(nwn)
      double precision            :: table_n0(nwn)
      double precision            :: table_k0(nwn)
      double precision            :: table_cn(nwn)
      double precision            :: table_ck(nwn)
      double precision            :: table_ncorr(nwn)
      double precision            :: table_kcorr(nwn)
      double precision            :: n_pure_water_SST(nwn)
      double precision            :: k_pure_water_SST(nwn)
      double precision            :: n_salt_water_SST_SSS(nwn)
      double precision            :: k_salt_water_SST_SSS(nwn)
      double precision            :: f_wn
      double precision, parameter :: c_light_cm_per_s = 29979245800.0
      integer                     :: index
      double precision            :: n_interp, k_interp, weight
      double complex              :: nk_complex
!
! Pinkley and Williams (1976) used standard sea water according to the
! prescription of Lyman and Fleming, J. Marine Res., 3, 134 (1940).
! Table VII in Lyman and Fleming gives a sum 34.482 g/kg
!
      double precision, parameter :: pw_salinity = 34.482  ! psu
!
! Read in tabulated refractive index data.
! Table columns, varying with wavenumber f_wn, are:
!  f_wn (cm^-1), n_T0, k_T0, coef_n, coef_k, n_salt_corr, k_salt_corr
!
      call getenv('OceanMod', Fort)
      hiresfile = Fort(1:lnblnk(Fort))//
     &'/Data/Refractive_Index/refindex_hifreq.dat'
      open(unit=42, file=hiresfile, status='old',
     &     action='read', iostat=stat)
      if (stat /= 0) stop 'Could not open file : '//hiresfile
! Skip header
      do i = 1, 25
        read (42,*) header_text
      end do
! Read table columns
      do i = 1, nwn
        read (42,*) table_wn(i), table_n0(i), table_k0(i),
     &    table_cn(i), table_ck(i), table_ncorr(i), table_kcorr(i)
      end do
      close(42)
!
! Calculate pure water refractive index (real, imag parts) for input SST
! refractive index (real) n = n_T0 + (T - 0.0)*coef_n for n_T0 at 0.0 degC
! refractive index (imag) k = k_T0 + (T - 0.0)*coef_k for k_T0 at 0.0 degC
!
      n_pure_water_SST(:) = table_n0 + SST*table_cn
      k_pure_water_SST(:) = table_k0 + SST*table_ck
!
! Apply salt water correction to real, imag parts of pure water refractive index
! refractive index (real) n = n_pure_water + n_salt_corr for input SSS
! refractive index (imag) k = k_pure_water + k_salt_corr for input SSS
!
      n_salt_water_SST_SSS(:) = n_pure_water_SST(:)
     &                        + table_ncorr*(SSS/pw_salinity)
      k_salt_water_SST_SSS(:) = k_pure_water_SST(:)
     &                        + table_kcorr*(SSS/pw_salinity)
!
! Determine complex refractive index at frequency f by interpolation -
! if f is outside the bounds of the tabulated data then abort
!
      f_wn = f/c_light_cm_per_s  ! wavenumbers in inverse cm
      if ((f_wn .lt. table_wn(1)) .or.
     &    (f_wn .gt. table_wn(nwn))) then
        print *, 'Supplied frequency is outside range of dataset: ',
     &             f_wn, 'cm^-1, ', f, 'Hz'
        print *, ''
        print *, 'Allowed range of frequency when using the',
     &           ' high-frequency dielectric model is: '
        print *, '         0.962 cm^-1 < - > 14999.6 cm^-1' 
        print *, '          28.8 GHz   < - > 449677  GHz'
        
        print *, 'Stopping'
        stop
      endif
      do index = 1, nwn-1
        if ((f_wn .ge. table_wn(index)) .and.
     &      (f_wn .le. table_wn(index+1))) then
          weight = (f_wn - table_wn(index))
     &           /(table_wn(index+1) - table_wn(index))
          n_interp = (1.0 - weight)*n_salt_water_SST_SSS(index)
     &             + weight*n_salt_water_SST_SSS(index+1)
          k_interp = (1.0 - weight)*k_salt_water_SST_SSS(index)
     &             + weight*k_salt_water_SST_SSS(index+1)
        endif
      end do

      nk_complex = dcmplx(n_interp, k_interp)
!
! Convert complex refractive index to complex dielectric constant
!
      epsi_hifreq = nk_complex**2

      end subroutine epsilon_hifreq


