c-----------------------------------------------------------------------
-                        SENSOR PARAMETERS                          -
c-----------------------------------------------------------------------
-- Electromagnetic Wavelength (meter) [e.g. L-band = 0.2123]
0.2123
-- Two-scale cutoff wavenumber(s) (rad/m) [Typical is 2*pi/wavelength/N with N = 3 -> 5]
7.4842
9999
-- Earth incidence angle(s) (degrees)
0.
5.
10.
15.
20.
25.
30.
35.
40.
45.
50.
55.
60.
9999 
c-----------------------------------------------------------------------
-                              WIND                                    -
c-----------------------------------------------------------------------
-- Wind Speed(s) (m/s)
3.
7.
12.
9999
-- Altitude for Wind Speeds provided (m)
10.
-- Model for drag coefficient(y : Cardone (69), c : Charnock (55), d:Donelan(93))
y
c-----------------------------------------------------------------------
-                              OCEAN                                   -
c-----------------------------------------------------------------------
-- Sea Surface Temperatures, SST (degree Celsius)
15.
9999
-- Sea Surface Salinity, SSS (practical salinity scale/unit)
36.
9999
-- Atmospheric Stability (T air - T water)
0.
9999
c-----------------------------------------------------------------------
-                              Output Files                            -
c-----------------------------------------------------------------------
-- Output file that includes foam impact ('none' => no output with foam)
test_foam.dat
-- Output file with no foam impact ('none' => no output without foam)
test_nofoam.dat
c-----------------------------------------------------------------------
-                               SEA STATE                              -
c-----------------------------------------------------------------------
-- Slope variance computation (c : Cox & Munk 55, s: sea spectrum)
s
-- Amplitude  coefficient for sea spectrum by Durden & Vesecky 1985 (e.g. 0.008 Yueh ; 0.004 original DV1985)
0.008
-- Sea Spectrum Model (e: Elfouhaily, d:Durden & Vesecky)
d
-- Inverse wave age
1.6
c-----------------------------------------------------------------------
-                               FOAM                                   -
c-----------------------------------------------------------------------
-- Model for Foam fraction (Monahan (86) : M1, M2 (least square); Monahan & Lu (90): M3 (active + passive), M4 (active); WISE2001 : M5; Yin et al. (2016), 5 versions: M-Du-E, M-Du-E1, M-Du-S, M-Ku-E, M-Ku-S)
M1
-- Model for Foam emissivity (Stogryn 72 : S; Yin et al. (2016), 5 versions: M-Du-E, M-Du-E1, M-Du-S, M-Ku-E, M-Ku-S; Anguelova et al. 2022 (multi-freq tuned): M-Du-Tune)
M-Du-Tune
c-----------------------------------------------------------------------
-                           PERMITTIVITY                               -
c-----------------------------------------------------------------------
-- Model for the sea water dielectric constant (e: Ellison et al. 1998, k: Klein & Swift 1977, m: Meissner et al. (2004, 2012, 2014), h: high frequency from tabulated data, 28.8-449677 GHz, see epsilon_hifreq.f)
m
-- Imaginary part of diecltric constant is (O: positive; N: negative)
o
c-----------------------------------------------------------------------
-                             SWELL                                    -
c-----------------------------------------------------------------------
-- Swell is prrsent (O: yes; N: no)
N
-- Number of points to compute integral of swell along the wind direction
100
-- Number of points to compute integral of swell across the wind direction
100
-- RMS for swell height (m)
1.
-- Half power width of swell Gaussian PDF along wind (1/m)
0.0025
-- Half power width of swell Gaussian PDF across wind (1/m)
0.0025
-- Peak PDF swell along wind (1/m)
0.
-- Peak PDF swell across wind (1/m)
0.0314
c-----------------------------------------------------------------------
-            PARAMETERS FOR SCATTERING LOOKUP TABLE                    -
c-----------------------------------------------------------------------
-- Maximum value for incidence angle and incremental step (degrees)
90. 
10.
-- lambdamin
0.
-- lambdamax
0.
-- iVal
0.
-- Model Type (2 : 2-scale model, p : small scales only, g : large scales only)
2
-- Number of intergation points on the slope distribution, upwind direction
200
-- Number of intergation points on the slope distribution, crosswind direction
200
-- Upper limit of integration on slope distribution (in number of slope variance)
5.
c-----------------------------------------------------------------------
-               COMPONENT FOR MULTI RELFECTION                         -
c-----------------------------------------------------------------------
--   Model to use (0 : Skip computation ; 1 : Masuda 2006 )
0
--   Orders of multi-reflection to run
2 
