c-----------------------------------------------------------------------
-                        SENSOR PARAMETERS                          -
c-----------------------------------------------------------------------
-- Electromagnetic Wavelength (meter) [e.g. L-band = 0.2123]
0.2381
-- Two-scale cutoff wavenumber(s) (rad/m) [Typical is 2*pi/wavelength/N with N = 3 -> 5] Value 9999 terminates the list
6.6
9999
-- Earth incidence angle(s) (degrees) Value 9999 terminates the list
0.
5.
10.
15.
20.
25.
30.
35.
38.
40.
45.
50.
52.
55.
57.
60.
65.
70.
75.
80.
82.
85.
88.
90.
9999 
c-----------------------------------------------------------------------
-                              WIND                                    -
c-----------------------------------------------------------------------
-- Wind Speed(s) (m/s)  Value 9999 terminates the list
1.
2.
3.
4.
5.
6.
7.
8.
9.
10.
11.
12.
13.
14.
15.
16.
17.
18.
19.
20.
21.
22.
23.
24.
25.
9999
-- Altitude for Wind Speeds provided (m)
10.
-- Model for drag coefficient(y : Cardone (69), c : Charnock (55), d:Donelan(93))
y
c-----------------------------------------------------------------------
-                              OCEAN                                   -
c-----------------------------------------------------------------------
-- Sea Surface Temperature(s), SST (degree Celsius) Value 9999 terminates the list
15.
9999
-- Sea Surface Salinity(ies), SSS (practical salinity scale/unit)  Value 9999 terminates the list
36.
9999
-- Atmospheric Stability (T air - T water)
0.
0.
9999
c-----------------------------------------------------------------------
-                              Output Files                            -
c-----------------------------------------------------------------------
-- Output file that includes foam impact ('none' => no output with foam)
aucun
-- Output file with no foam impact ('none' => no output without foam)
NRCS_Lband_DV2.dat
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
0.84
c-----------------------------------------------------------------------
-                               FOAM                                   -
c-----------------------------------------------------------------------
-- Model for Foam fraction (Monahan (86) : M1, M2 (least square), Monahan & Lu (90): M3 (active + passive),M4 (active), WISE2001 : M5)
M1
-- Model for Foam emissivity (S : Stogryn 72)
s
c-----------------------------------------------------------------------
-                           PERMITTIVITY                               -
c-----------------------------------------------------------------------
-- Model for the sea water dielectric constant (e: Ellison et al. 1998, k: Klein & Swift 1977, m: Meissner et al. (2004, 2012, 2014))
k
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
1
-- lambdamin
0.
-- lambdamax
0.
-- iVal
0.
-- Model Type (2 : 2-scale model, p : small scales only, g : large scales only)
2
-- Number of intergation points on the slope distribution, upwind direction
500
-- Number of intergation points on the slope distribution, crosswind direction
500
-- Upper limit of integration on slope distribution (in number of slope variance)
25.
