c----------TEMPLATE TO BE READ  AND MODIFIEDBY Tb_rough_vs_freq.m ------
-                        SENSOR PARAMETERS                          -
c-----------------------------------------------------------------------
-- Electromagnetic Wavelength (meter) [e.g. L-band = 0.2123]
XXlambda
-- Two-scale cutoff wavenumber(s) (rad/m) [Typical is 2*pi/wavelength/N with N = 3 -> 5]
XXKd
9999
-- Earth incidence angle(s) (degrees)
XXtht
9999 
c-----------------------------------------------------------------------
-                              WIND                                    -
c-----------------------------------------------------------------------
-- Wind Speed(s) (m/s)
XXU10
9999
-- Altitude for Wind Speeds provided (m)
XXUalt
-- Model for drag coefficient(y : Cardone (69), c : Charnock (55), d:Donelan(93))
y
c-----------------------------------------------------------------------
-                              OCEAN                                   -
c-----------------------------------------------------------------------
-- Sea Surface Temperatures, SST (degree Celsius)
XXSST
9999
-- Sea Surface Salinity, SSS (practical salinity scale/unit)
XXSSS
9999
-- Atmospheric Stability (T air - T water)
XXstab
9999
c-----------------------------------------------------------------------
-                              Output Files                            -
c-----------------------------------------------------------------------
-- Output file that includes foam impact ('none' => no output with foam)
XXFout1
-- Output file with no foam impact ('none' => no output without foam)
XXFout2
c-----------------------------------------------------------------------
-                               SEA STATE                              -
c-----------------------------------------------------------------------
-- Slope variance computation (c : Cox & Munk 55, s: sea spectrum)
s
-- Amplitude  coefficient for sea spectrum by Durden & Vesecky 1985 (e.g. 0.008 Yueh ; 0.004 original DV1985)
XXSpecAmp
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
-                           PERMITTIVITE                               -
c-----------------------------------------------------------------------
-- Modele de constante dielectrique (e: ellison, k:klein & swift)
k
-- Partie imaginaire de constante dielectrique > 0 (O/N)
o
c-----------------------------------------------------------------------
-                             HOULE                                    -
c-----------------------------------------------------------------------
-- Presence de houle (O/N)
N
-- Nombre de point de tabulation de la houle (direction du vent)
100
-- Nombre de point de tabulation de la houle (normale au vent)
100
-- RMS des hauteurs de la houle (m)
1.
-- Largeur à mi puissance de la houle (direction du vent) (m-1)
0.0025
-- Largeur à mi puissance de la houle (normale du vent) (m-1)
0.0025
-- pic densite spectrale houle (direction du vent) (m-1)
0.
-- pic densite spectrale houle (normale du vent) (m-1)
0.0314
c-----------------------------------------------------------------------
-            PARAMETERS FOR SCATTERING LOOKUP TABLE                    -
c-----------------------------------------------------------------------
-- Maximum value for incidence angle and incremental step (degrees)
90. 
5.
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
