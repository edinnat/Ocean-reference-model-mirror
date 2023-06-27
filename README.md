# PARMIO: A Reference Ocean Model for Emissivity and Reflectivity

The Passive and Active Reference Microwave to Infrared Ocean (PARMIO) model is a modular computer model coded in FORTRAN that computes the ocean surface emissivity and reflective properties from the microwave (~0.5 GHz) to the infrared (~400 nm). It was created to answer the need, identified in various reports and international workshops, for a community reference quality ocean emission and reflection model for use across a broad spectral range, as well as supporting passive and active remote sensing.

More information on the project and original team members can be found at the International Space Science Institute : https://www.issibern.ch/teams/oceansurfemiss/

Model details and validation effort are reported in

E. Dinnat et al., “PARMIO: A reference quality model for ocean surface emissivity and backscatter from the microwave to the infrared,” Bull. Am. Meteorol. Soc., vol. 104, no. 4, pp. E742–E748, Apr. 2023. https://journals.ametsoc.org/view/journals/bams/aop/BAMS-D-23-0023.1/BAMS-D-23-0023.1.xml


# Compiling the source code

The source code comes with a 'Makefile' located in the code's root directory that provides all the information needed for compiling exectuables. Also provided is a very basic shell script that uses the Makefile to compile the radiometer code and produces a "Tb" executable. The shell script is also in the root directory and named 'compile_Tb.sh'. Users may need to change the compiler definition in the Makefile to fit their installation. By default, environmental variables are used to define various compilation parameters (see below), but users can enter directly what they want to use as compiler and paths into the Makefile. 

Default environment variables used for compilation:

"compil" should point to the user's compiler (e.g. gfortran, using full path if necessary)
"OceanMod" should point to the root folder of the ocean model code

Both these variables can be redefined in the Makefile without defining the environment variables.

Three executables can be created:
* Tb 	: The main radiometer code that computes the brightness temperature of the ocean surface. The computation relies on azimuthal harmonic decomposition for the angle between the surface emission vector and the wind direction.
* TbTot (Obsolete): Same as previous code except the simulation is ran for mutliple azimuth angles and the  . This code is much slower than 'Tb'. It also has not be updated as often as 'Tb', and is most of the time far behind in terms of features and bug fixes.
* NRCS  : The code to compute the Normalized Radar Cross Section of the ocean surface.

Each exectuable above can be created from the root directory issuing one of the following  commands:

* make -f Makefile Tb
* make -f Makefile TbTot
* make -f Makefile NRCS

After compilation, the executable just created is placed in the 'Binaries' folder.

# Running the program

The program is run from one of the executables in the "Binaries" folder. The program will first look for a configuration (or parameters) file and some input data files at default paths. To define the paths properly for the user's installation, the user should define the following environment variables:


* inputFortran : where to find the configuration files that contain parameters to run the program, such as SST, SSS, wind speed, frequency .... 

* outputFortran: where to save the results.

* OceanMod : where the ocean emissivity/reflectivity model root is (the path where this README.md file is most likely located).

Users can override some default setting by passing the full path to the configuration file as a parameter to the executable. For exemple, from the root folder of the code, users can run the command:

	Binaries/Tb Config/Tb.p

 which will use the configuration file in "Config" folder named "Tb.p". If no configuration file is provided, the program will use the default 'Tb.p' from the default path in "inputFortran". 

The program will save the results in the folder defined by "outputFortran" with a filename defined in the configuration file. It can output 2 files, with or without foam effects included. 

# Format of output files
The Tb output files (both with and without foam)  will contain the following variables in columns (first number is column number):
 
01: freq, frequency in GHz \
02: SST, Sea surface temperature in degree C\
03: SSS, Sea Surface Salinity in psu\
04: U10, wind speed at 10 meter height in meter per second\
05: ustar, wind stress in cm/s, computed in program from U10 and drag coefficient model\
06: theta, incidence angle in degrees\
07: TvN, flat surface (Null wind speed) Tb at V-pol in Kelvin\
08: ThN, flat surface (Null wind speed) Tb at H-pol in Kelvin \
09: Tv0 , rough surface contribution to Tb at V-pol in Kelvin, omnidirectional component ("0-th harmonic")\
10: Th0, rough surface contribution to Tb at H-pol in Kelvin, omnidirectional component ("0-th harmonic")\
11: Tv1, first harmonic (upwind/downwind) of roughness contribution to Tb at V-pol in Kelvin\
12: Th1, first harmonic (upwind/downwind) of roughness contribution to Tb at H-pol in Kelvin\
13: U1, first harmonic (upwind/downwind) of roughness contribution to Third Stokes parameter in Kelvin\
14: V1, first harmonic (upwind/downwind) of roughness contribution to Forth Stokes parameter in Kelvin\
15: Tv2, second harmonic (crosswind) of roughness contribution to Tb at V-pol in Kelvin\
16: Th2, second harmonic (crosswind) of roughness contribution to Tb at H-pol in Kelvin\
17: U2, second harmonic (crosswind) of roughness contribution to Third Stokes parameter in Kelvin\
18: V2, second harmonic (crosswind) of roughness contribution to Forth Stokes parameter in Kelvin \
19: lam_d, cutoff wavelength between small and large scale roughness for the two-scale model (meters)\
20: Stab, atmospheric stability (degC)\
21: Re(epsi), real part of sea water dielectric constant\
22: Im(epsi), imaginary part of sea water dielectric constant\
23: Var_u, large-scale roughness slope variance in upwind direction\
24: Var_c, large-scale roughness slope variance in downwind direction \
25: Fr, foam fraction in percent (always 0 for the output file without foam)
26: dTb_v_MR: contribution to Tb in V-pol from  multiple reflections
27: dTb_h_MR: contribution to Tb in H-pol from  multiple reflections

The total Tb for a direction phi between propagation vector and wind vector is computed from:

	Tb = TpN + Tp0 + Tp1 * cos (phi) + Tp2 * cos ( 2*phi ) + dTb_p_MR

for p = v or h pol, and from:

	p = p1 * sin (phi) + p2 * sin ( 2*phi )

for p = U or V (for the 3rd and 4th Stokes parameters respectively).

Some intermediate results for the brightness temperatures are outputed to a file with the name derived from the "no foam" output file's name with the "Int_" prefix added :

01: freq, frequency in GHz \
02: SST, Sea surface temperature in degree C\
03: SSS, Sea Surface Salinity in psu\
04: U10, wind speed at 10 meter height in meter per second\
05: ustar, wind stress in cm/s, computed in program from U10 and drag coefficient model\
06: theta, incidence angle in degrees\
07-10 : Tv0_SPM, Th0_SPM, T30_SPM, T40_SPM = Small Perturbation Method Tb for v, h, Stokes 3 and Stokes 4 respectively, omnidirectional term.\
11-14 : Tv2_SPM, Th2_SPM, T32_SPM, T42_SPM = Small Perturbation Method Tb for v, h, Stokes 3 and Stokes 4 respectively, second harmonic in azimuth wrt wind direction.\
15-16 : Tv0_GO, Th0_GO = Geometric Optics Tb in v- and h-pol respectively, omnidirectional term.\
17-20 : Tv1_GO, Th1_GO, T31_GO, T41_GO  = Geometric Optics Tb in v-, h-pol, Stokes 3 and Stokes 4 respectively, first harmonic (upwind/downwind).\
21-24 : Tv2_GO, Th2_GO, T32_GO, T42_GO = Geometric Optics Tb in v-pol, h-pol, Stokes 3 and Stokes 4 respectively, second harmonic (crosswind).\

 
The NRCS output file contains:

01 : freq, frequency in GHz \
02 : SST, Sea surface temperature in degree C\
03 : SSS, Sea Surface Salinity in psu\
04 : U10, wind speed at 10 meter height in meter per second\
05 : ustar, wind stress in cm/s, computed in program from U10 and drag coefficient model\
06 : theta, incidence angle in degrees\
07 : sigmaVV0 , cross section from Bragg scattering on tilted large waves, omnidirectional\
08 : sigmaHH0 , cross section from Bragg scattering on tilted large waves, omnidirectional\
09 : sigmaVV1 , cross section from Bragg scattering on tilted large waves, 1st Harmonic (upwind/downwind)\
10 : sigmaHH1 , cross section from Bragg scattering on tilted large waves, 1st Harmonic (upwind/downwind)\
11 : sigmaHV1 , cross section from Bragg scattering on tilted large waves, 1st Harmonic (upwind/downwind)\
12 : sigmaVH1 , cross section from Bragg scattering on tilted large waves, 1st Harmonic (upwind/downwind)\
13 : sigmaVV2 , cross section from Bragg scattering on tilted large waves, 2nd Harmonic (crosswind)\
14 : sigmaHH2 , cross section from Bragg scattering on tilted large waves, 2nd Harmonic (crosswind)\
15 : sigmaHV2 , cross section from Bragg scattering on tilted large waves, 2nd Harmonic (crosswind)\
16 : sigmaVH2 , cross section from Bragg scattering on tilted large waves, 2nd Harmonic (crosswind)\
17 : sigmaQS , Quasi-Specular (QS) cross section\
18 : sigma_vv , cross section for Bragg scattering only\
19 : sigma_hh , cross section for Bragg scattering only\
20 : sigma_hv , cross section for Bragg scattering only\
21 : sigma_vh , cross section for Bragg scattering only\
22 : 0.D0, placeholder\
23 : 0.D0, placeholder\
24 : 0.D0, placeholder\
25 : 0.D0, placeholder\
26 : 0.D0, placeholder\
27 : 0.D0, placeholder\
28 : k0/kd, ratio between electromagnetic and cutoff wavenumbers\
29 : A0, sigmaQS factor, to be multiplied to QS component\
30 : Re(epsi), real part of sea water dielectric constant\
31 : Im(epsi), imaginary part of sea water dielectric constant\
32 : sigma_vv0 , cross section for Bragg scattering only, omnidirectional component\
33 : sigma_hh0 , cross section for Bragg scattering only, omnidirectional component\
34 : Fr, foam fraction in percent (always 0 for the output file without foam) /!\ NOT USED /!\

The total cross section is compute from the sum of (1) the Bragg scattering on tilted large waves component and the (2) modified Quasi-Specular component. Assuming the angle phi between the propagation vector and the wind vector, the VV cross section is computed from:

sigmaVV = sigmaVV_B [07,09,13] + A0 [29] * sigmaQS [17]

where the directional tilted Bragg scattering is derived from the harmonic terms and phi as:

sigmaVV_B = sigmaVV0 [07] + sigmaVV1 [09] * cos ( phi ) + sigmaVV2 [13] * cos ( 2*phi )


# References

Anguelova, M. D., Dinnat, E., Kilic, L., Bettenhausen, M. H., English, S., Prigent, C., Meissner, T., Boutin, J., Newman, S., Johnson, B., Yueh, S., Kazumori, M., Weng, F., Stoffelen, A., & Accadia, C. (2022, June). Foam emissivity modelling with foam properties tuned by frequency and polarization. International Geoscience and Remote Sensing Symposium (IGARSS).

Anguelova, M. D., & Gaiser, P. W. (2013). Microwave emissivity of sea foam layers with vertically inhomogeneous dielectric properties. Remote Sensing of Environment, 139, 81–96. https://doi.org/10.1016/j.rse.2013.07.017

Cardone, V. J. (1969). Specification of the wind distribution in the marine boundary layer for wave forecasting. https://doi.org/10.21236/AD0702490

Charnock, H. (1955). Wind stress on a water surface. Q. J. R. Meteorol. Soc., 81, 639–640.

Cox, C., & Munk, W. (1954). Measurement of the Roughness of the Sea Surface from Photographs of the Sun’s Glitter. Journal of the Optical Society of America, 44(11), 838–850. https://doi.org/10.1364/JOSA.44.000838

Dinnat, E. P., Boutin, J., Caudal, G., & Etcheto, J. (2003). Issues concerning the sea emissivity modeling at L band for retrieving surface salinity. Radio Science, 38(4), n/a-n/a. https://doi.org/10.1029/2002RS002637

Donelan, M. A., Dobson, F. W., Smith, S. D., & Anderson, R. J. (1993). On the dependence of sea surface roughness on wave development. Journal of Physical Oceanography, 23(9), 2143–2149.

Durden, S. L., & Vesecky, J. F. (1985). A Physical Radar Cross-Section Model for a Wind-Driven Sea with Swell. IEEE Journal of Oceanic Engineering, 10(4), 445–451. https://doi.org/10.1109/JOE.1985.1145133

Elfouhaily, T., Chapron, B., Katsaros, K., & Vandemark, D. (1997). A unified directional spectrum for long and short wind-driven waves. Journal of Geophysical Research, 102(C7), 15781–15796.

Ellison, W., Balana, A., Delbos, G., Lamkaouchi, K., Eymard, L., Guillou, C., & Prigent, C. (1998). New permittivity measurements of seawater. Radio Science, 33(3), 639. https://doi.org/10.1029/97RS02223

Klein, L., & Swift, C. (1977). An improved model for the dielectric constant of sea water at microwave frequencies. IEEE Transactions on Antennas and Propagation, OE-2(1), 104–111. https://doi.org/10.1109/JOE.1977.1145319

Meissner, T., & Wentz, F. J. (2004). The complex dielectric constant of pure and sea water from microwave satellite observations. IEEE Transactions on Geoscience and Remote Sensing, 42(9), 1836–1849. https://doi.org/10.1109/TGRS.2004.831888

Meissner, T., & Wentz, F. J. (2012). The Emissivity of the Ocean Surface Between 6 and 90 GHz Over a Large Range of Wind Speeds and Earth Incidence Angles. IEEE Transactions on Geoscience and Remote Sensing, 50(8), 3004–3026. https://doi.org/10.1109/TGRS.2011.2179662

Meissner, T., Wentz, F. J., & Ricciardulli, L. (2014). The emission and scattering of L-band microwave radiation from rough ocean surfaces and wind speed measurements from the Aquarius sensor. Journal of Geophysical Research C: Oceans, 119(9), 6499–6522. https://doi.org/10.1002/2014JC009837


Monahan, E. C., & Lu, M. (1990). Acoustically relevant bubble assemblages and their dependence on meteorological parameters. IEEE Journal of Oceanic Engineering, 15(4), 340–349.

Monahan, E. C., & O’Muircheartaigh, I. G. (1986). Whitecaps and the passive remote sensing of the ocean surface. International Journal of Remote Sensing, 7(5), 627–642.

Stogryn, A. (1972). The emissivity of sea foam at microwave frequencies. Journal of Geophysical Research, 77(9), 1658–1666.


Yin, X., Boutin, J., Dinnat, E., Song, Q., & Martin, A. (2016). Roughness and foam signature on SMOS-MIRAS brightness temperatures: A semi-theoretical approach. Remote Sensing of Environment, 180, 221–233. https://doi.org/10.1016/j.rse.2016.02.005


Yueh, S. H. (1997). Modeling of wind direction signals in polarimetric sea surface brightness temperatures. IEEE Transactions on Geoscience and Remote Sensing, 35(6), 1400–1418. https://doi.org/10.1109/36.649793
