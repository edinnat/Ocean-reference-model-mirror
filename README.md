# Ocean-model
Reference Quality Model For Ocean Surface Emissivity And Backscatter From The Microwave To The Infrared

# Compiling the source code

The source code comes with a 'Makefile' located in the code's root directory that provides all the information needed for compiling exectuables. Also provided is a very basic shell script that uses the Makefile to compile the radiometer code and produces a "Tb" executable. The shell script is also in the root directory and named 'compile_Tb.sh'. Users may need to change the compiler definition in the Makefile to fit their installation. By default, environmental variables are used to define various compilation parameters (see below), but users can enter directly what they want to use as compiler and paths into the Makefile. 

Default environment variables used for compilation:

"compil" should point to the user's compiler (e.g. gfortran, using full path if necessary)
"OceanMod" should point to the root folder of the ocean model code

Both these variables can be redefined in the Makefile without defining the environment variables.

Three executables can be created:
* Tb 	: The main radiometer code that computes the brightness temperature of the ocean surface. The computation relies on azimuthal harmonic decomposition for the angle between the surface emission vector and the wind direction.
* TbTot : Same as previous code except the simulation is ran for mutliple azimuth angles and the  . This code is much slower than 'Tb'. It is also likely to not be updated as often as 'Tb'.
* NRCS  : The code to compute the Normalized Radar Cross Section of the ocean surface.

Each exectuable above can be created from the root directory issuing one of the following  commands:
* make -f Makefile Tb
* make -f Makefile TbTot
* make -f Makefile NRCS

After compilation, the executable just created is placed in the 'Binaries' folder.

# Running the program

The program is run from one of the executables in the "Binaries" folder. The program will first look for a configuration (or parameters) file and some input data files at default paths. To define the paths properlyi for the user's installation, the user should define the following environment variables:


* inputFortran : Where to find the configuration files that contain parameters to run the program, such as SST, SSS, wind speed, frequency .... 

* outputFortran: Path where to save the results.

* OceanMod : Path to the the ocean emissivity/reflectivity model root (the path where this README.md file is most likely located).

Users can override some default setting by passing the full path to the configuration file as a parameter to the executable. For exemple, from the root folder of the code, users can run the command:

	Binaries/Tb Config/Tb.p

 which will use the config/parameters file in "Config" folder named "Tb.p". If no configuration file is provided, the program will use the default 'Tb.p' from the default path in "inputFortran". 

The program will save the results in the folder defined by "outputFortran" with a filename defined in the configuration file. It can output 2 files, with or without foam effects included. 

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

The Foam-free output file contains additional columns (26 - 43) with intermediate results for the brightness temperatures:

26-29 : Tv0_SPM, Th0_SPM, T30_SPM, T40_SPM = Small Perturbation Method Tb for v, h, Stokes 3 and Stokes 4 respectively, omnidirectional term.\
30-33 : Tv2_SPM, Th2_SPM, T32_SPM, T42_SPM = Small Perturbation Method Tb for v, h, Stokes 3 and Stokes 4 respectively, second harmonic in azimuth wrt wind direction.\
34-35 : Tv0_GO, Th0_GO = Geometric Optics Tb in v- and h-pol respectively, , omnidirectional term.\
36-39 : Tv1_GOTh1_GO, T31_GO, T41_GO  = Geometric Optics Tb in v- and h-pol respectively, first harmonic (upwind/downwind).\
40-43 : Tv2_GO, Th2_GO, T32_GO,T42_GO = Geometric Optics Tb in v- and h-pol respectively, second harmonic (crosswind).\

The total Tb for a direction phi between propagation vector and wind vector is computed from:

	Tb = TpN + Tp0 + Tp1 * cos (phi) + Tp2 * cos ( 2*phi )

for p = v or h pol, and from:

	p = p1 * sin (phi) + p2 * sin ( 2*phi )

for p = U or V (for the 3rd and 4th Stokes parameters respectively).
 
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
