Folder 'toGitHub' contains only the directories and the files that have been changed to add tuned foam emissivity.
These are:
	* 2 new F90 files for foam emissivity in dir Ocean/Foam/
	* Main progrma in modified file Tb.f in dir Radiometer
	* Input variables in modified file Tb.p in dir Config
	* Modified Makefile to compile the new F90 file. 

The changes in the main program (Tb.f) are:
	* Frequency band (nuBand = L,C,X,K,Ka,W) as input from the config file Tb.p (Line 461 in Tb.f)
	* Introduce case 7 for foam emissivity tuned for foam parameters (line 1757 in Tb.f)
	* Call to subroutine for tuned foam emissivity (line 1759 in Tb.f)
	* Stability parameter for foam coverage M1 (Monahan, 1986, eq.5) with minus sign "-" (Line 1278 in Tb.f)
		as atm. stability from the config file Tb.p is Tair-Twater
	* Added foam + atmosphere calculations (only foam without atmosphere in the original) (Lines 2082-2104 in Tb.f)

The foam emissivity uses the fast calculations with the semi-closed form of Yin et al (2016), but foam parameters are tuned by frequency and polarization. 


Maggie Anguelova
