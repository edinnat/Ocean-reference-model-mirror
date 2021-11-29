# Makefile to compile Fortran Code for the Ocean Reference Model Project
#

Fort = 'gfortran'
#Fort = $(compil) # Fortran Compiler (from environmnnt variable 'compil')
#				 # e.g. 'g77', 'gfortran'
#
# SELECT Compiler options here
# Opts = -Waliasing -Wall -Wampersand -Warray-bounds -Wc-binding-type -Wcharacter-truncation -Wconversion -Wfunction-elimination -Wimplicit-interface -Wimplicit-procedure -Wintrinsic-shadow -Wintrinsics-std -Wline-truncation -Wno-align-commons -Wno-tabs -Wreal-q-constant -Wsurprising -Wunderflow -Wunused-parameter -Wrealloc-lhs -Wrealloc-lhs-all -Wtarget-lifetime
Opts = -fPIC

# Do not write command lines to screen before execution
# 	Comment out to debug
.SILENT : 

#-----------------------------------------------------------------------
# Check existence of environment variable 'OceanMod'  -------------------
# that contains root path to code					 -------------------

ifndef OceanMod
.PHONY : ERROR
ERROR :
	
	echo "	ERROR : environment variable OceanMod is not defined"
	echo "  set it as the root path to FORTRAN code."
	exit

endif

#-----------------------------------------------------------------------
# Definition of sub-folder containing various componenet of model
#

FN = $(OceanMod) # Directory containing all Fortran source code
                # Environmnet variable OceanMod needs to be assigned

TR = $(join $(FN), /Atmosphere/Radiative_Transfer/)
VT = $(join $(FN), /Atmosphere/Wind/)
EP = $(join $(FN), /Ocean/Dielectric_Constant/)
EC = $(join $(FN), /Ocean/Foam/)
SP = $(join $(FN), /Ocean/Sea_Spectrum/)
HO = $(join $(FN), /Ocean/Swell/)
DI = $(join $(FN), /Scattering/)
RO = $(join $(FN), /Radiometer/)
RA = $(join $(FN), /Radar/)
NU = $(join $(FN), /Numerical/)
GE = $(join $(FN), /Geometry/)
IO = $(join $(FN), /IO/)
#-----------------------------------------------------------------------

# Group necessary files by topic

# 	Atmospheric Radiative Transfer
FilesTR = TbAtmo.f paramAtmo.f tauP.f convDpath.f

# 	Wind
FilesVT = wind.f funcUstar_y.f funcUstar_c.f 

#   Dielectric Constant for sea water
FilesEP = epsilon_KS.f epsilon_El.f epsilon_MW.f

#   Foam emissivity and foam fraction
FilesEC = foam.f esf.f WISE2001.f MonahanLu.f foam_fr_Yin16.f foam_emiss_Yin16.f
FilesEC += mod_foam_emiss.f90 foam_emiss.f90

# 	Sea Spectrum
FilesSP = c.f Su.f Sc.f P.f Sigma.f spectrum_DV_Base.f spectrum_DV.f spectrum_Lemaire.f 
FilesSP += spectrum_El.f spectrum_powerLaw.f fksup.f fkinf.f spectrum_Yueh.f SuScKud.f
FilesSP += spectrum_Kudryavtsev.f

#  Swell
FilesHO = infSwell.f TabSwell.f intSwellX.f intSwellY.f

# Scattering coefficients
FilesDI = Gvv2.f Gvv1.f Ghh2.f Ghh1.f Gvh1.f Ghv1.f Ghv2.f

# Numerical Routines
FilesNU = qromb2D1.f qromb2D2.f dichoto.f trapzd2D1.f polint.f qromb.f 
FilesNU += trapzd.f trapzd2D2.f zroots.f laguer.f mylib.f 
FilesNU += root.f extr_val.f
FilesNU += spline.f splint.f

# Geometric Transformations
FilesGE = Glob2Loc.f Loc2Glob.f

# 	Gestion Entrees / Sorties
FilesIO = readvec.f

# ---- Radiometer
#		For TbTot (to compute every azimuthal angle, without use of harmonics, TBC)
FilesRO1 = TbTot.f func2D1.f func2D21.f func2D22.f #func2D21_2.f func2D22_2.f

#		For Tb (compute azimuthal signal as Fourier series component)
FilesRO2 = Tb.f func2D1.f func2D21.f func2D22.f #func2D21_2.f func2D22_2.f

# ---- Radar
FilesRA = NRCS.f func2D1radar.f func2D21radar.f func2D22radar.f


# Concatenate paths and filenames
# 	1) Select file sources for radiometer program
#
SourcesRO = $(addprefix $(TR),$(FilesTR))
SourcesRO += $(addprefix $(VT),$(FilesVT))
SourcesRO += $(addprefix $(EP),$(FilesEP))
SourcesRO += $(addprefix $(EC),$(FilesEC))
SourcesRO += $(addprefix $(SP),$(FilesSP))
SourcesRO += $(addprefix $(HO),$(FilesHO))
SourcesRO += $(addprefix $(DI),$(FilesDI))
SourcesRO += $(addprefix $(NU),$(FilesNU))
SourcesRO += $(addprefix $(GE),$(FilesGE))
SourcesRO += $(addprefix $(IO),$(FilesIO))

# 	1.1) for Tb Tot
SourcesRO1 = $(SourcesRO) $(addprefix $(RO),$(FilesRO1))

# 	1.2) for Tb
SourcesRO2 = $(SourcesRO) $(addprefix $(RO),$(FilesRO2))

# 	2) for NRCS
SourcesRA = $(SourcesRO) $(addprefix $(RA),$(FilesRA))


#-----------------------------------------------------------------------

# Create objects *.o from source files *.f
ObjectsTbTot = $(SourcesRO1:.f=.o)
BasenamesTb = $(basename $(SourcesRO2))
ObjectsTb    = $(addsuffix .o, $(BasenamesTb))
#ObjectsTb    = $(SourcesRO2:.f=.o)
ObjectsNRCS   = $(SourcesRA:.f=.o)

#-----------------------------------------------------------------------

# Set up cleaning function : remove all object files

clean :
	echo Removing `find ${OceanMod}/* | grep '\.o'`
	rm -f `find ${OceanMod}/* | grep '\.o'`

#-----------------------------------------------------------------------

# Build Binaries		    --------------------------------------------

Tb : $(ObjectsTb)
	echo ======= linking Tb ============
	$(Fort) $(Opts) $^ -o $(join $(join $(FN), /Binaries/), $@)

TbTot : $(ObjectsTbTot)
	echo ======= linking TbTot ============
	$(Fort) $^ -o $(join $(join $(FN), /Binaries/), $@)

NRCS : $(ObjectsNRCS)
	echo ======= linking NRCS ============
	$(Fort) $^ -o $(join $(join $(FN), /Binaries/), $@) 

#-----------------------------------------------------------------------

# Rule to build all objects --------------------------------------------

%.o : %.f
	echo Compiling $*.f
	$(Fort) -c $*.f -o $*.o -O3 -funroll-loops
	echo '  -->' $*.o created

%.o : %.f90
	echo Compiling $*.f90
	$(Fort) -c $*.f90 -fPIC -o $*.o -O3 -funroll-loops
	echo '  -->' $*.o created
#-----------------------------------------------------------------------
