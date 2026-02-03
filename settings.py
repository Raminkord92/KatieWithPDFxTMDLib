
# Path to the AVHLIB directory
AVHLIBpath = '/root/physics/montecarlo/avhlib'

# Path to the directory where libLHAPDF.so is
LHAPDFpath = '/usr/local/lib'

# Fortran compiler with flags
FC = 'gfortran -fcheck=bounds'

########################################################################
# The following is necessary if you want to use TMDlib,
# and must be out-commented if you do not have TMDlib installed.
#
# Path to the directory where TMDlib-config is (for tmdlib-2.0.x or
# older you can still put the path where libTMDlib.so is)
TMDLIBpath = '/usr/local/lib'
#
# Path to the directory where libPDFxTMDLib.so is located.
# Leave empty if PDFxTMD is not available on this machine.
PDFXTMDpath = '/usr/local/lib'
#
# Do you want to use multiple TMD sets (available from version 2.0.2)?
#TMDLIBmulti = 'yes'
#
# Path to the directory where libgsl.so is
GSLpath =  '/usr/local/lib'
#
########################################################################

########################################################################
# No need to change the following if you want to keep the defaults:
#
# Default: 'dynamic'
#LINKMETHOD = 'static'
#
# The following only concerns dynamic linking.
# Default: 'gcc' for MacOS, and 'ld' for other platforms
#LINKER = 'gcc'

# If you want to combine KaTie with Mincas, put the path to
# the directory with the "src" and "include" of Mincas, and
# the directory with the source files of Merger here:
#MINCASpath = '/home/user0/physoft/MINCAS'
#MERGERpath = '/tmp/Merger'
#
########################################################################
