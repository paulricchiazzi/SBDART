# SBDART

SBDART (Santa Barbara DISORT Atmospheric Radiative Transfer) is a FORTRAN computer code designed for the analysis of a wide variety of plane-parallel radiative transfer problems encountered in satellite remote sensing and atmospheric energy budget studies. The program is based on a collection of highly developed and reliable physical models, which have been developed by the atmospheric science community over the past few decades. More information on the model is available from SBDART_Introduction.pdf

SBDART Installation

SBDART is written entirely in FORTRAN.  Some parts of the code were developed using old versions of FORTRAN dating back to the 1970s.   For this release, the most archaic parts of the code have been updated, and now the program conforms to the 2003 FORTRAN standard.  The gfortran compiler successfully compiles the code both on OSX and on Windows7.  The makefile included in the repository works on UNIX-like operating systems.  If you would like to use the makefile on Windows, you'll need to install something like Cygwin or MinGW to make it work.

RunRT

SBDART is a command line code.  It reads inputs from an input file (INPUT) and writes text output to stdout.  This mode of operation is
useful for generating output for graphics postprocessors or for generating databases for remote sensing retrievals, but it is not very
useful for students trying to understand radiative transfer. The RunRT graphical user interface is designed to fill this gap.
RunRT is a python (2.7) graphics code that makes it easy to run SBDART and display its output graphically in ways that bring out the essential features of the radiative transfer physics. More information on how to setup and run RunRT may be found in file runrtdoc.txt in the RunRT folder.  There should be only one editing step necessary to install set up RunRT to run with SBDART, and that is to change the value of the variable sbdartexe to the path to the SBDART executable on your system.  (see line 25 of RunRT.py).

  

