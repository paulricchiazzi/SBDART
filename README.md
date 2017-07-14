# SBDART

SBDART (Santa Barbara DISORT Atmospheric Radiative Transfer) is a FORTRAN computer code designed for the analysis of a wide variety of plane-parallel radiative transfer problems encountered in satellite remote sensing and atmospheric energy budget studies. The program is based on a collection of highly developed and reliable physical models, which have been developed by the atmospheric science community over the past few decades. More information on the model is available from SBDART_Introduction.pdf

Installation

SBDART is written entirely in FORTRAN. I have successfully compiled the code both on OSX and on Windows7 using the gfortran compiler using a switch to enforce FORTRAN 2003 compliance.  The makefile included in repository works on UNIX-like operating systems.  If you'ld like to use the makefile on Windows, you'll need to install something like Cygwin to make it work.
