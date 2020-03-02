# Simulation of foreground corresponding to EoR

The HI 21 cm signal is one of the most comprehensive tools to probe the early universe. The signal
that we receive on earth consists of the 21 cm signal along with noise and signals from various
astrophysical sources known as the foregrounds. The main component of foregrounds are Galactic
diffused sychrotron emission ( GDSE – 70 %) and the extragalactic point sources (27 %).
This project is aimed at simulating the different foreground components for a small patch of sky
using the flat sky approximation (Datta et al., 2007).
The simulated images are generated using C. For this purpose the simulated images are generated in
FITS format so that they can be analysed and different properties can be obtained from them.

## Foreground Component

**Galactic Diffused Synchrotron Emission (GDSE)**
There are 5 files for this which are:
1. *grf.c* – This is the main programme and almost all the parameters are defined here. It can
simulate the GDSE sky map for a variety of frequencies and the image can be analysed using
SAODS9.
2. *beam.c* – This file contains the expression for primary beam pattern of the radio telescopes which
has the form (J 1 (u)/u) 2 but can often be represented by a gaussian function.
3. *input.grf* - This file contains some of the parameters such as the number of grid points, number
of image slices corresponding to different frequencies etc.
4. *makefile* - Instead of typing lengthy commands in terminal to compile the programme , one can
simply use the makefile to compile and run it.
5. *fitsprog.c*
To compile and run the programme, one has to type the following

*make grf*
*./grf input.grf output.fits *

If you already have a file named ‘output.fits’, then the programme will ask you whether you want
to overwrite it. If you want to overwrite it, type ‘y’.

**Extragalactic point sources**
For simulating the extragalactic point sources at 150 MHz, an input differential source count has
been used (Ghosh et al., 2012).
There are 2 files, one generating the image in FITS format while the other is used for generating
only the dataset.
