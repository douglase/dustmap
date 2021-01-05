3 Jan 2020:
Added functionality to interface with a python wrapper [(dustmapPy)](https://github.com/douglase/dustmapPy) by Ewan Douglas (UArizona/MIT @douglase) and Yinzi Xin (Caltech/MIT @yinzi-xin).


dustmap v3.1.2

---

Last updated 4 Apr 2015

Written by Christopher Stark
STScI
cstark@stsci.edu

-------------------------

DESCRIPTION

dustmap is an IDL suite designed to synthesize images of debris disk models.  The disk models are assumed to be collections of 3D coordinates, each of which represents a dust particle.  Among other things, dustmap can produce density histograms, optical depth distributions, thermal emission images, and scattered light images.  Synthesized images can be viewed from an arbitrary viewing location, including within the dust disk.  The viewing direction is assumed to point toward the star.

Optical constants and scattering phase functions are calculated using Mie Theory.  dustmap includes a small library of dust compositions to use, or you can specifiy your own ASCII file containing the indices of refraction as a function of wavelength.  Alternatively, one can specify a generic power law for Qabs and Qsca and specifiy an albedo.  You can also override the scattering phase function with a Henyey-Greenstein scattering phase function.

-------------------------

INSTALLATION & USAGE

See included "Introduction_to_dustmap.pdf"

-------------------------

RELEASE NOTES

Please contact me if you find any bugs or have any comments/questions/suggestions.

What's new in version 3.1.2:

- Minor functional updates to dustmap_binary_writer.pro
- Improved the interpolation scheme for the Kurucz stellar models.
- Improved functionality of tsublimate keyword.  Previously setting tsublimate only removed the thermal emission from grains w/ t>tsublimate.  Now it also removes scattered light from those grains, too.
- The generic optical constants option is now functioning.  It was broken, unusable, and resulted in an error.
- Fixed a bug in the HG scattering phase function (only affected models with g > 0.9).
- Fixed a bug in the getKurucz.c file that resulted in buffer overflows on some systems.  If dustmap didn’t crash, this bug didn’t affect you.

What's new in version 3.1.1:

- Fixed a bug in the scaling factor.
- Fixed a bug in the HG scattering phase function.

What's new in version 3.1.0:

- Updated some of the syntax to make more intuitive (e.g. Rstar has been replaced with Lstar)
- Updated mie_single.c to return scattering phase function (scatt phase func is now calcualted via Mie Theory on the fly)
- Added COSTHETA and PFUNC keywords to return the scattering phase function calculated via Mie Theory.
- Added NCOSTHETA keyword which sets the resolution of the scattering phase function.  Default = 500.  Increase to better resolve scatt. phase func. at small and large scattering angles.  You can check the resolution of the scatt. phase func. by calculating delta_theta from the returned COSTHETA values.
- Added ability to load Kurucz stellar model (via KURUCZ keyword and necessary LOGG keyword)
- Added ability to ouput the stellar flux (via FSTAR keyword)
- Added the dustmap_binary_writer.pro code

What's new in version 3.0.3:

- Fixed a bug regarding the xshift keyword that prevented the disk from being illuminated by starlight properly

What's new in version 3.0.2:

- Fixed two minor memory leaks in dm2d.c
- Improved occultation-detection method in dm2d.c
- Fixed stellar extinction calculation in dm2d.c

What's new in version 3.0.1:

- Fixed bug in dmchecksyntax.pro that prevented DEGREES keyword from working properly (thanks to Alex Mustill for pointing this out)
- Fixed a comment in dustmap.pro that stated the output units were in erg cm^-2 per pixel when in fact they are in Jy per pixel

What's new in version 3.0.0:

- dustmap was completely rewritten in C for a significant improvement in run time
- input data format has changed:
      - identification of format more intuitive with DATATYPE keyword
      - input data format is now only binary; no more ASCII or .sav input files
      - no more need for VEL, BINARY, SAVFILE, or SWAPBIN keywords, which have been removed
- Multiple wavelengths now supported and coded for maximum efficiency
      - Output image can now be a 3D data cube of image vs. wavelength
- Histograms of the geometric optical depth can now be calculated via the OPTICALDEPTH keyword
- Optical constants now calculated from Mie Theory by default, requiring the COMPOSITION keyword
      - A small library of compositions have been included with the distribution
      - If Mie Theory is not desired, user can specify Qabs and Qsca power law (QEXP keyword) as well
         as the HG phase function g value (G keyword)
- NUMPIXELS keyword added and can be used in place of PIXELSIZE
- New ALLSKY keyword allows quick setting of an all-sky field of view
- New AITOFF keyword enables Aitoff projections (useful for all-sky images)
- Occultation of the star has be added and can be accessed using the OCCULT keyword
- TEMPSTAR keyword changed to TSTAR
- R_STAR keyword changed to RSTAR
- DUSTSIZE keyword changed to RDUST
- NODISP, NOPRINT, and REDUCEPRINT keywords allow user to reduce the amount of feedback
- COLLISIONRATE and IGNORECOLLISIONS keywords removed
- ITERATETDUST keyword removed; the dust grain temperature is iteratively calculated by default now
- IMAGE3D keyword removed; only 2D mode available now
