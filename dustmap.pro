pro dustmap, INPUTFILE, HISTO, $
             DISTANCE=distance, FOV=fov, PIXELSIZE=pixelsize, NUMPIXELS=numpixels, DATATYPE=datatype, $
             INCLINATION=inclination, LONGITUDE=longitude, PA=pa, $
             AU=au, DEGREES=degrees, $
             OPTICALDEPTH=opticaldepth, SCATTERED=scattered, THERMAL=thermal, LAMBDA=lambda, $
             TSTAR=Tstar, LSTAR=lstar, KURUCZ=kurucz, LOGG=logg, FSTAR=fstar, $
             SCALING=scaling, RDUST=rdust, COMPOSITION=composition, LNKFILE=LNKFILE, TSUBLIMATE=tsublimate, $
             QEXP=qexp, ALBEDO=albedo, OC_GENERIC_FLAG=oc_generic_flag, G=g, HG_flag=HG_flag, EXT=ext, $
             PFUNC=pfunc, COSTHETA=costheta, NCOSTHETA=NCOSTHETA, QABS=qabs, QSCA=qsca, $
             AITOFF=aitoff, ALLSKY=allsky, VGA=vga, HD=hd, $
             XSHIFT=xshift, AZAVG=azavg, DISTMASK=distmask, IWA=iwa, EFFRDUST=EFFRDUST, $
             MARKLOCATIONS=marklocations, MARKWEIGHT=markweight, MARKABS=markabs, $
             LOG=log, COLORBAR=colorbar, CBTITLE=cbtitle, $
             NODISP=nodisp, NOPRINT=noprint, REDUCEPRINT=reduceprint

; dustmap.pro
; Written by Christopher Stark
; NASA GSFC
; Last updated 10 June 2014

; Description: Uses discrete 3D locations (in AU) to create an image
; from an arbitrary viewing location of dust around a star.  Can
; create density histograms, optical depth images, scattered light
; images, and thermal emission images. The viewing direction is
; assumed to point toward the star.


; SYNTAX
; See dmchecksyntax.pro or just call dustmap without any parameters
;
; --- ARGUMENTS ---
; INPUTFILE: file path + name (e.g. '/user/data/file.dat') of postion data in
;            this format: particleID x y z (in AU)
;            NOTE: inputfile can be an array of input files
; HISTO: the name of the histogram in which dustmap returns the
;        results, with units of Jy per pixel 

; --- MANDATORY KEYWORDS ---
; DISTANCE: distance to the system (pc, unless /AU keyword is set)
; FOV: a 1 or 2-element array defining the field of view in mas
;      (FOV is mandatory unless /ALLSKY keyword is set)
; PIXELSIZE or NUMPIXELS: the size of each square pixel in mas, or
;                         the number of pixels
; DATATYPE: one of four data types that defines the format of the
;     binary input file (see dmgetdata.c for more info)
;     Note: you can use dustmap_binary_writer.pro to create a binary data file
;           1 = particle ID and 3D position
;           2 = particle ID, 3D position, 3D velocity
;           3 = particle ID, 3D position, intensity
;           4 = particle ID, 3D position, 3D velocity, intensity 

; --- OPTIONAL KEYWORDS ---
; INCLINATION: inclination of the system (degrees)
; LONGITUDE: rotates the system about its z-axis, i.e. controls the
;            longitude of the planets (degrees)
; PA: position angle of the system on the sky (degrees)
; AU: when set, the distance is treated as AU instead of pc
; DEGREES: when set, then fov and pixelsize are treated as degrees instead of mas
; OPTICALDEPTH: the output histogram will be weighted by the dust
;               grain cross-section
; THERMAL: the output histogram includes the flux from thermal emission
; SCATTERED: the output histogram includes the flux from scattered light
;      Must also specify...
;      LAMBDA: array of wavelength(s) (microns)
;      TSTAR: effective stellar temperature (degrees Kelvin)
;      LSTAR: stellar luminosity (Solar luminosities)
;      RDUST: radius of dust grains (microns). This can be a single value
;             (assumed to be same for all input files) or an array with
;             length equal to the # of input files (for a distribution of grain sizes).
;      COMPOSITION: the composition of the dust. Can be one of several strings:
;              'astrosil', 'aboliv', 'organics', 'tholin', troilite', or 'water'
;      ...or...
;      LNKFILE: if you have your own lnk file to use, you may bypass one of the above
;               compositions by using this keyword that contains the file's full path
;      Optional...
;      KURUCZ: set this flag to load a Kurucz stellar atmosphere model
;      LOGG: must specificy when requesting a Kurucz model
;      FSTAR: optional output stellar flux vs lambda as seen by the observer (Jy)
; SCALING: this scales all intensity values by the assigned amount.
;          If intensity values are set, this is a multiplicative factor.
;          If intensity values are not set, intensity = scaling.
;          This can be a single value (assumed to be same for all input
;          files) or an array w/ length equal to the # of input files.
; OC_GENERIC_FLAG: set this to a non-zero value to use generic optical constants
;      QEXP: Qabs and Qsca are set equal to 1 for lambda<2*PI*s, and (2*PI*s/lambda)^QEXP
;            otherwise; this should typically be a positive number if non-zero
;      ALBEDO: The Qsca values are all multiplied by albedo / (1 - albedo)
; HG_FLAG: set this to a non-zero value to use a Henyey-Greenstein
;          scattering phase function
;      G: g value of HG SPF instead of the SPF calculated from Mie Theory
; EXT: add negative flux to account for the dust that blocks the star
; PFUNC: the scattering phase function used vs wavelength (only for the last grain size)
; COSTHETA: the cos(scatt ang) values used to calculate PFUNC
; NCOSTHETA: # of scattering angles to calculate (default is 1001),
;            can be increased to improve resolution at small/large
;            scattering angles
; QABS: the returned absorption efficiency vs wavelength (only for the last grain size)
; QSCA: the returned scattering efficiency vs wavelength (only for the last grain size)
; ALLSKY: when set, the fov input is ignored (& not needed) and an
;         all-sky map is created
; AITOFF: the output histogram/image is an Aitoff projection (useful
;         for all-sky maps)
; VGA, HD: automatically sets the number of pixels and fov aspect
;          ratio for VGA or HD resolutions
; XSHIFT: shifts the observer off of the line of sight to star in the
;         x direction by this amount (in AU)
; AZAVG: azimuthally average the data points AZAVG number of times,
;         good for reducing Poisson noise in symmetric disks
; DISTMASK: defines the radius of a sphere centered on the observer,
;          within this sphere dust is ignored/deleted
; IWA: the angular radius of an occulting spot to place in the
;         image (in mas, unless degrees keyword is set)
; EFFRDUST: the (un)physical size of a grain in AU, used for Gaussian
;           smoothing. Useful for scattered light and thermal images
;           when some points are very bright because they are close
;           to the camera.  Really slows down code.
; MARKLOCATIONS:  This is an array of at least 3 parameters
;      (e.g. [x,y,z], or [[x1,y1,z1],[x2,y2,z2]]) of 3D positions to mark
;      in the 2D image.  Each supplied location is marked in the image
;      with the maximum brightness of the image.
;      Can also specify....
;      MARKWEIGHT: the "brightness" of the mark.  By default, it
;                  should be a fraction of the maximum brightness in histo
;                  and set to a value between 0 and 1.  If weightabs
;                  keyword is set, then it likely must be much larger.
;      MARKABS: if this keyword is set, the markweight value
;               is used literally, and is not treated as fraction of the
;               maximum value in histo
; LOG: if this keyword is set, the image displayed is in a
;      logarithmic scale.  The output histo is unaffected.
; COLORBAR: if this keyword is set, a color bar is displalyed along
;           with the image--default is no colorbar
;      CBTITLE: the title of the color bar
; NODISP: supress displaying histo when created
; NOPRINT: surpress all output text
; REDUCEPRINT: reduce the amount of output text



;Determine directory path that dustmap is currently located in
;All supporting files distributed with dustmap must be in the same
;directory
dustmaproutinename='DUSTMAP'
dustmappath=routine_info(dustmaproutinename,/source)
dustmappath=strmid(dustmappath.path,0,strlen(dustmappath.path)-strlen(dustmaproutinename)-4)




;%%%%%%%%%%
;CHECK SYNTAX
;%%%%%%%%%%

dmchecksyntax, inputfile, histo, $
               distance, inclination, longitude, pa, $
               fov, pixelsize, numpixels, $
               datatype, scaling, userscaling, $
               au, degrees, $
               opticaldepth, scattered, thermal, lambda, $
               tstar, lstar, logg, kurucz, $
               rdust, composition, Tsublimate, ncostheta, $
               iwa, log, colorbar, cbtitle, $
               aitoff, aitoff_flag, allsky, $
               nodisp, noprint, reduceprint, $
               lnkfile, densityhisto_flag, opticaldepth_flag, $
               scatteredlight_flag, thermalemission_flag, verbose_flag, $
               qexp, albedo, g, oc_generic_flag, HG_flag, ext, xshift, effrdust, $
               histoxsize, histoysize, vga, hd, $
               marklocations, markx, marky, markz, markweight, nmarks, markabs, markabs_flag, azavg, $
               dustmappath

if keyword_set(noprint) then reduceprint=1

userfov = fov
userpixelsize = pixelsize
userdistance = distance
userinclination = inclination
userlongitude = longitude
userpa = pa

kurucz_file = dustmappath+'fp00k2odfnew.pck'

;%%%%%%%%%%
;OPTIONAL STUFF
;%%%%%%%%%%


; Convert fov and pixelsize from degrees to mas if necessary
if keyword_set(degrees) then begin
   userfov *= 3600000.0
   userpixelsize *= 3600000.0
   iwa *= 3600000.0
endif


; Convert distance from AU to pc if necessary
if keyword_set(au) then userdistance /= 206265.0 ;convert from AU to parsec
if keyword_set(au) and keyword_set(distmask) then distmask /= 206265.0 ;convert from AU to parsec







; %%%%%%%%%%%%%%%%
; MAIN ROUTINE
; %%%%%%%%%%%%%%%%

numfiles = size(inputfile,/n_elements)
if (not keyword_set(noprint)) then print, strcompress(string(numfiles),/remove_all),' file(s) to process...'


;Make sure variables are defined to be compatible with dustmap_c.c
nlambda = n_elements(lambda)
histo = dblarr(histoxsize, histoysize, nlambda)
histoxsize = long(histoxsize)
histoysize = long(histoysize)
iflen=strlen(inputfile)
ifn=n_elements(inputfile)
ifn=ifn[0]
inputfile4c=bytarr(total(iflen) + ifn)
strloc=0
for i=0,ifn-1 do begin
   inputfile4c[strloc:strloc+iflen[i]-1] = byte(inputfile[i])
   inputfile4c[strloc+iflen[i]] = byte(0)
   strloc += (iflen[i] + 1)
endfor
numfiles = long(n_elements(inputfile))
datatype = long(datatype)
userdistance = float(userdistance)
fovx = float(userfov[0])
fovy = float(userfov[1])
userpixelsize = float(userpixelsize)
userinclination = float(userinclination)
userlongitude = float(userlongitude)
userpa = float(userpa)
lambda = float(lambda)
nlambda = long(nlambda)
lstar = float(lstar)
Tstar = float(Tstar)
logg = float(logg)
kurucz = long(kurucz)
fstar = fltarr(nlambda)
rdust = float(rdust)
Tsublimate = float(Tsublimate)
lnkfile4c = bytarr(strlen(lnkfile)+1)
lnkfile4c[0:strlen(lnkfile)-1] = byte(lnkfile)
lnkfile4c[strlen(lnkfile)] = byte(0)
kuruczfile4c = bytarr(strlen(kurucz_file)+1)
kuruczfile4c[0:strlen(kurucz_file)-1] = byte(kurucz_file)
kuruczfile4c[strlen(kurucz_file)] = byte(0)
userscaling = double(userscaling)
iwa = float(iwa)
aitoff_flag = long(aitoff_flag)
densityhisto_flag = long(densityhisto_flag)
opticaldepth_flag = long(opticaldepth_flag)
scatteredlight_flag = long(scatteredlight_flag)
thermalemission_flag = long(thermalemission_flag)
verbose_flag = long(verbose_flag)
generic_oc_exp = float(qexp)
generic_oc_albedo = float(albedo)
HG_g = float(g)
oc_generic_flag = long(oc_generic_flag)
HG_flag = long(HG_flag)
extinction_flag = long(ext)
xshift = float(xshift)
effrdust = float(effrdust)
markx = float(markx)
marky = float(marky)
markz = float(markz)
markweight = float(markweight)
nmarks = long(nmarks)
markabs_flag = long(markabs_flag)
azavg = long(azavg)
if not keyword_set(distmask) then distmask = 0.0
distmask = float(distmask)
ncostheta = long(ncostheta)
internal_pfunc = fltarr(nlambda,ncostheta)
costheta = fltarr(ncostheta)
nrdust = n_elements(uniq(rdust,sort(rdust)))
pfunc = fltarr(nrdust,nlambda,ncostheta)
qabs = fltarr(nrdust,nlambda)
qsca = fltarr(nrdust,nlambda)


;Now call dustmap_c.c and have it do the dirty work      
x = CALL_EXTERNAL(dustmappath+'dustmap_c.so', 'dustmap_c', $
                  histo, histoxsize, histoysize, inputfile4c, numfiles, datatype, $
                  userdistance, fovx, fovy, userpixelsize, userinclination, userlongitude, userpa, $
                  lambda, nlambda, Lstar, Tstar, logg, kurucz, kuruczfile4c, fstar, $
                  rdust, Tsublimate, lnkfile4c, userscaling, $
                  iwa, aitoff_flag, densityhisto_flag, opticaldepth_flag, $
                  scatteredlight_flag, thermalemission_flag, verbose_flag, $
                  oc_generic_flag, generic_oc_exp, generic_oc_albedo, HG_flag, HG_g, $
                  extinction_flag, xshift, effrdust, $
                  markx, marky, markz, markweight, nmarks, markabs_flag, azavg, distmask, $
                  ncostheta, internal_pfunc, costheta, pfunc, qabs, qsca, nrdust)

displayhisto = histo[*,*,0]
if keyword_set(log) then begin
   displayhisto += min(displayhisto[where(displayhisto gt 0)])
   displayhisto = alog10(displayhisto)
endif


if not keyword_set(nodisp) then begin
   ; Open window 0 for display (default)
   if keyword_set(colorbar) then window, 0, xsize=histoxsize, ysize=1.125*histoysize $
   else if (!d.x_vsize ne histoxsize or !d.y_vsize ne histoysize or !d.window ne 0) then window, 0, xsize=histoxsize, ysize=histoysize
   tvscl,displayhisto
endif



if keyword_set(colorbar) then begin
   if keyword_set(log) then begin
      dmcolorbar, cbtitle=cbtitle, xrange=[10.^min(displayhisto),10.^max(displayhisto)]
   endif else begin
      dmcolorbar, cbtitle=cbtitle, xrange=[min(displayhisto),max(displayhisto)]
   endelse
endif



if (not keyword_set(noprint)) then begin
   print,' '
   if not keyword_set(image3D) then begin
      if keyword_set(thermal) then print,'Thermal emission flux histogram created.'
      if keyword_set(scattered) then print,'Scattered light flux histogram created.'
      if not keyword_set(thermal) and not keyword_set(scattered) then print,'2D density distribution created.'
   endif
endif


histo=reform(histo) ;if nlambda=1, this removes the extra dimension so it's just a simple 2D array, not a datacube

end
