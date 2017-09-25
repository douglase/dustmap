pro dmchecksyntax, inputfile, histo, $
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
                   qexp, albedo, g, generic_oc_flag, HG_flag, ext, xshift, effrdust, $
                   histoxsize, histoysize, vga, hd, $
                   marklocations, markx, marky, markz, markweight, nmarks, markabs, markabs_flag, azavg, $
                   dustmappath


  ;dmchecksyntax.pro
  ;Checks to make sure the user input to dustmap.pro is correct.
  ;This isn't idiot-proof, but it's fairly robust.

  if (n_elements(inputfile) lt 1) then begin
     print, " "
     print, "Generalized syntax for dustmap.pro:"
     print, " "
     print, "dustmap, inputfile_list, image_output $" 
     print, "         , DISTANCE=distance in pc, FOV=field of view in mas"
     print, "         , PIXELSIZE=pixel size in mas, NUMPIXELS=number of pixels, DATATYPE=input data type code $"
     print, "         [, INCLINATION=inclination in degrees, LONGITUDE=longitude in degrees, PA=position angle in degrees] $"
     print, "         [, /AU, /DEGREES] $"
     print, "         [, /OPTICALDEPTH, /SCATTERED, /THERMAL] $"
     print, "         [, LAMBDA=wavelength(s), SCALING=scaling factor] $"
     print, "         [, TSTAR=stellar temperature in K, LSTAR=stellar luminosity in LSun, KURUCZ=Kurucz model flag, LOGG=surface gravity] $"
     print, "         [, RDUST=dust radius in microns, COMPOSITION=dust composition string, LNKFILE=composition file string] $"
     print, "         [, /OC_GENERIC_FLAG, QEXP=Qsca & Qabs exponent, ALBEDO=albedo, /HG_flag, G=HG scattering asymmetry parameter] $"
     print, "         [, PFUNC=phase function array, COSTHETA=cos(scattering angle), NCOSTHETA=# of dcostheta steps to resolve SPF] $"
     print, "         [, /AITOFF, /ALLSKY, /VGA, /HD] $"
     print, "         [, XSHIFT=horizontal shift of observer in AU, AZAVG=integer for azimuthal averaging] $"
     print, "         [, DISTMASK=radius from observer in AU within which dust is ignored, IWA=inner working angle] $"
     print, "         [, EFFRDUST=effective radius of dust in AU] $"
     print, "         [, MARKLOCATIONS=3D positions to mark, MARKWEIGHT=mark weight, /MARKABS] $"
     print, "         [, /LOG, /COLORBAR, CBTITLE=colorbar title] $"
     print, "         [, /NODISP, /NOPRINT, /REDUCEPRINT]"
     print, " "
     print, "See dustmap.pro file for more info."
     stop,' '
  endif

  ;Check for type
  if (size(inputfile,/type) ne 7) then stop, 'ERROR: inputfile must be an array of strings.'

  ;Check for necessary keywords
  if n_elements(inputfile) lt 1 then stop, 'ERROR: inputfile must be specified'
  if n_elements(distance) eq 0 then stop, 'ERROR: distance must be specified'
  if n_elements(fov) eq 0 and not keyword_set(allsky) then stop, 'ERROR: fov must be specified'
  if n_elements(pixelsize) eq 0 and n_elements(numpixels) eq 0 and not keyword_set(hd) and not keyword_set(vga) then stop, 'ERROR: pixelsize or numpixels must be specified'

  ;Check for confusion of keywords
  if n_elements(pixelsize) ne 0 and n_elements(numpixels) ne 0 then stop, 'ERROR: cannot set both pixelsize and numpixels'
  if n_elements(pixelsize) ne 0 and (keyword_set(vga) or keyword_set(hd)) then stop, 'ERROR: cannot set both pixelsize and HD or VGA keywords'
  if n_elements(numpixels) ne 0 and (keyword_set(vga) or keyword_set(hd)) then stop, 'ERROR: cannot set both pixelsize and HD or VGA keywords'
  if n_elements(datatype) eq 0 then begin
     print,'ERROR: datatype must be specified'
     print,'Possible datatypes:'
     print,' 1 = particle ID and 3D position'
     print,' 2 = particle ID, 3D position, 3D velocity'
     print,' 3 = particle ID, 3D position, intensity'
     print,' 4 = particle ID, 3D position, 3D velocity, intensity'
     print,'You can use the included dustmap_binary_writer to create the binary data file.'
     stop
  endif
  if keyword_set(generic_oc_flag) then begin
     if n_elements(qexp) eq 0 then stop,'ERROR: When using generic optical constants, you must specificy qexp.'
     if n_elements(albedo) eq 0 then stop,'ERROR: When using generic optical constants, you must specificy albedo.'
     if not keyword_set(HG_flag) then stop,'ERROR: When using generic optical constants, you must use a Henyey-Greenstein scattering phase function.'
  endif
  if (n_elements(qexp) ne 0 or n_elements(albedo) ne 0) and not keyword_set(generic_oc_flag) then stop,'ERROR: if you want to use generic optical constants, you must set oc_generic_flag=1'
  if keyword_set(HG_flag) then begin
     if n_elements(g) eq 0 then stop,'ERROR: When using generic optical constants, you must specificy qexp.'
  endif
  if n_elements(g) ne 0 and not keyword_set(HG_flag) then stop,'ERROR: if you want to use a Henyey-Greenstein scattering phase function, you must set HG_flag=1'

  ;Check for appropriate number of elements
  if n_elements(distance) gt 1 then stop,'ERROR: distance must be a scalar value.'
  if n_elements(inclination) gt 1 then stop,'ERROR: inclination must be a scalar value.'
  if n_elements(longitude) gt 1 then stop,'ERROR: longitude must be a scalar value.'
  if n_elements(pa) gt 1 then stop, 'ERROR: pa must be a scalar value.'
  if n_elements(tstar) gt 1 then stop,'ERROR: Tstar must be a scalar value.'
  if n_elements(lstar) gt 1 then stop, 'ERROR: Lstar must be a scalar value.'
  if n_elements(kurucz) gt 1 then stop, 'ERROR: kurucz must be a scalar value.'
  if n_elements(logg) gt 1 then stop, 'ERROR: logg must be a scalar value.'
  if n_elements(pixelsize) gt 1 then stop,'ERROR: pixelsize must be a scalar value.'
  if n_elements(numpixels) gt 1 then stop,'ERROR: numpixels must be a scalar value.'
  if n_elements(fov) gt 2 then stop,'ERROR: fov must be a scalar or 2-element vector.'
  if n_elements(datatype) gt 1 then stop,'ERROR: datatype must be a scalar value.'
  if keyword_set(scaling) then begin
     if n_elements(scaling) ne 1 and n_elements(scaling) ne n_elements(inputfile) and n_elements(scaling) ne 0 then stop, 'ERROR: scaling must be a single value or an array w/ length equal to the number of input files'
  endif
  if n_elements(marklocations) ne 0 then begin
     markdim = size(marklocations,/dim)
     if n_elements(markdim) gt 2 then stop, 'ERROR: marklocations must be an n by 3 array.'
     if n_elements(markdim) eq 1 then begin
        markdim = markdim[0]
        nmarks = 1
     endif else begin
        nmarks = markdim[1]
        markdim = markdim[0]
     endelse
     if nmarks eq 0 then stop, 'ERROR: must supply marklocations, an n by 3 array.'
     if markdim ne 3 then stop, 'ERROR: marklocations must be an n by 3 array.'
     if n_elements(markweight) eq 0 then markweight = fltarr(nmarks) + 1.0 ;weights are 1 by default
     if n_elements(markweight) ne nmarks then stop, 'ERROR: markweight must have same number of entries as marklocations.'
  endif
  if keyword_set(markabs) and (n_elements(marklocations) eq 0 or n_elements(markweight) eq 0) then stop, 'ERROR: /MARKABS specifies absolute mark weighting.  You must also supply markweight and marklocations.'

  ;Check for keyword-dependent options
  if keyword_set(scattered) or keyword_set(thermal) then begin
     if n_elements(lambda) eq 0 then stop,'ERROR: for a scattered light or thermal image, you must supply lambda (microns).'
     if n_elements(Lstar) eq 0 then stop,'ERROR: for a scattered light or thermal image, you must supply Lstar (L_solar).'
     if n_elements(Tstar) eq 0 then stop,'ERROR: for a scattered light or thermal image, you must supply Tstar (degrees K).'
     if keyword_set(kurucz) and n_elements(logg) eq 0 then stop,'ERROR: for a scattered light or thermal image with a Kurucz model, you must supply logg.'
     if n_elements(rdust) eq 0 then stop,'ERROR: for a scattered light or thermal image, you must supply rdust (microns).'
     if not keyword_set(composition) and not keyword_set(lnkfile) and not keyword_set(generic_oc_flag) then begin
        print,'ERROR: for a scattered light or thermal image, you must supply composition, lnkfile, or set generic_oc_flag=1.'
        stop,'Possible composition values: astrosil, iron, olivine, organics, orthopyr, troilite, waterice'
     endif
     if n_elements(rdust) ne 1 and n_elements(rdust) ne n_elements(inputfile) then stop,'ERROR: rdust must be a scalar or a vector with length equal to the number of input files'
     if n_elements(composition) gt 1 then stop,'ERROR: composition must be a single string.'
     if n_elements(lnkfile) gt 1 then stop, 'ERROR: lnkfile must be a single string.'
     if keyword_set(opticaldepth) then stop, 'ERROR: cannot set opticaldepth keyword simultaneously with scattered or thermal keywords'
     if n_elements(ncostheta) gt 1 then stop, 'ERROR: ncostheta must be a scalar integer.'
     if n_elements(ncostheta) gt 0 then begin
        if ncostheta lt 100 then print,'WARNING: ncostheta < 100'
     endif
  endif
  if n_elements(composition) gt 0 and (size(composition,/type) ne 7) then stop, 'ERROR: composition must be a string.'
  if keyword_set(opticaldepth) and n_elements(rdust) lt 1 then stop, 'ERROR: to compute optical depth, you must supply rdust (microns).'
  if keyword_set(effrdust) and not keyword_set(thermal) and not keyword_set(scattered) then stop, 'ERROR: effrdust is only used if thermal or scattered keywords are set.'


  ;Now, perform some minor fixes to the input variables
  ;This puts everything in the correct format
  if keyword_set(fov) then fov=float(fov)
  distance=float(distance)
  if keyword_set(inclination) then inclination=float(inclination)
  if keyword_set(longitude) then longitude=float(longitude)
  if keyword_set(pa) then pa=float(pa)
  if n_elements(fov) eq 1 then begin ;if fov is a single-element array, make it a 2-element array
     tempfov = fltarr(2)
     tempfov[0] = fov
     tempfov[1] = fov
     fov = tempfov
  endif
  if keyword_set(allsky) then begin
     fov = [360.0,180.] ;all sky fov in degrees
     if not keyword_set(degrees) then begin
        fov *= 60.*60.*1000. ;convert to mas
     endif
  endif
  if keyword_set(hd) then begin
     numpixels = 1920
     pixelsize = fov[0] / numpixels
     fov[1] = 1080 * pixelsize
  endif
  if keyword_set(vga) then begin
     numpixels = 640
     pixelsize = fov[0] / numpixels
     fov[1] = 480 * pixelsize
  endif
  if keyword_set(pixelsize) and n_elements(numpixels) lt 1 then begin
     numpixels = fix(fov[0] / pixelsize)
  endif
  if n_elements(pixelsize) lt 1 and keyword_set(numpixels)  then begin
     pixelsize = fov[0] / numpixels
  endif
  if keyword_set(scaling) then begin
     if n_elements(scaling) eq 1 then userscaling = fltarr(n_elements(inputfile)) + scaling[0] $
     else userscaling = scaling
  endif
  if n_elements(rdust) eq 1 then begin
     temprdust = fltarr(n_elements(inputfile)) + rdust[0]
     rdust = temprdust
  endif
  if n_elements(distance) eq 1 then distance=distance[0] ;convert to scalar if 1-element array
  if n_elements(datatype) eq 1 then datatype=datatype[0]
  if n_elements(numpixels) eq 1 then numpixels=numpixels[0]
  if n_elements(pixelsize) eq 1 then pixelsize=pixelsize[0] 
  if n_elements(Lstar) eq 1  then Lstar=Lstar[0]
  if n_elements(Tstar) eq 1 then Tstar=Tstar[0]
  if n_elements(kurucz) eq 1 then kurucz=kurucz[0]
  if n_elements(logg) eq 1 then logg=logg[0]



  ;Set some keywords that weren't set by the user
  aitoff_flag = 0
  if n_elements(generic_oc_flag) eq 0 then generic_oc_flag = 0
  if n_elements(albedo) eq 0 then albedo = 0.0
  if n_elements(HG_flag) eq 0 then HG_flag = 0
  if n_elements(iwa) eq 0 then iwa = 0.0
  if n_elements(ext) eq 0 then ext = 0
  if n_elements(qexp) eq 0 then qexp = 0.0
  if n_elements(g) eq 0 then g = 0.0
  if n_elements(xshift) eq 0 then xshift = 0.0
  if n_elements(inclination) eq 0 then inclination = 0.0
  if n_elements(longitude) eq 0 then longitude = 0.0
  if n_elements(pa) eq 0 then pa = 0.0
  if n_elements(scaling) eq 0 then userscaling = fltarr(n_elements(inputfile)) + 1.0
  if keyword_set(opticaldepth) then opticaldepth_flag = 1 else opticaldepth_flag = 0
  if keyword_set(scattered) then scatteredlight_flag = 1 else scatteredlight_flag = 0
  if keyword_set(thermal) then thermalemission_flag = 1 else thermalemission_flag = 0
  if not keyword_set(scattered) and not keyword_set(thermal) and not keyword_set(opticaldepth) then densityhisto_flag = 1 else densityhisto_flag = 0
  if not keyword_set(lnkfile) then begin
     lnkfile='lnkfile'
     if keyword_set(composition) then lnkfile=dustmappath+'lnkfiles/'+composition+'.lnk'
  endif
  if n_elements(Tsublimate) eq 0 then Tsublimate = 1e10 ;something impossibly large so that sublimation is ignored
  if n_elements(ncostheta) eq 0 then ncostheta = long(500) ;default value
  if n_elements(Lstar) eq 0 then Lstar = 0.0 ;doesn't get used if this happens, but it must exist
  if n_elements(Tstar) eq 0 then Tstar = 0.0 ;doesn't get used if this happens, but it must exist
  if n_elements(kurucz) eq 0 then kurucz = 0 ;doesn't get used if this happens, but it must exist
  if n_elements(logg) eq 0 then logg = 0.0 ;doesn't get used if this happens, but it must exist
  if n_elements(lambda) eq 0 then lambda = 0.0 ;doesn't get used if this happens, but it must exist
  if n_elements(rdust) eq 0 then rdust = 0.0 ;doesn't get used if this happens, but it must exist
  verbose_flag = 2
  if keyword_set(reduceprint) then verbose_flag = 1
  if keyword_set(noprint) then verbose_flag = 0
  if keyword_set(aitoff) then aitoff_flag = 1
  if n_elements(effrdust) eq 0 then effrdust = 0.0
  if keyword_set(markabs) then markabs_flag = 1 else markabs_flag = 0  
  if n_elements(nmarks) eq 0 then nmarks=0
  if nmarks gt 0 then begin
     markx = marklocations[0,*]
     marky = marklocations[1,*]
     markz = marklocations[2,*]
  endif else begin
     markx=[0]
     marky=[0]
     markz=[0]
     markweight=[0]
  endelse
  if not keyword_set(AU) and not keyword_set(degrees) then pixelsizeAU = (pixelsize/1000.) * distance
  if keyword_set(AU) and not keyword_set(degrees) then pixelsizeAU = (pixelsize/1000.) * distance * 4.84813681*10.^(-6.)
  if keyword_set(AU) and keyword_set(degrees) then pixelsizeAU = (pixelsize/(60.*60.)) * distance * 4.84813681*10.^(-6.)
  if effrdust gt 1 and not keyword_set(noprint) then begin
     print,'WARNING: Many points will have Gaussian smoothing.  This may take a while.  Press any key to continue or "q" to quit...'
     a=get_kbrd()
     if strlowcase(a) eq 'q' then return
  endif
  if not keyword_set(azavg) then azavg=1


  ;Print a summary of what's happening (if desired)
  histoxsize = long(2 * floor(fov[0]/pixelsize/2))
  histoysize = long(2 * floor(fov[1]/pixelsize/2))
  if (not keyword_set(noprint)) then begin
     print,' '
     print,'dustmap.pro: creating image...'
     print,' - Imaging '+strcompress(string(fov[0],format='(F0.1)'),/remove_all)+' x '+$
           strcompress(string(fov[1],format='(F0.1)'),/remove_all)+' mas of the sky using '+$
           strcompress(string(histoxsize,format='(I)'),/remove_all)+' x '+$
           strcompress(string(histoysize,format='(I)'),/remove_all)+' pixels.'
     datastring = ' - DATATYPE = '+strcompress(string(datatype,format='(I)'),/remove_all)
     datastring += ': processing binary data files with particle IDs, positions'
     if datatype eq 2 then datastring+=', and velocities.'
     if datatype eq 3 then datastring+=', and intensities.'
     if datatype eq 4 then datastring+=', velocities, and intensities.'
     print,datastring
  
     if densityhisto_flag eq 1 then print, format='($," - Creating surface density histogram")'
     if opticaldepth_flag eq 1 then print, format='($," - Creating geometric optical depth histogram")'
     if scatteredlight_flag eq 1 and thermalemission_flag eq 1 then print, format='($," - Creating image with scattered light and thermal emission flux")'
     if scatteredlight_flag eq 1 and thermalemission_flag eq 0 then print, format='($," - Creating image with scattered light flux")'
     if scatteredlight_flag eq 0 and thermalemission_flag eq 1 then print, format='($," - Creating image with thermal emission flux")'
     if keyword_set(scattered) or keyword_set(thermal) then begin
        print, format='($," at [", F0.3)',lambda[0]
        for i=1,n_elements(lambda)-1 do print, format='($,", ", F0.3)',lambda[i]
        print, format='($,"] microns")'
     endif
     if not keyword_set(log) and not keyword_set(colorbar) then print, ' '
     if keyword_set(log) and keyword_set(colorbar) then print, ', which will be displayed on a logarithmic scale with a colorbar'
     if keyword_set(log) and not keyword_set(colorbar) then print, ', which will be displayed on a logarithmic scale'
     if not keyword_set(log) and keyword_set(colorbar) then print, ', which will be displayed with a colorbar'
     if keyword_set(thermal) or keyword_set(scattered) then begin
        print, ' - Using the lnkfile: ',lnkfile
     endif
     print, ' '
  endif
  

  ;Check to see if input files actually exist
  file_exist = file_test(inputfile)
  if total(file_exist) ne n_elements(inputfile) then stop, 'ERROR: at least 1 of the specified input files are missing.'

  
end
