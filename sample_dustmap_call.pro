pro sample_dustmap_call


  
  inputfile = 'sample_inputfile.dat'
  distance = 10.                ;pc
  fov = 1000.                   ;mas
  numpixels = 200
  datatype = 1
  inclination = 60. ;degrees
  longitude = 90. ;degrees
  pa = 45. ;degrees
  ;opticaldepth = 1
  rdust = 1.
  thermal = 1
  lambda = 0.5
  Lstar = 1.
  Tstar = 5778.
  kurucz = 1
  logg=4.5
  composition = 'astrosil'
  scattered = 1
  scaling=1e30
  effrdust=0
  ;hg = 1
  ;g = 0.1
  ;HG_flag=1
  
!PATH=!PATH+':'+Expand_Path('+~/github/IDLAstro/')
!PATH=!PATH+':'+Expand_Path('+~/github/coyote/')

  
  
  dustmap, INPUTFILE, image, $	
           DISTANCE=distance, FOV=fov, PIXELSIZE=pixelsize, NUMPIXELS=numpixels, DATATYPE=datatype, $
           INCLINATION=inclination, LONGITUDE=longitude, PA=pa, $
           AU=au, DEGREES=degrees, $
           OPTICALDEPTH=opticaldepth, SCATTERED=scattered, THERMAL=thermal, LAMBDA=lambda, $
           TSTAR=Tstar, LSTAR=lstar, KURUCZ=kurucz, LOGG=logg, FSTAR=fstar, $
           SCALING=scaling, RDUST=rdust, COMPOSITION=composition, LNKFILE=lnkfile, TSUBLIMATE=tsublimate, $
           QEXP=qexp, ALBEDO=albedo, OC_GENERIC_FLAG=oc_generic_flag, G=g, HG_flag=HG_flag, EXT=ext, $
           PFUNC=pfunc, COSTHETA=costheta, NCOSTHETA=ncostheta, QABS=qabs, QSCA=qsca, $
           AITOFF=aitoff, ALLSKY=allsky, VGA=vga, HD=hd,effrdust=effrdust

  ;window,0,xsize=500,ysize=500
  ;loadct,39
  ;plot,lambda,fstar,/xlog,/ylog
  ;plot,lambda,total(total(image,1),1),/xlog,/ylog
  ;tvscl,image[*,*,0]^0.25

SAVE, /VARIABLES, FILENAME = 'sample_dustmap_call.sav'
 writefits,inputfile+".fits",image 
end
