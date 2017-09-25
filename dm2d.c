void dm2d(double *histo, int histoxsize, int histoysize, datatype3 *data, long ndatapoints, float fovx, float fovy, float pixelsize, float inclination, float longitude, float distance, float pa, float *lambda, long nlambda, float *Qabs_value, float *Qsca_value, float dustradius, double *Bnustar, float starradius, float *Tvsr_list, float delta_log_r, int nr, float log10minr, float Tsublimate, float *pfunc, int npfunc, float *Qwave, float *Qabs, int nQwave, float iwa, int aitoff_flag, int densityhisto_flag, int opticaldepth_flag, int scatteredlight_flag, int thermalemission_flag, int extinction_flag, float xshift, float effrdust, float distmask){



  // dm2d.c
  // Synthesizes an image of the data
  // 
  // Input:
  // dustradius: radius of the dust grains (in microns)
  // ndatapoints: the number of elements in data
  // fovx: The size of the field of view in the x-direction (in degrees)
  // fovy: The size of the field of view in the y-direction (in degrees)
  // pixelsize: The size of the pixels (in degrees)
  // inclination: (in degrees)
  // longitude: (in degrees)
  // distance:
  // pa: position angle of the disk (in degrees)
  // xshift: a translational shift in the x-direction (in AU)
  // effrdust: an effective radius for the dust for image smoothing (in AU)
  // distmask: dust within this distance of the observer is deleted (in AU)
  //
  // Output:
  // histo: the synthetic image histogram
  //
  // Created by Christopher Stark
  // NASA GSFC
  // Last updated 23 Apr 2014 by Christopher Stark


  long i, j, ilambda;
  long histosize, histoxsizeo2, histoysizeo2;
  long cosscattang_index, rindex;
  long ipx, ipxstart, ipxend, ipy, ipystart, ipyend;
  long rpx, tworpx;
  long gaussiansmoothed=0;

  float oneopixelsize;
  float pixelsizeAU, pixelsizecm, pixelareacm2;
  float fovxo2, fovyo2;
  float x, y, z;
  float xx, xx_stellocentric, yy, zz, dd;
  float r, l;
  float cosscattang;
  float pio180 = PI / 180.0;
  float one80opi = 180.0 / PI;
  float erg_per_AU2_to_Jy =  0.000446836; //conversion factor used later on
  float rotmatrix[3][3], irotmatrix[3][3], parotmatrix[3][3];
  float pfuncapplied;
  float starradiusAU, starradiusAU2;
  float dustradiuscm, dustareacm2;
  float Tdust;
  float npfuncm1o2;
  float sqrt2, oneof, alpha2, delta, cdec, oneodenom, tempxang, tempyang; //used if aitoff projection requested
  float rang2, starang, starang2;
  float oneotwosigma2;
  float totalamp, oneototalamp;
  float *amp=NULL;

  double Bnu;
  double power_received;
  double flux;
  double dustareacm2opixelareacm2;


  float *d;
  d = (float *) malloc (ndatapoints * sizeof(float));
  if (d==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating d array.\n"); exit (1);}

  float *dist;
  dist = (float *) malloc (ndatapoints * sizeof(float));
  if (dist==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating dist array.\n"); exit (1);}

  float *xang;
  xang = (float *) malloc (ndatapoints * sizeof(float));
  if (xang==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating xang array.\n"); exit (1);}

  float *yang;
  yang = (float *) malloc (ndatapoints * sizeof(float));
  if (yang==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating yang array.\n"); exit (1);}

  long *xindex;
  xindex = (long *) malloc (ndatapoints * sizeof(long));
  if (xindex==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating xindex array.\n"); exit (1);}

  long *yindex;
  yindex = (long *) malloc (ndatapoints * sizeof(long));
  if (yindex==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating yindex array.\n"); exit (1);}

  int *inlim;
  inlim = (int *) malloc (ndatapoints * sizeof(int));
  if (inlim==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating inlim array.\n"); exit (1);}

  double *Qabs_x_dustareacm2_x_convfactor;
  Qabs_x_dustareacm2_x_convfactor = (double *) malloc (nlambda * sizeof(double));
  if (Qabs_x_dustareacm2_x_convfactor==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Qabs_x_dustareacm2_x_convfactor array.\n"); exit (1);}

  double *Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor;
  Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor = (double *) malloc (nlambda * sizeof(double));
  if (Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor array.\n"); exit (1);}

  double *Bnu_x_Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor;
  Bnu_x_Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor = (double *) malloc (nlambda * sizeof(double));
  if (Bnu_x_Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Bnu_x_Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor array.\n"); exit (1);}

  double *Bnu_x_Qext_x_dustareacm2_x_oneodistance2_convfactor;
  Bnu_x_Qext_x_dustareacm2_x_oneodistance2_convfactor = (double *) malloc (nlambda * sizeof(double));
  if (Bnu_x_Qext_x_dustareacm2_x_oneodistance2_convfactor==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Bnu_x_Qext_x_dustareacm2_x_oneodistance2_convfactor array.\n"); exit (1);}


  histosize = histoxsize * histoysize;

  //Convert angles to radians
  inclination *= pio180;
  longitude *= -pio180; //longitude is a CW angle we want to rotate points CCW
  fovx *= pio180;
  fovy *= pio180;
  pixelsize *= pio180;
  //pa += 90.0;
  pa *= -pio180; //pa is a CW angle and we want to rotate points CCW
  iwa *= pio180/3600000.0; //occulted radius is in mas, we want radians
  
  //Here are some variables defined so that computations are performed more quickly later on
  //fovxo2 = fovx * 0.5;
  //fovyo2 = fovy * 0.5;
  oneopixelsize = 1. / pixelsize;
  histoxsizeo2 = histoxsize / 2;
  histoysizeo2 = histoysize / 2;
  npfuncm1o2 = 0.5*(npfunc-1);
  dustradiuscm = dustradius * 0.0001; //convert microns to cm
  dustareacm2 = PI * dustradiuscm * dustradiuscm; //dust area in cm^2
  pixelsizeAU = distance * tan(pixelsize);
  pixelsizecm = pixelsizeAU * 1.49598e13;
  pixelareacm2 = pixelsizecm * pixelsizecm;
  dustareacm2opixelareacm2 = (double) ((double) dustareacm2) / ((double) pixelareacm2);
  starradiusAU = starradius * 0.00464913; //convert from solar radii to AU
  starradiusAU2 = starradiusAU * starradiusAU; //simply the square of the above value
  starang = atan(starradiusAU / distance); //the angular radius of the star in radians
  starang2 = starang * starang; //the angular radius of the star in radians

  for(ilambda=0;ilambda<nlambda;ilambda++) {
    Qabs_x_dustareacm2_x_convfactor[ilambda] = ((double) Qabs_value[ilambda]) * dustareacm2  * erg_per_AU2_to_Jy; //removes some multiplication later on
    Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor[ilambda] = ((double) Qsca_value[ilambda]) * PI * dustareacm2 * starradiusAU2 * erg_per_AU2_to_Jy; //removes some multipl. later on
  }

  //Rotate axes so that z-axis of disk lies along the viewing vector and then rotate by the position angle.
  //This is written to be understandable, not fast.  This is small potatoes computationally compared to
  //the rest of the code.
  //First, create the initial rotation matrix
  irotmatrix[0][0] = cos(longitude);                      irotmatrix[0][1] = sin(longitude);                      irotmatrix[0][2] = 0;
  irotmatrix[1][0] = -sin(longitude) * cos(inclination);  irotmatrix[1][1] = cos(inclination) * cos(longitude);   irotmatrix[1][2] = sin(inclination);
  irotmatrix[2][0] = sin(longitude) * sin(inclination);   irotmatrix[2][1] = -sin(inclination) * cos(longitude);  irotmatrix[2][2] = cos(inclination);


  //Now create the position angle rotation matrix
  parotmatrix[0][0] = cos(pa);   parotmatrix[0][1] = sin(pa);  parotmatrix[0][2] = 0;
  parotmatrix[1][0] = -sin(pa);  parotmatrix[1][1] = cos(pa);  parotmatrix[1][2] = 0;
  parotmatrix[2][0] = 0;         parotmatrix[2][1] = 0;        parotmatrix[2][2] = 1;
  //Multiply the two rotation matrices
  rotmatrix[0][0] = irotmatrix[0][0]*parotmatrix[0][0] + irotmatrix[1][0]*parotmatrix[0][1];
  rotmatrix[0][1] = irotmatrix[0][1]*parotmatrix[0][0] + irotmatrix[1][1]*parotmatrix[0][1];
  rotmatrix[0][2] = irotmatrix[0][2]*parotmatrix[0][0] + irotmatrix[1][2]*parotmatrix[0][1];
  rotmatrix[1][0] = irotmatrix[0][0]*parotmatrix[1][0] + irotmatrix[1][0]*parotmatrix[1][1];
  rotmatrix[1][1] = irotmatrix[0][1]*parotmatrix[1][0] + irotmatrix[1][1]*parotmatrix[1][1];
  rotmatrix[1][2] = irotmatrix[0][2]*parotmatrix[1][0] + irotmatrix[1][2]*parotmatrix[1][1];
  rotmatrix[2][0] = irotmatrix[2][0];
  rotmatrix[2][1] = irotmatrix[2][1];
  rotmatrix[2][2] = irotmatrix[2][2];
  //Now rotate the data points
  for(i=0;i<ndatapoints;i++){
    x = rotmatrix[0][0] * data[i].x + rotmatrix[0][1] * data[i].y + rotmatrix[0][2] * data[i].z;
    y = rotmatrix[1][0] * data[i].x + rotmatrix[1][1] * data[i].y + rotmatrix[1][2] * data[i].z;
    z = rotmatrix[2][0] * data[i].x + rotmatrix[2][1] * data[i].y + rotmatrix[2][2] * data[i].z;
    data[i].x = x;
    data[i].y = y;
    data[i].z = z;
  }

  //Shift the data points if we are viewing from off the z-axis
  if(xshift != 0) {for(i=0;i<ndatapoints;i++) data[i].x -= xshift;}
  
  //Calculate the component of distance along the line of sight, d, to each point and the angle made with respect to the viewer
  for(i=0;i<ndatapoints;i++){
    d[i] = distance - data[i].z;
    dist[i] = sqrt(d[i]*d[i] + data[i].x*data[i].x + data[i].y*data[i].y);
    xang[i] = atan2(data[i].x, d[i]); //longitude
    yang[i] = atan2(data[i].y, sqrt(data[i].x*data[i].x + d[i]*d[i])); //latitude
  }


  
  //Make histo indeces
  for(i=0;i<ndatapoints;i++){
    xindex[i] = floor(floor(xang[i] * oneopixelsize) + histoxsizeo2);
    yindex[i] = floor(floor(yang[i] * oneopixelsize) + histoysizeo2);
    if( (yindex[i] < histoysize) && (xindex[i] < histoxsize) && (yindex[i] >= 0) && (xindex[i] >= 0) ) {
      inlim[i] = 1; //mark as "in bounds"
      if(iwa > 0) {  //remove those that are within occulting spot
	if(sqrt(xang[i] * xang[i] + yang[i] * yang[i]) < iwa) inlim[i] = 0;
      }
    }
    else inlim[i] = 0;
  }


  //If an Aitoff projection is requested, project xang and yang
  if(aitoff_flag == 1) {
    sqrt2 = sqrt(2.0);
    oneof = PI / (2.0 * sqrt2);
    for(i=0;i<ndatapoints;i++){
      if(inlim[i]) {
	alpha2 = xang[i] * 0.5;
	delta = yang[i];
	cdec = cos(delta);    
	oneodenom = 1.0 / sqrt(1.0 + cdec * cos(alpha2));
	tempxang = cdec * sin(alpha2) * 2. * sqrt2 * oneodenom;
	tempyang = sin(delta) * sqrt2 * oneodenom;
	tempxang = tempxang * oneof;
	tempyang = tempyang * oneof;
	
	xang[i] = tempxang; //assign new longitude for aitoff projection
	yang[i] = tempyang; //assign new latitude for aitoff projection
	xindex[i] = floor(floor(xang[i] * oneopixelsize) + histoxsizeo2); //remake indeces
	yindex[i] = floor(floor(yang[i] * oneopixelsize) + histoysizeo2);
      }
    }
  }
  

  //CASE 1: density histo, no radiative calulcations required
  if(densityhisto_flag){
    for(i=0;i<ndatapoints;i++) {
      if(inlim[i] && dist[i] >= distmask) histo[xindex[i] + yindex[i] * histoxsize] += data[i].intensity;
    }
  }
  else if(opticaldepth_flag){
    //CASE 2: optical depth histo
    for(i=0;i<ndatapoints;i++) {
      if(inlim[i] && dist[i] >= distmask) histo[xindex[i] + yindex[i] * histoxsize] += data[i].intensity * dustareacm2opixelareacm2;
    }
    
  }
  else {


    //CASE 3: scattered light and/or thermal emission requested

    for(ilambda=0;ilambda<nlambda;ilambda++) {
      Bnu_x_Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor[ilambda] = Bnustar[ilambda] * Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor[ilambda];
      Bnu_x_Qext_x_dustareacm2_x_oneodistance2_convfactor[ilambda] = Bnustar[ilambda] * ((double) Qabs_value[ilambda] + (double) Qsca_value[ilambda]) * dustareacm2 * erg_per_AU2_to_Jy / (distance * distance);
    }




    //PART 1: scattered light
    if(scatteredlight_flag) {
      for(i=0;i<ndatapoints;i++){
	if(inlim[i] && dist[i] >= distmask){
	  rang2 = xang[i] * xang[i] + yang[i] * yang[i];
	  if ( !( (rang2 < starang2) && (data[i].z < 0) ) ) { //If dust is not blocked by star...
	    xx_stellocentric = (data[i].x + xshift) * (data[i].x + xshift); //stellocentric x value (unshifted)
	    xx = data[i].x * data[i].x; //this has been shifted
	    yy = data[i].y * data[i].y;
	    zz = data[i].z * data[i].z;
	    dd = d[i] * d[i];
	    r = sqrt(xx_stellocentric + yy + zz);
	    rindex = floor((log10(r) - log10minr) / delta_log_r); //find the entry in Tvsr_list that corresponds to this distance
	    if(rindex < nr) Tdust = Tvsr_list[rindex]; //Use the previously tabulated values of temperature
	    else Tdust = Tvsr_list[nr-1];
	    //else dmtemperaturecalc(starradius, Tstar, &r, 1, Qwave, Qabs, nQwave, &Tdust);
	    if(Tdust < Tsublimate){
	      l = sqrt(xx + yy + dd);
	  
   	      //Here is something that looks messy, but is fast:
	      cosscattang = (data[i].z * d[i] - (data[i].x + xshift) * data[i].x - yy) / (r * l);
	  
	      //Here it is explained in comments:
	      //r_x = data[i].x + xshift; //x comp of vector joining star to dust grain
	      //r_y = data[i].y; //y comp of vector joining star to dust grain
	      //r_z = data[i].z; //z comp of vector joining star to dust grain
	      //l_x = -data[i].x; //x comp of vector joining dust grain to viewing position
	      //l_y = -data[i].y; //y comp of vector joining dust grain to viewing position
	      //l_z = d[i]; //z comp of vector joining dust grain to viewing position
	      //cosscattang = (r_x * l_x + r_y * l_y + r_z * l_z) / (r*l) //cos(scattering angle) = dot product of unit vector joining star to dust grain and unit vector joining dust grain to viewing position
	      cosscattang_index = floor(cosscattang * npfuncm1o2 + npfuncm1o2);


	      rpx = (long) ((effrdust / d[i]) * oneopixelsize); //effective imaging radius of dust grain in pixels
	      if(rpx > 0) {
		gaussiansmoothed=1;
		tworpx = 2 * rpx;
		oneotwosigma2 = 1.0 / (tworpx * rpx); //sigma value in pixels (based on distance and size)
		ipxstart = xindex[i] - tworpx;
		ipxend = xindex[i] + tworpx;
		ipystart = yindex[i] - tworpx;
		ipyend = yindex[i] + tworpx;
		ipxstart = (ipxstart > 0) ? ipxstart : 0; //ipxstart cannot be negative
		ipystart = (ipystart > 0) ? ipystart : 0; //ipystart cannot be negative
		ipxend = (ipxend < histoxsize-1) ? ipxend : histoxsize-1; //ipxend cannot be greater than histoxsize
		ipyend = (ipyend < histoysize-1) ? ipyend : histoysize-1; //ipxend cannot be greater than histoxsize
		float *tmp;
		tmp = (float *) realloc (amp, (ipxend-ipxstart+1)*(ipyend-ipystart+1) * sizeof(float));
		amp = tmp;
		if (amp==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating amp array.\n"); exit (1);}
		totalamp = 0.0;
		for(ipx=ipxstart;ipx<=ipxend;ipx++){
		  for(ipy=ipystart;ipy<=ipyend;ipy++){
		    amp[(ipx-ipxstart) * (ipyend-ipystart+1) + ipy-ipystart] = exp( - ((xindex[i]-ipx)*(xindex[i]-ipx) + (yindex[i]-ipy)*(yindex[i]-ipy)) * oneotwosigma2 );  
		    totalamp += amp[(ipx-ipxstart) * (ipyend-ipystart+1) + ipy-ipystart];
		  }
		}
		oneototalamp = 1.0 / totalamp;
		for(ipx=ipxstart;ipx<=ipxend;ipx++){
		for(ipy=ipystart;ipy<=ipyend;ipy++){
		  amp[(ipx-ipxstart) * (ipyend-ipystart+1) + ipy-ipystart] *= oneototalamp;
		}
		}
	      }
	      
	      for(ilambda=0;ilambda<nlambda;ilambda++) {
		//Note: for a 2D array entry c[ix,iy], the 1D index = ix + iy*xsize
		//Here, ix = ilambda, iy = cosscattang_index, and xsize = nlambda
		pfuncapplied = pfunc[ilambda+nlambda*cosscattang_index];
		power_received = data[i].intensity * Bnu_x_Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor[ilambda] / (xx_stellocentric + yy + zz);
		//The above looks messy, but it's just
		// power_received = (# of grains) * (PI * Bnu) * Qsca * (PI * rdust * rdust) * (rstar / distance)^2
		// There's also a conversion factor thrown in the above calculation to save computation time, so the above
		// is only power_received if you divide by convfactor.
		flux = (power_received * pfuncapplied) / (xx + yy + dd);
		if(rpx <= 0) { //Single pixel case, no gaussian weighting
		  histo[xindex[i] + yindex[i] * histoxsize + ilambda * histosize] += flux; //flux in Jy
		}
		else { //Multiple pixels to fill according to a gaussian weighting
		  for(ipx=ipxstart;ipx<=ipxend;ipx++){
		    for(ipy=ipystart;ipy<=ipyend;ipy++){
		      histo[ipx + ipy * histoxsize + ilambda * histosize] += flux * amp[(ipx-ipxstart) * (ipyend-ipystart+1) + ipy-ipystart];
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    } //end of part 1, the scattered light stuff





    //PART 2: thermal emission
    if(thermalemission_flag){
      for(i=0;i<ndatapoints;i++){
	if(inlim[i] && dist[i] >= distmask){
	  rang2 = xang[i] * xang[i] + yang[i] * yang[i];
	  if ( !( (rang2 < starang2) && (data[i].z < 0) ) ) { //If dust is not blocked by star...
	    xx_stellocentric = (data[i].x + xshift) * (data[i].x + xshift); //stellocentric x value (unshifted)
	    xx = data[i].x * data[i].x; //this has been shifted
	    yy = data[i].y * data[i].y;
	    dd = d[i] * d[i];
	    r = sqrt(xx_stellocentric + yy + data[i].z*data[i].z); //circumsmtellar distance of dust grain
	    rindex = floor((log10(r) - log10minr) / delta_log_r); //find the entry in Tvsr_list that corresponds to this distance
	    if(rindex < nr) Tdust = Tvsr_list[rindex]; //Use the previously tabulated values of temperature
	    else Tdust = Tvsr_list[nr-1];
	    //else dmtemperaturecalc(starradius, Tstar, &r, 1, Qwave, Qabs, nQwave, &Tdust);
	    if(Tdust < Tsublimate){
	      rpx = (long) ((effrdust / d[i]) * oneopixelsize); //effective imaging radius of dust grain in pixels
	      if(rpx > 0) {
		gaussiansmoothed=1;
		tworpx = 2 * rpx;
		oneotwosigma2 = 1.0 / (tworpx * rpx); //sigma value in pixels (based on distance and size)
		ipxstart = xindex[i] - tworpx;
		ipxend = xindex[i] + tworpx;
		ipystart = yindex[i] - tworpx;
		ipyend = yindex[i] + tworpx;
		ipxstart = (ipxstart > 0) ? ipxstart : 0; //ipxstart cannot be negative
		ipystart = (ipystart > 0) ? ipystart : 0; //ipystart cannot be negative
		ipxend = (ipxend < histoxsize-1) ? ipxend : histoxsize-1; //ipxend cannot be greater than histoxsize
		ipyend = (ipyend < histoysize-1) ? ipyend : histoysize-1; //ipxend cannot be greater than histoxsize
		float *tmp;
		tmp = (float *) realloc (amp, (ipxend-ipxstart+1)*(ipyend-ipystart+1) * sizeof(float));
		amp = tmp;
		if (amp==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating amp array.\n"); exit (1);}
		totalamp = 0.0;
		for(ipx=ipxstart;ipx<=ipxend;ipx++){
		  for(ipy=ipystart;ipy<=ipyend;ipy++){
		    amp[(ipx-ipxstart) * (ipyend-ipystart+1) + ipy-ipystart] = exp( - ((xindex[i]-ipx)*(xindex[i]-ipx) + (yindex[i]-ipy)*(yindex[i]-ipy)) * oneotwosigma2 );  
		    totalamp += amp[(ipx-ipxstart) * (ipyend-ipystart+1) + ipy-ipystart];
		  }
		}
		oneototalamp = 1.0 / totalamp;
		for(ipx=ipxstart;ipx<=ipxend;ipx++){
		  for(ipy=ipystart;ipy<=ipyend;ipy++){
		    amp[(ipx-ipxstart) * (ipyend-ipystart+1) + ipy-ipystart] *= oneototalamp;
		  }
		}
	      }


	      for(ilambda=0;ilambda<nlambda;ilambda++) {
		calcBnu(Tdust, lambda[ilambda], &Bnu);

		flux = (data[i].intensity * Bnu * Qabs_x_dustareacm2_x_convfactor[ilambda]) / (xx + yy + dd); //flux in Jy
		//The above looks messy, but it's just
		// flux = (# of grains) * (PI * Bnu) * (rdust / distance)^2, plus a conversion factor to convert erg AU^-2 to Jy
		if(effrdust <= 0) { //Single pixel case, no gaussian weighting
		  histo[xindex[i] + yindex[i] * histoxsize + ilambda * histosize] += flux; //flux in Jy
		}
		else { //Multiple pixels to fill according to a gaussian weighting
		  for(ipx=ipxstart;ipx<=ipxend;ipx++){
		    for(ipy=ipystart;ipy<=ipyend;ipy++){
		      histo[ipx + ipy * histoxsize + ilambda * histosize] += flux * amp[(ipx-ipxstart) * (ipyend-ipystart+1) + ipy-ipystart];
		    }
		  }
		}
	      


	      }
	    }
	  }
	}
      }
    } //end of part 2, the thermal emission stuff





    //PART 3: extinction
    //Note: This just adds a negative flux to the location of the star
    //based on the scattering/absorption of dust in front of the star.
    //Note: This is approximate.  It is not a good approximation when
    //dust is close to the star (r ~< 10 R_star), or when the observer
    //is very close to the star (distance ~< 10 R_star).
    //Note: There is no Gaussian smoothing on extinction--the negative flux
    //associated with extinction always gets placed at the projected
    //location of the dust.
    if(extinction_flag) {
      for(i=0;i<ndatapoints;i++){
	if(inlim[i] && dist[i] >= distmask){
	  rang2 = xang[i] * xang[i] + yang[i] * yang[i];
	  if ( (rang2 < starang2) && (data[i].z > 0) ) { //If dust is in front of star...

	    for(ilambda=0;ilambda<nlambda;ilambda++) {
	      //Here we have a grain that is blocking light from the star.  So, 
	      //to account for the dimming of the star, we need to subtract the
	      //flux intercepted by the grains.
	      if( data[i].z > d[i] * 0.01 ) { //circumstellar distance is not much less than distance to observer, do a more exact calculation
		flux = - data[i].intensity * Bnu_x_Qext_x_dustareacm2_x_oneodistance2_convfactor[ilambda] * (data[i].z + d[i]) * (data[i].z + d[i]) / (d[i] * d[i]);
	      } else {
		flux = - data[i].intensity * Bnu_x_Qext_x_dustareacm2_x_oneodistance2_convfactor[ilambda];
	      }
	      histo[xindex[i] + yindex[i] * histoxsize + ilambda * histosize] += flux; //flux in Jy
	    }

	  }
	}
      }
    } //end of part 3, the extinction stuff





  } //end of scattered light/thermal emission else command



  

  //Free some locally-allocated arrays that aren't used externally
  if (gaussiansmoothed==1) free(amp);
  free(d); free(dist); free(xang); free(yang); free(xindex); free(yindex); free(inlim);
  free(Qabs_x_dustareacm2_x_convfactor);
  free(Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor);
  free(Bnu_x_Qsca_x_PI_x_dustareacm2_x_starradiusAU2_x_convfactor);
  free(Bnu_x_Qext_x_dustareacm2_x_oneodistance2_convfactor);
}
