#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>

#define PI 3.141592653589793238

#include "dmdatatypes.h"

#include "dmgetdata.c"
#include "calcBnu.c"
#include "getKurucz.c"
#include "dmgetBnustar.c"
#include "hgphasefunction.c"
#include "dmtemperaturecalc.c"
#include "dmgetoc.c"
#include "dm2d.c"
#include "dmmarklocs.c"

int dustmap_c(int argc, void *argv[]) {


  //dustmap_c.c
  //Last updated 26 June 2014
  //Created by Christopher Stark
  //NASA GSFC

  //To compile this code via IDL:
  //IDL> a='dustmap_c'
  //IDL> b='DIRECTORY_WHERE_DUSTMAP_FILES_ARE_LOCATED'
  //IDL> make_dll,a,a,a,input_directory=b,output_directory=b,cc='gcc -c -fPIC -O3 %C -o %O'
  // Note: the -O3 is an optimization flag that will dramatically improve runtime of the code

  //-------------------- INPUT --------------------

  //Here are the variables passed by the IDL wrapper.
  //Note that IDL is not neccessary--you can create your own wrapper if you want.
  //Note: all scalars are interpreted below as passed-by-value, so dustmap.c won't change them
  int j = 0;
  double *histo = (double *) argv[j]; j++;     //The final image, passed by reference
  int histoxsize = *(int *) argv[j]; j++;
  int histoysize = *(int *) argv[j]; j++;
  char *inputfilelist = (char *) argv[j]; j++; //Concatenated strings of input files from IDL
  int numfiles = *(int *) argv[j]; j++;        //number of input file strings to expect
  int datatype = *(int *) argv[j]; j++;        //The data type (see dmgetdata.c)
  float distance = *(float *) argv[j]; j++;     //Distance to the star in pc
  float fovx = *(float *) argv[j]; j++;         //Field of view in the x-direction in mas
  float fovy = *(float *) argv[j]; j++;         //Field of view in the y-direction in mas
  float pixelsize = *(float *) argv[j]; j++;    //Pixel size in mas
  float inclination = *(float *) argv[j]; j++;  //Inclination of the system
  float longitude = *(float *) argv[j]; j++;    //Orbital phase of the system
  float pa = *(float *) argv[j]; j++;           //Position angle of the system in degrees
  float *lambda = (float *) argv[j]; j++;       //Wavelength of image in microns
  int nlambda = *(int *) argv[j]; j++;
  float Lstar = *(float *) argv[j]; j++;        //Luminosity of star in solar luminosities
  float Tstar = *(float *) argv[j]; j++;        //Effective temperature of star in degrees K
  float logg = *(float *) argv[j]; j++;         //Surface gravity of star
  int kurucz_flag = *(int *) argv[j]; j++;      //Flag signaling a Kurucz stellar atmosphere model is desired
  char *kurucz_file = (char *) argv[j]; j++;
  float *fstar = (float *) argv[j]; j++;       //Output stellar flux in Jy
  float *dustradius = (float *) argv[j]; j++;   //Radius of dust grain in microns, array with length equal to numfiles
  float Tsublimate = *(float *) argv[j]; j++;   //Sublimation temperature of the dust in degrees K
  char *lnkfile = (char *) argv[j]; j++;        //File name containing indeces of refraction vs. lambda for dust
  double *scaling = (double *) argv[j]; j++;    //Factor to scale up the amount of dust
  float iwa = *(float *) argv[j]; j++;          //Radius of an occulting spot to place on image
  int aitoff_flag = *(int *) argv[j]; j++;
  int densityhisto_flag = *(int *) argv[j]; j++;    //Signals a density histo is to be made, ignores image flags below
  int opticaldepth_flag = *(int *) argv[j]; j++;    //Signals a geometric optical depth histo is to be made
  int scatteredlight_flag = *(int *) argv[j]; j++;  //Signals image to include scattered light calculation
  int thermalemission_flag = *(int *) argv[j]; j++; //Signals image to include thermal emission calculation
  int verbose_flag = *(int *) argv[j]; j++;         //Print a lot of info
  int generic_oc_flag = *(int *) argv[j]; j++;         //Print a lot of info
  float generic_oc_exp = *(float *) argv[j]; j++;   //An input exponent for the Qabs/Qsca power law; non-zero value overides values from lnkfile
  float generic_oc_albedo = *(float *) argv[j]; j++;           //An input albedo to adjust generic Qsca
  int HG_flag = *(int *) argv[j]; j++;         //Print a lot of info
  float HG_g = *(float *) argv[j]; j++;     //An input g value; non-zero value overides values from lnkfile
  int extinction_flag = *(int *) argv[j]; j++;      //Signals that extinction of the star light should be taken into account
  float xshift = *(float *) argv[j]; j++;           //A left-right shift in the viewing direction (AU)
  float effrdust = *(float *) argv[j]; j++;         //The effective radius of the dust grains for gaussian image smoothing (AU)
  float *markx = (float *) argv[j]; j++;      //array of x coordinates of points to mark
  float *marky = (float *) argv[j]; j++;      //array of y coordinates of points to mark
  float *markz = (float *) argv[j]; j++;      //array of z coordinates of points to mark
  float *markweight = (float *) argv[j]; j++; //array of weights of the marks
  int nmarks = *(int *) argv[j]; j++;        //Number of entries in above mark arrays
  int markabs_flag = *(int *) argv[j]; j++;  //Whether weights are absolute (default is relative)
  int az_sym = *(int *) argv[j]; j++;       //A flag to force the disk to have azimuthal symmetry, value equals
                                            //the number of rotation angles w/in 360 deg. to average over
  float distmask = *(float *) argv[j]; j++; //dust within this radius from the observer is deleted
  int npfunc = *(int *) argv[j]; j++; //number of cos(theta) values
  float *pfunc = (float *) argv[j]; j++; //used internally for the phase function
  float *costheta = (float *) argv[j]; j++; //output costheta
  float *pfunc_out = (float *) argv[j]; j++; //output phase function
  float *qabs_out = (float *) argv[j]; j++; //output qabs
  float *qsca_out = (float *) argv[j]; j++; //output qsca
  int nrdust = *(int *) argv[j]; j++; //number of unique dust sizes
  //-----------------------------------------------





  //-------------------- DEFINITIONS --------------------

  //Definitions that actually control something, but are not input to the code.
  //These shouldn't need to be changed, but are available here if you need to.
  int nr=2000; //Number of entries in r_list, which controls the number of temperatures calcualted vs r

  //Definitions that shouldn't be touched
  int nQwave;
  int ifile;
  int irdust;
  long ndatapoints;
  int i, ilambda, deltaj, deltajmax, ipfunc;
  char **inputfiles;
  float prevdustradius;
  float delta_log_r, minr, maxr, log10minr; //delta log(r) for r_list (in log(AU))
  float twoPIs;
  float Rstar, kteff, klogg;

  float *r_list;
  r_list = (float *) malloc (nr * sizeof(float));
  if (r_list==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating r array.\n"); exit (1);}
  
  float *Tvsr_list;
  Tvsr_list = (float *) malloc (nr * sizeof(float));
  if (Tvsr_list==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating temperature array.\n"); exit (1);}
  
  float *Qabs_value; 
  Qabs_value = (float *) malloc (nlambda * sizeof(float));
  if (Qabs_value==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Qabs_value array.\n"); exit (1);}

  float *Qsca_value; 
  Qsca_value = (float *) malloc (nlambda * sizeof(float));
  if (Qsca_value==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Qsca_value array.\n"); exit (1);}

  double *Bnustar; //create an array to contain stellar spectrum
  Bnustar = (double *) malloc (nlambda * sizeof(double));
  if (Bnustar==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Bnustar array for lambda.\n"); exit (1);}

  //Create some pointers that are allocated later in other functions
  datatype3 *data=NULL;
  datatype3 *markdata;
  float *Qwave=NULL, *Qabs=NULL, *Qsca=NULL, *Qpfunc=NULL;
  
  for(i=0;i<npfunc;i++) costheta[i] = i*2.0/(npfunc-1) - 1.0; //calculate costheta array -- only has to be done once, needs to have constant dcostheta
  //-----------------------------------------------------








  //-------------------- REWORK SOME INPUT --------------------

  //Here we reinterpret the list of input files given by IDL.
  //IDL just gives us a series of concatenated strings, separated by a null terminating
  //character.  It's more elegant to deal with an array of strings, so we manipulate it a bit.

  //Find length of the longest filename, deltajmax
  j=0;
  deltajmax=0;
  for(i=0;i<numfiles;i++) {
    deltaj = (1 + strcspn (&inputfilelist[j],"\0"));
    j += deltaj;
    if(deltaj > deltajmax) deltajmax=deltaj;
  }
  deltajmax += 10;

  //Now create an array of inputfiles, each with length deltajmax
  inputfiles = (char **) malloc(numfiles * sizeof(char *));   // allocate storage for an array of pointers   
  if (inputfiles==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating inputfile array.\n"); exit (1);}
  for (i = 0; i < numfiles; i++) {
    inputfiles[i] = (char *) malloc(deltajmax * sizeof(char));   // for each pointer, allocate storage for an array of chars
    if (inputfiles[i] == NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating inputfile array.\n"); exit (1);} 
  }

  //Now print the concatenated strings into an array of strings
  j=0;
  for(i=0;i<numfiles;i++) {
    sprintf(inputfiles[i],"%s",&inputfilelist[j]);    
    j += (1 + strcspn (&inputfilelist[j],"\0"));
  }


  //Calculate delta_log_r
  maxr = 2. * distance * (fovx/1000.) * 0.5 * sqrt(3.0); //The maximum circumstellar distance in AU for temperature calculations
  if(distance < 0.1) maxr = 300.0;
  minr = 0.01; //calculations good down to 1/100th of an AU
  log10minr = log10(minr); //used later
  delta_log_r = (log10(maxr) - log10minr) / nr;

  //Do some conversions
  effrdust = effrdust * (pixelsize/1000.) * distance; //convert from fraction of pixel at r=distance to AU
  fovx /= 3600000.0; //convert to degrees
  fovy /= 3600000.0; //convert to degrees
  pixelsize /= 3600000.0; //convert to degrees
  distance *= 206265.0; //convert to AU
  distmask *= 206265.0; //convert to AU

  //-----------------------------------------------------------





  //-------------------- PRINT INFO IF DESIRED --------------------

  if (verbose_flag >= 2) {
    printf("\nCalling dustmap_c...\n");
    printf("Input values:\n");
    printf("  inputfiles  =\n");
    for(i=0;i<numfiles;i++) printf("    - %s\n",inputfiles[i]);
    printf("  numfiles    = %d\n",numfiles);
    printf("  datatype    = %d\n",datatype);
    printf("  distance    = %f pc\n",distance/206265.0);
    printf("  fovx        = %0.1f mas\n",fovx*3600000.0);
    printf("  fovy        = %0.1f mas\n",fovy*3600000.0);
    printf("  pixelsize   = %0.1f mas\n",pixelsize*3600000.0);
    printf("  inclination = %0.1f degrees\n",inclination);
    printf("  longitude   = %0.1f degrees\n",longitude);
    printf("  pa          = %0.1f degrees\n",pa);
    printf("  lambda      =");
    for(ilambda=0;ilambda<nlambda;ilambda++) printf(" %0.2f",lambda[ilambda]);
    printf(" microns\n");
    printf("  Lstar       = %0.2f Solar luminosities\n",Lstar);
    printf("  Tstar       = %0.1f degrees K\n",Tstar);
    printf("  logg        = %0.2f \n",logg);
    printf("  kurucz_flag = %d\n",kurucz_flag);
    printf("  kurucz_file = %s\n",kurucz_file);
    printf("  dustradius  =");
    for(i=0;i<numfiles;i++) printf(" %f",dustradius[i]);
    printf(" microns\n");
    printf("  Tsublimate  = %0.1f degrees K\n",Tsublimate);
    printf("  lnkfile     = %s\n",lnkfile);
    printf("  scaling     =");
    for(i=0;i<numfiles;i++) printf(" %e",scaling[i]);
    printf("\n");
    printf("  IWA         = %0.3f mas\n",iwa);
    printf("  aitoff flag = %d\n",aitoff_flag);
    printf("  dens. flag  = %d\n",densityhisto_flag);
    printf("  o.d. flag   = %d\n",opticaldepth_flag);
    printf("  scatt. flag = %d\n",scatteredlight_flag);
    printf("  therm. flag = %d\n",thermalemission_flag);
    printf("  ncostheta   = %d\n",npfunc);
    if(generic_oc_flag != 0){
      printf("  Using generic optical constants...\n");
      printf("  gen. oc exp = %0.3f\n",generic_oc_exp);
      printf("  gen. oc alb.= %0.3f\n\n",generic_oc_albedo);
    } else printf("  Using Mie Theory optical constants...\n");
    if(HG_flag != 0){
      printf("  Using Henyey Greenstein SPF...\n");
      printf("  HG g        = %0.3f\n",HG_g);
    } else printf("  Using Mie Theory SPF...\n");
  }

  //---------------------------------------------------------------



  //-------------------- STELLAR SPECTRUM -------------------------
  //If desired, load the stellar spectrum
  if (scatteredlight_flag == 1 || thermalemission_flag == 1) {
    dmgetBnustar(Tstar, logg, Lstar, kurucz_flag, kurucz_file, lambda, nlambda, &Rstar, &Bnustar[0], &kteff, &klogg); //calculate Bnustar for all lambda
    for(i=0;i<nlambda;i++) fstar[i] = PI*(Bnustar[i]*1e23)*(Rstar*0.00464913034/distance)*(Rstar*0.00464913034/distance); //in Jy
    if(verbose_flag >= 1) {
      if(kurucz_flag) printf("  Using a Kurucz model stellar atmosphere with Teff = %.0f, Lstar = %.1f, log(g) = %.1f, Rstar = %.2f\n", kteff, Lstar, klogg, Rstar);
      else printf("  Using a blackbody model stellar atmosphere with Teff = %.0f, Lstar = %.1f, Rstar = %.2f\n", Tstar, Lstar, Rstar);
    }
  }
  //---------------------------------------------------------------
  
  
  
  //-------------------- MAIN LOOP --------------------

  prevdustradius=0.0;
  irdust = -1;
  for(ifile=0;ifile<numfiles;ifile++){ //loop through all input files
    
    //Load the data
    dmgetdata(inputfiles[ifile], datatype, scaling[ifile], &data, &ndatapoints, az_sym);

    //------------ WARNINGS ----------
    if( (scaling[ifile] < 1e10) && (datatype == 1 || datatype == 2) && (thermalemission_flag==1 || scatteredlight_flag==1) ){
      printf("***** WARNING: No intensity values were detected and your scaling parameter is small.  Results may suffer from floating point underflow.\n");
    }

    //Print progress
    if(verbose_flag >= 1) {printf("\rProgress: %3.0f%% complete", (100.0 * ifile)/numfiles); fflush(stdout);}

    //Load optical constants
    if(scatteredlight_flag == 1 || thermalemission_flag == 1) {


      //Only load the optical constants if something has changed
      if(dustradius[ifile] != prevdustradius) { //must calculate OCs if this is true
	irdust++;
	prevdustradius=dustradius[ifile];
	  
	//If desired, calculate Qabs, Qsca, and scattering phase function using Mie theory
	if (generic_oc_flag <= 0) {

	  dmgetoc(lnkfile, dustradius[ifile], costheta, npfunc, &Qwave, &Qabs, &Qsca, &nQwave, &Qpfunc);

          //Interpolate optical constants to lambda microns
	  for(ilambda=0;ilambda<nlambda;ilambda++){
	    for(i=0;i<nQwave;i++) if(Qwave[i] > lambda[ilambda]) break; //find where lambda falls in the Qwave array
	    if(i == 0 || i == nQwave-1) { //use the endpoint values if reached
	      Qabs_value[ilambda]=Qabs[i];
	      Qsca_value[ilambda]=Qsca[i];
	      for(ipfunc=0;ipfunc<npfunc;ipfunc++) {
		//Note: for a 2D array entry c[ix,iy], the 1D index = ix + iy*xsize
		//Here, for pfunc: ix = ilambda, iy = ipfunc, and xsize = nlambda
		// and for Qpfunc: ix = i, iy = ipfunc, and xsize = nQwave
		pfunc[ilambda+nlambda*ipfunc] = Qpfunc[i+nQwave*ipfunc];
	      }
	    }
	    else { //otherwise, interpolate...
	      Qabs_value[ilambda] = Qabs[i-1] + (Qabs[i] - Qabs[i-1])/(Qwave[i] - Qwave[i-1]) * (lambda[ilambda] - Qwave[i-1]);
	      Qsca_value[ilambda] = Qsca[i-1] + (Qsca[i] - Qsca[i-1])/(Qwave[i] - Qwave[i-1]) * (lambda[ilambda] - Qwave[i-1]);
	      for(ipfunc=0;ipfunc<npfunc;ipfunc++) {
		//Note: for a 2D array entry c[ix,iy], the 1D index = ix + iy*xsize
		//Here, for pfunc: ix = ilambda, iy = ipfunc, and xsize = nlambda
		// and for Qpfunc: ix = i, iy = ipfunc, and xsize = nQwave
		pfunc[ilambda+nlambda*ipfunc] = Qpfunc[(i-1)+nQwave*ipfunc] + (Qpfunc[i+nQwave*ipfunc] - Qpfunc[(i-1)+nQwave*ipfunc])/(Qwave[i] - Qwave[i-1]) * (lambda[ilambda] - Qwave[i-1]);
	      }
	    }
	    qabs_out[irdust+ilambda*nrdust] = Qabs_value[ilambda];
	    qsca_out[irdust+ilambda*nrdust] = Qsca_value[ilambda];
	    for(ipfunc=0;ipfunc<npfunc;ipfunc++) pfunc_out[irdust+ilambda*nrdust+ipfunc*nlambda*nrdust] = pfunc[ilambda+nlambda*ipfunc];
	  }
	}


	//If desired, calculate Qabs and Qsca using generic optical constants
	if (generic_oc_flag > 0) {
	  twoPIs = 2.0 * PI * dustradius[ifile];
	  for(ilambda=0;ilambda<nlambda;ilambda++){
	    if(lambda[ilambda] < twoPIs) { 
	      Qabs_value[ilambda] = 1.;
	    }
	    else {
	      Qabs_value[ilambda] = pow(twoPIs / lambda[ilambda], generic_oc_exp);
	    }
	    Qsca_value[ilambda] = (generic_oc_albedo / (1. - generic_oc_albedo)) * Qabs_value[ilambda];
	    qabs_out[irdust+ilambda*nrdust] = Qabs_value[ilambda];
	    qsca_out[irdust+ilambda*nrdust] = Qsca_value[ilambda];
	  }
	  //also need to define nQwave, Qwave, and Qabs for temperature calculation
	  nQwave = 1000;
	  Qwave = (float *) malloc (nQwave * sizeof(float));
	  if (Qwave==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Qwave in generic_oc_flag section.\n"); exit (1);}
	  for(ilambda=0;ilambda<nQwave;ilambda++) Qwave[ilambda] = (float) pow(10.,(ilambda * 5.) / (nQwave-1) - 1.);
	  Qabs = (float *) malloc (nQwave * sizeof(float));
	  if (Qabs==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Qabs in generic_oc_flag section.\n"); exit (1);}
	  for(ilambda=0;ilambda<nQwave;ilambda++) {
	    if(Qwave[ilambda] < twoPIs) { 
	      Qabs[ilambda] = 1.;
	    }
	    else {
	      Qabs[ilambda] = pow(twoPIs / Qwave[ilambda], generic_oc_exp);
	    }
	  }
	}
	
	//If desired, calculate scattering phase function w/ a Henyey-Greenstein function
	if (scatteredlight_flag > 0 && HG_flag > 0) {
	  for(ilambda=0;ilambda<nlambda;ilambda++){
	    //Note: for a 2D array entry c[ix,iy], the 1D index = ix + iy*xsize
	    //Here, for pfunc: ix = ilambda, iy = 0, and xsize = nlambda
	    hgphasefunction(HG_g, costheta, npfunc, nlambda, &(pfunc[ilambda+nlambda*0]));
	    for(ipfunc=0;ipfunc<npfunc;ipfunc++) pfunc_out[irdust+ilambda*nrdust+ipfunc*nlambda*nrdust] = pfunc[ilambda+nlambda*ipfunc];
	  }
	}
	
	//Calculate temperature as a function of circumstellar distance
	//Need to do this for scattered light, too, so that tsublimate works
	if (thermalemission_flag == 1 || scatteredlight_flag == 1) {  
	  for(i=0;i<nr;i++) r_list[i] = minr * pow(10.,delta_log_r * i); //calculate r_list array

	  double *tempBnustar; //create an array to contain stellar spectrum
	  tempBnustar = (double *) malloc (nQwave * sizeof(double));
	  if (tempBnustar==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating tempBnustar array for Qwave.\n"); exit (1);}
	  dmgetBnustar(Tstar, logg, Lstar, kurucz_flag, kurucz_file, Qwave, nQwave, &Rstar, &tempBnustar[0], &kteff, &klogg);	  //calculate Bnustar for all wavelengths in Qwave

	  dmtemperaturecalc(Rstar, Tstar, tempBnustar, r_list, nr, Qwave, Qabs, nQwave, Tvsr_list);	  //calculate temperature
	  free(tempBnustar); //free stellar spectrum
	}

	
      } //end of prevdustradius if statement

      
    } //end of scattered light or thermal emission flag if statement


    //Create the image    
    dm2d(histo, histoxsize, histoysize, data, ndatapoints, fovx, fovy, pixelsize, inclination, longitude, distance, pa, lambda, nlambda, Qabs_value, Qsca_value, dustradius[ifile], Bnustar, Rstar, Tvsr_list, delta_log_r, nr, log10minr, Tsublimate, pfunc, npfunc, Qwave, Qabs, nQwave, iwa, aitoff_flag, densityhisto_flag, opticaldepth_flag, scatteredlight_flag, thermalemission_flag, extinction_flag, xshift, effrdust, distmask);

  }  
  
  if(verbose_flag >= 1) {printf("\rProgress: %3.0f%% complete\n", 100.0); fflush(stdout);}

  //Free memory
  free(data);
  free(Qabs_value);
  free(Qsca_value);
  free(r_list);
  free(Tvsr_list);
  free(Bnustar);
  for (i = 0; i < numfiles; i++) free(inputfiles[i]);
  free(inputfiles);
  if(scatteredlight_flag == 1 || thermalemission_flag == 1) {free(Qwave); free(Qabs); free(Qsca);}
  if(scatteredlight_flag == 1) free(Qpfunc);


  //Mark locations
  if(nmarks > 0) {
    if(verbose_flag >= 1) printf("Marked locations.\n");
    markdata = (datatype3 *) malloc (nmarks * sizeof(datatype3));
    if(verbose_flag >= 1) printf("Marked locations.\n");
    if (markdata==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating markdata array.\n"); exit (1);}
    if(verbose_flag >= 1) printf("Marked locations.\n");
    for(i=0;i<nmarks;i++){
      markdata[i].x = markx[i];
      markdata[i].y = marky[i];
      markdata[i].z = markz[i];
      markdata[i].intensity = markweight[i];
    }
    if(verbose_flag >= 1) printf("Marked locations.\n");
    dmmarklocs(histo, markdata, nmarks, fovx, fovy, pixelsize, inclination, longitude, distance, pa, aitoff_flag, markabs_flag, xshift);
    if(verbose_flag >= 1) printf("Marked locations.\n");
    free(markdata);
  }

  return 0;
}
