void getKurucz(char *kurucz_file, float Tstar, float logg, int *nkl, float **kl, float **knu, float **kdnu, double **kBnu, double **kBnucont, float *kteff, float *klogg){

  //This routine retrieves Kurucz model stellar atmosphere data
  
  //INPUTS
  //Tstar: stellar effective temperture (K)
  //logg: surface gravity

  //OUTPUTS
  //nkl: a integer specifying the length of the Kurucz vectors
  //kl: a vector of wavelengths in microns (Kurucz gives in nm, but we convert to microns)
  //knu: a vector of frequencies in Hz
  //kdnu: a vector of Delta_frequency in Hz
  //kBnu: surface brightness **SEE NOTE BELOW
  //kBnucont: surface brightness of the continuum **SEE NOTE BELOW

  //IMPORTANT NOTE: The Kurucz flux data is in units of erg s^-1 cm^-2 Hz^-1 steradian^-1.  This
  //is the same units as surface brightness, so we might think that Kurucz is providing Bnu, where
  //Fnu = pi * Bnu * (Rstar / d)^2 is the flux Fnu for a star of radius Rstar at distance d.  But in
  //fact, for whatever reason, Kurucz is providing the quantity (pi * Bnu) / (4 * pi).  This is the
  //surface flux per unit solid angle--this must be important in stellar atmospheres or something.
  //Anyway, to calculate Bnu, which is what this routine returns, we multiply Kurucz's numbers by 4.
  //If this sounds crazy to you, I note that this agrees with the SYNPHOT Data User's Guide, which
  //can be found at http://www.stsci.edu/hst/HST_overview/documents/synphot/AppA_Catalogs9.html

  int i;
  FILE *pFile;
  char header[500];
  char tempstr[20];
  char *pch;
  float lightspeed = 2.99792458e14; //microns s^-1

  *nkl = 1221; //number of wavelengths in Kurucz model files
  *kteff=0.0;
  *klogg=0.0;

  float *tmp;
  double *dtmp;
  tmp = (float *) realloc (*kl, *nkl * sizeof(float));
  *kl = tmp;
  if (*kl==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating kl array.\n"); exit (1);}
  tmp = (float *) realloc (*knu, *nkl * sizeof(float));
  *knu = tmp;
  if (*knu==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating knu array.\n"); exit (1);}
  tmp = (float *) realloc (*kdnu, *nkl * sizeof(float));
  *kdnu = tmp;
  if (*kdnu==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating kdnu array.\n"); exit (1);}
  dtmp = (double *) realloc (*kBnu, *nkl * sizeof(double));
  *kBnu = dtmp;
  if (*kBnu==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating kBnu array.\n"); exit (1);}
  dtmp = (double *) realloc (*kBnucont, *nkl * sizeof(double));
  *kBnucont = dtmp;
  if (*kBnucont==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating kBnucont array.\n"); exit (1);}



  // Open the Kurucz file and determine number of wavelength values
  pFile = fopen (kurucz_file, "r");
  if(pFile == NULL) {fprintf(stderr,"ERROR: cannot open file %s\n",kurucz_file); exit (1);}
  for(i=0;i<22;i++) {
    fgets(header,500,pFile); //First, read the header (must be 2 lines)
  }

  //Read in wavelength data
  for(i=0;i<*nkl;i++) {
    fscanf(pFile, "%f",&(*kl)[i]);
    (*kl)[i] *= 0.001; //convert from nm to microns
  }

  while(rint(*kteff) < rint(Tstar) || fabs(*klogg-logg) > 0.25) {
    while(strstr(header,"GRAVITY") == NULL) fgets(header,500,pFile);
    pch = strtok(header," "); //TEFF
    pch = strtok(NULL," "); //the value of TEFF
    *kteff = atof(pch);
    pch = strtok(NULL," "); //GRAVITY
    pch = strtok(NULL," "); //the value of GRAVITY
    *klogg = atof(pch);
    if (*kteff > 49999. && *klogg > 4.9) {fprintf(stderr,"ERROR: Tstar and logg beyond the scope of Kurucz models.  See http://wwwuser.oats.inaf.it/castelli/grids/gridp00k0odfnew/fp00k0tab.html for a grid of acceptable values.\n"); exit (1);}
  }

  for(i=0;i<*nkl;i++) {
    fgets(tempstr,11,pFile);
    if(strstr(tempstr,"\n") == NULL){
      (*kBnu)[i] = atof(tempstr);
      (*kBnu)[i] *= 4.0; //convert Kurucz's surface flux per unit solid angle (pi*Bnu)/(4*pi) to Bnu (see above note) 
    } else i--;
  }
  for(i=0;i<*nkl;i++) {
    fgets(tempstr,11,pFile);
    if(strstr(tempstr,"\n") == NULL){
      (*kBnucont)[i] = atof(tempstr);
      (*kBnucont)[i] *= 4.0; //convert Kurucz's surface flux per unit solid angle (pi*Bnu)/(4*pi) to Bnu (see above note) 
    } else i--;
  }

  fclose(pFile);


  //Calculate nu array
  for(i=0;i<*nkl;i++) (*knu)[i] = lightspeed / (*kl)[i];

  //Calculate dnu array
  for(i=0;i<*nkl-1;i++) (*kdnu)[i] = (*knu)[i+1] - (*knu)[i];
  (*kdnu)[*nkl-1] = (*kdnu)[*nkl-2];

}
