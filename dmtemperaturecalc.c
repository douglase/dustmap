void dmtemperaturecalc(float rstar, float tstar, double *bnustar, float *r, int nr, float *lambda, float *Qabs, int nlambda, float *temperature){

  // dmtemperaturecalc.c
  // Iteratively calculates the equilibrium temperature of a particle
  // at a given circumstellar distance using optical constants.
  // 
  // Input:
  // rstar: radius of the star in stellar radii
  // bnustar: the stellar surface brightness (either Planck law or Kurucz surf. brightness)
  // r: array of circumstellar distances in AU
  // nr: the nummber of elements in the r array
  // lambda: an array of wavelengths in microns
  // Qabs: an array of Qabs values for the dust grain as a function of lambda
  // nlambda: the number of elements in the lambda (and Qabs) array
  //
  // Output:
  // temperature: an array containing the equilibrium temperatures of the dust grain
  //
  // Created by Joannah Metz for IDL
  // Translated to C & updated by Christopher Stark
  // Carnegie Institution of Washington
  // cstark@dtm.ciw.edu
  // Last updated 20 Sep 2010 by Christopher Stark


  float tol = 0.001; // the fractional energy error that is acceptable
                     // -- this controls the precision of temperature 

  //*********************
  // Fundamental constants in CGS units
  double c = 2.9979e10; //in cm s^-1
  //*********************


  //*********************
  // Variable declarations

  int i, j;
  int ok, bad, oldbad;

  float tgrain, tstep;
  float *deltanu;

  deltanu = (float *) malloc ((nlambda-1) * sizeof(float));
  if (deltanu==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating deltanu array.\n"); exit (1);}

  double Eabs, Eemit, Eemit_ll, Eemit_ul;
  double rstarcm, rcm;
  double *bnu;
  bnu = (double *) malloc (nlambda * sizeof(double));
  if (bnu==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating bnu array.\n"); exit (1);}
  //********************



  


  for(i=0;i<nr;i++) if(r[i] <= 0) {fprintf(stderr,"ERROR: r must be non-zero.\n"); exit (1);}

  rstarcm = (double) rstar * 6.9550e10;

  // Make deltanu array
  for(i=0;i<nlambda-1;i++) deltanu[i] = fabs((c / (lambda[i+1] * 1e-4)) - (c / (lambda[i] * 1e-4)));


  // Loop through all r values
  for(j=0;j<nr;j++) {

    rcm = (double) r[j] * 1.495979e13; //distance in cm
    if(j==0) tgrain = tstar * sqrt(0.5 * rstarcm / rcm); //guess the local BB temp for a starting value
    else tgrain = temperature[j-1]; //every other starting value gets the temp at the neighboring r value



    // Calculate absorbed energy
    // Note: the 0.25 below is because two midpoints are taken during the sum, so
    // the values should be divided by 2 twice, but it's faster to do it at the end
    Eabs = 0.0;
    for(i=0;i<nlambda-1;i++) Eabs += (bnustar[i]+bnustar[i+1]) * deltanu[i] * (Qabs[i]+Qabs[i+1]);
    Eabs *= (rstarcm / rcm) * (rstarcm / rcm) * 0.25; //Note: this isn't really the energy absorbed -- we would need to 
                                                      //multiply by pi*pi*s^2 to get the right number.  But we're
                                                      //only interested in comparing it to Eemit, which is also
                                                      //off by the same factor -- done to save a little computation time.
    Eemit_ul = Eabs * (1.0 + tol); //upper limit given the tolerance, used later on
    Eemit_ll = Eabs * (1.0 - tol); //lower limit given the tolerance, used later on

    ok = 0;
    bad = 0;
    tstep = 20.0;
    
    while(ok == 0){

      // Calculate emitted energy
      for(i=0;i<nlambda;i++) calcBnu(tgrain,lambda[i],&bnu[i]); // bnu (erg s^-1 cm^-2 ster^-1 Hz^-1)
      Eemit = 0.0;
      for(i=0;i<nlambda-1;i++) Eemit += (bnu[i]+bnu[i+1]) * deltanu[i] * (Qabs[i]+Qabs[i+1]); //Note: this isn't really the
                                                                                              //energy emitted -- see above. 
      //Eemit *= 4.0 * 0.25;  //Normally we would have to multiply by 4, but b/c the midpoint method was used above,
                              //we also have to divide by 4.  So let's just not multiply by unity.

      
      oldbad = bad;
      
      if(Eemit > Eemit_ul) { //If emitted energy is too high, decrease the temperature
	bad = 1;
	tgrain -= tstep;
      }
      else if(Eemit < Eemit_ll) {  //If emitted energy is too low, increase the temperature
	bad = -1;
	tgrain += tstep;
      }
      else if(Eemit < Eemit_ul && Eemit > Eemit_ll) ok = 1; // If it's in the sweet spot, we've won!

      if(oldbad != bad) tstep *= 0.5; //If we've changed directions, reduce the step size so we don't keep overshooting it
      
    }

    temperature[j] = tgrain;
  }

  free(bnu);
  free(deltanu);

}
