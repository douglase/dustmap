
void dmgetBnustar(float Tstar, float logg, float Lstar, int kurucz_flag, char *kurucz_file, float *lambda, int nlambda, float *Rstar, double *Bnustar, float *kteff, float *klogg){

  // dmgetBnustar.c
  // Returns surface brightness and radius of a star of temperature Tstar, surface gravity logg,
  // and luminosity Lstar. Can use blackbody function or Kurucz model.
  // 
  // Input:
  // Tstar: stellar effective temperature (K)
  // Lstar: luminosity (L_sun)
  // lambda: wavelength vector at which to calculate Bnu
  // nlambda: the length of the lambda vector
  // kurucz_flag: if non-zero, then we load the appropriate kurucz stellar model
  //
  // Output:
  // Rstar: stellar radius required to achieve a luminosity of Lstar (solar radii)
  // Bnustar: the surface brightness as a function of lambda (in erg s^-1 Hz^-1 cm^-2 steradian^-1)

  // The flux per unit frequency from a star at distance D can be calculated from the above two
  // quantities using the expression Fnu = pi * Bnustar * (Rstar / D)^2
  // 
  //
  // Created by Christopher Stark
  // NASA Goddard
  // christopher.stark@nasa.gov
  // Last updated 14 May 2014 by Christopher Stark

  
  int i;
  double tempBnu;

  if (kurucz_flag == 0) { //Blackbody
    
    for(i=0;i<nlambda;i++) {
      calcBnu(Tstar, lambda[i], &tempBnu);
      Bnustar[i] = tempBnu;
    }
    *Rstar = pow(5778./Tstar,2.0) * sqrt(Lstar); //(in Rsun) comes from ratio of Stefan-Boltzmann Laws
  }
  else { //Kurucz model

    int ilambda;
    double Rstarcm;
    double LboloR2 = 0.0;
    float kB = 1.38065e-16;
    float h = 6.6261e-27;
    float klolambda, linterp, uinterp, w;

    // = 1.3806488e−16; //Boltzmann constant in CGS (erg K^-1)
    //float h = 6.6261e-27; //Planck constant in CGS (cm^2 g s^-1)

    int nkl;
    float *kl=NULL, *knu=NULL, *kdnu=NULL;
    double *kBnu=NULL, *kBnucont=NULL;

    //Retrieve the closest Kurucz model
    getKurucz(kurucz_file, Tstar, logg, &nkl, &kl, &knu, &kdnu, &kBnu, &kBnucont, kteff, klogg);
    
    //Calculate Bolometric luminosity
    for(i=0;i<nkl-1;i++) LboloR2 += 0.5 * (kBnu[i]+kBnu[i+1]) * fabs(kdnu[i]);
    LboloR2 *= (4.0 * PI * PI);
    
    //Now calculate Rstar for the desired bolometric luminosity
    Rstarcm = sqrt((((double) Lstar)*3.846e33) / LboloR2); //the 3.846e33 converts 1 Lsun to erg/s, Rstar is in cm
    *Rstar = (float) (Rstarcm / (6.955e10)); //convert to solar radii

    //Finally, interpolate the surface brightness to the desired wavelengths
    for(ilambda=0;ilambda<nlambda;ilambda++){


      for(i=0;i<nkl;i++) if(kl[i] > lambda[ilambda]) break; //find where lambda falls in the kl array
      if(i == 0) Bnustar[ilambda] = 0.0; //too far in the UV
      else { //otherwise, interpolate...
	if(h*knu[i-1]/(kB*Tstar) < 0.2) { //If this is true, we interpolate in the Rayleigh-Jean regime
	  klolambda = kl[i-1]/lambda[ilambda];
	  linterp = kBnu[i-1]*klolambda*klolambda; //Interpolation from the next shortest grided wavelength
	  if(i==nkl) {
	    Bnustar[ilambda] = linterp; //If we're beyond the longest gridded wavelength, just interpolate with lambda^-2 from there
	  } else {
	    klolambda = kl[i]/lambda[ilambda]; //If we're between two data points, trust the data and smoothly interpolate between them
	    uinterp = kBnu[i]*klolambda*klolambda;
	    w = (lambda[ilambda]-kl[i-1]) / (kl[i]-kl[i-1]);
	    Bnustar[ilambda] = (1.-w) * linterp + w * uinterp;
	  }
	}
	else { //If not in R-J regime, do a linear interpolation in log-log space b/c the spectrum could be really complex
	  Bnustar[ilambda] = log10(kBnu[i-1]) + (log10(kBnu[i]) - log10(kBnu[i-1]))/(log10(kl[i]) - log10(kl[i-1])) * (log10(lambda[ilambda]) - log10(kl[i-1]));
	  Bnustar[ilambda] = pow(10.0,Bnustar[ilambda]);
	}
      }


    }

    free(kl);
    free(knu);
    free(kdnu);
    free(kBnu);
    free(kBnucont);
  }


}
