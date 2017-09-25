void hgphasefunction(float g, float *costheta, int numthetas, int nlambda, float *pfunc){

  // hgphasefunction.c
  // Calculates the Henyey-Greenstein phase function
  //
  // Input:
  // g: the HG phase function asymmetry parameter
  // costheta: an array of costheta values, the cosine
  //           of the scattering angle
  // numthetas: the number of entries in costheta array
  //
  // Output:
  // pfunc: an array of scattering phase function amplitudes,
  //        with a 1:1 correspoding to the costheta values
  //
  // Created by Christopher Stark
  // NASA GSFC
  // christopher.c.stark@nasa.gov
  // Last updated 18 Nov 2014 by Christopher Stark

  int i;
  float tempval;
  float N;
 
  //calculate phase function as a function of costheta
  for(i=0;i<numthetas;i++){
    tempval = (1.0 - 2.0 * g * costheta[i] + g * g);
    pfunc[i*nlambda] = (1.0 - g * g) / sqrt(tempval*tempval*tempval) / (4 * PI);
  }

}
