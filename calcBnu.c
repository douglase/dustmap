void calcBnu(float temperature, float lambda, double *result){

  // calcBnu.c
  // Calculates the Planck Law per unit frequency in CGS units
  //
  // Input:
  // temperature = temperature in Kelvin
  // lambda = wavelength in microns
  //
  // Output:
  // result = Bnu in erg s^-1 Hz^-1 cm^-2 steradian^-1
  //
  // Created by Christopher Stark
  // Carnegie Institution of Washington
  // cstark@dtm.ciw.edu
  // Last updated 10 Sep 2010 by Christopher Stark

  double const1, const2;
  const1 = (double) 3.972895e-4;  //const1 = 2 * h * c * (1E6)^3 (in CGS w/ correction for microns)
  const2 = (double) 1.438769e4;   //const2 = h * c * 1E6 / k (in CGS w/ correction for microns)

  *result = const1 / (lambda*lambda*lambda * (exp(const2 / (lambda * temperature)) - 1.0));

}
