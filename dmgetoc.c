#include "mie_single.c"

void dmgetoc(char *lnkfile, float s, float *costheta, int numthetas, float **lambda, float **Qabs, float **Qsca, int *nlambda, float **Qpfunc){

  // dmgetoc.c
  // Returns the optical constants for dust grains of a given size
  // and composition as a function of wavelength
  // 
  // Input:
  // lnkfile: name of the lnk file to be opened
  //      - this is a ascii file, with each line containing l n k,
  //        where l is wavelength in microns, n and k are indices of refraction
  // s: dust grain size in microns
  // costheta: vector of cos(scattering angle) values
  // numthetas: the length of the costheta vector
  //
  // Output:
  // lambda: an array containing the wavelengths in microns
  // Qabs: the absorption efficiencies corresponding to each lambda
  // Qsca: the scattering efficiencies corresponding to each lambda
  // nlambda: the number of elements in the lambda array (as well as Qabs and Qsca arrays)
  // Qpfunc: the scattering phase function array for all wavlengths & thetas
  // 
  //
  // Created by Christopher Stark
  // Carnegie Institution of Washington
  // cstark@dtm.ciw.edu
  // Last updated 12 Sep 2010 by Christopher Stark

  
  FILE *pFile;
  char header[500];
  int i, ilambda;
  float *n, *k;
  float a, b, c;
  float twoPIs = 2.0 * PI * s;

  // Open the lnk file and determine number of lambda values
  pFile = fopen (lnkfile, "r");
  if(pFile == NULL) {fprintf(stderr,"ERROR: cannot open file %s\n",lnkfile); exit (1);}
  fgets(header,500,pFile); //First, read the header (must be 2 lines)
  fgets(header,500,pFile);
  *nlambda=0;
  while(!feof(pFile)){ //Now calculate number of lambda values
    fscanf(pFile, "%f %f %f\n",&a,&b,&c);
    (*nlambda)++;
  }

  //Now make the arrays the appropriate size
  float *tmp;
  tmp = (float *) realloc (*lambda, *nlambda * sizeof(float));
  *lambda = tmp;
  if (*lambda==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating lambda array.\n"); exit (1);}
  n = (float *) malloc (*nlambda * sizeof(float));
  if (n==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating n array.\n"); exit (1);}
  k = (float *) malloc (*nlambda * sizeof(float));
  if (k==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating k array.\n"); exit (1);}
  tmp = (float *) realloc (*Qabs, *nlambda * sizeof(float));
  *Qabs = tmp;
  if (*Qabs==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Qabs array.\n"); exit (1);}
  tmp = (float *) realloc (*Qsca, *nlambda * sizeof(float));
  *Qsca = tmp;
  if (*Qsca==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Qsca array.\n"); exit (1);}

  tmp = (float *) realloc (*Qpfunc, *nlambda * numthetas * sizeof(float));
  *Qpfunc = tmp;
  if (*Qpfunc==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating Qpfunc array.\n"); exit (1);}


  //Go back to the beginning and read in data
  fseek (pFile, 0, SEEK_SET);
  fgets(header,500,pFile); //First, read the header (must be 2 lines)
  fgets(header,500,pFile);
  for(i=0;i<*nlambda;i++) fscanf(pFile, "%f %f %f\n",&(*lambda)[i],&n[i],&k[i]);
  fclose(pFile);


  //Now we use lambda, n, and k to calculate Qabs, Qsca, and Qpfunc with some Mie code
  mie_single(s, *lambda, n, k, *nlambda, costheta, numthetas, *Qabs, *Qsca, *Qpfunc); //Use Mie theory

  free(n);
  free(k);


}
