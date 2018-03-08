typedef struct {
  int ID;
  float x;
  float y;
  float z;
} datatype1;

typedef struct {
  int ID;
  float x;
  float y;
  float z;
  float vx;
  float vy;
  float vz;
} datatype2;

typedef struct {
  int ID;
  float x;
  float y;
  float z;
  double intensity;
} datatype3;

typedef struct {
  int ID;
  float x;
  float y;
  float z;
  float vx;
  float vy;
  float vz;
  double intensity;
} datatype4;


int dustmap_c(int argc, void *argv[]);

//cython demo function:
void cos_doubles(double * in_array,
   double * out_array,
   int size);


void dustmap_call(double * hist,// = (double *) argv[j]; j++;     //The final image, passed by reference
  int histoxsize,// = *(int *) argv[j]; j++;
  int histoysize,// = *(int *) argv[j]; j++;
  char * inputfilelist ,//= (char *) argv[j]; j++; //Concatenated strings of input files from IDL
  int numfiles,// = *(int *) argv[j]; j++;        //number of input file strings to expect
  int datatype,// = *(int *) argv[j]; j++;        //The data type (see dmgetdata.c)
  float distance,// = *(float *) argv[j]; j++;     //Distance to the star in pc
  float fovx,// = *(float *) argv[j]; j++;         //Field of view in the x-direction in mas
  float fovy,// = *(float *) argv[j]; j++;         //Field of view in the y-direction in mas
  float pixelsize,// = *(float *) argv[j]; j++;    //Pixel size in mas
  float inclination,// = *(float *) argv[j]; j++;  //Inclination of the system
  float longitude,// = *(float *) argv[j]; j++;    //Orbital phase of the system
  float pa,// = *(float *) argv[j]; j++;           //Position angle of the system in degrees
  float * lambda,// = (float *) argv[j]; j++;       //Wavelength of image in microns
  int nlambda,// = *(int *) argv[j]; j++;
  float Lstar,// = *(float *) argv[j]; j++;        //Luminosity of star in solar luminosities
  float Tstar,// = *(float *) argv[j]; j++;        //Effective temperature of star in degrees K
  float logg,// = *(float *) argv[j]; j++;         //Surface gravity of star
  int kurucz_flag,// = *(int *) argv[j]; j++;      //Flag signaling a Kurucz stellar atmosphere model is desired
  char *kurucz_file,// = (char *) argv[j]; j++;
  float *fstar,// = (float *) argv[j]; j++;       //Output stellar flux in Jy
  float *dustradius,// = (float *) argv[j]; j++;   //Radius of dust grain in microns, array with length equal to numfiles
  float Tsublimate,// = *(float *) argv[j]; j++;   //Sublimation temperature of the dust in degrees K
  char *lnkfile,// = (char *) argv[j]; j++;        //File name containing indeces of refraction vs. lambda for dust
  double *scaling,// = (double *) argv[j]; j++;    //Factor to scale up the amount of dust
  float iwa,// = *(float *) argv[j]; j++;          //Radius of an occulting spot to place on image
  int aitoff_flag,// = *(int *) argv[j]; j++;
  int densityhisto_flag,// = *(int *) argv[j]; j++;    //Signals a density histo is to be made, ignores image flags below
  int opticaldepth_flag,// = *(int *) argv[j]; j++;    //Signals a geometric optical depth histo is to be made
  int scatteredlight_flag,// = *(int *) argv[j]; j++;  //Signals image to include scattered light calculation
  int thermalemission_flag,// = *(int *) argv[j]; j++; //Signals image to include thermal emission calculation
  int verbose_flag,// = *(int *) argv[j]; j++;         //Print a lot of info
  int generic_oc_flag,// = *(int *) argv[j]; j++;         //Print a lot of info
  float generic_oc_exp,// = *(float *) argv[j]; j++;   //An input exponent for the Qabs/Qsca power law; non-zero value overides values from lnkfile
  float generic_oc_albedo,// = *(float *) argv[j]; j++;           //An input albedo to adjust generic Qsca
  int HG_flag,// = *(int *) argv[j]; j++;         //Print a lot of info

  float HG_g,// = *(float *) argv[j]; j++;     //An input g value; non-zero value overides values from lnkfile
  int extinction_flag,// = *(int *) argv[j]; j++;      //Signals that extinction of the star light should be taken into account
  float xshift,// = *(float *) argv[j]; j++;           //A left-right shift in the viewing direction (AU)
  float effrdust,// = *(float *) argv[j]; j++;         //The effective radius of the dust grains for gaussian image smoothing (AU)
  float *markx,// = (float *) argv[j]; j++;      //array of x coordinates of points to mark
  float *marky,// = (float *) argv[j]; j++;      //array of y coordinates of points to mark
  float *markz,// = (float *) argv[j]; j++;      //array of z coordinates of points to mark
  float *markweight,// = (float *) argv[j]; j++; //array of weights of the marks
  int nmarks,// = *(int *) argv[j]; j++;        //Number of entries in above mark arrays
  int markabs_flag,// = *(int *) argv[j]; j++;  //Whether weights are absolute (default is relative)
  int az_sym,// = *(int *) argv[j]; j++;       //A flag to force the disk to have azimuthal symmetry, value equals
                                            //the number of rotation angles w/in 360 deg. to average over
  float distmask,// = *(float *) argv[j]; j++; //dust within this radius from the observer is deleted
  int npfunc,// = *(int *) argv[j]; j++; //number of cos(theta) values
  float *pfunc,// = (float *) argv[j]; j++; //used internally for the phase function
  float *costheta,// = (float *) argv[j]; j++; //output costheta
  float *pfunc_out,// = (float *) argv[j]; j++; //output phase function
  float *qabs_out,// = (float *) argv[j]; j++; //output qabs
  float *qsca_out,// = (float *) argv[j]; j++; //output qsca
  int nrdust// = *(int *) argv[j]; j++; //number of unique dust sizes
);
//https://cython.readthedocs.io/en/latest/src/userguide/external_C_code.html#c-api-declarations
