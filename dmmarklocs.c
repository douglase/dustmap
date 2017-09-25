void dmmarklocs(double *histo, datatype3 *data, long ndatapoints, float fovx, float fovy, float pixelsize, float inclination, float longitude, float distance, float pa, int aitoff_flag, int weightabs_flag, float xshift){

  // dmmarklocs.c
  // Adds weight to pixels at a specific location
  // 
  // This routine is very similar to dm2d.c.
  // See dm2d.c for more detail.
  //
  // Created by Christopher Stark
  // Carnegie Institution of Washington
  // cstark@dtm.ciw.edu
  // Last updated 26 Oct 2010 by Christopher Stark

  


  long i;
  long histoxsize, histoysize, histosize, histoxsizeo2, histoysizeo2;

  float oneopixelsize;
  float fovxo2, fovyo2;
  float x, y, z;
  float pio180 = PI / 180.0;
  float rotmatrix[3][3], irotmatrix[3][3], parotmatrix[3][3];
  float sqrt2, oneof, alpha2, delta, cdec, oneodenom, tempxang, tempyang; //used if aitoff projection requested

  double maxval;

  float *d;
  d = (float *) malloc (ndatapoints * sizeof(float));
  if (d==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating d array.\n"); exit (1);}

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


  histoxsize = (long) 2 * floor(fovx / pixelsize * 0.5);
  histoysize = (long) 2 * floor(fovy / pixelsize * 0.5);
  histosize = histoxsize * histoysize;

  //Convert angles to radians
  inclination *= pio180;
  longitude *= -pio180; //longitude is a CW angle we want to rotate points CCW
  fovx *= pio180;
  fovy *= pio180;
  pixelsize *= pio180;
  //pa += 90.0;
  pa *= -pio180; //pa is a CW angle and we want to rotate points CCW
  
  //Here are some variables defined so that computations are performed more quickly later on
  oneopixelsize = 1. / pixelsize;
  histoxsizeo2 = histoxsize / 2;
  histoysizeo2 = histoysize / 2;

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
    xang[i] = atan2(data[i].x, d[i]); //longitude
    yang[i] = atan2(data[i].y, sqrt(data[i].x*data[i].x + d[i]*d[i])); //latitude
  }


  
  //Make histo indeces
  for(i=0;i<ndatapoints;i++){
    xindex[i] = floor(floor(xang[i] * oneopixelsize) + histoxsizeo2);
    yindex[i] = floor(floor(yang[i] * oneopixelsize) + histoysizeo2);
    if( (yindex[i] < histoysize) && (xindex[i] < histoxsize) && (yindex[i] >= 0) && (xindex[i] >= 0) ) {
      inlim[i] = 1; //mark as "in bounds"
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
  
  //Find maximum value of histo
  maxval = 0.0;
  for(i=0;i<histosize;i++) if(histo[i] > maxval) maxval = histo[i];

  //Mark those points!
  for(i=0;i<ndatapoints;i++) {
    if(inlim[i]) {
      if(weightabs_flag == 0) {
	histo[xindex[i] + yindex[i] * histoxsize] = maxval * data[i].intensity;
      }
      else {
	histo[xindex[i] + yindex[i] * histoxsize] = data[i].intensity;
      }
    }
  } 



  free(d); free(xang); free(yang); free(xindex); free(yindex); free(inlim);

}
