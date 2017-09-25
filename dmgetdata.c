void dmgetdata(char *inputfilename, int datatype, double scaling, datatype3 **data, long *nentries, int az_sym){

  // Reads binary data from an input file with one of the following
  // four possible datatypes:

  // Type 1 -- Positions only
  //        - each entry is of the form: ID x y z
  //            - dust particle ID # (4-byte int)
  //            - x position (4-byte float)
  //            - y position (4-byte float)
  //            - z position (4-byte float)
  // Type 2 -- Positions & velocities
  //        - each entry is of the form: ID x y z vx vy vz
  //            - dust particle ID # (4-byte int)
  //            - x position (4-byte float)
  //            - y position (4-byte float)
  //            - z position (4-byte float)
  //            - vx velocity (4-byte float)
  //            - vy velocity (4-byte float)
  //            - vz velocity (4-byte float)
  // Type 3 -- Collision data, positions only
  //        - each entry is of the form: ID x y z intensity
  //            - dust particle ID # (4-byte int)
  //            - x position (4-byte float)
  //            - y position (4-byte float)
  //            - z position (4-byte float)
  //            - intensity (8-byte float)
  // Type 4 -- Collision data, positions & velocities
  //        - each entry is of the form: ID x y z vx vy vz
  //            - dust particle ID # (4-byte int)
  //            - x position (4-byte float)
  //            - y position (4-byte float)
  //            - z position (4-byte float)
  //            - vx velocity (4-byte float)
  //            - vy velocity (4-byte float)
  //            - vz velocity (4-byte float)
  //            - intensity (8-byte float)

  // Created by Christopher Stark
  // Carnegie Institution of Washington
  // cstark@dtm.ciw.edu
  // Last updated 7 Jun 2011 by Christopher Stark


  //------------ DEFINITIONS ---------------

  //General definitions
  FILE *pFile;
  long lSize;
  long i;
  long nentries_in_file;
  int bytesperentry, result;
  float rot;


  //Some calculations to make more definitions...
  //First, calculate how many entries are in the file
  switch(datatype)
    {
    case 1:
      bytesperentry = sizeof(datatype1);
      break;
    case 2:
      bytesperentry = sizeof(datatype2);
      break;
    case 3:
      bytesperentry = sizeof(datatype3);
      break;
    case 4:
      bytesperentry = sizeof(datatype4);
      break;
  }
  pFile = fopen ( inputfilename , "rb" );
  if (pFile==NULL) {fprintf(stderr,"ERROR: cannot open file %s\n",inputfilename); exit (1);}
  fseek (pFile , 0 , SEEK_END); //obtain file size
  lSize = ftell (pFile);
  rewind (pFile);
  nentries_in_file = lSize / bytesperentry;

  if(az_sym > 1) {
    *nentries = nentries_in_file * az_sym;
    rot = 2 * PI / az_sym; //rotational angle size for azimuthal averaging
  } else {
    *nentries = nentries_in_file;
  } 

  //Now define the permanent data structure to hold all of the entries
  datatype3 *tmp;
  tmp = (datatype3 *) realloc (*data,*nentries * sizeof(datatype3));
  *data = tmp;
  if (*data==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating data array.\n"); exit (1);}

  //Also define a temporary data structure to make file reading faster
  datatype1 *tempdata1;
  datatype2 *tempdata2;
  datatype3 *tempdata3;
  datatype4 *tempdata4;

  //--------------------------------------






  //------------ MAIN CODE ---------------

  //Read the data!
  switch(datatype)
    {
    case 1:
      tempdata1 = (datatype1 *) malloc (nentries_in_file * sizeof(datatype1));
      if (tempdata1==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating tempdata array.\n"); exit (1);}
      result = fread (tempdata1,1,lSize,pFile);
      if(az_sym > 1) {
	for(i=0;i<*nentries;i++) {
	  (*data)[i].ID=tempdata1[i % nentries_in_file].ID;
	  (*data)[i].x=tempdata1[i % nentries_in_file].x * cos((i / nentries_in_file) * rot) - tempdata1[i % nentries_in_file].y * sin((i / nentries_in_file) * rot);
	  (*data)[i].y=tempdata1[i % nentries_in_file].x * sin((i / nentries_in_file) * rot) + tempdata1[i % nentries_in_file].y * cos((i / nentries_in_file) * rot);
	  (*data)[i].z=tempdata1[i % nentries_in_file].z;
	  (*data)[i].intensity=1.0;
	}
      } else {
	for(i=0;i<*nentries;i++) {
	  (*data)[i].ID=tempdata1[i].ID;
	  (*data)[i].x=tempdata1[i].x;
	  (*data)[i].y=tempdata1[i].y;
	  (*data)[i].z=tempdata1[i].z;
	  (*data)[i].intensity=1.0;
	}
      }
      free(tempdata1);
      break;
    case 2:
      tempdata2 = (datatype2 *) malloc (nentries_in_file * sizeof(datatype2));
      if (tempdata2==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating tempdata array.\n"); exit (1);}
      result = fread (tempdata2,1,lSize,pFile);
      if(az_sym > 1) {
	for(i=0;i<*nentries;i++) {
	  (*data)[i].ID=tempdata2[i % nentries_in_file].ID;
	  (*data)[i].x=tempdata2[i % nentries_in_file].x * cos((i / nentries_in_file) * rot) - tempdata2[i % nentries_in_file].y * sin((i / nentries_in_file) * rot);
	  (*data)[i].y=tempdata2[i % nentries_in_file].x * sin((i / nentries_in_file) * rot) + tempdata2[i % nentries_in_file].y * cos((i / nentries_in_file) * rot);
	  (*data)[i].z=tempdata2[i % nentries_in_file].z;
	  (*data)[i].intensity=1.0;
	}
      } else {
	for(i=0;i<*nentries;i++) {
	  (*data)[i].ID=tempdata2[i].ID;
	  (*data)[i].x=tempdata2[i].x;
	  (*data)[i].y=tempdata2[i].y;
	  (*data)[i].z=tempdata2[i].z;
	  (*data)[i].intensity=1.0;
	}
      }
      free(tempdata2);
      break;
    case 3:
      tempdata3 = (datatype3 *) malloc (nentries_in_file * sizeof(datatype3));
      if (tempdata3==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating tempdata array.\n"); exit (1);}
      result = fread (tempdata3,1,lSize,pFile);
      if(az_sym > 1) {
	for(i=0;i<*nentries;i++) {
	  (*data)[i].ID=tempdata3[i % nentries_in_file].ID;
	  (*data)[i].x=tempdata3[i % nentries_in_file].x * cos((i / nentries_in_file) * rot) - tempdata3[i % nentries_in_file].y * sin((i / nentries_in_file) * rot);
	  (*data)[i].y=tempdata3[i % nentries_in_file].x * sin((i / nentries_in_file) * rot) + tempdata3[i % nentries_in_file].y * cos((i / nentries_in_file) * rot);
	  (*data)[i].z=tempdata3[i % nentries_in_file].z;
	  (*data)[i].intensity=tempdata3[i % nentries_in_file].intensity;
	}
      } else {
	for(i=0;i<*nentries;i++) {
	  (*data)[i].ID=tempdata3[i].ID;
	  (*data)[i].x=tempdata3[i].x;
	  (*data)[i].y=tempdata3[i].y;
	  (*data)[i].z=tempdata3[i].z;
	  (*data)[i].intensity=tempdata3[i].intensity;
	}
      }
      free(tempdata3);
      break;
    case 4:
      tempdata4 = (datatype4 *) malloc (nentries_in_file * sizeof(datatype4));
      if (tempdata4==NULL) {fprintf(stderr,"ERROR: memory acquisition failed when creating tempdata array.\n"); exit (1);}
      result = fread (tempdata4,1,lSize,pFile);
      if(az_sym > 1) {
	for(i=0;i<*nentries;i++) {
	  (*data)[i].ID=tempdata4[i % nentries_in_file].ID;
	  (*data)[i].x=tempdata4[i % nentries_in_file].x * cos((i / nentries_in_file) * rot) - tempdata4[i % nentries_in_file].y * sin((i / nentries_in_file) * rot);
	  (*data)[i].y=tempdata4[i % nentries_in_file].x * sin((i / nentries_in_file) * rot) + tempdata4[i % nentries_in_file].y * cos((i / nentries_in_file) * rot);
	  (*data)[i].z=tempdata4[i % nentries_in_file].z;
	  (*data)[i].intensity=tempdata4[i % nentries_in_file].intensity;
	}
      } else {
	for(i=0;i<*nentries;i++) {
	  (*data)[i].ID=tempdata4[i].ID;
	  (*data)[i].x=tempdata4[i].x;
	  (*data)[i].y=tempdata4[i].y;
	  (*data)[i].z=tempdata4[i].z;
	  (*data)[i].intensity=tempdata4[i].intensity;
	}
      }
      free(tempdata4);
      break;
  }
  fclose(pFile);

  //--------------------------------------




  //------------ OPTIONAL ITEMS ----------

  //Scale up the intensity values if desired
  if(scaling != 1.0) for(i=0;i<*nentries;i++) (*data)[i].intensity *= scaling;

  //--------------------------------------


}
