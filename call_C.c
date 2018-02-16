#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "dustmap_c.c"

//gcc call_C.c -g -o call_C

int main (){
  double balance[5];
  printf("(Hi0)\n" );
  char *args[55];
  balance[2]=9;
  int calc=sqrt(balance[2]);
  printf("(%f)\n",calc);
  printf("(leaving call_C)\n" );

    args[0] = balance;
    args[1] = &calc;
    args[2] = "9";
    args[3] = "5";
  //&balance = {1000.0, 2.0, 3.4, 7.0, 50.0};
  //dustmap_c(1, args);
  printf("(%f)\n",sqrt(PI) );

  return 0;


}
