
#include <stdio.h>
#include <math.h>
#include "solve.h"

void Ini_label(int *lab, int myid){
  /** correspond the rank to 2_bit number*/
  lab[0] = myid/4;
  lab[1] = (myid - lab[0]*4)/2;
  lab[2] = myid - lab[0]*4 - lab[1]*2;
}

