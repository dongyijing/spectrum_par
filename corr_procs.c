#include <stdio.h>
#include "solve.h"

int corr_procs (int myid, int direction){
  if (direction == 0){
    switch(myid){
	  case 0: return 4;
	  case 1: return 5;
	  case 2: return 6;
	  case 3: return 7;
	  case 4: return 0;
	  case 5: return 1;
	  case 6: return 2;
	  case 7: return 3;
    }
  }
  if (direction == 1){
    switch(myid){
	  case 0: return 2;
	  case 1: return 3;
	  case 2: return 0;
	  case 3: return 1;
	  case 4: return 6;
	  case 5: return 7;
	  case 6: return 4;
	  case 7: return 5;
	}
  }
  if (direction == 2){
    switch(myid){
	  case 0: return 1;
	  case 1: return 0;
	  case 2: return 3;
	  case 3: return 2;
	  case 4: return 5;
	  case 5: return 4;
	  case 6: return 7;
	  case 7: return 6;
	}
  }
}

