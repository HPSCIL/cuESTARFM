//#include"cuLayer.h"
typedef struct { 
 // char outpath[1000];
 // int A;
  int class_num;
 // int pf;
 // int pc;
 // float L_ERR;
  int NUM_PREDICTIONS;
  float uncertain;
  int NUM_PAIRS;             /* number of input data pair */ 
  int WIN_SIZE;              /* window size determined by the MAX_SEARCH_DIS */
} PARAMETER;
