#include<iostream>
#include"PARAMETER.h"
#include"cuLayer.h"
void Re_fusion2(const char * BufferIn0,const char * BufferIn1,const char * BufferIn2,const char * BufferIn3,const char * BufferIn4,const char * BufferOut,int win_size,int class_num,float M_err);
int parseParameters(char *fname, CuLayer *psensor,PARAMETER *par);
void usage(char *command);
void Re_fusion3(CuLayer *psensor,PARAMETER *par);