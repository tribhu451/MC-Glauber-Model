#ifndef SRC_DATA_H
#define SRC_DATA_H

#include<iostream>
#include<string>
using std::string;
 
typedef struct input_data
{

  double SNN;
  string projectile;
  string target;  

  double bmin,bmax;

  double xhard;
  double npp;

  string root_output_file_name ;
    
}InputData;

#endif
