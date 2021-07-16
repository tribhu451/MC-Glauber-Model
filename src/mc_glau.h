// *************************** //
// *  25 feb. 2021           * //
// *  tribhu.451@gmail.com   * //
// *  version 2.0            * //
// *************************** //


//  This class calculates npart, ncoll and eccentricity
//  in mc glauber model

#pragma once

#include<fstream>
#include<string>
#include<sstream>
#include<iostream>
#include<fstream>
#include <TRandom3.h>
#include <TF1.h>
#include "TMath.h"
#include "input_data.h"
	
using std::cout;
using std::endl;

class mc_glau
{
  
 public:
  mc_glau(InputData *InData1);
  ~mc_glau();
  void event();
  inline int get_npart(){return NPART;}
  inline int get_ncoll(){return NCOLL;}
  inline double get_impactf(){return IMPACT_PARAM;}
  void calculate_eccentricity(int aN_part,int aN_coll,double *Npart_x,
                     double *Npart_y,double *Ncoll_x,double *Ncoll_y);

  inline double eccen2(){return eccentricity[2];}
  inline double phi2(){return PhiN[2];}
  inline double eccen3(){return eccentricity[3];}
  inline double phi3(){return PhiN[3];}
  inline double eccen4(){return eccentricity[4];}
  inline double phi4(){return PhiN[4];}
  inline double eccen5(){return eccentricity[5];}
  inline double phi5(){return PhiN[5];}
  inline double eccen6(){return eccentricity[6];}
  inline double phi6(){return PhiN[6];}




  

  
 private:
  InputData *InData;
  
  int A; // mass no. of projrctile nucleus
  int B; // mass no. of target nucleus
  double sigma; // energy-> cross-section
  
// projectie wood-Saxon parameters
  double p_radius;        
  double p_dlt;        
  double p_beta2;
  double p_beta4;

// target wood-Saxon parameters
  double t_radius;        
  double t_dlt;        
  double t_beta2;
  double t_beta4;

// npart, ncoll and impact parameter
  int NPART;
  int NCOLL;
  double IMPACT_PARAM;

// two-component energy deposition
  double npp;
  double X_hard;

 //impact parameter range
  double bmin,bmax;

// eccentricity and partcipant plane angle
  double eccentricity[10];
  double PhiN[10];


  void generate_nucleus(double* X1, double* Y1,double* Z1,int A,
			double R, double dlt, double BETA2, double BETA4, double etaA, double psiA);
  
  void calculate_npart_ncoll(double* vxA,double* vyA,double* vxB,double* vyB, int &Npart, 
			     int &Ncoll, double* Npart_x, double* Npart_y, double* Ncoll_x, double* Ncoll_y);
  
  void shift_nucleus(double* X1, double* Y1, double* Z1,int A, double b,
                     double zhi,double* X2, double* Y2, double* Z2 );




void set_mc_glau_params()
{
  npp = InData->npp;
  X_hard = InData->xhard;

  // projectile nucleus
  if(InData->projectile == "Au") {A= 197; p_radius = 6.42; p_dlt=0.41; p_beta2 = -0.13; p_beta4 =0.0;} //arXiv: 1409.8375 [Table. 1]
  if(InData->projectile == "Au2") {A= 197; p_radius = 6.37; p_dlt=0.53; p_beta2 = 0.0; p_beta4 =0.0;} //arXiv: 1409.8375 [Table. 1]
  else if(InData->projectile == "Pb") {A= 208; p_radius = 6.66; p_dlt=0.45; p_beta2 = 0.0; p_beta4 =0.0;}
  else if(InData->projectile == "p") {A= 1; p_radius = 6.37; p_dlt=0.53; p_beta2 = 0.0; p_beta4 =0.0;} // woods-saxon params can be ignored-
  //because we will sample only one nucleon.
  else if(InData->projectile == "U") {A= 238; p_radius = 6.86; p_dlt=0.42; p_beta2 = 0.265; p_beta4 =0.093;}
  else {cout<<"projectile not recognized, it's : "<<InData->projectile<<endl; exit(1);}
  
  // target nucleus
  if(InData->target == "Au") {B= 197; t_radius = 6.42; t_dlt=0.41; t_beta2 = -0.13; t_beta4 =0.0;}
  if(InData->target == "Au2") {B= 197; t_radius = 6.37; t_dlt=0.53; t_beta2 = 0; t_beta4 =0.0;} 
  else if(InData->target == "Pb") {B= 208; t_radius = 6.66; t_dlt=0.45; t_beta2 = 0.0; t_beta4 =0.0;}
  else if(InData->target == "p") {B= 1; t_radius = 6.37; t_dlt=0.53; t_beta2 = 0.0; t_beta4 =0.0;}
  else if(InData->target == "U"){B =238;  t_radius = 6.86; t_dlt=0.42; t_beta2 = 0.265; t_beta4 =0.093;}
  else {cout<<"target not recognized, it's : "<<InData->target<<endl; exit(1);}
  
  // collision energy
  if(InData->SNN == 200.0){sigma =4.2;}
  else if(InData->SNN == 62.4){sigma =3.155;}
  else if(InData->SNN == 2760.0){sigma =6.4;}
  else if(InData->SNN == 5500.0){sigma =7.2;}
  else{cout<<"SNN(energy) not recognised"<<endl; exit(1);}
  
  // impact parameter range
  bmin = InData->bmin; bmax = InData->bmax;
  
  if(InData->projectile != "p")
    {
      cout<<"[Info] projectile mass number : "<<A<<endl;
      cout<<"[Info] wood-saxon | projectile Radius : "<<p_radius<<"   a: "<<p_dlt<<endl;
      cout<<"[Info] deformation parameter of projectile | Beta-2 : "<<p_beta2<<" and Beta-4 : "<<p_beta4<<"\n"<<endl; 
    }
  
  if(InData->target != "p")
    {
      cout<<"[Info] target mass number : "<<B<<endl;
      cout<<"[Info] wood-saxon | target Radius : "<<t_radius<<"   a: "<<t_dlt<<endl;
      cout<<"[Info] deformation parameter of target | Beta-2 : "<<t_beta2<<" and Beta-4 : "<<t_beta4<<"\n"<<endl;
    }
  
}



};
















