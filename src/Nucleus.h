#include <iostream>
#include <TF1.h>
#include "TMath.h"
#include <TRandom.h>
#include <TRandom3.h>
#include <TUUID.h>

using namespace std;

class Nucleus
  {

  public:
    Nucleus(int ,int,double,double,double,double);
    ~Nucleus();
    void generate_nucleons_A_position(double* ,double* ,double* );  // function to generate position of nucleus A
    void generate_nucleons_B_position(double* ,double* ,double* );  // function to generate position of nucleus B
    void Set_Impact_Parameter(double );                       
    double Get_ThetaA();                                            // orientation of nucleus during generation
    double Get_PhiA();
    double Get_ThetaB();
    double Get_PhiB();    
    double etaA;
    double psiA;
    double etaB;
    double psiB;
    
  private:
    TRandom3* tr1;
    TRandom3* tr2;
    TRandom3* tr3;
    TRandom3* tr4;
    TRandom3* tr5;
    TRandom3* tr6;
    TRandom3* tr7;
    TRandom3* tr8;
    TRandom3* tr9;    
    TRandom3* trA;
    TRandom3* trB;    
    TRandom3* rsduse;
    TF1* angl;

    double X[500];double Y[500];double Z[500];
    int A;             //Nucleus-1
    int B;             //Nucleus-2
    double b;          //Impact Parameter
    double R;          // Wood Saxon (Radius)
    double dlt;        // Wood Saxon (a)
    double ZETA;
    double BETA2;      // Wood Saxon (deformation parameter)
    double BETA4;
        
}; 
