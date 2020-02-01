#include <iostream>
#include <TMath.h>

using namespace std;

class nPart_nColl
   {

   public:     
     void Calculate_nPart_nColl(double* ,double* ,double* ,double*);
     void Set_Datas(int , int , double , double , double);
     double Get_Npart();
     double Get_Ncoll();
     double Get_Eccentricity2();
     double Get_PhiN2();
     double Get_Eccentricity3();
     double Get_PhiN3();
     double Get_Eccentricity4();
     double Get_PhiN4();
     double Get_Eccentricity5();
     double Get_PhiN5();
     double Get_Eccentricity6();
     double Get_PhiN6();
     
   private:
     double occA[500];double occB[500];
     double Ncoll_x[3000]; double Ncoll_y[3000];
     double Ncoll_r[3000]; double Ncoll_phi[3000];
     double Npart_x[1000]; double Npart_y[1000];
     double Npart_r[1000]; double Npart_phi[1000];
     double Nch_r[4000]; double Nch_phi[4000];
     double Nch_value[4000];
     int A;             //Nucleus-1
     int B;             //Nucleus-2
     double sigma;      //Energy
     int aN_coll;
     int aN_part;
     double PhiN[10];
     double eccentricity[10];
     double Ra; double Theta;
     double npp;
     double X_hard;
          
};
