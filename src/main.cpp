#include <iostream>
#include <fstream>
#include <TMath.h>
#include <TRandom.h>
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>
#include "Nucleus.h"
#include "nPart_nColl.h"

using namespace std;
int main()
  {

  
    // INPUTS
    int A=238;
    int B=238;
    double sigma=4.2;          
    int Event_No=50;
    double Radius=6.81;       
    double dlt=0.54;        
    double BETA2=0.28;
    double BETA4=0.093;
    double npp=2.49;
    double X_hard=0.87;


//     FILE setup
    ofstream myfile;
    myfile.open("Urandomb.dat");
    TFile* File2= new TFile ("Urandomb.root","recreate");




//     TREE setup
    double Impact;double E2; double Psi2; double E3; double Psi3;
    double E4; double Psi4; double E5; double Psi5;  double E6; double Psi6;
    double ota;double otb;double opa; double opb;double Nch;int Npart;int Ncoll;
    TTree* tree1= new TTree("TreeU","Tree for U+U at 200 Gev");
    tree1-> Branch("N_Part",&Npart,"NPart/I");
    tree1->Branch("N_Coll",&Ncoll,"NColl/I");
    tree1->Branch("NCh",&Nch,"Nch/D");
    tree1->Branch("Impact_parameter",&Impact,"b/D");
    tree1->Branch("Ori_TA",&ota,"ota/D");
    tree1->Branch("Ori_PA",&opa,"opa/D");
    tree1->Branch("Ori_TB",&otb,"otb/D");
    tree1->Branch("Ori_PB",&opb,"opb/D");
    tree1->Branch("E2",&E2,"E2/D");
    tree1->Branch("Psi_2",&Psi2,"Psi2/D");
    tree1->Branch("E3",&E3,"E3/D");
    tree1->Branch("Psi_3",&Psi3,"Psi3/D");
    tree1->Branch("E4",&E4,"E4/D");
    tree1->Branch("Psi_4",&Psi4,"Psi4/D");
    tree1->Branch("E5",&E5,"E5/D");
    tree1->Branch("Psi_5",&Psi5,"Psi5/D");
    tree1->Branch("E6",&E6,"E6/D");
    tree1->Branch("Psi_6",&Psi6,"Psi6/D");

    
    // Random number of Impact parameter
    TRandom* t1=new TRandom();
    t1->SetSeed(0);
    long kss=t1->GetSeed();
    gRandom->SetSeed(kss);
    TF1* f1= new TF1("f1","x",0.0,30.0);
    
    
    Nucleus* NN = new Nucleus(A,B,Radius,dlt,BETA2,BETA4);
    nPart_nColl* PC = new nPart_nColl(A,B,sigma,npp,X_hard);


    double NPart;double NColl;double NCh;double b;
    double Eccen2;double Eccen3;double Eccen4;double Eccen5;double Eccen6;
    double PhiN2;double PhiN3;double PhiN4;double PhiN5;double PhiN6;

    
    double XA[A];double YA[A];double ZA[A];
    double XB[B];double YB[B];double ZB[B];
    for(int j=0;j<=A;j++){XA[j]=0.0;YA[j]=0.0;ZA[j]=0.0;}
    for(int j=0;j<=B;j++){XB[j]=0.0;YB[j]=0.0;ZB[j]=0.0;}

    



    int event=0;
    do
       {
	 
	 if((event/1000)*1000 == event){cout<<event<<"\n";}
	 b=f1->GetRandom(0.0,20.0);                     //Randomly Impact parameter set
	 NN->Set_Impact_Parameter(b);                    //Impact parameter set to generate Nucleons
	 NN->generate_nucleons_A_position(XA,YA,ZA);     //Nucleons of A generated
	 double OTA=NN->Get_ThetaA();                    //Orientation Of Nucleus A (Theta)
	 double OPA=NN->Get_PhiA();                      //---do-----A (Phi)
	 NN->generate_nucleons_B_position(XB,YB,ZB);     //Nucleons of B generated
	 double OTB=NN->Get_ThetaB();                    //Orientation Of Nucleus B (Theta)
	 double OPB=NN->Get_PhiB();                      //---do-----B (Phi)
	 PC->Calculate_nPart_nColl(XA,YA,XB,YB);         //N_part and N-coll calculated
	 NColl=PC->Get_Ncoll();                          // Got the Ncoll
	 NPart= PC->Get_Npart();                         //Got the Npart
	 NCh=PC->Get_Multiplicity();                     // Got multiplicity
	 Eccen2=PC->Get_Eccentricity2();                 //Got Eccentricity
	 PhiN2=PC->Get_PhiN2();                          //Got Participant plane angle
	 Eccen3=PC->Get_Eccentricity3();                 //Got Eccentricity
	 PhiN3=PC->Get_PhiN3();                          //Got Participant plane angle
	 Eccen4=PC->Get_Eccentricity4();                 //Got Eccentricity
	 PhiN4=PC->Get_PhiN4();                          //Got Participant plane angle
	 Eccen5=PC->Get_Eccentricity5();                 //Got Eccentricity
	 PhiN5=PC->Get_PhiN5();                          //Got Participant plane angle
	 Eccen6=PC->Get_Eccentricity6();                 //Got Eccentricity
	 PhiN6=PC->Get_PhiN6();                          //Got Participant plane angle
	 
	 if(NPart != 0 )
	     { 
	   /*myfile<<b<<"  "<<OTA<<"  "<<OPA<<"  "<<OTB<<"  "<<OPB<<"  "<<NPart<<"  "<<NColl<<"  "<<NCh<<"  "<<
	     Eccen2<<"  "<<PhiN2<<" "<<Eccen3<<"  "<<PhiN3<<"  "<<Eccen4<<"  "<<PhiN4<<"  "<<Eccen5<<"  "<<PhiN5<<"  "<<Eccen6<<"  "<<PhiN6<<"\n";*/
	   
	       Npart=NPart;
	       Ncoll=NColl;
	       Nch=NCh;
	       Impact=b;
	       ota=OTA;
	       otb=OTB;
	       opa=OPA;
	       opb=OPB;
	       E2=Eccen2;
	       Psi2=PhiN2;
	       E3=Eccen3;
	       Psi3=PhiN3;
	       E4=Eccen4;
	       Psi4=PhiN4;
	       E5=Eccen5;
	       Psi5=PhiN5;
	       E6=Eccen6;
	       Psi6=PhiN6;
	       tree1->Fill();
	       
	       event=event+1;
	     }   // IF loop end
	 
       }
    while(event<Event_No);  // event loop
    
    
    myfile.close();
    tree1->Write();
    File2->Close();
    
    return 0;
  }

