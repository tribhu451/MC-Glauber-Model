#include <iostream>
#include <string>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "input_data.h"
#include "read_input_data.h"
#include "mc_glau.h"

using std::cout;
using std::endl;
using std::to_string;


int main(int argc, char **argv)
{
  // reading the input data from the file
  string input_file_name;char* event_no_s ;
  if(argc == 3){event_no_s = *(argv+1); input_file_name = *(argv+2);}
  else{cout<<"plz give 2 arguments\n1st argument - no of events you want to generate \n"
      "2nd argument - input filename"<<endl; exit(1);}
  ReadInputPars reader;
  InputData InData;
  reader.read_input_data_(&InData, input_file_name);
  int event_no = atof(event_no_s) ;
  cout<<"total no of events : "<<event_no<<endl;
  
  
  cout<<"      Monte Carlo Glauber Model     "<<endl;
  cout<<"      *************************     "<<endl;
  cout<<"[Info] "<<InData.projectile<<"+"<<InData.target<<" at "<<InData.SNN<<"GeV\n"<<endl;

  string collision_species = InData.projectile +"_"+ InData.target +"_" + to_string(int(InData.SNN)) ;
  
  mc_glau* MC = new mc_glau(&InData);


  // root file init
  double npart,ncoll,nch,b;
  TFile* root_file = new TFile (InData.root_output_file_name.c_str(),"recreate");
  TTree* root_tree= new TTree(collision_species.c_str(),collision_species.c_str());
  root_tree-> Branch("npart",&npart,"npart/D");
  root_tree->Branch("ncoll",&ncoll,"ncoll/D");
  root_tree->Branch("nch",&nch,"nch/D");
  root_tree->Branch("b",&b,"b/D");
  
  // event loop starts
  for(int i = 0; i<event_no; i++) 
    {
      MC->event();
      npart =  MC->get_npart();
      ncoll =  MC->get_ncoll();
      nch   =  InData.npp * ( (1.0-InData.xhard)*0.5*MC->get_npart() 
			      + InData.xhard*MC->get_ncoll() ) ;
      b = MC->get_impactf();
      double e2 = MC->eccen2();
      double phi2 = MC->phi2();
      if(npart >= 2) root_tree->Fill();
      cout<<i<<"\t"<<b<<"\t"<<npart<<"\t"<<ncoll<<"\t"<<nch <<"\t"<<e2<<"\t"<<phi2<<endl;
    } 
  
  root_tree->Write();
  root_file->Close(); 
  
  delete MC;
  
  return 0;
}

