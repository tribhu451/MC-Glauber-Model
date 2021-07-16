// root -b -q "centrality_determin.C(\"au1.root\")"

#include <iostream>
#include "TTree.h"
#include "TFile.h"
#include "TObject.h"
#include "TROOT.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TROOT.h"

using namespace std;

void centrality_determin_2()
 {
  
   
   double nch; double imp_b; double npart; double ncoll ;
   double max_b; double max_npart; double max_ncoll;
   double npart_should; double ncoll_should; double nch_should;
   
   
   // inputs [
   int bins=20000;
   double min_bin_val = 0.;
   double max_bin_val = 4000.;
   
   double xcut[20]={0.0};
   
   const int pr_f = 12;
   double percent[pr_f] = {0,  5,  10 , 20, 30, 40, 45 , 50, 55, 60, 65};
   //double percent[pr_f] = {0,  6,  15, 25 , 35,45,55};
   
   TFile* f1=new TFile("Pb_Pb_2760_min_bias_xhard_0.14.root","read");
   TH1F* h1 = new TH1F("h1","",bins,min_bin_val,max_bin_val);
   TTree* TreeAu =new TTree();
   f1->GetObject("Pb_Pb_2760",TreeAu);
   TreeAu->SetBranchAddress("nch",&nch);
   TreeAu->SetBranchAddress("b",&imp_b);
   TreeAu->SetBranchAddress("npart",&npart);
   TreeAu->SetBranchAddress("ncoll",&ncoll);
   // inputs ]
   
     
   double maxch = 0 ; 
   int Entries=TreeAu->GetEntries();        
   for(int j=0;j<Entries;j++)
     {
         TreeAu->GetEntry(j);
         if(nch > max_bin_val)
	  { cout<<"nch > maximum bin range : "<<nch<<endl; }
         if(npart < 2 ) continue ;
         h1->Fill(nch);
         if (nch > maxch){ maxch = nch ; max_b = imp_b ; max_npart=npart; max_ncoll=ncoll;} 
         if(imp_b < 0.5 ){npart_should = npart; ncoll_should = ncoll; nch_should = nch;}
     }
   

   cout << "MAXCH = "<< maxch << " MAX_B = "<<max_b << " MAX_NPART = " << max_npart << " MAX_NCOLL = " << max_ncoll << endl;
   cout << "MAX-NCH-SHOULD > "<< nch_should << " MAX-NPART-SHOULD > " << npart_should << " MAX-NCOLL-SHOULD > " << ncoll_should << endl;


   int n1=0;
   double Sum=0.0;
   int per_count = 0 ;
   xcut[n1] = maxch ;


   for(int j=bins; j>=1; j--)
     {
       double aSum = h1->GetBinContent(j);
       //if(aSum == 0){cout<<"Bin contetnt zero at "<<j<<"th bin"<<endl; }
       double perc=percent[per_count+1] - percent[per_count];   // perc=10 means 0-10% , 10-20%, 20-30%....                                     
       if(Sum>(Entries*perc*0.01))
	 {
           n1=n1+1; 
	   double p0=h1->GetBinCenter(j);                                       
	   xcut[n1]=p0;
           Sum=0.0;
	   
           per_count += 1;

           if(bins == 1){break;}  
         }
       else
	 {
	   Sum=Sum+aSum;
	 }
     }
   



 
   for(int i=0; i<pr_f-1; i++)
     {
       cout<<percent[i]<<"-"<<percent[i+1]<<" \%"<< endl ;
     }
   
   for(int i=0 ; i <= n1 ; i++)
     {
         cout<<"MULT : "<<2.25*xcut[i]<<" - "<<2.25*xcut[i+1]<<endl;
     }
   
   
   double bb[n1];
   int npartt[n1];
   for(int j=0;j<Entries;j++)
     {
       TreeAu->GetEntry(j);

       for (int i = 0; i <= n1 ; i ++){
	 if(nch < xcut[i]+1 && nch > xcut[i]-1)
	   { bb[i] = imp_b; npartt[i] = npart; }
	 
       }
     }
   

   for (int i = 0; i <= n1 ; i ++)
     {
       cout<<"multiplicity "<<2.25*xcut[i]<<"\t=>\tb : "<<bb[i]<<"\tnpart : "<<npartt[i]<<endl;
     }



   // for xhard determination
    double avg_b[n1];
    double avg_npart[n1];
    double avg_nch[n1];

   for(int i = 0; i <= n1-1 ; i++)
    {
     avg_b[i] = (2./3.) * ( pow(bb[i+1],3.0) - pow(bb[i],3.0) ) / ( pow(bb[i+1],2.0) - pow(bb[i],2.0) ) ;
     cout << "Avg b["<< i << "]  " << avg_b[i] << endl; 

    }


   for(int j=0;j<Entries;j++)
     {
       TreeAu->GetEntry(j);

       for (int i = 0; i <= n1-1 ; i ++){
	 if(imp_b < avg_b[i]+0.1 && imp_b > avg_b[i]-0.1)
	   { avg_npart[i] = npart; avg_nch[i] = 2.25*nch; } // 2.25 multiplied here
	 
       }
     }
   

   for (int i = 0; i <= n1-1 ; i ++)
     {
       cout<<"avg npart : "<<avg_npart[i]<<"\t avg_nch : "<<avg_nch[i] << "\t avgnch/(0.5*npart) : "<< avg_nch[i]/(0.5*avg_npart[i]) << endl;
     }


/*   std::ofstream Filex;
   Filex.open("cent_vs_avgnc_by_npart_pair.dat");

  for(int i=0; i<n1-1; i++)
     {

       Filex<<(percent[i]+percent[i+1])/2.0<<  "\t" << avg_nch[i]/(0.5*avg_npart[i]) << endl ;
     }
          cout<< "cent_vs_avgnc_by_npart_pair.dat created ..." << endl;
*/ 
   
 }
