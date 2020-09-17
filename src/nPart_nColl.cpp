#include "nPart_nColl.h"

using namespace std;

void nPart_nColl::Calculate_nPart_nColl(double* vxA,double* vyA,double* vxB,double* vyB)
     {
       
       aN_coll=0;
       aN_part=0;
       for(int k=2; k<=6; k++){PhiN[k]=0.0; eccentricity[k]=0.0;
       }
       
       for(int i=0;i<A;i++){occA[i]=0;}
       for(int i=0;i<B;i++){occB[i]=0;}
       
       for (int i=0; i<A; i++)
	 {
	   for (int j=0; j<B; j++)
	      {  
		double d=TMath::Sqrt( TMath::Power((vxB[j]-vxA[i]),2) + 
				      TMath::Power ( (vyB[j]-vyA[i]),2));
		double D=TMath::Sqrt(sigma/ (TMath::Pi())); 
	   
		if( d <= D)
		  { 
		    
		    aN_coll=aN_coll+1;
		    Ncoll_x[aN_coll]=(vxA[i]+vxB[j])/2;
		    Ncoll_y[aN_coll]=(vyA[i]+vyB[j])/2;
		    
		    if(occA[i]==0)
		      { //START1
			aN_part=aN_part+1;occA[i]=1;Npart_x[aN_part]=vxA[i]; Npart_y[aN_part]=vyA[i];
		      } //END1
		    
                    
		    if(occB[j]==0)
		      {//START2
			aN_part=aN_part+1;occB[j]=1;Npart_x[aN_part]=vxB[j]; Npart_y[aN_part]=vyB[j];
		      }//END2
		    
		    
		  }                                                           //End of if loop
	      }                                                          //End of j loop
	 }                                                         //End of i loop
       
       
       if(aN_part != 0 )
	 {  //IF 1
	   
	   //############shifting the energy distributions center to (0,0,0)#################//
	   double xref1=0.0;
	   double yref1=0.0;
	   double wref1=0.0;
	   for(int k=1;k<=aN_part;k++)
	     { 
	       xref1=xref1+(Npart_x[k]*(0.5*npp*X_hard));
	       yref1=yref1+(Npart_y[k]*(0.5*npp*X_hard));
	       wref1=wref1+(0.5*npp*X_hard);
	     }
	   
	   double xref2=0.0;
	   double yref2=0.0;
	   double wref2=0.0;
	   for(int k=1;k<=aN_coll;k++)
	     { 
	       xref2=xref2+(Ncoll_x[k]*(npp*(1-X_hard)));
	       yref2=yref2+(Ncoll_y[k]*(npp*(1-X_hard)));
	       wref2=wref2+(npp*(1-X_hard));
	     }
	   
	   double xAverage=((xref1+xref2)/(wref1+wref2));
	   double yAverage=((yref1+yref2)/(wref1+wref2));
	   
	   // cout<<xAverage<<"  "<<yAverage<<"\n";
	   
	   for(int k=1;k<=aN_part;k++)
	     { 
	       Ra= TMath::Sqrt( TMath::Power(Npart_x[k]-xAverage,2)+TMath::Power( Npart_y[k]-yAverage,2) );
	       Theta =(TMath::ATan2(Npart_y[k]-yAverage,Npart_x[k]-xAverage));
	       Npart_r[k]=Ra;
	       Npart_phi[k]=Theta;
	     }
	   for(int k=1;k<=aN_coll;k++)
	     {      
	       Ra=TMath::Sqrt( TMath::Power(Ncoll_x[k]-xAverage,2)+ TMath::Power(Ncoll_y[k]-yAverage,2));
	       Theta=(TMath::ATan2(Ncoll_y[k]-yAverage,Ncoll_x[k]-xAverage)); 
	       Ncoll_r[k]=Ra;
	       Ncoll_phi[k]=Theta;
	     }
	   //###########shifting done and the points set in polar co-ordinate#####################//
	   
	   
	   
	   
	   int Total_Nch=aN_part+aN_coll;
	   
	   for(int k=1;k<=aN_part;k++)
	     { //START3
	       Nch_r[k]=Npart_r[k];
	       Nch_phi[k]=Npart_phi[k];
	       Nch_value[k]=((0.5)*(npp)*(X_hard));
	     } //END3
	   
	   for(int k=1;k<=aN_coll;k++)
	     { //START4
	       Nch_r[k+aN_part]=Ncoll_r[k];
	       Nch_phi[k+aN_part]=Ncoll_phi[k];
	       Nch_value[k+aN_part]=((1-X_hard)*npp);
	     } //END4 
	   
	   
	   
	   
	   for(int N=2; N<=6; N++)
	     {   //N-loop Start
	       
	       double RXA=0.0;
	       double RXB=0.0;
	       double RXC=0.0;
	       
	       for(int k=1; k<=Total_Nch;k++)
		 { //START5
		   RXA=RXA+(Nch_value[k]* TMath::Power(Nch_r[k],(N)) );  
		   RXB=RXB+(Nch_value[k]* TMath::Power(Nch_r[k],(N))*TMath::Cos((N)*Nch_phi[k] )  );    
		   RXC=RXC+(Nch_value[k]* TMath::Power(Nch_r[k],(N))*TMath::Sin((N)*Nch_phi[k] )  );
		 } //END5
	       
	       double R1=-(RXC/RXA);
	       double R2=-(RXB/RXA);
	       eccentricity[N]= TMath::Sqrt((R1*R1)+(R2*R2));
	       PhiN[N]=((TMath::ATan2(R1,R2)))/ (N);
	       
	       
	     } //N-loop end
	   
	 }
       else
	 {
	   for(int g=2; g<=6; g++)
	     {
	       eccentricity[g]=0.0;PhiN[g]=0.0;
	     }
	 }  //ENDIF 1
       
     }


double nPart_nColl::Get_Npart(){return aN_part;}
double nPart_nColl::Get_Ncoll(){return aN_coll;}
double nPart_nColl::Get_Eccentricity2(){return eccentricity[2];}
double nPart_nColl::Get_PhiN2(){return PhiN[2];}
double nPart_nColl::Get_Eccentricity3(){return eccentricity[3];}
double nPart_nColl::Get_PhiN3(){return PhiN[3];}
double nPart_nColl::Get_Eccentricity4(){return eccentricity[4];}
double nPart_nColl::Get_PhiN4(){return PhiN[4];}
double nPart_nColl::Get_Eccentricity5(){return eccentricity[5];}
double nPart_nColl::Get_PhiN5(){return PhiN[5];}
double nPart_nColl::Get_Eccentricity6(){return eccentricity[6];}
double nPart_nColl::Get_PhiN6(){return PhiN[6];}
void nPart_nColl::Set_Datas(int aa,int ab, double ac,double anpp, double aX_hard){A=aa;B=ab;sigma=ac;npp=anpp;X_hard=aX_hard;}




