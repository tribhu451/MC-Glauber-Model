#include "mc_glau.h"

//#define ROTATED_SHIFT 

using std::cout;
using std::endl;

mc_glau::mc_glau(InputData *InData1)
{
  InData = InData1; 
  set_mc_glau_params();

}

mc_glau::~mc_glau(){}

// This is the master function
void mc_glau::event()
{
 

   NPART = 1E5;
   NCOLL = 1E5;
   IMPACT_PARAM = 1E5;

  
  TRandom* t1=new TRandom();
  t1->SetSeed(0);
  long kss=t1->GetSeed();
  gRandom->SetSeed(kss);
  TF1* f1= new TF1("f1","x",0.0,25.0);
  TF1* f2= new TF1("f2","sin(x)",0.0,TMath::Pi());
 
 
  double XA[A];double YA[A];double ZA[A];
  double XB[B];double YB[B];double ZB[B];
  double npart_x[500],npart_y[500];
  double ncoll_x[10000],ncoll_y[10000];
  for(int j=0;j<=A;j++){XA[j]=0.0;YA[j]=0.0;ZA[j]=0.0;}
  for(int j=0;j<=B;j++){XB[j]=0.0;YB[j]=0.0;ZB[j]=0.0;}
  for(int j=0;j<500;j++){npart_x[j]=0.0;npart_y[j]=0.0;}
  for(int j=0;j<10000;j++){ncoll_x[j]=0.0;ncoll_y[j]=0.0;}


  //generate orientation angles of target & projectile ...
  double p_ori_theta = f2->GetRandom(0.0,TMath::Pi());  
  double t_ori_theta = f2->GetRandom(0.0,TMath::Pi());  
  double p_ori_phi = (2.0*TMath::Pi())*(t1->Rndm());
  double t_ori_phi = (2.0*TMath::Pi())*(t1->Rndm());
  //cout<<"[Info] (projectile orientation) p_theta: "<<p_ori_theta<<" p_phi: "<<p_ori_phi<<endl;
  //cout<<"[Info] (target orientation) t_theta: "<<t_ori_theta<<" t_phi: "<<t_ori_phi<<endl;

  //generate nucleus
  generate_nucleus(XA,YA,ZA,A,p_radius,p_dlt,p_beta2,p_beta4,p_ori_theta,p_ori_phi);
  generate_nucleus(XB,YB,ZB,B,t_radius,t_dlt,t_beta2,t_beta4,t_ori_theta,t_ori_phi);

  // generate impact parameter between bmin-bmax ...
  double b=f1->GetRandom(bmin,bmax); 
  IMPACT_PARAM = b;                     
  //cout<<"[Info] b = "<<b<<" (fm)"<<endl;


   double zhi;
#ifdef ROTATED_SHIFT
   zhi=(2.0*TMath::Pi())*(t1->Rndm());
#else
   zhi =0.0;
#endif
   
   //shifting of nucleus 
   shift_nucleus( XA, YA,  ZA, A, +b/2.0, zhi, XA, YA, ZA );
   shift_nucleus( XB, YB,  ZB, B, -b/2.0, zhi, XB, YB, ZB);
   

   // calculating npart & ncoll ...
   calculate_npart_ncoll(XA,YA,XB,YB,NPART,NCOLL,npart_x,npart_y,ncoll_x, ncoll_y);                    
   //cout<<"[Info] No. of participants : "<<NPART<<endl;
   //cout<<"[Info] No. of binary collisions : "<<NCOLL<<endl;
   

   // calculating eccentricity ...
   calculate_eccentricity(NPART,NCOLL,npart_x,npart_y,ncoll_x,ncoll_y);
   
}


void mc_glau::generate_nucleus(double* X1, double* Y1,double* Z1,int A,
			       double R, double dlt, double BETA2, double BETA4, double etaA, double psiA)
{    
  double X[500];double Y[500];double Z[500];
  
  TRandom3* tr1 = new TRandom3();
  tr1->SetSeed(0);
  
  double CMx=0.0;double CMy=0.0;double CMz=0.0;
  int count=0;
  
  do
    {
      double r=(11.0)*(tr1->Rndm());
      double Theta=(TMath::Pi())*(tr1->Rndm());
      double Phi=((2.0)*TMath::Pi())*(tr1->Rndm());
      double test=tr1->Rndm();
      
      
      double Y20=0.25*TMath::Sqrt(5.0/TMath::Pi())*
	(3*TMath::Cos(Theta)*TMath::Cos(Theta)-1.0);
      double Y40=(3.0/(16.0*TMath::Sqrt(TMath::Pi())))* 
	((35*TMath::Power(TMath::Cos(Theta),4))-
	 (30*TMath::Power(TMath::Cos(Theta),2))+3);
      double RAT= R*(1+(BETA2*Y20)+(BETA4*Y40));
      double rho=(1.0/60.0)*(r*r*(TMath::Sin(Theta)))/(1.0+(TMath::Exp((r-RAT)/dlt)));
      
      
      if(test < rho )
	{      
	  
	  X[count]= (r*TMath::Sin(Theta)*TMath::Cos(Phi));
	  Y[count]=(r*TMath::Sin(Theta)*TMath::Sin(Phi));
	  Z[count]=(r*TMath::Cos(Theta));
	  CMx=CMx+X[count]; CMy=CMy+Y[count] ;CMz=CMz+Z[count];    
	  count=count+1;
	}   
    }   
  while(count<A);
  
  CMx=CMx/A;CMy=CMy/A;CMz=CMz/A;
  
  for(int j=0;j<A;j++){ X[j]=X[j]+(-CMx); Y[j]=Y[j]+(-CMy);Z[j]=Z[j]+(-CMz);}
  
  //etaA - nucleus orientaton angle (theta)
  //psiA - nucleus orientation angle (phi)
  
  for(int j=0;j<A;j++)
    {
      X1[j]=(TMath::Cos(psiA)*TMath::Cos(etaA)*X[j])+(-TMath::Sin(psiA)*Y[j])+(-TMath::Cos(psiA)*TMath::Sin(etaA)*Z[j]);
      Y1[j]=(TMath::Sin(psiA)*TMath::Cos(etaA)*X[j])+(TMath::Cos(psiA)*Y[j])+(-TMath::Sin(psiA)*TMath::Sin(etaA)*Z[j]);
      Z1[j]=(TMath::Sin(etaA)*X[j])+(TMath::Cos(etaA)*Z[j]);
    }
  
}



// this function shifts one of the nucleus
void mc_glau::shift_nucleus(double* X1, double* Y1, double* Z1,int A, double b,
                     double zhi,double* X2, double* Y2, double* Z2 )
{
  for(int j=0;j<A;j++){ X2[j]=X1[j]+((b)*TMath::Cos(zhi)); Y2[j]=Y1[j]+((b)*TMath::Sin(zhi));}
}



// this function calculates N_{part} & N_{coll}
void mc_glau::calculate_npart_ncoll(double* vxA,double* vyA,double* vxB,double* vyB, int &Npart, 
       int &Ncoll, double* Npart_x, double* Npart_y, double* Ncoll_x, double* Ncoll_y)
{
  
   Ncoll=0;
   Npart=0;

  double occA[1000];double occB[1000];         //flag during calc of Npart
  //double Ncoll_x[2000]; double Ncoll_y[2000];  // x & y co-ordinate of binary collision sources
  //double Npart_x[1000]; double Npart_y[1000];  // x & y co-ordinate of participant sources

  
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
	      Ncoll_x[Ncoll]=(vxA[i]+vxB[j])/2;
	      Ncoll_y[Ncoll]=(vyA[i]+vyB[j])/2;
	      Ncoll=Ncoll+1;
	      
	      if(occA[i]==0)
		{ 
		  occA[i]=1;Npart_x[Npart]=vxA[i]; Npart_y[Npart]=vyA[i];Npart=Npart+1;
		} 
	      
              
	      if(occB[j]==0)
		{
		  occB[j]=1;Npart_x[Npart]=vxB[j]; Npart_y[Npart]=vyB[j];Npart=Npart+1;
		}
	      
	      
	    }                                                           
	}                                                          
    }                                                         
      // [Info]  shifting the energy distributions center to (0,0,0)   
      double xref1=0.0;
      double yref1=0.0;
      double wref1=0.0;
      for(int k=0;k<Npart;k++)
	{ 
	  xref1=xref1+(Npart_x[k]*(0.5*npp*X_hard));
	  yref1=yref1+(Npart_y[k]*(0.5*npp*X_hard));
	  wref1=wref1+(0.5*npp*X_hard);
	}
      
      double xref2=0.0;
      double yref2=0.0;
      double wref2=0.0;
      for(int k=0;k<Ncoll;k++)
	{ 
	  xref2=xref2+(Ncoll_x[k]*(npp*(1-X_hard)));
	  yref2=yref2+(Ncoll_y[k]*(npp*(1-X_hard)));
	  wref2=wref2+(npp*(1-X_hard));
	}
      
      double xAverage=((xref1+xref2)/(wref1+wref2));
      double yAverage=((yref1+yref2)/(wref1+wref2));
      
      // cout<<xAverage<<"  "<<yAverage<<"\n";
      
      for(int k=0;k<Npart;k++)
	{ 
          Npart_x[k] = Npart_x[k]-xAverage;
          Npart_y[k] = Npart_y[k]-yAverage;
	}
     for(int k=0;k<Ncoll;k++)
	{ 
          Ncoll_x[k] = Ncoll_x[k]-xAverage;
          Ncoll_y[k] = Ncoll_y[k]-yAverage;
	}

}




// This function calculates eccentricity and participant plane angle
void mc_glau::calculate_eccentricity(int aN_part,int aN_coll,double *Npart_x,
                     double *Npart_y,double *Ncoll_x,double *Ncoll_y)
{

  for(int i=0; i<8; i++){ eccentricity[i] = 0.0; PhiN[i] = 0.0;}

	   int Total_Nch=aN_part+aN_coll;
           double Nch_r[Total_Nch];
           double Nch_phi[Total_Nch];
           double Nch_value[Total_Nch];

           double Npart_r[aN_part];
           double Npart_phi[aN_part];

           double Ncoll_r[aN_coll];
           double Ncoll_phi[aN_coll];

           double Ra,Theta;

	   for(int k=0;k<aN_part;k++)
	     { 
	       Ra= TMath::Sqrt( TMath::Power(Npart_x[k],2)+TMath::Power( Npart_y[k],2) );
	       Theta =(TMath::ATan2(Npart_y[k],Npart_x[k]));
	       Npart_r[k]=Ra;
	       Npart_phi[k]=Theta;
	     }
	   for(int k=0;k<aN_coll;k++)
	     {      
	       Ra=TMath::Sqrt( TMath::Power(Ncoll_x[k],2)+ TMath::Power(Ncoll_y[k],2));
	       Theta=(TMath::ATan2(Ncoll_y[k],Ncoll_x[k])); 
	       Ncoll_r[k]=Ra;
	       Ncoll_phi[k]=Theta;
	     }

	  if(aN_part != 0 && aN_coll !=0)
           {
	   
	   for(int k=0;k<aN_part;k++)
	     { //START3
	       Nch_r[k]=Npart_r[k];
	       Nch_phi[k]=Npart_phi[k];
	       Nch_value[k]=((0.5)*(npp)*(X_hard));
	     } //END3
	   
	   for(int k=0;k<aN_coll;k++)
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
	       
	       for(int k=0; k<Total_Nch;k++)
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


