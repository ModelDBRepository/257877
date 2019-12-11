#include <stdio.h>
#include <stdlib.h> 
#include <math.h> 
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>

using namespace std;


//---------Define the parameters---------------
#define   nKC    15000    // number of KC cells
#define   nAlpha 200  // number of Alpha cells
#define   iters  50000 //number of iterations
#define   g_ampa_AlphaKC 0.1//AMPA from KC to Alpha
//----------------------------------------------



//-- Regular spiking neuron cell---------------------------------------
class RS{
  double xp, xpp, mu, sigma_n, beta_n;
  double sigma_dc, beta_dc;
public:
  double x, y, alpha, sigma, sigma_e, beta_e, Idc, S_CX_DEND;
  int spike;

  RS(){
    mu=0.0005;
    spike=0;
    alpha=3.65;
    sigma=0.06;//0.17; //6.00E-2; //0.12; //6.00E-2;
    sigma_e=1.0;
    sigma_dc=1.0;
    beta_e=0.03; //0.133;
    beta_dc=0.133;
    Idc = 0.02;
    S_CX_DEND = 165.0e-6;
    xpp=-1+sigma;
    xp=-1+sigma;
    x=-1+sigma;
    y= x - alpha/(1-x);
  }

  void init(){
    mu=0.0005;
    spike=0;
    alpha=3.65;
    sigma=0.06;//0.17; //6.00E-2; //0.12; //6.00E-2;
    sigma_e=1.0;
    sigma_dc=1.0;
    beta_e=0.025; //0.133;
    beta_dc=0.133;
    Idc = 0.02;
    S_CX_DEND = 165.0e-6;
    xpp=-1+sigma;
    xp=-1+sigma;
    x=-1+sigma;
    y= x - alpha/(1-x);
  }
  void calc(double);
};
void RS::calc(double I){
  beta_n = beta_e * I + beta_dc * Idc;
  sigma_n = sigma_e * I + sigma_dc * Idc;

  //     if (beta_n < -0.0001) beta_n = -0.0001;
  if (beta_n < -1.0) beta_n = -1.0;
  if (beta_n > 1.0) beta_n = 1.0;
  if(xp <= 0.0) {
    x = alpha / (1.0 - xp) + y +beta_n;
    spike = 0;
  }
  else{
    if(xp <= alpha + y +beta_n && xpp <= 0.0) {
      x = alpha + y + beta_n;
      spike = 1;
    }
    else {

      x = -1;
      spike = 0;
    }
  }
  y = y - mu* (xp +1.0) + mu * sigma + mu * sigma_n;
  xpp = xp;
  xp = x;
  //     y=1;
}

//------------------------------------------------------------------

//------------------------------------------------------------------
class AMPAmap1 {
  static double E_AMPA;
  static double gamma;

public:
  double I, g;
  AMPAmap1() {
    I=0;
    g=0;
  }
  void calc(double, double, int);
};
double AMPAmap1::E_AMPA = -0.0, AMPAmap1::gamma = 0.4; //0.6;
void AMPAmap1::calc(double g_AMPA, double x_post, int spike) {

  if(spike > 0.1){
    g = gamma * g + g_AMPA;
  }
  else{
    g = gamma * g;
  }
  I = - g * (x_post - E_AMPA);
}


//------------------------------------------------------------------
//AMPA synapse with STDP
class AMPAmapS {
  static double E_AMPA;
  static double gamma;

public:
  double I, g,x,y,g_pre,g_post,a_pre,a_post,A_plus,A_minus,w,wmin,wmax,n_minus,n_plus;
  AMPAmapS() {
    I=0;
    g=0;
    x = 0.0;
    y = 0.0;
    g_pre = 0.95; //decay rate of x variable
    g_post = 0.95; //decay rate of y variable
    a_pre = 0.04; //in/decrement of x variable
    a_post = 0.04; //in/decrement of y variable
    A_minus = 1.0; //initial value of weight increase function
    A_plus = 1.0; //initial value of weight increase function
    w = 1.0; //initial weight
    wmax = 2.0; //max weight
    wmin = 0.0;//min weight
    n_minus = 1; //factor that multiplies max weight
    n_plus = 1;

  }
  void calc(double, double, int, int);
};
double AMPAmapS::E_AMPA = -0.0, AMPAmapS::gamma = 0.4; //0.6;
void AMPAmapS::calc(double g_AMPA, double x_post, int spike_pre, int spike_post) {

  //calculate the values of x and y variables
  if (spike_pre > 0.1){
    x = g_pre*x + a_pre;
    w = w - A_minus*y;
    
    if (w < wmin){
      w = wmin;
    }

    if (w > wmax){
      w = wmax;
    }

    g = gamma * g + g_AMPA*w;
  }

  if (spike_post > 0.1){
    y = g_post*y + a_post;
    w = w + A_plus*x;

    if (w < wmin){
      w = wmin;
    }

    if (w > wmax){
      w = wmax;
    }

    g = gamma*g;
  }

  if ((spike_pre <= 0.1) && (spike_post <= 0.1)){
    x = g_pre*x;
    y = g_post*y;
    g = gamma*g;
  }

  I = - g*(x_post - E_AMPA);
}

//-----------------------------------------------------------------





//----external functions--------------------------
void fun(int);

//--------external variables---------------------

AMPAmapS           *ampa_AlphaKC[nAlpha];
int              C_AlphaKC[nAlpha][nKC], spikesKC[nKC];
double           R;
FILE             *fspikesAlpha;


//-------Neurons----------
 
RS            Alpha[nAlpha]; // Alpha cells     



main(int argc, char **argv)
{
  int  t,ii,jj,xx,yy,zz,i,j;
  double kk, mm;
  ifstream fspikesKC("spikesKC2");
  string line;



  //***************initialization********//
  //---initialize connectivity matrix-----

  //1. KC to Alpha connectivity
  for (ii = 0; ii < nAlpha; ++ii){
    for (jj = 0; jj < nKC; ++jj){
      C_AlphaKC[ii][jj] = 0;
    }
  }
  //--- end initialize connectivity matrix-------------


  //-----read the connectivity matrix-----------------
  //randomly generated
  srand(1);
  for(i=0; i<nAlpha; i++)
    for(j=0; j<nKC; j++){
      R = 1.0 * rand()/(RAND_MAX + 1.0);
      
      if(R < 2) {
	C_AlphaKC[i][j]=1;  // all KCs connected to all Alphas initially
      }
    }
  //-----end read the connectivity matrix-----------------


  //----initialize synapses -------------------------
  ii = 0;

  //1. KC to Alpha --> AMPA 
  for(ii = 0; ii < nAlpha; ++ii){    
    ampa_AlphaKC[ii] = new AMPAmapS[nKC];
  }
  //---end initialize synapses--------------

  //**************end initialization*****************//

  //-----introduce some random variability--------

  srand(2);
  //1. variability in Alphas
  for(ii = 0; ii < nAlpha; ++ii){
    R = 2.0 * rand()/(RAND_MAX + 1.0) - 1.0;
    Alpha[ii].sigma = Alpha[ii].sigma + R*0.0005;

  }
  //--------------------------------------------





  //--------Start iteration--------------------

  t = 0;  

  fspikesAlpha = fopen("spikesAlphaSTDP","w");

  while(t < iters){
    

    //--initialize KC spike vector--
    
    for (ii = 0; ii < nKC; ++ii){
      spikesKC[ii] = 0;

    }
    //end initialize spike vector---
  
    // read KC spike vector
    getline(fspikesKC,line);
    int value;
    istringstream iss(line);
    while (iss >> value){
      spikesKC[value] = 1;
    }

    //evaluate rhs of all diffeqns
    fun(t);
    
    //     //write Alpha spikes to file
    //     for (ii = 0; ii < nAlpha; ++ii){
    //       fprintf(fspikesAlpha,"%lf \t",Alpha[ii].x);
    //     }
    //     fprintf(fspikesAlpha,"\n");


    for (ii = 0; ii < nAlpha; ++ii){
      if(Alpha[ii].spike == 1){
	fprintf(fspikesAlpha,"%d %d \n", t, ii);

      }
    }


    // printf("%d ",t); 
   ++t;
  }

  fclose(fspikesAlpha);
 

} //end main

  //-----------function to evaluate the Rhs of all the maps---


void fun(int t){
  int i, j, k, kmax;  
  double IsynAlphaKC[nAlpha];


  //------------ calculating synaptic currents-----------------------------------------


  //1. calculating AMPA type current from PN to KC-------------------------------------

  for(i = 0; i < nAlpha; ++i){
    k = 0;
    for(j = 0; j < nKC; ++j){
      if(C_AlphaKC[i][j] > 0){
	ampa_AlphaKC[i][k].calc(g_ampa_AlphaKC/Alpha[i].S_CX_DEND, Alpha[i].x, spikesKC[j],Alpha[i].spike);
	++k;
          
      } 
    } 


    kmax = k;
    IsynAlphaKC[i] = 0.0;
    for (k = 0; k < kmax; ++k) {

      if (kmax == 0) { } // no connections

      else{
	IsynAlphaKC[i] = IsynAlphaKC[i] + ampa_AlphaKC[i][k].I/kmax;
      } 

    }
  } 
  //---------end AMPA type current----------------------------------------------------------



  //----calculate the state of all Alpha at the next iteration based on synaptic input-------
    
  for(i = 0; i < nAlpha; ++i){
    Alpha[i].calc(IsynAlphaKC[i]);      
  }
  //---end calculate the state of all Alphas at the next iteration based on synaptic input-----

} //end fun
