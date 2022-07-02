#include <iostream>
#include <cmath>
#include <vector>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <fstream>
#include <TH1D.h>
#include "libreria.h"


using namespace std;

namespace dati{
  double kb = 1.38065e-23;   //in Joule/K
  double sigma = 3.4e-10;    //nel SI
  double epsilon = 119*kb;  //epsilon in Joule
  string a1,a2;
  double b1,b2;
  double T=300;    //in kelvin
  double T_quencing=10;
  double T_bh = 100;    //temperatura basin hopkins
  double tau = 5e-15;     //5 femptosecondi
  double m = 6.68e-26;    //è la massa in chili
  double eta = 5e11;

  //energia Joule -> ev devo moltiplicare per 6.241506e18
  //energia Joule -> epsilon=1 devo dividere per epsilon
}





//MAIN

int main(){


  TApplication app("app",0,NULL);




  vector<double> rx;
  vector<double> ry;
  vector<double> rz;
  vector<string> at;

  vector<double> vx;
  vector<double> vy;
  vector<double> vz;

  vector<double> fx;
  vector<double> fy;
  vector<double> fz;
  
  int N;
  
   readfile("initial_position.xyz", N, at, dati::a1, dati::a2, dati::b1, dati::b2 ,rx, ry,rz);
   // readfile("ar38_to.xyz", N, at, dati::a1, dati::a2, dati::b1, dati::b2 ,rx, ry,rz);
  ofstream file("/home/delle/Computazionale/ESAME/Dati1/output.xyz");
  ofstream file_quencing("/home/delle/Computazionale/ESAME/Dati1/output_quencing.xyz");
  writefile(file,N,at,dati::a1,dati::a2,dati::b1,dati::b2,rx,ry,rz);
  cout.precision(15);



  //GENERATORE VELOCITA INIZIALI
  
  TRandom3 rnd;
  rnd.SetSeed(time(0));
  TRandom3 rnd_uni;
  rnd_uni.SetSeed(time(0));
  vector<double> vx_star;
  vector<double> vy_star;
  vector<double> vz_star;
  double v_cm_x, v_cm_y, v_cm_z;
  for(int i=0; i<N; i++){
    
    vx_star.push_back(  sqrt(dati::kb*dati::T/dati::m)*rnd.Gaus(0,1));    //velocità casuali estratte da maxwelliana
    vy_star.push_back(  sqrt(dati::kb*dati::T/dati::m)*rnd.Gaus(0,1));
    vz_star.push_back(  sqrt(dati::kb*dati::T/dati::m)*rnd.Gaus(0,1));
    v_cm_x += vx_star[i];
    v_cm_y += vy_star[i];
    v_cm_z += vz_star[i];
  }
  v_cm_x = v_cm_x/N;
  v_cm_y = v_cm_y/N;
  v_cm_z = v_cm_z/N;


  for(int i=0; i<N; i++){
    vx.push_back(vx_star[i] - v_cm_x  );    //queste sono le velocità traslate
    vy.push_back(vy_star[i] - v_cm_y  );    //quelle che lui chiama con v barra
    vz.push_back(vz_star[i] - v_cm_z  );
  }

  double K_bar;
 
  K_bar = energia_cin(N,dati::m,vx,vy,vz);
 

  double K_0;
  K_0 = 3./2*N*dati::kb*dati::T ;
 

  for (int i=0; i<N; i++) {
    vx[i] = vx[i]*sqrt(K_0/K_bar);
    vy[i] = vy[i]*sqrt(K_0/K_bar);
    vz[i] = vz[i]*sqrt(K_0/K_bar);
  }


  

  //INIZIALIZZAZIONE DI CONFIGURAZIONE INIZIALE RANDOM

  
  for(int i=0; i<N; i++){
    double a = rnd.Rndm()*1;
    double b = rnd.Rndm()*1;
    double c = rnd.Rndm()*1;
    rx[i] += a*pow(10,-10);
    ry[i] += b*pow(10,-10);
    rz[i] += c*pow(10,-10);
  }
  writefile(file,N,at,dati::a1,dati::a2,dati::b1,dati::b2,rx,ry,rz);
  
  
  
  /*
  for(int i=0;i<N; i++){
    double dist = 0;
    bool test = false;
    int j=0;
    rx[i] = 15*rnd.Rndm()*pow(10,-10);
    ry[i] = 15*rnd.Rndm()*pow(10,-10);
    rz[i] = 15*rnd.Rndm()*pow(10,-10);
    dist = sqrt(pow(rx[i],2) + pow(ry[i],2) + pow(rz[i],2));
      while(test==false){
	while(j<N){
	  dist = sqrt(pow(rx[i]-rx[j],2) + pow(ry[i]-ry[j],2) + pow(rz[i]-rz[j],2));
	  if(dist < 4e-10 && j!=i){
	    rx[i] = 15*rnd.Rndm()*pow(10,-10);
	    ry[i] = 15*rnd.Rndm()*pow(10,-10);
	    rz[i] = 15*rnd.Rndm()*pow(10,-10);
	    j = 0;
	  }

	  if(j==N-1){
	    test = true;
	  }
	  j++;
	}
      }
    
  }
  writefile(file,N,at,dati::a1,dati::a2,dati::b1,dati::b2,rx,ry,rz);
  */
  


  
  //DINAMICA E QUENCING
  
  int passi = 300;
  int passi_quencing = 10000;    // numero passi che fa il quencing
  int N_finale = 1000;
  int num_quencing = 1;          // numero di volte che ho fatto il quencing per il file di output
  ofstream file_energie1("/home/delle/Computazionale/ESAME/Dati1/energie_prob.dat");
  ofstream file_energie2("/home/delle/Computazionale/ESAME/Dati1/energie_menoprob.dat");

  TGraph g1;
  TGraph g2, g3, g4;

 
  double Temp = dati::T;
 vector<double> rx_temp;
 vector<double> ry_temp;
 vector<double> rz_temp;
 double E_temp;
   double E_fin;
   double deltaE;
   int spintarella = 0;

 
 QUENCING(file_quencing,N,passi_quencing,num_quencing, dati::tau, dati::T_quencing,dati::eta, dati::epsilon,dati::sigma, dati::kb,dati::m,rx,ry,rz,vx,vy,vz,fx,fy,fz,at,dati::a1,dati::a2,dati::b1,dati::b2,g1,g2,g3,g4,E_fin);


 
 
  num_quencing += 1;

  double z;
  double P;

  //  DINAMICA(file,N,passi, dati::tau, dati::T,dati::eta, dati::epsilon,dati::sigma, dati::kb,dati::m,rx,ry,rz,vx,vy,vz,fx,fy,fz,at,dati::a1,dati::a2,dati::b1,dati::b2,g1,g2,g3,g4);



  
  while(num_quencing<N_finale){
    ofstream file_minimi("/home/delle/Computazionale/ESAME/Dati1/minimo"+ std::to_string(num_quencing) + ".xyz"); 
    file_minimi.precision(15);
    file_energie1.precision(15);
    file_energie2.precision(15);

    rx_temp = rx;
    ry_temp = ry;
    rz_temp = rz;
    E_temp = E_fin;
    
    DINAMICA(file,N,passi, dati::tau, Temp,dati::eta, dati::epsilon,dati::sigma, dati::kb,dati::m,rx,ry,rz,vx,vy,vz,fx,fy,fz,at,dati::a1,dati::a2,dati::b1,dati::b2,g1,g2,g3,g4);

    
  QUENCING(file_quencing,N,passi_quencing,num_quencing, dati::tau, dati::T_quencing,dati::eta, dati::epsilon,dati::sigma, dati::kb,dati::m,rx,ry,rz,vx,vy,vz,fx,fy,fz,at,dati::a1,dati::a2,dati::b1,dati::b2,g1,g2,g3,g4,E_fin);

 

  P = exp(-(E_fin - E_temp)/(dati::kb*dati::T_bh));
  z = rnd.Rndm();
  cout << P << "       ------------->        " << 1.*num_quencing/N_finale*100 << " % " << endl;
  if(z<=P){
    deltaE = abs(E_temp - E_fin);
    if(deltaE <= 1e-7){spintarella +=1;}
    if(deltaE > 1e-7){
      spintarella = 0;
    }
    E_temp = E_fin;
    writefile(file_minimi,N,at,dati::a1,dati::a2,E_fin, E_fin,rx,ry,rz);
    file_energie1 << num_quencing << '\t' << E_fin*6.241506e18 << "  eV "  << endl;
  
  }
  if(z>P){
    spintarella  = 0;
    file_energie2 << num_quencing << '\t' << E_fin*6.241506e18 << "  eV "  << endl;
    writefile(file_minimi,N,at,dati::a1,dati::a2, E_fin, E_fin,rx,ry,rz);
    rx = rx_temp;
    ry = ry_temp;
    rz = rz_temp;
    E_fin = E_temp;
  }

  if(spintarella >= 6){Temp=350;}
  if(spintarella < 6){Temp = dati::T;}
 
  
  num_quencing +=1;

  file_minimi.close();

  

  }
  
  
  file_energie1.close();
  file_energie2.close();
 
  
  //GRAFICI CONTROLLO


  /*
    TCanvas c("","",800,600);
  c.Divide(2,2);

  c.cd(1);
  // TCanvas c1;
  g1.SetTitle("Energia");
  g1.Draw("ap");
  c.cd(2);
  //TCanvas c2;
  g2.SetTitle("Energia potenziale");
  g2.Draw("ap");
  c.cd(3);
  // TCanvas c3;
  g3.SetTitle("Energia cinetica");
  g3.Draw("ap");
  c.cd(4);
  //TCanvas c4;
  g4.SetTitle("Temperatura");
  g4.Draw("ap");
  */





  

  //  app.Run(true);
  return 0;




}
