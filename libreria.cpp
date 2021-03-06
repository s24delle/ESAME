#include "libreria.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <TGraphErrors.h>
#include <fstream>

using namespace std;



void readfile(string nomefile, int& n, vector<string>& at,string& a1,string& a2, double& b1, double& b2, vector<double>& X, vector<double>& Y, vector<double>& Z){

  ifstream file(nomefile);
  if(!file.good()){
    cout << "Impossibile aprire il file" << endl;
  }
  file >>n>> a1 >> a2 >> b1 >> b2;
  double x,y,z;
  string Ar;
  for(int i=0; i<n; i++){
    file >> Ar >> x >> y >> z;
    at.push_back(Ar);
    X.push_back(x*pow(10,-10));    //converto i dati in sistema internazionale
    Y.push_back(y*pow(10,-10));
    Z.push_back(z*pow(10,-10));
  }
  file.close();
};


void writefile (ofstream& file, int& n, vector<string>& at, string& a1, string& a2, double& b1, double& b2, vector<double>& X, vector<double>& Y, vector<double>& Z){
  double eps = 119*1.38065e-23;
  file << n << endl;
  file.precision(15);
  file << a1 << '\t'  <<    a2 << '\t' <<  b1*6.241506e18 << '\t'   << b2/eps<< endl;
  file.precision(12);
								   
  for(int i=0; i<n; i++){
    file << at[i] << '\t' << X[i]*pow(10,10) << '\t' <<  Y[i]*pow(10,10) << '\t' << Z[i]*pow(10,10) << endl;     //ho riconvertito i dati nel file di output in armstrong
  }
};


double energia_pot(int& n,double& epsilon, double& sigma,double& kb, vector<double>& X, vector<double>& Y, vector<double>& Z){
  double U=0;
  double r_ij;

  for(int i=0; i<n-1; i++){
    for(int j=i+1; j<n; j++){
      r_ij = sqrt( pow(X[i] - X[j],2) +  pow(Y[i] - Y[j],2) +  pow(Z[i] - Z[j],2)  );
      U = U + 4*epsilon * ( pow(sigma / r_ij ,12) - pow(sigma / r_ij,6)  );
    }
  }
  return U;
};


double energia_cin(int& n, double& m, vector<double>& VX, vector<double>& VY, vector<double>& VZ){
  double K=0;
  for(int i=0; i<n; i++){
    K += 0.5*m*pow(VX[i],2) + 0.5*m*pow(VY[i],2) + 0.5*m*pow(VZ[i],2);
  }
  return K;
};



void forza(int& n, double& epsilon, double& sigma, vector<double>& X, vector<double>& Y, vector<double>& Z, vector<double> &fx, vector<double> &fy, vector<double> &fz){
  double r_ij;
  fx.clear();
  fy.clear();
  fz.clear();
  vector<double> fx_temp(n,0.0);
  vector<double> fy_temp(n,0.0);
  vector<double> fz_temp(n,0.0);
  double a,b,c;
  for(int i=0;i<n-1;i++){
    a=0; b=0; c=0;
    for(int j=i+1; j<n; j++){
      r_ij = sqrt( pow(X[i] - X[j],2) +  pow(Y[i] - Y[j],2) +  pow(Z[i] - Z[j],2)  );
      a = (4*epsilon*(- (6*pow(sigma,6))/(pow(r_ij,6)) + (12*pow(sigma,12))/(pow(r_ij,12)) )*(X[i]-X[j])/(r_ij*r_ij) );
      b = (4*epsilon*(- (6*pow(sigma,6))/(pow(r_ij,6)) + (12*pow(sigma,12))/(pow(r_ij,12)) )*(Y[i]-Y[j])/(r_ij*r_ij) );
      c = (4*epsilon*( -(6*pow(sigma,6))/(pow(r_ij,6)) + (12*pow(sigma,12))/(pow(r_ij,12)) )*(Z[i]-Z[j])/(r_ij*r_ij) );

      fx_temp[i] += a;
      fx_temp[j] += -a;
      fy_temp[i] += b;
      fy_temp[j] += -b;
      fz_temp[i] += c;
      fz_temp[j] += -c;
    }
  }
  fx = fx_temp;
  fy = fy_temp;
  fz = fz_temp;
};





void DINAMICA(ofstream& file, int& N, int& passi, double& tau, double& T, double& eta, double& epsilon, double& sigma,double& kb,double& m, vector<double>& rx, vector<double>& ry, vector<double>& rz,vector<double>& vx, vector<double>& vy, vector<double>& vz, vector<double> &fx, vector<double> &fy, vector<double> &fz, vector<string>& at, string& a1, string& a2, double& b1, double& b2,TGraph& g1,TGraph& g2,TGraph& g3,TGraph& g4){


  TRandom3 rnd;
  double t=0;
  double time =0;
  double dt = tau;
  double tmax = passi*dt;
  double ax, ay, az;
  vector<double> ax_temp,  ay_temp,  az_temp;
  int n_passi=0;
  int n_campioni=0;
  int num = 0;
  double E,K;
  double T_temp;
  double r1, r2, r3;
  double T_med = 0;
  double E_med, K_med, U_med;
  double x_cm, y_cm, z_cm;
  double U=0;


 
  forza(N,epsilon, sigma, rx,ry,rz,fx,fy,fz);      //prima del loop altrimenti le calcolo troppe volte
  // ofstream file1("/home/delle/Computazionale/ESAME/Dati/valori_medi.dat");
  // ofstream file2("/home/delle/Computazionale/ESAME/Dati/energia.dat");
  // file1.precision(12);
  //file2.precision(12);

 
    t = 0;
    T_med =0;
    K_med = 0;
    E_med = 0;
    U_med = 0;
    n_campioni = 0;
    num = 0;

    
    while(t<tmax){
      E=0;
      K=0;
      ax_temp.clear();
      ay_temp.clear();
      az_temp.clear();
      for(int i=0; i<N; i++){
	ax = fx[i]/m;     // questa ?? a(t) calcolata con f(t) per la i esima particella
	ax_temp.push_back(ax);
	ay = fy[i]/m;     // questa ?? a(t) calcolata con f(t) per la i esima particella
	ay_temp.push_back(ay);
	az = fz[i]/m;     // questa ?? a(t) calcolata con f(t) per la i esima particella
	az_temp.push_back(az);
	rx[i] = rx[i] + vx[i]*dt + 0.5 * ax*dt*dt;    //calcolo nuove posizioni
	ry[i] = ry[i] + vy[i]*dt + 0.5 * ay*dt*dt;
	rz[i] = rz[i] + vz[i]*dt + 0.5 * az*dt*dt;
      }

      if(n_passi%passi==0){
	x_cm = 0;
	y_cm = 0;
	z_cm = 0;
	for(int i=0; i<N; i++){
	  x_cm += rx[i];
	  y_cm += ry[i];
	  z_cm += rz[i];
	}
	x_cm = x_cm/N;
	y_cm = y_cm/N;
	z_cm = z_cm/N;
	for(int i=0; i<N; i++){
	  rx[i] = rx[i] - x_cm;   
	  ry[i] = ry[i] - y_cm;
	  rz[i] = rz[i] - z_cm;
	}
	writefile(file,N,at,a1,a2,b1,b2,rx,ry,rz);
	for(int i=0; i<N; i++){
	  rx[i] = rx[i] + x_cm;   
	  ry[i] = ry[i] + y_cm;
	  rz[i] = rz[i] + z_cm;
	}
      }	    
      forza(N,epsilon, sigma, rx,ry,rz,fx,fy,fz); 
    
      for(int i=0; i<N; i++){
     
	ax = fx[i]/m ;   //sono le nuove accelerazioni
	ay = fy[i]/m ;
	az = fz[i]/m ;

	r1 = rnd.Rndm();
	if(r1 > eta*dt){
	  vx[i] = vx[i] + 0.5*( ax_temp[i] + ax )*dt;      //calcolo nuove velocit??
	}
	if(r1 <= eta*dt){
	  vx[i] = sqrt(kb*T/m)*rnd.Gaus(0,1);
	}
      
	r2 = rnd.Rndm();
	if(r2 > eta*dt){
	  vy[i] = vy[i] + 0.5*( ay_temp[i] + ay )*dt;
	}
	if(r2 <= eta*dt){
	  vy[i] =  sqrt(kb*T/m)*rnd.Gaus(0,1);
	}
      
	r3 = rnd.Rndm();
	if(r3 > eta*dt){
	  vz[i] = vz[i] + 0.5*( az_temp[i] + az )*dt;
	}
	if(r3 <= eta*dt){
	  vz[i] =  sqrt(kb*T/m)*rnd.Gaus(0,1);
	} 
      }
  
      


      n_campioni += 1;         //conta il numero di passi da 1 a 100k per una sola temperatura
      n_passi += 1;            //?? il numero totale di passi di tutti i cicli
      t += dt;
      //  if(n_campioni>=5000){
      K = energia_cin(N,m,vx,vy,vz);
      //	K_med += K;
	U = energia_pot(N, epsilon, sigma, kb, rx, ry, rz);
	//	U_med += U;
	E = K+U;

	//	E_med += E;
	T_temp = 2./3 * K /(N*kb);
	//	T_med += T_temp;
	num += 1;
	if(n_campioni%10==0){
	  //	  file2 << time << '\t' << E << endl;
	}
	//  }
	//  g1.SetPoint(n_passi,time,E);
	// g2.SetPoint(n_passi,time,U);
	// g3.SetPoint(n_passi,time,K);
	// g4.SetPoint(n_passi,time,T_temp);
      
      time += dt;
      // file1 << T <<  '\t' << T_med <<  '\t' << K_med <<  '\t' << U_med <<   '\t'  << E_med << endl;
    }
    n_passi = n_passi -1;
    //   K_med = K_med/num;
    //   U_med = U_med/num;
    //  E_med = E_med/num;
    //  T_med = T_med/num;


  
    // file1.close();
    // file2.close();


};



void QUENCING(ofstream& file,int& N, int& passi, int& num_quencing, double& tau,double& T, double& eta, double& epsilon, double& sigma,double& kb,double& m, vector<double>& rx, vector<double>& ry, vector<double>& rz,vector<double>& vx, vector<double>& vy, vector<double>& vz, vector<double> &fx, vector<double> &fy, vector<double> &fz, vector<string>& at, string& a1, string& a2, double& b1, double& b2,TGraph& g1,TGraph& g2,TGraph& g3,TGraph& g4, double& E_fin){



 TRandom3 rnd;
  double t=0;
  double dt = tau;
  double tmax = passi*dt;
  double ax, ay, az;
  vector<double> ax_temp,  ay_temp,  az_temp;
  int n_passi=0;
  int n_campioni=0;
  double E,K;
  double T_temp;
  double E_med;
  double U=0;
  double dist = 0;

  // string name="/home/delle/Computazionale/ESAME/Dati/output_quencing" + std::to_string(num_quencing) + ".xyz"; 
  //   ofstream file1(name);
  forza(N,epsilon, sigma, rx,ry,rz,fx,fy,fz);      //prima del loop altrimenti le calcolo troppe volte


  while(t<tmax){

    E=0;
    K=0;
    ax_temp.clear();
    ay_temp.clear();
    az_temp.clear();
    for(int i=0; i<N; i++){
      ax = fx[i]/m;     // questa ?? a(t) calcolata con f(t) per la i esima particella
      ax_temp.push_back(ax);
      ay = fy[i]/m;     // questa ?? a(t) calcolata con f(t) per la i esima particella
      ay_temp.push_back(ay);
      az = fz[i]/m;     // questa ?? a(t) calcolata con f(t) per la i esima particella
      az_temp.push_back(az);
      rx[i] = rx[i] + vx[i]*dt + 0.5 * ax*dt*dt;    //calcolo nuove posizioni
      ry[i] = ry[i] + vy[i]*dt + 0.5 * ay*dt*dt;
      rz[i] = rz[i] + vz[i]*dt + 0.5 * az*dt*dt;
    }




    
    forza(N,epsilon, sigma, rx,ry,rz,fx,fy,fz);
    U = energia_pot(N,epsilon, sigma, kb, rx, ry, rz);
	
    for(int i=0; i<N; i++){
      
      if(fx[i]*vx[i] + fy[i]*vy[i] + fz[i]*vz[i] >= 0){
	ax = fx[i]/m ;   //sono le nuove accelerazioni
	ay = fy[i]/m ;
	az = fz[i]/m ; 
	vx[i] = vx[i] + 0.5*( ax_temp[i] + ax )*dt;       //calcolo nuove velocit??
	vy[i] = vy[i] + 0.5*( ay_temp[i] + ay )*dt;
	vz[i] = vz[i] + 0.5*( az_temp[i] + az )*dt;
      }

      if(fx[i]*vx[i] + fy[i]*vy[i] + fz[i]*vz[i] < 0){
	vx[i] = 0;
	vy[i] = 0;
	vz[i] = 0;
      }
    }
    
  
    K = energia_cin(N,m, vx,vy,vz);
    E = (K+U);
    T_temp = 2./3 * K /(N*kb);
    //  E_med += E;
    //  g1.SetPoint(n_campioni,t,E);
    //  g2.SetPoint(n_campioni,t,U);
    //  g3.SetPoint(n_campioni,t,K);
    //  g4.SetPoint(n_campioni,t,T_temp);
    n_campioni += 1;

    
    //  if(n_passi%100==0){
    //    writefile(file,N,at,a1, a2, b1, b2,rx,ry,rz);
 
    //  }
    
    n_passi += 1;
    t += dt;
  }
  n_passi -= 1;
  //  writefile(file1,N,at,a1,a2, E, E,rx,ry,rz);
  // file1.close();
  E_fin = E;

};
