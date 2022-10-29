/*
  Objective: To use the secant method to solve the transcendental equation
  for transmission through infinite barriers of alternating potiential and
  non-potienital regions.
*/


#include<iostream>
#include<cmath>
#include<sstream>
#include<fstream>
#include<vector>
// #include<complex>

using namespace std;

typedef vector<double> Vec;
typedef vector<Vec> Mat;

// const double m = .511; //
//const double h_bar = 4.; //
const double v = 13.60; // ev
const double b = 2; // angstroms
const double w = 3; //angstroms
const double a = b+w; // unit lattice length
const double A = 3.84; // h_bar/2m
// const float pahse_E = 1;
// const float beta_E = 1;
const int itr = 10000; // max iterations for secant method
const double tol = 1e-6; // tolerance level for root value

Mat energies_mn_mx(){
  double output1;
  double output2;
  double q;
  double beta;
  Vec func;
  Vec energies;
  Mat all_data;
  double E = 0.0001;

  while(E < v*6){
    if(E < v){
      q = sqrt(E/A);
      beta = sqrt((v-E)/A);
      func.push_back(cosh(beta*b)*cos(q*w) - ((pow(q,2) - pow(beta,2))/(2*q*beta))*sinh(beta*b)*sin(q*w));
    } else if(E >= v){
      q = sqrt(E/A);
      beta = sqrt(abs(v-E)/A);
      func.push_back(cos(beta*b)*cos(q*w) - ((pow(q,2) - pow(beta,2))/(2*q*beta))*sin(beta*b)*sin(q*w));
    }
    energies.push_back(E);
    E = E + 0.1;
    //cout<<energies.begin()<<endl;
  }
  all_data.push_back(func);
  all_data.push_back(energies);
  return all_data;
}

// When k = 0, max of cos(ka)
double func1(double E){
  double output;
  double q;
  double beta;
  q = sqrt(E/A);
  beta = sqrt((v-E)/A);
  output = cosh(beta*b)*cos(q*w) - ((pow(q,2) - pow(beta,2))/(2*q*beta))*sinh(beta*b)*sin(q*w) - 1;
  return output;
}
// When k = pi/2, rmin of cos(ka)
double func2(double E){
  double output;
  double q;
  double beta;
  q = sqrt(E/A);
  beta = sqrt((v-E)/A);
  output = cosh(beta*b)*cos(q*w) - ((pow(q,2) - pow(beta,2))/(2*q*beta))*sinh(beta*b)*sin(q*w) + 1;
  return output;
}

double func3(double E, double K){
  double output;
  double q;
  double beta;
  q = sqrt(E/A);
  beta = sqrt((v-E)/A);
  output = cosh(beta*b)*cos(q*w) - ((pow(q,2) - pow(beta,2))/(2*q*beta))*sinh(beta*b)*sin(q*w) - cos(K*a);
  return output;
}

double func4(double E){
  double output;
  double q;
  double beta;
  q = sqrt(E/A);
  beta = sqrt((v-E)/A);
  output = (cosh(beta*b)*cos(q*w) - ((pow(q,2) - pow(beta,2))/(2*q*beta))*sinh(beta*b)*sin(q*w));
  return output;
}


double root_finder(double x1, double x2, double(*function)(double x)){
  double x3;
  double root;
  for(int i = 0; i < itr; i++){
    x3 = x2;
    if(function(x3) != 0){
      x2 = x2 - ((x2-x1)*function(x2))/(function(x2)-function(x1));
      x1 = x3;
      root = x2;
    }
    if(abs(function(root)) < tol){
      break;
    }
  }
  //if gap region, return -1 to signify no real value
  if(isnan(root)){
    root = -1;
    //cout<<root<<endl;
  }
  // else{
  //   //cout<<root<<endl;
  // }
  cout<<root<<endl;
  return root;
}

double root_finder_k(double x1, double x2, double K, double(*function)(double x,double k)){
  double x3;
  double root;
  for(int i = 0; i < itr; i++){
    x3 = x2;
    if(function(x3,K) != 0){
      x2 = x2 - ((x2-x1)*function(x2,K))/(function(x2,K)-function(x1,K));
      x1 = x3;
      root = x2;
    }
    if(abs(function(root,K)) < tol){
      break;
    }
  }
  //cout<<root<<endl;
  // if gap region, return -1 to signify no real value
  // if(isnan(root)){
  //   root = -1;
  //   //cout<<root<<endl;
  // }
  // else{
  //   //cout<<root<<endl;
  // }
  return root;
}

int main(){
  Mat data;
  double E;
  Vec roots1;
  double dummy_root;
  Vec roots2;
  Mat roots_all;
  int k = 0;
  int m = 0;
  //int h = 0;
  // Find the max and min energies, then from there use the secant method to find the roots
  data = energies_mn_mx();
  ofstream myfile;
  string fileName = "energies.txt";
  myfile.open(fileName);
  for(int i = 0; i < data[0].size();i++){
    for(int j = 0; j < data.size(); j++){
      if(j != data.size()-1){
        myfile<<data[j][i]<<" , ";
      } else {
        myfile<<data[j][i];
      }
    }
    myfile<<endl;
  }
  myfile.close();

  while(E < v){
    E = E + 1;
    dummy_root = root_finder(E-1,E,&func1);
    if(abs(func4(dummy_root) - 1) < 1e-1){
      roots1.push_back(dummy_root);
    }
    dummy_root = root_finder(E-1,E,&func2);
    if(abs(func4(dummy_root) - 1) < 1e-1){
      roots2.push_back(dummy_root);
    }
    // roots1.push_back(root_finder(E-1,E,&func1));
    // roots2.push_back(root_finder(E-1,E,&func2));
    // cout<<"Top: "<<roots1[m]<<"    M:"<<m<<endl;
    // cout<<"Bottom: "<<roots2[m]<<"    M:"<<m<<endl;
    // m = m+1;
  }

  // Vec newRoot1;
  // for(int i=0;i<roots1.size()-1;i++){
  //   if(abs(roots1[i] - roots1[i+1]) > 1e-6){
  //     if(roots1[i] != -1){
  //       newRoot1.push_back(roots1[i]);
  //     } else{
  //       //cout<<"newRoot for loop"<<endl;
  //     }
  //   }
  // }
  //
  // Vec newRoot2;
  // for(int i=0;i<roots2.size()-1;i++){
  //   if(abs(roots2[i] - roots2[i+1]) > 1e-6){
  //     if(roots2[i] != -1){
  //       newRoot2.push_back(roots2[i]);
  //       //cout<<"newRoot2[i]""<<endl;
  //     }
  //   }
  // }

  roots_all.push_back(roots1);
  roots_all.push_back(roots2);

  fileName = "roots_mn_mx.txt";
  myfile.open(fileName);
  for(int i = 0; i < roots_all[0].size();i++){
    for(int j = 0; j < roots_all.size(); j++){
      if(j != roots_all.size()-1){
        myfile<<roots_all[j][i]<<" , ";
      } else {
        myfile<<roots_all[j][i];
      }
    }
    myfile<<endl;
  }
  myfile.close();

  Vec e_spec;
  Mat band;
  Vec k_spec;
  double K;
  int n = 0;
  for(int i=0;i < roots1.size();i++){
    K = - M_PI/a;
    while(K < M_PI/a){
      K = K + 0.01;
      if(n%2 == 0){
        e_spec.push_back(root_finder_k(roots1[n],roots2[n],K,&func3));
        k_spec.push_back(K);
      } else if(n%2 != 0){
          e_spec.push_back(root_finder_k(roots2[n],roots2[n],K,&func3));
          k_spec.push_back(K);
      }
    }
    //cout<<"N: "<<n<<endl;
    n = n+1;
  }
  band.push_back(e_spec);
  band.push_back(k_spec);

  fileName = "Bands.txt";
  myfile.open(fileName);
  for(int i = 0; i < band[0].size();i++){
    for(int j = 0; j < band.size(); j++){
      if(j != band.size()-1){
        myfile<<band[j][i]<<" , ";
      } else {
        myfile<<band[j][i];
      }
    }
    myfile<<endl;
  }


  return 0;
}
