#include<iostream>
#include<cmath>
#include<sstream>
#include<fstream>
#include<vector>
// #include<complex>

using namespace std;

typedef vector<double> Vec;
typedef vector<Vec> Mat;

const double v = 13.60; // ev
const double b = 2; // angstroms
const double w = 3*2; //angstroms
const double a = b+w; // unit lattice length
const double A = 3.84; // h_bar/2m

// These functions are used to find the max and min energies for each dispersion band
double func1(double E){}
double func2(double E){}
// Used to vary k between two energies
double func3(double E, double K)){}
// Maant as a function to plot the left hand side of Eq. Ev to K
Mat energies_mn_mx(){}
// Used to find the roots of max and min energies of a band
double root_finder(double x1, double x2, double(*function)(double x)){}
double root_finder_k(double x1, double x2, double K, double(*function)(double x,double k)){}

int main(){


  return 0;
}


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
    root = 0;
    //cout<<root<<endl;
  }
  else{
    //cout<<root<<endl;
  }
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
  // if gap region, return -1 to signify no real value
  if(isnan(root)){
    root = 0;
    //cout<<root<<endl;
  }
  else{
    //cout<<root<<endl;
  }
  return root;
}

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
