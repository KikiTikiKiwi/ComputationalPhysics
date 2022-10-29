#include<iostream>
#include<cmath>
#include<sstream>
#include<fstream>
#include<vector>
// #include<complex>

using namespace std;

typedef vector<double> Vec;
typedef vector<Vec> Mat;

const double v = 13.3; // ev
const double b = .9; // angstroms
const double w = 2; //angstroms
const double a = b+w; // unit lattice length
const double A = 3.84; // h_bar/2m
const double del_E = 0.01;

const int itr = 1e5; // max iterations for secant method
const double tol = 1e-6;

// These functions are used to find the max and min energies for each dispersion band
double func(double k,double E);
double func1(double E);
double func2(double E);
// Used to vary k between two energies
double func3(double E, double K);
// Maant as a function to plot the left hand side of Eq. Ev to K
Mat energies_mn_mx();
Mat mirror(Vec k, Vec e);
// Used to find the roots of max and min energies of a band
//double root_finder(double x1, double x2, double(*function)(double x));
double root_finder(double x1, double x2, double e, double(*function)(double x,double k));
double bisection_method(double x1, double x2, double(*function)(double x));

int main(){
  Mat data;
  double E;
  Vec roots1;
  Vec roots2;
  Mat roots_all;
  Mat bands;
  int k = 0;

  data = energies_mn_mx(); // Matrix where Energy is col 1 and FUNC is col 2
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
  return 0;
}

double func(double k,double E){
  double output;
  double q;
  double beta;
  q = sqrt(E/A);
  beta = sqrt((v-E)/A);
  output = cosh(beta*b)*cos(q*w) - ((pow(q,2) - pow(beta,2))/(2*q*beta))*sinh(beta*b)*sin(q*w) - cos(k*a);
  return output;
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

double func3(double K,double E){
  double output;
  double q;
  double beta;
  q = sqrt(E/A);
  beta = sqrt((v-E)/A);
  output = cosh(beta*b)*cos(q*w) - ((pow(q,2) - pow(beta,2))/(2*q*beta))*sinh(beta*b)*sin(q*w) - cos(K*a);
  return output;
}

double bisection_method(double x1, double x2, double(*function)(double x)){
  double x3;
  double root;
  for(int i = 0; i < itr; i++){
    x3 = (x2-x1)/2;
    if(function(x3) != 0){
      if(function(x3)*function(x1) < 0){
        x2 = x3;
        root = x2;
      }
      else if(function(x1)*function(x3) > 0){
        x1 = x3;
        root = x1;
      }
    }
    if(abs(function(root)) < tol){
      break;
    }
  }
  return root;
}

double root_finder(double x1, double x2, double e, double(*function)(double a,double b)){
  double x3;
  double root;
  for(int i = 0; i < itr; i++){
    x3 = x2;
    if(function(x3,e) != 0){
      x2 = x2 - ((x2-x1)*function(x2,e))/(function(x2,e)-function(x1,e));
      x1 = x3;
      root = x2;
    }
    if(abs(function(root,e)) < tol){
      break;
    }
  }
  //cout<<root<<endl;
  return root;
}

Mat mirror(Vec k, Vec e, Vec func){
  Vec dummy = k;
  Mat table;
  // dummy = -1*k;
  for(int i = 0; i < k.size();i++){
    k.push_back(-1*dummy[i]);
    e.push_back(e[i]);
    func.push_back(func[i]);
    cout<<k[i]<<endl;
    cout<<e[i]<<endl;
    cout<<func[i]<<endl;
  }
  table.push_back(k);
  table.push_back(e);
  table.push_back(func);

  return table;
}

Mat energies_mn_mx(){
  double output1;
  double output2;
  double q;
  double beta;
  Vec func;
  double func_value;
  Vec energies;
  Mat all_data;
  Vec k;
  double E = del_E;

  while(E < v){
    if(E < v){
      q = sqrt(E/A);
      beta = sqrt(abs(v-E)/A);
      func_value = (cosh(beta*b)*cos(q*w) - ((pow(q,2) - pow(beta,2))/(2*q*beta))*sinh(beta*b)*sin(q*w));
    } else if(E >= v){
      q = sqrt(E/A);
      beta = sqrt(abs(v-E)/A);
      func_value = (cos(beta*b)*cos(q*w) - ((pow(q,2) - pow(beta,2))/(2*q*beta))*sin(beta*b)*sin(q*w));
    }
    if(func_value <= 1 && func_value >=-1){
      func.push_back(func_value);
      energies.push_back(E);
    }
    E = E + 0.0001;
  }
  for(int i = 0; i < energies.size(); i++){
    k.push_back(root_finder(0,M_PI/a,energies[i],&func3));
  }
  // all_data = mirror(k,energies,func);
  all_data.push_back(func);
  all_data.push_back(energies);
  all_data.push_back(k);
  return all_data;
}
