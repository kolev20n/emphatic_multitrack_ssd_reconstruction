#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

#include "TFile.h"
#include "TTree.h"

#include "matej_minimum_data_structures_and_functions.h"

// using namespace TrackTools;

using namespace std;

MagneticField::MagneticField(string &filename) :
step(0), start{-50, -50, -50}
{
  
  
  TFile mfFile(filename.c_str(), "READ");
  cout << " ==> Opening file " << filename << " to read magnetic field..."
  << endl;
  
  if (mfFile.IsZombie()) {
    cerr << "Error opening file" << endl;
    exit(-1);
  }
  TTree *tree = (TTree*) mfFile.Get("magField");
  
  double x;
  double y;
  double z;
  double Bx;
  double By;
  double Bz;
  tree->SetBranchAddress("x", &x);
  tree->SetBranchAddress("y", &y);
  tree->SetBranchAddress("z", &z);
  tree->SetBranchAddress("Bx", &Bx);
  tree->SetBranchAddress("By", &By);
  tree->SetBranchAddress("Bz", &Bz);
  
  int nEntries = tree->GetEntries();
  
  tree->GetEntry(0);
  double xVal = x;
  double yVal = y;
  double zVal = z;
  
  //step = 0;
  tree->GetEntry(1);
  if(abs(xVal - x) > step) step = abs(xVal - x);
  else if(abs(yVal - y) > step) step = abs(yVal - y);
  else step = abs(zVal - z);
  
  start[0] = xVal;
  start[1] = yVal;
  start[2] = zVal;
  //start = {-60, -60, -60};
  
  for(int i = 0; i < nEntries; i++){
    tree->GetEntry(i);
    
    
    int indX = (x-start[0])/step;
    int indY = (y-start[1])/step;
    int indZ = (z-start[2])/step;
    
    std::vector<double> temp;
    temp.push_back(Bx);
    temp.push_back(By);
    temp.push_back(Bz);
    
    //if(By > 10) cout << x << " " << y << " " <<z << " " <<By << endl;
    
    field[indX][indY][indZ] = temp;
    /*if(By > 10) {cout << x << " " << y << " " <<z << " " <<By << endl;
     cout << indX << " " << indY << " " << indZ << " " << field[indX][indY][indZ].at(0) << endl;}*/
    if(indX == 50 && indY == 50 && indZ == 65) cout << By << endl;
  }
  
  ///////////////////////////////////////////////////////////////
  // Close the file
  cout << " ==> Closing file " << filename << endl;
  mfFile.Close();
  
  }
  
  
  MagneticField::~MagneticField() {
    
  }
  
  
  // Member functions
  
  void MagneticField::GetField(double x[3], double &Bx, double &By, double &Bz)
  {
    
    
    double indX = (x[0]-start[0])/step;
    double indY = (x[1]-start[1])/step;
    double indZ = (x[2]-start[2])/step;
    
    int ix[2] = {int(floor(indX)), int(ceil(indX))};
    int iy[2] = {int(floor(indY)), int(ceil(indY))};
    int iz[2] = {int(floor(indZ)), int(ceil(indZ))};
    
    bool skip = false;
    
    if(field.find(ix[0]) == field.end()) skip = true;
    else if(field.find(ix[1]) == field.end()) skip = true;
    else{
      if(field.at(ix[0]).find(iy[0]) == field.at(ix[0]).end()) skip = true;
      else if(field.at(ix[0]).find(iy[1]) == field.at(ix[0]).end()) skip = true;
      else if(field.at(ix[1]).find(iy[0]) == field.at(ix[1]).end()) skip = true;
      else if(field.at(ix[1]).find(iy[1]) == field.at(ix[1]).end()) skip = true;
      else{
        if(field.at(ix[0]).at(iy[0]).find(iz[0]) ==field.at(ix[0]).at(iy[0]).end()) skip = true;
        else if(field.at(ix[0]).at(iy[0]).find(iz[1]) ==field.at(ix[0]).at(iy[0]).end()) skip = true;
        else if(field.at(ix[0]).at(iy[1]).find(iz[0]) ==field.at(ix[0]).at(iy[1]).end()) skip = true;
        else if(field.at(ix[0]).at(iy[1]).find(iz[1]) ==field.at(ix[0]).at(iy[1]).end()) skip = true;
        else if(field.at(ix[1]).at(iy[0]).find(iz[0]) ==field.at(ix[1]).at(iy[0]).end()) skip = true;
        else if(field.at(ix[1]).at(iy[0]).find(iz[1]) ==field.at(ix[1]).at(iy[0]).end()) skip = true;
        else if(field.at(ix[1]).at(iy[1]).find(iz[0]) ==field.at(ix[1]).at(iy[1]).end()) skip = true;
        else if(field.at(ix[1]).at(iy[1]).find(iz[1]) ==field.at(ix[1]).at(iy[1]).end()) skip = true;
      }
    }
    
    
    if(skip){
      
      Bx = 0;
      By = 0;
      Bz = 0;
      return;
    }
    
    
    double sumx = 0;
    double sumy = 0;
    double sumz = 0;
    double norm = 0;
    
    for(int i = 0; i < 2; i++){
      for(int j = 0; j < 2; j++){
        for(int k = 0; k < 2; k++){
          double dist = sqrt((indX-ix[i])*(indX-ix[i]) + (indY-iy[j])*(indY-iy[j]) + (indZ-iz[k])*(indZ-iz[k]));
          sumx += field.at(ix[i]).at(iy[j]).at(iz[k]).at(0)*dist;
          sumy += field.at(ix[i]).at(iy[j]).at(iz[k]).at(1)*dist;
          sumz += field.at(ix[i]).at(iy[j]).at(iz[k]).at(2)*dist;
          norm += dist;
          
          /*sumx = field.at(ix[i]).at(iy[j]).at(iz[k]).at(0)*1;
           sumy = field.at(ix[i]).at(iy[j]).at(iz[k]).at(1)*1;
           sumz = field.at(ix[i]).at(iy[j]).at(iz[k]).at(2)*1;
           norm = 1;*/
          
        }
      }
    }
    //if(ix[50] == 0 && ix[50] == 0 && iz[0] == 65) cout << By << endl;
    //cout << field.at(50).at(50).at(80).at(1) << endl;
    Bx = (sumx/norm);
    By = (sumy/norm);
    Bz = (sumz/norm);
    /*cout << start[0] << " " << start[1] << start[2] << endl;
     cout << x[0] << " " << x[1] << " " << x[2] << endl;
     cout << ix[0] << " " << ix[0] << " " << iz[0] << " " << By << endl;
     cout << Bx << " " << By << " " << Bz << " " << endl;
     cout << "*************************************" << endl;*/
  }
  

TrackPar::TrackPar()
{
  fType = eUnknown;
  fCharge = 0;
  fZ = 0;
  fLength = 0;
  for (int i = 0; i < 15; i++){
    if(i < 5){
      fPar[i] = 0;
    }
    fCov[i] = 0;
  }
}
//**************************************************************************************************
TrackPar::TrackPar(TrackingType type)
{
  fType = type;
  fCharge = 0;
  fZ = 0;
  fLength = 0;
  for (int i = 0; i < 15; i++){
    if(i < 5){
      fPar[i] = 0;
    }
    fCov[i] = 0;
  }
}
//**************************************************************************************************
TrackPar::TrackPar(TrackingType type, double z, double *par)
{
  fType = type;
  fZ = z;
  fLength = 0;
  //ParameterCheck(par);
  ChargeCheck(type, par);
  
  for(int i = 0; i < 15; i++){
    if(i < 5){
      fPar[i] = par[i];
    }
    fCov[i] = 0;
  }
  
}
//**************************************************************************************************
TrackPar::TrackPar(TrackingType type, double z, double *par, double *cov)
{
  fType = type;
  fZ = z;
  fLength = 0;
  //ParameterCheck(par);
  //CovarianceCheck(cov);
  ChargeCheck(type, par);
  
  for(int i = 0; i < 15; i++){
    if(i < 5){
      fPar[i] = par[i];
    }
    fCov[i] = cov[i];
  }
}
//**************************************************************************************************

  void TrackPar::ChargeCheck(TrackingType type, double *par){
    
    switch(type){
      case eUnknown:
        cerr << "WARNING: " << __FUNCTION__ << ": TrackingType is unknown. Positive charge is assumed!" << endl;
        fCharge = 1;
        break;
      case eCartesian:
        cerr << "WARNING: " << __FUNCTION__ << ": Charge cannot be determined from eCartesian type. Positive charge is assumed!" << endl;
        fCharge = 1;
        break;
      case eNA61:
        if(par[0] > 0){
          fCharge = 1;
        }
        else fCharge = -1;
        break;
      case eKisel:
        if(par[4] > 0){
          fCharge = 1;
        }
        else fCharge = -1;
        break;
    }
  }
  
double TrackPar::GetPar(int i){
  //CheckIndex(i, fPar);
  return fPar[i];
}

//**************************************************************************************************
double TrackPar::GetCov(int i){
  //CheckIndex(i, fCov);
  return fCov[i];
}

//**************************************************************************************************
double TrackPar::GetStd(int i){
  //CheckIndex(i, fPar);
  
  int index = 7*(i+1)-6-(i+1)*(i+2)/2;
  return sqrt(fCov[index]);
}
//**************************************************************************************************
void TrackPar::SetPar(int i, double val){
  //CheckIndex(i, fPar);
  fPar[i] = val;
}
//**************************************************************************************************
void TrackPar::SetPar(double *par){
  //ParameterCheck(par);
  
  for(int i = 0; i < 5; i++){
    fPar[i] = par[i];
  }
}
//**************************************************************************************************
void TrackPar::SetCov(int i, double val){
  //CheckIndex(i, fCov);
  fCov[i] = val;
}
//**************************************************************************************************
void TrackPar::SetCov(double *cov)
{
  //ParameterCheck(cov);
  
  for(int i = 0; i < 15; i++){
    fCov[i] = cov[i];
  }
}

TrackPar& TrackPar::operator= (TrackPar &parIn)
{
  if(this != &parIn){
    fType = parIn.GetType();
    fCharge = parIn.GetCharge();
    fZ = parIn.GetZ();
    fLength = parIn.GetLength();
    
    for(int i = 0; i < 15; i++){
      if(i < 5){
        fPar[i] = parIn.GetPar(i);
      }
      fCov[i] = parIn.GetCov(i);
    }
  }
  return *this;
}
//**************************************************************************************************

void TrackPar::Print(){
  
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "Type = " << fType << endl;
  cout << "z = " << fZ << endl;
  cout << "q = " << fCharge << endl;
  cout << "*****************************TRACK PARAMETERS*****************************" << endl;
  for(int i = 0; i < 5; i++){
    cout << "p" << i+1 << " = " << fPar[i] << endl;
  }
  cout << "*****************************COVARIANCE MATRIX*****************************" << endl;
  int ind1 = 0;
  //int ind2 = 0;
  for(int i = 0; i < 5; i++){
    for(int j = 0; j < 5; j++){
      if(j >= i){
        cout << fCov[ind1] << "  ";
        ind1++;
      }
      else {
        cout << "    ";
        
      }
    }
    cout << endl;
  }
  cout << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
}

  void TransformCov(double (*J)[5], double *c){
    //double Jt[5][5];
    /*double cT[5][5] = {  {c[0][0] ,c[0][1] , c[0][2], c[0][3], c[0][4]},
     {c[1][0] ,c[1][1] , c[1][2], c[1][3], c[1][4]},
     {c[2][0] ,c[2][1] , c[2][2], c[2][3], c[2][4]},
     {c[3][0] ,c[3][1] , c[3][2], c[3][3], c[3][4]},
     {c[4][0] ,c[4][1] , c[4][2], c[4][3], c[4][4]}};*/
    
    double cT[15] = {c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], c[9], c[10], c[11], c[12], c[13], c[14]};
    /*for(int i = 0; i < 5; i++){
     for(int j = 0; j < 5; j++){
     Jt[i][j] = J[j][i];
     cTemp[i][j] = 0;
     }
     }*/
    
    /*for(int i = 0; i < 5; i++){
     for(int j = 0; j < 5; j++){
     for(int k = 0; k < 5; k++){
     cTemp[i][j] += cov[i][k]*J[j][k];
     }
     }
     }
     
     for(int i = 0; i < 5; i++){
     for(int j = 0; j < 5; j++){
     cov[i][j] = 0;
     for(int k = 0; k < 5; k++){
     cov[i][j] += J[i][k]*cTemp[k][j];
     }
     }
     }*/
    
    c[0] = J[0][0]*(cT[0]*J[0][0] + cT[1]*J[0][1] + cT[2]*J[0][2] + cT[3]*J[0][3] + cT[4]*J[0][4])
    + J[0][1]*(cT[1]*J[0][0] + cT[5]*J[0][1] + cT[6]*J[0][2] + cT[7]*J[0][3] + cT[8]*J[0][4])
    + J[0][2]*(cT[2]*J[0][0] + cT[6]*J[0][1] + cT[9]*J[0][2] + cT[10]*J[0][3] + cT[11]*J[0][4])
    + J[0][3]*(cT[3]*J[0][0] + cT[7]*J[0][1] + cT[10]*J[0][2] + cT[12]*J[0][3] + cT[13]*J[0][4])
    + J[0][4]*(cT[4]*J[0][0] + cT[8]*J[0][1] + cT[11]*J[0][2] + cT[13]*J[0][3] + cT[14]*J[0][4]);
    
    
    //cout << c[0] << endl;
    c[1] = J[0][0]*(cT[0]*J[1][0] + cT[1]*J[1][1] + cT[2]*J[1][2] + cT[3]*J[1][3] + cT[4]*J[1][4])
    + J[0][1]*(cT[1]*J[1][0] + cT[5]*J[1][1] + cT[6]*J[1][2] + cT[7]*J[1][3] + cT[8]*J[1][4])
    + J[0][2]*(cT[2]*J[1][0] + cT[6]*J[1][1] + cT[9]*J[1][2] + cT[10]*J[1][3] + cT[11]*J[1][4])
    + J[0][3]*(cT[3]*J[1][0] + cT[7]*J[1][1] + cT[10]*J[1][2] + cT[12]*J[1][3] + cT[13]*J[1][4])
    + J[0][4]*(cT[4]*J[1][0] + cT[8]*J[1][1] + cT[11]*J[1][2] + cT[13]*J[1][3] + cT[14]*J[1][4]);
    //cout << J[1][0] << " " << J[1][1] << " " << J[1][2] << " " << J[1][3] << " " << J[1][4] << endl;
    //cout << c[3] << " " << c[7] << " " << c[10] << " " << c[12] << " " << c[13] << endl;
    c[2] = J[0][0]*(cT[0]*J[2][0] + cT[1]*J[2][1] + cT[2]*J[2][2] + cT[3]*J[2][3] + cT[4]*J[2][4])
    + J[0][1]*(cT[1]*J[2][0] + cT[5]*J[2][1] + cT[6]*J[2][2] + cT[7]*J[2][3] + cT[8]*J[2][4])
    + J[0][2]*(cT[2]*J[2][0] + cT[6]*J[2][1] + cT[9]*J[2][2] + cT[10]*J[2][3] + cT[11]*J[2][4])
    + J[0][3]*(cT[3]*J[2][0] + cT[7]*J[2][1] + cT[10]*J[2][2] + cT[12]*J[2][3] + cT[13]*J[2][4])
    + J[0][4]*(cT[4]*J[2][0] + cT[8]*J[2][1] + cT[11]*J[2][2] + cT[13]*J[2][3] + cT[14]*J[2][4]);
    
    
    c[3] = J[0][0]*(cT[0]*J[3][0] + cT[1]*J[3][1] + cT[2]*J[3][2] + cT[3]*J[3][3] + cT[4]*J[3][4])
    + J[0][1]*(cT[1]*J[3][0] + cT[5]*J[3][1] + cT[6]*J[3][2] + cT[7]*J[3][3] + cT[8]*J[3][4])
    + J[0][2]*(cT[2]*J[3][0] + cT[6]*J[3][1] + cT[9]*J[3][2] + cT[10]*J[3][3] + cT[11]*J[3][4])
    + J[0][3]*(cT[3]*J[3][0] + cT[7]*J[3][1] + cT[10]*J[3][2] + cT[12]*J[3][3] + cT[13]*J[3][4])
    + J[0][4]*(cT[4]*J[3][0] + cT[8]*J[3][1] + cT[11]*J[3][2] + cT[13]*J[3][3] + cT[14]*J[3][4]);
    
    
    c[4] = J[0][0]*(cT[0]*J[4][0] + cT[1]*J[4][1] + cT[2]*J[4][2] + cT[3]*J[4][3] + cT[4]*J[4][4])
    + J[0][1]*(cT[1]*J[4][0] + cT[5]*J[4][1] + cT[6]*J[4][2] + cT[7]*J[4][3] + cT[8]*J[4][4])
    + J[0][2]*(cT[2]*J[4][0] + cT[6]*J[4][1] + cT[9]*J[4][2] + cT[10]*J[4][3] + cT[11]*J[4][4])
    + J[0][3]*(cT[3]*J[4][0] + cT[7]*J[4][1] + cT[10]*J[4][2] + cT[12]*J[4][3] + cT[13]*J[4][4])
    + J[0][4]*(cT[4]*J[4][0] + cT[8]*J[4][1] + cT[11]*J[4][2] + cT[13]*J[4][3] + cT[14]*J[4][4]);
    
    
    c[5] = J[1][0]*(cT[0]*J[1][0] + cT[1]*J[1][1] + cT[2]*J[1][2] + cT[3]*J[1][3] + cT[4]*J[1][4])
    + J[1][1]*(cT[1]*J[1][0] + cT[5]*J[1][1] + cT[6]*J[1][2] + cT[7]*J[1][3] + cT[8]*J[1][4])
    + J[1][2]*(cT[2]*J[1][0] + cT[6]*J[1][1] + cT[9]*J[1][2] + cT[10]*J[1][3] + cT[11]*J[1][4])
    + J[1][3]*(cT[3]*J[1][0] + cT[7]*J[1][1] + cT[10]*J[1][2] + cT[12]*J[1][3] + cT[13]*J[1][4])
    + J[1][4]*(cT[4]*J[1][0] + cT[8]*J[1][1] + cT[11]*J[1][2] + cT[13]*J[1][3] + cT[14]*J[1][4]);
    
    c[6] = J[1][0]*(cT[0]*J[2][0] + cT[1]*J[2][1] + cT[2]*J[2][2] + cT[3]*J[2][3] + cT[4]*J[2][4])
    + J[1][1]*(cT[1]*J[2][0] + cT[5]*J[2][1] + cT[6]*J[2][2] + cT[7]*J[2][3] + cT[8]*J[2][4])
    + J[1][2]*(cT[2]*J[2][0] + cT[6]*J[2][1] + cT[9]*J[2][2] + cT[10]*J[2][3] + cT[11]*J[2][4])
    + J[1][3]*(cT[3]*J[2][0] + cT[7]*J[2][1] + cT[10]*J[2][2] + cT[12]*J[2][3] + cT[13]*J[2][4])
    + J[1][4]*(cT[4]*J[2][0] + cT[8]*J[2][1] + cT[11]*J[2][2] + cT[13]*J[2][3] + cT[14]*J[2][4]);
    
    
    c[7] = J[1][0]*(cT[0]*J[3][0] + cT[1]*J[3][1] + cT[2]*J[3][2] + cT[3]*J[3][3] + cT[4]*J[3][4])
    + J[1][1]*(cT[1]*J[3][0] + cT[5]*J[3][1] + cT[6]*J[3][2] + cT[7]*J[3][3] + cT[8]*J[3][4])
    + J[1][2]*(cT[2]*J[3][0] + cT[6]*J[3][1] + cT[9]*J[3][2] + cT[10]*J[3][3] + cT[11]*J[3][4])
    + J[1][3]*(cT[3]*J[3][0] + cT[7]*J[3][1] + cT[10]*J[3][2] + cT[12]*J[3][3] + cT[13]*J[3][4])
    + J[1][4]*(cT[4]*J[3][0] + cT[8]*J[3][1] + cT[11]*J[3][2] + cT[13]*J[3][3] + cT[14]*J[3][4]);
    
    
    c[8] = J[1][0]*(cT[0]*J[4][0] + cT[1]*J[4][1] + cT[2]*J[4][2] + cT[3]*J[4][3] + cT[4]*J[4][4])
    + J[1][1]*(cT[1]*J[4][0] + cT[5]*J[4][1] + cT[6]*J[4][2] + cT[7]*J[4][3] + cT[8]*J[4][4])
    + J[1][2]*(cT[2]*J[4][0] + cT[6]*J[4][1] + cT[9]*J[4][2] + cT[10]*J[4][3] + cT[11]*J[4][4])
    + J[1][3]*(cT[3]*J[4][0] + cT[7]*J[4][1] + cT[10]*J[4][2] + cT[12]*J[4][3] + cT[13]*J[4][4])
    + J[1][4]*(cT[4]*J[4][0] + cT[8]*J[4][1] + cT[11]*J[4][2] + cT[13]*J[4][3] + cT[14]*J[4][4]);
    
    
    c[9] = J[2][0]*(cT[0]*J[2][0] + cT[1]*J[2][1] + cT[2]*J[2][2] + cT[3]*J[2][3] + cT[4]*J[2][4])
    + J[2][1]*(cT[1]*J[2][0] + cT[5]*J[2][1] + cT[6]*J[2][2] + cT[7]*J[2][3] + cT[8]*J[2][4])
    + J[2][2]*(cT[2]*J[2][0] + cT[6]*J[2][1] + cT[9]*J[2][2] + cT[10]*J[2][3] + cT[11]*J[2][4])
    + J[2][3]*(cT[3]*J[2][0] + cT[7]*J[2][1] + cT[10]*J[2][2] + cT[12]*J[2][3] + cT[13]*J[2][4])
    + J[2][4]*(cT[4]*J[2][0] + cT[8]*J[2][1] + cT[11]*J[2][2] + cT[13]*J[2][3] + cT[14]*J[2][4]);
    
    c[10] = J[2][0]*(cT[0]*J[3][0] + cT[1]*J[3][1] + cT[2]*J[3][2] + cT[3]*J[3][3] + cT[4]*J[3][4])
    + J[2][1]*(cT[1]*J[3][0] + cT[5]*J[3][1] + cT[6]*J[3][2] + cT[7]*J[3][3] + cT[8]*J[3][4])
    + J[2][2]*(cT[2]*J[3][0] + cT[6]*J[3][1] + cT[9]*J[3][2] + cT[10]*J[3][3] + cT[11]*J[3][4])
    + J[2][3]*(cT[3]*J[3][0] + cT[7]*J[3][1] + cT[10]*J[3][2] + cT[12]*J[3][3] + cT[13]*J[3][4])
    + J[2][4]*(cT[4]*J[3][0] + cT[8]*J[3][1] + cT[11]*J[3][2] + cT[13]*J[3][3] + cT[14]*J[3][4]);
    
    
    c[11] = J[2][0]*(cT[0]*J[4][0] + cT[1]*J[4][1] + cT[2]*J[4][2] + cT[3]*J[4][3] + cT[4]*J[4][4])
    + J[2][1]*(cT[1]*J[4][0] + cT[5]*J[4][1] + cT[6]*J[4][2] + cT[7]*J[4][3] + cT[8]*J[4][4])
    + J[2][2]*(cT[2]*J[4][0] + cT[6]*J[4][1] + cT[9]*J[4][2] + cT[10]*J[4][3] + cT[11]*J[4][4])
    + J[2][3]*(cT[3]*J[4][0] + cT[7]*J[4][1] + cT[10]*J[4][2] + cT[12]*J[4][3] + cT[13]*J[4][4])
    + J[2][4]*(cT[4]*J[4][0] + cT[8]*J[4][1] + cT[11]*J[4][2] + cT[13]*J[4][3] + cT[14]*J[4][4]);
    
    
    c[12] = J[3][0]*(cT[0]*J[3][0] + cT[1]*J[3][1] + cT[2]*J[3][2] + cT[3]*J[3][3] + cT[4]*J[3][4])
    + J[3][1]*(cT[1]*J[3][0] + cT[5]*J[3][1] + cT[6]*J[3][2] + cT[7]*J[3][3] + cT[8]*J[3][4])
    + J[3][2]*(cT[2]*J[3][0] + cT[6]*J[3][1] + cT[9]*J[3][2] + cT[10]*J[3][3] + cT[11]*J[3][4])
    + J[3][3]*(cT[3]*J[3][0] + cT[7]*J[3][1] + cT[10]*J[3][2] + cT[12]*J[3][3] + cT[13]*J[3][4])
    + J[3][4]*(cT[4]*J[3][0] + cT[8]*J[3][1] + cT[11]*J[3][2] + cT[13]*J[3][3] + cT[14]*J[3][4]);
    
    c[13] = J[3][0]*(cT[0]*J[4][0] + cT[1]*J[4][1] + cT[2]*J[4][2] + cT[3]*J[4][3] + cT[4]*J[4][4])
    + J[3][1]*(cT[1]*J[4][0] + cT[5]*J[4][1] + cT[6]*J[4][2] + cT[7]*J[4][3] + cT[8]*J[4][4])
    + J[3][2]*(cT[2]*J[4][0] + cT[6]*J[4][1] + cT[9]*J[4][2] + cT[10]*J[4][3] + cT[11]*J[4][4])
    + J[3][3]*(cT[3]*J[4][0] + cT[7]*J[4][1] + cT[10]*J[4][2] + cT[12]*J[4][3] + cT[13]*J[4][4])
    + J[3][4]*(cT[4]*J[4][0] + cT[8]*J[4][1] + cT[11]*J[4][2] + cT[13]*J[4][3] + cT[14]*J[4][4]);
    
    c[14] = J[4][0]*(cT[0]*J[4][0] + cT[1]*J[4][1] + cT[2]*J[4][2] + cT[3]*J[4][3] + cT[4]*J[4][4])
    + J[4][1]*(cT[1]*J[4][0] + cT[5]*J[4][1] + cT[6]*J[4][2] + cT[7]*J[4][3] + cT[8]*J[4][4])
    + J[4][2]*(cT[2]*J[4][0] + cT[6]*J[4][1] + cT[9]*J[4][2] + cT[10]*J[4][3] + cT[11]*J[4][4])
    + J[4][3]*(cT[3]*J[4][0] + cT[7]*J[4][1] + cT[10]*J[4][2] + cT[12]*J[4][3] + cT[13]*J[4][4])
    + J[4][4]*(cT[4]*J[4][0] + cT[8]*J[4][1] + cT[11]*J[4][2] + cT[13]*J[4][3] + cT[14]*J[4][4]);
  }
  
  void ConvertTrackPar(TrackPar &parIn, TrackingType typeOut){
    TrackingType typeIn = parIn.GetType();
    
    if(typeIn == typeOut){
      cerr << "WARNING: " << __FUNCTION__ << ": input and output tracking types are equal. Skiping conversion!" << endl;
      return;
    }
    else if(typeIn == eUnknown){
      cerr << "WARNING: " << __FUNCTION__ << ": unknown input tracking type. Skiping conversion!" << endl;
      return;
    }
    else if(typeOut == eUnknown){
      cerr << "WARNING: " << __FUNCTION__ << ": unknown output tracking type. Skiping conversion!" << endl;
      return;
    }
    
    
    double J[5][5] = {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};
    //double cov[5][5];
    
    
    /*for(int i = 0; i < 5; i++){
     //J[i] = new double[5];
     //cov[i] = new double[5];
     for(int j = 0; j < 5; j++){
     J[i][j] = 0;
     cov[i][j] = 0;
     }
     }  */
    //GetCovarianceMatrix(cov, parIn);
    
    switch(typeIn){
        
      case eCartesian:
      {
        switch(typeOut){
          case eNA61:
          {
            double pxz = sqrt(parIn.GetPar(2)*parIn.GetPar(2) + parIn.GetPar(4)*parIn.GetPar(4));
            double p1 = parIn.GetCharge()/pxz;
            double p2 = parIn.GetPar(3)/pxz;
            double p3 = atan(parIn.GetPar(4)/parIn.GetPar(2));
            double p4 = parIn.GetPar(0);
            double p5 = parIn.GetPar(1);
            
            /*J[0][0] = 1;
             J[1][1] = 1;
             J[2][2] = -parIn.GetCharge()*parIn.GetPar(2)/(pxz*pxz*pxz);
             J[2][4] = -parIn.GetCharge()*parIn.GetPar(4)/(pxz*pxz*pxz);
             J[3][2] = -parIn.GetPar(2)*parIn.GetPar(3)/(pxz*pxz*pxz);
             J[3][3] = 1/pxz;
             J[3][4] = -parIn.GetPar(3)*parIn.GetPar(4)/(pxz*pxz*pxz);
             J[4][2] = -2*parIn.GetPar(4)*parIn.GetPar(4)/(pxz*(pxz*pxz +(pxz-parIn.GetPar(2))*(pxz-parIn.GetPar(2))));
             J[4][4] = 2*parIn.GetPar(2)*parIn.GetPar(4)/(pxz*(pxz*pxz +(pxz-parIn.GetPar(2))*(pxz-parIn.GetPar(2))));*/
            cout << parIn.GetPar(0) << " " << parIn.GetPar(1) << endl;
            J[0][2] = -parIn.GetCharge()*parIn.GetPar(2)/(pxz*pxz*pxz);
            J[0][4] = -parIn.GetCharge()*parIn.GetPar(4)/(pxz*pxz*pxz);
            J[1][2] = -parIn.GetPar(2)*parIn.GetPar(3)/(pxz*pxz*pxz);
            J[1][3] = 1/pxz;
            J[1][4] = -parIn.GetPar(3)*parIn.GetPar(4)/(pxz*pxz*pxz);
            J[2][2] = -parIn.GetPar(4)/(pxz*pxz);
            J[2][4] = parIn.GetPar(2)/(pxz*pxz);
            J[3][0] = 1;
            J[4][1] = 1;
            
            parIn.SetPar(0, p1);
            parIn.SetPar(1, p2);
            parIn.SetPar(2, p3);
            parIn.SetPar(3, p4);
            parIn.SetPar(4, p5);
            parIn.SetType(typeOut);
            TransformCov(J, parIn.GetCov());
            break;
          }
          case eKisel:
          {
            double p = sqrt(parIn.GetPar(2)*parIn.GetPar(2) + parIn.GetPar(3)*parIn.GetPar(3) + parIn.GetPar(4)*parIn.GetPar(4));
            double p3 = parIn.GetPar(2)/parIn.GetPar(4);
            double p4 = parIn.GetPar(3)/parIn.GetPar(4);
            double p5 = parIn.GetCharge()/p;
            
            J[0][0] = 1;
            J[1][1] = 1;
            J[2][2] = 1/parIn.GetPar(4);
            J[2][4] = -parIn.GetPar(2)/parIn.GetPar(4)/parIn.GetPar(4);
            J[3][3] = 1/parIn.GetPar(4);
            J[3][4] = -parIn.GetPar(3)/parIn.GetPar(4)/parIn.GetPar(4);
            J[4][2] = -parIn.GetCharge()*parIn.GetPar(2)/p/p/p;
            J[4][3] = -parIn.GetCharge()*parIn.GetPar(3)/p/p/p;
            J[4][4] = -parIn.GetCharge()*parIn.GetPar(4)/p/p/p;
            parIn.SetPar(2, p3);
            parIn.SetPar(3, p4);
            parIn.SetPar(4, p5);
            parIn.SetType(typeOut);
            TransformCov(J, parIn.GetCov());
            break;
          }
            
        }
        break;
      }
      case eNA61:
      {
        switch(typeOut){
          case eCartesian:
          {
            //double ca = cos(parIn.GetPar(4)/2.);
            //double ta = tan(parIn.GetPar(4)/2.);
            double p1 = parIn.GetPar(3);
            double p2 = parIn.GetPar(4);
            double p3 = parIn.GetCharge()*cos(parIn.GetPar(2))/parIn.GetPar(0);
            double p4 = parIn.GetCharge()*parIn.GetPar(1)/parIn.GetPar(0);
            double p5 = parIn.GetCharge()*sin(parIn.GetPar(2))/parIn.GetPar(0);
            
            /*J[0][0] = 1;
             J[1][1] = 1;
             J[2][2] = -p3/parIn.GetPar(2);
             J[2][4] = -parIn.GetCharge()/(2*parIn.GetPar(2)*ca*ca);
             J[3][2] = -p4/parIn.GetPar(2);
             J[3][3] = parIn.GetCharge()/parIn.GetPar(2);
             J[4][2] = -p5/parIn.GetPar(2);
             J[4][4] = (1-ta)/(2*ca*ca*parIn.GetPar(2)*sqrt(ta*(2-ta)));*/
            
            J[0][3] = 1;
            J[1][4] = 1;
            J[2][0] = -cos(parIn.GetPar(2))/(parIn.GetPar(0)*parIn.GetPar(0));
            J[2][2] = -sin(parIn.GetPar(2))/parIn.GetPar(0);
            J[3][0] = -parIn.GetCharge()*parIn.GetPar(1)/(parIn.GetPar(0)*parIn.GetPar(0));
            J[3][1] = parIn.GetCharge()/parIn.GetPar(0);
            J[4][0] = -sin(parIn.GetPar(2))/(parIn.GetPar(0)*parIn.GetPar(0));
            J[4][2] = cos(parIn.GetPar(2))/parIn.GetPar(0);
            
            parIn.SetPar(0, p1);
            parIn.SetPar(1, p2);
            parIn.SetPar(2, p3);
            parIn.SetPar(3, p4);
            parIn.SetPar(4, p5);
            parIn.SetType(typeOut);
            TransformCov(J, parIn.GetCov());
            break;
          }
          case eKisel:
          {
            
            
            double b=sqrt(1+parIn.GetPar(1)*parIn.GetPar(1));
            double p1 = parIn.GetPar(3);
            double p2 = parIn.GetPar(4);
            double p3 = 1/tan(parIn.GetPar(2));
            double p4 = parIn.GetPar(1)/fabs(sin(parIn.GetPar(2)));
            double p5 = parIn.GetPar(0)/b;
            
            J[0][3] = 1;
            J[1][4] = 1;
            J[2][2] = -1/(sin(parIn.GetPar(2))*sin(parIn.GetPar(2)));
            J[3][1] = 1/sin(parIn.GetPar(2));
            J[3][2] = -parIn.GetPar(1)*cos(parIn.GetPar(2))/(sin(parIn.GetPar(2))*sin(parIn.GetPar(2)));
            J[4][0] = 1/b;
            J[4][1] = -parIn.GetPar(0)*parIn.GetPar(1)/(b*b*b);
            
            //cout << J[4][0] << endl;
            parIn.SetPar(0, p1);
            parIn.SetPar(1, p2);
            parIn.SetPar(2, p3);
            parIn.SetPar(3, p4);
            parIn.SetPar(4, p5);
            /*double ta = tan(parIn.GetPar(4)/2.);
             double d = sqrt(ta*(2-ta));
             double d1 = sqrt(1+parIn.GetPar(3)*parIn.GetPar(3));
             double ca = cos(parIn.GetPar(4)/2.)*cos(parIn.GetPar(4)/2.);
             //cout << ta << " " << parIn.GetPar(3) << " " << d << endl;
             double p3 = (1-ta)/d;
             double p4 = parIn.GetPar(3)/d;
             double p5 = parIn.GetPar(2)/d1;
             
             J[0][0] = 1;
             J[1][1] = 1;
             J[2][4] = -1/(2*ca*d*d*d);
             J[3][3] = 1/d;
             J[3][4] = -parIn.GetPar(3)*(1-ta)/(2*ca*d*d*d);
             J[4][2] = 1/d1;
             J[4][3] = -parIn.GetPar(2)*parIn.GetPar(3)/(d1*d1*d1);
             parIn.SetPar(2, p3);
             parIn.SetPar(3, p4);
             parIn.SetPar(4, p5);  */
            parIn.SetType(typeOut);
            TransformCov(J, parIn.GetCov());
            break;
          }
            
        }
        break;
      }
      case eKisel:
      {
        switch(typeOut){
          case eCartesian:
          {
            double a=sqrt(1+parIn.GetPar(2)*parIn.GetPar(2)+parIn.GetPar(3)*parIn.GetPar(3));
            double p3 = parIn.GetCharge()*parIn.GetPar(2)/parIn.GetPar(4)/a;
            double p4 = parIn.GetCharge()*parIn.GetPar(3)/parIn.GetPar(4)/a;
            double p5 = parIn.GetCharge()/parIn.GetPar(4)/a;
            
            J[0][0] = 1;
            J[1][1] = 1;
            J[2][2] = parIn.GetCharge()*(1+parIn.GetPar(3)*parIn.GetPar(3))/(a*a*a*parIn.GetPar(4));
            J[2][3] = -parIn.GetCharge()*parIn.GetPar(2)*parIn.GetPar(3)/(a*a*a*parIn.GetPar(4));
            J[2][4] = -parIn.GetCharge()*parIn.GetPar(2)/(a*parIn.GetPar(4)*parIn.GetPar(4));
            J[3][2] = -parIn.GetCharge()*parIn.GetPar(2)*parIn.GetPar(3)/(a*a*a*parIn.GetPar(4));
            J[3][3] = parIn.GetCharge()*(1+parIn.GetPar(2)*parIn.GetPar(2))/(a*a*a*parIn.GetPar(4));
            J[3][4] = -parIn.GetCharge()*parIn.GetPar(3)/(a*parIn.GetPar(4)*parIn.GetPar(4));
            J[4][2] = -parIn.GetCharge()*parIn.GetPar(2)/(a*a*a*parIn.GetPar(4));
            J[4][3] = -parIn.GetCharge()*parIn.GetPar(3)/(a*a*a*parIn.GetPar(4));
            J[4][4] = -parIn.GetCharge()/(a*parIn.GetPar(4)*parIn.GetPar(4));
            parIn.SetPar(2, p3);
            parIn.SetPar(3, p4);
            parIn.SetPar(4, p5);
            parIn.SetType(typeOut);
            TransformCov(J, parIn.GetCov());
            break;
          }
          case eNA61:
          {
            double a=sqrt(1+parIn.GetPar(2)*parIn.GetPar(2)+parIn.GetPar(3)*parIn.GetPar(3));
            double b=sqrt(1+parIn.GetPar(2)*parIn.GetPar(2));
            double p1 = parIn.GetPar(4)*a/b;
            double p2 = parIn.GetPar(3)/b;
            double p3 = atan(1/parIn.GetPar(2));
            double p4 = parIn.GetPar(0);
            double p5 = parIn.GetPar(1);
            
            J[0][2] = -parIn.GetPar(2)*parIn.GetPar(3)*parIn.GetPar(3)*parIn.GetPar(4)/(a*b*b*b);
            J[0][3] = parIn.GetPar(3)*parIn.GetPar(4)/(a*b);
            J[0][4] = a/b;
            J[1][2] = -parIn.GetPar(2)*parIn.GetPar(3)/(b*b*b);;
            J[1][3] = 1/b;
            J[2][2] = -1/(b*b);
            J[3][0] = 1;
            J[4][1] = 1;
            
            parIn.SetPar(0, p1);
            parIn.SetPar(1, p2);
            parIn.SetPar(2, p3);
            parIn.SetPar(3, p4);
            parIn.SetPar(4, p5);
            /*J[0][0] = 1;
             J[1][1] = 1;
             J[2][2] = -parIn.GetPar(2)*parIn.GetPar(3)*parIn.GetPar(3)*parIn.GetPar(4)/(a*b*b*b);
             J[2][3] = parIn.GetPar(3)*parIn.GetPar(4)/(a*b);
             J[2][4] = a/b;
             J[3][2] = -parIn.GetPar(2)*parIn.GetPar(3)/(b*b*b);
             J[3][3] = 1/b;
             J[4][2] = -2/(b*b+(b-parIn.GetPar(2))*(b-parIn.GetPar(2)))/b;
             parIn.SetPar(2, p3);
             parIn.SetPar(3, p4);
             parIn.SetPar(4, p5);  */
            parIn.SetType(typeOut);
            TransformCov(J, parIn.GetCov());
            break;
          }
            
        }
        break;
      }
        
    }
    
    
    
    //SetCovarianceMatrix(cov, parIn);
  }


double TrackExtrapolation::fKappa = 0.000299792458L;

TrackExtrapolation::TrackExtrapolation(bool cov, double step, MagneticField *a_magnetic_field, double radLength)
:fMagField(a_magnetic_field)
{
  fStep = step;
  fRadLength = radLength;
  TrackPar temp(eKisel);
  fStartPar = temp;
  fStopPar = temp;
  fErrorEstimation = cov;
  
  //fMagField = new MagneticField(field);
}
                   
TrackExtrapolation::TrackExtrapolation(bool cov, TrackPar trPar,  double step, MagneticField *a_magnetic_field, double radLength)
:fMagField(a_magnetic_field)
{
  SetTrackPar(trPar);
  fStep = step;
  fRadLength = radLength;
  fErrorEstimation = cov;
  fStartPar = trPar;
  fStopPar = trPar;
  
  //fMagField = new MagneticField(field);
}

TrackExtrapolation::~TrackExtrapolation()
{
  //delete fMagField;
  //fMagField = NULL;
}

/*void TrackExtrapolation:ExtrapolateToPlane(double step, double radLength, Plane& plane){
 fRadLength = radLength;
 fStep = step;
 
 ExtrapolateToPlane(zStop, plane);
 }
 
 void TrackExtrapolation::ExtrapolateToPlane(double radLength, Plane& plane){
 fRadLength = radLength;
 
 ExtrapolateToPlane(plane);
 }
 
 void TrackExtrapolation:ExtrapolateToPlane(Plane& plane){
 
 Cluster cPred;
 SVector3 vRot;
 vRot(0) = 0;
 vRot(1) = 0;
 vRot(2) = 1;
 
 vRot = plane.GetRotationMatrix()*vRot;
 
 
 
 }*/


void TrackExtrapolation::ExtrapolateToPlane(const double &a, const double&b, const double& c, const double& d)
{
#ifdef DEBUG
  cout << "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
  cout << "DEBUG OUTPUT" << endl;
  cout << __func__ << " in " << __FILE__ << endl;
  cout << "   Extrapolate from z = " << fStopPar.GetZ() << " cm" << endl;
  int iter = 0;
#endif
  double z = - (d + a*fStopPar.GetPar(0) + b*fStopPar.GetPar(1))/(a*fStopPar.GetPar(2) + b*fStopPar.GetPar(3) + c);
  
  double nStep = abs(z-fStopPar.GetZ())/fStep;
  
  
  while(nStep > 1){
    
#ifdef DEBUG
    cout << "   Iteration: " << iter << endl;
    cout << "   nSteps = " << nStep << ", z = " << z << " cm" << endl;
    cout << "   Extrapolate from z = " << fStopPar.GetZ() << " cm" << endl;
    cout << "   Extrapolate to z = " << z - 0.25*(z-fStopPar.GetZ()) << " cm " << endl;
    cout << a << " " << b << " " << c << " " << d << endl;
    iter++;
#endif
    Extrapolate(z - 0.25*(z-fStopPar.GetZ()));
    
    z = - (d + a*fStopPar.GetPar(0) + b*fStopPar.GetPar(1))/(a*fStopPar.GetPar(2) + b*fStopPar.GetPar(3) + c);
    nStep = abs(z-fStopPar.GetZ())/fStep;
  }
  
  Extrapolate(z);
  
#ifdef DEBUG
  cout << "   Final z = " << fStopPar.GetZ() << " cm" << endl;
#endif
}
  
  void TrackExtrapolation::Extrapolate(double zStop){
    // Calculate number of extrapolation sts
    //fStopPar = fStartPar;
    int nsteps = (int) fabs(zStop - fStopPar.GetZ())/fStep;
    //cout << nsteps << endl;
    // Check if the extrapolation is backward or forward
    double st = fStep;
    if(zStop < fStopPar.GetZ())
      st *= -1;
    
    
    //cout << "radi" << endl;
    double J[5][5] = {{0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}, {0, 0, 0, 0, 0}};
    
    
    double first = fStopPar.GetZ();
    // Do extrapolation in small steps
    for(int i = 0; i < nsteps+1; i++){
      //do{
      
      if(i == nsteps){
        st = fabs(zStop - fStopPar.GetZ());
        if(zStop < first)
          st *= -1;
      }
      
      
      // Calculate needed parameters
      double x = fStopPar.GetPar(0);
      double y = fStopPar.GetPar(1);
      double tx = fStopPar.GetPar(2);
      double ty = fStopPar.GetPar(3);
      double tx2 = tx*tx;
      double ty2 = ty*ty;
      double txy = tx*ty;
      double qdp = fStopPar.GetPar(4);
      double h = fKappa*qdp*sqrt(1+tx2+ty2);
      
      double dC = fKappa*qdp*fKappa*qdp;
      
      
      double magPoint[3] = {x + tx*st/2., y + ty*st/2., fStopPar.GetZ() + st/2.};
      
      double Bx = 0;
      double By = 0;
      double Bz = 0;
      fMagField->GetField(magPoint, Bx, By, Bz);
      //cout << magPoint[0] << " " << magPoint[1] << " " << magPoint[2] << " " << Bx << " " << By << " " << Bz << endl;
      
      fStopPar.SetZ(fStopPar.GetZ() + st);
      
      
      double dtx1 = txy*Bx - (1+tx2)*By + ty*Bz;
      double dtx2 = tx*(3*ty2+1)*Bx*Bx -2*ty*(3*tx2+1)*Bx*By + (3*ty2-tx2+1)*Bx*Bz + 3*tx*(tx2+1)*By*By - 4*txy*By*Bz - tx*Bz*Bz;
      double dtx3 = 3*txy*(5*ty2+3)*Bx*Bx*Bx - 3*(3*tx2+3*ty2+15*tx2*ty2+1)*Bx*Bx*By + ty*(-10*tx2+15*ty2+9)*Bx*Bx*Bz +  9*tx*ty*(5*tx2+3)*Bx*By*By
      + tx*(10*tx2-40*ty2-2)*Bx*By*Bz - 11*tx*ty*Bx*Bz*Bz - 3*(tx2+1)*(5*tx2+1)*By*By*By + ty*(25*tx2+7)*By*By*Bz
      + (7*tx2-4*ty2+1)*By*Bz*Bz - ty*Bz*Bz*Bz;
      
      
      double dty1 = (1+ty2)*Bx - txy*By - tx*Bz;
      double dty2 = 3*ty*(ty2+1)*Bx*Bx - 2*tx*(3*ty2+1)*Bx*By - 4*txy*Bx*Bz + ty*(3*tx2+1)*By*By + (3*tx2-ty2+1)*By*Bz - ty*Bz*Bz;
      double dty3 = 3*(ty2+1)*(5*ty2+1)*Bx*Bx*Bx -9*txy*(5*ty2+3)*Bx*Bx*By - tx*(25*ty2+7)*Bx*Bx*Bz + 3*(3*tx*tx+3*ty*ty+15*tx2*ty2+1)*Bx*By*By
      - ty*(-40*tx2+10*ty2-2)*Bx*By*Bz + (4*tx2-7*ty2-1)*Bx*Bz*Bz - 3*txy*(5*tx2+3)*By*By*By
      - tx*(15*tx2-10*ty2+9)*By*By*Bz + 11*txy*By*Bz*Bz + tx*Bz*Bz*Bz;
      
      x = x + tx*st + h*dtx1*st*st/2. + h*h*dtx2*st*st*st/6. + h*h*h*dtx3*st*st*st*st/24.;
      y = y + ty*st + h*dty1*st*st/2. + h*h*dty2*st*st*st/6. + h*h*h*dty3*st*st*st*st/24.;
      
      // fKappa = 2.99792e-4
      fStopPar.AddLength(sqrt(pow(tx*st + h*dtx1*st*st/2. + h*h*dtx2*st*st*st/6. + h*h*h*dtx3*st*st*st*st/24.,2) +
                              pow(ty*st + h*dty1*st*st/2. + h*h*dty2*st*st*st/6. + h*h*h*dty3*st*st*st*st/24.,2) + st*st));
      fStopPar.SetPar(0, x);
      fStopPar.SetPar(1, y);
      fStopPar.SetPar(2, tx + h*dtx1*st + h*h*dtx2*st*st/2. + h*h*h*dtx3*st*st*st/6.);
      fStopPar.SetPar(3, ty + h*dty1*st + h*h*dty2*st*st/2. + h*h*h*dty3*st*st*st/6.);
      
      
      //REMOVE THIS
      /*if(fStopPar.GetZ()<-567.51){
       if(qdp < 0){
       fStopPar.SetPar(4, -1/(-st*4.5/1000. - 1/qdp));
       }
       else{
       fStopPar.SetPar(4, 1/(-st*4.5/1000. + 1/qdp));
       }
       }*/
      //**************************************
      if(fErrorEstimation){
        
        double dtx1dtx = ty*Bx - 2*tx*By;
        double dtx1dty = tx*Bx + Bz;
        
        double dtx2dtx = (3*ty2+1)*Bx*Bx - 12*txy*Bx*By - 2*tx*Bx*Bz + 3*(3*tx2+1)*By*By - 4*ty*By*Bz - Bz*Bz;
        double dtx2dty = 6*txy*Bx*Bx - 2*(3*tx2+1)*Bx*By + 6*ty*Bx*Bz - 4*tx*By*Bz;
        
        double dtx3dtx = 3*ty*(5*ty2+3)*Bx*Bx*Bx - 18*tx*(5*ty2+1)*Bx*Bx*By - 20*txy*Bx*Bx*Bz + 27*ty*(5*tx2+1)*Bx*By*By + 2*(15*tx2-20*ty2-1)*Bx*By*Bz
        - 11*ty*Bx*Bz*Bz - 12*tx*(5*tx2+3)*By*By*By + 40*txy*By*By*Bz + 14*tx*By*Bz*Bz;
        double dtx3dty = 9*tx*(5*ty2+1)*Bx*Bx*Bx - 18*ty*(5*tx2+1)*Bx*Bx*By + (-10*tx2+45*ty2+9)*Bx*Bx*Bz + 9*tx*(5*tx2+3)*Bx*By*By -80*txy*Bx*By*Bz
        - 11*tx*Bx*Bz*Bz + (25*tx2+7)*By*By*Bz - 8*ty*By*Bz*Bz - Bz*Bz*Bz;
        
        
        double dty1dtx = -ty*By - Bz;
        double dty1dty = 2*ty*Bx + tx*By;
        
        double dty2dtx = -2*(3*ty2+1)*Bx*By - 4*ty*Bx*Bz + 6*txy*By*By + 6*tx*By*Bz;
        double dty2dty = 3*(3*ty2+1)*Bx*Bx - 12*txy*Bx*By - 4*tx*Bx*Bz + (3*tx2+1)*By*By - 2*ty*By*Bz - Bz*Bz;
        
        double dty3dtx = -9*ty*(5*ty2+3)*Bx*Bx*By - (25*ty2+7)*Bx*Bx*Bz + 18*tx*(5*ty2+1)*Bx*By*By + 80*txy*Bx*By*Bz + 8*tx*Bx*Bz*Bz
        - 9*ty*(5*tx2+1)*By*By*By - (45*tx2-10*ty2+9)*By*By*Bz + 11*ty*By*Bz*Bz + Bz*Bz*Bz;
        double dty3dty = -12*ty*(5*ty2+3)*Bx*Bx*Bx - 27*tx*(5*ty2+1)*Bx*Bx*By - 40*txy*Bx*Bx*Bz + 18*ty*(5*tx2+1)*Bx*By*By + 2*(20*tx2-15*ty2+1)*Bx*By*Bz
        - 14*ty*Bx*Bz*Bz - 3*tx*(5*tx2+3)*By*By*By + 20*txy*By*By*Bz + 11*By*Bz*Bz;
        
        
        J[0][0] = 1;
        J[0][2] = st + h*dtx1dtx*st*st/2. + dC*tx*dtx1*st*st/2./h + h*h*dtx2dtx*st*st*st/6. + dC*tx*dtx2*st*st*st/3.
        + h*h*h*dtx3dtx*st*st*st*st/24. + dC*tx*h*dtx3*st*st*st*st/12.;
        J[0][3] = h*dtx1dty*st*st/2. + dC*ty*dtx1*st*st/2./h + h*h*dtx2dty*st*st*st/6. + dC*ty*dtx2*st*st*st/3.
        + h*h*h*dtx3dty*st*st*st*st/24. + dC*ty*h*dtx3*st*st*st*st/12.;
        J[0][4] = h*dtx1*st*st/(2.*qdp) + h*h*dtx2*st*st*st/(3.*qdp) + h*h*h*dtx3*st*st*st*st/(8.*qdp);
        
        J[1][1] = 1;
        J[1][2] = h*dty1dtx*st*st/2. + dC*tx*dty1*st*st/2./h + h*h*dty2dtx*st*st*st/6. + dC*tx*dty2*st*st*st/3.
        + h*h*h*dty3dtx*st*st*st*st/24. + dC*tx*h*dty3*st*st*st*st/12.;
        J[1][3] = st + h*dty1dty*st*st/2. + dC*ty*dty1*st*st/2./h + h*h*dty2dty*st*st*st/6. + dC*ty*dty2*st*st*st/3.
        + h*h*h*dty3dty*st*st*st*st/24. + dC*ty*h*dty3*st*st*st*st/12.;
        J[1][4] = h*dty1*st*st/(2.*qdp) + h*h*dty2*st*st*st/(3.*qdp) + h*h*h*dty3*st*st*st*st/(8.*qdp);
        
        J[2][2] = 1 + h*dtx1dtx*st + dC*tx*dtx1*st/h + h*h*dtx2dtx*st*st/2. + dC*tx*dtx2*st*st
        + h*h*h*dtx3dtx*st*st*st/6. + dC*tx*h*dtx3*st*st*st/2.;
        J[2][3] = h*dtx1dty*st + dC*ty*dtx1*st/h + h*h*dtx2dty*st*st/2. + dC*ty*dtx2*st*st
        + h*h*h*dtx3dty*st*st*st/6. + dC*ty*h*dtx3*st*st*st/2.;
        J[2][4] = h*dtx1*st/qdp + h*h*dtx2*st*st/qdp + h*h*h*dtx3*st*st*st/(2.*qdp);
        
        J[3][2] = h*dty1dtx*st + dC*tx*dty1*st/h + h*h*dty2dtx*st*st/2. + dC*tx*dty2*st*st
        + h*h*h*dty3dtx*st*st*st/6. + dC*tx*h*dty3*st*st*st/2.;
        J[3][3] = 1 + h*dty1dty*st + dC*ty*dty1*st/h + h*h*dty2dty*st*st/2. + dC*ty*dty2*st*st
        + h*h*h*dty3dty*st*st*st/6. + dC*ty*h*dty3*st*st*st/2.;
        J[3][4] = h*dty1*st/qdp + h*h*dty2*st*st/qdp + h*h*h*dty3*st*st*st/(2.*qdp);
        
        J[4][4] = 1;
        
        //cout << J[0][0] << " " << J[0][1] << " " << J[0][2] << " " << J[0][3] << " " <<  J[0][4]  << endl;
        /*cout << "************************************************************************" <<endl;
         cout << fixed << setprecision(9);
         for(int k = 0; k < 5; k++){
         
         for(int l = 0; l < 5; l++){
         cout  << J[k][l] << "  ";
         }
         cout << endl;
         }*/
        
        TransformCov(J, fStopPar.GetCov());
        //fStopPar.Print();
        AddNoise(fStopPar, fabs(st));
      }
      
    }
    //SetCovarianceMatrix(cov, fStopPar);
    
    /*if(fabs(fStopPar.GetZ()-zStop) > 0){
     fStopPar.SetPar(0, fStopPar.GetPar(0)+fStopPar.GetPar(2)*(fStopPar.GetZ()-zStop));
     fStopPar.SetPar(1, fStopPar.GetPar(1)+fStopPar.GetPar(3)*(fStopPar.GetZ()-zStop));
     } */
  }
  
  
  
  
  
  
  
  
  void TrackExtrapolation::Extrapolate(double zStop, double st){
    fStep = st;
    Extrapolate(zStop);
  }
  
  void TrackExtrapolation::Extrapolate(double zStop, double st, double radLength){
    fStep = st;
    fRadLength = radLength;
    Extrapolate(zStop);
  }
  
  
  
  void TrackExtrapolation::AddNoise(TrackPar &trackPar, double X){
    
    //return;
    if (fRadLength == 0) return;
    
    double len = X/fRadLength;
    double SigTheta = 0.0136*fabs(trackPar.GetPar(4) * sqrt(len) * (1.+0.038*log(len)));
    
    double p3 = trackPar.GetPar(2);
    double p4 = trackPar.GetPar(3);
    
    double p3p3 = SigTheta*SigTheta * (1 + p3*p3) * (1 + p3*p3 + p4*p4);
    double p4p4 = SigTheta*SigTheta * (1 + p4*p4) * (1 + p3*p3 + p4*p4);
    double p3p4 = SigTheta*SigTheta * p3*p4       * (1 + p3*p3 + p4*p4);
    
    trackPar.SetCov(0, trackPar.GetCov(0)+fStep*fStep*p3p3);
    trackPar.SetCov(1, trackPar.GetCov(1)+fStep*fStep*p3p4);
    trackPar.SetCov(2, trackPar.GetCov(2)-fStep*p3p3);
    trackPar.SetCov(3, trackPar.GetCov(3)-fStep*p3p4);
    
    trackPar.SetCov(5, trackPar.GetCov(5)+fStep*fStep*p4p4);
    trackPar.SetCov(6, trackPar.GetCov(6)-fStep*p3p4);
    trackPar.SetCov(7, trackPar.GetCov(7)-fStep*p4p4);
    
    trackPar.SetCov(9, trackPar.GetCov(9)+p3p3);
    trackPar.SetCov(10, trackPar.GetCov(10)+p3p4);
    trackPar.SetCov(12, trackPar.GetCov(12)+p4p4);
    
  }

