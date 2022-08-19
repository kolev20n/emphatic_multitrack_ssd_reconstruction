// This file contains the minimum necessary to connect the general
// reconstruction code to Matej Pavin's single-track code.
// This code was taken from Matej's version with minimal editing.
//
// Blame: Nikolay Kolev, kolev20n@uregina.ca
//
// 2022
//

#ifndef matej_minimum_data_structures_and_functions_h
#define matej_minimum_data_structures_and_functions_h

class MagneticField
{
public:
  MagneticField(std::string &name);
  ~MagneticField();
  
  void GetField(double x[3], double &Bx, double &By, double &Bz) ;
  
private:
  std::map<int, std::map<int, std::map<int, std::vector<double> > > > field;
  double step;
  double start[3];
  
};

enum TrackingType
{
  eUnknown,  // It's not known! +
    // (p1,    p2,    p3,    p4,    p5)          +
  eCartesian,  // (x,     y,     px,   py,   pz)          +
  eNA61,    // (x,     y,     q/pxz,   py/pxz, 2atan((pxz-px)/pxz)  +
  eKisel     // (x,     y,     px/pz,   py/pz,   q/p)        +
};

class TrackPar
{
public:
  TrackPar();
  TrackPar(TrackingType type);
  TrackPar(TrackingType type, double z, double *par);
  TrackPar(TrackingType type, double z, double *par, double *cov);
  TrackPar(TrackingType type, int charge, double z, double *par, double *cov);
    
  TrackPar& operator= (TrackPar &parIn);
  TrackingType GetType(){return fType;}
  int GetCharge(){return fCharge;}
  double GetPar(int i);
  double GetCov(int i);
  double* GetCov(){return fCov;}
  double GetStd(int i);
  double GetZ(){return fZ;}
  double GetLength(){return fLength;}
    
  void SetType(TrackingType type){fType = type;}
  void SetCharge(int charge){fCharge = charge;}
  void SetPar(int i, double val);
  void SetPar(double *par);
  void SetCov(int i, double val);
  void SetCov(double *cov);
  void SetZ(double val){fZ = val;}
  void SetLength(double length){fLength = length;}
  void AddLength(double d){fLength += d;}
  void Print();

private:
  //void ParameterCheck(vector<double> &par);
  //void CovarianceCheck(vector<double> &cov);
  void ChargeCheck(TrackingType type, double *par);
  //void CheckIndex(int i, vector<double> vec);
    
  TrackingType fType;
  int fCharge;
    
  double fZ;
  
  double fPar[5];
  double fCov[15];
  
  double fLength;
};

void TransformCov(double (*J)[5], double *c);
void ConvertTrackPar(TrackPar &parIn, TrackingType type);



class TrackExtrapolation
{
public:
  TrackExtrapolation(bool cov, double step, MagneticField *a_magnetic_field, double radLength=0);
  TrackExtrapolation(bool cov, TrackPar trPar,  double step, MagneticField *a_magnetic_field, double radLength=0);
  ~TrackExtrapolation();
  //**********************************Get functions********************************
  double GetStep(){return fStep;}
  TrackPar& GetStartTrackParam(){return fStartPar;}
  TrackPar& GetStopTrackParam(){return fStopPar;}
  
  static double GetKappa(){return fKappa;}
  bool GetErrorEstimation(){return fErrorEstimation;}
  //**********************************Set functions********************************
  void SetTrackPar(TrackPar &par){fStartPar = par; fStopPar = par;}
  //void SetClusters(const evt::rec::Track& track, const evt::RecEvent& recEvent);
  void SetRadLength(double val){fRadLength = val;}
  void SetStep(double step){fStep = step;}
  void SetErrorEstimation(bool cov){fErrorEstimation = cov;};
  //**********************************Other functions******************************
  
  void Extrapolate(double zStop, double step, double radLength);
  void Extrapolate(double zStop, double radLength);
  void Extrapolate(double zStop);
  void ExtrapolateToPlane(const double &a, const double&b, const double& c, const double& d);
  //void ExtrapolateToPlane(double radLength, Plane& plane);
  //void ExtrapolateToPlane(Plane& plane);
  // void DoKalmanStep(Cluster &cluster, Plane &plane);
  
  
private:
  
  
  void AddNoise(TrackPar &trackPar, double X);
  
  
  TrackPar fStartPar;
  TrackPar fStopPar;
  
  bool fErrorEstimation;
  
  double fStep;
  double fLength;
  
  MagneticField *fMagField;
  
  static double fKappa;
  double fRadLength;
  
  //help variables;
  
};

#endif
