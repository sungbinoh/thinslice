#ifndef BETHEBLOCH_H
#define BETHEBLOCH_H

class TSpline3;

class BetheBloch {

 public:

  BetheBloch();
  BetheBloch(int pdg);

  void SetPdgCode(int pdg);

  int GetPdgCode(){ return pdgcode;};

  double meandEdx(double KE);

  double MPVdEdx(double KE, double pitch);

  double RangeFromKE(double KE, int n= 10000);

  double RangeFromKESpline(double KE);

  double KEFromRangeSpline(double range);

 private:

  int pdgcode;
  double mass;
  int charge;

  TSpline3 *sp_KE_range;
  TSpline3 *sp_range_KE;

  double densityEffect(double beta, double gamma);

  double betaGamma(double KE);

  void CreateSplines(int np = 1000, double minke = .01, double maxke = 2e5);

};

#endif
