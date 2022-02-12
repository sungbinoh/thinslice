#ifndef BETHEBLOCH_H
#define BETHEBLOCH_H

class BetheBloch {

 public:

  BetheBloch();
  BetheBloch(int pdg);

  double meandEdx(double KE);

  double MPVdEdx(double KE, double pitch);

  double RangeFromKE(double KE, int n= 10000);

 private:

  int pdgcode;
  double mass;
  int charge;

  void SetPdgCode(int pdg);

  double densityEffect(double beta, double gamma);

  double betaGamma(double KE);


};

#endif
