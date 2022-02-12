#include "BetheBloch.h"
#include "TF1.h"
#include "TSpline.h"
#include <iostream>
#include <cmath>

using namespace std;

BetheBloch::BetheBloch()
  : pdgcode(0)
  , mass(0)
  , charge(0)
  , sp_KE_range(0)
  , sp_range_KE(0){
}

BetheBloch::BetheBloch(int pdg)
  : pdgcode(0)
  , mass(0)
  , charge(0)
  , sp_KE_range(0)
  , sp_range_KE(0){
  SetPdgCode(pdg);
}

void BetheBloch::SetPdgCode(int pdg){

  pdgcode = pdg;

  if (abs(pdgcode) == 13){//muon
    mass = 105.6583755;
    charge = 1;
  }
  else if (abs(pdgcode) == 211){//pion
    mass = 139.57039;
    charge = 1;
  }
  else if (abs(pdgcode) == 321){//kaon
    mass = 493.677;
    charge = 1;
  }
  else if (pdgcode == 2212){//proton
    mass = 938.27208816;
    charge = 1;
  }
  else{
    cout<<"Unknown pdg code "<<pdgcode<<endl;
    exit(1);
  }

  CreateSplines();

}

double BetheBloch::densityEffect(double beta, double gamma){

  double lar_C = 5.215, lar_x0 = 0.201, lar_x1 = 3, lar_a = 0.196, lar_k = 3;
  double x = log10(beta * gamma);
  
  if( x >= lar_x1 ){
    return 2*log(10)*x - lar_C;
  }
  else if ( lar_x0 <= x && x < lar_x1){
    return 2*log(10)*x - lar_C + lar_a * pow(( lar_x1 - x ) , lar_k );
  }
  else{
    return 0; //if x < lar_x0
  }
}

double BetheBloch::betaGamma(double KE){

  double gamma, beta;
  gamma = (KE + mass) / mass;
  beta = sqrt( 1 - 1/pow(gamma,2));
   
  return beta*gamma;

}

double BetheBloch::meandEdx(double KE){

  //KE is kinetic energy in MeV
  
  double K = 0.307;
  double rho = 1.396;
  double Z = 18;
  double A = 39.948;
  double I = pow(10,-6)*10.5*18; //MeV
  double me = 0.511; //MeV me*c^2
  
  double gamma = (KE + mass) / mass;
  double beta = sqrt( 1 - 1/pow(gamma,2));
  
  double wmax = 2*me*pow(beta,2)*pow(gamma,2)/(1+2*gamma*me/mass + pow(me,2)/pow(mass,2));
  
  double dEdX = (rho*K*Z*pow(charge,2))/(A*pow(beta,2))*(0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2)) - pow(beta,2) - densityEffect( beta, gamma )/2 );
  
  return dEdX;
}

double BetheBloch::MPVdEdx(double KE, double pitch){

  //KE is kinetic energy in MeV
  //pitch is in cm

  double K = 0.307;
  double rho = 1.396;
  double Z = 18;
  double A = 39.948;
  double I = pow(10,-6)*10.5*18; //MeV
  double  me = 0.511; //MeV me*c^2

  double gamma = (KE + mass) / mass;
  double beta = sqrt( 1 - 1/pow(gamma,2));

  double xi = ( K/2 )*( Z/A )* ( pitch * rho / pow(beta,2));

  double eloss_mpv = xi*(log( 2*me*pow(gamma,2)*pow(beta,2) / I ) + log( xi / I ) + 0.2 - pow(beta,2) - densityEffect( beta, gamma ) )/pitch;

  return eloss_mpv;
}

double BetheBloch::RangeFromKE(double KE, int n){

  double step = KE/n;

  double area = 0.;

  for (int i = 0; i<n; ++i){
    area += 1/meandEdx((i+0.5)*step)*step;
  }
  return area;

}

void BetheBloch::CreateSplines(int np, double minke, double maxke){

  if (sp_KE_range) delete sp_KE_range;
  if (sp_range_KE) delete sp_range_KE;

  double *KE = new double[np];
  double *Range = new double[np];
//  vector<double> KE(np);
//  vector<double> Range(np);

  for (int i = 0; i<np; ++i){
    double ke = pow(10, log10(minke)+i*log10(maxke/minke)/np);
//    KE.push_back(ke);
//    Range.push_back(RangeFromKE(ke));
//    cout<<KE.back()<<" "<<Range.back()<<endl;
    KE[i] = ke;
    Range[i] = RangeFromKE(ke);
  }

  sp_KE_range = new TSpline3("sp_KE_range", KE, Range, np, "b2e2", 0, 0);
  sp_range_KE = new TSpline3("sp_range_KE", Range, KE, np, "b2e2", 0, 0);
  //cout<<sp_KE_range->Eval(10)<<endl;
  delete[] KE;
  delete[] Range;
  cout<<"Done creating splines for particle with pdgcode "<<pdgcode<<endl;
}

double BetheBloch::RangeFromKESpline(double KE){
  return sp_KE_range->Eval(KE);
}

double BetheBloch::KEFromRangeSpline(double range){
  return sp_range_KE->Eval(range);
}
