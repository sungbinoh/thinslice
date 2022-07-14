#include "BetheBloch.h"
#include "TSpline.h"
#include <iostream>
#include <cmath>
#include <algorithm>

ROOT::Math::VavilovAccurate vav;

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
  // == Estimate the density correction
  double density_y = TMath::Log10(beta * gamma);
  double ln10 = TMath::Log(10);
  double this_delta = 0.;
  if(density_y > density_y1){
    this_delta = 2.0 * ln10 * density_y - density_C;
  }
  else if (density_y < density_y0){
    this_delta = 0.;
  }
  else{
    this_delta = 2.0 * ln10 * density_y - density_C + density_a * pow(density_y1 - density_y, density_k);
  }

  return this_delta;
}

double BetheBloch::betaGamma(double KE){

  double gamma, beta;
  gamma = (KE + mass) / mass;
  beta = sqrt( 1 - 1/pow(gamma,2));
   
  return beta*gamma;

}

double BetheBloch::Landau_xi(double KE, double pitch){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double xi = rho * pitch * 0.5 * K * (Z / A) * pow(1. / beta, 2);
  return xi;
}

double BetheBloch::Get_Wmax(double KE){
  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double Wmax = (2.0 * me * pow(beta * gamma, 2)) / (1.0 + 2.0 * me * (gamma / mass) + pow((me / mass),2));

  return Wmax;
}

double BetheBloch::meandEdx(double KE){

  //KE is kinetic energy in MeV
  double gamma = (KE + mass) / mass;
  double beta = sqrt( 1 - 1/pow(gamma,2));
  
  double wmax = Get_Wmax(KE);

  //cout<<wmax<<" "<<rho*K*Z*pow(charge,2)<<" "<<A*pow(beta,2)<<" "<<0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2))<<" "<<pow(beta,2)<<" "<<densityEffect( beta, gamma )/2<<endl;
  double dEdX = (rho*K*Z*pow(charge,2))/(A*pow(beta,2))*(0.5*log(2*me*pow(gamma,2)*pow(beta,2)*wmax/pow(I,2)) - pow(beta,2) - densityEffect( beta, gamma )/2 );

  return dEdX;
}

double BetheBloch::MPVdEdx(double KE, double pitch){

  //KE is kinetic energy in MeV
  //pitch is in cm
  double gamma = (KE + mass) / mass;
  double beta = sqrt( 1 - 1/pow(gamma,2));

  double xi = Landau_xi(KE, pitch);
  
  double eloss_mpv = xi*(log( 2*me*pow(gamma,2)*pow(beta,2) / I ) + log( xi / I ) + 0.2 - pow(beta,2) - densityEffect( beta, gamma ) )/pitch;

  return eloss_mpv;
}

double BetheBloch::IntegratedEdx(double KE0, double KE1, int n){

  if (KE0>KE1) swap(KE0, KE1);

  double step = (KE1-KE0)/n;

  double area = 0;
  
  for (int i = 0; i<n; ++i){
    double dEdx = meandEdx(KE0 + (i+0.5)*step);
    if (dEdx)
      area += 1/dEdx*step;
  }
  return area;
}

double BetheBloch::RangeFromKE(double KE){

  return IntegratedEdx(0, KE);
}

void BetheBloch::CreateSplines(int np, double minke, double maxke){

  if (sp_KE_range) delete sp_KE_range;
  if (sp_range_KE) delete sp_range_KE;

  for (const auto & x : spmap){
    if (x.second) delete x.second;
  }
  spmap.clear();
  
  double *KE = new double[np];
  double *Range = new double[np];
//  vector<double> KE(np);
//  vector<double> Range(np);

  for (int i = 0; i<np; ++i){
    double ke = pow(10, log10(minke)+i*log10(maxke/minke)/np);
//    KE.push_back(ke);
//    Range.push_back(RangeFromKE(ke));
    KE[i] = ke;
    Range[i] = RangeFromKE(ke);
    //if (pdgcode ==2212) cout<<KE[i]<<" "<<Range[i]<<" "<<meandEdx(KE[i])<<endl;
  }

  sp_KE_range = new TSpline3("sp_KE_range", KE, Range, np, "b2e2", 0, 0);
  sp_range_KE = new TSpline3("sp_range_KE", Range, KE, np, "b2e2", 0, 0);
  //cout<<sp_KE_range->Eval(10)<<endl;
  delete[] KE;
  delete[] Range;
  cout<<"Done creating splines for particle with pdgcode "<<pdgcode<<endl;
}

double BetheBloch::RangeFromKESpline(double KE){
  if (!sp_KE_range){
    cout<<"Spline does not exist."<<endl;
    exit(1);
  }
  return sp_KE_range->Eval(KE);
}

double BetheBloch::KEFromRangeSpline(double range){
  if (!sp_range_KE){
    cout<<"Spline does not exit."<<endl;
    exit(1);
  }
  return sp_range_KE->Eval(range);
}

double BetheBloch::KEAtLength(double KE0, double tracklength){

  int iKE = int(KE0);

  if (spmap.find(iKE)==spmap.end()){
    CreateSplineAtKE(iKE);
  }

  double deltaE = spmap[iKE]->Eval(tracklength);

  if (deltaE < 0) return 0;//cout<<"Negative delta E: "<<deltaE<<endl;
  if (KE0 - deltaE < 0) return 0;//cout<<"Negative KE: "<<KE0 - deltaE<<endl;
  
  return KE0 - deltaE;

}

double BetheBloch::KEtoMomentum(double KE){
  return sqrt(pow(KE, 2) + 2.0 * KE * mass);
}

void BetheBloch::CreateSplineAtKE(int iKE){

  double KE0 = iKE;

  // Sample every 10 MeV
  int np = int(KE0/10);
  double *deltaE;
  double *trklength;
  if (np>1){
    deltaE = new double[np];
    trklength = new double[np];
    for (int i = 0; i<np; ++i){
      double KE = KE0 - i*10;
      deltaE[i] = KE0 - KE;
      trklength[i] = IntegratedEdx(KE, KE0);
      //cout<<iKE<<" "<<KE<<" "<<i<<"/"<<np<<" "<<trklength[i]<<endl;
    }
  }
  else{
    cout<<"KE too low: "<<iKE<<endl;
    np = 2;
    deltaE = new double[np];
    trklength = new double[np];
    deltaE[0] = 0;
    trklength[0] = 0;
    deltaE[1] = KE0;
    trklength[1] = RangeFromKE(KE0);
  }

  spmap[iKE] = new TSpline3(Form("KE %d",iKE), trklength, deltaE, np, "b2e2", 0, 0);
  delete[] trklength;
  delete[] deltaE;

}
double dEdx_PDF_fuction(double *x, double *par){
  // == par[5] = {kappa, beta^2, xi, <dE/dx>BB, width}
  double a = par[2] / par[4];
  double b = (0.422784 + par[1] + log(par[0])) * par[2] / par[4] + par[3];
  double y = (x[0] - b) / a;

  double this_vav = 0.;

  if(par[0] < 0.01){ // == Landau
    this_vav = TMath::Landau(y);
    this_vav =  this_vav / a;
  }
  else if(par[0] > 10.){ // == Gaussian
    double mu = vav.Mean(par[0], par[1]);
    double sigma = sqrt(vav.Variance(par[0], par[1]));
    this_vav =  TMath::Gaus(y, mu, sigma);
  }
  else{ // == Vavilov
    this_vav =  vav.Pdf(y, par[0], par[1]);
    this_vav =  this_vav / a;
  }

  return this_vav;
}


TF1 *BetheBloch::dEdx_PDF(double KE, double pitch){

  double gamma = (KE/mass)+1.0;
  double beta = TMath::Sqrt(1-(1.0/(gamma*gamma)));
  double this_xi = Landau_xi(KE, pitch);
  double this_Wmax = Get_Wmax(KE);
  double this_kappa = this_xi / this_Wmax;
  double this_dEdx_BB = meandEdx(KE);
  double par[5] = {this_kappa, beta * beta, this_xi, this_dEdx_BB, pitch};
  
  TF1 *out = new TF1("", dEdx_PDF_fuction, -100., 1000., 5);
  out -> SetParameters(par[0], par[1], par[2], par[3], par[4]);

  return out;
}
