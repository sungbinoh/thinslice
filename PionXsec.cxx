#include "PionXsec.h"
#include "HadAna.h"
#include "anavar.h"
#include "TH1D.h"
#include "util.h"
#include <iostream>

using namespace std;

PionXsec::PionXsec(){
  hadana.InitPi();
}

double PionXsec::Gaussian_Reweight(double mu1, double sigma1, double mu2, double sigma2, double x_min, double x_max, double x){
  // == Estimate y of gaussian1 / gaussian2
  TF1 *gaus1 = new TF1("gaus1", "gaus", x_min, x_max);
  gaus1 -> SetParameter(0, 1. / (sigma1 * sqrt(2.0 * TMath::Pi())));
  gaus1 -> SetParameter(1, mu1);
  gaus1 -> SetParameter(2, sigma1);

  TF1 *gaus2 = new TF1("gaus2","gaus",x_min, x_max);
  gaus2 -> SetParameter(0, 1. / (sigma2 * sqrt(2.0 * TMath::Pi())));
  gaus2 -> SetParameter(1, mu2);
  gaus2 -> SetParameter(2, sigma2);

  double y1 = gaus1 -> Eval(x);
  double y2 = gaus2 -> Eval(x);
  double out = y1/y2;
  //cout << "[PionXsec::Gaussian_Reweight] true beam P : " << x << ", y1 : " << y1 << ", y2 : " << y2 << ", y1/y2 : " << y1/y2 << endl;
  if(out > 10.) return 10.;
  return y1/y2;
}

void PionXsec::BookHistograms(){
  //cout << "[PionXsec::BookHistograms] Start" << endl;
  Hist.outfile = TFile::Open(fOutputFileName.c_str(), "recreate");
}

void PionXsec::ProcessEvent(const anavar & evt){
  hadana.ProcessEvent(evt);
}

void PionXsec::FillHistBeam(const anavar & evt, double weight, TString suffix){

  double beamP = evt.beam_inst_P*1000. * beamP_scale;
  double beamKE = sqrt(pow(beamP, 2) + pow(139.57, 2)) - 139.57;
  TString pitype_str = Form("%d", hadana.pitype);

  // ==== Beam
  // == Reco
  Hist.JSFillHist(suffix, "htrack_BeamP_" + suffix + "_" + pitype_str, beamP, weight, 5000., 0., 5000.);
  Hist.JSFillHist(suffix, "htrack_PandoraSlice_" + suffix + "_" + pitype_str, hadana.PassPandoraSliceCut(evt), weight, 2., 0., 2.);
  Hist.JSFillHist(suffix, "htrack_CaloSize_" + suffix + "_" + pitype_str, hadana.PassCaloSizeCut(evt),weight, 2., 0.,2.);
  Hist.JSFillHist(suffix, "htrack_beam_dx_" + suffix + "_" + pitype_str, hadana.beam_dx, weight, 200., -10., 10.);
  Hist.JSFillHist(suffix, "htrack_beam_dy_" + suffix + "_" + pitype_str, hadana.beam_dy, weight, 200., -10.,10.);
  Hist.JSFillHist(suffix, "htrack_beam_dz_" + suffix + "_" + pitype_str, hadana.beam_dz, weight, 200., -10.,10.);
  Hist.JSFillHist(suffix, "htrack_beam_dxy_" + suffix + "_" + pitype_str, hadana.beam_dxy, weight, 200., -10.,10.);
  Hist.JSFillHist(suffix, "htrack_beam_costh_" + suffix + "_" + pitype_str, hadana.beam_costh, weight, 200., -1.,1.);
  Hist.JSFillHist(suffix, "htrack_daughter_michel_score_" + suffix + "_" + pitype_str, hadana.daughter_michel_score, weight, 100., 0.,1.);
  Hist.JSFillHist(suffix, "htrack_chi2_proton_" + suffix + "_" + pitype_str, hadana.chi2_proton, weight, 1000., 0., 100.);
  Hist.JSFillHist(suffix, "htrack_beam_inst_XY_" + suffix + "_" + pitype_str, evt.beam_inst_X, evt.beam_inst_Y, weight, 400., -50., -10., 400., 400., 440.); 

  // == Truth
  if(evt.MC){
    double beamP_true = evt.true_beam_startP * 1000.;
    double true_KE_ff = hadana.true_ffKE;
    double true_P_ff = hadana.map_BB[abs(evt.true_beam_PDG)] -> KEtoMomentum(true_KE_ff);
    double beamP_res = (beamP_true - beamP) / beamP_true;
    Hist.JSFillHist(suffix, "htrack_BeamP_Res_" + suffix + "_" + pitype_str, beamP_res, weight, 4000, -2, 2);
    Hist.JSFillHist(suffix, "htrack_BeamP_true_" + suffix + "_" + pitype_str, beamP_true, weight, 5000, 0, 5000);
    Hist.JSFillHist(suffix, "htrack_KEffTruth_" + suffix + "_" + pitype_str, true_KE_ff, weight, 5000, 0, 5000);
    Hist.JSFillHist(suffix, "htrack_PffTruth_" + suffix + "_" + pitype_str, true_P_ff, weight, 5000, 0, 5000);
  }
}

void PionXsec::FillHistDaughters(const anavar & evt, double weight, TString suffix){
  
  // ==== Loop for Reco daughters
  for (size_t i = 0; i < evt.reco_daughter_allTrack_ID->size(); i++){
    if((*evt.reco_daughter_allTrack_dQdX_SCE)[i].empty()) continue;

    // ==== Truth level 
    if(!evt.MC) continue;
    // == proton
    if((*evt.reco_daughter_PFP_true_byHits_PDG).at(i) == 2212){
      Hist.JSFillHist(suffix, "hdaughter_proton_trackScore", (*evt.reco_daughter_PFP_trackScore).at(i), weight, 100., 0., 1.);
      Hist.JSFillHist(suffix, "hdaughter_proton_emScore", (*evt.reco_daughter_PFP_emScore).at(i), weight, 100., 0., 1.);
      double this_chi2 = (*evt.reco_daughter_allTrack_Chi2_proton).at(i) / (*evt.reco_daughter_allTrack_Chi2_ndof).at(i);
      Hist.JSFillHist(suffix, "hdaughter_proton_chi2_proton", this_chi2, weight, 1000., 0., 100.); 
      //cout << "proton this_chi2 : " << this_chi2 << endl;
    }
    // == pion
    else if(abs((*evt.reco_daughter_PFP_true_byHits_PDG).at(i)) == 211){
      Hist.JSFillHist(suffix, "hdaughter_pion_trackScore", (*evt.reco_daughter_PFP_trackScore).at(i), weight, 100., 0., 1.);
      Hist.JSFillHist(suffix, "hdaughter_pion_emScore", (*evt.reco_daughter_PFP_emScore).at(i), weight, 100., 0., 1.);
      double this_chi2 = (*evt.reco_daughter_allTrack_Chi2_proton).at(i) / (*evt.reco_daughter_allTrack_Chi2_ndof).at(i);
      Hist.JSFillHist(suffix, "hdaughter_pion_chi2_proton", this_chi2, weight, 1000., 0., 100.);
    }
    else{
      Hist.JSFillHist(suffix, "hdaughter_other_trackScore", (*evt.reco_daughter_PFP_trackScore).at(i), weight, 100., 0., 1.);
      Hist.JSFillHist(suffix, "hdaughter_other_emScore", (*evt.reco_daughter_PFP_emScore).at(i), weight, 100., 0., 1.);
      double this_chi2 = (*evt.reco_daughter_allTrack_Chi2_proton).at(i) / (*evt.reco_daughter_allTrack_Chi2_ndof).at(i);
      Hist.JSFillHist(suffix, "hdaughter_other_chi2_proton", this_chi2, weight, 1000., 0., 100.);
    }

    // == Truth matched
    for (size_t j = 0; j < evt.true_beam_daughter_ID->size(); j++){
      if (std::abs((*evt.true_beam_daughter_PDG)[j]) != 211 &&
          (*evt.true_beam_daughter_PDG)[j] != 2212) continue;

      if ((*evt.reco_daughter_PFP_true_byHits_ID)[i] == (*evt.true_beam_daughter_ID)[j]){

        double true_mom = (*evt.true_beam_daughter_startP)[j]*1000;
        double true_theta = GetTheta((*evt.true_beam_daughter_startPx)[j],
                                     (*evt.true_beam_daughter_startPy)[j],
                                     (*evt.true_beam_daughter_startPz)[j]);
        double true_phi = GetPhi((*evt.true_beam_daughter_startPx)[j],
                                 (*evt.true_beam_daughter_startPy)[j],
                                 (*evt.true_beam_daughter_startPz)[j]);
        bool foundreco = false;
      }

    }


  }


  // ==== Truth
  for (size_t i = 0; i<evt.true_beam_daughter_ID->size(); ++i){
    if (std::abs((*evt.true_beam_daughter_PDG)[i]) != 211 &&
        (*evt.true_beam_daughter_PDG)[i] != 2212) continue;

  }
  

}

void PionXsec::SaveHistograms(){
  
  outputFile->cd();
  outputFile->Write();
  Hist.WriteHist();
}

void PionXsec::Run(anavar & evt, Long64_t nentries=-1){

  BookHistograms();

  //Long64_t nentries = evt.fChain->GetEntriesFast();
  if (nentries == -1) nentries = evt.fChain->GetEntries();

  // == Beam P scale
  if(evt.MC) beamP_scale = 1.0;
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //if (jentry%100000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    if (jentry%1000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = evt.LoadTree(jentry);
    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //std::cout<<evt.run<<" "<<evt.event<<" "<<evt.MC<<" "<<evt.reco_beam_true_byE_matched<<" "<<evt.true_beam_PDG<<" "<<(*evt.true_beam_endProcess)<<std::endl;
    //std::cout<<GetParType(ana)<<std::endl;
    if (!hadana.isSelectedPart(evt)) continue;
    Hist.FillHist("Cutflow", 0.5, 1., 10., 0., 10.);
    ProcessEvent(evt);

    double weight = 1.;
    double weight_piOnly = 1.;
    if(evt.MC){
      double beamP_true = evt.true_beam_startP * 1000.;
      double data_mu = 1007.;
      double data_sigma = 68.17;
      double true_mu = 1007.;
      double true_sigma = 57.7;
      double P_reweight = Gaussian_Reweight(data_mu, data_sigma, true_mu, true_sigma, 0., 2000., beamP_true);
      weight = P_reweight;
      if(evt.true_beam_PDG == 211){
	weight_piOnly = P_reweight;
      }
    }
    FillHistBeam(evt, 1., "precut_noweight");
    FillHistBeam(evt, weight, "precut_Preweight");
    FillHistBeam(evt, weight_piOnly, "precut_Preweight_piOnly");

    if (!hadana.PassPiCuts(evt)) continue;
    FillHistBeam(evt, 1., "noweight");
    FillHistBeam(evt, weight, "Preweight");
    FillHistBeam(evt, weight_piOnly, "Preweight_piOnly");
    FillHistDaughters(evt, 1., "noweight");
    FillHistDaughters(evt, weight, "Preweight");
    FillHistDaughters(evt, weight_piOnly, "Preweight_piOnly");

    //cout << "(run:event)\t(" << evt.run << ":" << evt.event << ")\tdaughter_michel_score\t" << hadana.daughter_michel_score << endl;
  }
  SaveHistograms();
  //Make_dEdx_Range_Profile();
}

