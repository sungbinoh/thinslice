#include "ThinSlice.h"
#include "TGraphErrors.h"
#include "TVector3.h"
#include "RooUnfoldBayes.h"
#include "RooUnfoldSvd.h"
#include "util.h"

#include <iostream>

void ThinSlice::BookHistograms(){

  outputFile = TFile::Open(fOutputFileName.c_str(), "recreate");
  
  for (int i = 0; i<nthinslices; ++i){
    reco_incE[i] = new TH1D(Form("reco_incE_%d",i),Form("Reco incident energy, %.1f < z < %.1f (cm)",i*thinslicewidth, (i+1)*thinslicewidth), nbinse, 0, 1200.);
    true_incE[i] = new TH1D(Form("true_incE_%d",i),Form("True incident energy, %.1f < z < %.1f (cm)",i*thinslicewidth, (i+1)*thinslicewidth), nbinse, 0, 1200.);
    reco_incE[i]->Sumw2();
    true_incE[i]->Sumw2();
  }

  reco_AngCorr = new TH1D("reco_AngCorr","Reco angle correction", 100, 0, 1.);
  true_AngCorr = new TH1D("true_AngCorr","true angle correction", 100, 0, 1.);
  reco_AngCorr->Sumw2();
  true_AngCorr->Sumw2();

  h_truesliceid_pion_all = new TH1D("h_truesliceid_pion_all","h_truesliceid_pion_all;True SliceID", nthinslices + 2, -1, nthinslices + 1);
  h_truesliceid_pion_cuts = new TH1D("h_truesliceid_pion_cuts","h_truesliceid_pion_cuts;True SliceID", nthinslices + 2, -1, nthinslices + 1);
  h_truesliceid_pioninelastic_all = new TH1D("h_truesliceid_pioninelastic_all","h_truesliceid_pioninelastic_all;True SliceID", nthinslices + 2, -1, nthinslices + 1);
  h_truesliceid_pioninelastic_cuts = new TH1D("h_truesliceid_pioninelastic_cuts","h_truesliceid_pioninelastic_cuts;True SliceID", nthinslices + 2, -1, nthinslices + 1);
  h_recosliceid_allevts_cuts = new TH1D("h_recosliceid_allevts_cuts","h_recosliceid_allevts_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
  h_recosliceid_pion_cuts = new TH1D("h_recosliceid_pion_cuts","h_recosliceid_pion_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);
  h_recosliceid_pioninelastic_cuts = new TH1D("h_recosliceid_pioninelastic_cuts","h_recosliceid_pioninelastic_cuts;Reco SliceID", nthinslices + 2, -1, nthinslices + 1);

  h_truesliceid_pion_all->Sumw2();
  h_truesliceid_pion_cuts->Sumw2();
  h_truesliceid_pioninelastic_all->Sumw2();
  h_truesliceid_pioninelastic_cuts->Sumw2();
  h_recosliceid_allevts_cuts->Sumw2();
  h_recosliceid_pion_cuts->Sumw2();
  h_recosliceid_pioninelastic_cuts->Sumw2();

  for (int i = 0; i < nCuts; ++i){
    for (int j = 0; j < nParTypes+1; ++j){
      htrue_beam_endZ[i][j] = new TH1D(Form("htrue_beam_endZ_%d_%d",i,j),Form("true_beam_endZ, %s, %s;true_beam_endZ (cm)", cutName[i], parTypeName[j]), 70, -100, 600);
      htrue_beam_endZ[i][j]->Sumw2();
      hreco_beam_endZ[i][j] = new TH1D(Form("hreco_beam_endZ_%d_%d",i,j),Form("reco_beam_endZ, %s, %s;reco_beam_endZ (cm)", cutName[i], parTypeName[j]), 70, -100, 600);
      hreco_beam_endZ[i][j]->Sumw2();
      hreco_true_beam_endZ[i][j] = new TH1D(Form("hreco_true_beam_endZ_%d_%d",i,j), Form("reco_true_beam_endZ, %s, %s;reco_beam_endZ - true_beam_endZ (cm)", cutName[i], parTypeName[j]), 100, -100, 100);
      hreco_true_beam_endZ[i][j]->Sumw2();
      hreco_vs_true_beam_endZ[i][j]= new TH2D(Form("hreco_vs_true_beam_endZ_%d_%d",i,j), Form("%s, %s;true_beam_endZ (cm);reco_beam_endZ (cm)", cutName[i], parTypeName[j]), 70, -100, 600, 70, -100, 600);
      hreco_true_vs_true_beam_endZ[i][j]= new TH2D(Form("hreco_true_vs_true_beam_endZ_%d_%d",i,j), Form("%s, %s;true_beam_endZ (cm);reco - true_beam_endZ (cm)", cutName[i], parTypeName[j]), 70, -100, 600, 100, -100, 100);

      htrue_beam_endZ_SCE[i][j] = new TH1D(Form("htrue_beam_endZ_SCE_%d_%d",i,j),Form("true_beam_endZ_SCE, %s, %s;true_beam_endZ_SCE (cm)", cutName[i], parTypeName[j]), 70, -100, 600);
      htrue_beam_endZ_SCE[i][j]->Sumw2();
      hreco_beam_endZ_SCE[i][j] = new TH1D(Form("hreco_beam_endZ_SCE_%d_%d",i,j),Form("reco_beam_endZ_SCE, %s, %s;reco_beam_endZ_SCE (cm)", cutName[i], parTypeName[j]), 70, -100, 600);
      hreco_beam_endZ_SCE[i][j]->Sumw2();
      hreco_true_beam_endZ_SCE[i][j] = new TH1D(Form("hreco_true_beam_endZ_SCE_%d_%d",i,j), Form("reco_true_beam_endZ_SCE, %s, %s;reco_beam_endZ_SCE - true_beam_endZ_SCE (cm)", cutName[i], parTypeName[j]), 100, -100, 100);
      hreco_true_beam_endZ_SCE[i][j]->Sumw2();
      hreco_vs_true_beam_endZ_SCE[i][j]= new TH2D(Form("hreco_vs_true_beam_endZ_SCE_%d_%d",i,j), Form("%s, %s;true_beam_endZ_SCE (cm);reco_beam_endZ_SCE (cm)", cutName[i], parTypeName[j]), 70, -100, 600, 70, -100, 600);
      hreco_true_vs_true_beam_endZ_SCE[i][j]= new TH2D(Form("hreco_true_vs_true_beam_endZ_SCE_%d_%d",i,j), Form("%s, %s;true_beam_endZ_SCE (cm);reco - true_beam_endZ_SCE (cm)", cutName[i], parTypeName[j]), 70, -100, 600, 100, -100, 100);


      htrue_sliceID[i][j] = new TH1D(Form("htrue_sliceID_%d_%d",i,j),Form("true_sliceID, %s, %s;true_sliceID (cm)", cutName[i], parTypeName[j]), 26, -1, 25);
      htrue_sliceID[i][j]->Sumw2();
      hreco_sliceID[i][j] = new TH1D(Form("hreco_sliceID_%d_%d",i,j),Form("reco_sliceID, %s, %s;reco_sliceID", cutName[i], parTypeName[j]), 26, -1, 25);
      hreco_sliceID[i][j]->Sumw2();
      hreco_true_sliceID[i][j] = new TH1D(Form("hreco_true_sliceID_%d_%d",i,j), Form("reco_true_sliceID, %s, %s;reco_sliceID - true_sliceID", cutName[i], parTypeName[j]), 20, -10, 10);
      hreco_true_sliceID[i][j]->Sumw2();
      hreco_vs_true_sliceID[i][j]= new TH2D(Form("hreco_vs_true_sliceID_%d_%d",i,j), Form("%s, %s;true_sliceID;reco_sliceID", cutName[i], parTypeName[j]), 26, -1, 25, 26, -1, 25);
      hreco_true_vs_true_sliceID[i][j]= new TH2D(Form("hreco_true_vs_true_sliceID_%d_%d",i,j), Form("%s, %s;true_sliceID;reco_sliceID - true_sliceID", cutName[i], parTypeName[j]), 26, -1, 25, 20, -10, 10);

      hmediandEdx[i][j] = new TH1D(Form("hmediandEdx_%d_%d",i,j), Form("mediandEdx, %s, %s;Median dE/dx (MeV/cm)", cutName[i], parTypeName[j]), 100, 0, 5);
      hmediandEdx[i][j]->Sumw2();

      hdaughter_michel_score[i][j] = new TH1D(Form("hdaughter_michel_score_%d_%d",i,j), Form("daughter_michel_score, %s, %s;Michel score", cutName[i], parTypeName[j]), 110, -0.1, 1);
      hdaughter_michel_score[i][j]->Sumw2();

      htrackscore[i][j] = new TH1D(Form("htrackscore_%d_%d",i,j), Form("trackscore, %s, %s;Track score", cutName[i], parTypeName[j]), 110, -0.1, 1);
      htrackscore[i][j]->Sumw2();

      hdEdx_5cm[i][j] = new TH1D(Form("hdEdx_5cm_%d_%d",i,j), Form("dEdx_5cm, %s, %s;dE/dx_5cm (MeV/cm)", cutName[i], parTypeName[j]), 100, 0, 5);
      hdEdx_5cm[i][j]->Sumw2();

      hdeltax[i][j] = new TH1D(Form("hdeltax_%d_%d",i,j), Form("deltax, %s, %s;#Deltax/#sigma_{x}", cutName[i], parTypeName[j]), 100, -10, 10);
      hdeltax[i][j]->Sumw2();

      hdeltay[i][j] = new TH1D(Form("hdeltay_%d_%d",i,j), Form("deltay, %s, %s;#Deltay/#sigma_{y}", cutName[i], parTypeName[j]), 100, -10, 10);
      hdeltay[i][j]->Sumw2();

      hdeltaz[i][j] = new TH1D(Form("hdeltaz_%d_%d",i,j), Form("deltaz, %s, %s;#Deltaz/#sigma_{z}", cutName[i], parTypeName[j]), 100, -10, 10);
      hdeltaz[i][j]->Sumw2();

      hcostheta[i][j] = new TH1D(Form("hcostheta_%d_%d",i,j), Form("costheta, %s, %s;cos#theta", cutName[i], parTypeName[j]), 100, 0.9, 1);
      hcostheta[i][j]->Sumw2();

      htrklen[i][j] = new TH1D(Form("htrklen_%d_%d",i,j), Form("trklen, %s, %s;Track length (cm)", cutName[i], parTypeName[j]), 61, -10, 600);
      htrklen[i][j]->Sumw2();

      hreco_beam_startX_SCE[i][j] = new TH1D(Form("hreco_beam_startX_SCE_%d_%d",i,j), Form("reco_beam_startX_SCE, %s, %s; reco_beam_startX_SCE (cm)", cutName[i], parTypeName[j]), 100, -80, 20);
      hreco_beam_startX_SCE[i][j]->Sumw2();

      hreco_beam_startY_SCE[i][j] = new TH1D(Form("hreco_beam_startY_SCE_%d_%d",i,j), Form("reco_beam_startY_SCE, %s, %s; reco_beam_startY_SCE (cm)", cutName[i], parTypeName[j]), 100, 350, 500);
      hreco_beam_startY_SCE[i][j]->Sumw2();

      hreco_beam_startZ_SCE[i][j] = new TH1D(Form("hreco_beam_startZ_SCE_%d_%d",i,j), Form("reco_beam_startZ_SCE, %s, %s; reco_beam_startZ_SCE (cm)", cutName[i], parTypeName[j]), 100, -5, 10);
      hreco_beam_startZ_SCE[i][j]->Sumw2();

      hreco_beam_dcosX_SCE[i][j] = new TH1D(Form("hreco_beam_dcosX_SCE_%d_%d",i,j), Form("hreco_beam_dcosX_SCE, %s, %s; reco_beam_dcosX_SCE", cutName[i], parTypeName[j]), 100, -1, 1);
      hreco_beam_dcosX_SCE[i][j]->Sumw2();

      hreco_beam_dcosY_SCE[i][j] = new TH1D(Form("hreco_beam_dcosY_SCE_%d_%d",i,j), Form("hreco_beam_dcosY_SCE, %s, %s; reco_beam_dcosY_SCE", cutName[i], parTypeName[j]), 100, -1, 1);
      hreco_beam_dcosY_SCE[i][j]->Sumw2();

      hreco_beam_dcosZ_SCE[i][j] = new TH1D(Form("hreco_beam_dcosZ_SCE_%d_%d",i,j), Form("hreco_beam_dcosZ_SCE, %s, %s; reco_beam_dcosZ_SCE", cutName[i], parTypeName[j]), 100, -1, 1);
      hreco_beam_dcosZ_SCE[i][j]->Sumw2();

      hreco_beam_angleX_SCE[i][j] = new TH1D(Form("hreco_beam_angleX_SCE_%d_%d",i,j), Form("hreco_beam_angleX_SCE, %s, %s; #theta_{x} (deg)", cutName[i], parTypeName[j]), 180, 0, 180);
      hreco_beam_angleX_SCE[i][j]->Sumw2();

      hreco_beam_angleY_SCE[i][j] = new TH1D(Form("hreco_beam_angleY_SCE_%d_%d",i,j), Form("hreco_beam_angleY_SCE, %s, %s; #theta_{y} (deg)", cutName[i], parTypeName[j]), 180, 0, 180);
      hreco_beam_angleY_SCE[i][j]->Sumw2();

      hreco_beam_angleZ_SCE[i][j] = new TH1D(Form("hreco_beam_angleZ_SCE_%d_%d",i,j), Form("hreco_beam_angleZ_SCE, %s, %s; #theta_{z} (deg)", cutName[i], parTypeName[j]), 90, 0, 90);
      hreco_beam_angleZ_SCE[i][j]->Sumw2();

      hreco_beam_startXY_SCE[i][j] = new TH2D(Form("hreco_beam_startXY_SCE_%d_%d",i,j), Form("reco_beam_startXY_SCE, %s, %s;reco_beam_startX_SCE (cm);reco_beam_startY_SCE (cm)", cutName[i], parTypeName[j]), 1000, -360, 360, 1000, 0, 700);

    }
  }

   for (int i = 0; i<nthinslices; ++i){
     true_interactions[i] = 0;
     true_incidents[i] = 0;
   }
   
   //response_SliceID_Pion = new RooUnfoldResponse(nthinslices+2, -1, nthinslices+1, "response_SliceID_Pion");
   //response_SliceID_PionInEl = new RooUnfoldResponse(nthinslices+2, -1, nthinslices+1, "response_SliceID_PionInEl");

}

void ThinSlice::ProcessEvent(const HadAna & evt, Unfold & uf){

  reco_sliceID = -1;
  true_sliceID = -1;

  isTestSample = (evt.partype == kData);
  //if (evt.MC && evt.event%2 == 0) isTestSample = false;

  if (evt.MC){
    true_sliceID = int(evt.true_beam_endZ/thinslicewidth);
    if (true_sliceID < 0) true_sliceID = -1;
    if (evt.true_beam_endZ < 0) true_sliceID = -1;
    if (true_sliceID >= nthinslices) true_sliceID = nthinslices;
    if (evt.true_beam_PDG == 211){
      for (int i = 0; i<=true_sliceID; ++i){
        if (i<nthinslices) ++true_incidents[i];
      }
    }
    if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
      if (true_sliceID < nthinslices && true_sliceID>=0){
        ++true_interactions[true_sliceID];
      }
      // Reco info
      if (!(evt.reco_beam_calo_wire->empty()) && evt.reco_beam_true_byE_matched){
        std::vector<std::vector<double>> vincE(nthinslices);
        for (size_t i = 0; i<evt.reco_beam_calo_wire->size(); ++i){
          int this_sliceID = int((*evt.reco_beam_calo_Z)[i]/thinslicewidth);
          if (this_sliceID>=nthinslices) continue;
          if (this_sliceID<0) continue;
          double this_incE = (*evt.reco_beam_incidentEnergies)[i];
          vincE[this_sliceID].push_back(this_incE);
        }
        for (size_t i = 0; i<vincE.size(); ++i){
          if (!vincE[i].empty()){
            double sum_incE = 0;
            for (size_t j = 0; j<vincE[i].size(); ++j){
              sum_incE += vincE[i][j];
            }
            reco_incE[i]->Fill(sum_incE/vincE[i].size());
          }
        }
        TVector3 pt0(evt.reco_beam_calo_startX,
                     evt.reco_beam_calo_startY,
                     evt.reco_beam_calo_startZ);
        TVector3 pt1(evt.reco_beam_calo_endX,
                     evt.reco_beam_calo_endY,
                     evt.reco_beam_calo_endZ);
        TVector3 dir = pt1 - pt0;
        dir = dir.Unit();
        reco_AngCorr->Fill(dir.Z());
      }

      // True info
      if (!(evt.true_beam_traj_Z->empty())){
        std::vector<std::vector<double>> vincE(nthinslices);
        for (size_t i = 0; i<evt.true_beam_traj_Z->size()-1; ++i){//last point always has KE = 0
          int this_sliceID = int((*evt.true_beam_traj_Z)[i]/thinslicewidth);
          double this_incE = (*evt.true_beam_traj_KE)[i];
          if (this_sliceID>=nthinslices) continue;
          if (this_sliceID<0) continue;
          vincE[this_sliceID].push_back(this_incE);
        }
        for (size_t i = 0; i<vincE.size(); ++i){
          if (!vincE[i].empty()){
            double sum_incE = 0;
            for (size_t j = 0; j<vincE[i].size(); ++j){
              sum_incE += vincE[i][j];
            }
            true_incE[i]->Fill(sum_incE/vincE[i].size());
          }
        }
        TVector3 pt0(evt.true_beam_startX,
                     evt.true_beam_startY,
                     evt.true_beam_startZ);
        TVector3 pt1(evt.true_beam_endX,
                     evt.true_beam_endY,
                     evt.true_beam_endZ);
        TVector3 dir = pt1 - pt0;
        dir = dir.Unit();
        true_AngCorr->Fill(dir.Z());
      }
    }
  }

  if (!evt.reco_beam_calo_wire->empty()){
    reco_sliceID = int(evt.reco_beam_calo_endZ/thinslicewidth);
    if (reco_sliceID < 0) reco_sliceID = -1;
    if (evt.reco_beam_calo_endZ < 0) reco_sliceID = -1;
    if (reco_sliceID >= nthinslices) reco_sliceID = nthinslices;
  }

  if (evt.MC){
    if (evt.true_beam_PDG == 211){
      if (isTestSample){
        h_truesliceid_pion_all->Fill(true_sliceID);
      }
      else{
        uf.eff_den_Inc->Fill(true_sliceID);
      }
      if (evt.PassAllCuts() && evt.reco_beam_true_byE_matched){
        //if (evt.PassAllCuts()){
        if (isTestSample){
          h_recosliceid_pion_cuts->Fill(reco_sliceID);
          h_truesliceid_pion_cuts->Fill(true_sliceID);
        }
        else{
          uf.eff_num_Inc->Fill(true_sliceID);
          uf.pur_num_Inc->Fill(reco_sliceID);
          uf.response_SliceID_Inc.Fill(reco_sliceID, true_sliceID);
        }
      }
      else {
        if (!isTestSample){
          uf.response_SliceID_Inc.Miss(true_sliceID);
          //std::cout<<true_sliceID<<std::endl;
        }
      }
      
      if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
        if (isTestSample){
          h_truesliceid_pioninelastic_all->Fill(true_sliceID);
        }
        else{
          uf.eff_den_Int->Fill(true_sliceID);
        }
        if (evt.PassAllCuts() && evt.reco_beam_true_byE_matched){
          //if (evt.PassAllCuts()){
          if (isTestSample){
            h_recosliceid_pioninelastic_cuts->Fill(reco_sliceID);
            h_truesliceid_pioninelastic_cuts->Fill(true_sliceID);
          }
          else{
            uf.eff_num_Int->Fill(true_sliceID);
            uf.pur_num_Int->Fill(reco_sliceID);
            uf.response_SliceID_Int.Fill(reco_sliceID, true_sliceID);
          }
        }
        else{
          if (!isTestSample) uf.response_SliceID_Int.Miss(true_sliceID);
        }
      }
    }
    if (evt.PassAllCuts()){
      if (isTestSample){
        h_recosliceid_allevts_cuts->Fill(reco_sliceID);
      }
      else {
        uf.pur_den->Fill(reco_sliceID);
      }
    }
  }
}

void ThinSlice::FillHistograms(int cut, const HadAna & evt){
  
  if (cut>=0 && cut < nCuts){
    FillHistVec1D(htrue_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.partype);
    FillHistVec1D(htrue_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.partype);
    FillHistVec1D(htrue_sliceID[cut], true_sliceID, evt.partype);
    //    if (!evt.reco_beam_calo_wire->empty()){
    FillHistVec1D(hreco_beam_endZ[cut], evt.reco_beam_endZ, evt.partype);
    FillHistVec1D(hreco_true_beam_endZ[cut], evt.reco_beam_endZ - evt.true_beam_endZ_SCE, evt.partype);
    FillHistVec2D(hreco_vs_true_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.reco_beam_endZ, evt.partype);
    FillHistVec2D(hreco_true_vs_true_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.reco_beam_endZ - evt.true_beam_endZ_SCE, evt.partype);
    
    FillHistVec1D(hreco_beam_endZ_SCE[cut], evt.reco_beam_calo_endZ, evt.partype);
    FillHistVec1D(hreco_true_beam_endZ_SCE[cut], evt.reco_beam_calo_endZ - evt.true_beam_endZ, evt.partype);
    FillHistVec2D(hreco_vs_true_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.reco_beam_calo_endZ, evt.partype);
    FillHistVec2D(hreco_true_vs_true_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.reco_beam_calo_endZ - evt.true_beam_endZ, evt.partype);
    
    FillHistVec1D(hreco_sliceID[cut], reco_sliceID, evt.partype);
    FillHistVec1D(hreco_true_sliceID[cut], reco_sliceID - true_sliceID, evt.partype);
    FillHistVec2D(hreco_vs_true_sliceID[cut], true_sliceID, reco_sliceID, evt.partype);
    FillHistVec2D(hreco_true_vs_true_sliceID[cut], true_sliceID, reco_sliceID - true_sliceID, evt.partype);
    
    FillHistVec1D(hmediandEdx[cut], evt.median_dEdx, evt.partype);
    FillHistVec1D(hdaughter_michel_score[cut], evt.daughter_michel_score, evt.partype);
    FillHistVec1D(htrackscore[cut], evt.reco_beam_PFP_trackScore_collection, evt.partype);
    FillHistVec1D(hdEdx_5cm[cut], evt.dEdx_5cm, evt.partype);

    FillHistVec1D(hdeltax[cut], evt.beam_dx, evt.partype);
    FillHistVec1D(hdeltay[cut], evt.beam_dy, evt.partype);
    FillHistVec1D(hdeltaz[cut], evt.beam_dz, evt.partype);
    FillHistVec1D(hcostheta[cut], evt.beam_costh, evt.partype);

    FillHistVec1D(htrklen[cut], evt.reco_beam_alt_len, evt.partype);

    FillHistVec1D(hreco_beam_startX_SCE[cut], evt.reco_beam_calo_startX, evt.partype);
    FillHistVec1D(hreco_beam_startY_SCE[cut], evt.reco_beam_calo_startY, evt.partype);
    FillHistVec1D(hreco_beam_startZ_SCE[cut], evt.reco_beam_calo_startZ, evt.partype);

    if (!evt.reco_beam_calo_wire->empty()){    
      TVector3 pt0(evt.reco_beam_calo_startX,
                   evt.reco_beam_calo_startY,
                   evt.reco_beam_calo_startZ);
      TVector3 pt1(evt.reco_beam_calo_endX,
                   evt.reco_beam_calo_endY,
                   evt.reco_beam_calo_endZ);
      TVector3 dir = pt1 - pt0;
      dir = dir.Unit();
      FillHistVec1D(hreco_beam_dcosX_SCE[cut], dir.X(), evt.partype);
      FillHistVec1D(hreco_beam_dcosY_SCE[cut], dir.Y(), evt.partype);
      FillHistVec1D(hreco_beam_dcosZ_SCE[cut], dir.Z(), evt.partype);
      FillHistVec1D(hreco_beam_angleX_SCE[cut], acos(dir.X())*180/TMath::Pi(), evt.partype);
      FillHistVec1D(hreco_beam_angleY_SCE[cut], acos(dir.Y())*180/TMath::Pi(), evt.partype);
      FillHistVec1D(hreco_beam_angleZ_SCE[cut], acos(dir.Z())*180/TMath::Pi(), evt.partype);
    }

    FillHistVec2D(hreco_beam_startXY_SCE[cut], evt.reco_beam_calo_startX, evt.reco_beam_calo_startY, evt.partype);

  }
}

void ThinSlice::SaveHistograms(){
  outputFile->cd();
  outputFile->Write();
  h_truesliceid_pion_uf->Write("h_truesliceid_pion_uf");
  h_truesliceid_pioninelastic_uf->Write("h_truesliceid_pioninelastic_uf");
  //response_SliceID_Pion->Write("response_SliceID_Pion");
  //response_SliceID_PionInEl->Write("response_SliceID_PionInEl");
}

void ThinSlice::CalcXS(const Unfold & uf){

  double slcid[nthinslices] = {0};
  double avg_trueincE[nthinslices] = {0};
  double avg_recoincE[nthinslices] = {0};
  double err_trueincE[nthinslices] = {0};
  double err_recoincE[nthinslices] = {0};
  double reco_trueincE[nthinslices] = {0};
  double err_reco_trueincE[nthinslices] = {0};
  double truexs[nthinslices] = {0};
  double err_truexs[nthinslices] = {0};

  double NA=6.02214076e23;
  double MAr=39.95; //gmol
  double Density = 1.4; // g/cm^3

  for (int i = 0; i<nthinslices; ++i){
    
    slcid[i] = i;
    avg_trueincE[i] = true_incE[i]->GetMean();
    err_trueincE[i] = true_incE[i]->GetMeanError();
    avg_recoincE[i] = reco_incE[i]->GetMean();
    err_recoincE[i] = reco_incE[i]->GetMeanError();
    reco_trueincE[i] = avg_recoincE[i] - avg_trueincE[i];
    err_reco_trueincE[i] = sqrt(pow(err_trueincE[i],2)+pow(err_recoincE[i],2));
    //std::cout<<i<<" "<<avg_trueincE[i]<<std::endl;
    if (true_incidents[i] && true_interactions[i]){
      truexs[i] = MAr/(Density*NA*thinslicewidth/true_AngCorr->GetMean())*log(true_incidents[i]/(true_incidents[i]-true_interactions[i]))*1e27;
      err_truexs[i] = MAr/(Density*NA*thinslicewidth/true_AngCorr->GetMean())*1e27*sqrt(true_interactions[i]+pow(true_interactions[i],2)/true_incidents[i])/true_incidents[i];
    }
  }

  TGraphErrors *gr_trueincE = new TGraphErrors(nthinslices, &(slcid[0]), &(avg_trueincE[0]), 0, &(err_trueincE[0]));
  TGraphErrors *gr_recoincE = new TGraphErrors(nthinslices, &(slcid[0]), &(avg_recoincE[0]), 0, &(err_recoincE[0]));
  TGraphErrors *gr_reco_trueincE = new TGraphErrors(nthinslices, &(slcid[0]), &(reco_trueincE[0]), 0, &(err_reco_trueincE[0]));

  gr_trueincE->Write("gr_trueincE");
  gr_recoincE->Write("gr_recoincE");
  gr_reco_trueincE->Write("gr_reco_trueincE");

  TGraphErrors *gr_truexs = new TGraphErrors(nthinslices, &(avg_trueincE[0]), &(truexs[0]), 0, &(err_truexs[0]));
  
  gr_truexs->Write("gr_truexs");

  TH1D *hinc = (TH1D*)h_recosliceid_allevts_cuts->Clone("hinc");
  TH1D *hint = (TH1D*)h_recosliceid_allevts_cuts->Clone("hint");
  hinc->Multiply(uf.pur_Inc);
  hint->Multiply(uf.pur_Int);

//  RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, uf.pur_num_Inc, 4);
//  RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, uf.pur_num_Int, 4);

  RooUnfoldBayes   unfold_Inc (&uf.response_SliceID_Inc, hinc, 4);
  RooUnfoldBayes   unfold_Int (&uf.response_SliceID_Int, hint, 4);

//  RooUnfoldSvd     unfold_Inc (&uf.response_SliceID_Inc, hinc, 20);   // OR
//  RooUnfoldSvd     unfold_Int (&uf.response_SliceID_Int, hint, 20);   // OR

//  RooUnfoldSvd     unfold_Inc (&uf.response_SliceID_Inc, uf.pur_num_Inc, 20);   // OR
//  RooUnfoldSvd     unfold_Int (&uf.response_SliceID_Int, uf.pur_num_Int, 20);   // OR

  h_truesliceid_pion_uf = (TH1D*) unfold_Inc.Hreco();
  h_truesliceid_pioninelastic_uf = (TH1D*) unfold_Int.Hreco();

}

void ThinSlice::Run(HadAna & evt, Unfold & uf){

  BookHistograms();

  //Long64_t nentries = evt.fChain->GetEntriesFast();
  Long64_t nentries = evt.fChain->GetEntries();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%10000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = evt.LoadTree(jentry);
    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //std::cout<<evt.run<<" "<<evt.event<<" "<<evt.MC<<" "<<evt.reco_beam_true_byE_matched<<" "<<evt.true_beam_PDG<<" "<<(*evt.true_beam_endProcess)<<std::endl;
    //std::cout<<GetParType(ana)<<std::endl;
    if (!evt.isSelectedPart()) continue;
    evt.ProcessEvent();
    ProcessEvent(evt, uf);
    FillHistograms(kNocut, evt);
    if (evt.PassPandoraSliceCut()){
      FillHistograms(kPandoraSlice, evt);
      if (evt.PassCaloSizeCut()){
        FillHistograms(kCaloSize, evt);
        if (evt.PassBeamQualityCut()){
          FillHistograms(kBeamQuality, evt);
          if (evt.PassAPA3Cut()){
            FillHistograms(kAPA3, evt);
            /*
            if (evt.GetParType() == kPiInel && evt.daughter_michel_score > 0.99){
              std::cout<<evt.run<<" "<<evt.subrun<<" "<<evt.event<<std::endl;
              if (!evt.reco_daughter_PFP_michelScore_collection->empty()){
                for (size_t i = 0; i<evt.reco_daughter_PFP_michelScore_collection->size(); ++i){
                  std::cout<<i<<" "<<(*evt.reco_daughter_PFP_michelScore_collection)[i]<<" "<<(*evt.reco_daughter_PFP_nHits_collection)[i]<<std::endl;
                }
              }
            }
            */
            if (evt.PassMichelScoreCut()){
              FillHistograms(kMichelScore, evt);
              if (evt.PassMediandEdxCut()){
                FillHistograms(kMediandEdx, evt);
              }
            }
          }
        }
      }
    }
  }
  
  uf.SaveHistograms();
  CalcXS(uf);
  SaveHistograms();
}
