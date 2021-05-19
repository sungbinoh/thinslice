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
    reco_AngCorr[i] = new TH1D(Form("reco_AngCorr_%d",i),Form("Reco angle correction, %.1f < z < %.1f (cm)",i*thinslicewidth, (i+1)*thinslicewidth), 100, 0, 1.);
    true_AngCorr[i] = new TH1D(Form("true_AngCorr_%d",i),Form("True angle correction, %.1f < z < %.1f (cm)",i*thinslicewidth, (i+1)*thinslicewidth), 100, 0, 1.);
    reco_incE[i]->Sumw2();
    true_incE[i]->Sumw2();
    reco_AngCorr[i]->Sumw2();
    true_AngCorr[i]->Sumw2();
  }
  
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


      htrue_sliceID[i][j] = new TH1D(Form("htrue_sliceID_%d_%d",i,j),Form("true_sliceID, %s, %s;true_sliceID (cm)", cutName[i], parTypeName[j]), 50, -1, 49);
      htrue_sliceID[i][j]->Sumw2();
      hreco_sliceID[i][j] = new TH1D(Form("hreco_sliceID_%d_%d",i,j),Form("reco_sliceID, %s, %s;reco_sliceID", cutName[i], parTypeName[j]), 50, -1, 49);
      hreco_sliceID[i][j]->Sumw2();
      hreco_true_sliceID[i][j] = new TH1D(Form("hreco_true_sliceID_%d_%d",i,j), Form("reco_true_sliceID, %s, %s;reco_sliceID - true_sliceID", cutName[i], parTypeName[j]), 20, -10, 10);
      hreco_true_sliceID[i][j]->Sumw2();
      hreco_vs_true_sliceID[i][j]= new TH2D(Form("hreco_vs_true_sliceID_%d_%d",i,j), Form("%s, %s;true_sliceID;reco_sliceID", cutName[i], parTypeName[j]), 50, -1, 49, 50, -1, 49);
      hreco_true_vs_true_sliceID[i][j]= new TH2D(Form("hreco_true_vs_true_sliceID_%d_%d",i,j), Form("%s, %s;true_sliceID;reco_sliceID - true_sliceID", cutName[i], parTypeName[j]), 50, -1, 49, 20, -10, 10);

      hmediandEdx[i][j] = new TH1D(Form("hmediandEdx_%d_%d",i,j), Form("mediandEdx, %s, %s;Median dE/dx (MeV/cm)", cutName[i], parTypeName[j]), 100, 0, 5);
      hmediandEdx[i][j]->Sumw2();

      hdaughter_michel_score[i][j] = new TH1D(Form("hdaughter_michel_score_%d_%d",i,j), Form("daughter_michel_score, %s, %s;Michel score", cutName[i], parTypeName[j]), 100, 0, 1);
      hdaughter_michel_score[i][j]->Sumw2();

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

  if (evt.MC){
    true_sliceID = int(evt.true_beam_endZ/thinslicewidth);
    if (true_sliceID < 0) true_sliceID = -1;
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

        std::vector<std::vector<double>> vAngCorr(nthinslices);
        if (evt.reco_beam_calo_wire->size()>=2){
          for (size_t i = 0; i<evt.reco_beam_calo_wire->size(); ++i){
            int this_sliceID = int((*evt.reco_beam_calo_Z)[i]/thinslicewidth);
            if (this_sliceID>=nthinslices) continue;
            if (this_sliceID<0) continue;
            TVector3 pt0, pt1;
            if (i !=  evt.reco_beam_calo_wire->size()-2){
              pt0 = TVector3((*evt.reco_beam_calo_X)[i],
                             (*evt.reco_beam_calo_Y)[i],
                             (*evt.reco_beam_calo_Z)[i]);
              pt1 = TVector3((*evt.reco_beam_calo_X)[i+1],
                             (*evt.reco_beam_calo_Y)[i+1],
                             (*evt.reco_beam_calo_Z)[i+1]);
            }
            else{
              pt0 = TVector3((*evt.reco_beam_calo_X)[i],
                             (*evt.reco_beam_calo_Y)[i],
                             (*evt.reco_beam_calo_Z)[i]);
              pt1 = TVector3((*evt.reco_beam_calo_X)[i-1],
                             (*evt.reco_beam_calo_Y)[i-1],
                             (*evt.reco_beam_calo_Z)[i-1]);
            }
            TVector3 dir = pt1 - pt0;
            dir = dir.Unit();
            vAngCorr[this_sliceID].push_back(dir.Z());
          }
          for (size_t i = 0; i<vAngCorr.size(); ++i){
            if (!vAngCorr[i].empty()){
              double sum_AngCorr = 0;
              for (size_t j = 0; j<vAngCorr[i].size(); ++j){
                sum_AngCorr += vAngCorr[i][j];
              }
              reco_AngCorr[i]->Fill(sum_AngCorr/vAngCorr[i].size());
            }
          }
        }
      }

      // True info
      if (!(evt.true_beam_traj_Z_SCE->empty())){
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

        std::vector<std::vector<double>> vAngCorr(nthinslices);
        if (evt.true_beam_traj_Z->size()>=2){
          for (size_t i = 0; i<evt.true_beam_traj_Z->size(); ++i){
            int this_sliceID = int((*evt.true_beam_traj_Z)[i]/thinslicewidth);
            if (this_sliceID>=nthinslices) continue;
            if (this_sliceID<0) continue;
            TVector3 pt0, pt1;
            if (i !=  evt.true_beam_traj_Z->size()-2){
              pt0 = TVector3((*evt.true_beam_traj_X)[i],
                             (*evt.true_beam_traj_Y)[i],
                             (*evt.true_beam_traj_Z)[i]);
              pt1 = TVector3((*evt.true_beam_traj_X)[i+1],
                             (*evt.true_beam_traj_Y)[i+1],
                             (*evt.true_beam_traj_Z)[i+1]);
            }
            else{
              pt0 = TVector3((*evt.true_beam_traj_X)[i],
                             (*evt.true_beam_traj_Y)[i],
                             (*evt.true_beam_traj_Z)[i]);
              pt1 = TVector3((*evt.true_beam_traj_X)[i-1],
                             (*evt.true_beam_traj_Y)[i-1],
                             (*evt.true_beam_traj_Z)[i-1]);
            }
            TVector3 dir = pt1 - pt0;
            dir = dir.Unit();
            vAngCorr[this_sliceID].push_back(dir.Z());
          }
          for (size_t i = 0; i<vAngCorr.size(); ++i){
            if (!vAngCorr[i].empty()){
              double sum_AngCorr = 0;
              for (size_t j = 0; j<vAngCorr[i].size(); ++j){
                sum_AngCorr += vAngCorr[i][j];
              }
              true_AngCorr[i]->Fill(sum_AngCorr/vAngCorr[i].size());
            }
          }
        }
      }
    }
  }

  if (!evt.reco_beam_calo_wire->empty()){
    reco_sliceID = int(evt.reco_beam_calo_endZ/thinslicewidth);
    if (reco_sliceID < 0) reco_sliceID = -1;
    if (reco_sliceID >= nthinslices) reco_sliceID = nthinslices;
  }

  bool isTestSample = true;
  if (evt.MC && evt.event%2 == 0) isTestSample = false;

  if (evt.true_beam_PDG == 211){
    if (isTestSample){
      h_truesliceid_pion_all->Fill(true_sliceID);
    }
    else{
      uf.eff_den_Inc->Fill(true_sliceID);
    }
    //if (evt.PassAllCuts() && evt.reco_beam_true_byE_matched){
    if (evt.PassAllCuts()){
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
      //if (evt.PassAllCuts() && evt.reco_beam_true_byE_matched){
      if (evt.PassAllCuts()){
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

void ThinSlice::FillHistograms(int cut, const HadAna & evt){
  
  if (cut>=0 && cut < nCuts){
    FillHistVec1D(htrue_beam_endZ[cut], evt.true_beam_endZ_SCE, evt.partype);
    FillHistVec1D(htrue_beam_endZ_SCE[cut], evt.true_beam_endZ, evt.partype);
    FillHistVec1D(htrue_sliceID[cut], true_sliceID, evt.partype);
    if (!evt.reco_beam_calo_wire->empty()){
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
    }      
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
  double avg_trueAngCorr[nthinslices] = {0};
  double avg_recoAngCorr[nthinslices] = {0};
  double err_trueAngCorr[nthinslices] = {0};
  double err_recoAngCorr[nthinslices] = {0};
  double reco_trueAngCorr[nthinslices] = {0};
  double err_reco_trueAngCorr[nthinslices] = {0};
  double truexs[nthinslices] = {0};
  double err_truexs[nthinslices] = {0};

  double NA=6.02214076e23;
  double MAr=39.95; //gmol
  double Density = 1.39; // g/cm^3

  for (int i = 0; i<nthinslices; ++i){
    
    slcid[i] = i;
    avg_trueincE[i] = true_incE[i]->GetMean();
    err_trueincE[i] = true_incE[i]->GetMeanError();
    avg_recoincE[i] = reco_incE[i]->GetMean();
    err_recoincE[i] = reco_incE[i]->GetMeanError();
    reco_trueincE[i] = avg_recoincE[i] - avg_trueincE[i];
    err_reco_trueincE[i] = sqrt(pow(err_trueincE[i],2)+pow(err_recoincE[i],2));
    avg_trueAngCorr[i] = true_AngCorr[i]->GetMean();
    err_trueAngCorr[i] = true_AngCorr[i]->GetMeanError();
    avg_recoAngCorr[i] = reco_AngCorr[i]->GetMean();
    err_recoAngCorr[i] = reco_AngCorr[i]->GetMeanError();
    reco_trueAngCorr[i] = avg_recoAngCorr[i] - avg_trueAngCorr[i];
    err_reco_trueAngCorr[i] = sqrt(pow(err_trueAngCorr[i],2)+pow(err_recoAngCorr[i],2));
    //std::cout<<i<<" "<<avg_trueincE[i]<<std::endl;
    if (true_incidents[i] && true_interactions[i]){
      truexs[i] = MAr/(Density*NA*thinslicewidth/avg_trueAngCorr[i])*log(true_incidents[i]/(true_incidents[i]-true_interactions[i]))*1e27;
      err_truexs[i] = MAr/(Density*NA*thinslicewidth/avg_trueAngCorr[i])*1e27*sqrt(true_interactions[i]+pow(true_interactions[i],2)/true_incidents[i])/true_incidents[i];
    }
  }

  TGraphErrors *gr_trueincE = new TGraphErrors(nthinslices, &(slcid[0]), &(avg_trueincE[0]), 0, &(err_trueincE[0]));
  TGraphErrors *gr_recoincE = new TGraphErrors(nthinslices, &(slcid[0]), &(avg_recoincE[0]), 0, &(err_recoincE[0]));
  TGraphErrors *gr_reco_trueincE = new TGraphErrors(nthinslices, &(slcid[0]), &(reco_trueincE[0]), 0, &(err_reco_trueincE[0]));

  TGraphErrors *gr_trueAngCorr = new TGraphErrors(nthinslices, &(slcid[0]), &(avg_trueAngCorr[0]), 0, &(err_trueAngCorr[0]));
  TGraphErrors *gr_recoAngCorr = new TGraphErrors(nthinslices, &(slcid[0]), &(avg_recoAngCorr[0]), 0, &(err_recoAngCorr[0]));
  TGraphErrors *gr_reco_trueAngCorr = new TGraphErrors(nthinslices, &(slcid[0]), &(reco_trueAngCorr[0]), 0, &(err_reco_trueAngCorr[0]));

  gr_trueincE->Write("gr_trueincE");
  gr_recoincE->Write("gr_recoincE");
  gr_reco_trueincE->Write("gr_reco_trueincE");

  gr_trueAngCorr->Write("gr_trueAngCorr");
  gr_recoAngCorr->Write("gr_recoAngCorr");
  gr_reco_trueAngCorr->Write("gr_reco_trueAngCorr");

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
    if (!evt.isTrueSelectedPart()) continue;
    evt.ProcessEvent();
    ProcessEvent(evt, uf);
    FillHistograms(kNocut, evt);
    if (evt.PassPandoraSliceCut()){
      FillHistograms(kPandoraSlice, evt);
      if (evt.PassBeamQualityCut()){
        FillHistograms(kBeamQuality, evt);
        if (evt.PassAPA3Cut()){
          FillHistograms(kAPA3, evt);
          if (evt.PassCaloSizeCut()){
            FillHistograms(kCaloSize, evt);
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
