#include "HadAna.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <iostream>

bool HadAna::isTrueSelectedPart(){
  for (size_t i = 0; i<truepdglist.size(); ++i){
    if (true_beam_PDG == truepdglist[i]) return true;
  }
  return false;
}

int HadAna::GetParType(){

  if (!reco_beam_true_byE_matched){
    return kMisID;
  }
  else if (true_beam_PDG == -13){
    return kPrimMuP;
  }
  else if (true_beam_PDG == 13){
    return kPrimMuM;
  }
  else if (true_beam_PDG == 211){
    if ((*true_beam_endProcess) == "pi+Inelastic"){
      return kPrimPiPInEl;
    }
    else return kPrimPiPEl;
  }
  else if (true_beam_PDG == 2212){
    if ((*true_beam_endProcess) == "protonInelastic"){
      return kPrimProInEl;
    }
    else return kPrimProEl;
  }
  
  return kOther;
}


bool HadAna::PassPandoraSliceCut(){
  return reco_beam_type == pandora_slice_pdg;
}

bool HadAna::PassBeamQualityCut(){

  double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
  double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;

  double projectX = (true_beam_startX + -1*true_beam_startZ*(true_beam_startDirX/true_beam_startDirZ) );
  double projectY = (true_beam_startY + -1*true_beam_startZ*(true_beam_startDirY/true_beam_startDirZ) );
  double cos = true_beam_startDirX*reco_beam_trackDirX + true_beam_startDirY*reco_beam_trackDirY + true_beam_startDirZ*reco_beam_trackDirZ;
  
  if ( (reco_beam_startX - projectX) < xlow )
    return false;
  
  if ( (reco_beam_startX - projectX) > xhigh )
    return false;
  
  if ( (reco_beam_startY - projectY) < ylow )
    return false;
  
  if ( (reco_beam_startY - projectY) > yhigh )
    return false;
  
  if (reco_beam_startZ < zlow || zhigh < reco_beam_startZ)
    return false;
  
  if ( cos < coslow)
    return false;
  
  return true;
  
};

bool HadAna::PassAPA3Cut(){

  double cutAPA3_Z = 226.;

  return reco_beam_endZ < cutAPA3_Z;
}

bool HadAna::PassCaloSizeCut(){
  
  return !(reco_beam_calo_wire->empty());
}

bool HadAna::PassMichelScoreCut(){
  
  return daughter_michel_score < 0.55;
}

bool HadAna::PassMediandEdxCut(){

  return median_dEdx < 2.4;
}

void HadAna::BookHistograms(){

  outputFile = TFile::Open(fOutputFileName.c_str(), "recreate");

  for (int i = 0; i < nCuts; ++i){
    for (int j = 0; j < nParTypes+1; ++j){
      htrue_beam_endZ[i][j] = new TH1D(Form("htrue_beam_endZ_%d_%d",i,j),Form("true_beam_endZ, %s, %s;true_beam_endZ (cm)", cutName[i], parTypeName[j]), 100, -100, 900);
      htrue_beam_endZ[i][j]->Sumw2();
      hreco_beam_endZ[i][j] = new TH1D(Form("hreco_beam_endZ_%d_%d",i,j),Form("reco_beam_endZ, %s, %s;reco_beam_endZ (cm)", cutName[i], parTypeName[j]), 100, -100, 900);
      hreco_beam_endZ[i][j]->Sumw2();
      hreco_true_beam_endZ[i][j] = new TH1D(Form("hreco_true_beam_endZ_%d_%d",i,j), Form("reco_true_beam_endZ, %s, %s;reco_beam_endZ - true_beam_endZ (cm)", cutName[i], parTypeName[j]), 100, -100, 100);
      hreco_true_beam_endZ[i][j]->Sumw2();
      hreco_vs_true_beam_endZ[i][j]= new TH2D(Form("hreco_vs_true_beam_endZ_%d_%d",i,j), Form("reco_vs_true_beam_endZ, %s, %s;true_beam_endZ (cm);reco_beam_endZ (cm)", cutName[i], parTypeName[j]), 100, -100, 900, 100, -100, 900);
      hreco_true_vs_true_beam_endZ[i][j]= new TH2D(Form("hreco_true_vs_true_beam_endZ_%d_%d",i,j), Form("reco_true_vs_true_beam_endZ, %s, %s;true_beam_endZ (cm);reco - true_beam_endZ (cm)", cutName[i], parTypeName[j]), 100, -100, 900, 100, -100, 100);

      htrue_sliceID[i][j] = new TH1D(Form("htrue_sliceID_%d_%d",i,j),Form("true_sliceID, %s, %s;true_sliceID (cm)", cutName[i], parTypeName[j]), 50, -1, 49);
      htrue_sliceID[i][j]->Sumw2();
      hreco_sliceID[i][j] = new TH1D(Form("hreco_sliceID_%d_%d",i,j),Form("reco_sliceID, %s, %s;reco_sliceID", cutName[i], parTypeName[j]), 50, -1, 49);
      hreco_sliceID[i][j]->Sumw2();
      hreco_true_sliceID[i][j] = new TH1D(Form("hreco_true_sliceID_%d_%d",i,j), Form("reco_true_sliceID, %s, %s;reco_sliceID - true_sliceID", cutName[i], parTypeName[j]), 20, -10, 10);
      hreco_true_sliceID[i][j]->Sumw2();
      hreco_vs_true_sliceID[i][j]= new TH2D(Form("hreco_vs_true_sliceID_%d_%d",i,j), Form("reco_vs_true_sliceID, %s, %s;true_sliceID;reco_sliceID", cutName[i], parTypeName[j]), 50, -1, 49, 50, -1, 49);
      hreco_true_vs_true_sliceID[i][j]= new TH2D(Form("hreco_true_vs_true_sliceID_%d_%d",i,j), Form("reco_true_vs_true_sliceID, %s, %s;true_sliceID;reco_sliceID - true_sliceID", cutName[i], parTypeName[j]), 50, -1, 49, 20, -10, 10);

      hreco_interacting_Energy_vs_true_beam_endZ[i][j]= new TH2D(Form("hreco_interacting_Energy_vs_true_beam_endZ_%d_%d",i,j), Form("hreco_interacting_Energy_vs_true_beam_endZ, %s, %s;true_beam_endZ (cm);reco interacting energy (MeV)", cutName[i], parTypeName[j]), 100, -100, 900, 100, 0, 1000);
      htrue_interacting_Energy_vs_true_beam_endZ[i][j]= new TH2D(Form("htrue_interacting_Energy_vs_true_beam_endZ_%d_%d",i,j), Form("htrue_interacting_Energy_vs_true_beam_endZ, %s, %s;true_beam_endZ (cm);true interacting energy (MeV)", cutName[i], parTypeName[j]), 100, -100, 900, 100, 0, 1000);
      hreco_true_interacting_Energy_vs_true_beam_endZ[i][j] = new TH2D(Form("hreco_true_interacting_Energy_vs_true_beam_endZ_%d_%d",i,j), Form("hreco_true_interacting_Energy_vs_true_beam_endZ, %s, %s;true_beam_endZ (cm); reco - true interacting energy (MeV)", cutName[i], parTypeName[j]), 100, -100, 900, 100, -100, 100);

      hmediandEdx[i][j] = new TH1D(Form("hmediandEdx_%d_%d",i,j), Form("mediandEdx, %s, %s;Median dE/dx (MeV/cm)", cutName[i], parTypeName[j]), 100, 0, 5);
      hmediandEdx[i][j]->Sumw2();

      hdaughter_michel_score[i][j] = new TH1D(Form("hdaughter_michel_score_%d_%d",i,j), Form("daughter_michel_score, %s, %s;Michel score", cutName[i], parTypeName[j]), 100, 0, 1);
      hdaughter_michel_score[i][j]->Sumw2();

      hrecoZsce[i][j] = new TH1D(Form("hrecoZsce_%d_%d",i,j), Form("recoZsce, %s, %s;Reco z (cm)", cutName[i], parTypeName[j]), 100, -100, 900);
      hrecoZsce[i][j]->Sumw2();

    }
  }

   for (int i = 0; i<nslices; ++i){
     reco_incE[i] = new TH1D(Form("reco_incE_%d",i),Form("Reco incident energy, %d<=wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinse, 0, 1200.);
     true_incE[i] = new TH1D(Form("true_incE_%d",i),Form("True incident energy, %d<=wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinse, 0, 1200.);
     reco_pitch[i] = new TH1D(Form("reco_pitch_%d",i),Form("Slice thickness, %d<=wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinsthickness, 0, 20.);
     reco_dEdx[i] = new TH1D(Form("reco_dEdx_%d",i),Form("dE/dx, %d<=wire #<%d;dE/dx (MeV/cm)",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinsthickness, 0, 20.);
     true_X[i] = new TH1D(Form("true_X_%d",i),Form("true_X, %d<=wire #<%d;True x (cm)",i*nwires_in_slice, (i+1)*nwires_in_slice), 100, -360, 0);
     true_Y[i] = new TH1D(Form("true_Y_%d",i),Form("true_Y, %d<=wire #<%d;True y (cm)",i*nwires_in_slice, (i+1)*nwires_in_slice), 100, 0, 600);
     true_Z[i] = new TH1D(Form("true_Z_%d",i),Form("true_Z, %d<=wire #<%d;True z (cm)",i*nwires_in_slice, (i+1)*nwires_in_slice), 100, -100, 500);
   }

   for (int i = 0; i<nslices; ++i){
     true_interactions[i] = 0;
     true_incidents[i] = 0;
   }

}

void HadAna::FillHistVec1D(TH1D *hist[nParTypes+1], const double &value){
  hist[0]->Fill(value);
  if (partype>=1 && partype < nParTypes+1){
    hist[partype]->Fill(value);
  }
}

void HadAna::FillHistVec2D(TH2D *hist[nParTypes+1], const double &value1, const double &value2){
  hist[0]->Fill(value1, value2);
  if (partype>=1 && partype < nParTypes+1){
    hist[partype]->Fill(value1, value2);
  }
}

void HadAna::FillHistograms(int cut){
  
  if (cut>=0 && cut < nCuts){
    FillHistVec1D(htrue_beam_endZ[cut], true_beam_endZ_SCE);
    FillHistVec1D(htrue_sliceID[cut], true_sliceID);
    FillHistVec2D(htrue_interacting_Energy_vs_true_beam_endZ[cut], true_beam_interactingEnergy, true_beam_endZ_SCE);
    if (!reco_beam_calo_wire->empty()){
      double reco_beam_endZ_SCE = reco_beam_calo_wire->back()*0.479 + 0.5603500;
      FillHistVec1D(hreco_beam_endZ[cut], reco_beam_endZ_SCE);
      FillHistVec1D(hreco_true_beam_endZ[cut], reco_beam_endZ_SCE - true_beam_endZ_SCE);
      FillHistVec2D(hreco_vs_true_beam_endZ[cut], true_beam_endZ_SCE, reco_beam_endZ_SCE);
      FillHistVec2D(hreco_true_vs_true_beam_endZ[cut], true_beam_endZ_SCE, reco_beam_endZ_SCE - true_beam_endZ_SCE);
      FillHistVec1D(hreco_sliceID[cut], reco_sliceID);
      FillHistVec1D(hreco_true_sliceID[cut], reco_sliceID - true_sliceID);
      FillHistVec2D(hreco_vs_true_sliceID[cut], true_sliceID, reco_sliceID);
      FillHistVec2D(hreco_true_vs_true_sliceID[cut], true_sliceID, reco_sliceID - true_sliceID);
      FillHistVec2D(hreco_interacting_Energy_vs_true_beam_endZ[cut], true_beam_endZ_SCE, reco_beam_interactingEnergy);
      FillHistVec2D(hreco_true_interacting_Energy_vs_true_beam_endZ[cut], true_beam_endZ_SCE, reco_beam_interactingEnergy - true_beam_interactingEnergy);
      FillHistVec1D(hmediandEdx[cut], median_dEdx);
      FillHistVec1D(hdaughter_michel_score[cut], daughter_michel_score);
      FillHistVec1D(hrecoZsce[cut], reco_beam_calo_endZ);
    }      
  }

  // Fill reco incE and pitch for each slice
  if (cut == kNocut && partype == kPrimPiPInEl){
    // Reco info
    if (!(reco_beam_calo_wire->empty())){
      std::vector<std::vector<double>> vpitch(nslices);
      std::vector<std::vector<double>> vincE(nslices);
      std::vector<std::vector<double>> vdEdx(nslices);
      for (size_t i = 0; i<reco_beam_calo_wire->size(); ++i){
        int this_wire = (*reco_beam_calo_wire)[i];
        int this_sliceID = this_wire/nwires_in_slice;
        //ignore the last slice for pitch and incident energy calculations
        if (this_sliceID>=reco_sliceID) continue;
        if (this_sliceID>=nslices) continue;
        if (this_sliceID<0) continue;

        double this_incE = (*reco_beam_incidentEnergies)[i];
        double this_pitch = (*reco_beam_TrkPitch_SCE)[i];
        double this_dEdx = (*reco_beam_dEdX_SCE)[i];
        vpitch[this_sliceID].push_back(this_pitch);
        vincE[this_sliceID].push_back(this_incE);
        vdEdx[this_sliceID].push_back(this_dEdx);
      }
      for (size_t i = 0; i<vpitch.size(); ++i){
        if (!vpitch[i].empty()){
          double sum_pitch = 0;
          for (size_t j = 0; j<vpitch[i].size(); ++j){
            //std::cout<<vpitch[i][j]<<std::endl;
            sum_pitch += vpitch[i][j];
          }
          //std::cout<<sum_pitch<<" "<<vpitch[i].size()<<std::endl;
          reco_pitch[i]->Fill(sum_pitch/vpitch[i].size()*nwires_in_slice);
        }
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
      for (size_t i = 0; i<vdEdx.size(); ++i){
        if (!vdEdx[i].empty()){
          double sum_dEdx = 0;
          for (size_t j = 0; j<vdEdx[i].size(); ++j){
            sum_dEdx += vdEdx[i][j];
          }
          reco_dEdx[i]->Fill(sum_dEdx/vdEdx[i].size());
        }
      }
    }
    // True info
    if (!(true_beam_traj_Z_SCE->empty())){
      std::vector<std::vector<double>> vincE(nslices);
      std::vector<std::vector<double>> vX(nslices);
      std::vector<std::vector<double>> vY(nslices);
      std::vector<std::vector<double>> vZ(nslices);
      for (size_t i = 0; i<true_beam_traj_Z_SCE->size()-1; ++i){//last point always has KE = 0
        int this_sliceID = int(((*true_beam_traj_Z_SCE)[i]-0.5603500-0.479/2)/0.479/nwires_in_slice);
        double this_incE = (*true_beam_traj_KE)[i];
//        if (run == 17775907 && event == 924){
//          std::cout<<i<<" "<<this_sliceID<<" "<<(*true_beam_traj_Z_SCE)[i]<<" "<<(*true_beam_traj_KE)[i]<<std::endl;
//        }
        if (this_sliceID>=nslices) continue;
        if (this_sliceID<0) continue;

        vincE[this_sliceID].push_back(this_incE);
        vX[this_sliceID].push_back((*true_beam_traj_X)[i]);
        vY[this_sliceID].push_back((*true_beam_traj_Y)[i]);
        vZ[this_sliceID].push_back((*true_beam_traj_Z)[i]);
      }
      for (size_t i = 0; i<vincE.size(); ++i){
        if (!vincE[i].empty()){
          double sum_incE = 0;
          for (size_t j = 0; j<vincE[i].size(); ++j){
            sum_incE += vincE[i][j];
          }
          true_incE[i]->Fill(sum_incE/vincE[i].size());
        }
        if (!vX[i].empty()){
          double sum_X = 0;
          for (size_t j = 0; j<vX[i].size(); ++j){
            sum_X += vX[i][j];
          }
          true_X[i]->Fill(sum_X/vX[i].size());
        }
        if (!vY[i].empty()){
          double sum_Y = 0;
          for (size_t j = 0; j<vY[i].size(); ++j){
            sum_Y += vY[i][j];
          }
          true_Y[i]->Fill(sum_Y/vY[i].size());
        }
        if (!vZ[i].empty()){
          double sum_Z = 0;
          for (size_t j = 0; j<vZ[i].size(); ++j){
            sum_Z += vZ[i][j];
          }
          true_Z[i]->Fill(sum_Z/vZ[i].size());
        }
      }
    }
  }
  
}

void HadAna::ProcessEvent(){

  partype = -1;
  reco_sliceID = -100;
  true_sliceID = -100;
  median_dEdx = -1;
  daughter_michel_score = -1;

  if (!reco_beam_calo_wire->empty()){
    reco_sliceID = reco_beam_calo_wire->back()/nwires_in_slice;
    median_dEdx = TMath::Median(reco_beam_calibrated_dEdX_SCE->size(), &(*reco_beam_calibrated_dEdX_SCE)[0]);
    if (!reco_daughter_PFP_michelScore_collection->empty()){
      daughter_michel_score = TMath::MaxElement(reco_daughter_PFP_michelScore_collection->size(), &(*reco_daughter_PFP_michelScore_collection)[0]);
      //std::cout<<daughter_michel_score<<std::endl;
    }
    else{
      daughter_michel_score = 0;
    }
  }

  if (MC){
    partype = GetParType();
    true_sliceID = int((true_beam_endZ_SCE-0.5603500-0.479/2)/0.479/nwires_in_slice);
    if (true_beam_PDG == 211){
      for (int i = 0; i<=true_sliceID; ++i){
        if (i<nslices) ++true_incidents[i];
      }
    }
    if ((*true_beam_endProcess) == "pi+Inelastic"){
      if (true_sliceID < nslices && true_sliceID>=0){
        ++true_interactions[true_sliceID];
      }
    }
  }
}

void HadAna::SaveHistograms(){

  double slcid[nslices] = {0};
  double avg_recoincE[nslices] = {0};
  double avg_recotrueincE[nslices] = {0};
  double avg_recopitch[nslices] = {0};
  double avg_recodEdx[nslices] = {0};
  double avg_trueincE[nslices] = {0};
  double avg_trueX[nslices] = {0};
  double avg_trueY[nslices] = {0};
  double avg_trueZ[nslices] = {0};
  double avg_truepitch[nslices] = {0};
  double avg_truedEdx[nslices] = {0};
  double avg_truedE[nslices] = {0};
  double truexs_thinslice[nslices] = {0};
  double truexs_eslice[nslices] = {0};
  double err_recoincE[nslices] = {0};
  double err_recotrueincE[nslices] = {0};
  double err_recopitch[nslices] = {0};
  double err_recodEdx[nslices] = {0};
  double err_trueincE[nslices] = {0};
  double err_trueX[nslices] = {0};
  double err_trueY[nslices] = {0};
  double err_trueZ[nslices] = {0};
  double err_truepitch[nslices] = {0};
  double err_truedEdx[nslices] = {0};
  double err_truexs_thinslice[nslices] = {0};
  double err_truexs_eslice[nslices] = {0};

  for (int i = 0; i<nslices; ++i){
    slcid[i] = i;
    avg_recoincE[i] = reco_incE[i]->GetMean();
    avg_trueincE[i] = true_incE[i]->GetMean();
    avg_trueX[i] = true_X[i]->GetMean();
    avg_trueY[i] = true_Y[i]->GetMean();
    avg_trueZ[i] = true_Z[i]->GetMean();
    avg_recotrueincE[i] = reco_incE[i]->GetMean() - true_incE[i]->GetMean();
    avg_recopitch[i] = reco_pitch[i]->GetMean();
    avg_recodEdx[i] = reco_dEdx[i]->GetMean();
    err_recoincE[i] = reco_incE[i]->GetMeanError();
    err_trueincE[i] = true_incE[i]->GetMeanError();
    err_trueX[i] = true_X[i]->GetMeanError();
    err_trueY[i] = true_Y[i]->GetMeanError();
    err_trueZ[i] = true_Z[i]->GetMeanError();
    err_recotrueincE[i] = sqrt(pow(reco_incE[i]->GetMeanError(),2) + pow(true_incE[i]->GetMeanError(),2));
    err_recopitch[i] = reco_pitch[i]->GetMeanError();
    err_recopitch[i] = reco_pitch[i]->GetMeanError();
    err_recodEdx[i] = reco_dEdx[i]->GetMeanError();
  }

  for (int i = 0; i<nslices; ++i){
    if (i!=0 && i!=nslices-1){
      avg_truepitch[i] = sqrt(pow(avg_trueX[i+1]-avg_trueX[i-1],2)+
                              pow(avg_trueY[i+1]-avg_trueY[i-1],2)+
                              pow(avg_trueZ[i+1]-avg_trueZ[i-1],2))/2;
      avg_truedE[i] = (avg_trueincE[i-1]-avg_trueincE[i+1])/2;
      avg_truedEdx[i] = (avg_trueincE[i-1]-avg_trueincE[i+1])/2/avg_truepitch[i];
    }
    else if (i==0){
      avg_truepitch[i] = sqrt(pow(avg_trueX[i+1]-avg_trueX[i],2)+
                              pow(avg_trueY[i+1]-avg_trueY[i],2)+
                              pow(avg_trueZ[i+1]-avg_trueZ[i],2));
      avg_truedE[i] = avg_trueincE[i]-avg_trueincE[i+1];
      avg_truedEdx[i] = (avg_trueincE[i]-avg_trueincE[i+1])/avg_truepitch[i];
    }
    else if (i==nslices-1){
      avg_truepitch[i] = sqrt(pow(avg_trueX[i]-avg_trueX[i-1],2)+
                              pow(avg_trueY[i]-avg_trueY[i-1],2)+
                              pow(avg_trueZ[i]-avg_trueZ[i-1],2));
      avg_truedE[i] = avg_trueincE[i-1]-avg_trueincE[i];
      avg_truedEdx[i] = (avg_trueincE[i-1]-avg_trueincE[i])/avg_truepitch[i];
    }
  }

  double NA=6.02214076e23;
  double MAr=39.95; //gmol
  double Density = 1.39; // g/cm^3
  for (int i = 0; i<nslices; ++i){
    if (true_incidents[i] && true_interactions[i] && avg_recopitch[i]){
      truexs_thinslice[i] = MAr/(Density*NA*avg_recopitch[i])*log(true_incidents[i]/(true_incidents[i]-true_interactions[i]))*1e27;
      err_truexs_thinslice[i] = MAr/(Density*NA*avg_recopitch[i])*1e27*sqrt(true_interactions[i]+pow(true_interactions[i],2)/true_incidents[i])/true_incidents[i];
    }
    if (true_incidents[i] && true_interactions[i]){
      truexs_eslice[i] = MAr/(Density*NA*avg_truedE[i])*2.3*log(true_incidents[i]/(true_incidents[i]-true_interactions[i]))*1e27;
      err_truexs_eslice[i] = MAr/(Density*NA*avg_truedE[i])*2.3*1e27*sqrt(true_interactions[i]+pow(true_interactions[i],2)/true_incidents[i])/true_incidents[i];
    }
    //std::cout<<i<<" "<<avg_trueincE[i]<<" "<<truexs_thinslice[i]<<" "<<err_truexs_thinslice[i]<<std::endl;
  }

  TGraphErrors *gr_recoincE_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_recoincE[0]), 0, &(err_recoincE[0]));
  TGraphErrors *gr_trueincE_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_trueincE[0]), 0, &(err_trueincE[0]));
  TGraphErrors *gr_trueX_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_trueX[0]), 0, &(err_trueX[0]));
  TGraphErrors *gr_trueY_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_trueY[0]), 0, &(err_trueY[0]));
  TGraphErrors *gr_trueZ_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_trueZ[0]), 0, &(err_trueZ[0]));
  TGraphErrors *gr_recotrueincE_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_recotrueincE[0]), 0, &(err_recotrueincE[0]));
  TGraphErrors *gr_recopitch_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_recopitch[0]), 0, &(err_recopitch[0]));
  TGraphErrors *gr_recodEdx_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_recodEdx[0]), 0, &(err_recodEdx[0]));
  TGraphErrors *gr_truepitch_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_truepitch[0]), 0, &(err_truepitch[0]));
  TGraphErrors *gr_truedEdx_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_truedEdx[0]), 0, &(err_truedEdx[0]));
  TGraphErrors *gr_truexs_thinslice = new TGraphErrors(nslices, &(avg_trueincE[0]), &(truexs_thinslice[0]), 0, &(err_truexs_thinslice[0]));
  TGraphErrors *gr_truexs_eslice = new TGraphErrors(nslices, &(avg_trueincE[0]), &(truexs_eslice[0]), 0, &(err_truexs_eslice[0]));

  outputFile->cd();
  gr_recoincE_slc->Write("gr_recoincE_slc");
  gr_trueincE_slc->Write("gr_trueincE_slc");
  gr_recotrueincE_slc->Write("gr_recotrueincE_slc");
  gr_recopitch_slc->Write("gr_recopitch_slc");
  gr_recodEdx_slc->Write("gr_recodEdx_slc");
  gr_truedEdx_slc->Write("gr_truedEdx_slc");
  gr_trueX_slc->Write("gr_trueX_slc");
  gr_trueY_slc->Write("gr_trueY_slc");
  gr_trueZ_slc->Write("gr_trueZ_slc");
  gr_truepitch_slc->Write("gr_truepitch_slc");
  gr_truexs_thinslice->Write("gr_truexs_thinslice");
  gr_truexs_eslice->Write("gr_truexs_eslice");

  outputFile->Write();
}

