#include "HadAna.h"
#include "TGraphErrors.h"
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
  
  return kUnknown;
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
    }
  }

   for (int i = 0; i<nslices; ++i){
     reco_incE[i] = new TH1D(Form("reco_incE_%d",i),Form("Reco incident energy, %d<=wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinse, 0, 1200.);
     true_incE[i] = new TH1D(Form("true_incE_%d",i),Form("True incident energy, %d<=wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinse, 0, 1200.);
     reco_pitch[i] = new TH1D(Form("reco_pitch_%d",i),Form("Slice thickness, %d<=wire #<%d",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinsthickness, 0, 20.);
     reco_dEdx[i] = new TH1D(Form("reco_dEdx_%d",i),Form("dE/dx, %d<=wire #<%d;dE/dx (MeV/cm)",i*nwires_in_slice, (i+1)*nwires_in_slice), nbinsthickness, 0, 20.);
   }
}

void HadAna::FillHistograms(int cut){
  
  if (cut>=0 && cut < nCuts){
    htrue_beam_endZ[cut][0]->Fill(true_beam_endZ_SCE);
    htrue_sliceID[cut][0]->Fill(true_sliceID);
    htrue_interacting_Energy_vs_true_beam_endZ[cut][0]->Fill(true_beam_endZ_SCE, true_beam_interactingEnergy);
    if (!reco_beam_calo_wire->empty()){
      //std::cout<<true_beam_incidentEnergies->back()<<" "<<reco_beam_incidentEnergies->back()<<std::endl;
      double reco_beam_endZ_SCE = reco_beam_calo_wire->back()*0.479 + 0.5603500;
      hreco_beam_endZ[cut][0]->Fill(reco_beam_endZ_SCE);
      hreco_true_beam_endZ[cut][0]->Fill(reco_beam_endZ_SCE - true_beam_endZ_SCE);
      hreco_vs_true_beam_endZ[cut][0]->Fill(true_beam_endZ_SCE, reco_beam_endZ_SCE);
      hreco_true_vs_true_beam_endZ[cut][0]->Fill(true_beam_endZ_SCE, reco_beam_endZ_SCE - true_beam_endZ_SCE);

      hreco_sliceID[cut][0]->Fill(reco_sliceID);
      hreco_true_sliceID[cut][0]->Fill(reco_sliceID - true_sliceID);
      hreco_vs_true_sliceID[cut][0]->Fill(true_sliceID, reco_sliceID);
      hreco_true_vs_true_sliceID[cut][0]->Fill(true_sliceID, reco_sliceID - true_sliceID);

      hreco_interacting_Energy_vs_true_beam_endZ[cut][0]->Fill(true_beam_endZ_SCE, reco_beam_interactingEnergy);
      hreco_true_interacting_Energy_vs_true_beam_endZ[cut][0]->Fill(true_beam_endZ_SCE, reco_beam_interactingEnergy - true_beam_interactingEnergy);
      
    }
    if (partype>=1 && partype < nParTypes+1){
      htrue_beam_endZ[cut][partype]->Fill(true_beam_endZ_SCE);
      htrue_sliceID[cut][partype]->Fill(true_sliceID);
      htrue_interacting_Energy_vs_true_beam_endZ[cut][partype]->Fill(true_beam_endZ_SCE, true_beam_interactingEnergy);
      if (!reco_beam_calo_wire->empty()){
        double reco_beam_endZ_SCE = reco_beam_calo_wire->back()*0.479 + 0.5603500;
        hreco_beam_endZ[cut][partype]->Fill(reco_beam_endZ_SCE);
        hreco_true_beam_endZ[cut][partype]->Fill(reco_beam_endZ_SCE - true_beam_endZ_SCE);
        hreco_vs_true_beam_endZ[cut][partype]->Fill(true_beam_endZ_SCE, reco_beam_endZ_SCE);
        hreco_true_vs_true_beam_endZ[cut][partype]->Fill(true_beam_endZ_SCE, reco_beam_endZ_SCE - true_beam_endZ_SCE);

        hreco_sliceID[cut][partype]->Fill(reco_sliceID);
        hreco_true_sliceID[cut][partype]->Fill(reco_sliceID - true_sliceID);
        hreco_vs_true_sliceID[cut][partype]->Fill(true_sliceID, reco_sliceID);
        hreco_true_vs_true_sliceID[cut][partype]->Fill(true_sliceID, reco_sliceID - true_sliceID);

        hreco_interacting_Energy_vs_true_beam_endZ[cut][partype]->Fill(true_beam_endZ_SCE, reco_beam_interactingEnergy);
        hreco_true_interacting_Energy_vs_true_beam_endZ[cut][partype]->Fill(true_beam_endZ_SCE, reco_beam_interactingEnergy - true_beam_interactingEnergy);

      }
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
      for (size_t i = 0; i<true_beam_traj_Z_SCE->size()-1; ++i){//last point always has KE = 0
        int this_sliceID = int(((*true_beam_traj_Z_SCE)[i]-0.5603500-0.479/2)/0.479/nwires_in_slice);
        double this_incE = (*true_beam_traj_KE)[i];
//        if (run == 17775907 && event == 924){
//          std::cout<<i<<" "<<this_sliceID<<" "<<(*true_beam_traj_Z_SCE)[i]<<" "<<(*true_beam_traj_KE)[i]<<std::endl;
//        }
        if (this_sliceID>=nslices) continue;
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
//          if (i==9 && sum_incE/vincE[i].size() < 400){
//            //std::cout<<sum_incE<<" "<<vincE[i].size()<<std::endl;
//            //std::cout<<run<<" "<<event<<std::endl;
//          }
        }
      }
    }
  }
  
}

void HadAna::ProcessEvent(){

  partype = -1;
  reco_sliceID = -100;
  true_sliceID = -100;

  if (!reco_beam_calo_wire->empty()){
    reco_sliceID = reco_beam_calo_wire->back()/nwires_in_slice;
  }

  if (MC){
    partype = GetParType();
    true_sliceID = int((true_beam_endZ_SCE-0.5603500-0.479/2)/0.479/nwires_in_slice);
  }
}

void HadAna::SaveHistograms(){

  double slcid[nslices] = {0};
  double avg_recoincE[nslices] = {0};
  double avg_trueincE[nslices] = {0};
  double avg_recotrueincE[nslices] = {0};
  double avg_recopitch[nslices] = {0};
  double avg_recodEdx[nslices] = {0};
  double err_recoincE[nslices] = {0};
  double err_trueincE[nslices] = {0};
  double err_recotrueincE[nslices] = {0};
  double err_recopitch[nslices] = {0};
  double err_recodEdx[nslices] = {0};

  for (int i = 0; i<nslices; ++i){
    slcid[i] = i;
    avg_recoincE[i] = reco_incE[i]->GetMean();
    avg_trueincE[i] = true_incE[i]->GetMean();
    avg_recotrueincE[i] = reco_incE[i]->GetMean() - true_incE[i]->GetMean();
    avg_recopitch[i] = reco_pitch[i]->GetMean();
    avg_recodEdx[i] = reco_dEdx[i]->GetMean();
    err_recoincE[i] = reco_incE[i]->GetMeanError();
    err_trueincE[i] = true_incE[i]->GetMeanError();
    err_recotrueincE[i] = sqrt(pow(reco_incE[i]->GetMeanError(),2) + pow(true_incE[i]->GetMeanError(),2));
    err_recopitch[i] = reco_pitch[i]->GetMeanError();
    err_recopitch[i] = reco_pitch[i]->GetMeanError();
    err_recodEdx[i] = reco_dEdx[i]->GetMeanError();
  }

  TGraphErrors *gr_recoincE_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_recoincE[0]), 0, &(err_recoincE[0]));
  TGraphErrors *gr_trueincE_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_trueincE[0]), 0, &(err_trueincE[0]));
  TGraphErrors *gr_recotrueincE_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_recotrueincE[0]), 0, &(err_recotrueincE[0]));
  TGraphErrors *gr_recopitch_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_recopitch[0]), 0, &(err_recopitch[0]));
  TGraphErrors *gr_recodEdx_slc = new TGraphErrors(nslices, &(slcid[0]), &(avg_recodEdx[0]), 0, &(err_recodEdx[0]));

  outputFile->cd();
  gr_recoincE_slc->Write("gr_recoincE_slc");
  gr_trueincE_slc->Write("gr_trueincE_slc");
  gr_recotrueincE_slc->Write("gr_recotrueincE_slc");
  gr_recopitch_slc->Write("gr_recopitch_slc");
  gr_recodEdx_slc->Write("gr_recodEdx_slc");

  outputFile->Write();
}
