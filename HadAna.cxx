#include "HadAna.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TVector3.h"
#include "util.h"
#include <iostream>

HadAna::HadAna(TTree *tree) 
  : anavar(tree){

  AddTruePDG(-13);
  AddTruePDG(13);
  AddTruePDG(211);

  SetPandoraSlicePDG(13);

  SetBeamQualityCuts();

}

void HadAna::AddTruePDG(int pdg){
  truepdglist.push_back(pdg);
};

bool HadAna::isSelectedPart() const{
  if (MC){
    for (size_t i = 0; i<truepdglist.size(); ++i){
      if (true_beam_PDG == truepdglist[i]) return true;
    }
    return false;
  }
  else{
    if (!beam_inst_valid) return false;
    if (beam_inst_nMomenta != 1 || beam_inst_nTracks != 1) return false;
    for (size_t i = 0; i<truepdglist.size(); ++i){
      for (size_t j = 0; j<beam_inst_PDG_candidates->size(); ++j){
        if ((*beam_inst_PDG_candidates)[j] == truepdglist[i]) return true;
      }
    }
    return false;
  }
}

void HadAna::SetPandoraSlicePDG(int pdg){
  pandora_slice_pdg = pdg;
};

void HadAna::SetBeamQualityCuts(double dx_min, double dx_max,
                                double dy_min, double dy_max,
                                double dz_min, double dz_max,
                                double costh_min, double costh_max){
  beamcut_dx_min = dx_min; beamcut_dx_max = dx_max;
  beamcut_dy_min = dy_min; beamcut_dy_max = dy_max;
  beamcut_dz_min = dz_min; beamcut_dz_max = dz_max;
  beamcut_costh_min = costh_min; beamcut_costh_max = costh_max;
}

int HadAna::GetParType(){

  if (!MC){
    return kData;
  }
  else if (event%2){
    return kData;
  }
  else if (!reco_beam_true_byE_matched){
    if (reco_beam_true_byE_origin == 2) {
      return kMIDcosmic;
    }
    else if (std::abs(reco_beam_true_byE_PDG) == 211){
      return kMIDpi;
    }
    else if (reco_beam_true_byE_PDG == 2212){
      return kMIDp;
    }
    else if (std::abs(reco_beam_true_byE_PDG) == 13){
      return kMIDmu;
    }
    else if (std::abs(reco_beam_true_byE_PDG) == 11 ||
             reco_beam_true_byE_PDG == 22){
      return kMIDeg;
    }
    else {
      //std::cout<<reco_beam_true_byE_PDG<<std::endl;
      return kMIDother;
    }
  }
  else if (true_beam_PDG == -13){
    return kMuon;
  }
  else if (true_beam_PDG == 211){
    if ((*true_beam_endProcess) == "pi+Inelastic"){
      return kPiInel;
    }
    else return kPiElas;
  }
  
  return kMIDother;
}


bool HadAna::PassPandoraSliceCut() const{

  return (reco_beam_type == pandora_slice_pdg);
}

bool HadAna::PassBeamQualityCut() const{

  if (beam_dx<beamcut_dx_min) return false;
  if (beam_dx>beamcut_dx_max) return false;
  if (beam_dy<beamcut_dy_min) return false;
  if (beam_dy>beamcut_dy_max) return false;
  //if (beam_dx*beam_dx+beam_dy*beam_dy>1.5*1.5) return false;
  if (beam_dz<beamcut_dz_min) return false;
  if (beam_dz>beamcut_dz_max) return false;
  if (beam_costh<beamcut_costh_min) return false;
  if (beam_costh>beamcut_costh_max) return false;

  return true;
}

bool HadAna::PassAPA3Cut() const{

  double cutAPA3_Z = 226.;

  return reco_beam_endZ < cutAPA3_Z;
}

bool HadAna::PassCaloSizeCut() const{
  
  return !(reco_beam_calo_wire->empty());
}

bool HadAna::PassMichelScoreCut() const{
  
  return daughter_michel_score < 0.55;
}

bool HadAna::PassMediandEdxCut() const{

  return median_dEdx < 2.4;
}

bool HadAna::PassAllCuts() const{
  return PassPandoraSliceCut()&&
    PassCaloSizeCut()&&
    PassBeamQualityCut()&&
    PassAPA3Cut()&&
    PassMichelScoreCut()&&
    PassMediandEdxCut();
}

void HadAna::ProcessEvent(){

  partype = GetParType();

  median_dEdx = -1;
  daughter_michel_score = -999;

  if (!reco_beam_calo_wire->empty()){
    median_dEdx = TMath::Median(reco_beam_calibrated_dEdX_SCE->size(), &(*reco_beam_calibrated_dEdX_SCE)[0]);
//    daughter_michel_score = 0;
//    int nhits = 0;
//    for (size_t i = 0; i<reco_daughter_PFP_michelScore_collection->size(); ++i){
//      nhits += (*reco_daughter_PFP_nHits_collection)[i];
//      daughter_michel_score += (*reco_daughter_PFP_michelScore_collection)[0] * (*reco_daughter_PFP_nHits_collection)[i];
//    }
//    if (nhits) daughter_michel_score/=nhits;
//    else daughter_michel_score = -999;
    if (reco_beam_vertex_nHits) daughter_michel_score = reco_beam_vertex_michel_score/reco_beam_vertex_nHits;
  }

  beam_dx = -999;
  beam_dy = -999;
  beam_dz = -999;
  beam_costh = -999;

  if (!reco_beam_calo_wire->empty()){

    TVector3 pt0(reco_beam_calo_startX,
                 reco_beam_calo_startY,
                 reco_beam_calo_startZ);
    TVector3 pt1(reco_beam_calo_endX,
                 reco_beam_calo_endY,
                 reco_beam_calo_endZ);
    TVector3 dir = pt1 - pt0;
    dir = dir.Unit();

    if (MC){
      beam_dx = (reco_beam_calo_startX - beam_startX_mc)/beam_startX_rms_mc;
      beam_dy = (reco_beam_calo_startY - beam_startY_mc)/beam_startY_rms_mc;
      beam_dz = (reco_beam_calo_startZ - beam_startZ_mc)/beam_startZ_rms_mc;
      TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180),
                       cos(beam_angleY_mc*TMath::Pi()/180),
                       cos(beam_angleZ_mc*TMath::Pi()/180));
      beamdir = beamdir.Unit();
      beam_costh = dir.Dot(beamdir);
    }
    else{
      beam_dx = (reco_beam_calo_startX - beam_startX_data)/beam_startX_rms_data;
      beam_dy = (reco_beam_calo_startY - beam_startY_data)/beam_startY_rms_data;
      beam_dz = (reco_beam_calo_startZ - beam_startZ_data)/beam_startZ_rms_data;
      TVector3 beamdir(cos(beam_angleX_data*TMath::Pi()/180),
                       cos(beam_angleY_data*TMath::Pi()/180),
                       cos(beam_angleZ_data*TMath::Pi()/180));
      beamdir = beamdir.Unit();
      beam_costh = dir.Dot(beamdir);
    }
  }

  dEdx_5cm = -1;
  /*
  if (!reco_beam_allTrack_calibrated_dEdX->empty()){
    dEdx_5cm = 0;
    int nhits = 0;
    for (int i = 0; i<reco_beam_allTrack_calibrated_dEdX->size(); ++i){
      std::cout<<i<<" "<<reco_beam_allTrack_resRange->back()-(*reco_beam_allTrack_resRange)[i]<<" "<<(*reco_beam_allTrack_calibrated_dEdX)[i]<<endl;
      if (std::abs(reco_beam_allTrack_resRange->back()-(*reco_beam_allTrack_resRange)[i])<5){
        dEdx_5cm += (*reco_beam_allTrack_calibrated_dEdX)[i];
        ++nhits;
      }
    }
    if (nhits) dEdx_5cm/=nhits;
    else dEdx_5cm = -1;
  }
  */

  //if (event == 78467) cout<<reco_beam_calibrated_dEdX_SCE->size()<<endl;
  //cout<<reco_beam_allTrack_calibrated_dEdX->size()<<endl;
  if (!reco_beam_calibrated_dEdX_SCE->empty()){
    //dEdx_5cm = 0;
    //int nhits = 0;
    std::vector<double> vdEdx;
    for (int i = 0; i<reco_beam_calibrated_dEdX_SCE->size(); ++i){
      //std::cout<<i<<" "<<reco_beam_resRange_SCE->back()-(*reco_beam_resRange_SCE)[i]<<" "<<(*reco_beam_calibrated_dEdX_SCE)[i]<<endl;
      //if (event == 78467) cout<<(*reco_beam_resRange_SCE)[i]<<" "<<(*reco_beam_calo_Z)[i]<<" "<<(*reco_beam_calibrated_dEdX_SCE)[i]<<endl;
      if ((*reco_beam_resRange_SCE)[i]<5){
        vdEdx.push_back((*reco_beam_calibrated_dEdX_SCE)[i]);
        //dEdx_5cm += (*reco_beam_calibrated_dEdX_SCE)[i];
        //++nhits;
      }
    }
    //if (nhits) dEdx_5cm/=nhits;
    //else dEdx_5cm = -1;
    if (!vdEdx.empty()){
      dEdx_5cm = TMath::Median(vdEdx.size(), &vdEdx[0]);
    }
  }
//  if (!MC && reco_beam_PFP_trackScore_collection>=0 && reco_beam_PFP_trackScore_collection<0.01){
//    cout<<run<<" "<<event<<endl;
//  }

}

