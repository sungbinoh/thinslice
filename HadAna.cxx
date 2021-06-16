#include "anavar.h"
#include "HadAna.h"
#include "TMath.h"
#include "TVector3.h"
#include <iostream>

HadAna::HadAna(){
}

void HadAna::InitPi(){

  AddTruePDG(-13);
  AddTruePDG(13);
  AddTruePDG(211);

  SetPandoraSlicePDG(13);

  SetBeamQualityCuts();
}

void HadAna::InitP(){

  AddTruePDG(2212);

  SetPandoraSlicePDG(13);

  SetBeamQualityCuts();
}


void HadAna::AddTruePDG(int pdg){
  truepdglist.push_back(pdg);
};

bool HadAna::isSelectedPart(const anavar& evt) const{
  if (evt.MC){
    for (size_t i = 0; i<truepdglist.size(); ++i){
      if (evt.true_beam_PDG == truepdglist[i]) return true;
    }
    return false;
  }
  else{
    if (!evt.beam_inst_valid) return false;
    if (evt.beam_inst_nMomenta != 1 || evt.beam_inst_nTracks != 1) return false;
    for (size_t i = 0; i<truepdglist.size(); ++i){
      for (size_t j = 0; j<evt.beam_inst_PDG_candidates->size(); ++j){
        if ((*evt.beam_inst_PDG_candidates)[j] == truepdglist[i]) return true;
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
                                double dxy_min, double dxy_max,
                                double costh_min, double costh_max){
  beamcut_dx_min = dx_min; beamcut_dx_max = dx_max;
  beamcut_dy_min = dy_min; beamcut_dy_max = dy_max;
  beamcut_dz_min = dz_min; beamcut_dz_max = dz_max;
  beamcut_dxy_min = dxy_min; beamcut_dxy_max = dxy_max;
  beamcut_costh_min = costh_min; beamcut_costh_max = costh_max;
}

int HadAna::GetPiParType(const anavar& evt){

  if (!evt.MC){
    return pi::kData;
  }
  else if (evt.event%2){
    return pi::kData;
  }
  else if (!evt.reco_beam_true_byE_matched){
    if (evt.reco_beam_true_byE_origin == 2) {
      return pi::kMIDcosmic;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 211){
      return pi::kMIDpi;
    }
    else if (evt.reco_beam_true_byE_PDG == 2212){
      return pi::kMIDp;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 13){
      return pi::kMIDmu;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 11 ||
             evt.reco_beam_true_byE_PDG == 22){
      return pi::kMIDeg;
    }
    else {
      //std::cout<<reco_beam_true_byE_PDG<<std::endl;
      return pi::kMIDother;
    }
  }
  else if (evt.true_beam_PDG == -13){
    return pi::kMuon;
  }
  else if (evt.true_beam_PDG == 211){
    if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
      return pi::kPiInel;
    }
    else return pi::kPiElas;
  }
  
  return pi::kMIDother;
}

int HadAna::GetPParType(const anavar& evt){

  if (!evt.MC){
    return p::kData;
  }
  else if (evt.event%2){
    return p::kData;
  }
  else if (!evt.reco_beam_true_byE_matched){
    if (evt.reco_beam_true_byE_origin == 2) {
      return p::kMIDcosmic;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 211){
      return p::kMIDpi;
    }
    else if (evt.reco_beam_true_byE_PDG == 2212){
      return p::kMIDp;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 13){
      return p::kMIDmu;
    }
    else if (std::abs(evt.reco_beam_true_byE_PDG) == 11 ||
             evt.reco_beam_true_byE_PDG == 22){
      return p::kMIDeg;
    }
    else {
      //std::cout<<reco_beam_true_byE_PDG<<std::endl;
      return p::kMIDother;
    }
  }
  else if (evt.true_beam_PDG == 2212){
    if ((*evt.true_beam_endProcess) == "protonInelastic"){
      return p::kPInel;
    }
    else return p::kPElas;
  }
  
  return p::kMIDother;
}

bool HadAna::PassPandoraSliceCut(const anavar& evt) const{

  return (evt.reco_beam_type == pandora_slice_pdg);
}

bool HadAna::PassBeamQualityCut() const{

  if (beamcut_dx_min<beamcut_dx_max){
    if (beam_dx<beamcut_dx_min) return false;
    if (beam_dx>beamcut_dx_max) return false;
  }

  if (beamcut_dy_min<beamcut_dy_max){
    if (beam_dy<beamcut_dy_min) return false;
    if (beam_dy>beamcut_dy_max) return false;
  }

  if (beamcut_dz_min<beamcut_dz_max){
    if (beam_dz<beamcut_dz_min) return false;
    if (beam_dz>beamcut_dz_max) return false;
  }

  if (beamcut_dxy_min<beamcut_dxy_max){
    if (beam_dxy<beamcut_dxy_min) return false;
    if (beam_dxy>beamcut_dxy_max) return false;
  }

  if (beamcut_costh_min<beamcut_costh_max){
    if (beam_costh<beamcut_costh_min) return false;
    if (beam_costh>beamcut_costh_max) return false;
  }

  return true;
}

bool HadAna::PassAPA3Cut(const anavar& evt) const{

  double cutAPA3_Z = 220.;

  return evt.reco_beam_calo_endZ < cutAPA3_Z;
}

bool HadAna::PassCaloSizeCut(const anavar& evt) const{
  
  return !(evt.reco_beam_calo_wire->empty());
}

bool HadAna::PassMichelScoreCut() const{
  
  return daughter_michel_score < 0.55;
}

bool HadAna::PassMediandEdxCut() const{

  return median_dEdx < 2.4;
}

bool HadAna::PassPiCuts(const anavar& evt) const{
  return PassPandoraSliceCut(evt)&&
    PassCaloSizeCut(evt)&&
    PassBeamQualityCut()&&
    PassAPA3Cut(evt)&&
    PassMichelScoreCut()&&
    PassMediandEdxCut();
}

bool HadAna::PassPCuts(const anavar& evt) const{
  return PassPandoraSliceCut(evt)&&
    PassCaloSizeCut(evt)&&
    PassBeamQualityCut();
}

void HadAna::ProcessEvent(const anavar& evt){

  pitype = GetPiParType(evt);
  ptype = GetPParType(evt);

  median_dEdx = -1;
  daughter_michel_score = -999;

  if (!evt.reco_beam_calo_wire->empty()){
    median_dEdx = TMath::Median(evt.reco_beam_calibrated_dEdX_SCE->size(), &(*evt.reco_beam_calibrated_dEdX_SCE)[0]);
//    daughter_michel_score = 0;
//    int nhits = 0;
//    for (size_t i = 0; i<reco_daughter_PFP_michelScore_collection->size(); ++i){
//      nhits += (*reco_daughter_PFP_nHits_collection)[i];
//      daughter_michel_score += (*reco_daughter_PFP_michelScore_collection)[0] * (*reco_daughter_PFP_nHits_collection)[i];
//    }
//    if (nhits) daughter_michel_score/=nhits;
//    else daughter_michel_score = -999;
    if (evt.reco_beam_vertex_nHits) daughter_michel_score = evt.reco_beam_vertex_michel_score/evt.reco_beam_vertex_nHits;
  }

  beam_dx = -999;
  beam_dy = -999;
  beam_dz = -999;
  beam_dxy = -999;
  beam_costh = -999;

  if (!evt.reco_beam_calo_wire->empty()){

    TVector3 pt0(evt.reco_beam_calo_startX,
                 evt.reco_beam_calo_startY,
                 evt.reco_beam_calo_startZ);
    TVector3 pt1(evt.reco_beam_calo_endX,
                 evt.reco_beam_calo_endY,
                 evt.reco_beam_calo_endZ);
    TVector3 dir = pt1 - pt0;
    dir = dir.Unit();

    if (evt.MC){
      beam_dx = (evt.reco_beam_calo_startX - beam_startX_mc)/beam_startX_rms_mc;
      beam_dy = (evt.reco_beam_calo_startY - beam_startY_mc)/beam_startY_rms_mc;
      beam_dz = (evt.reco_beam_calo_startZ - beam_startZ_mc)/beam_startZ_rms_mc;
      beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
      TVector3 beamdir(cos(beam_angleX_mc*TMath::Pi()/180),
                       cos(beam_angleY_mc*TMath::Pi()/180),
                       cos(beam_angleZ_mc*TMath::Pi()/180));
      beamdir = beamdir.Unit();
      beam_costh = dir.Dot(beamdir);
    }
    else{
      beam_dx = (evt.reco_beam_calo_startX - beam_startX_data)/beam_startX_rms_data;
      beam_dy = (evt.reco_beam_calo_startY - beam_startY_data)/beam_startY_rms_data;
      beam_dz = (evt.reco_beam_calo_startZ - beam_startZ_data)/beam_startZ_rms_data;
      beam_dxy = sqrt(pow(beam_dx,2) + pow(beam_dy,2));
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
  if (!evt.reco_beam_calibrated_dEdX_SCE->empty()){
    //dEdx_5cm = 0;
    //int nhits = 0;
    std::vector<double> vdEdx;
    for (int i = 0; i<evt.reco_beam_calibrated_dEdX_SCE->size(); ++i){
      //std::cout<<i<<" "<<reco_beam_resRange_SCE->back()-(*reco_beam_resRange_SCE)[i]<<" "<<(*reco_beam_calibrated_dEdX_SCE)[i]<<endl;
      //if (event == 78467) cout<<(*reco_beam_resRange_SCE)[i]<<" "<<(*reco_beam_calo_Z)[i]<<" "<<(*reco_beam_calibrated_dEdX_SCE)[i]<<endl;
      if ((*evt.reco_beam_resRange_SCE)[i]<5){
        vdEdx.push_back((*evt.reco_beam_calibrated_dEdX_SCE)[i]);
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

