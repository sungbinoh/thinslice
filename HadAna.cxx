#include "HadAna.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "util.h"
#include <iostream>

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

  if (MC){
    double xlow = -3.,  xhigh = 7.,  ylow = -8.,  yhigh = 7.;
    double zlow = 27.5,  zhigh = 32.5,  coslow = 0.93;

    if ( beam_dx < xlow )
      return false;
    
    if ( beam_dx > xhigh )
      return false;
    
    if ( beam_dy < ylow )
      return false;
    
    if ( beam_dy > yhigh )
      return false;
    
    if (reco_beam_startZ < zlow || zhigh < reco_beam_startZ)
      return false;
    
    if ( beam_costh < coslow)
      return false;
    
    return true;
  }
  else{
    double data_xlow = 0., data_xhigh = 10., data_ylow= -5.;
    double data_yhigh= 10., data_zlow=30., data_zhigh=35., data_coslow=.93;

    if( (beam_dx < data_xlow) || (beam_dx > data_xhigh) )
      return false;
    
    if ( (beam_dy < data_ylow) || (beam_dy > data_yhigh) )
      return false;
    
    if ( (reco_beam_startZ < data_zlow) || (reco_beam_startZ > data_zhigh) )
      return false;
    
    if (beam_costh < data_coslow)
      return false;

    return true;
  }
    
};

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

  partype = -1;

  median_dEdx = -1;
  daughter_michel_score = 0;
  if (!reco_beam_calo_wire->empty()){

    median_dEdx = TMath::Median(reco_beam_calibrated_dEdX_SCE->size(), &(*reco_beam_calibrated_dEdX_SCE)[0]);
    int nhits = 0;
    for (size_t i = 0; i<reco_daughter_PFP_michelScore_collection->size(); ++i){
      nhits += (*reco_daughter_PFP_nHits_collection)[i];
      daughter_michel_score += (*reco_daughter_PFP_michelScore_collection)[0] * (*reco_daughter_PFP_nHits_collection)[i];
    }
    if (nhits) daughter_michel_score/=nhits;
    else daughter_michel_score = -999;
  }

  if (MC){
    double projectX = (true_beam_startX + -1*true_beam_startZ*(true_beam_startDirX/true_beam_startDirZ) );
    double projectY = (true_beam_startY + -1*true_beam_startZ*(true_beam_startDirY/true_beam_startDirZ) );
    beam_costh = true_beam_startDirX*reco_beam_trackDirX + true_beam_startDirY*reco_beam_trackDirY + true_beam_startDirZ*reco_beam_trackDirZ;
    beam_dx = reco_beam_startX - projectX;
    beam_dy = reco_beam_startY - projectY;
  }
  else{
    beam_dx = reco_beam_startX - beam_inst_X;
    beam_dy = reco_beam_startY - beam_inst_Y;
    beam_costh = beam_inst_dirX*reco_beam_trackDirX 
      + beam_inst_dirY*reco_beam_trackDirY 
      + beam_inst_dirZ*reco_beam_trackDirZ;
  }
  partype = GetParType();
}

