#include "HadAna.h"

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
      htrue_beam_endZ[i][j] = new TH1D(Form("htrue_beam_endZ_%d_%d",i,j),Form("true_beam_endZ,%s,%s", cutName[i], parTypeName[j]), 100, -100, 900);
      htrue_beam_endZ[i][j]->Sumw2();
    }
  }
}

void HadAna::FillHistograms(int cut){
  
  if (cut>=0 && cut < nCuts){
    htrue_beam_endZ[cut][0]->Fill(true_beam_endZ_SCE);
    int partype = GetParType();
    if (partype>=1 && partype < nParTypes+1){
      htrue_beam_endZ[cut][partype]->Fill(true_beam_endZ_SCE);
    }
  }
}

void HadAna::SaveHistograms(){
  outputFile->Write();
}
