#include "BeamSampleAna.h"
#include "BeamNtuple.h"

using namespace std;

BeamSampleAna::BeamSampleAna(){
  //hadana.InitPi();
}

void BeamSampleAna::BookHistograms(){
  //cout << "[BeamSampleAna::BookHistograms] Start" << endl;
  Hist.outfile = TFile::Open(fOutputFileName.c_str(), "recreate");

  // == Histograms for normalization
  h_cutflow = new TH1D("Cutflow", "Cutflow", 20, 0., 20.);

}

void BeamSampleAna::ProcessEvent(const BeamNtuple & evt){
}

void BeamSampleAna::FillHistograms(const BeamNtuple & evt){
  cout << "[BeamSampleAna::FillHistograms] Start" <<endl;
  cout << "AfterTarget_Pz : " << evt.AfterTarget_Pz << endl;
  double P_AfterTarget = sqrt(pow(evt.AfterTarget_Px, 2) + pow(evt.AfterTarget_Py, 2) + pow(evt.AfterTarget_Pz, 2));
  Hist.FillHist("P_AfterTarget", P_AfterTarget, 1., 5000., 0., 5000.); 
  double beamP_scale = 1.0;
  //if(evt.MC) beamP_scale = 0.5;
}

void BeamSampleAna::SaveHistograms(){
  Hist.WriteHist();
  //Hist.outfile->cd();
  //outputFile->Write();
}

void BeamSampleAna::Run(BeamNtuple & evt, Long64_t nentries=-1){

  BookHistograms();

  //Long64_t nentries = evt.fChain->GetEntriesFast();
  if (nentries == -1) nentries = evt.fChain->GetEntries();
  
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    //if (jentry%100000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    if (jentry%1000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = evt.LoadTree(jentry);
    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;
    h_cutflow -> Fill(0.5);
    ProcessEvent(evt);
    FillHistograms(evt);
  }
  SaveHistograms();
}

