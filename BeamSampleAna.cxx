#include "BeamSampleAna.h"
#include "BeamNtuple.h"
#include "BeamVirtualDetector.h"

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
  //cout << "[BeamSampleAna::FillHistograms] Start" <<endl;

  // == PIDs
  int PID_AfterTarget = evt.AfterTarget_PDGid;
  int PID_TOF1 = evt.TOF1_PDGid;
  int PID_COLL1 = evt.COLL1_PDGid;
  int PID_BPROF1 = evt.BPROF1_PDGid;
  int PID_BPROF2 = evt.BPROF2_PDGid;
  int PID_BPROF3 = evt.BPROF3_PDGid;
  int PID_TRIG1 = evt.TRIG1_PDGid;
  int PID_BPROFEXT = evt.BPROFEXT_PDGid;
  int PID_BPROF4 = evt.BPROF4_PDGid;
  int PID_TRIG2 = evt.TRIG2_PDGid;
  int PID_NP04front = evt.NP04front_PDGid;
  int PID_NP04FieldCage = evt.NP04FieldCage_PDGid;

  //if(abs(PID_AfterTarget) != 13 && PID_AfterTarget != 211) return; // == Select only pion and muon beam
  if(PID_AfterTarget != PID_NP04front){
    //cout << "PID_AfterTarget : " << PID_AfterTarget << ", PID_NP04front: " << PID_NP04front << endl;
  }
  
  TString PID_AfterTarget_str = Form("%d", abs(PID_AfterTarget));

  // == P
  double P_AfterTarget = sqrt(pow(evt.AfterTarget_Px, 2) + pow(evt.AfterTarget_Py, 2) + pow(evt.AfterTarget_Pz, 2));
  double P_TOF1 = sqrt(pow(evt.TOF1_Px, 2) + pow(evt.TOF1_Py, 2) + pow(evt.TOF1_Pz, 2));
  double P_COLL1 = sqrt(pow(evt.COLL1_Px, 2) + pow(evt.COLL1_Py, 2) + pow(evt.COLL1_Pz, 2));
  double P_BPROF1 = sqrt(pow(evt.BPROF1_Px, 2) + pow(evt.BPROF1_Py, 2) + pow(evt.BPROF1_Pz, 2));
  double P_BPROF2 = sqrt(pow(evt.BPROF2_Px, 2) + pow(evt.BPROF2_Py, 2) + pow(evt.BPROF2_Pz, 2));
  double P_BPROF3 = sqrt(pow(evt.BPROF3_Px, 2) + pow(evt.BPROF3_Py, 2) + pow(evt.BPROF3_Pz, 2));
  double P_TRIG1 = sqrt(pow(evt.TRIG1_Px, 2) + pow(evt.TRIG1_Py, 2) + pow(evt.TRIG1_Pz, 2));
  double P_BPROFEXT = sqrt(pow(evt.BPROFEXT_Px, 2) + pow(evt.BPROFEXT_Py, 2) + pow(evt.BPROFEXT_Pz, 2));
  double P_BPROF4 = sqrt(pow(evt.BPROF4_Px, 2) + pow(evt.BPROF4_Py, 2) + pow(evt.BPROF4_Pz, 2));
  double P_TRIG2 = sqrt(pow(evt.TRIG2_Px, 2) + pow(evt.TRIG2_Py, 2) + pow(evt.TRIG2_Pz, 2));
  double P_NP04front = sqrt(pow(evt.NP04front_Px, 2) + pow(evt.NP04front_Py, 2) + pow(evt.NP04front_Pz, 2));
  double P_NP04FieldCage = sqrt(pow(evt.NP04FieldCage_Px, 2) + pow(evt.NP04FieldCage_Py, 2) + pow(evt.NP04FieldCage_Pz, 2));

  // == After Target
  Hist.FillHist("P_AfterTarget", P_AfterTarget, 1., 5000., 0., 5000.); 
  Hist.FillHist("P_AfterTarget_PID" + PID_AfterTarget_str, P_AfterTarget, 1., 5000., 0., 5000.);
  /*
  cout << "===========" << endl;
  cout << "AfterTarget_x : " << evt.AfterTarget_x << endl;
  cout << "AfterTarget_y : " << evt.AfterTarget_y << endl;
  cout << "AfterTarget_z : " << evt.AfterTarget_z << endl;
  cout << "AfterTarget_InitX : " << evt.AfterTarget_InitX << endl;
  cout << "AfterTarget_InitY : " << evt.AfterTarget_InitY << endl;
  cout << "AfterTarget_InitZ : " << evt.AfterTarget_InitZ << endl;
  cout << "TRIG2_x : " << evt.TRIG2_x << endl;
  cout << "TRIG2_y : " << evt.TRIG2_y << endl;
  cout << "TRIG2_z : " << evt.TRIG2_z << endl;
  cout << "NP04front_x : " << evt.NP04front_x << endl;
  cout << "NP04front_y : " << evt.NP04front_y << endl;
  cout << "NP04front_z : " << evt.NP04front_z << endl;
  */
  // == TOF1
  Hist.FillHist("P_TOF1", P_TOF1, 1., 5000., 0., 5000.);
  Hist.FillHist("P_TOF1_PID" + PID_AfterTarget_str, P_TOF1, 1., 5000., 0., 5000.);

  // == COLL1
  Hist.FillHist("P_COLL1", P_COLL1, 1., 5000., 0., 5000.);
  Hist.FillHist("P_COLL1_PID" + PID_AfterTarget_str, P_COLL1, 1., 5000., 0., 5000.);

  // == BPROF1
  Hist.FillHist("P_BPROF1", P_BPROF1, 1., 5000., 0., 5000.);
  Hist.FillHist("P_BPROF1_PID" + PID_AfterTarget_str, P_BPROF1, 1., 5000., 0., 5000.);

  // == BPROF2
  Hist.FillHist("P_BPROF2", P_BPROF2, 1., 5000., 0., 5000.);
  Hist.FillHist("P_BPROF2_PID" + PID_AfterTarget_str, P_BPROF2, 1., 5000., 0., 5000.);

  // == BPROF3
  Hist.FillHist("P_BPROF3", P_BPROF3, 1., 5000., 0., 5000.);
  Hist.FillHist("P_BPROF3_PID" + PID_AfterTarget_str, P_BPROF3, 1., 5000., 0., 5000.);

  // == TRIG1
  Hist.FillHist("P_TRIG1", P_TRIG1, 1., 5000., 0., 5000.);
  Hist.FillHist("P_TRIG1_PID" + PID_AfterTarget_str, P_TRIG1, 1., 5000., 0., 5000.);

  // == BPROFEXT
  Hist.FillHist("P_BPROFEXT", P_BPROFEXT, 1., 5000., 0., 5000.);
  Hist.FillHist("P_BPROFEXT_PID" + PID_AfterTarget_str, P_BPROFEXT, 1., 5000., 0., 5000.);

  // == BPROF4
  Hist.FillHist("P_BPROF4", P_BPROF4, 1., 5000., 0., 5000.);
  Hist.FillHist("P_BPROF4_PID" + PID_AfterTarget_str, P_BPROF4, 1., 5000., 0., 5000.);

  // == TRIG2
  Hist.FillHist("P_TRIG2", P_TRIG2, 1., 5000., 0., 5000.);
  Hist.FillHist("P_TRIG2_PID" + PID_AfterTarget_str, P_TRIG2, 1., 5000., 0., 5000.);

  // == NP04front
  Hist.FillHist("P_NP04front", P_NP04front, 1., 5000., 0., 5000.);
  Hist.FillHist("P_NP04front_PID" + PID_AfterTarget_str, P_NP04front, 1., 5000., 0., 5000.);

  // == NP04FieldCage
  Hist.FillHist("P_NP04FieldCage", P_NP04FieldCage, 1., 5000., 0., 5000.);
  Hist.FillHist("P_NP04FieldCage_PID" + PID_AfterTarget_str, P_NP04FieldCage, 1., 5000., 0., 5000.);

}

void BeamSampleAna::FillHistograms(const BeamVirtualDetector & evt, TString detector_str){

  int PID = evt.PDGid;
  TString PID_str = Form("%d", abs(PID));

  double P = sqrt(pow(evt.Px, 2) + pow(evt.Py, 2) + pow(evt.Pz, 2));
  double InitKE = evt.InitKE;
  double z_in_meter = (evt.z - 681526.937500) / 1000.;
  //cout << "[BeamSampleAna::FillHistograms] " << detector_str << Form(", (x, y, z) =  (%f, %f, %f", evt.x, evt.y, evt.z - 681526.937500) << endl;
  /*  
  if(PID = -13){
    cout << "[BeamSampleAna::FillHistograms] " << endl;
    cout << "[BeamSampleAna::FillHistograms] " << detector_str << Form(" (PolX, PolY, PolZ) =  (%f, %f, %f", evt.PolX, evt.PolY, evt.PolZ) << endl;
    cout << "[BeamSampleAna::FillHistograms] " << detector_str << Form(" (Bx, By, Bz) =  (%f, %f, %f", evt.Bx, evt.By, evt.Bz) << endl;
    cout << "[BeamSampleAna::FillHistograms] " << detector_str << Form(" (Ex, Ey, Ez) =  (%f, %f, %f", evt.Ex, evt.Ey, evt.Ez) << endl;
    cout << "[BeamSampleAna::FillHistograms] " << detector_str << Form(", (x, y, z) =  (%f, %f, %f", evt.x, evt.y, evt.z) << endl;
    cout << "[BeamSampleAna::FillHistograms] " << detector_str << Form(" (InitX, InitY, InitZ) =  (%f, %f, %f", evt.InitX, evt.InitY, evt.InitZ) << endl;
    cout << "[BeamSampleAna::FillHistograms] " << detector_str << "PID : " << PID << ", InitKE : " << evt.InitKE << ", P : " << P << endl;
  }
  */
  double m_muon = 105.658;
  double m_pion = 139.57;
  double E_muon_from_P = sqrt(P*P + m_muon * m_muon);
  double E_pion_from_P = sqrt(P*P + m_pion * m_pion); 
  double E_muon_from_KE = InitKE + m_muon;
  double E_pion_from_KE = InitKE + m_pion;
  Hist.FillHist("VirtualDetector_P_" + detector_str, P, 1., 10000., 0., 100000.);
  Hist.FillHist("VirtualDetector_P_" + detector_str + "_" + PID_str, P, 1., 10000., 0., 100000.);
  Hist.FillHist("VirtualDetector_P_vs_InitKE_" + detector_str + "_" + PID_str, P, InitKE, 1., 400., 0., 2000., 1000., 0., 5000.);
  Hist.FillHist("VirtualDetector_E_muon_from_P_vs_E_muon_from_KE_" + detector_str + "_" + PID_str, E_muon_from_P, E_muon_from_KE, 1., 400., 0., 2000., 1000., 0., 5000.);
  Hist.FillHist("VirtualDetector_E_muon_from_P_vs_E_pion_from_KE_" + detector_str + "_" + PID_str, E_muon_from_P, E_pion_from_KE, 1., 400., 0., 2000., 1000., 0., 5000.);
  Hist.FillHist("VirtualDetector_E_pion_from_P_vs_E_pion_from_KE_" + detector_str + "_" + PID_str, E_pion_from_P, E_pion_from_KE, 1., 400., 0., 2000., 1000., 0., 5000.);
  Hist.FillHist("VirtualDetector_xz_" + detector_str + "_" + PID_str, z_in_meter, evt.x, 1., 400., 0., 40., 400., -200., 200.);
  Hist.FillHist("VirtualDetector_yz_" + detector_str + "_" + PID_str, z_in_meter, evt.y, 1., 400., 0.,40., 400., -200., 200.);
  
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
  
  //for (Long64_t jentry=0; jentry<nentries;jentry++) {
  for (Long64_t jentry=0; jentry<1;jentry++) { 
   //if (jentry%100000==0) std::cout<<jentry<<"/"<<nentries<<std::endl;
    if (jentry%1000==0) std::cout<<"[GoodParticle] " <<jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = evt.LoadTree(jentry);
    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;
    h_cutflow -> Fill(0.5);
    ProcessEvent(evt);
    //FillHistograms(evt);
  }
  SaveHistograms();
}

void BeamSampleAna::Run(BeamVirtualDetector & evt, Long64_t nentries=-1, TString detector_str = ""){

  BookHistograms();
  if (nentries == -1) nentries = evt.fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%1000==0) std::cout<<"[" << detector_str << "] " << jentry<<"/"<<nentries<<std::endl;
    Long64_t ientry = evt.LoadTree(jentry);
    if (ientry < 0) break;
    nb = evt.fChain->GetEntry(jentry);   nbytes += nb;
    h_cutflow -> Fill(0.5);
    FillHistograms(evt, detector_str);
  }

  SaveHistograms();
}
