#include "BeamTreeReproducer.h"
#include "BeamNtuple.h"
#include "BeamVirtualDetector.h"

using namespace std;

BeamTreeReproducer::BeamTreeReproducer(){
}

void BeamTreeReproducer::BookHistograms(){
  //cout << "[BeamTreeReproducer::BookHistograms] Start" << endl;
  Hist.outfile = TFile::Open(fOutputFileName.c_str(), "recreate");

  // == Histograms for normalization
  h_cutflow = new TH1D("Cutflow", "Cutflow", 20, 0., 20.);
}

void BeamTreeReproducer::ProcessEvent(const BeamNtuple & evt){
}

void BeamTreeReproducer::FillHistograms(const BeamNtuple & evt){
  //cout << "[BeamTreeReproducer::FillHistograms] Start" <<endl;

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

void BeamTreeReproducer::FillHistograms(const BeamVirtualDetector & evt, TString detector_str){

  int PID = evt.PDGid;
  TString PID_str = Form("%d", PID);

  double P = sqrt(pow(evt.Px, 2) + pow(evt.Py, 2) + pow(evt.Pz, 2));
  double InitKE = evt.InitKE;
  double z_in_meter = (evt.z - 681526.937500) / 1000.;
  //cout << "[BeamTreeReproducer::FillHistograms] " << detector_str << Form(", (x, y, z) =  (%f, %f, %f", evt.x, evt.y, evt.z - 681526.937500) << endl;
  /*
  if(PID = -13){
    cout << "[BeamTreeReproducer::FillHistograms] " << endl;
    cout << "[BeamTreeReproducer::FillHistograms] " << detector_str << Form(" (PolX, PolY, PolZ) =  (%f, %f, %f", evt.PolX, evt.PolY, evt.PolZ) << endl;
    cout << "[BeamTreeReproducer::FillHistograms] " << detector_str << Form(" (Bx, By, Bz) =  (%f, %f, %f", evt.Bx, evt.By, evt.Bz) << endl;
    cout << "[BeamTreeReproducer::FillHistograms] " << detector_str << Form(" (Ex, Ey, Ez) =  (%f, %f, %f", evt.Ex, evt.Ey, evt.Ez) << endl;
    cout << "[BeamTreeReproducer::FillHistograms] " << detector_str << Form(", (x, y, z) =  (%f, %f, %f", evt.x, evt.y, evt.z) << endl;
    cout << "[BeamTreeReproducer::FillHistograms] " << detector_str << Form(" (InitX, InitY, InitZ) =  (%f, %f, %f", evt.InitX, evt.InitY, evt.InitZ) << endl;
    cout << "[BeamTreeReproducer::FillHistograms] " << detector_str << "PID : " << PID << ", InitKE : " << evt.InitKE << ", P : " << P << endl;
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

  double xy = sqrt(evt.x *evt.x + evt.y * evt.y);
  Hist.FillHist("VirtualDetector_xz_" + detector_str + "_" + PID_str, z_in_meter, evt.x, 1., 400., 0., 40., 400., -200., 200.);
  Hist.FillHist("VirtualDetector_yz_" + detector_str + "_" + PID_str, z_in_meter, evt.y, 1., 400., 0.,40., 400., -200., 200.);
  Hist.FillHist("VirtualDetector_P_vs_xy_"  + detector_str + "_" + PID_str, P, xy, 1., 1500., 0., 1500., 200., 0., 200.);
}

void BeamTreeReproducer::SaveHistograms(){
  Hist.WriteHist();
  //Hist.outfile->cd();
  //outputFile->Write();
}

void BeamTreeReproducer::Find_Matched_Event(BeamVirtualDetector & detector_tree, int Evt_ID, int Trk_ID, int PDG_ID, TString detector_str){

  Long64_t this_nentries = detector_tree.fChain->GetEntries();
 
  int N_match_Evt = 0;
  int N_match_Evt_Trk = 0;
  int PDG_ID_match_Evt = 0;
  int PDG_ID_match_Evt_Trk = 0;
  for(Long64_t kentry = 0; kentry < this_nentries; kentry++){
    //if (kentry%10000==0) cout << "[" << detector_str << "] " << kentry << "/" << this_nentries << endl;
    Long64_t ientry = detector_tree.LoadTree(kentry);
    if (ientry < 0) break;
    detector_tree.fChain->GetEntry(kentry);
    int this_Evt_ID = detector_tree.EventID;
    int this_Trk_ID = detector_tree.TrackID;
    int this_PDG_ID = detector_tree.PDGid;
    if(Evt_ID == this_Evt_ID){
      N_match_Evt++;
      PDG_ID_match_Evt = this_PDG_ID;
    }
    if(Evt_ID == this_Evt_ID && Trk_ID == this_Trk_ID){
      N_match_Evt_Trk++;
      PDG_ID_match_Evt_Trk = this_PDG_ID;
    }
  }
  TString out = Form("N_match_Evt\t%d (%d)\tN_match_Evt_Trk\t%d (%d)", N_match_Evt, PDG_ID_match_Evt, N_match_Evt_Trk, PDG_ID_match_Evt_Trk);
  cout << "[" << detector_str << "]\t" << out << endl;

  return;
}

//void BeamTreeReproducer::Run(vector<BeamVirtualDetector> & branch_vector, vector<TString> branch_name_vector, Long64_t nentries = -1){
void BeamTreeReproducer::Run(BeamVirtualDetector & evt_AfterTarget, BeamVirtualDetector & evt_TOF1, BeamVirtualDetector & evt_COLL1, BeamVirtualDetector & evt_BPROF1,
			     BeamVirtualDetector & evt_BPROF2, BeamVirtualDetector & evt_BPROF3, BeamVirtualDetector & evt_TRIG1, BeamVirtualDetector & evt_BPROFEXT,
			     BeamVirtualDetector & evt_BPROF4, BeamVirtualDetector & evt_TRIG2, Long64_t nentries = -1){

  BookHistograms();
  int nentries_BPROF1 = evt_BPROF1.fChain->GetEntries();

  if (nentries == -1) nentries = evt_TRIG2.fChain->GetEntries();
  cout << "evt_AfterTarget nentries : " << evt_AfterTarget.fChain->GetEntries() << endl;
  cout << "evt_TOF1 nentries : " << evt_TOF1.fChain->GetEntries() << endl;
  cout << "evt_COLL1 nentries : " << evt_COLL1.fChain->GetEntries() << endl;
  cout << "evt_BPROF1 nentries : " << evt_BPROF1.fChain->GetEntries() << endl;
  cout << "evt_BPROF2 nentries : " << evt_BPROF2.fChain->GetEntries() << endl;
  cout << "evt_BPROF3 nentries : " << evt_BPROF3.fChain->GetEntries() << endl;
  cout << "evt_TRIG1 nentries : " << evt_TRIG1.fChain->GetEntries() << endl;
  cout << "evt_BPROFEXT nentries : " << evt_BPROFEXT.fChain->GetEntries() << endl;
  cout << "evt_BPROF4 nentries : " << evt_BPROF4.fChain->GetEntries() << endl;
  cout << "evt_TRIG2 nentries : " << evt_TRIG2.fChain->GetEntries() << endl;

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    if (jentry%1000==0) std::cout<<"[evt_BPROF4] " << jentry<<"/"<< nentries << " from total : " << evt_BPROF4.fChain->GetEntries() << std::endl;
    Long64_t ientry = evt_BPROF4.LoadTree(jentry);
    if (ientry < 0) break;
    
    nb = evt_BPROF4.fChain->GetEntry(jentry); nbytes += nb;
    int Evt_ID = evt_BPROF4.EventID;
    int Trk_ID = evt_BPROF4.TrackID;
    int PDG_ID = evt_BPROF4.PDGid;
    if(abs(PDG_ID) != 13) continue;

    cout << "--------------------" << endl;
    cout << "[BPROF4] " << jentry << "\tEvt_ID : " << Evt_ID << ", Trk_ID : " << Trk_ID << ", PDG_ID : " << PDG_ID << endl;
    
    // == Loop for upstream detectors
    Find_Matched_Event(evt_BPROF1, Evt_ID, Trk_ID, PDG_ID, "BPROF1");
    Find_Matched_Event(evt_BPROF2, Evt_ID, Trk_ID, PDG_ID, "BPROF2");
    Find_Matched_Event(evt_BPROF3, Evt_ID, Trk_ID, PDG_ID, "BPROF3");

    /*
    int N_match_Evt = 0;
    int N_match_Evt_Trk = 0;
    int PDG_ID_match_Evt = 0;
    int PDG_ID_match_Evt_Trk = 0;
    for(Long64_t kentry = 0; kentry < nentries_BPROF1; kentry++){
      //if (kentry%10000==0) cout << "[" << detector_str << "] " << kentry << "/" << this_nentries << endl;                                                                                  
      Long64_t ientry = evt_BPROF1.LoadTree(kentry);
      if (ientry < 0) break;
      evt_BPROF1.fChain->GetEntry(kentry);
      int this_Evt_ID = evt_BPROF1.EventID;
      int this_Trk_ID = evt_BPROF1.TrackID;
      int this_PDG_ID = evt_BPROF1.PDGid;
      if(Evt_ID == this_Evt_ID){
	N_match_Evt++;
	PDG_ID_match_Evt = this_PDG_ID;
      }
      if(Evt_ID == this_Evt_ID && Trk_ID == this_Trk_ID){
	N_match_Evt_Trk++;
	PDG_ID_match_Evt_Trk = this_PDG_ID;
      }
    }
    TString out = Form("N_match_Evt\t%d (%d)\tN_match_Evt_Trk\t%d (%d)", N_match_Evt, PDG_ID_match_Evt, N_match_Evt_Trk, PDG_ID_match_Evt_Trk);
    cout << "[BPROF1]\t" << out << endl;
    */

    h_cutflow -> Fill(0.5);
    //FillHistograms(evt_BPROF4, "BPROF4");
  }

  // == Save
  SaveHistograms();
}
