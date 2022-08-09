#include "RecoDaughter.h"

RecoDaughter::RecoDaughter(){

  j_IsEmpty = true;
  j_PFP_true_byHits_PDG = -1;
  j_PFP_true_byHits_ID = -1;
  j_PFP_true_byHits_origin = -1;
  j_PFP_true_byHits_parID = -1;
  j_PFP_true_byHits_parPDG = -1;
  j_PFP_true_byHits_process = "";
  j_PFP_true_byHits_sharedHits = 0;
  j_PFP_true_byHits_emHits = 0;
  j_PFP_true_byHits_len = -9999.;
  j_PFP_true_byHits_startX = -9999.;
  j_PFP_true_byHits_startY = -9999.;
  j_PFP_true_byHits_startZ = -9999.;
  j_PFP_true_byHits_endX = -9999.;
  j_PFP_true_byHits_endY = -9999.;
  j_PFP_true_byHits_endZ = -9999.;
  j_PFP_true_byHits_startPx = -9999.;
  j_PFP_true_byHits_startPy = -9999.;
  j_PFP_true_byHits_startPz = -9999.;
  j_PFP_true_byHits_startP = -9999.;
  j_PFP_true_byHits_startE = -9999.;
  j_PFP_true_byHits_endProcess = "";
  j_PFP_true_byHits_purity = -9999.;
  j_PFP_true_byHits_completeness = -9999.;
  j_PFP_true_byE_PDG = -1;
  j_PFP_true_byE_len = -9999.;
  j_PFP_ID = -1;
  j_PFP_nHits = -1;
  j_PFP_nHits_collection = -1;
  j_PFP_trackScore = -9999.;
  j_PFP_emScore = -9999.;
  j_PFP_michelScore = -9999.;
  j_PFP_trackScore_collection = -9999.;
  j_PFP_emScore_collection = -9999.;
  j_PFP_michelScore_collection = -9999.;
  j_allTrack_ID = -1;
  j_allTrack_EField_SCE.clear();
  j_allTrack_resRange_SCE.clear();
  j_allTrack_resRange_SCE_plane0.clear();
  j_allTrack_resRange_SCE_plane1.clear();
  j_allTrack_calibrated_dEdX_SCE.clear();
  j_allTrack_calibrated_dEdX_SCE_plane0.clear();
  j_allTrack_calibrated_dEdX_SCE_plane1.clear();
  j_allTrack_Chi2_proton = -9999.;
  j_allTrack_Chi2_pion = -9999.;
  j_allTrack_Chi2_muon = -9999.;
  j_allTrack_Chi2_ndof = -9999.;
  j_allTrack_Chi2_ndof_pion = -9999.;
  j_allTrack_Chi2_ndof_muon = -9999.;
  j_allTrack_Theta = -9999.;
  j_allTrack_Phi = -9999.;
  j_allTrack_startDirX = -9999.;
  j_allTrack_startDirY = -9999.;
  j_allTrack_startDirZ = -9999.;
  j_allTrack_alt_len = -9999.;
  j_allTrack_startX = -9999.;
  j_allTrack_startY = -9999.;
  j_allTrack_startZ = -9999.;
  j_allTrack_endX = -9999.;
  j_allTrack_endY = -9999.;
  j_allTrack_endZ = -9999.;
  j_allTrack_vertex_michel_score = -9999.;
  j_allTrack_vertex_nHits = -1;
  j_pandora_type = -1;
  j_Beam_Cos = -9999.;
}

RecoDaughter::~RecoDaughter(){}

void RecoDaughter::SetIsEmpty(bool i_IsEmpty){ j_IsEmpty = i_IsEmpty; }

void RecoDaughter::Set_PFP_true_byHits_PDG(int i_PFP_true_byHits_PDG){ j_PFP_true_byHits_PDG = i_PFP_true_byHits_PDG; }
void RecoDaughter::Set_PFP_true_byHits_ID(int i_PFP_true_byHits_ID){ j_PFP_true_byHits_ID = i_PFP_true_byHits_ID; }
void RecoDaughter::Set_PFP_true_byHits_origin(int i_PFP_true_byHits_origin){ j_PFP_true_byHits_origin = i_PFP_true_byHits_origin; }
void RecoDaughter::Set_PFP_true_byHits_parID(int i_PFP_true_byHits_parID){ j_PFP_true_byHits_parID = i_PFP_true_byHits_parID; }
void RecoDaughter::Set_PFP_true_byHits_parPDG(int i_PFP_true_byHits_parPDG){ j_PFP_true_byHits_parPDG = i_PFP_true_byHits_parPDG; }
void RecoDaughter::Set_PFP_true_byHits_process(string i_PFP_true_byHits_process){ j_PFP_true_byHits_process = i_PFP_true_byHits_process; }
void RecoDaughter::Set_PFP_true_byHits_sharedHits(unsigned long i_PFP_true_byHits_sharedHits){ j_PFP_true_byHits_sharedHits = i_PFP_true_byHits_sharedHits; }
void RecoDaughter::Set_PFP_true_byHits_emHits(unsigned long i_PFP_true_byHits_emHits){ j_PFP_true_byHits_emHits = i_PFP_true_byHits_emHits; }
void RecoDaughter::Set_PFP_true_byHits_len(double i_PFP_true_byHits_len){ j_PFP_true_byHits_len = i_PFP_true_byHits_len; }
void RecoDaughter::Set_PFP_true_byHits_startX(double i_PFP_true_byHits_startX){ j_PFP_true_byHits_startX = i_PFP_true_byHits_startX; }
void RecoDaughter::Set_PFP_true_byHits_startY(double i_PFP_true_byHits_startY){ j_PFP_true_byHits_startY = i_PFP_true_byHits_startY; }
void RecoDaughter::Set_PFP_true_byHits_startZ(double i_PFP_true_byHits_startZ){ j_PFP_true_byHits_startZ = i_PFP_true_byHits_startZ; }
void RecoDaughter::Set_PFP_true_byHits_endX(double i_PFP_true_byHits_endX){ j_PFP_true_byHits_endX = i_PFP_true_byHits_endX; }
void RecoDaughter::Set_PFP_true_byHits_endY(double i_PFP_true_byHits_endY){ j_PFP_true_byHits_endY = i_PFP_true_byHits_endY; }
void RecoDaughter::Set_PFP_true_byHits_endZ(double i_PFP_true_byHits_endZ){ j_PFP_true_byHits_endZ = i_PFP_true_byHits_endZ; }
void RecoDaughter::Set_PFP_true_byHits_startPx(double i_PFP_true_byHits_startPx){ j_PFP_true_byHits_startPx = i_PFP_true_byHits_startPx; }
void RecoDaughter::Set_PFP_true_byHits_startPy(double i_PFP_true_byHits_startPy){ j_PFP_true_byHits_startPy = i_PFP_true_byHits_startPy; }
void RecoDaughter::Set_PFP_true_byHits_startPz(double i_PFP_true_byHits_startPz){ j_PFP_true_byHits_startPz = i_PFP_true_byHits_startPz; }
void RecoDaughter::Set_PFP_true_byHits_startP(double i_PFP_true_byHits_startP){ j_PFP_true_byHits_startP = i_PFP_true_byHits_startP; }
void RecoDaughter::Set_PFP_true_byHits_startE(double i_PFP_true_byHits_startE){ j_PFP_true_byHits_startE = i_PFP_true_byHits_startE; }
void RecoDaughter::Set_PFP_true_byHits_endProcess(string i_PFP_true_byHits_endProcess){ j_PFP_true_byHits_endProcess = i_PFP_true_byHits_endProcess; }
void RecoDaughter::Set_PFP_true_byHits_purity(double i_PFP_true_byHits_purity){ j_PFP_true_byHits_purity = i_PFP_true_byHits_purity; }
void RecoDaughter::Set_PFP_true_byHits_completeness(double i_PFP_true_byHits_completeness){ j_PFP_true_byHits_completeness = i_PFP_true_byHits_completeness; }
void RecoDaughter::Set_PFP_true_byE_PDG(int i_PFP_true_byE_PDG){ j_PFP_true_byE_PDG = i_PFP_true_byE_PDG; }
void RecoDaughter::Set_PFP_true_byE_len(double i_PFP_true_byE_len){ j_PFP_true_byE_len = i_PFP_true_byE_len; }
void RecoDaughter::Set_PFP_ID(int i_PFP_ID){ j_PFP_ID = i_PFP_ID; }
void RecoDaughter::Set_PFP_nHits(int i_PFP_nHits){ j_PFP_nHits = i_PFP_nHits; }
void RecoDaughter::Set_PFP_nHits_collection(int i_PFP_nHits_collection){ j_PFP_nHits_collection = i_PFP_nHits_collection; }
void RecoDaughter::Set_PFP_trackScore(double i_PFP_trackScore){ j_PFP_trackScore = i_PFP_trackScore; }
void RecoDaughter::Set_PFP_emScore(double i_PFP_emScore){ j_PFP_emScore = i_PFP_emScore; }
void RecoDaughter::Set_PFP_michelScore(double i_PFP_michelScore){ j_PFP_michelScore = i_PFP_michelScore; }
void RecoDaughter::Set_PFP_trackScore_collection(double i_PFP_trackScore_collection){ j_PFP_trackScore_collection = i_PFP_trackScore_collection; }
void RecoDaughter::Set_PFP_emScore_collection(double i_PFP_emScore_collection){ j_PFP_emScore_collection = i_PFP_emScore_collection; }
void RecoDaughter::Set_PFP_michelScore_collection(double i_PFP_michelScore_collection){ j_PFP_michelScore_collection = i_PFP_michelScore_collection; }
void RecoDaughter::Set_allTrack_ID(int i_allTrack_ID){ j_allTrack_ID = i_allTrack_ID; }
void RecoDaughter::Set_allTrack_EField_SCE(vector<double> i_allTrack_EField_SCE){ j_allTrack_EField_SCE = i_allTrack_EField_SCE; }
void RecoDaughter::Set_allTrack_resRange_SCE(vector<double> i_allTrack_resRange_SCE){ j_allTrack_resRange_SCE = i_allTrack_resRange_SCE; }
void RecoDaughter::Set_allTrack_resRange_SCE_plane0(vector<double> i_allTrack_resRange_SCE_plane0){ j_allTrack_resRange_SCE_plane0 = i_allTrack_resRange_SCE_plane0; }
void RecoDaughter::Set_allTrack_resRange_SCE_plane1(vector<double> i_allTrack_resRange_SCE_plane1){ j_allTrack_resRange_SCE_plane1 = i_allTrack_resRange_SCE_plane1; }
void RecoDaughter::Set_allTrack_calibrated_dEdX_SCE(vector<double> i_allTrack_calibrated_dEdX_SCE){ j_allTrack_calibrated_dEdX_SCE = i_allTrack_calibrated_dEdX_SCE; }
void RecoDaughter::Set_allTrack_calibrated_dEdX_SCE_plane0(vector<double> i_allTrack_calibrated_dEdX_SCE_plane0){ j_allTrack_calibrated_dEdX_SCE_plane0 = i_allTrack_calibrated_dEdX_SCE_plane0; }
void RecoDaughter::Set_allTrack_calibrated_dEdX_SCE_plane1(vector<double> i_allTrack_calibrated_dEdX_SCE_plane1){ j_allTrack_calibrated_dEdX_SCE_plane1 = i_allTrack_calibrated_dEdX_SCE_plane1; }
void RecoDaughter::Set_allTrack_Chi2_proton(double i_allTrack_Chi2_proton){ j_allTrack_Chi2_proton = i_allTrack_Chi2_proton; }
void RecoDaughter::Set_allTrack_Chi2_pion(double i_allTrack_Chi2_pion){ j_allTrack_Chi2_pion = i_allTrack_Chi2_pion; }
void RecoDaughter::Set_allTrack_Chi2_muon(double i_allTrack_Chi2_muon){ j_allTrack_Chi2_muon = i_allTrack_Chi2_muon; }
void RecoDaughter::Set_allTrack_Chi2_ndof(double i_allTrack_Chi2_ndof){ j_allTrack_Chi2_ndof = i_allTrack_Chi2_ndof; }
void RecoDaughter::Set_allTrack_Chi2_ndof_pion(double i_allTrack_Chi2_ndof_pion){ j_allTrack_Chi2_ndof_pion = i_allTrack_Chi2_ndof_pion; }
void RecoDaughter::Set_allTrack_Chi2_ndof_muon(double i_allTrack_Chi2_ndof_muon){ j_allTrack_Chi2_ndof_muon = i_allTrack_Chi2_ndof_muon; }
void RecoDaughter::Set_allTrack_Theta(double i_allTrack_Theta){ j_allTrack_Theta = i_allTrack_Theta; }
void RecoDaughter::Set_allTrack_Phi(double i_allTrack_Phi){ j_allTrack_Phi = i_allTrack_Phi; }
void RecoDaughter::Set_allTrack_startDirX(double i_allTrack_startDirX){ j_allTrack_startDirX = i_allTrack_startDirX; }
void RecoDaughter::Set_allTrack_startDirY(double i_allTrack_startDirY){ j_allTrack_startDirY = i_allTrack_startDirY; }
void RecoDaughter::Set_allTrack_startDirZ(double i_allTrack_startDirZ){ j_allTrack_startDirZ = i_allTrack_startDirZ; }
void RecoDaughter::Set_allTrack_alt_len(double i_allTrack_alt_len){ j_allTrack_alt_len = i_allTrack_alt_len; }
void RecoDaughter::Set_allTrack_startX(double i_allTrack_startX){ j_allTrack_startX = i_allTrack_startX; }
void RecoDaughter::Set_allTrack_startY(double i_allTrack_startY){ j_allTrack_startY = i_allTrack_startY; }
void RecoDaughter::Set_allTrack_startZ(double i_allTrack_startZ){ j_allTrack_startZ = i_allTrack_startZ; }
void RecoDaughter::Set_allTrack_endX(double i_allTrack_endX){ j_allTrack_endX = i_allTrack_endX; }
void RecoDaughter::Set_allTrack_endY(double i_allTrack_endY){ j_allTrack_endY = i_allTrack_endY; }
void RecoDaughter::Set_allTrack_endZ(double i_allTrack_endZ){ j_allTrack_endZ = i_allTrack_endZ; }
void RecoDaughter::Set_allTrack_vertex_michel_score(double i_allTrack_vertex_michel_score){ j_allTrack_vertex_michel_score = i_allTrack_vertex_michel_score; }
void RecoDaughter::Set_allTrack_vertex_nHits(int i_allTrack_vertex_nHits){ j_allTrack_vertex_nHits = i_allTrack_vertex_nHits; }
void RecoDaughter::Set_pandora_type(int i_pandora_type){ j_pandora_type = i_pandora_type; }
void RecoDaughter::Set_Beam_Cos(double i_Beam_Cos) { j_Beam_Cos = i_Beam_Cos; }
void RecoDaughter::Set_Beam_Dist(double i_Beam_Dist) { j_Beam_Dist = i_Beam_Dist; }
