#include "ThinSlice.h"

void ThinSlice::BookHistograms(){

   for (int i = 0; i<nthinslices; ++i){
     reco_incE[i] = new TH1D(Form("reco_incE_%d",i),Form("Reco incident energy, %f<= z <%f (cm)",i*thinslicewidth, (i+1)*thinslicewidth), nbinse, 0, 1200.);
     true_incE[i] = new TH1D(Form("true_incE_%d",i),Form("True incident energy, %f<= z <%f (cm)",i*thinslicewidth, (i+1)*thinslicewidth), nbinse, 0, 1200.);
   }

   for (int i = 0; i<nthinslices; ++i){
     true_interactions[i] = 0;
     true_incidents[i] = 0;
   }
}

void ThinSlice::ProcessEvent(const HadAna & evt){

  if (evt.MC){
    true_sliceID = int(evt.true_beam_endZ/thinslicewidth);
    if (evt.true_beam_PDG == 211){
      for (int i = 0; i<=true_sliceID; ++i){
        if (i<nslices) ++true_incidents[i];
      }
    }
    if ((*evt.true_beam_endProcess) == "pi+Inelastic"){
      if (true_sliceID < nslices && true_sliceID>=0){
        ++true_interactions[true_sliceID];
      }
      // Reco info
      if (!(evt.reco_beam_calo_wire->empty())){
        std::vector<std::vector<double>> vincE(nslices);
        for (size_t i = 0; i<evt.reco_beam_calo_wire->size(); ++i){
          int this_sliceID = int((*evt.reco_beam_calo_wire_z)[i]/thinslicewidth);
          if (this_sliceID>=nthinslices) continue;
          if (this_sliceID<0) continue;
          double this_incE = (*evt.reco_beam_incidentEnergies)[i];
          vincE[this_sliceID].push_back(this_incE);
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
      }

      // True info
      if (!(evt.true_beam_traj_Z_SCE->empty())){
        std::vector<std::vector<double>> vincE(nslices);
        for (size_t i = 0; i<evt.true_beam_traj_Z->size()-1; ++i){//last point always has KE = 0
          int this_sliceID = int((*evt.true_beam_traj_Z)[i]/thinslicewidth);
          double this_incE = (*evt.true_beam_traj_KE)[i];
          if (this_sliceID>=nthinslices) continue;
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
          }
        }
      }
    }
  }
}

  
