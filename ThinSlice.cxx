#include "ThinSlice.h"
#include "TGraphErrors.h"
#include <iostream>

void ThinSlice::BookHistograms(){

   for (int i = 0; i<nthinslices; ++i){
     reco_incE[i] = new TH1D(Form("reco_incE_%d",i),Form("Reco incident energy, %.1f < z < %.1f (cm)",i*thinslicewidth, (i+1)*thinslicewidth), nbinse, 0, 1200.);
     true_incE[i] = new TH1D(Form("true_incE_%d",i),Form("True incident energy, %.1f < z < %.1f (cm)",i*thinslicewidth, (i+1)*thinslicewidth), nbinse, 0, 1200.);
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
      if (!(evt.reco_beam_calo_wire->empty()) && evt.reco_beam_true_byE_matched){
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

void ThinSlice::CalcXS(){

  double avg_trueincE[nthinslices] = {0};
  double truexs[nthinslices] = {0};
  double err_truexs[nthinslices] = {0};

  double NA=6.02214076e23;
  double MAr=39.95; //gmol
  double Density = 1.39; // g/cm^3

  for (int i = 0; i<nthinslices; ++i){
    
    avg_trueincE[i] = true_incE[i]->GetMean();
    //std::cout<<i<<" "<<avg_trueincE[i]<<std::endl;
    if (true_incidents[i] && true_interactions[i]){
      truexs[i] = MAr/(Density*NA*thinslicewidth)*log(true_incidents[i]/(true_incidents[i]-true_interactions[i]))*1e27;
      err_truexs[i] = MAr/(Density*NA*thinslicewidth)*1e27*sqrt(true_interactions[i]+pow(true_interactions[i],2)/true_incidents[i])/true_incidents[i];
    }
  }

  TGraphErrors *gr_truexs = new TGraphErrors(nthinslices, &(avg_trueincE[0]), &(truexs[0]), 0, &(err_truexs[0]));
  
  gr_truexs->Write("gr_truexs");
}
