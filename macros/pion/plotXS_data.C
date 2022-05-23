#include "../../SliceParams.h"

void plotXS_data(){

  gStyle->SetOptStat(0);

  TFile *file = TFile::Open("XS.root");
  
  TGraphErrors *gr_truexs = (TGraphErrors*)file->Get("gr_truexs");
  TGraphErrors *gr_recoxs = (TGraphErrors*)file->Get("gr_recoxs");
  
  TFile f2("../files/exclusive_xsec.root");
  TGraph *total_inel_KE = (TGraph*)f2.Get("total_inel_KE");
  
  /*for (int i: {0, 12, 12, 12, 12, 12, 12, 12}) {
    gr_truexs->RemovePoint(i);
    gr_recoxs->RemovePoint(i);
  }*/
  ofstream myfile;
  myfile.open ("toys/example.txt", ios::app);
  
  //Ninc
  for (int i=0; i<pi::nthinslices; ++i) {
    double ninc = gr_inc->GetPointY(i);
    myfile<<ninc<<"\t";
  }
  myfile<<endl;
  for (int i=0; i<pi::nthinslices; ++i) {
    double ninc_err = gr_inc->GetErrorY(i);
    myfile<<ninc_err<<"\t";
  }
  myfile<<endl;
  //Nint
  for (int i=0; i<pi::nthinslices; ++i) {
    double nint = gr_int->GetPointY(i);
    myfile<<nint<<"\t";
  }
  myfile<<endl;
  for (int i=0; i<pi::nthinslices; ++i) {
    double nint_err = gr_int->GetErrorY(i);
    myfile<<nint_err<<"\t";
  }
  myfile<<endl;
  /*double MCg4rw[20] = {
    1,1,1,1,1,1,1,1,1,1,
    1,1,1,1,1,1,1,1,1,1,
  };
  for (int i=0; i<pi::nthinslices; ++i) {
    double KE = gr_recoxs->GetPointX(i);
    double xs_true = total_inel_KE->Eval(KE)*MCg4rw[i];
    //bool cover = abs(xs - xs_true) < xs_err;
    cout<<xs_true<<"\t";
  }
  cout<<endl;*/
  //XS
  for (int i=0; i<pi::nthinslices; ++i) {
    double xs = gr_recoxs->GetPointY(i);
    myfile<<xs<<"\t";
  }
  myfile<<endl;
  for (int i=0; i<pi::nthinslices; ++i) {
    double xs_err = gr_recoxs->GetErrorY(i);
    myfile<<xs_err<<"\t";
  }
  myfile<<endl;
  
  myfile<<endl;
  myfile.close();
}
