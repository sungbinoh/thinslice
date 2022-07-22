#ifndef BeamNtuple_h
#define BeamNtuple_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

class BeamNtuple {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         AfterTarget_x;
   Float_t         AfterTarget_y;
   Float_t         AfterTarget_z;
   Float_t         AfterTarget_Px;
   Float_t         AfterTarget_Py;
   Float_t         AfterTarget_Pz;
   Float_t         AfterTarget_t;
   Float_t         AfterTarget_PDGid;
   Float_t         AfterTarget_EventID;
   Float_t         AfterTarget_TrackID;
   Float_t         AfterTarget_ParentID;
   Float_t         AfterTarget_Weight;
   Float_t         AfterTarget_Bx;
   Float_t         AfterTarget_By;
   Float_t         AfterTarget_Bz;
   Float_t         AfterTarget_Ex;
   Float_t         AfterTarget_Ey;
   Float_t         AfterTarget_Ez;
   Float_t         AfterTarget_ProperTime;
   Float_t         AfterTarget_PathLength;
   Float_t         AfterTarget_PolX;
   Float_t         AfterTarget_PolY;
   Float_t         AfterTarget_PolZ;
   Float_t         AfterTarget_InitX;
   Float_t         AfterTarget_InitY;
   Float_t         AfterTarget_InitZ;
   Float_t         AfterTarget_InitT;
   Float_t         AfterTarget_InitKE;
   Float_t         TOF1_x;
   Float_t         TOF1_y;
   Float_t         TOF1_z;
   Float_t         TOF1_Px;
   Float_t         TOF1_Py;
   Float_t         TOF1_Pz;
   Float_t         TOF1_t;
   Float_t         TOF1_PDGid;
   Float_t         TOF1_EventID;
   Float_t         TOF1_TrackID;
   Float_t         TOF1_ParentID;
   Float_t         TOF1_Weight;
   Float_t         TOF1_Bx;
   Float_t         TOF1_By;
   Float_t         TOF1_Bz;
   Float_t         TOF1_Ex;
   Float_t         TOF1_Ey;
   Float_t         TOF1_Ez;
   Float_t         TOF1_ProperTime;
   Float_t         TOF1_PathLength;
   Float_t         TOF1_PolX;
   Float_t         TOF1_PolY;
   Float_t         TOF1_PolZ;
   Float_t         TOF1_InitX;
   Float_t         TOF1_InitY;
   Float_t         TOF1_InitZ;
   Float_t         TOF1_InitT;
   Float_t         TOF1_InitKE;
   Float_t         COLL1_x;
   Float_t         COLL1_y;
   Float_t         COLL1_z;
   Float_t         COLL1_Px;
   Float_t         COLL1_Py;
   Float_t         COLL1_Pz;
   Float_t         COLL1_t;
   Float_t         COLL1_PDGid;
   Float_t         COLL1_EventID;
   Float_t         COLL1_TrackID;
   Float_t         COLL1_ParentID;
   Float_t         COLL1_Weight;
   Float_t         COLL1_Bx;
   Float_t         COLL1_By;
   Float_t         COLL1_Bz;
   Float_t         COLL1_Ex;
   Float_t         COLL1_Ey;
   Float_t         COLL1_Ez;
   Float_t         COLL1_ProperTime;
   Float_t         COLL1_PathLength;
   Float_t         COLL1_PolX;
   Float_t         COLL1_PolY;
   Float_t         COLL1_PolZ;
   Float_t         COLL1_InitX;
   Float_t         COLL1_InitY;
   Float_t         COLL1_InitZ;
   Float_t         COLL1_InitT;
   Float_t         COLL1_InitKE;
   Float_t         BPROF1_x;
   Float_t         BPROF1_y;
   Float_t         BPROF1_z;
   Float_t         BPROF1_Px;
   Float_t         BPROF1_Py;
   Float_t         BPROF1_Pz;
   Float_t         BPROF1_t;
   Float_t         BPROF1_PDGid;
   Float_t         BPROF1_EventID;
   Float_t         BPROF1_TrackID;
   Float_t         BPROF1_ParentID;
   Float_t         BPROF1_Weight;
   Float_t         BPROF1_Bx;
   Float_t         BPROF1_By;
   Float_t         BPROF1_Bz;
   Float_t         BPROF1_Ex;
   Float_t         BPROF1_Ey;
   Float_t         BPROF1_Ez;
   Float_t         BPROF1_ProperTime;
   Float_t         BPROF1_PathLength;
   Float_t         BPROF1_PolX;
   Float_t         BPROF1_PolY;
   Float_t         BPROF1_PolZ;
   Float_t         BPROF1_InitX;
   Float_t         BPROF1_InitY;
   Float_t         BPROF1_InitZ;
   Float_t         BPROF1_InitT;
   Float_t         BPROF1_InitKE;
   Float_t         BPROF2_x;
   Float_t         BPROF2_y;
   Float_t         BPROF2_z;
   Float_t         BPROF2_Px;
   Float_t         BPROF2_Py;
   Float_t         BPROF2_Pz;
   Float_t         BPROF2_t;
   Float_t         BPROF2_PDGid;
   Float_t         BPROF2_EventID;
   Float_t         BPROF2_TrackID;
   Float_t         BPROF2_ParentID;
   Float_t         BPROF2_Weight;
   Float_t         BPROF2_Bx;
   Float_t         BPROF2_By;
   Float_t         BPROF2_Bz;
   Float_t         BPROF2_Ex;
   Float_t         BPROF2_Ey;
   Float_t         BPROF2_Ez;
   Float_t         BPROF2_ProperTime;
   Float_t         BPROF2_PathLength;
   Float_t         BPROF2_PolX;
   Float_t         BPROF2_PolY;
   Float_t         BPROF2_PolZ;
   Float_t         BPROF2_InitX;
   Float_t         BPROF2_InitY;
   Float_t         BPROF2_InitZ;
   Float_t         BPROF2_InitT;
   Float_t         BPROF2_InitKE;
   Float_t         BPROF3_x;
   Float_t         BPROF3_y;
   Float_t         BPROF3_z;
   Float_t         BPROF3_Px;
   Float_t         BPROF3_Py;
   Float_t         BPROF3_Pz;
   Float_t         BPROF3_t;
   Float_t         BPROF3_PDGid;
   Float_t         BPROF3_EventID;
   Float_t         BPROF3_TrackID;
   Float_t         BPROF3_ParentID;
   Float_t         BPROF3_Weight;
   Float_t         BPROF3_Bx;
   Float_t         BPROF3_By;
   Float_t         BPROF3_Bz;
   Float_t         BPROF3_Ex;
   Float_t         BPROF3_Ey;
   Float_t         BPROF3_Ez;
   Float_t         BPROF3_ProperTime;
   Float_t         BPROF3_PathLength;
   Float_t         BPROF3_PolX;
   Float_t         BPROF3_PolY;
   Float_t         BPROF3_PolZ;
   Float_t         BPROF3_InitX;
   Float_t         BPROF3_InitY;
   Float_t         BPROF3_InitZ;
   Float_t         BPROF3_InitT;
   Float_t         BPROF3_InitKE;
   Float_t         TRIG1_x;
   Float_t         TRIG1_y;
   Float_t         TRIG1_z;
   Float_t         TRIG1_Px;
   Float_t         TRIG1_Py;
   Float_t         TRIG1_Pz;
   Float_t         TRIG1_t;
   Float_t         TRIG1_PDGid;
   Float_t         TRIG1_EventID;
   Float_t         TRIG1_TrackID;
   Float_t         TRIG1_ParentID;
   Float_t         TRIG1_Weight;
   Float_t         TRIG1_Bx;
   Float_t         TRIG1_By;
   Float_t         TRIG1_Bz;
   Float_t         TRIG1_Ex;
   Float_t         TRIG1_Ey;
   Float_t         TRIG1_Ez;
   Float_t         TRIG1_ProperTime;
   Float_t         TRIG1_PathLength;
   Float_t         TRIG1_PolX;
   Float_t         TRIG1_PolY;
   Float_t         TRIG1_PolZ;
   Float_t         TRIG1_InitX;
   Float_t         TRIG1_InitY;
   Float_t         TRIG1_InitZ;
   Float_t         TRIG1_InitT;
   Float_t         TRIG1_InitKE;
   Float_t         BPROFEXT_x;
   Float_t         BPROFEXT_y;
   Float_t         BPROFEXT_z;
   Float_t         BPROFEXT_Px;
   Float_t         BPROFEXT_Py;
   Float_t         BPROFEXT_Pz;
   Float_t         BPROFEXT_t;
   Float_t         BPROFEXT_PDGid;
   Float_t         BPROFEXT_EventID;
   Float_t         BPROFEXT_TrackID;
   Float_t         BPROFEXT_ParentID;
   Float_t         BPROFEXT_Weight;
   Float_t         BPROFEXT_Bx;
   Float_t         BPROFEXT_By;
   Float_t         BPROFEXT_Bz;
   Float_t         BPROFEXT_Ex;
   Float_t         BPROFEXT_Ey;
   Float_t         BPROFEXT_Ez;
   Float_t         BPROFEXT_ProperTime;
   Float_t         BPROFEXT_PathLength;
   Float_t         BPROFEXT_PolX;
   Float_t         BPROFEXT_PolY;
   Float_t         BPROFEXT_PolZ;
   Float_t         BPROFEXT_InitX;
   Float_t         BPROFEXT_InitY;
   Float_t         BPROFEXT_InitZ;
   Float_t         BPROFEXT_InitT;
   Float_t         BPROFEXT_InitKE;
   Float_t         BPROF4_x;
   Float_t         BPROF4_y;
   Float_t         BPROF4_z;
   Float_t         BPROF4_Px;
   Float_t         BPROF4_Py;
   Float_t         BPROF4_Pz;
   Float_t         BPROF4_t;
   Float_t         BPROF4_PDGid;
   Float_t         BPROF4_EventID;
   Float_t         BPROF4_TrackID;
   Float_t         BPROF4_ParentID;
   Float_t         BPROF4_Weight;
   Float_t         BPROF4_Bx;
   Float_t         BPROF4_By;
   Float_t         BPROF4_Bz;
   Float_t         BPROF4_Ex;
   Float_t         BPROF4_Ey;
   Float_t         BPROF4_Ez;
   Float_t         BPROF4_ProperTime;
   Float_t         BPROF4_PathLength;
   Float_t         BPROF4_PolX;
   Float_t         BPROF4_PolY;
   Float_t         BPROF4_PolZ;
   Float_t         BPROF4_InitX;
   Float_t         BPROF4_InitY;
   Float_t         BPROF4_InitZ;
   Float_t         BPROF4_InitT;
   Float_t         BPROF4_InitKE;
   Float_t         TRIG2_x;
   Float_t         TRIG2_y;
   Float_t         TRIG2_z;
   Float_t         TRIG2_Px;
   Float_t         TRIG2_Py;
   Float_t         TRIG2_Pz;
   Float_t         TRIG2_t;
   Float_t         TRIG2_PDGid;
   Float_t         TRIG2_EventID;
   Float_t         TRIG2_TrackID;
   Float_t         TRIG2_ParentID;
   Float_t         TRIG2_Weight;
   Float_t         TRIG2_Bx;
   Float_t         TRIG2_By;
   Float_t         TRIG2_Bz;
   Float_t         TRIG2_Ex;
   Float_t         TRIG2_Ey;
   Float_t         TRIG2_Ez;
   Float_t         TRIG2_ProperTime;
   Float_t         TRIG2_PathLength;
   Float_t         TRIG2_PolX;
   Float_t         TRIG2_PolY;
   Float_t         TRIG2_PolZ;
   Float_t         TRIG2_InitX;
   Float_t         TRIG2_InitY;
   Float_t         TRIG2_InitZ;
   Float_t         TRIG2_InitT;
   Float_t         TRIG2_InitKE;
   Float_t         NP04front_x;
   Float_t         NP04front_y;
   Float_t         NP04front_z;
   Float_t         NP04front_Px;
   Float_t         NP04front_Py;
   Float_t         NP04front_Pz;
   Float_t         NP04front_t;
   Float_t         NP04front_PDGid;
   Float_t         NP04front_EventID;
   Float_t         NP04front_TrackID;
   Float_t         NP04front_ParentID;
   Float_t         NP04front_Weight;
   Float_t         NP04front_Edep;
   Float_t         NP04front_VisibleEdep;
   Float_t         NP04front_Ntracks;
   Float_t         NP04FieldCage_x;
   Float_t         NP04FieldCage_y;
   Float_t         NP04FieldCage_z;
   Float_t         NP04FieldCage_Px;
   Float_t         NP04FieldCage_Py;
   Float_t         NP04FieldCage_Pz;
   Float_t         NP04FieldCage_t;
   Float_t         NP04FieldCage_PDGid;
   Float_t         NP04FieldCage_EventID;
   Float_t         NP04FieldCage_TrackID;
   Float_t         NP04FieldCage_ParentID;
   Float_t         NP04FieldCage_Weight;
   Float_t         NP04FieldCage_Edep;
   Float_t         NP04FieldCage_VisibleEdep;
   Float_t         NP04FieldCage_Ntracks;

   // List of branches
   TBranch        *b_AfterTarget_x;   //!
   TBranch        *b_AfterTarget_y;   //!
   TBranch        *b_AfterTarget_z;   //!
   TBranch        *b_AfterTarget_Px;   //!
   TBranch        *b_AfterTarget_Py;   //!
   TBranch        *b_AfterTarget_Pz;   //!
   TBranch        *b_AfterTarget_t;   //!
   TBranch        *b_AfterTarget_PDGid;   //!
   TBranch        *b_AfterTarget_EventID;   //!
   TBranch        *b_AfterTarget_TrackID;   //!
   TBranch        *b_AfterTarget_ParentID;   //!
   TBranch        *b_AfterTarget_Weight;   //!
   TBranch        *b_AfterTarget_Bx;   //!
   TBranch        *b_AfterTarget_By;   //!
   TBranch        *b_AfterTarget_Bz;   //!
   TBranch        *b_AfterTarget_Ex;   //!
   TBranch        *b_AfterTarget_Ey;   //!
   TBranch        *b_AfterTarget_Ez;   //!
   TBranch        *b_AfterTarget_ProperTime;   //!
   TBranch        *b_AfterTarget_PathLength;   //!
   TBranch        *b_AfterTarget_PolX;   //!
   TBranch        *b_AfterTarget_PolY;   //!
   TBranch        *b_AfterTarget_PolZ;   //!
   TBranch        *b_AfterTarget_InitX;   //!
   TBranch        *b_AfterTarget_InitY;   //!
   TBranch        *b_AfterTarget_InitZ;   //!
   TBranch        *b_AfterTarget_InitT;   //!
   TBranch        *b_AfterTarget_InitKE;   //!
   TBranch        *b_TOF1_x;   //!
   TBranch        *b_TOF1_y;   //!
   TBranch        *b_TOF1_z;   //!
   TBranch        *b_TOF1_Px;   //!
   TBranch        *b_TOF1_Py;   //!
   TBranch        *b_TOF1_Pz;   //!
   TBranch        *b_TOF1_t;   //!
   TBranch        *b_TOF1_PDGid;   //!
   TBranch        *b_TOF1_EventID;   //!
   TBranch        *b_TOF1_TrackID;   //!
   TBranch        *b_TOF1_ParentID;   //!
   TBranch        *b_TOF1_Weight;   //!
   TBranch        *b_TOF1_Bx;   //!
   TBranch        *b_TOF1_By;   //!
   TBranch        *b_TOF1_Bz;   //!
   TBranch        *b_TOF1_Ex;   //!
   TBranch        *b_TOF1_Ey;   //!
   TBranch        *b_TOF1_Ez;   //!
   TBranch        *b_TOF1_ProperTime;   //!
   TBranch        *b_TOF1_PathLength;   //!
   TBranch        *b_TOF1_PolX;   //!
   TBranch        *b_TOF1_PolY;   //!
   TBranch        *b_TOF1_PolZ;   //!
   TBranch        *b_TOF1_InitX;   //!
   TBranch        *b_TOF1_InitY;   //!
   TBranch        *b_TOF1_InitZ;   //!
   TBranch        *b_TOF1_InitT;   //!
   TBranch        *b_TOF1_InitKE;   //!
   TBranch        *b_COLL1_x;   //!
   TBranch        *b_COLL1_y;   //!
   TBranch        *b_COLL1_z;   //!
   TBranch        *b_COLL1_Px;   //!
   TBranch        *b_COLL1_Py;   //!
   TBranch        *b_COLL1_Pz;   //!
   TBranch        *b_COLL1_t;   //!
   TBranch        *b_COLL1_PDGid;   //!
   TBranch        *b_COLL1_EventID;   //!
   TBranch        *b_COLL1_TrackID;   //!
   TBranch        *b_COLL1_ParentID;   //!
   TBranch        *b_COLL1_Weight;   //!
   TBranch        *b_COLL1_Bx;   //!
   TBranch        *b_COLL1_By;   //!
   TBranch        *b_COLL1_Bz;   //!
   TBranch        *b_COLL1_Ex;   //!
   TBranch        *b_COLL1_Ey;   //!
   TBranch        *b_COLL1_Ez;   //!
   TBranch        *b_COLL1_ProperTime;   //!
   TBranch        *b_COLL1_PathLength;   //!
   TBranch        *b_COLL1_PolX;   //!
   TBranch        *b_COLL1_PolY;   //!
   TBranch        *b_COLL1_PolZ;   //!
   TBranch        *b_COLL1_InitX;   //!
   TBranch        *b_COLL1_InitY;   //!
   TBranch        *b_COLL1_InitZ;   //!
   TBranch        *b_COLL1_InitT;   //!
   TBranch        *b_COLL1_InitKE;   //!
   TBranch        *b_BPROF1_x;   //!
   TBranch        *b_BPROF1_y;   //!
   TBranch        *b_BPROF1_z;   //!
   TBranch        *b_BPROF1_Px;   //!
   TBranch        *b_BPROF1_Py;   //!
   TBranch        *b_BPROF1_Pz;   //!
   TBranch        *b_BPROF1_t;   //!
   TBranch        *b_BPROF1_PDGid;   //!
   TBranch        *b_BPROF1_EventID;   //!
   TBranch        *b_BPROF1_TrackID;   //!
   TBranch        *b_BPROF1_ParentID;   //!
   TBranch        *b_BPROF1_Weight;   //!
   TBranch        *b_BPROF1_Bx;   //!
   TBranch        *b_BPROF1_By;   //!
   TBranch        *b_BPROF1_Bz;   //!
   TBranch        *b_BPROF1_Ex;   //!
   TBranch        *b_BPROF1_Ey;   //!
   TBranch        *b_BPROF1_Ez;   //!
   TBranch        *b_BPROF1_ProperTime;   //!
   TBranch        *b_BPROF1_PathLength;   //!
   TBranch        *b_BPROF1_PolX;   //!
   TBranch        *b_BPROF1_PolY;   //!
   TBranch        *b_BPROF1_PolZ;   //!
   TBranch        *b_BPROF1_InitX;   //!
   TBranch        *b_BPROF1_InitY;   //!
   TBranch        *b_BPROF1_InitZ;   //!
   TBranch        *b_BPROF1_InitT;   //!
   TBranch        *b_BPROF1_InitKE;   //!
   TBranch        *b_BPROF2_x;   //!
   TBranch        *b_BPROF2_y;   //!
   TBranch        *b_BPROF2_z;   //!
   TBranch        *b_BPROF2_Px;   //!
   TBranch        *b_BPROF2_Py;   //!
   TBranch        *b_BPROF2_Pz;   //!
   TBranch        *b_BPROF2_t;   //!
   TBranch        *b_BPROF2_PDGid;   //!
   TBranch        *b_BPROF2_EventID;   //!
   TBranch        *b_BPROF2_TrackID;   //!
   TBranch        *b_BPROF2_ParentID;   //!
   TBranch        *b_BPROF2_Weight;   //!
   TBranch        *b_BPROF2_Bx;   //!
   TBranch        *b_BPROF2_By;   //!
   TBranch        *b_BPROF2_Bz;   //!
   TBranch        *b_BPROF2_Ex;   //!
   TBranch        *b_BPROF2_Ey;   //!
   TBranch        *b_BPROF2_Ez;   //!
   TBranch        *b_BPROF2_ProperTime;   //!
   TBranch        *b_BPROF2_PathLength;   //!
   TBranch        *b_BPROF2_PolX;   //!
   TBranch        *b_BPROF2_PolY;   //!
   TBranch        *b_BPROF2_PolZ;   //!
   TBranch        *b_BPROF2_InitX;   //!
   TBranch        *b_BPROF2_InitY;   //!
   TBranch        *b_BPROF2_InitZ;   //!
   TBranch        *b_BPROF2_InitT;   //!
   TBranch        *b_BPROF2_InitKE;   //!
   TBranch        *b_BPROF3_x;   //!
   TBranch        *b_BPROF3_y;   //!
   TBranch        *b_BPROF3_z;   //!
   TBranch        *b_BPROF3_Px;   //!
   TBranch        *b_BPROF3_Py;   //!
   TBranch        *b_BPROF3_Pz;   //!
   TBranch        *b_BPROF3_t;   //!
   TBranch        *b_BPROF3_PDGid;   //!
   TBranch        *b_BPROF3_EventID;   //!
   TBranch        *b_BPROF3_TrackID;   //!
   TBranch        *b_BPROF3_ParentID;   //!
   TBranch        *b_BPROF3_Weight;   //!
   TBranch        *b_BPROF3_Bx;   //!
   TBranch        *b_BPROF3_By;   //!
   TBranch        *b_BPROF3_Bz;   //!
   TBranch        *b_BPROF3_Ex;   //!
   TBranch        *b_BPROF3_Ey;   //!
   TBranch        *b_BPROF3_Ez;   //!
   TBranch        *b_BPROF3_ProperTime;   //!
   TBranch        *b_BPROF3_PathLength;   //!
   TBranch        *b_BPROF3_PolX;   //!
   TBranch        *b_BPROF3_PolY;   //!
   TBranch        *b_BPROF3_PolZ;   //!
   TBranch        *b_BPROF3_InitX;   //!
   TBranch        *b_BPROF3_InitY;   //!
   TBranch        *b_BPROF3_InitZ;   //!
   TBranch        *b_BPROF3_InitT;   //!
   TBranch        *b_BPROF3_InitKE;   //!
   TBranch        *b_TRIG1_x;   //!
   TBranch        *b_TRIG1_y;   //!
   TBranch        *b_TRIG1_z;   //!
   TBranch        *b_TRIG1_Px;   //!
   TBranch        *b_TRIG1_Py;   //!
   TBranch        *b_TRIG1_Pz;   //!
   TBranch        *b_TRIG1_t;   //!
   TBranch        *b_TRIG1_PDGid;   //!
   TBranch        *b_TRIG1_EventID;   //!
   TBranch        *b_TRIG1_TrackID;   //!
   TBranch        *b_TRIG1_ParentID;   //!
   TBranch        *b_TRIG1_Weight;   //!
   TBranch        *b_TRIG1_Bx;   //!
   TBranch        *b_TRIG1_By;   //!
   TBranch        *b_TRIG1_Bz;   //!
   TBranch        *b_TRIG1_Ex;   //!
   TBranch        *b_TRIG1_Ey;   //!
   TBranch        *b_TRIG1_Ez;   //!
   TBranch        *b_TRIG1_ProperTime;   //!
   TBranch        *b_TRIG1_PathLength;   //!
   TBranch        *b_TRIG1_PolX;   //!
   TBranch        *b_TRIG1_PolY;   //!
   TBranch        *b_TRIG1_PolZ;   //!
   TBranch        *b_TRIG1_InitX;   //!
   TBranch        *b_TRIG1_InitY;   //!
   TBranch        *b_TRIG1_InitZ;   //!
   TBranch        *b_TRIG1_InitT;   //!
   TBranch        *b_TRIG1_InitKE;   //!
   TBranch        *b_BPROFEXT_x;   //!
   TBranch        *b_BPROFEXT_y;   //!
   TBranch        *b_BPROFEXT_z;   //!
   TBranch        *b_BPROFEXT_Px;   //!
   TBranch        *b_BPROFEXT_Py;   //!
   TBranch        *b_BPROFEXT_Pz;   //!
   TBranch        *b_BPROFEXT_t;   //!
   TBranch        *b_BPROFEXT_PDGid;   //!
   TBranch        *b_BPROFEXT_EventID;   //!
   TBranch        *b_BPROFEXT_TrackID;   //!
   TBranch        *b_BPROFEXT_ParentID;   //!
   TBranch        *b_BPROFEXT_Weight;   //!
   TBranch        *b_BPROFEXT_Bx;   //!
   TBranch        *b_BPROFEXT_By;   //!
   TBranch        *b_BPROFEXT_Bz;   //!
   TBranch        *b_BPROFEXT_Ex;   //!
   TBranch        *b_BPROFEXT_Ey;   //!
   TBranch        *b_BPROFEXT_Ez;   //!
   TBranch        *b_BPROFEXT_ProperTime;   //!
   TBranch        *b_BPROFEXT_PathLength;   //!
   TBranch        *b_BPROFEXT_PolX;   //!
   TBranch        *b_BPROFEXT_PolY;   //!
   TBranch        *b_BPROFEXT_PolZ;   //!
   TBranch        *b_BPROFEXT_InitX;   //!
   TBranch        *b_BPROFEXT_InitY;   //!
   TBranch        *b_BPROFEXT_InitZ;   //!
   TBranch        *b_BPROFEXT_InitT;   //!
   TBranch        *b_BPROFEXT_InitKE;   //!
   TBranch        *b_BPROF4_x;   //!
   TBranch        *b_BPROF4_y;   //!
   TBranch        *b_BPROF4_z;   //!
   TBranch        *b_BPROF4_Px;   //!
   TBranch        *b_BPROF4_Py;   //!
   TBranch        *b_BPROF4_Pz;   //!
   TBranch        *b_BPROF4_t;   //!
   TBranch        *b_BPROF4_PDGid;   //!
   TBranch        *b_BPROF4_EventID;   //!
   TBranch        *b_BPROF4_TrackID;   //!
   TBranch        *b_BPROF4_ParentID;   //!
   TBranch        *b_BPROF4_Weight;   //!
   TBranch        *b_BPROF4_Bx;   //!
   TBranch        *b_BPROF4_By;   //!
   TBranch        *b_BPROF4_Bz;   //!
   TBranch        *b_BPROF4_Ex;   //!
   TBranch        *b_BPROF4_Ey;   //!
   TBranch        *b_BPROF4_Ez;   //!
   TBranch        *b_BPROF4_ProperTime;   //!
   TBranch        *b_BPROF4_PathLength;   //!
   TBranch        *b_BPROF4_PolX;   //!
   TBranch        *b_BPROF4_PolY;   //!
   TBranch        *b_BPROF4_PolZ;   //!
   TBranch        *b_BPROF4_InitX;   //!
   TBranch        *b_BPROF4_InitY;   //!
   TBranch        *b_BPROF4_InitZ;   //!
   TBranch        *b_BPROF4_InitT;   //!
   TBranch        *b_BPROF4_InitKE;   //!
   TBranch        *b_TRIG2_x;   //!
   TBranch        *b_TRIG2_y;   //!
   TBranch        *b_TRIG2_z;   //!
   TBranch        *b_TRIG2_Px;   //!
   TBranch        *b_TRIG2_Py;   //!
   TBranch        *b_TRIG2_Pz;   //!
   TBranch        *b_TRIG2_t;   //!
   TBranch        *b_TRIG2_PDGid;   //!
   TBranch        *b_TRIG2_EventID;   //!
   TBranch        *b_TRIG2_TrackID;   //!
   TBranch        *b_TRIG2_ParentID;   //!
   TBranch        *b_TRIG2_Weight;   //!
   TBranch        *b_TRIG2_Bx;   //!
   TBranch        *b_TRIG2_By;   //!
   TBranch        *b_TRIG2_Bz;   //!
   TBranch        *b_TRIG2_Ex;   //!
   TBranch        *b_TRIG2_Ey;   //!
   TBranch        *b_TRIG2_Ez;   //!
   TBranch        *b_TRIG2_ProperTime;   //!
   TBranch        *b_TRIG2_PathLength;   //!
   TBranch        *b_TRIG2_PolX;   //!
   TBranch        *b_TRIG2_PolY;   //!
   TBranch        *b_TRIG2_PolZ;   //!
   TBranch        *b_TRIG2_InitX;   //!
   TBranch        *b_TRIG2_InitY;   //!
   TBranch        *b_TRIG2_InitZ;   //!
   TBranch        *b_TRIG2_InitT;   //!
   TBranch        *b_TRIG2_InitKE;   //!
   TBranch        *b_NP04front_x;   //!
   TBranch        *b_NP04front_y;   //!
   TBranch        *b_NP04front_z;   //!
   TBranch        *b_NP04front_Px;   //!
   TBranch        *b_NP04front_Py;   //!
   TBranch        *b_NP04front_Pz;   //!
   TBranch        *b_NP04front_t;   //!
   TBranch        *b_NP04front_PDGid;   //!
   TBranch        *b_NP04front_EventID;   //!
   TBranch        *b_NP04front_TrackID;   //!
   TBranch        *b_NP04front_ParentID;   //!
   TBranch        *b_NP04front_Weight;   //!
   TBranch        *b_NP04front_Edep;   //!
   TBranch        *b_NP04front_VisibleEdep;   //!
   TBranch        *b_NP04front_Ntracks;   //!
   TBranch        *b_NP04FieldCage_x;   //!
   TBranch        *b_NP04FieldCage_y;   //!
   TBranch        *b_NP04FieldCage_z;   //!
   TBranch        *b_NP04FieldCage_Px;   //!
   TBranch        *b_NP04FieldCage_Py;   //!
   TBranch        *b_NP04FieldCage_Pz;   //!
   TBranch        *b_NP04FieldCage_t;   //!
   TBranch        *b_NP04FieldCage_PDGid;   //!
   TBranch        *b_NP04FieldCage_EventID;   //!
   TBranch        *b_NP04FieldCage_TrackID;   //!
   TBranch        *b_NP04FieldCage_ParentID;   //!
   TBranch        *b_NP04FieldCage_Weight;   //!
   TBranch        *b_NP04FieldCage_Edep;   //!
   TBranch        *b_NP04FieldCage_VisibleEdep;   //!
   TBranch        *b_NP04FieldCage_Ntracks;   //!

   BeamNtuple(TTree *tree=0);
   virtual ~BeamNtuple();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef BeamNtuple_cxx
BeamNtuple::BeamNtuple(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("H4_v34b_1GeV_-27.7_10M_30.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("H4_v34b_1GeV_-27.7_10M_30.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("H4_v34b_1GeV_-27.7_10M_30.root:/NTuples");
      dir->GetObject("GoodParticle",tree);

   }
   Init(tree);
}

BeamNtuple::~BeamNtuple()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t BeamNtuple::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t BeamNtuple::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void BeamNtuple::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("AfterTarget_x", &AfterTarget_x, &b_AfterTarget_x);
   fChain->SetBranchAddress("AfterTarget_y", &AfterTarget_y, &b_AfterTarget_y);
   fChain->SetBranchAddress("AfterTarget_z", &AfterTarget_z, &b_AfterTarget_z);
   fChain->SetBranchAddress("AfterTarget_Px", &AfterTarget_Px, &b_AfterTarget_Px);
   fChain->SetBranchAddress("AfterTarget_Py", &AfterTarget_Py, &b_AfterTarget_Py);
   fChain->SetBranchAddress("AfterTarget_Pz", &AfterTarget_Pz, &b_AfterTarget_Pz);
   fChain->SetBranchAddress("AfterTarget_t", &AfterTarget_t, &b_AfterTarget_t);
   fChain->SetBranchAddress("AfterTarget_PDGid", &AfterTarget_PDGid, &b_AfterTarget_PDGid);
   fChain->SetBranchAddress("AfterTarget_EventID", &AfterTarget_EventID, &b_AfterTarget_EventID);
   fChain->SetBranchAddress("AfterTarget_TrackID", &AfterTarget_TrackID, &b_AfterTarget_TrackID);
   fChain->SetBranchAddress("AfterTarget_ParentID", &AfterTarget_ParentID, &b_AfterTarget_ParentID);
   fChain->SetBranchAddress("AfterTarget_Weight", &AfterTarget_Weight, &b_AfterTarget_Weight);
   fChain->SetBranchAddress("AfterTarget_Bx", &AfterTarget_Bx, &b_AfterTarget_Bx);
   fChain->SetBranchAddress("AfterTarget_By", &AfterTarget_By, &b_AfterTarget_By);
   fChain->SetBranchAddress("AfterTarget_Bz", &AfterTarget_Bz, &b_AfterTarget_Bz);
   fChain->SetBranchAddress("AfterTarget_Ex", &AfterTarget_Ex, &b_AfterTarget_Ex);
   fChain->SetBranchAddress("AfterTarget_Ey", &AfterTarget_Ey, &b_AfterTarget_Ey);
   fChain->SetBranchAddress("AfterTarget_Ez", &AfterTarget_Ez, &b_AfterTarget_Ez);
   fChain->SetBranchAddress("AfterTarget_ProperTime", &AfterTarget_ProperTime, &b_AfterTarget_ProperTime);
   fChain->SetBranchAddress("AfterTarget_PathLength", &AfterTarget_PathLength, &b_AfterTarget_PathLength);
   fChain->SetBranchAddress("AfterTarget_PolX", &AfterTarget_PolX, &b_AfterTarget_PolX);
   fChain->SetBranchAddress("AfterTarget_PolY", &AfterTarget_PolY, &b_AfterTarget_PolY);
   fChain->SetBranchAddress("AfterTarget_PolZ", &AfterTarget_PolZ, &b_AfterTarget_PolZ);
   fChain->SetBranchAddress("AfterTarget_InitX", &AfterTarget_InitX, &b_AfterTarget_InitX);
   fChain->SetBranchAddress("AfterTarget_InitY", &AfterTarget_InitY, &b_AfterTarget_InitY);
   fChain->SetBranchAddress("AfterTarget_InitZ", &AfterTarget_InitZ, &b_AfterTarget_InitZ);
   fChain->SetBranchAddress("AfterTarget_InitT", &AfterTarget_InitT, &b_AfterTarget_InitT);
   fChain->SetBranchAddress("AfterTarget_InitKE", &AfterTarget_InitKE, &b_AfterTarget_InitKE);
   fChain->SetBranchAddress("TOF1_x", &TOF1_x, &b_TOF1_x);
   fChain->SetBranchAddress("TOF1_y", &TOF1_y, &b_TOF1_y);
   fChain->SetBranchAddress("TOF1_z", &TOF1_z, &b_TOF1_z);
   fChain->SetBranchAddress("TOF1_Px", &TOF1_Px, &b_TOF1_Px);
   fChain->SetBranchAddress("TOF1_Py", &TOF1_Py, &b_TOF1_Py);
   fChain->SetBranchAddress("TOF1_Pz", &TOF1_Pz, &b_TOF1_Pz);
   fChain->SetBranchAddress("TOF1_t", &TOF1_t, &b_TOF1_t);
   fChain->SetBranchAddress("TOF1_PDGid", &TOF1_PDGid, &b_TOF1_PDGid);
   fChain->SetBranchAddress("TOF1_EventID", &TOF1_EventID, &b_TOF1_EventID);
   fChain->SetBranchAddress("TOF1_TrackID", &TOF1_TrackID, &b_TOF1_TrackID);
   fChain->SetBranchAddress("TOF1_ParentID", &TOF1_ParentID, &b_TOF1_ParentID);
   fChain->SetBranchAddress("TOF1_Weight", &TOF1_Weight, &b_TOF1_Weight);
   fChain->SetBranchAddress("TOF1_Bx", &TOF1_Bx, &b_TOF1_Bx);
   fChain->SetBranchAddress("TOF1_By", &TOF1_By, &b_TOF1_By);
   fChain->SetBranchAddress("TOF1_Bz", &TOF1_Bz, &b_TOF1_Bz);
   fChain->SetBranchAddress("TOF1_Ex", &TOF1_Ex, &b_TOF1_Ex);
   fChain->SetBranchAddress("TOF1_Ey", &TOF1_Ey, &b_TOF1_Ey);
   fChain->SetBranchAddress("TOF1_Ez", &TOF1_Ez, &b_TOF1_Ez);
   fChain->SetBranchAddress("TOF1_ProperTime", &TOF1_ProperTime, &b_TOF1_ProperTime);
   fChain->SetBranchAddress("TOF1_PathLength", &TOF1_PathLength, &b_TOF1_PathLength);
   fChain->SetBranchAddress("TOF1_PolX", &TOF1_PolX, &b_TOF1_PolX);
   fChain->SetBranchAddress("TOF1_PolY", &TOF1_PolY, &b_TOF1_PolY);
   fChain->SetBranchAddress("TOF1_PolZ", &TOF1_PolZ, &b_TOF1_PolZ);
   fChain->SetBranchAddress("TOF1_InitX", &TOF1_InitX, &b_TOF1_InitX);
   fChain->SetBranchAddress("TOF1_InitY", &TOF1_InitY, &b_TOF1_InitY);
   fChain->SetBranchAddress("TOF1_InitZ", &TOF1_InitZ, &b_TOF1_InitZ);
   fChain->SetBranchAddress("TOF1_InitT", &TOF1_InitT, &b_TOF1_InitT);
   fChain->SetBranchAddress("TOF1_InitKE", &TOF1_InitKE, &b_TOF1_InitKE);
   fChain->SetBranchAddress("COLL1_x", &COLL1_x, &b_COLL1_x);
   fChain->SetBranchAddress("COLL1_y", &COLL1_y, &b_COLL1_y);
   fChain->SetBranchAddress("COLL1_z", &COLL1_z, &b_COLL1_z);
   fChain->SetBranchAddress("COLL1_Px", &COLL1_Px, &b_COLL1_Px);
   fChain->SetBranchAddress("COLL1_Py", &COLL1_Py, &b_COLL1_Py);
   fChain->SetBranchAddress("COLL1_Pz", &COLL1_Pz, &b_COLL1_Pz);
   fChain->SetBranchAddress("COLL1_t", &COLL1_t, &b_COLL1_t);
   fChain->SetBranchAddress("COLL1_PDGid", &COLL1_PDGid, &b_COLL1_PDGid);
   fChain->SetBranchAddress("COLL1_EventID", &COLL1_EventID, &b_COLL1_EventID);
   fChain->SetBranchAddress("COLL1_TrackID", &COLL1_TrackID, &b_COLL1_TrackID);
   fChain->SetBranchAddress("COLL1_ParentID", &COLL1_ParentID, &b_COLL1_ParentID);
   fChain->SetBranchAddress("COLL1_Weight", &COLL1_Weight, &b_COLL1_Weight);
   fChain->SetBranchAddress("COLL1_Bx", &COLL1_Bx, &b_COLL1_Bx);
   fChain->SetBranchAddress("COLL1_By", &COLL1_By, &b_COLL1_By);
   fChain->SetBranchAddress("COLL1_Bz", &COLL1_Bz, &b_COLL1_Bz);
   fChain->SetBranchAddress("COLL1_Ex", &COLL1_Ex, &b_COLL1_Ex);
   fChain->SetBranchAddress("COLL1_Ey", &COLL1_Ey, &b_COLL1_Ey);
   fChain->SetBranchAddress("COLL1_Ez", &COLL1_Ez, &b_COLL1_Ez);
   fChain->SetBranchAddress("COLL1_ProperTime", &COLL1_ProperTime, &b_COLL1_ProperTime);
   fChain->SetBranchAddress("COLL1_PathLength", &COLL1_PathLength, &b_COLL1_PathLength);
   fChain->SetBranchAddress("COLL1_PolX", &COLL1_PolX, &b_COLL1_PolX);
   fChain->SetBranchAddress("COLL1_PolY", &COLL1_PolY, &b_COLL1_PolY);
   fChain->SetBranchAddress("COLL1_PolZ", &COLL1_PolZ, &b_COLL1_PolZ);
   fChain->SetBranchAddress("COLL1_InitX", &COLL1_InitX, &b_COLL1_InitX);
   fChain->SetBranchAddress("COLL1_InitY", &COLL1_InitY, &b_COLL1_InitY);
   fChain->SetBranchAddress("COLL1_InitZ", &COLL1_InitZ, &b_COLL1_InitZ);
   fChain->SetBranchAddress("COLL1_InitT", &COLL1_InitT, &b_COLL1_InitT);
   fChain->SetBranchAddress("COLL1_InitKE", &COLL1_InitKE, &b_COLL1_InitKE);
   fChain->SetBranchAddress("BPROF1_x", &BPROF1_x, &b_BPROF1_x);
   fChain->SetBranchAddress("BPROF1_y", &BPROF1_y, &b_BPROF1_y);
   fChain->SetBranchAddress("BPROF1_z", &BPROF1_z, &b_BPROF1_z);
   fChain->SetBranchAddress("BPROF1_Px", &BPROF1_Px, &b_BPROF1_Px);
   fChain->SetBranchAddress("BPROF1_Py", &BPROF1_Py, &b_BPROF1_Py);
   fChain->SetBranchAddress("BPROF1_Pz", &BPROF1_Pz, &b_BPROF1_Pz);
   fChain->SetBranchAddress("BPROF1_t", &BPROF1_t, &b_BPROF1_t);
   fChain->SetBranchAddress("BPROF1_PDGid", &BPROF1_PDGid, &b_BPROF1_PDGid);
   fChain->SetBranchAddress("BPROF1_EventID", &BPROF1_EventID, &b_BPROF1_EventID);
   fChain->SetBranchAddress("BPROF1_TrackID", &BPROF1_TrackID, &b_BPROF1_TrackID);
   fChain->SetBranchAddress("BPROF1_ParentID", &BPROF1_ParentID, &b_BPROF1_ParentID);
   fChain->SetBranchAddress("BPROF1_Weight", &BPROF1_Weight, &b_BPROF1_Weight);
   fChain->SetBranchAddress("BPROF1_Bx", &BPROF1_Bx, &b_BPROF1_Bx);
   fChain->SetBranchAddress("BPROF1_By", &BPROF1_By, &b_BPROF1_By);
   fChain->SetBranchAddress("BPROF1_Bz", &BPROF1_Bz, &b_BPROF1_Bz);
   fChain->SetBranchAddress("BPROF1_Ex", &BPROF1_Ex, &b_BPROF1_Ex);
   fChain->SetBranchAddress("BPROF1_Ey", &BPROF1_Ey, &b_BPROF1_Ey);
   fChain->SetBranchAddress("BPROF1_Ez", &BPROF1_Ez, &b_BPROF1_Ez);
   fChain->SetBranchAddress("BPROF1_ProperTime", &BPROF1_ProperTime, &b_BPROF1_ProperTime);
   fChain->SetBranchAddress("BPROF1_PathLength", &BPROF1_PathLength, &b_BPROF1_PathLength);
   fChain->SetBranchAddress("BPROF1_PolX", &BPROF1_PolX, &b_BPROF1_PolX);
   fChain->SetBranchAddress("BPROF1_PolY", &BPROF1_PolY, &b_BPROF1_PolY);
   fChain->SetBranchAddress("BPROF1_PolZ", &BPROF1_PolZ, &b_BPROF1_PolZ);
   fChain->SetBranchAddress("BPROF1_InitX", &BPROF1_InitX, &b_BPROF1_InitX);
   fChain->SetBranchAddress("BPROF1_InitY", &BPROF1_InitY, &b_BPROF1_InitY);
   fChain->SetBranchAddress("BPROF1_InitZ", &BPROF1_InitZ, &b_BPROF1_InitZ);
   fChain->SetBranchAddress("BPROF1_InitT", &BPROF1_InitT, &b_BPROF1_InitT);
   fChain->SetBranchAddress("BPROF1_InitKE", &BPROF1_InitKE, &b_BPROF1_InitKE);
   fChain->SetBranchAddress("BPROF2_x", &BPROF2_x, &b_BPROF2_x);
   fChain->SetBranchAddress("BPROF2_y", &BPROF2_y, &b_BPROF2_y);
   fChain->SetBranchAddress("BPROF2_z", &BPROF2_z, &b_BPROF2_z);
   fChain->SetBranchAddress("BPROF2_Px", &BPROF2_Px, &b_BPROF2_Px);
   fChain->SetBranchAddress("BPROF2_Py", &BPROF2_Py, &b_BPROF2_Py);
   fChain->SetBranchAddress("BPROF2_Pz", &BPROF2_Pz, &b_BPROF2_Pz);
   fChain->SetBranchAddress("BPROF2_t", &BPROF2_t, &b_BPROF2_t);
   fChain->SetBranchAddress("BPROF2_PDGid", &BPROF2_PDGid, &b_BPROF2_PDGid);
   fChain->SetBranchAddress("BPROF2_EventID", &BPROF2_EventID, &b_BPROF2_EventID);
   fChain->SetBranchAddress("BPROF2_TrackID", &BPROF2_TrackID, &b_BPROF2_TrackID);
   fChain->SetBranchAddress("BPROF2_ParentID", &BPROF2_ParentID, &b_BPROF2_ParentID);
   fChain->SetBranchAddress("BPROF2_Weight", &BPROF2_Weight, &b_BPROF2_Weight);
   fChain->SetBranchAddress("BPROF2_Bx", &BPROF2_Bx, &b_BPROF2_Bx);
   fChain->SetBranchAddress("BPROF2_By", &BPROF2_By, &b_BPROF2_By);
   fChain->SetBranchAddress("BPROF2_Bz", &BPROF2_Bz, &b_BPROF2_Bz);
   fChain->SetBranchAddress("BPROF2_Ex", &BPROF2_Ex, &b_BPROF2_Ex);
   fChain->SetBranchAddress("BPROF2_Ey", &BPROF2_Ey, &b_BPROF2_Ey);
   fChain->SetBranchAddress("BPROF2_Ez", &BPROF2_Ez, &b_BPROF2_Ez);
   fChain->SetBranchAddress("BPROF2_ProperTime", &BPROF2_ProperTime, &b_BPROF2_ProperTime);
   fChain->SetBranchAddress("BPROF2_PathLength", &BPROF2_PathLength, &b_BPROF2_PathLength);
   fChain->SetBranchAddress("BPROF2_PolX", &BPROF2_PolX, &b_BPROF2_PolX);
   fChain->SetBranchAddress("BPROF2_PolY", &BPROF2_PolY, &b_BPROF2_PolY);
   fChain->SetBranchAddress("BPROF2_PolZ", &BPROF2_PolZ, &b_BPROF2_PolZ);
   fChain->SetBranchAddress("BPROF2_InitX", &BPROF2_InitX, &b_BPROF2_InitX);
   fChain->SetBranchAddress("BPROF2_InitY", &BPROF2_InitY, &b_BPROF2_InitY);
   fChain->SetBranchAddress("BPROF2_InitZ", &BPROF2_InitZ, &b_BPROF2_InitZ);
   fChain->SetBranchAddress("BPROF2_InitT", &BPROF2_InitT, &b_BPROF2_InitT);
   fChain->SetBranchAddress("BPROF2_InitKE", &BPROF2_InitKE, &b_BPROF2_InitKE);
   fChain->SetBranchAddress("BPROF3_x", &BPROF3_x, &b_BPROF3_x);
   fChain->SetBranchAddress("BPROF3_y", &BPROF3_y, &b_BPROF3_y);
   fChain->SetBranchAddress("BPROF3_z", &BPROF3_z, &b_BPROF3_z);
   fChain->SetBranchAddress("BPROF3_Px", &BPROF3_Px, &b_BPROF3_Px);
   fChain->SetBranchAddress("BPROF3_Py", &BPROF3_Py, &b_BPROF3_Py);
   fChain->SetBranchAddress("BPROF3_Pz", &BPROF3_Pz, &b_BPROF3_Pz);
   fChain->SetBranchAddress("BPROF3_t", &BPROF3_t, &b_BPROF3_t);
   fChain->SetBranchAddress("BPROF3_PDGid", &BPROF3_PDGid, &b_BPROF3_PDGid);
   fChain->SetBranchAddress("BPROF3_EventID", &BPROF3_EventID, &b_BPROF3_EventID);
   fChain->SetBranchAddress("BPROF3_TrackID", &BPROF3_TrackID, &b_BPROF3_TrackID);
   fChain->SetBranchAddress("BPROF3_ParentID", &BPROF3_ParentID, &b_BPROF3_ParentID);
   fChain->SetBranchAddress("BPROF3_Weight", &BPROF3_Weight, &b_BPROF3_Weight);
   fChain->SetBranchAddress("BPROF3_Bx", &BPROF3_Bx, &b_BPROF3_Bx);
   fChain->SetBranchAddress("BPROF3_By", &BPROF3_By, &b_BPROF3_By);
   fChain->SetBranchAddress("BPROF3_Bz", &BPROF3_Bz, &b_BPROF3_Bz);
   fChain->SetBranchAddress("BPROF3_Ex", &BPROF3_Ex, &b_BPROF3_Ex);
   fChain->SetBranchAddress("BPROF3_Ey", &BPROF3_Ey, &b_BPROF3_Ey);
   fChain->SetBranchAddress("BPROF3_Ez", &BPROF3_Ez, &b_BPROF3_Ez);
   fChain->SetBranchAddress("BPROF3_ProperTime", &BPROF3_ProperTime, &b_BPROF3_ProperTime);
   fChain->SetBranchAddress("BPROF3_PathLength", &BPROF3_PathLength, &b_BPROF3_PathLength);
   fChain->SetBranchAddress("BPROF3_PolX", &BPROF3_PolX, &b_BPROF3_PolX);
   fChain->SetBranchAddress("BPROF3_PolY", &BPROF3_PolY, &b_BPROF3_PolY);
   fChain->SetBranchAddress("BPROF3_PolZ", &BPROF3_PolZ, &b_BPROF3_PolZ);
   fChain->SetBranchAddress("BPROF3_InitX", &BPROF3_InitX, &b_BPROF3_InitX);
   fChain->SetBranchAddress("BPROF3_InitY", &BPROF3_InitY, &b_BPROF3_InitY);
   fChain->SetBranchAddress("BPROF3_InitZ", &BPROF3_InitZ, &b_BPROF3_InitZ);
   fChain->SetBranchAddress("BPROF3_InitT", &BPROF3_InitT, &b_BPROF3_InitT);
   fChain->SetBranchAddress("BPROF3_InitKE", &BPROF3_InitKE, &b_BPROF3_InitKE);
   fChain->SetBranchAddress("TRIG1_x", &TRIG1_x, &b_TRIG1_x);
   fChain->SetBranchAddress("TRIG1_y", &TRIG1_y, &b_TRIG1_y);
   fChain->SetBranchAddress("TRIG1_z", &TRIG1_z, &b_TRIG1_z);
   fChain->SetBranchAddress("TRIG1_Px", &TRIG1_Px, &b_TRIG1_Px);
   fChain->SetBranchAddress("TRIG1_Py", &TRIG1_Py, &b_TRIG1_Py);
   fChain->SetBranchAddress("TRIG1_Pz", &TRIG1_Pz, &b_TRIG1_Pz);
   fChain->SetBranchAddress("TRIG1_t", &TRIG1_t, &b_TRIG1_t);
   fChain->SetBranchAddress("TRIG1_PDGid", &TRIG1_PDGid, &b_TRIG1_PDGid);
   fChain->SetBranchAddress("TRIG1_EventID", &TRIG1_EventID, &b_TRIG1_EventID);
   fChain->SetBranchAddress("TRIG1_TrackID", &TRIG1_TrackID, &b_TRIG1_TrackID);
   fChain->SetBranchAddress("TRIG1_ParentID", &TRIG1_ParentID, &b_TRIG1_ParentID);
   fChain->SetBranchAddress("TRIG1_Weight", &TRIG1_Weight, &b_TRIG1_Weight);
   fChain->SetBranchAddress("TRIG1_Bx", &TRIG1_Bx, &b_TRIG1_Bx);
   fChain->SetBranchAddress("TRIG1_By", &TRIG1_By, &b_TRIG1_By);
   fChain->SetBranchAddress("TRIG1_Bz", &TRIG1_Bz, &b_TRIG1_Bz);
   fChain->SetBranchAddress("TRIG1_Ex", &TRIG1_Ex, &b_TRIG1_Ex);
   fChain->SetBranchAddress("TRIG1_Ey", &TRIG1_Ey, &b_TRIG1_Ey);
   fChain->SetBranchAddress("TRIG1_Ez", &TRIG1_Ez, &b_TRIG1_Ez);
   fChain->SetBranchAddress("TRIG1_ProperTime", &TRIG1_ProperTime, &b_TRIG1_ProperTime);
   fChain->SetBranchAddress("TRIG1_PathLength", &TRIG1_PathLength, &b_TRIG1_PathLength);
   fChain->SetBranchAddress("TRIG1_PolX", &TRIG1_PolX, &b_TRIG1_PolX);
   fChain->SetBranchAddress("TRIG1_PolY", &TRIG1_PolY, &b_TRIG1_PolY);
   fChain->SetBranchAddress("TRIG1_PolZ", &TRIG1_PolZ, &b_TRIG1_PolZ);
   fChain->SetBranchAddress("TRIG1_InitX", &TRIG1_InitX, &b_TRIG1_InitX);
   fChain->SetBranchAddress("TRIG1_InitY", &TRIG1_InitY, &b_TRIG1_InitY);
   fChain->SetBranchAddress("TRIG1_InitZ", &TRIG1_InitZ, &b_TRIG1_InitZ);
   fChain->SetBranchAddress("TRIG1_InitT", &TRIG1_InitT, &b_TRIG1_InitT);
   fChain->SetBranchAddress("TRIG1_InitKE", &TRIG1_InitKE, &b_TRIG1_InitKE);
   fChain->SetBranchAddress("BPROFEXT_x", &BPROFEXT_x, &b_BPROFEXT_x);
   fChain->SetBranchAddress("BPROFEXT_y", &BPROFEXT_y, &b_BPROFEXT_y);
   fChain->SetBranchAddress("BPROFEXT_z", &BPROFEXT_z, &b_BPROFEXT_z);
   fChain->SetBranchAddress("BPROFEXT_Px", &BPROFEXT_Px, &b_BPROFEXT_Px);
   fChain->SetBranchAddress("BPROFEXT_Py", &BPROFEXT_Py, &b_BPROFEXT_Py);
   fChain->SetBranchAddress("BPROFEXT_Pz", &BPROFEXT_Pz, &b_BPROFEXT_Pz);
   fChain->SetBranchAddress("BPROFEXT_t", &BPROFEXT_t, &b_BPROFEXT_t);
   fChain->SetBranchAddress("BPROFEXT_PDGid", &BPROFEXT_PDGid, &b_BPROFEXT_PDGid);
   fChain->SetBranchAddress("BPROFEXT_EventID", &BPROFEXT_EventID, &b_BPROFEXT_EventID);
   fChain->SetBranchAddress("BPROFEXT_TrackID", &BPROFEXT_TrackID, &b_BPROFEXT_TrackID);
   fChain->SetBranchAddress("BPROFEXT_ParentID", &BPROFEXT_ParentID, &b_BPROFEXT_ParentID);
   fChain->SetBranchAddress("BPROFEXT_Weight", &BPROFEXT_Weight, &b_BPROFEXT_Weight);
   fChain->SetBranchAddress("BPROFEXT_Bx", &BPROFEXT_Bx, &b_BPROFEXT_Bx);
   fChain->SetBranchAddress("BPROFEXT_By", &BPROFEXT_By, &b_BPROFEXT_By);
   fChain->SetBranchAddress("BPROFEXT_Bz", &BPROFEXT_Bz, &b_BPROFEXT_Bz);
   fChain->SetBranchAddress("BPROFEXT_Ex", &BPROFEXT_Ex, &b_BPROFEXT_Ex);
   fChain->SetBranchAddress("BPROFEXT_Ey", &BPROFEXT_Ey, &b_BPROFEXT_Ey);
   fChain->SetBranchAddress("BPROFEXT_Ez", &BPROFEXT_Ez, &b_BPROFEXT_Ez);
   fChain->SetBranchAddress("BPROFEXT_ProperTime", &BPROFEXT_ProperTime, &b_BPROFEXT_ProperTime);
   fChain->SetBranchAddress("BPROFEXT_PathLength", &BPROFEXT_PathLength, &b_BPROFEXT_PathLength);
   fChain->SetBranchAddress("BPROFEXT_PolX", &BPROFEXT_PolX, &b_BPROFEXT_PolX);
   fChain->SetBranchAddress("BPROFEXT_PolY", &BPROFEXT_PolY, &b_BPROFEXT_PolY);
   fChain->SetBranchAddress("BPROFEXT_PolZ", &BPROFEXT_PolZ, &b_BPROFEXT_PolZ);
   fChain->SetBranchAddress("BPROFEXT_InitX", &BPROFEXT_InitX, &b_BPROFEXT_InitX);
   fChain->SetBranchAddress("BPROFEXT_InitY", &BPROFEXT_InitY, &b_BPROFEXT_InitY);
   fChain->SetBranchAddress("BPROFEXT_InitZ", &BPROFEXT_InitZ, &b_BPROFEXT_InitZ);
   fChain->SetBranchAddress("BPROFEXT_InitT", &BPROFEXT_InitT, &b_BPROFEXT_InitT);
   fChain->SetBranchAddress("BPROFEXT_InitKE", &BPROFEXT_InitKE, &b_BPROFEXT_InitKE);
   fChain->SetBranchAddress("BPROF4_x", &BPROF4_x, &b_BPROF4_x);
   fChain->SetBranchAddress("BPROF4_y", &BPROF4_y, &b_BPROF4_y);
   fChain->SetBranchAddress("BPROF4_z", &BPROF4_z, &b_BPROF4_z);
   fChain->SetBranchAddress("BPROF4_Px", &BPROF4_Px, &b_BPROF4_Px);
   fChain->SetBranchAddress("BPROF4_Py", &BPROF4_Py, &b_BPROF4_Py);
   fChain->SetBranchAddress("BPROF4_Pz", &BPROF4_Pz, &b_BPROF4_Pz);
   fChain->SetBranchAddress("BPROF4_t", &BPROF4_t, &b_BPROF4_t);
   fChain->SetBranchAddress("BPROF4_PDGid", &BPROF4_PDGid, &b_BPROF4_PDGid);
   fChain->SetBranchAddress("BPROF4_EventID", &BPROF4_EventID, &b_BPROF4_EventID);
   fChain->SetBranchAddress("BPROF4_TrackID", &BPROF4_TrackID, &b_BPROF4_TrackID);
   fChain->SetBranchAddress("BPROF4_ParentID", &BPROF4_ParentID, &b_BPROF4_ParentID);
   fChain->SetBranchAddress("BPROF4_Weight", &BPROF4_Weight, &b_BPROF4_Weight);
   fChain->SetBranchAddress("BPROF4_Bx", &BPROF4_Bx, &b_BPROF4_Bx);
   fChain->SetBranchAddress("BPROF4_By", &BPROF4_By, &b_BPROF4_By);
   fChain->SetBranchAddress("BPROF4_Bz", &BPROF4_Bz, &b_BPROF4_Bz);
   fChain->SetBranchAddress("BPROF4_Ex", &BPROF4_Ex, &b_BPROF4_Ex);
   fChain->SetBranchAddress("BPROF4_Ey", &BPROF4_Ey, &b_BPROF4_Ey);
   fChain->SetBranchAddress("BPROF4_Ez", &BPROF4_Ez, &b_BPROF4_Ez);
   fChain->SetBranchAddress("BPROF4_ProperTime", &BPROF4_ProperTime, &b_BPROF4_ProperTime);
   fChain->SetBranchAddress("BPROF4_PathLength", &BPROF4_PathLength, &b_BPROF4_PathLength);
   fChain->SetBranchAddress("BPROF4_PolX", &BPROF4_PolX, &b_BPROF4_PolX);
   fChain->SetBranchAddress("BPROF4_PolY", &BPROF4_PolY, &b_BPROF4_PolY);
   fChain->SetBranchAddress("BPROF4_PolZ", &BPROF4_PolZ, &b_BPROF4_PolZ);
   fChain->SetBranchAddress("BPROF4_InitX", &BPROF4_InitX, &b_BPROF4_InitX);
   fChain->SetBranchAddress("BPROF4_InitY", &BPROF4_InitY, &b_BPROF4_InitY);
   fChain->SetBranchAddress("BPROF4_InitZ", &BPROF4_InitZ, &b_BPROF4_InitZ);
   fChain->SetBranchAddress("BPROF4_InitT", &BPROF4_InitT, &b_BPROF4_InitT);
   fChain->SetBranchAddress("BPROF4_InitKE", &BPROF4_InitKE, &b_BPROF4_InitKE);
   fChain->SetBranchAddress("TRIG2_x", &TRIG2_x, &b_TRIG2_x);
   fChain->SetBranchAddress("TRIG2_y", &TRIG2_y, &b_TRIG2_y);
   fChain->SetBranchAddress("TRIG2_z", &TRIG2_z, &b_TRIG2_z);
   fChain->SetBranchAddress("TRIG2_Px", &TRIG2_Px, &b_TRIG2_Px);
   fChain->SetBranchAddress("TRIG2_Py", &TRIG2_Py, &b_TRIG2_Py);
   fChain->SetBranchAddress("TRIG2_Pz", &TRIG2_Pz, &b_TRIG2_Pz);
   fChain->SetBranchAddress("TRIG2_t", &TRIG2_t, &b_TRIG2_t);
   fChain->SetBranchAddress("TRIG2_PDGid", &TRIG2_PDGid, &b_TRIG2_PDGid);
   fChain->SetBranchAddress("TRIG2_EventID", &TRIG2_EventID, &b_TRIG2_EventID);
   fChain->SetBranchAddress("TRIG2_TrackID", &TRIG2_TrackID, &b_TRIG2_TrackID);
   fChain->SetBranchAddress("TRIG2_ParentID", &TRIG2_ParentID, &b_TRIG2_ParentID);
   fChain->SetBranchAddress("TRIG2_Weight", &TRIG2_Weight, &b_TRIG2_Weight);
   fChain->SetBranchAddress("TRIG2_Bx", &TRIG2_Bx, &b_TRIG2_Bx);
   fChain->SetBranchAddress("TRIG2_By", &TRIG2_By, &b_TRIG2_By);
   fChain->SetBranchAddress("TRIG2_Bz", &TRIG2_Bz, &b_TRIG2_Bz);
   fChain->SetBranchAddress("TRIG2_Ex", &TRIG2_Ex, &b_TRIG2_Ex);
   fChain->SetBranchAddress("TRIG2_Ey", &TRIG2_Ey, &b_TRIG2_Ey);
   fChain->SetBranchAddress("TRIG2_Ez", &TRIG2_Ez, &b_TRIG2_Ez);
   fChain->SetBranchAddress("TRIG2_ProperTime", &TRIG2_ProperTime, &b_TRIG2_ProperTime);
   fChain->SetBranchAddress("TRIG2_PathLength", &TRIG2_PathLength, &b_TRIG2_PathLength);
   fChain->SetBranchAddress("TRIG2_PolX", &TRIG2_PolX, &b_TRIG2_PolX);
   fChain->SetBranchAddress("TRIG2_PolY", &TRIG2_PolY, &b_TRIG2_PolY);
   fChain->SetBranchAddress("TRIG2_PolZ", &TRIG2_PolZ, &b_TRIG2_PolZ);
   fChain->SetBranchAddress("TRIG2_InitX", &TRIG2_InitX, &b_TRIG2_InitX);
   fChain->SetBranchAddress("TRIG2_InitY", &TRIG2_InitY, &b_TRIG2_InitY);
   fChain->SetBranchAddress("TRIG2_InitZ", &TRIG2_InitZ, &b_TRIG2_InitZ);
   fChain->SetBranchAddress("TRIG2_InitT", &TRIG2_InitT, &b_TRIG2_InitT);
   fChain->SetBranchAddress("TRIG2_InitKE", &TRIG2_InitKE, &b_TRIG2_InitKE);
   fChain->SetBranchAddress("NP04front_x", &NP04front_x, &b_NP04front_x);
   fChain->SetBranchAddress("NP04front_y", &NP04front_y, &b_NP04front_y);
   fChain->SetBranchAddress("NP04front_z", &NP04front_z, &b_NP04front_z);
   fChain->SetBranchAddress("NP04front_Px", &NP04front_Px, &b_NP04front_Px);
   fChain->SetBranchAddress("NP04front_Py", &NP04front_Py, &b_NP04front_Py);
   fChain->SetBranchAddress("NP04front_Pz", &NP04front_Pz, &b_NP04front_Pz);
   fChain->SetBranchAddress("NP04front_t", &NP04front_t, &b_NP04front_t);
   fChain->SetBranchAddress("NP04front_PDGid", &NP04front_PDGid, &b_NP04front_PDGid);
   fChain->SetBranchAddress("NP04front_EventID", &NP04front_EventID, &b_NP04front_EventID);
   fChain->SetBranchAddress("NP04front_TrackID", &NP04front_TrackID, &b_NP04front_TrackID);
   fChain->SetBranchAddress("NP04front_ParentID", &NP04front_ParentID, &b_NP04front_ParentID);
   fChain->SetBranchAddress("NP04front_Weight", &NP04front_Weight, &b_NP04front_Weight);
   fChain->SetBranchAddress("NP04front_Edep", &NP04front_Edep, &b_NP04front_Edep);
   fChain->SetBranchAddress("NP04front_VisibleEdep", &NP04front_VisibleEdep, &b_NP04front_VisibleEdep);
   fChain->SetBranchAddress("NP04front_Ntracks", &NP04front_Ntracks, &b_NP04front_Ntracks);
   fChain->SetBranchAddress("NP04FieldCage_x", &NP04FieldCage_x, &b_NP04FieldCage_x);
   fChain->SetBranchAddress("NP04FieldCage_y", &NP04FieldCage_y, &b_NP04FieldCage_y);
   fChain->SetBranchAddress("NP04FieldCage_z", &NP04FieldCage_z, &b_NP04FieldCage_z);
   fChain->SetBranchAddress("NP04FieldCage_Px", &NP04FieldCage_Px, &b_NP04FieldCage_Px);
   fChain->SetBranchAddress("NP04FieldCage_Py", &NP04FieldCage_Py, &b_NP04FieldCage_Py);
   fChain->SetBranchAddress("NP04FieldCage_Pz", &NP04FieldCage_Pz, &b_NP04FieldCage_Pz);
   fChain->SetBranchAddress("NP04FieldCage_t", &NP04FieldCage_t, &b_NP04FieldCage_t);
   fChain->SetBranchAddress("NP04FieldCage_PDGid", &NP04FieldCage_PDGid, &b_NP04FieldCage_PDGid);
   fChain->SetBranchAddress("NP04FieldCage_EventID", &NP04FieldCage_EventID, &b_NP04FieldCage_EventID);
   fChain->SetBranchAddress("NP04FieldCage_TrackID", &NP04FieldCage_TrackID, &b_NP04FieldCage_TrackID);
   fChain->SetBranchAddress("NP04FieldCage_ParentID", &NP04FieldCage_ParentID, &b_NP04FieldCage_ParentID);
   fChain->SetBranchAddress("NP04FieldCage_Weight", &NP04FieldCage_Weight, &b_NP04FieldCage_Weight);
   fChain->SetBranchAddress("NP04FieldCage_Edep", &NP04FieldCage_Edep, &b_NP04FieldCage_Edep);
   fChain->SetBranchAddress("NP04FieldCage_VisibleEdep", &NP04FieldCage_VisibleEdep, &b_NP04FieldCage_VisibleEdep);
   fChain->SetBranchAddress("NP04FieldCage_Ntracks", &NP04FieldCage_Ntracks, &b_NP04FieldCage_Ntracks);
   Notify();
}

Bool_t BeamNtuple::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void BeamNtuple::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t BeamNtuple::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef BeamNtuple_cxx
