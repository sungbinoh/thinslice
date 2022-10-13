#include "BeamTreeReproducer.h"
#include "EventType.h"
#include "EventSelection.h"
#include "TChain.h"
#include "json/json.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

int main(int argc, char ** argv){

  bool found_config = false;

  string config_file;

  for (int iArg = 1; iArg < argc; iArg++) {
    if (!strcasecmp(argv[iArg],"-c")) {
     config_file = argv[++iArg];
     found_config = true;
    }
    if (!strcasecmp(argv[iArg],"-h")) {
      std::cout << "Usage: RunBeamTreeReproducer " <<
                   "-c config.json " << std::endl;
      return 1;
    }
  }

  if (!found_config){
    cout<<"Error: no configuration file was provided! Please provide with '-i'" <<endl;
    return 1;
  }

  Json::Value root;
  ifstream file(config_file);
  file >> root;
  cout<<root<<endl;
  //cout<<root["mcfile"].asString()<<endl;

  TChain *mcchain = new TChain();
  mcchain->Add(Form("%s/NTuples/GoodParticle", root["mcfile"].asString().c_str()));
  
  BeamNtuple mcevt(mcchain);

  BeamTreeReproducer mcpe;
  mcpe.SetOutputFileName(root["mcoutfile"].asString());
  mcpe.Run(mcevt, root["nevents"].asInt());

  TChain *Chain_VirtualDetector_AfterTarget = new TChain();
  Chain_VirtualDetector_AfterTarget -> Add(Form("%s/VirtualDetector/AfterTarget", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_AfterTarget(Chain_VirtualDetector_AfterTarget);
  mcpe.Run(evt_AfterTarget, root["nevents"].asInt(), "AfterTarget");

  TChain *Chain_VirtualDetector_COLL1 = new TChain();
  Chain_VirtualDetector_COLL1 -> Add(Form("%s/VirtualDetector/COLL1", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_COLL1(Chain_VirtualDetector_COLL1);
  mcpe.Run(evt_COLL1, root["nevents"].asInt(), "COLL1");

  TChain *Chain_VirtualDetector_BPROF1 = new TChain();
  Chain_VirtualDetector_BPROF1 -> Add(Form("%s/VirtualDetector/BPROF1", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_BPROF1(Chain_VirtualDetector_BPROF1);
  mcpe.Run(evt_BPROF1, root["nevents"].asInt(), "BPROF1");

  TChain *Chain_VirtualDetector_BPROF3 = new TChain();
  Chain_VirtualDetector_BPROF3 -> Add(Form("%s/VirtualDetector/BPROF3", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_BPROF3(Chain_VirtualDetector_BPROF3);
  mcpe.Run(evt_BPROF3, root["nevents"].asInt(), "BPROF3");

  TChain *Chain_VirtualDetector_TRIG1 = new TChain();
  Chain_VirtualDetector_TRIG1 -> Add(Form("%s/VirtualDetector/TRIG1", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_TRIG1(Chain_VirtualDetector_TRIG1);
  mcpe.Run(evt_TRIG1, root["nevents"].asInt(), "TRIG1");

  TChain *Chain_VirtualDetector_BPROF4 = new TChain();
  Chain_VirtualDetector_BPROF4 -> Add(Form("%s/VirtualDetector/BPROF4", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_BPROF4(Chain_VirtualDetector_BPROF4);
  mcpe.Run(evt_BPROF4, root["nevents"].asInt(), "BPROF4");

  TChain *Chain_VirtualDetector_TRIG2 = new TChain();
  Chain_VirtualDetector_TRIG2 -> Add(Form("%s/VirtualDetector/TRIG2", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_TRIG2(Chain_VirtualDetector_TRIG2);
  mcpe.Run(evt_TRIG2, root["nevents"].asInt(), "TRIG2");

  return 0;
}
