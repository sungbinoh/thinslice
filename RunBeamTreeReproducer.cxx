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
  BeamTreeReproducer mcpe;
  mcpe.SetOutputFileName(root["mcoutfile"].asString());
 
  TChain *Chain_VirtualDetector_AfterTarget = new TChain();
  Chain_VirtualDetector_AfterTarget -> Add(Form("%s/VirtualDetector/AfterTarget", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_AfterTarget(Chain_VirtualDetector_AfterTarget);

  TChain *Chain_VirtualDetector_TOF1 = new TChain();
  Chain_VirtualDetector_TOF1 -> Add(Form("%s/VirtualDetector/TOF1", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_TOF1(Chain_VirtualDetector_TOF1);

  TChain *Chain_VirtualDetector_COLL1 = new TChain();
  Chain_VirtualDetector_COLL1 -> Add(Form("%s/VirtualDetector/COLL1", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_COLL1(Chain_VirtualDetector_COLL1);

  TChain *Chain_VirtualDetector_BPROF1 = new TChain();
  Chain_VirtualDetector_BPROF1 -> Add(Form("%s/VirtualDetector/BPROF1", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_BPROF1(Chain_VirtualDetector_BPROF1);

  TChain *Chain_VirtualDetector_BPROF2 = new TChain();
  Chain_VirtualDetector_BPROF2 -> Add(Form("%s/VirtualDetector/BPROF2", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_BPROF2(Chain_VirtualDetector_BPROF2);

  TChain *Chain_VirtualDetector_BPROF3 = new TChain();
  Chain_VirtualDetector_BPROF3 -> Add(Form("%s/VirtualDetector/BPROF3", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_BPROF3(Chain_VirtualDetector_BPROF3);

  TChain *Chain_VirtualDetector_TRIG1 = new TChain();
  Chain_VirtualDetector_TRIG1 -> Add(Form("%s/VirtualDetector/TRIG1", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_TRIG1(Chain_VirtualDetector_TRIG1);

  TChain *Chain_VirtualDetector_BPROFEXT = new TChain();
  Chain_VirtualDetector_BPROFEXT -> Add(Form("%s/VirtualDetector/BPROFEXT", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_BPROFEXT(Chain_VirtualDetector_BPROFEXT);

  TChain *Chain_VirtualDetector_BPROF4 = new TChain();
  Chain_VirtualDetector_BPROF4 -> Add(Form("%s/VirtualDetector/BPROF4", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_BPROF4(Chain_VirtualDetector_BPROF4);

  TChain *Chain_VirtualDetector_TRIG2 = new TChain();
  Chain_VirtualDetector_TRIG2 -> Add(Form("%s/VirtualDetector/TRIG2", root["mcfile"].asString().c_str()));
  BeamVirtualDetector evt_TRIG2(Chain_VirtualDetector_TRIG2);

  mcpe.Run(evt_AfterTarget, evt_TOF1, evt_COLL1, evt_BPROF1, evt_BPROF2, evt_BPROF3, evt_TRIG1, evt_BPROFEXT, evt_BPROF4, evt_TRIG2, root["nevents"].asInt());

  return 0;
}
