#include "PionXsec.h"
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
      std::cout << "Usage: RunPionXsec " <<
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
  mcchain->Add(Form("%s/pduneana/beamana", root["mcfile"].asString().c_str()));
  
  anavar mcevt(mcchain);

  PionXsec mcpe;
  mcpe.SetOutputFileName(root["mcoutfile"].asString());
  mcpe.Run(mcevt, root["nevents"].asInt());
  /*
  TChain *datachain = new TChain();
  datachain->Add(Form("%s/pduneana/beamana", root["datafile"].asString().c_str()));

  anavar dataevt(datachain);

  PionXsec datape;
  datape.SetOutputFileName(root["dataoutfile"].asString());
  datape.Run(dataevt, root["nevents"].asInt());
  */
  return 0;

}
