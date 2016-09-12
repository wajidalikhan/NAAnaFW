#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include <vector>
#include <assert.h>
#include <TMVA/Reader.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <cassert>
#include <sstream>
#include <string>
#include "TFileCollection.h"
#include "THashList.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Math/interface/LorentzVector.h"
//#include "Tprime/TprimeAnalysis/interface/Weights.h"
#include "DataFormats/Math/interface/deltaR.h"

using namespace std;

typedef vector<double> vdouble;
typedef vector<float> vfloat;
typedef vector<int> vint;
typedef vector<bool> vbool;
typedef vector<string> vstring;

void callme(){
  std::cout<<" NaN value"<<std::endl;
}

int main(int argc, char **argv) {

  std::cout<<"Info: Starting to run ..."<<endl;
  std::cout<<"--------------------"<<endl; 
  std::cout<<"Usage:\t "<<argv[0]<<" "<<argv[1]<<" "<<argv[2]<<" "<<argv[3]<<" "<<argv[4]<<" "<<argv[5]<<" "<<argv[6]<<" "<<argv[7]<<endl;
  

  string sample(argv[1]) ;
  std::cout<<"\t Sample: "<<sample<<endl;

  string path(argv[2]);
  std::cout<<"\t File list to open: "<<path<<endl;

  string channel(argv[3]);
  std::cout<<"\t Channel: "<<channel<<endl;
  
  string cat(argv[4]);
  std::cout<<"\t Category: " <<cat<<endl;
  
  string sys(argv[5]);
  std::cout<<"\t Systematics: "<<sys<<endl;

  string sync(argv[6]);
  std::cout<<"\t Synchro: "<<sync<<endl;

  string isData(argv[7]);
  std::cout<<"\t isData: "<<isData<<endl;

  TString path_ = path ; 
  std::cout<<"\t File to open: "<<path_<<endl;
  std::cout<<"--------------------"<<endl; 
  
  std::cout << "Info: Loading file collection from " << path << std::endl;
  TFileCollection fc(sample.c_str(),sample.c_str(),path.c_str());
  std::cout << "Info: Files found : " << fc.GetNFiles() << std::endl;

  string reportName = "SelectedEvents_"+channel+"_"+cat+"_"+sample+".txt";
  ofstream fileout;
  fileout.open(reportName.c_str(),ios::in | ios::out | ios::trunc);
  fileout<<"RunNumber EvtNumber Lumi "<<std::endl;

  TString outfile = "test/"+sample + "_" +cat+"_"+channel+".root";
  TFile fout(outfile, "recreate");
  
  //std::cout<<"File to open: "<<path_<<endl;
  TString treePath = "DMTreesDumper/ttDM__noSyst";

  TChain chain(treePath);
  chain.AddFileInfoList(fc.GetList());

  Int_t nEvents = (Int_t)chain.GetEntries();
  std::cout<<"Info: Number of Events: "<<nEvents<< endl;
  //nEvents = std::min(nEvents, 1000000);
  
  int sizeMax=50;
  int jetSize;
  //float passTrigHT(0.), Ht(0.);
  float Ht(0.);
  float runNumber(0.), lumiSec(0.);
  double evtNumber(0.);
  float jete[sizeMax], jetpt[sizeMax], jetphi[sizeMax], jeteta[sizeMax];

  chain.SetBranchAddress("jetsAK4CHS_E",    &jete);
  chain.SetBranchAddress("jetsAK4CHS_Pt",   &jetpt);
  chain.SetBranchAddress("jetsAK4CHS_Phi",  &jetphi);
  chain.SetBranchAddress("jetsAK4CHS_Eta",  &jeteta);
  chain.SetBranchAddress("jetsAK4CHS_size", &jetSize);
  
  chain.SetBranchAddress("Event_Ht", &Ht);
  chain.SetBranchAddress("Event_RunNumber", &runNumber);
  chain.SetBranchAddress("Event_LumiBlock", &lumiSec);
  chain.SetBranchAddress("Event_EventNumber", &evtNumber);

  /********************************************************************************/
  /**************                    Histogram booking              ***************/
  /********************************************************************************/

  
  TH1F *h_jete    = new TH1F("h_jete",    "AK4 jet e",    100, 0, 1500);
  TH1F *h_jetpt   = new TH1F("h_jetpt",   "AK4 jet pt",   100, 0, 1500);
  TH1F *h_jetphi  = new TH1F("h_jetphi",  "AK4 jet phi",  100, -3.5, 3.5 );
  TH1F *h_jeteta  = new TH1F("h_jeteta",  "AK4 jet eta",  100, -5.2, 5.2);
  TH1F *h_jetsn   = new TH1F("h_jetsn",   "AK4 njets",    15, 0, 15);
  
  TH1F *h_ht      = new TH1F("h_ht",      "ht",           100, 0, 1500);

  TH1F *h_cutFlow = new TH1F("h_cutFlow","cutFlow",7,-0.5,6.5);
  
  for(Int_t i=0; i<nEvents; i++ ){
    
    if(i%100000==1 ){
      cout<<"Info: Running on event: "<<i<<endl; 
      }
    chain.GetEntry(i);

    int maxJetLoop = min(10, jetSize);
	    
    for (int j = 0; j <maxJetLoop;++j){	 
      if(jetpt[j]>30. &&  fabs(jeteta[j])<4) {
	      h_jete->Fill(jete[j]);
	      h_jetpt->Fill(jetpt[j]);
	      h_jeteta->Fill(jeteta[j]);
	      h_jetphi->Fill(jetphi[j]);
	      h_jetsn->Fill(jetSize);
        }	    
    }//end of the jet loop
    
    h_ht->Fill(Ht);      
    fileout
      <<std::fixed<<std::setprecision(0)
      <<runNumber<<"   "
      <<evtNumber<<"   "
      <<lumiSec<<"   "
      <<std::endl;
    
  }//end of loop over events 
  
  // Cut Flow
  h_cutFlow->SetBinContent(1,nEvents);
  h_cutFlow->GetXaxis()->SetBinLabel(1,"no selection");
 
  fout.cd();
  
  h_jete->Write();
  h_jetpt->Write();
  h_jeteta->Write();
  h_jetphi->Write();
  h_jetsn->Write();
  h_ht->Write();

  h_cutFlow->Write();
 
  fout.Close();
  fileout.close();
  
  std::cout<< "Info: Total No. of events for sample["<<sample<<"] : "<<nEvents<<std::endl;
  //return h
  
}//end of main
