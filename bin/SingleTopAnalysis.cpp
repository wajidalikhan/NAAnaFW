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
#include "Analysis/NAAnaFW/src/DMTopVariables.h"

using namespace std;

typedef vector<double> vdouble;
typedef vector<float> vfloat;
typedef vector<int> vint;
typedef vector<bool> vbool;
typedef vector<string> vstring;

enum weightedSysts { NOSYST=0, BTAGUP = 1,BTAGDOWN=2,MISTAGUP=3,MISTAGDOWN=4, PUUP=5, PUDOWN=6, LEPUP=7, LEPDOWN=8, MAXSYSTS=9};

struct systWeights{
  void initHistogramsSysts(TH1F** histo, TString name, TString, int, float, float);
  void createFilesSysts(TFile ** allFiles, TString basename, TString opt="RECREATE");
  void fillHistogramsSysts(TH1F** histo, float v, float W, double *wcats= NULL,bool verbose=false);
  void fillHistogramsSysts(TH1F** histo, float v, float W,  float *systWeights, int nFirstSysts=(int)MAXSYSTS, double *wcats=NULL, bool verbose=false);

  void closeFilesSysts(TFile ** allFiles);
  void writeHistogramsSysts(TH1F** histo, TFile ** allFiles );
  void writeSingleHistogramSysts(TH1F* histo,TFile ** allMyFiles);
  void setMax(int max);
  void setMaxNonPDF(int max);
  //  TFile** initFilesSysts();
  void setSystValue(string name, double value, bool mult=false);
  void setSystValue(int systPlace, double value, bool mult=false);
  void setPDFWeights(float * wpdfs, int numPDFs, float wzero=1.0,bool mult=true);
  void setQ2Weights(float q2up, float q2down, float wzero=1.0,bool mult=true);
  void setTWeight(float tweight, float wtotsample=1.0,bool mult=true);
  void setVHFWeight(int vhf,bool mult=true, double shiftval=0.65);
  void setPDFValue(int numPDF, double value);
  double getPDFValue(int numPDF);
  void setWeight(string name, double value, bool mult=false);
  void setWeight(int systPlace, double value, bool mult=false);
  void prepareDefault(bool addDefault, bool addPDF, bool addQ2, bool addTopPt, bool addVHF, bool addTTSplit, int numPDF=102);
  void addSyst(string name);
  void addSystNonPDF(string name);
  void setWCats(double *wcats);

  void addkFact(string name);
  void setkFact(string name,float kfact_nom, float kfact_up,float kfact_down,  bool mult=true);

  void copySysts(systWeights sys);
  void calcPDFHisto(TH1F** histos, TH1F* singleHisto, double scalefactor=1.0, int c = 0);
  void setOnlyNominal(bool useOnlyNominal=false);
  bool onlyNominal;
  bool addPDF, addQ2, addTopPt, addVHF, addTTSplit;
  int maxSysts, maxSystsNonPDF;
  int nPDF;
  int nCategories;
  float weightedSysts[150];
  double wCats[10];
  string weightedNames[150];
  string categoriesNames[10];
  };

struct lept{
  TLorentzVector p4;
  int flavour;
  };

struct by_pt{
    bool operator()(lept const &l1, lept const &l2){
    return l1.p4.Pt() > l2.p4.Pt();
    }
  };

bool pt_ordered (double i,double j) { return (j<i); }

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
  bool doSynch=false;
  if(doSynch){
    fileout.open(reportName.c_str(),ios::in | ios::out | ios::trunc);
    fileout<<"RunNumber EvtNumber Lumi "<<std::endl;
  }
  TString outfile = "test/"+sample + "_" +cat+"_"+channel+".root";
  
  TString treePath = "DMTreesDumper/ttDM__noSyst";
  TString treePathNEvents = "DMTreesDumper/WeightHistory";

  TChain chain(treePath);
  chain.AddFileInfoList(fc.GetList());
  TChain chainNEvents(treePathNEvents);
  chainNEvents.AddFileInfoList(fc.GetList());
  Int_t nEventsTot = (Int_t)chainNEvents.GetEntries();

  TH1D totalWeightTop("w_top_total","Top pt reweighting: overall sample weight",2000,0,2.0);
  chainNEvents.Project("w_top_total","Event_T_Weight","Event_T_Weight!=1.00");
  double topWeight=totalWeightTop.GetMean();
  cout << "totaltopweight is "<< topWeight<<endl;
  if(topWeight==0)topWeight=1;

  Int_t nEvents = (Int_t)chain.GetEntries();
  std::cout<<"Info: Number of Events: "<<nEvents<< endl;
  //nEvents = std::min(nEvents, 1000);
  
  TString weightedSystsNames (weightedSysts sy);
  
  systWeights systZero,syst0BM,syst1BM,syst2BM; 
  int maxSysts=0; 
  int sizeMax=50;
  int muSize, elSize, jetSize;
  //float passTrigHT(0.), Ht(0.);
  float Ht(0.), mt(0.);
  float runNumber(0.), lumiSec(0.);
  double evtNumber(0.);
  float jetE[sizeMax], jetPt[sizeMax], jetPhi[sizeMax], jetEta[sizeMax];
  float jet_iscsvl[sizeMax], jetIsCSVM[sizeMax], jet_iscsvt[sizeMax],jetIsLoose[sizeMax],jetIsTight[sizeMax],jetak4chs_csvv2[sizeMax];
  float nTightMuons, nTightElectrons, nVetoElectrons, nLooseMuons, nJets,nCSVJets;//, nCSVLJets;
  float  nGoodPV, nPV, numTrueInt;
  float metPt[1],metPhi[1],metPx[1],metPy[1];
  
  double met,metpx,metpy;
  float muE[sizeMax], muPt[sizeMax], muEta[sizeMax], muIso[sizeMax], muIsTight[sizeMax],muIsLoose[sizeMax], muPhi[sizeMax],elPt[sizeMax];
  float muCharge[sizeMax], elCharge[sizeMax];

  float elEta[sizeMax], elIso[sizeMax], elIsTight[sizeMax], elIsVeto[sizeMax], elE[sizeMax], elPhi[sizeMax], elPassesDRmu[sizeMax];    
  float lep1Pt(0), lep1Eta(0), lep1Phi(0), lep1E(0), lep2Pt(0), lep2Eta(0), lep2Phi(0),  lep2E(0), lep1Flavour(0), lep1Charge(0), lep2Flavour(0), lep2Charge(0);
  int nPDF=102;

  float LHEWeightSign[1] = {1.};
  float w(1.),w_q2up(1.),w_q2down(1.),w_zero(1.),w_top(1.);
  float w_pdfs[nPDF];

  float bWeight0CSVM, bWeight0CSVMBTagUp,bWeight0CSVMBTagDown,bWeight0CSVMMisTagUp,bWeight0CSVMMisTagDown;
  float bWeight1CSVM, bWeight1CSVMBTagUp,bWeight1CSVMBTagDown,bWeight1CSVMMisTagUp,bWeight1CSVMMisTagDown;
  float bWeight2CSVM, bWeight2CSVMBTagUp,bWeight2CSVMBTagDown,bWeight2CSVMMisTagUp,bWeight2CSVMMisTagDown;

  //  float bWeight0CSVT, bWeight0CSVTBTagUp,bWeight0CSVTBTagDown,bWeight0CSVTMisTagUp,bWeight0CSVTMisTagDown;
  //  float bWeight1CSVT, bWeight1CSVTBTagUp,bWeight1CSVTBTagDown,bWeight1CSVTMisTagUp,bWeight1CSVTMisTagDown;
  //  float bWeight2CSVT, bWeight2CSVTBTagUp,bWeight2CSVTBTagDown,bWeight2CSVTMisTagUp,bWeight2CSVTMisTagDown;
  
  bool onlyNominal=false;
  systZero.setOnlyNominal(onlyNominal);

  bool addPDF=false,addQ2=false,addTopPt=false,addVHF=false,addTTSplit=false;
  if(sample=="TT"){
    addTTSplit=false;
    }


  systZero.prepareDefault(true, addQ2, addPDF, addTopPt,addVHF,addTTSplit);
  //Here must add all systematics we want to put there so that the size of the syst vector is set

  //this maxSysts is to set the number of systematics we do want IMPORTANT FOR ALL INITIALIZATIONS
  maxSysts= systZero.maxSysts;
  TFile * allMyFiles[maxSysts];
  
  systZero.createFilesSysts(allMyFiles,"./res/"+sample + "_" +channel);
  //addWZNLO=false, 
  //addPDF=true;
  //addQ2=true;
  //addTopPt=true;


  //addVHF=false;
  //addWZNLO=false;
  //if(isData=="DATA"){addPDF=false, addQ2=false;addVHF=false;addTTSplit=false;} 

  if(isData=="MC"){chain.SetBranchAddress("Event_LHEWeight0", &w_zero); 
  }
  
  if(addQ2){
    chain.SetBranchAddress("Event_LHEWeight4", &w_q2up);
    chain.SetBranchAddress("Event_LHEWeight8", &w_q2down);
  }
  if(addPDF){
    for (int p = 1; p <= nPDF;++p){
      stringstream pdfss;
      pdfss<<(p+8);
      string pstr =(pdfss.str());
      chain.SetBranchAddress(("Event_LHEWeight"+pstr).c_str(), &w_pdfs[p-1]);
    //    cout << "Event_LHEWeight"+pstr<<w_pdfs<<endl;
    }
    chain.SetBranchAddress("Event_T_Weight",&w_top);
    //    chain.SetBranchAddress("Event_T_Weight",&w_top);
  }

  chain.SetBranchAddress("electronsTight_E", elE);
  chain.SetBranchAddress("electronsTight_Phi", elPhi);
  chain.SetBranchAddress("electronsTight_Eta", elEta);
  chain.SetBranchAddress("electronsTight_Pt", elPt);
  chain.SetBranchAddress("electronsTight_Charge", elCharge);
  chain.SetBranchAddress("electronsTight_Iso03", elIso);
  chain.SetBranchAddress("electronsTight_vidTight", elIsTight);
  chain.SetBranchAddress("electronsTight_vidVeto", elIsVeto);
  chain.SetBranchAddress("electronsTight_PassesDRmu", elPassesDRmu);
  chain.SetBranchAddress("electronsTight_size", &elSize);

  chain.SetBranchAddress("muonsTight_E", muE);
  chain.SetBranchAddress("muonsTight_Phi", muPhi);
  chain.SetBranchAddress("muonsTight_Eta", muEta);
  chain.SetBranchAddress("muonsTight_Pt", muPt);
  chain.SetBranchAddress("muonsTight_Charge", muCharge);
  chain.SetBranchAddress("muonsTight_Iso04", muIso);
  chain.SetBranchAddress("muonsTight_IsTightMuon", muIsTight);
  chain.SetBranchAddress("muonsTight_IsLooseMuon", muIsLoose);
  chain.SetBranchAddress("muonsTight_size", &muSize);

  chain.SetBranchAddress("jetsAK4CHSTight_E",      &jetE);
  chain.SetBranchAddress("jetsAK4CHSTight_Pt",     &jetPt);
  chain.SetBranchAddress("jetsAK4CHSTight_Phi",    &jetPhi);
  chain.SetBranchAddress("jetsAK4CHSTight_Eta",    &jetEta);
  chain.SetBranchAddress("jetsAK4CHSTight_size",   &jetSize);
  chain.SetBranchAddress("jetsAK4CHSTight_IsCSVL", &jet_iscsvl);
  chain.SetBranchAddress("jetsAK4CHSTight_IsCSVM", &jetIsCSVM);
  chain.SetBranchAddress("jetsAK4CHSTight_IsCSVT", &jet_iscsvt);
  chain.SetBranchAddress("jetsAK4CHSTight_IsLoose",&jetIsLoose);
  chain.SetBranchAddress("jetsAK4CHSTight_IsTight",&jetIsTight);
  chain.SetBranchAddress("jetsAK4CHSTight_CSVv2",  &jetak4chs_csvv2);
  
  chain.SetBranchAddress("Event_RunNumber", &runNumber);
  chain.SetBranchAddress("Event_LumiBlock", &lumiSec);
  chain.SetBranchAddress("Event_EventNumber", &evtNumber);
  
  chain.SetBranchAddress("Event_Ht", &Ht);
  chain.SetBranchAddress("Event_mt", &mt);
  chain.SetBranchAddress("Event_nGoodPV",&nGoodPV);
  chain.SetBranchAddress("Event_nPV",&nPV);
  chain.SetBranchAddress("Event_nTruePV",&numTrueInt);
 
  chain.SetBranchAddress("Event_Lepton1_Pt", &lep1Pt);
  chain.SetBranchAddress("Event_Lepton1_Eta", &lep1Eta);
  chain.SetBranchAddress("Event_Lepton1_Phi", &lep1Phi);
  chain.SetBranchAddress("Event_Lepton1_E", &lep1E);
  chain.SetBranchAddress("Event_Lepton1_Charge", &lep1Charge);
  chain.SetBranchAddress("Event_Lepton1_Flavour", &lep1Flavour);

  chain.SetBranchAddress("Event_Lepton2_Pt", &lep2Pt);
  chain.SetBranchAddress("Event_Lepton2_Eta", &lep2Eta);
  chain.SetBranchAddress("Event_Lepton2_Phi", &lep2Phi);
  chain.SetBranchAddress("Event_Lepton2_E", &lep2E);
  chain.SetBranchAddress("Event_Lepton2_Charge", &lep2Charge);
  chain.SetBranchAddress("Event_Lepton2_Flavour", &lep2Flavour);

  chain.SetBranchAddress("Event_nTightElectrons",&nTightElectrons);
  chain.SetBranchAddress("Event_nVetoElectrons",&nVetoElectrons);
  chain.SetBranchAddress("Event_nLooseMuons",&nLooseMuons);
  chain.SetBranchAddress("Event_nTightMuons",&nTightMuons);


  //  float bWeight0CSVM, bWeight0CSVMBTagUp,bWeight0CSVMBTagDown,bWeight0CSVMMisTagUp,bWeight0CSVMMisTagDown;
  //  float bWeight1CSVM, bWeight1CSVMBTagUp,bWeight1CSVMBTagDown,bWeight1CSVMMisTagUp,bWeight1CSVMMisTagDown;
  //  float bWeight2CSVM, bWeight2CSVMBTagUp,bWeight2CSVMBTagDown,bWeight2CSVMMisTagUp,bWeight2CSVMMisTagDown;
  //b-weight settings
  chain.SetBranchAddress("Event_bWeight2CSVMTight", &bWeight2CSVM);
  chain.SetBranchAddress("Event_bWeightMisTagUp2CSVMTight", &bWeight2CSVMMisTagUp);
  chain.SetBranchAddress("Event_bWeightMisTagDown2CSVMTight", &bWeight2CSVMMisTagDown);
  chain.SetBranchAddress("Event_bWeightBTagUp2CSVMTight", &bWeight2CSVMBTagUp);
  chain.SetBranchAddress("Event_bWeightBTagDown2CSVMTight", &bWeight2CSVMBTagDown);

  chain.SetBranchAddress("Event_bWeight0CSVMTight", &bWeight0CSVM);

  chain.SetBranchAddress("Event_bWeightMisTagUp0CSVMTight", &bWeight0CSVMMisTagUp);
  chain.SetBranchAddress("Event_bWeightMisTagDown0CSVMTight", &bWeight0CSVMMisTagDown);
  chain.SetBranchAddress("Event_bWeightBTagUp0CSVMTight", &bWeight0CSVMBTagUp);
  chain.SetBranchAddress("Event_bWeightBTagDown0CSVMTight", &bWeight0CSVMBTagDown);

  chain.SetBranchAddress("Event_bWeight1CSVMTight", &bWeight1CSVM);
  chain.SetBranchAddress("Event_bWeightMisTagUp1CSVMTight", &bWeight1CSVMMisTagUp);
  chain.SetBranchAddress("Event_bWeightMisTagDown1CSVMTight", &bWeight1CSVMMisTagDown);
  chain.SetBranchAddress("Event_bWeightBTagUp1CSVMTight", &bWeight1CSVMBTagUp);
  chain.SetBranchAddress("Event_bWeightBTagDown1CSVMTight", &bWeight1CSVMBTagDown);


  if(isData=="MC"){
    chain.SetBranchAddress("metFull_CorrPt",metPt);
    chain.SetBranchAddress("metFull_CorrPhi",metPhi);}
  else {
    chain.SetBranchAddress("metFull_CorrPt",metPt);
    chain.SetBranchAddress("metFull_CorrPhi",metPhi);}

  chain.SetBranchAddress("metFull_Px",metPx);
  chain.SetBranchAddress("metFull_Py",metPy);

  /********************************************************************************/
  /**************                    Histogram booking              ***************/
  /********************************************************************************/
  TH1F *h_cutFlow = new TH1F("h_cutFlow","cutflow",10,-0.5,9.5);
  TH1F *h_weight_sign = new TH1F("h_weight_sign","weight correction factor before all selections, gives the effective number of events",2,-2.0,2.0);
  TH1F *h_weightZero[maxSysts];      systZero.initHistogramsSysts(h_weightZero,"h_weightZero","weight before all selections, normalized by the total number of events",2000,0,10.0);
  TH1F *h_nPV[maxSysts];      systZero.initHistogramsSysts(h_nPV,"h_nPV","nPV",102,0,51);
  TH1F *h_nGoodPV[maxSysts];  systZero.initHistogramsSysts(h_nGoodPV,"h_nGoodPV","nGoodPV",92,0,46);
  TH1F *h_nTruePV[maxSysts];  systZero.initHistogramsSysts(h_nTruePV,"h_nTruePV","nTruePV",80,0,40);
  
  TH1F *h_nJets[maxSysts];    systZero.initHistogramsSysts(h_nJets,"h_nJets","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets[maxSysts];   systZero.initHistogramsSysts(h_nbJets,"h_nbJets","Number of tight b-jets",11,-0.5,10.5);

  //2j0t 
  TH1F *h_2j0t_jetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j0t_jetpt40_1st,   "h_2j0t_jetpt40_leading","2j0t leading jet Pt ",100,0,500);
  TH1F *h_2j0t_jetpt40_2nd[maxSysts];    systZero.initHistogramsSysts(h_2j0t_jetpt40_2nd,   "h_2j0t_jetpt40_subleading","2j0t Sub. leading jet Pt ",100,0,500);
  TH1F *h_2j0t_mtw[maxSysts];            systZero.initHistogramsSysts(h_2j0t_mtw,           "h_2j0t_mtw",     "2j0t mtw ",100,0,500);
  
  //2j1t Category
  TH1F *h_2j1t_bjetpt[maxSysts];   systZero.initHistogramsSysts(h_2j1t_bjetpt,  "h_2j1t_bjetpt",  "2j1t b jet pt",100,0,500);
  TH1F *h_2j1t_jetpt[maxSysts];    systZero.initHistogramsSysts(h_2j1t_jetpt,   "h_2j1t_jetpt",   "2j1t light jet pt" ,100,0,500);
  TH1F *h_2j1t_jeteta[maxSysts];   systZero.initHistogramsSysts(h_2j1t_jeteta,  "h_2j1t_jeteta",  "2j1t light jet eta",64,-3.2,3.2);
  TH1F *h_2j1t_nMu[maxSysts];      systZero.initHistogramsSysts(h_2j1t_nMu,     "h_2j1t_nMu",     "2j1t Number of tight Muons",13,-0.5,12.5);
  TH1F *h_2j1t_MuPt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_MuPt,    "h_2j1t_MuPt",    "2j1t muon pt ",100,0,500);
  TH1F *h_2j1t_MuEta[maxSysts];    systZero.initHistogramsSysts(h_2j1t_MuEta,   "h_2j1t_MuEta",   "2j1t muon eta ",64,-3.2,3.2);
  TH1F *h_2j1t_MuPhi[maxSysts];    systZero.initHistogramsSysts(h_2j1t_MuPhi,   "h_2j1t_MuPhi",   "2j1t muon phi ",16, -3.2, 3.2);
  TH1F *h_2j1t_MuE[maxSysts];      systZero.initHistogramsSysts(h_2j1t_MuE,     "h_2j1t_MuE",     "2j1t muon e ",100,0,500);
  TH1F *h_2j1t_MuCharge[maxSysts]; systZero.initHistogramsSysts(h_2j1t_MuCharge,"h_2j1t_MuCharge","2j1t muon charge ",2,-1,1);
  TH1F *h_2j1t_mtw[maxSysts];      systZero.initHistogramsSysts(h_2j1t_mtw,     "h_2j1t_mtw",     "2j1t mtw ",100,0,500);
  TH1F *h_2j1t_topMass[maxSysts];  systZero.initHistogramsSysts(h_2j1t_topMass, "h_2j1t_topMass", "2j1t top mass",200,100,500);
  TH1F *h_2j1t_topPt[maxSysts];    systZero.initHistogramsSysts(h_2j1t_topPt,   "h_2j1t_topPt",   "2j1t top pt",200,0,500);
  
  //3j1t 
  TH1F *h_3j1t_bjetpt[maxSysts];   systZero.initHistogramsSysts(h_3j1t_bjetpt,  "h_3j1t_bjetpt",  "3j1t b jet pt ",100,0,500);
  TH1F *h_3j1t_MuPt[maxSysts];     systZero.initHistogramsSysts(h_3j1t_MuPt,    "h_3j1t_MuPt",    "muon pt ",100,0,500);
  TH1F *h_3j1t_MuEta[maxSysts];    systZero.initHistogramsSysts(h_3j1t_MuEta,   "h_3j1t_MuEta",   "muon eta ",64,-3.2,3.2);
  TH1F *h_3j1t_MuPhi[maxSysts];    systZero.initHistogramsSysts(h_3j1t_MuPhi,   "h_3j1t_MuPhi",   "muon phi ",16, -3.2, 3.2);
  TH1F *h_3j1t_MuE[maxSysts];      systZero.initHistogramsSysts(h_3j1t_MuE,     "h_3j1t_MuE",     "muon e ",100,0,500);
  TH1F *h_3j1t_mtw[maxSysts];      systZero.initHistogramsSysts(h_3j1t_mtw,     "h_3j1t_mtw",     "3j1t mtw ",100,0,500);
  
  //3j2t
  TH1F *h_3j2t_bjetpt[maxSysts];   systZero.initHistogramsSysts(h_3j2t_bjetpt,    "h_3j2t_bjetpt",     "3j2t b jet pt ",100,0,500);
  TH1F *h_3j2t_2ndbjetpt[maxSysts];systZero.initHistogramsSysts(h_3j2t_2ndbjetpt, "h_3j2t_2ndbjetpt",  "3j2t sub lead. b jet pt",100,0,500);
  TH1F *h_3j2t_MuPt[maxSysts];     systZero.initHistogramsSysts(h_3j2t_MuPt,      "h_3j2t_MuPt",       "3j2t muon pt ",100,0,500);
  TH1F *h_3j2t_MuEta[maxSysts];    systZero.initHistogramsSysts(h_3j2t_MuEta,     "h_3j2t_MuEta",      "3j2t muon eta ",64,-3.2,3.2);
  TH1F *h_3j2t_MuPhi[maxSysts];    systZero.initHistogramsSysts(h_3j2t_MuPhi,     "h_3j2t_MuPhi",      "3j2t muon phi ",16, -3.2, 3.2);
  TH1F *h_3j2t_MuE[maxSysts];      systZero.initHistogramsSysts(h_3j2t_MuE,       "h_3j2t_MuE",        "3j2t muon e ",100,0,500);
  TH1F *h_3j2t_MuCharge[maxSysts]; systZero.initHistogramsSysts(h_3j2t_MuCharge,  "h_3j2t_MuCharge",   "3j2t muon charge ",2,-1,1);
  TH1F *h_3j2t_mtw[maxSysts];      systZero.initHistogramsSysts(h_3j2t_mtw,       "h_3j2t_mtw",        "3j2t mtw ",100,0,500);

  TopUtilities topUtils;


for(Int_t evt=0; evt<nEvents; evt++ ){
    if(evt%10000==1 ){
    cout<<"Info: Running on event: "<<evt<<endl; 
    }
  chain.GetEntry(evt);
  int maxJetLoop = min(15, jetSize);
  int maxMuLoop = min(6, muSize);
  int maxElLoop = min(6, elSize);



  if(isData=="MC"){
    LHEWeightSign[0]=w_zero/fabs(w_zero);
    w = LHEWeightSign[0];
    
    //       w=1.;
    //    w_pu = LumiWeights_.weight(numTrueInt);
    
    //    w = w * w_pu;
    //    w*=k_fact;
    //    cout << " --> bef top pt " << w << endl;
    
    if(sample=="TTNoTT"){
      w*=w_top/topWeight;
    }

    //    cout << "     after top pt " << w << endl; 
       //w*= ((METturnon13TeV(metPt[0],1,true))/(METturnon13TeV(metPt[0],1,false)));
       //       cout << "vhf "<<vhf << " eventFlavour "<<eventFlavour << endl;
    //vhf=(int)eventFlavour;
       //FIXXX
       //       w_top=1.0;
       //       topWeight=1.0;
       //std::cout<<"K-factor: "<<k_fact<<std::endl;
     }
     
     //added for sync

  systZero.setWeight(0,1.);
  systZero.setWeight("btagUp",1.);
  systZero.setWeight("btagDown",1.);
  systZero.setWeight("mistagUp",1.);
  systZero.setWeight("mistagDown",1.);
  systZero.setWeight("puDown",1.);
  systZero.setWeight("puUp",1.);
  systZero.setWeight("lepDown",1.);
  systZero.setWeight("lepUp",1.);

  
  if(addPDF)systZero.setPDFWeights(w_pdfs,nPDF,w_zero,true);
  if(addQ2)systZero.setQ2Weights(w_q2up,w_q2down,w_zero,true);
  if(addTopPt)systZero.setTWeight(w_top,topWeight,true);
  

  syst0BM.copySysts(systZero);
  syst0BM.setWeight(0,bWeight0CSVM);
  syst0BM.setWeight("btagUp",bWeight0CSVMBTagUp);
  syst0BM.setWeight("btagDown",bWeight0CSVMBTagDown);
  syst0BM.setWeight("mistagUp",bWeight0CSVMMisTagUp);
  syst0BM.setWeight("mistagDown",bWeight0CSVMMisTagDown);

  syst1BM.copySysts(systZero);
  syst1BM.setWeight(0,bWeight1CSVM);
  syst1BM.setWeight("btagUp",bWeight1CSVMBTagUp);
  syst1BM.setWeight("btagDown",bWeight1CSVMBTagDown);
  syst1BM.setWeight("mistagUp",bWeight1CSVMMisTagUp);
  syst1BM.setWeight("mistagDown",bWeight1CSVMMisTagDown);


  syst2BM.copySysts(systZero);
  syst2BM.setWeight(0,bWeight2CSVM);
  syst2BM.setWeight("btagUp",bWeight2CSVMBTagUp);
  syst2BM.setWeight("btagDown",bWeight2CSVMBTagDown);
  syst2BM.setWeight("mistagUp",bWeight2CSVMMisTagUp);
  syst2BM.setWeight("mistagDown",bWeight2CSVMMisTagDown);

  
  met = metPt[0];
  metpx = metPx[0];
  metpy = metPy[0];
  
  TLorentzVector lep1;
  TLorentzVector lep2;
  TLorentzVector lep, mu, el;

  vector<TLorentzVector> tightEl, tightMu;
  int nMu(0.), nEl(0.);//, nVetoEl(0.), nLooseMu(0.);
  for(int e= 0; e<maxElLoop;++e ){
      if(elPt[e]>0){
      el.SetPtEtaPhiE(elPt[e], elEta[e], elPhi[e],elE[e]);
      tightEl.push_back(el);
    
    }
  }//end loop on electrons     
 
  for(int m= 0; m<maxMuLoop;++m ){
      if(muPt[m]>0){
      mu.SetPtEtaPhiE(muPt[m], muEta[m], muPhi[m],muE[m]);
      tightMu.push_back(mu);

    }
  }//end loop on muons
  
  nMu = tightMu.size();
  nEl = tightEl.size();
  
  vector<lept> leptons;
  
  for (size_t e = 0; e < (size_t)tightEl.size(); ++e){
    lept lepton;
    lepton.p4 = tightEl[e];
    lepton.flavour = 11;
    leptons.push_back(lepton);
    }
  
  for (size_t m = 0; m < (size_t)tightMu.size(); ++m){
    lept lepton;
    lepton.p4 = tightMu[m];
    lepton.flavour = 13;
    leptons.push_back(lepton);
    }

  std::sort(leptons.begin(), leptons.end(), by_pt());

  vector<float> jetsPhi;
  vector<TLorentzVector> bjets, jets_nob, jets;    
  
  struct btag{
      TLorentzVector vect;
      float csv;
  };
  
  std::vector<btag> bvects;
    
  struct by_csv{
      bool operator()(btag const &b1, btag const &b2){
      return b1.csv > b2.csv;
      }
  };

  struct by_pt_jet{
      bool operator()(btag const &jet1, btag const &jet2){
      return jet2.vect.Pt() < jet2.vect.Pt();
      }
  };
  
  for (int j = 0; j <maxJetLoop;++j){	
      TLorentzVector jet, all_jets;
      
      if(jetPt[j]>40.0 && abs(jetEta[j])<2.4){ 
      all_jets.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
      jets.push_back(all_jets);
      }   

  if(jetIsCSVM[j]<=0.){
      jet.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
      jets_nob.push_back(jet);
      }
  
  jetsPhi.push_back(jetPhi[j]);
  TLorentzVector bjet;
 
  if( jetPt[j]>40.0 && jetIsTight[j] && jetIsCSVM[j] && abs(jetEta[j])<2.4){
     bjet.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
     bjets.push_back(bjet);
     btag b;
     b.vect = bjet;
     b.csv = jetak4chs_csvv2[j];
     bvects.push_back(b);
     
   }  
}// End of jet loop 
 
  nJets = jetsPhi.size();
  nCSVJets=bjets.size();

  std::sort(bvects.begin(), bvects.end(), by_pt_jet()); 
  
  //2j
  if((jets.size() == 2 && bjets.size() == 0 && nMu == 1 && nEl==0)){
  for(size_t i= 0; i< (size_t)jets.size();++i ){
    //cout <<runNumber<<"   "<<evtNumber<<"   "<<"jet["<<i<< "] "<<jets[i].Pt()<<"   bjetsize "<< bjets.size() <<"   Mu["<<i<< "] "<<tightMu[i].Pt()<< std::endl;
    if(i==0)syst0BM.fillHistogramsSysts(h_2j0t_jetpt40_1st,jets[i].Pt(),w);  
    if(i==1)syst0BM.fillHistogramsSysts(h_2j0t_jetpt40_2nd,jets[i].Pt(),w);  
   }
    
  for(size_t i = 0; i < (size_t)tightMu.size();++i ){
        if((tightMu.size())< 2 && tightMu[i].Pt()>20 ){
          TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
          float phi_lmet = fabs(deltaPhi(tightMu[i].Phi(), metPhi[0]) );
          mt = sqrt(2* tightMu[i].Pt() * met* ( 1- cos(phi_lmet)));
          syst0BM.fillHistogramsSysts(h_2j0t_mtw,mt,w);
          }
      }
  }
  //2j1t
  if((jets.size() == 2 && bjets.size() == 1 && nMu == 1 && nEl==0)){
    for(size_t i= 0; i< (size_t)bjets.size();++i ){
      //      cout << " w "<< w <<" b weight 1 "<< bWeight1CSVM << endl;
      if(i==0){
	syst1BM.fillHistogramsSysts(h_2j1t_bjetpt,bjets[i].Pt(),w,NULL,false);
	//	cout << bjets[i].Pt()<< endl;
      }
    }
  
  for (size_t j= 0; j< (size_t)jets.size();++j ){
    if(!(jetIsCSVM[j]) && jetPt[j]>40.0 && abs(jetEta[j])<2.4){
      syst1BM.fillHistogramsSysts(h_2j1t_jetpt,jets[j].Pt(),w);
      syst1BM.fillHistogramsSysts(h_2j1t_jeteta,jets[j].Eta(),w);
    }  
  }
  
  for(size_t i = 0; i < (size_t)tightMu.size();++i ){
   if(tightMu[i].Pt()>20){
   syst1BM.fillHistogramsSysts(h_2j1t_MuPt,tightMu[i].Pt(),w);
   syst1BM.fillHistogramsSysts(h_2j1t_MuEta,tightMu[i].Eta(),w);
   syst1BM.fillHistogramsSysts(h_2j1t_MuPhi,tightMu[i].Phi(),w);
   syst1BM.fillHistogramsSysts(h_2j1t_MuE,tightMu[i].E(),w);
   
   if((tightMu.size())<2 ){
        TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
        float phi_lmet = fabs(deltaPhi(tightMu[i].Phi(), metPhi[0]) );
        mt = sqrt(2* tightMu[i].Pt() * met* ( 1- cos(phi_lmet)));
        syst1BM.fillHistogramsSysts(h_2j1t_mtw,mt,w);
        }
      }
    }    
  
  //Adding top related variables
  vector<TLorentzVector> tops= topUtils.top4Momenta(tightMu, bjets ,metpx, metpy);
  if(tops.size()!=0){
    double mass = tops.at(0).M();
    double pt = tops.at(0).Pt();
    syst1BM.fillHistogramsSysts(h_2j1t_topMass,mass,w);
    syst1BM.fillHistogramsSysts(h_2j1t_topPt,pt,w);
    }
  }
  
  //3j1t   
  if((jets.size() == 3 && bjets.size() ==1 && nMu==1 && nEl==0)) {
  for(size_t i= 0; i< (size_t)bjets.size();++i ){
   if(i==0)systZero.fillHistogramsSysts(h_3j1t_bjetpt,bjets[i].Pt(),w);
   } 
  
  for(size_t i = 0; i < (size_t)tightMu.size();++i ){
        syst1BM.fillHistogramsSysts(h_3j1t_MuPt,tightMu[i].Pt(),w);
        syst1BM.fillHistogramsSysts(h_3j1t_MuEta,tightMu[i].Eta(),w);
        syst1BM.fillHistogramsSysts(h_3j1t_MuPhi,tightMu[i].Phi(),w);
        syst1BM.fillHistogramsSysts(h_3j1t_MuE,tightMu[i].E(),w);
        if((tightMu.size())<2 ){
          TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
          float phi_lmet = fabs(deltaPhi(tightMu[i].Phi(), metPhi[0]) );
          mt = sqrt(2* tightMu[i].Pt() * met* ( 1- cos(phi_lmet)));
          systZero.fillHistogramsSysts(h_3j1t_mtw,mt,w);
          }
      }
  }
  
  //3j2t
  if((jets.size() == 3 && bjets.size()==2 && nMu==1 && nEl==0)){
    for(size_t i=0; i< (size_t)bjets.size();++i ){
      //      cout << " ,bWeight2CSVM "<< <<",bWeight2CSVMBTagUp "<<bWeight2CSVMBTagUp<<endl;

      if(i==0)syst2BM.fillHistogramsSysts(h_3j2t_bjetpt,bjets[i].Pt(),w,NULL,false);
      if(i==1)syst2BM.fillHistogramsSysts(h_3j2t_2ndbjetpt,bjets[i].Pt(),w);
      }
      for(size_t i = 0; i < (size_t)tightMu.size();++i ){
        if(tightMu[i].Pt()>20){
        syst2BM.fillHistogramsSysts(h_3j2t_MuPt,tightMu[i].Pt(),w);
        syst2BM.fillHistogramsSysts(h_3j2t_MuEta,tightMu[i].Eta(),w);
        syst2BM.fillHistogramsSysts(h_3j2t_MuPhi,tightMu[i].Phi(),w);
        syst2BM.fillHistogramsSysts(h_3j2t_MuE,tightMu[i].E(),w);
        if((tightMu.size())<2 ){
          TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
          float phi_lmet = fabs(deltaPhi(tightMu[i].Phi(), metPhi[0]) );
          mt = sqrt(2* tightMu[i].Pt() * met* ( 1- cos(phi_lmet)));
          systZero.fillHistogramsSysts(h_3j2t_mtw,mt,w);
          }
        }
      }
    } 
  
  systZero.fillHistogramsSysts(h_nJets,nJets,w); 
  systZero.fillHistogramsSysts(h_nbJets,nCSVJets,w); 
  systZero.fillHistogramsSysts(h_nPV,nPV,w);
  systZero.fillHistogramsSysts(h_nGoodPV,nPV,w);
  systZero.fillHistogramsSysts(h_nTruePV,numTrueInt,w);
  
  if(nTightElectrons != nEl)cout << "warning! problem with tight el"<<endl;
  if(nTightMuons != nMu)cout << "warning! problem with tight mu"<<endl;
  
  nTightElectrons = nEl;
  nTightMuons = nMu;    


  if(doSynch)fileout<<std::fixed<<std::setprecision(0)<<runNumber<<"   "<<evtNumber<<"   "<<lumiSec<<"   "<<std::endl;
  
  }//end of loop over events 
  
 if(doSynch)fileout.close();  //return h

 //Step 0:get only the negative weights version, 
 chainNEvents.Project("h_weight_sign","Event_LHEWeight10/abs(Event_LHEWeight10)"); // I messed up and weight zero is not  saved in the full chain... but this one should do the trick as well! Just do h_weight_sign->GetMean() and multiply this to get the correct number :)

 // TO FIX! ADD PROJECTION OF ALL PDFS TO THE EVENT!
 // systZero.projectAllPDF(h_weightZero,chainNEvents);

  
  // Cut Flow
 
   
 h_cutFlow->SetBinContent(0,nEventsTot*h_weight_sign->GetMean());//Underflow: number of events pre-preselection.
 h_cutFlow->SetBinContent(1,nEvents);
  h_cutFlow->GetXaxis()->SetBinLabel(1,"no selection");
  h_cutFlow->SetBinContent(2, 0);
  h_cutFlow->GetXaxis()->SetBinLabel(2, "trigger");
  
  //Write the Histogramms here  
  systZero.writeSingleHistogramSysts(h_cutFlow, allMyFiles); 
  systZero.writeSingleHistogramSysts(h_weight_sign, allMyFiles); 
  systZero.writeHistogramsSysts(h_nPV, allMyFiles); 
  systZero.writeHistogramsSysts(h_nGoodPV, allMyFiles); 
  systZero.writeHistogramsSysts(h_nTruePV, allMyFiles); 
  
  systZero.writeHistogramsSysts(h_2j0t_jetpt40_1st, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_jetpt40_2nd, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtw,   allMyFiles); 
  
  systZero.writeHistogramsSysts(h_2j1t_bjetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_jetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_jeteta,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuPt,    allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuEta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuPhi,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuE,     allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuCharge,allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtw,     allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_topMass, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_topPt,   allMyFiles); 
  
  systZero.writeHistogramsSysts(h_3j1t_bjetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_MuPt,    allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_MuEta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_MuPhi,   allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_MuE,     allMyFiles); 
  systZero.writeHistogramsSysts(h_3j1t_mtw,     allMyFiles); 
  
  systZero.writeHistogramsSysts(h_3j2t_bjetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_2ndbjetpt, allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuPt,    allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuEta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuPhi,   allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuE,     allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_MuCharge,allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_mtw,     allMyFiles); 
  
  
  systZero.writeHistogramsSysts(h_nJets, allMyFiles); 
  systZero.writeHistogramsSysts(h_nbJets, allMyFiles); 


  
  
  systZero.closeFilesSysts(allMyFiles);
  
  std::cout<< "Info: Total No. of events for sample["<<sample<<"] : "<<nEvents<<std::endl;
}

void systWeights::closeFilesSysts(TFile *filesout[(int)MAXSYSTS]){
  int MAX= this->maxSystsNonPDF;
  //  int MAXTOT= this->maxSysts;
  
  bool useOnlyNominal = this->onlyNominal;
  for (int c = 0; c < this->nCategories; c++){
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      //cout << "c is now "<< c << " sy "<< sy << " location "<< sy+(MAXTOT+1)*c <<" is histo there? " << histo[sy+(MAXTOT+1)*c] << " file location "<<sy+(MAX+1)*c << " is file there "<< filesout[sy+(MAX+1)*c]<< endl;
      if(!(!useOnlyNominal || sy==0)) continue;
      filesout[(int)sy+(MAX+1)*(c)]->Close();
    }
    
    if(this->addPDF){
      if(!useOnlyNominal){
	filesout[MAX+(MAX+1)*(c)]->Close();
      }
    }
  }
}
void systWeights::setOnlyNominal(bool useOnlyNominal){
  this->onlyNominal=useOnlyNominal;
}
void systWeights::setSystValue(string name, double value, bool mult){
  float zerofact=1.0;
  if(mult)zerofact=this->weightedSysts[0];
  int MAX = this->maxSysts;
  for(int sy=0;sy<(int)MAX;++sy){
    if(this->weightedNames[(int)sy] ==name){
      this->weightedSysts[(int)sy] =value*zerofact;
    }
  }
}


void systWeights::setSystValue(int place, double value, bool mult){
  float zerofact=1.0;
  if(mult)zerofact=this->weightedSysts[0];
  this->weightedSysts[place] =value*zerofact;
}

void systWeights::setWeight(string name, double value, bool mult){
  this->setSystValue(name, value, mult);
}

void systWeights::setWeight(int place, double value, bool mult){
  this->setSystValue(place, value, mult);
}

void systWeights::initHistogramsSysts(TH1F** histo,TString name, TString title, int nbins, float min, float max){
  for (int c = 0; c < this->nCategories; c++){
    int MAX = this->maxSysts;
    bool useOnlyNominal = this->onlyNominal;
    TString cname= (this->categoriesNames[c]).c_str();
    for(int sy=0;sy<(int)MAX;++sy){
      TString ns= (this->weightedNames[sy]).c_str();
      if(sy==0){
          if(c==0) histo[sy+((MAX+1)*(c))]=new TH1F(name,title,nbins,min,max);
          else histo[sy+((MAX+1)*(c))]=new TH1F(name+"_"+cname,title,nbins,min,max);
          }
      
      if(sy!=0 && !useOnlyNominal) {
          if(c==0)histo[sy+(MAX+1)*c]=new TH1F(name+"_"+ns,title,nbins,min,max);
          else histo[sy+(MAX+1)*c]=new TH1F(name+"_"+ns+"_"+cname,title,nbins,min,max);
          }
          //cout << " initialized histogram "<< histo[sy+(MAX+1)*c]->GetName() <<" sy " << sy << " c  "<< c <<" location " << sy+(MAX+1)*c << endl;
    }
  }
}
void writeHistogramsSysts(TH1F* histo[(int)MAXSYSTS], TFile *filesout[(int)MAXSYSTS], bool useOnlyNominal=false){
  for(int sy=0;sy<(int)MAXSYSTS;++sy){
    //cout << " writing histo "<< histo[(int)sy]->GetName()<< " in file "<< filesout[(int)sy]->GetName()<<endl;;
    //TString ns= weightedSystsNames((weightedSysts)sy);
    filesout[(int)sy]->cd();
    histo[sy]->Write(histo[0]->GetName());
    //    histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
  }
}
void systWeights::writeHistogramsSysts(TH1F** histo, TFile **filesout){
  int MAX= this->maxSystsNonPDF;
  int MAXTOT= this->maxSysts;
  bool useOnlyNominal = this->onlyNominal;
  for (int c = 0; c < this->nCategories; c++){
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      //cout << "c is now "<< c << " sy "<< sy << " location "<< sy+(MAXTOT+1)*c <<" is histo there? " << histo[sy+(MAXTOT+1)*c] << " file location "<<sy+(MAX+1)*c << " is file there "<< filesout[sy+(MAX+1)*c]<< endl;
      
      //cout << " writing histo "<< histo[sy+(MAXTOT+1)*c]->GetName()<< " in file "<< filesout[sy+(MAX+1)*c]->GetName()<<endl;;
      
      //TString ns= weightedSystsNames((weightedSysts)sy);
      if(!(!useOnlyNominal || sy==0)) continue;
      filesout[(int)sy+(MAX+1)*(c)]->cd();
      if(this->addPDF){
      //if(this->weightedNames[sy]=="pdf_totalUp")calcPDFHisto(histo, histo[sy+(MAXTOT+1)*(c)],1.0,c);
      //if(this->weightedNames[sy]=="pdf_totalDown")calcPDFHisto(histo, histo[sy+(MAXTOT+1)*(c)],-1.0,c);
      ;      //this->
      }

      histo[sy+(MAXTOT+1)*c]->Write(histo[0]->GetName());
      //histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
      }
      
      if(this->addPDF){
      if(!useOnlyNominal){
        filesout[MAX+(MAX+1)*(c)]->cd();
        //cout << " file max is "<< filesout[MAX+(MAX+1)*c]->GetName()<<endl;
        //int npdf=this->maxSysts-this->maxSystsNonPdf;
        int MAXPDF=this->maxSysts;
        for(int sy=MAX;sy<MAXPDF;++sy){
        //    cout << " writing sy "<<sy+(MAXTOT+1)*c<<endl;
        //    cout << " histo is there? "<< histo[sy+(MAXTOT+1)*c]<<endl;
        histo[sy+(MAXTOT+1)*(c)]->Write();
        //    cout << " written sy "<< histo[sy+(MAXTOT+1)*c]->GetName()<<endl;
        }
      }
    }

  }

}

void systWeights::createFilesSysts(  TFile ** allFiles, TString basename, TString opt){
  for (int c = 0; c < this->nCategories; c++){
    int MAX = this->maxSystsNonPDF;
    int MAXTOT = this->maxSystsNonPDF;
    bool useOnlyNominal = this->onlyNominal;
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      TString ns= (this->weightedNames[(int)sy]);
      cout << " creating file for syst "<< ns<<endl;
      if (c!=0)     cout << " category is "<< c<<endl;
      cout << "onlynominal is "<<useOnlyNominal<<endl;

      if(sy==0){
	      cout<<" filename is "<< basename+ns+cname+".root"<<endl;
	      allFiles[sy+(MAX+1)*c]= TFile::Open((basename+ns+cname+".root"), opt);
        }
      
      else{
	      if(!useOnlyNominal){
	          //if((ns!="1lep") && (ns!="2lep")&& (ns!="0lep")){
	          cout<<" filename is "<< basename+ns+cname+".root"<<endl;
	          allFiles[sy+(MAX+1)*c]= TFile::Open((basename+"_"+ns+cname+".root"), opt);
	          }
          }
        //TFile *outTree = TFile::Open(("trees/tree_"+outFileName).c_str(), "RECREATE");
        //cout << " created file at c "<< c << " s "<< sy << " location "<< sy+(MAXTOT+1)*c<< " fname "<<allFiles[sy+(MAXTOT+1)*c]->GetName()<<endl;   
        }
      
      if(this->addPDF){
      if(!useOnlyNominal)allFiles[MAX+((MAX+1)*c)]= TFile::Open((basename+"_pdf"+cname+".root"), opt);
      //cout << " created file at c "<< c << " s "<< MAX+(MAX+1)*c << " location "<< MAX+(MAX+1)*c<<endl;
      cout<< " fname "<<allFiles[MAX+(MAXTOT+1)*c]->GetName()<<endl;
        }
      }
    //return allFiles;
}



 
void systWeights::fillHistogramsSysts(TH1F** histo, float v, float w, double * wcats, bool verbose ){
  if(wcats==NULL){
    wcats=this->wCats;
  }
  for (int c = 0; c < this->nCategories; c++){
    if(verbose){
      cout<< " in cat loop "<< c<<endl;
      cout<< " value "<< wcats[c] <<endl;
    }
    int MAX = this->maxSysts;
    bool useOnlyNominal = this->onlyNominal;
    for(int sy=0;sy<(int)MAX;++sy){
      if(verbose){
	cout<< " in syst loop "<< sy<< endl;
	cout<<" value "<< this->weightedSysts[(int)sy] <<endl ;
      }
      if(sy!=0&& useOnlyNominal)continue;
      float ws = (this->weightedSysts[(int)sy])*wcats[c];
      if(verbose)cout << " filling histogram "<< histo[(int)sy]->GetName() << " with value "<< v <<" and weight "<< w <<" ws "<< ws<<endl;
      histo[sy+(MAX+1)*(c)]->Fill(v, w*ws);
    }
  }
}

//void  initHistogramsSysts (TH1F* histo[(int)MAXSYSTS],TString name, TString title, int nbins, float min, float max, bool useOnlyNominal=false){
//    for(int sy=0;sy<(int)MAXSYSTS;++sy){
//          //TString ns= weightedSystsNames((weightedSysts)sy);
//          histo[sy]=new TH1F(name,title,nbins,min,max);
//          }
//        }

void systWeights::setWCats(double * wcats){
  for(int i =0;i<this->nCategories;++i){
    //    cout << "setting wcat #"<< i << " to be "<<wcats[i]<<endl;
    this->wCats[i]=wcats[i];
  }
 
}

void systWeights::copySysts(systWeights sys){
  for(int i =0; i < sys.maxSysts;++i){
    this->weightedNames[i]=sys.weightedNames[i];
    this->weightedSysts[i]=sys.weightedSysts[i];

  }
  this->setOnlyNominal(sys.onlyNominal);
  this->setMax(sys.maxSysts);
  this->setMaxNonPDF(sys.maxSystsNonPDF);
  this->nPDF=sys.nPDF;
  this->nCategories=sys.nCategories;  
  this->addQ2=sys.addQ2;
  this->addPDF=sys.addPDF;
  this->addTopPt=sys.addTopPt;
  this->addVHF=sys.addVHF;
  this->addTTSplit=sys.addTTSplit;
  this->setWCats(sys.wCats);

}

void systWeights::setMax(int max){
  this->maxSysts =  max;
  }
void systWeights::setMaxNonPDF(int max){
  this->maxSystsNonPDF =  max;
  }

void systWeights::prepareDefault(bool addDefault, bool addQ2, bool addPDF, bool addTopPt,bool addVHF, bool addTTSplit, int numPDF){ 
  this->addPDF=addPDF;
  this->addQ2=addQ2;
  this->addTopPt=addTopPt;
  this->addVHF=addVHF;
  this->addTTSplit=addTTSplit;
  this->nPDF=numPDF;
  this->nCategories=1;
  categoriesNames[0]="";
  this->wCats[0]=1.0;
  if(addDefault){
    this->weightedNames[0]="";
    this->weightedNames[1]="btagUp";
    this->weightedNames[2]="btagDown";
    this->weightedNames[3]="mistagUp";
    this->weightedNames[4]="mistagDown";
    this->weightedNames[5]="puUp";
    this->weightedNames[6]="puDown";
    this->weightedNames[7]="lepUp";
    this->weightedNames[8]="lepDown";
    this->setMax(9);
    this->setMaxNonPDF(9);
    this->weightedNames[this->maxSysts]="";
  }
  if(addQ2){
    this->weightedNames[this->maxSysts]= "q2Up";
    this->weightedNames[this->maxSysts+1]= "q2Down";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    this->weightedNames[this->maxSysts]= "";
  }

  if(addTopPt){
    this->weightedNames[this->maxSysts]="topPtWeightUp";
    this->weightedNames[this->maxSysts+1]="topPtWeightDown";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    this->weightedNames[this->maxSysts]= "";
  }

  if(addVHF){
    this->weightedNames[this->maxSysts]="VHFWeightUp";
    this->weightedNames[this->maxSysts+1]="VHFWeightDown";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    this->weightedNames[this->maxSysts]= "";
  }

  if(addTTSplit){
    this->nCategories=4;
    categoriesNames[1]="TT0lep";
    categoriesNames[2]="TT1lep";
    categoriesNames[3]="TT2lep";
    this->wCats[1]=1.0;
    this->wCats[2]=1.0;
    this->wCats[3]=1.0;
    }

  if(addPDF){
    this->weightedNames[this->maxSysts]= "pdf_totalUp";
    this->weightedNames[this->maxSysts+1]= "pdf_totalDown";
    this->weightedNames[this->maxSysts+2]= "pdf_asUp";
    this->weightedNames[this->maxSysts+3]= "pdf_asDown";
    this->weightedNames[this->maxSysts+4]= "pdf_zmUp";
    this->weightedNames[this->maxSysts+5]= "pdf_zmDown";
    this->setMax(this->maxSysts+6);
    this->setMaxNonPDF(this->maxSystsNonPDF+6);
    int nPDF=this->nPDF;
    for(int i =0; i < nPDF;++i){
      stringstream ss;
      ss<< i+1;
      this->weightedNames[i+this->maxSysts]= "pdf"+ss.str();
      }
    this->setMax(maxSysts+nPDF);
    this->weightedNames[this->maxSysts]= "";
    }
}

TString  weightedSystsNames (weightedSysts sy){
  switch(sy){
  case NOSYST : return "";
  case BTAGUP : return "btagUp";
  case BTAGDOWN : return "btagDown";
  case MISTAGUP : return "mistagUp";
  case MISTAGDOWN : return "mistagDown";
  case PUUP : return "puUp";
  case PUDOWN : return "puDown";
  case LEPUP : return "lepUp";
  case LEPDOWN : return "lepDown";
  //case ISOUP : return "isoUp";
  //case ISODOWN : return "isoDown";
  //case TRIGUP : return "trigUp";
  //case TRIGDOWN : return "trigDown";
  case MAXSYSTS : return "";
  }
  return "noSyst";
}

void systWeights::writeSingleHistogramSysts(TH1F* histo, TFile **filesout){  
  int MAX= this->maxSystsNonPDF;
  bool useOnlyNominal = this->onlyNominal;
  for (int c = 0; c < this->nCategories; c++){
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      if(!(!useOnlyNominal || sy==0)) continue;
      //cout << " writing histo "<< histo->GetName()<< " in file "<< filesout[(int)sy]->GetName()<<endl;;
      filesout[(int)sy+(MAX+1)*c]->cd();
      histo->Write();
      //    histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
    }
    if(this->addPDF){
      if(!useOnlyNominal){
        filesout[MAX+(MAX+1)*c]->cd();
        int MAXPDF=this->maxSysts;
        for(int sy=MAX;sy<MAXPDF;++sy){
        //cout << " writing sy "<< histo[sy]->GetName()<<endl;
        histo->Write();
        //cout << " written sy "<< histo[sy]->GetName()<<endl;
        }
      }
    }
  }
}


void systWeights::setPDFWeights(float * wpdfs, int numPDFs, float wzero,bool mult){
  float zerofact=1.0;
  if(mult)zerofact=this->weightedSysts[0];
  for (int i = 0; i < numPDFs; ++i){
    this->setPDFValue(i,zerofact*wpdfs[i]/wzero);
  }
  this->setSystValue("pdf_asUp", this->getPDFValue(this->nPDF-2)/wzero);
  this->setSystValue("pdf_asDown", zerofact);
  this->setSystValue("pdf_zmUp", this->getPDFValue(this->nPDF-1)/wzero);
  this->setSystValue("pdf_zmDown", zerofact);
  this->setSystValue("pdf_totalUp", zerofact);
  this->setSystValue("pdf_totalDown", zerofact);
}

//void systWeights::setTWeight(float tweight, float totalweight){
void systWeights::setTWeight(float tweight, float wtotsample,bool mult){
  float zerofact=1.0;
  //  cout << " weighted syst 0 is "<< weightedSysts[0]<<endl;
  if(mult)zerofact=this->weightedSysts[0];
  this->setSystValue("topPtWeightUp", zerofact*tweight/wtotsample);
  this->setSystValue("topPtWeightDown", zerofact/tweight*wtotsample);
}

void systWeights::setVHFWeight(int vhf,bool mult,double shiftval){
  float zerofact=1.0;
  double w_shift=0.0;
  //  cout << "vhf is "<<vhf<<endl;
  if (vhf>1)w_shift=shiftval;
  //  cout << " weighted syst 0 is "<< weightedSysts[0]<<endl;
  if(mult)zerofact=this->weightedSysts[0];
  this->setSystValue("VHFWeightUp", zerofact*(1+w_shift));
  this->setSystValue("VHFWeightDown", zerofact*(1-w_shift));
}


void systWeights::setQ2Weights(float q2up, float q2down, float wzero, bool mult){
  float zerofact=1.0;
  if(mult){
    zerofact=this->weightedSysts[0];
    //    cout <<  "zerofact "<< zerofact << endl;
  }
  //  cout <<  "zerofact "<< zerofact << " q2up weight "<< q2up/wzero << " tot to fill "<< zerofact*q2up/wzero<<endl;
  //  cout <<  "zerofact "<< zerofact << " q2down weight "<< q2down/wzero << " tot to fill "<< zerofact*q2down/wzero<<endl;
  this->setSystValue("q2Up", zerofact*q2up/wzero);
  this->setSystValue("q2Down", zerofact*q2down/wzero);
}

double systWeights::getPDFValue(int numPDF){
  if(!addPDF){ cout << "error! No PDF used, this will do nothing."<<endl;return 0.;}
  int MIN = this->maxSystsNonPDF;
  return (double)this->weightedSysts[numPDF+MIN];

}
void systWeights::setPDFValue(int numPDF, double w){
  if(!addPDF){ cout << "error! No PDF used, this will do nothing."<<endl;return;}
  int MIN = this->maxSystsNonPDF;
  this->weightedSysts[numPDF+MIN]=w;

}

void systWeights::calcPDFHisto(TH1F** histo, TH1F* singleHisto, double scalefactor, int c){//EXPERIMENTAL
  if(!addPDF){ cout << "error! No PDF used, this will do nothing."<<endl;return;}
  int MAX = this->maxSysts;
  //  for (int c = 0; c < this->nCategories; c++){
    int MIN = this->maxSystsNonPDF+(MAX+1)*c;
    for(int b = 0; b< singleHisto->GetNbinsX();++b){
      float val = singleHisto->GetBinContent(b);
      //      cout << "bin # "<<b << " val "<<val<<endl;
      float mean = 0, devst=0;
      //      cout << "name is "<< singleHisto->GetName()<<endl;
      for(int i = 0; i<this->nPDF;++i ){
	//cout<< " now looking at pdf # "<<i<< " coordinate is "<< MIN+i<<endl;
	//	cout << "is histo there? "<< histo[i+MIN]<<endl;
	//	cout << " histo should be "<< (histo[i+MIN])->GetName()<<endl;
	mean = mean+ histo[i+MIN]->GetBinContent(b);
      }
      mean = mean/this->nPDF;
      //mean = val;//mean/this->nPDF;
      for(int i = 0; i<this->nPDF;++i ){
	devst+=(mean-histo[i+MIN]->GetBinContent(b))*(mean-histo[i+MIN]->GetBinContent(b));
      }
      devst= sqrt(devst/this->nPDF);
      singleHisto->SetBinContent(b,val+devst*scalefactor);
      //      singleHisto->SetBinContent(b,mean+devst*scalefactor);
    }
    //}
}

