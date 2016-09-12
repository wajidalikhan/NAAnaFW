 /**
 *\ttDManalysis.cpp
 *
 * Applies the ttDM analysis workflow for semileptonic and fullhadronic channels. 
 * 
 *
 * \Authors: A. Orso M. Iorio, A. de Cosa
 *
 *
 *\version  $Id:
 *
 **************
 **/


#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TBranch.h"
#include "TH1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLorentzVector.h"
#include <TMVA/Reader.h>
#include "TFileCollection.h"
#include "THashList.h"
#include "TBenchmark.h"
#include "TStopwatch.h"

#include <cstdlib> 
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

#include "ttDM/TopTagResolved/interface/KinematicFitter.hh"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "ttDM/newTTDManalysis/interface/Weights.h"
#include "ttDM/newTTDManalysis/interface/kFactors.h"
#include "METTriggerTurnonEfficiency.C"
//#include "./interface/Weights.h"
#include "ttDM/newTTDManalysis/interface/MT2Utility.h"
#include "ttDM/newTTDManalysis/interface/mt2w_bisect.h"
#include "ttDM/newTTDManalysis/interface/mt2bl_bisect.h"
#include "ttDM/newTTDManalysis/interface/Mt2Com_bisect.h"
#include "DataFormats/Math/interface/deltaR.h"

//#include "ttDM/newTTDManalysis/METTriggerTurnonEfficiency.C"

using namespace std;

typedef vector<double> vdouble;
typedef vector<float> vfloat;
typedef vector<int> vint;
typedef vector<bool> vbool;
typedef vector<string> vstring;

//enum weightedSysts { NOSYST=0, BTAGUP = 1,BTAGDOWN=2,MISTAGUP=3,MISTAGDOWN=4, PUUP=5, PUDOWN=6, LEPIDUP=7, LEPIDDOWN=8, ISOUP=9, ISODOWN=10, TRIGUP=11, TRIGDOWN=12, MAXSYSTS=13};
enum weightedSysts { NOSYST=0, BTAGUP = 1,BTAGDOWN=2,MISTAGUP=3,MISTAGDOWN=4, PUUP=5, PUDOWN=6, LEPUP=7, LEPDOWN=8, MAXSYSTS=9};
enum theoSysts {SCALEUP=101,SCALEDOWN=102, NNPDF1=100, NNPDF2=102};
int wLimit =150;


struct systWeights{
  void initHistogramsSysts(TH1F** histo, TString name, TString, int, float, float);
  void createFilesSysts(TFile ** allFiles, TString basename, TString opt="RECREATE");
  void fillHistogramsSysts(TH1F** histo, float v, float W, double *wcats= NULL,bool verbose=false);
  void fillHistogramsSysts(TH1F** histo, float v, float W,  float *systWeights, int nFirstSysts=(int)MAXSYSTS, double *wcats=NULL, bool verbose=false);
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


struct TopFitResults {
  int idxB, idxJ1, idxJ2;
  float mva;
  float bdt_qgid1;
  float bdt_qgid2;
  float bdt_dphij1b;
  float bdt_dphij2b;
  float bdt_drj1b;
  float bdt_drj2b;
  float bdt_bjcsv;
  float bdt_jet1csv; 
  float bdt_jet2csv;
  float bdt_prob; 

  TLorentzVector topPreFit, WPreFit;

  float WMPhiPreFit;
  float TMPhiPreFit;
  float BMPhiPreFit;
  float WBPhiPreFit;
  float LMPhiPreFit;
		 
  TLorentzVector topPostFit, WPostFit;

  float WMPhiPostFit;
  float TMPhiPostFit;
  float BMPhiPostFit;
  float WBPhiPostFit;
  float LMPhiPostFit;
     
};


void callme(){
  std::cout<<" NaN value"<<std::endl;
}

int main(int argc, char **argv) {

  std::cout<<"Let's start"<<endl;

  // Create a benchmark object to track execution time 
  TBenchmark bench;

  string sample(argv[1]) ;
  std::cout<<"sample: "<<sample<<endl;

  //  TString path(argv[2]);
  //  std::cout<<"File to open: "<<path<<endl;


  string path(argv[2]);
  std::cout<<"File list to open: "<<path<<endl;


  string channel(argv[3]);
  std::cout<<"channel: "<<channel<<endl;

  string sys(argv[4]);
  std::cout<<"systematics: "<<sys<<endl;

  string sync(argv[5]);
  std::cout<<"synchro: "<<sync<<endl;

  string isData(argv[6]);
  std::cout<<"isData: "<<isData<<endl;

  // TString lepton(argv[7]);
  // std::cout<<"Divide leptons: "<<lepton<<endl;

  std::string outdir(argv[7]);
  std::cout<<"Output directory: "<<outdir<<endl;

  // All histograms have weights enabled
  TH1::SetDefaultSumw2 (kTRUE);
  

  std::cout << "Loading file collection from " << path << std::endl;
  TFileCollection fc(sample.c_str(),sample.c_str(),path.c_str());
  std::cout << "Files found : " << fc.GetNFiles() << std::endl;


  // TString path_ = path + "/trees*.root";
  // std::cout<<"File to open: "<<path_<<endl;
  // //TString path_ = path + "/*.root";

  //  bool saveTree=false;
  //  saveTree=true;
  TH1F * initproduct(TH1F * h_A,TH1F* h_B, int rebinA = 1, int rebinB=1,double integral = -1.);
  TH1F * makeproduct(TH1F * h_A,TH1F* h_B, int rebinA = 1, int rebinB=1,double integral = -1.);
  TH1F * makeproduct(TH2F * h);

  TString syststr = "";
  string syststrname = "";
  // if (lepton != "" && lepton != "any")lepstr= "_"+ lepton;
  if (sys != "noSys"){syststr= "_"+ sys; syststrname= "_"+ sys;}
  
  string reportDir = outdir+"/txt";
  string reportName = reportDir+"/SelectedEvents_"+channel+"_"+sample+syststrname+".txt";


  ofstream fileout;
  fileout.open(reportName.c_str(),ios::in | ios::out | ios::trunc);
  //fileout<<"RunNumber  EvtNumber Lumi  nLooseLep   MET   NJets  nbjets  mindphi  PassMETFilter  HadronicTrigger"<<std::endl;
  fileout<<"RunNumber  EvtNumber Lumi  MET   NJets  nbjets w_pu k_fact b_weight el_weight mu_weight met_trig "<<std::endl;
  
  string reportName_presel = reportDir+"/SelectedEvents_Preselection_"+channel+"_"+sample+syststrname+".txt";
  ofstream fileout_presel;
  fileout_presel.open(reportName_presel.c_str(),ios::in | ios::out  | ios::trunc);
  fileout_presel<<"RunNumber  EvtNumber Lumi  nLooseLep   MET   NJets  nbjets  mindphi  PassMETFilter  HadronicTrigger"<<std::endl;
  
  string reportName_CR3 = reportDir+"/SelectedEvents_CR3_"+channel+"_"+sample+syststrname+".txt";
  ofstream fileout_CR3;
  fileout_CR3.open(reportName_CR3.c_str(),ios::in | ios::out  | ios::trunc);
  fileout_CR3<<"RunNumber  EvtNumber Lumi  nLooseLep   MET   NJets  nbjets  mindphi  PassMETFilter  HadronicTrigger"<<std::endl;
  
  string reportName_CR5 = reportDir+"/SelectedEvents_CR5_"+channel+"_"+sample+syststrname+".txt";
  ofstream fileout_CR5;
  fileout_CR5.open(reportName_CR5.c_str(),ios::in | ios::out  | ios::trunc);
  //  fileout_CR5<<"RunNumber  EvtNumber Lumi  nLooseLep   MET   NJets  nbjets  mindphi  PassMETFilter  HadronicTrigger"<<std::endl;
  fileout_CR5<<"RunNumber   EvtNumber Lumi  nLooseLep   MET   NJets  nbjets  mindphi b_weight b_weightUp b_weightDown b_mistagUp b_mistagDown"<<std::endl;
  
  string reportName_CR6 = reportDir+"/SelectedEvents_CR6_"+channel+"_"+sample+syststrname+".txt";
  ofstream fileout_CR6;
  fileout_CR6.open(reportName_CR6.c_str(),ios::in | ios::out  | ios::trunc);
  fileout_CR6<<"RunNumber  EvtNumber Lumi  nLooseLep   MET   NJets  nbjets  mindphi  PassMETFilter  HadronicTrigger"<<std::endl;
  
  string reportName_CR7 = reportDir+"/SelectedEvents_CR7_"+channel+"_"+sample+syststrname+".txt";
  ofstream fileout_CR7;
  fileout_CR7.open(reportName_CR7.c_str(),ios::in | ios::out  | ios::trunc);
  fileout_CR7<<"RunNumber  EvtNumber Lumi  nLooseLep   MET   NJets  nbjets  mindphi  PassMETFilter  HadronicTrigger"<<std::endl;
  
  string reportName_CRouttop = reportDir+"/SelectedEvents_CRouttop_"+channel+"_"+sample+syststrname+".txt";
  ofstream fileout_CRouttop;
  fileout_CRouttop.open(reportName_CRouttop.c_str(),ios::in | ios::out  | ios::trunc);
  fileout_CRouttop<<"RunNumber  EvtNumber Lumi  nLooseLep   MET   NJets  nbjets  mindphi  PassMETFilter  HadronicTrigger"<<std::endl;
  
  string reportName_CR1lep = reportDir+"/SelectedEvents_CR1lep_"+channel+"_"+sample+syststrname+".txt";
  ofstream fileout_CR1lep;
  fileout_CR1lep.open(reportName_CR1lep.c_str(),ios::in | ios::out  | ios::trunc);
  fileout_CR1lep<<"RunNumber  EvtNumber Lumi  nLooseLep   MET   NJets  nbjets  mindphi  PassMETFilter  HadronicTrigger"<<std::endl;
  
  string reportName_CR2lep = reportDir+"/SelectedEvents_CR2lep_"+channel+"_"+sample+syststrname+".txt";
  ofstream fileout_CR2lep;
  fileout_CR2lep.open(reportName_CR2lep.c_str(),ios::in | ios::out  | ios::trunc);
  fileout_CR2lep<<"RunNumber  EvtNumber Lumi  nLooseLep   MET   NJets  nbjets j0pt j1pt j2pt  mindphi kfact PassMETFilter  HadronicTrigger"<<std::endl;
  
  //  bool useAdditionalWeightedSysts=false;

  TString weightedSystsNames (weightedSysts sy);
  void initHistogramsSysts(TH1F* histo[(int)MAXSYSTS],TString, TString, int, float, float , bool useOnlyNominal=false);
  void createFilesSysts(  TFile * allFiles[(int)MAXSYSTS], TString basename, bool useOnlyNominal=false, TString opt="RECREATE");
  void fillHistogramsSysts(TH1F* histo[(int)MAXSYSTS], float v, float W, float systWeight[(int)MAXSYSTS] , bool useOnlyNominal=false);
  void writeHistogramsSysts(TH1F* histo[(int)MAXSYSTS], TFile * allFiles[(int)MAXSYSTS] , bool useOnlyNominal=false);
  void writeSingleHistogramSysts(TH1F* histo,TFile * allMyFiles[(int)MAXSYSTS],bool useOnlyNominal=false);

  //  void fillHistogramSysts();
  //  TFile * allMyFiles[(int)MAXSYSTS];
  systWeights systZero,syst12B,syst12BL,syst0B,syst2B;
  int maxSysts=0;
  bool addPDF=false,addQ2=false,addTopPt=false,addVHF=false,addWZNLO=false, addTTSplit=false;
  addPDF=true;
  addQ2=true;
  addTopPt=true;

  if(sample=="TT"){
    addTTSplit=true;
  }

  addVHF=true;
  addWZNLO=true;
  int nPDF=102;
  if(isData=="DATA"){addPDF=false, addQ2=false;addVHF=false;addTTSplit=false;}
  systZero.prepareDefault(true, addQ2, addPDF, addTopPt,addVHF,addTTSplit);

  if(addWZNLO){
    systZero.addkFact("QCDRen");
    systZero.addkFact("QCDFac");
    systZero.addkFact("EWK");
  }
  maxSysts= systZero.maxSysts;
  if(sample=="TT"){
    maxSysts= (systZero.maxSysts+1)*4;
    cout << "maxSysts "<< maxSysts<<endl;
  }
  

  //  int MAX = syst;
  TFile * allMyFiles[maxSysts];
  
  //  TString initdir = "/afs/cern.ch/work/o/oiorio/public/xTTDM/Feb2016/Freezing/";
  // TString initdir = "";
  //  TString outfile = "/afs/cern.ch/user/o/oiorio/public/xAnnapaola/Nov10/res/"+sample + "_" +channel+".root";
  TString outfile = outdir+"/res/"+sample + "_" +channel+".root";
  //  if(sys!="noSys") outfile = "/afs/cern.ch/user/o/oiorio/public/xAnnapaola/Nov10/res/"+sample + "_" +channel+"_"+sys+".root";

  // TString lepstr = "";


  bool onlyNominal=false;
  //if (sys != "noSys")  onlyNominal=true;
  systZero.setOnlyNominal(onlyNominal);
  syst12B.setOnlyNominal(onlyNominal);
  syst12BL.setOnlyNominal(onlyNominal);
  syst0B.setOnlyNominal(onlyNominal);
  syst2B.setOnlyNominal(onlyNominal);

  systZero.createFilesSysts(allMyFiles,outdir+"/res/"+sample + "_" +channel+syststr);
  
  //  systZero.createFilesSysts(allMyFiles0Lep,outdir+"/res/"+sample + "_0lep_" +channel+syststr);
  //  systZero.createFilesSysts(allMyFiles1Lep,outdir+"/res/"+sample + "_1lep_" +channel+syststr);
  //  systZero.createFilesSysts(allMyFiles2Lep,outdir+"/res/"+sample + "_2lep_" +channel+syststr);


  //  createFilesSysts(allMyFiles,"/afs/cern.ch/user/o/oiorio/public/xAnnapaola/Nov10/res/"+sample + "_" +channel+syststr+lepstr);
  // systZero.createFilesSysts(allMyFiles,outdir+"res/"+sample + "_" +channel+syststr+lepstr);
  
  //  if(sys=="noSys"){//outfile = "/afs/cern.ch/user/o/oiorio/public/xAnnapaola/Nov10/res/"+sample + "_" +channel+"_"+sys+".root";
  //    createFilesSysts(allMyFiles,"/afs/cern.ch/user/o/oiorio/public/xAnnapaola/Nov10/res/"+sample + "_" +channel+lepstr);
  //  }
  //  if(sys!="noSys"){//outfile = "/afs/cern.ch/user/o/oiorio/public/xAnnapaola/Nov10/res/"+sample + "_" +channel+"_"+sys+".root";
  //    createFilesSysts(allMyFiles,"/afs/cern.ch/user/o/oiorio/public/xAnnapaola/Nov10/res/"+sample + "_" +channel+"_"+sys+lepstr);
  //  }

  //  TFile fout(outfile, "RECREATE");

  //  std::cout<<"File to open: "<<path_<<endl;
  TString treePath = "DMTreesDumper/ttDM__noSyst";
  TString treePathNEvents = "DMTreesDumper/WeightHistory";
  //  TString treePathNEvents = "DMTreesDumper/EventHistory";
  if(sys=="jesUp") treePath = "DMTreesDumper/ttDM__jes__up";
  else if(sys=="jesDown") treePath = "DMTreesDumper/ttDM__jes__down";
  else if(sys=="jerUp") treePath = "DMTreesDumper/ttDM__jer__up";
  else if(sys=="jerDown")treePath = "DMTreesDumper/ttDM__jer__down";
  else if(sys=="metUnclUp") treePath = "DMTreesDumper/ttDM__unclusteredMet__up";
  else if(sys=="metUnclDown")treePath = "DMTreesDumper/ttDM__unclusteredMet__down";

  std::cout<<"Tree: "<<treePath<<endl;
    
  //TChain chainNEvents(treePathNEvents);
  //chainNEvents.Add(path_);

  bench.Start("NEvents");
  TChain chainNEvents(treePathNEvents);
  chainNEvents.AddFileInfoList(fc.GetList());
  Int_t nEvents = (Int_t)chainNEvents.GetEntries();
  bench.Stop("NEvents");
  bench.Print("NEvents");

  bench.Start("NEventsPrePres");


  //  Int_t nEvents = (Int_t)chainNEvents.GetEntries();

  

  TH1D totalWeightTop("w_top_total","Top pt reweighting: overall sample weight",2000,0,2.0);
  chainNEvents.Project("w_top_total","Event_T_Weight","Event_T_Weight!=1.00");
  double topWeight=totalWeightTop.GetMean();
  cout << "totaltopweight is "<< topWeight<<endl;
  if(topWeight==0)topWeight=1;

  //  TChain chain(treePath);
  //  chain.Add(path_);
  TChain chain(treePath);
  chain.AddFileInfoList(fc.GetList());
  Int_t nEventsPrePres = (Int_t)chain.GetEntries();
  std::cout<<"--> --> Number of Events: "<<nEvents<< " after preselection "<< nEventsPrePres << endl;
  bench.Stop("NEventsPrePres");
  bench.Print("NEventsPrePres");

  // Int_t nEventsPrePres = (Int_t)chain.GetEntries();
  //std::cout<<"Number of Events: "<<nEvents<< " after preselection "<< nEventsPrePres << endl;
  //nEvents = 10;
  // Initiating variables 

  //cout << "nevent " << nEvents << " nevent pre-sel " << nEventsPrePres << endl;

  //
  // Set up MVA reader
  //
  // spectator variables, not used for MVA evaluation
  int isSig, b_mis, w_mis, wb_mis;
  float mtop;
  // MVA input variables
  float bdt_qgid1, bdt_qgid2;
  float bdt_dphij1b, bdt_dphij2b, bdt_drj1b, bdt_drj2b;
  float bdt_bjcsv, bdt_jet1csv, bdt_jet2csv;
  float bdt_prob;
  TMVA::Reader res_topmvaReader("");
  res_topmvaReader.AddSpectator("isSig",  &isSig);
  res_topmvaReader.AddSpectator("b_mis",  &b_mis);
  res_topmvaReader.AddSpectator("w_mis",  &w_mis);
  res_topmvaReader.AddSpectator("wb_mis", &wb_mis);
  res_topmvaReader.AddSpectator("mtop",   &mtop);
  res_topmvaReader.AddVariable("qgid1",   &bdt_qgid1);    // QGL for one of the W-jet candidates
  res_topmvaReader.AddVariable("qgid2",   &bdt_qgid2);    // QGL for the other W-jet candidate
  res_topmvaReader.AddVariable("dphij1b", &bdt_dphij1b);  // |deltaPhi| between W-jet #1 and b-jet
  res_topmvaReader.AddVariable("dphij2b", &bdt_dphij2b);  // |deltaPhi| between W-jet #2 and b-jet
  res_topmvaReader.AddVariable("drj1b",   &bdt_drj1b);    // |deltaR| between W-jet #1 and b-jet
  res_topmvaReader.AddVariable("drj2b",   &bdt_drj2b);    // |deltaR| between W-jet #2 and b-jet
  res_topmvaReader.AddVariable("bjcsv",   &bdt_bjcsv);    // CSVv2+IVF value of b-jet
  res_topmvaReader.AddVariable("jet1csv", &bdt_jet1csv);  // CSVv2+IVF value of W-jet #1
  res_topmvaReader.AddVariable("jet2csv", &bdt_jet2csv);  // CSVv2+IVF value of W-jet #2
  res_topmvaReader.AddVariable("prob",    &bdt_prob);     // probability of kinematic fit
  
  const std::string cmssw_base = getenv("CMSSW_BASE");
  const std::string weightsfile = cmssw_base + std::string("/src/ttDM/TopTagResolved/data/toptrainingbits_prob.root_BDTG.weights.xml");
  res_topmvaReader.BookMVA("BDTG", weightsfile.c_str()); 




  
  int sizeMax = 100, sizeMaxTopHad = 100;
  Float_t nTightMuons, nTightElectrons, nVetoElectrons, nLooseMuons, nJets, nCSVJets, nCSVLJets;
  Int_t nResolvedSemiLep, nResolvedHad;
  Int_t muSize, elSize, jetSize;
  //  Float_t evtNumber = 0.;
  Float_t  nType1, nType2;//, nBoostedT2, nBoostedT1; // number of lepton    Heroes of the Storm BlizzCon Brawl Match Game 1 s/jets/tops
  //Float_t evtNumber, runNumber;
  Float_t  mt, mt2w;
  vector<float> mva;
  //vector<int> idxB, idxJ1, idxJ2;
  vector<TopFitResults> topFitRes; 
  vector<TopFitResults> topFitResFull; 
  Float_t  nGoodPV, nPV;

  Float_t numTrueInt;

  Float_t w_pu;
  //for bkg estimation
  //int sizeMaxLept = 4;


  float LHEWeightSign[1] = {1.};
  float w(1.);
  float passElTrig(0.), passMuTrig(0.);
  float passHadTrig(0.) ;
  float bWeightZero = 1.,bWeightZeroBTagUp= 1., bWeightZeroBTagDown=1.0, bWeightZeroMisTagUp=1.0,bWeightZeroMisTagDown=1.0;
  float bWeight2 = 1.,bWeight2BTagUp= 1., bWeight2BTagDown=1.0, bWeight2MisTagUp=1.0,bWeight2MisTagDown=1.0;
  float bWeight12 = 1.,bWeight12BTagUp= 1., bWeight12BTagDown=1.0, bWeight12MisTagUp=1.0,bWeight12MisTagDown=1.0;
  float bWeight12L = 1.,bWeight12LBTagUp= 1., bWeight12LBTagDown=1.0, bWeight12LMisTagUp=1.0,bWeight12LMisTagDown=1.0;
  //float passLepTrig(0.), passHadTrig(0.) ;
  float topSLPt[sizeMax], topSLE[sizeMax], topSLPhi[sizeMax],topSLEta[sizeMax], topSLMT[sizeMax], topSLMass[sizeMax], topSLLeptonFlavour[sizeMax],topSLIndexL[sizeMax],topSLIndexB[sizeMax],metPt[1],metPhi[1],topSLTMPhi[sizeMax],topSLBMPhi[sizeMax], topSLLMPhi[sizeMax], topSLLBMPhi[sizeMax], jetPt[sizeMax], jetEta[sizeMax], jetE[sizeMax], jetCSV[sizeMax], jetIsCSVM[sizeMax],jetIsTight[sizeMax], jetPhi[sizeMax], jetPassID[sizeMax],jetPassDR[sizeMax], muE[sizeMax], muPt[sizeMax], muEta[sizeMax], muIso[sizeMax], muIsTight[sizeMax],muIsLoose[sizeMax], muPhi[sizeMax],elPt[sizeMax], elEta[sizeMax], scEta[sizeMax],  elIso[sizeMax], elIsTight[sizeMax], elIsVeto[sizeMax], elE[sizeMax], elPhi[sizeMax], elPassesDRmu[sizeMax], topHMass[sizeMaxTopHad], topHPt[sizeMaxTopHad], topHE[sizeMaxTopHad],topHEta[sizeMaxTopHad], topHPhi[sizeMaxTopHad], topHBMPhi[sizeMaxTopHad], topHTMPhi[sizeMaxTopHad], topHWMPhi[sizeMaxTopHad] , topHWBPhi[sizeMaxTopHad], topHIndexJ1[sizeMaxTopHad], topHIndexJ2[sizeMaxTopHad], topHIndexB[sizeMaxTopHad], topHWMass[sizeMaxTopHad]; //electronsEta[sizeMax],electronsPhi[sizeMax] ,muonsEta[sizeMax],muonsPhi[sizeMax];
  
  float topHmassDrop[sizeMaxTopHad];
  float jetQGL[sizeMax];
  
  float w_pdfs[nPDF];
  float w_zero,w_q2up,w_q2down,w_top,eventFlavour=0.,NMCLeptons=-1;
  int vhf=0;
  //for bkg estimation
  //float elPt2[sizeMaxLept], elEta2[sizeMaxLept], elE2[sizeMaxLept],  elPhi2[sizeMaxLept],  muPt2[sizeMaxLept], muEta2[sizeMaxLept], muE2[sizeMaxLept],  muPhi2[sizeMaxLept]; 
  


  float n_binlep(0),n_binjet(0), n_binbjet(0), n_binmet160(0), n_binmet320(0), n_binmt2w(0), n_binmt(0), n_phi(0),  n_trig(0), n_noBoost(0), n_boost11(0),n_boost22(0), n_boost12(0), n_boost1Res(0),n_boost2Res(0),n_boostFullRes(0), lep1Pt(0), lep1Eta(0), lep1Phi(0), lep1E(0), lep2Pt(0), lep2Eta(0), lep2Phi(0),  lep2E(0), hadTrigNoiseCleaned_v1(0), hadTrigNoiseCleaned_v2(0), hadTrigJetCleaned_v1(0), hadTrigJetCleaned_v2(0), k_fact(1.), lep1Flavour(0), lep1Charge(0), lep2Flavour(0), lep2Charge(0), runNumber(0.), lumiSec(0.), slTrigEle_v1(0.), slTrigEle_v2(0.), slTrigIsoMu20_v1(0.), slTrigIsoMu20_v2(0.), slTrigIsoMu20_v3(0.), slTrigIsoMuTk20_v1(0.), slTrigIsoMuTk20_v2(0.), slTrigIsoMuTk20_v3(0.), slTrigIsoMuTk20_v4(0.), ZPt(-1.), WPt(-1);
double evtNumber(0.);

  chain.SetBranchAddress("resolvedTopHad_massDrop", topHmassDrop);
  chain.SetBranchAddress("jetsAK4Tight_QGL", jetQGL);
  chain.SetBranchAddress("resolvedTopSemiLep_Pt", topSLPt);
  chain.SetBranchAddress("resolvedTopSemiLep_Phi", topSLPhi);
  chain.SetBranchAddress("resolvedTopSemiLep_Eta", topSLEta);
  chain.SetBranchAddress("resolvedTopSemiLep_E", topSLE);
  chain.SetBranchAddress("resolvedTopSemiLep_Mass", topSLMass);
  chain.SetBranchAddress("resolvedTopSemiLep_MT", topSLMT);
  chain.SetBranchAddress("resolvedTopSemiLep_LBMPhi", topSLLBMPhi);
  chain.SetBranchAddress("resolvedTopSemiLep_BMPhi", topSLBMPhi);// Inverted in the current iteration! Remember for >v8!
  chain.SetBranchAddress("resolvedTopSemiLep_LMPhi", topSLLMPhi);
  chain.SetBranchAddress("resolvedTopSemiLep_TMPhi", topSLTMPhi);

  chain.SetBranchAddress("resolvedTopSemiLep_IndexL", topSLIndexL);
  chain.SetBranchAddress("resolvedTopSemiLep_IndexB", topSLIndexB);
  chain.SetBranchAddress("resolvedTopSemiLep_LeptonFlavour", topSLLeptonFlavour);
  
  chain.SetBranchAddress("resolvedTopHad_IndexB", topHIndexB);
  chain.SetBranchAddress("resolvedTopHad_IndexJ1", topHIndexJ1);
  chain.SetBranchAddress("resolvedTopHad_IndexJ2", topHIndexJ2);

  chain.SetBranchAddress("resolvedTopHad_WMass", topHWMass);
  chain.SetBranchAddress("resolvedTopHad_Mass", topHMass);
  chain.SetBranchAddress("resolvedTopHad_Pt", topHPt);
  chain.SetBranchAddress("resolvedTopHad_E", topHE);
  chain.SetBranchAddress("resolvedTopHad_Phi", topHPhi);
  chain.SetBranchAddress("resolvedTopHad_Eta", topHEta);
  chain.SetBranchAddress("resolvedTopHad_BMPhi", topHBMPhi);
  chain.SetBranchAddress("resolvedTopHad_WMPhi", topHWMPhi);
  chain.SetBranchAddress("resolvedTopHad_TMPhi", topHTMPhi);
  chain.SetBranchAddress("resolvedTopHad_WBPhi", topHWBPhi);  
      
  chain.SetBranchAddress("jetsAK4Tight_IsTight", jetIsTight);

  if(isData=="MC")  chain.SetBranchAddress("jetsAK4Tight_CorrPt", jetPt);
  else  chain.SetBranchAddress("jetsAK4Tight_CorrPt", jetPt);

  chain.SetBranchAddress("jetsAK4Tight_Phi", jetPhi);
  chain.SetBranchAddress("jetsAK4Tight_Eta", jetEta);
  chain.SetBranchAddress("jetsAK4Tight_CSV", jetCSV);

  if(isData=="MC") chain.SetBranchAddress("jetsAK4Tight_CorrE", jetE);
  else  chain.SetBranchAddress("jetsAK4Tight_CorrE", jetE);

  chain.SetBranchAddress("jetsAK4Tight_IsCSVM", jetIsCSVM);
  chain.SetBranchAddress("jetsAK4Tight_PassesID", jetPassID);
  chain.SetBranchAddress("jetsAK4Tight_PassesDR", jetPassDR);
  chain.SetBranchAddress("jetsAK4Tight_size", &jetSize);

  float passHBHE(0.), passHBHEIso(0.), passMETFilters(0.) ;
  chain.SetBranchAddress("Event_passesHBHE", &passHBHE);
  chain.SetBranchAddress("Event_passesHBHEIso", &passHBHEIso);
  chain.SetBranchAddress("Event_passesMETFilters", &passMETFilters);

  chain.SetBranchAddress("electronsTight_E", elE);
  chain.SetBranchAddress("electronsTight_Phi", elPhi);
  chain.SetBranchAddress("electronsTight_Eta", elEta);
  chain.SetBranchAddress("electronsTight_scEta", scEta);
  chain.SetBranchAddress("electronsTight_Pt", elPt);
  chain.SetBranchAddress("electronsTight_Iso03", elIso);
  // chain.SetBranchAddress("electrons_IsEB", elIsEB);
  //  chain.SetBranchAddress("electrons_isTight", elIsTight);
  //  chain.SetBranchAddress("electrons_isVeto", elIsVeto);
  chain.SetBranchAddress("electronsTight_vidTight", elIsTight);
  chain.SetBranchAddress("electronsTight_vidVeto", elIsVeto);
  chain.SetBranchAddress("electronsTight_PassesDRmu", elPassesDRmu);
  chain.SetBranchAddress("electronsTight_size", &elSize);
  chain.SetBranchAddress("muonsTight_E", muE);
  chain.SetBranchAddress("muonsTight_Phi", muPhi);
  chain.SetBranchAddress("muonsTight_Eta", muEta);
  chain.SetBranchAddress("muonsTight_Pt", muPt);
  chain.SetBranchAddress("muonsTight_Iso04", muIso);
  chain.SetBranchAddress("muonsTight_IsTightMuon", muIsTight);
  chain.SetBranchAddress("muonsTight_IsLooseMuon", muIsLoose);
  chain.SetBranchAddress("muonsTight_size", &muSize);
  
  chain.SetBranchAddress("Event_nTightElectrons",&nTightElectrons);
  chain.SetBranchAddress("Event_nVetoElectrons",&nVetoElectrons);
  chain.SetBranchAddress("Event_nLooseMuons",&nLooseMuons);
  chain.SetBranchAddress("Event_nTightMuons",&nTightMuons);
  chain.SetBranchAddress("Event_nJetsCut30",&nJets);
  chain.SetBranchAddress("Event_nCSVMJetsCut30",&nCSVJets);
  chain.SetBranchAddress("Event_nCSVLJetsCut30",&nCSVLJets);
  chain.SetBranchAddress("resolvedTopSemiLep_size",&nResolvedSemiLep);
  chain.SetBranchAddress("resolvedTopHad_size",&nResolvedHad);
  chain.SetBranchAddress("Event_nType1TopJets", &nType1);
  chain.SetBranchAddress("Event_nType2TopJets", &nType2);

  chain.SetBranchAddress("Event_passesHLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v1", &hadTrigNoiseCleaned_v1);
  chain.SetBranchAddress("Event_passesHLT_PFMETNoMu120_NoiseCleaned_PFMHTNoMu120_IDTight_v2", &hadTrigNoiseCleaned_v2);
  chain.SetBranchAddress("Event_passesHLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v1", &hadTrigJetCleaned_v1);
  chain.SetBranchAddress("Event_passesHLT_PFMETNoMu120_JetIdCleaned_PFMHTNoMu120_IDTight_v2", &hadTrigJetCleaned_v2);
  chain.SetBranchAddress("Event_passesHLT_IsoMu20_v1", &slTrigIsoMu20_v1);
  chain.SetBranchAddress("Event_passesHLT_IsoMu20_v2", &slTrigIsoMu20_v2);
  chain.SetBranchAddress("Event_passesHLT_IsoMu20_v3", &slTrigIsoMu20_v3);
  chain.SetBranchAddress("Event_passesHLT_IsoTkMu20_v1", &slTrigIsoMuTk20_v1);
  chain.SetBranchAddress("Event_passesHLT_IsoTkMu20_v2", &slTrigIsoMuTk20_v2);
  chain.SetBranchAddress("Event_passesHLT_IsoTkMu20_v3", &slTrigIsoMuTk20_v3);
  chain.SetBranchAddress("Event_passesHLT_IsoTkMu20_v4", &slTrigIsoMuTk20_v4);

  if(isData=="MC"){
    chain.SetBranchAddress("Event_passesHLT_Ele27_eta2p1_WP75_Gsf_v1", &slTrigEle_v1);
    chain.SetBranchAddress("Event_passesHLT_Ele27_eta2p1_WP75_Gsf_v2", &slTrigEle_v2);
  }
  else{
    chain.SetBranchAddress("Event_passesHLT_Ele27_eta2p1_WPLoose_Gsf_v1", &slTrigEle_v1);
    chain.SetBranchAddress("Event_passesHLT_Ele27_eta2p1_WPLoose_Gsf_v2", &slTrigEle_v2);
  }

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
  chain.SetBranchAddress("Event_RunNumber", &runNumber);
  chain.SetBranchAddress("Event_LumiBlock", &lumiSec);
  chain.SetBranchAddress("Event_EventNumber", &evtNumber);

  chain.SetBranchAddress("Event_Z_Pt", &ZPt);
  chain.SetBranchAddress("Event_W_Pt", &WPt);
 
  const char* Wlabel = "WJets";
  const char* Zlabel = "ZJets";
  const char* DYlabel = "DY";

  bool iswzjets=false;
  double shiftval = 0.0;
  if( !strncmp(sample.c_str(), Wlabel , strlen(Wlabel))) {chain.SetBranchAddress("Event_W_Weight", &k_fact);iswzjets=true;shiftval=0.2;}
  if( !strncmp(sample.c_str(), Zlabel , strlen(Zlabel))) {chain.SetBranchAddress("Event_Z_Weight", &k_fact);iswzjets=true;shiftval=0.2;}
  if( !strncmp(sample.c_str(), DYlabel , strlen(DYlabel))) {chain.SetBranchAddress("Event_Z_Weight", &k_fact);iswzjets=true;shiftval=0.2;}

  std::cout<<"k-factor "<<k_fact<<std::endl;


 // if(channel.find("WJet")!=std::string::npos && channel.find("HT")!=std::string::npos)iswzjets=true;
 // if((channel.find("ZJet")!=std::string::npos && channel.find("HT")!=std::string::npos)iswzjets=true;
 // if(channel.find("DY")!=std::string::npos && channel.find("HT")!=std::string::npos)iswzjets=true;
 


  //  chain.SetBranchAddress("Event_bWeight1_2CSVL", &nType1);
  //  chain.SetBranchAddress("Event_bWeight1_2CSVT", &nType1);

  //  chain.SetBranchAddress("Event_bWeight2CSVL", &nType1);
  //  chain.SetBranchAddress("Event_bWeight2CSVT", &nType1);

  //  float bWeight2 = 1.,bWeight2BTagUp= 1., bWeight2BTagDown=1.0, bWeight2MisTagUp=1.0,bWeight2MisTagDown=1.0;
  //  float bWeight12 = 1.,bWeight12BTagUp= 1., bWeight12BTagDown=1.0, bWeight12MisTagUp=1.0,bWeight12MisTagDown=1.0;

  if(isData=="MC") chain.SetBranchAddress("Event_LHEWeight0", &w_zero);
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
  if(iswzjets && addVHF){
    chain.SetBranchAddress("Event_eventFlavour",&eventFlavour);
  }
  if(addTTSplit) chain.SetBranchAddress("Event_NMCLeptons",&NMCLeptons);
  

  chain.SetBranchAddress("Event_bWeight2CSVM", &bWeight2);

  chain.SetBranchAddress("Event_bWeightMisTagUp2CSVM", &bWeight2MisTagUp);
  chain.SetBranchAddress("Event_bWeightMisTagDown2CSVM", &bWeight2MisTagDown);
  chain.SetBranchAddress("Event_bWeightBTagUp2CSVM", &bWeight2BTagUp);
  chain.SetBranchAddress("Event_bWeightBTagDown2CSVM", &bWeight2BTagDown);

  chain.SetBranchAddress("Event_bWeight0CSVM", &bWeightZero);

  chain.SetBranchAddress("Event_bWeightMisTagUp0CSVM", &bWeightZeroMisTagUp);
  chain.SetBranchAddress("Event_bWeightMisTagDown0CSVM", &bWeightZeroMisTagDown);
  chain.SetBranchAddress("Event_bWeightBTagUp0CSVM", &bWeightZeroBTagUp);
  chain.SetBranchAddress("Event_bWeightBTagDown0CSVM", &bWeightZeroBTagDown);

  chain.SetBranchAddress("Event_bWeight1_2CSVM", &bWeight12);

  chain.SetBranchAddress("Event_bWeightMisTagUp1_2CSVM", &bWeight12MisTagUp);
  chain.SetBranchAddress("Event_bWeightMisTagDown1_2CSVM", &bWeight12MisTagDown);
  chain.SetBranchAddress("Event_bWeightBTagUp1_2CSVM", &bWeight12BTagUp);
  chain.SetBranchAddress("Event_bWeightBTagDown1_2CSVM", &bWeight12BTagDown);

  chain.SetBranchAddress("Event_bWeight1_2CSVL", &bWeight12L);
  chain.SetBranchAddress("Event_bWeightMisTagUp1_2CSVL", &bWeight12LMisTagUp);
  chain.SetBranchAddress("Event_bWeightMisTagDown1_2CSVL", &bWeight12LMisTagDown);
  chain.SetBranchAddress("Event_bWeightBTagUp1_2CSVL", &bWeight12LBTagUp);
  chain.SetBranchAddress("Event_bWeightBTagDown1_2CSVL", &bWeight12LBTagDown);

  chain.SetBranchAddress("Event_mt",&mt);
  chain.SetBranchAddress("Event_Mt2w",&mt2w);

  if(isData=="MC"){
    chain.SetBranchAddress("metFull_CorrT1Pt",metPt);  
    chain.SetBranchAddress("metFull_CorrT1Phi",metPhi);}
  else {
    chain.SetBranchAddress("metFull_CorrT1Pt",metPt);
    chain.SetBranchAddress("metFull_CorrT1Phi",metPhi);}
  
  if(sample =="SingleTop_T_tchan" || sample == "SingleTop_Tbar_tchan" || sample == "WJets" || sample == "DY")chain.SetBranchAddress("Event_LHEWeightSign", LHEWeightSign);

  //cout << "test " << endl;

  ////chain.SetBranchAddress("Event_passesLeptonicTriggers", &passLepTrig);
  chain.SetBranchAddress("Event_passesSingleElTriggers", &passElTrig);
  chain.SetBranchAddress("Event_passesSingleMuTriggers", &passMuTrig);
  chain.SetBranchAddress("Event_passesHadronicTriggers", &passHadTrig);
  
  chain.SetBranchAddress("Event_nGoodPV",&nGoodPV);
  chain.SetBranchAddress("Event_nPV",&nPV);

  chain.SetBranchAddress("Event_nTruePV",&numTrueInt);
  

  bool elePD(0), muPD(0);
  
  const char* Elelabel = "SingleEl";
  const char* Mulabel = "SingleMu";
  
  if( !strncmp(sample.c_str(), Elelabel , strlen(Elelabel))) elePD = 1;
  if( !strncmp(sample.c_str(), Mulabel , strlen(Mulabel))) muPD = 1;

  std::cout<<"Is electron PD? "<<elePD<<std::endl;
  std::cout<<"Is muon PD? "<<muPD<<std::endl;


  /********************************************************************************/
  /**************                    Histogram booking              ***************/
  /********************************************************************************/


  
  //After trigger
  TH1F *h_nPV   [maxSysts] ;systZero.initHistogramsSysts(h_nPV,"nPV","nPV",500,0,500);
  TH1F *h_nGoodPV   [maxSysts] ;systZero.initHistogramsSysts(h_nGoodPV,"n_GoodPV","nGoodPV",500,0,500);
  TH1F *h_nPV_w   [maxSysts] ;systZero.initHistogramsSysts(h_nPV_w,"nPV_w","nPV_w",500,0,500);
  TH1F *h_nGoodPV_w   [maxSysts] ;systZero.initHistogramsSysts(h_nGoodPV_w,"n_GoodPV_w","nGoodPV_w",500,0,500);
  TH1F *h_nJets[maxSysts] ;systZero.initHistogramsSysts(h_nJets,"h_nJets","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets[maxSysts] ;systZero.initHistogramsSysts(h_nbJets,"h_nbJets","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_jet1Pt[maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt,"h_jet1Pt","Leading jet Pt distribution",100,0,500);
  TH1F *h_jet2Pt[maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt,"h_jet2Pt","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_jet3Pt[maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt,"h_jet3Pt","Third Jet Pt distribution",100,0,500);

  TH1F *h_bjet1Pt[maxSysts] ;systZero.initHistogramsSysts(h_bjet1Pt,"h_bjet1Pt","Leading jet Pt distribution",100,0,500);
  TH1F *h_bjet2Pt[maxSysts] ;systZero.initHistogramsSysts(h_bjet2Pt,"h_bjet2Pt","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_bjet3Pt[maxSysts] ;systZero.initHistogramsSysts(h_bjet3Pt,"h_bjet3Pt","Third Jet Pt distribution",100,0,500);

  TH1F *h_bjetsPt[maxSysts] ;systZero.initHistogramsSysts(h_bjetsPt,"h_bjetsPt","B-Jets Pt distribution",100,0,500);

  //N-1 
  TH1F *h_dphi[maxSysts] ;systZero.initHistogramsSysts(h_dphi,"dphi","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muondphi[maxSysts] ;systZero.initHistogramsSysts(h_muondphi,"muondphi","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electrondphi[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi,"electrondphi","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_dphi_6j[maxSysts] ;systZero.initHistogramsSysts(h_dphi_6j,"dphi_6j","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muondphi_6j[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_6j,"muondphi_6j","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electrondphi_6j[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_6j,"electrondphi_6j","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);

  TH1F *h_mt[maxSysts] ;systZero.initHistogramsSysts(h_mt,"mt","M_{T}",100,0,500);
  TH1F *h_muonmt[maxSysts] ;systZero.initHistogramsSysts(h_muonmt,"muonmt","M_{T}",100,0,500);
  TH1F *h_electronmt[maxSysts] ;systZero.initHistogramsSysts(h_electronmt,"electronmt","M_{T}",100,0,500);

  TH1F *h_mt2w[maxSysts] ;systZero.initHistogramsSysts(h_mt2w,"mt2w","MT2W",100,50,500);
  TH1F *h_muonmt2w[maxSysts] ;systZero.initHistogramsSysts(h_muonmt2w,"muonmt2w","MT2W",100,50,500);
  TH1F *h_electronmt2w[maxSysts] ;systZero.initHistogramsSysts(h_electronmt2w,"electronmt2w","MT2W",100,50,500);

  TH1F *h_met[maxSysts] ;systZero.initHistogramsSysts(h_met,"met","MET",200,0,2000);
  TH1F *h_muonmet[maxSysts] ;systZero.initHistogramsSysts(h_muonmet,"muonmet","MET",200,0,2000);
  TH1F *h_electronmet[maxSysts] ;systZero.initHistogramsSysts(h_electronmet,"electronmet","MET",200,0,2000);

  //SR shape analysis
  TH1F *h_metFinal[maxSysts] ;systZero.initHistogramsSysts(h_metFinal,"metFinal","MET",200,0,2000);

  TH1F *h_metFinal_4j[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_4j,"metFinal_4j","MET",200,0,2000);
  TH1F *h_metFinal_6j[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_6j,"metFinal_6j","MET",200,0,2000);

  TH1F *h_metFinal_noMT2W[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_noMT2W,"metFinal_noMT2W","MET",200,0,2000);
  TH1F *h_muonmetFinal[maxSysts] ;systZero.initHistogramsSysts(h_muonmetFinal,"muonmetFinal","MET",200,0,2000);
  TH1F *h_electronmetFinal[maxSysts] ;systZero.initHistogramsSysts(h_electronmetFinal,"electronmetFinal","MET",200,0,2000);
  
  TH1F *h_nJets_fin[maxSysts] ;systZero.initHistogramsSysts(h_nJets_fin,"h_nJets_fin","Number of tight jets",12,-0.5,11.5);
  TH1F *h_nbJets_fin[maxSysts] ;systZero.initHistogramsSysts(h_nbJets_fin,"h_nbJets_fin","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_jet1Pt_fin[maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_fin,"h_jet1Pt_fin","Leading jet Pt distribution",100,0,500);
  TH1F *h_jet2Pt_fin[maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_fin,"h_jet2Pt_fin","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_jet3Pt_fin[maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_fin,"h_jet3Pt_fin","Third Jet Pt distribution",100,0,500);

  //pre-selection
  TH1F *h_jetQGL1_preS  [maxSysts] ;systZero.initHistogramsSysts(h_jetQGL1_preS,"h_jetQGL1_preS","QGL 1st jet",50,0,1);
  TH1F *h_jetQGL2_preS  [maxSysts] ;systZero.initHistogramsSysts(h_jetQGL2_preS,"h_jetQGL2_preS","QGL 2nd jet",50,0,1);
  
  TH1F *h_nPV_preS   [maxSysts] ;systZero.initHistogramsSysts(h_nPV_preS,"nPV_preS","nPV",500,0,500);
  TH1F *h_nGoodPV_preS [maxSysts] ;systZero.initHistogramsSysts(h_nGoodPV_preS,"n_GoodPV_preS","nGoodPV",500,0,500);
  TH1F *h_nPV_w_preS   [maxSysts] ;systZero.initHistogramsSysts(h_nPV_w_preS,"nPV_w_preS","nPV_w",500,0,500);
  TH1F *h_nGoodPV_w_preS [maxSysts] ;systZero.initHistogramsSysts(h_nGoodPV_w_preS,"n_GoodPV_w_preS","nGoodPV_w",500,0,500);
  
  TH1F *h_met_preS[maxSysts] ;systZero.initHistogramsSysts(h_met_preS,"met_preS","MET",200,0,2000);
  TH1F *h_dphi_preS[maxSysts] ;systZero.initHistogramsSysts(h_dphi_preS,"dphi_preS","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_dphi_6j_preS[maxSysts] ;systZero.initHistogramsSysts(h_dphi_6j_preS,"dphi_6j_preS","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  
  TH1F *h_mt_preS  [maxSysts] ;systZero.initHistogramsSysts(h_mt_preS,"mt_preaS","M_{T}",500,0,500);
  TH1F *h_mt_hadpreS  [maxSysts] ;systZero.initHistogramsSysts(h_mt_hadpreS,"mt_hadpreS","M_{T}",500,0,500);
  TH1F *h_mt2w_preS[maxSysts] ;systZero.initHistogramsSysts(h_mt2w_preS,"mt2w_preS","MT2W",45,50,500);
  
  TH1F *h_nJets_preS [maxSysts] ;systZero.initHistogramsSysts(h_nJets_preS,"h_nJets_preS","Number of tight jets",12,-0.5,11.5);
  TH1F *h_nbJets_preS[maxSysts] ;systZero.initHistogramsSysts(h_nbJets_preS,"h_nbJets_preS","Number of tight b-jets",11,-0.5,10.5);
  
  TH1F *h_jet1Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_preS,"h_jet1Pt_preS","Leading jet Pt distribution",500,0,500);
  TH1F *h_jet2Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_preS,"h_jet2Pt_preS","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_jet3Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_preS,"h_jet3Pt_preS","Third Jet Pt distribution",500,0,500);
  
  TH1F *h_jet1Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_jet1Eta_preS,"h_jet1Eta_preS","Leading jet Eta distribution",64,-3.2,3.2);
  TH1F *h_jet2Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_jet2Eta_preS,"h_jet2Eta_preS","Trailing Jet Eta distribution",64,-3.2,3.2);
  TH1F *h_jet3Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_jet3Eta_preS,"h_jet3Eta_preS","Third Jet Eta distribution",64,-3.2,3.2);
  
  TH1F *h_jet1Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_jet1Phi_preS,"h_jet1Phi_preS","Leading jet Phi distribution",64,-3.2,3.2);
  TH1F *h_jet2Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_jet2Phi_preS,"h_jet2Phi_preS","Trailing Jet Phi distribution",64,-3.2,3.2);
  TH1F *h_jet3Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_jet3Phi_preS,"h_jet3Phi_preS","Third Jet Phi distribution",64,-3.2,3.2);
  
  TH1F *h_jet1CSV_preS[maxSysts] ;systZero.initHistogramsSysts(h_jet1CSV_preS,"h_jet1CSV_preS","Leading jet CSV distribution",100,0,1);
  TH1F *h_jet2CSV_preS[maxSysts] ;systZero.initHistogramsSysts(h_jet2CSV_preS,"h_jet2CSV_preS","Trailing Jet CSV distribution",100,0,1);
  TH1F *h_jet3CSV_preS[maxSysts] ;systZero.initHistogramsSysts(h_jet3CSV_preS,"h_jet3CSV_preS","Third Jet CSV distribution",100,0,1);
  
  TH1F *h_bjet1Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_bjet1Pt_preS,"h_bjet1Pt_preS","Leading bjet Pt distribution",500,0,500);
  TH1F *h_bjet2Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_bjet2Pt_preS,"h_bjet2Pt_preS","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_bjet3Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_bjet3Pt_preS,"h_bjet3Pt_preS","Third Bjet Pt distribution",500,0,500);
  
  TH1F *h_bjet1Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_bjet1Eta_preS,"h_bjet1Eta_preS","Leading bjet Eta distribution",64,-3.2,3.2);
  TH1F *h_bjet2Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_bjet2Eta_preS,"h_bjet2Eta_preS","Trailing Bjet Eta distribution",64,-3.2,3.2);
  TH1F *h_bjet3Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_bjet3Eta_preS,"h_bjet3Eta_preS","Third Bjet Eta distribution",64,-3.2,3.2);
  
  TH1F *h_bjet1Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_bjet1Phi_preS,"h_bjet1Phi_preS","Leading bjet Phi distribution",64,-3.2,3.2);
  TH1F *h_bjet2Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_bjet2Phi_preS,"h_bjet2Phi_preS","Trailing Bjet Phi distribution",64,-3.2,3.2);
  TH1F *h_bjet3Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_bjet3Phi_preS,"h_bjet3Phi_preS","Third Bjet Phi distribution",64,-3.2,3.2);
  
  
  TH1F *h_lep1Pt_preS [maxSysts] ;systZero.initHistogramsSysts(h_lep1Pt_preS,"h_lep1Pt_preS","Lepton pt distribution",400, 0.,200);
  TH1F *h_lep1Eta_preS [maxSysts] ;systZero.initHistogramsSysts(h_lep1Eta_preS,"h_lep1Eta_preS","Lepton pt distribution",48, -2.4,2.4);
  TH1F *h_lep1Phi_preS [maxSysts] ;systZero.initHistogramsSysts(h_lep1Phi_preS,"h_lep1Phi_preS","Lepton pt distribution",64, -3.2,3.2);
  TH1F *h_lep1E_preS [maxSysts] ;systZero.initHistogramsSysts(h_lep1E_preS,"h_lep1E_preS","Lepton pt distribution",800, 0.,400);
  
  TH1F *h_mu1Pt_preS [maxSysts] ;systZero.initHistogramsSysts(h_mu1Pt_preS,"h_mu1Pt_preS","Lepton pt distribution",400, 0.,200);
  TH1F *h_mu1Eta_preS [maxSysts] ;systZero.initHistogramsSysts(h_mu1Eta_preS,"h_mu1Eta_preS","Lepton pt distribution",48, -2.4,2.4);
  TH1F *h_mu1Phi_preS [maxSysts] ;systZero.initHistogramsSysts(h_mu1Phi_preS,"h_mu1Phi_preS","Lepton pt distribution",64, -3.2,3.2);
  TH1F *h_mu1E_preS [maxSysts] ;systZero.initHistogramsSysts(h_mu1E_preS,"h_mu1E_preS","Lepton pt distribution",800, 0.,400);
  
  TH1F *h_el1Pt_preS [maxSysts] ;systZero.initHistogramsSysts(h_el1Pt_preS,"h_el1Pt_preS","Lepton pt distribution",400, 0.,200);
  TH1F *h_el1Eta_preS [maxSysts] ;systZero.initHistogramsSysts(h_el1Eta_preS,"h_el1Eta_preS","Lepton pt distribution",48, -2.4,2.4);
  TH1F *h_el1Phi_preS [maxSysts] ;systZero.initHistogramsSysts(h_el1Phi_preS,"h_el1Phi_preS","Lepton pt distribution",64, -3.2,3.2);
  TH1F *h_el1E_preS [maxSysts] ;systZero.initHistogramsSysts(h_el1E_preS,"h_el1E_preS","Lepton pt distribution",800, 0.,400);
  
  TH1F *h_muonmet_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonmet_preS,"muonmet_preS","MET",200,0,2000);
  TH1F *h_muondphi_preS[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_preS,"muondphi_preS","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muondphi_6j_preS[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_6j_preS,"muondphi_6j_preS","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muonmt_preS  [maxSysts] ;systZero.initHistogramsSysts(h_muonmt_preS,"muonmt_preS","M_{T}",100,0,500);
  TH1F *h_muonmt2w_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonmt2w_preS,"muonmt2w_preS","MT2W",100,50,500);
  
  TH1F *h_muonnJets_preS  [maxSysts] ;systZero.initHistogramsSysts(h_muonnJets_preS,"h_muonnJets_preS","Number of tight jets",12,-0.5,11.5);
  TH1F *h_muonnbJets_preS [maxSysts] ;systZero.initHistogramsSysts(h_muonnbJets_preS,"h_muonnbJets_preS","Number of tight b-jets",11,-0.5,10.5);
  
  TH1F *h_muonjet1Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonjet1Pt_preS,"h_muonjet1Pt_preS","Leading jet Pt distribution",500,0,500);
  TH1F *h_muonjet2Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonjet2Pt_preS,"h_muonjet2Pt_preS","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_muonjet3Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonjet3Pt_preS,"h_muonjet3Pt_preS","Third Jet Pt distribution",500,0,500);
  
  TH1F *h_muonjet1Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonjet1Eta_preS,"h_muonjet1Eta_preS","Leading jet Eta distribution",64,-3.2,3.2);
  TH1F *h_muonjet2Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonjet2Eta_preS,"h_muonjet2Eta_preS","Trailing Jet Eta distribution",64,-3.2,3.2);
  TH1F *h_muonjet3Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonjet3Eta_preS,"h_muonjet3Eta_preS","Third Jet Eta distribution",64,-3.2,3.2);
  
  TH1F *h_muonjet1Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonjet1Phi_preS,"h_muonjet1Phi_preS","Leading jet Phi distribution",64,-3.2,3.2);
  TH1F *h_muonjet2Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonjet2Phi_preS,"h_muonjet2Phi_preS","Trailing Jet Phi distribution",64,-3.2,3.2);
  TH1F *h_muonjet3Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonjet3Phi_preS,"h_muonjet3Phi_preS","Third Jet Phi distribution",64,-3.2,3.2);
  
  TH1F *h_muonjet1CSV_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonjet1CSV_preS,"h_muonjet1CSV_preS","Leading jet CSV distribution",100,0,1);
  TH1F *h_muonjet2CSV_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonjet2CSV_preS,"h_muonjet2CSV_preS","Trailing Jet CSV distribution",100,0,1);
  TH1F *h_muonjet3CSV_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonjet3CSV_preS,"h_muonjet3CSV_preS","Third Jet CSV distribution",100,0,1);
  
  TH1F *h_muonbjet1Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet1Pt_preS,"h_muonbjet1Pt_preS","Leading bjet Pt distribution",500,0,500);
  TH1F *h_muonbjet2Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet2Pt_preS,"h_muonbjet2Pt_preS","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_muonbjet3Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet3Pt_preS,"h_muonbjet3Pt_preS","Third Bjet Pt distribution",500,0,500);
  
  TH1F *h_muonbjet1Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet1Eta_preS,"h_muonbjet1Eta_preS","Leading bjet Eta distribution",64,-3.2,3.2);
  TH1F *h_muonbjet2Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet2Eta_preS,"h_muonbjet2Eta_preS","Trailing Bjet Eta distribution",64,-3.2,3.2);
  TH1F *h_muonbjet3Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet3Eta_preS,"h_muonbjet3Eta_preS","Third Bjet Eta distribution",64,-3.2,3.2);
  
  TH1F *h_muonbjet1Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet1Phi_preS,"h_muonbjet1Phi_preS","Leading bjet Phi distribution",64,-3.2,3.2);
  TH1F *h_muonbjet2Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet2Phi_preS,"h_muonbjet2Phi_preS","Trailing Bjet Phi distribution",64,-3.2,3.2);
  TH1F *h_muonbjet3Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet3Phi_preS,"h_muonbjet3Phi_preS","Third Bjet Phi distribution",64,-3.2,3.2);
  
  TH1F *h_electronmet_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronmet_preS,"electronmet_preS","MET",200,0,2000);
  TH1F *h_electrondphi_preS[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_preS,"electrondphi_preS","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electrondphi_6j_preS[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_6j_preS,"electrondphi_6j_preS","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  
  TH1F *h_electronmt_preS  [maxSysts] ;systZero.initHistogramsSysts(h_electronmt_preS,"electronmt_preS","M_{T}",100,0,500);
  TH1F *h_mt_EC_preS  [maxSysts] ;systZero.initHistogramsSysts(h_mt_EC_preS,"mt_EC_preS","M_{T}",500,0,500);
  TH1F *h_mt_BC_preS  [maxSysts] ;systZero.initHistogramsSysts(h_mt_BC_preS,"mt_BC_preS","M_{T}",500,0,500);
  TH1F *h_mt_EC_hadpreS  [maxSysts] ;systZero.initHistogramsSysts(h_mt_EC_hadpreS,"mt_EC_hadpreS","M_{T}",500,0,500);
  TH1F *h_mt_BC_hadpreS  [maxSysts] ;systZero.initHistogramsSysts(h_mt_BC_hadpreS,"mt_BC_hadpreS","M_{T}",500,0,500);
  
  TH1F *h_electronmt2w_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronmt2w_preS,"electronmt2w_preS","MT2W",100,50,500);
  
  TH1F *h_electronnJets_preS  [maxSysts] ;systZero.initHistogramsSysts(h_electronnJets_preS,"h_electronnJets_preS","Number of tight jets",12,-0.5,11.5);
  TH1F *h_electronnbJets_preS [maxSysts] ;systZero.initHistogramsSysts(h_electronnbJets_preS,"h_electronnbJets_preS","Number of tight b-jets",11,-0.5,10.5);
  
  TH1F *h_electronjet1Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronjet1Pt_preS,"h_electronjet1Pt_preS","Leading jet Pt distribution",500,0,500);
  TH1F *h_electronjet2Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronjet2Pt_preS,"h_electronjet2Pt_preS","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_electronjet3Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronjet3Pt_preS,"h_electronjet3Pt_preS","Third Jet Pt distribution",500,0,500);
  
  TH1F *h_electronjet1Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronjet1Eta_preS,"h_electronjet1Eta_preS","Leading jet Eta distribution",64,-3.2,3.2);
  TH1F *h_electronjet2Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronjet2Eta_preS,"h_electronjet2Eta_preS","Trailing Jet Eta distribution",64,-3.2,3.2);
  TH1F *h_electronjet3Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronjet3Eta_preS,"h_electronjet3Eta_preS","Third Jet Eta distribution",64,-3.2,3.2);
  
  TH1F *h_electronjet1Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronjet1Phi_preS,"h_electronjet1Phi_preS","Leading jet Phi distribution",64,-3.2,3.2);
  TH1F *h_electronjet2Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronjet2Phi_preS,"h_electronjet2Phi_preS","Trailing Jet Phi distribution",64,-3.2,3.2);
  TH1F *h_electronjet3Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronjet3Phi_preS,"h_electronjet3Phi_preS","Third Jet Phi distribution",64,-3.2,3.2);
  
  TH1F *h_electronjet1CSV_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronjet1CSV_preS,"h_electronjet1CSV_preS","Leading jet CSV distribution",100,0,1);
  TH1F *h_electronjet2CSV_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronjet2CSV_preS,"h_electronjet2CSV_preS","Trailing Jet CSV distribution",100,0,1);
  TH1F *h_electronjet3CSV_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronjet3CSV_preS,"h_electronjet3CSV_preS","Third Jet CSV distribution",100,0,1);
  
  TH1F *h_electronbjet1Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet1Pt_preS,"h_electronbjet1Pt_preS","Leading bjet Pt distribution",500,0,500);
  TH1F *h_electronbjet2Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet2Pt_preS,"h_electronbjet2Pt_preS","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_electronbjet3Pt_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet3Pt_preS,"h_electronbjet3Pt_preS","Third Bjet Pt distribution",500,0,500);
  
  TH1F *h_electronbjet1Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet1Eta_preS,"h_electronbjet1Eta_preS","Leading bjet Eta distribution",64,-3.2,3.2);
  TH1F *h_electronbjet2Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet2Eta_preS,"h_electronbjet2Eta_preS","Trailing Bjet Eta distribution",64,-3.2,3.2);
  TH1F *h_electronbjet3Eta_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet3Eta_preS,"h_electronbjet3Eta_preS","Third Bjet Eta distribution",64,-3.2,3.2);
  
  TH1F *h_electronbjet1Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet1Phi_preS,"h_electronbjet1Phi_preS","Leading bjet Phi distribution",64,-3.2,3.2);
  TH1F *h_electronbjet2Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet2Phi_preS,"h_electronbjet2Phi_preS","Trailing Bjet Phi distribution",64,-3.2,3.2);
  TH1F *h_electronbjet3Phi_preS[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet3Phi_preS,"h_electronbjet3Phi_preS","Third Bjet Phi distribution",64,-3.2,3.2);  

  //CR tt(2lep) (CR1)
  TH1F *h_metFinal_2lep   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_2lep,"metFinal_2lep","MET",200,0,2000);

  TH1F *h_metFinal_2lep_4j   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_2lep_4j,"metFinal_2lep_4j","MET",200,0,2000);
  TH1F *h_metFinal_2lep_6j   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_2lep_6j,"metFinal_2lep_6j","MET",200,0,2000);

  //TH1F *h_metFinal_2lep_Full   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_2lep_Full,"metFinal_2lep_Full","MET",200,0,2000);
  TH1F *h_electronmetFinal_2lep   [maxSysts] ;systZero.initHistogramsSysts(h_electronmetFinal_2lep,"electronmetFinal_2lep","MET",200,0,2000);
  TH1F *h_muonmetFinal_2lep   [maxSysts] ;systZero.initHistogramsSysts(h_muonmetFinal_2lep,"muonmetFinal_2lep","MET",200,0,2000);
  TH1F *h_mixmetFinal_2lep   [maxSysts] ;systZero.initHistogramsSysts(h_mixmetFinal_2lep,"mixmetFinal_2lep","MET",200,0,2000);
  
  TH1F *h_jet1Pt_CRtt2_lep   [maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_CRtt2_lep,"h_jet1Pt_CRtt2_lep","Leading jet Pt distribution",100,0,500);
  TH1F *h_jet2Pt_CRtt2_lep   [maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_CRtt2_lep,"h_jet2Pt_CRtt2_lep","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_jet3Pt_CRtt2_lep   [maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_CRtt2_lep,"h_jet3Pt_CRtt2_lep","Third Jet Pt distribution",100,0,500);
  TH1F *h_nJets_CRtt2_lep   [maxSysts] ;systZero.initHistogramsSysts(h_nJets_CRtt2_lep,"h_nJets_CRtt2_lep","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets_CRtt2_lep   [maxSysts] ;systZero.initHistogramsSysts(h_nbJets_CRtt2_lep,"h_nbJets_CRtt2_lep","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_dphi_CRtt2_lep   [maxSysts] ;systZero.initHistogramsSysts(h_dphi_CRtt2_lep,"h_dphi_CRtt2_lep","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_lepPt_CRtt2_lep   [maxSysts] ;systZero.initHistogramsSysts(h_lepPt_CRtt2_lep,"h_lepPt_CRtt2_lep","Lepton pt distribution",200,0,100);
  //TH1F *h_topSLMass_CRtt2_lep   [maxSysts] ;systZero.initHistogramsSysts(h_topSLMass_CRtt2_lep,"h_topSLMass_CRtt2_lep","Top Mass distribution",100,0,500);
  TH1F *h_mt_CR1  [maxSysts] ;systZero.initHistogramsSysts(h_mt_CR1,"mt_CR1","M_{T}",100,0,500);
  TH1F *h_mt2w_CR1[maxSysts] ;systZero.initHistogramsSysts(h_mt2w_CR1,"mt2w_CR1","MT2W",100,50,500);

  TH1F *h_lep1Pt_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_lep1Pt_CR1,"h_lep1Pt_CR1","Lepton pt distribution",400, 0.,200);
  TH1F *h_lep1Eta_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_lep1Eta_CR1,"h_lep1Eta_CR1","Lepton pt distribution",48, -2.4,2.4);
  TH1F *h_lep1Phi_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_lep1Phi_CR1,"h_lep1Phi_CR1","Lepton pt distribution",64, -3.2,3.2);
  TH1F *h_lep1E_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_lep1E_CR1,"h_lep1E_CR1","Lepton pt distribution",800, 0.,400);
  
  TH1F *h_lep2Pt_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_lep2Pt_CR1,"h_lep2Pt_CR1","Lepton pt distribution",400, 0.,200);
  TH1F *h_lep2Eta_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_lep2Eta_CR1,"h_lep2Eta_CR1","Lepton pt distribution",48, -2.4,2.4);
  TH1F *h_lep2Phi_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_lep2Phi_CR1,"h_lep2Phi_CR1","Lepton pt distribution",64, -3.2,3.2);
  TH1F *h_lep2E_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_lep2E_CR1,"h_lep2E_CR1","Lepton pt distribution",800, 0.,400);
  
  TH1F *h_mu1Pt_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_mu1Pt_CR1,"h_mu1Pt_CR1","Lepton pt distribution",400, 0.,200);
  TH1F *h_mu2Pt_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_mu2Pt_CR1,"h_mu2Pt_CR1","Lepton pt distribution",400, 0.,200);
  
  TH1F *h_el1Pt_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_el1Pt_CR1,"h_el1Pt_CR1","Lepton pt distribution",400, 0.,200);
  TH1F *h_el2Pt_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_el2Pt_CR1,"h_el2Pt_CR1","Lepton pt distribution",400, 0.,200);
  
  TH1F *h_muondphi_CR1[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_CR1,"muondphi_CR1","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muondphi_6j_CR1[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_6j_CR1,"muondphi_6j_CR1","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muonmt_CR1  [maxSysts] ;systZero.initHistogramsSysts(h_muonmt_CR1,"muonmt_CR1","M_{T}",100,0,500);
  TH1F *h_muonmt2w_CR1[maxSysts] ;systZero.initHistogramsSysts(h_muonmt2w_CR1,"muonmt2w_CR1","MT2W",100,50,500);
  TH1F *h_muonnJets_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_muonnJets_CR1,"h_muonnJets_CR1","Number of tight jets",12,-0.5,11.5);
  TH1F *h_muonnbJets_CR1[maxSysts] ;systZero.initHistogramsSysts(h_muonnbJets_CR1,"h_muonnbJets_CR1","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_muonjet1Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_muonjet1Pt_CR1,"h_muonjet1Pt_CR1","Leading jet Pt distribution",500,0,500);
  TH1F *h_muonjet2Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_muonjet2Pt_CR1,"h_muonjet2Pt_CR1","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_muonjet3Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_muonjet3Pt_CR1,"h_muonjet3Pt_CR1","Third Jet Pt distribution",500,0,500);
  TH1F *h_muonbjet1Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet1Pt_CR1,"h_muonbjet1Pt_CR1","Leading bjet Pt distribution",500,0,500);
  TH1F *h_muonbjet2Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet2Pt_CR1,"h_muonbjet2Pt_CR1","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_muonbjet3Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet3Pt_CR1,"h_muonbjet3Pt_CR1","Third Bjet Pt distribution",500,0,500);
  
  TH1F *h_electrondphi_CR1[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_CR1,"electrondphi_CR1","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electrondphi_6j_CR1[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_6j_CR1,"electrondphi_6j_CR1","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electronmt_CR1  [maxSysts] ;systZero.initHistogramsSysts(h_electronmt_CR1,"electronmt_CR1","M_{T}",100,0,500);
  TH1F *h_electronmt2w_CR1[maxSysts] ;systZero.initHistogramsSysts(h_electronmt2w_CR1,"electronmt2w_CR1","MT2W",100,50,500);
  TH1F *h_electronnJets_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_electronnJets_CR1,"h_electronnJets_CR1","Number of tight jets",12,-0.5,11.5);
  TH1F *h_electronnbJets_CR1[maxSysts] ;systZero.initHistogramsSysts(h_electronnbJets_CR1,"h_electronnbJets_CR1","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_electronjet1Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_electronjet1Pt_CR1,"h_electronjet1Pt_CR1","Leading jet Pt distribution",500,0,500);
  TH1F *h_electronjet2Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_electronjet2Pt_CR1,"h_electronjet2Pt_CR1","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_electronjet3Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_electronjet3Pt_CR1,"h_electronjet3Pt_CR1","Third Jet Pt distribution",500,0,500);
  TH1F *h_electronbjet1Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet1Pt_CR1,"h_electronbjet1Pt_CR1","Leading bjet Pt distribution",500,0,500);
  TH1F *h_electronbjet2Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet2Pt_CR1,"h_electronbjet2Pt_CR1","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_electronbjet3Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet3Pt_CR1,"h_electronbjet3Pt_CR1","Third Bjet Pt distribution",500,0,500);
  
  TH1F *h_mixdphi_CR1[maxSysts] ;systZero.initHistogramsSysts(h_mixdphi_CR1,"mixdphi_CR1","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_mixdphi_6j_CR1[maxSysts] ;systZero.initHistogramsSysts(h_mixdphi_6j_CR1,"mixdphi_6j_CR1","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_mixmt_CR1  [maxSysts] ;systZero.initHistogramsSysts(h_mixmt_CR1,"mixmt_CR1","M_{T}",100,0,500);
  TH1F *h_mixmt2w_CR1[maxSysts] ;systZero.initHistogramsSysts(h_mixmt2w_CR1,"mixmt2w_CR1","MT2W",100,50,500);
  TH1F *h_mixnJets_CR1 [maxSysts] ;systZero.initHistogramsSysts(h_mixnJets_CR1,"h_mixnJets_CR1","Number of tight jets",12,-0.5,11.5);
  TH1F *h_mixnbJets_CR1[maxSysts] ;systZero.initHistogramsSysts(h_mixnbJets_CR1,"h_mixnbJets_CR1","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_mixjet1Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_mixjet1Pt_CR1,"h_mixjet1Pt_CR1","Leading jet Pt distribution",500,0,500);
  TH1F *h_mixjet2Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_mixjet2Pt_CR1,"h_mixjet2Pt_CR1","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_mixjet3Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_mixjet3Pt_CR1,"h_mixjet3Pt_CR1","Third Jet Pt distribution",500,0,500);
  TH1F *h_mixbjet1Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_mixbjet1Pt_CR1,"h_mixbjet1Pt_CR1","Leading bjet Pt distribution",500,0,500);
  TH1F *h_mixbjet2Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_mixbjet2Pt_CR1,"h_mixbjet2Pt_CR1","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_mixbjet3Pt_CR1[maxSysts] ;systZero.initHistogramsSysts(h_mixbjet3Pt_CR1,"h_mixbjet3Pt_CR1","Third Bjet Pt distribution",500,0,500);
  
  //with additonal Z mass window
  TH1F *h_metFinal_2lep_Z_nobtag   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_2lep_Z_nobtag,"metFinal_2lep_Z_nobtag","MET",200,0,2000);
  TH1F *h_muonmetFinal_2lep_Z_nobtag   [maxSysts] ;systZero.initHistogramsSysts(h_muonmetFinal_2lep_Z_nobtag,"muonmetFinal_2lep_Z_nobtag","MET",200,0,2000);
  TH1F *h_electronmetFinal_2lep_Z_nobtag   [maxSysts] ;systZero.initHistogramsSysts(h_electronmetFinal_2lep_Z_nobtag,"electronmetFinal_2lep_Z_nobtag","MET",200,0,2000);
  TH1F *h_mixmetFinal_2lep_Z_nobtag   [maxSysts] ;systZero.initHistogramsSysts(h_mixmetFinal_2lep_Z_nobtag,"mixmetFinal_2lep_Z_nobtag","MET",200,0,2000);
  
  TH1F *h_jet1Pt_CRZ_lep   [maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_CRZ_lep,"h_jet1Pt_CRZ_lep","Leading jet Pt distribution",100,0,500);
  TH1F *h_jet2Pt_CRZ_lep   [maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_CRZ_lep,"h_jet2Pt_CRZ_lep","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_jet3Pt_CRZ_lep   [maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_CRZ_lep,"h_jet3Pt_CRZ_lep","Third Jet Pt distribution",100,0,500);
  TH1F *h_nJets_CRZ_lep   [maxSysts] ;systZero.initHistogramsSysts(h_nJets_CRZ_lep,"h_nJets_CRZ_lep","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets_CRZ_lep   [maxSysts] ;systZero.initHistogramsSysts(h_nbJets_CRZ_lep,"h_nbJets_CRZ_lep","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_dphi_CRZ_lep   [maxSysts] ;systZero.initHistogramsSysts(h_dphi_CRZ_lep,"h_dphi_CRZ_lep","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  //TH1F *h_mt_CRZ_lep     [maxSysts] ;systZero.initHistogramsSysts(h_mt_CRZ_lep,"h_mt_CRZ_lep","M_{T}",100,0,500);
  TH1F *h_mt2w_CRZ_lep   [maxSysts] ;systZero.initHistogramsSysts(h_mt2w_CRZ_lep,"h_mt2w_CRZ_lep","MT2W",100,50,500);
  TH1F *h_lepPt_CRZ_lep   [maxSysts] ;systZero.initHistogramsSysts(h_lepPt_CRZ_lep,"h_lepPt_CRZ_lep","Lepton pt distribution",200,0,100);
  //  TH1F *h_topSLMass_CRZ_lep   [maxSysts] ;systZero.initHistogramsSysts(h_topSLMass_CRZ_lep,"h_topSLMass_CRZ_lep","Top Mass distribution",100,0,500);
  
  
  //CR W+jets 1lep (CR2)
  TH1F *h_metFinal_met_0btag   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_met_0btag,"metFinal_met_0btag","MET",200,0,2000);

  TH1F *h_metFinal_met_0btag_4j   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_met_0btag_4j,"metFinal_met_0btag_4j","MET",200,0,2000);
  TH1F *h_metFinal_met_0btag_6j   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_met_0btag_6j,"metFinal_met_0btag_6j","MET",200,0,2000);

  TH1F *h_muonmetFinal_met_0btag   [maxSysts] ;systZero.initHistogramsSysts(h_muonmetFinal_met_0btag,"muonmetFinal_met_0btag","MET",200,0,2000);
  TH1F *h_electronmetFinal_met_0btag   [maxSysts] ;systZero.initHistogramsSysts(h_electronmetFinal_met_0btag,"electronmetFinal_met_0btag","MET",200,0,2000);
  
  TH1F *h_mu1Pt_CR2 [maxSysts] ;systZero.initHistogramsSysts(h_mu1Pt_CR2,"h_mu1Pt_CR2","Lepton pt distribution",400, 0.,200);
  TH1F *h_muondphi_CR2[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_CR2,"muondphi_CR2","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muondphi_6j_CR2[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_6j_CR2,"muondphi_6j_CR2","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muonnJets_CR2 [maxSysts] ;systZero.initHistogramsSysts(h_muonnJets_CR2,"h_muonnJets_CR2","Number of tight jets",12,-0.5,11.5);
  TH1F *h_muonnbJets_CR2[maxSysts] ;systZero.initHistogramsSysts(h_muonnbJets_CR2,"h_muonnbJets_CR2","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_muonjet1Pt_CR2[maxSysts] ;systZero.initHistogramsSysts(h_muonjet1Pt_CR2,"h_muonjet1Pt_CR2","Leading jet Pt distribution",500,0,500);
  TH1F *h_muonjet2Pt_CR2[maxSysts] ;systZero.initHistogramsSysts(h_muonjet2Pt_CR2,"h_muonjet2Pt_CR2","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_muonjet3Pt_CR2[maxSysts] ;systZero.initHistogramsSysts(h_muonjet3Pt_CR2,"h_muonjet3Pt_CR2","Third Jet Pt distribution",500,0,500);
  
  TH1F *h_el1Pt_CR2 [maxSysts] ;systZero.initHistogramsSysts(h_el1Pt_CR2,"h_el1Pt_CR2","Lepton pt distribution",400, 0.,200);
  TH1F *h_electrondphi_CR2[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_CR2,"electrondphi_CR2","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electrondphi_6j_CR2[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_6j_CR2,"electrondphi_6j_CR2","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electronnJets_CR2 [maxSysts] ;systZero.initHistogramsSysts(h_electronnJets_CR2,"h_electronnJets_CR2","Number of tight jets",12,-0.5,11.5);
  TH1F *h_electronnbJets_CR2[maxSysts] ;systZero.initHistogramsSysts(h_electronnbJets_CR2,"h_electronnbJets_CR2","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_electronjet1Pt_CR2[maxSysts] ;systZero.initHistogramsSysts(h_electronjet1Pt_CR2,"h_electronjet1Pt_CR2","Leading jet Pt distribution",500,0,500);
  TH1F *h_electronjet2Pt_CR2[maxSysts] ;systZero.initHistogramsSysts(h_electronjet2Pt_CR2,"h_electronjet2Pt_CR2","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_electronjet3Pt_CR2[maxSysts] ;systZero.initHistogramsSysts(h_electronjet3Pt_CR2,"h_electronjet3Pt_CR2","Third Jet Pt distribution",500,0,500);
  
  TH1F *h_dphi_CRwj_lep   [maxSysts] ;systZero.initHistogramsSysts(h_dphi_CRwj_lep,"h_dphi_CRwj_lep","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_dphi_6j_CRwj_lep[maxSysts] ;systZero.initHistogramsSysts(h_dphi_6j_CRwj_lep,"dphi_6j_CRwj_lep","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_jet1Pt_CRwj_lep   [maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_CRwj_lep,"h_jet1Pt_CRwj_lep","Leading jet Pt distribution",100,0,500);
  TH1F *h_jet2Pt_CRwj_lep   [maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_CRwj_lep,"h_jet2Pt_CRwj_lep","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_jet3Pt_CRwj_lep   [maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_CRwj_lep,"h_jet3Pt_CRwj_lep","Third Jet Pt distribution",100,0,500);
  TH1F *h_nJets_CRwj_lep   [maxSysts] ;systZero.initHistogramsSysts(h_nJets_CRwj_lep,"h_nJets_CRwj_lep","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets_CRwj_lep   [maxSysts] ;systZero.initHistogramsSysts(h_nbJets_CRwj_lep,"h_nbJets_CRwj_lep","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_mt_CRwj_lep     [maxSysts] ;systZero.initHistogramsSysts(h_mt_CRwj_lep,"h_mt_CRwj_lep","M_{T}",500,0,500);
  TH1F *h_mt2w_CRwj_lep   [maxSysts] ;systZero.initHistogramsSysts(h_mt2w_CRwj_lep,"h_mt2w_CRwj_lep","MT2W",100,50,500);
  TH1F *h_lepPt_CRwj_lep   [maxSysts] ;systZero.initHistogramsSysts(h_lepPt_CRwj_lep,"h_lepPt_CRwj_lep","MT2W",200,0,100);
  
  
  //CR tt(1lep) (CR3)
  TH1F *h_metFinal_SR_1lep   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_SR_1lep,"metFinal_SR_1lep","MET",200,0,2000);
  TH1F *h_metFinal_SR_1lep_4j   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_SR_1lep_4j,"metFinal_SR_1lep_4j","MET",200,0,2000);
  TH1F *h_metFinal_SR_1lep_6j   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_SR_1lep_6j,"metFinal_SR_1lep_6j","MET",200,0,2000);

  TH1F *h_muonmetFinal_SR_1lep   [maxSysts] ;systZero.initHistogramsSysts(h_muonmetFinal_SR_1lep,"muonmetFinal_SR_1lep","MET",200,0,2000);
  TH1F *h_electronmetFinal_SR_1lep   [maxSysts] ;systZero.initHistogramsSysts(h_electronmetFinal_SR_1lep,"electronmetFinal_SR_1lep","MET",200,0,2000);
  
  TH1F *h_lep1Pt_CR3 [maxSysts] ;systZero.initHistogramsSysts(h_lep1Pt_CR3,"h_lep1Pt_CR3","Lepton pt distribution",400, 0.,200);
  TH1F *h_lep1Eta_CR3 [maxSysts] ;systZero.initHistogramsSysts(h_lep1Eta_CR3,"h_lep1Eta_CR3","Lepton pt distribution",48, -2.4,2.4);
  TH1F *h_lep1Phi_CR3 [maxSysts] ;systZero.initHistogramsSysts(h_lep1Phi_CR3,"h_lep1Phi_CR3","Lepton pt distribution",64, -3.2,3.2);
  TH1F *h_lep1E_CR3 [maxSysts] ;systZero.initHistogramsSysts(h_lep1E_CR3,"h_lep1E_CR3","Lepton pt distribution",800, 0.,400);
  
  TH1F *h_mu1Pt_CR3 [maxSysts] ;systZero.initHistogramsSysts(h_mu1Pt_CR3,"h_mu1Pt_CR3","Lepton pt distribution",400, 0.,200);
  TH1F *h_mu1Eta_CR3 [maxSysts] ;systZero.initHistogramsSysts(h_mu1Eta_CR3,"h_mu1Eta_CR3","Lepton pt distribution",48, -2.4,2.4);
  TH1F *h_mu1Phi_CR3 [maxSysts] ;systZero.initHistogramsSysts(h_mu1Phi_CR3,"h_mu1Phi_CR3","Lepton pt distribution",64, -3.2,3.2);
  TH1F *h_mu1E_CR3 [maxSysts] ;systZero.initHistogramsSysts(h_mu1E_CR3,"h_mu1E_CR3","Lepton pt distribution",800, 0.,400);
  
  TH1F *h_el1Pt_CR3 [maxSysts] ;systZero.initHistogramsSysts(h_el1Pt_CR3,"h_el1Pt_CR3","Lepton pt distribution",400, 0.,200);
  TH1F *h_el1Eta_CR3 [maxSysts] ;systZero.initHistogramsSysts(h_el1Eta_CR3,"h_el1Eta_CR3","Lepton pt distribution",48, -2.4,2.4);
  TH1F *h_el1Phi_CR3 [maxSysts] ;systZero.initHistogramsSysts(h_el1Phi_CR3,"h_el1Phi_CR3","Lepton pt distribution",64, -3.2,3.2);
  TH1F *h_el1E_CR3 [maxSysts] ;systZero.initHistogramsSysts(h_el1E_CR3,"h_el1E_CR3","Lepton pt distribution",800, 0.,400);

  TH1F *h_jet1Pt_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_CRtt_had,"h_jet1Pt_CRtt_had","Leading jet Pt distribution",100,0,500);
  TH1F *h_jet2Pt_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_CRtt_had,"h_jet2Pt_CRtt_had","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_jet3Pt_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_CRtt_had,"h_jet3Pt_CRtt_had","Third Jet Pt distribution",100,0,500);
  TH1F *h_nJets_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_nJets_CRtt_had,"h_nJets_CRtt_had","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_nbJets_CRtt_had,"h_nbJets_CRtt_had","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_lepPt_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_lepPt_CRtt_had,"h_lepPt_CRtt_had","Lepton pt distribution",200, 0., 100);
  
  TH1F *h_dphi_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_dphi_CRtt_had,"dphi_CRtt_had","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_dphi_6j_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_dphi_6j_CRtt_had,"dphi_6j_CRtt_had","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_mt_CRtt_had  [maxSysts] ;systZero.initHistogramsSysts(h_mt_CRtt_had,"mt_CRtt_had","M_{T}",100,0,500);
  TH1F *h_mt2w_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_mt2w_CRtt_had,"mt2w_CRtt_had","MT2W",100,50,500);
  
  TH1F *h_bjet1Pt_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_bjet1Pt_CRtt_had,"h_bjet1Pt_CRtt_had","Leading bjet Pt distribution",500,0,500);
  TH1F *h_bjet2Pt_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_bjet2Pt_CRtt_had,"h_bjet2Pt_CRtt_had","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_bjet3Pt_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_bjet3Pt_CRtt_had,"h_bjet3Pt_CRtt_had","Third Bjet Pt distribution",500,0,500);
  
  TH1F *h_muonjet1Pt_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_muonjet1Pt_CRtt_had,"h_muonjet1Pt_CRtt_had","Leading jet Pt distribution",100,0,500);
  TH1F *h_muonjet2Pt_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_muonjet2Pt_CRtt_had,"h_muonjet2Pt_CRtt_had","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_muonjet3Pt_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_muonjet3Pt_CRtt_had,"h_muonjet3Pt_CRtt_had","Third Jet Pt distribution",100,0,500);
  TH1F *h_muonnJets_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_muonnJets_CRtt_had,"h_muonnJets_CRtt_had","Number of tight jets",13,-0.5,12.5);
  TH1F *h_muonnbJets_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_muonnbJets_CRtt_had,"h_muonnbJets_CRtt_had","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_muonlepPt_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_muonlepPt_CRtt_had,"h_muonlepPt_CRtt_had","Lepton pt distribution",200, 0., 100);
  TH1F *h_muondphi_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_CRtt_had,"muondphi_CRtt_had","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muondphi_6j_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_6j_CRtt_had,"muondphi_6j_CRtt_had","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muonmt_CRtt_had  [maxSysts] ;systZero.initHistogramsSysts(h_muonmt_CRtt_had,"muonmt_CRtt_had","M_{T}",100,0,500);
  TH1F *h_muonmt2w_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_muonmt2w_CRtt_had,"muonmt2w_CRtt_had","MT2W",100,50,500);
  TH1F *h_muonbjet1Pt_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet1Pt_CRtt_had,"h_muonbjet1Pt_CRtt_had","Leading bjet Pt distribution",500,0,500);
  TH1F *h_muonbjet2Pt_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet2Pt_CRtt_had,"h_muonbjet2Pt_CRtt_had","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_muonbjet3Pt_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet3Pt_CRtt_had,"h_muonbjet3Pt_CRtt_had","Third Bjet Pt distribution",500,0,500);
  
  TH1F *h_electronjet1Pt_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_electronjet1Pt_CRtt_had,"h_electronjet1Pt_CRtt_had","Leading jet Pt distribution",100,0,500);
  TH1F *h_electronjet2Pt_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_electronjet2Pt_CRtt_had,"h_electronjet2Pt_CRtt_had","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_electronjet3Pt_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_electronjet3Pt_CRtt_had,"h_electronjet3Pt_CRtt_had","Third Jet Pt distribution",100,0,500);
  TH1F *h_electronnJets_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_electronnJets_CRtt_had,"h_electronnJets_CRtt_had","Number of tight jets",13,-0.5,12.5);
  TH1F *h_electronnbJets_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_electronnbJets_CRtt_had,"h_electronnbJets_CRtt_had","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_electronlepPt_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_electronlepPt_CRtt_had,"h_electronlepPt_CRtt_had","Lepton pt distribution",200, 0., 100);
  TH1F *h_electrondphi_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_CRtt_had,"electrondphi_CRtt_had","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electrondphi_6j_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_6j_CRtt_had,"electrondphi_6j_CRtt_had","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electronmt_CRtt_had  [maxSysts] ;systZero.initHistogramsSysts(h_electronmt_CRtt_had,"electronmt_CRtt_had","M_{T}",100,0,500);
  TH1F *h_electronmt2w_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_electronmt2w_CRtt_had,"electronmt2w_CRtt_had","MT2W",100,50,500);
  TH1F *h_electronbjet1Pt_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet1Pt_CRtt_had,"h_electronbjet1Pt_CRtt_had","Leading bjet Pt distribution",500,0,500);
  TH1F *h_electronbjet2Pt_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet2Pt_CRtt_had,"h_electronbjet2Pt_CRtt_had","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_electronbjet3Pt_CRtt_had[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet3Pt_CRtt_had,"h_electronbjet3Pt_CRtt_had","Third Bjet Pt distribution",500,0,500);
  
  TH1F *h_lepMetCosDPhi_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_lepMetCosDPhi_CRtt_had,"h_lepMetCosDPhi_CRtt_had","Cosinus of #Delta #phi between lepton and MET",200, -1., 1.);
  TH1F *h_topSLMass_CRtt_had   [maxSysts] ;systZero.initHistogramsSysts(h_topSLMass_CRtt_had,"h_topSLMass_CRtt_had","Top Mass distribution",100,0,500);


  //CR tt(1lep) (CR3_nw)
  TH1F *h_metFinal_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_CR3nw,"metFinal_CR3nw","MET",200,0,2000);
  TH1F *h_muonmetFinal_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_muonmetFinal_CR3nw,"muonmetFinal_CR3nw","MET",200,0,2000);
  TH1F *h_electronmetFinal_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_electronmetFinal_CR3nw,"electronmetFinal_CR3nw","MET",200,0,2000);
  
  TH1F *h_lep1Pt_CR3nw [maxSysts] ;systZero.initHistogramsSysts(h_lep1Pt_CR3nw,"h_lep1Pt_CR3nw","Lepton pt distribution",400, 0.,200);
  TH1F *h_lep1Eta_CR3nw [maxSysts] ;systZero.initHistogramsSysts(h_lep1Eta_CR3nw,"h_lep1Eta_CR3nw","Lepton pt distribution",48, -2.4,2.4);
  TH1F *h_lep1Phi_CR3nw [maxSysts] ;systZero.initHistogramsSysts(h_lep1Phi_CR3nw,"h_lep1Phi_CR3nw","Lepton pt distribution",64, -3.2,3.2);
  TH1F *h_lep1E_CR3nw [maxSysts] ;systZero.initHistogramsSysts(h_lep1E_CR3nw,"h_lep1E_CR3nw","Lepton pt distribution",800, 0.,400);
  
  TH1F *h_mu1Pt_CR3nw [maxSysts] ;systZero.initHistogramsSysts(h_mu1Pt_CR3nw,"h_mu1Pt_CR3nw","Lepton pt distribution",400, 0.,200);
  TH1F *h_mu1Eta_CR3nw [maxSysts] ;systZero.initHistogramsSysts(h_mu1Eta_CR3nw,"h_mu1Eta_CR3nw","Lepton pt distribution",48, -2.4,2.4);
  TH1F *h_mu1Phi_CR3nw [maxSysts] ;systZero.initHistogramsSysts(h_mu1Phi_CR3nw,"h_mu1Phi_CR3nw","Lepton pt distribution",64, -3.2,3.2);
  TH1F *h_mu1E_CR3nw [maxSysts] ;systZero.initHistogramsSysts(h_mu1E_CR3nw,"h_mu1E_CR3nw","Lepton pt distribution",800, 0.,400);
  
  TH1F *h_el1Pt_CR3nw [maxSysts] ;systZero.initHistogramsSysts(h_el1Pt_CR3nw,"h_el1Pt_CR3nw","Lepton pt distribution",400, 0.,200);
  TH1F *h_el1Eta_CR3nw [maxSysts] ;systZero.initHistogramsSysts(h_el1Eta_CR3nw,"h_el1Eta_CR3nw","Lepton pt distribution",48, -2.4,2.4);
  TH1F *h_el1Phi_CR3nw [maxSysts] ;systZero.initHistogramsSysts(h_el1Phi_CR3nw,"h_el1Phi_CR3nw","Lepton pt distribution",64, -3.2,3.2);
  TH1F *h_el1E_CR3nw [maxSysts] ;systZero.initHistogramsSysts(h_el1E_CR3nw,"h_el1E_CR3nw","Lepton pt distribution",800, 0.,400);

  TH1F *h_jet1Pt_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_CR3nw,"h_jet1Pt_CR3nw","Leading jet Pt distribution",100,0,500);
  TH1F *h_jet2Pt_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_CR3nw,"h_jet2Pt_CR3nw","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_jet3Pt_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_CR3nw,"h_jet3Pt_CR3nw","Third Jet Pt distribution",100,0,500);
  TH1F *h_nJets_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_nJets_CR3nw,"h_nJets_CR3nw","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_nbJets_CR3nw,"h_nbJets_CR3nw","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_lepPt_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_lepPt_CR3nw,"h_lepPt_CR3nw","Lepton pt distribution",200, 0., 100);
  
  TH1F *h_dphi_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_dphi_CR3nw,"dphi_CR3nw","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_dphi_6j_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_dphi_6j_CR3nw,"dphi_6j_CR3nw","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_mt_CR3nw  [maxSysts] ;systZero.initHistogramsSysts(h_mt_CR3nw,"mt_CR3nw","M_{T}",100,0,500);
  TH1F *h_mt2w_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_mt2w_CR3nw,"mt2w_CR3nw","MT2W",100,50,500);
  
  TH1F *h_bjet1Pt_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_bjet1Pt_CR3nw,"h_bjet1Pt_CR3nw","Leading bjet Pt distribution",500,0,500);
  TH1F *h_bjet2Pt_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_bjet2Pt_CR3nw,"h_bjet2Pt_CR3nw","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_bjet3Pt_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_bjet3Pt_CR3nw,"h_bjet3Pt_CR3nw","Third Bjet Pt distribution",500,0,500);
  
  TH1F *h_muonjet1Pt_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_muonjet1Pt_CR3nw,"h_muonjet1Pt_CR3nw","Leading jet Pt distribution",100,0,500);
  TH1F *h_muonjet2Pt_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_muonjet2Pt_CR3nw,"h_muonjet2Pt_CR3nw","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_muonjet3Pt_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_muonjet3Pt_CR3nw,"h_muonjet3Pt_CR3nw","Third Jet Pt distribution",100,0,500);
  TH1F *h_muonnJets_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_muonnJets_CR3nw,"h_muonnJets_CR3nw","Number of tight jets",13,-0.5,12.5);
  TH1F *h_muonnbJets_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_muonnbJets_CR3nw,"h_muonnbJets_CR3nw","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_muonlepPt_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_muonlepPt_CR3nw,"h_muonlepPt_CR3nw","Lepton pt distribution",200, 0., 100);
  TH1F *h_muondphi_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_CR3nw,"muondphi_CR3nw","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muondphi_6j_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_6j_CR3nw,"muondphi_6j_CR3nw","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muonmt_CR3nw  [maxSysts] ;systZero.initHistogramsSysts(h_muonmt_CR3nw,"muonmt_CR3nw","M_{T}",100,0,500);
  TH1F *h_muonmt2w_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_muonmt2w_CR3nw,"muonmt2w_CR3nw","MT2W",100,50,500);
  TH1F *h_muonbjet1Pt_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet1Pt_CR3nw,"h_muonbjet1Pt_CR3nw","Leading bjet Pt distribution",500,0,500);
  TH1F *h_muonbjet2Pt_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet2Pt_CR3nw,"h_muonbjet2Pt_CR3nw","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_muonbjet3Pt_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_muonbjet3Pt_CR3nw,"h_muonbjet3Pt_CR3nw","Third Bjet Pt distribution",500,0,500);
  
  TH1F *h_electronjet1Pt_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_electronjet1Pt_CR3nw,"h_electronjet1Pt_CR3nw","Leading jet Pt distribution",100,0,500);
  TH1F *h_electronjet2Pt_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_electronjet2Pt_CR3nw,"h_electronjet2Pt_CR3nw","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_electronjet3Pt_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_electronjet3Pt_CR3nw,"h_electronjet3Pt_CR3nw","Third Jet Pt distribution",100,0,500);
  TH1F *h_electronnJets_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_electronnJets_CR3nw,"h_electronnJets_CR3nw","Number of tight jets",13,-0.5,12.5);
  TH1F *h_electronnbJets_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_electronnbJets_CR3nw,"h_electronnbJets_CR3nw","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_electronlepPt_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_electronlepPt_CR3nw,"h_electronlepPt_CR3nw","Lepton pt distribution",200, 0., 100);
  TH1F *h_electrondphi_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_CR3nw,"electrondphi_CR3nw","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electrondphi_6j_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_6j_CR3nw,"electrondphi_6j_CR3nw","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electronmt_CR3nw  [maxSysts] ;systZero.initHistogramsSysts(h_electronmt_CR3nw,"electronmt_CR3nw","M_{T}",100,0,500);
  TH1F *h_electronmt2w_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_electronmt2w_CR3nw,"electronmt2w_CR3nw","MT2W",100,50,500);
  TH1F *h_electronbjet1Pt_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet1Pt_CR3nw,"h_electronbjet1Pt_CR3nw","Leading bjet Pt distribution",500,0,500);
  TH1F *h_electronbjet2Pt_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet2Pt_CR3nw,"h_electronbjet2Pt_CR3nw","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_electronbjet3Pt_CR3nw[maxSysts] ;systZero.initHistogramsSysts(h_electronbjet3Pt_CR3nw,"h_electronbjet3Pt_CR3nw","Third Bjet Pt distribution",500,0,500);
  
  TH1F *h_lepMetCosDPhi_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_lepMetCosDPhi_CR3nw,"h_lepMetCosDPhi_CR3nw","Cosinus of #Delta #phi between lepton and MET",200, -1., 1.);
  TH1F *h_topSLMass_CR3nw   [maxSysts] ;systZero.initHistogramsSysts(h_topSLMass_CR3nw,"h_topSLMass_CR3nw","Top Mass distribution",100,0,500);
  
  
  //CR QCD (CR4)
  TH1F *h_metFinal_outtop   [maxSysts] ;systZero.initHistogramsSysts(h_metFinal_outtop,"metFinal_outtop","MET",200,0,2000);
  
  TH1F *h_nJets_CRqcd_had   [maxSysts] ;systZero.initHistogramsSysts(h_nJets_CRqcd_had,"h_nJets_CRqcd_had","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets_CRqcd_had   [maxSysts] ;systZero.initHistogramsSysts(h_nbJets_CRqcd_had,"h_nbJets_CRqcd_had","Number of tight b-jets",11,-0.5,10.5);
  
  TH1F *h_jet1Pt_CRqcd_had   [maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_CRqcd_had,"h_jet1Pt_CRqcd_had","Leading jet Pt distribution",100,0,500);
  TH1F *h_jet2Pt_CRqcd_had   [maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_CRqcd_had,"h_jet2Pt_CRqcd_had","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_jet3Pt_CRqcd_had   [maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_CRqcd_had,"h_jet3Pt_CRqcd_had","Third Jet Pt distribution",100,0,500);
  
  TH1F *h_bjet1Pt_CRqcd_had[maxSysts] ;systZero.initHistogramsSysts(h_bjet1Pt_CRqcd_had,"h_bjet1Pt_CRqcd_had","Leading bjet Pt distribution",500,0,500);
  TH1F *h_bjet2Pt_CRqcd_had[maxSysts] ;systZero.initHistogramsSysts(h_bjet2Pt_CRqcd_had,"h_bjet2Pt_CRqcd_had","Trailing Bjet Pt distribution",500,0,500);
  TH1F *h_bjet3Pt_CRqcd_had[maxSysts] ;systZero.initHistogramsSysts(h_bjet3Pt_CRqcd_had,"h_bjet3Pt_CRqcd_had","Third Bjet Pt distribution",500,0,500);
  
  TH1F *h_dphi_CRqcd_had[maxSysts] ;systZero.initHistogramsSysts(h_dphi_CRqcd_had,"dphi_CRqcd_had","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_dphi_6j_CRqcd_had[maxSysts] ;systZero.initHistogramsSysts(h_dphi_6j_CRqcd_had,"dphi_6j_CRqcd_had","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  
  
  //CR W/Z+jets 0 lep (CR5)
  TH1F *h_metFinal_CR5[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_CR5,"metFinal_CR5","MET",200,0,2000);
  TH1F *h_metFinal_CR5_4j[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_CR5_4j,"metFinal_CR5_4j","MET",200,0,2000);
  TH1F *h_metFinal_CR5_6j[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_CR5_6j,"metFinal_CR5_6j","MET",200,0,2000);

  TH1F *h_nJets_CR5   [maxSysts] ;systZero.initHistogramsSysts(h_nJets_CR5,"h_nJets_CR5","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets_CR5   [maxSysts] ;systZero.initHistogramsSysts(h_nbJets_CR5,"h_nbJets_CR5","Number of tight b-jets",11,-0.5,10.5);
  
  TH1F *h_jet1Pt_CR5   [maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_CR5,"h_jet1Pt_CR5","Leading jet Pt distribution",100,0,500);
  TH1F *h_jet2Pt_CR5   [maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_CR5,"h_jet2Pt_CR5","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_jet3Pt_CR5   [maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_CR5,"h_jet3Pt_CR5","Third Jet Pt distribution",100,0,500);
  
  TH1F *h_dphi_CR5[maxSysts] ;systZero.initHistogramsSysts(h_dphi_CR5,"dphi_CR5","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_dphi_6j_CR5[maxSysts] ;systZero.initHistogramsSysts(h_dphi_6j_CR5,"dphi_6j_CR5","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  
  
  //CR W+jets 1 lep (CR6)
  TH1F *h_metFinal_CR6[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_CR6,"metFinal_CR6","MET",200,0,2000);
  TH1F *h_muonmetFinal_CR6[maxSysts] ;systZero.initHistogramsSysts(h_muonmetFinal_CR6,"muonmetFinal_CR6","MET",200,0,2000);
  TH1F *h_electronmetFinal_CR6[maxSysts] ;systZero.initHistogramsSysts(h_electronmetFinal_CR6,"electronmetFinal_CR6","MET",200,0,2000);
  
  TH1F *h_nJets_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_nJets_CR6,"h_nJets_CR6","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_nbJets_CR6,"h_nbJets_CR6","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_muonnJets_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_muonnJets_CR6,"h_muonnJets_CR6","Number of tight jets",13,-0.5,12.5);
  TH1F *h_muonnbJets_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_muonnbJets_CR6,"h_muonnbJets_CR6","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_electronnJets_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_electronnJets_CR6,"h_electronnJets_CR6","Number of tight jets",13,-0.5,12.5);
  TH1F *h_electronnbJets_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_electronnbJets_CR6,"h_electronnbJets_CR6","Number of tight b-jets",11,-0.5,10.5);
  
  TH1F *h_jet1Pt_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_CR6,"h_jet1Pt_CR6","Leading jet Pt distribution",100,0,500);
  TH1F *h_jet2Pt_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_CR6,"h_jet2Pt_CR6","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_jet3Pt_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_CR6,"h_jet3Pt_CR6","Third Jet Pt distribution",100,0,500);
  
  TH1F *h_dphi_CR6[maxSysts] ;systZero.initHistogramsSysts(h_dphi_CR6,"dphi_CR6","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_dphi_6j_CR6[maxSysts] ;systZero.initHistogramsSysts(h_dphi_6j_CR6,"dphi_6j_CR6","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  
  TH1F *h_mu1Pt_CR6 [maxSysts] ;systZero.initHistogramsSysts(h_mu1Pt_CR6,"h_mu1Pt_CR6","Lepton pt distribution",400, 0.,200);
  TH1F *h_muonjet1Pt_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_muonjet1Pt_CR6,"h_muonjet1Pt_CR6","Leading jet Pt distribution",100,0,500);
  TH1F *h_muonjet2Pt_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_muonjet2Pt_CR6,"h_muonjet2Pt_CR6","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_muonjet3Pt_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_muonjet3Pt_CR6,"h_muonjet3Pt_CR6","Third Jet Pt distribution",100,0,500);
  TH1F *h_muondphi_CR6[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_CR6,"muondphi_CR6","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muondphi_6j_CR6[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_6j_CR6,"muondphi_6j_CR6","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  
  TH1F *h_el1Pt_CR6 [maxSysts] ;systZero.initHistogramsSysts(h_el1Pt_CR6,"h_el1Pt_CR6","Lepton pt distribution",400, 0.,200);
  TH1F *h_electronjet1Pt_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_electronjet1Pt_CR6,"h_electronjet1Pt_CR6","Leading jet Pt distribution",100,0,500);
  TH1F *h_electronjet2Pt_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_electronjet2Pt_CR6,"h_electronjet2Pt_CR6","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_electronjet3Pt_CR6   [maxSysts] ;systZero.initHistogramsSysts(h_electronjet3Pt_CR6,"h_electronjet3Pt_CR6","Third Jet Pt distribution",100,0,500);
  TH1F *h_electrondphi_CR6[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_CR6,"lectrondphi_CR6","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electrondphi_6j_CR6[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_6j_CR6,"lectrondphi_6j_CR6","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  
  //CR W+jets 1 lep (CR6nw)
  TH1F *h_metFinal_CR6nw[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_CR6nw,"metFinal_CR6nw","MET",200,0,2000);
  TH1F *h_metFinal_CR6nw_4j[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_CR6nw_4j,"metFinal_CR6nw_4j","MET",200,0,2000);
  TH1F *h_metFinal_CR6nw_6j[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_CR6nw_6j,"metFinal_CR6nw_6j","MET",200,0,2000);

  TH1F *h_muonmetFinal_CR6nw[maxSysts] ;systZero.initHistogramsSysts(h_muonmetFinal_CR6nw,"muonmetFinal_CR6nw","MET",200,0,2000);
  TH1F *h_electronmetFinal_CR6nw[maxSysts] ;systZero.initHistogramsSysts(h_electronmetFinal_CR6nw,"electronmetFinal_CR6nw","MET",200,0,2000);
  
  TH1F *h_nJets_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_nJets_CR6nw,"h_nJets_CR6nw","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_nbJets_CR6nw,"h_nbJets_CR6nw","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_muonnJets_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_muonnJets_CR6nw,"h_muonnJets_CR6nw","Number of tight jets",13,-0.5,12.5);
  TH1F *h_muonnbJets_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_muonnbJets_CR6nw,"h_muonnbJets_CR6nw","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_electronnJets_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_electronnJets_CR6nw,"h_electronnJets_CR6nw","Number of tight jets",13,-0.5,12.5);
  TH1F *h_electronnbJets_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_electronnbJets_CR6nw,"h_electronnbJets_CR6nw","Number of tight b-jets",11,-0.5,10.5);
  
  TH1F *h_jet1Pt_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_CR6nw,"h_jet1Pt_CR6nw","Leading jet Pt distribution",100,0,500);
  TH1F *h_jet2Pt_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_CR6nw,"h_jet2Pt_CR6nw","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_jet3Pt_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_CR6nw,"h_jet3Pt_CR6nw","Third Jet Pt distribution",100,0,500);
  
  TH1F *h_dphi_CR6nw[maxSysts] ;systZero.initHistogramsSysts(h_dphi_CR6nw,"dphi_CR6nw","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_dphi_6j_CR6nw[maxSysts] ;systZero.initHistogramsSysts(h_dphi_6j_CR6nw,"dphi_6j_CR6nw","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  
  TH1F *h_mu1Pt_CR6nw [maxSysts] ;systZero.initHistogramsSysts(h_mu1Pt_CR6nw,"h_mu1Pt_CR6nw","Lepton pt distribution",400, 0.,200);
  TH1F *h_muonjet1Pt_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_muonjet1Pt_CR6nw,"h_muonjet1Pt_CR6nw","Leading jet Pt distribution",100,0,500);
  TH1F *h_muonjet2Pt_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_muonjet2Pt_CR6nw,"h_muonjet2Pt_CR6nw","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_muonjet3Pt_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_muonjet3Pt_CR6nw,"h_muonjet3Pt_CR6nw","Third Jet Pt distribution",100,0,500);
  TH1F *h_muondphi_CR6nw[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_CR6nw,"muondphi_CR6nw","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muondphi_6j_CR6nw[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_6j_CR6nw,"muondphi_6j_CR6nw","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  
  TH1F *h_el1Pt_CR6nw [maxSysts] ;systZero.initHistogramsSysts(h_el1Pt_CR6nw,"h_el1Pt_CR6nw","Lepton pt distribution",400, 0.,200);
  TH1F *h_electronjet1Pt_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_electronjet1Pt_CR6nw,"h_electronjet1Pt_CR6nw","Leading jet Pt distribution",100,0,500);
  TH1F *h_electronjet2Pt_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_electronjet2Pt_CR6nw,"h_electronjet2Pt_CR6nw","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_electronjet3Pt_CR6nw   [maxSysts] ;systZero.initHistogramsSysts(h_electronjet3Pt_CR6nw,"h_electronjet3Pt_CR6nw","Third Jet Pt distribution",100,0,500);
  TH1F *h_electrondphi_CR6nw[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_CR6nw,"lectrondphi_CR6nw","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electrondphi_6j_CR6nw[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_6j_CR6nw,"lectrondphi_6j_CR6nw","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  
  //CR Z+jets 2 lep (CR7)
  TH1F *h_metFinal_CR7[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_CR7,"metFinal_CR7","MET",200,0,2000);
  TH1F *h_muonmetFinal_CR7[maxSysts] ;systZero.initHistogramsSysts(h_muonmetFinal_CR7,"muonmetFinal_CR7","MET",200,0,2000);
  TH1F *h_electronmetFinal_CR7[maxSysts] ;systZero.initHistogramsSysts(h_electronmetFinal_CR7,"electronmetFinal_CR7","MET",200,0,2000);
  TH1F *h_mixmetFinal_CR7[maxSysts] ;systZero.initHistogramsSysts(h_mixmetFinal_CR7,"mixmetFinal_CR7","MET",200,0,2000);
  
  TH1F *h_dphi_CR7[maxSysts] ;systZero.initHistogramsSysts(h_dphi_CR7,"dphi_CR7","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_dphi_6j_CR7[maxSysts] ;systZero.initHistogramsSysts(h_dphi_6j_CR7,"dphi_6j_CR7","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  
  TH1F *h_nJets_CR7 [maxSysts] ;systZero.initHistogramsSysts(h_nJets_CR7,"h_nJets_CR7","Number of tight jets",12,-0.5,11.5);
  TH1F *h_nbJets_CR7[maxSysts] ;systZero.initHistogramsSysts(h_nbJets_CR7,"h_nbJets_CR7","Number of tight b-jets",11,-0.5,10.5);
  
  TH1F *h_jet1Pt_CR7[maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_CR7,"h_jet1Pt_CR7","Leading jet Pt distribution",500,0,500);
  TH1F *h_jet2Pt_CR7[maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_CR7,"h_jet2Pt_CR7","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_jet3Pt_CR7[maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_CR7,"h_jet3Pt_CR7","Third Jet Pt distribution",500,0,500);
  
  TH1F *h_mu1Pt_CR7 [maxSysts] ;systZero.initHistogramsSysts(h_mu1Pt_CR7,"h_mu1Pt_CR7","Lepton pt distribution",400, 0.,200);
  TH1F *h_mu2Pt_CR7 [maxSysts] ;systZero.initHistogramsSysts(h_mu2Pt_CR7,"h_mu2Pt_CR7","Lepton pt distribution",400, 0.,200);
  TH1F *h_el1Pt_CR7 [maxSysts] ;systZero.initHistogramsSysts(h_el1Pt_CR7,"h_el1Pt_CR7","Lepton pt distribution",400, 0.,200);
  TH1F *h_el2Pt_CR7 [maxSysts] ;systZero.initHistogramsSysts(h_el2Pt_CR7,"h_el2Pt_CR7","Lepton pt distribution",400, 0.,200);
  
  TH1F *h_muondphi_CR7[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_CR7,"muondphi_CR7","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muondphi_6j_CR7[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_6j_CR7,"muondphi_6j_CR7","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muonnJets_CR7 [maxSysts] ;systZero.initHistogramsSysts(h_muonnJets_CR7,"h_muonnJets_CR7","Number of tight jets",12,-0.5,11.5);
  TH1F *h_muonnbJets_CR7[maxSysts] ;systZero.initHistogramsSysts(h_muonnbJets_CR7,"h_muonnbJets_CR7","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_muonjet1Pt_CR7[maxSysts] ;systZero.initHistogramsSysts(h_muonjet1Pt_CR7,"h_muonjet1Pt_CR7","Leading jet Pt distribution",500,0,500);
  TH1F *h_muonjet2Pt_CR7[maxSysts] ;systZero.initHistogramsSysts(h_muonjet2Pt_CR7,"h_muonjet2Pt_CR7","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_muonjet3Pt_CR7[maxSysts] ;systZero.initHistogramsSysts(h_muonjet3Pt_CR7,"h_muonjet3Pt_CR7","Third Jet Pt distribution",500,0,500);
  
  TH1F *h_electrondphi_CR7[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_CR7,"electrondphi_CR7","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electrondphi_6j_CR7[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_6j_CR7,"electrondphi_6j_CR7","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electronnJets_CR7 [maxSysts] ;systZero.initHistogramsSysts(h_electronnJets_CR7,"h_electronnJets_CR7","Number of tight jets",12,-0.5,11.5);
  TH1F *h_electronnbJets_CR7[maxSysts] ;systZero.initHistogramsSysts(h_electronnbJets_CR7,"h_electronnbJets_CR7","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_electronjet1Pt_CR7[maxSysts] ;systZero.initHistogramsSysts(h_electronjet1Pt_CR7,"h_electronjet1Pt_CR7","Leading jet Pt distribution",500,0,500);
  TH1F *h_electronjet2Pt_CR7[maxSysts] ;systZero.initHistogramsSysts(h_electronjet2Pt_CR7,"h_electronjet2Pt_CR7","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_electronjet3Pt_CR7[maxSysts] ;systZero.initHistogramsSysts(h_electronjet3Pt_CR7,"h_electronjet3Pt_CR7","Third Jet Pt distribution",500,0,500);
  
  TH1F *h_mixdphi_CR7[maxSysts] ;systZero.initHistogramsSysts(h_mixdphi_CR7,"mixdphi_CR7","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_mixdphi_6j_CR7[maxSysts] ;systZero.initHistogramsSysts(h_mixdphi_6j_CR7,"mixdphi_6j_CR7","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_mixnJets_CR7 [maxSysts] ;systZero.initHistogramsSysts(h_mixnJets_CR7,"h_mixnJets_CR7","Number of tight jets",12,-0.5,11.5);
  TH1F *h_mixnbJets_CR7[maxSysts] ;systZero.initHistogramsSysts(h_mixnbJets_CR7,"h_mixnbJets_CR7","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_mixjet1Pt_CR7[maxSysts] ;systZero.initHistogramsSysts(h_mixjet1Pt_CR7,"h_mixjet1Pt_CR7","Leading jet Pt distribution",500,0,500);
  TH1F *h_mixjet2Pt_CR7[maxSysts] ;systZero.initHistogramsSysts(h_mixjet2Pt_CR7,"h_mixjet2Pt_CR7","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_mixjet3Pt_CR7[maxSysts] ;systZero.initHistogramsSysts(h_mixjet3Pt_CR7,"h_mixjet3Pt_CR7","Third Jet Pt distribution",500,0,500);



  //CR Z+jets 2 lep (CR7nw)
  TH1F *h_metFinal_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_CR7nw,"metFinal_CR7nw","MET",200,0,2000);
  TH1F *h_metFinal_CR7nw_4j[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_CR7nw_4j,"metFinal_CR7nw_4j","MET",200,0,2000);
  TH1F *h_metFinal_CR7nw_6j[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_CR7nw_6j,"metFinal_CR7nw_6j","MET",200,0,2000);

  TH1F *h_muonmetFinal_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_muonmetFinal_CR7nw,"muonmetFinal_CR7nw","MET",200,0,2000);
  TH1F *h_electronmetFinal_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_electronmetFinal_CR7nw,"electronmetFinal_CR7nw","MET",200,0,2000);
  TH1F *h_mixmetFinal_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_mixmetFinal_CR7nw,"mixmetFinal_CR7nw","MET",200,0,2000);
  
  TH1F *h_dphi_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_dphi_CR7nw,"dphi_CR7nw","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_dphi_6j_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_dphi_6j_CR7nw,"dphi_6j_CR7nw","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  
  TH1F *h_nJets_CR7nw [maxSysts] ;systZero.initHistogramsSysts(h_nJets_CR7nw,"h_nJets_CR7nw","Number of tight jets",12,-0.5,11.5);
  TH1F *h_nbJets_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_nbJets_CR7nw,"h_nbJets_CR7nw","Number of tight b-jets",11,-0.5,10.5);
  
  TH1F *h_jet1Pt_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt_CR7nw,"h_jet1Pt_CR7nw","Leading jet Pt distribution",500,0,500);
  TH1F *h_jet2Pt_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt_CR7nw,"h_jet2Pt_CR7nw","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_jet3Pt_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt_CR7nw,"h_jet3Pt_CR7nw","Third Jet Pt distribution",500,0,500);
  
  TH1F *h_mu1Pt_CR7nw [maxSysts] ;systZero.initHistogramsSysts(h_mu1Pt_CR7nw,"h_mu1Pt_CR7nw","Lepton pt distribution",400, 0.,200);
  TH1F *h_mu2Pt_CR7nw [maxSysts] ;systZero.initHistogramsSysts(h_mu2Pt_CR7nw,"h_mu2Pt_CR7nw","Lepton pt distribution",400, 0.,200);
  TH1F *h_el1Pt_CR7nw [maxSysts] ;systZero.initHistogramsSysts(h_el1Pt_CR7nw,"h_el1Pt_CR7nw","Lepton pt distribution",400, 0.,200);
  TH1F *h_el2Pt_CR7nw [maxSysts] ;systZero.initHistogramsSysts(h_el2Pt_CR7nw,"h_el2Pt_CR7nw","Lepton pt distribution",400, 0.,200);
  
  TH1F *h_muondphi_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_CR7nw,"muondphi_CR7nw","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muondphi_6j_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_muondphi_6j_CR7nw,"muondphi_6j_CR7nw","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_muonnJets_CR7nw [maxSysts] ;systZero.initHistogramsSysts(h_muonnJets_CR7nw,"h_muonnJets_CR7nw","Number of tight jets",12,-0.5,11.5);
  TH1F *h_muonnbJets_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_muonnbJets_CR7nw,"h_muonnbJets_CR7nw","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_muonjet1Pt_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_muonjet1Pt_CR7nw,"h_muonjet1Pt_CR7nw","Leading jet Pt distribution",500,0,500);
  TH1F *h_muonjet2Pt_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_muonjet2Pt_CR7nw,"h_muonjet2Pt_CR7nw","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_muonjet3Pt_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_muonjet3Pt_CR7nw,"h_muonjet3Pt_CR7nw","Third Jet Pt distribution",500,0,500);
  
  TH1F *h_electrondphi_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_CR7nw,"electrondphi_CR7nw","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electrondphi_6j_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_electrondphi_6j_CR7nw,"electrondphi_6j_CR7nw","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_electronnJets_CR7nw [maxSysts] ;systZero.initHistogramsSysts(h_electronnJets_CR7nw,"h_electronnJets_CR7nw","Number of tight jets",12,-0.5,11.5);
  TH1F *h_electronnbJets_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_electronnbJets_CR7nw,"h_electronnbJets_CR7nw","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_electronjet1Pt_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_electronjet1Pt_CR7nw,"h_electronjet1Pt_CR7nw","Leading jet Pt distribution",500,0,500);
  TH1F *h_electronjet2Pt_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_electronjet2Pt_CR7nw,"h_electronjet2Pt_CR7nw","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_electronjet3Pt_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_electronjet3Pt_CR7nw,"h_electronjet3Pt_CR7nw","Third Jet Pt distribution",500,0,500);
  
  TH1F *h_mixdphi_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_mixdphi_CR7nw,"mixdphi_CR7nw","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_mixdphi_6j_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_mixdphi_6j_CR7nw,"mixdphi_6j_CR7nw","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_mixnJets_CR7nw [maxSysts] ;systZero.initHistogramsSysts(h_mixnJets_CR7nw,"h_mixnJets_CR7nw","Number of tight jets",12,-0.5,11.5);
  TH1F *h_mixnbJets_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_mixnbJets_CR7nw,"h_mixnbJets_CR7nw","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_mixjet1Pt_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_mixjet1Pt_CR7nw,"h_mixjet1Pt_CR7nw","Leading jet Pt distribution",500,0,500);
  TH1F *h_mixjet2Pt_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_mixjet2Pt_CR7nw,"h_mixjet2Pt_CR7nw","Trailing Jet Pt distribution",500,0,500);
  TH1F *h_mixjet3Pt_CR7nw[maxSysts] ;systZero.initHistogramsSysts(h_mixjet3Pt_CR7nw,"h_mixjet3Pt_CR7nw","Third Jet Pt distribution",500,0,500);


  //**************************
 
  
  //Full selection histograms
  TH1F* h_topHadMass_allCutsMS[maxSysts] ;
  //  float systWeights[maxSysts];
  float systWeights2B[maxSysts];
  float systWeights12B[maxSysts];
  float systWeights12BL[maxSysts];
  float systWeightsZeroB[maxSysts];
  float systWeightsNoSyst[maxSysts];
  systZero.initHistogramsSysts(h_topHadMass_allCutsMS,"h_topHadMass_allCutsMS","Top Mass distribution",100,0,500);//EXAMPLE
  
  TH1F *h_topHadMass_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topHadMass_allCuts,"h_topHadMass_allCuts","Top Mass distribution",100,0,500);
  TH1F *h_topHadPt_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topHadPt_allCuts,"h_topHadPt_allCuts","Top Pt distribution",100,0,500);
  TH1F *h_topHadCosBM_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topHadCosBM_allCuts,"h_topHadCosBM_allCuts","Best top candidate phi(b-jet, met)",100,-1,1);
  TH1F *h_topHadCosWM_allCuts  [maxSysts] ;systZero.initHistogramsSysts(h_topHadCosWM_allCuts,"h_topHadCosWM_allCuts","Best top candidate cos(W, met)",100,-1,1);
  TH1F *h_topHadCosTM_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topHadCosTM_allCuts,"h_topHadCosTM_allCuts","Bbest top candidate cos(t, met)",100,-1,1);
  TH1F *h_topHadCosWB_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topHadCosWB_allCuts,"h_topHadCosWB_allCuts","Best top candidate cos(W, b-jet)",100,-1,1);
  TH1F *h_topCosBB_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topCosBB_allCuts,"h_topHadCosBB_allCuts","Best top candidate cos(b-jet, b-jet)",100,-1,1);
  TH1F *h_topEtaBB_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topEtaBB_allCuts,"h_topHadEtaBB_allCuts","Best top candidate #Delta #eta(b-jet, b-jet)",100,0.,6);
  TH1F *h_topHadCosTT_MET_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topHadCosTT_MET_allCuts,"h_topHadCosTT_MET_allCuts","Best top candidate cos(top+top, MET)",100,-1,1);
  
  TH1F *h_topCosBB_MET_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topCosBB_MET_allCuts,"h_topHadCosBB_MET_allCuts","Best top candidate cos(b-jet + b-jet, MET)",100,-1,1);
  
  TH1F *h_topSLMass_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topSLMass_allCuts,"h_topSLMass_allCuts","Top Mass distribution",100,0,500);
  TH1F *h_topSLPt_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topSLPt_allCuts,"h_topSLPt_allCuts","Top Pt distribution",100,0,500);
  TH1F *h_topSLMT_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topSLMT_allCuts,"h_topSLMT_allCuts","Top MT distribution",100,0,500);
  TH1F *h_topSLCosBM_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topSLCosBM_allCuts,"h_topSLCosBM_allCuts","Best top candidate phi(b-jet, met)",100,-1,1);
  TH1F *h_topSLCosLM_allCuts  [maxSysts] ;systZero.initHistogramsSysts(h_topSLCosLM_allCuts,"h_topSLCosLM_allCuts","Best top candidate cos(l, met)",100,-1,1);
  TH1F *h_topSLCosTM_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topSLCosTM_allCuts,"h_topSLCosTM_allCuts","Bbest top candidate cos(t, met)",100,-1,1);
  TH1F *h_topSLCosLBM_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_topSLCosLBM_allCuts,"h_topSLCosLBM_allCuts","Best top candidate cos(l+b-jet, met)",100,-1,1);
  
  TH1F *h_minDPhi[maxSysts] ;systZero.initHistogramsSysts(h_minDPhi,"h_minDPhi","#Delta #phi (j_{1,2},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_minDPhi_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_minDPhi_allCuts,"h_minDPhi_allCuts","min dphi",100,0,3.1416);
  
  TH1F *h_minDPhi_6j[maxSysts] ;systZero.initHistogramsSysts(h_minDPhi_6j,"h_minDPhi_6j","#Delta #phi (j_{1,..,6},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_minDPhi_10j[maxSysts] ;systZero.initHistogramsSysts(h_minDPhi_10j,"h_minDPhi_10j","#Delta #phi (j_{1,..,10},E^{miss}_{T})",16, 0., 3.20);
  TH1F *h_minDPhi_6j_allCuts[maxSysts] ;systZero.initHistogramsSysts(h_minDPhi_6j_allCuts,"h_minDPhi_6j_allCuts","min dphi 6j",100,0,3.1416);
  
  TH1F *h_cutFlow = new TH1F("cutFlow","",10,-0.5,9.5);
  TH1F *h_boostCat = new TH1F("boostCat","",10,-0.5,9.5);

  TH1F *h_metFinal_Angular[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_Angular,"metFinal_Angular","MET",200,0,2000);
  TH1F *h_metFinal_Angular_tag[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_Angular_tag,"metFinal_Angular_tag","MET",200,0,2000);
  TH1F *h_metFinal_Angular_untag[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_Angular_untag,"metFinal_Angular_untag","MET",200,0,2000);
  TH1F *h_metFinal_Angular_tag_CR3[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_Angular_tag_CR3,"metFinal_Angular_tag_CR3","MET",200,0,2000);
  TH1F *h_metFinal_Angular_1tag_CR3[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_Angular_1tag_CR3,"metFinal_Angular_1tag_CR3","MET",200,0,2000);
  TH1F *h_metFinal_Angular_untag_CR3[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_Angular_untag_CR3,"metFinal_Angular_untag_CR3","MET",200,0,2000);
  //  TH1F *h_metFinal_Angular_tag_lep[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_Angular_tag_lep,"metFinal_Angular_tag_lep","MET",200,0,2000);
  //  TH1F *h_metFinal_Angular_untag_lep[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_Angular_untag_lep,"metFinal_Angular_untag_lep","MET",200,0,2000);
  
  TH1F *h_metFinal_tag[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_tag,"metFinal_tag","MET",200,0,2000);
  TH1F *h_metFinal_1tag[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_1tag,"metFinal_1tag","MET",200,0,2000);
  TH1F *h_metFinal_untag[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_untag,"metFinal_untag","MET",200,0,2000);
  TH1F *h_metFinal_tag_CR1[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_tag_CR1,"metFinal_tag_CR1","MET",200,0,2000);
  TH1F *h_metFinal_untag_CR1[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_untag_CR1,"metFinal_untag_CR1","MET",200,0,2000);
  TH1F *h_metFinal_tag_CR2[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_tag_CR2,"metFinal_tag_CR2","MET",200,0,2000);
  TH1F *h_metFinal_untag_CR2[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_untag_CR2,"metFinal_untag_CR2","MET",200,0,2000);
  TH1F *h_metFinal_tag_CR3[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_tag_CR3,"metFinal_tag_CR3","MET",200,0,2000);
  TH1F *h_metFinal_1tag_CR3[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_1tag_CR3,"metFinal_1tag_CR3","MET",200,0,2000);
  TH1F *h_metFinal_untag_CR3[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_untag_CR3,"metFinal_untag_CR3","MET",200,0,2000);
  TH1F *h_metFinal_tag_CR4[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_tag_CR4,"metFinal_tag_CR4","MET",200,0,2000);
  TH1F *h_metFinal_1tag_CR4[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_1tag_CR4,"metFinal_1tag_CR4","MET",200,0,2000);
  TH1F *h_metFinal_untag_CR4[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_untag_CR4,"metFinal_untag_CR4","MET",200,0,2000);
  TH1F *h_metFinal_tag_CR5[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_tag_CR5,"metFinal_tag_CR5","MET",200,0,2000);
  TH1F *h_metFinal_1tag_CR5[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_1tag_CR5,"metFinal_1tag_CR5","MET",200,0,2000);
  TH1F *h_metFinal_untag_CR5[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_untag_CR5,"metFinal_untag_CR5","MET",200,0,2000);
  TH1F *h_metFinal_tag_CR6[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_tag_CR6,"metFinal_tag_CR6","MET",200,0,2000);
  TH1F *h_metFinal_1tag_CR6[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_1tag_CR6,"metFinal_1tag_CR6","MET",200,0,2000);
  TH1F *h_metFinal_untag_CR6[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_untag_CR6,"metFinal_untag_CR6","MET",200,0,2000);
  TH1F *h_metFinal_tag_CR7[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_tag_CR7,"metFinal_tag_CR7","MET",200,0,2000);
  TH1F *h_metFinal_1tag_CR7[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_1tag_CR7,"metFinal_1tag_CR7","MET",200,0,2000);
  TH1F *h_metFinal_untag_CR7[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_untag_CR7,"metFinal_untag_CR7","MET",200,0,2000);
  TH1F *h_metFinal_had6Jets[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_had6Jets,"metFinal_had6Jets","MET",200,0,2000);
  TH1F *h_metFinal_had45Jets[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_had45Jets,"metFinal_had45Jets","MET",200,0,2000);
  TH1F *h_metFinal_sl3Jets[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_sl3Jets,"metFinal_sl3Jets","MET",200,0,2000);
  TH1F *h_metFinal_sl4Jets[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_sl4Jets,"metFinal_sl4Jets","MET",200,0,2000);
  TH1F *h_metFinal_noBoost[maxSysts] ;systZero.initHistogramsSysts(h_metFinal_noBoost,"metFinal_noBoost","MET",200,0,2000);
  TH1F *h_metFinal11[maxSysts] ;systZero.initHistogramsSysts(h_metFinal11,"metFinal11","MET",200,0,2000);
  TH1F *h_metFinal12[maxSysts] ;systZero.initHistogramsSysts(h_metFinal12,"metFinal12","MET",200,0,2000);
  TH1F *h_metFinal22[maxSysts] ;systZero.initHistogramsSysts(h_metFinal22,"metFinal22","MET",200,0,2000);
  TH1F *h_metFinal1Res[maxSysts] ;systZero.initHistogramsSysts(h_metFinal1Res,"metFinal1Res","MET",200,0,2000);
  TH1F *h_metFinal2Res[maxSysts] ;systZero.initHistogramsSysts(h_metFinal2Res,"metFinal2res","MET",200,0,2000);  
  
  TH1F *h_topMassPreFit[maxSysts] ;systZero.initHistogramsSysts(h_topMassPreFit,"h_topMassPreFit","Top Mass distribution",100,0,500);
  TH1F *h_topMassPostFit[maxSysts] ;systZero.initHistogramsSysts(h_topMassPostFit,"h_topMassPostFit","Top Mass distribution",100,0,500);
  TH1F *h_topMassPreFit_1tag[maxSysts] ;systZero.initHistogramsSysts(h_topMassPreFit_1tag,"h_topMassPreFit_1tag","Top Mass distribution",100,0,500);
  TH1F *h_topMassPostFit_1tag[maxSysts] ;systZero.initHistogramsSysts(h_topMassPostFit_1tag,"h_topMassPostFit_1tag","Top Mass distribution",100,0,500);
  TH1F *h_topMassPostFit_preS[maxSysts] ;systZero.initHistogramsSysts(h_topMassPostFit_preS,"h_topMassPostFit_preS","Top Mass distribution",100,0,500);
  TH1F *h_topMassPreFit_preS[maxSysts] ;systZero.initHistogramsSysts(h_topMassPreFit_preS,"h_topMassPreFit_preS","Top Mass distribution",100,0,500);
  
  TH1F *h_mva[maxSysts] ;systZero.initHistogramsSysts(h_mva,"mva","MVA",50,-1,1);
  TH1F *h_mva_preS[maxSysts] ;systZero.initHistogramsSysts(h_mva_preS,"mva_preS","MVA",50,-1,1);
  TH1F *h_mva_first[maxSysts] ;systZero.initHistogramsSysts(h_mva_first,"mva_first","MVA",50,-1,1);
  TH1F *h_mva_second[maxSysts] ;systZero.initHistogramsSysts(h_mva_second,"mva_second","MVA",50,-1,1);
  
  TH1F *h_bdt_qgid1[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid1,"h_bdt_qgl1","QGL 1st jet",50,0,1);
  TH1F *h_bdt_qgid2[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid2,"h_bdt_qgl2","QGL 2nd jet",50,0,1);
  TH1F *h_bdt_dphij1b[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij1b,"h_bdt_dphij1b","#Delta #phi (j_{1},b)",16, 0., 3.20);
  TH1F *h_bdt_dphij2b[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij2b,"h_bdt_dphij2b","#Delta #phi (j_{2},b)",16, 0., 3.20);
  TH1F *h_bdt_drj1b[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj1b,"h_bdt_drj1b","Best top candidate #Delta #eta(jet_{1}, b-jet)",100,0.,6);  
  TH1F *h_bdt_drj2b[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj2b,"h_bdt_drj2b","Best top candidate #Delta #eta(jet_{2}, b-jet)",100,0.,6);
  TH1F *h_bdt_jet1csv[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet1csv,"h_bdt_jet1csv","CSV 1st jet",100,0.,1);
  TH1F *h_bdt_jet2csv[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet2csv,"h_bdt_jet2csv","CSV 2nd jet",100,0.,1);
  TH1F *h_bdt_bjcsv[maxSysts] ;systZero.initHistogramsSysts(h_bdt_bjcsv,"h_bdt_bjcsv","CSV b jet",100,0.,1);
  TH1F *h_bdt_prob[maxSysts] ;systZero.initHistogramsSysts(h_bdt_prob,"h_bdt_prob","Fit probability",100,0.,1);
  
  TH1F *h_topHadTagCosTM  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosTM,"h_topHadTagCosTM","Best top candidate cos(t, met), vtag",100,-1,1);
  TH1F *h_topHadTagCosBM  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosBM,"h_topHadTagCosBM","Best top candidate cos(b, met), vtag",100,-1,1);
  TH1F *h_topHadTagCosWM  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosWM,"h_topHadTagCosWM","Best top candidate cos(W, met), vtag",100,-1,1);
  TH1F *h_topHadTagCosWB  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosWB,"h_topHadTagCosWB","Best top candidate cos(W, b), vtag",100,-1,1); 
  TH1F *h_topHadTagCosBB  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosBB,"h_topHadTagCosBB","Top tag, cos(b, b)",100,-1,1);
  
  TH1F *h_topHadTagCosTM_CR3  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosTM_CR3,"h_topHadTagCosTM_CR3","Best top candidate cos(t, met), vtag",100,-1,1);
  TH1F *h_topHadTagCosBM_CR3  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosBM_CR3,"h_topHadTagCosBM_CR3","Best top candidate cos(b, met), vtag",100,-1,1);
  TH1F *h_topHadTagCosWM_CR3  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosWM_CR3,"h_topHadTagCosWM_CR3","Best top candidate cos(W, met), vtag",100,-1,1);
  TH1F *h_topHadTagCosWB_CR3  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosWB_CR3,"h_topHadTagCosWB_CR3","Best top candidate cos(W, b), vtag",100,-1,1); 
  TH1F *h_topHadTagCosBB_CR3  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosBB_CR3,"h_topHadTagCosBB_CR3","Top tag, cos(b, b)",100,-1,1);
  
  TH1F *h_topHadTagCosTM_lep  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosTM_lep,"h_topHadTagCosTM_lep","Best top candidate cos(t, met), vtag",100,-1,1);
  TH1F *h_topHadTagCosBM_lep  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosBM_lep,"h_topHadTagCosBM_lep","Best top candidate cos(b, met), vtag",100,-1,1);
  TH1F *h_topHadTagCosWM_lep  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosWM_lep,"h_topHadTagCosWM_lep","Best top candidate cos(W, met), vtag",100,-1,1);
  TH1F *h_topHadTagCosWB_lep  [maxSysts] ;systZero.initHistogramsSysts(h_topHadTagCosWB_lep,"h_topHadTagCosWB_lep","Best top candidate cos(W, b), vtag",100,-1,1); 
  
  // ======================= TOP TAG ===================================
  
  TH1F *h_mva_first_CR1[maxSysts] ;systZero.initHistogramsSysts(h_mva_first_CR1,"mva_first_CR1","MVA",50,-1,1);
  
  TH1F *h_bdt_qgid1_CR1[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid1_CR1,"h_bdt_qgl1_CR1","QGL 1st jet",50,0,1);
  TH1F *h_bdt_qgid2_CR1[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid2_CR1,"h_bdt_qgl2_CR1","QGL 2nd jet",50,0,1);
  TH1F *h_bdt_dphij1b_CR1[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij1b_CR1,"h_bdt_dphij1b_CR1","#Delta #phi (j_{1},b)",16, 0., 3.20);
  TH1F *h_bdt_dphij2b_CR1[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij2b_CR1,"h_bdt_dphij2b_CR1","#Delta #phi (j_{2},b)",16, 0., 3.20);
  TH1F *h_bdt_drj1b_CR1[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj1b_CR1,"h_bdt_drj1b_CR1","Best top candidate #Delta #eta(jet_{1}, b-jet)",100,0.,6);  
  TH1F *h_bdt_drj2b_CR1[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj2b_CR1,"h_bdt_drj2b_CR1","Best top candidate #Delta #eta(jet_{2}, b-jet)",100,0.,6);
  TH1F *h_bdt_jet1csv_CR1[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet1csv_CR1,"h_bdt_jet1csv_CR1","CSV 1st jet",100,0.,1);
  TH1F *h_bdt_jet2csv_CR1[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet2csv_CR1,"h_bdt_jet2csv_CR1","CSV 2nd jet",100,0.,1);
  TH1F *h_bdt_bjcsv_CR1[maxSysts] ;systZero.initHistogramsSysts(h_bdt_bjcsv_CR1,"h_bdt_bjcsv_CR1","CSV b jet",100,0.,1);
  TH1F *h_bdt_prob_CR1[maxSysts] ;systZero.initHistogramsSysts(h_bdt_prob_CR1,"h_bdt_prob_CR1","Fit probability",100,0.,1);
  
  TH1F *h_mva_first_CR2[maxSysts] ;systZero.initHistogramsSysts(h_mva_first_CR2,"mva_first_CR2","MVA",50,-1,1);
  
  TH1F *h_bdt_qgid1_CR2[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid1_CR2,"h_bdt_qgl1_CR2","QGL 1st jet",50,0,1);
  TH1F *h_bdt_qgid2_CR2[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid2_CR2,"h_bdt_qgl2_CR2","QGL 2nd jet",50,0,1);
  TH1F *h_bdt_dphij1b_CR2[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij1b_CR2,"h_bdt_dphij1b_CR2","#Delta #phi (j_{1},b)",16, 0., 3.20);
  TH1F *h_bdt_dphij2b_CR2[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij2b_CR2,"h_bdt_dphij2b_CR2","#Delta #phi (j_{2},b)",16, 0., 3.20);
  TH1F *h_bdt_drj1b_CR2[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj1b_CR2,"h_bdt_drj1b_CR2","Best top candidate #Delta #eta(jet_{1}, b-jet)",100,0.,6);  
  TH1F *h_bdt_drj2b_CR2[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj2b_CR2,"h_bdt_drj2b_CR2","Best top candidate #Delta #eta(jet_{2}, b-jet)",100,0.,6);
  TH1F *h_bdt_jet1csv_CR2[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet1csv_CR2,"h_bdt_jet1csv_CR2","CSV 1st jet",100,0.,1);
  TH1F *h_bdt_jet2csv_CR2[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet2csv_CR2,"h_bdt_jet2csv_CR2","CSV 2nd jet",100,0.,1);
  TH1F *h_bdt_bjcsv_CR2[maxSysts] ;systZero.initHistogramsSysts(h_bdt_bjcsv_CR2,"h_bdt_bjcsv_CR2","CSV b jet",100,0.,1);
  TH1F *h_bdt_prob_CR2[maxSysts] ;systZero.initHistogramsSysts(h_bdt_prob_CR2,"h_bdt_prob_CR2","Fit probability",100,0.,1);
  
  TH1F *h_mva_first_CR3[maxSysts] ;systZero.initHistogramsSysts(h_mva_first_CR3,"mva_first_CR3","MVA",50,-1,1);
  TH1F *h_mva_second_CR3[maxSysts] ;systZero.initHistogramsSysts(h_mva_second_CR3,"mva_second_CR3","MVA",50,-1,1);
  
  TH1F *h_bdt_qgid1_CR3[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid1_CR3,"h_bdt_qgl1_CR3","QGL 1st jet",50,0,1);
  TH1F *h_bdt_qgid2_CR3[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid2_CR3,"h_bdt_qgl2_CR3","QGL 2nd jet",50,0,1);
  TH1F *h_bdt_dphij1b_CR3[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij1b_CR3,"h_bdt_dphij1b_CR3","#Delta #phi (j_{1},b)",16, 0., 3.20);
  TH1F *h_bdt_dphij2b_CR3[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij2b_CR3,"h_bdt_dphij2b_CR3","#Delta #phi (j_{2},b)",16, 0., 3.20);
  TH1F *h_bdt_drj1b_CR3[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj1b_CR3,"h_bdt_drj1b_CR3","Best top candidate #Delta #eta(jet_{1}, b-jet)",100,0.,6);  
  TH1F *h_bdt_drj2b_CR3[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj2b_CR3,"h_bdt_drj2b_CR3","Best top candidate #Delta #eta(jet_{2}, b-jet)",100,0.,6);
  TH1F *h_bdt_jet1csv_CR3[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet1csv_CR3,"h_bdt_jet1csv_CR3","CSV 1st jet",100,0.,1);
  TH1F *h_bdt_jet2csv_CR3[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet2csv_CR3,"h_bdt_jet2csv_CR3","CSV 2nd jet",100,0.,1);
  TH1F *h_bdt_bjcsv_CR3[maxSysts] ;systZero.initHistogramsSysts(h_bdt_bjcsv_CR3,"h_bdt_bjcsv_CR3","CSV b jet",100,0.,1);
  TH1F *h_bdt_prob_CR3[maxSysts] ;systZero.initHistogramsSysts(h_bdt_prob_CR3,"h_bdt_prob_CR3","Fit probability",100,0.,1);
  
  TH1F *h_mva_first_CR4[maxSysts] ;systZero.initHistogramsSysts(h_mva_first_CR4,"mva_first_CR4","MVA",50,-1,1);
  TH1F *h_mva_second_CR4[maxSysts] ;systZero.initHistogramsSysts(h_mva_second_CR4,"mva_second_CR4","MVA",50,-1,1);
  
  TH1F *h_bdt_qgid1_CR4[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid1_CR4,"h_bdt_qgl1_CR4","QGL 1st jet",50,0,1);
  TH1F *h_bdt_qgid2_CR4[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid2_CR4,"h_bdt_qgl2_CR4","QGL 2nd jet",50,0,1);
  TH1F *h_bdt_dphij1b_CR4[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij1b_CR4,"h_bdt_dphij1b_CR4","#Delta #phi (j_{1},b)",16, 0., 3.20);
  TH1F *h_bdt_dphij2b_CR4[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij2b_CR4,"h_bdt_dphij2b_CR4","#Delta #phi (j_{2},b)",16, 0., 3.20);
  TH1F *h_bdt_drj1b_CR4[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj1b_CR4,"h_bdt_drj1b_CR4","Best top candidate #Delta #eta(jet_{1}, b-jet)",100,0.,6);  
  TH1F *h_bdt_drj2b_CR4[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj2b_CR4,"h_bdt_drj2b_CR4","Best top candidate #Delta #eta(jet_{2}, b-jet)",100,0.,6);
  TH1F *h_bdt_jet1csv_CR4[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet1csv_CR4,"h_bdt_jet1csv_CR4","CSV 1st jet",100,0.,1);
  TH1F *h_bdt_jet2csv_CR4[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet2csv_CR4,"h_bdt_jet2csv_CR4","CSV 2nd jet",100,0.,1);
  TH1F *h_bdt_bjcsv_CR4[maxSysts] ;systZero.initHistogramsSysts(h_bdt_bjcsv_CR4,"h_bdt_bjcsv_CR4","CSV b jet",100,0.,1);
  TH1F *h_bdt_prob_CR4[maxSysts] ;systZero.initHistogramsSysts(h_bdt_prob_CR4,"h_bdt_prob_CR4","Fit probability",100,0.,1);
  
  TH1F *h_mva_first_CR5[maxSysts] ;systZero.initHistogramsSysts(h_mva_first_CR5,"mva_first_CR5","MVA",50,-1,1);
  TH1F *h_mva_second_CR5[maxSysts] ;systZero.initHistogramsSysts(h_mva_second_CR5,"mva_second_CR5","MVA",50,-1,1);
  
  TH1F *h_bdt_qgid1_CR5[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid1_CR5,"h_bdt_qgl1_CR5","QGL 1st jet",50,0,1);
  TH1F *h_bdt_qgid2_CR5[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid2_CR5,"h_bdt_qgl2_CR5","QGL 2nd jet",50,0,1);
  TH1F *h_bdt_dphij1b_CR5[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij1b_CR5,"h_bdt_dphij1b_CR5","#Delta #phi (j_{1},b)",16, 0., 3.20);
  TH1F *h_bdt_dphij2b_CR5[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij2b_CR5,"h_bdt_dphij2b_CR5","#Delta #phi (j_{2},b)",16, 0., 3.20);
  TH1F *h_bdt_drj1b_CR5[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj1b_CR5,"h_bdt_drj1b_CR5","Best top candidate #Delta #eta(jet_{1}, b-jet)",100,0.,6);  
  TH1F *h_bdt_drj2b_CR5[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj2b_CR5,"h_bdt_drj2b_CR5","Best top candidate #Delta #eta(jet_{2}, b-jet)",100,0.,6);
  TH1F *h_bdt_jet1csv_CR5[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet1csv_CR5,"h_bdt_jet1csv_CR5","CSV 1st jet",100,0.,1);
  TH1F *h_bdt_jet2csv_CR5[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet2csv_CR5,"h_bdt_jet2csv_CR5","CSV 2nd jet",100,0.,1);
  TH1F *h_bdt_bjcsv_CR5[maxSysts] ;systZero.initHistogramsSysts(h_bdt_bjcsv_CR5,"h_bdt_bjcsv_CR5","CSV b jet",100,0.,1);
  TH1F *h_bdt_prob_CR5[maxSysts] ;systZero.initHistogramsSysts(h_bdt_prob_CR5,"h_bdt_prob_CR5","Fit probability",100,0.,1);
  
  TH1F *h_mva_first_CR6[maxSysts] ;systZero.initHistogramsSysts(h_mva_first_CR6,"mva_first_CR6","MVA",50,-1,1);
  TH1F *h_mva_second_CR6[maxSysts] ;systZero.initHistogramsSysts(h_mva_second_CR6,"mva_second_CR6","MVA",50,-1,1);
  
  TH1F *h_bdt_qgid1_CR6[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid1_CR6,"h_bdt_qgl1_CR6","QGL 1st jet",50,0,1);
  TH1F *h_bdt_qgid2_CR6[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid2_CR6,"h_bdt_qgl2_CR6","QGL 2nd jet",50,0,1);
  TH1F *h_bdt_dphij1b_CR6[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij1b_CR6,"h_bdt_dphij1b_CR6","#Delta #phi (j_{1},b)",16, 0., 3.20);
  TH1F *h_bdt_dphij2b_CR6[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij2b_CR6,"h_bdt_dphij2b_CR6","#Delta #phi (j_{2},b)",16, 0., 3.20);
  TH1F *h_bdt_drj1b_CR6[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj1b_CR6,"h_bdt_drj1b_CR6","Best top candidate #Delta #eta(jet_{1}, b-jet)",100,0.,6);  
  TH1F *h_bdt_drj2b_CR6[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj2b_CR6,"h_bdt_drj2b_CR6","Best top candidate #Delta #eta(jet_{2}, b-jet)",100,0.,6);
  TH1F *h_bdt_jet1csv_CR6[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet1csv_CR6,"h_bdt_jet1csv_CR6","CSV 1st jet",100,0.,1);
  TH1F *h_bdt_jet2csv_CR6[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet2csv_CR6,"h_bdt_jet2csv_CR6","CSV 2nd jet",100,0.,1);
  TH1F *h_bdt_bjcsv_CR6[maxSysts] ;systZero.initHistogramsSysts(h_bdt_bjcsv_CR6,"h_bdt_bjcsv_CR6","CSV b jet",100,0.,1);
  TH1F *h_bdt_prob_CR6[maxSysts] ;systZero.initHistogramsSysts(h_bdt_prob_CR6,"h_bdt_prob_CR6","Fit probability",100,0.,1);
  
  TH1F *h_mva_first_CR7[maxSysts] ;systZero.initHistogramsSysts(h_mva_first_CR7,"mva_first_CR7","MVA",50,-1,1);
  TH1F *h_mva_second_CR7[maxSysts] ;systZero.initHistogramsSysts(h_mva_second_CR7,"mva_second_CR7","MVA",50,-1,1);
  
  TH1F *h_bdt_qgid1_CR7[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid1_CR7,"h_bdt_qgl1_CR7","QGL 1st jet",50,0,1);
  TH1F *h_bdt_qgid2_CR7[maxSysts] ;systZero.initHistogramsSysts(h_bdt_qgid2_CR7,"h_bdt_qgl2_CR7","QGL 2nd jet",50,0,1);
  TH1F *h_bdt_dphij1b_CR7[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij1b_CR7,"h_bdt_dphij1b_CR7","#Delta #phi (j_{1},b)",16, 0., 3.20);
  TH1F *h_bdt_dphij2b_CR7[maxSysts] ;systZero.initHistogramsSysts(h_bdt_dphij2b_CR7,"h_bdt_dphij2b_CR7","#Delta #phi (j_{2},b)",16, 0., 3.20);
  TH1F *h_bdt_drj1b_CR7[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj1b_CR7,"h_bdt_drj1b_CR7","Best top candidate #Delta #eta(jet_{1}, b-jet)",100,0.,6);  
  TH1F *h_bdt_drj2b_CR7[maxSysts] ;systZero.initHistogramsSysts(h_bdt_drj2b_CR7,"h_bdt_drj2b_CR7","Best top candidate #Delta #eta(jet_{2}, b-jet)",100,0.,6);
  TH1F *h_bdt_jet1csv_CR7[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet1csv_CR7,"h_bdt_jet1csv_CR7","CSV 1st jet",100,0.,1);
  TH1F *h_bdt_jet2csv_CR7[maxSysts] ;systZero.initHistogramsSysts(h_bdt_jet2csv_CR7,"h_bdt_jet2csv_CR7","CSV 2nd jet",100,0.,1);
  TH1F *h_bdt_bjcsv_CR7[maxSysts] ;systZero.initHistogramsSysts(h_bdt_bjcsv_CR7,"h_bdt_bjcsv_CR7","CSV b jet",100,0.,1);
  TH1F *h_bdt_prob_CR7[maxSysts] ;systZero.initHistogramsSysts(h_bdt_prob_CR7,"h_bdt_prob_CR7","Fit probability",100,0.,1);
  
  
  //pre-sel: top variables
  TH1F *h_topCosBB_preS [maxSysts] ;systZero.initHistogramsSysts(h_topCosBB_preS,"h_topHadCosBB_preS","Best top candidate cos(b-jet, b-jet)",100,-1,1);
  TH1F *h_topEtaBB_preS [maxSysts] ;systZero.initHistogramsSysts(h_topEtaBB_preS,"h_topHadEtaBB_preS","Best top candidate #Delta #eta(b-jet, b-jet)",100,0.,6);
  TH1F *h_topCosBB_MET_preS [maxSysts] ;systZero.initHistogramsSysts(h_topCosBB_MET_preS,"h_topHadCosBB_MET_preS","Best top candidate cos(b-jet + b-jet, MET)",100,-1,1);
  
  TH1F *h_topHadCosBM_preS [maxSysts] ;systZero.initHistogramsSysts(h_topHadCosBM_preS,"h_topHadCosBM_preS","Best top candidate phi(b-jet, met)",100,-1,1);
  TH1F *h_topHadCosWM_preS [maxSysts] ;systZero.initHistogramsSysts(h_topHadCosWM_preS,"h_topHadCosWM_preS","Best top candidate cos(W, met)",100,-1,1);
  TH1F *h_topHadCosTM_preS  [maxSysts] ;systZero.initHistogramsSysts(h_topHadCosTM_preS,"h_topHadCosTM_preS","Bbest top candidate cos(t, met)",100,-1,1);
  TH1F *h_topHadCosWB_preS  [maxSysts] ;systZero.initHistogramsSysts(h_topHadCosWB_preS,"h_topHadCosWB_preS","Best top candidate cos(W, b-jet)",100,-1,1);
  TH1F *h_topHadPt_preS [maxSysts] ;systZero.initHistogramsSysts(h_topHadPt_preS,"h_topHadPt_preS","Top Pt distribution",100,0,500);
  TH1F *h_topHadMass_preMET  [maxSysts] ;systZero.initHistogramsSysts(h_topHadMass_preMET,"h_topHadMass_preMET","Top Mass distribution",100,0,500);
  TH1F *h_topHadMass_preS  [maxSysts] ;systZero.initHistogramsSysts(h_topHadMass_preS,"h_topHadMass_preS","Top Mass distribution",100,0,500);
  TH1F *h_topHadPhiBM_preS [maxSysts] ;systZero.initHistogramsSysts(h_topHadPhiBM_preS,"h_topHadPhiBM_preS","Best top candidate phi(b-jet, met)",100,-1,1);
  TH1F *h_topHadPhiWM_preS  [maxSysts] ;systZero.initHistogramsSysts(h_topHadPhiWM_preS,"h_topHadPhiWM_preS","Best top candidate phi(W, met)",100,-1,1);
  TH1F *h_topHadPhiTM_preS  [maxSysts] ;systZero.initHistogramsSysts(h_topHadPhiTM_preS,"h_topHadPhiTM_preS","Bbest top candidate phi(t, met)",100,-1,1);
  TH1F *h_topHadPhiWB_preS   [maxSysts] ;systZero.initHistogramsSysts(h_topHadPhiWB_preS,"h_topHadPhiWB_preS","Best top candidate phi(W, b-jet)",100,-1,1);
  TH1F *h_topHadCosTT_MET_preS  [maxSysts] ;systZero.initHistogramsSysts(h_topHadCosTT_MET_preS,"h_topHadCosTT_MET_preS","Best top candidate cos(top+top, MET)",100,-1,1);
  TH1F *h_topHadChi_preS    [maxSysts] ;systZero.initHistogramsSysts(h_topHadChi_preS,"h_topHadChi_preS","Chi distribution",100,0,15);
  
  TH1F *h_topSLMT_preS [maxSysts] ;systZero.initHistogramsSysts(h_topSLMT_preS,"h_topSLMT_preS","Top MT distribution",100,0,500);
  TH1F *h_topSLCosBM_preS [maxSysts] ;systZero.initHistogramsSysts(h_topSLCosBM_preS,"h_topSLCosBM_preS","Best top candidate phi(b-jet, met)",100,-1,1);
  TH1F *h_topSLCosLM_preS  [maxSysts] ;systZero.initHistogramsSysts(h_topSLCosLM_preS,"h_topSLCosLM_preS","Best top candidate cos(l, met)",100,-1,1);
  TH1F *h_topSLCosTM_preS  [maxSysts] ;systZero.initHistogramsSysts(h_topSLCosTM_preS,"h_topSLCosTM_preS","Bbest top candidate cos(t, met)",100,-1,1);
  TH1F *h_topSLCosLBM_preS [maxSysts] ;systZero.initHistogramsSysts(h_topSLCosLBM_preS,"h_topSLCosLBM_preS","Best top candidate cos(l+b-jet, met)",100,-1,1);
  
  TH1F *h_topSLPt_preS  [maxSysts] ;systZero.initHistogramsSysts(h_topSLPt_preS,"h_topSLPt_preS","Top Pt distribution",100,0,500);
  TH1F *h_topSLMass_preS   [maxSysts] ;systZero.initHistogramsSysts(h_topSLMass_preS,"h_topSLMass_preS","Top Mass distribution",100,0,500);
  TH1F *h_topSLMass_preS_NS   [maxSysts] ;systZero.initHistogramsSysts(h_topSLMass_preS_NS,"h_topSLMass_preS_NS","Top Mass distribution",100,0,500);
  TH1F *h_topSLMass_preMET   [maxSysts] ;systZero.initHistogramsSysts(h_topSLMass_preMET,"h_topSLMass_preMET","Top Mass distribution",100,0,500);
  
  TH1F *h_topSLPhiBM_preS [maxSysts] ;systZero.initHistogramsSysts(h_topSLPhiBM_preS,"h_topSLPhiBM_preS","Best top candidate phi(b-jet, met)",100,-1,1);
  TH1F *h_topSLPhiLM_preS [maxSysts] ;systZero.initHistogramsSysts(h_topSLPhiLM_preS,"h_topSLPhiLM_preS","Best top candidate phi(l, met)",100,-1,1);
  TH1F *h_topSLPhiTM_preS [maxSysts] ;systZero.initHistogramsSysts(h_topSLPhiTM_preS,"h_topSLPhiTM_preS","Bbest top candidate phi(t, met)",100,-1,1);
  TH1F *h_topSLPhiLBM_preS [maxSysts] ;systZero.initHistogramsSysts(h_topSLPhiLBM_preS,"h_topSLPhiLBM_preS","Best top candidate phi(l+b-jet, met)",100,-1,1);
  TH1F *h_topSLChi_preS    [maxSysts] ;systZero.initHistogramsSysts(h_topSLChi_preS,"h_topSLChi_preS","Chi distribution",100,0,15);
  
  TH1F *h_topHadmassDrop_preS    [maxSysts] ;systZero.initHistogramsSysts(h_topHadmassDrop_preS,"h_topHadmassDrop_preS","Best top candidate mass drop",40,0,40);
  //TH1F *h_jetQGL_preS  [maxSysts] ;systZero.initHistogramsSysts(h_jetQGL_preS,"h_jetQGL_preS","QGL jets",50,0,1);
  

 // **********************************
 // Definition of output tree
 // **********************************

 string outFileName =sample + "_" +channel+".root";;
 
 TFile *outTree = TFile::Open((outdir+"/trees/tree_"+outFileName+syststrname).c_str(), "RECREATE");
 TTree ttDMTree("ttDMTree","ttDM tree");
 
 //for after-cut variables
 // double d_topHadMass_allCuts, d_topHadPt_allCuts, d_topHadCosBM_allCuts, d_topHadCosWM_allCuts, d_topHadCosTM_allCuts, d_topHadCosWB_allCuts, d_topHadCosBB_allCuts, d_topHadEtaBB_allCuts, d_topHadSecondTopMass_allCuts, d_topHadSecondTopPt_allCuts, d_topHadCosTT_allCuts;

 // double d_topSLMass_allCuts, d_topSLPt_allCuts  , d_topSLMT_allCuts, d_topSLCosBM_allCuts, d_topSLCosLM_allCuts, d_topSLCosTM_allCuts, d_topSLCosLBM_allCuts,  d_topSLSecondTopMass_allCuts, d_topSLSecondTopPt_allCuts, d_topSLCosTT_allCuts;
 
 //for preselection-level variables
 // double d_topHadMass_presel, d_topHadPt_presel, d_topHadCosBM_presel, d_topHadCosWM_presel, d_topHadCosTM_presel, d_topHadCosWB_presel, d_topHadCosBB_presel, d_topHadEtaBB_presel, d_topHadSecondTopMass_presel, d_topHadSecondTopPt_presel, d_topHadCosTT_presel;

 // double d_topSLMass_presel, d_topSLPt_presel  , d_topSLMT_presel, d_topSLCosBM_presel, d_topSLCosLM_presel, d_topSLCosTM_presel, d_topSLCosLBM_presel,d_topSLSecondTopMass_presel, d_topSLSecondTopPt_presel, d_topSLCosTT_presel;

 
 double d_topHadMass, d_topHadWMass, d_topHadPt, d_topHadChi, d_topHadCosBM, d_topHadCosWM, d_topHadCosTM, d_topHadCosWB, d_topCosBB, d_topCosBB_MET,d_topHadCosTT_MET,  d_topEtaBB, d_topHadSecondTopMass, d_topHadSecondTopPt, d_topHadCosTT, d_topHadEtaTT, d_topHadSecondTopWMass;

 double d_topSLMass, d_topSLPt  , d_topSLMT, d_topSLChi, d_topSLCosBM, d_topSLCosLM, d_topSLCosTM, d_topSLCosLBM, d_topSLSecondTopMass, d_topSLSecondTopPt, d_topSLCosTT;//,d_topHadCosTL,d_topHadEtaTL;

 double metFinal, d_njets,d_cutLevel, met;
 double minDPhi(0.), minDPhi_6j(0.),  minDPhi_10j(0.); 
 double d_mt,d_mt2w;
   //d_dphi, d_mt, d_mt2w, d_njets, d_cutLevel;

 double d_nGoodPV,d_nPV;
 
 ttDMTree.Branch("metFinal", &metFinal, "metFinal/D");
 ttDMTree.Branch("cutLevel", &d_cutLevel, "cutLevel/D");
 
 ttDMTree.Branch("njets", &d_njets, "njets/D");
 ttDMTree.Branch("met", &met, "met/D");
 ttDMTree.Branch("mt", &d_mt, "mt/D");
 ttDMTree.Branch("mt2w", &d_mt2w, "mt2w/D");
 ttDMTree.Branch("minDPhi", &minDPhi, "minDPhi/D");
 ttDMTree.Branch("minDPhi_6j", &minDPhi_6j, "minDPhi_6j/D");


 ttDMTree.Branch("nGoodPV", &d_nGoodPV, "nGoodPV/D");
 ttDMTree.Branch("nPV", &d_nPV, "nPV/D");
 
 // Fullhad
 //FullCuts
 ttDMTree.Branch("topHadMass", &d_topHadMass, "topHadMass/D");
 ttDMTree.Branch("topHadWMass", &d_topHadWMass, "topHadWMass/D");
 ttDMTree.Branch("topHadPt", &d_topHadPt, "topHadPt/D");
 ttDMTree.Branch("topHadChi", &d_topHadChi, "topHadChi/D");

 ttDMTree.Branch("topHadCosBM", &d_topHadCosBM, "topHadCosBM/D");
 ttDMTree.Branch("topHadCosWM", &d_topHadCosWM, "topHadCosWM/D");
 ttDMTree.Branch("topHadCosTM", &d_topHadCosTM, "topHadCosTM/D");
 
 ttDMTree.Branch("topHadCosWB", &d_topHadCosWB, "topHadCosWB/D");
 ttDMTree.Branch("topCosBB", &d_topCosBB, "topCosBB/D");
 ttDMTree.Branch("topEtaBB", &d_topEtaBB, "topEtaBB/D");
 ttDMTree.Branch("topCosBB_MET", &d_topCosBB_MET, "topCosBB_MET/D");
 ttDMTree.Branch("topHadCosTT_MET", &d_topHadCosTT_MET, "topHadCosTT_MET/D");


 ttDMTree.Branch("topHadSecondTopMass", &d_topHadSecondTopMass, "topHadSecondTopMass/D");
 ttDMTree.Branch("topHadSecondTopWMass", &d_topHadSecondTopWMass, "topHadSecondTopWMass/D");
 ttDMTree.Branch("topHadSecondTopPt", &d_topHadSecondTopPt, "topHadSecondTopPt/D");
 ttDMTree.Branch("topHadCosTT", &d_topHadCosTT, "topHadCosTT/D");
 ttDMTree.Branch("topHadEtaTT", &d_topHadEtaTT, "topHadEtaTT/D");
 
 //Presel
 /* ttDMTree.Branch("topHadMass_presel", &d_topHadMass_presel, "topHadMass_presel/D");
 ttDMTree.Branch("topHadPt_presel", &d_topHadPt_presel, "topHadPt_presel/D");

 ttDMTree.Branch("topHadCosBM_presel", &d_topHadCosBM_presel, "topHadCosBM_presel/D");
 ttDMTree.Branch("topHadCosWM_presel", &d_topHadCosWM_presel, "topHadCosWM_presel/D");
 ttDMTree.Branch("topHadCosTM_presel", &d_topHadCosTM_presel, "topHadCosTM_presel/D");

 ttDMTree.Branch("topHadCosWB_presel", &d_topHadCosWB_presel, "topHadCosWB_presel/D");
 ttDMTree.Branch("topHadCosBB_presel", &d_topHadCosBB_presel, "topHadCosBB_presel/D");
 ttDMTree.Branch("topHadEtaBB_presel", &d_topHadEtaBB_presel, "topHadEtaBB_presel/D");

 ttDMTree.Branch("topHadSecondTopMass_presel", &d_topHadSecondTopMass_presel, "topHadSecondTopMass_presel/D");
 ttDMTree.Branch("topHadSecondTopPt_presel", &d_topHadSecondTopPt_presel, "topHadSecondTopPt_presel/D");
 ttDMTree.Branch("topHadCosTT_presel", &d_topHadCosTT_presel, "topHadCosTT_presel/D");
 */
 
 //SemiLep
 //allCuts
 ttDMTree.Branch( "topSLMass", &d_topSLMass ,"topSLMass/D");
 ttDMTree.Branch( "topSLPt", &d_topSLPt ,"topSLPt/D");
 ttDMTree.Branch( "topSLMT", &d_topSLMT,"topSLMT/D");
 ttDMTree.Branch( "topSLChi", &d_topSLChi ,"topSLChi/D");

 ttDMTree.Branch( "topSLCosBM", &d_topSLCosBM,"topSLCosBM/D");
 ttDMTree.Branch( "topSLCosLM", &d_topSLCosLM,"topSLCosLM/D");
 ttDMTree.Branch( "topSLCosTM", &d_topSLCosTM,"topSLCosTM/D");
 ttDMTree.Branch( "topSLCosLBM", &d_topSLCosLBM,"topSLCosLBM/D");
 ttDMTree.Branch( "topSLSecondTopMass", &d_topSLSecondTopMass,"topSLSecondTopMass/D");
 ttDMTree.Branch( "topSLSecondTopPt", &d_topSLSecondTopPt,"topSLSecondTopPt/D");
 ttDMTree.Branch( "topSLCosTT", &d_topSLCosTT,"topSLCosTT/D");
 //ttDMTree.Branch("topHadCosTL", &d_topHadCosTL, "topHadCosTL/D");
 //ttDMTree.Branch("topHadEtaTL", &d_topHadEtaTL, "topHadEtaTL/D");

 
 //Presel
 /* ttDMTree.Branch( "topSLMass_presel", d_topSLMass_presel ,"topSLMass_presel");
 ttDMTree.Branch( "topSLPt_presel", d_topSLPt_presel ,"topSLPt_presel");
 ttDMTree.Branch( "topSLMT_presel",d_topSLMT_presel,"topSLMT_presel");
 ttDMTree.Branch( "topSLCosBM_presel",d_topSLCosBM_presel,"topSLCosBM_presel");
 ttDMTree.Branch( "topSLCosLM_presel",d_topSLCosLM_presel,"topSLCosLM_presel");
 ttDMTree.Branch( "topSLCosTM_presel",d_topSLCosTM_presel,"topSLCosTM_presel");
 ttDMTree.Branch( "topSLCosLBM_presel",d_topSLCosLBM_presel,"topSLCosLBM_presel");
 ttDMTree.Branch( "topSLSecondTopMass_presel",d_topSLSecondTopMass_presel,"topSLSecondTopMass_presel");
 ttDMTree.Branch( "topSLSecondTopPt_presel",d_topSLSecondTopPt_presel,"topSLSecondTopPt_presel");
 ttDMTree.Branch( "topSLCosTT_presel",d_topSLCosTT_presel,"topSLCosTT_presel");
 */
 
 //PU
 //**** 
 edm::LumiReWeighting LumiWeights_, LumiWeightsUp_, LumiWeightsDown_;
 
 if(isData=="MC"){
   LumiWeights_ = edm::LumiReWeighting("data/puMC.root", "data/MyDataPileupHistogram.root","MC_pu","pileup");
   LumiWeightsUp_ = edm::LumiReWeighting("data/puMC.root", "data/MyDataPileupHistogramUP.root","MC_pu","pileup");
   LumiWeightsDown_ = edm::LumiReWeighting("data/puMC.root", "data/MyDataPileupHistogramDOWN.root","MC_pu","pileup");
 }


  // File with lepton eff weights
 TFile* file_eff = TFile::Open("data/Weights.root");
 Weights elEff( file_eff , "EleIDEffTight" );
 Weights elTrigEff( file_eff , "Ele27TriggerEff" );
 Weights muEff( file_eff , "MuIDEffTight" );
 Weights muTrigEff( file_eff , "IsoMu20ORIsoTkMu20Eff" );
 Weights muIsoEff( file_eff , "MuTightRelIsoTightID" );

 kFactor wQCD("WQCD");
 kFactor wEWK("WEWK");
 kFactor wQCDrenUp("WQCDrenUp");
 kFactor wQCDrenDown("WQCDrenDown");
 kFactor wQCDfacUp("WQCDfacUp");
 kFactor wQCDfacDown("WQCDfacDown");

 kFactor zQCD("ZQCD");
 kFactor zEWK("ZEWK");
 kFactor zQCDrenUp("ZQCDrenUp");
 kFactor zQCDrenDown("ZQCDrenDown");
 kFactor zQCDfacUp("ZQCDfacUp");
 kFactor zQCDfacDown("ZQCDfacDown");

 
 float kfact(1.);
 float kfact_qcd(1.),  kfact_ewk(1.);
 float kfact_qcdrenUp(1.), kfact_qcdrenDown(1.), kfact_qcdfacUp(1.), kfact_qcdfacDown(1.), kfact_ewkUp(1.) ,kfact_ewkDown(1.);


 cout << "LOADING EFFICIENCY TABLES" << endl;


 // **********************************
 // EVENT LOOP
 // **********************************

 //nEventsPrePres = 1000;
 // nEventsPrePres = 10;


 TStopwatch watch;

 for(Int_t i=0; i<nEventsPrePres; i++ )
   {

     mva.clear(); 
     topFitRes.clear();
     topFitResFull.clear();

     // if(i%100000==1 ){
     //   cout<<"Event: "<<i<<endl; 
     //   cout << "ok " << endl;
     // }


     if (i % 100000 == 0 ) {
       cout << "--- Event: " << i << endl;

       if ( i != 0 ) {
         // cout << "ok " << endl;
         watch.Stop();
         watch.Print();
         Double_t real = watch.RealTime();

         Double_t leftReal = real / ((double)i) * (nEventsPrePres - i);

         std::cout << "--- Estimated Time to Completion "
         << std::fixed
         << setprecision(0) << setw(2) << setfill('0') << TMath::Floor(leftReal / 60.) << ":"
         << setprecision(2) << setw(5) << setfill('0') <<fmod(leftReal, 60.)
         // << " (estimated)"
         << std::endl;
         watch.Continue();
       }
     }



     chain.GetEntry(i);
     w = LHEWeightSign[0];

     int maxJetLoop = min(15, jetSize);
     int maxMuLoop = min(6, muSize);
     int maxElLoop = min(6, elSize);



     if( !strncmp(sample.c_str(), Wlabel , strlen(Wlabel))) {
        kfact_qcd = wQCD.getkFact(WPt);
        kfact_ewk = wEWK.getkFact(WPt);

	kfact_qcdrenUp = wQCDrenUp.getSyst(WPt); 
	kfact_qcdrenDown = wQCDrenDown.getSyst(WPt); 
	kfact_qcdfacUp = wQCDfacUp.getSyst(WPt); 
	kfact_qcdfacDown = wQCDfacDown.getSyst(WPt); 

	kfact_ewkUp=1.;
	kfact_ewkDown=kfact_ewk;

	kfact =  kfact_qcd *  kfact_ewk;
	k_fact = kfact;
	//std::cout<< "Boson Pt: "<< WPt<<std::endl;
     }

     else if( !strncmp(sample.c_str(), Zlabel , strlen(Zlabel))  or !strncmp(sample.c_str(), Zlabel , strlen(DYlabel)) ) {
        kfact_qcd = zQCD.getkFact(ZPt);
        kfact_ewk = zEWK.getkFact(ZPt);

	kfact_qcdrenUp = zQCDrenUp.getSyst(ZPt); 
	kfact_qcdrenDown = zQCDrenDown.getSyst(ZPt); 
	kfact_qcdfacUp = zQCDfacUp.getSyst(ZPt); 
	kfact_qcdfacDown = zQCDfacDown.getSyst(ZPt); 

	kfact_ewkUp=1.;
	kfact_ewkDown=kfact_ewk;
	
	kfact =  kfact_qcd *  kfact_ewk;
	k_fact = kfact;
	//std::cout<< "Boson Pt: "<< ZPt<<std::endl;
     }
     else{

	kfact =  1.;
     }

     //     std::cout <<" Factor: "<<kfact<<std::endl;
     //     std::cout<<"k-factor "<<k_fact<<std::endl;


     bool skip = false;    
     if(isData=="DATA"){
       if(runNumber==259626 || runNumber== 259637 || runNumber== 259681 || runNumber== 259682 || runNumber== 259683 || runNumber== 259685){
     	 skip = true;
       }
     }

     if(skip==true && isData=="DATA") {
       cout << "skipping event " << runNumber << " " << lumiSec << " " << evtNumber << endl;
       continue;
     }



     //if (isData=="DATA" && runNumber>runNumbBlinding) continue;

     //if(evtNumber!=13230 and evtNumber!=17394 and evtNumber!=26088 and evtNumber!=29912  and evtNumber!=40514)continue;
     //else{std::cout<<"Before entering the selection"<<std::endl;}

     if(isData=="MC"){
       w = LHEWeightSign[0];
       //       w=1.;
       w_pu = LumiWeights_.weight(numTrueInt);
       
       w = w * w_pu;
       w*=k_fact;
       //cout << " --> bef top pt " << w << endl;
       
       if(sample=="TT"){
	 w*=w_top/topWeight;
       }

       //cout << "     after top pt " << w << endl; 
       //w*= ((METturnon13TeV(metPt[0],1,true))/(METturnon13TeV(metPt[0],1,false)));
       //       cout << "vhf "<<vhf << " eventFlavour "<<eventFlavour << endl;
       vhf=(int)eventFlavour;
       //FIXXX
       //       w_top=1.0;
       //       topWeight=1.0;
       //std::cout<<"K-factor: "<<k_fact<<std::endl;
     }
     
     //added for sync
     //if(isData=="DATA"){
     //w=1.;
     //w*=k_fact;
     //w*= ((METturnon13TeV(metPt[0],2,true))/(METturnon13TeV(metPt[0],2,false)));
     //cout << "--> met trig " << (METturnon13TeV(metPt[0],2,true)/METturnon13TeV(metPt[0],2,false)) << endl;
     //}
     

     //initialization tree variables
     d_topHadMass=-999.; d_topHadPt=-999.; d_topHadChi=-999.; d_topHadCosBM=-999.; d_topHadCosWM=-999.; d_topHadCosTM=-999.; d_topHadCosWB=-999.; d_topCosBB=-999.;  d_topCosBB_MET=-999.; d_topHadCosTT_MET=-999.;d_topEtaBB=-999.; d_topHadSecondTopMass=-999.; d_topHadSecondTopPt=-999.; d_topHadCosTT=-999., d_topHadEtaTT=-999., d_topHadSecondTopMass=-999., d_topHadSecondTopPt=-999., d_topHadSecondTopWMass=-999.;


     d_topSLMass=-999.; d_topSLPt  =-999.; d_topSLChi  =-999.; d_topSLMT=-999.; d_topSLCosBM=-999.; d_topSLCosLM=-999.; d_topSLCosTM=-999.; d_topSLCosLBM=999.; d_topSLSecondTopMass=-999.; d_topSLSecondTopPt=-999.; d_topSLCosTT=-999.;

     metFinal=-999.;  d_njets=0; d_cutLevel=-999.;
     //     std::cout<<"We are in the loop and met filters are: "<<(passHBHE > 0.0 &&  passMETFilters > 0.0)<<std::endl;
     if( isData=="DATA" && !(passHBHE > 0.0 && passHBHEIso > 0.0 &&  passMETFilters > 0.0))continue ;

     passElTrig = slTrigEle_v1>0. || slTrigEle_v2>0.;
     passMuTrig = slTrigIsoMu20_v1>0. || slTrigIsoMu20_v2>0. || slTrigIsoMu20_v3>0. || slTrigIsoMuTk20_v1>0. || slTrigIsoMuTk20_v2>0. ||  slTrigIsoMuTk20_v3>0. || slTrigIsoMuTk20_v4>0.;

     bool dataTriggerSL = ( ((elePD || isData=="MC") && passElTrig > 0.) || ((muPD || isData=="MC") && passMuTrig > 0. && !(passElTrig>0.)) );
     bool dataTriggerHad = ( (hadTrigNoiseCleaned_v1>0. or hadTrigNoiseCleaned_v2>0. or hadTrigJetCleaned_v1>0. or hadTrigJetCleaned_v2>0.) );
     //&& !(passElTrig>0. || passMuTrig>0.) );

     //if(channel == "semileptonic" && lepton == "Mu" && not(passMuTrig>0.)) continue;
     //if(channel == "semileptonic" && lepton == "El" && not(passElTrig>0.)) continue;

       
     // Removing connection of channel and trigger for fullhadronic case
     if(dataTriggerHad || dataTriggerSL){ 

       d_cutLevel =0;
              
       // Compute minDPhi     
       met = metPt[0];

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
       }
       //end loop on electrons     

       for(int m= 0; m<maxMuLoop;++m ){
	 if(muPt[m]>0){
	   mu.SetPtEtaPhiE(muPt[m], muEta[m], muPhi[m],muE[m]);
	   tightMu.push_back(mu);
	 }
       }
       //end loop on muons 
       
       nMu = tightMu.size();
       nEl = tightEl.size();
       //       nLooseMu=1.;//nLooseMuons;
       //       nVetoEl=1.;nVetoElectrons;
              
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

       /*Define lepton 1 and lepton 2*/

       float muIDweight(1.), muTrigweight(1), muIsoweight(1.) ;
       float elIDweight(1.), elTrigweight(1.);

       float muIDweightUp(1.), muIDweightDown(1.), muTrigweightUp(1), muTrigweightDown(1), muIsoweightUp(1.), muIsoweightDown(1.) ;
       float elIDweightUp(1.), elIDweightDown(1.), elTrigweightUp(1.), elTrigweightDown(1.);

       float muID2weightUp(1.), muID2weightDown(1.), muIso2weightUp(1.), muIso2weightDown(1.) ;
       float elID2weightUp(1.), elID2weightDown(1.);
       
       if( leptons.size()>0){
	 lep1 = leptons[0].p4;
	 lep1Flavour = leptons[0].flavour;
	 /* Apply trigger plus ID/Iso efficiencies */
	 if(lep1Flavour==11){
	   elIDweight = elEff.getEff(lep1.Eta(), lep1.Pt()); 
	   elIDweightUp = (elIDweight + elEff.getErr(lep1.Eta(), lep1.Pt()))/elIDweight;
	   elIDweightDown = (elIDweight - elEff.getErr(lep1.Eta(), lep1.Pt()))/elIDweight;

       	   elTrigweight = elTrigEff.getEff( fabs(lep1.Eta()), lep1.Pt());	 
       	   elTrigweightUp = (elTrigweight + (elTrigweight*0.02) )/elTrigweight;	 
       	   elTrigweightDown = (elTrigweight - (elTrigweight*0.02))/elTrigweight;	 

       	   if(isData=="MC")w = w*elTrigweight * elIDweight;
       	 }
       	 else{
       	   muTrigweight = muTrigEff.getEff( fabs(lep1.Eta()), lep1.Pt() );
	   muTrigweightUp = (muTrigweight + elEff.getErr(fabs(lep1.Eta()), lep1.Pt()))/muTrigweight;
           muTrigweightDown = (muTrigweight - elEff.getErr(fabs(lep1.Eta()), lep1.Pt()))/muTrigweight;

	   muIDweight = muEff.getEff(fabs(lep1.Eta()), lep1.Pt() );
	   muIDweightUp = (muIDweight + muEff.getErr(fabs(lep1.Eta()), lep1.Pt()))/muIDweight;
           muIDweightDown = (muIDweight - muEff.getErr(fabs(lep1.Eta()), lep1.Pt()))/muIDweight;

	   muIsoweight = muIsoEff.getEff( fabs(lep1.Eta()), lep1.Pt() );   
	   muIsoweightUp = (muIsoweight + muEff.getErr(fabs(lep1.Eta()), lep1.Pt()))/muIsoweight;
           muIsoweightDown = (muIsoweight - muEff.getErr(fabs(lep1.Eta()), lep1.Pt()))/muIsoweight;
    	 
	   if(isData=="MC")w = w* muTrigweight * muIsoweight * muIDweight ;
       	 }
	 
	 if( leptons.size()>1){
	   lep2 = leptons[1].p4;
	   lep2Flavour = leptons[1].flavour;
	   if(lep2Flavour==11){
	     elIDweight = elEff.getEff(lep2.Eta(), lep2.Pt()); 
	     elID2weightUp = (elIDweight + elEff.getErr(lep2.Eta(), lep2.Pt()))/elIDweight;
	     elID2weightDown = (elIDweight - elEff.getErr(lep2.Eta(), lep2.Pt()))/elIDweight;

	     if(isData=="MC")w = w* elIDweight;
	   }
	   else{
	     muIDweight = muEff.getEff(fabs(lep2.Eta()), lep2.Pt() );
	     muID2weightUp = (muIDweight + muEff.getErr(fabs(lep2.Eta()), lep2.Pt()))/muIDweight;
	     muID2weightDown = (muIDweight - muEff.getErr(fabs(lep2.Eta()), lep2.Pt()))/muIDweight;

	     muIsoweight = muIsoEff.getEff( fabs(lep2.Eta()), lep2.Pt() );       	 
	     muIso2weightUp = (muIsoweight + muEff.getErr(fabs(lep2.Eta()), lep2.Pt()))/muIsoweight;
	     muIso2weightDown = (muIsoweight - muEff.getErr(fabs(lep2.Eta()), lep2.Pt()))/muIsoweight;

	     if(isData=="MC")w = w* muIsoweight * muIDweight ;
	   }
	 }
       }

       if(isData=="MC" and ((nTightElectrons+nTightMuons ==0) &&  (nLooseMuons+nVetoElectrons ==0))){
         w*= ((METturnon13TeV(metPt[0],1,true))/(METturnon13TeV(metPt[0],1,false)));
       }
       
       double  wcats[10]={1.,1.,1.,1.,1.,1.,1.,1.,1.,1.};

       if(isData=="DATA"){
	 /*	 systWeights[NOSYST]=1;
		 systWeights[BTAGUP]=1;
		 systWeights[BTAGDOWN]=1;
		 systWeights[MISTAGUP]=1;
		 systWeights[MISTAGDOWN]=1;
		 systWeights[PUUP]=1;
		 systWeights[PUDOWN]=1;
	 */

       
	 systZero.setWeight(0,1.);
	 systZero.setWeight("btagUp",1.);
	 systZero.setWeight("btagDown",1.);
	 systZero.setWeight("mistagUp",1.);
	 systZero.setWeight("mistagDown",1.);
	 systZero.setWeight("puDown",1.);
	 systZero.setWeight("puUp",1.);
	 systZero.setWeight("lepDown",1.);
	 systZero.setWeight("lepUp",1.);

	 if(addTTSplit){
	   wcats[0]=1.0;
	   wcats[1]=double(NMCLeptons==0);
	   wcats[2]=double(NMCLeptons==1);
	   wcats[3]=double(NMCLeptons==2);
	   systZero.setWCats(wcats);
	 }

	 //systZero.setWeight("isoUp",1.);
	 //systZero.setWeight("trigDown",1.);
	 //systZero.setWeight("trigUp",1.);

	 syst12B.copySysts(systZero);
	 syst12B.setWeight("btagUp",1.);
	 syst12B.setWeight("btagDown",1.);

	 syst12BL.copySysts(systZero);
	 syst12BL.setWeight("btagUp",1.);
	 syst12BL.setWeight("btagDown",1.);

	 syst2B.copySysts(systZero);
	 syst2B.setWeight("btagUp",1.);
	 syst2B.setWeight("btagDown",1.);

	 syst0B.copySysts(systZero);
	 syst0B.setWeight("btagUp",1.);
	 syst0B.setWeight("btagDown",1.);

	 systWeights12B[NOSYST]=1.;
	 systWeights12B[BTAGUP]=1.;
	 systWeights12B[BTAGDOWN]=1.;
	 systWeights12B[MISTAGUP]=1.;
	 systWeights12B[MISTAGDOWN]=1.;
	 systWeights12B[PUUP]=1.;
	 systWeights12B[PUDOWN]=1.;
	 systWeights12B[LEPUP]=1.;
         systWeights12B[LEPDOWN]=1.;
	 //systWeights12B[ISOUP]=1.;
         //systWeights12B[ISODOWN]=1.;
	 //systWeights12B[TRIGUP]=1.;
         //systWeights12B[TRIGDOWN]=1.;

	 systWeights12BL[NOSYST]=1.;
	 systWeights12BL[BTAGUP]=1.;
	 systWeights12BL[BTAGDOWN]=1.;
	 systWeights12BL[MISTAGUP]=1.;
	 systWeights12BL[MISTAGDOWN]=1.;
	 systWeights12BL[PUUP]=1.;
	 systWeights12BL[PUDOWN]=1.;
	 systWeights12BL[LEPUP]=1.;
         systWeights12BL[LEPDOWN]=1.;
	 //systWeights12BL[ISOUP]=1.;
         //systWeights12BL[ISODOWN]=1.;
	 //systWeights12BL[TRIGUP]=1.;
         //systWeights12BL[TRIGDOWN]=1.;

	 systWeights2B[NOSYST]=1.;
	 systWeights2B[BTAGUP]=1.;
	 systWeights2B[BTAGDOWN]=1.;
	 systWeights2B[MISTAGUP]=1.;
	 systWeights2B[MISTAGDOWN]=1.;
	 systWeights2B[PUUP]=1.;
	 systWeights2B[PUDOWN]=1.;
	 systWeights2B[LEPUP]=1.;
         systWeights2B[LEPDOWN]=1.;
	 //systWeights2B[ISOUP]=1.;
         //systWeights2B[ISODOWN]=1.;
	 //systWeights2B[TRIGUP]=1.;
         //systWeights2B[TRIGDOWN]=1.;

	 systWeightsZeroB[NOSYST]=1.;
	 systWeightsZeroB[BTAGUP]=1.;
	 systWeightsZeroB[BTAGDOWN]=1.;
	 systWeightsZeroB[MISTAGUP]=1.;
	 systWeightsZeroB[MISTAGDOWN]=1.;
	 systWeightsZeroB[PUUP]=1.;
	 systWeightsZeroB[PUDOWN]=1.;
	 systWeightsZeroB[LEPUP]=1.;
         systWeightsZeroB[LEPDOWN]=1.;
	 //systWeightsZeroB[ISOUP]=1.;
         //systWeightsZeroB[ISODOWN]=1.;
	 //systWeightsZeroB[TRIGUP]=1.;
         //systWeightsZeroB[TRIGDOWN]=1.;

       }

       if(isData=="MC"){

//	 systWeights[NOSYST]=bWeight2;
//	 systWeights[BTAGUP]=bWeight2BTagUp;
//	 systWeights[BTAGDOWN]=bWeight2BTagDown;
//	 systWeights[MISTAGUP]=bWeight2MisTagUp;
//	 systWeights[MISTAGDOWN]=bWeight2MisTagDown;
//	 systWeights[PUUP]=1.;
//	 systWeights[PUDOWN]=1.;
	 //const char* QCDlabel = "QCD";
	 //const char* DMttlabel = "DMtt";
	 
	 double puUpFact=(LumiWeightsUp_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));
	 double puDownFact=(LumiWeightsDown_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));

	 if(numTrueInt>49){
	   cout << " --> numTrueInt very high!!" << endl;
	   puUpFact =0;
	   puDownFact=0;
	 }

	 systZero.setWeight(0,1.);
	 systZero.setWeight("btagUp",1.);
	 systZero.setWeight("btagDown",1.);
	 systZero.setWeight("mistagUp",1.);
	 systZero.setWeight("mistagDown",1.);
	 systZero.setWeight("puUp",1.);
	 systZero.setWeight("puDown",1.);
 	 systZero.setWeight("lepUp",1.);
	 systZero.setWeight("lepDown",1.);
	 //systZero.setWeight("isoUp",1.);
	 //systZero.setWeight("isoDown",1.);
	 //systZero.setWeight("trigUp",1.);
	 //systZero.setWeight("trigDown",1.);

	   //	   cout << " number of MC leptons " << NMCLeptons<< " weight 0 lep " << float(NMCLeptons==0)<< " weight 1 lep " <<float(NMCLeptons==1)<< " weight 2 lep " <<float(NMCLeptons==2)<<endl;

	 if(addTTSplit){
	   wcats[0]=1.0;
	   wcats[1]=double(NMCLeptons==0);
	   wcats[2]=double(NMCLeptons==1);
	   wcats[3]=double(NMCLeptons==2);
	   systZero.setWCats(wcats);
	   //	   systZero.setWeight("0lep",float(NMCLeptons==0),true);	   systZero.setWeight("1lep",float(NMCLeptons==1),true);	   systZero.setWeight("2lep",float(NMCLeptons==2),true);
	 }
	 
	 if(addPDF)systZero.setPDFWeights(w_pdfs,nPDF,w_zero,true);
	 if(addQ2)systZero.setQ2Weights(w_q2up,w_q2down,w_zero,true);
	 if(addTopPt)systZero.setTWeight(w_top,topWeight,true);
	 if(addVHF)systZero.setVHFWeight(vhf,true,shiftval);
	 if(addWZNLO){
	   systZero.setkFact("QCDRen",kfact_qcd,kfact_qcdrenUp,kfact_qcdrenDown);
	   systZero.setkFact("QCDFac",kfact_qcd,kfact_qcdfacUp,kfact_qcdfacDown);
	   systZero.setkFact("EWK",kfact_ewk,kfact_ewkUp,kfact_ewkDown);
	 }


	 syst12B.copySysts(systZero);
	 syst12B.setWeight(0,bWeight12);
	 syst12B.setWeight("btagUp",bWeight12BTagUp);
	 syst12B.setWeight("btagDown",bWeight12BTagDown);
	 syst12B.setWeight("mistagUp",bWeight12MisTagUp);
	 syst12B.setWeight("mistagDown",bWeight12MisTagDown);
	 syst12B.setWeight("puUp",bWeight12 * puUpFact);
	 syst12B.setWeight("puDown",bWeight12 * puDownFact);
	 syst12B.setWeight("lepUp",bWeight12 * elIDweightUp * muIDweightUp * elID2weightUp * muID2weightUp* muIsoweightUp * muIso2weightUp* elTrigweightUp * muTrigweightUp);
	 syst12B.setWeight("lepDown",bWeight12 * elIDweightDown * muIDweightDown * elID2weightDown * muID2weightDown* muIsoweightDown * muIso2weightDown* elTrigweightDown * muTrigweightDown);
	 //syst12B.setWeight("lepidUp",bWeight12 * elIDweightUp * muIDweightUp * elID2weightUp * muID2weightUp);
	 //syst12B.setWeight("lepidDown",bWeight12 * elIDweightDown * muIDweightDown * elID2weightDown * muID2weightDown);
	 //syst12B.setWeight("isoUp",bWeight12 * muIsoweightUp * muIso2weightUp);
	 //syst12B.setWeight("isoDown",bWeight12 * muIsoweightDown * muIso2weightDown);
	 //syst12B.setWeight("trigUp",bWeight12 * elTrigweightUp * muTrigweightUp);
	 //syst12B.setWeight("trigDown",bWeight12 * elTrigweightDown * muTrigweightDown);
	 if(addTTSplit){
	   syst12B.setWCats(wcats);
	   //	   syst12B.setWeight("0lep",float(NMCLeptons==0),true);	   syst12B.setWeight("1lep",float(NMCLeptons==1),true);	   syst12B.setWeight("2lep",float(NMCLeptons==2),true);
	 }
	 if(addPDF)syst12B.setPDFWeights(w_pdfs,nPDF,w_zero,true);
	 if(addQ2)syst12B.setQ2Weights(w_q2up,w_q2down,w_zero,true);
	 if(addTopPt)syst12B.setTWeight(w_top,topWeight,true);
	 if(addVHF)syst12B.setVHFWeight(vhf,true,shiftval);
	 if(addWZNLO){
	   syst12B.setkFact("QCDRen",kfact_qcd,kfact_qcdrenUp,kfact_qcdrenDown);
	   syst12B.setkFact("QCDFac",kfact_qcd,kfact_qcdfacUp,kfact_qcdfacDown);
	   syst12B.setkFact("EWK",kfact_ewk,kfact_ewkUp,kfact_ewkDown);
	 }

	 syst2B.copySysts(systZero);
	 syst2B.setWeight(0,bWeight2);
	 syst2B.setWeight("btagUp",bWeight2BTagUp);
	 syst2B.setWeight("btagDown",bWeight2BTagDown);
	 syst2B.setWeight("mistagUp",bWeight2MisTagUp);
	 syst2B.setWeight("mistagDown",bWeight2MisTagDown);
	 syst2B.setWeight("puUp",bWeight2 * puUpFact);
	 syst2B.setWeight("puDown",bWeight2 * puDownFact);
	 syst2B.setWeight("lepUp",bWeight2 * elIDweightUp * muIDweightUp * elID2weightUp * muID2weightUp* muIsoweightUp * muIso2weightUp* elTrigweightUp * muTrigweightUp);
	 syst2B.setWeight("lepDown",bWeight2 * elIDweightDown * muIDweightDown * elID2weightDown * muID2weightDown* muIsoweightDown * muIso2weightDown* elTrigweightDown * muTrigweightDown);
    	 //syst2B.setWeight("isoUp",bWeight2 * muIsoweightUp * muIso2weightUp);
	 //syst2B.setWeight("isoDown",bWeight2 * muIsoweightDown * muIso2weightDown);
	 //syst2B.setWeight("trigUp",bWeight2 * elTrigweightUp * muTrigweightUp);
	 //syst2B.setWeight("trigDown",bWeight2 * elTrigweightDown * muTrigweightDown);
	 if(addTTSplit){
	   syst2B.setWCats(wcats);
	   //	   syst2B.setWeight("0lep",float(NMCLeptons==0),true);	   syst2B.setWeight("1lep",float(NMCLeptons==1),true);	   syst2B.setWeight("2lep",float(NMCLeptons==2),true);
	 }
	 if(addPDF)syst2B.setPDFWeights(w_pdfs,nPDF,w_zero,true);
	 if(addQ2)syst2B.setQ2Weights(w_q2up,w_q2down,w_zero,true);
	 if(addTopPt)syst2B.setTWeight(w_top,topWeight,true);
	 if(addVHF)syst2B.setVHFWeight(vhf,true,shiftval);
	 if(addWZNLO){
	   syst2B.setkFact("QCDRen",kfact_qcd,kfact_qcdrenUp,kfact_qcdrenDown);
	   syst2B.setkFact("QCDFac",kfact_qcd,kfact_qcdfacUp,kfact_qcdfacDown);
	   syst2B.setkFact("EWK",kfact_ewk,kfact_ewkUp,kfact_ewkDown);
	 }

	 syst12BL.copySysts(systZero);
	 syst12BL.setWeight(0,bWeight12L);
	 syst12BL.setWeight("btagUp",bWeight12LBTagUp);
	 syst12BL.setWeight("btagDown",bWeight12LBTagDown);
	 syst12BL.setWeight("mistagUp",bWeight12LMisTagUp);
	 syst12BL.setWeight("mistagDown",bWeight12LMisTagDown);
	 syst12BL.setWeight("puUp",bWeight12L * puUpFact);
	 syst12BL.setWeight("puDown",bWeight12L * puDownFact);
	 syst12BL.setWeight("lepUp",bWeight12L * elIDweightUp * muIDweightUp * elID2weightUp * muID2weightUp* muIsoweightUp * muIso2weightUp* elTrigweightUp * muTrigweightUp);
       syst12BL.setWeight("lepDown",bWeight12L * elIDweightDown * muIDweightDown * elID2weightDown * muID2weightDown* muIsoweightDown * muIso2weightDown* elTrigweightDown * muTrigweightDown);
	 //syst12BL.setWeight("trigUp",bWeight12L * elTrigweightUp * muTrigweightUp);
	 //syst12BL.setWeight("trigDown",bWeight12L * elTrigweightDown * muTrigweightDown);
	 if(addTTSplit){
	   syst12BL.setWCats(wcats);
	   //	   syst12BL.setWeight("0lep",float(NMCLeptons==0),true);	   syst12BL.setWeight("1lep",float(NMCLeptons==1),true);	   syst12BL.setWeight("2lep",float(NMCLeptons==2),true);
	 }
	 if(addPDF)syst12BL.setPDFWeights(w_pdfs,nPDF,w_zero,true);
	 if(addQ2)syst12BL.setQ2Weights(w_q2up,w_q2down,w_zero,true);
	 if(addTopPt)syst12BL.setTWeight(w_top,topWeight,true);
	 if(addVHF)syst12BL.setVHFWeight(vhf,true,shiftval);
	 if(addWZNLO){
	   syst12BL.setkFact("QCDRen",kfact_qcd,kfact_qcdrenUp,kfact_qcdrenDown);
	   syst12BL.setkFact("QCDFac",kfact_qcd,kfact_qcdfacUp,kfact_qcdfacDown);
	   syst12BL.setkFact("EWK",kfact_ewk,kfact_ewkUp,kfact_ewkDown);
	 }

	 syst0B.copySysts(systZero);
	 syst0B.setWeight(0,bWeightZero);
	 syst0B.setWeight("btagUp",bWeightZeroBTagUp);
	 syst0B.setWeight("btagDown",bWeightZeroBTagDown);
	 syst0B.setWeight("mistagUp",bWeightZeroMisTagUp);
	 syst0B.setWeight("mistagDown",bWeightZeroMisTagDown);
	 syst0B.setWeight("puUp",bWeightZero * puUpFact);
	 syst0B.setWeight("puDown",bWeightZero * puDownFact);
	 syst0B.setWeight("lepUp",bWeightZero * elIDweightUp * muIDweightUp * elID2weightUp * muID2weightUp* muIsoweightUp * muIso2weightUp* elTrigweightUp * muTrigweightUp);
	 syst0B.setWeight("lepDown",bWeightZero * elIDweightDown * muIDweightDown * elID2weightDown * muID2weightDown* muIsoweightDown * muIso2weightDown* elTrigweightDown * muTrigweightDown);
	 //syst0B.setWeight("isoUp",bWeightZero * muIsoweightUp * muIso2weightUp);
	 //syst0B.setWeight("isoDown",bWeightZero * muIsoweightDown * muIso2weightDown);
	 //syst0B.setWeight("trigUp",bWeightZero * elTrigweightUp * muTrigweightUp);
	 //syst0B.setWeight("trigDown",bWeightZero * elTrigweightDown * muTrigweightDown);
	 if(addTTSplit){
	   syst0B.setWCats(wcats);
	   //	   syst0B.setWeight("0lep",float(NMCLeptons==0),true);	   syst0B.setWeight("1lep",float(NMCLeptons==1),true);	   syst0B.setWeight("2lep",float(NMCLeptons==2),true);
	 }
	 if(addPDF)syst0B.setPDFWeights(w_pdfs,nPDF,w_zero,true);
	 if(addQ2)syst0B.setQ2Weights(w_q2up,w_q2down,w_zero,true);
	 if(addTopPt)syst0B.setTWeight(w_top,topWeight,true);
  	 if(addVHF)syst0B.setVHFWeight(vhf,true,shiftval);
	 if(addWZNLO){
	   syst0B.setkFact("QCDRen",kfact_qcd,kfact_qcdrenUp,kfact_qcdrenDown);
	   syst0B.setkFact("QCDFac",kfact_qcd,kfact_qcdfacUp,kfact_qcdfacDown);
	   syst0B.setkFact("EWK",kfact_ewk,kfact_ewkUp,kfact_ewkDown);
	 }
 

	 systWeights12B[PUUP]= puUpFact * bWeight12;
	 systWeights12B[PUDOWN]=puDownFact * bWeight12;
	 systWeights12B[LEPUP]= bWeight12 * elIDweightUp * muIDweightUp * elID2weightUp * muID2weightUp* muIsoweightUp * muIso2weightUp* elTrigweightUp * muTrigweightUp;
         systWeights12B[LEPDOWN]= bWeight12 * elIDweightDown * muIDweightDown * elID2weightDown * muID2weightDown* muIsoweightDown * muIso2weightDown* elTrigweightDown * muTrigweightDown;
	 //systWeights12B[ISOUP]= bWeight12 * muIsoweightUp * muIso2weightUp;
         //systWeights12B[ISODOWN]= bWeight12 * muIsoweightDown * muIso2weightDown;
	 //systWeights12B[TRIGUP]= bWeight12 * elTrigweightUp * muTrigweightUp;
         //systWeights12B[TRIGDOWN]= bWeight12 * elTrigweightDown * muTrigweightDown;
	 systWeights12B[NOSYST]=bWeight12;
	 systWeights12B[BTAGUP]=bWeight12BTagUp;
	 systWeights12B[BTAGDOWN]=bWeight12BTagDown;
	 systWeights12B[MISTAGUP]=bWeight12MisTagUp;
	 systWeights12B[MISTAGDOWN]=bWeight12MisTagDown;

	 systWeights12BL[PUUP]=puUpFact * bWeight12L;
	 systWeights12BL[PUDOWN]=puDownFact * bWeight12L;
	 systWeights12BL[LEPUP]= bWeight12L * elIDweightUp * muIDweightUp * elID2weightUp * muID2weightUp* muIsoweightUp * muIso2weightUp* elTrigweightUp * muTrigweightUp;
         systWeights12BL[LEPDOWN]= bWeight12L * elIDweightDown * muIDweightDown * elID2weightDown * muID2weightDown* muIsoweightDown * muIso2weightDown* elTrigweightDown * muTrigweightDown;
	 //systWeights12BL[ISOUP]= bWeight12L * muIsoweightUp * muIso2weightUp;
         //systWeights12BL[ISODOWN]= bWeight12L * muIsoweightDown * muIso2weightDown;
	 //systWeights12BL[TRIGUP]= bWeight12L * elTrigweightUp * muTrigweightUp;
         //systWeights12BL[TRIGDOWN]= bWeight12L * elTrigweightDown * muTrigweightDown;
	 systWeights12BL[NOSYST]=bWeight12L;
	 systWeights12BL[BTAGUP]=bWeight12LBTagUp;
	 systWeights12BL[BTAGDOWN]=bWeight12LBTagDown;
	 systWeights12BL[MISTAGUP]=bWeight12LMisTagUp;
	 systWeights12BL[MISTAGDOWN]=bWeight12LMisTagDown;
	 
	 systWeights2B[PUUP]=puUpFact * bWeight2;
	 systWeights2B[PUDOWN]=puDownFact * bWeight2;
	 systWeights2B[LEPUP]= bWeight2 * elIDweightUp * muIDweightUp * elID2weightUp * muID2weightUp* muIsoweightUp * muIso2weightUp* elTrigweightUp * muTrigweightUp;
         systWeights2B[LEPDOWN]= bWeight2 * elIDweightDown * muIDweightDown * elID2weightDown * muID2weightDown* muIsoweightDown * muIso2weightDown* elTrigweightDown * muTrigweightDown;
	 //systWeights2B[ISOUP]= bWeight2 * muIsoweightUp * muIso2weightUp;
         //systWeights2B[ISODOWN]= bWeight2 * muIsoweightDown * muIso2weightDown;
	 //systWeights2B[TRIGUP]= bWeight2 * elTrigweightUp * muTrigweightUp;
         //systWeights2B[TRIGDOWN]= bWeight2 * elTrigweightDown * muTrigweightDown;
       	 systWeights2B[NOSYST]=bWeight2;
	 systWeights2B[BTAGUP]=bWeight2BTagUp;
	 systWeights2B[BTAGDOWN]=bWeight2BTagDown;
	 systWeights2B[MISTAGUP]=bWeight2MisTagUp;
	 systWeights2B[MISTAGDOWN]=bWeight2MisTagDown;
	 
	 //bWeightZero = 1.,bWeightZeroBTagUp= 1., bWeightZeroBTagDown=1.0, bWeightZeroMisTagUp=1.0,bWeightZeroMisTagDown=1.0;
	 systWeightsZeroB[PUUP]=puUpFact * bWeightZero;
         systWeightsZeroB[PUDOWN]= puDownFact * bWeightZero;
      	 systWeightsZeroB[LEPUP]= bWeightZero * elIDweightUp * muIDweightUp * elID2weightUp * muID2weightUp* muIsoweightUp * muIso2weightUp* elTrigweightUp * muTrigweightUp;
         systWeightsZeroB[LEPDOWN]= bWeightZero * elIDweightDown * muIDweightDown * elID2weightDown * muID2weightDown* muIsoweightDown * muIso2weightDown* elTrigweightDown * muTrigweightDown;
      	 //systWeightsZeroB[ISOUP]= bWeightZero * muIsoweightUp * muIso2weightUp;
         //systWeightsZeroB[ISODOWN]= bWeightZero * muIsoweightDown * muIso2weightDown;
	 //systWeightsZeroB[TRIGUP]= bWeightZero * elTrigweightUp * muTrigweightUp;
         //systWeightsZeroB[TRIGDOWN]= bWeightZero * elTrigweightDown * muTrigweightDown;
      	 systWeightsZeroB[NOSYST]=bWeightZero;
	 systWeightsZeroB[BTAGUP]=bWeightZeroBTagUp;
	 systWeightsZeroB[BTAGDOWN]=bWeightZeroBTagDown;
	 systWeightsZeroB[MISTAGUP]=bWeightZeroMisTagUp;
	 systWeightsZeroB[MISTAGDOWN]=bWeightZeroMisTagDown;
       }
       
       systWeightsNoSyst[PUUP]=1.0;
       systWeightsNoSyst[PUDOWN]=1.0;
       systWeightsNoSyst[LEPUP]=1.0;
       systWeightsNoSyst[LEPDOWN]=1.0;
       //systWeightsNoSyst[ISOUP]=1.0;
       //systWeightsNoSyst[ISODOWN]=1.0;
       //systWeightsNoSyst[TRIGUP]=1.0;
       //systWeightsNoSyst[TRIGDOWN]=1.0;
       systWeightsNoSyst[NOSYST]=1.0;
       systWeightsNoSyst[BTAGUP]=1.0;
       systWeightsNoSyst[BTAGDOWN]=1.0;
       systWeightsNoSyst[MISTAGUP]=1.0;
       systWeightsNoSyst[MISTAGDOWN]=1.0;
 	 
       n_trig += w;
       
       //       cout<< " prefirstfill "<<endl;
       syst0B.fillHistogramsSysts(h_nPV,nPV,1.0,systWeightsZeroB);
       //       cout<< " postfirstfill "<<endl;
       syst0B.fillHistogramsSysts(h_nGoodPV,nPV,1.0,systWeightsZeroB);
       
       syst0B.fillHistogramsSysts(h_nPV_w,nPV,w,systWeightsZeroB);
       syst0B.fillHistogramsSysts(h_nGoodPV_w,nPV,w,systWeightsZeroB);

       syst0B.fillHistogramsSysts(h_nJets,nJets,w,systWeightsZeroB);       
       syst0B.fillHistogramsSysts(h_nbJets,nCSVJets,w,systWeightsZeroB);


       //fillHistogramsSysts(h_topHadMass_allCutsMS, 170.0, w, systWeights);


       if(nTightElectrons != nEl)cout << "warning! problem with tight el"<<endl;
       if(nTightMuons != nMu)cout << "warning! problem with tight mu"<<endl;
       nTightElectrons = nEl;

       //nLooseMuons= nLooseMu;
       nTightMuons = nMu;
 
       float dphiBB =-999.;
       float dphiBB_MET =-999.;
       //float dphiTT_MET =-999.;
       float detaBB =-999.;
       // float detaTT_MET =-999.;
       //float dphiTL =-999.;
       //float detaTL =-999.;
       vector<float> jetsPhi;
       vector<TLorentzVector> bjets, jets_nob;

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

       //       vector<bool> passDR_vect;
       for (int j = 0; j <maxJetLoop;++j){
	 TLorentzVector jet;
	 if(jetIsCSVM[j]<=0.){
	   jet.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
	   jets_nob.push_back(jet);
	 }
	 jetsPhi.push_back(jetPhi[j]);
	 d_njets+=1;

	 if(jetsPhi.size()==1)syst0B.fillHistogramsSysts(h_jet1Pt,jetPt[j],w,systWeightsZeroB);
	 if(jetsPhi.size()==2)syst0B.fillHistogramsSysts(h_jet2Pt,jetPt[j],w,systWeightsZeroB);
	 if(jetsPhi.size()==3)syst0B.fillHistogramsSysts(h_jet3Pt,jetPt[j],w,systWeightsZeroB);
	 
	 TLorentzVector bjet;
	 //	   if(jetIsTight[j] && jetIsCSVM[j] && abs(jetEta[j])<2.4){
	 if( jetIsCSVM[j] && abs(jetEta[j])<2.4){
	   bjet.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
	   bjets.push_back(bjet);
	   btag b;
	   b.vect = bjet;
	   b.csv = jetCSV[j];
	   bvects.push_back(b);
	   
	   syst0B.fillHistogramsSysts(h_bjetsPt,jetPt[j],w ,systWeightsZeroB);
	   if(bjets.size()==1)syst0B.fillHistogramsSysts(h_bjet1Pt,jetPt[j],w,systWeightsZeroB);
	   if(bjets.size()==2)syst0B.fillHistogramsSysts(h_bjet2Pt,jetPt[j],w,systWeightsZeroB);
	   if(bjets.size()==3)syst0B.fillHistogramsSysts(h_bjet3Pt,jetPt[j],w,systWeightsZeroB);
	 }
       	 //passDR_vect.push_back(passesDR);
       } //end of the for loop

       //assert(nCSVJets==bjets.size());
       
       nJets = jetsPhi.size();
       nCSVJets=bjets.size();
       
       if(nJets>2){
	 //	 cout<< " test1 "<<endl;
	 // vector<float> minDphi_6jMET, minDphi_10jMET;
	 float dphi1 = std::fabs(deltaPhi(metPhi[0],jetsPhi[0]));
	 float dphi2 = std::fabs(deltaPhi(metPhi[0],jetsPhi[1]));	
	 minDPhi =  std::min(dphi1,dphi2);

	 float dphi3 = 90;
	 float dphi4 = 90;	
	 float dphi5 = 90;
	 float dphi6 = 90;
	 float dphi7 = 90;
	 float dphi8 = 90;
	 float dphi9 = 90;
	 float dphi10 = 90;
 
	 if(nJets>2) dphi3 = std::fabs(deltaPhi(metPhi[0],jetsPhi[2]));
	 if(nJets>3) dphi4 = std::fabs(deltaPhi(metPhi[0],jetsPhi[3]));
	 if(nJets>4) dphi5 = std::fabs(deltaPhi(metPhi[0],jetsPhi[4]));
	 if(nJets>5) dphi6 = std::fabs(deltaPhi(metPhi[0],jetsPhi[5]));	
	 if(nJets>6) dphi7 = std::fabs(deltaPhi(metPhi[0],jetsPhi[6]));	
	 if(nJets>7) dphi8 = std::fabs(deltaPhi(metPhi[0],jetsPhi[7]));
	 if(nJets>8) dphi9 = std::fabs(deltaPhi(metPhi[0],jetsPhi[8]));	
	 if(nJets>9) dphi10 = std::fabs(deltaPhi(metPhi[0],jetsPhi[9]));	

	 vector<float> minDPhi_10jMET{dphi1, dphi2, dphi3, dphi4, dphi5, dphi6, dphi7, dphi8, dphi9, dphi10};
	 vector<float> minDPhi_6jMET{dphi1, dphi2, dphi3, dphi4, dphi5, dphi6};
	 
	 minDPhi_6j = *std::min_element(std::begin(minDPhi_6jMET), std::end(minDPhi_6jMET));
	 minDPhi_10j = *std::min_element(std::begin(minDPhi_10jMET), std::end(minDPhi_10jMET));
	 //std::cout << "min element new method: " <<  mDPhi_j6;
	 
	 // float min_1_2 = std::min(dphi1,dphi2);
	 // float min_3_4 = std::min(dphi3,dphi4);
	 // float min_5_6 = std::min(dphi5,dphi6);
	 // float min_34_56 = std::min(min_3_4,min_5_6);
	 // minDPhi_6j = std::min(min_1_2, min_34_56);
	 // std::cout << "min element old method: " <<  minDPhi_6j;

	 // assert(minDPhi_6j==mDPhi_j6);
	 
	 if(nCSVJets>1){
	   dphiBB = std::fabs(deltaPhi(bjets[0].Phi(),bjets[1].Phi()));
	   detaBB = std::fabs(bjets[0].Eta()-bjets[1].Eta());
	   std::sort(bvects.begin(), bvects.end(), by_csv());
	   dphiBB_MET = fabs(deltaPhi( (bvects[0].vect + bvects[1].vect).Phi(), metPhi[0] ) );
	 }
       }// end if njets > 2
       
       //int nLooseMuons(0);//, nVetoElectrons(0);
       //nTightElectrons =0;
       //nTightMuons = 0;

       bool semileptonic = ( (nTightElectrons+nTightMuons ==1) && (nLooseMuons+nVetoElectrons <= 1)) && channel == "semileptonic";
       bool fullhadronic = ( (nTightElectrons+nTightMuons ==0) &&  (nLooseMuons+nVetoElectrons ==0)) && channel == "fullhadronic";

       //for bkg estimation
       bool fullhadronic_1lep = (nTightElectrons+nTightMuons ==1) && (nLooseMuons+nVetoElectrons <= 1) && channel == "fullhadronic";
       bool fullhadronic_2lep = (nTightElectrons+nTightMuons ==2) && (nLooseMuons+nVetoElectrons <= 2) && channel == "fullhadronic";
       bool semileptonic_2lep = (nTightElectrons+nTightMuons ==2) && (nLooseMuons+nVetoElectrons <= 2) && channel == "semileptonic";

       // metNoLep is evaluating removing the 1st lepton in case only 1 lepton is present in the event 
       // or the 2nd lepton in case two or more leptons are present in the event, this can be used for tt(1l) tt(2l) CRs for example
       // metNoTwoLep is evaluating removing from met computation the momemtum of 2 leptons, and it's is used for Z+Jets CR in hadronic channel 
       float metNoLep = met;
       float metNoTwoLep = met;
       // float metNoLepPhi = metPhi[0];
       float newMetPx = met*TMath::Cos(metPhi[0]);
       float newMetPy = met*TMath::Sin(metPhi[0]);
       float newMetPx_2 = met*TMath::Cos(metPhi[0]);
       float newMetPy_2 = met*TMath::Sin(metPhi[0]);
       //  TLorentzVector lep1;
       // TLorentzVector lep2;
       float mll(0.);
       if(nTightElectrons+nTightMuons==1){
	 //	 lep1.SetPtEtaPhiE(lep1Pt,lep1Eta,lep1Phi,lep1E);
	 newMetPx += lep1.Px();
	 newMetPy += lep1.Py();
       }
       else if((nTightElectrons==2 or nTightMuons==2) and (nTightElectrons+nTightMuons==2)){
	 //	 lep1.SetPtEtaPhiE(lep1Pt,lep1Eta,lep1Phi,lep1E);
	 //	 lep2.SetPtEtaPhiE(lep2Pt,lep2Eta,lep2Phi,lep2E);
	 newMetPx += lep2.Px();
	 newMetPy += lep2.Py();
	 newMetPx_2 = newMetPx + lep1.Px();
	 newMetPy_2 = newMetPy + lep1.Py();
	 mll = (lep1 + lep2).M();
	 // I leptoni al momento sono vuoti.
	 //std::cout<<"Mll "<<mll<<std::endl;		
       }
       else if(nTightElectrons+nTightMuons>=2){
	 // lep1.SetPtEtaPhiE(lep1Pt,lep1Eta,lep1Phi,lep1E);
	 //	 lep2.SetPtEtaPhiE(lep2Pt,lep2Eta,lep2Phi,lep2E);
	 newMetPx += lep2.Px();
	 newMetPy += lep2.Py();
	 newMetPx_2 = newMetPx + lep1.Px();
	 newMetPy_2 = newMetPy + lep1.Py();
       }

       // if(newMetPx<0){
       // 	 if(newMetPy>0)metNoLepPhi = atan(newMetPy/newMetPx)+3.141592;
       // 	 if(newMetPy<0)metNoLepPhi = atan(newMetPy/newMetPx)-3.141592;
       // }
       // else  metNoLepPhi = (atan(newMetPy/newMetPx));
       
       metNoLep = TMath::Sqrt( (pow(newMetPx,2) + pow(newMetPy,2)));
       metNoTwoLep = TMath::Sqrt( (pow(newMetPx_2,2) + pow(newMetPy_2,2)));
       //std::cout<<"met "<<metNoTwoLep<<std::endl;

       // =========================================================================
       // KINEMATIC FITTER + MVA FOR BOTH singleLepton and full hadronic, SR and CR
       // =========================================================================

       
       if( ( channel == "fullhadronic" && nJets>3) or ( channel == "semileptonic" && nJets>2 && (nTightElectrons+nTightMuons)>0) ){
	 float  csv_l, csv_j, csv_k, jet1csv, jet2csv, jet3csv, jet1qgid, jet2qgid, mva_;
	 vector<float> csv, topMassPreFit, topPtPreFit, topEtaPreFit, topPhiPreFit, WMassPreFit, WMPhiPreFit , TMPhiPreFit, BMPhiPreFit, WBPhiPreFit;
	 vector<float> topMassPostFit, topPtPostFit, topEtaPostFit, topPhiPostFit, WMassPostFit,  WMPhiPostFit , TMPhiPostFit, BMPhiPostFit, WBPhiPostFit, LMPhiPostFit ;
	 vector<int> idx;
	 int idxCSV=-1;
	 TLorentzVector jet1, jet2, jet3;
	 KinematicFitter fitter;


	 int countcombos=0;
	 for (int l = 0; l <maxJetLoop;++l){
	   for (int j = l+1; j <maxJetLoop;++j){
	     for (int k = j+1; k <maxJetLoop;++k){
	       //	       bool pass_l = (jetPt[l]>30. && fabs(jetEta[l])<4. && jetPassID[l]>0. && passDR_vect[l]>0.);
	       //	       bool pass_j = (jetPt[j]>30. && fabs(jetEta[j])<4. && jetPassID[j]>0. &&  passDR_vect[j]>0.);
	       //	       bool pass_k = (jetPt[k]>30. && fabs(jetEta[k])<4. && jetPassID[k]>0. &&  passDR_vect[k]>0.);

	       bool pass_l = true;
	       bool pass_j = true;
	       bool pass_k = true;
	       if(pass_l & pass_j && pass_k){
	       //	       if(jetIsTight[l] and jetPt[l]>30. and jetIsTight[j] and jetPt[j]>30. and jetIsTight[k] and jetPt[k]>30.) {

                 csv.clear(); idx.clear();
		 TopFitResults topKinFit;
		 //std::cout<<"Performing kinematic fitter"<<std::endl;
		 csv_l =  jetCSV[l]; csv_j = jetCSV[j]; csv_k =  jetCSV[k];
		 csv.push_back(!TMath::IsNaN(csv_l) ? csv_l : -10. );
		 csv.push_back(!TMath::IsNaN(csv_j) ? csv_j : -10.);
		 csv.push_back(!TMath::IsNaN(csv_k) ? csv_k : -10.);
		 idx.push_back(l);  idx.push_back(j); idx.push_back(k);
	

		 idxCSV = distance(csv.begin(), max_element(csv.begin(), csv.end()) );

		 //std::cout<<"B index as max_element "<<distance(csv.begin(), max_element(csv.begin(), csv.end()) );

		 // Setting 3rd jet-> the b jet
		 jet3.SetPtEtaPhiE(jetPt[idx[idxCSV]],jetEta[idx[idxCSV]],jetPhi[idx[idxCSV]],jetE[idx[idxCSV]]);
		 jet3csv = jetCSV[idx[idxCSV]];
		 //idxB.push_back(idx[idxCSV]);
		 topKinFit.idxB = idx[idxCSV];
		 idx.erase(idx.begin() + idxCSV);
		 jet1.SetPtEtaPhiE(jetPt[idx[0]],jetEta[idx[0]],jetPhi[idx[0]],jetE[idx[0]]);
		 jet1csv = jetCSV[idx[0]]; jet1qgid = jetQGL[idx[0]];

		 jet2.SetPtEtaPhiE(jetPt[idx[1]],jetEta[idx[1]],jetPhi[idx[1]],jetE[idx[1]]);
		 jet2csv = jetCSV[idx[1]]; jet2qgid = jetQGL[idx[1]];

		 // Patch temporanea
		 if(TMath::IsNaN(jetCSV[idx[0]]))jet1csv = -10;
		 if(TMath::IsNaN(jetCSV[idx[1]]))jet2csv = -10;

		 //idxJ1.push_back(idx[0]);idxJ2.push_back(idx[1]);
		 topKinFit.idxJ1 = idx[0];topKinFit.idxJ2 = idx[1];
		 // initialize object that stores fit results
		 FitResults fitres;
		 fitres.converged = false;
		 fitres.prob      = 0.;
		 fitres.chisq     = 999.;
		 fitres.cost      = 999.;
		 fitres.fitmass   = 0.;
		 fitres.fitmassW  = 0.;
		 
		 //                 std::cout << "-- jets ---" << std::endl;
		 //                 std::cout << jet1.Pt() << " " << jet2.Pt() << " " << jet3.Pt() << "Bjet csv: "<<jet3csv<<std::endl;
		 TopCandidate::TopCandidateParticle wjet1(jet1, std::string("unmatched"), 3, 0);
		 TopCandidate::TopCandidateParticle wjet2(jet2, std::string("unmatched"), 3, 0);
		 TopCandidate::TopCandidateParticle bjet (jet3, std::string("unmatched"), 3, 0);
		 
		 // The constructor, TopCandidate(j1,j2,j3), assumes "j3" corresponds to the b-jet while "j1" and "j2" are the W-jets
		 TopCandidate combo(wjet1, wjet2, bjet);
		 //Prefit
		 //std::cout<< "- Indexes as they are: "<<std::endl;
		 //std::cout<< "- " << i << ", " << j << ", " << k << " : "<<std::endl;
		 //std::cout<< "- Indexes as I take them: "<<std::endl;
		 //std::cout<< "- " <<  topKinFit.idxB << ", " <<  topKinFit.idxJ1 << ", " <<  topKinFit.idxJ2 << " : "<<std::endl;
		 //std::cout<<combo.topvec.M()<<" "<<combo.topvec.Pt()<<" "<<combo.topvec.Phi()<<" "<<combo.topvec.Eta()<<std::endl;


		 //==========PRE FIT===================
		 topKinFit.topPreFit = combo.topvec;
		 topKinFit.WPreFit = combo.Wvec;
		 
		 topKinFit.WMPhiPreFit = deltaPhi(combo.Wvec.Phi(), met);
		 topKinFit.TMPhiPreFit = deltaPhi(combo.topvec.Phi(), met);
		 topKinFit.BMPhiPreFit = deltaPhi(combo.particles[2].vec.Phi(), met);
		 topKinFit.WBPhiPreFit = deltaPhi(combo.Wvec.Phi(), combo.particles[2].vec.Phi());
		 
		 // FIT
		 //		 combo.reset();
		 //		 fitter.fit(combo, fitres);

		 //==========POST FIT===================
		 topKinFit.topPostFit = combo.topvec;
		 topKinFit.WPostFit = combo.Wvec;
		 
		 topKinFit.WMPhiPostFit = deltaPhi(combo.Wvec.Phi(), met);
		 topKinFit.TMPhiPostFit = deltaPhi(combo.topvec.Phi(), met);
		 topKinFit.BMPhiPostFit = deltaPhi(combo.particles[2].vec.Phi(), met);
		 topKinFit.WBPhiPostFit = deltaPhi(combo.Wvec.Phi(), combo.particles[2].vec.Phi());

		 // Set the inputs for the MVA
		 bdt_qgid1   = jet1qgid;
		 bdt_qgid2   = jet2qgid;
		 bdt_dphij1b = fabs(jet1.DeltaPhi(jet3));
		 bdt_dphij2b = fabs(jet2.DeltaPhi(jet3));
		 bdt_drj1b   = jet1.DeltaR(jet3);
		 bdt_drj2b   = jet2.DeltaR(jet3);
		 bdt_bjcsv   = jet3csv;
		 bdt_jet1csv = jet1csv;
		 bdt_jet2csv = jet2csv;
		 bdt_prob    = fitres.prob;
		 
		 // =======Filling up TopFitResults===============
		 topKinFit.bdt_qgid1   = jet1qgid;
		 topKinFit.bdt_qgid2   = jet2qgid;
		 topKinFit.bdt_dphij1b = fabs(jet1.DeltaPhi(jet3));
		 topKinFit.bdt_dphij2b = fabs(jet2.DeltaPhi(jet3));
		 topKinFit.bdt_drj1b   = jet1.DeltaR(jet3);
		 topKinFit.bdt_drj2b   = jet2.DeltaR(jet3);
		 topKinFit.bdt_bjcsv   = jet3csv;
		 topKinFit.bdt_jet1csv = jet1csv;
		 topKinFit.bdt_jet2csv = jet2csv;
		 topKinFit.bdt_prob    = fitres.prob;
		 
		 //
		 // Compute the MVA value
		 //

		 if(TMath::IsNaN(bdt_prob))bdt_prob=0;
		 if(TMath::IsNaN(bdt_bjcsv))callme();

		 mva_ = res_topmvaReader.EvaluateMVA("BDTG");
		 mva.push_back(mva_);
		 
		 if(mva_ == -999){
		   std::cout << "-- jets ---" << std::endl;
		   std::cout << jetPt[l] << " " << jetPt[j] << " " << jetPt[k] << std::endl;
		   std::cout << "-- csv ---" << std::endl;
		   std::cout <<jetCSV[l] << " " << jetCSV[j] << " " << jetCSV[k] <<std::endl;}

		 topKinFit.mva = mva_;
		 //topFitRes.push_back(topKinFit);
		 topFitResFull.push_back(topKinFit);
		 //std::cout << "MVA value = " << mva << std::endl;
		 countcombos++;
		 if(countcombos>0)break;
		 //		 if(l==0 && j==1 && k==2)break;
	       } // end of tight selection on jets
	       if(countcombos>0)break;
	     }// loop over jets -k 
	     if(countcombos>0)break;
	   }// loop over jets -j 
	   if(countcombos>0)break;
	 }// loop over jets -i 
      
	 //std::cout << "-- Raw MVAs --- size "  << topFitRes.size()<< std::endl;
	 //for ( const TopFitResults& res : topFitRes) {
	 //  std::cout << "idxs: b = " << res.idxB << " j1 = " << res.idxJ1 << " j2 = " << res.idxJ2 << " -> mva = " << res.mva << std::endl;
	 //}    
	 //std::cout << "--------" << std::endl;
	 
       // Sorting the vector topFitRes according to the value of mva
	 // sort(topFitRes.begin(), topFitRes.end(), 
	 //      [](const TopFitResults & a, const TopFitResults & b) -> bool
	 //      { 
	 // 	return a.mva > b.mva; 
	 //      });
	 //	 cout<< " test2 "<<endl;

	 sort(topFitResFull.begin(), topFitResFull.end(), 
	      [](const TopFitResults & a, const TopFitResults & b) -> bool
	      { 
		return a.mva > b.mva; 
	      });

	 
	 const TopFitResults& res0 = topFitResFull[0];	 
	 std::set<int> idx0 = { res0.idxB, res0.idxJ1, res0.idxJ2 };
	 topFitRes.push_back(res0);
	 
	 for (auto it = (topFitResFull.begin()+1); it!=topFitResFull.end();++it){
	   std::set<int> idxs  = { it->idxB, it->idxJ1, it->idxJ2 };
	   std::set<int> common;
	   std::set_intersection(idx0.begin(), idx0.end(), 
				 idxs.begin(), idxs.end(), 
				 std::inserter(common, common.begin()));
	   if(common.empty())topFitRes.push_back(*it);
	 }

	 /*
	 std::cout << "-- Original tri-jet combinations--- size "  << topFitResFull.size()<< std::endl;
	 for ( const TopFitResults& res : topFitResFull) {
	   std::cout << "idxs: b = " << res.idxB << " j1 = " << res.idxJ1 << " j2 = " << res.idxJ2 << " -> mva = " << res.mva << std::endl;
	 }    
	 std::cout << "--------" << std::endl;
	 
	 std::cout << "-- Cleand tri-jet combinations--- size "  << topFitRes.size()<< std::endl;
	 
	 for ( const TopFitResults& res : topFitRes) {
	   std::cout << "idxs: b = " << res.idxB << " j1 = " << res.idxJ1 << " j2 = " << res.idxJ2 << " -> mva = " << res.mva << std::endl;
	 }    
	 std::cout << "--------" << std::endl;
	 */

	 // Reference to first element


	 // std::remove_if( (topFitRes.begin()+1), 
	 // 	      topFitRes.end(),
	 // 	      []( TopFitResults& r ){
	 // 		std::set<int> common;
	 // 		std::set<int> idxs  = { r.idxB, r.idxJ1, r.idxJ2 };
	 // 		std::set_intersection(idx0.begin(), idx0.end(), 
	 // 				      idxs.begin(), idxs.end(), 
	 // 				      std::inserter(common, common.begin()));
	 // 		return !common.empty();
	 // 	      } 
	 // )

	 
	 //std::cout << "-- Sorted MVAs --- size "  << topFitRes.size()<< std::endl;
	 //for ( const TopFitResults& res : topFitRes) {
	 //  std::cout << "idxs: b = " << res.idxB << " j1 = " << res.idxJ1 << " j2 = " << res.idxJ2 << " -> mva = " << res.mva << std::endl;
	 //}    
	 //std::cout << "--------" << std::endl;

       }//END KINEMATIC FITTER + MVA
        
       // ***************************************************
       // *************BOOSTED CATEGORIES********************
       // ***************************************************
       //       cout<< " test3 "<<endl;
       if(semileptonic && met > 300){
	 if(nType1==1 && nType2==0){
	   n_boost1Res += w;
	   syst0B.fillHistogramsSysts(h_metFinal1Res,met,w,systWeightsZeroB);
	 }	     
	 else if(nType1==0 && nType2==1){
	   n_boost2Res += w;
	   syst0B.fillHistogramsSysts(h_metFinal2Res,met,w,systWeightsZeroB);
	 }
	 else if(nType1==0 && nType2==0) n_boostFullRes += w;
       }
       if(fullhadronic && met > 300){
	 if(nType1>1 && nCSVJets>0){
	   n_boost11 += w;
	   syst12B.fillHistogramsSysts(h_metFinal11,met,w,systWeights12B);
	 }
	 else if(nType2>1){
	   n_boost22 += w;
	   syst0B.fillHistogramsSysts(h_metFinal22,met,w,systWeightsZeroB);
	 }
	 else if(nType1==1 && nType2==1){
	   n_boost12 += w;
	   syst0B.fillHistogramsSysts(h_metFinal12,met,w,systWeightsZeroB);
	 }
	 else if(nType1==1 && nType2==0){
	   n_boost1Res += w;
	   syst0B.fillHistogramsSysts(h_metFinal1Res,met,w,systWeightsZeroB);
	 }
	 else if(nType1==0 && nType2==1){
	   n_boost2Res += w;
	   syst0B.fillHistogramsSysts(h_metFinal2Res,met,w,systWeightsZeroB);
	 }
	 else if(nType1==0 && nType2==0){
	   n_boostFullRes += w;
	 }
       }                               
       
       if( (nJets-nCSVJets)<2 ){
         TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
         float phi_lmet = fabs(deltaPhi(lep1.Phi(), metPhi[0]) );
         mt = sqrt(2* lep1.Pt() * met* ( 1- cos(phi_lmet)));
         Mt2Com_bisect *Mt2cal = new Mt2Com_bisect();
         mt2w = Mt2cal->calculateMT2w(jets_nob,bjets,lep1, met_,"MT2w");
       }

       
       // ***************************************************
       // *************FULLHADRONIC CHANNEL******************
       // ***************************************************
       
       
       //for bkg estimation
       float mt_full1 = 0;
       
       if(fullhadronic_1lep or semileptonic){
	 mt_full1 = sqrt(2* lep1.Pt() * metPt[0] * ( 1- cos(fabs(deltaPhi(lep1.Phi(), metPhi[0])))));
       }
       
       if(fullhadronic || semileptonic) d_cutLevel = 1.0;
       //       cout<< " test4 "<<endl;

       bool CR3     = fullhadronic_1lep && dataTriggerSL && nJets>3 && nCSVJets >1 && met > 200 && mt_full1<160 && minDPhi_6j > 1.; 
       bool CR3_nw  = fullhadronic_1lep && dataTriggerSL && nJets>3 && nCSVJets >1 && met > 200 && mt_full1<160; 
       bool CR3_tag = fullhadronic_1lep && dataTriggerSL && nJets>3 && met > 200 && minDPhi_6j >  1. && mt_full1<160 && nCSVLJets>0;
       bool CR4_tag = fullhadronic && dataTriggerHad && nJets>3 && met > 200  && abs(topFitRes.at(0).topPostFit.M()-172.5)>20 && nCSVLJets>0;
       bool CR5_tag = fullhadronic && dataTriggerHad && nJets>3 && met > 200 && minDPhi_6j >  1. && nCSVJets==0;
       bool CR6_tag = fullhadronic_1lep && dataTriggerSL && nJets>3 && met > 200 && minDPhi_6j >  1. && nCSVJets==0  && mt_full1<160.;     
       bool CR6_nw  = fullhadronic_1lep && dataTriggerSL && nJets>3 && met > 200 && nCSVJets==0  && mt_full1<160.;     
       bool CR7_tag = fullhadronic_2lep && dataTriggerSL && nJets>3 && metNoTwoLep > 200. && minDPhi_6j > 1. && nCSVJets==0 && mll> 60. && mll<120. && (lep1Flavour==lep2Flavour) && (lep1Charge != lep2Charge); //add charge requirement (30Jan16)
       bool CR7_nw = fullhadronic_2lep && dataTriggerSL && nJets>3 && metNoTwoLep > 200. && nCSVJets==0 && mll> 60. && mll<120. && (lep1Flavour==lep2Flavour) && (lep1Charge != lep2Charge); //add charge requirement (30Jan16)
       
       if(fullhadronic && dataTriggerHad && nJets>3 &&  nCSVJets>1 ){
	 
       	 bool SR_had = met > 200 && minDPhi_6j > 1.;
	 bool preS_had = met > 200;
	 
	 double residTopMassMin = 99999.;
	 // double residTopMassMinSecond = 99999.;//,residTopMassMinThird=99999.;
	 double chi_tt = 99999.;
	 double chi_min = 99999.;
	 int besttopidx = -1,  secondtopidx=-1;//, thirdtopidx=-1;
	 d_cutLevel =2.0;
	 //if (nResolvedHad==0)cout <<" nresolvedHad zero! "<<endl;
	 for (int t = 0; t < nResolvedHad;++t){
	   if (abs( topHMass[t]-172.5)<residTopMassMin){
	     residTopMassMin = abs(topHMass[t]-172.5);
	     besttopidx = t;	
	   }//end loop over resolved Had tops for the event j-th
	 }//end loop over resolved Had tops for the event i-th
	 //	 cout << " best top idx is "<< besttopidx << " res is "<< residTopMassMin<<endl;
	 
	 if(nJets > 5){
	   for (int t = 0; t < nResolvedHad;++t){
	     for (int tt = t; tt < nResolvedHad;++tt){
	       if(topHIndexB[t]==topHIndexB[tt])continue;
	       
	       else if(topHIndexJ1[t]==topHIndexJ1[tt] ||
		       topHIndexJ1[t]==topHIndexJ2[tt] ||
		       topHIndexJ2[t]==topHIndexJ1[tt] ||
		       topHIndexJ2[t]==topHIndexJ2[tt] 
		       ){continue;}
	       else{
		 chi_tt = abs( topHMass[t]-172.5) + abs( topHMass[tt]-172.5) + abs( topHWMass[t]-80.4) + abs( topHWMass[tt]-80.4);
		 if (chi_tt < chi_min){
		   chi_min =  chi_tt;
		   besttopidx = t;
		   secondtopidx = tt;
		   //cout << "best t1 " <<besttopidx<<endl;
		   //cout << "best t2 " << secondtopidx<<endl;
		 }	
		 //		  cout << "after " << residTopMassMinSecond<<endl;
	       }
	     } // 2nd for on top quarks
	   } // 1st for on top quarks
	   d_topHadChi = chi_min;
	 } // end if (nJets>5)
	 
	 //************** Filling Histograms *******************	
	 if( minDPhi_6j > 1. and nCSVJets>1)syst2B.fillHistogramsSysts(h_met,met,w,systWeights2B);
	 if(met > 160 and nCSVJets>1){
	   syst2B.fillHistogramsSysts(h_dphi,minDPhi,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_dphi_6j,minDPhi_6j,w,systWeights2B);
	   //syst2B.fillHistogramsSysts(h_dphi_10j,minDPhi_10j,w,systWeights2B);
	 }

	 //===================================================	  
	 // PRESELECTION ALL-HADRONIC
	 //	 std::cout << "preS_had: " << preS_had << " nCSVJets: " << nCSVJets << std::endl;
	 //       	 if(preS_had &&  nCSVJets>1 && !(isData=="DATA" && runNumber>runNumbBlinding) ) {
	 if(preS_had &&  nCSVJets>1  ) {
	   //	   std::cout << "met: " << met << " sync: " << sync << std::endl;
	   
	   if(met>200 && sync=="sync"){
	     /* Synchro exercise */
	     fileout_presel
	       <<std::fixed<<std::setprecision(0)
	       <<runNumber<<"   "
	       <<evtNumber<<"   "
	       <<lumiSec<<"   "
	       <<(nLooseMuons+nVetoElectrons)<<"   "
	       <<std::setprecision(3)
	       <<met<<"   "
	       <<std::setprecision(0)
	       <<nJets<<"   "
	       <<nCSVJets<<"   "
	       <<std::setprecision(3)
	       <<minDPhi_6j<<"   "
	       <<std::setprecision(0)
	       <<passMETFilters<<"    "
	       <<dataTriggerHad
	       <<std::endl;
	     }
	   
	   syst2B.fillHistogramsSysts(h_mva_first,topFitRes.at(0).mva,w,systWeights2B);
	   if(topFitRes.size()>1)syst2B.fillHistogramsSysts(h_mva_second,topFitRes.at(1).mva,w,systWeights2B);

	   //	   cout<< " test5 "<<endl;
	   syst2B.fillHistogramsSysts(h_topMassPostFit_preS,topFitRes.at(0).topPostFit.M(),w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_topMassPreFit_preS,topFitRes.at(0).topPreFit.M(),w,systWeights2B);
	   for ( const TopFitResults& res : topFitRes) {
	     syst2B.fillHistogramsSysts(h_bdt_qgid1,res.bdt_qgid1,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_bdt_qgid2,res.bdt_qgid2,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_bdt_dphij1b,res.bdt_dphij1b,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_bdt_dphij2b,res.bdt_dphij2b,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_bdt_drj1b,res.bdt_drj1b,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_bdt_drj1b,res.bdt_drj2b,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_bdt_bjcsv,res.bdt_bjcsv,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_bdt_jet1csv,res.bdt_jet1csv,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_bdt_jet2csv,res.bdt_jet2csv,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_bdt_prob,res.bdt_prob,w,systWeights2B);
	   }
	 }//end preselection without btagging
	 // BEGING SIGNAL SELECTION
	 //	 if(SR_had  && minDPhi_10j > 0.4 and !(isData=="DATA" && runNumber>runNumbBlinding) ){
	 //	 if(SR_had  &&  !(isData=="DATA" && runNumber>runNumbBlinding) ){
	 if(SR_had   ){
	   //if(SR_had  && minDPhi_10j > 0.4){
	   metFinal = met;
	   if(nJets> 3 and nCSVJets>1){
	     syst2B.fillHistogramsSysts(h_metFinal,met,w,systWeights2B);
	     if(nJets == 4) syst2B.fillHistogramsSysts(h_metFinal_4j,met,w,systWeights2B);
	     if(nJets == 6) syst2B.fillHistogramsSysts(h_metFinal_6j,met,w,systWeights2B);

	     if(met>200 && sync=="sync"){
	       //   /* Synchro exercise */
	       //std::cout<<runNumber<<" "<<evtNumber<<" "<<(nLooseMuons+nVetoElectrons)<<" "<<met<<" "<<nJets<<" "<<nCSVJets<<" "<<minDPhi_6j<<" "<<passMETFilters<<std::endl;
	       // fileout<<"RunNumber  EvtNumber Lumi  MET   NJets  nbjets w_pu k_fact w_top b_weight el_weight mu_weight met_trig "<<std::endl;
	       fileout
		 <<std::fixed<<std::setprecision(0)
		 <<runNumber<<"   "
		 <<evtNumber<<"   "
		 <<lumiSec<<"   "
		 <<std::setprecision(3)
		 <<met<<"   "
		 <<std::setprecision(0)
		 <<nJets<<"   "
		 <<nCSVJets<<"   "
		 <<std::setprecision(3)
		 <<w_pu<<"    "
		 <<k_fact<<"    "
		 <<bWeight2<<"    "
		 <<" ---   "
		 <<" ---   "
		 << ((METturnon13TeV(metPt[0],1,true))/(METturnon13TeV(metPt[0],1,false)))
		 <<std::endl;

	       //fileout
	       //<<std::fixed<<std::setprecision(0)
	       //<<runNumber<<"   "
	       //<<evtNumber<<"   "
	       //<<lumiSec<<"   "
	       //<<(nLooseMuons+nVetoElectrons)<<"   "
	       //<<std::setprecision(3)
	       //<<met<<"   "
	       //<<std::setprecision(0)
	       //<<nJets<<"   "
	       //<<nCSVJets<<"   "
	       //<<std::setprecision(3)
	       //<<minDPhi_6j<<"   "
	       //<<std::setprecision(0)
	       //<<passMETFilters<<"    "
	       //<<dataTriggerHad
	       //<<std::endl;

	     }
	   }
	 
	   // ***********************************************
	   // *********** VTAGGING CATEGORY *****************
	   // ***********************************************
	   float maxMVA = topFitRes.at(0).mva;
	   float maxMVA_2 = - 10.;
	   if(topFitRes.size()>1) maxMVA_2 = topFitRes.at(1).mva;
	   //float maxMVA_2 = topFitResFull.at(1).mva;
	   if(maxMVA>0.1 and maxMVA_2>-0.4 and minDPhi_10j>0.4){
	     //std::cout<<"Tagged Category"<<std::endl;
	     syst2B.fillHistogramsSysts(h_metFinal_tag,met,w,systWeights2B);	     
	     syst2B.fillHistogramsSysts(h_topMassPreFit,topFitRes.at(0).topPreFit.M(),w,systWeights2B);
	      syst2B.fillHistogramsSysts(h_topMassPostFit,topFitRes.at(0).topPostFit.M(),w,systWeights2B);
	    
	      double d_topHadTagCosTM = cos(topFitRes.at(0).TMPhiPostFit);//cos(deltaPhi(topFitRes.at(0).topPostFit.Phi(),metPhi[0]));
	      double d_topHadTagCosBM = cos(topFitRes.at(0).BMPhiPostFit);//deltaPhi(topFitRes.at(0).topPostFit.Phi(),metPhi[0]));
	      double d_topHadTagCosWM= cos(topFitRes.at(0).WMPhiPostFit);
	      double d_topHadTagCosWB= cos(topFitRes.at(0).WBPhiPostFit);
	      double d_topHadTagCosBB = cos(dphiBB);
	      
	      syst2B.fillHistogramsSysts(h_topHadTagCosTM,d_topHadTagCosTM,w,systWeights2B);
	      syst2B.fillHistogramsSysts(h_topHadTagCosBM,d_topHadTagCosBM,w,systWeights2B);
	      syst2B.fillHistogramsSysts(h_topHadTagCosWM,d_topHadTagCosWM,w,systWeights2B);
	      syst2B.fillHistogramsSysts(h_topHadTagCosWB,d_topHadTagCosWB,w,systWeights2B);
	      syst2B.fillHistogramsSysts(h_topHadTagCosBB,d_topHadTagCosBB,w,systWeights2B);
	      
	      if( d_topHadTagCosTM <0.8 
		  && d_topHadTagCosBM < 0.8
		  && d_topHadTagCosWM < 0.8
		  && d_topHadTagCosWB > -0.8
		  //		  && d_topHadTagCosBB >- 0.8
		  ){
		syst2B.fillHistogramsSysts(h_metFinal_Angular_tag,met,w,systWeights2B);
	      }
	   }//end of 2-tag cat
	   else if( (maxMVA>0.1 and maxMVA_2<-0.4) and  minDPhi_10j> 0.4 and nCSVJets>1){
	     //std::cout<<"1 tag: "<<maxMVA<<" "<<maxMVA_2<<std::endl;
	     syst2B.fillHistogramsSysts(h_metFinal_1tag,met,w,systWeights2B);	     
	     syst2B.fillHistogramsSysts(h_topMassPreFit_1tag,topFitRes.at(0).topPreFit.M(),w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_topMassPostFit_1tag,topFitRes.at(0).topPostFit.M(),w,systWeights2B);
	   }//end of 1-tag cat
	   else if((maxMVA<0.1 and maxMVA_2<-0.4) and minDPhi_10j> 0.6 and nCSVJets>1){
	     d_topCosBB= cos(dphiBB);
	     syst2B.fillHistogramsSysts(h_metFinal_untag,met,w,systWeights2B);	     
	     if(cos(topHBMPhi[besttopidx])<0.8 
		&& cos(topHTMPhi[besttopidx])<0.8 
		&& cos(topHWMPhi[besttopidx])<0.8 
		&& cos(topHWBPhi[besttopidx])>-0.8 
		&& d_topCosBB > -0.8
		){
	       syst2B.fillHistogramsSysts(h_metFinal_Angular_untag,met,w,systWeights2B);
	     }	
	   }//end of 0-tag cat
	   // vatgger selection 
	 
	   if( (nJets == 4|| nJets==5) and nCSVJets>1)syst2B.fillHistogramsSysts(h_metFinal_had45Jets,met,w,systWeights2B);
	   if( nJets>5 and nCSVJets>1)syst2B.fillHistogramsSysts(h_metFinal_had6Jets,met,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_nJets_fin,nJets,w,systWeights2B);       
	   syst2B.fillHistogramsSysts(h_nbJets_fin,nCSVJets,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_jet1Pt_fin,jetPt[0],w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_jet2Pt_fin,jetPt[1],w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_jet3Pt_fin,jetPt[2],w,systWeights2B);
	   if(nType1 ==0 && nType2 == 0  and nCSVJets>1) syst2B.fillHistogramsSysts(h_metFinal_noBoost,met,w,systWeights2B);
	 }// END SIGNAL REGION FOR HADRONIC CHANNEL
	 
	 // Plots for resolved hadronic top
	 //	 cout<< " test6 "<<endl;
	 if(besttopidx>=0 && nCSVJets>1){

	   d_cutLevel = 3.0;
	   d_topHadMass=topHMass[besttopidx]; 
	   d_topHadWMass=topHWMass[besttopidx]; 
	   d_topHadPt=topHPt[besttopidx];
	   d_topHadCosBM=cos(topHBMPhi[besttopidx]);
	   d_topHadCosWM=cos(topHWMPhi[besttopidx]);
	   d_topHadCosTM=cos(topHTMPhi[besttopidx]); 
	   d_topHadCosWB=cos(topHWBPhi[besttopidx]);
	   d_topCosBB=cos(dphiBB); 
	   d_topEtaBB=detaBB; 
	   d_topCosBB_MET=cos(dphiBB_MET);
	   if(nJets>5 && secondtopidx>=0){
	     d_topHadSecondTopMass=topHMass[secondtopidx];
	     d_topHadSecondTopWMass=topHWMass[secondtopidx];
	     d_topHadSecondTopPt=topHPt[secondtopidx];
	     d_topHadCosTT=cos(topHPhi[besttopidx]-topHPhi[secondtopidx]);
	     d_topHadEtaTT=abs(topHEta[besttopidx]-topHEta[secondtopidx]);
	     TLorentzVector topone, toptwo, ditop;	    
	     topone.SetPtEtaPhiE(topHPt[besttopidx],topHEta[besttopidx],topHPhi[besttopidx],topHE[besttopidx]);
	     toptwo.SetPtEtaPhiE(topHPt[secondtopidx],topHEta[secondtopidx],topHPhi[secondtopidx],topHE[secondtopidx]);
	     ditop = topone + toptwo;
	     d_topHadCosTT_MET=cos(deltaPhi(ditop.Phi(), metPhi[0]));
	   }//
	   else{;
	   }
	   //d_topHadSecondTopMass=-999.; d_topHadSecondTopPt=-999.; d_topHadCosTT=-999.;
	   //d_topSLMass=-999.; d_topSLPt  =-999.; d_topSLMT=-999.; d_topSLCosBM=-999.; d_topSLCosLM=-999.; d_topSLCosTM=-999.; d_topSLCosLBM=999.; d_topSLSecondTopMass=-999.; d_topSLSecondTopPt=-999.; d_topSLCosTT=-999.;
	   //     metFinal=-999.; d_minDPhi=-999.; d_dphi=-999.; d_mt=-999.; d_mt2w=-999.; d_njets=-999.; d_cutLevel=-999.;
	   d_mt =  mt; d_mt2w = mt2w;
	   if(metPt[0]>200 and nCSVJets>1){
	     d_cutLevel = 4.0;
	     // syst2B.fillHistogramsSysts(h_topHadCosBM,cos(topHBMPhi[besttopidx]),w,systWeights);
	     // syst2B.fillHistogramsSysts(h_topHadCosWM,cos(topHWMPhi[besttopidx]),w,systWeights);
	     // syst2B.fillHistogramsSysts(h_topHadCosTM,cos(topHTMPhi[besttopidx]),w,systWeights);
	     // syst2B.fillHistogramsSysts(h_topHadCosWB,cos(topHWBPhi[besttopidx]),w,systWeights);
	     syst2B.fillHistogramsSysts(h_minDPhi,minDPhi,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_minDPhi_6j,minDPhi_6j,w,systWeights2B);
	     //syst2B.fillHistogramsSysts(h_minDPhi_10j,minDPhi_10j,w,systWeights2B);
	     // syst2B.fillHistogramsSysts(h_topHadPt,topHPt[besttopidx],w,systWeights);
	     // syst2B.fillHistogramsSysts(h_topHadMass,topHMass[besttopidx],w,systWeights);
	     if(minDPhi_6j>1.2 and nCSVJets>1){
	       d_cutLevel = 5.0;
	       if(cos(topHBMPhi[besttopidx])<0.8 
		  && cos(topHTMPhi[besttopidx])<0.8 
		  && cos(topHWMPhi[besttopidx])<0.8 
		  && cos(topHWBPhi[besttopidx])>-0.8 
		  && d_topCosBB > -0.8
		  ){
		 syst2B.fillHistogramsSysts(h_metFinal_Angular,met,w,systWeights2B);
	      }
	     }
	     if(metPt[0]> 320 && minDPhi_6j>1. and nCSVJets>1){
	       d_cutLevel = 5.0;
	       // syst2B.fillHistogramsSysts(h_topHadPhiBM,topHBMPhi[besttopidx],w,systWeights);
	       // syst2B.fillHistogramsSysts(h_topHadPhiWM,topHWMPhi[besttopidx],w,systWeights);
	       // syst2B.fillHistogramsSysts(h_topHadPhiTM,topHTMPhi[besttopidx],w,systWeights);
	       // syst2B.fillHistogramsSysts(h_topHadPhiWB,topHWBPhi[besttopidx],w,systWeights);
	       syst2B.fillHistogramsSysts(h_topHadCosBM_allCuts,cos(topHBMPhi[besttopidx]),w,systWeights2B);
	       syst2B.fillHistogramsSysts(h_topHadCosWM_allCuts,cos(topHWMPhi[besttopidx]),w,systWeights2B);
	       syst2B.fillHistogramsSysts(h_topHadCosTM_allCuts,cos(topHTMPhi[besttopidx]),w,systWeights2B);
	       syst2B.fillHistogramsSysts(h_topHadCosWB_allCuts,cos(topHWBPhi[besttopidx]),w,systWeights2B);
	       syst2B.fillHistogramsSysts(h_minDPhi_allCuts,minDPhi,w,systWeights2B);
	       syst2B.fillHistogramsSysts(h_minDPhi_6j_allCuts,minDPhi_6j,w,systWeights2B);
	       syst2B.fillHistogramsSysts(h_topHadPt_allCuts,topHPt[besttopidx],w,systWeights2B);
	       syst2B.fillHistogramsSysts(h_topHadMass_allCuts,topHMass[besttopidx],w,systWeights2B);
	       syst2B.fillHistogramsSysts(h_topCosBB_allCuts,cos(dphiBB),w,systWeights2B);
	       syst2B.fillHistogramsSysts(h_topEtaBB_allCuts,detaBB,w,systWeights2B);
	       syst2B.fillHistogramsSysts(h_topCosBB_MET_allCuts,dphiBB_MET,w,systWeights2B);
	      

	      //	      fillSystHisto(h_topHadMass_allCutsMS,w,)
	      if(nJets>5 and nCSVJets>1){
		syst2B.fillHistogramsSysts(h_topHadCosTT_MET_allCuts,d_topHadCosTT_MET,w,systWeights2B);
		//syst2B.fillHistogramsSysts(h_topHadChi,d_topSLChi,w,systWeights2B);

	      }
	    }	    //end final selection
	  } // end preselection on met
	 } // end if bestCandidate found
	 
	 //for bkg estimation
	 //	 cout<< " test7 "<<endl;
	 if(abs(topHMass[besttopidx]-172.5)>20 and met>200){
	   syst2B.fillHistogramsSysts(h_metFinal_outtop,met,w,systWeights2B);

	  syst2B.fillHistogramsSysts(h_jet1Pt_CRqcd_had,jetPt[0],w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_jet2Pt_CRqcd_had,jetPt[1],w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_jet3Pt_CRqcd_had,jetPt[2],w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_nJets_CRqcd_had,nJets,w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_nbJets_CRqcd_had,nCSVJets,w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_dphi_CRqcd_had,minDPhi,w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_dphi_6j_CRqcd_had,minDPhi_6j,w,systWeights2B);
	  //syst2B.fillHistogramsSysts(h_dphi_10j_CRqcd_had,minDPhi_10j,w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_bjet1Pt_CRqcd_had,bjets[0].Pt(),w,systWeights2B);
	  if(nCSVJets>1) syst2B.fillHistogramsSysts(h_bjet2Pt_CRqcd_had,bjets[1].Pt(),w,systWeights2B);
	  if(nCSVJets>2) syst2B.fillHistogramsSysts(h_bjet3Pt_CRqcd_had,bjets[2].Pt(),w,systWeights2B);

	  if(met>200 && sync=="sync"){
	     /* Synchro exercise */
	    fileout_CRouttop
	      <<std::fixed<<std::setprecision(0)
	      <<runNumber<<"   "
	      <<evtNumber<<"   "
	      <<lumiSec<<"   "
	      <<(nLooseMuons+nVetoElectrons)<<"   "
	      <<std::setprecision(3)
	      <<met<<"   "
	      <<std::setprecision(0)
	      <<nJets<<"   "
	      <<nCSVJets<<"   "
	      <<std::setprecision(3)
	      <<minDPhi_6j<<"   "
	      <<std::setprecision(0)
	      <<passMETFilters<<"    "
	      <<dataTriggerHad
	      <<std::endl;
	  }
	}
              
	syst2B.fillHistogramsSysts(h_topHadMass_preMET,topHMass[besttopidx],w,systWeights2B);

	if(met>200  and nCSVJets>1){

	  syst2B.fillHistogramsSysts(h_jetQGL1_preS,jetQGL[1],w,systWeights2B); 
	  syst2B.fillHistogramsSysts(h_jetQGL2_preS,jetQGL[2],w,systWeights2B); 

	  syst2B.fillHistogramsSysts(h_jet1Pt_preS,jetPt[0],w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_jet2Pt_preS,jetPt[1],w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_jet3Pt_preS,jetPt[2],w,systWeights2B);

	  syst2B.fillHistogramsSysts(h_jet1Eta_preS,jetEta[0],w,systWeights2B);
          syst2B.fillHistogramsSysts(h_jet2Eta_preS,jetEta[1],w,systWeights2B);
          syst2B.fillHistogramsSysts(h_jet3Eta_preS,jetEta[2],w,systWeights2B);

	  syst2B.fillHistogramsSysts(h_jet1Phi_preS,jetPhi[0],w,systWeights2B);
          syst2B.fillHistogramsSysts(h_jet2Phi_preS,jetPhi[1],w,systWeights2B);
          syst2B.fillHistogramsSysts(h_jet3Phi_preS,jetPhi[2],w,systWeights2B);
	  
	  syst2B.fillHistogramsSysts(h_jet1CSV_preS,jetCSV[0],w,systWeights2B);
          syst2B.fillHistogramsSysts(h_jet2CSV_preS,jetCSV[1],w,systWeights2B);
          syst2B.fillHistogramsSysts(h_jet3CSV_preS,jetCSV[2],w,systWeights2B);

	  syst2B.fillHistogramsSysts(h_bjet1Pt_preS,bjets[0].Pt(),w,systWeights2B);
          if(nCSVJets>1) syst2B.fillHistogramsSysts(h_bjet2Pt_preS,bjets[1].Pt(),w,systWeights2B);
          if(nCSVJets>2) syst2B.fillHistogramsSysts(h_bjet3Pt_preS,bjets[2].Pt(),w,systWeights2B);

	  syst2B.fillHistogramsSysts(h_bjet1Eta_preS,bjets[0].Eta(),w,systWeights2B);
          if(nCSVJets>1) syst2B.fillHistogramsSysts(h_bjet2Eta_preS,bjets[1].Eta(),w,systWeights2B);
          if(nCSVJets>2) syst2B.fillHistogramsSysts(h_bjet3Eta_preS,bjets[2].Eta(),w,systWeights2B);

	  syst2B.fillHistogramsSysts(h_bjet1Phi_preS,bjets[0].Phi(),w,systWeights2B);
          if(nCSVJets>1) syst2B.fillHistogramsSysts(h_bjet2Phi_preS,bjets[1].Phi(),w,systWeights2B);
          if(nCSVJets>2) syst2B.fillHistogramsSysts(h_bjet3Phi_preS,bjets[2].Phi(),w,systWeights2B);

	  syst2B.fillHistogramsSysts(h_nPV_preS,nPV,1.0,systWeights2B);
	  syst2B.fillHistogramsSysts(h_nGoodPV_preS,nPV,1.0,systWeights2B);
	  syst2B.fillHistogramsSysts(h_nPV_w_preS,nPV,w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_nGoodPV_w_preS,nPV,w,systWeights2B);

	  //	  syst2B.fillHistogramsSysts(h_met_preS,met,w,systWeights2B,7,true, true);
	  syst2B.fillHistogramsSysts(h_met_preS,met,w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_dphi_preS,minDPhi,w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_dphi_6j_preS,minDPhi_6j,w,systWeights2B);
	  //syst2B.fillHistogramsSysts(h_dphi_10j_preS,minDPhi_10j,w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_nJets_preS,nJets,w,systWeights2B);
	  systZero.fillHistogramsSysts(h_nbJets_preS,nCSVJets,w,systWeightsNoSyst);//use 0-b weight
	
	  syst2B.fillHistogramsSysts(h_topHadCosBM_preS,cos(topHBMPhi[besttopidx]),w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_topHadCosWM_preS,cos(topHWMPhi[besttopidx]),w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_topHadCosTM_preS,cos(topHTMPhi[besttopidx]),w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_topHadCosWB_preS,cos(topHWBPhi[besttopidx]),w,systWeights2B);
	  	 
	  syst2B.fillHistogramsSysts(h_topHadPt_preS,topHPt[besttopidx],w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_topHadMass_preS,topHMass[besttopidx],w,systWeights2B);
	  
	  syst2B.fillHistogramsSysts(h_topHadPhiBM_preS,topHBMPhi[besttopidx],w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_topHadPhiWM_preS,topHWMPhi[besttopidx],w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_topHadPhiTM_preS,topHTMPhi[besttopidx],w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_topHadPhiWB_preS,topHWBPhi[besttopidx],w,systWeights2B);
	  
	  syst2B.fillHistogramsSysts(h_topCosBB_preS,cos(dphiBB),w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_topEtaBB_preS,detaBB,w,systWeights2B);
	  syst2B.fillHistogramsSysts(h_topCosBB_MET_preS,dphiBB_MET,w,systWeights2B);

	  syst2B.fillHistogramsSysts(h_topHadmassDrop_preS,topHmassDrop[besttopidx],w,systWeights2B);

	  if(nJets>5){
	    syst2B.fillHistogramsSysts(h_topHadCosTT_MET_preS,d_topHadCosTT_MET,w,systWeights2B);
	    syst2B.fillHistogramsSysts(h_topHadChi_preS,d_topSLChi,w,systWeights2B);
	  }
	}
       }//end of HADRONIC RESOLVED CAT
       
       //       cout<< " test8 "<<endl;
       if(CR3_tag){
	 double  residTopMassMin = 99999.;
	 int besttopidx = -1;
	 for (int t = 0; t < nResolvedHad;++t){
	     if (abs( topHMass[t]-172.5)<residTopMassMin){
	       residTopMassMin = abs(topHMass[t]-172.5);
	       besttopidx = t;	
	     }//end loop over resolved Had tops for the event j-th
	 }//end loop over resolved Had tops for the event i-th
	 syst12BL.fillHistogramsSysts(h_mva_first_CR3,topFitRes.at(0).mva,w,systWeights12BL);
	 if(topFitRes.size()>1)syst12BL.fillHistogramsSysts(h_mva_second_CR3,topFitRes.at(1).mva,w,systWeights12BL);
	 
	 for ( const TopFitResults& res : topFitRes) {
	   syst12BL.fillHistogramsSysts(h_bdt_qgid1_CR3,res.bdt_qgid1,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_qgid2_CR3,res.bdt_qgid2,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_dphij1b_CR3,res.bdt_dphij1b,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_dphij2b_CR3,res.bdt_dphij2b,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_drj1b_CR3,res.bdt_drj1b,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_drj2b_CR3,res.bdt_drj2b,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_bjcsv_CR3,res.bdt_bjcsv,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_jet1csv_CR3,res.bdt_jet1csv,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_jet2csv_CR3,res.bdt_jet2csv,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_prob_CR3,res.bdt_prob,w,systWeights12BL);
	 }   
	 
	 float maxMVA = topFitRes.at(0).mva;
	 float maxMVA_2 = -10.;
	 if(topFitRes.size()>1) maxMVA_2 = topFitRes.at(1).mva;
	 if( (maxMVA>0.1 and maxMVA_2>-0.4)  and minDPhi_10j>0.4 ){
	   syst12BL.fillHistogramsSysts(h_metFinal_tag_CR3,met,w,systWeights12BL);	     
	   double d_topHadTagCosTM = cos(topFitRes.at(0).TMPhiPostFit);//cos(deltaPhi(topFitRes.at(0).topPostFit.Phi(),metPhi[0]));
	   double d_topHadTagCosBM = cos(topFitRes.at(0).BMPhiPostFit);//deltaPhi(topFitRes.at(0).topPostFit.Phi(),metPhi[0]));
	   double d_topHadTagCosWM= cos(topFitRes.at(0).WMPhiPostFit);
	   double d_topHadTagCosWB= cos(topFitRes.at(0).WBPhiPostFit);
	   //	   double d_topHadTagCosBB = cos(dphiBB);
	   
	   syst12BL.fillHistogramsSysts(h_topHadTagCosTM_CR3,d_topHadTagCosTM,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_topHadTagCosBM_CR3,d_topHadTagCosBM,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_topHadTagCosWM_CR3,d_topHadTagCosWM,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_topHadTagCosWB_CR3,d_topHadTagCosWB,w,systWeights12BL);
	   //syst12BL.fillHistogramsSysts(h_topHadTagCosBB_CR3,d_topHadTagCosBB,w,systWeights12BL);
	   
	   //if(d_topHadTagCosTM)
	   if( d_topHadTagCosTM <0.8 
	       && d_topHadTagCosBM < 0.8
	       && d_topHadTagCosWM < 0.8
	       && d_topHadTagCosWB > -0.8
	       //	       && d_topHadTagCosBB >- 0.8
	       ){
	     //		cout << " tight fill "<<endl;
	     
	     syst12BL.fillHistogramsSysts(h_metFinal_Angular_tag_CR3,met,w,systWeights12BL);
	   }
	 }
	 else if(( maxMVA>0.1 and maxMVA_2<-0.4) and  minDPhi_10j> 0.4 and nCSVJets>1){
	   
	   syst2B.fillHistogramsSysts(h_metFinal_1tag_CR3,met,w,systWeights2B);	     
	   if(cos(topHBMPhi[besttopidx])<0.8 
	      && cos(topHTMPhi[besttopidx])<0.8 
	      && cos(topHWMPhi[besttopidx])<0.8 
	      && cos(topHWBPhi[besttopidx])>-0.8 
	      && d_topCosBB > -0.8
	      ){
	     syst2B.fillHistogramsSysts(h_metFinal_Angular_1tag_CR3,met,w,systWeights2B);	     
	   }
	 }
	 else if((maxMVA<0.1 and maxMVA_2<-0.4) and minDPhi_10j> 0.6 and nCSVJets>1){
	   syst2B.fillHistogramsSysts(h_metFinal_untag_CR3,met,w,systWeights2B);	     
	   if(cos(topHBMPhi[besttopidx])<0.8 
	      && cos(topHTMPhi[besttopidx])<0.8 
	      && cos(topHWMPhi[besttopidx])<0.8 
	      && cos(topHWBPhi[besttopidx])>-0.8 
	      && d_topCosBB > -0.8
	      ){syst2B.fillHistogramsSysts(h_metFinal_Angular_untag_CR3,met,w,systWeights2B);}
	 }

	 
       }//end of CR_3
       //       cout<< " test8.5 "<<endl;
              // bool CR4_had = fullhadronic && nJets>3 && nCSVJets >1 && met > 200 &&  abs(topFitRes.at(0).topPostFit.M()-172.5)>20.;
       //std::cout<<"CR4_tag: "<<CR4_tag<<std::endl;
       if(CR4_tag){
	 syst12BL.fillHistogramsSysts(h_mva_first_CR4,topFitRes.at(0).mva,w,systWeights12BL);
	 if(topFitRes.size()>1)syst12BL.fillHistogramsSysts(h_mva_second_CR4,topFitRes.at(1).mva,w,systWeights12BL);
	 for ( const TopFitResults& res : topFitRes) {
	   syst12BL.fillHistogramsSysts(h_bdt_qgid1_CR4,res.bdt_qgid1,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_qgid2_CR4,res.bdt_qgid2,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_dphij1b_CR4,res.bdt_dphij1b,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_dphij2b_CR4,res.bdt_dphij2b,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_drj1b_CR4,res.bdt_drj1b,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_drj2b_CR4,res.bdt_drj2b,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_bjcsv_CR4,res.bdt_bjcsv,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_jet1csv_CR4,res.bdt_jet1csv,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_jet2csv_CR4,res.bdt_jet2csv,w,systWeights12BL);
	   syst12BL.fillHistogramsSysts(h_bdt_prob_CR4,res.bdt_prob,w,systWeights12BL);
	 }   
	 float maxMVA = topFitRes.at(0).mva;
	 float maxMVA_2 = -10.;
	 if(topFitRes.size()>1) maxMVA_2 = topFitRes.at(1).mva;
	 
	 if((maxMVA>0.1 and maxMVA_2>-0.4)  and minDPhi_10j>0.4){
	     syst12BL.fillHistogramsSysts(h_metFinal_tag_CR4,met,w,systWeights12BL);	     
	   }
	   else if((maxMVA>0.1 and maxMVA_2<-0.4)  and  minDPhi_10j> 0.4 and nCSVJets>1){
	      syst2B.fillHistogramsSysts(h_metFinal_1tag_CR4,met,w,systWeights2B);	     
	    }
	   else if((maxMVA<0.1 and maxMVA_2<-0.4) and minDPhi_10j> 0.6 and nCSVJets>1){
	      syst2B.fillHistogramsSysts(h_metFinal_untag_CR4,met,w,systWeights2B);	     
	   }

       }//end of CR4_tag

       //       cout<< " test8.75 "<<endl;
       if(CR5_tag){
	 syst0B.fillHistogramsSysts(h_metFinal_CR5,met,w,systWeightsZeroB);	     
	 if(nJets == 4) syst0B.fillHistogramsSysts(h_metFinal_CR5_4j,met,w,systWeightsZeroB);     
	 if(nJets == 6) syst0B.fillHistogramsSysts(h_metFinal_CR5_6j,met,w,systWeightsZeroB);     

	 syst0B.fillHistogramsSysts(h_jet1Pt_CR5,jetPt[0],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_jet2Pt_CR5,jetPt[1],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_jet3Pt_CR5,jetPt[2],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_nJets_CR5,nJets,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_nbJets_CR5,nCSVJets,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_dphi_CR5,minDPhi,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_dphi_6j_CR5,minDPhi_6j,w,systWeightsZeroB);
	 //syst0B.fillHistogramsSysts(h_dphi_10j_CR5,minDPhi_10j,w,systWeightsZeroB);

	 //std::cout<<"In CR_category"<<std::endl;
	 syst0B.fillHistogramsSysts(h_mva_first_CR5,topFitRes.at(0).mva,w,systWeightsZeroB);
	 if(topFitRes.size()>1)syst0B.fillHistogramsSysts(h_mva_second_CR5,topFitRes.at(1).mva,w,systWeightsZeroB);
	 for ( const TopFitResults& res : topFitRes) {
	   syst0B.fillHistogramsSysts(h_bdt_qgid1_CR5,res.bdt_qgid1,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_qgid2_CR5,res.bdt_qgid2,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_dphij1b_CR5,res.bdt_dphij1b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_dphij2b_CR5,res.bdt_dphij2b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_drj1b_CR5,res.bdt_drj1b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_drj2b_CR5,res.bdt_drj2b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_bjcsv_CR5,res.bdt_bjcsv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_jet1csv_CR5,res.bdt_jet1csv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_jet2csv_CR5,res.bdt_jet2csv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_prob_CR5,res.bdt_prob,w,systWeightsZeroB);
	 }   
	 float maxMVA = topFitRes.at(0).mva;
	 float maxMVA_2 = -10.;
	 if(topFitRes.size()>1) maxMVA_2 = topFitRes.at(1).mva;
	 if((maxMVA>0.1 and maxMVA_2>-0.4)  and minDPhi_10j > 0.4){
	   syst0B.fillHistogramsSysts(h_metFinal_tag_CR5,met,w,systWeightsZeroB);	     
	 }
	 else if( (maxMVA>0.1 and maxMVA_2<-0.4) and  minDPhi_10j> 0.4){
	   syst0B.fillHistogramsSysts(h_metFinal_1tag_CR5,met,w,systWeightsZeroB);	     
	 }
	 else if((maxMVA<0.1 and maxMVA_2<-0.4) and minDPhi_10j> 0.6){
	   syst0B.fillHistogramsSysts(h_metFinal_untag_CR5,met,w,systWeightsZeroB);	     
	 }

	 //	 cout<< " test9 "<<endl;
	 if(met>200 && sync=="sync"){
	   /* Synchro exercise */
	   fileout_CR5
	     <<std::fixed<<std::setprecision(0)
	     <<runNumber<<"   "
	     <<evtNumber<<"   "
	     <<lumiSec<<"   "
	     <<(nLooseMuons+nVetoElectrons)<<"   "
	     <<std::setprecision(3)
	     <<met<<"   "
	     <<std::setprecision(0)
	     <<nJets<<"   "
	     <<nCSVJets<<"   "
	     <<std::setprecision(3)
	     <<minDPhi_6j<<"   "
	     <<bWeightZero<<"    "
	     <<bWeightZeroBTagUp<<"    "
	     <<bWeightZeroBTagDown<<"    "
	     <<bWeightZeroMisTagUp<<"    "
	     <<bWeightZeroMisTagDown<<"    "
	     <<std::setprecision(0)
	     //<<passMETFilters<<"    "
	     //<<dataTriggerHad
	     <<std::endl;
	 }
       }//end of CR5_tag

       if(CR6_tag){
	 syst0B.fillHistogramsSysts(h_metFinal_CR6,met,w,systWeightsZeroB);	     
	 
	 syst0B.fillHistogramsSysts(h_jet1Pt_CR6,jetPt[0],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_jet2Pt_CR6,jetPt[1],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_jet3Pt_CR6,jetPt[2],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_nJets_CR6,nJets,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_nbJets_CR6,nCSVJets,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_dphi_CR6,minDPhi,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_dphi_6j_CR6,minDPhi_6j,w,systWeightsZeroB);
	 //syst0B.fillHistogramsSysts(h_dphi_10j_CR6,minDPhi_10j,w,systWeightsZeroB);

	 if(lep1Flavour==13){
	   syst0B.fillHistogramsSysts(h_muonmetFinal_CR6,met,w,systWeightsZeroB);

	   syst0B.fillHistogramsSysts(h_mu1Pt_CR6,lep1Pt,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet1Pt_CR6,jetPt[0],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet2Pt_CR6,jetPt[1],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet3Pt_CR6,jetPt[2],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonnJets_CR6,nJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonnbJets_CR6,nCSVJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muondphi_CR6,minDPhi,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muondphi_6j_CR6,minDPhi_6j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_muondphi_10j_CR6,minDPhi_10j,w,systWeightsZeroB);
	 }

	 if(lep1Flavour==11){
           syst0B.fillHistogramsSysts(h_electronmetFinal_CR6,met,w,systWeightsZeroB);

           syst0B.fillHistogramsSysts(h_el1Pt_CR6,lep1Pt,w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electronjet1Pt_CR6,jetPt[0],w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electronjet2Pt_CR6,jetPt[1],w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electronjet3Pt_CR6,jetPt[2],w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electronnJets_CR6,nJets,w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electronnbJets_CR6,nCSVJets,w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electrondphi_CR6,minDPhi,w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electrondphi_6j_CR6,minDPhi_6j,w,systWeightsZeroB);
           //syst0B.fillHistogramsSysts(h_electrondphi_10j_CR6,minDPhi_10j,w,systWeightsZeroB);
	 }

	 syst0B.fillHistogramsSysts(h_mva_first_CR6,topFitRes.at(0).mva,w,systWeightsZeroB);
	  if(topFitRes.size()>1)syst0B.fillHistogramsSysts(h_mva_second_CR6,topFitRes.at(1).mva,w,systWeightsZeroB);
	 for ( const TopFitResults& res : topFitRes) {
	   syst0B.fillHistogramsSysts(h_bdt_qgid1_CR6,res.bdt_qgid1,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_qgid2_CR6,res.bdt_qgid2,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_dphij1b_CR6,res.bdt_dphij1b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_dphij2b_CR6,res.bdt_dphij2b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_drj1b_CR6,res.bdt_drj1b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_drj2b_CR6,res.bdt_drj2b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_bjcsv_CR6,res.bdt_bjcsv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_jet1csv_CR6,res.bdt_jet1csv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_jet2csv_CR6,res.bdt_jet2csv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_prob_CR6,res.bdt_prob,w,systWeightsZeroB);
	 }   
	   float maxMVA = topFitRes.at(0).mva;
	   float maxMVA_2 = -10.;
	   if(topFitRes.size()>1) maxMVA_2 = topFitRes.at(1).mva;
	 
	   if((maxMVA>0.1 and maxMVA_2>-0.4)  and minDPhi_10j>0.4){
	     syst0B.fillHistogramsSysts(h_metFinal_tag_CR6,met,w,systWeightsZeroB);	     
	   }
	   else if((maxMVA>0.1 and maxMVA_2<-0.4)  and  minDPhi_10j> 0.4 ){
	      syst0B.fillHistogramsSysts(h_metFinal_1tag_CR6,met,w,systWeightsZeroB);	     
	    }
	   else if((maxMVA<0.1 and maxMVA_2<-0.4) and minDPhi_10j> 0.6 ){
	      syst0B.fillHistogramsSysts(h_metFinal_untag_CR6,met,w,systWeightsZeroB);	     
	      
	   }

       }// end of CR6_tag

       if(CR6_nw){
	 syst0B.fillHistogramsSysts(h_metFinal_CR6nw,met,w,systWeightsZeroB);	     
	 if(nJets == 4) syst0B.fillHistogramsSysts(h_metFinal_CR6nw_4j,met,w,systWeightsZeroB);     
	 if(nJets == 6) syst0B.fillHistogramsSysts(h_metFinal_CR6nw_6j,met,w,systWeightsZeroB);     

	 syst0B.fillHistogramsSysts(h_jet1Pt_CR6nw,jetPt[0],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_jet2Pt_CR6nw,jetPt[1],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_jet3Pt_CR6nw,jetPt[2],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_nJets_CR6nw,nJets,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_nbJets_CR6nw,nCSVJets,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_dphi_CR6nw,minDPhi,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_dphi_6j_CR6nw,minDPhi_6j,w,systWeightsZeroB);
	 //syst0B.fillHistogramsSysts(h_dphi_10j_CR6nw,minDPhi_10j,w,systWeightsZeroB);

	 if(lep1Flavour==13){
	   syst0B.fillHistogramsSysts(h_muonmetFinal_CR6nw,met,w,systWeightsZeroB);

	   syst0B.fillHistogramsSysts(h_mu1Pt_CR6nw,lep1Pt,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet1Pt_CR6nw,jetPt[0],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet2Pt_CR6nw,jetPt[1],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet3Pt_CR6nw,jetPt[2],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonnJets_CR6nw,nJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonnbJets_CR6nw,nCSVJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muondphi_CR6nw,minDPhi,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muondphi_6j_CR6nw,minDPhi_6j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_muondphi_10j_CR6nw,minDPhi_10j,w,systWeightsZeroB);
	 }

	 if(lep1Flavour==11){
           syst0B.fillHistogramsSysts(h_electronmetFinal_CR6nw,met,w,systWeightsZeroB);

           syst0B.fillHistogramsSysts(h_el1Pt_CR6nw,lep1Pt,w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electronjet1Pt_CR6nw,jetPt[0],w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electronjet2Pt_CR6nw,jetPt[1],w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electronjet3Pt_CR6nw,jetPt[2],w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electronnJets_CR6nw,nJets,w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electronnbJets_CR6nw,nCSVJets,w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electrondphi_CR6nw,minDPhi,w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_electrondphi_6j_CR6nw,minDPhi_6j,w,systWeightsZeroB);
           //syst0B.fillHistogramsSysts(h_electrondphi_10j_CR6nw,minDPhi_10j,w,systWeightsZeroB);
	 }

	 if(met>200 && sync=="sync"){
	   /* Synchro exercise */
	   fileout_CR6
	     <<std::fixed<<std::setprecision(0)
	     <<runNumber<<"   "
	     <<evtNumber<<"   "
	     <<lumiSec<<"   "
	     <<(nLooseMuons+nVetoElectrons)<<"   "
	     <<std::setprecision(3)
	     <<met<<"   "
	     <<std::setprecision(0)
	     <<nJets<<"   "
	     <<nCSVJets<<"   "
	     <<std::setprecision(3)
	     <<minDPhi_6j<<"   "
	     <<std::setprecision(0)
	       <<passMETFilters<<"    "
	     <<dataTriggerHad
	     <<std::endl;
	 }
       }// end of CR6_nw
    
       if(CR7_tag){
	 syst0B.fillHistogramsSysts(h_metFinal_CR7,metNoTwoLep,w,systWeightsZeroB);	     

	 syst0B.fillHistogramsSysts(h_dphi_CR7,minDPhi,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_dphi_6j_CR7,minDPhi_6j,w,systWeightsZeroB);
	 //syst0B.fillHistogramsSysts(h_dphi_10j_CR7,minDPhi_10j,w,systWeightsZeroB);
	 //syst0B.fillHistogramsSysts(h_mt_CR7,mt,w,systWeightsZeroB);

	 syst0B.fillHistogramsSysts(h_nJets_CR7,nJets,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_nbJets_CR7,nCSVJets,w,systWeightsZeroB);

	 syst0B.fillHistogramsSysts(h_jet1Pt_CR7,jetPt[0],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_jet2Pt_CR7,jetPt[1],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_jet3Pt_CR7,jetPt[2],w,systWeightsZeroB);

         if(lep1Flavour==11 and lep2Flavour==11){
	   syst0B.fillHistogramsSysts(h_el1Pt_CR7,lep1Pt,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_el2Pt_CR7,lep2Pt,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronmetFinal_CR7,met,w,systWeightsZeroB);

	   syst0B.fillHistogramsSysts(h_electrondphi_CR7,minDPhi,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electrondphi_6j_CR7,minDPhi_6j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_electrondphi_10j_CR7,minDPhi_10j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_electronmt_CR7,mt,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_electronnJets_CR7,nJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronnbJets_CR7,nCSVJets,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_electronjet1Pt_CR7,jetPt[0],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronjet2Pt_CR7,jetPt[1],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronjet3Pt_CR7,jetPt[2],w,systWeightsZeroB);
	 }

         if(lep1Flavour==13 and lep2Flavour==13){
	   syst0B.fillHistogramsSysts(h_mu1Pt_CR7,lep1Pt,w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_mu2Pt_CR7,lep2Pt,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonmetFinal_CR7,met,w,systWeightsZeroB);

	   syst0B.fillHistogramsSysts(h_muondphi_CR7,minDPhi,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muondphi_6j_CR7,minDPhi_6j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_muondphi_10j_CR7,minDPhi_10j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_muonmt_CR7,mt,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_muonnJets_CR7,nJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonnbJets_CR7,nCSVJets,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_muonjet1Pt_CR7,jetPt[0],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet2Pt_CR7,jetPt[1],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet3Pt_CR7,jetPt[2],w,systWeightsZeroB);
	 }

	 if((lep1Flavour==13 and lep2Flavour==11) or (lep1Flavour==11 and lep2Flavour==13)) {
	   syst0B.fillHistogramsSysts(h_mixmetFinal_CR7,met,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_mixdphi_CR7,minDPhi,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_mixdphi_6j_CR7,minDPhi_6j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_mixdphi_10j_CR7,minDPhi_10j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_mixmt_CR7,mt,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_mixnJets_CR7,nJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_mixnbJets_CR7,nCSVJets,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_mixjet1Pt_CR7,jetPt[0],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_mixjet2Pt_CR7,jetPt[1],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_mixjet3Pt_CR7,jetPt[2],w,systWeightsZeroB);
	 }

	 syst0B.fillHistogramsSysts(h_mva_first_CR7,topFitRes.at(0).mva,w,systWeightsZeroB);
	 if(topFitRes.size()>1) syst0B.fillHistogramsSysts(h_mva_second_CR7,topFitRes.at(1).mva,w,systWeightsZeroB);
	 for ( const TopFitResults& res : topFitRes) {
	   syst0B.fillHistogramsSysts(h_bdt_qgid1_CR7,res.bdt_qgid1,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_qgid2_CR7,res.bdt_qgid2,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_dphij1b_CR7,res.bdt_dphij1b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_dphij2b_CR7,res.bdt_dphij2b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_drj1b_CR7,res.bdt_drj1b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_drj2b_CR7,res.bdt_drj2b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_bjcsv_CR7,res.bdt_bjcsv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_jet1csv_CR7,res.bdt_jet1csv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_jet2csv_CR7,res.bdt_jet2csv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_prob_CR7,res.bdt_prob,w,systWeightsZeroB);
	 }   
	   float maxMVA = topFitRes.at(0).mva;
	   float maxMVA_2 = -10.;
	   if(topFitRes.size()>1) maxMVA_2 = topFitRes.at(1).mva;

	   if((maxMVA>0.1 and maxMVA_2>-0.4)  and minDPhi_10j>0.4){
	     syst0B.fillHistogramsSysts(h_metFinal_tag_CR7,metNoTwoLep,w,systWeightsZeroB);	     
	   }
	   else if( (maxMVA>0.1 and maxMVA_2<-0.4) and  minDPhi_10j> 0.4 ){
	      syst0B.fillHistogramsSysts(h_metFinal_1tag_CR7,metNoTwoLep,w,systWeightsZeroB);	     
	    }
	   else if((maxMVA<0.1 and maxMVA_2<-0.4) and minDPhi_10j> 0.6 ){
	      syst0B.fillHistogramsSysts(h_metFinal_untag_CR7,metNoTwoLep,w,systWeightsZeroB);	     
	      
	   }

       }// end of CR7_tag


       if(CR7_nw){
	 syst0B.fillHistogramsSysts(h_metFinal_CR7nw,metNoTwoLep,w,systWeightsZeroB);	     
	 if(nJets == 4) syst0B.fillHistogramsSysts(h_metFinal_CR7nw_4j,metNoTwoLep,w,systWeightsZeroB);     
	 if(nJets == 6) syst0B.fillHistogramsSysts(h_metFinal_CR7nw_6j,metNoTwoLep,w,systWeightsZeroB);     

	 syst0B.fillHistogramsSysts(h_dphi_CR7nw,minDPhi,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_dphi_6j_CR7nw,minDPhi_6j,w,systWeightsZeroB);
	 //syst0B.fillHistogramsSysts(h_dphi_10j_CR7nw,minDPhi_10j,w,systWeightsZeroB);
	 //syst0B.fillHistogramsSysts(h_mt_CR7nw,mt,w,systWeightsZeroB);

	 syst0B.fillHistogramsSysts(h_nJets_CR7nw,nJets,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_nbJets_CR7nw,nCSVJets,w,systWeightsZeroB);

	 syst0B.fillHistogramsSysts(h_jet1Pt_CR7nw,jetPt[0],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_jet2Pt_CR7nw,jetPt[1],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_jet3Pt_CR7nw,jetPt[2],w,systWeightsZeroB);

         if(lep1Flavour==11 and lep2Flavour==11){
	   syst0B.fillHistogramsSysts(h_el1Pt_CR7nw,lep1Pt,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_el2Pt_CR7nw,lep2Pt,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronmetFinal_CR7nw,met,w,systWeightsZeroB);

	   syst0B.fillHistogramsSysts(h_electrondphi_CR7nw,minDPhi,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electrondphi_6j_CR7nw,minDPhi_6j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_electrondphi_10j_CR7nw,minDPhi_10j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_electronmt_CR7nw,mt,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_electronnJets_CR7nw,nJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronnbJets_CR7nw,nCSVJets,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_electronjet1Pt_CR7nw,jetPt[0],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronjet2Pt_CR7nw,jetPt[1],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronjet3Pt_CR7nw,jetPt[2],w,systWeightsZeroB);
	 }

         if(lep1Flavour==13 and lep2Flavour==13){
	   syst0B.fillHistogramsSysts(h_mu1Pt_CR7nw,lep1Pt,w,systWeightsZeroB);
           syst0B.fillHistogramsSysts(h_mu2Pt_CR7nw,lep2Pt,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonmetFinal_CR7nw,met,w,systWeightsZeroB);

	   syst0B.fillHistogramsSysts(h_muondphi_CR7nw,minDPhi,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muondphi_6j_CR7nw,minDPhi_6j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_muondphi_10j_CR7nw,minDPhi_10j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_muonmt_CR7nw,mt,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_muonnJets_CR7nw,nJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonnbJets_CR7nw,nCSVJets,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_muonjet1Pt_CR7nw,jetPt[0],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet2Pt_CR7nw,jetPt[1],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet3Pt_CR7nw,jetPt[2],w,systWeightsZeroB);
	 }

	 if((lep1Flavour==13 and lep2Flavour==11) or (lep1Flavour==11 and lep2Flavour==13)) {
	   syst0B.fillHistogramsSysts(h_mixmetFinal_CR7nw,met,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_mixdphi_CR7nw,minDPhi,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_mixdphi_6j_CR7nw,minDPhi_6j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_mixdphi_10j_CR7nw,minDPhi_10j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_mixmt_CR7nw,mt,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_mixnJets_CR7nw,nJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_mixnbJets_CR7nw,nCSVJets,w,systWeightsZeroB);
	   
	   syst0B.fillHistogramsSysts(h_mixjet1Pt_CR7nw,jetPt[0],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_mixjet2Pt_CR7nw,jetPt[1],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_mixjet3Pt_CR7nw,jetPt[2],w,systWeightsZeroB);
	 }

	 if(met>200 && sync=="sync"){
	   /* Synchro exercise */
	   fileout_CR7
	     <<std::fixed<<std::setprecision(0)
	     <<runNumber<<"   "
	     <<evtNumber<<"   "
	     <<lumiSec<<"   "
	     <<(nLooseMuons+nVetoElectrons)<<"   "
	     <<std::setprecision(3)
	     <<met<<"   "
	     <<std::setprecision(0)
	     <<nJets<<"   "
	     <<nCSVJets<<"   "
	     <<std::setprecision(3)
	     <<minDPhi_6j<<"   "
	     <<std::setprecision(0)
	     <<passMETFilters<<"    "
	     <<dataTriggerHad
	     <<std::endl;
	 }
       }// end of CR7nw_tag


       //for bkg estimation
       if(fullhadronic_1lep && nJets>3 ){
	 double  residTopMassMin_1lep = 99999.;
	 //double  residTopMassMin = 99999.;
         int topid_1lep = -1;

         for (int t = 0; t < nResolvedSemiLep;++t){
	   if (abs( topSLMass[t]-172.5)<residTopMassMin_1lep){
             residTopMassMin_1lep = abs(topSLMass[t]-172.5);
             topid_1lep  = t;
           }//end loop over resolved Had tops for the event j-th            
	 }//end loop over resolved Had tops for the event i-th              
         //cos(deltaPhi(ditop.Phi(), metPhi[0]))
       
	 // int besttopidx = -1;
	 // for (int t = 0; t < nResolvedHad;++t){
	 //   if (abs( topHMass[t]-172.5)<residTopMassMin){
	 //     residTopMassMin = abs(topHMass[t]-172.5);
	 //     besttopidx = t;	
	 //   }//end loop over resolved Had tops for the event j-th
	 // }//end loop over resolved Had tops for the event i-th

	 //float dlepmet = cos(fabs(deltaPhi(lep1.Phi(),metPhi[0])));

	 if(CR3){
	   //fullhadronic_1lep && nJets>3 + fullhadronic_1lep && dataTriggerSL && nJets>3 && nCSVJets >1 && met > 200 && minDPhi_6j > 1. && mt_full1<160
	   syst2B.fillHistogramsSysts(h_metFinal_SR_1lep,met,w,systWeights2B);
	 
	   syst2B.fillHistogramsSysts(h_lep1Pt_CR3,lep1Pt,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_lep1Eta_CR3,lep1Eta,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_lep1Phi_CR3,lep1Phi,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_lep1E_CR3,lep1E,w,systWeights2B);

	   syst2B.fillHistogramsSysts(h_jet1Pt_CRtt_had,jetPt[0],w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_jet2Pt_CRtt_had,jetPt[1],w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_jet3Pt_CRtt_had,jetPt[2],w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_nJets_CRtt_had,nJets,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_nbJets_CRtt_had,nCSVJets,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_lepPt_CRtt_had,lep1.Pt() ,w,systWeights2B);
	   
	   syst2B.fillHistogramsSysts(h_dphi_CRtt_had,minDPhi,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_dphi_6j_CRtt_had,minDPhi_6j,w,systWeights2B);
	   //syst2B.fillHistogramsSysts(h_dphi_10j_CRtt_had,minDPhi_10j,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_mt_CRtt_had,mt_full1,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_mt2w_CRtt_had,mt2w,w,systWeights2B);
	   
	   syst2B.fillHistogramsSysts(h_bjet1Pt_CRtt_had,bjets[0].Pt(),w,systWeights2B);
	   if(nCSVJets>1) syst2B.fillHistogramsSysts(h_bjet2Pt_CRtt_had,bjets[1].Pt(),w,systWeights2B);
	   if(nCSVJets>2) syst2B.fillHistogramsSysts(h_bjet3Pt_CRtt_had,bjets[2].Pt(),w,systWeights2B);

	   if(lep1Flavour==13){
	     syst2B.fillHistogramsSysts(h_mu1Pt_CR3,lep1Pt,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_mu1Eta_CR3,lep1Eta,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_mu1Phi_CR3,lep1Phi,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_mu1E_CR3,lep1E,w,systWeights2B);

	     syst2B.fillHistogramsSysts(h_muonjet1Pt_CRtt_had,jetPt[0],w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonjet2Pt_CRtt_had,jetPt[1],w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonjet3Pt_CRtt_had,jetPt[2],w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonnJets_CRtt_had,nJets,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonnbJets_CRtt_had,nCSVJets,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonlepPt_CRtt_had,lep1.Pt() ,w,systWeights2B);

	     syst2B.fillHistogramsSysts(h_muondphi_CRtt_had,minDPhi,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muondphi_6j_CRtt_had,minDPhi_6j,w,systWeights2B);
	     //syst2B.fillHistogramsSysts(h_muondphi_10j_CRtt_had,minDPhi_10j,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonmt_CRtt_had,mt_full1,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonmt2w_CRtt_had,mt2w,w,systWeights2B);
	   
	     syst2B.fillHistogramsSysts(h_muonbjet1Pt_CRtt_had,bjets[0].Pt(),w,systWeights2B);
	     if(nCSVJets>1) syst2B.fillHistogramsSysts(h_muonbjet2Pt_CRtt_had,bjets[1].Pt(),w,systWeights2B);
	     if(nCSVJets>2) syst2B.fillHistogramsSysts(h_muonbjet3Pt_CRtt_had,bjets[2].Pt(),w,systWeights2B);

	     syst2B.fillHistogramsSysts(h_muonmetFinal_SR_1lep,met,w,systWeights2B);	     
	   }

	   if(lep1Flavour==11){
	     syst2B.fillHistogramsSysts(h_el1Pt_CR3,lep1Pt,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_el1Eta_CR3,lep1Eta,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_el1Phi_CR3,lep1Phi,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_el1E_CR3,lep1E,w,systWeights2B);

	     syst2B.fillHistogramsSysts(h_electronjet1Pt_CRtt_had,jetPt[0],w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronjet2Pt_CRtt_had,jetPt[1],w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronjet3Pt_CRtt_had,jetPt[2],w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronnJets_CRtt_had,nJets,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronnbJets_CRtt_had,nCSVJets,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronlepPt_CRtt_had,lep1.Pt() ,w,systWeights2B);

	     syst2B.fillHistogramsSysts(h_electrondphi_CRtt_had,minDPhi,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electrondphi_6j_CRtt_had,minDPhi_6j,w,systWeights2B);
	     //syst2B.fillHistogramsSysts(h_electrondphi_10j_CRtt_had,minDPhi_10j,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronmt_CRtt_had,mt_full1,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronmt2w_CRtt_had,mt2w,w,systWeights2B);
	   
	     syst2B.fillHistogramsSysts(h_electronbjet1Pt_CRtt_had,bjets[0].Pt(),w,systWeights2B);
	     if(nCSVJets>1) syst2B.fillHistogramsSysts(h_electronbjet2Pt_CRtt_had,bjets[1].Pt(),w,systWeights2B);
	     if(nCSVJets>2) syst2B.fillHistogramsSysts(h_electronbjet3Pt_CRtt_had,bjets[2].Pt(),w,systWeights2B);

	     syst2B.fillHistogramsSysts(h_electronmetFinal_SR_1lep,met,w,systWeights2B);
	   }
	   
	   if(topid_1lep >=0) syst2B.fillHistogramsSysts(h_topSLMass_CRtt_had,topSLMass[topid_1lep],w,systWeights2B);
           else syst2B.fillHistogramsSysts(h_topSLMass_CRtt_had,-10.,w,systWeights2B);
	 
	 }
	 
	 if(CR3_nw){
	   //fullhadronic_1lep && nJets>3 + fullhadronic_1lep && dataTriggerSL && nJets>3 && nCSVJets >1 && met > 200 && minDPhi_6j > 1. && mt_full1<160
	   syst2B.fillHistogramsSysts(h_metFinal_CR3nw,met,w,systWeights2B);
	   
	   syst2B.fillHistogramsSysts(h_lep1Pt_CR3nw,lep1Pt,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_lep1Eta_CR3nw,lep1Eta,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_lep1Phi_CR3nw,lep1Phi,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_lep1E_CR3nw,lep1E,w,systWeights2B);
	   
	   syst2B.fillHistogramsSysts(h_jet1Pt_CR3nw,jetPt[0],w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_jet2Pt_CR3nw,jetPt[1],w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_jet3Pt_CR3nw,jetPt[2],w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_nJets_CR3nw,nJets,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_nbJets_CR3nw,nCSVJets,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_lepPt_CR3nw,lep1.Pt() ,w,systWeights2B);
	   
	   syst2B.fillHistogramsSysts(h_dphi_CR3nw,minDPhi,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_dphi_6j_CR3nw,minDPhi_6j,w,systWeights2B);
	   //syst2B.fillHistogramsSysts(h_dphi_10j_CR3nw,minDPhi_10j,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_mt_CR3nw,mt_full1,w,systWeights2B);
	   syst2B.fillHistogramsSysts(h_mt2w_CR3nw,mt2w,w,systWeights2B);
	   
	   syst2B.fillHistogramsSysts(h_bjet1Pt_CR3nw,bjets[0].Pt(),w,systWeights2B);
	   if(nCSVJets>1) syst2B.fillHistogramsSysts(h_bjet2Pt_CR3nw,bjets[1].Pt(),w,systWeights2B);
	   if(nCSVJets>2) syst2B.fillHistogramsSysts(h_bjet3Pt_CR3nw,bjets[2].Pt(),w,systWeights2B);
	   
	   if(lep1Flavour==13){
	     syst2B.fillHistogramsSysts(h_mu1Pt_CR3nw,lep1Pt,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_mu1Eta_CR3nw,lep1Eta,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_mu1Phi_CR3nw,lep1Phi,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_mu1E_CR3nw,lep1E,w,systWeights2B);
	     
	     syst2B.fillHistogramsSysts(h_muonjet1Pt_CR3nw,jetPt[0],w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonjet2Pt_CR3nw,jetPt[1],w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonjet3Pt_CR3nw,jetPt[2],w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonnJets_CR3nw,nJets,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonnbJets_CR3nw,nCSVJets,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonlepPt_CR3nw,lep1.Pt() ,w,systWeights2B);
	     
	     syst2B.fillHistogramsSysts(h_muondphi_CR3nw,minDPhi,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muondphi_6j_CR3nw,minDPhi_6j,w,systWeights2B);
	     //syst2B.fillHistogramsSysts(h_muondphi_10j_CR3nw,minDPhi_10j,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonmt_CR3nw,mt_full1,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_muonmt2w_CR3nw,mt2w,w,systWeights2B);
	     
	     syst2B.fillHistogramsSysts(h_muonbjet1Pt_CR3nw,bjets[0].Pt(),w,systWeights2B);
	     if(nCSVJets>1) syst2B.fillHistogramsSysts(h_muonbjet2Pt_CR3nw,bjets[1].Pt(),w,systWeights2B);
	     if(nCSVJets>2) syst2B.fillHistogramsSysts(h_muonbjet3Pt_CR3nw,bjets[2].Pt(),w,systWeights2B);
	     
	     syst2B.fillHistogramsSysts(h_muonmetFinal_CR3nw,met,w,systWeights2B);	     
	   }
	   
	   if(lep1Flavour==11){
	     syst2B.fillHistogramsSysts(h_el1Pt_CR3nw,lep1Pt,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_el1Eta_CR3nw,lep1Eta,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_el1Phi_CR3nw,lep1Phi,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_el1E_CR3nw,lep1E,w,systWeights2B);
	     
	     syst2B.fillHistogramsSysts(h_electronjet1Pt_CR3nw,jetPt[0],w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronjet2Pt_CR3nw,jetPt[1],w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronjet3Pt_CR3nw,jetPt[2],w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronnJets_CR3nw,nJets,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronnbJets_CR3nw,nCSVJets,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronlepPt_CR3nw,lep1.Pt() ,w,systWeights2B);
	     
	     syst2B.fillHistogramsSysts(h_electrondphi_CR3nw,minDPhi,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electrondphi_6j_CR3nw,minDPhi_6j,w,systWeights2B);
	     //syst2B.fillHistogramsSysts(h_electrondphi_10j_CR3nw,minDPhi_10j,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronmt_CR3nw,mt_full1,w,systWeights2B);
	     syst2B.fillHistogramsSysts(h_electronmt2w_CR3nw,mt2w,w,systWeights2B);
	     
	     syst2B.fillHistogramsSysts(h_electronbjet1Pt_CR3nw,bjets[0].Pt(),w,systWeights2B);
	     if(nCSVJets>1) syst2B.fillHistogramsSysts(h_electronbjet2Pt_CR3nw,bjets[1].Pt(),w,systWeights2B);
	     if(nCSVJets>2) syst2B.fillHistogramsSysts(h_electronbjet3Pt_CR3nw,bjets[2].Pt(),w,systWeights2B);
	     
	     syst2B.fillHistogramsSysts(h_electronmetFinal_CR3nw,met,w,systWeights2B);
	   }
	   
	   if(topid_1lep >=0) syst2B.fillHistogramsSysts(h_topSLMass_CR3nw,topSLMass[topid_1lep],w,systWeights2B);
	   else syst2B.fillHistogramsSysts(h_topSLMass_CR3nw,-10.,w,systWeights2B);
	   
	   
	   if(met>200 && sync=="sync"){
	     /* Synchro exercise */
	     fileout_CR3
	       <<std::fixed<<std::setprecision(0)
	       <<runNumber<<"   "
	       <<evtNumber<<"   "
	       <<lumiSec<<"   "
	       <<(nLooseMuons+nVetoElectrons)<<"   "
	       <<std::setprecision(3)
	       <<met<<"   "
	       <<std::setprecision(0)
	       <<nJets<<"   "
	       <<nCSVJets<<"   "
	       <<std::setprecision(3)
	       <<minDPhi_6j<<"   "
	       <<std::setprecision(0)
	       <<passMETFilters<<"    "
	       <<dataTriggerHad
	       <<std::endl;
	   }
	 }
       }

       // **********************************
       // Cut flow for fullhadronic analysis
       // **********************************
       
       if(fullhadronic){ 
	 n_binlep += w;
	 if (nJets>3){ 
	   n_binjet += w;
	   if(nCSVJets >1){ 
	     n_binbjet += w;
	     if(metPt[0]>200){ 
	       n_binmet160 += w;
	       if(met>320){ 
		 n_binmet320 += w;
		 if(minDPhi > 1.2){
		   n_phi += w;
		   if(nType1 ==0 && nType2 == 0)n_noBoost += w;
		 }
	       }
	     }
	   }
	 }
       }

       // ***************************************************
       // *************SEMILEPTONIC CHANNEL******************
       // ***************************************************
     
       bool CR1_lep = semileptonic_2lep && dataTriggerSL && nJets>2 && nCSVJets >0 && met>160;
       bool CR1_tag = semileptonic_2lep && dataTriggerSL && nJets>2 && met>160;
       bool CR2_lep = semileptonic  && dataTriggerSL && nJets>2 && nCSVJets== 0 && met> 160 &&  mt_full1>160.; 
       	 
       /* ==== Computing mt and mt2w =========*/
       if(((nTightMuons+nTightElectrons)>=1 && mt==0. && mt2w==0) ){
	 //	 if(((nTightMuons+nTightElectrons)>=1 && mt==0. && mt2w==0) or (evtNumber==25485) or (evtNumber==61147) or (evtNumber==57055) or  (evtNumber==74755)){
	 TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
	 float phi_lmet = fabs(deltaPhi(lep1.Phi(), metPhi[0]) );
	 mt = sqrt(2* lep1.Pt() * met* ( 1- cos(phi_lmet)));
	 Mt2Com_bisect *Mt2cal = new Mt2Com_bisect();
	 mt2w = Mt2cal->calculateMT2w(jets_nob,bjets,lep1, met_,"MT2w");
       }    

       if(semileptonic && nJets>2 && dataTriggerSL && nCSVJets >0 ){
	 bool SR_lep = (met > 160) && (mt > 160) &&  (mt2w > 200 and mt2w<9999) && minDPhi > 1.2  && nCSVJets >0;
	 bool SR_tag = (met > 160) && (mt > 160) &&  (mt2w > 200 and mt2w<9999) && minDPhi > 1.2;
	 bool preS_tag = met > 160;
	 
	 double residTopMassMin = 99999.,residTopMassMinHad=99999.;
	 int besttopidx = -1,besttophadidx=-1;
	 double chi_tt = 99999.;
	 double chi_min = 99999.;
	 
	 d_cutLevel = 2.0;
	 for (int t = 0; t < nResolvedSemiLep;++t){
	   if (abs( topSLMass[t]-172.5)<residTopMassMin){
	     residTopMassMin = (topSLMass[t]-172.5);
	     besttopidx = t;
	   }
	 }//end loop over resolved SL tops for the event i-th
	 
	 for (int t = 0; t < nResolvedHad;++t){
	   if (abs( topHMass[t]-172.5)<residTopMassMinHad){
	     residTopMassMinHad = abs(topHMass[t]-172.5);
	     besttophadidx = t;
	   }
	 }//end loop over resolved SL tops for the event i-th

	 if(nJets > 3 && nCSVJets >1){
	  for (int t = 0; t < nResolvedSemiLep;++t){
	   for (int tt = t; tt < nResolvedHad;++tt){

	     if(topSLIndexB[t]==topHIndexB[tt])continue;

	     else{
	       chi_tt = abs( topSLMass[t]-172.5) + abs( topHMass[tt]-172.5) + abs( topHWMass[tt]-80.4);
	       if (chi_tt < chi_min){
		 chi_min =  chi_tt;
		 besttopidx = t;
		 besttophadidx = tt;
	       }	
	       //		  cout << "after " << residTopMassMinSecond<<endl;
	     }
	   //} //
	   } // 2nd for on top quarks
	  } // 1st for on top quarks
	  d_topSLChi = chi_min;
	  //fillHistogramsSysts(h_topSLChi,chi_min,w,systWeights);
	 } // end if (nJets>5)

	 
	//************** Filling Histograms *******************	
	 d_mt= mt;
	 d_mt2w = mt2w;

	 if(met>160 && sync=="sync"){
	   //   /* Synchro exercise */
	   //std::cout<<runNumber<<" "<<evtNumber<<" "<<(nLooseMuons+nVetoElectrons)<<" "<<met<<" "<<nJets<<" "<<nCSVJets<<" "<<minDPhi_6j<<" "<<passMETFilters<<std::endl;
	   fileout_presel
	     <<std::fixed<<std::setprecision(0)
	     <<runNumber<<"   "
	     <<evtNumber<<"   "
	     <<lumiSec<<"   "
	     <<(nLooseMuons+nVetoElectrons)<<"   "
	     <<std::setprecision(3)
	     <<met<<"   "
	     <<std::setprecision(0)
	     <<nJets<<"   "
	     <<nCSVJets<<"   "
	     <<std::setprecision(3)
	     <<minDPhi<<"   "
	     <<std::setprecision(0)
	     <<passMETFilters<<"    "
	     <<dataTriggerSL
	     <<std::endl;
	 }
       
	 if((met > 160) && (mt > 160)  && minDPhi > 1.2  && nCSVJets >0){
	   syst12B.fillHistogramsSysts(h_mt2w,mt2w,w,systWeights12B);
	   if(lep1Flavour==13) syst12B.fillHistogramsSysts(h_muonmt2w,mt2w,w,systWeights12B);
	   if(lep1Flavour==11) syst12B.fillHistogramsSysts(h_electronmt2w,mt2w,w,systWeights12B);
	 }
       
	 if((met > 160) && (mt2w > 200 and mt2w<9999) && minDPhi > 1.2  && nCSVJets >0){
	   syst12B.fillHistogramsSysts(h_mt,mt,w,systWeights12B);
	   if(lep1Flavour==13) syst12B.fillHistogramsSysts(h_muonmt,mt,w,systWeights12B);
           if(lep1Flavour==11) syst12B.fillHistogramsSysts(h_electronmt,mt,w,systWeights12B);
	 }

	 if((mt > 160) &&  (mt2w > 200 and mt2w<9999) && minDPhi > 1.2  && nCSVJets >0){
	   syst12B.fillHistogramsSysts(h_met,met,w,systWeights12B);
	   if(lep1Flavour==13) syst12B.fillHistogramsSysts(h_muonmet,met,w,systWeights12B);
           if(lep1Flavour==11) syst12B.fillHistogramsSysts(h_electronmet,met,w,systWeights12B);
	 }

	 if((met > 160) &&  (mt2w > 200 and mt2w<9999) && (mt > 160)  && nCSVJets >0){
	   syst12B.fillHistogramsSysts(h_dphi,minDPhi,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_dphi_6j,minDPhi_6j,w,systWeights12B);
	   //syst12B.fillHistogramsSysts(h_dphi_10j,minDPhi_10j,w,systWeights12B);
	   if(lep1Flavour==13){
	     syst12B.fillHistogramsSysts(h_muondphi,minDPhi,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muondphi_6j,minDPhi_6j,w,systWeights12B);
	     //syst12B.fillHistogramsSysts(h_muondphi_10j,minDPhi_10j,w,systWeights12B);
	   }
	   if(lep1Flavour==11){
	     syst12B.fillHistogramsSysts(h_electrondphi,minDPhi,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_electrondphi_6j,minDPhi_6j,w,systWeights12B);
	     //syst12B.fillHistogramsSysts(h_electrondphi_10j,minDPhi_10j,w,systWeights12B);
	   }
	 }

	 //	 if((met > 160) && (mt > 160) && (mt2w<9999) && minDPhi > 1.2  && nCSVJets >0 && !(isData=="DATA" && runNumber>runNumbBlinding)){
	 if((met > 160) && (mt > 160) && (mt2w<9999) && minDPhi > 1.2  && nCSVJets >0 ){
	   syst12B.fillHistogramsSysts(h_metFinal_noMT2W,met,w,systWeights12B);
	 }

	 //if(SR_lep){
	 //	 if(SR_lep && !(isData=="DATA" && runNumber>runNumbBlinding)){
	 if(SR_lep ){
	   syst12B.fillHistogramsSysts(h_metFinal,met,w,systWeights12B);
	   cout << "filling " << endl;
	   if(nJets == 4) syst12B.fillHistogramsSysts(h_metFinal_4j,met,w,systWeights12B);
	   if(nJets == 6) syst12B.fillHistogramsSysts(h_metFinal_6j,met,w,systWeights12B);

	   if(lep1Flavour==13) syst12B.fillHistogramsSysts(h_muonmetFinal,met,w,systWeights12B);
           if(lep1Flavour==11) syst12B.fillHistogramsSysts(h_electronmetFinal,met,w,systWeights12B);
	   if(sync=="sync"){
	     /* Synchro exercise */
	     fileout
	       <<std::fixed<<std::setprecision(0)
	       <<runNumber<<"   "
	       <<evtNumber<<"   "
	       <<lumiSec<<"   "
	       <<std::setprecision(3)
	       <<met<<"   "
	       <<std::setprecision(0)
	       <<nJets<<"   "
	       <<nCSVJets<<"   "
	       <<std::setprecision(3)
	       <<w_pu<<"    "
	       <<k_fact<<"    "
	       <<bWeight12<<"    "
	       <<elTrigweight * elIDweight<<"    " 
	       <<muTrigweight * muIsoweight * muIDweight <<"    "
	       <<" ---   "
	       <<std::endl;

	     //fileout
	     //<<std::fixed<<std::setprecision(0)
	     ///<<runNumber<<"   "
	     ///<<evtNumber<<"   "
	     //<<lumiSec<<"   "
	     //<<(nLooseMuons+nVetoElectrons)<<"   "
	     //<<std::setprecision(3)
	     //<<met<<"   "
	     //<<std::setprecision(0)
	     //<<nJets<<"   "
	     //<<nCSVJets<<"   "
	     //<<std::setprecision(3)
	     //<<minDPhi<<"   "
	     //<<std::setprecision(0)
	     //<<passMETFilters<<"    "
	     //<<dataTriggerSL
	     //<<std::endl;
	   }
	  
	   metFinal = met;
	   if( nJets == 3  && nCSVJets >0)syst12B.fillHistogramsSysts(h_metFinal_sl3Jets,met,w,systWeights12B);
	   if( nJets>3  && nCSVJets >0)syst12B.fillHistogramsSysts(h_metFinal_sl4Jets,met,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_nJets_fin,nJets,w,systWeights12B);       
	   syst12B.fillHistogramsSysts(h_nbJets_fin,nCSVJets,w,systWeights12B);
	   if(nType1 ==0 && nType2 == 0 && nCSVJets >0) syst12B.fillHistogramsSysts(h_metFinal_noBoost,met,w,systWeights12B);
	 }// END FINAL SEMILEPTONIC BASELINE SELECTION	 

	 // ***********************************************
	 // *********** VTAGGING CATEGORY *****************
	 // ***********************************************
	 
	 if(preS_tag){
	   syst0B.fillHistogramsSysts(h_topMassPostFit_preS,topFitRes.at(0).topPostFit.M(),w,systWeightsZeroB); 
	   syst0B.fillHistogramsSysts(h_mva_first,topFitRes.at(0).mva,w,systWeightsZeroB);
	 }
	 //if(SR_tag){
	 //	 if(SR_tag && !(isData=="DATA" && runNumber>runNumbBlinding)){
	 if(SR_tag ){
	   float maxMVA = topFitRes.at(0).mva;
	   for ( const TopFitResults& res : topFitRes) {
	     syst0B.fillHistogramsSysts(h_bdt_qgid1,res.bdt_qgid1,w,systWeightsZeroB);
	     syst0B.fillHistogramsSysts(h_bdt_qgid2,res.bdt_qgid2,w,systWeightsZeroB);
	     syst0B.fillHistogramsSysts(h_bdt_dphij1b,res.bdt_dphij1b,w,systWeightsZeroB);
	     syst0B.fillHistogramsSysts(h_bdt_dphij2b,res.bdt_dphij2b,w,systWeightsZeroB);
	     syst0B.fillHistogramsSysts(h_bdt_drj1b,res.bdt_drj1b,w,systWeightsZeroB);
	     syst0B.fillHistogramsSysts(h_bdt_drj2b,res.bdt_drj2b,w,systWeightsZeroB);
	     syst0B.fillHistogramsSysts(h_bdt_bjcsv,res.bdt_bjcsv,w,systWeightsZeroB);
	     syst0B.fillHistogramsSysts(h_bdt_jet1csv,res.bdt_jet1csv,w,systWeightsZeroB);
	     syst0B.fillHistogramsSysts(h_bdt_jet2csv,res.bdt_jet2csv,w,systWeightsZeroB);
	     syst0B.fillHistogramsSysts(h_bdt_prob,res.bdt_prob,w,systWeightsZeroB);
	   }   
	   syst0B.fillHistogramsSysts(h_topMassPreFit,topFitRes.at(0).topPreFit.M(),w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_topMassPostFit,topFitRes.at(0).topPostFit.M(),w,systWeightsZeroB);
	   for(unsigned int i=0.;i<mva.size(); ++i){syst0B.fillHistogramsSysts(h_mva,mva[i],w,systWeightsZeroB);}
	   if(maxMVA>0.2  && nCSVLJets >0){
	     //td::cout<<"Tagged Category"<<std::endl;
	     syst12B.fillHistogramsSysts(h_metFinal_tag,met,w,systWeights12B);	     
	     double d_topHadTagCosTM = cos(topFitRes.at(0).TMPhiPostFit);//cos(deltaPhi(topFitRes.at(0).topPostFit.Phi(),metPhi[0]));
	     double d_topHadTagCosBM = cos(topFitRes.at(0).BMPhiPostFit);//deltaPhi(topFitRes.at(0).topPostFit.Phi(),metPhi[0]));
	     double d_topHadTagCosWM = cos(topFitRes.at(0).WMPhiPostFit);
	     double d_topHadTagCosWB = cos(topFitRes.at(0).WBPhiPostFit);
	     
	     syst12B.fillHistogramsSysts(h_topHadTagCosTM_lep,d_topHadTagCosTM,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_topHadTagCosBM_lep,d_topHadTagCosBM,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_topHadTagCosWM_lep,d_topHadTagCosWM,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_topHadTagCosWB_lep,d_topHadTagCosWB,w,systWeights12B);
	     
	     if ( d_topHadTagCosTM <0.8 && 
		  d_topHadTagCosBM < 0.8 && 
		  d_topHadTagCosWM < 0.8 && 
		  d_topHadTagCosWB > -0.8 ) { 
	       syst12B.fillHistogramsSysts(h_metFinal_Angular_tag,met,w,systWeights12B);        
	       //fillHistogramsSysts(h_metFinal_Angular,met,w,systWeights12B);        
	     }
	   }//end of 1-tag cat
	   else if( nCSVJets >0){
	     //std::cout<<"UnTagged Category"<<std::endl;
	     syst12B.fillHistogramsSysts(h_metFinal_untag,met,w,systWeights12B);	     
	     if (besttopidx>=0){
	       if(topSLMT[besttopidx]  > 175
		  && cos(topSLLMPhi[besttopidx]) < 0.8
		  && cos(topSLTMPhi[besttopidx]) < 0.8
		  && topSLMass[besttopidx] < 250.0) {
		 //		 cout <<" fuffa "<< " besttopidx "<< besttopidx << " topSLMT[besttopidx] "<< topSLMT[besttopidx] <<endl;
		 syst12B.fillHistogramsSysts(h_metFinal_Angular_untag,met,w,systWeights12B);        
		 //fillHistogramsSysts(h_metFinal_Angular,met,w,systWeights12B);        
	       }
	     }
	   }//end 0-tag cat
	 }// SR vatgger selection
	 
	 // Plots for resolved semileptonic top
	 if(besttopidx>=0 and nCSVJets >0){
	   d_cutLevel = 3.0;
	   d_topSLMass = topSLMass[besttopidx];
	   d_topSLPt = topSLPt[besttopidx];
	   d_topSLMT = topSLMT[besttopidx];
	   d_topSLCosBM = cos(topSLBMPhi[besttopidx]);
	   d_topSLCosLM = cos(topSLLMPhi[besttopidx]);
	   d_topSLCosTM = cos(topSLTMPhi[besttopidx]);
	   d_topSLCosLBM = cos(topSLLBMPhi[besttopidx]);
	   
	   if(besttophadidx>=0){
	     d_topHadMass = topHMass[besttophadidx];
	     d_topHadWMass = topHWMass[besttophadidx];
	     d_topHadPt = topHPt[besttophadidx];
	     d_topHadCosBM=cos(topHBMPhi[besttophadidx]);
	     d_topHadCosWM=cos(topHWMPhi[besttophadidx]);
	     d_topHadCosTM=cos(topHTMPhi[besttophadidx]); 
	     d_topHadCosWB=cos(topHWBPhi[besttophadidx]);

	     if(nJets>3 && nCSVJets >1 )
	       {
		 d_topHadCosTT=cos(topSLPhi[besttopidx]-topHPhi[besttophadidx]);
		 d_topHadEtaTT=abs(topSLEta[besttopidx]-topHEta[besttophadidx]);
		 
		 TLorentzVector topone, toptwo, ditop;	    
		 topone.SetPtEtaPhiE(topSLPt[besttopidx],topSLEta[besttopidx],topSLPhi[besttopidx],topSLE[besttopidx]);
		 toptwo.SetPtEtaPhiE(topHPt[besttophadidx],topHEta[besttophadidx],topHPhi[besttophadidx],topHE[besttophadidx]);
		 ditop = topone + toptwo;
		 d_topHadCosTT_MET=cos(deltaPhi(ditop.Phi(), metPhi[0]));
	       }
	     
	   }
	   
	   if(nCSVJets >1){
	     d_topCosBB=cos(dphiBB); 
	     d_topEtaBB=detaBB; 
	     d_topCosBB_MET=dphiBB_MET; 
	   }
	 	   
	   if(met>160 && topSLMT[besttopidx] > 0.  && minDPhi>1.2 && nCSVJets >0 ){
	     d_cutLevel = 4.0;
	     //fillHistogramsSysts(h_topSLMT,topSLMT[besttopidx],w,systWeights12B);
	     if(topSLMT[besttopidx]>172.5){
	       // fillHistogramsSysts(h_topSLCosBM,cos(topSLBMPhi[besttopidx]),w,systWeights12B);
	       // fillHistogramsSysts(h_topSLCosLM,cos(topSLLMPhi[besttopidx]),w,systWeights12B);
	       // fillHistogramsSysts(h_topSLCosTM,cos(topSLTMPhi[besttopidx]),w,systWeights12B);
	       // fillHistogramsSysts(h_topSLCosLBM,cos(topSLLBMPhi[besttopidx]),w,systWeights12B);
	       syst12B.fillHistogramsSysts(h_minDPhi,minDPhi,w,systWeights12B);
	       syst12B.fillHistogramsSysts(h_minDPhi_6j,minDPhi_6j,w,systWeights12B);
	       syst12B.fillHistogramsSysts(h_minDPhi_10j,minDPhi_10j,w,systWeights12B);
	       // fillHistogramsSysts(h_topSLPt,topSLPt[besttopidx],w,systWeights12B);
	       // fillHistogramsSysts(h_topSLMass,topSLMass[besttopidx],w,systWeights12B);
	       
	       // fillHistogramsSysts(h_topSLPhiBM,topSLBMPhi[besttopidx],w,systWeights12B);
	       // fillHistogramsSysts(h_topSLPhiLM,topSLLMPhi[besttopidx],w,systWeights12B);
	       // fillHistogramsSysts(h_topSLPhiTM,topSLTMPhi[besttopidx],w,systWeights12B);
	       // fillHistogramsSysts(h_topSLPhiLBM,topSLLBMPhi[besttopidx],w,systWeights12B);
	       if(mt > 160 && mt2w > 200 and mt2w<9999){
		 d_cutLevel = 5.0;
                 if(topSLMT[besttopidx]  > 175
                    && cos(topSLLMPhi[besttopidx]) < 0.8
                    && cos(topSLTMPhi[besttopidx]) < 0.8
                    && topSLMass[besttopidx] < 250.0
		    ){
		   syst12B.fillHistogramsSysts(h_metFinal_Angular,met, w,systWeights12B);
                 };
	       }

	       if(met> 160 && mt > 160 && mt2w > 200 && mt2w<9999){
		 d_cutLevel = 5.0;
		 syst12B.fillHistogramsSysts(h_topSLCosBM_allCuts,cos(topSLBMPhi[besttopidx]),w,systWeights12B);
		 syst12B.fillHistogramsSysts(h_topSLCosLM_allCuts,cos(topSLLMPhi[besttopidx]),w,systWeights12B);
		 syst12B.fillHistogramsSysts(h_topSLCosTM_allCuts,cos(topSLTMPhi[besttopidx]),w,systWeights12B);
		 syst12B.fillHistogramsSysts(h_topSLCosLBM_allCuts,cos(topSLLBMPhi[besttopidx]),w,systWeights12B);
		 syst12B.fillHistogramsSysts(h_minDPhi_allCuts,minDPhi,w,systWeights12B);
		 syst12B.fillHistogramsSysts(h_minDPhi_6j_allCuts,minDPhi_6j,w,systWeights12B);
		 syst12B.fillHistogramsSysts(h_topSLMT_allCuts,topSLMT[besttopidx],w,systWeights12B);
		 syst12B.fillHistogramsSysts(h_topSLPt_allCuts,topSLPt[besttopidx],w,systWeights12B);
		 syst12B.fillHistogramsSysts(h_topSLMass_allCuts,topSLMass[besttopidx],w,systWeights12B);
		 if(nCSVJets >1){
		   syst12B.fillHistogramsSysts(h_topCosBB_allCuts,cos(dphiBB),w,systWeights12B);
		   syst12B.fillHistogramsSysts(h_topEtaBB_allCuts,detaBB,w,systWeights12B);
		   syst12B.fillHistogramsSysts(h_topCosBB_MET_allCuts,dphiBB_MET,w,systWeights12B);
		 }
	       }//end final selection
	     }//end of cut on MT top
	   } // end preselection on met
	 } // end if bestCandidate found and btags > 0
	 
	 if(nCSVJets >0)syst12B.fillHistogramsSysts(h_topSLMass_preMET,topSLMass[besttopidx],w,systWeights12B);
	
	 if(met>160 and nCSVJets >0){
	   syst12B.fillHistogramsSysts(h_met_preS,met,w,systWeights12B);
           syst12B.fillHistogramsSysts(h_dphi_preS,minDPhi,w,systWeights12B);
           syst12B.fillHistogramsSysts(h_dphi_6j_preS,minDPhi_6j,w,systWeights12B);
           //syst12B.fillHistogramsSysts(h_dphi_10j_preS,minDPhi_10j,w,systWeights12B);
           syst12B.fillHistogramsSysts(h_mt_preS,mt,w,systWeights12B);

	   if(met>200) syst12B.fillHistogramsSysts(h_mt_hadpreS,mt,w,systWeights12B);

           syst12B.fillHistogramsSysts(h_mt2w_preS,mt2w,w,systWeights12B);
           syst12B.fillHistogramsSysts(h_nJets_preS,nJets,w,systWeights12B);
           systZero.fillHistogramsSysts(h_nbJets_preS,nCSVJets,w,systWeightsNoSyst);

	   syst12B.fillHistogramsSysts(h_topSLMT_preS,topSLMT[besttopidx],w,systWeights12B);

	   syst12B.fillHistogramsSysts(h_jetQGL1_preS,jetQGL[1],w,systWeights12B); 
	   syst12B.fillHistogramsSysts(h_jetQGL2_preS,jetQGL[2],w,systWeights12B); 

	   syst12B.fillHistogramsSysts(h_jet1Pt_preS,jetPt[0],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_jet2Pt_preS,jetPt[1],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_jet3Pt_preS,jetPt[2],w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_jet1Eta_preS,jetEta[0],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_jet2Eta_preS,jetEta[1],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_jet3Eta_preS,jetEta[2],w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_jet1Phi_preS,jetPhi[0],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_jet2Phi_preS,jetPhi[1],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_jet3Phi_preS,jetPhi[2],w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_jet1CSV_preS,jetCSV[0],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_jet2CSV_preS,jetCSV[1],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_jet3CSV_preS,jetCSV[2],w,systWeights12B);

	   syst12B.fillHistogramsSysts(h_bjet1Pt_preS,bjets[0].Pt(),w,systWeights12B);
	   if(nCSVJets>1) syst12B.fillHistogramsSysts(h_bjet2Pt_preS,bjets[1].Pt(),w,systWeights12B);
	   if(nCSVJets>2) syst12B.fillHistogramsSysts(h_bjet3Pt_preS,bjets[2].Pt(),w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_bjet1Eta_preS,bjets[0].Eta(),w,systWeights12B);
	   if(nCSVJets>1) syst12B.fillHistogramsSysts(h_bjet2Eta_preS,bjets[1].Eta(),w,systWeights12B);
	   if(nCSVJets>2) syst12B.fillHistogramsSysts(h_bjet3Eta_preS,bjets[2].Eta(),w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_bjet1Phi_preS,bjets[0].Phi(),w,systWeights12B);
	   if(nCSVJets>1) syst12B.fillHistogramsSysts(h_bjet2Phi_preS,bjets[1].Phi(),w,systWeights12B);
	   if(nCSVJets>2) syst12B.fillHistogramsSysts(h_bjet3Phi_preS,bjets[2].Phi(),w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_nPV_preS,nPV,1.0,systWeights12B);
	   syst12B.fillHistogramsSysts(h_nGoodPV_preS,nPV,1.0,systWeights12B);
	   syst12B.fillHistogramsSysts(h_nPV_w_preS,nPV,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_nGoodPV_w_preS,nPV,w,systWeights12B);

	   syst12B.fillHistogramsSysts(h_lep1Pt_preS,lep1Pt,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_lep1Eta_preS,lep1Eta,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_lep1Phi_preS,lep1Phi,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_lep1E_preS,lep1E,w,systWeights12B);

	   if(lep1Flavour==13){
	     syst12B.fillHistogramsSysts(h_mu1Pt_preS,lep1Pt,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_mu1Eta_preS,lep1Eta,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_mu1Phi_preS,lep1Phi,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_mu1E_preS,lep1E,w,systWeights12B);

	     syst12B.fillHistogramsSysts(h_muonmet_preS,met,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muondphi_preS,minDPhi,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muondphi_6j_preS,minDPhi_6j,w,systWeights12B);
	     //syst12B.fillHistogramsSysts(h_muondphi_10j_preS,minDPhi_10j,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muonmt_preS,mt,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muonmt2w_preS,mt2w,w,systWeights12B);

	     syst12B.fillHistogramsSysts(h_muonnJets_preS,nJets,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muonnbJets_preS,nCSVJets,w,systWeights12B);

	     syst12B.fillHistogramsSysts(h_muonjet1Pt_preS,jetPt[0],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muonjet2Pt_preS,jetPt[1],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muonjet3Pt_preS,jetPt[2],w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_muonjet1Eta_preS,jetEta[0],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muonjet2Eta_preS,jetEta[1],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muonjet3Eta_preS,jetEta[2],w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_muonjet1Phi_preS,jetPhi[0],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muonjet2Phi_preS,jetPhi[1],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muonjet3Phi_preS,jetPhi[2],w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_muonjet1CSV_preS,jetCSV[0],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muonjet2CSV_preS,jetCSV[1],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_muonjet3CSV_preS,jetCSV[2],w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_muonbjet1Pt_preS,bjets[0].Pt(),w,systWeights12B);
	     if(nCSVJets>1) syst12B.fillHistogramsSysts(h_muonbjet2Pt_preS,bjets[1].Pt(),w,systWeights12B);
	     if(nCSVJets>2) syst12B.fillHistogramsSysts(h_muonbjet3Pt_preS,bjets[2].Pt(),w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_muonbjet1Eta_preS,bjets[0].Eta(),w,systWeights12B);
	     if(nCSVJets>1) syst12B.fillHistogramsSysts(h_muonbjet2Eta_preS,bjets[1].Eta(),w,systWeights12B);
	     if(nCSVJets>2) syst12B.fillHistogramsSysts(h_muonbjet3Eta_preS,bjets[2].Eta(),w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_muonbjet1Phi_preS,bjets[0].Phi(),w,systWeights12B);
	     if(nCSVJets>1) syst12B.fillHistogramsSysts(h_muonbjet2Phi_preS,bjets[1].Phi(),w,systWeights12B);
	     if(nCSVJets>2) syst12B.fillHistogramsSysts(h_muonbjet3Phi_preS,bjets[2].Phi(),w,systWeights12B);
	   }

	   if(lep1Flavour==11){
	     syst12B.fillHistogramsSysts(h_el1Pt_preS,lep1Pt,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_el1Eta_preS,lep1Eta,w,systWeights12B);
             syst12B.fillHistogramsSysts(h_el1Phi_preS,lep1Phi,w,systWeights12B);
             syst12B.fillHistogramsSysts(h_el1E_preS,lep1E,w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_electronmet_preS,met,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_electrondphi_preS,minDPhi,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_electrondphi_6j_preS,minDPhi_6j,w,systWeights12B);
	     //syst12B.fillHistogramsSysts(h_electrondphi_10j_preS,minDPhi_10j,w,systWeights12B);

	     syst12B.fillHistogramsSysts(h_electronmt_preS,mt,w,systWeights12B);

	     if(abs(scEta[0])> 1.479 and abs(scEta[0])<2.5 ) {
	       syst12B.fillHistogramsSysts(h_mt_EC_preS,mt,w,systWeights12B);
	       if(met>200) syst12B.fillHistogramsSysts(h_mt_EC_hadpreS,mt,w,systWeights12B);
	     }
	     if(abs(scEta[0])<= 1.479) {
	       syst12B.fillHistogramsSysts(h_mt_BC_preS,mt,w,systWeights12B);
	       if(met>200) syst12B.fillHistogramsSysts(h_mt_BC_hadpreS,mt,w,systWeights12B);
	     }
	     syst12B.fillHistogramsSysts(h_electronmt2w_preS,mt2w,w,systWeights12B);

	     syst12B.fillHistogramsSysts(h_electronnJets_preS,nJets,w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_electronnbJets_preS,nCSVJets,w,systWeights12B);

	     syst12B.fillHistogramsSysts(h_electronjet1Pt_preS,jetPt[0],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_electronjet2Pt_preS,jetPt[1],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_electronjet3Pt_preS,jetPt[2],w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_electronjet1Eta_preS,jetEta[0],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_electronjet2Eta_preS,jetEta[1],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_electronjet3Eta_preS,jetEta[2],w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_electronjet1Phi_preS,jetPhi[0],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_electronjet2Phi_preS,jetPhi[1],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_electronjet3Phi_preS,jetPhi[2],w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_electronjet1CSV_preS,jetCSV[0],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_electronjet2CSV_preS,jetCSV[1],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_electronjet3CSV_preS,jetCSV[2],w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_electronbjet1Pt_preS,bjets[0].Pt(),w,systWeights12B);
	     if(nCSVJets>1) syst12B.fillHistogramsSysts(h_electronbjet2Pt_preS,bjets[1].Pt(),w,systWeights12B);
	     if(nCSVJets>2) syst12B.fillHistogramsSysts(h_electronbjet3Pt_preS,bjets[2].Pt(),w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_electronbjet1Eta_preS,bjets[0].Eta(),w,systWeights12B);
	     if(nCSVJets>1) syst12B.fillHistogramsSysts(h_electronbjet2Eta_preS,bjets[1].Eta(),w,systWeights12B);
	     if(nCSVJets>2) syst12B.fillHistogramsSysts(h_electronbjet3Eta_preS,bjets[2].Eta(),w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_electronbjet1Phi_preS,bjets[0].Phi(),w,systWeights12B);
	     if(nCSVJets>1) syst12B.fillHistogramsSysts(h_electronbjet2Phi_preS,bjets[1].Phi(),w,systWeights12B);
	     if(nCSVJets>2) syst12B.fillHistogramsSysts(h_electronbjet3Phi_preS,bjets[2].Phi(),w,systWeights12B);
	   }

	   if(topSLMT[besttopidx]>172.5){
	     syst12B.fillHistogramsSysts(h_topSLCosBM_preS,cos(topSLBMPhi[besttopidx]),w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_topSLCosLM_preS,cos(topSLLMPhi[besttopidx]),w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_topSLCosTM_preS,cos(topSLTMPhi[besttopidx]),w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_topSLCosLBM_preS,cos(topSLLBMPhi[besttopidx]),w,systWeights12B);
	     
	     syst12B.fillHistogramsSysts(h_topSLPt_preS,topSLPt[besttopidx],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_topSLMass_preS,topSLMass[besttopidx],w,systWeights12B);
	     systZero.fillHistogramsSysts(h_topSLMass_preS_NS,topSLMass[besttopidx],w,systWeightsNoSyst);
	     
	     syst12B.fillHistogramsSysts(h_topSLPhiBM_preS,topSLBMPhi[besttopidx],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_topSLPhiLM_preS,topSLLMPhi[besttopidx],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_topSLPhiTM_preS,topSLTMPhi[besttopidx],w,systWeights12B);
	     syst12B.fillHistogramsSysts(h_topSLPhiLBM_preS,topSLLBMPhi[besttopidx],w,systWeights12B);
	     
	     if(nCSVJets >1){
	       syst2B.fillHistogramsSysts(h_topCosBB_preS,cos(dphiBB),w,systWeights2B);
	       syst2B.fillHistogramsSysts(h_topEtaBB_preS,detaBB,w,systWeights2B);
	       syst2B.fillHistogramsSysts(h_topCosBB_MET_preS,dphiBB_MET,w,systWeights2B);
	     }
	   }
	 }//end of plotting at preSelection 
	 
       }//end of SEMILEPTONIC RESOLVED CAT 
       
       // ==============Starting Bkgs==============================

       //for bkg estimation
       if(CR1_tag){
	 syst0B.fillHistogramsSysts(h_mva_first_CR1,topFitRes.at(0).mva,w,systWeightsZeroB);
	 for ( const TopFitResults& res : topFitRes) {
	   //	   std::cout<<res.bdt_qgid2<<std::endl;
	   syst0B.fillHistogramsSysts(h_bdt_qgid1_CR1,res.bdt_qgid1,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_qgid2_CR1,res.bdt_qgid2,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_dphij1b_CR1,res.bdt_dphij1b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_dphij2b_CR1,res.bdt_dphij2b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_drj1b_CR1,res.bdt_drj1b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_drj2b_CR1,res.bdt_drj2b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_bjcsv_CR1,res.bdt_bjcsv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_jet1csv_CR1,res.bdt_jet1csv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_jet2csv_CR1,res.bdt_jet2csv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_prob_CR1,res.bdt_prob,w,systWeightsZeroB);
	 }   
	 float maxMVA = topFitRes.at(0).mva;
	
	 if(maxMVA>0.2 &&  nCSVLJets >0 ){
	     syst12B.fillHistogramsSysts(h_metFinal_tag_CR1,metNoLep,w,systWeights12B);
	 
	 }// end of 1-tag cat
	 else if(maxMVA<0.2 &&  nCSVJets >0 ){
	   syst12B.fillHistogramsSysts(h_metFinal_untag_CR1,metNoLep,w,systWeights12B);
	 }//end of 0-tag cat
       }// end of CR1_tag


       if(CR1_lep){ 

	 syst12B.fillHistogramsSysts(h_metFinal_2lep,met,w,systWeights12B);

	 if(nJets == 4) syst12B.fillHistogramsSysts(h_metFinal_2lep_4j,met,w,systWeights12B);
	 if(nJets == 6) syst12B.fillHistogramsSysts(h_metFinal_2lep_6j,met,w,systWeights12B);
	   
	 syst12B.fillHistogramsSysts(h_lep1Pt_CR1,lep1Pt,w,systWeights12B);
	 syst12B.fillHistogramsSysts(h_lep1Eta_CR1,lep1Eta,w,systWeights12B);
	 syst12B.fillHistogramsSysts(h_lep1Phi_CR1,lep1Phi,w,systWeights12B);
	 syst12B.fillHistogramsSysts(h_lep1E_CR1,lep1E,w,systWeights12B);
	 
	 syst12B.fillHistogramsSysts(h_lep2Pt_CR1,lep2Pt,w,systWeights12B);
         syst12B.fillHistogramsSysts(h_lep2Eta_CR1,lep2Eta,w,systWeights12B);
         syst12B.fillHistogramsSysts(h_lep2Phi_CR1,lep2Phi,w,systWeights12B);
         syst12B.fillHistogramsSysts(h_lep2E_CR1,lep2E,w,systWeights12B);
         
	 syst12B.fillHistogramsSysts(h_jet1Pt_CRtt2_lep,jetPt[0],w,systWeights12B);
	 syst12B.fillHistogramsSysts(h_jet2Pt_CRtt2_lep,jetPt[1],w,systWeights12B);
	 syst12B.fillHistogramsSysts(h_jet3Pt_CRtt2_lep,jetPt[2],w,systWeights12B);
	 syst12B.fillHistogramsSysts(h_nJets_CRtt2_lep,nJets,w,systWeights12B);
	 syst12B.fillHistogramsSysts(h_nbJets_CRtt2_lep,nCSVJets,w,systWeights12B);
	 syst12B.fillHistogramsSysts(h_mt_CR1,mt,w,systWeights12B);
	 syst12B.fillHistogramsSysts(h_mt2w_CR1,mt2w,w,systWeights12B);
	 
	 //if(topid_2lep>=0) fillHistogramsSysts(h_topSLMass_CRtt2_lep,topSLMass[topid_2lep],w,systWeights12B);
	 //else fillHistogramsSysts(h_topSLMass_CRtt2_lep,-10,w,systWeights12B);

         if(lep1Flavour==11 and lep2Flavour==11){
	   syst12B.fillHistogramsSysts(h_el1Pt_CR1,lep1Pt,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_el2Pt_CR1,lep2Pt,w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_electronmetFinal_2lep,met,w,systWeights12B);

	   syst12B.fillHistogramsSysts(h_electrondphi_CR1,minDPhi,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_electrondphi_6j_CR1,minDPhi_6j,w,systWeights12B);
	   //syst12B.fillHistogramsSysts(h_electrondphi_10j_CR1,minDPhi_10j,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_electronmt_CR1,mt,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_electronmt2w_CR1,mt2w,w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_electronnJets_CR1,nJets,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_electronnbJets_CR1,nCSVJets,w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_electronjet1Pt_CR1,jetPt[0],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_electronjet2Pt_CR1,jetPt[1],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_electronjet3Pt_CR1,jetPt[2],w,systWeights12B);

	   syst12B.fillHistogramsSysts(h_electronbjet1Pt_CR1,bjets[0].Pt(),w,systWeights12B);
	   if(nCSVJets>1) syst12B.fillHistogramsSysts(h_electronbjet2Pt_CR1,bjets[1].Pt(),w,systWeights12B);
	   if(nCSVJets>2) syst12B.fillHistogramsSysts(h_electronbjet3Pt_CR1,bjets[2].Pt(),w,systWeights12B);
	 }

         if(lep1Flavour==13 and lep2Flavour==13){
	   syst12B.fillHistogramsSysts(h_mu1Pt_CR1,lep1Pt,w,systWeights12B);
           syst12B.fillHistogramsSysts(h_mu2Pt_CR1,lep2Pt,w,systWeights12B);

	   syst12B.fillHistogramsSysts(h_muonmetFinal_2lep,met,w,systWeights12B);

	   syst12B.fillHistogramsSysts(h_muondphi_CR1,minDPhi,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_muondphi_6j_CR1,minDPhi_6j,w,systWeights12B);
	   //syst12B.fillHistogramsSysts(h_muondphi_10j_CR1,minDPhi_10j,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_muonmt_CR1,mt,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_muonmt2w_CR1,mt2w,w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_muonnJets_CR1,nJets,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_muonnbJets_CR1,nCSVJets,w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_muonjet1Pt_CR1,jetPt[0],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_muonjet2Pt_CR1,jetPt[1],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_muonjet3Pt_CR1,jetPt[2],w,systWeights12B);

	   syst12B.fillHistogramsSysts(h_muonbjet1Pt_CR1,bjets[0].Pt(),w,systWeights12B);
	   if(nCSVJets>1) syst12B.fillHistogramsSysts(h_muonbjet2Pt_CR1,bjets[1].Pt(),w,systWeights12B);
	   if(nCSVJets>2) syst12B.fillHistogramsSysts(h_muonbjet3Pt_CR1,bjets[2].Pt(),w,systWeights12B);
	 }

	 if((lep1Flavour==13 and lep2Flavour==11) or (lep1Flavour==11 and lep2Flavour==13)) {
	   syst12B.fillHistogramsSysts(h_mixmetFinal_2lep,met,w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_mixdphi_CR1,minDPhi,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_mixdphi_6j_CR1,minDPhi_6j,w,systWeights12B);
	   //syst12B.fillHistogramsSysts(h_mixdphi_10j_CR1,minDPhi_10j,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_mixmt_CR1,mt,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_mixmt2w_CR1,mt2w,w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_mixnJets_CR1,nJets,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_mixnbJets_CR1,nCSVJets,w,systWeights12B);
	   
	   syst12B.fillHistogramsSysts(h_mixjet1Pt_CR1,jetPt[0],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_mixjet2Pt_CR1,jetPt[1],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_mixjet3Pt_CR1,jetPt[2],w,systWeights12B);

	   syst12B.fillHistogramsSysts(h_mixbjet1Pt_CR1,bjets[0].Pt(),w,systWeights12B);
	   if(nCSVJets>1) syst12B.fillHistogramsSysts(h_mixbjet2Pt_CR1,bjets[1].Pt(),w,systWeights12B);
	   if(nCSVJets>2) syst12B.fillHistogramsSysts(h_mixbjet3Pt_CR1,bjets[2].Pt(),w,systWeights12B);
	 }
      
	 if( (abs(mll-90)<20)){
	   syst12B.fillHistogramsSysts(h_metFinal_2lep_Z_nobtag,met,w,systWeights12B);
	   
	   if(lep1Flavour==13 and lep2Flavour==13) syst12B.fillHistogramsSysts(h_muonmetFinal_2lep_Z_nobtag,met,w,systWeights12B);
	   if(lep1Flavour==11 and lep2Flavour==11) syst12B.fillHistogramsSysts(h_electronmetFinal_2lep_Z_nobtag,met,w,systWeights12B);
	   if((lep1Flavour==13 and lep2Flavour==11) or (lep1Flavour==11 and lep2Flavour==13)) syst12B.fillHistogramsSysts(h_mixmetFinal_2lep_Z_nobtag,met,w,systWeights12B);

	   syst12B.fillHistogramsSysts(h_jet1Pt_CRZ_lep,jetPt[0],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_jet2Pt_CRZ_lep,jetPt[1],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_jet3Pt_CRZ_lep,jetPt[2],w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_nJets_CRZ_lep,nJets,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_nbJets_CRZ_lep,nCSVJets,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_dphi_CRZ_lep,minDPhi,w,systWeights12B);
	   //syst12B.fillHistogramsSysts(h_mt_CRZ_lep,mt,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_mt2w_CRZ_lep,mt2w,w,systWeights12B);
	   syst12B.fillHistogramsSysts(h_lepPt_CRZ_lep,lep1.Pt(),w,systWeights12B);
	   //if(topid_2lep >=0) syst12B.fillHistogramsSysts(h_topSLMass_CRZ_lep,topHMass[topid_2lep ],w,systWeights12B);
           //else syst12B.fillHistogramsSysts(h_topSLMass_CRZ_lep,-10.,w,systWeights12B);
	 }
	 // syst12B.fillHistogramsSysts(h_metFinal_2lep_Full,met,w,systWeights12B);
	 // syst12B.fillHistogramsSysts(h_jet1Pt_2lep_Full_lep,jetPt[0],w,systWeights12B);
	 // syst12B.fillHistogramsSysts(h_jet2Pt_2lep_Full_lep,jetPt[1],w,systWeights12B);
	 // syst12B.fillHistogramsSysts(h_jet3Pt_2lep_Full_lep,jetPt[2],w,systWeights12B);
	 // syst12B.fillHistogramsSysts(h_nJets_2lep_Full_lep,nJets,w,systWeights12B);
	 // syst12B.fillHistogramsSysts(h_nbJets_2lep_Full_lep,nCSVJets,w,systWeights12BNoSyst);
	 // syst12B.fillHistogramsSysts(h_dphi_2lep_Full_lep,minDPhi,w,systWeights12B);
	 // syst12B.fillHistogramsSysts(h_mt_CRttFull_lep,mt,w,systWeights12B);
	 // //syst12B.fillHistogramsSysts(h_mt_CRZ_lep,mt,w,systWeights12B);
	 // syst12B.fillHistogramsSysts(h_mt2w_2lep_Full_lep,mt2w,w,systWeights12B);
	 // syst12B.fillHistogramsSysts(h_lepPt_2lep_Full_lep,lept_ord[0].Pt(),w,systWeights12B);
	 //if(topid_2lep >=0) syst12B.fillHistogramsSysts(h_topSLMass_CRZ_lep,topHMass[topid_2lep ],w,systWeights12B);
	 //else syst12B.fillHistogramsSysts(h_topSLMass_CRZ_lep,-10.,w,systWeights12B);
       	 
	 if(met>160 && sync=="sync"){
	   /* Synchro exercise */
	   fileout_CR1lep
	     <<std::fixed<<std::setprecision(0)
	     <<runNumber<<"   "
	     <<evtNumber<<"   "
	     <<lumiSec<<"   "
	     <<(nLooseMuons+nVetoElectrons)<<"   "
	     <<std::setprecision(3)
	     <<met<<"   "
	     <<std::setprecision(0)
	     <<nJets<<"   "
	     <<nCSVJets<<"   "
	     <<std::setprecision(3)
	     <<minDPhi<<"   "
	     <<std::setprecision(0)
	     <<passMETFilters<<"    "
	     <<dataTriggerHad
	     <<std::endl;
	 }

       }//end of CR1_lep
     
       if(CR2_lep){
	 syst0B.fillHistogramsSysts(h_mva_first_CR2,topFitRes.at(0).mva,w,systWeightsZeroB);
	 for ( const TopFitResults& res : topFitRes) {
	   syst0B.fillHistogramsSysts(h_bdt_qgid1_CR2,res.bdt_qgid1,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_qgid2_CR2,res.bdt_qgid2,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_dphij1b_CR2,res.bdt_dphij1b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_dphij2b_CR2,res.bdt_dphij2b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_drj1b_CR2,res.bdt_drj1b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_drj2b_CR2,res.bdt_drj2b,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_bjcsv_CR2,res.bdt_bjcsv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_jet1csv_CR2,res.bdt_jet1csv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_jet2csv_CR2,res.bdt_jet2csv,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_bdt_prob_CR2,res.bdt_prob,w,systWeightsZeroB);
	 }   
	 float maxMVA = topFitRes.at(0).mva;
	 if(maxMVA>0.2){
	   syst0B.fillHistogramsSysts(h_metFinal_tag_CR2,met,w,systWeightsZeroB);	     
	 }//end 1-tag cat
	 else{
	   syst0B.fillHistogramsSysts(h_metFinal_untag_CR2,met,w,systWeightsZeroB);	     
	 }//end 0-tag cat
	 
	 syst0B.fillHistogramsSysts(h_metFinal_met_0btag,met,w,systWeightsZeroB);

	 if(nJets == 4) syst0B.fillHistogramsSysts(h_metFinal_met_0btag_4j,met,w,systWeightsZeroB);
	 if(nJets == 6) syst0B.fillHistogramsSysts(h_metFinal_met_0btag_6j,met,w,systWeightsZeroB);
	 
	 double mt_f = sqrt(2* lep1.Pt() * metPt[0] * ( 1- cos(fabs(deltaPhi(lep1.Phi(), metPhi[0])))));
	 syst0B.fillHistogramsSysts(h_jet1Pt_CRwj_lep,jetPt[0],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_jet2Pt_CRwj_lep,jetPt[1],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_jet3Pt_CRwj_lep,jetPt[2],w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_nJets_CRwj_lep,nJets,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_nbJets_CRwj_lep,nCSVJets,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_dphi_6j_CRwj_lep,minDPhi_6j,w,systWeightsZeroB);
	 //syst0B.fillHistogramsSysts(h_dphi_10j_CRwj_lep,minDPhi_10j,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_dphi_CRwj_lep,minDPhi,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_mt_CRwj_lep,mt_f,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_mt2w_CRwj_lep,mt2w,w,systWeightsZeroB);
	 syst0B.fillHistogramsSysts(h_lepPt_CRwj_lep,lep1.Pt(),w,systWeightsZeroB);

	 if(lep1Flavour==13) {
	   syst0B.fillHistogramsSysts(h_muonmetFinal_met_0btag,met,w,systWeightsZeroB);

	   syst0B.fillHistogramsSysts(h_mu1Pt_CR2,lep1Pt,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet1Pt_CR2,jetPt[0],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet2Pt_CR2,jetPt[1],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonjet3Pt_CR2,jetPt[2],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonnJets_CR2,nJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muonnbJets_CR2,nCSVJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muondphi_6j_CR2,minDPhi_6j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_muondphi_10j_CR2,minDPhi_10j,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_muondphi_CR2,minDPhi,w,systWeightsZeroB);
	 }

	 if(lep1Flavour==11) {
	   syst0B.fillHistogramsSysts(h_electronmetFinal_met_0btag,met,w,systWeightsZeroB);

	   syst0B.fillHistogramsSysts(h_el1Pt_CR2,lep1Pt,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronjet1Pt_CR2,jetPt[0],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronjet2Pt_CR2,jetPt[1],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronjet3Pt_CR2,jetPt[2],w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronnJets_CR2,nJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electronnbJets_CR2,nCSVJets,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electrondphi_6j_CR2,minDPhi_6j,w,systWeightsZeroB);
	   //syst0B.fillHistogramsSysts(h_electrondphi_10j_CR2,minDPhi_10j,w,systWeightsZeroB);
	   syst0B.fillHistogramsSysts(h_electrondphi_CR2,minDPhi,w,systWeightsZeroB);
	 }

	 if(met>160 && sync=="sync"){
	   /* Synchro exercise */
	   fileout_CR2lep
	     <<std::fixed<<std::setprecision(0)
	     <<runNumber<<"   "
	     <<evtNumber<<"   "
	     <<lumiSec<<"   "
	     <<(nLooseMuons+nVetoElectrons)<<"   "
	     <<std::setprecision(3)
	     <<met<<"   "
	     <<std::setprecision(0)
	     <<nJets<<"   "
	     <<nCSVJets<<"   "
	     <<std::setprecision(3)
	     <<jetPt[0]<<"   "
	     <<jetPt[1]<<"   "
	     <<jetPt[2]<<"   "
	     <<minDPhi<<"   "
	     <<k_fact<<"   "
	     <<std::setprecision(0)
	     <<passMETFilters<<"    "
	     <<dataTriggerSL
	     <<std::endl;
	 }
       }
       
       // **********************************
       // Cut flow for semileptonic analysis
       // **********************************
       
       if(semileptonic){ 
	 n_binlep += w;
	 if(nJets>2){ 
	   n_binjet += w;
	   if(nCSVJets >0){ 
	     n_binbjet += w;
	     if(met>160){ 
	       n_binmet160 += w;
	       if(met>320){ 
		 n_binmet320 += w;
		 if(mt > 160){
		   n_binmt += w;
		   if(mt2w > 200 and mt2w<9999){
		     n_binmt2w += w;
		     if(minDPhi > 1.2){
		       n_phi += w;
		       if(nType1 ==0 && nType2 == 0)n_noBoost += w;
		     }
		   }
		 }
	       }
	     }
	   }
	 }
       }//end of Counting events for semileptonic channel
       
     } // end of if pass trigger
     
     // Filling up the tree per each event
     ttDMTree.Fill();
     
   }//end of loop over events 
 // cout << "pass " << pass << endl;
 // Cut Flow
 h_cutFlow->SetBinContent(0, nEvents);
 h_cutFlow->SetBinContent(1, n_trig);
 h_cutFlow->GetXaxis()->SetBinLabel(1, "trigger");
 h_cutFlow->SetBinContent(2, n_binlep);
 h_cutFlow->GetXaxis()->SetBinLabel(2, "lepton cuts");
 h_cutFlow->SetBinContent(3, n_binjet);
 h_cutFlow->GetXaxis()->SetBinLabel(3, "jet cuts");
 h_cutFlow->SetBinContent(4, n_binbjet);
 h_cutFlow->GetXaxis()->SetBinLabel(4, "b jet cuts");
 h_cutFlow->SetBinContent(5, n_binmet160);
 h_cutFlow->GetXaxis()->SetBinLabel(5, "MET > 160 GeV");
 h_cutFlow->SetBinContent(6, n_binmet320);
 h_cutFlow->GetXaxis()->SetBinLabel(6, "MET > 320 GeV");
 if(channel == "semileptonic"){
   h_cutFlow->SetBinContent(7, n_binmt);
   h_cutFlow->GetXaxis()->SetBinLabel(7, "M_{T} >  160 GeV");
   h_cutFlow->SetBinContent(8, n_binmt2w);
   h_cutFlow->GetXaxis()->SetBinLabel(8, "M_{T2}^{W} >  200 GeV");
   h_cutFlow->SetBinContent(9, n_phi);
   h_cutFlow->GetXaxis()->SetBinLabel(9, "#Delta #phi (j_{1,2},E^{miss}_{T}) >  1.2");
   h_cutFlow->SetBinContent(10, n_noBoost);
   h_cutFlow->GetXaxis()->SetBinLabel(10, "No Type1 or Type2");
 }
 else if(channel == "fullhadronic"){
   h_cutFlow->SetBinContent(7, n_phi);
   h_cutFlow->GetXaxis()->SetBinLabel(7, "#Delta #phi (j_{1,2},E^{miss}_{T}) >  1.2");
   h_cutFlow->SetBinContent(8, n_noBoost);
   h_cutFlow->GetXaxis()->SetBinLabel(8, "No Type1 or Type2");
   }     
 

 


 h_boostCat->SetBinContent(1,n_boostFullRes );
 h_boostCat->GetXaxis()->SetBinLabel(1, "FullRes");
 h_boostCat->SetBinContent(2,n_boost2Res );
 h_boostCat->GetXaxis()->SetBinLabel(2, "2Res");
 h_boostCat->SetBinContent(3,n_boost1Res );
 h_boostCat->GetXaxis()->SetBinLabel(3, "1Res");
 if(channel == "fullhadronic"){
   h_boostCat->SetBinContent(4,n_boost22 );
   h_boostCat->GetXaxis()->SetBinLabel(4, "22");
   h_boostCat->SetBinContent(5,n_boost12 );
   h_boostCat->GetXaxis()->SetBinLabel(5, "12");
   h_boostCat->SetBinContent(6,n_boost11 );
   h_boostCat->GetXaxis()->SetBinLabel(6, "11");
 }

 //cout << "test"<<endl;;
 // Saving tree

 outTree->cd();
 // cout << "test"<<saveTree<<endl;;

 // if(saveTree){
 ttDMTree.Write();
   // }
 
 outTree->Close();

 // Saving histograms
 // fout.cd();
 // writeHistogramsSysts(h_topHadPhiBM, allMyFiles);
 // writeHistogramsSysts(h_topHadPhiWM, allMyFiles);
 // writeHistogramsSysts(h_topHadPhiTM, allMyFiles);
 // writeHistogramsSysts(h_topHadPhiWB, allMyFiles);
 
 // writeHistogramsSysts(h_topHadCosBM, allMyFiles);
 // writeHistogramsSysts(h_topHadCosWM, allMyFiles);
 // writeHistogramsSysts(h_topHadCosTM, allMyFiles);
 // writeHistogramsSysts(h_topHadCosWB, allMyFiles);
 cout << " now writing "<<endl;
 systZero.writeHistogramsSysts(h_minDPhi, allMyFiles);
 systZero.writeHistogramsSysts(h_minDPhi_6j, allMyFiles);
 //systZero.writeHistogramsSysts(h_minDPhi_10j, allMyFiles);
		      
 // systZero.writeHistogramsSysts(h_topHadPt, allMyFiles);
 // systZero.writeHistogramsSysts(h_topHadMass, allMyFiles);

 systZero.writeHistogramsSysts(h_topHadTagCosTM, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadTagCosBM, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadTagCosWM, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadTagCosWB, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadTagCosBB, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_Angular_tag, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_Angular_untag, allMyFiles);

 systZero.writeHistogramsSysts(h_topHadTagCosTM_CR3, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadTagCosBM_CR3, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadTagCosWM_CR3, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadTagCosWB_CR3, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadTagCosBB_CR3, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_Angular_tag_CR3, allMyFiles);

 systZero.writeHistogramsSysts(h_topHadTagCosTM_lep, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadTagCosBM_lep, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadTagCosWM_lep, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadTagCosWB_lep, allMyFiles);
 // systZero.writeHistogramsSysts(h_topHadTagCosBB_lep, allMyFiles);
 // systZero.writeHistogramsSysts(h_metFinal_Angular_tag_lep, allMyFiles);
 // systZero.writeHistogramsSysts(h_metFinal_Angular_untag_lep, allMyFiles);
   
 systZero.writeHistogramsSysts(h_topHadCosBM_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadCosWM_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadCosTM_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadCosWB_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_minDPhi_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_minDPhi_6j_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadPt_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadMass_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topCosBB_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topEtaBB_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topCosBB_MET_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadCosTT_MET_allCuts, allMyFiles);
 //systZero.writeHistogramsSysts(h_topHadChi, allMyFiles);
 // systZero.writeHistogramsSysts(h_topSLPhiBM, allMyFiles);
 // systZero.writeHistogramsSysts(h_topSLPhiLM, allMyFiles);
 // systZero.writeHistogramsSysts(h_topSLPhiTM, allMyFiles);
 // systZero.writeHistogramsSysts(h_topSLPhiLBM, allMyFiles);
 // systZero.writeHistogramsSysts(h_topSLCosBM, allMyFiles);
 // systZero.writeHistogramsSysts(h_topSLCosLM, allMyFiles);
 // systZero.writeHistogramsSysts(h_topSLCosTM, allMyFiles);
 // systZero.writeHistogramsSysts(h_topSLCosLBM, allMyFiles);
 // systZero.writeHistogramsSysts(h_topSLPt, allMyFiles);
 // systZero.writeHistogramsSysts(h_topSLMass, allMyFiles);
 // systZero.writeHistogramsSysts(h_topSLMT, allMyFiles);

 systZero.writeHistogramsSysts(h_topSLCosBM_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLCosLM_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLCosTM_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLCosLBM_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLPt_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLMass_allCuts, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLMT_allCuts, allMyFiles);
 //systZero.writeHistogramsSysts(h_topSLChi, allMyFiles);
 
 systZero.writeHistogramsSysts(h_met, allMyFiles);

 systZero.writeHistogramsSysts(h_metFinal, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_4j, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_6j, allMyFiles);

 systZero.writeHistogramsSysts(h_metFinal_tag, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_1tag, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_untag, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_tag_CR1, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_untag_CR1, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_tag_CR2, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_untag_CR2, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_tag_CR3, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_1tag_CR3, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_untag_CR3, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_tag_CR4, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_1tag_CR4, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_untag_CR4, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_tag_CR5, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_1tag_CR5, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_untag_CR5, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_tag_CR6, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_1tag_CR6, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_untag_CR6, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_tag_CR7, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_1tag_CR7, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_untag_CR7, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_had45Jets, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_had6Jets, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_sl3Jets, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_sl4Jets, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal11, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal12, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal22, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal1Res, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal2Res, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_noBoost, allMyFiles);
 
 systZero.writeHistogramsSysts(h_mva, allMyFiles);
 systZero.writeHistogramsSysts(h_mva_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topMassPreFit, allMyFiles);
 systZero.writeHistogramsSysts(h_topMassPostFit, allMyFiles);
 systZero.writeHistogramsSysts(h_topMassPreFit_1tag, allMyFiles);
 systZero.writeHistogramsSysts(h_topMassPostFit_1tag, allMyFiles);
 systZero.writeHistogramsSysts(h_topMassPostFit_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topMassPreFit_preS, allMyFiles);

 systZero.writeHistogramsSysts(h_dphi, allMyFiles);
 systZero.writeHistogramsSysts(h_dphi_6j, allMyFiles);
 //systZero.writeHistogramsSysts(h_dphi_10j, allMyFiles);
 
 systZero.writeSingleHistogramSysts(h_cutFlow, allMyFiles);
 systZero.writeSingleHistogramSysts(h_boostCat, allMyFiles);
 // h_cutFlow->Write();
 // h_boostCat->Write();
 systZero.writeHistogramsSysts(h_metFinal_Angular, allMyFiles);
 systZero.writeHistogramsSysts(h_mt, allMyFiles);
 systZero.writeHistogramsSysts(h_mt2w, allMyFiles);
 systZero.writeHistogramsSysts(h_jet1Pt, allMyFiles);
 systZero.writeHistogramsSysts(h_jet2Pt, allMyFiles);
 systZero.writeHistogramsSysts(h_jet3Pt, allMyFiles);
 
 systZero.writeHistogramsSysts(h_bjet1Pt, allMyFiles);
 systZero.writeHistogramsSysts(h_bjet2Pt, allMyFiles);
 systZero.writeHistogramsSysts(h_bjet3Pt, allMyFiles);
 
 systZero.writeHistogramsSysts(h_jet1Pt_fin, allMyFiles);
 systZero.writeHistogramsSysts(h_jet2Pt_fin, allMyFiles);
 systZero.writeHistogramsSysts(h_jet3Pt_fin, allMyFiles);
 systZero.writeHistogramsSysts(h_bjetsPt, allMyFiles);
 systZero.writeHistogramsSysts(h_nJets, allMyFiles);
 systZero.writeHistogramsSysts(h_nbJets, allMyFiles);
 systZero.writeHistogramsSysts(h_nJets_fin, allMyFiles);
 systZero.writeHistogramsSysts(h_nbJets_fin, allMyFiles);

 //pre-selection
 //systZero.writeHistogramsSysts(h_jetQGL_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jetQGL1_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jetQGL2_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadmassDrop_preS, allMyFiles); 

 systZero.writeHistogramsSysts(h_nPV_preS, allMyFiles);   
 systZero.writeHistogramsSysts(h_nGoodPV_preS, allMyFiles); 
 systZero.writeHistogramsSysts(h_nPV_w_preS, allMyFiles);   
 systZero.writeHistogramsSysts(h_nGoodPV_w_preS, allMyFiles); 
 
 systZero.writeHistogramsSysts(h_mt_EC_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_mt_BC_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_mt_EC_hadpreS, allMyFiles);
 systZero.writeHistogramsSysts(h_mt_BC_hadpreS, allMyFiles);

 systZero.writeHistogramsSysts(h_nJets_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_nbJets_preS, allMyFiles);
 
 systZero.writeHistogramsSysts(h_jet1Pt_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jet2Pt_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jet3Pt_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jet1Eta_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jet2Eta_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jet3Eta_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jet1Phi_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jet2Phi_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jet3Phi_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jet1CSV_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jet2CSV_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_jet3CSV_preS, allMyFiles);
 
 systZero.writeHistogramsSysts(h_bjet1Pt_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_bjet2Pt_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_bjet3Pt_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_bjet1Eta_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_bjet2Eta_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_bjet3Eta_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_bjet1Phi_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_bjet2Phi_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_bjet3Phi_preS, allMyFiles);

 systZero.writeHistogramsSysts(h_met_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_dphi_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_dphi_6j_preS, allMyFiles);
 //systZero.writeHistogramsSysts(h_dphi_10j_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_mt_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_mt_hadpreS, allMyFiles);
 systZero.writeHistogramsSysts(h_mt2w_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_nJets_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_nbJets_preS, allMyFiles);

 systZero.writeHistogramsSysts(h_topSLMT_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLCosBM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLCosLM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLCosTM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLCosLBM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLPt_preS, allMyFiles); 
 systZero.writeHistogramsSysts(h_topSLMass_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLMass_preS_NS, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLMass_preMET, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLPhiBM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLPhiLM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLPhiTM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topSLPhiLBM_preS, allMyFiles);
 //systZero.writeHistogramsSysts(h_topSLChi_preS, allMyFiles);   

 systZero.writeHistogramsSysts(h_topHadCosBM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadCosWM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadCosTM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadCosWB_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadPt_preS, allMyFiles); 
 systZero.writeHistogramsSysts(h_topHadMass_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadMass_preMET, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadPhiBM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadPhiWM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadPhiTM_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadPhiWB_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topCosBB_preS, allMyFiles);   
 systZero.writeHistogramsSysts(h_topEtaBB_preS, allMyFiles);   
 systZero.writeHistogramsSysts(h_topCosBB_MET_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadCosTT_MET_preS, allMyFiles);
 systZero.writeHistogramsSysts(h_topHadChi_preS, allMyFiles);    

 systZero.writeHistogramsSysts(h_nPV, allMyFiles);
 systZero.writeHistogramsSysts(h_nPV_w, allMyFiles);
 
 systZero.writeHistogramsSysts(h_nGoodPV, allMyFiles);
 systZero.writeHistogramsSysts(h_nGoodPV_w, allMyFiles);
 
 systZero.writeHistogramsSysts(h_mva_first, allMyFiles);
 systZero.writeHistogramsSysts(h_bdt_qgid1, allMyFiles);
 systZero.writeHistogramsSysts(h_bdt_qgid2, allMyFiles);
 systZero.writeHistogramsSysts(h_bdt_dphij1b, allMyFiles);
 systZero.writeHistogramsSysts(h_bdt_dphij2b, allMyFiles);
 systZero.writeHistogramsSysts(h_bdt_jet1csv, allMyFiles);
 systZero.writeHistogramsSysts(h_bdt_jet2csv, allMyFiles);
 systZero.writeHistogramsSysts(h_bdt_bjcsv, allMyFiles);
 systZero.writeHistogramsSysts(h_bdt_prob, allMyFiles);

 systZero.writeHistogramsSysts(h_metFinal_SR_1lep_4j, allMyFiles);
 systZero.writeHistogramsSysts(h_metFinal_SR_1lep_6j, allMyFiles);

  
 //for bkg estimation
 if(channel == "fullhadronic"){

   systZero.writeHistogramsSysts(h_lep1Pt_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_lep1Eta_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_lep1Phi_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_lep1E_CR3, allMyFiles);
  
   systZero.writeHistogramsSysts(h_mu1Pt_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_mu1Eta_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_mu1Phi_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_mu1E_CR3, allMyFiles);
   
   systZero.writeHistogramsSysts(h_el1Pt_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_el1Eta_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_el1Phi_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_el1E_CR3, allMyFiles);
  
   systZero.writeHistogramsSysts(h_metFinal_outtop, allMyFiles);

   systZero.writeHistogramsSysts(h_metFinal_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_CR5_4j, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_CR5_6j, allMyFiles);

   systZero.writeHistogramsSysts(h_metFinal_CR6, allMyFiles);

   systZero.writeHistogramsSysts(h_metFinal_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_CR6nw_4j, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_CR6nw_6j, allMyFiles);

   systZero.writeHistogramsSysts(h_metFinal_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_CR7nw_4j, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_CR7nw_6j, allMyFiles);

   systZero.writeHistogramsSysts(h_muonmetFinal_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmetFinal_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_mixmetFinal_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmetFinal_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmetFinal_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mixmetFinal_CR7nw, allMyFiles);

   systZero.writeHistogramsSysts(h_metFinal_SR_1lep, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_SR_1lep_4j, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_SR_1lep_6j, allMyFiles);

   systZero.writeHistogramsSysts(h_muonmetFinal_SR_1lep, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmetFinal_SR_1lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_metFinal_SR_1lep_mt, allMyFiles);

   //CR3
   systZero.writeHistogramsSysts(h_jet1Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_jet2Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_jet3Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_nJets_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_nbJets_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_lepPt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_topSLMass_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_6j_CRtt_had, allMyFiles);
   //systZero.writeHistogramsSysts(h_dphi_10j_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_mt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_mt2w_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_bjet1Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_bjet2Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_bjet3Pt_CRtt_had, allMyFiles);

   systZero.writeHistogramsSysts(h_muonjet1Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet2Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet3Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnJets_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnbJets_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_muonlepPt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_6j_CRtt_had, allMyFiles);
   //systZero.writeHistogramsSysts(h_muondphi_10j_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmt2w_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet1Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet2Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet3Pt_CRtt_had, allMyFiles);

   systZero.writeHistogramsSysts(h_electronjet1Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet2Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet3Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnJets_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnbJets_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_electronlepPt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_6j_CRtt_had, allMyFiles);
   //systZero.writeHistogramsSysts(h_electrondphi_10j_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmt2w_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet1Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet2Pt_CRtt_had, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet3Pt_CRtt_had, allMyFiles);

   //CR3_nw
   systZero.writeHistogramsSysts(h_metFinal_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmetFinal_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmetFinal_CR3nw, allMyFiles);

   systZero.writeHistogramsSysts(h_lep1Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_lep1Eta_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_lep1Phi_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_lep1E_CR3nw, allMyFiles);
  
   systZero.writeHistogramsSysts(h_mu1Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mu1Eta_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mu1Phi_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mu1E_CR3nw, allMyFiles);
   
   systZero.writeHistogramsSysts(h_el1Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_el1Eta_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_el1Phi_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_el1E_CR3nw, allMyFiles);

   systZero.writeHistogramsSysts(h_jet1Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_jet2Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_jet3Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_nJets_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_nbJets_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_lepPt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_topSLMass_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_6j_CR3nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_dphi_10j_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mt2w_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_bjet1Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_bjet2Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_bjet3Pt_CR3nw, allMyFiles);

   systZero.writeHistogramsSysts(h_muonjet1Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet2Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet3Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnJets_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnbJets_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonlepPt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_6j_CR3nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_muondphi_10j_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmt2w_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet1Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet2Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet3Pt_CR3nw, allMyFiles);

   systZero.writeHistogramsSysts(h_electronjet1Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet2Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet3Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnJets_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnbJets_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_6j_CR3nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_electrondphi_10j_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmt2w_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet1Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet2Pt_CR3nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet3Pt_CR3nw, allMyFiles);

   //CR4
   systZero.writeHistogramsSysts(h_jet1Pt_CRqcd_had, allMyFiles);
   systZero.writeHistogramsSysts(h_jet2Pt_CRqcd_had, allMyFiles);
   systZero.writeHistogramsSysts(h_jet3Pt_CRqcd_had, allMyFiles);
   systZero.writeHistogramsSysts(h_nJets_CRqcd_had, allMyFiles);
   systZero.writeHistogramsSysts(h_nbJets_CRqcd_had, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_CRqcd_had, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_6j_CRqcd_had, allMyFiles);
   //systZero.writeHistogramsSysts(h_dphi_10j_CRqcd_had, allMyFiles);
   systZero.writeHistogramsSysts(h_bjet1Pt_CRqcd_had, allMyFiles);
   systZero.writeHistogramsSysts(h_bjet2Pt_CRqcd_had, allMyFiles);
   systZero.writeHistogramsSysts(h_bjet3Pt_CRqcd_had, allMyFiles);

   //CR5
   systZero.writeHistogramsSysts(h_jet1Pt_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_jet2Pt_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_jet3Pt_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_nJets_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_nbJets_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_6j_CR5, allMyFiles);
   //systZero.writeHistogramsSysts(h_dphi_10j_CR5, allMyFiles);

   //CR6
   systZero.writeHistogramsSysts(h_muonmetFinal_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_mu1Pt_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet1Pt_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet2Pt_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet3Pt_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnJets_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnbJets_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_6j_CR6, allMyFiles);
   //systZero.writeHistogramsSysts(h_muondphi_10j_CR6, allMyFiles);

   systZero.writeHistogramsSysts(h_electronmetFinal_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_el1Pt_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet1Pt_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet2Pt_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet3Pt_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnJets_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnbJets_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_6j_CR6, allMyFiles);
   //systZero.writeHistogramsSysts(h_electrondphi_10j_CR6, allMyFiles);

   systZero.writeHistogramsSysts(h_jet1Pt_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_jet2Pt_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_jet3Pt_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_nJets_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_nbJets_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_6j_CR6, allMyFiles);
   //systZero.writeHistogramsSysts(h_dphi_10j_CR6, allMyFiles);

   //CR6nw
   systZero.writeHistogramsSysts(h_muonmetFinal_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mu1Pt_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet1Pt_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet2Pt_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet3Pt_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnJets_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnbJets_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_6j_CR6nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_muondphi_10j_CR6nw, allMyFiles);

   systZero.writeHistogramsSysts(h_electronmetFinal_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_el1Pt_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet1Pt_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet2Pt_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet3Pt_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnJets_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnbJets_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_6j_CR6nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_electrondphi_10j_CR6nw, allMyFiles);

   systZero.writeHistogramsSysts(h_jet1Pt_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_jet2Pt_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_jet3Pt_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_nJets_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_nbJets_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_CR6nw, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_6j_CR6nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_dphi_10j_CR6nw, allMyFiles);

   //CR7
   systZero.writeHistogramsSysts(h_jet1Pt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_jet2Pt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_jet3Pt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_nJets_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_nbJets_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_6j_CR7, allMyFiles);
   //systZero.writeHistogramsSysts(h_dphi_10j_CR7, allMyFiles);
   //systZero.writeHistogramsSysts(h_mt_CR7, allMyFiles);

   systZero.writeHistogramsSysts(h_mu1Pt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_mu2Pt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_el1Pt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_el2Pt_CR7, allMyFiles);

   systZero.writeHistogramsSysts(h_muondphi_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_6j_CR7, allMyFiles);
   //systZero.writeHistogramsSysts(h_muondphi_10j_CR7, allMyFiles);
   //systZero.writeHistogramsSysts(h_muonmt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnJets_CR7 , allMyFiles);
   systZero.writeHistogramsSysts(h_muonnbJets_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet1Pt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet2Pt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet3Pt_CR7, allMyFiles);

   systZero.writeHistogramsSysts(h_electrondphi_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_6j_CR7, allMyFiles);
   //systZero.writeHistogramsSysts(h_electrondphi_10j_CR7, allMyFiles);
   //systZero.writeHistogramsSysts(h_electronmt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnJets_CR7 , allMyFiles);
   systZero.writeHistogramsSysts(h_electronnbJets_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet1Pt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet2Pt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet3Pt_CR7, allMyFiles);

   systZero.writeHistogramsSysts(h_mixdphi_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_mixdphi_6j_CR7, allMyFiles);
   //systZero.writeHistogramsSysts(h_mixdphi_10j_CR7, allMyFiles);
   //systZero.writeHistogramsSysts(h_mixmt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_mixnJets_CR7 , allMyFiles);
   systZero.writeHistogramsSysts(h_mixnbJets_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_mixjet1Pt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_mixjet2Pt_CR7, allMyFiles);
   systZero.writeHistogramsSysts(h_mixjet3Pt_CR7, allMyFiles);

   //CR7nw
   systZero.writeHistogramsSysts(h_jet1Pt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_jet2Pt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_jet3Pt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_nJets_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_nbJets_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_6j_CR7nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_dphi_10j_CR7nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_mt_CR7nw, allMyFiles);

   systZero.writeHistogramsSysts(h_mu1Pt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mu2Pt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_el1Pt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_el2Pt_CR7nw, allMyFiles);

   systZero.writeHistogramsSysts(h_muondphi_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_6j_CR7nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_muondphi_10j_CR7nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_muonmt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnJets_CR7nw , allMyFiles);
   systZero.writeHistogramsSysts(h_muonnbJets_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet1Pt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet2Pt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet3Pt_CR7nw, allMyFiles);

   systZero.writeHistogramsSysts(h_electrondphi_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_6j_CR7nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_electrondphi_10j_CR7nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_electronmt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnJets_CR7nw , allMyFiles);
   systZero.writeHistogramsSysts(h_electronnbJets_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet1Pt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet2Pt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet3Pt_CR7nw, allMyFiles);

   systZero.writeHistogramsSysts(h_mixdphi_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mixdphi_6j_CR7nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_mixdphi_10j_CR7nw, allMyFiles);
   //systZero.writeHistogramsSysts(h_mixmt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mixnJets_CR7nw , allMyFiles);
   systZero.writeHistogramsSysts(h_mixnbJets_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mixjet1Pt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mixjet2Pt_CR7nw, allMyFiles);
   systZero.writeHistogramsSysts(h_mixjet3Pt_CR7nw, allMyFiles);
   
   // QCD CR in fullhadronic
   systZero.writeHistogramsSysts(h_mva_first_CR4, allMyFiles);
   systZero.writeHistogramsSysts(h_mva_second_CR4, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid1_CR4, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid2_CR4, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij1b_CR4, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij2b_CR4, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet1csv_CR4, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet2csv_CR4, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_bjcsv_CR4, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_prob_CR4, allMyFiles); 
   // tt CR in fullhadronic
   systZero.writeHistogramsSysts(h_mva_first_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_mva_second_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid1_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid2_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij1b_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij2b_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet1csv_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet2csv_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_bjcsv_CR3, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_prob_CR3, allMyFiles); 
   // V+Jets CR in fullhadronic
   systZero.writeHistogramsSysts(h_mva_first_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_mva_second_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid1_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid2_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij1b_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij2b_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet1csv_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet2csv_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_bjcsv_CR5, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_prob_CR5, allMyFiles); 
   // W+Jets CR in fullhadronic
   systZero.writeHistogramsSysts(h_mva_first_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_mva_second_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid1_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid2_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij1b_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij2b_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet1csv_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet2csv_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_bjcsv_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_prob_CR6, allMyFiles); 
   // Z+Jets CR in fullhadronic
   systZero.writeHistogramsSysts(h_mva_first_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_mva_second_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid1_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid2_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij1b_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij2b_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet1csv_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet2csv_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_bjcsv_CR6, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_prob_CR6, allMyFiles);    
 }
  
 if(channel == "semileptonic"){

   systZero.writeHistogramsSysts(h_metFinal_noMT2W, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmet, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmetFinal, allMyFiles);

   systZero.writeHistogramsSysts(h_muonmet_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_6j_preS, allMyFiles);
   //systZero.writeHistogramsSysts(h_muondphi_10j_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmt2w_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_muonnJets_preS , allMyFiles);
   systZero.writeHistogramsSysts(h_muonnbJets_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_muonjet1Pt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet2Pt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet3Pt_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_muonjet1Eta_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet2Eta_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet3Eta_preS, allMyFiles);
   
   systZero.writeHistogramsSysts(h_muonjet1Phi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet2Phi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet3Phi_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_muonjet1CSV_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet2CSV_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet3CSV_preS, allMyFiles);
 
   systZero.writeHistogramsSysts(h_muonbjet1Pt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet2Pt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet3Pt_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_muonbjet1Eta_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet2Eta_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet3Eta_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_muonbjet1Phi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet2Phi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet3Phi_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_muondphi, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_6j, allMyFiles);
   //systZero.writeHistogramsSysts(h_muondphi_10j, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmt, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmt2w, allMyFiles);

   systZero.writeHistogramsSysts(h_muonmetFinal_2lep, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmetFinal_2lep_Z_nobtag, allMyFiles);
   systZero.writeHistogramsSysts(h_mixmetFinal_2lep, allMyFiles);
   systZero.writeHistogramsSysts(h_mixmetFinal_2lep_Z_nobtag, allMyFiles);

   systZero.writeHistogramsSysts(h_electronmet, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmetFinal, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmetFinal_2lep, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmetFinal_2lep_Z_nobtag, allMyFiles);

   systZero.writeHistogramsSysts(h_electronmet_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_6j_preS, allMyFiles);
   //systZero.writeHistogramsSysts(h_electrondphi_10j_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmt2w_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_electronnJets_preS , allMyFiles);
   systZero.writeHistogramsSysts(h_electronnbJets_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_electronjet1Pt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet2Pt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet3Pt_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_electronjet1Eta_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet2Eta_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet3Eta_preS, allMyFiles);
   
   systZero.writeHistogramsSysts(h_electronjet1Phi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet2Phi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet3Phi_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_electronjet1CSV_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet2CSV_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet3CSV_preS, allMyFiles);
 
   systZero.writeHistogramsSysts(h_electronbjet1Pt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet2Pt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet3Pt_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_electronbjet1Eta_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet2Eta_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet3Eta_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_electronbjet1Phi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet2Phi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet3Phi_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_lep1Pt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_lep1Eta_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_lep1Phi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_lep1E_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_mu1Pt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_mu1Eta_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_mu1Phi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_mu1E_preS, allMyFiles);

   systZero.writeHistogramsSysts(h_el1Pt_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_el1Eta_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_el1Phi_preS, allMyFiles);
   systZero.writeHistogramsSysts(h_el1E_preS, allMyFiles);

   //CR1
   systZero.writeHistogramsSysts(h_lep1Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_lep1Eta_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_lep1Phi_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_lep1E_CR1, allMyFiles);

   systZero.writeHistogramsSysts(h_lep2Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_lep2Eta_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_lep2Phi_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_lep2E_CR1, allMyFiles);

   systZero.writeHistogramsSysts(h_mu1Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mu2Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_el1Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_el2Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mt2w_CR1, allMyFiles);
   
   systZero.writeHistogramsSysts(h_muondphi_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_6j_CR1, allMyFiles);
   //systZero.writeHistogramsSysts(h_muondphi_10j_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_muonmt2w_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnJets_CR1 , allMyFiles);
   systZero.writeHistogramsSysts(h_muonnbJets_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet1Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet2Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet3Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet1Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet2Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_muonbjet3Pt_CR1, allMyFiles);

   systZero.writeHistogramsSysts(h_electrondphi_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_6j_CR1, allMyFiles);
   //systZero.writeHistogramsSysts(h_electrondphi_10j_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmt2w_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnJets_CR1 , allMyFiles);
   systZero.writeHistogramsSysts(h_electronnbJets_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet1Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet2Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet3Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet1Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet2Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_electronbjet3Pt_CR1, allMyFiles);

   systZero.writeHistogramsSysts(h_mixdphi_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mixdphi_6j_CR1, allMyFiles);
   //systZero.writeHistogramsSysts(h_mixdphi_10j_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mixmt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mixmt2w_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mixnJets_CR1 , allMyFiles);
   systZero.writeHistogramsSysts(h_mixnbJets_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mixjet1Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mixjet2Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mixjet3Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mixbjet1Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mixbjet2Pt_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_mixbjet3Pt_CR1, allMyFiles);
   
   //CR2
   systZero.writeHistogramsSysts(h_metFinal_met_0btag, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_met_0btag_4j, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_met_0btag_6j, allMyFiles);

   systZero.writeHistogramsSysts(h_muonmetFinal_met_0btag, allMyFiles);
   systZero.writeHistogramsSysts(h_electronmetFinal_met_0btag, allMyFiles);

   systZero.writeHistogramsSysts(h_mu1Pt_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_muondphi_6j_CR2, allMyFiles);
   //systZero.writeHistogramsSysts(h_muondphi_10j_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnJets_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_muonnbJets_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet1Pt_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet2Pt_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_muonjet3Pt_CR2, allMyFiles);

   systZero.writeHistogramsSysts(h_el1Pt_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_electrondphi_6j_CR2, allMyFiles);
   //systZero.writeHistogramsSysts(h_electrondphi_10j_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnJets_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_electronnbJets_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet1Pt_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet2Pt_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_electronjet3Pt_CR2, allMyFiles);

   systZero.writeHistogramsSysts(h_metFinal_2lep, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_2lep_4j, allMyFiles);
   systZero.writeHistogramsSysts(h_metFinal_2lep_6j, allMyFiles);

   systZero.writeHistogramsSysts(h_jet1Pt_CRtt2_lep , allMyFiles);
   systZero.writeHistogramsSysts(h_jet2Pt_CRtt2_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_jet3Pt_CRtt2_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_nJets_CRtt2_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_nbJets_CRtt2_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_CRtt2_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_mt_CRtt2_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_mt2w_CRtt2_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_lepPt_CRtt2_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_topSLMass_CRtt2_lep, allMyFiles);

   systZero.writeHistogramsSysts(h_metFinal_2lep_Z_nobtag, allMyFiles);

   systZero.writeHistogramsSysts(h_jet1Pt_CRZ_lep , allMyFiles);
   systZero.writeHistogramsSysts(h_jet2Pt_CRZ_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_jet3Pt_CRZ_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_nJets_CRZ_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_nbJets_CRZ_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_CRZ_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_mt_CRZ_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_mt2w_CRZ_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_lepPt_CRZ_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_topSLMass_CRZ_lep, allMyFiles);

   systZero.writeHistogramsSysts(h_jet1Pt_CRwj_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_jet2Pt_CRwj_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_jet3Pt_CRwj_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_nJets_CRwj_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_nbJets_CRwj_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_CRwj_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_dphi_6j_CRwj_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_dphi_10j_CRwj_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_mt_CRwj_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_mt2w_CRwj_lep, allMyFiles);
   systZero.writeHistogramsSysts(h_lepPt_CRwj_lep, allMyFiles);

   //systZero.writeHistogramsSysts(h_jet1Pt_CRttFull_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_jet2Pt_CRttFull_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_jet3Pt_CRttFull_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_nJets_CRttFull_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_nbJets_CRttFull_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_dphi_CRttFull_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_mt_CRttFull_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_lepPt_CRttFull_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_lepPt_CRttFull_lep2, allMyFiles);
   //systZero.writeHistogramsSysts(h_topSLMass_CRttFull_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_jet1Pt_2lep_Full_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_jet2Pt_2lep_Full_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_jet3Pt_2lep_Full_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_nJets_2lep_Full_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_nbJets_2lep_Full_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_dphi_2lep_Full_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_mt2w_2lep_Full_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_lepPt_2lep_Full_lep, allMyFiles);
   //systZero.writeHistogramsSysts(h_topSLMass_2lep_Full_lep, allMyFiles);
   
   systZero.writeHistogramsSysts(h_mva_first_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid1_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid2_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij1b_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij2b_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet1csv_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet2csv_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_bjcsv_CR1, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_prob_CR1, allMyFiles);
   
   systZero.writeHistogramsSysts(h_mva_first_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid1_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_qgid2_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij1b_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_dphij2b_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet1csv_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_jet2csv_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_bjcsv_CR2, allMyFiles);
   systZero.writeHistogramsSysts(h_bdt_prob_CR2, allMyFiles); 
  
 }
 systZero.writeHistogramsSysts(h_topHadMass_allCutsMS, allMyFiles);

 // std::cout<<"Integral: "<<h_topHadPt->Integral()<<std::endl;
 // fout.Close();
 
 fileout.close();
 fileout_presel.close();
 fileout_CR3.close();
 fileout_CR5.close();
 fileout_CR6.close();
 fileout_CR7.close();
 //fileout_CRouttop.close();
 fileout_CRouttop.close();
 fileout_CR1lep.close();
 fileout_CR2lep.close();

   std::cout<< "---> "<<sample<<std::endl;
   std::cout<< "Number of events           : "<<nEvents<<std::endl;
   std::cout<< "Events after pre-sel.      : "<<nEventsPrePres<<std::endl;
   std::cout<< "Events after trigger cut   : "<<n_trig<<std::endl; 
   std::cout<< "Events after lepton cut    : " <<n_binlep<<std::endl;
   std::cout<< "Events after jet cut       : " <<n_binjet<<std::endl;
   std::cout<< "Events after bjet cut      : " <<n_binbjet<<std::endl;
   std::cout<< "Events after met160 cut    : " <<n_binmet160<<std::endl;
   std::cout<< "Events after met320 cut    : " <<n_binmet320<<std::endl;
   if(channel == "semileptonic"){
     std::cout<< "Events after mt cut        : " <<n_binmt<<std::endl; 
     std::cout<< "Events after mt2wcut       : " <<n_binmt2w<<std::endl; 
   }  
   std::cout<< "Events after minDPhi cut   : " <<n_phi<<std::endl;
   std::cout<< "Events after noBoost cut   : " <<n_noBoost<<std::endl;
   
   // std::cout<< "-" * 20<<std::endl;
   std::cout<< "Boosted Categories"<<std::endl;
   if(channel == "fullhadronic"){
     std::cout<< "Events in cat 11     : " <<n_boost11<<std::endl;
     std::cout<< "Events in cat 12     : " <<n_boost12<<std::endl;
     std::cout<< "Events in cat 22     : " <<n_boost22<<std::endl;
   }
   std::cout<< "Events in cat 1Res   : " <<n_boost1Res<<std::endl;
   std::cout<< "Events in cat 2Res   : " <<n_boost2Res<<std::endl;
   std::cout<< "Events in cat Res    : " <<n_boostFullRes<<std::endl;
   

   
   //return h

}//end of main

TH1F * initproduct(TH1F * hA,TH1F* hB, int rebinA = 1, int rebinB=1,double integral = -1.){
  int nbinsA = hA->GetNbinsX();
  int nbinsB = hA->GetNbinsX();
  double min = hA->GetBinLowEdge(1)*hB->GetBinLowEdge(1);
  double max = hA->GetBinLowEdge(nbinsA+1)*hB->GetBinLowEdge(nbinsB+1);
  //Get the actual name from the original histograms 
  string name =(string)(hA->GetName()) +"_vs_"+ (string)(hB->GetName());
  
  //Initialize histogram 
  TH1F * result = new TH1F(name.c_str(),name.c_str(),nbinsA*nbinsB,min,max);
  return result;
}

TH1F * makeproduct(TH1F * hA,TH1F* hB, int rebinA = 1, int rebinB=1,double integral = -1.){

  //Make temporary histos to rebin
  //  TH1F *hA = (TH1F*)h_A->Clone("hA");
  // TH1F *hB = (TH1F*)h_B->Clone("hB");

  //  hA->Rebin(rebinA);
  // hB->Rebin(rebinB);
  
  //get nbins from new histos
  int nbinsA = hA->GetNbinsX();
  int nbinsB = hA->GetNbinsX();
  double min = hA->GetBinLowEdge(1)*hB->GetBinLowEdge(1);
  double max = hA->GetBinLowEdge(nbinsA+1)*hB->GetBinLowEdge(nbinsB+1);
  //Get the actual name from the original histograms 
  string name =(string)(hA->GetName()) +"_vs_"+ (string)(hB->GetName());
  
  //Initialize histogram 
  TH1F * result = new TH1F(name.c_str(),name.c_str(),nbinsA*nbinsB,min,max);
  //Fill histogram
  for(int i =1; i<= nbinsA;++i){
    for(int j =1; j<= nbinsB;++j){
      double value = hA->GetBinContent(i)*hB->GetBinContent(j);
      int k = ((i-1)*nbinsB)+j;
      result->SetBinContent(k,value);
    }
  }
  if( integral <= 0.)integral = hB->Integral()/result->Integral();
  else integral = integral / result->Integral();
  result->Scale(integral);
  return result;

}

//void initHistogramsSysts(TH1F* histo, TString name, TString, int, float, float , bool useOnlyNominal=false);

void systWeights::copySysts(systWeights sys){
  for(int i =0; i < sys.maxSysts;++i){
    this->weightedNames[i]=sys.weightedNames[i];
    this->weightedSysts[i]=sys.weightedSysts[i];

  }
  this->setMax(sys.maxSysts);
  this->setMaxNonPDF(sys.maxSystsNonPDF);
  this->nPDF=sys.nPDF;
  this->nCategories=sys.nCategories;  
  this->addQ2=sys.addQ2;
  this->addPDF=sys.addPDF;
  this->addTopPt=sys.addTopPt;
  this->addVHF=sys.addVHF;
  this->addTTSplit=sys.addTTSplit;
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
    //    int MAX = this->maxSysts;
    this->weightedNames[0]="";
    this->weightedNames[1]="btagUp";
    this->weightedNames[2]="btagDown";
    this->weightedNames[3]="mistagUp";
    this->weightedNames[4]="mistagDown";
    this->weightedNames[5]="puUp";
    this->weightedNames[6]="puDown";
    this->weightedNames[7]="lepUp";
    this->weightedNames[8]="lepDown";
    //this->weightedNames[9]="isoUp";
    //this->weightedNames[10]="isoDown";
    //this->weightedNames[11]="trigUp";
    //this->weightedNames[12]="trigDown";
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
    //    this->weightedNames[this->maxSysts]="2lep";
    //    this->weightedNames[this->maxSysts+1]="1lep";
    //    this->weightedNames[this->maxSysts+2]="0lep";
    //    this->setMax(this->maxSysts+3);
    //    this->setMaxNonPDF(this->maxSystsNonPDF+3);
    //    this->weightedNames[this->maxSysts]= "";
    this->nCategories=4;
    categoriesNames[1]="TT0lep";
    categoriesNames[2]="TT1lep";
    categoriesNames[3]="TT2lep";
    this->wCats[1]=1.0;
    this->wCats[2]=1.0;
    this->wCats[3]=1.0;

  }


  /*  if(addkFact){
    this->weightedNames[this->maxSysts]="VHFWeightUp";
    this->weightedNames[this->maxSysts+1]="VHFWeightDown";
    this->setMax(this->maxSysts+2);
    this->setMaxNonPDF(this->maxSystsNonPDF+2);
    this->weightedNames[this->maxSysts]= "";
    }*/

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
void systWeights::addSyst(string name){
  this->weightedNames[this->maxSysts]= name;
  this->setMax(maxSysts+1);
  if(name.find("pdf")!=std::string::npos)this->setMaxNonPDF(maxSysts+1);
  this->weightedNames[this->maxSysts]= "";
}

void systWeights::addSystNonPDF(string name){
  this->weightedNames[this->maxSystsNonPDF]= name;
  this->setMaxNonPDF(maxSystsNonPDF+1);
  int nPDF=this->nPDF;
  for(int i =0; i < nPDF;++i){
    stringstream ss;
    ss<< i+1;
    this->weightedNames[i+this->maxSystsNonPDF]= "pdf"+ss.str();
  }
  this->setMax(maxSystsNonPDF+nPDF);
  this->weightedNames[this->maxSysts]= "";
}

void systWeights::addkFact(string name){
  string up=name+"Up";
  string down=name+"Down";
  cout << " adding syst "<< up<<endl;
  this->addSystNonPDF(up);
  this->addSystNonPDF(down);
}

void systWeights::setkFact(string name, float kfact_nom, float kfact_up,float kfact_down, bool mult){
  //  void setkFact(string name,float kfact_nom, float kfact_up,float kfact_down, float w_zero=1.0, mult=true);
  float zerofact=1.0;
  if(mult)zerofact=this->weightedSysts[0];
  string up = name+"Up";
  string down = name+"Down";
  float valueup=kfact_up/kfact_nom;
  float valuedown=kfact_down/kfact_nom;
  //  cout << "setting syst "<< up<<endl;
  //  cout << "values nom "<<kfact_nom<< " up "<< kfact_up << " down "<< kfact_down << " valup "<< valueup<< " valdown "<< valuedown <<" zerofact "<< zerofact<<endl;
  this->setSystValue(up, valueup*zerofact);
  this->setSystValue(down, valuedown*zerofact);

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
      //      cout << " initialized histogram "<< histo[sy+(MAX+1)*c]->GetName() <<" sy " << sy << " c  "<< c <<" location " << sy+(MAX+1)*c << endl;

    }
  }
}

void systWeights::setOnlyNominal(bool useOnlyNominal){
  this->onlyNominal=useOnlyNominal;
}

void systWeights::setWCats(double * wcats){
  for(int i =0;i<this->nCategories;++i){
    //    cout << "setting wcat #"<< i << " to be "<<wcats[i]<<endl;
    this->wCats[i]=wcats[i];
  }
 
}
void systWeights::fillHistogramsSysts(TH1F** histo, float v, float w,  float *systWeights, int nFirstSysts, double * wcats, bool verbose){
  if(wcats== NULL){
    wcats = this->wCats;
  }
	     
  for (int c = 0; c < this->nCategories; c++){
    int MAX = this->maxSysts;
    bool useOnlyNominal = this->onlyNominal;
    //cout << " filling histo " << histo[0+(MAX+1)*(c)]->GetName()<< " MAX "<<MAX*(1+c)<<" nFirstSysts"<< nFirstSysts<< endl;
    //    cout << "weight 0 "<< systWeights[0]<< " weighted syst 0 "<< this->weightedSysts[0]<<endl;
    for(int sy=0;sy<(int)MAX;++sy){
      if(sy!=0 && useOnlyNominal)continue;
      float ws=1.0;
      if(sy<nFirstSysts){
	ws=systWeights[sy]*wcats[c];
      }
      else {
	ws = (this->weightedSysts[(int)sy]*wcats[c]);
      }
      if(verbose) cout << "filling histo "<< histo[sy+((MAX+1)*(c))]->GetName()<<" value "<< v << " wevt "<< w << " syst number "<< sy<< " name "<< weightedNames[sy]<<" ws value " <<ws<<endl;
      histo[sy+((MAX+1)*(c))]->Fill(v, w*ws);
    }
  }
}

void systWeights::fillHistogramsSysts(TH1F** histo, float v, float w, double * wcats, bool verbose ){
  if(wcats==NULL){
    wcats=this->wCats;
  }
  for (int c = 0; c < this->nCategories; c++){
    int MAX = this->maxSysts;
    bool useOnlyNominal = this->onlyNominal;
    for(int sy=0;sy<(int)MAX;++sy){
      if(sy!=0&& useOnlyNominal)continue;
      float ws = (this->weightedSysts[(int)sy])*wcats[c];
      // cout << " filling histogram "<< histo[(int)sy]->GetName() << " with value "<< v <<" and weight "<< w <<" ws "<< ws<<endl;
      histo[sy+(MAX+1)*(c)]->Fill(v, w*ws);
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

      //    "/afs/cern.ch/user/o/oiorio/public/xAnnapaola/Nov10/res/"+sample + "_" +channel+".root";
      if(sy==0){
	//cout<<" filename is "<< basename+ns+cname+".root"<<endl;
	allFiles[sy+(MAX+1)*c]= TFile::Open((basename+ns+cname+".root"), opt);
      }
      else{
	if(!useOnlyNominal){
	  //if((ns!="1lep") && (ns!="2lep")&& (ns!="0lep")){
	  //	  cout<<" filename is "<< basename+ns+cname+".root"<<endl;
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

void systWeights::writeHistogramsSysts(TH1F** histo, TFile **filesout){  
  int MAX= this->maxSystsNonPDF;
  int MAXTOT= this->maxSysts;
  bool useOnlyNominal = this->onlyNominal;
  for (int c = 0; c < this->nCategories; c++){
    TString cname= (this->categoriesNames[c]).c_str();
    if (c!=0) cname= "_"+cname;
    for(int sy=0;sy<(int)MAX;++sy){
      //      cout << "c is now "<< c << " sy "<< sy << " location "<< sy+(MAXTOT+1)*c <<" is histo there? " << histo[sy+(MAXTOT+1)*c] << " file location "<<sy+(MAX+1)*c << " is file there "<< filesout[sy+(MAX+1)*c]<< endl;
      //      cout << " writing histo "<< histo[sy+(MAXTOT+1)*c]->GetName()<< " in file "<< filesout[sy+(MAX+1)*c]->GetName()<<endl;;
      //TString ns= weightedSystsNames((weightedSysts)sy);
      if(!(!useOnlyNominal || sy==0)) continue;
      
      filesout[(int)sy+(MAX+1)*(c)]->cd();
      if(this->addPDF){
	if(this->weightedNames[sy]=="pdf_totalUp")calcPDFHisto(histo, histo[sy+(MAXTOT+1)*(c)],1.0,c);
	if(this->weightedNames[sy]=="pdf_totalDown")calcPDFHisto(histo, histo[sy+(MAXTOT+1)*(c)],-1.0,c);
	;      //this->
      }
      
      histo[sy+(MAXTOT+1)*c]->Write(histo[0]->GetName());
      //    histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
      //    filesout[(int)sy]->Close();
    }
    if(this->addPDF){
      if(!useOnlyNominal){
	filesout[MAX+(MAX+1)*(c)]->cd();
	//	cout << " file max is "<< filesout[MAX+(MAX+1)*c]->GetName()<<endl;
	//	int npdf=this->maxSysts-this->maxSystsNonPdf;
	int MAXPDF=this->maxSysts;
	for(int sy=MAX;sy<MAXPDF;++sy){
	  //	  cout << " writing sy "<<sy+(MAXTOT+1)*c<<endl;
	  //	  cout << " histo is there? "<< histo[sy+(MAXTOT+1)*c]<<endl;
	  histo[sy+(MAXTOT+1)*(c)]->Write();
	  //	  cout << " written sy "<< histo[sy+(MAXTOT+1)*c]->GetName()<<endl;
	}
      }
    }
  }
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
	  //      cout << " writing sy "<< histo[sy]->GetName()<<endl;
	  histo->Write();
	  //      cout << " written sy "<< histo[sy]->GetName()<<endl;
	}
      }
    }
  }
}


void systWeights::setMax(int max){
  this->maxSysts =  max;
}
void systWeights::setMaxNonPDF(int max){
  this->maxSystsNonPDF =  max;
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
//void systWeights::setSystValue(weightedSysts systPlace, double value);

/*
void createFilesSysts(  TFile * allFiles[maxSysts], TString basename, bool useOnlyNominal=false,TString opt = "RECREATE"){
  for(int sy=0;sy<maxSysts;++sy){
    TString ns= weightedSystsNames((weightedSysts)sy);
    //    "/afs/cern.ch/user/o/oiorio/public/xAnnapaola/Nov10/res/"+sample + "_" +channel+".root";
    if(sy==0){
	allFiles[sy]= TFile::Open((basename+ns+".root"), opt);}
    else{
      allFiles[sy]= TFile::Open((basename+"_"+ns+".root"), opt);}
    //TFile *outTree = TFile::Open(("trees/tree_"+outFileName).c_str(), "RECREATE");
    
  }
  //return allFiles;
}

void systZero.writeHistogramsSysts(TH1F* histo[maxSysts], TFile *filesout[(int)MAXSYSTS], bool useOnlyNominal=false){  
  for(int sy=0;sy<(int)MAXSYSTS;++sy){
    //cout << " writing histo "<< histo[(int)sy]->GetName()<< " in file "<< filesout[(int)sy]->GetName()<<endl;;
    //TString ns= weightedSystsNames((weightedSysts)sy);
    filesout[(int)sy]->cd();
    histo[sy]->Write(histo[0]->GetName());
    //    histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
  }
}

void writeSingleHistogramSysts(TH1F* histo, TFile *filesout[(int)MAXSYSTS], bool useOnlyNominal=false){  
  for(int sy=0;sy<(int)MAXSYSTS;++sy){
    cout << " writing histo "<< histo->GetName()<< " in file "<< filesout[(int)sy]->GetName()<<endl;;
    //TString ns= weightedSystsNames((weightedSysts)sy);
    filesout[(int)sy]->cd();
    histo->Write();
    //    histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
  }
}
*/
//
//{
//  void createFileSyst();
//  void fillHistogramSyst();
//  void writeHistogramSyst();
//  void setMax(int max);
//  void setSystValue(string name, double value);
//  void setSystValue(weightedSysts systPlace, double value);
//  int maxSysts;
//  int[wLimit] weightedSysts;
//  string[wLimit] weightedNames;
//}


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

void  initHistogramsSysts (TH1F* histo[(int)MAXSYSTS],TString name, TString title, int nbins, float min, float max, bool useOnlyNominal=false){
  for(int sy=0;sy<(int)MAXSYSTS;++sy){
    TString ns= weightedSystsNames((weightedSysts)sy);
    histo[sy]=new TH1F(name+ns,title,nbins,min,max);
  }
}

void fillHistogramsSysts(TH1F* histo[(int)MAXSYSTS], float v, float w, float systWeight[(int)MAXSYSTS] , bool useOnlyNominal=false){
  for(int sy=0;sy<(int)MAXSYSTS;++sy){
    float ws = systWeight[(int)sy];
    //    cout << " filling histogram "<< histo[(int)sy]->GetName() << " with value "<< v <<" and weight "<< w <<" ws "<< ws<<endl;
    histo[sy]->Fill(v, w*ws);
  }
}

void createFilesSysts(  TFile * allFiles[(int)MAXSYSTS], TString basename, bool useOnlyNominal=false,TString opt = "RECREATE"){
  for(int sy=0;sy<(int)MAXSYSTS;++sy){
    TString ns= weightedSystsNames((weightedSysts)sy);
    //    "/afs/cern.ch/user/o/oiorio/public/xAnnapaola/Nov10/res/"+sample + "_" +channel+".root";
    if(sy==0){
	allFiles[sy]= TFile::Open((basename+ns+".root"), opt);}
    else{
      if(!useOnlyNominal) allFiles[sy]= TFile::Open((basename+"_"+ns+".root"), opt);}
    //TFile *outTree = TFile::Open(("trees/tree_"+outFileName).c_str(), "RECREATE");
    
  }
  //return allFiles;
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

void writeSingleHistogramSysts(TH1F* histo, TFile *filesout[(int)MAXSYSTS], bool useOnlyNominal=false){  
  for(int sy=0;sy<(int)MAXSYSTS;++sy){
    cout << " writing histo "<< histo->GetName()<< " in file "<< filesout[(int)sy]->GetName()<<endl;;
    //TString ns= weightedSystsNames((weightedSysts)sy);
    filesout[(int)sy]->cd();
    histo->Write();
    //    histo[sy]=new TH1F(name+ns,name+ns,nbins,min,max);
  }
}

TH1F * makeproduct(TH2F * h){

  //Make temporary histos to rebin
  //  TH1F *hA = (TH1F*)h_A->Clone("hA");
  // TH1F *hB = (TH1F*)h_B->Clone("hB");

  //  hA->Rebin(rebinA);
  // hB->Rebin(rebinB);
  
  //get nbins from new histos
  int nbinsA = h->GetNbinsX();
  int nbinsB = h->GetNbinsY();
  double min = h->GetXaxis()->GetBinLowEdge(1)*h->GetYaxis()->GetBinLowEdge(1);
  double max = h->GetXaxis()->GetBinLowEdge(nbinsA+1)*h->GetYaxis()->GetBinLowEdge(nbinsB+1);
  //Get the actual name from the original histograms 
  string name = (string)(h->GetName()) + "_1D";
  
  //Initialize histogram 
  TH1F * result = new TH1F(name.c_str(),name.c_str(),nbinsA*nbinsB,min,max);
  //Fill histogram
  for(int i =1; i<= nbinsA;++i){
    for(int j =1; j<= nbinsB;++j){
      double value = h->GetBinContent(i,j);
      int k = ((i-1)*nbinsB)+j;
      result->SetBinContent(k,value);
    }
  }
  //  if( integral <= 0.)integral = hA->Integral()/result->Integral();
  //else integral = integral / result->Integral();
  //result->Scale(integral);
  return result;

}
