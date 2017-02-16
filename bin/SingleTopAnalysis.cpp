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
#include "Analysis/NAAnaFW/bin/Weights.h"
//#include "Tprime/TprimeAnalysis/interface/Weights.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "Analysis/NAAnaFW/src/DMTopVariables.h"
#include "SystematicsUtilities.h"

using namespace std;

typedef vector<double> vdouble;
typedef vector<float> vfloat;
typedef vector<int> vint;
typedef vector<bool> vbool;
typedef vector<string> vstring;



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

    string treesDir(argv[8]);
    if(treesDir != "noTrees")std::cout<<"\t treesDir: "<<treesDir<<endl;
    
    TString path_ = path ; 
    std::cout<<"\t File to open: "<<path_<<endl;
    std::cout<<"--------------------"<<endl; 
    
    std::cout << "Info: Loading file collection from " << path << std::endl;
    TFileCollection fc(sample.c_str(),sample.c_str(),path.c_str());
    std::cout << "Info: Files found : " << fc.GetNFiles() << std::endl;

    bool useCSVMSelection = true;
    bool useCSVTSelection = false;
    bool useCMVATSelection = false;
    //  cout << "cat is "<<endl;
    if(cat=="CSVT"){
        useCSVMSelection = false;
        useCSVTSelection = true;
        useCMVATSelection = false;
        }

    if(cat=="CMVAT"){
        useCSVMSelection = false;
        useCSVTSelection = false;
        useCMVATSelection = true;
        }
    
    string btagalgo="CSVM";
    
    if(useCSVMSelection) btagalgo= "CSVM";
    if(useCSVTSelection) btagalgo= "CSVT";
    if(useCMVATSelection) btagalgo= "CMVAT";

    string reportName = "SelectedEvents_"+channel+"_"+cat+"_"+sample+".txt";
    string reportName_step0 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step0.txt";
    string reportName_step1 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step1.txt";
    string reportName_step2 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step2.txt";
    string reportName_step3 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step3.txt";
    string reportName_step4 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step4.txt";
    string reportName_step5 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step5.txt";
    string reportName_step6 = "SelectedEvents_"+channel+"_"+cat+"_"+sample+"_step6.txt";
    //  string reportName = "SelectedEvents_"+channel+""cat"_"+sample+".txt";
    ofstream fileout;
    ofstream fileout_step0;
    ofstream fileout_step1;
    ofstream fileout_step2;
    ofstream fileout_step3;
    ofstream fileout_step4;
    ofstream fileout_step5;
    ofstream fileout_step6;

    bool doSynch=sync.find("sync")!=std::string::npos;
    if(doSynch){
        fileout.open(reportName.c_str(),ios::in | ios::out | ios::trunc);
        fileout<<"RunNumber EvtNumber Lumi "<<std::endl;
        fileout_step0.open(reportName_step0.c_str(),ios::in | ios::out | ios::trunc);
        fileout_step1.open(reportName_step1.c_str(),ios::in | ios::out | ios::trunc);
        //fileout_step1<<"RunNumber | EvtNumber | Lumi | #TightMu | #LooseMu | #LooseEl | #Jets | #BJets | TightMu Pt | TightMuRelIso | AddLoose.Mu Pt | Add.LooseEl Pt | 1stJetPt | 2ndJetPt | 1stB Dis.| 2ndB Dis.| MET | MTW | "<<std::endl;
        fileout_step2.open(reportName_step2.c_str(),ios::in | ios::out | ios::trunc);
        fileout_step3.open(reportName_step3.c_str(),ios::in | ios::out | ios::trunc);
        fileout_step4.open(reportName_step4.c_str(),ios::in | ios::out | ios::trunc);
        fileout_step5.open(reportName_step5.c_str(),ios::in | ios::out | ios::trunc);
        fileout_step6.open(reportName_step6.c_str(),ios::in | ios::out | ios::trunc);
        }
    
    TString outfile = "test/"+sample + "_"+channel+".root";
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
    
    bool addTrees=false;
    if (treesDir!="noTrees") addTrees=true;
    
    TFile *outTreeFile = new TFile;
    if(addTrees){
      outTreeFile = TFile::Open((treesDir+"/trees_"+sample+"_"+channel+".root").c_str(), "RECREATE");
    }
    
    
    Int_t nEvents = (Int_t)chain.GetEntries();
    std::cout<<"Info: Number of Events: "<<nEvents<< endl;
    //nEvents = std::min(nEvents, 1000);
    
    TString weightedSystsNames (weightedSysts sy);
    
    systWeights systZero,syst0BM,syst1BM,syst2BM; 
    int maxSysts=0; 
    int sizeMax=50;
    int muSize, elSize, jetSize,jet20Size, muLooseSize,elLooseSize;
    //float passTrigHT(0.), Ht(0.);
    float Ht(0.), mt(0.);
    float runNumber(0.), lumiSec(0.);
    double evtNumber(0.);
    float jetE[sizeMax], jetPt[sizeMax], jetPhi[sizeMax], jetEta[sizeMax];
    float jetReshapeFactorCSV[sizeMax],jetReshapeFactorCSV_SD[sizeMax];
    float jetReshapeFactorCMVA[sizeMax],jetReshapeFactorCMVA_SD[sizeMax];
    float jetHadronFlavour[sizeMax],jetPartonFlavour[sizeMax];

    float jet20E[sizeMax], jet20Pt[sizeMax], jet20Phi[sizeMax], jet20Eta[sizeMax];
    float jetIsCSVL[sizeMax], jetIsCSVM[sizeMax], jetIsCSVT[sizeMax],jetIsLoose[sizeMax],jetIsTight[sizeMax],jetak4chs_csvv2[sizeMax];
    float jetIsCMVAL[sizeMax], jetIsCMVAM[sizeMax], jetIsCMVAT[sizeMax],jetak4chs_cmvav2[sizeMax];
    float jetPassesB[sizeMax], jet20PassesB[sizeMax];
    float jet20IsCSVL[sizeMax], jet20IsCSVM[sizeMax], jet20IsCSVT[sizeMax],jet20IsLoose[sizeMax],jet20IsTight[sizeMax],jet20ak4chs_csvv2[sizeMax];
    float jet20ReshapeFactorCSV[sizeMax],jet20ReshapeFactorCSV_SD[sizeMax];
    float jet20ReshapeFactorCMVA[sizeMax],jet20ReshapeFactorCMVA_SD[sizeMax];
    float jet20HadronFlavour[sizeMax],jet20PartonFlavour[sizeMax];



    float nTightMuons, nTightElectrons, nVetoElectrons, nLooseMuons, nJets,nCSVJets;//, nCSVLJets;
    float nGoodPV, nPV, numTrueInt, w_pu;
    float metPt[1],metPhi[1],metPx[1],metPy[1];
    
    //double met,metpx,metpy;
    float met,metpx,metpy;
    float muE[sizeMax], muPt[sizeMax], muEta[sizeMax], muIso[sizeMax], muIsTight[sizeMax],muIsLoose[sizeMax], muPhi[sizeMax],elPt[sizeMax],muLoosePhi[sizeMax],muLooseEta[sizeMax];
    float muAntiIsoE[sizeMax], muAntiIsoPt[sizeMax], muAntiIsoEta[sizeMax], muAntiIsoIso[sizeMax], muAntiIsoIsTight[sizeMax],muAntiIsoIsLoose[sizeMax], muAntiIsoPhi[sizeMax];
    int   muAntiIsoSize;
    float muCharge[sizeMax], elCharge[sizeMax];
    float muAntiIsoCharge[sizeMax];
    float topcharge,topsize,antitopsize;

    float elEta[sizeMax], elIso[sizeMax], elIsTight[sizeMax], elIsVeto[sizeMax], elE[sizeMax], elPhi[sizeMax], elPassesDRmu[sizeMax];    
    float lep1Pt(0), lep1Eta(0), lep1Phi(0), lep1E(0), lep2Pt(0), lep2Eta(0), lep2Phi(0),  lep2E(0), lep1Flavour(0), lep1Charge(0), lep2Flavour(0), lep2Charge(0);
    int nPDF=102;

    float slTrigIsoMu20_v1(0.), slTrigIsoMu20_v2(0.), slTrigIsoMu20_v3(0.);
    float slTrigIsoMu22_v1(0.), slTrigIsoMu22_v2(0.), slTrigIsoMu22_v3(0.);
    float slTrigIsoMu24_v1(0.), slTrigIsoMu24_v2(0.), slTrigIsoMu24_v3(0.),slTrigIsoMu24_v4(0.);

    float slTrigIsoTkMu20_v1(0.), slTrigIsoTkMu20_v2(0.), slTrigIsoTkMu20_v3(0.);
    float slTrigIsoTkMu22_v1(0.), slTrigIsoTkMu22_v2(0.), slTrigIsoTkMu22_v3(0.);
    float slTrigIsoTkMu24_v1(0.), slTrigIsoTkMu24_v2(0.), slTrigIsoTkMu24_v3(0.),slTrigIsoTkMu24_v4(0.);

    float slTrigEle_v1(0.), slTrigEle_v2(0.);
    
    bool  TrigIsoMu20 = false;
    bool  TrigIsoMu22 = false;
    bool  TrigIsoMu24 = false;
    
    float n_trig(0.),n_lepton(0.),n_loose_veto(0.),n_lepton_cross_veto(0),n_2j(0), n_2j1t(0), n_2j1tmtw(0.);
    int nev_trig(0),nev_lepton(0),nev_loose_veto(0),nev_lepton_cross_veto(0),nev_2j(0),nev_2j1t(0), nev_2j1tmtw(0.);

    float LHEWeightSign[1] = {1.};
    float w(1.),w_q2up(1.),w_q2down(1.),w_zero(1.),w_top(1.);
    float w_pdfs[nPDF];
    float passElTrig(0.), passMuTrig(0.);
    float passHadTrig(0.) ;

    float bWeight0CSVM(1.), bWeight0CSVMBTagUp(1.),bWeight0CSVMBTagDown(1.),bWeight0CSVMMisTagUp(1.),bWeight0CSVMMisTagDown(1.);
    float bWeight1CSVM(1.), bWeight1CSVMBTagUp(1.),bWeight1CSVMBTagDown(1.),bWeight1CSVMMisTagUp(1.),bWeight1CSVMMisTagDown(1.);
    float bWeight2CSVM(1.), bWeight2CSVMBTagUp(1.),bWeight2CSVMBTagDown(1.),bWeight2CSVMMisTagUp(1.),bWeight2CSVMMisTagDown(1.);
    

    float lepWeight1Mu(1.), lepWeight1MuLepUp(1.), lepWeight1MuLepDown(1.);
    float lepWeight1Ele(1.), lepWeight1EleLepUp(1.), lepWeight1EleLepDown(1.);
    
    //  float bWeight0CSVT, bWeight0CSVTBTagUp,bWeight0CSVTBTagDown,bWeight0CSVTMisTagUp,bWeight0CSVTMisTagDown;
    //  float bWeight1CSVT, bWeight1CSVTBTagUp,bWeight1CSVTBTagDown,bWeight1CSVTMisTagUp,bWeight1CSVTMisTagDown;
    //  float bWeight2CSVT, bWeight2CSVTBTagUp,bWeight2CSVTBTagDown,bWeight2CSVTMisTagUp,bWeight2CSVTMisTagDown;
    
    bool doTtoSDDecay=false;
    bool reweightall=false;

    bool onlyNominal=false;
    systZero.setOnlyNominal(onlyNominal);

    bool addPDF=false,addQ2=false,addTopPt=false,addVHF=false,addTTSplit=false;
    if(sample=="TT"){
        addTTSplit=false;
        }

    if(sample.find("_sd")!=std::string::npos){
        doTtoSDDecay=true;
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

    if(isData=="MC"){
        chain.SetBranchAddress("Event_LHEWeight0", &w_zero); 
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
            //cout << "Event_LHEWeight"+pstr<<w_pdfs<<endl;
            }
        chain.SetBranchAddress("Event_T_Weight",&w_top);
      //    chain.SetBranchAddress("Event_T_Weight",&w_top);
        }
    
    chain.SetBranchAddress("Event_T_size",&topsize);
    chain.SetBranchAddress("Event_Tbar_size",&antitopsize);

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
    chain.SetBranchAddress("electronsVeto_size", &elLooseSize);
    
    chain.SetBranchAddress("muonsTight_E", muE);
    chain.SetBranchAddress("muonsTight_Phi", muPhi);
    chain.SetBranchAddress("muonsTight_Eta", muEta);
    chain.SetBranchAddress("muonsTight_Pt", muPt);
    chain.SetBranchAddress("muonsTight_Charge", muCharge);
    chain.SetBranchAddress("muonsTight_Iso04", muIso);
    chain.SetBranchAddress("muonsTight_IsTightMuon", muIsTight);
    chain.SetBranchAddress("muonsTight_IsLooseMuon", muIsLoose);
    chain.SetBranchAddress("muonsTight_size", &muSize);
    chain.SetBranchAddress("muonsLoose_size", &muLooseSize);
    chain.SetBranchAddress("muonsLoose_Phi", muLoosePhi);
    chain.SetBranchAddress("muonsLoose_Eta", muLooseEta);
    
    chain.SetBranchAddress("muonsTightAntiIso_E", muAntiIsoE);
    chain.SetBranchAddress("muonsTightAntiIso_Phi", muAntiIsoPhi);
    chain.SetBranchAddress("muonsTightAntiIso_Eta", muAntiIsoEta);
    chain.SetBranchAddress("muonsTightAntiIso_Pt", muAntiIsoPt);
    chain.SetBranchAddress("muonsTightAntiIso_Charge", muAntiIsoCharge);
    chain.SetBranchAddress("muonsTightAntiIso_Iso04", muAntiIsoIso);
    chain.SetBranchAddress("muonsTightAntiIso_IsTightMuon", muAntiIsoIsTight);
    chain.SetBranchAddress("muonsTightAntiIso_IsLooseMuon", muAntiIsoIsLoose);
    chain.SetBranchAddress("muonsTightAntiIso_size", &muAntiIsoSize);


    chain.SetBranchAddress("jetsAK4CHSTight_CorrE",      &jetE);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt",     &jetPt);
    chain.SetBranchAddress("jetsAK4CHSTight_Phi",    &jetPhi);
    chain.SetBranchAddress("jetsAK4CHSTight_Eta",    &jetEta);
    chain.SetBranchAddress("jetsAK4CHSTight_size",   &jetSize);
    chain.SetBranchAddress("jetsAK4CHSTight_IsCSVL", &jetIsCSVL);
    chain.SetBranchAddress("jetsAK4CHSTight_IsCSVM", &jetIsCSVM);
    chain.SetBranchAddress("jetsAK4CHSTight_IsCSVT", &jetIsCSVT);
    chain.SetBranchAddress("jetsAK4CHSTight_IsLoose",&jetIsLoose);
    chain.SetBranchAddress("jetsAK4CHSTight_IsTight",&jetIsTight);
    chain.SetBranchAddress("jetsAK4CHSTight_CSVv2",  &jetak4chs_csvv2);
    chain.SetBranchAddress("jetsAK4CHSTight_reshapeFactorCSV",  &jetReshapeFactorCSV);
    chain.SetBranchAddress("jetsAK4CHSTight_reshapeFactorCSV_SD",  &jetReshapeFactorCSV_SD);
    chain.SetBranchAddress("jetsAK4CHSTight_reshapeFactorCMVA",  &jetReshapeFactorCMVA);
    chain.SetBranchAddress("jetsAK4CHSTight_reshapeFactorCMVA_SD",  &jetReshapeFactorCMVA_SD);
    
    chain.SetBranchAddress("jetsAK4CHSTight_HadronFlavour", &jetHadronFlavour);
    chain.SetBranchAddress("jetsAK4CHSTight_PartonFlavour",  &jetPartonFlavour);

    //  float jetReshapeFactorCSV[sizeMax],jetReshapeFactorCSV_SD[sizeMax];
     chain.SetBranchAddress("jetsAK4CHSTight_IsCMVAL", &jetIsCMVAL);
     chain.SetBranchAddress("jetsAK4CHSTight_IsCMVAM", &jetIsCMVAM);
     chain.SetBranchAddress("jetsAK4CHSTight_IsCMVAT", &jetIsCMVAT);
     chain.SetBranchAddress("jetsAK4CHSTight_CMVAv2",  &jetak4chs_cmvav2);

     chain.SetBranchAddress(("jetsAK4CHSTight_Is"+btagalgo).c_str(),  &jetPassesB);
     chain.SetBranchAddress(("jetsAK4CHSTight_CorrPt_20_Is"+btagalgo).c_str(),  &jet20PassesB);
        

    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_CorrE",      &jet20E);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_CorrPt",     &jet20Pt);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_Phi",    &jet20Phi);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_Eta",    &jet20Eta);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_size",   &jet20Size);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_IsCSVL", &jet20IsCSVL);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_IsCSVM", &jet20IsCSVM);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_IsCSVT", &jet20IsCSVT);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_IsLoose",&jet20IsLoose);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_IsTight",&jet20IsTight);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_CSVv2",  &jet20ak4chs_csvv2);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_reshapeFactorCSV",  &jet20ReshapeFactorCSV);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_reshapeFactorCSV_SD",  &jet20ReshapeFactorCSV_SD);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_reshapeFactorCSV",  &jet20ReshapeFactorCMVA);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_reshapeFactorCSV_SD",  &jet20ReshapeFactorCMVA_SD);

    

    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_HadronFlavour", &jet20HadronFlavour);
    chain.SetBranchAddress("jetsAK4CHSTight_CorrPt_20_PartonFlavour",  &jet20PartonFlavour);  


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
    

    if(isData!="DATA"){
      
      
     
      chain.SetBranchAddress(("Event_bWeight2"+btagalgo+"Tight").c_str(), &bWeight2CSVM);
      chain.SetBranchAddress(("Event_bWeightMisTagUp2"+btagalgo+"Tight").c_str(), &bWeight2CSVMMisTagUp);
      chain.SetBranchAddress(("Event_bWeightMisTagDown2"+btagalgo+"Tight").c_str(), &bWeight2CSVMMisTagDown);
      chain.SetBranchAddress(("Event_bWeightBTagUp2"+btagalgo+"Tight").c_str(), &bWeight2CSVMBTagUp);
      chain.SetBranchAddress(("Event_bWeightBTagDown2"+btagalgo+"Tight").c_str(), &bWeight2CSVMBTagDown);
      
      chain.SetBranchAddress(("Event_bWeight0"+btagalgo+"Tight").c_str(), &bWeight0CSVM);
      
      chain.SetBranchAddress(("Event_bWeightMisTagUp0"+btagalgo+"Tight").c_str(), &bWeight0CSVMMisTagUp);
      chain.SetBranchAddress(("Event_bWeightMisTagDown0"+btagalgo+"Tight").c_str(), &bWeight0CSVMMisTagDown);
      chain.SetBranchAddress(("Event_bWeightBTagUp0"+btagalgo+"Tight").c_str(), &bWeight0CSVMBTagUp);
      chain.SetBranchAddress(("Event_bWeightBTagDown0"+btagalgo+"Tight").c_str(), &bWeight0CSVMBTagDown);
      
      chain.SetBranchAddress(("Event_bWeight1"+btagalgo+"Tight").c_str(), &bWeight1CSVM);
      chain.SetBranchAddress(("Event_bWeightMisTagUp1"+btagalgo+"Tight").c_str(), &bWeight1CSVMMisTagUp);
      chain.SetBranchAddress(("Event_bWeightMisTagDown1"+btagalgo+"Tight").c_str(), &bWeight1CSVMMisTagDown);
      chain.SetBranchAddress(("Event_bWeightBTagUp1"+btagalgo+"Tight").c_str(), &bWeight1CSVMBTagUp);
      chain.SetBranchAddress(("Event_bWeightBTagDown1"+btagalgo+"Tight").c_str(), &bWeight1CSVMBTagDown);
      
    }

    if(isData=="MC"){
        chain.SetBranchAddress("metFull_Pt",metPt);
        chain.SetBranchAddress("metFull_CorrPhi",metPhi);}
    else {
        chain.SetBranchAddress("metFull_Pt",metPt);
        chain.SetBranchAddress("metFull_CorrPhi",metPhi);
        }

    chain.SetBranchAddress("metFull_Px",metPx);
    chain.SetBranchAddress("metFull_Py",metPy);
    
    //Muon trigger
    //string ptmu="20";
    //if(runNumber>=272023 && runNumber<=274443) ptmu="20";
    //if(runNumber>=274954 && runNumber<=276811) ptmu="22";
    //if(runNumber>=276824 && runNumber<=99999999) ptmu="24";
    //chain.SetBranchAddress(("Event_passesHLT_IsoMu"+ptmu+"_v1").c_str(), &slTrigIsoMu22_v1);
    //chain.SetBranchAddress(("Event_passesHLT_IsoMu"+ptmu+"_v2").c_str(), &slTrigIsoMu22_v2);
    //chain.SetBranchAddress(("Event_passesHLT_IsoMu"+ptmu+"_v3").c_str(), &slTrigIsoMu22_v3);
    //Event_passesHLT_IsoTkMu22_v3
    
    chain.SetBranchAddress("Event_passesHLT_IsoMu20_v1", &slTrigIsoMu20_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoMu20_v2", &slTrigIsoMu20_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoMu20_v3", &slTrigIsoMu20_v3);
    
    chain.SetBranchAddress("Event_passesHLT_IsoMu22_v1", &slTrigIsoMu22_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoMu22_v2", &slTrigIsoMu22_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoMu22_v3", &slTrigIsoMu22_v3);
    
    chain.SetBranchAddress("Event_passesHLT_IsoMu24_v1", &slTrigIsoMu24_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoMu24_v2", &slTrigIsoMu24_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoMu24_v3", &slTrigIsoMu24_v3);
    //
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu20_v1", &slTrigIsoTkMu20_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu20_v2", &slTrigIsoTkMu20_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu20_v3", &slTrigIsoTkMu20_v3);
    
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu22_v1", &slTrigIsoTkMu22_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu22_v2", &slTrigIsoTkMu22_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu22_v3", &slTrigIsoTkMu22_v3);
    
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu24_v1", &slTrigIsoTkMu24_v1);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu24_v2", &slTrigIsoTkMu24_v2);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu24_v3", &slTrigIsoTkMu24_v3);
    
    chain.SetBranchAddress("Event_passesHLT_IsoMu24_v4", &slTrigIsoMu24_v4);
    chain.SetBranchAddress("Event_passesHLT_IsoTkMu24_v4", &slTrigIsoTkMu24_v4);
        

  
    //Electron Triger
    if(isData=="MC"){
        chain.SetBranchAddress("Event_passesHLT_Ele27_eta2p1_WP75_Gsf_v1", &slTrigEle_v1);
        chain.SetBranchAddress("Event_passesHLT_Ele27_eta2p1_WP75_Gsf_v2", &slTrigEle_v2);
        }
    else{
        chain.SetBranchAddress("Event_passesHLT_Ele27_eta2p1_WPLoose_Gsf_v1", &slTrigEle_v1);
        chain.SetBranchAddress("Event_passesHLT_Ele27_eta2p1_WPLoose_Gsf_v2", &slTrigEle_v2);
        }
    chain.SetBranchAddress("Event_passesSingleElTriggers", &passElTrig);
    chain.SetBranchAddress("Event_passesSingleMuTriggers", &passMuTrig);
    chain.SetBranchAddress("Event_passesHadronicTriggers", &passHadTrig);
  


    /********************************************************************************/
    /**************                    MiniTree branching              ***************/
    /********************************************************************************/
    if(addTrees){
    outTreeFile->cd();
    }
    TTree * trees1T[6];
    TTree * trees2T[3];
    //    float * flo = new float(1.0); 
    //    float * mtw_2j1t = new float(-999.0); 
    
    //
    //Standard selection: 2j1t + mtw cut
    float * w_2j1t_mtwcut = new float(-999.0); 
    float * etajprime_2j1t_mtwcut = new float(-999.0); 
    float * topMass_2j1t_mtwcut = new float(-999.0); 
    float * nextrajets_2j1t_mtwcut = new float(-999.0); 
    float * leadingextrajetpt_2j1t_mtwcut = new float(-999.0); 
    float * leadingextrajetcsv_2j1t_mtwcut = new float(-999.0); 
    float * leadingextrajetcsvweight_2j1t_mtwcut = new float(-999.0); 
    float * leadingextrajetcsvweight_sd_2j1t_mtwcut = new float(-999.0); 
    float * topMassExtra_2j1t_mtwcut = new float(-999.0); 
    //Tight t-channel selection: 2j1t + mtw cut + etajprime > 2.5
    float * w_2j1t_mtwcut_sr = new float(-999.0); 
    float * topMass_2j1t_mtwcut_sr = new float(-999.0); 
    float * nextrajets_2j1t_mtwcut_sr = new float(-999.0); 
    float * leadingextrajetpt_2j1t_mtwcut_sr = new float(-999.0); 
    float * leadingextrajetcsv_2j1t_mtwcut_sr = new float(-999.0); 
    float * leadingextrajetcsvweight_2j1t_mtwcut_sr = new float(-999.0); 
    float * leadingextrajetcsvweight_sd_2j1t_mtwcut_sr = new float(-999.0); 
    float * topMassExtra_2j1t_mtwcut_sr = new float(-999.0); 
    //3j2t selection
    float * w_3j2t = new float(-999.0); 
    float * mtw_3j2t = new float(-999.0); 
    float * etajprime_3j2t = new float(-999.0); 
    float * topMassLeading_3j2t = new float(-999.0); 
    float * topMassSecond_3j2t = new float(-999.0); 
    float * nextrajets_3j2t = new float(-999.0); 

    syst0BM.copySysts(systZero);
    syst1BM.copySysts(systZero);
    syst2BM.copySysts(systZero);
    if(addTrees){
      cout << "adding trees "<<endl;
      //syst1BM.addSelection("3j1t_mtwcut");
      syst1BM.addSelection("2j1t_mtwcut");
      syst1BM.addSelection("2j1t_mtwcut_sr");
      syst1BM.initTreesSysts(trees1T, outTreeFile);

      syst2BM.addSelection("3j2t");
      syst2BM.initTreesSysts(trees2T, outTreeFile);
      cout << " bef branch trees" <<endl;

      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut","w", outTreeFile, w_2j1t_mtwcut);

      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut","etajprime", outTreeFile, etajprime_2j1t_mtwcut);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut","topMass", outTreeFile, topMass_2j1t_mtwcut);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut","topMassExtra", outTreeFile, topMassExtra_2j1t_mtwcut);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut","nextrajets", outTreeFile, nextrajets_2j1t_mtwcut);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut","leadingextrajetcsv", outTreeFile, leadingextrajetcsv_2j1t_mtwcut);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut","leadingextrajetpt", outTreeFile, leadingextrajetpt_2j1t_mtwcut);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut","leadingextrajetcsvweight", outTreeFile, leadingextrajetcsvweight_2j1t_mtwcut);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut","leadingextrajetcsvweight_sd", outTreeFile, leadingextrajetcsvweight_sd_2j1t_mtwcut);
      
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut_sr","w", outTreeFile, w_2j1t_mtwcut_sr);

      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut_sr","topMass", outTreeFile, topMass_2j1t_mtwcut_sr);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut_sr","topMassExtra", outTreeFile, topMassExtra_2j1t_mtwcut_sr);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut_sr","nextrajets", outTreeFile, nextrajets_2j1t_mtwcut_sr);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut_sr","leadingexrajetcsv", outTreeFile, leadingextrajetcsv_2j1t_mtwcut_sr);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut_sr","leadingextrajetpt", outTreeFile, leadingextrajetpt_2j1t_mtwcut_sr);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut","leadingextrajetcsvweight", outTreeFile, leadingextrajetcsvweight_2j1t_mtwcut_sr);
      syst1BM.branchTreesSysts(trees1T,"2j1t_mtwcut","leadingextrajetcsvweight_sd", outTreeFile, leadingextrajetcsvweight_sd_2j1t_mtwcut_sr);

      
      syst2BM.branchTreesSysts(trees2T,"3j2t","w", outTreeFile, w_3j2t);

      syst2BM.branchTreesSysts(trees2T,"3j2t","mtw", outTreeFile, mtw_3j2t);
      syst2BM.branchTreesSysts(trees2T,"3j2t","etajprime", outTreeFile, etajprime_3j2t);
      syst2BM.branchTreesSysts(trees2T,"3j2t","topMassLeading", outTreeFile, topMassLeading_3j2t);
      syst2BM.branchTreesSysts(trees2T,"3j2t","topMassSecond", outTreeFile, topMassSecond_3j2t);
    }

    
    /********************************************************************************/
    /**************                    Histogram booking              ***************/
    /********************************************************************************/
    TH1F *h_cutFlow = new TH1F("h_cutFlow","cutflow",10,-0.5,9.5);
    TH1F *h_weight_sign = new TH1F("h_weight_sign","weight correction factor before all selections, gives the effective number of events",2,-2.0,2.0);
    chainNEvents.Project("h_weight_sign","Event_LHEWeight10/abs(Event_LHEWeight10)"); // I messed up and weight zero is not  saved in the full chain... but this one should do the trick as well! Just do h_weight_sign->GetMean() and multiply this to get the correct number :)
    double zeroWeight=h_weight_sign->GetMean(); // Will have to add a zeroWeight vector in the systematics for a proper treatment of the acceptance in systematics. For now we do it via hardcoding outside.
    cout << " mean  "<< zeroWeight<<endl;
    
  
    TH1F *h_weightZero[maxSysts]; systZero.initHistogramsSysts(h_weightZero,"h_weightZero","weight before all selections, normalized by the total number of events",2000,0,10.0);
    TH1F *h_nPV[maxSysts];        systZero.initHistogramsSysts(h_nPV,"h_nPV","nPV",102,0,51);
    TH1F *h_nGoodPV[maxSysts];    systZero.initHistogramsSysts(h_nGoodPV,"h_nGoodPV","nGoodPV",92,0,46);
    TH1F *h_nTruePV[maxSysts];    systZero.initHistogramsSysts(h_nTruePV,"h_nTruePV","nTruePV",80,0,40);
    TH1F *h_nJets[maxSysts];      systZero.initHistogramsSysts(h_nJets,"h_nJets","Number of tight jets",13,-0.5,12.5);
    TH1F *h_nbJets[maxSysts];     systZero.initHistogramsSysts(h_nbJets,"h_nbJets","Number of tight b-jets",11,-0.5,10.5);

    //2j0t 
    TH1F *h_2j0t_dR_lepjetpt40_1st[maxSysts];      systZero.initHistogramsSysts(h_2j0t_dR_lepjetpt40_1st,   "h_2j0t_dR_lepjetpt40_1st",  "dR lep-jet ",50,0,10);
    TH1F *h_2j0t_dPhi_lepjetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j0t_dPhi_lepjetpt40_1st,   "h_2j0t_dPhi_lepjetpt40_1st",  "dPhi lep-jet ",50,-3.2,3.2);
    TH1F *h_2j0t_dEta_lepjetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j0t_dEta_lepjetpt40_1st,   "h_2j0t_dEta_lepjetpt40_1st",  "dPhi lep-jet ",50,-3.2,3.2);
    TH1F *h_2j0t_MuEta[maxSysts];          systZero.initHistogramsSysts(h_2j0t_MuEta,       "h_2j0t_MuEta",      "muon eta ",100,-2.1,2.1); 
    TH1F *h_2j0t_MuPhi[maxSysts];          systZero.initHistogramsSysts(h_2j0t_MuPhi,       "h_2j0t_MuPhi",      "muon phi ",100,-3.2,3.2); 
    TH1F *h_2j0t_jetpt40_1st[maxSysts];    systZero.initHistogramsSysts(h_2j0t_jetpt40_1st,   "h_2j0t_jetpt40_leading","2j0t leading jet Pt ",250,0,500);
    TH1F *h_2j0t_jetpt40_2nd[maxSysts];    systZero.initHistogramsSysts(h_2j0t_jetpt40_2nd,   "h_2j0t_jetpt40_subleading","2j0t Sub. leading jet Pt ",250,0,500);
    TH1F *h_2j0t_mtw[maxSysts];            systZero.initHistogramsSysts(h_2j0t_mtw,           "h_2j0t_mtw",     "2j0t mtw ",250,0,500);

    //After mtw>50 cut
    TH1F *h_2j0t_mtwcut_jetpt40_1st[maxSysts];   systZero.initHistogramsSysts(h_2j0t_mtwcut_jetpt40_1st,   "h_2j0t_mtwcut_jetpt40_leading","2j0t leading jet Pt ",250,0,500);
    TH1F *h_2j0t_mtwcut_jetpt40_2nd[maxSysts];   systZero.initHistogramsSysts(h_2j0t_mtwcut_jetpt40_2nd,   "h_2j0t_mtwcut_jetpt40_subleading","2j0t Sub. leading jet Pt ",250,0,500);
    TH1F *h_2j0t_mtwcut_jeteta40_1st[maxSysts];  systZero.initHistogramsSysts(h_2j0t_mtwcut_jeteta40_1st,   "h_2j0t_mtwcut_jeteta40_leading","2j0t leading jet eta ",200,-4.7,4.7);
    TH1F *h_2j0t_mtwcut_incljeteta[maxSysts];    systZero.initHistogramsSysts(h_2j0t_mtwcut_incljeteta,   "h_2j0t_mtwcut_incljeteta","2j0t leading inclusive jet eta ",200,-4.7,4.7);

    TH1F *h_2j0t_mtwcut_jetabseta40_1st[maxSysts]; systZero.initHistogramsSysts(h_2j0t_mtwcut_jetabseta40_1st,   "h_2j0t_mtwcut_jetabseta40_leading","2j0t leading jet eta ",200,0,4.7);
    TH1F *h_2j0t_mtwcut_incljetabseta[maxSysts];   systZero.initHistogramsSysts(h_2j0t_mtwcut_incljetabseta,   "h_2j0t_mtwcut_incljetabseta","2j0t leading inclusive jet eta ",200,0,4.7);
    TH1F *h_2j0t_incljetabseta[maxSysts];          systZero.initHistogramsSysts(h_2j0t_incljetabseta,   "h_2j0t_incljetabseta","2j0t leading inclusive jet eta ",200,0,4.7);

  
  //2j1t Category
  TH1F *h_2j1t_bjetpt[maxSysts];    systZero.initHistogramsSysts(h_2j1t_bjetpt,  "h_2j1t_bjetpt",  "2j1t b jet pt",250,0,500);
  TH1F *h_2j1t_ljetpt[maxSysts];    systZero.initHistogramsSysts(h_2j1t_ljetpt,  "h_2j1t_ljetpt",   "2j1t light jet pt" ,250,0,500);
  TH1F *h_2j1t_ljeteta[maxSysts];   systZero.initHistogramsSysts(h_2j1t_ljeteta, "h_2j1t_ljeteta",  "2j1t light jet eta",200,-4.7,4.7);
  TH1F *h_2j1t_nMu[maxSysts];       systZero.initHistogramsSysts(h_2j1t_nMu,     "h_2j1t_nMu",     "2j1t Number of tight Muons",13,-0.5,12.5);
  TH1F *h_2j1t_MuPt[maxSysts];      systZero.initHistogramsSysts(h_2j1t_MuPt,    "h_2j1t_MuPt",    "2j1t muon pt ",250,0,500);
  TH1F *h_2j1t_MuEta[maxSysts];     systZero.initHistogramsSysts(h_2j1t_MuEta,   "h_2j1t_MuEta",   "2j1t muon eta ",100,-2.1,2.1);
  TH1F *h_2j1t_MuPhi[maxSysts];     systZero.initHistogramsSysts(h_2j1t_MuPhi,   "h_2j1t_MuPhi",   "2j1t muon phi ",100, -3.2, 3.2);
  TH1F *h_2j1t_MuE[maxSysts];       systZero.initHistogramsSysts(h_2j1t_MuE,     "h_2j1t_MuE",     "2j1t muon e ",250,0,500);
  TH1F *h_2j1t_MuCharge[maxSysts];  systZero.initHistogramsSysts(h_2j1t_MuCharge,"h_2j1t_MuCharge","2j1t muon charge ",2,-1,1);
  TH1F *h_2j1t_mtw[maxSysts];       systZero.initHistogramsSysts(h_2j1t_mtw,     "h_2j1t_mtw",     "2j1t mtw ",250,0,500);
  TH1F *h_2j1t_topMass[maxSysts];   systZero.initHistogramsSysts(h_2j1t_topMass, "h_2j1t_topMass", "2j1t top mass",200,100,500);
  TH1F *h_2j1t_topPt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_topPt,   "h_2j1t_topPt",   "2j1t top pt",250,0,500);
  TH1F *h_2j1t_etajprime[maxSysts];   systZero.initHistogramsSysts(h_2j1t_etajprime,  "h_2j1t_etajprime",  "2j1t light jet abs(eta)",100,0.,4.7);

  //Jets pt 20 spectra.
  TH1F *h_2j1t_nextrajets[maxSysts];   systZero.initHistogramsSysts(h_2j1t_nextrajets,  "h_2j1t_nextrajets",  "2j1t number of jets with 20<pt<40",5,-0.5,4.5);
  TH1F *h_2j1t_leadingextrajetpt[maxSysts];   systZero.initHistogramsSysts(h_2j1t_leadingextrajetpt,  "h_2j1t_leadingextrajetpt",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_leadingextrajetcsv[maxSysts];   systZero.initHistogramsSysts(h_2j1t_leadingextrajetcsv,  "h_2j1t_leadingextrajetcsv",  "2j1t leading extra jet csv value",100,0.,1.0);
  TH1F *h_2j1t_leadingextrajetpt_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_leadingextrajetpt_sd_b,  "h_2j1t_leadingextrajetpt_sd_b",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_leadingextrajetcsv_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_leadingextrajetcsv_sd_b,  "h_2j1t_leadingextrajetcsv_sd_b",  "2j1t leading extra jet csv value for b from top decays",100,0.,1.0);
  
  TH1F *h_2j1t_leadingextrajetcsv_reshape_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_leadingextrajetcsv_reshape_sd_b,  "h_2j1t_leadingextrajetcsv_reshape_sd_b",  "2j1t leading extra jet csv reshaping from b to sd for b from top decays",1000,0.,10.0);

  
  //After mtw>50 cut, label mtwcut
  TH1F *h_2j1t_mtwcut_topMass[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_topMass, "h_2j1t_mtwcut_topMass", "2j1t top mass",200,100,500);
  TH1F *h_2j1t_mtwcut_mtw[maxSysts];systZero.initHistogramsSysts(h_2j1t_mtwcut_mtw,     "h_2j1t_mtwcut_mtw",     "2j1t mtw ",250,0,500);
  TH1F *h_2j1t_mtwcut_bjetpt[maxSysts];    systZero.initHistogramsSysts(h_2j1t_mtwcut_bjetpt,  "h_2j1t_mtwcut_bjetpt",  "2j1t b jet pt",250,0,500);
  TH1F *h_2j1t_mtwcut_ljetpt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_mtwcut_ljetpt,   "h_2j1t_mtwcut_ljetpt",   "2j1t light jet pt" ,250,0,500);
  TH1F *h_2j1t_mtwcut_ljeteta[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_ljeteta,  "h_2j1t_mtwcut_ljeteta",  "2j1t light jet eta",200,-4.7,4.7);
  TH1F *h_2j1t_mtwcut_topPt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_mtwcut_topPt,   "h_2j1t_mtwcut_topPt",   "2j1t top pt",250,0,500);
  TH1F *h_2j1t_mtwcut_etajprime[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_etajprime,  "h_2j1t_mtwcut_etajprime",  "2j1t light jet abs(eta)",100,0.,4.7);
  TH1F *h_2j1t_mtwcut_nextrajets[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_nextrajets,  "h_2j1t_mtwcut_nextrajets",  "2j1t number of jets with 20<pt<40",5,-0.5,4.5);
  TH1F *h_2j1t_mtwcut_leadingextrajetpt[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt,  "h_2j1t_mtwcut_leadingextrajetpt",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_mtwcut_leadingextrajetcsv[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv,  "h_2j1t_mtwcut_leadingextrajetcsv",  "2j1t csv of the leading extra jet",100,0.,1.0);
  TH1F *h_2j1t_mtwcut_topMassExtra[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_topMassExtra, "h_2j1t_mtwcut_topMassExtra", "2j1t top mass with the first extra jet as b",200,100,500);

  TH1F *h_2j1t_mtwcut_leadingextrajetpt_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt_sd_b,  "h_2j1t_mtwcut_leadingextrajetpt_sd_b",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_mtwcut_leadingextrajetcsv_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv_sd_b,  "h_2j1t_mtwcut_leadingextrajetcsv_sd_b",  "2j1t csv of the leading extra jet  for b from top decays",100,0.,1.0);
  TH1F *h_2j1t_mtwcut_topMassExtra_sd_b[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_topMassExtra_sd_b, "h_2j1t_mtwcut_topMassExtra_sd_b", "2j1t top mass with the first extra jet as b, for b-quarks from top",200,100,500);

  TH1F *h_2j1t_mtwcut_leadingextrajetcsv_reshape_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv_reshape_sd_b,  "h_2j1t_mtwcut_leadingextrajetcsv_reshape_sd_b",  "2j1t leading extra jet csv reshaping from b to sd for b from top decays",1000,0.,10.0);

  
  //mtw > 50; |eta_j'|>2.5, label mtwcut_sr
  TH1F *h_2j1t_mtwcut_sr_topMass[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_topMass, "h_2j1t_mtwcut_sr_topMass", "2j1t top mass",200,100,500);
  TH1F *h_2j1t_mtwcut_sr_topPt[maxSysts];     systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_topPt,   "h_2j1t_mtwcut_sr_topPt",   "2j1t top pt",250,0,500);
  TH1F *h_2j1t_mtwcut_sr_nextrajets[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_nextrajets,  "h_2j1t_mtwcut_sr_nextrajets",  "2j1t number of jets with 20<pt<40",5,-0.5,4.5);

  TH1F *h_2j1t_mtwcut_sr_leadingextrajetpt[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt,  "h_2j1t_mtwcut_sr_leadingextrajetpt",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_mtwcut_sr_leadingextrajetcsv[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv,  "h_2j1t_mtwcut_sr_leadingextrajetcsv",  "2j1t csv of the leading extra jet",100,0.,1.0);
  TH1F *h_2j1t_mtwcut_sr_topMassExtra[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra, "h_2j1t_mtwcut_sr_topMassExtra", "2j1t top mass with the first extra jet as b",200,100,500);
  TH1F *h_2j1t_mtwcut_sr_leadingextrajetcsv_reshape_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv_reshape_sd_b,  "h_2j1t_mtwcut_sr_leadingextrajetcsv_reshape_sd_b",  "2j1t leading extra jet csv reshaping from b to sd for b from top decays",1000,0.,10.0);

  
  TH1F *h_2j1t_mtwcut_sr_leadingextrajetpt_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt_sd_b,  "h_2j1t_mtwcut_sr_leadingextrajetpt_sd_b",  "2j1t pt of the leading extra jet",100,0.,200);
  TH1F *h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b[maxSysts];   systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b,  "h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b",  "2j1t csv of the leading extra jet",100,0.,1.0);
  TH1F *h_2j1t_mtwcut_sr_topMassExtra_sd_b[maxSysts];  systZero.initHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra_sd_b, "h_2j1t_mtwcut_sr_topMassExtra_sd_b", "2j1t top mass with the first extra jet as b",200,100,500);
  


  //Other possible cr:
  //Features of top quarks reconstructed with pt 20.
  
  ////
  //Additional t-channel enriched regions from TOP-12-038
  //Top mass SR 130<mtop<220, cut
  ////Top mass SR 130<mtop<220, mtw > 50 cuts
  //Top mass SB (130<mtop || mtop>220)
  ////

  //3j1t 
  TH1F *h_3j1t_bjetpt[maxSysts];   systZero.initHistogramsSysts(h_3j1t_bjetpt,  "h_3j1t_bjetpt",  "3j1t b jet pt ",250,0,500);
  TH1F *h_3j1t_MuPt[maxSysts];     systZero.initHistogramsSysts(h_3j1t_MuPt,    "h_3j1t_MuPt",    "muon pt ",250,0,500);
  TH1F *h_3j1t_MuEta[maxSysts];    systZero.initHistogramsSysts(h_3j1t_MuEta,   "h_3j1t_MuEta",   "muon eta ",100,-2.1,2.1);
  TH1F *h_3j1t_MuPhi[maxSysts];    systZero.initHistogramsSysts(h_3j1t_MuPhi,   "h_3j1t_MuPhi",   "muon phi ",100, -3.2, 3.2);
  TH1F *h_3j1t_MuE[maxSysts];      systZero.initHistogramsSysts(h_3j1t_MuE,     "h_3j1t_MuE",     "muon e ",250,0,500);
  TH1F *h_3j1t_mtw[maxSysts];      systZero.initHistogramsSysts(h_3j1t_mtw,     "h_3j1t_mtw",     "3j1t mtw ",250,0,500);
  
  //3j2t
  TH1F *h_3j2t_bjetpt[maxSysts];   systZero.initHistogramsSysts(h_3j2t_bjetpt,    "h_3j2t_bjetpt",     "3j2t b jet pt ",100,0,500);
  TH1F *h_3j2t_2ndbjetpt[maxSysts];systZero.initHistogramsSysts(h_3j2t_2ndbjetpt, "h_3j2t_2ndbjetpt",  "3j2t sub lead. b jet pt",100,0,500);
  TH1F *h_3j2t_ljetpt[maxSysts];   systZero.initHistogramsSysts(h_3j2t_ljetpt,    "h_3j2t_ljetpt",     "3j2t light jet pt ",100,0,500);
  TH1F *h_3j2t_ljeteta[maxSysts];   systZero.initHistogramsSysts(h_3j2t_ljeteta,    "h_3j2t_ljeteta",     "3j2t light jet eta ",100,0,500);
  TH1F *h_3j2t_etajprime[maxSysts];   systZero.initHistogramsSysts(h_3j2t_etajprime,    "h_3j2t_etajprime",     "3j2t eta j' ",100,0.,4.7);
  TH1F *h_3j2t_nextrajets[maxSysts];   systZero.initHistogramsSysts(h_3j2t_nextrajets,    "h_3j2t_nextrajets",     "3j2t number of extra jets ",5,-0.5,4.5);
  TH1F *h_3j2t_MuPt[maxSysts];     systZero.initHistogramsSysts(h_3j2t_MuPt,      "h_3j2t_MuPt",       "3j2t muon pt ",250,0,500);
  TH1F *h_3j2t_MuEta[maxSysts];    systZero.initHistogramsSysts(h_3j2t_MuEta,     "h_3j2t_MuEta",      "3j2t muon eta ",100,-2.1,2.1);
  TH1F *h_3j2t_MuPhi[maxSysts];    systZero.initHistogramsSysts(h_3j2t_MuPhi,     "h_3j2t_MuPhi",      "3j2t muon phi ",100, -3.2, 3.2);
  TH1F *h_3j2t_MuE[maxSysts];      systZero.initHistogramsSysts(h_3j2t_MuE,       "h_3j2t_MuE",        "3j2t muon e ",250,0,500);
  TH1F *h_3j2t_MuCharge[maxSysts]; systZero.initHistogramsSysts(h_3j2t_MuCharge,  "h_3j2t_MuCharge",   "3j2t muon charge ",2,-1,1);
  TH1F *h_3j2t_mtw[maxSysts];      systZero.initHistogramsSysts(h_3j2t_mtw,       "h_3j2t_mtw",        "3j2t mtw ",250,0,500);
  
  //After mtw>50 cut

  TopUtilities topUtils;
  
  //Pileup Reweighting
  edm::LumiReWeighting LumiWeights_, LumiWeightsUp_, LumiWeightsDown_;

  TFile* file_eff_trig = TFile::Open("data/SingleMuonTrigger_Z_RunBCD_prompt80X_7p65.root");
  TFile* file_sf_id = TFile::Open("data/MuonID_Z_RunBCD_prompt80X_7p65.root");
  TFile* file_sf_iso = TFile::Open("data/MuonIso_Z_RunBCD_prompt80X_7p65.root");
  
  //cout << "before getting weights "<<endl;

  Weights muTightTrigEff(file_eff_trig,"IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093/efficienciesDATA/abseta_pt_DATA");
  Weights muTightIdSF(file_sf_id,"MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio");
  Weights muTightIsoSF(file_sf_iso,"MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio");

  //  cout << "got weights" <<endl;

    if(isData=="MC"){
        //LumiWeights_ = edm::LumiReWeighting("pu/MyDataPileupHistogram.root", "pu/MyDataPileupHistogram.root","pileup","pileup");
        //LumiWeightsUp_ = edm::LumiReWeighting("pu/MyDataPileupHistogram.root", "pu/MyDataPileupHistogram.root","pileup","pileup");
        //LumiWeightsDown_ = edm::LumiReWeighting("pu/MyDataPileupHistogram.root", "pu/MyDataPileupHistogram.root","pileup","pileup");
        LumiWeights_ = edm::LumiReWeighting("pu/MCPU.root", "pu/MyDataPileupHistogram.root","pileup","pileup");
        LumiWeightsUp_ = edm::LumiReWeighting("pu/MCPU.root", "pu/MyDataPileupHistogram_up.root","pileup","pileup");
        LumiWeightsDown_ = edm::LumiReWeighting("pu/MCPU.root", "pu/MyDataPileupHistogram_down.root","pileup","pileup");
        //LumiWeightsUp_ = edm::LumiReWeighting("pu/puMC.root", "pu/MyDataPileupHistogram_up.root","MC_pu","pileup");
        //LumiWeightsDown_ = edm::LumiReWeighting("pu/puMC.root", "pu/MyDataPileupHistogram_down.root","MC_pu","pileup");
        }
  
    for(Int_t evt=0; evt<nEvents; evt++ ){
        if(evt%10000==1 ){
        cout<<"Info: Running on event: "<<evt<<endl; 
        }
    chain.GetEntry(evt);
    int maxJetLoop = min(15, jetSize);
    int maxJet20Loop = min(15, jet20Size);

    //  cout << " jet 20 loop size"<< maxJet20Loop <<endl;
    int maxMuLoop = min(6, muSize);
    if(channel == "muonantiiso") maxMuLoop = min(6, muAntiIsoSize); 
    int maxElLoop = min(6, elSize);
  
    //step 1 Trigger
    if(isData=="DATA"){
        TrigIsoMu20= false;
        //TrigIsoMu20= (runNumber>=272023 && runNumber<=274443 && (slTrigIsoMu20_v1 || slTrigIsoMu20_v2 || slTrigIsoMu20_v3));
        //TrigIsoMu22= (runNumber>=274954 && runNumber<=276811 && (slTrigIsoMu22_v1 || slTrigIsoMu22_v2 || slTrigIsoMu22_v3));  
	//        TrigIsoMu22=(runNumber>=0.0&&runNumber<=276811&&(slTrigIsoMu22_v1 || slTrigIsoMu22_v2 || slTrigIsoMu22_v3 || slTrigIsoTkMu22_v1 || slTrigIsoTkMu22_v2 || slTrigIsoTkMu22_v3));  
	//        TrigIsoMu24=(runNumber>=276824&&runNumber<=999999&&(slTrigIsoMu24_v1 || slTrigIsoMu24_v2 || slTrigIsoMu24_v3 ||slTrigIsoTkMu24_v1 || slTrigIsoTkMu24_v2 || slTrigIsoTkMu24_v3));  
	//        TrigIsoMu24=((slTrigIsoMu24_v1 || slTrigIsoMu24_v2 || slTrigIsoMu24_v3 || slTrigIsoMu24_v4 ||slTrigIsoTkMu24_v1 || slTrigIsoTkMu24_v2 || slTrigIsoTkMu24_v3|| slTrigIsoTkMu24_v4));  
        TrigIsoMu24 = slTrigIsoMu24_v1 || slTrigIsoMu24_v2|| slTrigIsoMu24_v3 || slTrigIsoMu24_v4 || slTrigIsoTkMu24_v1 || slTrigIsoTkMu24_v2|| slTrigIsoTkMu24_v3|| slTrigIsoTkMu24_v4;
        }

    if(isData=="MC"){
        TrigIsoMu20=false;
        TrigIsoMu24=false;
        TrigIsoMu22= true;//Trigger efficiency:use the measured one atm --> always pass the trigger in MC// (slTrigIsoMu22_v1 || slTrigIsoMu22_v2 || slTrigIsoMu22_v3);  
        if (doSynch ){TrigIsoMu24 = slTrigIsoMu24_v1 || slTrigIsoMu24_v2|| slTrigIsoMu24_v3 || slTrigIsoMu24_v4 || slTrigIsoTkMu24_v1 || slTrigIsoTkMu24_v2|| slTrigIsoTkMu24_v3|| slTrigIsoTkMu24_v4; TrigIsoMu22=false;}
	//TrigIsoMu22= slTrigIsoMu22_v2 || slTrigIsoMu22_v3; //ckecking only for QCD if it works 
        }

    //if (TrigIsoMu20) cout << "TrigIsoMu20 fired in Event = "<<evt<<" Run Number = "<<runNumber<<endl;
    //if (TrigIsoMu22) cout << "TrigIsoMu22 fired in Event = "<<evt<<" Run Number = "<<runNumber<<endl;
    //if (TrigIsoMu24) cout << "TrigIsoMu24 fired in Event = "<<evt<<" Run Number = "<<runNumber<<endl;
  
    bool muonTrigger = (TrigIsoMu20 || TrigIsoMu22 || TrigIsoMu24); 
    bool passesAnyTrigger = muonTrigger || muonTrigger;//REMINDER: ADD ALL TRIGGERS HERE 
    if(!passesAnyTrigger)continue;
    if(doSynch && passesAnyTrigger){
      fileout_step0<<std::fixed<<std::setprecision(0)<<evtNumber<<std::endl;
    }
    
    n_trig += w;
    nev_trig +=1;

    double puUpFact=1.0;//(LumiWeightsUp_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));
    double puDownFact=1.0;//(LumiWeightsDown_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));

    if(isData=="MC"){
        LHEWeightSign[0]=w_zero/fabs(w_zero);
        w = LHEWeightSign[0];
        //    w=1.;
        w_pu = LumiWeights_.weight(numTrueInt);
        //cout << " ntrueint "<< numTrueInt <<" puweight "<<w_pu<<endl;
        //w_pu =1;

        w = w * w_pu;
        //w = w  w_pu;
    if(sample=="TTNoTT"){
        w*=w_top/topWeight;
        }
    
    if(w_pu!=0){
        puUpFact=(LumiWeightsUp_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));
        puDownFact=(LumiWeightsDown_.weight(numTrueInt))/(LumiWeights_.weight(numTrueInt));
        }
        //cout << " ntrueint "<< numTrueInt <<" puweight "<<w_pu<<" wup "<< puUpFact<<" wdown "<<puDownFact<<endl;
        } 

    //if(numTrueInt>49){
    //  cout << " --> numTrueInt very high!!" << endl;
    //  puUpFact =0;
    //  puDownFact=0;
    //  }
    met = metPt[0];
    metpx = metPx[0];
    metpy = metPy[0];
  
    TLorentzVector lep1;
    TLorentzVector lep2;
    TLorentzVector lep, mu, el;

    vector<TLorentzVector> tightEl, tightMu;
    int nMu(0.), nEl(0.);//, nVetoEl(0.), nLooseMu(0.);
    for(int e= 0; e<maxElLoop;++e ){
        if(elPt[e]>20){
            el.SetPtEtaPhiE(elPt[e], elEta[e], elPhi[e],elE[e]);
            tightEl.push_back(el);
            }
    }//end loop on electrons     
  
    vector<float> selectedIso;

    vector<float> mindeltaRmuonloop;
    if(channel!="muonantiiso"){
        for(int m= 0; m<maxMuLoop;++m ){
            if(muPt[m]>26 && muEta[m]<2.4 && muIso[m]<0.06){
                selectedIso.push_back(muIso[m]);
                mu.SetPtEtaPhiE(muPt[m], muEta[m], muPhi[m],muE[m]);
	            tightMu.push_back(mu);
                //------------
                for (int j = 0; j <maxJetLoop;++j){
                   TLorentzVector all_jets;
                   all_jets.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
                   float dRtemp = deltaR(all_jets.Eta(),all_jets.Phi(),mu.Eta(),mu.Phi());
                   mindeltaRmuonloop.push_back(dRtemp);
                   } 
                 //-----------   
                }
            }
        }
    else{
        for(int m= 0; m<maxMuLoop;++m ){
            bool passesDR=true;
            double dR=1111111;
            if(muAntiIsoPt[m]>20){
            //if(muAntiIsoPt[m]>20 && muAntiIsoIso[m]>=0.15){
                selectedIso.push_back(muAntiIsoIso[m]);
                mu.SetPtEtaPhiE(muAntiIsoPt[m], muAntiIsoEta[m], muAntiIsoPhi[m],muAntiIsoE[m]);
                for (int j = 0; j <maxJetLoop;++j){
                    TLorentzVector all_jets;
                    all_jets.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
                    float dRtemp = deltaR(all_jets.Eta(),all_jets.Phi(),mu.Eta(),mu.Phi());
                    mindeltaRmuonloop.push_back(dRtemp);
                    if(dRtemp<dR){
                        dR=dRtemp;
                        }
                    if(dR<0.4){
                        passesDR=false; //throw away the muon if deltaR btw jet and tightmu in less than 0.4
                        continue;
                        }
                    }
            if(!passesDR)continue;
                //cout <<"Evt No. after "<<evt<<"   "<<"muAntiIsoPt["<<m<< "] "<<muAntiIsoPt[m]<< " dR(jet,mu) "<< dR <<" passesDR "<< passesDR<<endl;
                tightMu.push_back(mu);
                }
            }
        }//end loop on muons
  
    nMu = tightMu.size();
    nEl = tightEl.size();
    vector<lept> leptons;
    //lepWeight1Mu(1.), lepWeight1MuLepUp(1.), lepWeight1MuLepDown(1.);
    //  lepWeight1Ele(1.), lepWeight1EleLepUp(1.), lepWeight1EleLepDown(1.);

    //Weights muTightTrigEff(file_eff_trig,"IsoMu22_OR_IsoTkMu22_PtEtaBins_Run273158_to_274093/efficienciesDATA/abseta_pt");
    //Weights muTightIdSF(file_eff_trig,"MC_NUM_TightIDandIPCut_DEN_genTracks_PAR_pt_spliteta_bin1/abseta_pt_ratio");
    //Weights muTightIsoSF(file_eff_trig,"MC_NUM_TightRelIso_DEN_TightID_PAR_pt_spliteta_bin1/abseta_pt_ratio");
    if(isData=="MC"){
        if(tightMu.size()==1){
            /*      cout << "event with 1 muon, "<< evt<< " getting trigg eff" <<endl;
            cout << muTightTrigEff.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())<<endl;
            cout << "getting id sf"<< endl;
            cout << muTightIdSF.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())<<endl;
            cout << "getting iso sf"<< endl;
            cout << muTightIsoSF.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())<<endl;
            cout << "lep sf retrieved successfully!"<<endl;*/
            lepWeight1Mu = muTightTrigEff.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())*
	        muTightIdSF.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())*
	        muTightIsoSF.getEff(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt());
            float errMu= muTightTrigEff.getErr(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())*
	        muTightIdSF.getErr(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt())*
	        muTightIsoSF.getErr(fabs(tightMu.at(0).Eta()),tightMu.at(0).Pt());
            lepWeight1MuLepUp=lepWeight1Mu+errMu;
            lepWeight1MuLepDown=lepWeight1Mu-errMu;
            }
        }
    
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
  
    float w_lep(1.0),sf_lepUp(1.0),sf_lepDown(1.0);
    if(isData=="MC"){
        if(channel=="muon"){
            w_lep=lepWeight1Mu;
            if(w_lep!=0){
	            sf_lepUp=lepWeight1MuLepUp/w_lep;
	            sf_lepDown=lepWeight1MuLepDown/w_lep;
                }
            }
        if(channel=="electron"){
            w_lep=lepWeight1Ele;
            if(w_lep!=0){
                sf_lepUp=lepWeight1EleLepUp/w_lep;
	            sf_lepDown=lepWeight1EleLepDown/w_lep;
                }
            }
            w*=w_lep;
        }
  
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
    syst0BM.setWeight("puUp",bWeight0CSVM*puUpFact);
    syst0BM.setWeight("puDown",bWeight0CSVM*puDownFact);
    syst0BM.setWeight("lepUp",bWeight0CSVM*sf_lepUp);
    syst0BM.setWeight("lepDown",bWeight0CSVM*sf_lepDown);
    
    syst1BM.copySysts(systZero);
    syst1BM.setWeight(0,bWeight1CSVM);
    syst1BM.setWeight("btagUp",bWeight1CSVMBTagUp);
    syst1BM.setWeight("btagDown",bWeight1CSVMBTagDown);
    syst1BM.setWeight("mistagUp",bWeight1CSVMMisTagUp);
    syst1BM.setWeight("mistagDown",bWeight1CSVMMisTagDown);
    syst1BM.setWeight("lepUp",bWeight1CSVM*sf_lepUp);
    syst1BM.setWeight("lepDown",bWeight1CSVM*sf_lepDown);
    
    syst2BM.copySysts(systZero);
    syst2BM.setWeight(0,bWeight2CSVM);
    syst2BM.setWeight("btagUp",bWeight2CSVMBTagUp);
    syst2BM.setWeight("btagDown",bWeight2CSVMBTagDown);
    syst2BM.setWeight("mistagUp",bWeight2CSVMMisTagUp);
    syst2BM.setWeight("mistagDown",bWeight2CSVMMisTagDown);
    syst2BM.setWeight("puUp",bWeight2CSVM*puUpFact);
    syst2BM.setWeight("puDown",bWeight2CSVM*puDownFact);
    syst2BM.setWeight("lepUp",bWeight2CSVM*sf_lepUp);
    syst2BM.setWeight("lepDown",bWeight2CSVM*sf_lepDown);
    
    vector<float> jetsPhi;
    
    vector<TLorentzVector> bjets, jets_nob, jets, jets20;    

    struct btag{
        TLorentzVector vect;
        float csv;
        float csvreweight;
        float csvreweightsd;
        float partonFlavour;
        float hadronFlavour;
        };
  
    std::vector<btag> bvects,extrajets;
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
  
    for (int j = 0; j <maxJet20Loop;++j){
        //     cout <<" in j20 loop"<<endl;
        TLorentzVector jet20;
        if(jet20Pt[j]<40){
            //cout <<" # j= "<<j<< " pt "<< jet20Pt[j] <<endl;
            jet20.SetPtEtaPhiE(jet20Pt[j], jet20Eta[j], jet20Phi[j], jet20E[j]);
            jets20.push_back(jet20);
            btag bj;
            bj.vect = jet20;
            bj.csv = jet20ak4chs_csvv2[j];
            bj.csvreweight = jet20ReshapeFactorCSV[j];
            bj.csvreweightsd = 1.;//jet20ak4chsReweightFactorCSV_SD[j];
            if( doTtoSDDecay){
                bj.csvreweightsd = jet20ReshapeFactorCSV_SD[j];
	            bj.partonFlavour=jet20PartonFlavour[j];
                }
            extrajets.push_back(bj);
            }
        }
   
    for (int j = 0; j <maxJetLoop;++j){
        TLorentzVector jet, all_jets;
        all_jets.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);

        bool passesDR=true;
        double dR=0; 
        //--------------------------------------------------------------------------------
      
        if(channel=="muonantiiso"){
            for (size_t i = 0; i<(size_t)tightMu.size();++i){
                //dR = deltaR(jetEta[j],jetPhi[j],tightMu[i].Eta(),tightMu[i].Phi());   
                dR= deltaR(all_jets.Eta(),all_jets.Phi(),tightMu[i].Eta(),tightMu[i].Phi());
                if (dR<0.4) passesDR=false; //throw away the jet if deltaR btw jet and tightmu in less than 0.4
                }
            for (size_t i = 0; i<(size_t)muLooseSize;++i){ 
                dR= deltaR(all_jets.Eta(),all_jets.Phi(),muLooseEta[i],muLoosePhi[i]);   
	            if (dR<0.4) passesDR=false; //throw away the jet if deltaR btw jet and tightmu in less than 0.4
                }
            } 
     
     if(!passesDR)continue;
     //cout <<"Evt No. after "<<evt<<"   "<<"jet["<<j<< "] "<<jetPt[j]<< " dR(jet,mu) "<< dR <<endl;
     
     jets.push_back(all_jets);
     //     bool btagcond = jetIsCSVM[j]>0.;
     //if(useCSVTSelection){
     bool btagcond = jetPassesB[j]>0.;
	// }
     
     if(!btagcond){
         jet.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
         jets_nob.push_back(jet);
     }
     
     jetsPhi.push_back(jetPhi[j]);
     TLorentzVector bjet;
     
     if( btagcond && abs(jetEta[j])<2.4){
       bjet.SetPtEtaPhiE(jetPt[j], jetEta[j], jetPhi[j], jetE[j]);
       bjets.push_back(bjet);
       btag b;
       b.vect = bjet;
       b.csv = jetak4chs_csvv2[j];
       bvects.push_back(b);
     }
    }//End of jet loop
    
    nJets = jetsPhi.size();
    nCSVJets=bjets.size();
    //  cout <<"passes trig "<<endl; 
    //  std::sort(jets20.begin(), jets20.end(), by_pt()); 
    std::sort(bvects.begin(), bvects.end(), by_pt_jet()); 
    std::sort(extrajets.begin(), extrajets.end(), by_pt_jet()); 
    //  bool passmuon = muonTrigger && nMu == 1 && muLooseSize==1 && nEl==0 && elLooseSize == 0;
    
    bool passmuon = muonTrigger && nMu == 1 ;
    if(passmuon){n_lepton+=w;nev_lepton+=1;}
    if(doSynch && passmuon){
      //RunNumber | EvtNumber | Lumi | #TightMu | #LooseMu | #LooseEl | #Jets | #BJets | TightMu Pt | TightMuRelIso | AddLoose.Mu Pt | Add.LooseEl Pt | 1stJetPt | 2ndJetPt | 1stB Dis.| 2ndB Dis.| MET | MTW | 
      //string spaces ="   ";
      //fileout_step1<<std::fixed<<std::setprecision(0)<<runNumber<<"   "<<evtNumber<<"   "<<lumiSec<<"   "<<std::setprecision(4)<<muSize<<"   "<<muLooseSize<<"   "<<elLooseSize
      //<<"   "<<jetSize<<"   "<<bjets.size()<<"   "<<muPt[0]<<"   "<<selectedIso[0]
      //<<std::endl;
      fileout_step1<<std::fixed<<std::setprecision(0)<<evtNumber<<endl;
    }
    
    passmuon = passmuon && muLooseSize==1;
    if(passmuon){ n_loose_veto+=w; nev_loose_veto+=1;}
    if(doSynch && passmuon){
      fileout_step2<<std::fixed<<std::setprecision(0)<<evtNumber<<std::endl;
    }
    
    passmuon = passmuon && nEl==0 && elLooseSize == 0;
    if(passmuon){n_lepton_cross_veto+=w;nev_lepton_cross_veto+=1;}
    if(doSynch && passmuon){
      fileout_step3<<std::fixed<<std::setprecision(0)<<evtNumber<<std::endl;
    }
    
    if(topsize==1)topcharge=1.;
    else if(antitopsize==1)topcharge=-1.;
    //bool passantiisomuon = muonTrigger && nMu == 1 && muLooseSize == 0 && nEl==0 && elLooseSize == 0;
    //bool passantiisomuon = muonTrigger && muAntiIsoSize==1 && nMu==1 && nEl==0 && elLooseSize == 0;
    bool passantiisomuon = muonTrigger && muAntiIsoSize==1 && nMu==1 && nEl==0 && elLooseSize == 0;
    
    bool passelectron = false;
    bool passsinglelepton = passmuon || passelectron;

    if(channel=="muonantiiso" && !passantiisomuon) continue;
    if(channel=="muon" && !passmuon)continue; 
    //cout <<"Evt No. after "<<evt<<" tightMuons "<<nMu<<" Loose Muons "<<muLooseSize<<endl; 
    
    if(channel=="electron" && !passelectron)continue;
    if(channel=="singlelepton" && !passsinglelepton)continue;
  
    
    if(jets.size() ==2){
      n_2j+=w; nev_2j+=1;
      if(doSynch){
	fileout_step4<<std::fixed<<std::setprecision(0)<<evtNumber<<std::endl;
      }
    }
    
    if(jets.size() ==2 && bjets.size()==1){ n_2j1t+=w; nev_2j1t+=1;  
      
      if(doSynch){
	fileout_step5<<std::fixed<<std::setprecision(0)<<evtNumber<<std::endl;
      }
      if(doSynch){
	if((tightMu.size())==1){
	  TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
	  float phi_lmet = fabs(deltaPhi(tightMu[0].Phi(), metPhi[0]) );
	  mt = sqrt(2* tightMu[0].Pt() * met* ( 1- cos(phi_lmet)));
	}
	if(mt> 50){
	  n_2j1tmtw+=w;
	  nev_2j1tmtw+=1.;
	  
	  if(doSynch){
	    fileout_step6<<std::fixed<<std::setprecision(0)<<evtNumber<<std::endl;
	  }
	  cout <<std::fixed<<std::setprecision(0)<<"---Evt No: "<<evtNumber<<endl;
	  //cout <<"Jet[0] (pt,eta,E) ("<<jets[0].Pt()<<","<<jets[0].Eta()<<","<<jets[0].E()<<")" 
	  //<<" | Jet[1] (pt,eta,E) ("<<jets[1].Pt()<<","<<jets[1].Eta()<<","<<jets[1].E()<<")" 
	  //<<" | Mu[0] (pt,eta,E) "<<tightMu[0].Pt()<<","<<tightMu[0].Eta()<<","<<tightMu[0].E() <<" | RelIso. "<<selectedIso[0]
	  //<<" | bjet Pt "<<bjets[0].Pt()
	//<<" | MET "<<met<< " | metPhi "<<metPhi[0]
	//<<" | mt "<< sqrt(2* tightMu[0].Pt() * met* ( 1- cos(fabs(deltaPhi(tightMu[0].Phi(), metPhi[0]))))) 
	//<<endl;
	}
      }
    }
    
  
  
    //2j0t 
    if((jets.size() == 2 && bjets.size() == 0)){
      for(size_t i = 0; i < (size_t)tightMu.size();++i ){
	if((tightMu.size()) == 1){
            TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
            float phi_lmet = fabs(deltaPhi(tightMu[i].Phi(), metPhi[0]) );
            mt = sqrt(2* tightMu[i].Pt() * met* ( 1- cos(phi_lmet)));
	        //cout << "muon iso "<< selectedIso[i]<< " mtw "<< mt <<endl;
            syst0BM.fillHistogramsSysts(h_2j0t_mtw,mt,w);
            }
        }
        for(size_t j= 0; j< (size_t)jets.size();++j ){
        //cout <<"Evt No. "<<evt<<"   "<<"   "<<"jet["<<j<< "] "<<jets[j].Pt()<<"   bjetsize "<< bjets.size() <<"   Mu[0] "<<tightMu[0].Pt()<< std::endl;
        //cout << "weight 0 "  << bWeight0CSVM <<" w "<< w <<endl;
            if(j==0){
                syst0BM.fillHistogramsSysts(h_2j0t_jetpt40_1st,jets[j].Pt(),w);
                syst0BM.fillHistogramsSysts(h_2j0t_MuEta,tightMu[0].Eta(),w);
                syst0BM.fillHistogramsSysts(h_2j0t_MuPhi,tightMu[0].Phi(),w);     
                syst0BM.fillHistogramsSysts(h_2j0t_dR_lepjetpt40_1st,deltaR(jets[j].Eta(),jets[j].Phi(),tightMu[0].Eta(),tightMu[0].Phi()),w);
                syst0BM.fillHistogramsSysts(h_2j0t_dPhi_lepjetpt40_1st,fabs(deltaPhi(jets[j].Phi(),tightMu[0].Phi())),w);
                syst0BM.fillHistogramsSysts(h_2j0t_dEta_lepjetpt40_1st,jets[j].Eta()-tightMu[0].Eta(),w);
                }
      
            if(j==1)syst0BM.fillHistogramsSysts(h_2j0t_jetpt40_2nd,jets[j].Pt(),w);
            syst0BM.fillHistogramsSysts(h_2j0t_incljetabseta,fabs(jets[j].Eta()),w);  

            if(mt>50.0){
                if(j==0){
                syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_jetpt40_1st,jets[j].Pt(),w);  
                syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_jeteta40_1st,jets[j].Pt(),w);  
                syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_jetabseta40_1st,fabs(jets[j].Eta()),w);  
                syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_incljeteta,jets[j].Eta(),w);  
                syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_incljetabseta,fabs(jets[j].Eta()),w);  
	            }
        
            if(j==1){
                syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_jetpt40_2nd,jets[j].Pt(),w);  
                syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_incljeteta,jets[j].Eta(),w);  
                syst0BM.fillHistogramsSysts(h_2j0t_mtwcut_incljetabseta,fabs(jets[j].Eta()),w);  
                }
            }
        }
    }


    double etathreshold=2.5;
    //2j1t
    if((jets.size() == 2 && bjets.size() == 1)){
    for(size_t b= 0; b< (size_t)bjets.size();++b ){
      //      cout << " w "<< w <<" b weight 1 "<< bWeight1CSVM << endl;
      //      cout <<" pu weight = "<< w_pu<<endl; 
      //      cout <<" lep weight = "<< w_lep<<endl; 
      
      if(b==0){
          syst1BM.fillHistogramsSysts(h_2j1t_bjetpt,bjets[b].Pt(),w,NULL,false);
          //cout << bjets[i].Pt()<< endl;
          }
    }
    
    bool qcddepleted=false, signalenriched=false;
    //define QCD depletion condition:
    for(size_t i = 0; i < (size_t)tightMu.size();++i ){
        syst1BM.fillHistogramsSysts(h_2j1t_MuPt,tightMu[i].Pt(),w);
        syst1BM.fillHistogramsSysts(h_2j1t_MuEta,tightMu[i].Eta(),w);
        syst1BM.fillHistogramsSysts(h_2j1t_MuPhi,tightMu[i].Phi(),w);
        syst1BM.fillHistogramsSysts(h_2j1t_MuE,tightMu[i].E(),w);
	if((tightMu.size()) == 1 ){
            TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
            float phi_lmet = fabs(deltaPhi(tightMu[i].Phi(), metPhi[0]) );
            mt = sqrt(2* tightMu[i].Pt() * met* ( 1 - cos(phi_lmet)));
            //cout <<" mt = "<<mt<<endl;
            //cout <<" calculate_mtw() = "<<calculate_mtw(metPt, metPhi,tightMu)<<endl;     
            syst1BM.fillHistogramsSysts(h_2j1t_mtw,mt,w);
            if (calculate_mtw(metPt, metPhi,tightMu) > 50.0){
                qcddepleted=true;
                syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_mtw,mt,w);
                } 
            }
        }

    //define signal enriching condition:
    for (size_t j= 0; j< (size_t)jets.size();++j ){
      bool btagcond = jetPassesB[j]>0.;
      //      if(useCSVTSelection){
      //          btagcond = jetIsCSVT[j]>0.;
      //          }
      if(!(btagcond)){
	syst1BM.fillHistogramsSysts(h_2j1t_ljetpt,jets[j].Pt(),w);
	syst1BM.fillHistogramsSysts(h_2j1t_ljeteta,jets[j].Eta(),w);
	syst1BM.fillHistogramsSysts(h_2j1t_etajprime,fabs(jets[j].Eta()),w);
	signalenriched=signalenriched||fabs(jets[j].Eta())>etathreshold;
	if(qcddepleted) {
	  *etajprime_2j1t_mtwcut = fabs(jets[j].Eta());
	  syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_ljetpt,jets[j].Pt(),w);
	  syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_ljeteta,jets[j].Eta(),w);
	  syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_etajprime,fabs(jets[j].Eta()),w);
	}
      } 
    }
    //    cout << " filling extrajets "<<endl;
    syst1BM.fillHistogramsSysts(h_2j1t_nextrajets,extrajets.size(),w);
    //    cout << " filled nextrajets extrajets "<<endl;
    if(extrajets.size()){
      
      //      syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetpt,jets20.at(0).Pt(),w);
      //      cout << "method 1 ok "<<endl;
      syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);
      //      cout << "method 2 ok "<<endl;
      float weightreshape=extrajets.at(0).csvreweight;
      bool fillcond=true;

      if( doTtoSDDecay){
	fillcond = !(extrajets.at(0).partonFlavour ==5 && extrajets.at(0).partonFlavour*topcharge>0&&fabs(extrajets.at(0).vect.Eta())<2.4);
      }
      //      cout << " fillcond, parton flavour "<< extrajets.at(0).partonFlavour<< " charge "  <<topcharge<<  " eta "   <<  fabs(extrajets.at(0).vect.Eta())<<endl;
      if(fillcond){
	syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);
	syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetcsv,extrajets.at(0).csv,w*weightreshape);
      }
      else{
	//	cout << " not fillcond, filling top from b"<<endl;
	float weight_sd_b=extrajets.at(0).csvreweightsd;
	syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetcsv_sd_b,extrajets.at(0).csv,w*weightreshape*weight_sd_b);
	syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetcsv_reshape_sd_b, weight_sd_b,w);
	if(reweightall){syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetpt_sd_b,extrajets.at(0).vect.Pt(),w*weight_sd_b);}
	else {syst1BM.fillHistogramsSysts(h_2j1t_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);}
      }
      
      //      cout << "csv ok"<<endl;
    }
    //    cout << "extrajets ok"<<endl;
    
    //Adding top related variables
    vector<TLorentzVector> tops= topUtils.top4Momenta(tightMu, bjets ,metpx, metpy);
    
    if(tops.size()!=0){
      double mass = tops.at(0).M();
      double pt = tops.at(0).Pt();
      *topMass_2j1t_mtwcut=mass;*topMass_2j1t_mtwcut_sr=mass;
      syst1BM.fillHistogramsSysts(h_2j1t_topMass,mass,w);
      syst1BM.fillHistogramsSysts(h_2j1t_topPt,pt,w);
      //      if(calculate_mtw(metPt, metPhi,tightMu) > 50.0) {
      if(qcddepleted){ 

	syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_topMass,mass,w);
        syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_topPt,pt,w);
      }
      if(signalenriched && qcddepleted){
	syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_topMass,mass,w);
        syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_topPt,pt,w);
      }
    }
    //Adding extra-jet related variables: only for QCD and signal enriched regions.
    //    cout << "tops ok"<<endl;
    if(qcddepleted){
       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_nextrajets,extrajets.size(),w);
       *nextrajets_2j1t_mtwcut=extrajets.size();

       *leadingextrajetpt_2j1t_mtwcut=-99;
       *leadingextrajetcsv_2j1t_mtwcut=-99;
       *topMassExtra_2j1t_mtwcut=-99;

       *leadingextrajetcsvweight_2j1t_mtwcut=0;
       *leadingextrajetcsvweight_sd_2j1t_mtwcut=0;

       if(signalenriched){
	 syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_nextrajets,extrajets.size(),w);
	 *nextrajets_2j1t_mtwcut_sr=extrajets.size();
	 
	 *leadingextrajetpt_2j1t_mtwcut_sr=-99;
	 *leadingextrajetcsv_2j1t_mtwcut_sr=-99;
	 *topMassExtra_2j1t_mtwcut_sr=-99;

	 *leadingextrajetcsvweight_2j1t_mtwcut_sr=0;
	 *leadingextrajetcsvweight_sd_2j1t_mtwcut_sr=0;

       }
       
       if(extrajets.size()){
	 vector<TLorentzVector> topsextra= topUtils.top4Momenta(tightMu, jets20 ,metpx, metpy);

	 float weightreshape=extrajets.at(0).csvreweight;
	 bool fillcond=true;
	 if( doTtoSDDecay){
	   fillcond = (!(extrajets.at(0).partonFlavour ==5
			&& (extrajets.at(0).partonFlavour*topcharge>0)
			 && fabs(extrajets.at(0).vect.Eta())<2.4));		 
	 }
	 
	 if(topsextra.size()!=0){
	   *leadingextrajetcsvweight_2j1t_mtwcut=extrajets.at(0).csvreweight;
	   if(fillcond){
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_topMassExtra,topsextra[0].M(),w);
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv,extrajets.at(0).csv,w*weightreshape);
	     *leadingextrajetcsvweight_sd_2j1t_mtwcut=0;
	   }
	   else{
	     float weight_sd_b=extrajets.at(0).csvreweightsd;
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv_sd_b,extrajets.at(0).csv,w*weightreshape*weight_sd_b);
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv_reshape_sd_b, weight_sd_b,w);
	     *leadingextrajetcsvweight_sd_2j1t_mtwcut=weight_sd_b;
	     if(reweightall){
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_topMassExtra_sd_b,topsextra[0].M(),w*weight_sd_b);
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt_sd_b,extrajets.at(0).vect.Pt(),w*weight_sd_b);
	     }
	     else{
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_topMassExtra,topsextra[0].M(),w);
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);
	     }

	   }
	   //	   cout << "extrajet csvs ok"<<endl;
	   
	   *leadingextrajetpt_2j1t_mtwcut=extrajets.at(0).vect.Pt();
	   *leadingextrajetcsv_2j1t_mtwcut=extrajets.at(0).csv;
	   *topMassExtra_2j1t_mtwcut=topsextra[0].M();
	 }
	 if(signalenriched){
	   *leadingextrajetcsvweight_2j1t_mtwcut_sr=extrajets.at(0).csvreweight;
	   if(fillcond){
	   //	   cout<< " filling histograms for the signal enriched region"<<endl;
	   syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra,topsextra[0].M(),w);
	   syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt,extrajets.at(0).vect.Pt(),w);
	   syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv,extrajets.at(0).csv,w*weightreshape);
	   *leadingextrajetpt_2j1t_mtwcut_sr=extrajets.at(0).vect.Pt();
	   *leadingextrajetcsv_2j1t_mtwcut_sr=extrajets.at(0).csv;
	   *topMassExtra_2j1t_mtwcut_sr=topsextra[0].M();
	   *leadingextrajetcsvweight_sd_2j1t_mtwcut_sr=0;
	   }
	   else{
	     float weight_sd_b=extrajets.at(0).csvreweightsd;
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b,extrajets.at(0).csv,w*weightreshape*weight_sd_b);
	     syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv_reshape_sd_b, weight_sd_b,w);
	     *leadingextrajetcsvweight_sd_2j1t_mtwcut_sr=weight_sd_b;
	     if(reweightall){
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra_sd_b,topsextra[0].M(),w*weight_sd_b);
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt_sd_b,extrajets.at(0).vect.Pt(),w*weight_sd_b);
	     }
	     else{
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra_sd_b,topsextra[0].M(),w);
	       syst1BM.fillHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt_sd_b,extrajets.at(0).vect.Pt(),w);
	     }
	   }
	   //	   cout << "extrajet csvs sr ok"<<endl;
	 }
       }
    }
    if(qcddepleted){
      *w_2j1t_mtwcut=w;
      if(addTrees)syst1BM.fillTreesSysts(trees1T,"2j1t_mtwcut");
    }
    if(signalenriched && qcddepleted) {
      *w_2j1t_mtwcut_sr=w;
      if(addTrees)syst1BM.fillTreesSysts(trees1T,"2j1t_mtwcut_sr");}
    //AOB for 2j1t?
    }
  //End 2j1t

  //3j1t   
  if((jets.size() == 3 && bjets.size() ==1)) {
  for(size_t b= 0; b< (size_t)bjets.size();++b ){
  if(b==0)systZero.fillHistogramsSysts(h_3j1t_bjetpt,bjets[b].Pt(),w);
   } 
  for(size_t i = 0; i < (size_t)tightMu.size();++i ){
        syst1BM.fillHistogramsSysts(h_3j1t_MuPt,tightMu[i].Pt(),w);
        syst1BM.fillHistogramsSysts(h_3j1t_MuEta,tightMu[i].Eta(),w);
        syst1BM.fillHistogramsSysts(h_3j1t_MuPhi,tightMu[i].Phi(),w);
        syst1BM.fillHistogramsSysts(h_3j1t_MuE,tightMu[i].E(),w);
        if((tightMu.size()) == 2 ){
          TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
          float phi_lmet = fabs(deltaPhi(tightMu[i].Phi(), metPhi[0]) );
          mt = sqrt(2* tightMu[i].Pt() * met* ( 1- cos(phi_lmet)));
          systZero.fillHistogramsSysts(h_3j1t_mtw,mt,w);
          }
      }
  //AOB for 3j1t?
  }
  //End 3j1t


  //3j2t
  if((jets.size() == 3 && bjets.size()==2)){
    //      cout << " w "<< w <<" b weight 2 "<< bWeight2CSVM << endl;
    //      cout <<" pu weight = "<< w_pu<<endl; 
    //      cout <<" lep weight = "<< w_lep<<endl; 

    //    cout << " event "<<evt<<endl;
    for(size_t b=0; b< (size_t)bjets.size();++b ){
      //      cout << " ,bWeight2CSV "<< bWeight2CSVM <<",bWeight2CSVMBTagUp "<<bWeight2CSVMBTagUp<<endl;
      // cout << " i is "<< i << endl;
      if(b==0)syst2BM.fillHistogramsSysts(h_3j2t_bjetpt,bjets[b].Pt(),w,NULL,false);
      if(b==1)syst2BM.fillHistogramsSysts(h_3j2t_2ndbjetpt,bjets[b].Pt(),w);
    }
    for(size_t i = 0; i < (size_t)tightMu.size();++i ){
      if(tightMu[i].Pt()>20){
        syst2BM.fillHistogramsSysts(h_3j2t_MuPt,tightMu[i].Pt(),w);
        syst2BM.fillHistogramsSysts(h_3j2t_MuEta,tightMu[i].Eta(),w);
        syst2BM.fillHistogramsSysts(h_3j2t_MuPhi,tightMu[i].Phi(),w);
        syst2BM.fillHistogramsSysts(h_3j2t_MuE,tightMu[i].E(),w);
        if((tightMu.size()) == 1 ){
          TVector2 met_( met*cos(metPhi[0]), met*sin(metPhi[0]));
          float phi_lmet = fabs(deltaPhi(tightMu[i].Phi(), metPhi[0]) );
          mt = sqrt(2* tightMu[i].Pt() * met* ( 1- cos(phi_lmet)));
          syst2BM.fillHistogramsSysts(h_3j2t_mtw,mt,w);
	  *mtw_3j2t=mt;
	}
      }
    }
    for(size_t j=0; j< (size_t)jets.size();++j ){
      //      bool btagcond=false;
      //      bool btagcond = jetIsCSVM[j]>0.;

      bool btagcond = jetPassesB[j]>0.;
     
      if(!btagcond){
	syst2BM.fillHistogramsSysts(h_3j2t_ljetpt,(jets.at(j).Pt()),w);
	syst2BM.fillHistogramsSysts(h_3j2t_ljeteta,(jets.at(j).Eta()),w);
	syst2BM.fillHistogramsSysts(h_3j2t_etajprime,(fabs(jets.at(j).Eta())),w);
	*etajprime_3j2t=fabs(jets.at(j).Eta());
      }
    }

    vector<TLorentzVector> tops= topUtils.top4Momenta(tightMu, bjets ,metpx, metpy);
    if(tops.size()!=0){
      double massLeading = tops.at(0).M();
      double massSecond = tops.at(1).M();
      *topMassLeading_3j2t=massLeading;
      *topMassSecond_3j2t=massSecond;
    }    
    syst2BM.fillHistogramsSysts(h_3j2t_nextrajets,extrajets.size(),w);
    *nextrajets_3j2t= extrajets.size();
    
    *w_3j2t=w;
    if(addTrees)syst2BM.fillTreesSysts(trees2T,"3j2t");
    //AOB for 3j2t?
  }
  //End 3j2t
  
  //To add here: 2l3j2t: 2 leptons, 3 jets, 2 tags for ttbar control

  //End signal and control regions.

  systZero.fillHistogramsSysts(h_nJets,nJets,w); 
  systZero.fillHistogramsSysts(h_nbJets,nCSVJets,w); 
  systZero.fillHistogramsSysts(h_nPV,nPV,w);
  systZero.fillHistogramsSysts(h_nGoodPV,nPV,w);
  systZero.fillHistogramsSysts(h_nTruePV,numTrueInt,w);
  
  //if(nTightElectrons != nEl)cout << "warning! problem with tight el"<<endl;
  //if(nTightMuons != nMu)cout << "warning! problem with tight mu"<<endl;
  
  nTightElectrons = nEl;
  nTightMuons = nMu;    


  if(doSynch)fileout<<std::fixed<<std::setprecision(0)<<runNumber<<"   "<<evtNumber<<"   "<<lumiSec<<"   "<<std::endl;
  
  }//end of loop over events 
  
 
 if(doSynch){
   //  int nev_trig(0),nev_lepton(0),nev_loose_veto(0),nev_lepton_cross_veto(0),nev_2j(0),nev_2j1t(0);
   fileout<<" ntrig "<< " nlep "<< " nlooseveto "<< " nloosevetoele "<< " n2jets "<< " n2jets1tag "<< " n2jets1tagmtw "<<endl;
   fileout<< " " << nev_trig << " " << nev_lepton << " " << nev_loose_veto << " " << nev_lepton_cross_veto << " " << nev_2j << " " << nev_2j1t <<" "<< nev_2j1tmtw <<endl;
   fileout<<" wtrig "<< " wlep "<< " wlooseveto "<<" wloosevetoele "<< " w2jets "<< " w2jets1tag "<<" w2jets1tagmtw "<<endl;
   fileout<< " " << n_trig << " " << n_lepton << " " << n_loose_veto <<n_lepton_cross_veto << " " << n_2j << " " << n_2j1t <<" "<< n_2j1tmtw<<endl;

  fileout.close();  
 }//return h

 //Step 0:get only the negative weights version, 
 if(isData=="MC"){
   //chainNEvents.Project("h_weight_sign","Event_LHEWeight10/abs(Event_LHEWeight10)"); // I messed up and weight zero is not  saved in the full chain... but this one should do the trick as well! Just do h_weight_sign->GetMean() and multiply this to get the correct number :)
   cout << " mean  "<< zeroWeight<<endl;
 }
 // TO FIX! ADD PROJECTION OF ALL PDFS TO THE EVENT!
 // systZero.projectAllPDF(h_weightZero,chainNEvents);

 bool hasLHEWeights=true; 
 cout << "haslheweights?" << hasLHEWeights<<endl;
 if(sample.find("QCDMu")!=std::string::npos || sample.find("_tW_")!=std::string::npos || isData=="DATA")hasLHEWeights=false;
 // Cut Flow
 
 if(hasLHEWeights) {
   h_cutFlow->SetBinContent(0,nEventsTot*zeroWeight);//Underflow: number of events pre-preselection.}
   cout << " nEvents tot "<< nEventsTot << " * weight"<< zeroWeight<<" = "<< h_cutFlow->GetBinContent(0)<<endl; 
 }
 else{ h_cutFlow->SetBinContent(0,nEventsTot);//Underflow: number of events pre-preselection.
   cout << " nEvents tot "<< nEventsTot <<" = "<< h_cutFlow->GetBinContent(0)<<endl; 
}

  h_cutFlow->SetBinContent(1,nEvents);
  h_cutFlow->GetXaxis()->SetBinLabel(1,"no selection");
  h_cutFlow->SetBinContent(2, n_trig);
  h_cutFlow->GetXaxis()->SetBinLabel(2, "trigger");
  h_cutFlow->SetBinContent(3,n_lepton);
  h_cutFlow->GetXaxis()->SetBinLabel(3, "=1 Iso. Muon");
  h_cutFlow->SetBinContent(4,n_lepton_cross_veto);
  h_cutFlow->GetXaxis()->SetBinLabel(4, "=0 Loose Ele.");
  h_cutFlow->SetBinContent(5,n_2j);
  h_cutFlow->GetXaxis()->SetBinLabel(5, "2-jet");
  h_cutFlow->SetBinContent(6,n_2j1t);
  h_cutFlow->GetXaxis()->SetBinLabel(6, "2-jet 1-tag");

  //Rescale the extra jet csv by the average of the sf
  systZero.rescaleHistograms(h_2j1t_leadingextrajetcsv_sd_b,h_2j1t_leadingextrajetcsv_reshape_sd_b);
  systZero.addHistograms(h_2j1t_leadingextrajetcsv,h_2j1t_leadingextrajetcsv_sd_b);

  systZero.rescaleHistograms(h_2j1t_mtwcut_leadingextrajetcsv_sd_b,h_2j1t_mtwcut_leadingextrajetcsv_reshape_sd_b);
  systZero.addHistograms(h_2j1t_mtwcut_leadingextrajetcsv,h_2j1t_mtwcut_leadingextrajetcsv_sd_b);

  systZero.rescaleHistograms(h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b,h_2j1t_mtwcut_sr_leadingextrajetcsv_reshape_sd_b);
  systZero.addHistograms(h_2j1t_mtwcut_sr_leadingextrajetcsv,h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b);


  //Write the Histogramms here  
  systZero.writeSingleHistogramSysts(h_cutFlow, allMyFiles); 
  systZero.writeSingleHistogramSysts(h_weight_sign, allMyFiles); 
  systZero.writeHistogramsSysts(h_nPV, allMyFiles); 
  systZero.writeHistogramsSysts(h_nGoodPV, allMyFiles); 
  systZero.writeHistogramsSysts(h_nTruePV, allMyFiles); 
  
  systZero.writeHistogramsSysts(h_2j0t_dR_lepjetpt40_1st,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_dPhi_lepjetpt40_1st,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_dEta_lepjetpt40_1st,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_MuEta,      allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_MuPhi,      allMyFiles); 
  
  systZero.writeHistogramsSysts(h_2j0t_jetpt40_1st, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_jetpt40_2nd, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtw,   allMyFiles); 

  systZero.writeHistogramsSysts(h_2j0t_mtwcut_jeteta40_1st,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_jetabseta40_1st,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_jetpt40_1st,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_jetpt40_2nd,   allMyFiles); 
  
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_incljeteta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_mtwcut_incljetabseta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j0t_incljetabseta,   allMyFiles); 
  
  //2j1t
  systZero.writeHistogramsSysts(h_2j1t_bjetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_ljetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_ljeteta,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_etajprime,  allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuPt,    allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuEta,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuPhi,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuE,     allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_MuCharge,allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtw,     allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_topMass, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_topPt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_nextrajets,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_leadingextrajetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_leadingextrajetcsv,   allMyFiles); 
  if(doTtoSDDecay){
    systZero.writeHistogramsSysts(h_2j1t_leadingextrajetpt_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_leadingextrajetcsv_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_leadingextrajetcsv_reshape_sd_b,   allMyFiles); 
  }

  systZero.writeHistogramsSysts(h_2j1t_mtwcut_mtw,     allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_topMass, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_topPt, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_bjetpt, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_ljetpt, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_ljeteta, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_etajprime, allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_nextrajets,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_topMassExtra,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_topMassExtra,   allMyFiles); 
  if(doTtoSDDecay){
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetpt_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_topMassExtra_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_leadingextrajetcsv_reshape_sd_b,   allMyFiles); 
  }
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_topMass,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_topPt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_nextrajets,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv,   allMyFiles); 
  systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra,   allMyFiles); 
    if(doTtoSDDecay){
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetpt_sd_b,   allMyFiles); //Part for the reweighted csv
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_topMassExtra_sd_b,   allMyFiles); 
    systZero.writeHistogramsSysts(h_2j1t_mtwcut_sr_leadingextrajetcsv_reshape_sd_b,   allMyFiles); 
  }

  
 


  //3j1t
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

  systZero.writeHistogramsSysts(h_3j2t_ljetpt,  allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_ljeteta,allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_etajprime,allMyFiles); 
  systZero.writeHistogramsSysts(h_3j2t_nextrajets,allMyFiles); 
  
  systZero.writeHistogramsSysts(h_nJets, allMyFiles); 
  systZero.writeHistogramsSysts(h_nbJets, allMyFiles); 


  if(addTrees){
    syst1BM.writeTreesSysts(trees1T,outTreeFile);
    syst2BM.writeTreesSysts(trees2T,outTreeFile);
  }
  
  systZero.closeFilesSysts(allMyFiles);
  
  std::cout<< "Info: Total No. of events for sample["<<sample<<"] : "<<nEvents<<std::endl;
}

