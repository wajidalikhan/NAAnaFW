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
  
  TString treePath = "DMTreesDumper/ttDM__noSyst";

  TChain chain(treePath);
  chain.AddFileInfoList(fc.GetList());

  Int_t nEvents = (Int_t)chain.GetEntries();
  //std::cout<<"Info: Number of Events: "<<nEvents<< endl;
  nEvents = std::min(nEvents, 1000);
  
  TString weightedSystsNames (weightedSysts sy);
  
  systWeights systZero,syst0B; 
  int maxSysts=0; 
  int sizeMax=50;
  int jet_size;
  //float passTrigHT(0.), Ht(0.);
  float Ht(0.);
  float runNumber(0.), lumiSec(0.);
  double evtNumber(0.);
  float jet_e[sizeMax], jet_pt[sizeMax], jet_phi[sizeMax], jet_eta[sizeMax];
  float jet_iscsvl[sizeMax], jet_iscsvm[sizeMax], jet_iscsvt[sizeMax],jet_isloose[sizeMax],jet_istight[sizeMax];

  bool onlyNominal=false;
  systZero.setOnlyNominal(onlyNominal);

  bool addPDF=false,addQ2=false,addTopPt=false,addVHF=false,addTTSplit=false;
  if(sample=="TT"){
    addTTSplit=true;
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
  //int nPDF=102;
  //if(isData=="DATA"){addPDF=false, addQ2=false;addVHF=false;addTTSplit=false;} 
  
  chain.SetBranchAddress("jetsAK4CHS_E",      &jet_e);
  chain.SetBranchAddress("jetsAK4CHS_Pt",     &jet_pt);
  chain.SetBranchAddress("jetsAK4CHS_Phi",    &jet_phi);
  chain.SetBranchAddress("jetsAK4CHS_Eta",    &jet_eta);
  chain.SetBranchAddress("jetsAK4CHS_size",   &jet_size);
  chain.SetBranchAddress("jetsAK4CHS_IsCSVL", &jet_iscsvl);
  chain.SetBranchAddress("jetsAK4CHS_IsCSVM", &jet_iscsvm);
  chain.SetBranchAddress("jetsAK4CHS_IsCSVT", &jet_iscsvt);
  chain.SetBranchAddress("jetsAK4CHS_IsLoose",&jet_isloose);
  chain.SetBranchAddress("jetsAK4CHS_IsTight",&jet_istight);
  chain.SetBranchAddress("Event_Ht", &Ht);
  chain.SetBranchAddress("Event_RunNumber", &runNumber);
  chain.SetBranchAddress("Event_LumiBlock", &lumiSec);
  chain.SetBranchAddress("Event_EventNumber", &evtNumber);

  /********************************************************************************/
  /**************                    Histogram booking              ***************/
  /********************************************************************************/
  TH1F *h_nJets[maxSysts] ;systZero.initHistogramsSysts(h_nJets,"h_nJets","Number of tight jets",13,-0.5,12.5);
  TH1F *h_nbJets[maxSysts] ;systZero.initHistogramsSysts(h_nbJets,"h_nbJets","Number of tight b-jets",11,-0.5,10.5);
  TH1F *h_jet1Pt[maxSysts] ;systZero.initHistogramsSysts(h_jet1Pt,"h_jet1Pt","Leading jet Pt distribution",100,0,500);
  TH1F *h_jet2Pt[maxSysts] ;systZero.initHistogramsSysts(h_jet2Pt,"h_jet2Pt","Trailing Jet Pt distribution",100,0,500);
  TH1F *h_jet3Pt[maxSysts] ;systZero.initHistogramsSysts(h_jet3Pt,"h_jet3Pt","Third Jet Pt distribution",100,0,500);

  TH1F *h_bjet1Pt[maxSysts] ;systZero.initHistogramsSysts(h_bjet1Pt,"h_bjet1Pt","Leading b-jet Pt distribution",100,0,500);
  TH1F *h_bjet2Pt[maxSysts] ;systZero.initHistogramsSysts(h_bjet2Pt,"h_bjet2Pt","Trailing b-jet Pt distribution",100,0,500);
  TH1F *h_bjet3Pt[maxSysts] ;systZero.initHistogramsSysts(h_bjet3Pt,"h_bjet3Pt","Third b-jet Pt distribution",100,0,500);

  TH1F *h_bjetsPt[maxSysts] ;systZero.initHistogramsSysts(h_bjetsPt,"h_bjetsPt","B-Jets Pt distribution",100,0,500);

  for(Int_t i=0; i<nEvents; i++ ){
    if(i%100000==1 ){
    cout<<"Info: Running on event: "<<i<<endl; 
    }
  chain.GetEntry(i);

  systZero.setWeight(0,1.);
  systZero.setWeight("btagUp",1.);
  systZero.setWeight("btagDown",1.);
  systZero.setWeight("mistagUp",1.);
  systZero.setWeight("mistagDown",1.);
  systZero.setWeight("puDown",1.);
  systZero.setWeight("puUp",1.);
  systZero.setWeight("lepDown",1.);
  systZero.setWeight("lepUp",1.);


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


  int maxJetLoop = min(10, jet_size);
	    
  for (int j = 0; j <maxJetLoop;++j){	
      TLorentzVector jet;
      if(jet_iscsvm[j]<=0.){
          jet.SetPtEtaPhiE(jet_pt[j], jet_eta[j], jet_phi[j], jet_e[j]);
          jets_nob.push_back(jet);
          }
  
  jetsPhi.push_back(jet_phi[j]);
    
  if(jetsPhi.size()==1)systZero.fillHistogramsSysts(h_jet1Pt,jet_pt[j],1.);
  if(jetsPhi.size()==2)systZero.fillHistogramsSysts(h_jet2Pt,jet_pt[j],1.);
  if(jetsPhi.size()==3)systZero.fillHistogramsSysts(h_jet3Pt,jet_pt[j],1.);
    
  }//end of the jet loop
    
  fileout<<std::fixed<<std::setprecision(0)<<runNumber<<"   "<<evtNumber<<"   "<<lumiSec<<"   "<<std::endl;
  }//end of loop over events 
  
  fileout.close();  //return h
  
  //Write the Histogramms here  
  //systZero.writeHistogramsSysts(h_nJets, allMyFiles); 
  //systZero.writeHistogramsSysts(h_nbJets, allMyFiles); 

  systZero.writeHistogramsSysts(h_jet1Pt, allMyFiles); 
  systZero.writeHistogramsSysts(h_jet2Pt, allMyFiles); 
  systZero.writeHistogramsSysts(h_jet3Pt, allMyFiles); 

  //systZero.writeHistogramsSysts(h_bjet1Pt, allMyFiles); 
  //systZero.writeHistogramsSysts(h_bjet2Pt, allMyFiles); 
  //systZero.writeHistogramsSysts(h_bjet3Pt, allMyFiles); 
  
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
      
      cout << " writing histo "<< histo[sy+(MAXTOT+1)*c]->GetName()<< " in file "<< filesout[sy+(MAX+1)*c]->GetName()<<endl;;
      
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
      
      //fixes the error at the end 
      //filesout[(int)sy]->Close();
      //filesout[sy+(MAX+1)*c]->Close();
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
    //cout<< " in cat loop "<< c<<endl;
    //cout<< " value "<< wcats[c] <<endl;
    int MAX = this->maxSysts;
    bool useOnlyNominal = this->onlyNominal;
    for(int sy=0;sy<(int)MAX;++sy){
      //cout<< " in syst loop "<< sy<< endl;
      //cout<<" value "<< this->weightedSysts[(int)sy] <<endl ;
      if(sy!=0&& useOnlyNominal)continue;
      float ws = (this->weightedSysts[(int)sy])*wcats[c];
      //cout << " filling histogram "<< histo[(int)sy]->GetName() << " with value "<< v <<" and weight "<< w <<" ws "<< ws<<endl;
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

