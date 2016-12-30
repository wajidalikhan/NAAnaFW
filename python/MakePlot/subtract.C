#include <vector>
#include <iostream>
#include <vector>

    TString ipath="/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin";
    TString opath="/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin";
    
    TFile *f1 = TFile::Open("/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin/DD-QCD_muonantiiso.root");
    TFile *f2 = TFile::Open("/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin/ST_tch_muonantiiso.root");
    TFile *f3 = TFile::Open("/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin/ST_sch_muonantiiso.root");
    TFile *f4 = TFile::Open("/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin/ST_tW_muonantiiso.root");
    TFile *f5 = TFile::Open("/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin/TT_muonantiiso.root");
    TFile *f6 = TFile::Open("/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin/WJets_muonantiiso.root");
    TFile *f7 = TFile::Open("/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin/DYJets_muonantiiso.root");
    TFile *f8 = TFile::Open("/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin/VV_muonantiiso.root");

void getobservable(TString obs){
    
    TH1D *h1 = (TH1D*) f1->FindObjectAny(obs);
    TH1D *h2 = (TH1D*) f2->FindObjectAny(obs);
    TH1D *h3 = (TH1D*) f3->FindObjectAny(obs);
    TH1D *h4 = (TH1D*) f4->FindObjectAny(obs);
    TH1D *h5 = (TH1D*) f5->FindObjectAny(obs);
    TH1D *h6 = (TH1D*) f6->FindObjectAny(obs);
    TH1D *h7 = (TH1D*) f7->FindObjectAny(obs);
    TH1D *h8 = (TH1D*) f8->FindObjectAny(obs);
    
    h2->Add(h3);
    h2->Add(h4);
    h2->Add(h5);
    h2->Add(h6);
    h2->Add(h7);
    h2->Add(h8);
   
    //Subtracting all the Non-QCD from the QCD 
    
    if(obs=="h_2j0t_mtw"|| obs=="h_2j1t_mtw" || obs=="h_3j1t_mtw" || obs=="h_3j2t_mtw"){
    cout <<"Entries in QCD = "<<h1->Integral()<<" Entries in non-QCD = "<<h2->Integral()<<" Purity = "<<((h1->Integral()-h2->Integral())/h1->Integral())*100<<"%"<<endl;
    }
    h1->Add(h2,-1);
    
    h1->Write();

}
void subtract(){
    vector<TString> variables;
    TFile *file = TFile::Open("/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin/DD-QCD_muonantiiso.root");
    for (auto&& keyAsObj : *file->GetListOfKeys()){
        auto key = (TKey*) keyAsObj;
        //cout << key->GetName() << " " << key->GetClassName() << endl;
        variables.push_back(key->GetName());
        }
    
    //TFile *fileout = new TFile("/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muon/histos_lin/DD-QCD_muon.root","RECREATE");
    TFile *fileout = new TFile("/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin/DD-QCD_muon.root","RECREATE");
    fileout->cd();    
   
    for (int i; i<variables.size();i++){
      cout<<"Variable written : "<<variables[i]<<endl;
      getobservable(variables[i]);
      }
     
    fileout->Close();

}


