#include "TFile.h" 
#include "TString.h" 
#include "TKey.h" 
#include <iostream>

void subtraction(
    const char *idirname="/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muonantiiso/histos_lin/", 
    //const char *odirname="/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muon/histos_lin/", 
    const char *odirname="./", 
    const char *type=".root"){
   
    TSystemDirectory idir(idirname, idirname);
    TSystemDirectory odir(odirname, odirname);
    TList *files = idir.GetListOfFiles();
    cout <<" In Path : "<<idirname<<endl;
    cout <<" Out Path : "<<odirname<<endl;
    
    TH1F *hsig, *hbkg, *htmp, *histsig, *histbkg, *histtmp;
    std::vector<TH1F*> listSig, listBkg, listTemp;
    std::vector<string> keys;

    if(files){
      TSystemFile *file;
      TString fname;
      TIter next(files);
      while ((file=(TSystemFile*)next())){
        fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(type)){

	      if (fname.Contains("DDQCD_muonantiiso")){
	        cout <<" File names : " <<fname.Data() << endl;
	        TFile f(idirname+fname);
	        for(auto && keyAsObj : *f.GetListOfKeys()){
	          auto key = (TKey*) keyAsObj;
              keys.push_back(key->GetName());
              histsig = (TH1F*) f.Get(key->GetName())->Clone(); 
              listSig.push_back((TH1F*)histsig);
              histsig->SetDirectory(0);

              histbkg = (TH1F*) f.Get(key->GetName())->Clone();    
              histbkg->Reset("ICES");
              listBkg.push_back((TH1F*)histbkg);
              histbkg->SetDirectory(0);           
          
              histtmp = (TH1F*) f.Get(key->GetName());    
              histtmp->Reset("ICES");
              histtmp->SetDirectory(0);           
          
              hsig = (TH1F*)f.Get("h_2j1t_mtw")->Clone();
              hbkg = (TH1F*)f.Get("h_2j1t_mtw")->Clone();
              hsig->SetDirectory(0);
	          hbkg->Reset("ICES");
	          hbkg->SetDirectory(0);
              }
            }
          }
        }
    next = TIter(files);
    while ((file=(TSystemFile*)next())){
      fname = file->GetName();
      if (!file->IsDirectory() && fname.EndsWith(type)){
        if (!fname.Contains("DDQCD_muonantiiso")&& !fname.Contains("Data_muonantiiso")){
          cout <<" File names : " <<fname.Data() << endl;
          Char_t filename[0];
          sprintf(filename,"%s%s",idirname,fname.Data());
          TFile g(filename);
          for(int i=0; i<keys.size();i++){
            TString name;
            name=keys.at(i);
            //cout <<filename<<":"<<name<<endl;
            TH1F *h =(TH1F*)g.Get(name);
            h->SetDirectory(0);
            listTemp.push_back(h);
            listSig.at(i)->Add(h,-1);
            listBkg.at(i)->Add(h);
            }
          htmp = (TH1F*)g.Get("h_2j1t_mtw")->Clone();
          cout <<" htmp : "<<htmp->Integral()<<endl;    
          htmp->SetDirectory(0);
          hsig->Add(htmp,-1);
          hbkg->Add(htmp);
          g.Close();
          }
        }
      }//while-loop
    }
 
 for(int i=0; i<listSig.size();i++){
      if(keys.at(i)=="h_2j0t_mtw" || keys.at(i)=="h_2j1t_mtw" || keys.at(i)=="h_3j1t_mtw" || keys.at(i)=="h_3j2t_mtw"){
      cout <<" key Name : "<<keys.at(i)<<" | Entries in MC-Subtracted QCD "<<listSig.at(i)->Integral()
      <<" | Total bkg MC "<<listBkg.at(i)->Integral()
      <<" | QCD Purity "<<listSig.at(i)->Integral()*100/(listSig.at(i)->Integral()+listBkg.at(i)->Integral())<<" %"<<endl;
      }
   }
    cout << "Entries in MC-Subtracted : "<<hsig->Integral()<<" Entries in MC : "<<hbkg->Integral()<< " QCD Purity : "<<hsig->Integral()*100/(hsig->Integral()+hbkg->Integral())<<"%"<<endl;
    
    Char_t outfilename[0];
    sprintf(outfilename,"%s%s",odirname,"DDQCD_muon.root");
    TFile outfile(outfilename,"RECREATE");
    for(int i=0; i<listSig.size();i++){
      listSig.at(i)->Write();
      }
    
    outfile.Close();

delete files;
delete hsig;
delete hbkg;
delete htmp;
delete histsig;
delete histbkg;
delete histtmp;
}
