{
  string basestring = "/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin/";
  string endstring  = "/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/bin/pu/";
  string txtlocation = "files/final/ST_T_tch";

  TChain *t=new TChain("DMTreesDumper/WeightHistory");
  //TChain *t=new TChain("DMTreesDumper/ttDM__noSyst");
    
  TH1D  *pileup = new TH1D("pileup", "pileup", 74, -0.0, 74.);
  TString rootFile;
  ifstream is((basestring + txtlocation + string(".txt")).c_str());
    
  while(rootFile.ReadLine(is) && (!rootFile.IsNull())){
    printf(" Adding file: %s\n", rootFile.Data());
    string currfile = rootFile.Data();
    currfile = t->AddFile(currfile.c_str());
    cout << "entries " << t->GetEntries() << endl;
  }
  
  t->Draw("Event_nTruePV>>pileup", "");

  string outname = endstring + "MCPU.root";
  TFile *outfile = new TFile(outname.c_str(), "RECREATE");
  pileup->Write();
  outfile->Close();
  // cleanup
  delete pileup;
  delete outfile;
  delete t;
  cout << endl << "finished " << endl;
}
