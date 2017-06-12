#include <algorithm>
#include <vector>
#include <iostream>     
#include <iomanip>      
TCanvas *c1 = new TCanvas("c1","Singletop",800,600);
void qcdFoM(){
  TString path = "/afs/cern.ch/work/w/wajid/NapoliFW/CMSSW_8_0_20/src/Analysis/NAAnaFW/python/MakePlot/output/muon2p4/histos_lin";
  TString var = "h_2j1t_mtwcut_sr_MuPt";
  
  vector<string> listOfVariable;
  listOfVariable.push_back("h_2j1t_mtwcut_sr_MuPt");
  listOfVariable.push_back("h_2j1t_mtwcut_sr_MuPt");
  listOfVariable.push_back("h_2j1t_mtwcut_sr_MuPt");
  listOfVariable.push_back("h_2j1t_mtwcut_sr_MuPt");
 
  for(auto i : listOfVariable) {
  // process i
     cout << i << " ";
  }
  gROOT->SetStyle("Plain");       // Switches off the ROOT default style
  gPad->UseCurrentStyle();        // this makes everything black and white,
  gROOT->ForceStyle();            // forces the style chosen above to be used,
                                  // not the style the rootfile was made with

  gStyle->SetCanvasColor(-1);
  gStyle->SetPadColor(-1);
  gStyle->SetFrameFillColor(-1);
  gStyle->SetHistFillColor(-1);
  gStyle->SetTitleFillColor(-1);
  gStyle->SetFillColor(-1);
  gStyle->SetFillStyle(4000);
  gStyle->SetStatStyle(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetCanvasBorderSize(0);
  gStyle->SetFrameBorderSize(0);
  gStyle->SetLegendBorderSize(0);
  gStyle->SetStatBorderSize(0);
  gStyle->SetTitleBorderSize(0);
 
  
  TH1::SetDefaultSumw2(true);
  gStyle->SetOptStat(0);
  c1->Range(-1.335008,-1154.067,1.331658,9580.143);
  c1->SetFillColor(0);
  c1->SetBorderMode(0);
  c1->SetBorderSize(2);
  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->SetLeftMargin(0.15);
  c1->SetFrameBorderMode(0);
  c1->SetFrameBorderMode(0);

  TFile *sig0 = TFile::Open(path+"/"+"ST_tch_muon.root");
  TFile *bkg1 = TFile::Open(path+"/"+"QCDMuPt20toInf_muon.root");
  TFile *bkg2 = TFile::Open(path+"/"+"DYJets_muon.root");
  TFile *bkg3 = TFile::Open(path+"/"+"ST_sch_muon.root");
  TFile *bkg4 = TFile::Open(path+"/"+"ST_tW_muon.root");
  TFile *bkg5 = TFile::Open(path+"/"+"TT_muon.root");
  TFile *bkg6 = TFile::Open(path+"/"+"VV_muon.root");
  TFile *bkg7 = TFile::Open(path+"/"+"WJets_muon.root");

  TH1F *h0  = (TH1F*)sig0-> Get(var);
  TH1F *h1  = (TH1F*)bkg1-> Get(var);
  TH1F *h2  = (TH1F*)bkg2-> Get(var);
  TH1F *h3  = (TH1F*)bkg3-> Get(var);
  TH1F *h4  = (TH1F*)bkg4-> Get(var);
  TH1F *h5  = (TH1F*)bkg5-> Get(var);
  TH1F *h6  = (TH1F*)bkg6-> Get(var);
  TH1F *h7  = (TH1F*)bkg7-> Get(var);

  TH1F *stack = (TH1F*)h1->Clone();
  stack->Add(h2); 
  stack->Add(h3); 
  stack->Add(h4); 
  stack->Add(h5); 
  stack->Add(h6); 
  stack->Add(h7); 

  cout <<"Total Signal : "<<h0->Integral()<<" Total Bkg : "<<stack->Integral()<<endl; 
  cout <<" Get Bins = "<< stack->GetNbinsX()<<endl;
  //int nbins = stack->GetNbinsX();
  int nbins = 25; 

  vector<double> x,x2,x3,x4;
  vector<double> y,y2,y3,y4;
  vector<float_t> cutvalue;
  vector<float_t> cutvalue_bkg;
  vector<float_t> max;

  for (int i =1;i<=nbins;i++){
    cout <<" "+var+" value = "<< stack->GetBinLowEdge(i) << endl;
    
    double Bkg =    h1->Integral(i,nbins) 
                + h2->Integral(i,nbins)
                + h3->Integral(i,nbins)
                + h4->Integral(i,nbins)
                + h5->Integral(i,nbins)
                + h6->Integral(i,nbins)
                + h7->Integral(i,nbins);
 
    double Sig = h0->Integral(i,nbins);
    if((stack->GetBinContent(i)<=0) &&  (h0->GetBinContent(i)<=0))continue;
    double fog = Sig/TMath::Sqrt(Sig+Bkg);
    double fog_bkgunc2 = Sig/TMath::Sqrt(Sig+Bkg+(pow(0.01*Bkg,2)));
    double fog_bkgunc3 = Sig/TMath::Sqrt((pow(0.01*Bkg,2)));
    double fog_bkgunc4 = Sig/TMath::Sqrt((Bkg));
    
    cutvalue.push_back(fog);
    cutvalue_bkg.push_back(fog_bkgunc2);
    max.push_back(fog_bkgunc3);
    

    x.push_back(stack->GetBinLowEdge(i));    
    y.push_back(fog);

    x2.push_back(stack->GetBinLowEdge(i));    
    y2.push_back(fog_bkgunc2);
    
    x3.push_back(stack->GetBinLowEdge(i));    
    y3.push_back(fog_bkgunc3);
    
    x4.push_back(stack->GetBinLowEdge(i));    
    y4.push_back(fog_bkgunc4);
    
    cout <<" Bin["<<i<<"]"<<endl;
    cout <<" Sig = "<<Sig<<endl;
    cout <<" Bkg = "<<Bkg<<endl;
    cout <<" S/B = "<<Sig/Bkg<<endl;
    cout <<" S/SQRT(S+B) = "<<fog<<endl;
    cout <<" S/SQRT(S+B+deltaB^2) = "<<fog_bkgunc2<<endl;
    cout <<" S/SQRT(deltaB^2) = "<<fog_bkgunc3<<endl;
    cout <<" S/SQRT(B) = "<<fog_bkgunc4<<endl;
    cout <<endl;
    }
 
    cout <<endl;
    cout <<" Maximum s/sqrt(s+b) = "<<*max_element(cutvalue.begin(), cutvalue.end())<<" | Minimum s/sqrt(s+b) = " << *min_element(cutvalue.begin(), cutvalue.end())<<endl;
    cout <<" Maximum s/sqrt(s+b+dB2) = "<<*max_element(cutvalue_bkg.begin(),cutvalue_bkg.end())<<" | Minimum s/sqrt(s+b+dB2) = " <<*min_element(cutvalue_bkg.begin(),cutvalue_bkg.end())<<endl;
    
    TGraph *gr = new TGraph(x.size(), &(x[0]), &(y[0]));
    gr->SetTitle("");
    gr->SetName("S/#sqrt{S+B}");
    gr->SetFillStyle(0);
    gr->SetLineColor(2);
    gr->SetLineWidth(3);

    gr->SetMarkerColor(kBlack);
    gr->SetMarkerStyle(0); 
    gr->GetXaxis()->SetTitle(var);
    gr->GetXaxis()->SetLabelFont(42);
    gr->GetXaxis()->SetLabelSize(0.05);
    gr->GetXaxis()->SetTitleSize(0.05);
    gr->GetXaxis()->SetTitleFont(42);
    gr->GetYaxis()->SetTitle("Arbitrary Units");
    gr->GetYaxis()->SetLabelFont(42);
    gr->GetYaxis()->SetLabelOffset(0.003);
    gr->GetYaxis()->SetLabelSize(0.05);
    gr->GetYaxis()->SetTitleSize(0.05);
    gr->GetYaxis()->SetTitleOffset(1.4);
    gr->GetYaxis()->SetTitleFont(42);
    if(var.Contains("h_2j1t_mtwcut_sr_MuPt") || var.Contains("h_3j2t_mtwcut_MuPt") || var.Contains("h_3j2t_MuPt"))gr->SetMaximum(*max_element(max.begin(), max.end())*1.4);
    else gr->SetMaximum(*max_element(cutvalue.begin(), cutvalue.end())*1.3);

    gr->Draw("AC");

    TGraph *gr2 = new TGraph(x2.size(), &(x2[0]), &(y2[0]));
    gr2->SetTitle("");
    gr2->SetName("S/#sqrt{S+B+#deltaB^{2}}");
    gr2->SetFillStyle(0);
    gr2->SetLineColor(3); 
    gr2->SetLineWidth(3);
    gr2->Draw("SAME");
 

    TGraph *gr3 = new TGraph(x3.size(), &(x3[0]), &(y3[0]));
    gr3->SetTitle("");
    gr3->SetName("S/#sqrt{#deltaB^{2}}");
    gr3->SetFillStyle(0);
    gr3->SetLineColor(1); 
    gr3->SetLineWidth(3);
    gr3->Draw("SAME");
    
    TGraph *gr4 = new TGraph(x4.size(), &(x4[0]), &(y4[0]));
    gr4->SetTitle("");
    gr4->SetName("S/#sqrt{B}");
    gr4->SetFillStyle(0);
    gr4->SetLineColor(4); 
    gr4->SetLineWidth(3);
    gr4->Draw("SAME");
    TLegend *leg= c1->BuildLegend();
    
    leg->SetFillColor(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.05);
    leg->SetLineColor(0);
    leg->SetLineStyle(1);
    leg->SetLineWidth(1);
    
    TLine *vline = new TLine(*max_element(cutvalue.begin(), cutvalue.end()),gr->GetMinimum(),*max_element(cutvalue.begin(), cutvalue.end()),gr->GetMaximum()); 
    vline->SetLineColor(kBlue);
    vline->SetLineWidth(2);
    vline->SetLineStyle(7); 
    //vline->Draw("SAME"); 
    
    TLine *hline = new TLine(gr->GetXaxis()->GetXmin(),*max_element(cutvalue.begin(), cutvalue.end()),gr->GetXaxis()->GetXmax(),*max_element(cutvalue.begin(), cutvalue.end()));
    hline->SetLineColor(kBlue);
    hline->SetLineWidth(2);
    hline->SetLineStyle(7);
    //hline->Draw("SAME");

    c1->SaveAs("SoverB_"+var+".png");

}
