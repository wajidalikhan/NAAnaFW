#!/usr/bin/env python
from ROOT import TFile, TCanvas, TColor, gStyle, TLegend, TLatex, TH1F, TTree, TH2F
from math import sqrt
from array import array

# style
gStyle.SetOptStat(0)
gStyle.SetLabelSize(0.06,"xy")
gStyle.SetTitleSize(0.06,"xy")
gStyle.SetTitleOffset(1.6,"x")
gStyle.SetTitleOffset(0.4,"y")
gStyle.SetPadTopMargin(0.1)
gStyle.SetPadRightMargin(0.1)
gStyle.SetPadBottomMargin(0.16)
gStyle.SetPadLeftMargin(0.14)

# create canvas
canvas1 = TCanvas("c1","c1", 800, 800)
canvas1.SetTickx()
canvas1.SetTicky()

variables=[]

# open file
muon = TFile("res/QCDMuEPt20toInf_muon.root")
muonantiiso = TFile("res/QCDMuEPt20toInf_muonantiiso.root")

for key in muon.GetListOfKeys():
    kname = key.GetName()
    #print kname
    variables.append(kname)

#muon = TFile("res/test_muon.root")
#muonantiiso = TFile("res/test_muonantiiso.root")

for var in variables:
    print "---- ", var
    
    # get histograms 
    hist_muon = muon.Get(var)
    hist_muon.Rebin(5)
    hist_muon.SetLineColor(TColor.GetColor(1))
    hist_muon.SetLineWidth(3)
    hist_muon.SetTitle('')
    hist_muon.GetYaxis().SetTitle('Normalized to Unity')
    hist_muon.GetYaxis().SetTitleOffset(1.6)
    hist_muon.GetXaxis().SetTitleOffset(1.3)
    hist_muon.GetXaxis().SetTitle(var)
    hist_muon.SetMaximum(hist_muon.GetMaximum()*1.5);
    hist_muon.DrawNormalized("Hist")
    
    hist_muonantiiso = muonantiiso.Get(var)
    hist_muonantiiso.Rebin(5)
    hist_muonantiiso.SetLineColor(TColor.GetColor(200))
    hist_muonantiiso.SetLineWidth(3)
    hist_muonantiiso.DrawNormalized("Histsame")
    
    # set legend for 1D
    l = TLegend(0.65, 0.65, 0.89, 0.85)
    l.SetTextSize(0.04)
    l.SetTextFont(42)
    l.SetFillColor(10)
    l.SetLineColor(10)
    l.SetBorderSize(0)
    l.AddEntry(hist_muon,"mu","l")
    l.AddEntry(hist_muonantiiso,"mu-antiiso","l")
    
    #save canvas
    canvas1.RedrawAxis()
    l.Draw("same")
    canvas1.Print(var+".png")

muon.Close()
muonantiiso.Close()


