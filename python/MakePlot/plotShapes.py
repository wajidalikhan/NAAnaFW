#!/usr/bin/env python
import ROOT
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
ddqcdmcsub = TFile("DDQCD-Purity/DDQCD_muon.root")
ddqcd = TFile("output/muonantiiso/histos_lin/DDQCD_muonantiiso.root")

for key in ddqcd.GetListOfKeys():
    kname = key.GetName()
    #print kname
    #variables.append(kname)
    variables.append('h_2j1t_mtw')


# get histograms 
hist_ddqcd = ddqcd.Get('h_2j1t_mtw')
#hist_ddqcd.Rebin(5)
#hist_ddqcd.SetLineColor(TColor.GetColor(1))
hist_ddqcd.SetLineColor(1)
hist_ddqcd.SetFillColor(0)
hist_ddqcd.SetLineWidth(3)
hist_ddqcd.SetTitle('')
hist_ddqcd.GetYaxis().SetTitle('Normalized to Unity')
hist_ddqcd.GetYaxis().SetTitleOffset(1.6)
hist_ddqcd.GetXaxis().SetTitleOffset(1.3)
#hist_ddqcd.GetXaxis().SetTitle(var)
hist_ddqcd.GetXaxis().SetTitle('MTW [GeV]')
hist_ddqcd.SetMaximum(hist_ddqcd.GetMaximum()*1.5);
hist_ddqcd.DrawNormalized("Hist")

hist_ddqcdmcsub = ddqcdmcsub.Get('h_2j1t_mtw')
#hist_ddqcdmcsub.Rebin(5)
hist_ddqcdmcsub.SetLineColor(2)
hist_ddqcdmcsub.SetFillColor(0)
hist_ddqcdmcsub.SetLineWidth(3)
hist_ddqcdmcsub.DrawNormalized("Histsame")

# set legend for 1D
l = TLegend(0.4, 0.6, 0.8, 0.8)
l.SetTextSize(0.04)
l.SetTextFont(42)
l.SetFillColor(10)
l.SetLineColor(10)
l.SetBorderSize(0)
l.AddEntry(hist_ddqcd,"DDQCD","l")
l.AddEntry(hist_ddqcdmcsub,"DDQCD MC SUB","l")

#save canvas
canvas1.RedrawAxis()
l.Draw("same")
canvas1.Print("h_2j1t_mtw.png")

ddqcd.Close()
ddqcdmcsub.Close()


