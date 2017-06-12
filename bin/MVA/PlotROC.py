#!/usr/bin/env python
import ROOT
from ROOT import TFile, TCanvas, TColor, gStyle, TLegend, TLatex, TH1F, TTree, TH2F, gROOT
from math import sqrt
from array import array

label=['ST_vs_QCDMu','ST_vs_TT','ST_vs_VJ']
for i in label:
  print i

gStyle.SetOptStat(0)
gStyle.SetLabelSize(0.06,"xy")
gStyle.SetTitleSize(0.06,"xy")
gStyle.SetTitleOffset(1.2,"xy")
gStyle.SetTitleOffset(1.2,"xy")
gStyle.SetPadTopMargin(0.03)
gStyle.SetPadRightMargin(0.03)
gStyle.SetPadBottomMargin(0.08)
gStyle.SetPadLeftMargin(0.1)

file1 = TFile("tmvaST_vs_QCDMu.root")
file2 = TFile("tmvaST_vs_QCDMu_1.root")
file3 = TFile("tmvaST_vs_QCDMu_2.root")
file4 = TFile("tmvaST_vs_QCDMu_3.root")

canvas = TCanvas("c1","c1", 1000, 800)
canvas.SetTickx()
canvas.SetTicky()

h1 = file1.Get("Method_BDT/BDT/MVA_BDT_rejBvsS")
h2 = file2.Get("Method_BDT/BDT/MVA_BDT_rejBvsS")
h3 = file3.Get("Method_BDT/BDT/MVA_BDT_rejBvsS")
h4 = file4.Get("Method_BDT/BDT/MVA_BDT_rejBvsS")

h1.SetTitle('')
h1.GetYaxis().SetTitle('Background Rejection')
h1.GetYaxis().SetTitleOffset(1.4);
h1.GetXaxis().SetTitle('Signal Efficiency')

h1.SetLineColor(1)
h1.SetLineWidth(3)
h1.SetFillStyle(0)

h2.SetLineColor(2)
h2.SetLineWidth(3)
h2.SetFillStyle(0)

h3.SetLineColor(3)
h3.SetLineWidth(3)
h3.SetFillStyle(0)

h4.SetLineColor(4)
h4.SetLineWidth(3)
h4.SetFillStyle(0)

l = TLegend(0.15, 0.25, 0.35, 0.45)
l.SetTextSize(0.03)
l.SetTextFont(42)
l.SetFillColor(10)
l.SetLineColor(10)
l.SetBorderSize(0)

print "Integral of ROC = %4.2f"%(h1.Integral()/100);
s = "%8.2f" %h1.Integral()
l.AddEntry(h1,"MinNodeSize = 0.1 : ROC Integral = "+"%4.2f"%(h1.Integral()/100),"l")
l.AddEntry(h2,"MinNodeSize = 0.3 : ROC Integral = "+"%4.2f"%(h2.Integral()/100),"l")
l.AddEntry(h3,"MinNodeSize = 0.5 : ROC Integral = "+"%4.2f"%(h3.Integral()/100),"l")
l.AddEntry(h4,"MinNodeSize = 0.7 : ROC Integral = "+"%4.2f"%(h4.Integral()/100),"l")

h1.Draw("HIST")
h2.Draw("HISTsame")
h3.Draw("HISTsame")
h4.Draw("HISTsame")

l.Draw("same")



canvas.RedrawAxis()
canvas.Print("plots/ROC.png")
canvas.Print("plots/ROC.pdf")

file1.Close()


