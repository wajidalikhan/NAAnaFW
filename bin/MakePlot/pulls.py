#! /usr/bin/env python

import os, multiprocessing
import copy
import math
from array import array
from ROOT import ROOT, gROOT, gStyle, gRandom, TSystemDirectory
from ROOT import TFile, TChain, TTree, TCut, TH1F, TH2F, THStack, TGraph
from ROOT import TStyle, TCanvas, TPad
from ROOT import TLegend, TLatex, TText, TLine

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-f", "--fileName", action="store", type="string", dest="fileName", default="")
parser.add_option("-b", "--binbybin",    action="store_true", dest="binbybin", default=False)
parser.add_option("-n", "--no-binbybin", action="store_false",dest="binbybin", default=False)
(options, args) = parser.parse_args()

fileName = options.fileName
binbybin = options.binbybin
# Format:
#0 : sys name
#1 : b-only DX
#2 : b-only dDX
#3 : s+b DX
#4 : s+b dDX
#5 : rho

gStyle.SetOptStat(0)


def pulls(filename):
    
    with open(filename) as f:
        fullcontent = f.readlines()
    
#    content_nostat = [x for x in fullcontent if not 'CMS_stat' in x]
    content_nostat = [x for x in fullcontent if not 'mcstat_' in x]

    content = fullcontent if binbybin else content_nostat

    nbins = len(content)
    b_pulls = TH1F("b_pulls", ";;pulls", nbins, 0., nbins)
    s_pulls = TH1F("s_pulls", ";;pulls", nbins, 0., nbins)
    rho     = TH1F("rho",     ";;rho",   nbins, 0., nbins)
    
    for i, s in enumerate(content):
        l = s.split()

        b_pulls.GetXaxis().SetBinLabel(i+1, l[0])
        s_pulls.GetXaxis().SetBinLabel(i+1, l[0])
        rho.GetXaxis().SetBinLabel(    i+1, l[0])

        b_pulls.SetBinContent(i+1, float(l[1]))
        b_pulls.SetBinError  (i+1, float(l[2]))
        s_pulls.SetBinContent(i+1, float(l[3]))
        s_pulls.SetBinError  (i+1, float(l[4]))
        rho.SetBinContent    (i+1, float(l[5]))
    
    b_pulls.SetFillStyle(3003)
    b_pulls.SetFillColor(ROOT.kBlue+2)
    b_pulls.SetLineColor(0)
    b_pulls.SetMarkerColor(ROOT.kBlue+2)
    b_pulls.SetLineWidth(0)
    b_pulls.SetMarkerStyle(20)
    b_pulls.SetMarkerSize(1.2)

    s_pulls.SetLineColor(ROOT.kRed+2)
    s_pulls.SetMarkerColor(ROOT.kRed+2)
    s_pulls.SetMarkerStyle(24)
    s_pulls.SetLineWidth(2)

    rho.SetLineColor(ROOT.kRed+2)
    rho.SetMarkerColor(ROOT.kRed+2)
    rho.SetMarkerStyle(24)
    rho.SetLineWidth(2)

    b_pulls.GetYaxis().SetRangeUser(-3., 3.)
    s_pulls.GetYaxis().SetRangeUser(-3., 3.)
    rho.GetYaxis().SetRangeUser(-1.1, 1.1)
    
    
    c1 = TCanvas("c1", "pulls", 1600, 800)
    c1.cd()
    c1.SetTopMargin(0.06)
    c1.SetRightMargin(0.05)
    c1.SetTicks(1, 1)
    
    # Draw
    b_pulls.Draw("E2")
    s_pulls.Draw("SAME, PE1")
    
    leg = TLegend(0.35, 0.95, 0.65, 0.995)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetNColumns(2)
    leg.AddEntry(b_pulls, "bkg-only fit", "fp")
    leg.AddEntry(s_pulls, "sgn+bkg fit",  "lp")
    
    line0 = TLine()
    line0.SetLineColor(ROOT.kGray+2)
    line0.SetLineStyle(6)
    line0.SetLineWidth(2)
    line0.DrawLine(0.,  0., nbins,  0.)

    line1 = TLine()
    line1.SetLineColor(ROOT.kGreen+2)
    line1.SetLineStyle(6)
    line1.SetLineWidth(2)
    line1.DrawLine(0.,  1., nbins,  1.)
    line1.DrawLine(0., -1., nbins, -1.)
    
    line2 = TLine()
    line2.SetLineColor(ROOT.kYellow+2)
    line2.SetLineStyle(6)
    line2.SetLineWidth(2)
    line2.DrawLine(0.,  2., nbins,  2.)
    line2.DrawLine(0., -2., nbins, -2.)
    
    leg.Draw()

    suffix = "all" if binbybin else "no_binbybin"
    print suffix, binbybin
    outputname = "pulls_"+ suffix
    c1.Print(outputname+".png")
    c1.Print(outputname+".pdf")
    
    c2 = TCanvas("c2", "rho", 1600, 800)
    c2.cd()
    c2.SetTopMargin(0.06)
    c2.SetRightMargin(0.05)
    c2.SetTicks(1, 1)
    
    # Draw
    rho.Draw("E0")
    
    leg2 = TLegend(0.35, 0.95, 0.65, 0.995)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetFillColor(0)
    leg2.AddEntry(rho, "rho", "flp")
    
    outputname = "rho_"+ suffix
    c2.Print(outputname+".png")
    c2.Print(outputname+".pdf")
    
    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")

pulls(fileName)
