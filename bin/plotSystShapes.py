#!/usr/bin/env python
import ROOT
from ROOT import TFile, TCanvas, TColor, gStyle, TLegend, TLatex, TH1F
from math import sqrt
from array import array
import os, optparse, os.path, optparse, subprocess, sys, glob
from os.path import join, exists
ROOT.gROOT.SetBatch(True)

usage = 'python plotSystShapes.py -o /home/student -c muon -S ST_tch,TT,ST_tW,ST_sch,VV,WJets,QCDMuPt20toInf,DYJets -s pu,jes,btag,mistag,lep'
parser = optparse.OptionParser(usage)

parser.add_option('-v', '--variable',  dest='variable', type='string', default = 'none', help ='list of input variable')
parser.add_option('-o', '--outdir',  dest='outdir', type='string', default = 'plots', help ='output for pdf, png')
parser.add_option('-c', '--channel', dest='channel', type='string', default = 'muon', help ='physics channel')
parser.add_option('-s', '--syst', dest='syst', type='string', default = 'none', help ='List of systematics comma seperated pu,jes,jer PS space between the commas are not supported at the moment')
parser.add_option('-S', '--sample', dest='sample', type='string', default = 'none', help ='List of samples comma seperated ST_tch,TT,ST_sch PS space between the commas are not supported at the moment')

(opt, args) = parser.parse_args()

if opt.channel not in ['muon','electron']:
    parser.error('Channel available too choose: muon, electron')


# create a output directory if it dosen't exists
ndir = 'plots'
if opt.outdir:
    if not exists(opt.outdir + '/' + ndir):
        os.makedirs(opt.outdir + '/' + ndir)
else:
    if not exists(ndir):
        os.makedirs(ndir)
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

def makeList(inlist):
    tmplist = []
    systs = inlist.split(",")
    print systs
    for s in systs:
        tmplist.append(s)
        if s.startswith("_"):
            tmplist.append(s[1:])
    return tmplist

variables=[]
samples = []
systematic = []
systematic = makeList(opt.syst)
samples = makeList(opt.sample)

for s in samples:
    # open file
    nominal = TFile("histos_lin/"+s+"_"+opt.channel+".root")
    for key in nominal.GetListOfKeys():
        kname = key.GetName()
        print kname
        variables.append(kname)

for var in variables:
    if var == 'h_2j1t_mtw':
        for s in samples:
            f_nominal = TFile("histos_lin/"+s+"_"+opt.channel+".root")
            for sy in systematic:
                f_up = TFile("histos_lin/"+s+"_"+opt.channel+"_"+sy+"Up.root")
                f_down = TFile("histos_lin/"+s+"_"+opt.channel+"_"+sy+"Down.root")
                print "Info: sample + systematic : " + s + "_" + sy
                hist = f_nominal.Get(var)
                hist.SetLineColor(2)
                hist.SetMarkerSize(0.9)
                hist.SetMarkerColor(2)
                hist.SetLineWidth(3)
                hist.SetTitle('')
                hist.GetYaxis().SetTitle('norm. to unity')
                hist.GetYaxis().SetTitleOffset(1.6)
                hist.GetXaxis().SetTitleOffset(1.3)
                hist.DrawNormalized()


                # taking the up-variation
                histUp = f_up.Get(var)
                histUp.SetLineColor(4)
                histUp.SetLineWidth(3)
                histUp.SetMarkerSize(0.9)
                histUp.SetMarkerColor(4)
                histUp.SetTitle('')
                histUp.DrawNormalized("SAME")

                # taking the down-variation
                histDown = f_down.Get(var)
                histDown.SetLineColor(3)
                histDown.SetLineWidth(3)
                histDown.SetMarkerSize(0.9)
                histDown.SetMarkerColor(3)
                histDown.SetTitle('')
                histDown.DrawNormalized("SAME")

                tex1 = TLatex(0.13,0.903,"#bf{CMS} Preliminary");
                tex1.SetNDC();
                tex1.SetTextAngle(0);
                tex1.SetTextFont(42);
                tex1.SetTextSize(0.04);
                tex1.SetTextAlign(11);
                #tex1.Draw();

                if s == 'QCDMuPt20toInf':
                    tex2 = TLatex(0.15,0.905,"Sample: "+s+ " Syst: "+sy+ "13 TeV")
                else:
                    tex2 = TLatex(0.2,0.905,"Sample: "+s+ " Syst: "+sy+ "13 TeV")
                tex2.SetNDC();
                tex2.SetTextAngle(0);
                tex2.SetTextFont(42);
                #tex2.SetTextAlign(11);
                tex2.SetTextSize(0.04);
                tex2.Draw();


                # set legend for 1D
                l = TLegend(0.63, 0.7, 0.93, 0.85)
                l.SetTextSize(0.04)
                l.SetTextFont(42)
                l.SetFillColor(10)
                l.SetLineColor(10)
                l.SetBorderSize(0)

                #l.AddEntry(hist,s+" "+sy,"p")
                l.AddEntry(hist,"nominal","l")
                l.AddEntry(histUp,"up","l")
                l.AddEntry(histDown,"down","l")
                l.Draw("SAME")
                canvas1.RedrawAxis()
                canvas1.Print(opt.outdir + '/' + ndir +'/' +s+"_"+var+"_"+sy+".png")
                canvas1.Print(opt.outdir + '/' + ndir +'/' +s+"_"+var+"_"+sy+".pdf")

        f_up.Close()
        f_nominal.Close()
        f_down.Close()
