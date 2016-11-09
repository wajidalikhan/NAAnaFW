#!/bin/env python


########################################
###  Author: Annapaola de Cosa             
###  Code to produce pre e post fit plots  
###  from combine output (mlfit.root)
########################################


import sys
import os
import ROOT

from plots.services import Histo, Stack, Legend
from samples.toPlot import samples

import tdrstyle, CMS_lumi



ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)
from ROOT import SetOwnership


    
def raiseExep(s):
    raise RuntimeError('Failed to retrieve %s' %s )

def getHisto(samp):
    h = ROOT.gDirectory.Get(samp)
    if not isinstance(h, ROOT.TH1):
        raise RuntimeError('Failed to retrieve %s' %samp )
    return h

                  
def getSamples(fitPhase, region, samps):
    ifilename= "mlfit.root"
    ifile = ROOT.TFile.Open(ifilename)
    #print "-----> Opening input file ", ifilename
    ROOT.gDirectory.cd("shapes_"+fitPhase+"/"+region)
    #print "--- Looking into directory: " + "shapes_"+fitPhase+"/"+region
    histos = { s: getHisto(s) for s in samps}
    [h.Scale(h.GetBinWidth(1)) for h in histos.itervalues()]
    #print "Histo type after getting: ", histos["TT"]
    return histos

def getData(channel, var):
    ifilename= "histos/Data_"+channel+".root"
    print " file name data ", ifilename
    ifile = ROOT.TFile.Open(ifilename)
    h = ROOT.gDirectory.Get(var)
    if not isinstance(h, ROOT.TH1):
        raise RuntimeError('Failed to retrieve Data')
    h_data = Histo.fromTH1(h)
    h_data.SetStyle(samples["Data"].color, samples["Data"].style, samples["Data"].fill)
    return h_data
    
def setSampStyle(histos):
    histos_new = { s: Histo.fromTH1(h) for (s,h) in histos.iteritems()}
    [ h.SetStyle(samples[s].color, samples[s].style, samples[s].fill) if s in samples.keys() else  raiseExep(s) for (s,h) in histos_new.iteritems() ]
    return histos_new

    
def doStack(histos):
    stack = Stack("met","met")
    leg = Legend()
    for s,h in histos.iteritems():
        stack.Add(h)
        leg.AddEntry(h, samples[s].leglabel, "f")

    return stack,leg

def getStack(fitPhase, region, samps):
    histos =  setSampStyle(getSamples(fitPhase, region, samps))
    print "+ TT integral prefit: ", histos["TT"].Integral("width")
    stack,leg = doStack(histos)
    return stack,leg


def plot(fitPhase, region, samps, lumi = "1.26"):

    ### Setting canvas and pad style
    tdrstyle.setTDRStyle();
    H=600
    W=700
    #T = 0.08*H
    #B = 0.12*H
    L = 0.12*W
    R = 0.08*W
    c1 = ROOT.TCanvas("c1","c1",50,50,H,W)
    #c1 = ROOT.TCanvas()
    c1.SetFillColor(0)
    c1.SetBorderMode(0)
    c1.SetFrameFillStyle(0)
    c1.SetFrameBorderMode(0)
    c1.SetLeftMargin( L/W )
    c1.SetRightMargin( R/W )
    #c1.SetTopMargin( T/H )
    #c1.SetBottomMargin(B/H)
    c1.SetTickx(1)
    c1.SetTicky(1)
    c1.cd()
    unit = "fb^{-1}"
    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_sqrtS = str(lumi)+ unit +" (13 TeV)"
    iPeriod = 0
    iPos = 11
    # writing the lumi information and the CMS "logo"

    pad= ROOT.TPad("pad", "pad", 0, 0 , 1, 1)
    SetOwnership(pad, 0)
    pad.SetBorderMode(0);
    pad.SetTickx(0);
    pad.SetTicky(0);
    pad.Draw()
    pad.cd()

    ### Creating and drawing TStack and Legend
      
    stack,leg = getStack(fitPhase, region[1], samps)
    
    print "region " , region[0]

    stack.SetMaximum(stack.GetMaximum()*1.2)
    stack.DrawStack()
    tmp = ROOT.TH1F()
    h_data = Histo.fromTH1(tmp)
    if(region[1].startswith("semilep")):
        h_data=getData("semileptonic",region[0])
    if(region[1].startswith("fullhad")):
        h_data=getData("fullhadronic",region[0])
    h_data.Draw("same")
    leg.AddEntry(h_data, samples["Data"].leglabel, "lp")
    leg.Draw("same")
    
    CMS_lumi.CMS_lumi(pad, iPeriod, iPos)
    
    c1.cd()
    c1.Update()
    c1.RedrawAxis()
    path = "./pdf"
    if not os.path.isdir(path):
        os.makedirs(path)
    pdfname = "pdf/"+ fitPhase+"_" + region[1]+".pdf"   
    c1.Print(pdfname)
  



    
    




