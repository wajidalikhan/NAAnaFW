#! /usr/bin/env python

import os, multiprocessing
import copy
import math
from array import array

from utils import *

########## ######## ##########

import optparse
usage = "usage: %prog [options]"
parser = optparse.OptionParser(usage)
parser.add_option("-r", "--region", action="store", type="string", dest="region", default="fullhadronic")
parser.add_option("-c", "--channel", action="store", type="string", dest="channel", default="fullhadronic")
parser.add_option("-s", "--signal", action="store", type="string", dest="signal", default="")
parser.add_option("-e", "--pre", action="store_true", default=False, dest="pre")
parser.add_option("-o", "--post", action="store_true", default=False, dest="post")
parser.add_option("-a", "--all", action="store_true", default=False, dest="all")
parser.add_option("-b", "--batch", action="store_true", default=False, dest="batch")
parser.add_option("-B", "--blind", action="store_true", default=False, dest="blind")
(options, args) = parser.parse_args()
if options.batch: gROOT.SetBatch(True)

########## SAMPLES ##########
regions_list = {}
regions_list['fullhadronic'] = [ 'fullhadronic_CR_TT', 'fullhadronic_CR_VJets', 'fullhadronic_CR_WJets', 'fullhadronic']
regions_list['semileptonic'] = [ 'semileptonic_CR_TT', 'semileptonic_CR_WJets', 'semileptonic' ]

data = ["data_obs"]
back = {}
back['fullhadronic_CR_TT']    = ["otherBkg", "SingleTop", "WJets", "DY", "TT"]
back['fullhadronic_CR_VJets'] = ["otherBkg","QCD", "SingleTop", "WJets", "DY", "ZToNuNu", "TT"]
back['fullhadronic_CR_WJets'] = ["otherBkg","QCD", "SingleTop", "WJets", "DY", "ZToNuNu", "TT"]
back['fullhadronic']          = ["otherBkg","QCD", "SingleTop", "WJets", "DY", "ZToNuNu", "TT"]
back['semileptonic_CR_TT']    = ["otherBkg", "SingleTop", "WJets", "DY", "TT"]
back['semileptonic_CR_WJets'] = ["otherBkg", "SingleTop", "WJets", "DY", "ZToNuNu", "TT"]
back['semileptonic']          = ["otherBkg", "SingleTop", "WJets", "DY", "ZToNuNu", "TT"]

sign = ['DMtt_sc_Mchi1Mphi10']
exps = {'DMtt_sc_Mchi1Mphi10' : 1. }

#FILE = "mlfit_oldShort.root"
FILE = "mlfit.root"
#FILE = "mlfit_mia.root"
HFILEDIR = "histos/" #Data_fullhadronic.root

########## ######## ##########


def getData(dfile,channel,region):
#    print "channel",channel, "region", region
    hname = "metFinal"
    if channel == "fullhadronic":
        if region == "fullhadronic":
            hname = "metFinal"
        elif region == "fullhadronic_CR_TT":
            hname = "metFinal_SR_1lep"
#            hname = "metFinal_CR7"
#            hname = "metFinal_CR3nw"
        elif region == "fullhadronic_CR_VJets":
            hname = "metFinal_CR5"
        elif region == "fullhadronic_CR_WJets":
            hname = "metFinal_CR6nw"
#            hname = "metFinal_CR6"
    elif channel == "semileptonic":
        if region == "semileptonic":
            hname = "metFinal"
        elif region == "semileptonic_CR_TT":
            hname = "metFinal_2lep"
        elif region == "semileptonic_CR_WJets":
            hname = "metFinal_met_0btag"

    h = dfile.Get(hname)    
    return h

def plotfit(region):
    print 'plotting ', region
    if options.pre: hdir = 'shapes_prefit'
    elif options.post: hdir = 'shapes_fit_b'
    else:
        print "Please indicate if prefit (--pre) or postfit (--post)"
        exit()
    
    file = {}
    hist = {}
    tmphist = {}
    
    # get background
    err = array('d', [1.0])
    fileFit = TFile.Open(FILE, "READ")
    for i, s in enumerate(back[region]):
        hname = hdir + "/" + region + "/" + s
        tmphist[s] = fileFit.Get(hname)
    hname = hdir + "/" + region + "/total_background"
    tmphist['BkgSum'] = fileFit.Get(hname)
    [h.Scale(h.GetBinWidth(1)) for h in tmphist.itervalues()]

    # get data
    hfile = TFile.Open(HFILEDIR+'Data_'+options.channel+'.root',"READ")
    tmphist['data_obs'] = getData(hfile,options.channel,region)

    # Set the proper binning
    for s, h in tmphist.iteritems():
        hist[s] = tmphist['data_obs'].Clone(s)
        hist[s].Reset()
        for i in range(hist[s].GetNbinsX()+1):
            try:
                hist[s].SetBinContent(i, tmphist[s].GetBinContent(i))
                hist[s].SetBinError(i, tmphist[s].GetBinError(i))
            except:
                continue
    
    applyStyle(hist)

    # Fix fill
    fit = "prefit" if options.pre else "postfit" 
    out = draw(hist, fit, region, data, back[region], sign, 1, 1260, 4, False)
    
    if gROOT.IsBatch():
        outputname = ("prefit" if options.pre else "postfit")+"_"+region
        out[0].Print(outputname+".png")
        out[0].Print(outputname+".pdf")

    samples = [x for x in tmphist.keys() if not 'data' in x and not 'BkgSum' in x and not x in sign]

    print "-"*50
    print " --- Region:", region
    print "Sample                  Events"
    print "-"*50
    for i, s in enumerate(['data_obs']+samples+['BkgSum']):
        if i==1 or i==len(samples)+1: print "-"*40
        option = ""
        print "%-20s" % s, "\t%-10.2f +- %-10.2f" % (tmphist[s].IntegralAndError(1, tmphist[s].GetNbinsX(), err,option), err[0])
    print "-"*50
    for i, s in enumerate(sign):
        if not sample[s]['plot']: continue
        print "%-20s" % s, "\t%-10.2f +- %-10.2f" % (hist[s].IntegralAndError(1, hist[s].GetNbinsX(), err,option), err[0])
        pass
    print "-"*50
    
    if not gROOT.IsBatch(): raw_input("Press Enter to continue...")


def plotNevents(log=False):

    file = {}
    hist = {}
    tmphist = {} 
    err = array('d', [1.0])
    option = ""

    markers = {}
    markers['data']   = ROOT.kFullCircle
    markers['prefit'] = 0
    markers['fit_b']  = 0

    colors = {}
    colors['data']   = ROOT.kBlack
    colors['prefit'] = ROOT.kOrange-2
    colors['fit_b']  = ROOT.kBlue
    fill = {}
    fill['data']   = 0
    fill['prefit'] = 3001
    fill['fit_b']  = 0
    style = {}
    style['data']   = "pe1"
    style['prefit'] = "histo"
    style['fit_b']  = "histo"

    number_of_region = len(regions_list['fullhadronic']) + len(regions_list['semileptonic'])

    nevents_in_region = {}
    for i, s in enumerate(['data']+['prefit']+['fit_b']):
        nevents_in_region[s] = TH1F("nevents_in_region_"+s, "; regions; events",number_of_region,0.,number_of_region)
        nevents_in_region[s].GetYaxis().SetTitle("events")
        for j, r in enumerate(regions_list['fullhadronic']+regions_list['semileptonic']):
            nevents_in_region[s].GetXaxis().SetBinLabel(j+1,r)

    ibin=0
    for channel in (['fullhadronic']+['semileptonic']):
        for j, r in enumerate(regions_list[channel]):
            ibin=ibin+1
            region = r

            # get data
            hfile = TFile.Open(HFILEDIR+'Data_'+channel+'.root',"READ")
            h = getData(hfile,channel,region)
            nevents = h.IntegralAndError(1, h.GetNbinsX(), err,option)
            nevents_in_region['data'].SetBinContent(ibin,nevents)
            nevents_in_region['data'].SetBinError(ibin,err[0])

            # get background
            fileFit = TFile.Open(FILE, "READ")
            for fit in ['prefit','fit_b']:
                hdir = 'shapes_'+fit        
                hname = hdir + "/" + region + "/total_background"
                h = fileFit.Get(hname).Clone("bkg_"+fit)
                h.Scale(h.GetBinWidth(1))
                nevents = h.IntegralAndError(1, h.GetNbinsX(), err,option)
                nevents_in_region[fit].SetBinContent(ibin,nevents)
                nevents_in_region[fit].SetBinError(ibin,err[0])
            
    for i, s in enumerate(['prefit']+['fit_b']+['data']):
        nevents_in_region[s].SetMarkerStyle(markers[s])
        nevents_in_region[s].SetLineColor(colors[s])
        nevents_in_region[s].SetLineWidth(2)
        nevents_in_region[s].SetFillColor(colors[s])
        nevents_in_region[s].SetFillStyle(fill[s])
        setHistStyle(nevents_in_region[s],1.15)

    err = nevents_in_region['fit_b'].Clone("BkgErr;")
    err.SetTitle("")
    err.GetYaxis().SetTitle("data / bkg")
    err.SetFillColor(ROOT.kBlue)
    err.SetLineColor(ROOT.kBlue)
    err.SetFillStyle(3001)
    for i in range(1, err.GetNbinsX()+1):
        err.SetBinContent(i, 1)
        if nevents_in_region['fit_b'].GetBinContent(i) > 0:
            err.SetBinError(i, nevents_in_region['fit_b'].GetBinError(i)/nevents_in_region['fit_b'].GetBinContent(i))

    errLine = err.Clone("errLine")
    errLine.SetLineWidth(1)
    errLine.SetFillStyle(0)

    res = {}
    for s in (['fit_b']+['prefit']):
        res[s] = nevents_in_region['data'].Clone("residues")
        for i in range(0, res[s].GetNbinsX()+1):
            if nevents_in_region[s].GetBinContent(i) > 0: 
                res[s].SetBinContent(i, res[s].GetBinContent(i)/nevents_in_region[s].GetBinContent(i))
                res[s].SetBinError(i, res[s].GetBinError(i)/nevents_in_region[s].GetBinContent(i))
        res[s].SetLineColor(colors[s]+2)
        res[s].SetMarkerColor(colors[s]+2)

    # Legend
    n = 3
    leg = TLegend(0.54, 0.86-0.05*n, 0.8, 0.86)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetTextSize(0.04)
    leg.AddEntry(nevents_in_region['data'],  "data",                      "pe")
    leg.AddEntry(nevents_in_region['prefit'],"background (pre-fit)",      "f")
    leg.AddEntry(nevents_in_region['fit_b'], "background (b-only fit hp)","l")
        

    c1 = TCanvas("c1", "c1", 800, 800)
    gStyle.SetOptStat(0)
    c1.Divide(1, 2)
    setTopPad(c1.GetPad(1))
    setBotPad(c1.GetPad(2))
    c1.cd(1)
    c1.GetPad(True).SetTopMargin(0.08)
    c1.GetPad(True).SetRightMargin(0.05)
    c1.GetPad(True).SetTicks(1, 1)
    if log:
        c1.GetPad(True).SetLogy()
    
    for i, s in enumerate(['prefit']+['fit_b']+['data']):
        option = style[s] if i==0 else style[s]+"same"
        nevents_in_region[s].Draw(option)
    leg.Draw()

    lumi = 1260
    drawCMS(lumi, "Private", True)
        
    if gROOT.IsBatch():
        c1.Print("nevents_in_region.png")
        c1.Print("nevents_in_region.pdf")

    c1.cd(2)
    err.GetYaxis().SetRangeUser(0.5,1.5)
    setBotStyle(err,3,False)
    for s in (['fit_b']+['prefit']):
        setBotStyle(res[s])
    err.Draw("E2")
    errLine.Draw("SAME, HIST")
    for s in (['prefit']+['fit_b']):
        res[s].Draw("SAME, PE0")
    drawRatio(nevents_in_region['data'], nevents_in_region['prefit'],0.12,0.85,colors['prefit'])
    drawKolmogorov(nevents_in_region['data'], nevents_in_region['prefit'],0.73,0.85,colors['prefit'])
    drawRatio(nevents_in_region['data'], nevents_in_region['fit_b'],0.12,0.45,colors['fit_b'])
    drawKolmogorov(nevents_in_region['data'], nevents_in_region['fit_b'],0.73,0.45,colors['fit_b'])

    c1.Update()

    if gROOT.IsBatch():
        output_name="nevents_in_region" if not log else "nevents_in_region_log"
        c1.Print(output_name+".png")
        c1.Print(output_name+".pdf")


channel = options.channel
region = options.region

if options.all:
    regions = regions_list[channel]
    for r in regions :
        plotfit(r)
else:
    plotfit(region)

plotNevents(False)    
#plotNevents(True)    


