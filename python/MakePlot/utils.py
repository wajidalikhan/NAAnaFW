import os, multiprocessing
import copy
import math
from array import array
from ROOT import ROOT, gROOT, gStyle, gRandom, TSystemDirectory
from ROOT import TFile, TChain, TTree, TCut, TH1F, TH2F, THStack, TGraph
from ROOT import TStyle, TCanvas, TPad
from ROOT import TLegend, TLatex, TText, TLine

from samples import *

def draw(hist, fit, channel, data, back, sign, snorm=1, lumi=-1, ratio=0, log=False):
    # If not present, create BkgSum
    if not 'BkgSum' in hist.keys():
        hist['BkgSum'] = hist['data_obs'].Clone("BkgSum") if 'data_obs' in hist else hist[back[0]].Clone("BkgSum")
        hist['BkgSum'].Reset("MICES")
        for i, s in enumerate(back): hist['BkgSum'].Add(hist[s])
    hist['BkgSum'].SetMarkerStyle(0)
    
    setHistStyle(hist['BkgSum'], 1.1 if ratio else 1.0)
    # Create stack
    bkg = THStack("Bkg", ";"+hist['BkgSum'].GetXaxis().GetTitle()+";"+hist['BkgSum'].GetYaxis().GetTitle())
    for i, s in enumerate(back): bkg.Add(hist[s])
    
    # Legend
    n = len([x for x in data+back+['BkgSum']+sign if sample[x]['plot']])
    leg = TLegend(0.69, 0.86-0.04*n, 0.95, 0.86)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    if len(data) > 0:
        leg.AddEntry(hist[data[0]], sample[data[0]]['label'], "pl")
    for i, s in reversed(list(enumerate(['BkgSum']+back))):
        leg.AddEntry(hist[s], sample[s]['label'], "f")
    for i, s in enumerate(sign):
        if sample[s]['plot']: leg.AddEntry(hist[s], sample[s]['label'].replace("m_{#Chi}=1 GeV", ""), "fl")
        
    ### data/bkg ratio and systematics
    err = hist['BkgSum'].Clone("BkgErr;")
    err.SetTitle("")
    err.GetYaxis().SetTitle("data / bkg")
    if fit == "prefit":
        err.SetFillColor(ROOT.kOrange-2)
    elif fit == "postfit":
        err.SetFillColor(ROOT.kBlue)
    err.SetFillStyle(3001)
    for i in range(1, err.GetNbinsX()+1):
        err.SetBinContent(i, 1)
        if hist['BkgSum'].GetBinContent(i) > 0:
            err.SetBinError(i, hist['BkgSum'].GetBinError(i)/hist['BkgSum'].GetBinContent(i))
    errLine = err.Clone("errLine")
    errLine.SetLineWidth(1)
    errLine.SetFillStyle(0)
    res = hist['data_obs'].Clone("residues")
    for i in range(0, res.GetNbinsX()+1):
        if hist['BkgSum'].GetBinContent(i) > 0: 
            res.SetBinContent(i, res.GetBinContent(i)/hist['BkgSum'].GetBinContent(i))
            res.SetBinError(i, res.GetBinError(i)/hist['BkgSum'].GetBinContent(i))

    # Legend
    leg1 = TLegend(0.12, 0.45, 0.25, 0.5)
    leg1.SetBorderSize(0)
    leg1.SetFillStyle(0)
    leg1.SetFillColor(0)
    leg1.SetTextSize(0.05)
    leg1.AddEntry(err, "systematic uncertainty ("+fit+")", "f")
    
    # --- Display ---
    c1 = TCanvas("c1", hist.values()[0].GetXaxis().GetTitle(), 800, 800 if ratio else 600)
    gStyle.SetOptStat(0)
    if ratio:
        c1.Divide(1, 2)
        setTopPad(c1.GetPad(1), ratio)
        setBotPad(c1.GetPad(2), ratio)
    c1.cd(1)
    c1.GetPad(bool(ratio)).SetTopMargin(0.08)
    c1.GetPad(bool(ratio)).SetRightMargin(0.05)
    c1.GetPad(bool(ratio)).SetTicks(1, 1)
    if log:
        c1.GetPad(bool(ratio)).SetLogy()
        
    # Draw
    bkg.Draw("HIST") # stack
#    hist['BkgSum'].Draw("SAME, E2") # sum of bkg
    if len(data) > 0: hist['data_obs'].Draw("SAME, PE") # data
    for i, s in enumerate(sign):
        if sample[s]['plot']:
            hist[s].DrawNormalized("SAME, HIST", hist[s].Integral()*snorm) # signals
        pass
  
    bkg.SetMaximum((2. if log else 1.1)*max(bkg.GetMaximum(), hist['data_obs'].GetMaximum()+hist['data_obs'].GetBinError(hist['data_obs'].GetMaximumBin())))
    bkg.SetMinimum(max(min(hist['BkgSum'].GetBinContent(hist['BkgSum'].GetMinimumBin()), hist['data_obs'].GetMinimum()), 1.e-1)  if log else 0.)
    if log:
        bkg.GetYaxis().SetNoExponent(bkg.GetMaximum() < 1.e4)
        bkg.GetYaxis().SetMoreLogLabels(True)
    
    setHistStyle(bkg, 1.1 if ratio else 1.0)
       
    leg.Draw()
#    drawCMS(lumi, "Preliminary",True)
    drawCMS(lumi, "Private", True)
    drawRegion(channel)    
    drawAnalysis(channel)
    drawFit(fit)
   
    if ratio:
        c1.cd(2)
        setBotStyle(err,3,True)
        setBotStyle(res)
        err.Draw("E2")
        errLine.Draw("SAME, HIST")
        if len(data) > 0:
            res.Draw("SAME, PE0")
            drawRatio(hist['data_obs'], hist['BkgSum'])
            drawKolmogorov(hist['data_obs'], hist['BkgSum'])

        leg1.Draw()


    c1.Update()
    
    # return list of objects created by the draw() function
    return [c1, bkg, leg, leg1, err, errLine, res]




##################################################

def applyStyle(hist):
    for s, h in hist.iteritems():
        if s=="Data" or s=="data_obs":
            h.SetMarkerStyle(20)
            h.SetMarkerSize(1.25)
        elif s=="BkgSum":
            h.SetFillStyle(3003)
            h.SetFillColor(1)
        else:
            h.SetFillColor(sample[s]['fillcolor'])
            h.SetFillStyle(sample[s]['fillstyle'])
            h.SetLineColor(sample[s]['linecolor'])
            h.SetLineWidth(sample[s]['linewidth'])
            h.SetLineStyle(sample[s]['linestyle'])


def printTable(hist, sign=[]):
    samples = [x for x in hist.keys() if not 'data' in x and not 'BkgSum' in x and not x in sign]
    print "Sample                  Events          Entries         %"
    print "-"*80
    for i, s in enumerate(['data_obs']+samples+['BkgSum']):
        if i==1 or i==len(samples)+1: print "-"*80
        print "%-20s" % s, "\t%-10.2f" % hist[s].Integral(), "\t%-10.0f" % (hist[s].GetEntries()-2), "\t%-10.2f" % (100.*hist[s].Integral()/hist['BkgSum'].Integral()) if hist['BkgSum'].Integral() > 0 else 0, "%"
    print "-"*80
    for i, s in enumerate(sign):
        if not sample[s]['plot']: continue
        print "%-20s" % s, "\t%-10.2f" % hist[s].Integral(), "\t%-10.0f" % (hist[s].GetEntries()-2), "\t%-10.2f" % (100.*hist[s].GetEntries()/float(hist[s].GetOption())) if float(hist[s].GetOption()) > 0 else 0, "%"    
    print "-"*80




##################
### DRAW UTILS ###
##################

def drawCMS(lumi, text, onTop=False):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.05)
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(33)
    if float(lumi) > 0: latex.DrawLatex(0.95, 0.98, "%.2f fb^{-1}  (13 TeV)" % (float(lumi)/1000.))
    if not onTop: latex.SetTextAlign(11)
    latex.SetTextFont(62)
    latex.SetTextSize(0.055)
    if not onTop: latex.DrawLatex(0.2, 0.83, "CMS")
    else: latex.DrawLatex(0.2, 0.98, "CMS")
    latex.SetTextSize(0.055)
    latex.SetTextFont(52)
    if not onTop: latex.DrawLatex(0.33, 0.83, text)
    else: latex.DrawLatex(0.33, 0.98, text) # 0.40 preliminary; 0.33 private
    latex.SetTextSize(0.04)
    latex.SetTextFont(62)
    latex.SetTextAlign(13)

def drawAnalysis(s):
    text = ""
    if "DM" in s: text = "DM + t#bar{t}"
    else: return True
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)
    latex.SetTextFont(42)
    #latex.SetTextAlign(33)
    latex.DrawLatex(0.15, 0.95, text)

def drawFit(s):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextFont(72)
    latex.DrawLatex(0.46, 0.8, s)

def drawRegion(region, left=False):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextFont(72) #52
    latex.SetTextSize(0.035)
    if left: latex.DrawLatex(0.15, 0.75, region)
    else:
        latex.SetTextAlign(22)
        latex.DrawLatex(0.5, 0.85, region)

def drawMass(m):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextAlign(22)
    latex.SetTextSize(0.04)
    latex.DrawLatex(0.75, 0.85, "m_{X} = %.0f GeV" % m)

def drawChi2(f):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextColor(1)
    latex.SetTextFont(62)
    latex.SetTextSize(0.08)
    latex.DrawLatex(0.75, 0.85, "#chi^{2}/ndf = %.3f" % f.chiSquare())

def drawMediator(med):
    text = "Pseudoscalar" if 'ps' in med else "Scalar"
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.035)
    latex.SetTextFont(42)
    latex.DrawLatex(0.15, 0.70, text+" mediator")
    latex.DrawLatex(0.15, 0.65, "m_{#chi}=1 GeV")


def drawNorm(y, text, secondline=""):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextColor(1)
    latex.SetTextFont(42)
    latex.SetTextSize(0.025)
    latex.DrawLatex(0.72, y, text)
    if len(secondline) > 0: latex.DrawLatex(0.72, y-0.035, secondline)
    
def drawRatio(data, bkg,x=0.12,y=0.85,color=ROOT.kBlack):
    errData = array('d', [1.0])
    errBkg = array('d', [1.0])
    intData = data.IntegralAndError(1, data.GetNbinsX(), errData)
    intBkg = bkg.IntegralAndError(1, bkg.GetNbinsX(), errBkg)
    ratio = intData / intBkg if intBkg!=0 else 0.
    error = math.hypot(errData[0]*ratio/intData,  errBkg[0]*ratio/intBkg) if intData>0 and intBkg>0 else 0
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextColor(color)
    latex.SetTextFont(62)
    latex.SetTextSize(0.08)
    latex.DrawLatex(x, y, "data/bkg = %.2f #pm %.2f" % (ratio, error))
    print "  ratio:\t%.2f +- %.2f" % (ratio, error)
    #return [ratio, error]

def drawKolmogorov(data, bkg,x=0.73,y=0.85,color=ROOT.kBlack):
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextColor(color)
    latex.SetTextFont(62)
    latex.SetTextSize(0.08)
    latex.DrawLatex(x, y, "#chi^{2}/ndf = %.2f,   KS = %.2f" % (data.Chi2Test(bkg, "CHI2/NDF"), data.KolmogorovTest(bkg)))
    print "  CHI2/NDF:\t%.2f" % (data.Chi2Test(bkg, "CHI2/NDF"))
    print "  KS:\t\t%.2f" % (data.KolmogorovTest(bkg))

def drawCut(cut, ymin, ymax):
    line = TLine()
    line.SetLineWidth(2)
    line.SetLineStyle(7)
    line.SetLineColor(1)
    line.PaintLineNDC(cut, ymin, cut, ymax)

def setHistStyle(hist, r=1.1):
    hist.GetXaxis().SetTitleSize(hist.GetXaxis().GetTitleSize()*r)
    hist.GetYaxis().SetTitleSize(hist.GetYaxis().GetTitleSize()*r)
    hist.GetXaxis().SetLabelOffset(hist.GetXaxis().GetLabelOffset()*r*r*r*r*2)
    hist.GetXaxis().SetTitleOffset(hist.GetXaxis().GetTitleOffset()*r*r)
    hist.GetYaxis().SetTitleOffset(hist.GetYaxis().GetTitleOffset()*r)
    if hist.GetXaxis().GetTitle().find("MET") != -1: # and not hist.GetXaxis().IsVariableBinSize()
        div = (hist.GetXaxis().GetXmax() - hist.GetXaxis().GetXmin()) / hist.GetXaxis().GetNbins()
        hist.GetYaxis().SetTitle("events / %.2f GeV" % div)

def addOverflow(hist, addUnder=True):
    n = hist.GetNbinsX()
    hist.SetBinContent(n, hist.GetBinContent(n) + hist.GetBinContent(n+1))
    hist.SetBinError(n, math.sqrt( hist.GetBinError(n)**2 + hist.GetBinError(n+1)**2 ) )
    hist.SetBinContent(n+1, 0.)
    hist.SetBinError(n+1, 0.)
    if addUnder:
        hist.SetBinContent(1, hist.GetBinContent(0) + hist.GetBinContent(1))
        hist.SetBinError(1, math.sqrt( hist.GetBinError(0)**2 + hist.GetBinError(1)**2 ) )
        hist.SetBinContent(0, 0.)
        hist.SetBinError(0, 0.)

def setTopPad(TopPad, r=4):
    print TopPad
    TopPad.SetPad("TopPad", "", 0., 1./r, 1.0, 1.0, 0, -1, 0)
    TopPad.SetTopMargin(0.24/r)
    TopPad.SetBottomMargin(0.04/r)
    TopPad.SetRightMargin(0.05)
    TopPad.SetTicks(1, 1)

def setBotPad(BotPad, r=4):
    BotPad.SetPad("BotPad", "", 0., 0., 1.0, 1./r, 0, -1, 0)
    BotPad.SetTopMargin(r/100.)
    BotPad.SetBottomMargin(r/10.)
    BotPad.SetRightMargin(0.05)
    BotPad.SetTicks(1, 1)

def setBotStyle(h, r=4, fixRange=True):
    h.GetXaxis().SetLabelSize(h.GetXaxis().GetLabelSize()*(r-1));
    h.GetXaxis().SetLabelOffset(h.GetXaxis().GetLabelOffset()*(r-1));
    h.GetXaxis().SetTitleSize(h.GetXaxis().GetTitleSize()*(r-1));
    h.GetYaxis().SetLabelSize(h.GetYaxis().GetLabelSize()*(r-1));
    h.GetYaxis().SetNdivisions(505);
    h.GetYaxis().SetTitleSize(h.GetYaxis().GetTitleSize()*(r-1));
    h.GetYaxis().SetTitleOffset(h.GetYaxis().GetTitleOffset()/(r+1));
    if fixRange:
        h.GetYaxis().SetRangeUser(0., 2.)
        for i in range(1, h.GetNbinsX()+1):
            if h.GetBinContent(i)<1.e-6:
                h.SetBinContent(i, -1.e-6)
    
def drawCut(hist):
    #drawCut(80, 140, 0., hist['BkgSum'].GetMaximum())
    line1 = TLine(80, 0, 80, hist['BkgSum'].GetMaximum())
    line1.SetLineWidth(2)
    line1.SetLineStyle(7)
    line1.SetLineColor(1)
    line1.Draw()
    
    line2 = TLine(140, 0, 140, hist['BkgSum'].GetMaximum())
    line2.SetLineWidth(2)
    line2.SetLineStyle(7)
    line2.SetLineColor(1)
    line2.Draw()
    
    line1 = TLine(0.841, 0, 0.841, 15)
    line1.SetLineWidth(2)
    line1.SetLineStyle(7)
    line1.SetLineColor(1)
    line1.Draw()
    
    line1 = TLine(100, 0, 100, hist['BkgSum'].GetMaximum())
    line1.SetLineWidth(2)
    line1.SetLineStyle(7)
    line1.SetLineColor(1)
    line1.Draw()


