######################################
#
# Annapaola de Cosa, January 2015
#
######################################

import time

import sys
#sys.argv.append('-b')
import ROOT
import commands, os
import numpy

ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)
#ROOT.gSystem.Load("libFWCoreFWLite.so");
#ROOT.AutoLibraryLoader.enable();
#ROOT.gSystem.Load("libDataFormatsFWLite.so");

c = ROOT.TCanvas()

debug = False

# ====================================================
# Class Histo: to manage histograms and setting style
# ====================================================

class Histo(object):

    def __init__(self, name="", title="", nBins=100, xmin=0, xmax=100):
        self._name = name
        self._title = title
        self._nBins = nBins
        self._xmin = xmin
        self._xmax = xmax
        self._h = ROOT.TH1F( self._name, "",  self._nBins, self._xmin, self._xmax)

    def SetStyle(self, color, style = 0, fill = 0):
        self._h.GetXaxis().SetLabelFont(42);
        self._h.GetYaxis().SetLabelFont(42);
        self._h.GetXaxis().SetTitleFont(42);
        self._h.GetYaxis().SetTitleFont(42);
        self._h.GetXaxis().SetTitleOffset(0.9);
        self._h.GetYaxis().SetTitleOffset(1.2);
        self._h.SetTitleFont(42);
        self._h.SetTitle("");
        if(color != ROOT.kRed): self._h.SetLineColor(1)
        else: self._h.SetLineColor(color)
        self._h.SetLineWidth(1)
        self._h.SetLineStyle(style)
        self._h.SetFillColor(color)
        self._h.SetFillStyle(fill)
        self._h.GetXaxis().SetTitle(self._title)
        nEvts = (self._h.GetXaxis().GetXmax() - self._h.GetXaxis().GetXmin()) / self._h.GetNbinsX()
        self._h.GetYaxis().SetTitle("Events/"+ str.format("{0:.2f}", nEvts) );
    
    def Fill(self, value, weight = "1"):
        self._h.Fill(value, weight)

    def Rebin(self, value):
        self._h.Rebin(value)

    def Scale(self, value):
        self._h.Scale(value)

    def Integral(self):
        return self._h.Integral()

    def Draw(self, options=""):
        self._h.Draw(options)

    def Normalize(self):
        if(self._h.Integral()>0):self._h.Scale(1./self._h.Integral())
        color = self._h.GetFillColor()
        self._h.SetLineColor(color)        
        self._h.SetLineWidth(2)        
        self._h.SetFillStyle(0)

    def SetBinContent(self, bin, content):
        self._h.SetBinContent(bin, content)

    def GetBinContent(self, bin):
        return self._h.GetBinContent(bin)

    def GetHisto(self):
        return self._h



class Histo1(Histo):

    def __init__(self, histo):
        self._h = histo.Clone()
        self._title = histo.GetTitle()
        self._name = histo.GetName()
    
# ====================================================
# Class Stack: to manage THStacks
# ====================================================

class Stack(object):

    def __init__(self, name, title):
        self._name = name
        self._title = title
        self._hs = ROOT.THStack(self._name,"")
        self._latex = ROOT.TLatex()
        self._latex.SetNDC()
        self._latex.SetTextSize(0.04)
        self._latex.SetTextFont(42)
        self._latex.SetTextAlign(11)


    def SetStyle(self, options = ""):
        self._hs.Draw(options)
        self._hs.GetHistogram().GetXaxis().SetLabelFont(42);
        self._hs.GetHistogram().GetYaxis().SetLabelFont(42);
        self._hs.GetHistogram().GetXaxis().SetTitleFont(42);
        self._hs.GetHistogram().GetYaxis().SetTitleFont(42);
        self._hs.GetHistogram().GetXaxis().SetTitleOffset(0.9);
        self._hs.GetHistogram().GetYaxis().SetTitleOffset(1.);
        self._hs.GetHistogram().SetTitleFont(42);
        self._hs.GetHistogram().SetTitle("");

        self._hs.GetHistogram().GetXaxis().SetLabelSize(0.05);
        self._hs.GetHistogram().GetYaxis().SetLabelSize(0.05);
        self._hs.GetHistogram().GetXaxis().SetTitleSize(0.06);
        self._hs.GetHistogram().GetYaxis().SetTitleSize(0.06);

        self._hs.GetHistogram().GetXaxis().SetTitle(self._title)
        nEvts = (self._hs.GetHistogram().GetXaxis().GetXmax() - self._hs.GetHistogram().GetXaxis().GetXmin()) / self._hs.GetHistogram().GetNbinsX()
        self._hs.GetHistogram().GetYaxis().SetTitle("Number of events / "+ str.format("{0:.2f}", nEvts) );
 #       self._hs.Draw(options)


    def Add(self, h, options="hist"):
        self._hs.Add(h.GetHisto(), options)

    def Draw(self, options = ""):
        self._hs.Draw(options)

    def GetHistogram(self):
        return self._hs.GetHistogram()

    def GetMaximum(self):
        return self._hs.GetMaximum()

    def GetMinimum(self):
        return self._hs.GetMinimum()

    def SetMaximum(self, maximum):
        self._hs.SetMaximum(maximum)

    def SetMinimum(self, minimum):
        self._hs.SetMinimum(minimum)

    def MakeStack(self, h, options = "hist"):
        self._hs.Add(h.GetHisto(), options)
        return self._hs


    def DrawStack(self, lumi = 1, opt = ''):
        self.SetStyle(opt)
#        self.Draw(opt)
#        self._latex.DrawLatex(0.1, 0.92, "CMS preliminary");
#        self._latex.DrawLatex(0.6,0.92, " 40.1 pb^{-1} at #sqrt{s} = 13 TeV");



    def SetRangeUser(self, xmin, xmax):
        self._hs.Draw()
        self._hs.GetHistogram().GetXaxis().SetRangeUser(xmin, xmax)

# ====================================================
# Class Legend: to draw a legend
# ====================================================

class Legend(object):

    def __init__(self):
        #self._leg = ROOT.TLegend(.61, .58, .94, .89)
        self._leg = ROOT.TLegend(0.47, 0.65, 0.94, 0.89)
        self._leg.SetNColumns(3)
        self._leg.SetFillColor(0)
        self._leg.SetTextSize(0.045)
        self._leg.SetTextFont(42)
#        self._leg.SetBorderSize(0)
        
    def AddEntry(self, h, label, option = "f"):
        self._leg.AddEntry(h.GetHisto(), label, option)

    def Draw(self, options = "SAME"):
        self._leg.Draw(options)


# ====================================================
# Auxiliary functions
# ====================================================


def deltaPhi(phi1, phi2):
    delta_phi = numpy.abs(phi1 - phi2)
    if delta_phi > numpy.pi:
        delta_phi = 2 * numpy.pi - delta_phi
    return delta_phi



# ====================================================
# Class Tree: to manage trees
# ====================================================


class Tree(object):

    def __init__(self, filename, treename):
        self._filename = filename
        self._treename = treename
        self._tree = None
        self._file = None
        self._outfile = None
        self._newtree = None
        self._cutLabel = ''
        
    def _openFile(self, filename = '', option = ''):
        if filename == '': filename = self._filename
        file = ROOT.TFile.Open(filename, option)
        if not file.__nonzero__() or not file.IsOpen():
            raise NameError('File '+filename+' not open')
        return file
    
    def initTree(self):
        self._file = self._openFile()
        self._tree = self._file.Get(self._treename)

    def GetTree(self):
        return self._tree

    def SetLabel(self, label):
        self.label = label
        
    def GetEntries(self, cut = ''):
        return self._tree.GetEntries(cut)

    def SetCutLabel(self, label):
        self._cutLabel = label

    def GetCutLabel():
        return self._cutLabel
        
    def disconnect(self):
        self._newtree.Write()
        self._outfile.Close()
        self._file.Close()
        
        self._file = None
        self._tree = None
        self._outfile = None
        self._newtree = None
        
    def cloneTree(self, outfilename, branches = []):

        #self._tree.SetBranchStatus('*'  ,0)        
        #for b in branches:
         #   print"Branch Alias:", b
          #  bName = self._tree.GetAlias(b)
           # print "Branch Name", bName
            #self._tree.SetBranchStatus(bName,1)
            #print "brach turned on"

        self._outfile = self._openFile(outfilename, 'recreate')            
        self._newtree = self._tree.CloneTree()
        print "tree cloned"
            ## BUT keep all branches "active" in the old tree
        self._tree.SetBranchStatus('*'  ,1)

            
    def filterTree(self, cut):
        self._tree.Draw(">>elist", cut, "entrylist")
        self._list = ROOT.TEntryList(ROOT.gDirectory.Get("elist") )
        self._tree.SetEntryList(self._list)
        return self._list

    def saveEvtList(self):
        tp = ROOT.TTreePlayer()
        tp = self._tree.GetPlayer()
        tp.SetScanRedirect(ROOT.kTRUE) # You may have to cast t.GetPlayer() to a TTreePlayer*
        tp.SetScanFileName("output_"+self._cutLabel+".txt")
        self._tree.Scan("LumiBlock:RunNum:EvtNum")



    def pruneClonedTree(self):

        for i in xrange(self._list.GetN()):
            
            self._tree.GetEntry(self._list.GetEntry(i))
            
            self._newtree.Fill()

        self.disconnect()
            
                                                        
# ====================================================
# Code to create a cutFlow
# ====================================================

    
cutNames = ["Trigger", "OneTightLepton", "LooseMuonVeto", "LooseElectronVeto", "nJets#ge 1", "nJets#ge 2", "nJets#ge 3", "nJets#ge 4"]
elcuts = ["(nTightEl == 1 && tightElectronsMvaTrigV0 > 0.5 && tightElectronsTrigPresel)",
        "nVetoMu == 0 ",
        "nVetoEl == 1 ",
        "(nJets >= 1 && topJetsPFPtNoJER > 30.)",
        "(nJets >= 2 && topJetsPFPtNoJER > 30.)",
        "(nJets >= 3 && topJetsPFPtNoJER > 30.)",
        "(nJets >= 4 && topJetsPFPtNoJER > 30.)"
        ]


mucuts = ["(nTightMu == 1 )",
        "nVetoMu == 1 ",
        "nVetoEl == 0 ",
        "(nJets >= 1 && topJetsPFPtNoJER > 30.)",
        "(nJets >= 2 && topJetsPFPtNoJER > 30.)",
        "(nJets >= 3 && topJetsPFPtNoJER > 30.)",
        "(nJets >= 4 && topJetsPFPtNoJER > 30.)"
        ]


def _cutFlow(tree,cuts, ch):

    histo = Histo("CutFlow", "CutFlow", len(cutNames), 0., len(cutNames))
    histo.SetStyle(ROOT.kBlue, "")
    h = histo.GetHisto()
    h.SetBinContent(1, tree.GetEntries())
    h.GetXaxis().SetBinLabel(1, cutNames[0])
    print 'Step 0: ', h.GetBinContent(1)
    tree.SetCutLabel(ch +"_Step0")
    tree.saveEvtList()
    
    list = tree.filterTree("topJetsPFPtNoJER > 30.")
    print 'cut on jet pt' , tree.GetEntries("topJetsPFPtNoJER > 30.")
    
    cut = ''
    strcut = ''
    for i in range(len(cuts)):
        cut = cuts[i]
        strcut += cuts[i]
        strcut += " && "
        if(debug): print "CUT ----- ", strcut
        print "Step " + str(i+1) + ": " + str(tree.GetEntries(cut) )
        list = tree.filterTree(cut)
        h.SetBinContent(i+2, list.GetN())
        h.GetXaxis().SetBinLabel(i+2, cutNames[i+1])
        tree.SetCutLabel(ch +"_Step"+str(i+1))
        tree.saveEvtList()

        
    h.Draw()
    c.Print("cutFlow.pdf")
    return cut



    
def setAliases(t):
    t.SetAlias("nTightEl","floats_nTupleElectrons_tightElectronsE_SingleTop.@obj.size()")
    t.SetAlias("nTightMu","floats_nTupleMuons_tightMuonsE_SingleTop.@obj.size()")
    t.SetAlias("nVetoEl","floats_nTupleVetoElectrons_vetoElectronsE_SingleTop.@obj.size()")
    t.SetAlias("nVetoMu","floats_nTupleVetoMuons_vetoMuonsE_SingleTop.@obj.size()")
    t.SetAlias("nJets","floats_nTupleTopJetsPF_topJetsPFE_SingleTop.@obj.size()")
    t.SetAlias("nJets","floats_nTupleTopJetsPF_topJetsPFENoJER_SingleTop.@obj.size()")
    t.SetAlias("LumiBlock", "topJetsPFLumiblock")
    t.SetAlias("RunNum", "topJetsPFRunNumber")
    t.SetAlias("EvtNum", "topJetsPFEventNumber")

    
def doCutFlow(filename, cuts, ch):

    tree = Tree(filename, "Events")
    tree.initTree()
    t = tree.GetTree()

    ### CLONE DOES NOT WORK
    #tree.cloneTree("newtree.root", ["vetoMuonsE"])    
    #    file = ROOT.TFile(filename)
    #   t = file.Get("Events")
    setAliases(t)
    t.Draw("nJets")
    t.ls()
    c.Print("nJets.pdf")
    cut = _cutFlow(tree, cuts, ch)
    #t.Draw("nJets", "(nTightMu == 1 ) && nVetoMu == 1  && nVetoEl == 0  ")
    #c.Print("nJets_mu.pdf")
    #tree.pruneClonedTree()
    #    t.Draw("nJets")
    #t.ls()

    #c.Print("nJetsEntryList.pdf")

    #list.Print("all")


def plot(var, t1, t2, nBins = 100, xmin = 0,xmax = 100, cut = ''):

    setAliases(t1)
    setAliases(t2)
    def setH(sample, t, color):
        hist = Histo( sample, sample, nBins, xmin, xmax)
        hist.SetStyle(color, var)
        h = hist.GetHisto()
        h.SetDirectory(ROOT.gDirectory.func())
        #print type(h)
        #print "check1 ", h.GetHisto().ls()
        # print ROOT.gDirectory.ls()
        t.Draw(var+">>"+sample, cut)
        
        if(h.Integral()!=0):h.Scale(1./h.Integral())
#        h.SetLineColor(color)
        #print ROOT.gDirectory.ls()
        
        #h = t1.GetHistogram()
        #print ROOT.gDirectory.ls()
        #print h.Integral()

        #h.SetDirectory(0)
        return h
    
    
    b = setH("background", t1, ROOT.kBlue +2)
    #print 'b integral ', b.Integral()
    #print ROOT.gDirectory.ls()
    s = setH("signal", t2, ROOT.kRed +2)
    #assert(False)
    cplot = ROOT.TCanvas()
    s.SetMaximum(max([s.GetMaximum(), b.GetMaximum()]))
    #h.SetMinimum(minim)
    s.Draw("HIST")
    b.Draw("HISTSAME")
    cplot.Print(var+".pdf")
    cplot.Print(var+".png")
    #cplot.Print(var+".root")
    del s, b
