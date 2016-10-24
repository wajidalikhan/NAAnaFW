#!/bin/env python

import optparse

# Partse command line before ANY action
usage = 'usage: %prog -l lumi'
parser = optparse.OptionParser(usage)

#parser.add_option('-l', '--lumi', dest='lumi', type='float', default = '0.20923', help='Luminosity')
parser.add_option('-c', '--channel', dest='channel', type='string', default = 'semileptonic', help='Channel to analyze: semileptonic or fullhadronic')
parser.add_option('-s', '--sys', dest='sys', type='string', default = 'noSys', help='Systematics: noSys, jesUp, jesDown')
#parser.add_option('-n', '--normData', dest='normData', type='int', default = '0', help='Normalise to data?')
#parser.add_option('-r', '--resdir', dest='resdir', type='string', default = './', help='res directory')
#parser.add_option('-d', '--no-data', dest='data', action='store_false', default = True, help='Hide data')
#parser.add_option('--focusOn', dest='focus', default = None, help='Focus on a single plot')

(opt, args) = parser.parse_args()
 
import time
import sys
import ROOT
import copy
import commands, os.path
import numpy as n

# from plots.services import Histo, Stack, Legend, deltaPhi, Histo1
#from plots.services import deltaPhi
from samples.toPlot import samples
from plots.services import Histo
#import plots.common, plots.fullhadronic, plots.semileptonic


import tdrstyle, CMS_lumi

def sumerrors(h):
    return sum([h.GetBinError(ib) for ib in xrange(0,h.GetNbinsX()+2)])


ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)

tdrstyle.setTDRStyle();

#settings = {}
#store = []

#settings.update(plots.common.settings)
#store += plots.common.store

if opt.channel == 'semileptonic':
    # Add semileptonic specific settings and plots to store
#    settings.update(plots.semileptonic.settings)
#    store += plots.semileptonic.store

    outhistos = 'output/sl/histos_lin'
    outpdfs = 'output/sl/pdfs_lin'
    outtxt = 'output/sl/txt_lin'

elif opt.channel == 'fullhadronic':
    # Add fullhadronic specific settings and plots to store
#    settings.update(plots.fullhadronic.settings)
#    store += plots.fullhadronic.store

    outhistos = 'output/fh/histos_lin'
    outpdfs = 'output/fh/pdfs_lin'
    outtxt = 'output/fh/txt_lin'

else:
    print 'What?!?!'
    import sys
    sys.exit(0)


vars = {"metFinal":"", "metFinal_2lep":"CR_TT", "metFinal_met_0btag":"CR_WJets"}
if(opt.channel == "fullhadronic"):
    vars = {"metFinal":"", "metFinal_SR_1lep":"CR_TT", "metFinal_CR5":"CR_VJets", "metFinal_CR6nw":"CR_WJets", "metFinal_CR7nw":"CR_ZJets"}

    
def shift(h, ib, shift, slabel, channel, var):
    corr = 1
    if(shift == "Down"): corr = -1
    oh = 0.
    eoh = 0.
   # nhName = ""
    oh = h.GetBinContent(ib)
    eoh = h.GetBinError(ib)
    oh += eoh * corr
    #print "++ Bin 1: ", h.GetAt(ib)
    #print "++ W Bin 1: ", h.GetBinError(ib)
    
    nh = copy.deepcopy(h)
    nh.SetAt(oh, ib)
    #print "++ nh Bin 1: ", nh.GetBinContent(ib)
    #print "---- Old name: ",h.GetName()
    if(vars[var]!=""): nhName = h.GetName() + "_mcstat_"+ channel + "_" + vars[var] + "_"+ slabel+"_bin"+str(ib)+"_"+shift
    else: nhName = h.GetName() + "_mcstat_"+ channel + "_" + slabel+"_bin"+str(ib)+"_"+shift
    nh.SetName(nhName)
    print "---- New name: ",nh.GetName()
    #nh.GetSumw2().SetAt(X, ib)
    return nh
    

def makeBinHistos(var):
    for s in samples.itervalues():
        hfilename = "%s/%s_%s.root" % (outhistos,s.label,opt.channel)
    
        hfile = ROOT.TFile(hfilename, "UPDATE")

        print "file ", hfile

        hin = hfile.Get(var)
        if not isinstance(hin, ROOT.TH1):
            raise RuntimeError('Failed to load histogram %s from %s' % (hin, var))
        else:
            ibin = 7
            #print "Sample: ",s.label
            #print "Variable: ",hin.GetName()
            #print "Integral: ",hin.Integral()
            #print "Number of bins: ",hin.GetXaxis().GetNbins()
            #        nh = shift(hin, ib, "up")
            #nh = shift(hin, ibin, "down")
            mcStat_down = [ shift(hin, ib, "Down",s.label, opt.channel, var) for ib in xrange(1, hin.GetXaxis().GetNbins()+1) ]
            mcStat_up = [ shift(hin, ib, "Up", s.label, opt.channel, var) for ib in xrange(1, hin.GetXaxis().GetNbins()+1) ]
                    
            hfile.cd()
            [nh.Write() for nh in mcStat_down]
            [nh.Write() for nh in mcStat_up]
            hfile.Close()

        

            
[makeBinHistos(v) for v in vars.keys()]
