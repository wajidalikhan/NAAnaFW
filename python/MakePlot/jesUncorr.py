### Decorrelate Q2 syst among the processes

import sys
import os, commands
import shutil
import ROOT
import copy
from samples.toPlot import samples
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)
import optparse


usage = 'usage: %prog -d histosDir'
parser = optparse.OptionParser(usage)

parser.add_option('-d', '--dir', dest='path', type='string', default = './histos/', help='Histos folder')
parser.add_option('', '--debug', dest='debug', type='int', default = '1', help='Increment level of verbosity for debugging purposes')
parser.add_option('-s', '--sys', dest='sys', type='string', default = 'jes', help='Systematics to de-correlate')

(opt, args) = parser.parse_args()


channels = ["fullhadronic", "semileptonic"]

varsSl = {"metFinal":"SR_SL", "metFinal_2lep":"CR_TT_SL", "metFinal_met_0btag":"CR_WJets_SL"}

varsFh = {"metFinal":"SR_FH", "metFinal_SR_1lep":"CR_TT_FH", "metFinal_CR5":"CR_VJets_FH", "metFinal_CR6nw":"CR_WJets_FH"}#,  "metFinal_outtop":"CR_QCD","metFinal_CR7nw":"CR_ZJets"}

regsFh_SlUnc = ["SR_SL", "CR_TT_SL", "CR_WJets_SL"]

regsSl_FhUnc = ["SR_FH", "CR_TT_FH", "CR_WJets_FH","CR_VJets_FH"]


 

def raiseError(histo, var, filename): 
    if not isinstance(histo, ROOT.TH1):
        raise RuntimeError('Failed to load histogram  %s from file  %s' % ( var, filename))

                
def symmetrize(syst):
    for s in samples.itervalues():
        for ch in channels:

            nom = opt.path + s.label + "_" + ch + ".root"
            up = opt.path + s.label + "_" + ch + "_"+syst+"Up.root"
            down = opt.path + s.label + "_" + ch + "_"+syst+"Down.root"
           
            f_nom = ROOT.TFile(nom)
            f_up = ROOT.TFile(up)
            f_down = ROOT.TFile(down)

            
            vars = varsFh
            regsOff=regsFh_SlUnc
            if (ch=="semileptonic"): 
                vars = varsSl
                regsOff=regsSl_FhUnc

            
            for  v, region in vars.iteritems():
                sysUp = opt.path + s.label + "_" + ch +"_"+ syst  + region + "Up.root"
                sysDown =opt.path + s.label + "_" + ch +"_"+ syst  + region + "Down.root"
                #print "++ New File Up: ", sysUp
                #print "++ New File Down: ", sysDown

                print "+ Copying " + nom + " to "+ sysUp
                shutil.copyfile(nom, sysUp)
                print "+ Copying " + nom + " to "+ sysDown
                shutil.copyfile(nom, sysDown)


                h_up = f_up.Get(v)
                raiseError(h_up, v, up)
                h_down = f_down.Get(v)
                raiseError(h_down, v, down)
                
                nf_up = ROOT.TFile(sysUp, "UPDATE")
                nf_up.cd()
                nh_up = copy.deepcopy(h_up)
                nh_up.Write()
                nf_up.Close()
                
                nf_down = ROOT.TFile(sysDown, "UPDATE")
                nf_down.cd()
                nh_down = copy.deepcopy(h_down)
                nh_down.Write()
                nf_down.Close()
            
            for region in regsOff:
                sysUp = opt.path + s.label + "_" + ch +"_"+ syst  + region + "Up.root"
                sysDown =opt.path + s.label + "_" + ch +"_"+ syst  + region + "Down.root"
                print "+ Copying " + nom + " to "+ sysUp
                shutil.copyfile(nom, sysUp)
                print "+ Copying " + nom + " to "+ sysDown
                shutil.copyfile(nom, sysDown)
                
                

symmetrize(opt.sys)

