### Decorrelate Q2 syst among the processes

import sys
import os, commands
import shutil

from samples.toPlot import samples

import optparse


usage = 'usage: %prog -d histosDir'
parser = optparse.OptionParser(usage)

parser.add_option('-d', '--dir', dest='path', type='string', default = './histos/', help='Histos folder')
(opt, args) = parser.parse_args()


channels = ["fullhadronic", "semileptonic"]
syst = "q2"
vs = ["Up", "Down"]

newsyst = ["q2TT", "q2SingleTop", "q2VV", "q2QCD", "q2DY", "q2WJets", "q2ZToNuNu", "q2otherBkg", "q2DMtt", "q2DMj", "q2ZJets"]


def createSyst(versus, verbose = True):

    for s in samples.itervalues():
        for ch in channels:
            src = s.label + "_" + ch + ".root"
            dst =  [s.label + "_" + ch + "_" + nSyst + versus + ".root" for nSyst in newsyst]

            if(verbose):
                for d in dst:
                    print "Copying " + opt.path + src + " to "+ opt.path + d
    
            [shutil.copyfile(opt.path + src, opt.path + d) for d in dst]
    
            srcSys = s.label + "_" + ch + "_" + syst + versus + ".root"  
            dstSys = s.label + "_" + ch + "_" + syst + s.label + versus + ".root" 
            if s.label.startswith("DM"):
                dstSys = s.label + "_" + ch + "_" + syst + "DMtt"+ versus +".root"
            if (s.label.startswith("ZToNuNu") or s.label.startswith("DY")):
                
                dstSys = s.label + "_" + ch + "_" + syst + "ZJets"+ versus +".root"

    
            if verbose: print "Copying " + opt.path + srcSys + " to "+ opt.path + dstSys
            shutil.copyfile(opt.path + srcSys, opt.path + dstSys)


createSyst("Up")
createSyst("Down")
