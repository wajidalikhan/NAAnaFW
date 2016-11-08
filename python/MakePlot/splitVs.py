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



def createSyst(syst, newsyst, versus, verbose = False):

    verbose = True 
    
    for s in samples.itervalues():
#        if(not s.label.startswith("Z")):continue

        for ch in channels:
            src = s.label + "_" + ch + ".root"
            dst =  [s.label + "_" + ch + "_" + nSyst + versus + ".root" for nSyst in newsyst]

            if(verbose):
                for d in dst:
                    print "Copying 1 " + opt.path + src + " to "+ opt.path + d
    
            [shutil.copyfile(opt.path + src, opt.path + d) for d in dst]
    

            srcSys = s.label + "_" + ch + "_" + syst + versus + ".root"  

            if(s.label.startswith("WJets")):
                dstSys = s.label + "_" + ch + "_" + newsyst[0]+ versus + ".root" 
                if verbose: print "Copying W ++ " + opt.path + srcSys + " to "+ opt.path + dstSys
                shutil.copyfile(opt.path + srcSys, opt.path + dstSys)
            elif(s.label.startswith("ZToNuNu") or s.label.startswith("DY")):
                dstSys = s.label + "_" + ch + "_" + newsyst[1]+ versus + ".root" 
                if verbose: print "Copying Z ++ " + opt.path + srcSys + " to "+ opt.path + dstSys
                shutil.copyfile(opt.path + srcSys, opt.path + dstSys)



#syst = "VHFWeight"
#vs = ["Up", "Down"]

#newsyst = ["WHF", "ZHF"]

#createSyst("VHFWeight", ["WHF", "ZHF"], "Up")
#createSyst("VHFWeight", ["WHF", "ZHF"],"Down")
#createSyst("QCDRen", ["WQCDRen", "ZQCDRen"], "Up")
#createSyst("QCDRen", ["WQCDRen", "ZQCDRen"], "Down")
#createSyst("QCDFac", ["WQCDFac", "ZQCDFac"], "Up")
#createSyst("QCDFac", ["WQCDFac", "ZQCDFac"], "Down")
#createSyst("EWK", ["WEWK", "ZEWK"], "Up")
#createSyst("EWK", ["WEWK", "ZEWK"], "Down")
