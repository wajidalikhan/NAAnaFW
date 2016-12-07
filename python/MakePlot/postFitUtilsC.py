#from utils import *
#import ROOT
#Reads a rootfile and prints out a scale factor, the upgrade would be to print also the nuisance parameter. For now it does it with Combine only
from array import array
#import ROOT
import time

import sys
#sys.argv.append('-b')
import ROOT
import commands, os
import numpy

class postFitUtilities:
	_rootfile="./mlfit.root",
	_regions=[]
	_channel="muon"
	_sys=""
	_histodir="./histos/"
	_format="combine"

        def __init__(self,rootfile,regions, channel="muon", sys="",histodir="./histos/",format="combine"):
            self._rootfile=rootfile
            self._regions=regions
            self._channel=channel
            self._sys=sys
            self._histodir=histodir
            self._format=format
	
	def readSF(self,rootfile,  process, region, channel='muon', sys='', histodir="./histos/", format="combine"):
#	    fileFit = ROOT.TFile.Open(rootfile)
	    fileFit = ROOT.TFile.Open(str(rootfile), "READ")
	    option = ""
	    err = array('d', [1.0])
	    if region == "no_region":
	        #        print "returning 1.0 for region not in the fit."
	        return 1.0
	    if format=="combine":
	        hdirprefit= 'shapes_prefit'
                #	        print hdirprefit+'/'+channel+"_"+region+'/'+process
                hprefit = fileFit.Get(hdirprefit+'/'+channel+"_"+region+'/'+process).Clone("prefit_"+region+"_process")
	        hdirpostfit= 'shapes_fit_s'
	        hpostfit = fileFit.Get(hdirpostfit+'/'+channel+"_"+region+'/'+process).Clone("postfit_"+region+"_process")
	        neventspre = hprefit.IntegralAndError(1, hprefit.GetNbinsX(), err,option)
	        neventspost = hpostfit.IntegralAndError(1, hpostfit.GetNbinsX(), err,option)
	        sf = neventspost/neventspre
	        return sf
	
	#Tests: uncomment to check this part up to the above function:
	#sftest = readSF(rootfile="mlfit.root",process="WJets",region="2j1t")
	#print sftest
	
	
	#Prints out scale factors for region and sample:
	def listScaleFactors(self,rootfile,samples,regions,channel="muon",sys="", histodir="./histos/",format="combine"):
	    scalefactors ={}
	
	    for s in samples:
	        scalefactors[s]={}
	        for r in regions:
	            sf = self.readSF(rootfile=rootfile,process = s,region=r,channel= channel, sys=sys,histodir=histodir,format=format)
	            scalefactors[s][r]=sf
	
	    return scalefactors
	
	#Tests: uncomment to check this part up to the above function:
	#samples = ["VV","WJets","QCDMuEPt20toInf"]
	#regions = ["2j0t","2j1t"]
	#listtest = listScaleFactors(rootfile="./mlfit.root",samples=samples,regions = regions, channel="muon")
	#print listtest
	
	#Reads maps variables names to regions readable in the fit
	def mapVarsToRegions(self,variables,regions,defaultSF="no_region"):
	    mapvars = {}
	    for var in variables:
	        mapvars[var]=defaultSF
	        for r in regions:
	            if r in var:
	                mapvars[var]=r
	
	    return mapvars
	
	#This function imports the scale factors onto a format readable by the makePlot:
	def importScaleFactors(self,samples,variables,defaultSF="no_region"):
	#,regions,rootfile,channel="muon",sys="",histodir="./histos/",format="combine",defaultSF="2j0t"):
	    samplesvarsf={}
	    mapvars=self.mapVarsToRegions(variables=variables,regions=self._regions,defaultSF=defaultSF)
	    samplesregionssf={}
	    regionsnew=mapvars.values()
	    samplesregionssf=self.listScaleFactors(rootfile=self._rootfile,samples=samples, regions=regionsnew,channel=self._channel,sys=self._sys, histodir=self._histodir,format=self._format)
	    for s in samples:
	        samplesvarsf[s]={}
	        for v in variables:
	            samplesvarsf[s][v]=samplesregionssf[s][mapvars[v]]
	    return samplesvarsf
	
	
	#samples = ["VV","WJets","QCDMuEPt20toInf"]
	#regions = ["2j0t","2j1t"]
	#variables = ["2j0t_mtw","2j1t_etajprime","herpderp"]
	#listtest = importScaleFactors(variables= variables,rootfile="./mlfit.root",samples=samples,regions = regions, channel="muon")
	#print listtest
