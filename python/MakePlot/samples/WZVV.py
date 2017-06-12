######################################
#
# Annapaola de Cosa, January 2015
#
######################################
from utils import *

#---------------WJets Samples 
WToLNu0J = sample()
WToLNu0J.files = outlist (d,"WToLNu0J")
WToLNu0J.skimEff = 1.
WToLNu0J.sigma = 49670.0 
WToLNu0J.jpref = jetLabel 
WToLNu0J.jp = jetLabel
WToLNu0J.color = ROOT.kGreen+2
WToLNu0J.style = 1
WToLNu0J.fill = 1001
WToLNu0J.leglabel = "WToLNu0J"
WToLNu0J.label = "WToLNu0J"

WToLNu1J = sample()
WToLNu1J.files = outlist (d,"WToLNu1J")
WToLNu1J.skimEff = 1.
WToLNu1J.sigma = 8264.0 
WToLNu1J.jpref = jetLabel 
WToLNu1J.jp = jetLabel
WToLNu1J.color = ROOT.kGreen+2
WToLNu1J.style = 1
WToLNu1J.fill = 1001
WToLNu1J.leglabel = "WToLNu1J"
WToLNu1J.label = "WToLNu1J"

WToLNu2J = sample()
WToLNu2J.files = outlist (d,"WToLNu2J")
WToLNu2J.skimEff = 1.
WToLNu2J.sigma = 2544.0 
WToLNu2J.jpref = jetLabel 
WToLNu2J.jp = jetLabel
WToLNu2J.color = ROOT.kGreen+2
WToLNu2J.style = 1
WToLNu2J.fill = 1001
WToLNu2J.leglabel = "WToLNu2J"
WToLNu2J.label = "WToLNu2J"
#---------------WJets Samples 

#---------------ZJets Samples 
DYJetsToLL = sample()
DYJetsToLL.files = outlist (d,"DYJetsToLL")
DYJetsToLL.skimEff = 1.
DYJetsToLL.sigma =  1921.8*3.
DYJetsToLL.jpref = jetLabel 
DYJetsToLL.jp = jetLabel
DYJetsToLL.color = ROOT.kViolet
DYJetsToLL.style = 1
DYJetsToLL.fill = 1001
DYJetsToLL.leglabel = "DYJets"
DYJetsToLL.label = "DYJetsToLL"

#---------------DiBoson Samples 
WWTo1L1Nu2Q = sample()
WWTo1L1Nu2Q.files = outlist (d,"WWTo1L1Nu2Q")
WWTo1L1Nu2Q.jpref = jetLabel
WWTo1L1Nu2Q.jp = jetLabel
WWTo1L1Nu2Q.skimEff = 1.
WWTo1L1Nu2Q.sigma = 45.85
WWTo1L1Nu2Q.style = 1
WWTo1L1Nu2Q.fill = 1001
WWTo1L1Nu2Q.leglabel = "WWTo1L1Nu2Q"
WWTo1L1Nu2Q.label = "WWTo1L1Nu2Q"

WWTo2L2Nu = sample()
WWTo2L2Nu.files = outlist (d,"WWTo2L2Nu")
WWTo2L2Nu.jpref = jetLabel
WWTo2L2Nu.jp = jetLabel
WWTo2L2Nu.skimEff = 1.
WWTo2L2Nu.sigma = 12.178
WWTo2L2Nu.style = 1
WWTo2L2Nu.fill = 1001
WWTo2L2Nu.leglabel = "WWTo2L2Nu"
WWTo2L2Nu.label = "WWTo2L2Nu"

WZTo1L1Nu2Q = sample()
WZTo1L1Nu2Q.files = outlist (d,"WZTo1L1Nu2Q")
WZTo1L1Nu2Q.jpref = jetLabel
WZTo1L1Nu2Q.jp = jetLabel
WZTo1L1Nu2Q.skimEff = 1.
WZTo1L1Nu2Q.sigma = 10.71
WZTo1L1Nu2Q.style = 1
WZTo1L1Nu2Q.fill = 1001
WZTo1L1Nu2Q.leglabel = "WZTo1L1Nu2Q"
WZTo1L1Nu2Q.label = "WZTo1L1Nu2Q"

WZTo2L2Q = sample()
WZTo2L2Q.files = outlist (d,"WZTo2L2Q")
WZTo2L2Q.jpref = jetLabel 
WZTo2L2Q.jp = jetLabel 
WZTo2L2Q.skimEff = 1.
WZTo2L2Q.sigma = 5.595
WZTo2L2Q.style = 1
WZTo2L2Q.fill = 1001
WZTo2L2Q.leglabel = "WZTo2L2Q"
WZTo2L2Q.label = "WZTo2L2Q"

ZZTo2L2Q = sample()
ZZTo2L2Q.files = outlist (d,"ZZTo2L2Q")
ZZTo2L2Q.jpref = jetLabel
ZZTo2L2Q.jp = jetLabel
ZZTo2L2Q.skimEff = 1.
ZZTo2L2Q.sigma = 3.22
ZZTo2L2Q.style = 1
ZZTo2L2Q.fill = 1001
ZZTo2L2Q.leglabel = "ZZTo2L2Q"
ZZTo2L2Q.label = "ZZTo2L2Q"

#Grouping up the components: ZJets + Diboson
WZJetsVV = sample()
WZJetsVV.color = ROOT.kGreen+2 
WZJetsVV.style = 1
WZJetsVV.fill = 1001
WZJetsVV.leglabel = "W/Z-Jets + VV"
WZJetsVV.label = "WZJetsVV"
WZJetsVV.components = [WToLNu0J,WToLNu1J,WToLNu2J,DYJetsToLL,WWTo2L2Nu, WWTo1L1Nu2Q,  WZTo2L2Q, WZTo1L1Nu2Q, ZZTo2L2Q]
