from utils import *

TT = sample()
TT.files = outlist (d,"TT")
TT.files = TT.files[:1]
TT.jpref = jetLabel 
TT.jp = jetLabel 
TT.skimEff = 1.
TT.sigma = 831.76
TT.color = ROOT.kOrange
TT.style = 1
TT.fill = 1001
TT.leglabel = "t#bar{t}"
TT.label = "TT"

ST_T_tW = sample()
ST_T_tW.files = outlist (d,"ST_T_tW")
ST_T_tW.skimEff = 1.0
ST_T_tW.sigma = 35.6
ST_T_tW.jpref = jetLabel 
ST_T_tW.jp = jetLabel
ST_T_tW.color = ROOT.kYellow - 7
ST_T_tW.style = 1
ST_T_tW.fill = 1001
ST_T_tW.leglabel = "ST"
ST_T_tW.label = "ST_T_tW"

ST_Tbar_tW = sample()
ST_Tbar_tW.files = outlist (d,"ST_Tbar_tW")
ST_Tbar_tW.skimEff = 1.
ST_Tbar_tW.sigma = 35.6
ST_Tbar_tW.jpref = jetLabel 
ST_Tbar_tW.jp = jetLabel
#ST_Tbar_tW.color = ROOT.kYellow -7
ST_Tbar_tW.style = 1
ST_Tbar_tW.fill = 1001
ST_Tbar_tW.leglabel = "ST"
ST_Tbar_tW.label = "ST_Tbar_tW"

ST_T_sch = sample()
ST_T_sch.files = outlist (d,"ST_T_sch")
ST_T_sch.skimEff = 1.0
ST_T_sch.sigma = 10.32#Actually it's 6.35, 10.32 is the total top+antitop, which is represented by this sample at the moment
ST_T_sch.jpref = jetLabel 
ST_T_sch.jp = jetLabel
ST_T_sch.style = 1
ST_T_sch.fill = 1001
ST_T_sch.leglabel = "ST"
ST_T_sch.label = "ST_T_sch"

ST_Tbar_sch = sample()
ST_Tbar_sch.files = outlist (d,"ST_Tbar_sch")
ST_Tbar_sch.skimEff = 1.0
ST_Tbar_sch.sigma = 3.97
ST_Tbar_sch.jpref = jetLabel 
ST_Tbar_sch.jp = jetLabel
ST_Tbar_sch.style = 1
ST_Tbar_sch.fill = 1001
ST_Tbar_sch.leglabel = "ST"
ST_Tbar_sch.label = "ST_Tbar_sch"

TT_tWST = sample()
TT_tWST.color = ROOT.kOrange
TT_tWST.style = 1
TT_tWST.fill = 1001
TT_tWST.leglabel = "t#bar{t} + tW"
TT_tWST.label = "TT_tWST"
TT_tWST.components = [TT,ST_T_tW,ST_Tbar_tW,ST_T_sch] # summing sch+tW
