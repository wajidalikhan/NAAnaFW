######################################
#
# Annapaola de Cosa, January 2015
#
######################################


from utils import *

ST_T_tch = sample()
ST_T_tch.files = outlist (d,"ST_T_tch")
ST_T_tch.skimEff = 1
#ST_T_tch.sigma =  70.29828*4.65068
#ST_T_tch.sigma =  216.99 #* 0.324
ST_T_tch.sigma =  136.02 #* 0.324
ST_T_tch.jpref = jetLabel 
ST_T_tch.jp = jetLabel
ST_T_tch.style = 1
ST_T_tch.fill = 1001
ST_T_tch.leglabel = "ST"
ST_T_tch.label = "ST_T_tch"

ST_Tbar_tch = sample()
ST_Tbar_tch.files = outlist (d,"ST_Tbar_tch")
ST_Tbar_tch.skimEff = 1.
ST_Tbar_tch.sigma = 80.95 #* 0.324
ST_Tbar_tch.jpref = jetLabel 
ST_Tbar_tch.jp = jetLabel
ST_Tbar_tch.style = 1
ST_Tbar_tch.fill = 1001
ST_Tbar_tch.leglabel = "ST"
ST_Tbar_tch.label = "ST_Tbar_tch"

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

ST_T_tW = sample()
ST_T_tW.files = outlist (d,"ST_T_tW")
ST_T_tW.skimEff = 1.0
ST_T_tW.sigma = 35.6
ST_T_tW.jpref = jetLabel 
ST_T_tW.jp = jetLabel
#ST_T_tW.color = ROOT.kYellow - 7
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


ST_T_tch_sd = sample()
ST_T_tch_sd.files = outlist (d,"ST_T_tch_sd")
ST_T_tch_sd.skimEff = 1
#ST_T_tch_sd.sigma =  70.29828*4.65068
#ST_T_tch_sd.sigma =  216.99 #* 0.324
ST_T_tch_sd.sigma =  136.02 *0.005 #* 0.324
ST_T_tch_sd.jpref = jetLabel 
ST_T_tch_sd.jp = jetLabel
ST_T_tch_sd.style = 1
ST_T_tch_sd.fill = 1001
ST_T_tch_sd.leglabel = "ST"
ST_T_tch_sd.label = "ST_T_tch_sd"

ST_Tbar_tch_sd = sample()
ST_Tbar_tch_sd.files = outlist (d,"ST_Tbar_tch_sd")
ST_Tbar_tch_sd.skimEff = 1.
ST_Tbar_tch_sd.sigma = 80.95 *0.005 #* 0.324
ST_Tbar_tch_sd.jpref = jetLabel 
ST_Tbar_tch_sd.jp = jetLabel
ST_Tbar_tch_sd.style = 1
ST_Tbar_tch_sd.fill = 1001
ST_Tbar_tch_sd.leglabel = "ST"
ST_Tbar_tch_sd.label = "ST_Tbar_tch_sd"


#Grouping up the components:
ST_tch = sample()
ST_tch.color = ROOT.kRed
ST_tch.style = 1
ST_tch.fill = 1001
ST_tch.leglabel = "t, t-ch"
ST_tch.label = "ST_tch"
ST_tch.components = [ST_T_tch,ST_Tbar_tch]

ST_tch_sd = sample()
ST_tch_sd.color = ROOT.kBlue
ST_tch_sd.style = 1
ST_tch_sd.fill = 1001
ST_tch_sd.leglabel = "t, t-ch sd"
ST_tch_sd.label = "ST_tch_sd"
ST_tch_sd.components = [ST_T_tch_sd,ST_Tbar_tch_sd]

ST_sch = sample()
ST_sch.color = 95 
ST_sch.style = 1
ST_sch.fill = 1001
ST_sch.leglabel = "s-ch"
ST_sch.label = "ST_sch"
ST_sch.components = [ST_T_sch]

ST_tW = sample()
ST_tW.color = ROOT.kOrange+9
ST_tW.style = 1
ST_tW.fill = 1001
ST_tW.leglabel = "tW"
ST_tW.label = "ST_tW"
ST_tW.components = [ST_T_tW,ST_Tbar_tW]

ST_stW = sample()
ST_stW.color = ROOT.kOrange+9
ST_stW.style = 1
ST_stW.fill = 1001
ST_stW.leglabel = "tW+s-ch"
ST_stW.label = "ST_tW"
ST_stW.components = [ST_T_tW,ST_Tbar_tW,ST_T_sch] # summing sch+tW





