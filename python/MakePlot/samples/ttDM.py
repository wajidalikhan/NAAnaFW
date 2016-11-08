from utils import *

DMEFT_M100 = sample()
DMEFT_M100.files = outlist (d,"DMEFT_M100")
DMEFT_M100.skimEff = 1.
DMEFT_M100.sigma = 19.59 #to be changed with actula one, for now only shape comparison so no problem
DMEFT_M100.jpref = jetLabel
DMEFT_M100.jp = jetLabel
DMEFT_M100.color = ROOT.kBlack
DMEFT_M100.style = 1
DMEFT_M100.fill = 0
DMEFT_M100.leglabel = "EFT: M_{*} = 100 GeV"
DMEFT_M100.label = "DMEFT_M100"


DMtt_sc_Mchi1Mphi10 = sample()
DMtt_sc_Mchi1Mphi10.files = outlist (d,"DMtt_sc_Mchi1Mphi10")
DMtt_sc_Mchi1Mphi10.skimEff = 1.
DMtt_sc_Mchi1Mphi10.sigma = 19.59
DMtt_sc_Mchi1Mphi10.jpref = jetLabel
DMtt_sc_Mchi1Mphi10.jp = jetLabel
DMtt_sc_Mchi1Mphi10.color = ROOT.kRed
DMtt_sc_Mchi1Mphi10.style = 1
DMtt_sc_Mchi1Mphi10.fill = 0
DMtt_sc_Mchi1Mphi10.leglabel = "#splitline{sc: M_{#chi} 1 GeV,}{M_{#phi} 10 GeV}"
DMtt_sc_Mchi1Mphi10.label = "DMtt_sc_Mchi1Mphi10"

DMtt_sc_Mchi10Mphi10 = sample()
DMtt_sc_Mchi10Mphi10.files = outlist (d,"DMtt_sc_Mchi10Mphi10")
DMtt_sc_Mchi10Mphi10.skimEff = 1.
DMtt_sc_Mchi10Mphi10.sigma = 0.09487
DMtt_sc_Mchi10Mphi10.jpref = jetLabel
DMtt_sc_Mchi10Mphi10.jp = jetLabel
DMtt_sc_Mchi10Mphi10.color = ROOT.kRed
DMtt_sc_Mchi10Mphi10.style = 1
DMtt_sc_Mchi10Mphi10.fill = 0
DMtt_sc_Mchi10Mphi10.leglabel = "#splitline{sc: M_{#chi} 10 GeV,}{M_{#phi} 10 GeV}"
DMtt_sc_Mchi10Mphi10.label = "DMtt_sc_Mchi10Mphi10"

DMtt_sc_Mchi1Mphi20 = sample()
DMtt_sc_Mchi1Mphi20.files = outlist (d,"DMtt_sc_Mchi1Mphi20")
DMtt_sc_Mchi1Mphi20.skimEff = 1.
DMtt_sc_Mchi1Mphi20.sigma = 10.48
DMtt_sc_Mchi1Mphi20.jpref = jetLabel
DMtt_sc_Mchi1Mphi20.jp = jetLabel
DMtt_sc_Mchi1Mphi20.color = ROOT.kRed
DMtt_sc_Mchi1Mphi20.style = 1
DMtt_sc_Mchi1Mphi20.fill = 0
DMtt_sc_Mchi1Mphi20.leglabel = "#splitline{sc: M_{#chi} 1 GeV,}{M_{#phi} 20 GeV}"
DMtt_sc_Mchi1Mphi20.label = "DMtt_sc_Mchi1Mphi20"

DMtt_sc_Mchi1Mphi50 = sample()
DMtt_sc_Mchi1Mphi50.files = outlist (d,"DMtt_sc_Mchi1Mphi50")
DMtt_sc_Mchi1Mphi50.skimEff = 1.
DMtt_sc_Mchi1Mphi50.sigma = 2.941
DMtt_sc_Mchi1Mphi50.jpref = jetLabel
DMtt_sc_Mchi1Mphi50.jp = jetLabel
DMtt_sc_Mchi1Mphi50.color = ROOT.kBlue
DMtt_sc_Mchi1Mphi50.style = ROOT.kDashed
DMtt_sc_Mchi1Mphi50.fill = 0
DMtt_sc_Mchi1Mphi50.leglabel = "M_{#chi} = 1 GeV, M_{#phi} = 50 GeV"
DMtt_sc_Mchi1Mphi50.label = "DMtt_sc_Mchi1Mphi50"

DMtt_sc_Mchi10Mphi50 = sample()
DMtt_sc_Mchi10Mphi50.files = outlist (d,"DMtt_sc_Mchi10Mphi50")
DMtt_sc_Mchi10Mphi50.skimEff = 1.
DMtt_sc_Mchi10Mphi50.sigma = 2.942
DMtt_sc_Mchi10Mphi50.jpref = jetLabel
DMtt_sc_Mchi10Mphi50.jp = jetLabel
DMtt_sc_Mchi10Mphi50.color = ROOT.kBlue
DMtt_sc_Mchi10Mphi50.style = ROOT.kDashed
DMtt_sc_Mchi10Mphi50.fill = 0
DMtt_sc_Mchi10Mphi50.leglabel = "M_{#chi} = 10 GeV, M_{#phi} = 50 GeV"
DMtt_sc_Mchi10Mphi50.label = "DMtt_sc_Mchi10Mphi50"

DMtt_sc_Mchi50Mphi50 = sample()
DMtt_sc_Mchi50Mphi50.files = outlist (d,"DMtt_sc_Mchi50Mphi50")
DMtt_sc_Mchi50Mphi50.skimEff = 1.
DMtt_sc_Mchi50Mphi50.sigma = 0.002329
DMtt_sc_Mchi50Mphi50.jpref = jetLabel
DMtt_sc_Mchi50Mphi50.jp = jetLabel
DMtt_sc_Mchi50Mphi50.color = ROOT.kBlue
DMtt_sc_Mchi50Mphi50.style = ROOT.kDashed
DMtt_sc_Mchi50Mphi50.fill = 0
DMtt_sc_Mchi50Mphi50.leglabel = "M_{#chi} = 50 GeV, M_{#phi} = 50 GeV"
DMtt_sc_Mchi50Mphi50.label = "DMtt_sc_Mchi50Mphi50"

DMtt_sc_Mchi1Mphi100 = sample()
DMtt_sc_Mchi1Mphi100.files = outlist (d,"DMtt_sc_Mchi1Mphi100")
DMtt_sc_Mchi1Mphi100.skimEff = 1.
DMtt_sc_Mchi1Mphi100.sigma = 0.6723
DMtt_sc_Mchi1Mphi100.jpref = jetLabel
DMtt_sc_Mchi1Mphi100.jp = jetLabel
DMtt_sc_Mchi1Mphi100.color = ROOT.kRed
#DMtt_sc_Mchi1Mphi100.style = ROOT.kDotted
DMtt_sc_Mchi1Mphi100.style = 1
DMtt_sc_Mchi1Mphi100.fill = 0
DMtt_sc_Mchi1Mphi100.leglabel = "#splitline{sc: M_{#chi} 1 GeV,}{M_{#phi} 100 GeV}"
DMtt_sc_Mchi1Mphi100.label = "DMtt_sc_Mchi1Mphi100"

DMtt_sc_Mchi10Mphi100 = sample()
DMtt_sc_Mchi10Mphi100.files = outlist (d,"DMtt_sc_Mchi10Mphi100")
DMtt_sc_Mchi10Mphi100.skimEff = 1.
DMtt_sc_Mchi10Mphi100.sigma = 0.6732
DMtt_sc_Mchi10Mphi100.jpref = jetLabel
DMtt_sc_Mchi10Mphi100.jp = jetLabel
DMtt_sc_Mchi10Mphi100.color = ROOT.kRed
DMtt_sc_Mchi10Mphi100.style = ROOT.kDotted
DMtt_sc_Mchi10Mphi100.fill = 0
DMtt_sc_Mchi10Mphi100.leglabel = "#splitline{sc: M_{#chi} 10 GeV,}{M_{#phi} 100 GeV}"
DMtt_sc_Mchi10Mphi100.label = "DMtt_sc_Mchi10Mphi100"

DMtt_sc_Mchi1Mphi200 = sample()
DMtt_sc_Mchi1Mphi200.files = outlist (d,"DMtt_sc_Mchi1Mphi200")
DMtt_sc_Mchi1Mphi200.skimEff = 1.
DMtt_sc_Mchi1Mphi200.sigma = 0.09327
DMtt_sc_Mchi1Mphi200.jpref = jetLabel
DMtt_sc_Mchi1Mphi200.jp = jetLabel
DMtt_sc_Mchi1Mphi200.color = ROOT.kRed
DMtt_sc_Mchi1Mphi200.style = ROOT.kDotted
DMtt_sc_Mchi1Mphi200.fill = 0
DMtt_sc_Mchi1Mphi200.leglabel = "#splitline{sc: M_{#chi} 1 GeV,}{M_{#phi} 200 GeV}"
DMtt_sc_Mchi1Mphi200.label = "DMtt_sc_Mchi1Mphi200"

DMtt_sc_Mchi50Mphi200 = sample()
DMtt_sc_Mchi50Mphi200.files = outlist (d,"DMtt_sc_Mchi50Mphi200")
DMtt_sc_Mchi50Mphi200.skimEff = 1.
DMtt_sc_Mchi50Mphi200.sigma = 0.09224
DMtt_sc_Mchi50Mphi200.jpref = jetLabel
DMtt_sc_Mchi50Mphi200.jp = jetLabel
DMtt_sc_Mchi50Mphi200.color = ROOT.kRed
DMtt_sc_Mchi50Mphi200.style = ROOT.kDotted
DMtt_sc_Mchi50Mphi200.fill = 0
DMtt_sc_Mchi50Mphi200.leglabel = "#splitline{sc: M_{#chi} 50 GeV,}{M_{#phi} 200 GeV}"
DMtt_sc_Mchi50Mphi200.label = "DMtt_sc_Mchi50Mphi200"

DMtt_sc_Mchi150Mphi200 = sample()
DMtt_sc_Mchi150Mphi200.files = outlist (d,"DMtt_sc_Mchi150Mphi200")
DMtt_sc_Mchi150Mphi200.skimEff = 1.
DMtt_sc_Mchi150Mphi200.sigma = 0.00013
DMtt_sc_Mchi150Mphi200.jpref = jetLabel
DMtt_sc_Mchi150Mphi200.jp = jetLabel
DMtt_sc_Mchi150Mphi200.color = ROOT.kRed
DMtt_sc_Mchi150Mphi200.style = ROOT.kDotted
DMtt_sc_Mchi150Mphi200.fill = 0
DMtt_sc_Mchi150Mphi200.leglabel = "#splitline{sc: M_{#chi} 150 GeV,}{M_{#phi} 200 GeV}"
DMtt_sc_Mchi150Mphi200.label = "DMtt_sc_Mchi150Mphi200"

DMtt_sc_Mchi1Mphi300 = sample()
DMtt_sc_Mchi1Mphi300.files = outlist (d,"DMtt_sc_Mchi1Mphi300")
DMtt_sc_Mchi1Mphi300.skimEff = 1.
DMtt_sc_Mchi1Mphi300.sigma = 0.0295
DMtt_sc_Mchi1Mphi300.jpref = jetLabel
DMtt_sc_Mchi1Mphi300.jp = jetLabel
DMtt_sc_Mchi1Mphi300.color = ROOT.kRed
DMtt_sc_Mchi1Mphi300.style = ROOT.kDotted
DMtt_sc_Mchi1Mphi300.fill = 0
DMtt_sc_Mchi1Mphi300.leglabel = "#splitline{scalar: M_{#chi} 1 GeV,}{M_{#phi} 300 GeV}"
DMtt_sc_Mchi1Mphi300.label = "DMtt_sc_Mchi1Mphi300"

DMtt_sc_Mchi50Mphi300 = sample()
DMtt_sc_Mchi50Mphi300.files = outlist (d,"DMtt_sc_Mchi50Mphi300")
DMtt_sc_Mchi50Mphi300.skimEff = 1.
DMtt_sc_Mchi50Mphi300.sigma = 0.02901
DMtt_sc_Mchi50Mphi300.jpref = jetLabel
DMtt_sc_Mchi50Mphi300.jp = jetLabel
DMtt_sc_Mchi50Mphi300.color = ROOT.kRed
#DMtt_sc_Mchi50Mphi300.style = ROOT.kDotted
DMtt_sc_Mchi50Mphi300.style = ROOT.kDashed
DMtt_sc_Mchi50Mphi300.fill = 0
DMtt_sc_Mchi50Mphi300.leglabel = "#splitline{sc: M_{#chi} 50 GeV,}{M_{#phi} 300 GeV}"
DMtt_sc_Mchi50Mphi300.label = "DMtt_sc_Mchi50Mphi300"

DMtt_sc_Mchi1Mphi500 = sample()
DMtt_sc_Mchi1Mphi500.files = outlist (d,"DMtt_sc_Mchi1Mphi500")
DMtt_sc_Mchi1Mphi500.skimEff = 1.
DMtt_sc_Mchi1Mphi500.sigma = 0.00518
DMtt_sc_Mchi1Mphi500.jpref = jetLabel
DMtt_sc_Mchi1Mphi500.jp = jetLabel
DMtt_sc_Mchi1Mphi500.color = ROOT.kRed
DMtt_sc_Mchi1Mphi500.style = ROOT.kDotted
DMtt_sc_Mchi1Mphi500.fill = 0
DMtt_sc_Mchi1Mphi500.leglabel = "#splitline{sc: M_{#chi} 1 GeV,}{M_{#phi} 500 GeV}"
DMtt_sc_Mchi1Mphi500.label = "DMtt_sc_Mchi1Mphi500"

DMtt_sc_Mchi150Mphi500 = sample()
DMtt_sc_Mchi150Mphi500.files = outlist (d,"DMtt_sc_Mchi150Mphi500")
DMtt_sc_Mchi150Mphi500.skimEff = 1.
DMtt_sc_Mchi150Mphi500.sigma = 0.003754
DMtt_sc_Mchi150Mphi500.jpref = jetLabel
DMtt_sc_Mchi150Mphi500.jp = jetLabel
DMtt_sc_Mchi150Mphi500.color = ROOT.kRed
DMtt_sc_Mchi150Mphi500.style = ROOT.kDotted
DMtt_sc_Mchi150Mphi500.fill = 0
DMtt_sc_Mchi150Mphi500.leglabel = "#splitline{sc: M_{#chi} 150 GeV,}{M_{#phi} 500 GeV}"
DMtt_sc_Mchi150Mphi500.label = "DMtt_sc_Mchi150Mphi500"

DMtt_sc_Mchi500Mphi500 = sample()
DMtt_sc_Mchi500Mphi500.files = outlist (d,"DMtt_sc_Mchi500Mphi500")
DMtt_sc_Mchi500Mphi500.skimEff = 1.
DMtt_sc_Mchi500Mphi500.sigma = 0.0000009894
DMtt_sc_Mchi500Mphi500.jpref = jetLabel
DMtt_sc_Mchi500Mphi500.jp = jetLabel
DMtt_sc_Mchi500Mphi500.color = ROOT.kRed
DMtt_sc_Mchi500Mphi500.style = ROOT.kDotted
DMtt_sc_Mchi500Mphi500.fill = 0
DMtt_sc_Mchi500Mphi500.leglabel = "#splitline{sc: M_{#chi} 500 GeV,}{M_{#phi} 200 GeV}"
DMtt_sc_Mchi500Mphi500.label = "DMtt_sc_Mchi500Mphi500"

DMtt_sc_Mchi1Mphi1000 = sample()
DMtt_sc_Mchi1Mphi1000.files = outlist (d,"DMtt_sc_Mchi1Mphi1000")
DMtt_sc_Mchi1Mphi1000.skimEff = 1.
DMtt_sc_Mchi1Mphi1000.sigma = 0.0003687
DMtt_sc_Mchi1Mphi1000.jpref = jetLabel
DMtt_sc_Mchi1Mphi1000.jp = jetLabel
DMtt_sc_Mchi1Mphi1000.color = ROOT.kRed
DMtt_sc_Mchi1Mphi1000.style = ROOT.kDashed
DMtt_sc_Mchi1Mphi1000.fill = 0
DMtt_sc_Mchi1Mphi1000.leglabel = "#splitline{sc: M_{#chi} 1 GeV,}{M_{#phi} 1000 GeV}"
DMtt_sc_Mchi1Mphi1000.label = "DMtt_sc_Mchi1Mphi1000"

DMtt_sc_Mchi150Mphi1000 = sample()
DMtt_sc_Mchi150Mphi1000.files = outlist (d,"DMtt_sc_Mchi150Mphi1000")
DMtt_sc_Mchi150Mphi1000.skimEff = 1.
DMtt_sc_Mchi150Mphi1000.sigma = 0. #?
DMtt_sc_Mchi150Mphi1000.jpref = jetLabel
DMtt_sc_Mchi150Mphi1000.jp = jetLabel
DMtt_sc_Mchi150Mphi1000.color = ROOT.kRed
DMtt_sc_Mchi150Mphi1000.style = ROOT.kDashed
DMtt_sc_Mchi150Mphi1000.fill = 0
DMtt_sc_Mchi150Mphi1000.leglabel = "#splitline{sc: M_{#chi} 150 GeV,}{M_{#phi} 1000 GeV}"
DMtt_sc_Mchi150Mphi1000.label = "DMtt_sc_Mchi150Mphi1000"


#pseudo-scalar

DMtt_ps_Mchi1Mphi10 = sample()
DMtt_ps_Mchi1Mphi10.files = outlist (d,"DMtt_ps_Mchi1Mphi10")
DMtt_ps_Mchi1Mphi10.skimEff = 1.
DMtt_ps_Mchi1Mphi10.sigma = 0.4409
DMtt_ps_Mchi1Mphi10.jpref = jetLabel
DMtt_ps_Mchi1Mphi10.jp = jetLabel
DMtt_ps_Mchi1Mphi10.color = ROOT.kRed
DMtt_ps_Mchi1Mphi10.style = 1
DMtt_ps_Mchi1Mphi10.fill = 0
DMtt_ps_Mchi1Mphi10.leglabel = "#splitline{#splitline{pseudo-scalar:}{ M_{#chi} 1 GeV,}}{M_{#phi} 10 GeV}"
DMtt_ps_Mchi1Mphi10.label = "DMtt_ps_Mchi1Mphi10"

DMtt_ps_Mchi10Mphi10 = sample()
DMtt_ps_Mchi10Mphi10.files = outlist (d,"DMtt_ps_Mchi10Mphi10")
DMtt_ps_Mchi10Mphi10.skimEff = 1.
DMtt_ps_Mchi10Mphi10.sigma = 0.01499
DMtt_ps_Mchi10Mphi10.jpref = jetLabel
DMtt_ps_Mchi10Mphi10.jp = jetLabel
DMtt_ps_Mchi10Mphi10.color = ROOT.kRed
DMtt_ps_Mchi10Mphi10.style = 1
DMtt_ps_Mchi10Mphi10.fill = 0
DMtt_ps_Mchi10Mphi10.leglabel = "#splitline{ps: M_{#chi} 10 GeV,}{M_{#phi} 10 GeV}"
DMtt_ps_Mchi10Mphi10.label = "DMtt_ps_Mchi10Mphi10"

DMtt_ps_Mchi1Mphi20 = sample()
DMtt_ps_Mchi1Mphi20.files = outlist (d,"DMtt_ps_Mchi1Mphi20")
DMtt_ps_Mchi1Mphi20.skimEff = 1.
DMtt_ps_Mchi1Mphi20.sigma = 0.3992
DMtt_ps_Mchi1Mphi20.jpref = jetLabel
DMtt_ps_Mchi1Mphi20.jp = jetLabel
DMtt_ps_Mchi1Mphi20.color = ROOT.kRed
DMtt_ps_Mchi1Mphi20.style = 1
DMtt_ps_Mchi1Mphi20.fill = 0
DMtt_ps_Mchi1Mphi20.leglabel = "#splitline{ps: M_{#chi} 1 GeV,}{M_{#phi} 20 GeV}"
DMtt_ps_Mchi1Mphi20.label = "DMtt_ps_Mchi1Mphi20"

DMtt_ps_Mchi1Mphi50 = sample()
DMtt_ps_Mchi1Mphi50.files = outlist (d,"DMtt_ps_Mchi1Mphi50")
DMtt_ps_Mchi1Mphi50.skimEff = 1.
DMtt_ps_Mchi1Mphi50.sigma = 0.3032
DMtt_ps_Mchi1Mphi50.jpref = jetLabel
DMtt_ps_Mchi1Mphi50.jp = jetLabel
DMtt_ps_Mchi1Mphi50.color = ROOT.kBlue
DMtt_ps_Mchi1Mphi50.style = ROOT.kDashed
DMtt_ps_Mchi1Mphi50.fill = 0
DMtt_ps_Mchi1Mphi50.leglabel = "M_{#chi} = 1 GeV, M_{#phi} = 50 GeV"
DMtt_ps_Mchi1Mphi50.label = "DMtt_ps_Mchi1Mphi50"

DMtt_ps_Mchi10Mphi50 = sample()
DMtt_ps_Mchi10Mphi50.files = outlist (d,"DMtt_ps_Mchi10Mphi50")
DMtt_ps_Mchi10Mphi50.skimEff = 1.
DMtt_ps_Mchi10Mphi50.sigma = 0.3034
DMtt_ps_Mchi10Mphi50.jpref = jetLabel
DMtt_ps_Mchi10Mphi50.jp = jetLabel
DMtt_ps_Mchi10Mphi50.color = ROOT.kBlue
DMtt_ps_Mchi10Mphi50.style = ROOT.kDashed
DMtt_ps_Mchi10Mphi50.fill = 0
DMtt_ps_Mchi10Mphi50.leglabel = "M_{#chi} = 10 GeV, M_{#phi} = 50 GeV"
DMtt_ps_Mchi10Mphi50.label = "DMtt_ps_Mchi10Mphi50"

DMtt_ps_Mchi50Mphi50 = sample()
DMtt_ps_Mchi50Mphi50.files = outlist (d,"DMtt_ps_Mchi50Mphi50")
DMtt_ps_Mchi50Mphi50.skimEff = 1.
DMtt_ps_Mchi50Mphi50.sigma = 0.002979
DMtt_ps_Mchi50Mphi50.jpref = jetLabel
DMtt_ps_Mchi50Mphi50.jp = jetLabel
DMtt_ps_Mchi50Mphi50.color = ROOT.kBlue
DMtt_ps_Mchi50Mphi50.style = ROOT.kDashed
DMtt_ps_Mchi50Mphi50.fill = 0
DMtt_ps_Mchi50Mphi50.leglabel = "M_{#chi} = 50 GeV, M_{#phi} = 50 GeV"
DMtt_ps_Mchi50Mphi50.label = "DMtt_ps_Mchi50Mphi50"

DMtt_ps_Mchi1Mphi100 = sample()
DMtt_ps_Mchi1Mphi100.files = outlist (d,"DMtt_ps_Mchi1Mphi100")
DMtt_ps_Mchi1Mphi100.skimEff = 1.
DMtt_ps_Mchi1Mphi100.sigma = 0.1909
DMtt_ps_Mchi1Mphi100.jpref = jetLabel
DMtt_ps_Mchi1Mphi100.jp = jetLabel
DMtt_ps_Mchi1Mphi100.color = ROOT.kRed+1
DMtt_ps_Mchi1Mphi100.style = ROOT.kSolid
DMtt_ps_Mchi1Mphi100.fill = 0
DMtt_ps_Mchi1Mphi100.leglabel = "#splitline{pseudoscalar:}{m_{MED}=100 GeV, m_{DM} 1 GeV}"
DMtt_ps_Mchi1Mphi100.label = "DMtt_ps_Mchi1Mphi100"

DMtt_ps_Mchi10Mphi100 = sample()
DMtt_ps_Mchi10Mphi100.files = outlist (d,"DMtt_ps_Mchi10Mphi100")
DMtt_ps_Mchi10Mphi100.skimEff = 1.
DMtt_ps_Mchi10Mphi100.sigma = 0.1901
DMtt_ps_Mchi10Mphi100.jpref = jetLabel
DMtt_ps_Mchi10Mphi100.jp = jetLabel
DMtt_ps_Mchi10Mphi100.color = ROOT.kRed
DMtt_ps_Mchi10Mphi100.style = ROOT.kDotted
DMtt_ps_Mchi10Mphi100.fill = 0
DMtt_ps_Mchi10Mphi100.leglabel = "#splitline{ps: M_{#chi} 10 GeV,}{M_{#phi} 100 GeV}"
DMtt_ps_Mchi10Mphi100.label = "DMtt_ps_Mchi10Mphi100"

DMtt_ps_Mchi1Mphi200 = sample()
DMtt_ps_Mchi1Mphi200.files = outlist (d,"DMtt_ps_Mchi1Mphi200")
DMtt_ps_Mchi1Mphi200.skimEff = 1.
DMtt_ps_Mchi1Mphi200.sigma = 0.0836
DMtt_ps_Mchi1Mphi200.jpref = jetLabel
DMtt_ps_Mchi1Mphi200.jp = jetLabel
DMtt_ps_Mchi1Mphi200.color = ROOT.kRed
DMtt_ps_Mchi1Mphi200.style = ROOT.kDotted
DMtt_ps_Mchi1Mphi200.fill = 0
DMtt_ps_Mchi1Mphi200.leglabel = "#splitline{ps: M_{#chi} 1 GeV,}{M_{#phi} 200 GeV}"
DMtt_ps_Mchi1Mphi200.label = "DMtt_ps_Mchi1Mphi200"

DMtt_ps_Mchi50Mphi200 = sample()
DMtt_ps_Mchi50Mphi200.files = outlist (d,"DMtt_ps_Mchi50Mphi200")
DMtt_ps_Mchi50Mphi200.skimEff = 1.
DMtt_ps_Mchi50Mphi200.sigma = 0.08382
DMtt_ps_Mchi50Mphi200.jpref = jetLabel
DMtt_ps_Mchi50Mphi200.jp = jetLabel
DMtt_ps_Mchi50Mphi200.color = ROOT.kRed
DMtt_ps_Mchi50Mphi200.style = ROOT.kDotted
DMtt_ps_Mchi50Mphi200.fill = 0
DMtt_ps_Mchi50Mphi200.leglabel = "#splitline{ps: M_{#chi} 50 GeV,}{M_{#phi} 200 GeV}"
DMtt_ps_Mchi50Mphi200.label = "DMtt_ps_Mchi50Mphi200"

DMtt_ps_Mchi150Mphi200 = sample()
DMtt_ps_Mchi150Mphi200.files = outlist (d,"DMtt_ps_Mchi150Mphi200")
DMtt_ps_Mchi150Mphi200.skimEff = 1.
DMtt_ps_Mchi150Mphi200.sigma = 0.0004124
DMtt_ps_Mchi150Mphi200.jpref = jetLabel
DMtt_ps_Mchi150Mphi200.jp = jetLabel
DMtt_ps_Mchi150Mphi200.color = ROOT.kRed
DMtt_ps_Mchi150Mphi200.style = ROOT.kDotted
DMtt_ps_Mchi150Mphi200.fill = 0
DMtt_ps_Mchi150Mphi200.leglabel = "#splitline{ps: M_{#chi} 150 GeV,}{M_{#phi} 200 GeV}"
DMtt_ps_Mchi150Mphi200.label = "DMtt_ps_Mchi150Mphi200"


DMtt_ps_Mchi1Mphi300 = sample()
DMtt_ps_Mchi1Mphi300.files = outlist (d,"DMtt_ps_Mchi1Mphi300")
DMtt_ps_Mchi1Mphi300.skimEff = 1.
DMtt_ps_Mchi1Mphi300.sigma = 0.03999
DMtt_ps_Mchi1Mphi300.jpref = jetLabel
DMtt_ps_Mchi1Mphi300.jp = jetLabel
DMtt_ps_Mchi1Mphi300.color = ROOT.kRed
DMtt_ps_Mchi1Mphi300.style = ROOT.kDotted
DMtt_ps_Mchi1Mphi300.fill = 0
DMtt_ps_Mchi1Mphi300.leglabel = "#splitline{ps: M_{#chi} 1 GeV,}{M_{#phi} 300 GeV}"
DMtt_ps_Mchi1Mphi300.label = "DMtt_ps_Mchi1Mphi300"

DMtt_ps_Mchi50Mphi300 = sample()
DMtt_ps_Mchi50Mphi300.files = outlist (d,"DMtt_ps_Mchi50Mphi300")
DMtt_ps_Mchi50Mphi300.skimEff = 1.
DMtt_ps_Mchi50Mphi300.sigma = 0.03989
DMtt_ps_Mchi50Mphi300.jpref = jetLabel
DMtt_ps_Mchi50Mphi300.jp = jetLabel
DMtt_ps_Mchi50Mphi300.color = ROOT.kRed
DMtt_ps_Mchi50Mphi300.style = ROOT.kDotted
DMtt_ps_Mchi50Mphi300.fill = 0
DMtt_ps_Mchi50Mphi300.leglabel = "#splitline{ps: M_{#chi} 50 GeV,}{M_{#phi} 300 GeV}"
DMtt_ps_Mchi50Mphi300.label = "DMtt_ps_Mchi50Mphi300"

DMtt_ps_Mchi1Mphi500 = sample()
DMtt_ps_Mchi1Mphi500.files = outlist (d,"DMtt_ps_Mchi1Mphi500")
DMtt_ps_Mchi1Mphi500.skimEff = 1.
DMtt_ps_Mchi1Mphi500.sigma = 0.005408
DMtt_ps_Mchi1Mphi500.jpref = jetLabel
DMtt_ps_Mchi1Mphi500.jp = jetLabel
DMtt_ps_Mchi1Mphi500.color = ROOT.kRed
DMtt_ps_Mchi1Mphi500.style = ROOT.kDotted
DMtt_ps_Mchi1Mphi500.fill = 0
DMtt_ps_Mchi1Mphi500.leglabel = "#splitline{ps: M_{#chi} 1 GeV,}{M_{#phi} 500 GeV}"
DMtt_ps_Mchi1Mphi500.label = "DMtt_ps_Mchi1Mphi500"

DMtt_ps_Mchi150Mphi500 = sample()
DMtt_ps_Mchi150Mphi500.files = outlist (d,"DMtt_ps_Mchi150Mphi500")
DMtt_ps_Mchi150Mphi500.skimEff = 1.
DMtt_ps_Mchi150Mphi500.sigma = 0.004611
DMtt_ps_Mchi150Mphi500.jpref = jetLabel
DMtt_ps_Mchi150Mphi500.jp = jetLabel
DMtt_ps_Mchi150Mphi500.color = ROOT.kRed
DMtt_ps_Mchi150Mphi500.style = ROOT.kDotted
DMtt_ps_Mchi150Mphi500.fill = 0
DMtt_ps_Mchi150Mphi500.leglabel = "#splitline{ps: M_{#chi} 150 GeV,}{M_{#phi} 500 GeV}"
DMtt_ps_Mchi150Mphi500.label = "DMtt_ps_Mchi150Mphi500"

DMtt_ps_Mchi500Mphi500 = sample()
DMtt_ps_Mchi500Mphi500.files = outlist (d,"DMtt_ps_Mchi500Mphi500")
DMtt_ps_Mchi500Mphi500.skimEff = 1.
DMtt_ps_Mchi500Mphi500.sigma = 0.000003275 
DMtt_ps_Mchi500Mphi500.jpref = jetLabel
DMtt_ps_Mchi500Mphi500.jp = jetLabel
DMtt_ps_Mchi500Mphi500.color = ROOT.kRed
DMtt_ps_Mchi500Mphi500.style = ROOT.kDotted
DMtt_ps_Mchi500Mphi500.fill = 0
DMtt_ps_Mchi500Mphi500.leglabel = "#splitline{ps: M_{#chi} 500 GeV,}{M_{#phi} 200 GeV}"
DMtt_ps_Mchi500Mphi500.label = "DMtt_ps_Mchi500Mphi500"

DMtt_ps_Mchi150Mphi1000 = sample()
DMtt_ps_Mchi150Mphi1000.files = outlist (d,"DMtt_ps_Mchi150Mphi1000")
DMtt_ps_Mchi150Mphi1000.skimEff = 1.
DMtt_ps_Mchi150Mphi1000.sigma = 0. #?
DMtt_ps_Mchi150Mphi1000.jp = jetLabel
DMtt_ps_Mchi150Mphi1000.color = ROOT.kRed
DMtt_ps_Mchi150Mphi1000.style = ROOT.kDashed
DMtt_ps_Mchi150Mphi1000.fill = 0
DMtt_ps_Mchi150Mphi1000.leglabel = "#splitline{ps: M_{#chi} 150 GeV,}{M_{#phi} 1000 GeV}"
DMtt_ps_Mchi150Mphi1000.label = "DMtt_ps_Mchi150Mphi1000"




ttDM1 = sample()
ttDM1.files = outlist (d,"ttDM1")
ttDM1.skimEff = 1.
ttDM1.sigma = 1.32
ttDM1.jpref = jetLabel 
ttDM1.jp = jetLabel 
ttDM1.color = ROOT.kRed
ttDM1.style = 1
ttDM1.fill = 0
ttDM1.leglabel = "M_{#chi} = 1 GeV"
ttDM1.label = "ttDM1"

ttDM10 = sample()
ttDM10.files = outlist (d,"ttDM10")
ttDM10.skimEff = 1.
ttDM10.sigma = 1.322
ttDM10.jpref = jetLabel 
ttDM10.jp = jetLabel 
ttDM10.sf = 0.006651036898
ttDM10.color = ROOT.kRed
ttDM10.style = 4
ttDM10.fill = 0
ttDM10.leglabel = "M_{#chi} = 10 GeV"
ttDM10.label = "ttDM10"

ttDM50 = sample()
ttDM50.files = outlist (d,"ttDM50")
ttDM50.skimEff = 1.
ttDM50.sigma = 1.187
ttDM50.jpref = jetLabel 
ttDM50.jp = jetLabel 
ttDM50.color = ROOT.kRed
ttDM50.style = 9
ttDM50.fill = 0
ttDM50.leglabel = "M_{#chi} = 50 GeV"
ttDM50.label = "ttDM50"

ttDM100 = sample()
ttDM100.files = outlist (d,"ttDM100")
ttDM100.jpref = jetLabel 
ttDM100.jp = jetLabel 
ttDM100.skimEff = 1.
ttDM100.sigma = 0.9553
ttDM100.color = ROOT.kRed
ttDM100.style = 2
ttDM100.fill = 0
ttDM100.leglabel = "M_{#chi} = 100 GeV"
ttDM100.label = "ttDM100"

ttDM200 = sample()
ttDM200.files = outlist (d,"ttDM200")
ttDM200.jpref = jetLabel 
ttDM200.jp = jetLabel 
ttDM200.skimEff = 1.
ttDM200.sigma = 0.6301
ttDM200.color = ROOT.kRed
ttDM200.style = 8
ttDM200.fill = 0
ttDM200.leglabel = "M_{#chi} = 200 GeV"
ttDM200.label = "ttDM200"

ttDM600 = sample()
ttDM600.files = outlist (d,"ttDM600")
ttDM600.skimEff = 1.
ttDM600.sigma = 0.1038
ttDM600.jpref = jetLabel 
ttDM600.jp = jetLabel 
ttDM600.color = ROOT.kRed
ttDM600.style = 5
ttDM600.fill = 0
ttDM600.leglabel = "M_{#chi} = 600 GeV"
ttDM600.label = "ttDM600"

ttDM1000 = sample()
ttDM1000.files = outlist (d,"ttDM1000")
ttDM1000.skimEff = 1.
ttDM1000.sigma = 0.01585
ttDM1000.jpref = jetLabel 
ttDM1000.jp = jetLabel 
ttDM1000.color = ROOT.kRed
ttDM1000.style = 3
ttDM1000.fill = 0
ttDM1000.leglabel = "M_{#chi} = 1000 GeV"
ttDM1000.label = "ttDM1000"
