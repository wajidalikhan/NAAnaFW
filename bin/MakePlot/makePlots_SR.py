#!/bin/env python
#########################################################
############ Macro for plots with ratio pad #############
#########################################################

import optparse

# Partse command line before ANY action
usage = 'usage: %prog -l lumi'
parser = optparse.OptionParser(usage)

parser.add_option('-l', '--lumi', dest='lumi', type='float', default = '0.20923', help='Luminosity')
parser.add_option('-c', '--channel', dest='channel', type='string', default = 'semileptonic', help='Channel to analyze: semileptonic or fullhadronic')
parser.add_option('-s', '--sys', dest='sys', type='string', default = 'noSys', help='Systematics: noSys, jesUp, jesDown')
parser.add_option('-n', '--normData', dest='normData', type='int', default = '0', help='Normalise to data?')
parser.add_option('-r', '--resdir', dest='resdir', type='string', default = './', help='res directory')
parser.add_option('-d', '--no-data', dest='data', action='store_false', default = True, help='Hide data')

(opt, args) = parser.parse_args()

import time
import sys
import ROOT
import copy
import commands, os.path
import numpy as n

# from plots.services import Histo, Stack, Legend, deltaPhi, Histo1
from plots.services import deltaPhi
from samples.toPlot import samples
import plots.common_SR, plots.fullhadronic_SR, plots.semileptonic_SR


import tdrstyle, CMS_lumi


def sumerrors(h):
    return sum([h.GetBinError(ib) for ib in xrange(0,h.GetNbinsX()+2)])

def writeSummary(label, channel, sys, h, lumi):
    if sys=="noSys":
        filename = outtxt + '/' + label + "_" + channel + ".txt"
    else: 
        filename = outtxt + '/' + label + "_" + channel + "_" + sys + ".txt"

    ofile = open( filename, 'w')
    #print "Opening file: ", filename
    ofile.write("Lumi: %.1f\n"%(lumi))
    ofile.write("---> %s\n"%(label))
    ofile.write("Events before selection   : %.2f\n" %(h.GetBinContent(0)) )
    ofile.write("Events after trigger      : %.2f\n" %(h.GetBinContent(1)) )
    ofile.write("Events after lepton cut   : %.2f\n" %(h.GetBinContent(2)) )
    ofile.write("Events after jet cut      : %.2f\n" %(h.GetBinContent(3)) )
    ofile.write("Events after bjet cut     : %.2f\n" %(h.GetBinContent(4)) )
    ofile.write("Events after met160 cut   : %.2f\n" %(h.GetBinContent(5)) )
    ofile.write("Events after met320 cut   : %.2f\n" %(h.GetBinContent(6)) )
    if(channel == "semileptonic"):
        ofile.write("Events after mt cut       : %.2f\n" %(h.GetBinContent(7)) )
        ofile.write("Events after mt2w cut     : %.2f\n" %(h.GetBinContent(8)) )
        ofile.write("Events after minDPhi cut  : %.2f\n" %(h.GetBinContent(9)) )
        if(h.GetBinContent(0)!=0):ofile.write("Selection efficiency      : %.2f\n" %(h.GetBinContent(9)/h.GetBinContent(0)) )
    elif(channel == "fullhadronic"):
        ofile.write("Events after minDPhi cut  : %.2f\n" %(h.GetBinContent(7)) )
        if(h.GetBinContent(0)!=0):ofile.write("Selection efficiency      : %.2f\n" %(h.GetBinContent(7)/h.GetBinContent(0)) )


ROOT.gROOT.Reset();
ROOT.gROOT.SetStyle('Plain')
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.TH1.SetDefaultSumw2()
ROOT.TH1.AddDirectory(False)

tdrstyle.setTDRStyle();

settings = {}
store = []

settings.update(plots.common_SR.settings)
store += plots.common_SR.store

if opt.channel == 'semileptonic':
    # Add semileptonic specific settings and plots to store
    settings.update(plots.semileptonic_SR.settings)
    store += plots.semileptonic_SR.store

    # Define sl output paths
    outhistos = 'output/sl/histos_SR'
    outpdfs = 'output/sl/pdfs_SR'
    outtxt = 'output/sl/txt_SR'

elif opt.channel == 'fullhadronic':
    # Add fullhadronic specific settings and plots to store
    settings.update(plots.fullhadronic_SR.settings)
    store += plots.fullhadronic_SR.store

    # Define fh output paths
    outhistos = 'output/fh/histos_SR'
    outpdfs = 'output/fh/pdfs_SR'
    outtxt = 'output/fh/txt_SR'

else:
    print 'What?!?!'
    import sys
    sys.exit(0)


# subset = ['metFinal']
# settings = { k:settings[k] for k in subset }
# settings = { k:settings[k] for k in settings if k.startswith('met') }


# Create output direcotries
for d in [outhistos, outpdfs, outtxt]:
    if os.path.exists(d): continue
    print 'Creating',d
    os.makedirs(d)

# import sys
# sys.exit(0)

histoSig = {}


if opt.data:
    print 'Generating TDR-ratio plots'
    from plots.services import HistoRatio as Histo
    from plots.services import StackRatio as Stack
    from plots.services import LegendRatio as Legend
else:
    print 'Generating TDR plots'
    from plots.services import Histo, Stack, Legend


# for var, title in vars.iteritems():
for var,(title,scale,rebin, usrrng) in settings.iteritems():

    if not title:
        title = ''

    print "Variable: ",var
    print "Title: ", title
    print "Scale: ", scale
    print "Rebin: ", rebin

    #THStack for plotting
    stack_bkg = Stack(var,title)
    stack_sig = Stack(var,title)
    stack_data = Stack(var,title)
    stack_sh  = Stack(var,title)

    leg = Legend()
    leg_sh = Legend()
    # 
    leg_sign = ROOT.TLegend(0.7,0.7,0.9,0.9)

    h1 = None

    for s in samples.itervalues():
        nEntries = 0
        #print "+ sample" , s.label


        if(hasattr(s, "components")):
            #samples with components
            histos = []
            notFound = []
            for c in s.components:
                #print "comp ", c.label

                if c.label.startswith("MET") and opt.channel == "semileptonic": continue

                elif ((c.label.startswith("SingleMu") or c.label.startswith("SingleEl")) and (("CR3" not in var) and ("CR6" not in var) and ("CR7" not in var) and ("CRtt" not in var) and ("SR_1lep" not in var))  and (opt.channel == "fullhadronic")):
                    print "Fullhadronic channel, skipping dataset: " + c.label +" for variable " + var
                    continue

                elif ( c.label.startswith("MET") and (("CR3" in var) or ("CR6" in var) or ("CR7" in var) or ("CRtt" in var) or ("SR_1lep" in var)) and opt.channel == "fullhadronic"):
                    print "Fullhadronic channel, skipping dataset: " + c.label +" for variable " + var
                    continue

                #print "c.label " , c.label

                if(opt.sys=="noSys"): 
                    filename = opt.resdir+'/res_SR/'+c.label+"_" + opt.channel +".root"
                elif (
                        ( opt.sys=="jesUp" or opt.sys=="jesDown" or opt.sys=="jerUp" or opt.sys=="jerDown" ) and 
                        ( c.label.startswith("SingleMu") or c.label.startswith("SingleEl") or c.label.startswith("MET"))
                    ):
                    filename = opt.resdir+'/res_SR/'+c.label+"_" + opt.channel +".root"
                else: 
                    filename = opt.resdir+'/res_SR/'+c.label+"_" + opt.channel + "_"+ opt.sys +".root"

                #cutFlow histo saved only in nominal root file
                filename_nEvt = opt.resdir+'/res_SR/'+c.label+"_" + opt.channel +".root"
                infile_nEvt = ROOT.TFile.Open(filename_nEvt)
                nEntries = infile_nEvt.Get("cutFlow").GetBinContent(0)

                infile = ROOT.TFile.Open(filename)
                # htmp = infile.Get(var)
                hin = infile.Get(var)
                if not isinstance(hin, ROOT.TH1):
                    notFound.append( (c.label,filename) )
                    # raise RuntimeError('Failed to load histogram %s from %s' % (var,filename))
                    continue

                htmp = hin.Clone()
                # print '--> Sums before scaling: w2:',htmp.GetSumw2().GetSize(),' I(h):',htmp.Integral(), ' E(h):', sumerrors(htmp)

                # Ensure the error structure is activated before doing anything else
                htmp.Sumw2()


                # Applu lumi scale if not data
                if( nEntries != 0 and not(c.label.startswith("SingleMu") or c.label.startswith("SingleEl") or c.label.startswith("JetHT") or c.label.startswith("MET"))):
                    htmp.Scale((1./nEntries) * c.sigma * 1000.* float(opt.lumi) )

                # print '--> Sums after scaling: I(h):',htmp.Integral(), ' E(h):', sumerrors(htmp)

                # If a cutflow print a nice summary!
                if(var == "cutFlow"): writeSummary(c.label, opt.channel, opt.sys, htmp, opt.lumi)

                # To remove
                if(htmp.Integral()!=0):
                    htmp.SetMarkerStyle(20)
                    htmp.SetMarkerSize(1.2)
                    htmp.SetLineWidth(2)
                    htmp.SetLineColor(ROOT.kBlack)
                histos.append(htmp)

            if notFound:
                raise RuntimeError('Failed to retrieve %s' % notFound)

            # Sum components histos
            h1 = histos[0]
            [h1.Add(hi) for hi in histos[1:]]

            # h = Histo1(h1)
            # print 'h1 has sumw2', h1.GetSumw2().GetSize()
            
            h = Histo.fromTH1(h1)

        else:
            #sample does not have components
            #print "sample label ", s.label
            if opt.sys=="noSys":
                filename = opt.resdir+'/res_SR/'+s.label+"_" + opt.channel +".root"
            else:
                filename = opt.resdir+'/res_SR/'+s.label+"_" + opt.channel +"_"+ opt.sys +".root"

            #cutFlow histo saved only in nominal root file
            filename_nEvt = opt.resdir+'/res_SR/'+s.label+"_" + opt.channel +".root"
            infile_nEvt = ROOT.TFile.Open(filename_nEvt)
            nEntries =  infile_nEvt.Get("cutFlow").GetBinContent(0)

            infile = ROOT.TFile.Open(filename)
            hin = infile.Get(var)
            if not isinstance(hin, ROOT.TH1):
                raise RuntimeError('Failed to load histogram %s from %s' % (filename, var))

            htmp = hin.Clone()


            # Ensure the error structure is activated before doing anything else
            htmp.Sumw2()

            htmp.SetMarkerStyle(20)
            htmp.SetMarkerSize(1.2)

            if(s.label.startswith("DM")):
               htmp.SetLineWidth(3)

            # Create the Histo wrapper
            h = Histo.fromTH1(htmp)

            scaleFactor = 1.

            if (nEntries != 0  and not(s.label.startswith("SingleMu") or s.label.startswith("SingleEl") or s.label.startswith("JetHT") or s.label.startswith("MET"))):
                scaleFactor = 1./nEntries * s.skimEff * s.sigma * 1000.* float(opt.lumi)
                htmp.Scale(scaleFactor)

            h.Scale(scaleFactor)

        ### Apply rebin
        # if(var in rebin.keys()):
        if rebin:
            h.Rebin(rebin)

        ### Re-range
        if usrrng is not None:

            # Don't chop by default
            if len(usrrng) == 2:
                # Don't chop if the third parameter is omitted
                pass
            elif len(usrrng) == 3 and usrrng[2]:
                uFirst, uLast = usrrng[0], usrrng[1]

                # Pull out the histogram 
                hx = h.GetHisto()

                # Extract the old histogra,
                oldint = hx.Integral(0, hx.GetNbinsX()+1)
                xax = hx.GetXaxis()
                # print '+ Usrrange:',usrrng
                # print '+ First, Last:',xax.GetFirst(), xax.GetLast()

                # Use SetRangeUser to calculate first and last bib
                xax.SetRangeUser(uFirst, uLast)
                # print '+ First, Last:',xax.GetFirst(), xax.GetLast()

                # Here they are
                first, last = xax.GetFirst(), xax.GetLast()

                # new overflow
                ofw = 0
                # and new overflow weight
                w2ofw = 0

                # loop from last to overflow to calculate the new values
                for ib in xrange(last, hx.GetNbinsX()+2):
                    ofw += hx.GetAt(ib)
                    w2ofw += hx.GetSumw2().GetAt(ib)

                    # Zero the bin
                    hx.SetAt(0, ib)
                    hx.GetSumw2().SetAt(0, ib)

                # print '+ ofw, w2ofw: ',ofw, w2ofw

                # Re-insert overflows into the last bin
                hx.SetAt(ofw, last)
                hx.GetSumw2().SetAt(w2ofw, last)

                # Now check we didn't screw anything up
                newint = hx.Integral(0, hx.GetNbinsX()+1)

                # print '+ old->new', oldint, newint

        ### Set style and create stack plots
        h.SetStyle(s.color, s.style, s.fill)

        # Writing outfile: only 8 bin for met  Plot, overflow bin included
        if opt.sys=="noSys":
            outfilename = "%s/%s_%s.root" % (outhistos,s.label,opt.channel)
        else:
            outfilename = "%s/%s_%s_%s.root" % (outhistos,s.label,opt.channel, opt.sys)

        outfile = ROOT.TFile(outfilename, "UPDATE")

        outfile.cd()
        if var in store:
            h.GetHisto().SetName(var)
            h.GetHisto().Write()
        outfile.Close()

        #adding histo in bkg, sign and data categories
        if (s.label=="DMtt_sc_Mchi1Mphi10" or s.label=="DMtt_sc_Mchi1Mphi500"):
            stack_sig.Add(h)
            label = s.leglabel if scale == 1. else '%s (x%d)' % (s.leglabel,scale)
            leg_sign.AddEntry(h.GetHisto(), label , "l")
            h.Scale(scale) # make stack plots of bkground and signal x10 (for few variables)
        elif (s.label.startswith("Data")):
            stack_data.Add(h)
            if opt.data:
                leg.AddEntry(h, s.leglabel, "lp")
        elif ( not(s.label.startswith("DM") or s.label.startswith("Data"))) :
            stack_bkg.Add(h)
            leg.AddEntry(h, s.leglabel, "f")

        ### Make a summary
        if(var == "cutFlow"): writeSummary(s.label, opt.channel,opt.sys, h, opt.lumi)

        # make stack plot of histos normalized to unity
        h_norm = copy.deepcopy(h)
        h_norm.Normalize()
        leg_sh.AddEntry(h_norm, s.leglabel, "l")
        stack_sh.Add(h_norm)


    ### End of samples loop
    
    H=600
    W=700

    H_ref = 600
    W_ref = 700

    T = 0.08*H_ref
    B = 0.12*H_ref
    L = 0.12*W_ref
    R = 0.08*W_ref

    tdrstyle.setTDRStyle();

    # Update the legend's stype
    leg_sign.SetNColumns(2)
    leg_sign.SetFillColor(0)
    leg_sign.SetFillStyle(0)
    leg_sign.SetTextFont(42)
    leg_sign.SetBorderSize(0)

    ### Adjust ranges and create/save plots: some settings are defined in tdrstyle.py and service_ratio.py
    c1 = ROOT.TCanvas("c1","c1",50,50,W,H)
    c1.SetFillColor(0);
    c1.SetBorderMode(0);
    c1.SetFrameFillStyle(0);
    c1.SetFrameBorderMode(0);
    c1.SetLeftMargin( L/W );
    c1.SetRightMargin( R/W );
    c1.SetTopMargin( 1 );
    c1.SetBottomMargin(0);
    c1.SetTickx(0);
    c1.SetTicky(0);
    c1.cd()

    if opt.data:
        pad1= ROOT.TPad("pad1", "pad1", 0, 0.30 , 1, 1)
        pad1.SetTopMargin(0.1)
        pad1.SetBottomMargin(0)
        pad1.SetLeftMargin(0.12)
        pad1.SetRightMargin(0.05)


        # print 'x1ndc',leg_sign.GetX1NDC(), leg_sign.GetX1()
        leg_sign.SetTextSize(0.04)

        leg_sign.SetX1(.40)
        leg_sign.SetY1(.59)
        leg_sign.SetX2(.91)
        leg_sign.SetY2(.74)

    else:
        pad1= ROOT.TPad("pad1", "pad1", 0, 0 , 1, 1)

        leg_sign.SetTextSize(0.032)

        leg_sign.SetX1(.40)
        leg_sign.SetY1(.69)
        leg_sign.SetX2(.94)
        leg_sign.SetY2(.84)

    pad1.SetBorderMode(0);
    #if(not(var.startswith("h_topHadMass") or var.startswith("h_topSLMass"))):
    #pad1.SetLogy() # log scale, can be asked to not be applied for some variables
    pad1.SetTickx(0);
    pad1.SetTicky(0);
    pad1.Draw()
    pad1.cd()

    # Normalizing to data
    if(opt.normData>0):
        h_data_ = stack_data._hs.GetStack().Last()
        nData = h_data_.Integral()
        h_bkg_ = stack_bkg._hs.GetStack().Last()
        nBkg = h_bkg_.Integral()
        sf_norm = nData/nBkg

        for i in xrange(stack_bkg._hs.GetNhists()):
            h_i=stack_bkg._hs.GetHists().At(i)
            h_i.Scale(sf_norm)

        stack_bkg._hs.Modified()
        #stack_data.GetStack().At(stack_data.GetNhists()-1)
        print " NORMALIZING TO DATA"
        h_bkg.Scale(sf_norm)


    #Set Range user for a selection of plot: use same format for other variables if needed
    if usrrng is not None:
        stack_bkg.SetRangeUser( usrrng[0], usrrng[1] )
        stack_sig.SetRangeUser( usrrng[0], usrrng[1] )
        stack_data.SetRangeUser( usrrng[0], usrrng[1] )
        stack_sh.SetRangeUser( usrrng[0], usrrng[1] )

    h_data = stack_data.GetLast()
    h_bkg = stack_bkg.GetLast()

    # Print some yields
    # stack_bkg.PrintYields()

    maximum = max([stack_bkg.GetMaximum(),stack_data.GetMaximum(), h_data.GetMaximum(), stack_sig.GetMaximum()])
    minimum = min([stack_bkg.GetMinimum(),stack_data.GetMinimum(), h_data.GetMinimum(), stack_sig.GetMinimum()])

    #Maximum for log scale
    stack_bkg.SetMaximum(maximum*( 5000 if opt.data else 10000) )
    #Maximum for non log scale
    stack_bkg.SetMaximum(maximum*1.8)
    if(minimum > 0):
        stack_bkg.SetMinimum(0.01*minimum)
    else: 
        stack_bkg.SetMinimum(0.01)

    #Drawing THStacks
    stack_bkg.DrawStack(lumi = opt.lumi)
    stack_sig.DrawStack(lumi = opt.lumi, opt = "nostack, SAME")

    ROOT.gStyle.SetHatchesSpacing(2)
    ROOT.gStyle.SetHatchesLineWidth(2)

    h_err = h_bkg.Clone("h_err")
    h_err.SetLineWidth(100)
    h_err.SetFillStyle(3154)
    h_err.SetMarkerSize(0)
    h_err.SetFillColor(ROOT.kGray+2)
    h_err.Draw("e2same")


    if opt.data:
        # h_data.Sumw2()
        h_data.Draw("eSAMEpx0")
        stack_bkg.GetHistogram().GetXaxis().SetLabelOffset(1.8)


    leg.Draw("SAME")
    leg_sign.Draw("same")

    #***
    if opt.data:
        c1.cd()
        pad2 = ROOT.TPad("pad2", "pad2", 0, 0.03, 1, 0.29)
        pad2.SetTopMargin(0.05)#4
        pad2.SetBottomMargin(0.36)
        pad2.SetLeftMargin(0.12)
        pad2.SetRightMargin(0.05)

        c1.cd()
        pad2.Draw()
        pad2.cd()

        ratio  = h_data.Clone("ratio")
        ratio.SetLineColor(ROOT.kBlack)
        ratio.SetMinimum(0.1)
        ratio.SetMaximum(3.0)
        ratio.Sumw2()
        ratio.SetStats(0)
        if usrrng is not None:
            ratio.GetXaxis().SetRangeUser( usrrng[0], usrrng[1] )
            
        
        denom = h_bkg.Clone("denom")
        denom.Sumw2()

        if(denom.Integral() !=0):
            ratio.Divide(denom)
            ratio.SetMarkerStyle(20)
            ratio.SetMarkerSize(1.2)
            ratio.Draw("epx0")
            ratio.SetTitle("")

        h_bkg_err = h_bkg.Clone("h_err")
        h_bkg_err.Reset()
        h_bkg_err.Sumw2()

        for i in range(1,h_bkg.GetNbinsX()):
            h_bkg_err.SetBinContent(i,1)
            if(h_bkg.GetBinContent(i)):
                h_bkg_err.SetBinError(i, (h_bkg.GetBinError(i)/h_bkg.GetBinContent(i)))
            else:
                h_bkg_err.SetBinError(i, 0)
        #h_bkg_err = h_bkg.Clone("h_bkg_err")
        #h_bkg_err.Sumw2()
        #h_bkg_err.Divide(h_bkg);
        #h_bkg_err.Draw("E2same");
        #h_bkg_err.SetMaximum(2.);
        #h_bkg_err.SetMinimum(0);
        h_bkg_err.SetLineWidth(100)
        h_bkg_err.SetFillStyle(3154)
        h_bkg_err.SetMarkerSize(0)
        h_bkg_err.SetFillColor(ROOT.kGray+2)

        #for i in range(1,h_bkg.GetNbinsX()):
        #    if(var=="met_preS"):
        #        print "data ", h_data.GetBinContent(i), " error " , h_data.GetBinError(i)
        #        print "bkg ", h_bkg.GetBinContent(i), " error " , h_bkg.GetBinError(i)
        #        print "bkg_err ", h_bkg_err.GetBinContent(i), " error " , h_bkg_err.GetBinError(i)


        f1 = ROOT.TF1("myfunc","[0]",-100000,10000);
        f1.SetLineColor(ROOT.kBlack)
        f1.SetLineStyle(ROOT.kDashed)
        f1.SetParameter(0,1);
        f1.Draw("same")

        #some settings for histo are defined in tdrstyle.py and service_ratio.py
        ratio.Draw("epx0same")
        ratio.Sumw2()
        ratio.GetYaxis().SetTitle("#frac{Data}{MC}")
        ratio.GetYaxis().SetNdivisions(503)

        ratio.GetXaxis().SetLabelFont(42);
        ratio.GetYaxis().SetLabelFont(42);
        ratio.GetXaxis().SetTitleFont(42);
        ratio.GetYaxis().SetTitleFont(42);

        ratio.GetXaxis().SetTitleOffset(1.2);
        ratio.GetYaxis().SetTitleOffset(0.29);#0.31

        ratio.GetXaxis().SetLabelSize(0.12);
        ratio.GetYaxis().SetLabelSize(0.12);
        ratio.GetXaxis().SetTitleSize(0.14);
        ratio.GetYaxis().SetTitleSize(0.16);

        ratio.GetYaxis().SetRangeUser(0.5,1.5);#-2.3,3.7

        ratio.GetXaxis().SetTitle(title)
        ratio.GetXaxis().SetLabelOffset(0.06)
        ratio.GetYaxis().SetLabelOffset(0.015)

    CMS_lumi.writeExtraText = 1
    CMS_lumi.extraText = "Preliminary"
    lumi = opt.lumi
    if(lumi<1.):
        lumi = lumi*1000
        unit = " pb^{-1}"
    else: unit = " fb^{-1}"

    CMS_lumi.lumi_sqrtS = str(lumi)+ unit +" (13 TeV)" # used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

    iPeriod = 0
    iPos = 11

    # writing the lumi information and the CMS "logo"
    # Ratio Check HERE
    CMS_lumi.CMS_lumi(pad1, iPeriod, iPos)

    c1.cd()
    c1.Update();
    c1.RedrawAxis();
    #c1.GetFrame().Draw();

    pdfbasename = var+'_'+opt.channel

    syst_label = '' if opt.sys == 'noSys' else ('_'+opt.sys)
    data_label = '' if opt.data else '_nodata'

    pdfname = outpdfs+'/'+pdfbasename+syst_label+data_label+'.pdf'
    rootname = outpdfs+'/'+pdfbasename+syst_label+data_label+'.root'
    pngname = outpdfs+'/'+pdfbasename+syst_label+data_label+'.png'

    c1.Print(pdfname)
    c1.Print(rootname)
    c1.Print(pngname)

