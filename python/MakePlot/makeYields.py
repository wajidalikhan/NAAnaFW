#!/bin/env python

from extern.matrix2latex import matrix2latex

import os.path
import collections
import re

###########################################################
# Regions, Samples and Labels
###########################################################

regions = collections.OrderedDict( [
    ('SR met>160', ('metFinal', None) ),
    ('SR met>320', ('metFinal', (320, None)) ),
    # Fullhadronic specific
    ('FH-CR tt(1l)' , ('metFinal_SR_1lep', None) ),
    #('FH-CR QCD'    , ('metFinal_outtop', None) ),
    ('FH-CR 0b'     , ('metFinal_CR5', None) ),
    ('FH-CR 0b 1l'  , ('metFinal_CR6nw', None) ),
    ('FH-CR 0b 2l'  , ('metFinal_CR7nw', None) ),
    # Semileptonic specific
    ('SL-CR tt(2l)'  , ('metFinal_2lep', None) ),
    ('SL-CR 0b'  , ('metFinal_met_0btag', None) ),
    ]
)

data = [
    'Data'
    ]
signal = [
    'DMtt_ps_Mchi1Mphi100'#,
    #'DMtt_sc_Mchi50Mphi300'
    ]
backgrounds = [
    'TT_1lep',
    'TT_2lep',
    'WJets',
    'DY',
    'SingleTop',
    'QCD',
    'ZToNuNu',
    'otherBkg',
    'VV'
]

labels = {
    # Region labels
    'SR met>160': r'SR ${E\!\!\!/}_{\!\mathrm{T}}>$160',
    'SR met>320': r'SR ${E\!\!\!/}_{\!\mathrm{T}}>$320',
    # Fullhadronic
    'FH-CR tt(1l)' : 'CR tt(1l)',
    #'FH-CR QCD'    : 'CR QCD',
    'FH-CR 0b'     : 'CR 0b',
    'FH-CR 0b 1l'  : 'CR 0b 1l',
    'FH-CR 0b 2l'  : 'CR 0b 2l',
    # Semileptonic
    'SL-CR tt(2l)': 'CR tt(2l)',
    'SL-CR 0b' : 'CR 0b',
    # Sample labels
    'DMtt_sc_Mchi1Mphi10': r'sc: $M_{\chi} 1\;\mathrm{GeV}, M_{\phi} 10\;\mathrm{GeV}$',
    'DMtt_sc_Mchi50Mphi300': r'sc: $M_{\chi} 1\;\mathrm{GeV}, M_{\phi} 500\;\mathrm{GeV}$',
    'DMtt_ps_Mchi1Mphi100': r'ps: $M_{\chi} 1\;\mathrm{GeV}, M_{\phi} 100\;\mathrm{GeV}$',
    'ZToNuNu': r'Z$\rightarrow\nu\nu$',
}


###########################################################
# Yield helper class
###########################################################


class Yields(object):
    """docstring for Yields"""

    _backgrounds_label = 'Backgrounds'

    def __init__(self,regions, sigs, bkgs, data, labels={}):
        super(Yields, self).__init__()
        self.signals = sigs
        self.backgrounds = bkgs
        self.data = data
        self.regions = regions
        self.labels = labels
        self.values = { r:{} for r in regions }

    @property
    def samples(self):
        return self.signals+self.backgrounds+self.data

    def backgroundSums(self):

        sums = { r:{} for r in regions }

        for r in self.regions:
            nbkg = 0
            e2bkg = 0

            for s in self.backgrounds:
                nev,err = self.values[r][s]

                nbkg += nev
                e2bkg += err**2

            sums[r][self._backgrounds_label] = (nbkg, e2bkg**0.5)

        return sums

    def ascii(self):
        for r in self.regions:
            print 'Region:',r
            for s in self.samples:
                print '   ',s,'(%.1f+/-%.1f)' % self.values[r][s]

    def latex(self, bkgsum=False):
        matrix = []
        hdr = ['']+[ self.labels[r] if r in self.labels else r for r in self.regions]

        import copy
        # Copy values, in case we need to mangle them
        values = copy.deepcopy(self.values)

        if bkgsum:
            for r in self.regions:
                sums = self.backgroundSums()
                values[r].update(sums[r])

            samples = self.backgrounds+[self._backgrounds_label]+self.signals+self.data
        else:
            samples = self.backgrounds+self.signals+self.data

        for s in samples:

            # row starts with sample name
            row = [self.labels[s] if s in self.labels else s]

            # followed by events and errors in order
            row += [ '$%.2f \pm %.2f$' % values[r][s] for r in self.regions ]

            matrix.append(row)

        # Create the table
        table =  matrix2latex(matrix,headerRow=hdr)

        rows = table.splitlines()
        # create a regex that matches all rules
        reRules = re.compile(r'(\\toprule|\\midrule|\\bottomrule)')

        # replace lines with rules (hoping doens't break anything)
        rows = [ reRules.sub(r'\hline',r) for r in rows ]


        if bkgsum:
            for i,r in enumerate(rows):
                # search for the table block
                try:
                    j = r.index('&')
                except:
                    # not found
                    continue

                if not self._backgrounds_label in r[:j]: continue

                # count how many tabs to re-insert
                m = re.search('[^\t]', r)
                # m.start is the first non-tab
                # '\t'*m.start() is a string with n tabs
                tabs = '\t'*m.start() if m else ''
                # insert additional horizontal lines : order matters!
                rows.insert(i+1,tabs+'\hline')
                rows.insert(i,tabs+'\hline')
                break

        table = '\n'.join(rows) 
                    
        return table




###########################################################
# Class to fill a Yeild obj based on a region definiton set
###########################################################
class ChannelYieldFiller(object):
    """docstring for ChannelYieldFiller"""
    def __init__(self, regions, histopath, channel):
        super(ChannelYieldFiller, self).__init__()
        self.regions = regions
        self.histopath = histopath
        self.channel = channel

    def fill( self, yields ):
         # Import ROOT only here, to avoid it interfering with parsing
        from ROOT import TFile, TH1, Double

        for s in yields.samples:
            for r in yields.regions:
                filename = s+'_'+self.channel+'.root'
                filepath = os.path.join(self.histopath,filename)
                hname,rng = self.regions[r]

                rf = TFile(filepath)
                if rf.IsZombie():
                    raise RuntimeError('Failed to open %s',filepath)
                htmp = rf.Get(hname)
                # print htmp
                if not isinstance(htmp,TH1):
                    raise RuntimeError('Histogram %s not found in %s' % (hname,filepath))

                h = htmp.Clone()
                # Adapt the range, if required
                xax = h.GetXaxis()
                first, last = xax.GetFirst(), xax.GetLast()
                xmin, xmax = xax.GetBinLowEdge(first), xax.GetBinUpEdge(last)

                if rng is not None:
                    umin,umax = rng
                    umin = umin if umin is not None else xmin
                    umax = umax if umax is not None else xmax

                    # print umin,umax
                    # Apply new range (use ROOT to do the difficult bit :P)
                    xax.SetRangeUser(umin,umax)
                
                    # Update limits
                    first, last = xax.GetFirst(), xax.GetLast()
                    xmin, xmax = xax.GetBinLowEdge(first), xax.GetBinUpEdge(last)

                err = Double(0.)
                nev = h.IntegralAndError(first,last,err)
                rf.Close()

                print filename, hname,rng,
                print ' | ','%.2f(%d)'%(xmin,first),'-', last,'%.2f(%d)'%(xmax,last), ' -- ', nev,err

                del h

                yields.values[r][s] = (nev,err)       


###########################################################
# Main body
###########################################################
if __name__ == '__main__':

    import optparse

    usage = 'usage: %prog -l lumi'
    parser = optparse.OptionParser(usage)
    parser.add_option('-c', '--channel', dest='channel', type='string', default = 'semileptonic', help='Channel to analyze: semileptonic or fullhadronic')
    (opt, args) = parser.parse_args()

    basepath = 'output'
    histodir = 'histos_lin'

    if opt.channel == 'fullhadronic':
        chanLabel = 'fh'

        # And the yields to fill
        yields = {
            'signalRegions': Yields( ['SR met>160', 'SR met>320'], signal, backgrounds, data, labels=labels ),
            'controlRegions': Yields( ['FH-CR tt(1l)','FH-CR 0b', 'FH-CR 0b 1l','FH-CR 0b 2l'], signal, backgrounds, data, labels=labels ),
        }

    elif opt.channel == 'semileptonic':
        chanLabel = 'sl'

        # And the yields to fill
        yields = {
            'signalRegions': Yields( ['SR met>160', 'SR met>320'], signal, backgrounds, data, labels=labels ),
            'controlRegions': Yields( ['SL-CR tt(2l)', 'SL-CR 0b' ], signal, backgrounds, data, labels=labels ),
        }

    else:
        parser.error('What channel?')

    # Where are the histograms?
    histopath = os.path.join(basepath,chanLabel,histodir)

    # Create a filler object
    filler = ChannelYieldFiller(regions, histopath, opt.channel)

    # destination directory
    outtables = 'output/%s/tables' % (chanLabel,)
    if not os.path.exists(outtables):
        os.makedirs(outtables)
    
    for name,y in yields.iteritems():
        print 'Calculating yeilds:',name
        filler.fill(y)

        # print some details
        y.ascii()
        
        table = y.latex(bkgsum=True)

        # more debug stuff
        print table


        with open(outtables+'/'+name+'.tex','w') as dottex:
            dottex.write(table)

