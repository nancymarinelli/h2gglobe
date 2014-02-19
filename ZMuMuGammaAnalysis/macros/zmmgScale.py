#!/bin/env python

from optparse import OptionParser, make_option
import sys, os
import array
from math import sqrt, fabs

objs = []

def trMean(hist,alpha):

    beta = 1.-alpha
    probs = array.array('d',[beta*0.5,1.-beta*0.5])
    quantiles = array.array('d',[0.,0.])

    hist.GetQuantiles(2,quantiles,probs)

    print quantiles, probs
    
    limits = [ hist.GetXaxis().FindBin(q) for q in quantiles ]
    hist.GetXaxis().SetRange(*limits)
    
    mean = hist.GetMean()
    err  = hist.GetRMS() / ( ( 1. - beta )*sqrt( hist.Integral() ) )
        
    hist.GetXaxis().SetRange()
    return mean,err

def winMean(hist,rng):

    limits = [ hist.GetXaxis().FindBin(q) for q in rng ]
    hist.GetXaxis().SetRange(*limits)
    
    mean = hist.GetMean()
    err  = hist.GetRMS() / ( sqrt( hist.Integral(*limits) ) )
        
    hist.GetXaxis().SetRange()
    return mean,err

def getMean(h,alpha):
    if type(alpha) == float:
        print alpha
        x,xerr = trMean(h,alpha)
    elif type(alpha) == tuple:
        x,xerr = winMean(h,alpha)
    else:
        x,xerr = alpha(h)

    return x,xerr

def recFit(h):
    f = ROOT.TF1("f","gaus",80.,100.)
    oldmean = 90.
    oldsigma = 10.
    
    h.Fit(f,"QRN")
    errmean  = f.GetParError(1)
    mean  = f.GetParameter(1)
    sigma = f.GetParameter(2)*1
    iter = 0
    while iter < 10 or fabs(1. - mean/oldmean) > 0.005:
        g = ROOT.TF1("g","gaus",mean-sigma,mean+sigma)
        h.Fit(g,"QRN")
        oldmean = float(mean)
        oldsigma = float(sigma)
        errmean  = g.GetParError(1)
        mean  = g.GetParameter(1)
        sigma = g.GetParameter(2)*1
        iter+=1
        ## print iter, mean, oldmean, sigma

    g = ROOT.TF1("final_fit_%s" % h.GetName(),"gaus",mean-sigma,mean+sigma)
    g.SetLineColor(h.GetLineColor())
    h.Fit(g,"QRO+")
    ## h.GetListOfFunctions().Add(g.Clone())
    objs.append(g)
    
    return mean,errmean
        

def toStr(alpha):
    if type(alpha) == float:
        return "%d" % (100.*alpha)
    elif type(alpha) == tuple:
        return "%1.2g_%1.2g" % alpha
    else:
        return alpha.__name__

def buildCalib(hmcs,name,title,alpha=0.683):

    
    gr = ROOT.TGraphErrors()
    gr.SetName(name)
    gr.SetTitle(title)
    
    for y in sorted(hmcs.keys()):
        x,xerr = getMean(hmcs[y],alpha)
        ip = gr.GetN()
        gr.SetPoint(ip,x,y)
        gr.SetPointError(ip,xerr,0.)
        
    return gr
    
def main(options,args):
    objs = []

    infile = args[0]
    fin = ROOT.TFile.Open(infile)

    results = []
    ## alphas = [1,0.95,0.8,0.683,0.5,0.4,0.3,0.2]
    ## alphas = [recFit,0.683,(85.,100.),(80.,100.),(86.,96.),(88.,94.)]
    alphas = [recFit]
    for icat in range(options.ncat):
        hdata = fin.Get("th1f_data_mass_cat%d" % icat)

        hmcs = {}
        hmcs[0] = fin.Get("th1f_sig_dymm_mass_m90_cat%d" % icat)
        
        for isig in range(1,options.nsigma+1):
            hmcs[options.step*float(isig)] = fin.Get("th1f_sig_dymm_mass_m90_cat%d_E_scaleUp0%d_sigma" % (icat,isig) )
            hmcs[-options.step*float(isig)] = fin.Get("th1f_sig_dymm_mass_m90_cat%d_E_scaleDown0%d_sigma" % (icat,isig) )

        hdata.SetLineColor(ROOT.kBlack)
        hdata.SetMarkerColor(ROOT.kBlack)
        for s,h in hmcs.iteritems():
            h.SetLineColor(ROOT.kRed)
        hmcs[0].SetLineColor(ROOT.kBlue)
        
        canv = ROOT.TCanvas("cat_%d" % icat, "cat_%d" % icat )
        canv.cd()
        hmcs[0].GetXaxis().SetRangeUser(80.,100.)
        res = {}
        for alpha in alphas:
            name = "mVsDeltaE_%s" % toStr(alpha)
                
            data = getMean(hdata,alpha)
            res[alpha] = ( buildCalib(hmcs,name,";m_{ll#gamma} (GeV/c^{2});#Delta E_{#gamma} / E_{#gamma} (#times 10^{-2})", alpha ),
                           data )

        hmcs[0].Draw("hist")
        hdata.SetLineColor(ROOT.kBlack)
        hdata.SetMarkerColor(ROOT.kBlack)
        hdata.Draw("e sames")
        for s,h in hmcs.iteritems():
            if s == 0.:
                continue
            h.SetLineColor(ROOT.kRed)
            h.Draw("histsames")
        objs.append( (canv, hdata, hmcs ))
        
        for fmt in "C","png","pdf":
            canv.SaveAs("%s.%s" % ( canv.GetName(), fmt) )
        

        results.append(res)
        
    objs.append(results)

    meas = {}
    for alpha in alphas:
        scales = []
        name = "results_%s" % toStr(alpha)
        canv = ROOT.TCanvas(name,name)
        canv.Divide(options.ncat/2,2)
        objs.append(canv)
        for icat,res in enumerate(results):
            canv.cd(icat+1)
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()
            calib, data = res[alpha]
            calib.Fit("pol1","W+")
            func = calib.GetListOfFunctions().At(0)
            scale = func.Eval(data[0])
            scaleP = func.Eval(data[0]+data[1])
            scaleM = func.Eval(data[0]-data[1])
            scales.append( (scale,scaleM,scaleP) )
            ip = calib.GetN()
            calib.SetPoint( ip, data[0], scale )
            calib.SetPointError( ip, data[1], sqrt(2.)*0.5*(scaleP-scaleM) )
            calib.Draw("ap")
        for fmt in "C","png","pdf":
            canv.SaveAs("%s.%s" % (canv.GetName(),fmt) )
        print "Alpha: ", alpha
        print "Scales: ", scales
        
    return objs

if __name__ == "__main__":
    parser = OptionParser(option_list=[
        make_option("-s", "--step",
                    action="store", type="float", dest="step",
                    default=0.2,
                    ),
        make_option("-N", "--nsigma",
                    action="store", type="int", dest="nsigma",
                    default=2,
                    ),
        make_option("-n", "--ncat",
                    action="store", type="int", dest="ncat",
                    default=8,
                    ),
        make_option("-l", "--label",
                    action="append", dest="labels",
                    default=[],
                    ),
        make_option("-o", "--outdir",
                    action="store", dest="outdir", type="string",
                    default="zmmgScale",
                    ),
        ])
    
    (options, args) = parser.parse_args()

    print options
    
    sys.argv.append("-b")
    import ROOT
    
    objs = main( options, args )
