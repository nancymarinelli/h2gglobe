#!/bin/env python

from optparse import OptionParser, make_option
import sys, os
import array
from math import sqrt, fabs

objs = []

def setStat(h,prev=None,vert=-1,horiz=0.):
    ROOT.gPad.Update()
    st = h.FindObject('stats')
    print st, h, prev
    st.SetLineColor(h.GetLineColor())
    st.SetTextColor(h.GetLineColor())

    if prev:
        shiftx = (prev.GetX2NDC() - st.GetX1NDC())*horiz
        shifty = (prev.GetY2NDC() - st.GetY1NDC())*vert

        st.SetX1NDC(st.GetX1NDC()+shiftx)
        st.SetX2NDC(st.GetX2NDC()+shiftx)

        st.SetY1NDC(st.GetY1NDC()+shifty)
        st.SetY2NDC(st.GetY2NDC()+shifty)

    ROOT.gPad.Update()
    return st

def trMean(hist,alpha):

    beta = 1.-alpha
    probs = array.array('d',[beta*0.5,1.-beta*0.5])
    quantiles = array.array('d',[0.,0.])

    hist.GetQuantiles(2,quantiles,probs)

    print quantiles, probs
    
    limits = [ hist.GetXaxis().FindBin(q) for q in quantiles ]
    hist.GetXaxis().SetRange(*limits)
    
    mean = hist.GetMean()
    err  = hist.GetMeanError() ## hist.GetRMS() / ( ( 1. - beta )*sqrt( hist.Integral() ) )
        
    hist.GetXaxis().SetRange()
    return mean,err

def winMean(hist,rng):

    limits = [ hist.GetXaxis().FindBin(q) for q in rng ]
    hist.GetXaxis().SetRange(*limits)
    
    mean = hist.GetMean()
    err  = hist.GetMeanError()
        
    hist.GetXaxis().SetRange()
    return mean,err

def getMean(h,method):
    if type(method) == float:
        print method
        x,xerr = trMean(h,method)
    elif type(method) == tuple:
        x,xerr = winMean(h,method)
    else:
        x,xerr = method(h)

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
    while iter < 5 or fabs(1. - mean/oldmean) > 0.005:
        g = ROOT.TF1("g","gaus",mean-sigma,mean+sigma)
        h.Fit(g,"QRN")
        oldmean = float(mean)
        oldsigma = float(sigma)
        errmean  = g.GetParError(1)
        mean  = g.GetParameter(1)
        sigma = g.GetParameter(2)*1
        iter+=1
        ## print iter, mean, oldmean, sigma

    fcanv = ROOT.TCanvas("fit_%s" % h.GetName(),"fit_%s" % h.GetName())
    fcanv.cd()
    g = ROOT.TF1("final_fit_%s" % h.GetName(),"gaus",mean-sigma,mean+sigma)
    g.SetLineColor(h.GetLineColor())
    h.Fit(g,"QRO")
    ## h.GetListOfFunctions().Add(g.Clone())
    h.Draw()
    g.Draw("same")
    objs.append(g)
    objs.append(fcanv)
    objs.append(h)
    
    for fmt in ["png"]:
        fcanv.SaveAs( "%s.%s" % ( fcanv.GetName(), fmt ) )
                     
    return mean,errmean
        

def toStr(method):
    if type(method) == float:
        return "%d" % (100.*method)
    elif type(method) == tuple:
        return "%1.2g_%1.2g" % method
    else:
        return method.__name__

def buildCalib(hmcs,name,title,method):
    
    gr = ROOT.TGraphErrors()
    gr.SetName(name)
    gr.SetTitle(title)
    
    for y in sorted(hmcs.keys()):
        x,xerr = getMean(hmcs[y],method)
        ip = gr.GetN()
        gr.SetPoint(ip,x,y)
        gr.SetPointError(ip,xerr,0.)
        
    return gr
    
def main(options,args):

    infile = args[0]
    fin = ROOT.TFile.Open(infile)
    os.chdir(options.outdir)

    ROOT.gStyle.SetOptFit(1)

    results = []
    ## methods = [1,0.95,0.8,0.683,0.5,0.4,0.3,0.2]
    ## methods = [recFit,0.683,(85.,100.),(80.,100.),(86.,96.),(88.,94.)]
    methods = [recFit,1.,0.683,(85.,97.),(86.,96.),(88,94),(89,93),(90,92)]
    if len(options.groups) > 0:
        categories = [ [int(t) for t in g.split(",")] for g in options.groups ]
    else:
        categories = [ [i] for i in range(options.ncat) ]
    for igroup,group in enumerate(categories):

        hdata = fin.Get("th1f_data_mass_cat%d" % group[0]).Clone("data_cat%d" % igroup)

        hmcs = {}
        hmcs[0] = fin.Get("th1f_sig_dymm_mass_m90_cat%d" % group[0]).Clone("mc_cat%d" % igroup)
        
        for isig in range(1,options.nsigma+1):
            hmcs[options.step*float(isig)] = fin.Get("th1f_sig_dymm_mass_m90_cat%d_E_scaleUp0%d_sigma" % (group[0],isig) ).Clone("mc_cat%d_Up0%d" % (igroup,options.step*float(isig)*10.))
            hmcs[-options.step*float(isig)] = fin.Get("th1f_sig_dymm_mass_m90_cat%d_E_scaleDown0%d_sigma" % (group[0],isig) ).Clone("mc_cat%d_Down0%d" % (igroup,options.step*float(isig)*10.))

        for icat in group[1:]:
            hdata.Add( fin.Get("th1f_data_mass_cat%d" % icat) )

            hmcs[0].Add(fin.Get("th1f_sig_dymm_mass_m90_cat%d" % icat))
        
            for isig in range(1,options.nsigma+1):
                hmcs[options.step*float(isig)].Add(fin.Get("th1f_sig_dymm_mass_m90_cat%d_E_scaleUp0%d_sigma" % (icat,isig) ))
                hmcs[-options.step*float(isig)].Add(fin.Get("th1f_sig_dymm_mass_m90_cat%d_E_scaleDown0%d_sigma" % (icat,isig) ))
                

        icat = igroup
        hdata.SetLineColor(ROOT.kBlack)
        hdata.SetMarkerColor(ROOT.kBlack)
        for s,h in hmcs.iteritems():
            h.SetLineColor(ROOT.kRed)
        hmcs[0].SetLineColor(ROOT.kBlue)
        
        ## hmcs[0].GetXaxis().SetRangeUser(80.,100.)
        ROOT.gStyle.SetOptFit(1)
        res = {}
        for method in methods:
            name = "mVsDeltaE_%s" % toStr(method)
                
            data = ( getMean(hdata,method), getMean(hmcs[0],method) )
            res[method] = ( buildCalib(hmcs,name,";m_{ll#gamma} (GeV/c^{2});#Delta E_{#gamma} / E_{#gamma} (#times 10^{-2})", method ),
                            data )

        ROOT.gStyle.SetOptFit(0)
        hdata.GetListOfFunctions().Clear()
        for s,h in hmcs.iteritems():
            h.GetListOfFunctions().Clear()

        canv = ROOT.TCanvas("cat_%d" % icat, "cat_%d" % icat )
        canv.cd()

        ndata = hdata.Integral()
        hmcs[0].Rebin(4)
        hdata.Rebin(4)
        h = hmcs[0].DrawNormalized("hist",ndata)
        st = setStat(h)
        objs.extend((st,h))
        
        hdata.SetLineColor(ROOT.kBlack)
        hdata.SetMarkerColor(ROOT.kBlack)
        hdata.Draw("e sames")
        objs.append(setStat(hdata,st,0.,-1))

        ip = 1
        for s,h in hmcs.iteritems():
            if fabs(s)  != 0.4:
                continue
            h.Rebin(4)
            h.SetLineColor(ROOT.kRed)
            if s < 0:
                h.SetLineStyle(ROOT.kDashed)
            h = h.DrawNormalized("histsames",ndata)
            q = setStat(h,st,-1.*ip)
            objs.extend([q,h])
            ip += 1
        objs.append( (canv, hdata, hmcs ))
        
        for fmt in "C","png","pdf":
            canv.SaveAs("%s.%s" % ( canv.GetName(), fmt) )
        
        results.append(res)
        
    objs.append(results)

    meas = {}
    out = open("results.txt","w+")
    for method in methods:
        scales = []
        name = "results_%s" % toStr(method)
        canv = ROOT.TCanvas(name,name,2800,2000)
        canv.Divide(options.ncat/2,2)
        objs.append(canv)
        for icat,res in enumerate(results):
            canv.cd(icat+1)
            ROOT.gPad.SetGridx()
            ROOT.gPad.SetGridy()
            calib, data = res[method]
            data,mc = data
            calib.Fit("pol1","W+")
            func = calib.GetListOfFunctions().At(0)
            scale = func.Eval(data[0])
            err = sqrt(data[1]**2 + mc[1]**2)
            scaleP = func.Eval(data[0]+err)
            scaleM = func.Eval(data[0]-err)
            scales.append( (data[0],mc[0],100.*(1.-data[0]/mc[0]),100.*(err/mc[0]),0.5*(scaleP-scaleM),(data[0]-mc[0])/err,scale,scaleM,scaleP) )
            ip = calib.GetN()
            calib.SetPoint( ip, data[0], scale )
            calib.SetPointError( ip, err, 0.5*(scaleP-scaleM) )
            calib.SetMarkerStyle(10)
            calib.Draw("ap")
        for fmt in "C","png","pdf":
            canv.SaveAs("%s.%s" % (canv.GetName(),fmt) )
        out.write( "Method: %s\n" % toStr(method ) )
        for icat, scale in enumerate(scales):
            out.write("cat: %d " % icat )
            for num in scale:
                out.write(("%.4g" % num).ljust(10) )
            out.write("\n")
    out.close()
    out = open("results.txt")
    print out.read()
    
        
if __name__ == "__main__":
    parser = OptionParser(option_list=[
        make_option("-s", "--step",
                    action="store", type="float", dest="step",
                    default=0.2,
                    ),
        make_option("-N", "--nsigma",
                    action="store", type="int", dest="nsigma",
                    default=4,
                    ),
        make_option("-n", "--ncat",
                    action="store", type="int", dest="ncat",
                    default=8,
                    ),
        make_option("-g", "--group", 
                    action="append",  dest="groups",
                    default=[],
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
    ## options.groups.extend( [ "0,1", "2,3", "4,5", "6,7" ] )
    ## options.groups.extend( [ "0,4", "1,5", "2,6", "3,7" ] )
    ## options.groups.extend( [ "0,1,4,5", "2,3,6,7" ] )

    try:
        os.mkdir(options.outdir)
    except:
        pass

    print options
    
    sys.argv.append("-b")
    import ROOT
    
    main( options, args )
