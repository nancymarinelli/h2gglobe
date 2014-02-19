#!/bin/env python

import sys
import os
import ROOT

from math import fabs

def tdrGrid(tdrStyle, gridOn):
#def tdrGrid(TStyle tdrStyle, bool gridOn):
    tdrStyle.SetPadGridX(gridOn)
    tdrStyle.SetPadGridY(gridOn)
    return

# fixOverlay: Redraws the axis

def fixOverlay():
  gPad.RedrawAxis()
  return

def setTDRStyle():
  tdrStyle=ROOT.TStyle("tdrStyle","Style for P-TDR")

# For the canvas:
  tdrStyle.SetCanvasBorderMode(0) 
  tdrStyle.SetCanvasColor(ROOT.kWhite) 
  tdrStyle.SetCanvasDefH(600)  #Height of canvas
  tdrStyle.SetCanvasDefW(800)  #Width of canvas
  tdrStyle.SetCanvasDefX(0)    #POsition on screen
  tdrStyle.SetCanvasDefY(0) 

# For the Pad:
  tdrStyle.SetPadBorderMode(0) 
  # tdrStyle.SetPadBorderSize(Width_t size = 1) 
  tdrStyle.SetPadColor(ROOT.kWhite) 
  tdrStyle.SetPadGridX(False) 
  tdrStyle.SetPadGridY(False) 
  tdrStyle.SetGridColor(0) 
  tdrStyle.SetGridStyle(3) 
  tdrStyle.SetGridWidth(1) 

# For the frame:
  tdrStyle.SetFrameBorderMode(0) 
  tdrStyle.SetFrameBorderSize(1) 
  tdrStyle.SetFrameFillColor(0) 
  tdrStyle.SetFrameFillStyle(0) 
  tdrStyle.SetFrameLineColor(1) 
  tdrStyle.SetFrameLineStyle(1) 
  tdrStyle.SetFrameLineWidth(1) 

# For the histo:
  # tdrStyle.SetHistFillColor(1) 
  # tdrStyle.SetHistFillStyle(0) 
  tdrStyle.SetHistLineColor(1) 
  tdrStyle.SetHistLineStyle(0) 
  tdrStyle.SetHistLineWidth(1) 
  # tdrStyle.SetLegoInnerR(Float_t rad = 0.5) 
  # tdrStyle.SetNumberContours(Int_t number = 20) 

  tdrStyle.SetEndErrorSize(2) 
  #tdrStyle.SetErrorMarker(20)   # Seems to give an error
  tdrStyle.SetErrorX(0.) 
  
  tdrStyle.SetMarkerStyle(20) 

#For the fit/function:
  tdrStyle.SetOptFit(0) 
  tdrStyle.SetFitFormat("5.4g") 
  tdrStyle.SetFuncColor(2) 
  tdrStyle.SetFuncStyle(1) 
  tdrStyle.SetFuncWidth(1) 

#For the date:
  tdrStyle.SetOptDate(0) 
  # tdrStyle.SetDateX(Float_t x = 0.01) 
  # tdrStyle.SetDateY(Float_t y = 0.01) 

# For the statistics box:
  tdrStyle.SetOptFile(0) 
  tdrStyle.SetOptStat(0)  # To display the mean and RMS:   SetOptStat("mr") 
  tdrStyle.SetStatColor(ROOT.kWhite) 
  tdrStyle.SetStatFont(42) 
  tdrStyle.SetStatFontSize(0.025) 
  tdrStyle.SetStatTextColor(1) 
  tdrStyle.SetStatFormat("6.4g") 
  tdrStyle.SetStatBorderSize(1) 
  tdrStyle.SetStatH(0.1) 
  tdrStyle.SetStatW(0.15) 
  # tdrStyle.SetStatStyle(Style_t style = 1001) 
  # tdrStyle.SetStatX(Float_t x = 0) 
  # tdrStyle.SetStatY(Float_t y = 0) 

# Margins:
  tdrStyle.SetPadTopMargin(0.05) 
  tdrStyle.SetPadBottomMargin(0.13) 
  tdrStyle.SetPadLeftMargin(0.16) 
  tdrStyle.SetPadRightMargin(0.10) 

# For the Global title:
  tdrStyle.SetOptTitle(0)     # 0=No Title
  tdrStyle.SetTitleFont(42) 
  tdrStyle.SetTitleColor(1) 
  tdrStyle.SetTitleTextColor(1) 
  tdrStyle.SetTitleFillColor(10) 
  tdrStyle.SetTitleFontSize(0.07) 
  # tdrStyle.SetTitleH(0)  # Set the height of the title box
  # tdrStyle.SetTitleW(0)  # Set the width of the title box
  # tdrStyle.SetTitleX(0)  # Set the position of the title box
  # tdrStyle.SetTitleY(0.985)  # Set the position of the title box
  # tdrStyle.SetTitleStyle(Style_t style = 1001) 
  # tdrStyle.SetTitleBorderSize(2) 

# For the axis titles:
  tdrStyle.SetTitleColor(1, "XYZ") 
  tdrStyle.SetTitleFont(42, "XYZ") 
  tdrStyle.SetTitleSize(0.08, "XYZ") 
  # tdrStyle.SetTitleXSize(Float_t size = 0.02)  # Another way to set the size?
  # tdrStyle.SetTitleYSize(Float_t size = 0.02) 
  tdrStyle.SetTitleXOffset(1.5) 
  tdrStyle.SetTitleYOffset(1.5)
  # tdrStyle.SetTitleOffset(1.1, "Y")  # Another way to set the Offset

# For the axis labels:
  tdrStyle.SetLabelColor(1, "XYZ") 
  tdrStyle.SetLabelFont(42, "XYZ") 
  tdrStyle.SetLabelOffset(0.007, "XYZ") 
  tdrStyle.SetLabelSize(0.05, "XYZ") 

# For the axis:
  tdrStyle.SetAxisColor(1, "XYZ") 
  tdrStyle.SetStripDecimals(ROOT.kTRUE) 
  tdrStyle.SetTickLength(0.03, "XYZ") 
  tdrStyle.SetNdivisions(510, "XYZ") 
  tdrStyle.SetPadTickX(0)   # 0=Text labels (and tics) only on bottom, 1=Text labels on top and bottom
  tdrStyle.SetPadTickY(1) 

# Change for log plots:
  tdrStyle.SetOptLogx(0) 
  tdrStyle.SetOptLogy(0) 
  tdrStyle.SetOptLogz(0) 

# Postscript options:
  tdrStyle.SetPaperSize(20.,20.) 
  # tdrStyle.SetLineScalePS(Float_t scale = 3) 
  # tdrStyle.SetLineStyleString(Int_t i, const char* text) 
  # tdrStyle.SetHeaderPS(const char* header) 
  # tdrStyle.SetTitlePS(const char* pstitle) 

  # tdrStyle.SetBarOffset(Float_t baroff = 0.5) 
  # tdrStyle.SetBarWidth(Float_t barwidth = 0.5) 
  # tdrStyle.SetPaintTextFormat(const char* format = "g") 
  # tdrStyle.SetPalette(Int_t ncolors = 0, Int_t* colors = 0) 
  # tdrStyle.SetTimeOffset(Double_t toffset) 
  # tdrStyle.SetHistMinimumZero(kTRUE) 

  #gROOT.ForceStyle()   # Try this if stuff doesn't work right
  
  # Find RooFit include
  #gSystem.SetIncludePath("-I$ROOTSYS/include -I$ROOFITSYS/include") 

  tdrStyle.cd() 
  return


#
# Handy histogram manipulation
# 
def applyModifs(h,modifs):
    for method in modifs:
        args = None
        ret = None
        if type(method) == tuple:
            method, args = method
        if type(method) == str:
            if hasattr(h,method):
                method = getattr(h,method)
            else:
                method = globals()[method]
        if args == None:            
            try:
                ret = method(h)
            except:
                ret = method()
        else:            
            if not ( type(args) == tuple or type(args) == list ):
                args = tuple([args])
            try:
                ret = method(h,*args)
            except:
                ret = method(*args)
        if ret and ret != h:
            print "Replacing ", h, ret, method
            h = ret
            ret = None
    return h
# 
def setcolors(h,color):
    h.SetLineColor(color)
    h.SetFillColor(color)
    h.SetMarkerColor(color)
# 
def legopt(h,opt):
    h.legopt = opt
# 
def xrange(h,a,b):
    h.GetXaxis().SetRangeUser(a,b)
# 
def xtitle(h,tit):
    h.GetXaxis().SetTitle(tit)
# 
def ytitle(h,tit):
    h.GetYaxis().SetTitle( (tit % { "binw" : h.GetBinWidth(0), "xtitle" : h.GetXaxis().GetTitle() }).replace("(GeV)","") )

def rename(h,cmap):
    toks = h.GetName().split("_")
    cat = None
    for t in toks:
        if t in cmap:
            cat = t
            break
    if cat:
        newname = h.GetName().replace(cat,cmap[cat])
        h.SetName(newname)
#
# Read plots from globe histogram files
#
def readPlot(fin, name, cat=0, which=["_sequential","_nminus1"], samples=["diphojet_8TeV","ggh_m125_8TeV","vbf_m125_8TeV",]):

    ret = []
    for p in which:
        hists = []
        for s in samples:
            nam = "%s%s_cat%d_%s" % ( name, p, cat, s )
            h = fin.Get(nam)
            print nam
            h.GetXaxis().SetTitle(  h.GetXaxis().GetTitle().replace("@"," ") )
            ## print nam, h
            hists.append(h)
        ret.append(hists)
            
    return ret


#
# Dump event yields for different samples
#
def eventYield(filenames,categories=[5,6],procs=["ggh_m125_8TeV","vbf_m125_8TeV","diphojet_8TeV"]):
    files = [ ROOT.TFile.Open(f) for f in filenames ]

    
    for fin in files:
        print fin.GetName()
        for cat in categories:
            plots=readPlot(fin, "mass", cat=cat, which=[""], samples=procs)[0]
            print "cat%d " % cat-1,
            for iproc in range(len(procs)):
                ## print "%1.3g" % plots[iproc].Integral(18,28), 115-135
                print "%1.3g" % plots[iproc].Integral(18,28),
            print

#
# Read histograms for a given process, applying manipulators
#
def readProc(fin,name,title,style,subproc,plot,plotmodifs,category):

    names = subproc.keys()
    print fin, plot, names, name
    histos = readPlot(fin,plot,which=[""],samples=names,cat=category)[0]
    print histos
    for iplot in range(len(histos)):
        h = histos[iplot]
        hname = names[iplot]
        print h,hname
        h = applyModifs(h,subproc[hname])

    print len(histos)
    sum = histos[0].Clone(name)
    sum.SetTitle(title)
    
    for h in histos[1:]:
        sum.Add(h)
        
    sum = applyModifs(sum,plotmodifs)
    sum = applyModifs(sum,style)
    
    return sum

#
# Prepare canvas and legend
#
def makeCanvAndLeg(name,legPos):
    canv = ROOT.TCanvas(name)
    leg  = ROOT.TLegend(*legPos)

    leg.SetFillStyle(0), leg.SetLineColor(ROOT.kWhite)## , leg.SetShadowColor(ROOT.kWhite)

    return canv, leg

def makeLegend(legPos):
    leg  = ROOT.TLegend(*legPos)
    leg.SetFillStyle(0), leg.SetLineColor(ROOT.kWhite)## , leg.SetShadowColor(ROOT.kWhite)

    return leg

#
# Make THStack out of python list
#
def makeStack(name,histos):
    stk = ROOT.THStack()
    for h in histos:
        stk.Add(h)
    return stk

def stackTitles(stk):
    stk.GetHistogram().GetXaxis().SetTitle( stk.GetHists()[0].GetXaxis().GetTitle() )
    stk.GetHistogram().GetYaxis().SetTitle( stk.GetHists()[0].GetYaxis().GetTitle() )

#
# Make THStack out of python list
#
def makeEnvelope(name,histos,stPlus=None,stMinus=None):
    nominal = histos[0]
    errPlus  = nominal.Clone( "%s_ErrPlus" % name )
    errMinus = nominal.Clone( "%s_ErrMinus" % name )
    if stPlus:
        applyModifs( errPlus, stPlus )
        applyModifs( errMinus, stMinus )
    for ibin in range(nominal.GetNbinsX()):
        hist = ROOT.TH1F("hist","hist",11,-5.,5.)
        hist.Reset("ICE")
        hist.Sumw2()
        points = []
        
        plus  = nominal.GetBinContent(ibin+1)
        minus = nominal.GetBinContent(ibin+1)
        nom   = nominal.GetBinContent(ibin+1)
        err   = nominal.GetBinError(ibin+1)
        points.append( [nom,err] )
        hist.Fill(nom)
        for h in histos[1:]:
            content =  h.GetBinContent(ibin+1)
            err     =  h.GetBinError(ibin+1)
            hist.Fill(content)
            points.append( [content,err] )
            if content < minus:
                minus = content
            if content > plus:
                plus = content

        if hist.GetRMS() == 0.:
            errPlus.SetBinContent(ibin+1,plus)
            errMinus.SetBinContent(ibin+1,minus)
            continue
            
        hist2 = ROOT.TH1F("hist2","hist2",11,hist.GetMean()-5.*hist.GetRMS(),hist.GetMean()+5.*hist.GetRMS())
        hist2.Sumw2()
        for p,e in points:
            hist2.Fill(p)
            
        func = ROOT.TF1("func","[0]*exp( -0.5*pow( (x-[1])/( (x>=0)*[2] + (x<=0)*[3] ) ,2.) )",hist2.GetMean()-5.*hist2.GetRMS(),hist2.GetMean()+5.*hist2.GetRMS())
        func.SetParameters(len(histos),hist2.GetMean(),hist2.GetRMS(), hist2.GetRMS())
        ## func.SetParLimits(2,-0.5*hist.GetMean(),0.5*hist.GetMean())
        ## func.SetParLimits(3,-0.1*hist.GetMean(),0.1*hist.GetMean())
        stat = int( hist2.Fit( func, "L" ) )
        if stat == 0:
            errPlus.SetBinContent(ibin+1, max(nom,func.GetParameter(1)+fabs(func.GetParameter(2))))
            errMinus.SetBinContent(ibin+1,min(nom,func.GetParameter(1)-fabs(func.GetParameter(3))))
        else:
            errPlus.SetBinContent(ibin+1,plus)
            errMinus.SetBinContent(ibin+1,minus)
        
    return makeStack(name,[errPlus,errMinus,nominal])


def shiftHisto(offset,nominal,shifted):
    for ibin in range(offset.GetNbinsX()):
        shifted.SetBinContent(ibin+1, shifted.GetBinContent(ibin+1) + offset.GetBinContent(ibin+1) - nominal.GetBinContent(ibin+1) )
        shifted.SetBinError(ibin+1, ROOT.TMath.Sqrt( shifted.GetBinError(ibin+1)*shifted.GetBinError(ibin+1) - 
                                                     nominal.GetBinError(ibin+1)*nominal.GetBinError(ibin+1) ) )
        ## shifted.SetBinContent(ibin+1, shifted.GetBinContent(ibin+1) - nominal.GetBinContent(ibin+1) )
    return shifted

def ratioHisto(num,den,ytit):
    num.Divide(den)
    ytitle(num,ytit)
    return num

#
# Draw a THStack
#
def drawStack(stk, method, option):
    ymax = 0.
    if "DrawNormalized" in method:
        rng = None
        if "[" in method:
            rng = [ float(f) for f in method.split("DrawNormalized")[1].split("[")[1].split("]")[0].split(",") ]
            print rng
        histos = [ stk.GetStack().At(0) ]
        if "nostack" in option:
            histos = stk.GetHists()
            option = option.replace("nostack","")
        for h in histos:
            h.SetFillStyle(0)
            h.SetLineWidth(2)
            bmin = -1
            bmax = -1
            if rng:
                bmin = h.FindBin(rng[0])
                bmax = h.FindBin(rng[1])
            h.Scale(1./h.Integral(bmin,bmax))
            h.Draw("%s SAME" % option)
            ymax = max(ymax, h.GetMaximum()) 
    else:
        getattr(stk,method.split(",")[0])("%s" % option)
        ymax = stk.GetMaximum(option)
        print stk.GetName(), stk.GetHists()[0].Integral()
        
    return ymax

#
# Perform data/MC comparison for many plots and categories
#
def dataMcComparison(data, bkg, sig, plots, categories=[0], savefmts=["C","png","pdf"], catlabels={}):

    objs = []
    canvs = []
    # loop over categories
    for cat in categories:
        print cat
        # loop over plots
        for plot in plots:

            plotname, plotmodifs, drawopts, legPos = plot
            dm, dataopt, bkgopt, sigopt = drawopts
            
            bkghists = []
            sighists = []
            datahists = []
            
            # read background MC
            if bkg != None:
                bkgfile, bkgprocs = bkg
                bkghists = [ readProc(bkgfile,*subprocs,plot=plotname,plotmodifs=plotmodifs,category=cat) for subprocs in bkgprocs ]
                
            # read signal MC
            if sig != None:
                sigfile, sigprocs = sig
                sighists = [ readProc(sigfile,*subprocs,plot=plotname,plotmodifs=plotmodifs,category=cat) for subprocs in sigprocs ]
                
            # read data
            if data != None:
                datafile, dataprocs = data
                datahists = [ readProc(datafile,*subprocs,plot=plotname,plotmodifs=plotmodifs,category=cat) for subprocs in dataprocs ]
                
            # collect histograms
            allhists = datahists+bkghists+sighists
            objs += allhists
            
            # make empty frame histogram for convenience
            frame = allhists[0].Clone()
            frame.Reset("ICE")
            frame.SetEntries(0)
            objs.append(frame)
            ymax = 0.
            ymin = 0.

            # allocate canvas and legend and draw frame
            catname = "cat%d" % cat
            if catname in catlabels:
                catname = catlabels[catname]
            canv,leg = makeCanvAndLeg("%s_%s" % ( plotname, catname), legPos )
            objs.append(canv)
            objs.append(leg)
            canvs.append(canv)
            frame.Draw()

            # draw background first
            if len(bkghists) > 0:
                bkgstk = makeStack("bkg_%s_cat%d" % ( plotname, cat), bkghists)
                ### getattr(bkgstk,dm)("%s SAME" % bkgopt)
                ### ymax = max(ymax,bkgstk.GetMaximum(bkgopt))
                ymax = max(ymax,drawStack(bkgstk,dm,"%s SAME"%bkgopt))
                objs.append(bkgstk)
                
            # then data
            if len(datahists) > 0:
                datastk = makeStack("data_%s_cat%d" % ( plotname, cat),datahists)
                ### getattr(datastk,dm)("%s SAME" % dataopt)
                ### ymax = max(ymax,datastk.GetMaximum())
                ymax = max(ymax,drawStack(datastk,dm,"%s SAME"%dataopt))
                objs.append(datastk)

            # and finally signal
            if len(sighists) > 0:
                sigstk = makeStack("sig_%s_cat%d" % ( plotname, cat),sighists)
                ### getattr(sigstk,dm)("%s SAME" % sigopt)
                ### ymax = max(ymax,sigstk.GetMaximum(sigopt))
                ymax = max(ymax,drawStack(sigstk,dm,"%s SAME"%sigopt))
                objs.append(sigstk)

            # make legend
            for h in allhists:
                legopt = "f"
                if hasattr(h,"legopt"):
                    legopt = h.legopt
                leg.AddEntry(h,"",legopt)

            # adjust yaxis
            frame.GetYaxis().SetRangeUser(ymin,ymax*1.2)
            leg.Draw("same")
            canv.RedrawAxis()

            # if needed draw inset with zoom-in
            if "DrawInset" in dm:
                inset =  [ float(f) for f in dm.split("DrawInset")[1].split("[")[1].split("]")[0].split(",") ]
                rng = inset[0:2]
                pos = inset[2:]
                
                padname = "%s_cat%d_inset" % ( plotname, cat)
                pad = ROOT.TPad(padname, padname, *pos)
                objs.append(pad)
                pad.Draw("")
                pad.SetFillStyle(0)
                
                pad.cd()
                padframe = frame.Clone()
                padframe.GetXaxis().SetRangeUser(*rng)
                padframe.GetYaxis().SetRangeUser(ymin,ymax*1.2)
                padframe.Draw()
                objs.append(padframe)
                
                if len(bkghists) > 0:
                    drawStack(bkgstk,"Draw",bkgopt+" same")
                
                if len(datahists) > 0:
                    drawStack(datastk,"Draw",dataopt+" same")
                    
                if len(sighists) > 0:
                    drawStack(sigstk,"Draw",sigopt+" same")

                pad.RedrawAxis()
                
    # save plots
    for c in canvs:
        for fmt in savefmts:
            try:
                c.SaveAs("%s.%s" % (c.GetName(),fmt))
            except:
                print c
                
            
    return objs

if __name__ == "__main__":

    #
    # open input files
    #
    fdata = ROOT.TFile.Open(sys.argv[1])
    fsig  = fdata
    
    #
    # prepare output folder
    #
    outdir = sys.argv[2]
    try:
        os.mkdir(outdir)
    except:
        pass
    os.chdir(outdir)

    #
    # draw options and style
    # 
    defdrawopt = ("Draw","e","hist","hist nostack") ## drawing method(s), drawing option data, MC background, MC signal
    mvadrawopt = ("Draw,DrawInset[0.85,1,0.45,0.48,0.9,0.92]","e","hist","hist nostack")
    deflegPos  = (0.137,0.725,0.374,0.888)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetOptStat(0)
    
    #
    # datasets to be plot
    #

    # data
    data = [ fdata, [ ("data", "data",
                       [("SetMarkerStyle",(ROOT.kFullCircle)),(legopt,"pl")],
                       { "Data" : [] }
                       )
                      ]
             ]
    # background
    sig  = None
    # signal
    bkg = [ fsig,  [ ("Zmmg","Z #rightarrow #mu #mu #gamma",
                        [(setcolors,ROOT.kBlue),("SetLineWidth",2),("SetFillStyle",1),("Scale",1),(legopt,"l")],
                        { "dymm_m90" : [] }
                        )
                       ]
              ]

    catlabels = {"cat0": "all",
                 "cat1": "eb",
                 "cat2": "ee",
                 "cat3": "eb_highpt",
                 "cat4": "eb_lowpt",
                 "cat5": "ee_highpt",
                 "cat6": "ee_lowpt"
                 }

    
    # Make data/MC comparison
    objs=dataMcComparison( data = data,
                           sig = sig,
                           bkg = bkg,
                           categories = [0,1,2,3,4,5,6],
                           catlabels = catlabels,
                           plots = [ ### "vbf_mva",[("SetBinContent",(1,0.)),("Rebin",4),(xrange,(-1.,1)),
                                     ###            ("SetBinError",(1,0.)),(xtitle,"MVA"),## (ytitle,"A.U.")
                                     ###            (ytitle,"Events/%(binw)1.2g")
                                     ###            ],
                                     ### mvadrawopt,deflegPos),
                                     
                                     ("mmg_pt"   ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mmg_eta"  ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mmg_phi"  ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mmg_mass" ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("nvtx"     ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("pho_n"    ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("pho_e"    ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("pho_pt"   ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("pho_eta"  ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("pho_phi"  ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("pho_sce"  ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("pho_sceta",[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("pho_r9"   ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mu1_pt"   ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mu1_eta"  ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mu1_phi"  ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mu2_pt"   ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mu2_eta"  ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mu2_phi"  ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mumu_pt"  ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mumu_eta" ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mumu_phi" ,[(rename,catlabels)],defdrawopt,deflegPos),
                                     ("mumu_mass",[(rename,catlabels)],defdrawopt,deflegPos),
                                     ]
                           )

    
    
    ## eventYield(sys.argv[1:])
