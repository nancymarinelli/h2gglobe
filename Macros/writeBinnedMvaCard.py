import ROOT
import os,numpy,sys,math
ROOT.gROOT.SetStyle("Plain")
ROOT.gROOT.SetBatch(True)

ROOT.gROOT.ProcessLine(".L quadInterpolate.C+")
from ROOT import quadInterpolate

from optparse import OptionParser

r=ROOT.TRandom3(0)

# Some "Global" Variables
# PLOT OPS ----------------
lumistring = "4.76 fb^{-1}"
sigscale = 5.
# THEORY SYSTEMATICS ------
lumi 		= "1.045"
QCDscale_ggH  = "0.918/1.125"
PDF_gg_1      = "0.923/1.079"
PDF_gg_2      = "0.915/1.085"
QCDscale_qqH  = "0.997/1.005"
PDF_qqbar_1   = "0.979/1.027"
PDF_qqbar_2   = "0.958/1.042" 
QCDscale_VH   = "0.982/1.018"
QCDscale_ttH  = "0.905/1.036"
# SHAPE SYSTEMATICS -------
systematics = [
	       "E_res"
	      ,"E_scale"
	      ,"idEff"
	      ,"phoIdMva"
	      ,"kFactor"
	      ,"triggerEff"
	      ,"vtxEff"
	      ]
# ADDITIONAL SYSTEMATICS --
JetID_vbf = 0.1
JetID_ggh = 0.7
# -------------------------

def plainBin(hist):
	nb = hist.GetNbinsX()
	h2 = ROOT.TH1F(hist.GetName()+"new","",nb,0,nb)
	for i in range (1,nb+1):
		h2.SetBinContent(i,hist.GetBinContent(i))
		h2.SetBinError(i,hist.GetBinError(i))
		if (options.includeVBF):
			#h2.GetXaxis().SetLabelSize(0.03)
			print hist.GetBinLowEdge(i)
			if hist.GetBinLowEdge(i+1) <= 1.:
			 # h2.GetXaxis().SetBinLabel(i,"%.4f < BDT < %.4f "%(hist.GetBinLowEdge(i),hist.GetBinLowEdge(i+1)))
			  h2.GetXaxis().SetBinLabel(i,"BDT Bin %d "%(i))
			else: 
			  h2.GetXaxis().SetBinLabel(i," Di-jet ")
	h2.GetXaxis().SetNdivisions(nb)
	return h2

def plotDistributions(mass,data,signals,bkg,errors):

	if options.splitSignal: # last signal is separated off
	  for i in range(1,len(signals)-1):
		signals[0].Add(signals[i])
	else:
	  for i in range(1,len(signals)):
		signals[0].Add(signals[i])

	nbins = data.GetNbinsX()

	flatdata = plainBin(data)
	flatsignal  = plainBin(signals[0])
	flatsignal1  = plainBin(signals[-1])

	flatbkg  = plainBin(bkg)
	flatbkg.SetLineColor(4)
	flatbkg.SetLineWidth(2)

	fNew  = flatbkg.Clone()
	fNew2 = flatbkg.Clone()

	#flatdata.Sumw2();
	flatdata.SetMarkerStyle(20)
	flatdata.SetMarkerSize(1.0)
		
	fNew.SetLineColor(1)
	fNew.SetLineWidth(2)
	fNew.SetFillStyle(1001)
	fNew.SetFillColor(3)
	
	fNew2.SetFillColor(5)
	fNew2.SetLineColor(1)
	fNew2.SetLineWidth(2)
	fNew2.SetFillStyle(1001)
	fNewT = fNew.Clone()
	fNew2T = fNew2.Clone()
	fNewT.SetFillStyle(1001)
	fNew2T.SetFillStyle(1001)

	flatsignal.SetLineWidth(2)
	flatsignal.SetLineColor(ROOT.kRed)
#	flatsignal.SetFillColor(ROOT.kRed-10)
#	flatsignal.SetFillStyle(1001)
	flatsignal.Scale(sigscale)
	flatsignal1.SetLineWidth(2)
	flatsignal1.SetLineColor(ROOT.kGreen+4)
	flatsignal1.Scale(sigscale)
		
	leg = ROOT.TLegend(0.56,0.56,0.88,0.88)
	leg.SetFillColor(0)
	leg.SetBorderSize(0)

	for b in range(1,nbins+1):
		additional = errors[b-1]
  		fNew.SetBinError(b,((fNew.GetBinError(b)**2)+((fNew.GetBinContent(b)*additional)**2))**0.5)
  		fNew2.SetBinError(b,2*(((fNew2.GetBinError(b)**2)+((fNew2.GetBinContent(b)*additional)**2))**0.5))
  	#	flatbkg.SetBinError(b,((flatbkg.GetBinError(b)**2)+((flatbkg.GetBinContent(b)*additional)**2))**0.5)
	
	#c = ROOT.TCanvas("c","c",900,600)
	c = ROOT.TCanvas()
	c.SetLogy()
	if (not options.includeVBF): flatdata.GetXaxis().SetTitle("Category")
	flatdata.Draw("9")
	fNew2.Draw("9sameE2")
	fNew.Draw("9sameE2")
	flatbkg.Draw("9samehist")
	if options.splitSignal: 
	  sigst = ROOT.THStack()
	  sigst.Add(flatsignal)
	  sigst.Add(flatsignal1)
	  sigst.Draw("9samehist")
	else:
	  flatsignal.Draw("9samehist")
	flatdata.Draw("9sameP")
	flatdata.SetMinimum(1)

	leg.AddEntry(flatdata,"Data","PLE")
	if options.splitSignal:
	  leg.AddEntry(flatsignal,"Higgs (GG,WZ,TT), m_{H}=%3.1f GeV (x%d)"%(mass,int(sigscale)) ,"L")
	  leg.AddEntry(flatsignal1,"Higgs, m_{H}=%3.1f GeV (x%d)"%(mass,int(sigscale)) ,"L")

	else: leg.AddEntry(flatsignal,"Higgs, m_{H}=%3.1f GeV (x%d)"%(mass,int(sigscale)) ,"L")
	leg.AddEntry(flatbkg,"Background","L")
	leg.AddEntry(fNewT,"\pm 1\sigma","F")
	leg.AddEntry(fNew2T,"\pm 2\sigma","F")
	leg.Draw()
	mytext = ROOT.TLatex()
	mytext.SetTextSize(0.03)
	mytext.SetNDC()
#	mytext.DrawLatex(0.28,0.8,"#splitline{CMS preliminary}{#sqrt{s} = 7 TeV L = %s}"%(lumistring))
	mytext.DrawLatex(0.1,0.92,"CMS preliminary,  #sqrt{s} = 7 TeV ")
	
	mytext.SetTextSize(0.04)
	mytext.DrawLatex(0.2,0.2,"#int L = %s"%(lumistring))
	leg.Draw()
	c.SaveAs(plotOutDir+"/model_m%3.1f.pdf"%mass)
	c.SaveAs(plotOutDir+"/model_m%3.1f.png"%mass)
	
	d = ROOT.TCanvas()
	leg2 = ROOT.TLegend(0.56,0.56,0.88,0.88)
	leg2.SetFillColor(0)
	leg2.SetBorderSize(0)
	#d.SetLogy()
	if (not options.includeVBF): flatdata.GetXaxis().SetTitle("Category")
	flatdata.GetYaxis().SetTitle("Data - Background");
	datErrs = []
	for b in range(1,nbins+1): datErrs.append((flatdata.GetBinContent(b))**0.5);
	flatdata.Add(flatbkg,-1)
	for b in range(1,nbins+1): 
		flatdata.SetBinError(b,((datErrs[b-1]*datErrs[b-1]) +(fNew.GetBinError(b)*fNew.GetBinError(b)))**0.5 )
	flatbkg.Add(flatbkg,-1)

	flatdata.Draw("9")
	#fNew2.Draw("9sameE2")
	#fNew.Draw("9sameE2")
	flatbkg.Draw("9samehist")
	if options.splitSignal: 
	  sigst = ROOT.THStack()
	  sigst.Add(flatsignal)
	  sigst.Add(flatsignal1)
	  sigst.Draw("9samehist")
	else:
	  flatsignal.Draw("9samehist")
	flatdata.Draw("9sameP")
	flatdata.SetMaximum(250)
	flatdata.SetMinimum(-100)

	leg2.AddEntry(flatdata,"Data","PLE")
	if options.splitSignal:
	  leg2.AddEntry(flatsignal,"Higgs (GG,WZ,TT), m_{H}=%3.0f GeV (x%d)"%(mass,int(sigscale)) ,"L")
	  leg2.AddEntry(flatsignal1,"Higgs, m_{H}=%3.0f GeV (x%d)"%(mass,int(sigscale)) ,"L")

	else: leg2.AddEntry(flatsignal,"Higgs, m_{H}=%3.0f GeV (x%d)"%(mass,int(sigscale)) ,"L")
	leg2.AddEntry(flatbkg,"Background","L")

	mytext.SetTextSize(0.04)

	leg2.Draw()
	mytext = ROOT.TLatex()
	mytext.SetTextSize(0.03)
	mytext.SetNDC()
#	mytext.DrawLatex(0.28,0.8,"#splitline{CMS preliminary}{#sqrt{s} = 7 TeV L = %s}"%(lumistring))
	mytext.DrawLatex(0.1,0.92,"CMS preliminary,  #sqrt{s} = 7 TeV ")
	mytext.SetTextSize(0.04)
	mytext.DrawLatex(0.2,0.8,"#int L = %s"%(lumistring))
	d.SaveAs(plotOutDir+"/diff_model_m%3.1f.pdf"%mass)
	d.SaveAs(plotOutDir+"/diff_model_m%3.1f.png"%mass)
	


#def getBinningMass(mass):

#	if mass >= 115.0 and mass <= 117.0: return "115"
#	if mass >= 117.5 and mass <= 122.0: return "120"
#	if mass >= 122.5 and mass <= 127.0: return "125"
#	if mass >= 127.5 and mass <= 132.0: return "130"
#	if mass >= 132.5 and mass <= 137.0: return "135"
#	if mass >= 137.5 and mass <= 144.5: return "140"
#	if mass >= 145.0 and mass <= 150.0: return "150"

def py_quadInterpolate(C,X1,X2,X3,Y1,Y2,Y3):
	resL = quadInterpolate(-1*C,X1,X2,X3,Y1,Y2,Y3)
	resH = quadInterpolate(C,X1,X2,X3,Y1,Y2,Y3)
	if math.isnan(resL) or math.isinf(resL) or  math.isnan(resH) or math.isinf(resL): return " - "
	if abs(resL - 1) < 0.00001 or abs(resL - 1) > 1: return " - "
	if abs(resH - 1) < 0.00001 or abs(resH - 1) > 1: return " - "
	return " %.8f/%.8f "%(resL,resH) 

def getBinContent(hist,b):
  
	res = hist.GetBinContent(b)
	if res==0: return 0.00000001
	else: return res

def getPoissonBinContent(hist,b,exp):
  
	res = exp*(hist.GetBinContent(b))
	if res==0: return 0.00000001
	else: return r.Poisson(res)

def writeCard(tfile,mass,scaleErr):

  print "Writing Datacard for mass -> ", mass
  outPut = open(cardOutDir+"/mva-datacard_"+type+"_%3.1f.txt"%mass,"w")

  # Get All of the histograms we are going to use
  # Data ->
  dataHist = tfile.Get("th1f_data_"+type+"_%3.1f_cat0"%mass)
  nBins    = dataHist.GetNbinsX()
  print "Number of Channels -> ", nBins
  # bkg model ->
  bkgHist  	= tfile.Get("th1f_bkg_"+type+"_%3.1f_cat0"%mass)
  
  if options.Bias: bkgHistCorr   = tfile.Get("th1f_bkg_"+type+"_%3.1f_cat0_fitsb_biascorr"%mass)
  # 4 signal channels ->
  gghHist  = tfile.Get("th1f_sig_"+type+"_ggh_%3.1f_cat0"%mass)
  vbfHist  = tfile.Get("th1f_sig_"+type+"_vbf_%3.1f_cat0"%mass)
  wzhHist  = tfile.Get("th1f_sig_"+type+"_wzh_%3.1f_cat0"%mass)
  tthHist  = tfile.Get("th1f_sig_"+type+"_tth_%3.1f_cat0"%mass)

 
  ###############################################################################
  # Write the basics
  outPut.write("Mva Binned Analysis DataCard (mH=%3.1f) \n"%mass)
  outPut.write("DataCard extracted from %s \n"%tfile.GetName())
  outPut.write("--------------------------------------------------------------\n")
  outPut.write("imax *\n")
  outPut.write("jmax *\n")
  outPut.write("kmax *\n")
  outPut.write("--------------------------------------------------------------\n")
  ## Now its the observation
  outPut.write("bin 	     ")
  for b in range(1,nBins+1): outPut.write(" cat%d "%b)
  outPut.write("\nobservation")

  backgroundContents = [bkgHist.GetBinContent(b) for b in range(1,nBins+1)]
  if options.Bias:
	print "Using Bkg Model Corrected for mass bias"
	backgroundContents = [bkgHistCorr.GetBinContent(b) for b in range(1,nBins+1)]

  if options.throwToy:
	print "Throwing toy dataset"
	for b in range(1,nBins+1): 
	  nd = r.Poisson(backgroundContents[b-1])
	  ns = 0
	  if options.expSig>0:
		print "Injecting %.f x SM"%expSig
		ns+=getPoissonBinContent(gghHist,b,options.expSig)
		ns+=getPoissonBinContent(vbfHist,b,options.expSig)
		ns+=getPoissonBinContent(wzhHist,b,options.expSig)
		ns+=getPoissonBinContent(tthHist,b,options.expSig)

	  outPut.write(" %.2f "%(nd+ns))
	  dataHist.SetBinContent(b,nd+ns)
	  dataHist.SetBinError(b,(nd+ns)**0.5)

  else:
	for b in range(1,nBins+1): outPut.write(" %d "%dataHist.GetBinContent(b))
  outPut.write("\n--------------------------------------------------------------\n")
  ## Now we do the signal and background parts
  outPut.write("bin 	     ")
  for b in range(1,nBins+1): outPut.write(" cat%d  cat%d  cat%d  cat%d  cat%d "%(b,b,b,b,b))
  outPut.write("\nprocess    ")
  for b in range(1,nBins+1): outPut.write("  ggh    vbf    wzh    tth    bkg  ")
  outPut.write("\nprocess    ")
  for b in range(1,nBins+1): outPut.write("   0      0      0      0    1    ")
  outPut.write("\nrate       ")


  for b in range(1,nBins+1): outPut.write(" %.12f   %.12f   %.12f   %.12f   %.12f "\
    %(getBinContent(gghHist,b),getBinContent(vbfHist,b),getBinContent(wzhHist,b),getBinContent(tthHist,b)\
     ,backgroundContents[b-1]))
  outPut.write("\n--------------------------------------------------------------\n")


  # This next bit is for the signal systematics, first lets do the easy ones, lumi and theory
  outPut.write("\nlumi          lnN ")

  if options.theorySys:
    for b in range(1,nBins+1): outPut.write(" %s  %s  %s  %s  -  "%(lumi,lumi,lumi,lumi))
    outPut.write("\nQCDscale_ggH  lnN ")
    for b in range(1,nBins+1): outPut.write(" %s  -   -   -   -  "%(QCDscale_ggH))
    outPut.write("\nPDF_gg        lnN ")
    for b in range(1,nBins+1): outPut.write(" %s  -   -   %s  -  "%(PDF_gg_1,PDF_gg_2))
    outPut.write("\nQCDscale_qqH  lnN ")
    for b in range(1,nBins+1): outPut.write(" -   %s  -   -   -  "%(QCDscale_qqH))
    outPut.write("\nPDF_qqbar     lnN ")
    for b in range(1,nBins+1): outPut.write(" -   %s  %s  -   -  "%(PDF_qqbar_1,PDF_qqbar_2))
    outPut.write("\nQCDscale_VH   lnN ")
    for b in range(1,nBins+1): outPut.write(" -   -   %s  -   -  "%(QCDscale_VH))
    outPut.write("\nQCDscale_ttH  lnN ")
    for b in range(1,nBins+1): outPut.write(" -   -   -   %s  -  "%(QCDscale_ttH))

  outPut.write("\n")

  # includeVBF means the last bin is the VBF tagged bin and we apply and additional 
  # 70% GGH(TTH) and 10% on the VBF(WZH) part of that category (configurable above)
  if options.includeVBF:
    print "Including VBF (last Bin) Systematics"
    # calculate the effect on each bin, eg 70%, always assume last bin is VBF tag
    numberOfGGH_dijet = gghHist.GetBinContent(nBins)*JetID_ggh
    numberOfTTH_dijet = tthHist.GetBinContent(nBins)*JetID_ggh
    numberOfVBF_dijet = vbfHist.GetBinContent(nBins)*JetID_vbf
    numberOfWZH_dijet = wzhHist.GetBinContent(nBins)*JetID_vbf
    
    numberOfGGH_incl  = sum([gghHist.GetBinContent(b) for b in range(1,nBins)])
    numberOfTTH_incl  = sum([tthHist.GetBinContent(b) for b in range(1,nBins)])
    numberOfVBF_incl  = sum([vbfHist.GetBinContent(b) for b in range(1,nBins)])
    numberOfWZH_incl  = sum([wzhHist.GetBinContent(b) for b in range(1,nBins)])

    outPut.write("\nJetID_ggh  lnN ")
    for b in range(1,nBins): outPut.write(" %.3f/%.3f   -   -   %.3f/%.3f   -  "%\
		    (1.-(numberOfGGH_dijet/numberOfGGH_incl),1.+(numberOfGGH_dijet/numberOfGGH_incl),\
		     1.-(numberOfTTH_dijet/numberOfTTH_incl),1.+(numberOfTTH_dijet/numberOfTTH_incl)))
    outPut.write(" %.3f/%.3f   -   -   %.3f/%.3f   -  "%(1+JetID_ggh,1-JetID_ggh,1+JetID_ggh,1-JetID_ggh))
    outPut.write("\nJetID_vbf  lnN ")
    for b in range(1,nBins): outPut.write(" -  %.3f/%.3f  %.3f/%.3f  -   -  "%\
		    (1.-(numberOfVBF_dijet/numberOfVBF_incl),1.+(numberOfVBF_dijet/numberOfVBF_incl),\
		     1.-(numberOfWZH_dijet/numberOfWZH_incl),1.+(numberOfWZH_dijet/numberOfWZH_incl)))
    outPut.write(" -  %.3f/%.3f   %.3f/%.3f  -   -  "%(1+JetID_vbf,1-JetID_vbf,1+JetID_vbf,1-JetID_vbf))
    outPut.write("\n")

    
   
  # Now is the very tedious part of the signal shape systematics, for each shape, simply do -/+ sigma
  
  if options.signalSys:
   print "Writing Systematics Part (coule be slow)"
   for sys in systematics:

    gghHistU  = tfile.Get("th1f_sig_"+type+"_ggh_%3.1f_cat0_%sUp01_sigma"%(mass,sys))
    vbfHistU  = tfile.Get("th1f_sig_"+type+"_vbf_%3.1f_cat0_%sUp01_sigma"%(mass,sys))
    wzhHistU  = tfile.Get("th1f_sig_"+type+"_wzh_%3.1f_cat0_%sUp01_sigma"%(mass,sys))
    tthHistU  = tfile.Get("th1f_sig_"+type+"_tth_%3.1f_cat0_%sUp01_sigma"%(mass,sys))
    gghHistD  = tfile.Get("th1f_sig_"+type+"_ggh_%3.1f_cat0_%sDown01_sigma"%(mass,sys))
    vbfHistD  = tfile.Get("th1f_sig_"+type+"_vbf_%3.1f_cat0_%sDown01_sigma"%(mass,sys))
    wzhHistD  = tfile.Get("th1f_sig_"+type+"_wzh_%3.1f_cat0_%sDown01_sigma"%(mass,sys))
    tthHistD  = tfile.Get("th1f_sig_"+type+"_tth_%3.1f_cat0_%sDown01_sigma"%(mass,sys))

    outPut.write("\n%s lnN "%sys)
    for b in range(1,nBins+1): 
	 outPut.write(" %s %s %s %s - "%(\
				     py_quadInterpolate(1.,-3.,0.,3.,gghHistD.GetBinContent(b)  \
				        			  ,gghHist.GetBinContent(b)  \
                                        			  ,gghHistU.GetBinContent(b)) \
				    ,py_quadInterpolate(1.,-3.,0.,3.,vbfHistD.GetBinContent(b)  \
				        			  ,vbfHist.GetBinContent(b)  \
                                        			  ,vbfHistU.GetBinContent(b)) \
				    ,py_quadInterpolate(1.,-3.,0.,3.,wzhHistD.GetBinContent(b)  \
				        			  ,wzhHist.GetBinContent(b)  \
                                        			  ,wzhHistU.GetBinContent(b)) \
				    ,py_quadInterpolate(1.,-3.,0.,3.,tthHistD.GetBinContent(b)  \
				        			  ,tthHist.GetBinContent(b)  \
                                        			  ,tthHistU.GetBinContent(b)) \
 				    ))
  outPut.write("\n")
  # Finally the background errors, these are realtively simple
  outPut.write("\nbkg_norm lnN ")
  for b in range(1,nBins+1): outPut.write(" -   -   -   -  %.8f "%(scaleErr))

  ## now for the David errors
  if options.Bias:
	print "Including Mass Bias nuisances (bin-to-bin stat error included)"

	# Input Signed Error Matrix from Fit 
	th2f_errmatrix = tfile.Get("fUncorrErr_%s_%3.1f"%(options.bdtType,mass))
	for b in range(1,nBins):  # This error matrix is nBins-1 X nBins-1 due to constraint on sum on fractions
           outPut.write("\nmassBias%d lnN"%b)
	   for q in range(1,nBins+1):
	   	f_errentry = th2f_errmatrix.GetBinContent(b,q)
		bkgC = backgroundContents[q-1]/sum(backgroundContents)
		bias_nuis  = py_quadInterpolate(1.,-1.,0.,1.,bkgC-f_errentry,bkgC,bkgC+f_errentry)
		outPut.write(" - - - - %s "%(bias_nuis))
        	
        outPut.write("\n")
	

  if options.B2B and (not options.Bias):
   # bkg bins will be gmN errors instead 
   for b in range(1,nBins+1):
        bkgScale = bkgHist.Integral()/bkgHist.GetEntries()
        outPut.write("\nbkg_stat%d gmN %d "%(b,int(backgroundContents[b-1]/bkgScale)))
	for q in range(1,nBins+1):
		if q==b: outPut.write(" - - - - %.8f "%bkgScale)
		else:    outPut.write(" - - - - - ")

  # Finally make a plot of what will go into the limit
  if options.makePlot: 
	print "Plotting Limit Setting inputs to m%3.1f.pdf"%mass
	if options.Bias : plotBKG = bkgHistCorr.Clone()
	else: plotBKG = bkgHist.Clone()
	plotDistributions(mass,dataHist.Clone() \
			,[gghHist.Clone(),wzhHist.Clone(),tthHist.Clone(),vbfHist.Clone()]\
			,plotBKG \
			,[scaleErr-1. for b in range(nBins)])# put VBF at the end for plot option

  outPut.write("\n")  # done
  outPut.close()
    
#################################################################################  

parser = OptionParser()
parser.add_option("-i","--input",dest="tfileName")
parser.add_option("","--noB2B",action="store_false",dest="B2B",default=True)
#parser.add_option("","--addBias",dest="biasFile",default=None)
parser.add_option("","--noBias",dest="Bias",default=True,action="store_false")
parser.add_option("","--noSignalSys",action="store_false",dest="signalSys",default=True)
parser.add_option("","--noTheorySys",action="store_false",dest="theorySys",default=True)
parser.add_option("","--throwToy",action="store_true",dest="throwToy",default=False)
parser.add_option("","--expSig",dest="expSig",default=-1.,type="float")
parser.add_option("","--makePlot",dest="makePlot",default=False,action="store_true")
parser.add_option("","--noVbfTag",dest="includeVBF",default=True,action="store_false")
parser.add_option("","--plotCombineVbf",dest="splitSignal",default=True,action="store_false")
parser.add_option("-m","--mass",dest="singleMass",default=-1.,type="float")
parser.add_option("-t","--type",dest="bdtType",default="grad");


(options,args)=parser.parse_args()
print "Creating Binned Datacards from workspace -> ", options.tfileName
if options.throwToy: print ("Throwing Toy dataset from BKG")
if options.expSig > 0: print ("(Also throwing signal SMx%f)"%options.expSig)

if not options.includeVBF: options.splitSignal=False
if (not options.Bias) and options.includeVBF: sys.exit("Cannot use summed sideband (ie assume no mH dependance) for VBF (currently unsupported)")

type=options.bdtType

# create output folders
cardOutDir="mva-datacards-"+type
if not os.path.isdir(cardOutDir):
  os.makedirs(cardOutDir)
if options.makePlot:
  plotOutDir="mva-plots-"+type
  if not os.path.isdir(plotOutDir):
    os.makedirs(plotOutDir)

#if options.biasFile:
#	biasROOTFile = ROOT.TFile(options.biasFile)

genMasses     = [110,115,120,125,130,135,140,145,150]
#scalingErrors = [1.008,1.008,1.008,1.008,1.008,1.009,1.01,1.011] # Takes from P.Dauncey studies -> 7% window
#scalingErrors =  [1.007,1.007,1.006,1.008,1.007,1.008,1.009,1.01] # Takes from P.Dauncey studies -> 2% window
#scalingErrors = [1.013,1.013,1.012,1.012,1.014,1.015,1.016,1.016] # Takes from P.Dauncey studies -> 7% window (100-180)
#scalingErrors = [1.011,1.01,1.009,1.011,1.011,1.013,1.014,1.014] 	  # Takes from P.Dauncey studies -> 2% window (100-180)
#scalingErrors = [1.00815,1.01024,1.01076,1.01197,1.0099,1.009,1.00928,1.01054 ] 	  # Takes from P.Dauncey studies -> 2% window (100-180) / MIT Preselection
#scalingErrors = [ 1.008,1.008,1.008,1.008,1.01,1.010,1.011,1.012,1.012] # P.Dauncey 100-180 2% window /MIT preselction +BDT>-0.5
scalingErrors = [1.0136,1.0152,1.01425,1.01102,1.01283,1.01568,1.02157,1.02467,1.02466] # P.Dauncey 100-180, 2% window, MIT presel + BDT > 0.05 (Pow2 Fit)
#scalingErrors = [1.025,1.025,1.025,1.025,1.025,1.025,1.025,1.025,1.025] # FLAT 25%
#scalingErrors = [ 1.01185,1.01292,1.01378,1.01378,1.01594,1.01539,1.01814,1.02052,1.02257] # P.Dauncey 100-180, 2% window, MIT presel + BDT > 0.05 (Pol5 Fit)


evalMasses    = numpy.arange(110,150.5,0.5)
normG = ROOT.TGraph(len(genMasses))

# Fill the errors graph
can = ROOT.TCanvas()
for i,ne in enumerate(scalingErrors):
  normG.SetPoint(i,genMasses[i],ne)
normG.SetMarkerStyle(20)
normG.GetXaxis().SetTitle("mH")
normG.GetYaxis().SetTitle("(N+dN)/N")
#normG.Draw("ALP")
print "Check the Errors Look Sensible -> plot saved to normErrors_%s"%(options.tfileName)
can.SaveAs("normErrors_%s.pdf"%options.tfileName)

# Now we can write the cards
#tfileName = sys.argv[1]
tfile = ROOT.TFile(options.tfileName)
if options.singleMass>0: evalMasses=[float(options.singleMass)]
for m in evalMasses: writeCard(tfile,m,normG.Eval(m))
print "Done Writing Cards"

