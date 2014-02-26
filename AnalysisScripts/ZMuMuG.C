{
  TFile* f1=new 
TFile("/afs/cern.ch/work/n/nancy/private/Higgs_paper/zmumuMassStudies/CMSSW_6_1_2/src/h2gglobe_zmmg/AnalysisScripts/fromPasquale/histograms_CMS-HGG.root");
  f1->cd();
  
  // Photons plots
  string var[8]={"n","eta","phi","pt","r9","e","sce","sceta"};
  for(int iVar=0;iVar<8;iVar++){
    for(int j=0;j<15;j++){
      TCanvas *c1 = new TCanvas("c1", "c1",430, 10, 600,600);    
      TString h_mcname;
      TString h_dataname;
      TString fileName;
      h_mcname.Form("pho_%s_cat%d_dymm_m90",var[iVar].c_str(),j);
      h_dataname.Form("pho_%s_cat%d_Data",var[iVar].c_str(),j);
      fileName.Form("pho_%s_cat%d",var[iVar].c_str(),j);
      cout << h_mcname << " " << h_dataname << endl;
      TH1F *h_mc = (TH1F*)gDirectory->Get(h_mcname);
      TH1F *h_data = (TH1F*)gDirectory->Get(h_dataname);
      if ( h_mc && h_data) {
	float mcEntries = h_mc->Integral();
	float dataEntries = h_data->Integral();
	float scaleF= dataEntries/mcEntries;
	//	h_mc->Scale(scaleF);
	h_mc->SetLineWidth(2);
	h_mc->SetFillColor(3);
        h_mc->SetMinimum(0.);
	h_data->SetMarkerStyle(20);
	h_data->SetMarkerSize(0.8);
	h_mc->Draw("hist"); 
	h_data->Draw("same");
	//
        int nBins=h_mc->GetNbinsX();
        float xMin=h_mc->GetBinLowEdge(1);
	float xMax=h_mc->GetBinLowEdge(nBins)+h_mc->GetBinWidth(nBins);
	TH1F* h_ratio= new TH1F("ratio"," ",nBins,xMin,xMax);
	h_ratio->Divide(h_data,h_mc);
	h_ratio->SetStats(0);
        for ( int i=1; i<=ratio->GetNbinsX(); i++ ) {
	  float num=h_data->GetBinContent(i);
	  float den=h_mc->GetBinContent(i);
	  float dNum=h_data->GetBinError(i);
	  float dDen=h_mc->GetBinError(i);
	  float erro=0;
	  if ( num!=0 && den!=0) {
	    erro= ((1./den)*(1./den)*dNum*dNum) + ((num*num)/(den*den*den*den) * (dDen*dDen));
	    erro=sqrt(erro);
	  }
	  h_ratio->SetBinError(i, erro);
	}
	h_ratio->SetLineColor(1);
	h_ratio->SetLineWidth(2);
	h_ratio->SetMinimum(0.);
	h_ratio->SetMaximum(2.);
	
	TCanvas *c2 = new TCanvas("c2", "c2",430, 10, 600,600);    
	h_ratio->Draw("e");
	TLine *l = new TLine(xMin,1.,xMax,1.);
	l->Draw();
	
	//
      }
      c1->SaveAs("fromPasquale/plots/"+fileName+".png");
      c2->SaveAs("fromPasquale/plots/"+fileName+"_ratio.png");
    }
  }
  
  // Muon Plots 
  string varMu[3]={"eta","phi","pt"};
  for (int iObj=1; iObj<3; iObj++) {
    for(int iVar=0;iVar<3;iVar++){
      for(int j=0;j<15;j++){
	TCanvas *c1 = new TCanvas("c1", "c1",430, 10, 600,600);    
	TString h_mcname;
	TString h_dataname;
	TString fileName;
	h_mcname.Form("mu%d_%s_cat%d_dymm_m90",iObj,varMu[iVar].c_str(),j);
	h_dataname.Form("mu%d_%s_cat%d_Data",iObj,varMu[iVar].c_str(),j);
	fileName.Form("mu%d_%s_cat%d",iObj,varMu[iVar].c_str(),j);
	cout << h_mcname << " " << h_dataname << endl;
	TH1F *h_mc = (TH1F*)gDirectory->Get(h_mcname);
	TH1F *h_data = (TH1F*)gDirectory->Get(h_dataname);
	if ( h_mc && h_data) {
	  float mcEntries = h_mc->Integral();
	  float dataEntries = h_data->Integral();
	  float scaleF= dataEntries/mcEntries;
	  // h_mc->Scale(scaleF);
	  h_mc->SetLineWidth(2);
	  h_mc->SetFillColor(3);
	  h_mc->SetMinimum(0.);
	  h_data->SetMarkerStyle(20);
	  h_data->SetMarkerSize(0.8);
	  h_mc->Draw("hist"); 
	  h_data->Draw("same");
	//
        int nBins=h_mc->GetNbinsX();
        float xMin=h_mc->GetBinLowEdge(1);
	float xMax=h_mc->GetBinLowEdge(nBins)+h_mc->GetBinWidth(nBins);
	TH1F* h_ratio= new TH1F("ratio"," ",nBins,xMin,xMax);
	h_ratio->Divide(h_data,h_mc);
	h_ratio->SetStats(0);
        for ( int i=1; i<=ratio->GetNbinsX(); i++ ) {
	  float num=h_data->GetBinContent(i);
	  float den=h_mc->GetBinContent(i);
	  float dNum=h_data->GetBinError(i);
	  float dDen=h_mc->GetBinError(i);
	  float erro=0;
	  if ( num!=0 && den!=0) {
	    erro= ((1./den)*(1./den)*dNum*dNum) + ((num*num)/(den*den*den*den) * (dDen*dDen));
	    erro=sqrt(erro);
	  }
	  h_ratio->SetBinError(i, erro);
	}
	h_ratio->SetLineColor(1);
	h_ratio->SetLineWidth(2);
	h_ratio->SetMinimum(0.);
	h_ratio->SetMaximum(2.);
	
	TCanvas *c2 = new TCanvas("c2", "c2",430, 10, 600,600);    
	h_ratio->Draw("e");
	TLine *l = new TLine(xMin,1.,xMax,1.);
	l->Draw();
	
	//

	}
	c1->SaveAs("fromPasquale/plots/"+fileName+".png");
	c2->SaveAs("fromPasquale/plots/"+fileName+"_ratio.png");
      }
    }
  }

  // mumu and mumug system plots  
  // Muon Plots
  string obj[2]={"mumu","mmg"};
  string varM[4]={"eta","phi","pt","mass"};
  for (int iObj=0; iObj<2; iObj++) {
    for(int iVar=0;iVar<4;iVar++){
      for(int j=0;j<15;j++){
	TCanvas *c1 = new TCanvas("c1", "c1",430, 10, 600,600);    
	TString h_mcname;
	TString h_dataname;
	TString fileName;
	h_mcname.Form("%s_%s_cat%d_dymm_m90",obj[iObj].c_str(),varM[iVar].c_str(),j);
	h_dataname.Form("%s_%s_cat%d_Data",obj[iObj].c_str(),varM[iVar].c_str(),j);
	fileName.Form("%s_%s_cat%d",obj[iObj].c_str(),varM[iVar].c_str(),j);
	cout << h_mcname << " " << h_dataname << endl;
	TH1F *h_mc = (TH1F*)gDirectory->Get(h_mcname);
	TH1F *h_data = (TH1F*)gDirectory->Get(h_dataname);
	if ( h_mc && h_data) {
	  float mcEntries = h_mc->Integral();
	  float dataEntries = h_data->Integral();
	  float scaleF= dataEntries/mcEntries;
	  //h_mc->Scale(scaleF);
	  h_mc->SetLineWidth(2);
	  h_mc->SetFillColor(3);
	  h_mc->SetMinimum(0.);
	  h_data->SetMarkerStyle(20);
	  h_data->SetMarkerSize(0.8);
	  h_mc->Draw("hist"); 
	  h_data->Draw("same");
	//
        int nBins=h_mc->GetNbinsX();
        float xMin=h_mc->GetBinLowEdge(1);
	float xMax=h_mc->GetBinLowEdge(nBins)+h_mc->GetBinWidth(nBins);
	TH1F* h_ratio= new TH1F("ratio"," ",nBins,xMin,xMax);
	h_ratio->Divide(h_data,h_mc);
	h_ratio->SetStats(0);
        for ( int i=1; i<=ratio->GetNbinsX(); i++ ) {
	  float num=h_data->GetBinContent(i);
	  float den=h_mc->GetBinContent(i);
	  float dNum=h_data->GetBinError(i);
	  float dDen=h_mc->GetBinError(i);
	  float erro=0;
	  if ( num!=0 && den!=0) {
	    erro= ((1./den)*(1./den)*dNum*dNum) + ((num*num)/(den*den*den*den) * (dDen*dDen));
	    erro=sqrt(erro);
	  }
	  h_ratio->SetBinError(i, erro);
	}
	h_ratio->SetLineColor(1);
	h_ratio->SetLineWidth(2);
	h_ratio->SetMinimum(0.);
	h_ratio->SetMaximum(2.);
	
	TCanvas *c2 = new TCanvas("c2", "c2",430, 10, 600,600);    
	h_ratio->Draw("e");
	TLine *l = new TLine(xMin,1.,xMax,1.);
	l->Draw();
	
	//

	}
	c1->SaveAs("fromPasquale/plots/"+fileName+".png");
	c2->SaveAs("fromPasquale/plots/"+fileName+"_ratio.png");
      }
    }
  }





  
}
