#include "../interface/ZMuMuGammaAnalysis.h"

#include "PhotonReducedInfo.h"
#include "Sorters.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <algorithm>

#include "JetAnalysis/interface/JetHandler.h"
#include "CMGTools/External/interface/PileupJetIdentifier.h"

#define PADEBUG 0

using namespace std;

ZMuMuGammaAnalysis::ZMuMuGammaAnalysis(){}
ZMuMuGammaAnalysis::~ZMuMuGammaAnalysis(){}

ZMuMuGammaAnalysis::TreeVariables::TreeVariables() : leadMu(0), subMu(0), photon(0)
{}

// ----------------------------------------------------------------------------------------------------
void ZMuMuGammaAnalysis::Init(LoopAll& l)
{
    l.InitTrees("zmmgAnalysis");
    l.BookExternalTreeBranch( "run",         &treevars_.run, "zmmgAnalysis" );         
    l.BookExternalTreeBranch( "event",       &treevars_.event, "zmmgAnalysis" );         
    l.BookExternalTreeBranch( "lumi",        &treevars_.lumi, "zmmgAnalysis" );         
    l.BookExternalTreeBranch( "mass",        &treevars_.mass, "zmmgAnalysis" );         
    l.BookExternalTreeBranch( "weight",      &treevars_.weight, "zmmgAnalysis" );         
    l.BookExternalTreeBranch( "nvtx",        &treevars_.nvtx, "zmmgAnalysis" );         
    l.BookExternalTreeBranch( "leadMu",      &treevars_.leadMu, "zmmgAnalysis" );         
    l.BookExternalTreeBranch( "subMu",       &treevars_.subMu, "zmmgAnalysis" );         
    l.BookExternalTreeBranch( "photon",      &treevars_.photon, "zmmgAnalysis" );         
    l.BookExternalTreeBranch( "idmva",       &treevars_.idmva, "zmmgAnalysis" );         
    l.BookExternalTreeBranch( "ciclevel",    &treevars_.ciclevel, "zmmgAnalysis" );         

    l.rooContainer->BlindData(false);
    
    // initialize Hgg machinery, forcing     
    StatAnalysis::Init(l);

    // book photon ID MVA
    l.SetAllMVA();
    if( photonLevel2013IDMVA_EB != "" && photonLevel2013IDMVA_EE != "" ) {
	l.tmvaReaderID_2013_Barrel->BookMVA("AdaBoost",photonLevel2013IDMVA_EB.c_str());
	l.tmvaReaderID_2013_Endcap->BookMVA("AdaBoost",photonLevel2013IDMVA_EE.c_str());
    } else if( photonLevel2012IDMVA_EB != "" && photonLevel2012IDMVA_EE != "" ) {
    	l.tmvaReaderID_Single_Barrel->BookMVA("AdaBoost",photonLevel2012IDMVA_EB.c_str());
    	l.tmvaReaderID_Single_Endcap->BookMVA("AdaBoost",photonLevel2012IDMVA_EE.c_str());
	assert( bdtTrainingType == "Moriond2013" ); 
    } else if (photonLevel2013_7TeV_IDMVA_EB != "" && photonLevel2013_7TeV_IDMVA_EE != "" ) {
    	l.tmvaReaderID_2013_7TeV_MIT_Barrel->BookMVA("AdaBoost",photonLevel2013_7TeV_IDMVA_EB.c_str());
    	l.tmvaReaderID_2013_7TeV_MIT_Endcap->BookMVA("AdaBoost",photonLevel2013_7TeV_IDMVA_EE.c_str());
    } else { 
    	assert( run7TeV4Xanalysis );
    }

    // replace trigger selection
    triggerSelections.clear();
    triggerSelections.push_back(TriggerSelection(1,-1));
    triggerSelections.back().addpath("HLT_Mu17_TkMu8");
}

// ----------------------------------------------------------------------------------------------------
bool ZMuMuGammaAnalysis::AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4,
        float & mass, float & evweight, int & category, int & diphoton_id, bool & isCorrectVertex,
        float &kinematic_bdtout,
        bool isSyst,
        float syst_shift, bool skipSelection,
        BaseGenLevelSmearer *genSys, BaseSmearer *phoSys, BaseDiPhotonSmearer * diPhoSys)
{
    assert( ! skipSelection );
    
    int cur_type = l.itype[l.current];
    float sampleweight = l.sampleContainer[l.current_sample_index].weight();
    
    // do gen-level dependent first (e.g. k-factor); only for signal
    genLevWeight=1.;
    
    // event selection
    int muVtx=-1;
    int mu_ind=-1;
    int elVtx=-1;
    int el_ind=-1;
    
    int leadpho_ind=-1;
    int subleadpho_ind=-1;
    
    // Muon selection
    // -------------------------------------------------------------------------------------------------
    std::vector<int> sorted_mus;
    for(int imu=0; imu<l.mu_glo_n; ++imu) { 
    	    TLorentzVector * p4 = (TLorentzVector*)l.mu_glo_p4->At(imu);
    	    bool passSelection = true;
    	    if( passSelection ) {
    		    sorted_mus.push_back(imu);
    	    }
    }
    if( sorted_mus.size() < 2 ) { return false; }
    std::sort(sorted_mus.begin(),sorted_mus.end(),
    	      ClonesSorter<TLorentzVector,double,std::greater<double> >(l.mu_glo_p4,&TLorentzVector::Pt));
    
    int ileadMu = sorted_mus[0];
    int isubMu  = sorted_mus[1];
    
    TLorentzVector & leadMu =  *( (TLorentzVector*)l.mu_glo_p4->At(ileadMu));
    TLorentzVector & subMu  =  *( (TLorentzVector*)l.mu_glo_p4->At(isubMu) );
    
    TLorentzVector diMu = leadMu + subMu;
    
    // Photon selection
    // -------------------------------------------------------------------------------------------------
    
    // First apply corrections and smearing on the single photons
    smeared_pho_energy.clear(); smeared_pho_energy.resize(l.pho_n,0.);
    smeared_pho_r9.clear();     smeared_pho_r9.resize(l.pho_n,0.);
    smeared_pho_weight.clear(); smeared_pho_weight.resize(l.pho_n,1.);
    applySinglePhotonSmearings(smeared_pho_energy, smeared_pho_r9, smeared_pho_weight, cur_type, l, energyCorrected, energyCorrectedError,
    			       phoSys, syst_shift);
        
    std::vector<int> sorted_phos;
    TClonesArray phos_p4(TLorentzVector::Class(),l.pho_n);
    for(int ipho=0; ipho<l.pho_n; ++ipho) { 
    	    TLorentzVector p4 = l.get_pho_p4( ipho, 0, &smeared_pho_energy[0]);
    	    // Fill TClonesArray with corrected 4-vectors
    	    new(phos_p4[ipho]) TLorentzVector(p4);
    	    bool passSelection = true;
    	    if( passSelection ) {
    		    sorted_phos.push_back(ipho);
    	    }
    }
    if( sorted_phos.size() < 1 ) { return false; }
    std::sort(sorted_phos.begin(),sorted_phos.end(),
    	      ClonesSorter<TLorentzVector,double,std::greater<double> >(&phos_p4,&TLorentzVector::Pt));
    
    int iselPho = sorted_phos[0];
    TLorentzVector & selPho =  *((TLorentzVector*)phos_p4.At(iselPho));
    
    // Three body system
    // -------------------------------------------------------------------------------------------------
    TLorentzVector mmg = diMu + selPho;
    
    int etacat = (l.pho_isEB[iselPho]);
    int ptcat  = (selPho.Pt()>30.);
    int r9cat  = (l.pho_r9[iselPho] < 0.94);
    category = r9cat  + 2*etacat + 4*ptcat;
    mass = mmg.M();
    
    evweight = weight * smeared_pho_weight[iselPho] * genLevWeight;
    if( ! isSyst ) {
    	    l.countersred[diPhoCounter_]++;
    }
    
    // fill control plots
    fillPlots(0,evweight,l,leadMu,subMu,diMu,iselPho,selPho,mmg);
    fillPlots(1+etacat,evweight,l,leadMu,subMu,diMu,iselPho,selPho,mmg);
    fillPlots(3+2*etacat+ptcat,evweight,l,leadMu,subMu,diMu,iselPho,selPho,mmg);
    fillPlots(7+category,evweight,l,leadMu,subMu,diMu,iselPho,selPho,mmg);
    
    return (category >= 0 && mass>=massMin && mass<=massMax);
}

void ZMuMuGammaAnalysis::fillPlots(int cat, float evweight, 
				   LoopAll &l, TLorentzVector & leadMu, TLorentzVector & subMu, TLorentzVector & diMu, 
				   int iselPho, TLorentzVector & selPho, TLorentzVector & mmg)
{
	TVector3 & selSc = *((TVector3*)l.sc_xyz->At(l.pho_scind[iselPho]));
	double selScE =  l.sc_raw[l.pho_scind[iselPho]];
        l.FillHist("pt"       ,cat,evweight,mmg.Pt());
        l.FillHist("eta"      ,cat,evweight,mmg.Eta());
        l.FillHist("phi"      ,cat,evweight,mmg.Phi());
        l.FillHist("mass"     ,cat,evweight,mmg.M());
        l.FillHist("nvtx"     ,cat,evweight,l.vtx_std_n);
        l.FillHist("pho_n"    ,cat,evweight,l.pho_n);
        l.FillHist("pho_e"    ,cat,evweight,selPho.E());
        l.FillHist("pho_pt"   ,cat,evweight,selPho.Pt());
        l.FillHist("pho_eta"  ,cat,evweight,selPho.Eta());
        l.FillHist("pho_phi"  ,cat,evweight,selPho.Phi());
        l.FillHist("pho_sce"  ,cat,evweight,selScE);
        l.FillHist("pho_sceta",cat,evweight,selSc.Eta());
        l.FillHist("pho_r9"   ,cat,evweight,l.pho_r9[iselPho]);
        l.FillHist("mu1_pt"   ,cat,evweight,leadMu.Pt());
        l.FillHist("mu1_eta"  ,cat,evweight,leadMu.Eta());
        l.FillHist("mu1_phi"  ,cat,evweight,leadMu.Phi());
        l.FillHist("mu2_pt"   ,cat,evweight,subMu.Pt());
        l.FillHist("mu2_eta"  ,cat,evweight,subMu.Eta());
        l.FillHist("mu2_phi"  ,cat,evweight,subMu.Phi());
        l.FillHist("mumu_pt"  ,cat,evweight,diMu.Pt());
        l.FillHist("mumu_eta" ,cat,evweight,diMu.Eta());
        l.FillHist("mumu_phi" ,cat,evweight,diMu.Phi());
        l.FillHist("mumu_mass",cat,evweight,diMu.M());
	
	if( cat == 0 ) { 
                treevars_.run       = l.run;
                treevars_.event     = l.event;
                treevars_.lumi      = l.lumis;
                treevars_.mass      = mmg.M();
                treevars_.weight    = evweight;
                treevars_.nvtx      = l.vtx_std_n;
                *(treevars_.leadMu) = leadMu;
                *(treevars_.subMu)  = subMu;
                *(treevars_.photon) = selPho;
		treevars_.idmva     = l.photonIDMVA(iselPho,0,selPho,bdtTrainingType.c_str());
		std::vector<std::vector<bool> > ph_passcut;
		treevars_.ciclevel  = l.PhotonCiCPFSelectionLevel(iselPho, 0, ph_passcut, 4, 0, &smeared_pho_energy[0]);
	}
}
