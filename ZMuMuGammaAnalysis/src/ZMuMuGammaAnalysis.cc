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

// ----------------------------------------------------------------------------------------------------
void ZMuMuGammaAnalysis::Init(LoopAll& l)
{
    // initialize Hgg machinery, forcing 
    nEtaCategories = 2;
    nR9Categories = 2;
    nPtCategories = 2;

    includeVBF = false;
    includeVHhad = false;
    includeVHhadBtag = false;
    includeTTHhad = false;
    includeTTHlep = false;
    includeVHlep = false;
    includeVHlepPlusMet = false;
    nVHmetCategories = false;
    
    doVtxEffSmear = false;
    doTriggerEffSmear = false;
    doKFactorSmear = false;
    doPtSpinSmear = false;
    doPdfWeightSmear = false;
    doInterferenceSmear = false;
    doCosThetaDependentInterferenceSmear = false;
    
    skimOnDiphoN = false;
    nPhoMin = 1;

    massMin = 60;
    massMax = 120;
    nDataBins = 240;
    
    StatAnalysis::Init(l);
    
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
    TClonesArray mus_p4("mus_p4",l.mu_n);
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
    
    TLorentzVector * leadMu =  (TLorentzVector*)l.mu_glo_p4->At(ileadMu);
    TLorentzVector * subMu  =  (TLorentzVector*)l.mu_glo_p4->At(isubMu);
    
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
    TClonesArray phos_p4("phos_p4",l.pho_n);
    for(int ipho=0; ipho<l.pho_n; ++ipho) { 
    	    TLorentzVector * p4 = l.get_pho_p4( ipho, 0, &smeared_pho_energy[0]);
    	    // Fill TClonesArray with corrected 4-vectors
    	    new(phos_pt[sorted_phos.size()]) TLorentzVector( *p4 );
    	    bool passSelection = true;
    	    if( passSelection ) {
    		    sorted_phos.push_back(ipho);
    	    }
    	    delete p4;
    }
    if( sorted_phos.size() < 1 ) { return false; }
    std::sort(sorted_jets.begin(),sorted_jets.end(),
    	      ClonesSorter<TLorentzVector,double,std::greater<double> >(phos_p4,&TLorentzVector::Pt));
    
    int iselPho = sorted_phos[0];
    TLorentzVector * iselPho_p4 =  (TLorentzVector*)l.mu_glo_p4->At(iselPho);
    
    // Three body system
    // -------------------------------------------------------------------------------------------------
    TLorentzVector mmg = diMu + iselPho_p4;
    
    category = l.pho_r9[iselPho] < 0.94  + 2*(iselPho_p4.Pt()>30.) + 4*(l.pho_isEB[iselPho]);
    mass = mmg.M();
    
    evweight = weight * smeared_pho_weight[iselPho] * genLevWeight;
    if( ! isSyst ) {
    	    l.countersred[diPhoCounter_]++;
    }
    
    return (category >= 0 && mass>=massMin && mass<=massMax);
}
