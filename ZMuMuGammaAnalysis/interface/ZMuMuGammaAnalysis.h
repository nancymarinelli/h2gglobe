#ifndef __ZMUMUGAMMAANALYSIS__
#define __ZMUMUGAMMAANALYSIS__

#include "PhotonAnalysis/interface/StatAnalysis.h"

// ------------------------------------------------------------------------------------
class ZMuMuGammaAnalysis : public  StatAnalysis
{
 public:

  void Init(LoopAll& l);

  ZMuMuGammaAnalysis();
  virtual ~ZMuMuGammaAnalysis();
  
  struct TreeVariables 
  {
	  TreeVariables();
	  
	  unsigned int run, event, lumi;
	  TLorentzVector *leadMu, *subMu, *photon, *mm, *mmg;
	  int category;
	  double mass;
	  double weight;
	  int nvtx;
	  double idmva;
	  int ciclevel;
  };
  


private:
  virtual bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, 
			    float & evweight, int & category, int & diphoton_id,
			    bool & isCorrectVertex, float &kinematic_bdtout,
			    bool isSyst=false, 
			    float syst_shift=0., bool skipSelection=false,
			    BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0); 
  
  void fillPlots(int cat, float weight, LoopAll &l, TLorentzVector & leadMu, TLorentzVector & subMu, TLorentzVector & diMu, 
		 int iselPho, TLorentzVector & selPho, TLorentzVector & mmg);

  bool muonSelection (LoopAll& l, int iMu); 
  bool photonSelection (TLorentzVector& p4, TVector3 & sc);
  bool FSRselection ( LoopAll& l, int iMu1, int iMu2, int iPho, TClonesArray& p4Vec ); 
  TreeVariables treevars_;
  TLorentzVector&  getMumugP4() {return mumugMass_;}

 private:
  TLorentzVector mumugMass_;

 public:
  float muPtMin;
  int muTkLayers;
  bool muInnerHits;
  int muPixelHits;
  int muValidChambers;
  int muNmatches;
  float muNormChi2;
  float muD0Vtx;
  float muDZVtx;
  
  float phoPtMin;
  float dEtaMin;
  bool applyPhoPresel;

};

#endif
