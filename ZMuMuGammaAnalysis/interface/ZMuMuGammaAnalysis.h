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
	  TLorentzVector *leadMu, *subMu, *photon;
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
	  
  TreeVariables treevars_;

};

#endif
