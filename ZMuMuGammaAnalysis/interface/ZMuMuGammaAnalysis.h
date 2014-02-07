#ifndef __ZMUMUGAMMAANALYSIS__
#define __ZMUMUGAMMAANALYSIS__

#include "PhotonAnalysis/interface/StatAnalysis.h"

// ------------------------------------------------------------------------------------
class ZMuMuGammaAnalysis : public StatAnalysis 
{
 public:

  void Init(LoopAll& l);
  bool SkimEvents(LoopAll&, int);

  ZMuMuGammaAnalysis();
  virtual ~ZMuMuGammaAnalysis();
  
private:
  virtual bool AnalyseEvent(LoopAll& l, Int_t jentry, float weight, TLorentzVector & gP4, float & mass, float & evweight, int & category, int & diphoton_id,
			    bool & isCorrectVertex, float &kinematic_bdtout,
			    bool isSyst=false, 
			    float syst_shift=0., bool skipSelection=false,
			    BaseGenLevelSmearer *genSys=0, BaseSmearer *phoSys=0, BaseDiPhotonSmearer * diPhoSys=0); 

};

#endif
