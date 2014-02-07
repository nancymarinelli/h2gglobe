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
    PhotonAnalysis::Init(l);

    l.SetAllMVA();
    cout << "Weights file is: " << photonLevel2012IDMVA_EB.c_str() << endl;
    l.tmvaReaderID_2013_Barrel->BookMVA("AdaBoost",photonLevel2012IDMVA_EB.c_str());
    l.tmvaReaderID_2013_Endcap->BookMVA("AdaBoost",photonLevel2012IDMVA_EE.c_str());
    
}

//----------------------------------------------------------------------------------------------------
bool ZMuMuGammaAnalysis::SkimEvents(LoopAll& l, int jentry)
{
    if( run7TeV4Xanalysis ) { l.version=12; }
    else { l.version=13; }

    return true;
}
