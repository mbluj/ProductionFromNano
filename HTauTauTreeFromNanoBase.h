/*****************************
* M. Bluj, NCBJ, Poland
* 6/12/2017
*****************************/

#ifndef HTauTauTreeFromNanoBase_h
#define HTauTauTreeFromNanoBase_h

/*Include header of base class providing interface to the NanoAOD Events tree 
  as produced by nanoEventsChain->MakeClass("NanoEventsSkeleton") */
#undef NanoEventsSkeleton_cxx //undefine to protect againist problems with multile implementation
#include "NanoEventsSkeleton.h"
//#include <TROOT.h>
//#include <TChain.h>
//#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TMatrixD.h>
#include "Math/PtEtaPhiE4D.h"
#include "Math/PtEtaPhiM4D.h"
#include "Math/LorentzVector.h"

#include "HTTEvent.h"
#include <vector>
#include <iostream>

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfit.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/svFitHistogramAdapter.h"

#include "HTT-utilities/RecoilCorrections/interface/RecoilCorrector.h"

#include "CondFormats/JetMETObjects/interface/JetCorrectorParameters.h"
#include "CondFormats/JetMETObjects/interface/JetCorrectionUncertainty.h"

#include "DataFormats/Provenance/interface/LuminosityBlockRange.h"

class HTauTauTreeFromNanoBase : public NanoEventsSkeleton {
public :

  /////////////////////////////////////////////////
  /// Lorentz vector
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;
  typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PolarLorentzVector;
  /// Trigger data struct
  struct TriggerData {
    std::string path_name;
    unsigned int leg1Id, leg2Id;//0-undefined, 11-electron, 13-muon, 15-tau
    int leg1BitMask, leg2BitMask;//definition depends on Id, cf. PhysicsTools/NanoAOD/python/triggerObjects_cff.py
    float leg1Pt, leg2Pt, leg1L1Pt, leg2L1Pt;
    float leg1Eta, leg2Eta;
  };

  virtual void initHTTTree(const TTree *tree, std::string prefix="HTT");
  void initJecUnc(std::string correctionFile);

  void fillEvent();
  virtual bool buildPairs();
  virtual void fillPairs(unsigned int bestPairIndex);
  virtual void fillJets(unsigned int bestPairIndex);
  virtual void fillLeptons();
  virtual void fillGenLeptons();
  void applyMetRecoilCorrections();
  virtual bool thirdLeptonVeto(unsigned int signalLeg1Index, unsigned int signalLeg2Index, int leptonPdg, double dRmin=-1);
  virtual bool extraMuonVeto(unsigned int signalLeg1Index, unsigned int signalLeg2Index, double dRmin=-1);
  virtual bool extraElectronVeto(unsigned int signalLeg1Index, unsigned int signalLeg2Index, double dRmin=-1);
  bool muonSelection(unsigned int index);
  bool electronSelection(unsigned int index);
  virtual bool pairSelection(unsigned int index);
  virtual unsigned int bestPair(std::vector<unsigned int> &pairIndexes);
  void computeSvFit(HTTPair &aPair, HTTAnalysis::sysEffects type=HTTAnalysis::NOMINAL);
  TLorentzVector runSVFitAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons,
			      const TVector2 &aMET, const TMatrixD &covMET);
  bool jetSelection(unsigned int index, unsigned int bestPairIndex);
  int getMCMatching(unsigned int index, std::string colType="");
  int getTriggerMatchning(unsigned int index, bool checkBit=false, std::string colType="");
  int getMetFilterBits();
  double getPtReweight(const TLorentzVector &genBosonP4, bool doSUSY=false);
  bool isGoodToMatch(unsigned int ind);
  TLorentzVector getGenComponentP4(std::vector<unsigned int> &indexes, unsigned int iAbsCharge);
  bool eventInJson();

  template<typename T> T getBranchValue(const char *branchAddress, unsigned int index);
  Double_t getProperty(std::string name, unsigned int index, std::string colType="");
  std::vector<Double_t> getProperties(const std::vector<std::string> & propertiesList, unsigned int index, std::string colType="");
  void writePropertiesHeader(const std::vector<std::string> & propertiesList);
  void writeTriggersHeader(const std::vector<TriggerData> &triggerBits);
  double getJecUnc(unsigned int index, std::string name="Total", bool up=true);
  static bool compareLeptons(const HTTParticle& i, const HTTParticle& j);
  static bool comparePairs(const HTTPair& i, const HTTPair& j);
  //int isGenPartDaughterPdgId(int index, unsigned int aPdgId);
  //bool isGenPartDaughterIdx(int index, int mother);
  bool getDirectDaughterIndexes(std::vector<unsigned int> & indexes, unsigned int motherIndex, bool ignoreNeutrinos=true);
  unsigned int findFinalCopy(unsigned int index);
  unsigned int findFirstCopy(unsigned int index);
  int genTauDecayMode(std::vector<unsigned int> &daughterIndexes);
  bool findBosonP4(TLorentzVector &bosonP4, TLorentzVector &visBosonP4);
  bool findTopP4(TLorentzVector &topP4, TLorentzVector &antiTopP4);

  std::vector<HTTPair> httPairCollection, httPairs_;
  std::vector<HTTParticle> httJetCollection;
  std::vector<HTTParticle> httLeptonCollection;
  std::vector<HTTParticle> httGenLeptonCollection;
  std::vector<TriggerData> triggerBits_;
  std::vector<std::string> filterBits_;
  TTree *httTree;
  TFile *httFile;
  HTTEvent *httEvent;
  TH1F* hStats;
  TH2F* zptmass_histo, *zptmass_histo_SUSY;
  
  unsigned int bestPairIndex_;

  ClassicSVfit *svFitAlgo_;
  RecoilCorrector* recoilCorrector_;
  TFile* zPtReweightFile, *zPtReweightSUSYFile;
  TLorentzVector p4SVFit, p4Leg1SVFit, p4Leg2SVFit;   

  std::vector<edm::LuminosityBlockRange> jsonVector;

  bool firstWarningOccurence_; // used to print warnings only at first occurnece in the event loop

  std::vector<std::string> leptonPropertiesList, genLeptonPropertiesList, jecUncertList;
  std::vector<JetCorrectionUncertainty*> jecUncerts;

  HTauTauTreeFromNanoBase(TTree *tree=0, bool doSvFit=false, bool correctRecoil=false, std::vector<std::string> lumis = std::vector<std::string>(), string prefix="HTT");
  virtual ~HTauTauTreeFromNanoBase();
  virtual Int_t    Cut(Long64_t entry);
  virtual void     Loop();
};

#endif

