#define HTauTauTreeFromNanoBase_cxx

#include "HTauTauTreeFromNanoBase.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TLeaf.h>

#include <iostream>
#include <fstream>
#include <algorithm>

HTauTauTreeFromNanoBase::HTauTauTreeFromNanoBase(TTree *tree, bool doSvFit, bool correctRecoil, std::vector<std::string> lumis, std::string prefix) : NanoEventsSkeleton(tree)
{

  ///Init HTT ntuple
  initHTTTree(tree, prefix);

  ///Parse lumis to be processed
  for(unsigned int iL=0; iL<lumis.size(); ++iL){
    std::size_t pos = lumis[iL].find("-");
    const char *block1 = (lumis[iL].substr(0,pos)).c_str();
    const char *block2 = (lumis[iL].substr(pos+1)).c_str();
    char *pEnd = NULL;
    long int r1, l1, r2, l2;
    r1 = strtol(block1,&pEnd,10);
    l1 = strtol(pEnd+1,NULL,10);
    pEnd = NULL;
    r2 = strtol(block2,&pEnd,10);
    l2 = strtol(pEnd+1,NULL,10);
    edm::LuminosityBlockRange lumiRange(edm::LuminosityBlockID(r1,l1),
					edm::LuminosityBlockID(r2,l2));
    jsonVector.push_back(lumiRange);
  }
  std::cout<<"[HTauTauTreeFromNanoBase]: Size of jsonVec: "<<jsonVector.size()<<std::endl;
  //for(unsigned int iL=0; iL<jsonVector.size(); ++iL){
  //  std::cout<<"\t lumi range: "<<jsonVector[iL]<<std::endl;
  //}

  ///Initialization of SvFit
  if(doSvFit){
    std::cout<<"[HTauTauTreeFromNanoBase]: Run w/ SVFit"<<std::endl;
    unsigned int verbosity = 0;//Set the debug level to 3 for testing
    svFitAlgo_ = new ClassicSVfit(verbosity);
    //svFitAlgo_->setMaxObjFunctionCalls(100000); // CV: default is 100000 evaluations of integrand per event
    svFitAlgo_->setHistogramAdapter(new classic_svFit::DiTauSystemHistogramAdapter());//needed?
    //svFitAlgo_->setLikelihoodFileName("testClassicSVfit.root");//needed?
    svFitAlgo_->setDiTauMassConstraint(-1);//argument>0 constraints di-tau mass to its value
  } else {
    std::cout<<"[HTauTauTreeFromNanoBase]: Run w/o SVFit"<<std::endl;
    svFitAlgo_=nullptr;
  }

  ///Initialization of RecoilCorrector
  if(correctRecoil){
    std::cout<<"[HTauTauTreeFromNanoBase]: Apply MET recoil corrections"<<std::endl;
    //std::string correctionFile = std::string(getenv("CMSSW_BASE"))+"/src/";
    //correctionFile += "HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root";
    std::string correctionFile = "HTT-utilities/RecoilCorrections/data/TypeI-PFMet_Run2016BtoH.root";
    recoilCorrector_= new RecoilCorrector(correctionFile);
  } else{
    std::cout<<"[HTauTauTreeFromNanoBase]: Do not apply MET recoil corrections"<<std::endl;
    recoilCorrector_=nullptr;
  }

  ///Get files with weights
  zPtReweightFile = new TFile("zpt_weights_2016_BtoH.root");  
  if(!zPtReweightFile) std::cout<<"Z pt reweight file zpt_weights.root is missing."<<std::endl;
  zptmass_histo = (TH2F*)zPtReweightFile->Get("zptmass_histo");

  zPtReweightSUSYFile = new TFile("zpt_weights_summer2016.root");  
  if(!zPtReweightSUSYFile) std::cout<<"SUSY Z pt reweight file zpt_weights.root is missing."<<std::endl;
  zptmass_histo_SUSY = (TH2F*)zPtReweightSUSYFile->Get("zptmass_histo");

  ///Instantiate JEC uncertainty sources
  ///https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECDataMC
  initJecUnc("Summer16_23Sep2016V4_MC_UncertaintySources_AK4PFchs.txt");//need to data file to process

  firstWarningOccurence_=true;
}

HTauTauTreeFromNanoBase::~HTauTauTreeFromNanoBase()
{

  if(httFile){
    httFile->Write();
    delete httFile;
  }
  if(svFitAlgo_) delete svFitAlgo_;
  if(recoilCorrector_) delete recoilCorrector_;
  if(zPtReweightFile) delete zPtReweightFile;
  if(zPtReweightSUSYFile) delete zPtReweightSUSYFile;
}

/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::initHTTTree(const TTree *tree, std::string prefix){

  if(prefix=="") prefix="HTT";
  prefix += "_";
  std::string filePath(tree->GetCurrentFile()->GetName());
  size_t location = filePath.find_last_of("/");
  if(location==std::string::npos) location = 0;
  else location+=1;
  std::string fileName = prefix+filePath.substr(location,filePath.size());
  httFile = new TFile(fileName.c_str(),"RECREATE");
  httEvent = new HTTEvent();
  httTree = new TTree("HTauTauTree","");
  httTree->SetDirectory(httFile);
  TBranch *eventBranch = httTree->Branch("HTTEvent.",&httEvent);
  TBranch *pairBranch = httTree->Branch("HTTPairCollection",&httPairCollection);
  TBranch *jetBranch = httTree->Branch("HTTJetCollection",&httJetCollection);
  TBranch *leptonBranch = httTree->Branch("HTTLeptonCollection",&httLeptonCollection);
  TBranch *genLeptonBranch = httTree->Branch("HTTGenLeptonCollection",&httGenLeptonCollection);
  hStats = new TH1F("hStats","Bookkeeping histogram",11,-0.5,10.5);
  hStats->SetDirectory(httFile);

  leptonPropertiesList.push_back("pdgId");
  leptonPropertiesList.push_back("charge");
  leptonPropertiesList.push_back("Tau_decayMode");
  //leptonPropertiesList.push_back("discriminator");
  //leptonPropertiesList.push_back("daughters_muonID");
  leptonPropertiesList.push_back("Muon_mediumId");
  leptonPropertiesList.push_back("Electron_mvaSpring16GP_WP80");
  leptonPropertiesList.push_back("Electron_mvaSpring16GP_WP90");
  leptonPropertiesList.push_back("Electron_lostHits");
  leptonPropertiesList.push_back("Electron_convVeto");
  //leptonPropertiesList.push_back("daughters_typeOfMuon");
  leptonPropertiesList.push_back("Tau_rawIso");
  leptonPropertiesList.push_back("Tau_photonsOutsideSignalCone");
  leptonPropertiesList.push_back("Tau_rawMVAoldDM");
  leptonPropertiesList.push_back("Tau_rawAntiEleCat");
  leptonPropertiesList.push_back("dxy");
  leptonPropertiesList.push_back("dz");
  leptonPropertiesList.push_back("sip3d");
  leptonPropertiesList.push_back("Tau_idDecayMode");
  leptonPropertiesList.push_back("Tau_idMVAoldDM");//bits: 1-VL, 2-L, 4-M, 8-T, 16-VT, 32-VVT 
  leptonPropertiesList.push_back("Tau_idAntiEle");//bits: 1-VL, 2-L, 4-M, 8-T, 16-VT 
  leptonPropertiesList.push_back("Tau_idAntiMu");//bits: 1-L, 2-T
  leptonPropertiesList.push_back("pfRelIso03_all");//R=0.4 used for mu?
  leptonPropertiesList.push_back("Muon_pfRelIso04_all");//R=0.4 used for mu?
  //leptonPropertiesList.push_back("Tau_chargedIso");
  leptonPropertiesList.push_back("Tau_leadTkPtOverTauPt");
  //leptonPropertiesList.push_back("daughters_isGoodTriggerType");
  //leptonPropertiesList.push_back("daughters_FilterFired");
  leptonPropertiesList.push_back("isGoodTriggerType");
  leptonPropertiesList.push_back("FilterFired");
  //leptonPropertiesList.push_back("daughters_L3FilterFired");
  //leptonPropertiesList.push_back("daughters_L3FilterFiredLast");
  leptonPropertiesList.push_back("mc_match");

  leptonPropertiesList.push_back("Jet_rawFactor");//was rawPt is rawFactor=rawPt/corrPt
  leptonPropertiesList.push_back("Jet_area");
  leptonPropertiesList.push_back("Jet_puId");
  leptonPropertiesList.push_back("Jet_partonFlavour");
  //leptonPropertiesList.push_back("bDiscriminator");?
  leptonPropertiesList.push_back("Jet_btagCSVV2");
  leptonPropertiesList.push_back("Jet_jetId");//bit1=L,bit2=T
  ////////////////////////////////////////////////////////////
  ///Gen Lepton properties MUST be synchronized with lepton properties
  ///since the branches name are not uniform, we need a second names vector.
  genLeptonPropertiesList.push_back("genpart_pdg");//needed?
  genLeptonPropertiesList.push_back("genpart_TauGenDetailedDecayMode");//needed?

  ////////////////////////////////////////////////////////////
  ///Trigger bits to check
  ///FIXME: is there a nicer way to define trigger list, e.g. a cfg file?
  TriggerData aTrgData;
  // Single muon triggers
  ///2nd bit (1<<1) for IsoMuon (not sure if correctly encoded in NanoAOD for 80X)
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_IsoMu22";
  triggerBits_.back().leg1Id=13;
  triggerBits_.back().leg1BitMask=(1<<1);
  triggerBits_.back().leg1Pt=22;
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_IsoTkMu22";
  triggerBits_.back().leg1Id=13;
  triggerBits_.back().leg1BitMask=(1<<3);
  triggerBits_.back().leg1Pt=22;
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_IsoTkMu22_eta2p1";
  triggerBits_.back().leg1Id=13;
  triggerBits_.back().leg1BitMask=(1<<3);
  triggerBits_.back().leg1Pt=22;
  triggerBits_.back().leg1Eta=2.1;
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_IsoMu22_eta2p1";
  triggerBits_.back().leg1Id=13;
  triggerBits_.back().leg1BitMask=(1<<1);
  triggerBits_.back().leg1Pt=22;
  triggerBits_.back().leg1Eta=2.1;
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_IsoMu24";
  triggerBits_.back().leg1Id=13;
  triggerBits_.back().leg1BitMask=(1<<1);
  triggerBits_.back().leg1Pt=24;
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_IsoTkMu24";
  triggerBits_.back().leg1Id=13;
  triggerBits_.back().leg1BitMask=(1<<3);
  triggerBits_.back().leg1Pt=24;
  // Mu tauh triggers
  ///2nd mu bit (1<<1) for IsoMuon (not sure if correctly encoded in NanoAOD for 80X)
  ///3rd mu bit (1<<2) for mu-tau overlap (should be OK)
  ///6th tau bit(1<<5) for mu-tau overlap  (is OK for all 80X triggers?)
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_IsoMu19_eta2p1_LooseIsoPFTau20";
  triggerBits_.back().leg1Id=13;
  triggerBits_.back().leg1BitMask=(1<<1)+(1<<2);
  triggerBits_.back().leg1Pt=19;
  triggerBits_.back().leg1Eta=2.1;
  triggerBits_.back().leg2Id=15;
  triggerBits_.back().leg2BitMask=(1<<5);
  triggerBits_.back().leg2Pt=20;
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_IsoMu19_eta2p1_LooseIsoPFTau20_SingleL1";
  triggerBits_.back().leg1Id=13;
  triggerBits_.back().leg1BitMask=(1<<1)+(1<<2);
  triggerBits_.back().leg1Pt=19;
  triggerBits_.back().leg1Eta=2.1;
  triggerBits_.back().leg2Id=15;
  triggerBits_.back().leg2BitMask=(1<<5);
  triggerBits_.back().leg2Pt=20;
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_IsoMu21_eta2p1_LooseIsoPFTau20_SingleL1";
  triggerBits_.back().leg1Id=13;
  triggerBits_.back().leg1BitMask=(1<<1)+(1<<2);
  triggerBits_.back().leg1Pt=21;
  triggerBits_.back().leg1Eta=2.1;
  triggerBits_.back().leg2Id=15;
  triggerBits_.back().leg2BitMask=(1<<5);
  triggerBits_.back().leg2Pt=20;
  // Single tauh triggers
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_VLooseIsoPFTau120_Trk50_eta2p1";
  triggerBits_.back().leg1Id=15;
  triggerBits_.back().leg1BitMask=(1<<2);//nothing good for 80X
  triggerBits_.back().leg1Pt=120;
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_VLooseIsoPFTau140_Trk50_eta2p1";
  triggerBits_.back().leg1Id=15;
  triggerBits_.back().leg1BitMask=(1<<2);//nothing good for 80X
  triggerBits_.back().leg1Pt=140;
  // tauh tauh triggers
  ///9th tau bit(1<<8) for di-tau dz filter  (should be OK 80X triggers)
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg";
  triggerBits_.back().leg1Id=15;
  triggerBits_.back().leg1BitMask=(1<<8);
  triggerBits_.back().leg1Pt=35;
  triggerBits_.back().leg1Eta=2.1;
  triggerBits_.back().leg2Id=15;
  triggerBits_.back().leg2BitMask=(1<<8);
  triggerBits_.back().leg2Pt=35;
  triggerBits_.back().leg2Eta=2.1;
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_DoubleMediumCombinedIsoPFTau35_Trk1_eta2p1_Reg";
  triggerBits_.back().leg1Id=15;
  triggerBits_.back().leg1BitMask=(1<<8);
  triggerBits_.back().leg1Pt=35;
  triggerBits_.back().leg1Eta=2.1;
  triggerBits_.back().leg2Id=15;
  triggerBits_.back().leg2BitMask=(1<<8);
  triggerBits_.back().leg2Pt=35;
  triggerBits_.back().leg2Eta=2.1;
  triggerBits_.push_back(aTrgData);
  triggerBits_.back().path_name="HLT_DoubleMediumCombinedIsoPFTau40_Trk1_eta2p1_Reg";
  triggerBits_.back().leg1Id=15;
  triggerBits_.back().leg1BitMask=(1<<8);
  triggerBits_.back().leg1Pt=40;
  triggerBits_.back().leg1Eta=2.1;
  triggerBits_.back().leg2Id=15;
  triggerBits_.back().leg2BitMask=(1<<8);
  triggerBits_.back().leg2Pt=40;
  triggerBits_.back().leg2Eta=2.1;

  ////////////////////////////////////////////////////////////
  ///Filter bits to check
  filterBits_.push_back("Flag_HBHENoiseFilter");
  filterBits_.push_back("Flag_HBHENoiseIsoFilter");
  filterBits_.push_back("Flag_EcalDeadCellTriggerPrimitiveFilter");
  filterBits_.push_back("Flag_goodVertices");
  filterBits_.push_back("Flag_eeBadScFilter");
  filterBits_.push_back("Flag_globalTightHalo2016Filter");
  filterBits_.push_back("Flag_chargedHadronTrackResolutionFilter");//bad ChargedCand filter?
  filterBits_.push_back("Flag_muonBadTrackFilter");//bad pfMuonFilter filter?

  return;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::Loop(){

   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntries();   
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      httEvent->clear();

      //if(jentry%1000==0) std::cout<<"Processing "<<jentry<<"th event"<<std::endl;//FIXME
      //Check if event is contained in JSon
      if( !eventInJson() ) continue;
      //if(jentry%1000==0) std::cout<<"\t"<<jentry<<"th event in JSon"<<std::endl;//FIXME

      unsigned int bestPairIndex = Cut(ientry);

      fillEvent();

      hStats->Fill(0);//Number of events analyzed
      hStats->Fill(1,httEvent->getMCWeight());//Sum of weights

      bestPairIndex_ = bestPairIndex;

      if(bestPairIndex<9999){
	//if(jentry%1000==0) std::cout<<"\t"<<jentry<<"th event with good pair"<<std::endl;//FIXME

	///Call pairSelection again to set selection bits for the selected pair.
        pairSelection(bestPairIndex);

	fillJets(bestPairIndex);
	//fillLeptons();//moved
	fillGenLeptons();
	fillPairs(bestPairIndex);
	applyMetRecoilCorrections();//should be done after the best pair is found and thus full event (jets) is defined. Therefore, corrected Met (and releted eg. mT) cannot be used to select the best pair

	HTTPair & bestPair = httPairCollection[0];
        for(unsigned int sysType = (unsigned int)HTTAnalysis::NOMINAL;
	    sysType<(unsigned int)HTTAnalysis::DUMMY_SYS;++sysType){
	  HTTAnalysis::sysEffects type = static_cast<HTTAnalysis::sysEffects>(sysType);
	  computeSvFit(bestPair, type);
	  //break; ///TEST for synch. ntuple
	}
	httTree->Fill();
	hStats->Fill(2);//Number of events saved to ntuple
	hStats->Fill(3,httEvent->getMCWeight());//Sum of weights saved to ntuple
	if(firstWarningOccurence_)//stop to warn once the first pair is found and filled
	  firstWarningOccurence_ = false;
      }
   }

   writePropertiesHeader(leptonPropertiesList);

   writeTriggersHeader(triggerBits_);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
Int_t HTauTauTreeFromNanoBase::Cut(Long64_t entry){

  fillLeptons();
  if( !(httLeptonCollection.size()>1) ) return 9999;
  //std::cout<<"leptons: "<<httLeptonCollection.size()<<std::endl;

  //build pairs
  if(!buildPairs()) return 9999;
  //std::cout<<"pairs: "<<httPairs_.size()<<std::endl;
  std::vector<unsigned int> pairIndexes;
  for(unsigned int iPair=0;iPair<httPairs_.size();++iPair){
    if(pairSelection(iPair)) pairIndexes.push_back(iPair);
  }
  //std::cout<<"passed pairs: "<<pairIndexes.size()<<std::endl;
  
  return bestPair(pairIndexes);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
unsigned int HTauTauTreeFromNanoBase::bestPair(std::vector<unsigned int> &pairIndexes){

  ///Pair are already sorted during the ntuple creation
  if(!pairIndexes.empty()) return pairIndexes[0];
  else return 9999;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::pairSelection(unsigned int iPair){

  ///Requires channel specific implementation based on the following
  ///Baseline+post sync selection as on
  ///https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2015#Baseline_mu_tau_h_AN1
  ///Indexes for multiplexed ID variables taken from   LLRHiggsTauTau/NtupleProducer/plugins/
  ///HTauTauNtuplizer.cc, MuFiller.cc, TauFiller.cc, EleFiller.cc

  return false;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::thirdLeptonVeto(unsigned int signalLeg1Index, unsigned int signalLeg2Index, int leptonPdg, double dRmin){

  TLorentzVector leg1P4 = httLeptonCollection[signalLeg1Index].getP4();
  TLorentzVector leg2P4 = httLeptonCollection[signalLeg2Index].getP4();

  for(unsigned int iLepton=0;iLepton<httLeptonCollection.size();++iLepton){
    if(iLepton==signalLeg1Index || iLepton==signalLeg2Index) continue;
    TLorentzVector leptonP4 = httLeptonCollection[iLepton].getP4();
    double dr = std::min(leg1P4.DeltaR(leptonP4),leg2P4.DeltaR(leptonP4));
    if(dr<dRmin) continue;
    if(leptonPdg == 13 && std::abs(httLeptonCollection[iLepton].getPDGid())==leptonPdg && muonSelection(iLepton) && true) return true;
    /*       (std::abs(httLeptonCollection[signalLeg1].getPDGid())!=leptonPdg || muonSelection(signalLeg1Index)) &&
	     (std::abs(httLeptonCollection[signalLeg2].getPDGid())!=leptonPdg || muonSelection(signalLeg2Index)) ) return true; */ //AP 25/05/17: this part is outdated
    if(leptonPdg == 11 && std::abs(httLeptonCollection[iLepton].getPDGid())==leptonPdg && electronSelection(iLepton) && true) return true;
    /*       (std::abs(httLeptonCollection[signalLeg1].getPDGid())!=leptonPdg || electronSelection(signalLeg1Index)) &&
	     (std::abs(httLeptonCollection[signalLeg2].getPDGid())!=leptonPdg || electronSelection(signalLeg2Index)) ) return true; */ //AP 25/05/17: this part is outdated
  }
  return false;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::extraMuonVeto(unsigned int signalLeg1Index, unsigned int signalLeg2Index, double dRmin){
  return thirdLeptonVeto(signalLeg1Index,signalLeg2Index,13,dRmin);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::extraElectronVeto(unsigned int signalLeg1Index, unsigned int signalLeg2Index, double dRmin){
  return thirdLeptonVeto(signalLeg1Index,signalLeg2Index,11,dRmin);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::muonSelection(unsigned int index){

  TLorentzVector aP4 = httLeptonCollection[index].getP4();

  //Should be fine with NanoAOD?
  // int muonIdBit = 7;//Standard Medium ID
  // if(RunNumber<278808 && RunNumber>100000) muonIdBit = 6;//ICHEP Medium MuonID
  
  bool passSelection = aP4.Pt()>10 && std::abs(aP4.Eta())<2.4 &&
                       std::abs(httLeptonCollection[index].getProperty(PropertyEnum::dz))<0.2 &&
                       std::abs(httLeptonCollection[index].getProperty(PropertyEnum::dxy))<0.045 &&
   		       httLeptonCollection[index].getProperty(PropertyEnum::mediumId)>0 &&
   		       httLeptonCollection[index].getProperty(PropertyEnum::pfRelIso04_all)<0.3;

  return passSelection;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::electronSelection(unsigned int index){

  TLorentzVector aP4 = httLeptonCollection[index].getP4();

  bool passSelection = aP4.Pt()>10 && std::abs(aP4.Eta())<2.5 &&
                       std::abs(httLeptonCollection[index].getProperty(PropertyEnum::dz))<0.2 &&
                       std::abs(httLeptonCollection[index].getProperty(PropertyEnum::dxy))<0.045 &&
                       httLeptonCollection[index].getProperty(PropertyEnum::mvaSpring16GP_WP90)>0.5 &&
                       httLeptonCollection[index].getProperty(PropertyEnum::convVeto)>0.5 &&
                       httLeptonCollection[index].getProperty(PropertyEnum::lostHits)<=1 &&
   		       httLeptonCollection[index].getProperty(PropertyEnum::pfRelIso03_all)<0.3;

  return passSelection;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::jetSelection(unsigned int index, unsigned int bestPairIndex){

  TLorentzVector aP4;
  aP4.SetPtEtaPhiM(Jet_pt[index],
		   Jet_eta[index],
		   Jet_phi[index],
		   Jet_mass[index]);

  bool passSelection = aP4.Pt()>20 && std::abs(aP4.Eta())<4.7 &&
                       Jet_jetId[index]>=1;//it means at least loose
 
  if(bestPairIndex<9999){
    TLorentzVector leg1P4 = httPairs_[bestPairIndex].getLeg1().getP4();
    TLorentzVector leg2P4 = httPairs_[bestPairIndex].getLeg2().getP4();

    passSelection &= aP4.DeltaR(leg1P4) > 0.5 &&
                     aP4.DeltaR(leg2P4) > 0.5;
  }

  return passSelection;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::fillEvent(){


  httEvent->setRun(run);
  httEvent->setEvent(event);
  httEvent->setLS(luminosityBlock);
  httEvent->setNPV(PV_npvs);
  httEvent->setRho(fixedGridRhoFastjetAll);

  httEvent->setAODPV(TVector3(PV_x,PV_y,PV_z));
  httEvent->setRefittedPV(TVector3(PV_x,PV_y,PV_z));//FIXME
  httEvent->setIsRefit(false);//FIXME

  TVector2 metPF;
  metPF.SetMagPhi(MET_pt, MET_phi);
  httEvent->setMET(metPF);

  httEvent->setMETFilterDecision(getMetFilterBits());

  if(b_Pileup_nTrueInt!=nullptr){//Assume that all those are filled for MC
    httEvent->setNPU(Pileup_nTrueInt); //??Pileup_nPU or Pileup_nTrueInt
    httEvent->setMCWeight(genWeight);
    httEvent->setMCatNLOWeight(LHEWeight_originalXWGTUP);//??
    httEvent->setLHE_Ht(LHE_HT);
    httEvent->setLHEnOutPartons(LHE_Njets);
    //FIXMEhttEvent->setGenPV(TVector3(pvGen_x,pvGen_y,pvGen_z));

    TLorentzVector genBosonP4, genBosonVisP4;
    double ptReWeight = 1., ptReWeightSusy = 1.;
    if( findBosonP4(genBosonP4,genBosonVisP4) ){
      //std::cout<<"GenBos found! M="<<genBosonP4.M()<<", visM="<<genBosonVisP4.M()<<std::endl;
      httEvent->setGenBosonP4(genBosonP4,genBosonVisP4);
      ptReWeight = getPtReweight(genBosonP4);
      bool doSUSY = true;
      ptReWeightSusy = getPtReweight(genBosonP4,doSUSY);
    }
    else {
      TLorentzVector topP4, antitopP4;
      findTopP4(topP4, antitopP4);
      ///TT reweighting according to
      ///https://twiki.cern.ch/twiki/bin/view/CMS/TopSystematics#pt_top_Reweighting
      if(topP4.M()>1E-3 && antitopP4.M()>1E-3){
	double topPt = topP4.Perp();
	double antitopPt = antitopP4.Perp();
	double weightTop = exp(0.0615-0.0005*topPt);
	double weightAntitop= exp(0.0615-0.0005*antitopPt);
	ptReWeight = sqrt(weightTop*weightAntitop);
      }
    }
    httEvent->setPtReWeight(ptReWeight);
    httEvent->setPtReWeightSUSY(ptReWeightSusy);

    for(unsigned int iGenPart=0;iGenPart<nGenPart;++iGenPart){
      int absPDGId = std::abs(GenPart_pdgId[iGenPart]);
      if(absPDGId == 25 || absPDGId == 23 || absPDGId == 35 || absPDGId == 36){
	std::vector<unsigned int> daughterIndexes;
	if(!getDirectDaughterIndexes(daughterIndexes,(int)iGenPart)) continue;
	int ntau = 0, nele = 0, nmu = 0;
	for(unsigned int idx=0; idx<daughterIndexes.size(); ++idx) {
	  int pdg_id = std::abs(GenPart_pdgId[daughterIndexes[idx]]);
	  if(pdg_id == 11) nele++;
	  else if(pdg_id == 13) nmu++;
	  else if(pdg_id == 15) {
	    ntau++;
	    unsigned int tauIdx = findFinalCopy(daughterIndexes[idx]);
	    std::vector<unsigned int> tauDaughterIndexes;
	    if(!getDirectDaughterIndexes(tauDaughterIndexes,tauIdx))
	      std::cout<<"isNotFinal, pt1="<<GenPart_pt[(int)daughterIndexes[idx]]
		       <<",  pt2="<<GenPart_pt[(int)daughterIndexes[tauIdx]]
		       <<", #dau1="<<tauDaughterIndexes.size()<<std::endl;
	    //get DM and then translate it on the basics DM
	    int dm = genTauDecayMode(tauDaughterIndexes);
	    switch ( dm ){
	    case HTTAnalysis::tauDecayMuon :
	      nmu++;
	      break;
	    case HTTAnalysis::tauDecaysElectron :
	      nele++;
	      break;
	    default :
	      break;
	    }
	  }
	}
	// determine H/Z decay mode
	int hzDecay = 8;//other
	if (ntau == 0) {
	    if(nele == 0 && nmu == 2) hzDecay = 7;
	    else if(nele == 2 && nmu == 0) hzDecay = 6;
	    else hzDecay = 8;
	}
	else if(ntau == 2) {
	    if(nmu == 0 && nele == 0) hzDecay = 2;    
	    else if(nmu == 0 && nele == 1) hzDecay = 1;    
	    else if(nmu == 0 && nele == 2) hzDecay = 4;    
	    else if(nmu == 1 && nele == 0) hzDecay = 0;    
	    else if(nmu == 1 && nele == 1) hzDecay = 5;    
	    else if(nmu == 2 && nele == 0) hzDecay = 3;    
	}
	httEvent->setDecayModeBoson(hzDecay);
      }
      else if(absPDGId == 24) {
	std::vector<unsigned int> daughterIndexes;
	if(!getDirectDaughterIndexes(daughterIndexes,(int)iGenPart)) continue;
	int ntau = 0, nele = 0, nmu = 0, nquark = 0;
	for(unsigned int idx=0; idx<daughterIndexes.size(); ++idx) {
	  int pdg_id = std::abs(GenPart_pdgId[daughterIndexes[idx]]);
	  if(pdg_id < 5 ) nquark++;
	  else if(pdg_id == 11) nele++;
	  else if(pdg_id == 13) nmu++;
	  else if(pdg_id == 15) {
	    ntau++;
	    unsigned int tauIdx = findFinalCopy(daughterIndexes[idx]);
	    std::vector<unsigned int> tauDaughterIndexes;
	    if(!getDirectDaughterIndexes(tauDaughterIndexes,tauIdx))
	      std::cout<<"isNotFinal, pt1="<<GenPart_pt[(int)daughterIndexes[idx]]
		       <<",  pt2="<<GenPart_pt[(int)daughterIndexes[tauIdx]]
		       <<", #dau1="<<tauDaughterIndexes.size()<<std::endl;
	    //get DM and then translate it on the basics DM
	    int dm = genTauDecayMode(tauDaughterIndexes);
	    switch ( dm ){
	    case HTTAnalysis::tauDecayMuon :
	      nmu++;
	      break;
	    case HTTAnalysis::tauDecaysElectron :
	      nele++;
	      break;
	    }
	  }
	}
	// determine W decay mode
	int wDecay = 6;//other
	if(nquark == 2 && (nmu+nele+ntau)==0) wDecay = 0;
	else if(nquark ==0 && ntau == 0) {
	  if(nele == 0 && nmu == 1) wDecay = 1;
	  else if(nele == 1 && nmu == 0) wDecay = 2;
	  else wDecay = 6;
	}
	else if(nquark==0 && ntau == 1) {
	  if(nmu == 0 && nele == 0) wDecay = 5;    
	  if(nmu == 0 && nele == 1) wDecay = 4;    
	  if(nmu == 1 && nele == 0) wDecay = 3;
	  else wDecay = 6;
	}
	httEvent->setDecayModeBoson(10+wDecay);
      }
      else if(GenPart_pdgId[iGenPart]==15){
	TLorentzVector p4;
	p4.SetPtEtaPhiM(GenPart_pt[iGenPart],
			GenPart_eta[iGenPart],
			GenPart_phi[iGenPart],
			1.777);//should use pdg mass as masses below 10GeV are zeroed
	//do not consider low momentum candidates??
	if( !(p4.P()>10) ) continue;
	//find direct daughters
	std::vector<unsigned int> daughterIndexes;
	if(!getDirectDaughterIndexes(daughterIndexes,(int)iGenPart)) continue;
	//get DM and then translate it on the basics DM
	int dm = genTauDecayMode(daughterIndexes);
	switch ( dm ){
	case HTTAnalysis::tauDecayMuon ://0
	  httEvent->setDecayModeMinus(0);
	  break;
	case HTTAnalysis::tauDecaysElectron ://1
	  httEvent->setDecayModeMinus(1);
	  break;
	default: //2
	  httEvent->setDecayModeMinus(2);
	  break;
	}
      }
      if(GenPart_pdgId[iGenPart]==-15){
	TLorentzVector p4;
	p4.SetPtEtaPhiM(GenPart_pt[iGenPart],
			GenPart_eta[iGenPart],
			GenPart_phi[iGenPart],
			1.777);//should use pdg mass as masses below 10GeV are zeroed
	//do not consider low momentum candidates??
	if( !(p4.P()>10) ) continue;
	//find direct daughters
	std::vector<unsigned int> daughterIndexes;
	if(!getDirectDaughterIndexes(daughterIndexes,(int)iGenPart)) continue;
	//get DM and then translate it on the basics DM
	int dm = genTauDecayMode(daughterIndexes);
	switch ( dm ){
	case HTTAnalysis::tauDecayMuon ://0
	  httEvent->setDecayModePlus(0);
	  break;
	case HTTAnalysis::tauDecaysElectron ://1
	  httEvent->setDecayModePlus(1);
	  break;
	default: //2
	  httEvent->setDecayModePlus(2);
	  break;
	}
      }      
    }
  }

  std::string fileName(fChain->GetCurrentFile()->GetName());
  HTTEvent::sampleTypeEnum aType = HTTEvent::DUMMY;
  httEvent->setSampleType(aType);

}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::fillJets(unsigned int bestPairIndex){

  httJetCollection.clear();

  for(unsigned int iJet=0;iJet<nJet;++iJet){

    if(!jetSelection(iJet, bestPairIndex)) continue;
    HTTParticle aJet;

    TLorentzVector p4;
    p4.SetPtEtaPhiM(Jet_pt[iJet],
		    Jet_eta[iJet],
		    Jet_phi[iJet],
		    Jet_mass[iJet]);
    std::vector<Double_t> aProperties = getProperties(leptonPropertiesList, iJet,"Jet");
    ///Set jet PDG id by hand
    aProperties[(unsigned int)PropertyEnum::pdgId] = 98.0;
    ///JEC uncertaintes
    // Up
    for(unsigned int iUnc=0; iUnc<jecUncertList.size(); ++iUnc)
      aProperties.push_back(getJecUnc(iJet, jecUncertList[iUnc] ,true));
    // Down
    for(unsigned int iUnc=0; iUnc<jecUncertList.size(); ++iUnc)
      aProperties.push_back(getJecUnc(iJet, jecUncertList[iUnc] ,false));

    aJet.setProperties(aProperties);

    aJet.setP4(p4);
    aJet.setProperties(aProperties);
    httJetCollection.push_back(aJet);
  }
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::fillLeptons(){

  httLeptonCollection.clear();

  //Muons
  for(unsigned int iMu=0; iMu<nMuon; ++iMu){
    if( !(Muon_pt[iMu]>5) ) continue;
    HTTParticle aLepton;
    TLorentzVector p4;
    p4.SetPtEtaPhiM(Muon_pt[iMu],
		    Muon_eta[iMu],
		    Muon_phi[iMu],
		    0.10566); //muon mass
    TVector3 pca;//FIXME: can partly recover with ip3d and momentum?
    aLepton.setP4(p4);
    aLepton.setChargedP4(p4);//same as p4 for muon
    //aLepton.setNeutralP4(p4Neutral); not defined for muon
    aLepton.setPCA(pca);
    std::vector<Double_t> aProperties = getProperties(leptonPropertiesList, iMu, "Muon");
    aLepton.setProperties(aProperties);
    httLeptonCollection.push_back(aLepton);
  }//Muons
  //Electrons
  for(unsigned int iEl=0; iEl<nElectron; ++iEl){
    if( !(Electron_pt[iEl]>7) ) continue;
    HTTParticle aLepton;
    TLorentzVector p4;
    p4.SetPtEtaPhiM(Electron_pt[iEl],
		    Electron_eta[iEl],
		    Electron_phi[iEl],
		    0.51100e-3); //electron mass
    TVector3 pca;//FIXME: can partly recover with ip3d and momentum?
    aLepton.setP4(p4);
    aLepton.setChargedP4(p4);//same as p4 for electron
    //aLepton.setNeutralP4(p4Neutral); not defined for electron
    aLepton.setPCA(pca);
    std::vector<Double_t> aProperties = getProperties(leptonPropertiesList, iEl, "Electron");
    aLepton.setProperties(aProperties);
    httLeptonCollection.push_back(aLepton);
  }//Electrons
  //Taus
  for(unsigned int iTau=0; iTau<nTau; ++iTau){
    if( !(Tau_pt[iTau]>18) ) continue;
    if( !Tau_idDecayMode[iTau] ) continue; //oldDMs
    HTTParticle aLepton;
    TLorentzVector p4;
    p4.SetPtEtaPhiM(Tau_pt[iTau],
		    Tau_eta[iTau],
		    Tau_phi[iTau],
		    Tau_mass[iTau]);
    TVector3 pca;//FIXME: can partly recover with dxy,dz and momentum?
    aLepton.setP4(p4);
    TLorentzVector chargedP4;//approximate by leadTrack p4
    double leadTkPhi = Tau_phi[iTau]+Tau_leadTkDeltaPhi[iTau];
    if(std::abs(leadTkPhi)>TMath::Pi())//limit phi to the [-pi,pi] range
      leadTkPhi -= 2.*TMath::Pi()*std::round(leadTkPhi/(2.*TMath::Pi()));
    chargedP4.SetPtEtaPhiM(Tau_pt[iTau]*Tau_leadTkPtOverTauPt[iTau],
			   Tau_eta[iTau]+Tau_leadTkDeltaEta[iTau],
			   leadTkPhi,
			   0.13957); //pi+/- mass
    aLepton.setChargedP4(chargedP4);
    aLepton.setNeutralP4(p4-chargedP4);
    aLepton.setPCA(pca);
    std::vector<Double_t> aProperties = getProperties(leptonPropertiesList, iTau, "Tau");
    aLepton.setProperties(aProperties);
    //FIXME: for synch tests, should be removed(?) -->
    bool vlMVA=(((int)Tau_idMVAoldDM[iTau] & (1<<0))==(1<<0));
    bool vlMVANew = (((int)Tau_idMVAnewDM[iTau] & (1<<0))==(1<<0));
    //bool vlMVAR03 = (((int)Tau_idMVAoldDMdR03[iTau] & (1<<0))==(1<<0));
    bool lComb = (Tau_rawIso[iTau]<2.5 && Tau_photonsOutsideSignalCone[iTau]<0.1*Tau_pt[iTau]);
    if( !(vlMVA||vlMVANew||lComb) ) continue;
    //if( !(vlMVA) ) continue;
    //<--
    httLeptonCollection.push_back(aLepton);
  }//Taus
  //Sort leptons
  std::sort(httLeptonCollection.begin(),httLeptonCollection.end(),compareLeptons);

}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::fillGenLeptons(){

  httGenLeptonCollection.clear();

  if(!fChain->FindBranch("nGenPart")) return;

  for(unsigned int iGenPart=0;iGenPart<nGenPart;++iGenPart){
    if(std::abs(GenPart_pdgId[iGenPart])!=15) continue;
    TLorentzVector p4;
    p4.SetPtEtaPhiM(GenPart_pt[iGenPart],
		    GenPart_eta[iGenPart],
		    GenPart_phi[iGenPart],
		    1.777);//should use pdg mass as masses below 10GeV are zeroed
    //do not consider low momentum candidates??
    if( !(p4.P()>10) ) continue;
    //find direct daughters
    std::vector<unsigned int> daughterIndexes;
    bool isFinalTau=getDirectDaughterIndexes(daughterIndexes,(int)iGenPart);
    if(!isFinalTau) continue;

    HTTParticle aLepton;
    aLepton.setP4(p4);
    aLepton.setChargedP4(getGenComponentP4(daughterIndexes,1));
    aLepton.setNeutralP4(getGenComponentP4(daughterIndexes,0));
    //TVector3 pca(genpart_pca_x->at(iGenPart), genpart_pca_y->at(iGenPart), genpart_pca_z->at(iGenPart));
    //aLepton.setPCA(pca);

    //std::vector<Double_t> aProperties = getProperties(genLeptonPropertiesList, iGenPart);
    //set properties by hand (keep correct order)
    std::vector<Double_t> aProperties;
    aProperties.push_back(GenPart_pdgId[iGenPart]);
    aProperties.push_back(genTauDecayMode(daughterIndexes));    
    aLepton.setProperties(aProperties);

    httGenLeptonCollection.push_back(aLepton);

  }
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
TLorentzVector HTauTauTreeFromNanoBase::getGenComponentP4(std::vector<unsigned int> &indexes, unsigned int iAbsCharge){

  TLorentzVector aNeutralP4, aChargedP4, aLeptonP4;

  for(unsigned int idx=0;idx<indexes.size();++idx){
    unsigned int iGenPart = indexes[idx];
    unsigned int pdg_id = std::abs(GenPart_pdgId[iGenPart]);
    if(pdg_id == 11 || pdg_id == 13) 
      aLeptonP4.SetPtEtaPhiM(GenPart_pt[iGenPart],
			     GenPart_eta[iGenPart],
			     GenPart_phi[iGenPart],
			     (pdg_id==11?0.51100e-3:0.10566));//set mass
    else if(pdg_id == 211 || pdg_id == 321 )
      aChargedP4.SetPtEtaPhiM(GenPart_pt[iGenPart],
			      GenPart_eta[iGenPart],
			      GenPart_phi[iGenPart],
			      (pdg_id==211?0.1396:0.4937));//set mass
    else if(pdg_id == 111 || pdg_id == 130 || pdg_id == 310 || pdg_id == 311 )
      aNeutralP4.SetPtEtaPhiM(GenPart_pt[iGenPart],
			      GenPart_eta[iGenPart],
			      GenPart_phi[iGenPart],
			      (pdg_id==111?0.1350:0.4976));//set mass
  }

  TLorentzVector aP4;
  if(iAbsCharge==0) aP4 = aNeutralP4;
  else if(aChargedP4.E()>0) aP4 = aChargedP4;
  else if(aLeptonP4.E()>0) aP4 = aLeptonP4;

  return aP4;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::fillPairs(unsigned int bestPairIndex){

  httPairCollection.clear();
  if(bestPairIndex<9999 &&  bestPairIndex<httPairs_.size())
    httPairCollection.push_back(httPairs_[bestPairIndex]);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::buildPairs(){

  httPairs_.clear();
  for(unsigned int iL1=0; iL1<httLeptonCollection.size()-1; ++iL1){
    for(unsigned int iL2=iL1+1; iL2<httLeptonCollection.size(); ++iL2){
      if( !(httLeptonCollection[iL1].getP4().DeltaR(httLeptonCollection[iL2].getP4())>0.3) ) continue;
      TLorentzVector p4 = httLeptonCollection[iL1].getP4()+httLeptonCollection[iL2].getP4();
      if( !(p4.M()>0) ) continue;
      TVector2 met; met.SetMagPhi(MET_pt, MET_phi);
      double mTLeg1 = TMath::Sqrt(2.*httLeptonCollection[iL1].getP4().Pt()*MET_pt*(1.-TMath::Cos(httLeptonCollection[iL1].getP4().Phi()-MET_phi)));
      double mTLeg2 = TMath::Sqrt(2.*httLeptonCollection[iL2].getP4().Pt()*MET_pt*(1.-TMath::Cos(httLeptonCollection[iL2].getP4().Phi()-MET_phi)));
      HTTPair aHTTpair;
      aHTTpair.setP4(p4);
      aHTTpair.setMET(met);
      aHTTpair.setMETMatrix(MET_covXX, MET_covXY, MET_covXY, MET_covYY);
      aHTTpair.setMTLeg1(mTLeg1);
      aHTTpair.setMTLeg2(mTLeg2);
      aHTTpair.setLeg1(httLeptonCollection.at(iL1),iL1);
      aHTTpair.setLeg2(httLeptonCollection.at(iL2),iL2);
      httPairs_.push_back(aHTTpair); 
    }
  }
  std::sort(httPairs_.begin(),httPairs_.end(),comparePairs);

  return !httPairs_.empty();
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
template<class T> T HTauTauTreeFromNanoBase::getBranchValue(const char *branchAddress, unsigned int index){

  std::vector<T> *aVector = *(std::vector<T> **)(branchAddress);

  if(aVector->size()<=index){
    //std::cout<<"Index - size mismatch "<<std::endl;
    return 0;
  }
  return aVector->at(index);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
Double_t  HTauTauTreeFromNanoBase::getProperty(std::string name, unsigned int index, std::string colType){

  if(name=="mc_match") return getMCMatching(index,colType);
  
  if(name=="isGoodTriggerType") return getTriggerMatchning(index,false,colType);
  if(name=="FilterFired") return getTriggerMatchning(index,true,colType);//some overhead due to calling it again with a different option, but kept for backward compatibility (and debug)

  if(colType=="Electron"){
    if(name.find("Muon_")!=std::string::npos ||
       name.find("Tau_")!=std::string::npos ||
       name.find("Jet_")!=std::string::npos ||
       false)
      return 0;
  }
  else if(colType=="Muon"){
    if(name.find("Electron_")!=std::string::npos ||
       name.find("Tau_")!=std::string::npos ||
       name.find("Jet_")!=std::string::npos ||
       false)
      return 0;
  }
  else if(colType=="Tau"){
    if(name.find("Electron_")!=std::string::npos ||
       name.find("Muon_")!=std::string::npos ||
       name.find("Jet_")!=std::string::npos ||
       name.find("pfRelIso03_all")!=std::string::npos ||
       name.find("sip3d")!=std::string::npos ||
       false)
      return 0;
    if(name.find("pdgId")!=std::string::npos){
      TBranch *branch = fChain->GetBranch("Tau_charge");
      if(!branch){
	if(firstWarningOccurence_)
	  std::cout<<"Branch: Tau_charge not found in the TTree, return pdgId=-15"<<std::endl;
	return -15;
      }
      TLeaf *leaf = branch->FindLeaf("Tau_charge");
      int charge = leaf!=nullptr ? leaf->GetValue(index) : 0;
      return charge>0 ? -15 : 15; 
    }      
  }
  else if(colType=="Jet"){
    if(name.find("Jet_")==std::string::npos)
      return 0;
  }

  if(colType!="" && name.find(colType+"_")==std::string::npos){
    name = colType+"_"+name;
  }
  TBranch *branch = fChain->GetBranch(name.c_str());
  if(!branch){
    if(firstWarningOccurence_)
      std::cout<<"Branch: "<<name<<" not found in the TTree."<<std::endl;
    return 0;
  }

  const char *branchAddress = branch->GetAddress();
  std::string branchClass(branch->GetClassName());
  if(branchClass=="vector<double>") return getBranchValue<double>(branchAddress, index);
  if(branchClass=="vector<int>") return getBranchValue<int>(branchAddress, index);
  if(branchClass=="vector<bool>") return getBranchValue<bool>(branchAddress, index);
  if(branchClass=="vector<Long64_t>") return getBranchValue<Long64_t>(branchAddress, index);
  //otherwise assume that index value makes sense
  if(branchClass==""){//check leafs
    TLeaf *leaf = branch->FindLeaf(name.c_str());
    return leaf!=nullptr ? leaf->GetValue(index) : 0;
  }

  return 0;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::initJecUnc(std::string correctionFile){
  ///uncertainties
  const int nsrc = 27;//+6;
  const std::string srcnames[nsrc] =
    {"AbsoluteStat", "AbsoluteScale", "AbsoluteMPFBias",
     "Fragmentation",
     "SinglePionECAL", "SinglePionHCAL",
     "FlavorQCD",
     "TimePtEta",
     "RelativeJEREC1", "RelativeJEREC2", "RelativeJERHF",
     "RelativePtBB","RelativePtEC1", "RelativePtEC2", "RelativePtHF",
     "RelativeBal",
     "RelativeFSR",
     "RelativeStatFSR","RelativeStatEC", "RelativeStatHF",
     "PileUpDataMC",
     "PileUpPtRef", "PileUpPtBB", "PileUpPtEC1", "PileUpPtEC2", "PileUpPtHF",
     //"SubTotalPileUp", "SubTotalRelative", "SubTotalPt", "SubTotalScale", "SubTotalAbsolute", "SubTotalMC", //MB: 6 quadratic sums of subsets of uncertainties
     "Total"};//MB: sum of all uncertainties in quadrature

  ofstream outputFile("JecUncEnum.h");
  outputFile<<"enum class JecUncEnum { ";
  for(unsigned int isrc = 0; isrc < nsrc; isrc++) {
    JetCorrectorParameters *p = new JetCorrectorParameters(correctionFile, srcnames[isrc]);
    JetCorrectionUncertainty *unc = new JetCorrectionUncertainty(*p);
    jecUncertList.push_back(srcnames[isrc]);
    jecUncerts.push_back(unc);
    outputFile<<srcnames[isrc]<<" = "<<isrc<<", "<<std::endl;
  }
  outputFile<<"NONE"<<" = "<<nsrc<<std::endl;
  outputFile<<"};"<<std::endl;
  outputFile.close();
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
double HTauTauTreeFromNanoBase::getJecUnc(unsigned int index, std::string name,bool up){
  if(b_nGenPart==nullptr) return 0;//MB: do not check it for data
  double result = 0;
  double jetpt = Jet_pt[index];
  double jeteta =  Jet_eta[index];
  for(unsigned int isrc = 0; isrc < jecUncertList.size(); isrc++) {
    if(jecUncertList[isrc]==name){
      JetCorrectionUncertainty *unc = jecUncerts[isrc];
      unc->setJetPt(jetpt);
      unc->setJetEta(jeteta);
      result = unc->getUncertainty(up);
      break;
    }
  }
  return result;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void  HTauTauTreeFromNanoBase::writeTriggersHeader(const std::vector<TriggerData> &triggerBits){

  ofstream outputFile("TriggerEnum.h");

  outputFile<<"enum class TriggerEnum { ";
  for(unsigned int iBit=0;iBit<triggerBits.size();++iBit){
    std::string name = triggerBits[iBit].path_name;
    outputFile<<name<<" = "<<iBit<<", "<<std::endl;
  }
  outputFile<<"NONE"<<" = "<<triggerBits.size()<<std::endl;
  outputFile<<"};"<<std::endl;

  outputFile.close();
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void  HTauTauTreeFromNanoBase::writePropertiesHeader(const std::vector<std::string> & propertiesList){

  ofstream outputFile("PropertyEnum.h");

  outputFile<<"enum class PropertyEnum { ";
  for(unsigned int iItem=0;iItem<propertiesList.size();++iItem){
    std::string name = propertiesList[iItem];
    std::string pattern = "daughters_";
    if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
    pattern = "Daughters";
    if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
    pattern = "jets_";
    if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
    pattern = "Jet_";
    if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
    pattern = "Tau_";
    if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
    pattern = "Muon_";
    if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
    pattern = "Electron_";
    if(name.find(pattern)!=std::string::npos) name.erase(name.find(pattern), pattern.size());
    outputFile<<name<<" = "<<iItem<<", "<<std::endl;
  }
  outputFile<<"NONE"<<" = "<<propertiesList.size()<<std::endl;
  outputFile<<"};"<<std::endl;
  outputFile.close();
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
std::vector<Double_t>  HTauTauTreeFromNanoBase::getProperties(const std::vector<std::string> & propertiesList,
							      unsigned int index,
							      std::string colType){

  std::vector<Double_t> aProperties;

  for(auto propertyName:propertiesList){
    aProperties.push_back(getProperty(propertyName,index,colType));
  }

  return aProperties;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
int HTauTauTreeFromNanoBase::getMCMatching(unsigned int index, std::string colType){

  
  if(b_nGenPart==nullptr) return -999;
  if(nGenPart==0) return -999;

  double dR = 100;
  int genFlags = 0;
  int gen_ind = -1;
  unsigned int absPdgId_2 = 0;
  TLorentzVector p4_1, p4_2;
  if(colType=="Muon"){
    if(index>=nMuon) return -999;
    p4_1.SetPtEtaPhiM(Muon_pt[index],
		      Muon_eta[index],
		      Muon_phi[index],
		      0.10566); //muon mass
  }
  else if(colType=="Electron"){
    if(index>=nElectron) return -999;
    p4_1.SetPtEtaPhiM(Electron_pt[index],
		      Electron_eta[index],
		      Electron_phi[index],
		      0.51100e-3); //electron mass
  }
  else if(colType=="Tau"){
    if(index>=nTau) return -999;
    p4_1.SetPtEtaPhiM(Tau_pt[index],
		      Tau_eta[index],
		      Tau_phi[index],
		      Tau_mass[index]);
  }
  else
    return -999;
  
  for(unsigned int iGen = 0; iGen < nGenPart; ++iGen) {
    //consider only leptons with Pt>8GeV
    if( !(GenPart_pt[iGen]>8) ) continue;
    unsigned int absPdgId = std::abs(GenPart_pdgId[iGen]);
    if(absPdgId!=11 && absPdgId!=13 && absPdgId!=15 ) continue;
    //consider only final copies
    std::vector<unsigned int> daughterIndexes;
    bool isFinal=getDirectDaughterIndexes(daughterIndexes,(int)iGen);
    if(!isFinal) continue;
    int genFlags_tmp=0;
    TLorentzVector p4_tmp;
    if(absPdgId==11 || absPdgId==13) {
      p4_tmp.SetPtEtaPhiM(GenPart_pt[iGen],
			GenPart_eta[iGen],
			GenPart_phi[iGen],
			(absPdgId==11?0.51100e-3:0.10566));
      /*MB: Stauts flags in NanoAOD:
	"statusFlags().isLastCopyBeforeFSR()                  * 16384 +" <<14
	"statusFlags().isLastCopy()                           * 8192  +" <<12
	"statusFlags().isFirstCopy()                          * 4096  +" <<12
	"statusFlags().fromHardProcessBeforeFSR()             * 2048  +" <<11
	"statusFlags().isDirectHardProcessTauDecayProduct()   * 1024  +" <<10
	"statusFlags().isHardProcessTauDecayProduct()         * 512   +" <<9
	"statusFlags().fromHardProcess()                      * 256   +" <<8
	"statusFlags().isHardProcess()                        * 128   +" <<7
	"statusFlags().isDirectHadronDecayProduct()           * 64    +" <<6
	"statusFlags().isDirectPromptTauDecayProduct()        * 32    +" <<5
	"statusFlags().isDirectTauDecayProduct()              * 16    +" <<4
	"statusFlags().isPromptTauDecayProduct()              * 8     +" <<3
	"statusFlags().isTauDecayProduct()                    * 4     +" <<2
	"statusFlags().isDecayedLeptonHadron()                * 2     +" <<1
	"statusFlags().isPrompt()                             * 1      ",<<0
      */
      genFlags_tmp = GenPart_statusFlags[iGen];

      /*MB: old code from times when statusFlags were not available
      int mother_idx = GenPart_genPartIdxMother[findFirstCopy(iGen)];
      mother_idx = mother_idx<0?iGen:mother_idx;//protection
      unsigned int mother_absPdgId = std::abs(GenPart_pdgId[mother_idx]);
      bool mother_fromME = (GenPart_status[mother_idx]==3 || (GenPart_status[mother_idx]>20 && GenPart_status[mother_idx]<30)); // from matrix element
      if(mother_absPdgId==15) {//daugther of tau
	genFlags_tmp |= (1 << 2); //best possible approx of what is done with MiniAOD
	genFlags_tmp |= (1 << 5); //best possible approx of what is done with MiniAOD
      }
      else if(mother_absPdgId==23 || //assume that it is prompt if from Z
	      mother_absPdgId==24 || //W+-
	      mother_absPdgId==25 || //h
	      mother_absPdgId==35 || //H
	      mother_absPdgId==36 || //A
	      mother_absPdgId==32 || //Z'
	      mother_absPdgId==33 || //Z''
	      mother_absPdgId==34 || //W'
	      mother_absPdgId==37 || //H+-
	      mother_fromME ||  //or comes from ME
	      false) {
	genFlags_tmp |= (1 << 0);
      }
      */
    }
    else {//tau
      int dm = genTauDecayMode(daughterIndexes);
      if(dm!=HTTAnalysis::tauDecayMuon && dm!=HTTAnalysis::tauDecaysElectron){
	genFlags_tmp = GenPart_statusFlags[iGen];

	/*MB: old code from times when statusFlags were not available
	int mother_idx = GenPart_genPartIdxMother[findFirstCopy(iGen)];
	mother_idx = mother_idx<0?iGen:mother_idx;//protection
	unsigned int mother_absPdgId = std::abs(GenPart_pdgId[mother_idx]);
	bool mother_fromME = (GenPart_status[mother_idx]==3 || (GenPart_status[mother_idx]>20 && GenPart_status[mother_idx]<30)); // from matrix element
	if(mother_absPdgId==23 || //assume that it is prompt if from Z
	   mother_absPdgId==24 || //W+-
	   mother_absPdgId==25 || //h
	   mother_absPdgId==35 || //H
	   mother_absPdgId==36 || //A
	   mother_absPdgId==32 || //Z'
	   mother_absPdgId==33 || //Z''
	   mother_absPdgId==34 || //W'
	   mother_absPdgId==37 || //H+-
	   mother_fromME || //or comes from ME
	   false) {
	  genFlags_tmp |= (1 << 0);
	}
	*/
	for(unsigned int idau=0; idau<daughterIndexes.size(); ++idau){
	  //it assumes that neutrinos are excluded
	  TLorentzVector p4_dau;
	  p4_dau.SetPtEtaPhiM(GenPart_pt[daughterIndexes[idau]],
			      GenPart_eta[daughterIndexes[idau]],
			      GenPart_phi[daughterIndexes[idau]],
			      (GenPart_pdgId[daughterIndexes[idau]]==111?0.1350:0.1396));//use pi mass as simple appox
	  p4_tmp+=p4_dau;
	}
      }
    }
    if( !(p4_tmp.P()>0) ) continue;
    if(dR > p4_1.DeltaR(p4_tmp)) {
      dR = p4_1.DeltaR(p4_tmp); 
      p4_2 = p4_tmp;
      genFlags = genFlags_tmp;
      gen_ind = (int)iGen;
      absPdgId_2 = absPdgId==15?66615:absPdgId;
    }
  }

  if( gen_ind < 0 || !(p4_2.P()>0) ) return 6;
  if(dR > 0.2) return 6;

  if(absPdgId_2 == 11 && p4_2.Pt() > 8 && (genFlags & (1<<0)) == (1<<0)) return 1;
  if(absPdgId_2 == 13 && p4_2.Pt() > 8 && (genFlags & (1<<0)) == (1<<0)) return 2;
  if(absPdgId_2 == 11 && p4_2.Pt() > 8 && (genFlags & (1<<5)) == (1<<5)) return 3;
  if(absPdgId_2 == 13 && p4_2.Pt() > 8 && (genFlags & (1<<5)) == (1<<5)) return 4;
  if(absPdgId_2 == 66615 && p4_2.Pt() > 15 && (genFlags & (1<<0)) == (1<<0)) return 5;
  return 6;

}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::isGoodToMatch(unsigned int ind){
  
  //MB: This method not mandatory as selection of good to match is applied upper in the flow. Anyway, an trivial implementation kept.
  return true;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
int HTauTauTreeFromNanoBase::getTriggerMatchning(unsigned int index, bool checkBit, std::string colType){

  TLorentzVector p4_1;
  unsigned int particleId=0;
  double dRmax=0.5;//was 0.25
  if(colType=="Muon"){
    if(index>=nMuon) return 0;
    particleId=13;
    p4_1.SetPtEtaPhiM(Muon_pt[index],
		      Muon_eta[index],
		      Muon_phi[index],
		      0.10566); //muon mass
  }
  else if(colType=="Electron"){
    if(index>=nElectron) return 0;
    particleId=11;
    p4_1.SetPtEtaPhiM(Electron_pt[index],
		      Electron_eta[index],
		      Electron_phi[index],
		      0.51100e-3); //electron mass
  }
  else if(colType=="Tau"){
    if(index>=nTau) return 0;
    particleId=15;
    p4_1.SetPtEtaPhiM(Tau_pt[index],
		      Tau_eta[index],
		      Tau_phi[index],
		      Tau_mass[index]);
  }
  else
    return 0;

  int firedBits = 0;
  for(unsigned int iTrg=0; iTrg<triggerBits_.size(); ++iTrg){
    bool decision = false;
    //check if trigger is fired
    TBranch *branch = fChain->GetBranch(triggerBits_[iTrg].path_name.c_str());
    if(branch!=nullptr){
      TLeaf *leaf = branch->FindLeaf(triggerBits_[iTrg].path_name.c_str());
      decision = leaf!=nullptr ? leaf->GetValue() : false;
    }
    if(!decision) continue; // do not check rest if trigger is not fired
    //check legs
    decision = false;
    //first leg
    if(particleId==triggerBits_[iTrg].leg1Id){
      for(unsigned int iObj=0; iObj<nTrigObj; ++iObj){
	if(TrigObj_id[iObj]!=(int)particleId) continue;
	TLorentzVector p4_trg;
	p4_trg.SetPtEtaPhiM(TrigObj_pt[iObj],
			    TrigObj_eta[iObj],
			    TrigObj_phi[iObj],
			    0.);
	if( !(p4_1.DeltaR(p4_trg)<dRmax) ) continue;
	if( triggerBits_[iTrg].leg1Pt>0 && !(TrigObj_pt[iObj]>triggerBits_[iTrg].leg1Pt) ) continue;
	if( triggerBits_[iTrg].leg1Eta>0 && !(std::abs(TrigObj_eta[iObj])<triggerBits_[iTrg].leg1Eta) ) continue;
	if( triggerBits_[iTrg].leg1L1Pt>0 && !(TrigObj_l1pt[iObj]>triggerBits_[iTrg].leg1L1Pt) ) continue;
	if( checkBit && !( ((int)TrigObj_filterBits[iObj] & triggerBits_[iTrg].leg1BitMask)==triggerBits_[iTrg].leg1BitMask) ) continue;
	decision = true;
	if(decision) break;
      }
    }
    if(decision){
      firedBits |= (1<<iTrg);
      continue;
    }
    decision = false;
    //second leg
    if(particleId==triggerBits_[iTrg].leg2Id){
      for(unsigned int iObj=0; iObj<nTrigObj; ++iObj){
	if(TrigObj_id[iObj]!=(int)particleId) continue;
	TLorentzVector p4_trg;
	p4_trg.SetPtEtaPhiM(TrigObj_pt[iObj],
			    TrigObj_eta[iObj],
			    TrigObj_phi[iObj],
			    0.);
	if( !(p4_1.DeltaR(p4_trg)<dRmax) ) continue;
	if( triggerBits_[iTrg].leg2Pt>0 && !(TrigObj_pt[iObj]>triggerBits_[iTrg].leg2Pt) ) continue;
	if( triggerBits_[iTrg].leg2Eta>0 && !(std::abs(TrigObj_eta[iObj])<triggerBits_[iTrg].leg2Eta) ) continue;
	if( triggerBits_[iTrg].leg2L1Pt>0 && !(TrigObj_l1pt[iObj]>triggerBits_[iTrg].leg2L1Pt) ) continue;
	if( checkBit && !( ((int)TrigObj_filterBits[iObj] & triggerBits_[iTrg].leg2BitMask)==triggerBits_[iTrg].leg2BitMask) ) continue;
	decision = true;
	if(decision) break;
      }
    }
    if(decision){
      firedBits |= (1<<iTrg);
      continue;
    }
  }

  return firedBits;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
int HTauTauTreeFromNanoBase::getMetFilterBits(){
  int firedBits=0;
  for(unsigned int iFlt=0; iFlt<filterBits_.size(); ++iFlt){
    bool decision = false;
    //check if trigger is fired
    TBranch *branch = fChain->GetBranch(filterBits_[iFlt].c_str());
    if(branch!=nullptr){
      TLeaf *leaf = branch->FindLeaf(filterBits_[iFlt].c_str());
      decision = leaf!=nullptr ? leaf->GetValue() : false;
    }
    if(decision) firedBits |= (1<<iFlt);
  }
  return firedBits;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::findBosonP4(TLorentzVector &bosonP4, TLorentzVector &visBosonP4){

  bosonP4.SetXYZM(0,0,0,0);
  visBosonP4.SetXYZM(0,0,0,0);

  if(b_nGenPart==nullptr) return false;

  for(unsigned int iGen = 0; iGen < nGenPart; ++iGen) {
    unsigned int absPdgId = std::abs(GenPart_pdgId[iGen]);
    if(absPdgId != 11 && absPdgId != 13 && absPdgId != 15 && //charged leptons
       absPdgId != 12 && absPdgId != 14 && absPdgId != 16 ) //neutrinos
      continue;
    //mass is stored only m>10GeV and photons m>1GeV
    double mass = GenPart_mass[iGen]>0 ? GenPart_mass[iGen] :
      absPdgId==15 ? 1.7776 :
      absPdgId==13 ? 0.10566 :
      absPdgId==11 ? 0.51100e-3 : 0.;

    TLorentzVector p4;
    p4.SetPtEtaPhiM(GenPart_pt[iGen],
		    GenPart_eta[iGen],
		    GenPart_phi[iGen],
		    mass);

    std::vector<unsigned int> daughterIndexes;
    bool isFinal=getDirectDaughterIndexes(daughterIndexes,(int)iGen,false);//store neutrinos for further use
    if(!isFinal) continue;

    /*MB: old code from times when statusFlags were not available
    bool fromHardProcessFinalState = false; //MB: note that it is actually not checked if the particles are final state ones!
    unsigned int firstCp_idx = findFirstCopy(iGen);
    if(GenPart_status[firstCp_idx]==3 || (GenPart_status[firstCp_idx]>20 && GenPart_status[firstCp_idx]<30)) fromHardProcessFinalState |= true;
    if(!fromHardProcessFinalState && GenPart_status[firstCp_idx]==22){
      std::vector<unsigned int> fsrIndexes;
      getDirectDaughterIndexes(fsrIndexes,(int)firstCp_idx);
      bool fsrBranching = false;
      for(unsigned iFsr=0; iFsr<fsrIndexes.size(); ++iFsr){
	if(GenPart_pdgId[fsrIndexes[iFsr]]==GenPart_pdgId[firstCp_idx] && GenPart_status[fsrIndexes[iFsr]]>50 && GenPart_status[fsrIndexes[iFsr]]<60) fsrBranching |= true;
      }
      if(!fsrBranching) fromHardProcessFinalState |= true;
    }
    unsigned int mother_idx = GenPart_genPartIdxMother[firstCp_idx]<0?firstCp_idx:GenPart_genPartIdxMother[firstCp_idx];
    unsigned int mother_absPdgId = std::abs(GenPart_pdgId[mother_idx]);
    bool mother_fromME = (GenPart_status[mother_idx]==3 || (GenPart_status[mother_idx]>20 && GenPart_status[mother_idx]<30)); // from matrix element
    if(mother_absPdgId==23 || //assume that it is prompt if from Z
       mother_absPdgId==24 || //W+-
       mother_absPdgId==25 || //h
       mother_absPdgId==35 || //H
       mother_absPdgId==36 || //A
       mother_absPdgId==32 || //Z'
       mother_absPdgId==33 || //Z''
       mother_absPdgId==34 || //W'
       mother_absPdgId==37 || //H+-
       mother_fromME || //or comes from ME
       false) {
      fromHardProcessFinalState |= true;
    }*/
    /*MB: "statusFlags().fromHardProcess()                      * 256   +" <<8 */
    bool fromHardProcessFinalState = ( (GenPart_statusFlags[iGen] & (1<<8)) == (1<<8) );
    bool isElectron = (absPdgId == 11);
    bool isMuon = (absPdgId == 13);
    bool isNeutrino = (absPdgId == 12 || absPdgId == 14 || absPdgId == 16);
    bool isTau = (absPdgId == 15);
    /*
    bool isDirectHardProcessTauDecayProduct = true;//FIXME(genFlags & (1<<10)) == (1<<10);
    if ( (fromHardProcessFinalState && (isMuon || isElectron || isNeutrino)) || isDirectHardProcessTauDecayProduct){
      genBosonP4 += p4;
    }
    ///This GenParticle list is missing pions, so we have
    ///to add hadronic tau, subtract neutral component
    if(absPdgId == 66615) genBosonP4 += p4;
    if(absPdgId == 77715) genBosonP4 -= p4;
    ///
    */
    //std::cout<<"\tpdgId="<<absPdgId<<", fromHardProcess:"<<fromHardProcessFinalState<<std::endl;
   //is the following correct? When undecayed particles are used (incl. taus) there should not be double counting
    if ( fromHardProcessFinalState && (isMuon || isElectron || isNeutrino || isTau) ){
      bosonP4 += p4;
      if(!isNeutrino)
	visBosonP4 += p4;
      if(isTau){
	//subtract p4 of tau neutrinos from visBosonP4
	for(unsigned int iDau=0; iDau<daughterIndexes.size(); ++iDau){
	  unsigned int absPdgIdDau = std::abs(GenPart_pdgId[daughterIndexes[iDau]]);
	  if(absPdgIdDau != 12 && absPdgIdDau != 14 && absPdgIdDau != 16 ) //neutrinos
	    continue;
	  TLorentzVector p4Dau;
	  p4Dau.SetPtEtaPhiM(GenPart_pt[daughterIndexes[iDau]],
			     GenPart_eta[daughterIndexes[iDau]],
			     GenPart_phi[daughterIndexes[iDau]],
			     0.);
	  visBosonP4 -= p4Dau;
	}
      }
    }
  }
  //std::cout<<"boson mass="<<bosonP4.M()<<std::endl;

 return ( bosonP4.M()>1E-3 || bosonP4.P()>1E-3 );
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::findTopP4(TLorentzVector &topP4, TLorentzVector &antiTopP4){

  topP4.SetXYZM(0,0,0,0);
  antiTopP4.SetXYZM(0,0,0,0);

  if(b_nGenPart==nullptr) return false;

  for(unsigned int iGen = 0; iGen < nGenPart; ++iGen) {
    if(std::abs(GenPart_pdgId[iGen])!=6)
      continue;

    //mass is stored only m>10GeV and photons m>1GeV (fine for top)
    TLorentzVector p4;
    p4.SetPtEtaPhiM(GenPart_pt[iGen],
		    GenPart_eta[iGen],
		    GenPart_phi[iGen],
		    GenPart_mass[iGen]);

    std::vector<unsigned int> daughterIndexes;
    bool isFinal=getDirectDaughterIndexes(daughterIndexes,(int)iGen);
    if(!isFinal) continue;

    if(GenPart_pdgId[iGen]==6) topP4 = p4;
    if(GenPart_pdgId[iGen]==-6) antiTopP4 = p4;
  }

  return (topP4.M()>1E-3 && antiTopP4.M()>1E-3);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
double HTauTauTreeFromNanoBase::getPtReweight(const TLorentzVector &genBosonP4, bool doSUSY){

  double weight = 1.0;

  //Z pt reweighting
  TH2F *hWeight = zptmass_histo;
  if(doSUSY) hWeight = zptmass_histo_SUSY;
  
  if(genBosonP4.M()>1E-3){
    double mass = genBosonP4.M();
    double pt = genBosonP4.Perp();    
    int massBin = hWeight->GetXaxis()->FindBin(mass);
    int ptBin = hWeight->GetYaxis()->FindBin(pt);
    weight = hWeight->GetBinContent(massBin,ptBin);
  }

  return weight;
  }
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::computeSvFit(HTTPair &aPair,
					   HTTAnalysis::sysEffects type){

  if(svFitAlgo_==nullptr) return;

  //Legs
  HTTParticle leg1 = aPair.getLeg1();
  double mass1;
  int decay1 = -1;
  classic_svFit::MeasuredTauLepton::kDecayType type1;
  if(std::abs(leg1.getPDGid())==11){
    mass1 = 0.51100e-3; //electron mass
    type1 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
  }
  else if(std::abs(leg1.getPDGid())==13){
    mass1 = 0.10566; //muon mass
    type1 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
  }
  else{//tau->hadrs.
    decay1 = leg1.getProperty(PropertyEnum::decayMode);
    mass1 = leg1.getP4().M();
    if(decay1==0)
      mass1 = 0.13957; //pi+/- mass
    type1 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
  }
  HTTParticle leg2 = aPair.getLeg2();
  double mass2;
  int decay2 = -1;
  classic_svFit::MeasuredTauLepton::kDecayType type2;
  if(std::abs(leg2.getPDGid())==11){
    mass2 = 0.51100e-3; //electron mass
    type2 = classic_svFit::MeasuredTauLepton::kTauToElecDecay;
  }
  else if(std::abs(leg2.getPDGid())==13){
    mass2 = 0.10566; //muon mass
    type2 = classic_svFit::MeasuredTauLepton::kTauToMuDecay;
  }
  else{//tau->hadrs.
    decay2 = leg2.getProperty(PropertyEnum::decayMode);
    mass2 = leg2.getP4().M();
    if(decay2==0)
      mass2 = 0.13957; //pi+/- mass
    type2 = classic_svFit::MeasuredTauLepton::kTauToHadDecay;
  }
  //Leptons for SvFit
  std::vector<classic_svFit::MeasuredTauLepton> measuredTauLeptons;
  measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(type1, leg1.getP4(type).Pt(), leg1.getP4(type).Eta(),
								leg1.getP4(type).Phi(), mass1, decay1) );
  measuredTauLeptons.push_back(classic_svFit::MeasuredTauLepton(type2, leg2.getP4(type).Pt(), leg2.getP4(type).Eta(),
								leg2.getP4(type).Phi(), mass2, decay2) );
  //MET
  TVector2 aMET = aPair.getMET(type);
  TMatrixD covMET(2, 2);
  covMET[0][0] = aPair.getMETMatrix().at(0);
  covMET[0][1] = aPair.getMETMatrix().at(1);
  covMET[1][0] = aPair.getMETMatrix().at(2);
  covMET[1][1] = aPair.getMETMatrix().at(3);

  if(covMET[0][0]==0 && covMET[1][0]==0 && covMET[0][1]==0 && covMET[1][1]==0) return; //singular covariance matrix

  TLorentzVector p4SVFit = aPair.getP4(HTTAnalysis::NOMINAL);
  TLorentzVector leg1P4Nominal = leg1.getP4(HTTAnalysis::NOMINAL);
  TLorentzVector leg2P4Nominal = leg2.getP4(HTTAnalysis::NOMINAL);

  
  if(type==HTTAnalysis::NOMINAL ||
     leg1.getP4(type)!=leg1P4Nominal ||
     leg2.getP4(type)!=leg2P4Nominal){
       p4SVFit = runSVFitAlgo(measuredTauLeptons, aMET, covMET);
     }
  aPair.setP4(p4SVFit,type);
  aPair.setLeg1P4(p4Leg1SVFit,type);
  aPair.setLeg2P4(p4Leg2SVFit,type);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
TLorentzVector HTauTauTreeFromNanoBase::runSVFitAlgo(const std::vector<classic_svFit::MeasuredTauLepton> & measuredTauLeptons,
						     const TVector2 &aMET, const TMatrixD &covMET){

  TLorentzVector p4SVFit;
  if(measuredTauLeptons.size()!=2 || svFitAlgo_==nullptr) return p4SVFit;

  //set logM regularization term which is final state dependent
  double kappa = 4;
  if(measuredTauLeptons[0].type()==classic_svFit::MeasuredTauLepton::kTauToElecDecay || measuredTauLeptons[0].type()==classic_svFit::MeasuredTauLepton::kTauToMuDecay) { //1st tau is lepton
    if(measuredTauLeptons[1].type()==classic_svFit::MeasuredTauLepton::kTauToElecDecay || measuredTauLeptons[1].type()==classic_svFit::MeasuredTauLepton::kTauToMuDecay)
      kappa = 3; //ll decay
    else
      kappa = 4; //lt decay
  }
  else {//1st tau is hadron
    if(measuredTauLeptons[1].type()==classic_svFit::MeasuredTauLepton::kTauToElecDecay || measuredTauLeptons[1].type()==classic_svFit::MeasuredTauLepton::kTauToMuDecay)
      kappa = 4; //ltau decay
    else
      kappa = 5; //tt decay
  }    
  svFitAlgo_->addLogM_fixed(true, kappa);
  svFitAlgo_->integrate(measuredTauLeptons, aMET.X(), aMET.Y(), covMET);
  if(svFitAlgo_->isValidSolution() ){//Get solution

    classic_svFit::DiTauSystemHistogramAdapter* aHistogramAdapter = static_cast< classic_svFit::DiTauSystemHistogramAdapter*>(svFitAlgo_->getHistogramAdapter());
    p4SVFit.SetPtEtaPhiM(aHistogramAdapter->getPt(),
			 aHistogramAdapter->getEta(),
			 aHistogramAdapter->getPhi(),
			 aHistogramAdapter->getMass());

    /*not available with official version
    double tauMass = 1.77686; //GeV, PDG value
    p4Leg1SVFit.SetPtEtaPhiM(aHistogramAdapter->getLeg1Pt(),
			     aHistogramAdapter->getLeg1Eta(),
			     aHistogramAdapter->getLeg1Phi(),
			     tauMass);

    p4Leg2SVFit.SetPtEtaPhiM(aHistogramAdapter->getLeg2Pt(),
			     aHistogramAdapter->getLeg2Eta(),
			     aHistogramAdapter->getLeg2Phi(),
			     tauMass);
    */
  }
  else{
    p4SVFit.SetPtEtaPhiM(0,0,0,0);
    p4Leg1SVFit.SetPtEtaPhiM(0,0,0,0);
    p4Leg2SVFit.SetPtEtaPhiM(0,0,0,0);
  }

  return p4SVFit;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
void HTauTauTreeFromNanoBase::applyMetRecoilCorrections(){

  // Do nothing if there is not best pair or recoilCorrector is not initialized
  if( recoilCorrector_==nullptr || httPairCollection.empty() )
    return;
  TVector2 theUncorrMEt;
  float corrMEtPx, corrMEtPy;
  int nJets = httJetCollection.size();
  if(httEvent->getDecayModeBoson()>=10)//W, add jet for rake tau
    nJets++;
  TLorentzVector genBosonP4 = httEvent->getGenBosonP4();
  TLorentzVector genBosonVisP4 = httEvent->getGenBosonP4(true);
  /* Do not correct Met in the event, keep it as it is
  // Correct Met in the event
  theUncorrMEt = httEvent->getMET();
  recoilCorrector_->CorrectByMeanResolution(
  //recoilCorrector_->Correct( //Quantile correction works better for MVA MET
      theUncorrMEt.Px(),
      theUncorrMEt.Py(),
      genBosonP4.Px(),
      genBosonP4.Py(),
      genBosonVisP4.Px(),
      genBosonVisP4.Py(),
      nJets,
      corrMEtPx,
      corrMEtPy
  );
  httEvent->setMET( TVector2(corrMEtPx,corrMEtPy) );
  */
  // Correct Met in pairs (a priori it can be by-pair Met)
  for(unsigned int iPair=0; iPair<httPairCollection.size(); ++iPair){
    //theUncorrMEt = httEvent->getMET();
    theUncorrMEt = httPairCollection[iPair].getMET();//TES corrected, fine??
    recoilCorrector_->CorrectByMeanResolution(
    //recoilCorrector_->Correct( //Quantile correction works better for MVA MET
        theUncorrMEt.Px(),
        theUncorrMEt.Py(),
        genBosonP4.Px(),
        genBosonP4.Py(),
        genBosonVisP4.Px(),
        genBosonVisP4.Py(),
        nJets,
        corrMEtPx,
        corrMEtPy
    );
    //Remove TES correction to not have it twice
    if( (std::abs(httPairCollection[iPair].getLeg1().getPDGid())==15 && httPairCollection[iPair].getLeg1().getProperty(PropertyEnum::mc_match)==5) ||
        (std::abs(httPairCollection[iPair].getLeg2().getPDGid())==15 && httPairCollection[iPair].getLeg2().getProperty(PropertyEnum::mc_match)==5) ){
      corrMEtPx-=httPairCollection[iPair].getLeg1().getP4(HTTAnalysis::DUMMY_SYS).X(); //uncor
      corrMEtPx-=httPairCollection[iPair].getLeg2().getP4(HTTAnalysis::DUMMY_SYS).X(); //uncor
      corrMEtPx+=httPairCollection[iPair].getLeg1().getP4(HTTAnalysis::NOMINAL).X();
      corrMEtPx+=httPairCollection[iPair].getLeg2().getP4(HTTAnalysis::NOMINAL).X();
      corrMEtPy-=httPairCollection[iPair].getLeg1().getP4(HTTAnalysis::DUMMY_SYS).Y(); //uncor
      corrMEtPy-=httPairCollection[iPair].getLeg2().getP4(HTTAnalysis::DUMMY_SYS).Y(); //uncor
      corrMEtPy+=httPairCollection[iPair].getLeg1().getP4(HTTAnalysis::NOMINAL).Y();
      corrMEtPy+=httPairCollection[iPair].getLeg2().getP4(HTTAnalysis::NOMINAL).Y();
    }
    httPairCollection[iPair].setMET( TVector2(corrMEtPx,corrMEtPy) );

    //recompute mT's using consistently TES corrected MEt and Pt
    double mTLeg1 = TMath::Sqrt(2.*httPairCollection[iPair].getLeg1().getP4().Pt()*httPairCollection[iPair].getMET().Mod()*(1.-TMath::Cos(httPairCollection[iPair].getLeg1().getP4().Phi()-httPairCollection[iPair].getMET().Phi())));
    double mTLeg2 = TMath::Sqrt(2.*httPairCollection[iPair].getLeg2().getP4().Pt()*httPairCollection[iPair].getMET().Mod()*(1.-TMath::Cos(httPairCollection[iPair].getLeg2().getP4().Phi()-httPairCollection[iPair].getMET().Phi())));

    httPairCollection[iPair].setMTLeg1(mTLeg1);
    httPairCollection[iPair].setMTLeg2(mTLeg2);
  }

  return;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::compareLeptons(const HTTParticle& i, const HTTParticle& j){
  unsigned int i_type=2, j_type=2;
  if(std::abs(i.getPDGid())==11)//ele
    i_type=1;
  else if(std::abs(i.getPDGid())==13)//mu
    i_type=0;
  if(std::abs(j.getPDGid())==11)//ele
    j_type=1;
  else if(std::abs(j.getPDGid())==13)//mu
    j_type=0;
  if(i_type > j_type) return false;
  if(i_type == j_type && i.getP4().Pt() < j.getP4().Pt() ) return false;

  return true;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::comparePairs(const HTTPair& i, const HTTPair& j){
  double i_iso=999, j_iso=999;
  unsigned int i_type=2, j_type=2;

  //step 0.5, leg 1 type
  i_type = std::abs(i.getLeg1().getPDGid())==15 ? 2: std::abs(i.getLeg1().getPDGid())==11 ? 1: 0;
  j_type = std::abs(j.getLeg1().getPDGid())==15 ? 2: std::abs(j.getLeg1().getPDGid())==11 ? 1: 0;
  if (i_type<j_type) return true;
  else if(i_type>j_type) return false;

  //step 1, leg 1 ISO
  i_iso = std::abs(i.getLeg1().getPDGid())==15 ? -i.getLeg1().getProperty(PropertyEnum::rawMVAoldDM) : std::abs(i.getLeg1().getPDGid())==11 ? i.getLeg1().getProperty(PropertyEnum::pfRelIso03_all): i.getLeg1().getProperty(PropertyEnum::pfRelIso04_all);
  if(i_iso<-1) i_iso=999; //something went wrong
  j_iso = std::abs(j.getLeg1().getPDGid())==15 ? -j.getLeg1().getProperty(PropertyEnum::rawMVAoldDM) : std::abs(j.getLeg1().getPDGid())==11 ? j.getLeg1().getProperty(PropertyEnum::pfRelIso03_all): j.getLeg1().getProperty(PropertyEnum::pfRelIso04_all);
  if(j_iso<-1) j_iso=999; //something went wrong
  if (i_iso<j_iso) return true;
  else if(i_iso>j_iso) return false;

  //step 2, leg 1 Pt
  if(i.getLeg1().getP4().Pt()>j.getLeg1().getP4().Pt()) return true;
  else if(i.getLeg1().getP4().Pt()<j.getLeg1().getP4().Pt()) return false;

  //step 2.5, leg 2 type
  i_type = std::abs(i.getLeg2().getPDGid())==15 ? 2: std::abs(i.getLeg2().getPDGid())==11 ? 1: 0;
  j_type = std::abs(j.getLeg2().getPDGid())==15 ? 2: std::abs(j.getLeg2().getPDGid())==11 ? 1: 0;
  if (i_type<j_type) return true;
  else if(i_type>j_type) return false;

  //step 3, leg 2 ISO
  i_iso = std::abs(i.getLeg2().getPDGid())==15 ? -i.getLeg2().getProperty(PropertyEnum::rawMVAoldDM) : std::abs(i.getLeg2().getPDGid())==11 ? i.getLeg2().getProperty(PropertyEnum::pfRelIso03_all): i.getLeg2().getProperty(PropertyEnum::pfRelIso04_all);
  if(i_iso<-1) i_iso=999; //something went wrong
  j_iso = std::abs(j.getLeg2().getPDGid())==15 ? -j.getLeg2().getProperty(PropertyEnum::rawMVAoldDM) : std::abs(j.getLeg2().getPDGid())==11 ? j.getLeg2().getProperty(PropertyEnum::pfRelIso03_all): j.getLeg2().getProperty(PropertyEnum::pfRelIso04_all);
  if(j_iso<-1) j_iso=999; //something went wrong
  if (i_iso<j_iso) return true;
  else if(i_iso>j_iso) return false;

  //step 4, leg 2 Pt
  if(i.getLeg2().getP4().Pt()>j.getLeg2().getP4().Pt()) return true;

  return false;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
// int HTauTauTreeFromNanoBase::isGenPartDaughterPdgId(int index, unsigned int absPdgId){
//   if(index>nGenPart || index<0 || GenPart_genPartIdxMother[index]<0) return -1;
//   if(std::abs(GenPart_pdgId[GenPart_genPartIdxMother[index]])==absPdgId)
//     return GenPart_genPartIdxMother[index];
//   else
//     return isGenPartDaughterPdgId(GenPart_genPartIdxMother[index],absPdgId);
// }
/////////////////////////////////////////////////
/////////////////////////////////////////////////
// bool HTauTauTreeFromNanoBase::isGenPartDaughterIdx(int index, int mother){
//   if(index>nGenPart || index<0 || GenPart_genPartIdxMother[index]<0) return false;
//   if(GenPart_genPartIdxMother[index]==mother)
//     return true;
//   else
//     return isGenPartDaughterIdx(GenPart_genPartIdxMother[index],mother);
// }
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::getDirectDaughterIndexes(std::vector<unsigned int> &indexes, unsigned int motherIndex, bool ignoreNeutrinos){  
  indexes.clear();
  bool isFinal=true;
  for(unsigned int iDau=0;iDau<nGenPart;++iDau){
    if(GenPart_genPartIdxMother[iDau]==(int)motherIndex){
      unsigned int aPdgId = std::abs(GenPart_pdgId[iDau]);
      if(aPdgId==std::abs(GenPart_pdgId[motherIndex])){
	isFinal=false;
	//break;
      }
      if((aPdgId==12 || aPdgId==14 || aPdgId==16)&&ignoreNeutrinos)
	continue;
      indexes.push_back(iDau);
    }
    //if(!isFinal) break;
  }
  return isFinal;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
unsigned int HTauTauTreeFromNanoBase::findFinalCopy(unsigned int index){  
  std::vector<unsigned int> daughterIndexes;
  if(getDirectDaughterIndexes(daughterIndexes,index))
     return index;
  unsigned int newIndex=index;
  for(unsigned int iDau=0;iDau<daughterIndexes.size();++iDau){
    int pdgId = GenPart_pdgId[daughterIndexes[iDau]];
    if(pdgId==GenPart_pdgId[index]){
	newIndex=daughterIndexes[iDau];
	break;
    }
  }
  //check if found copy is decaying (should not happen for unstable particles)
  getDirectDaughterIndexes(daughterIndexes,newIndex,false);
  if(daughterIndexes.empty()) return index;

  return findFinalCopy(newIndex);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
unsigned int HTauTauTreeFromNanoBase::findFirstCopy(unsigned int index){
  int mo_idx = GenPart_genPartIdxMother[index];
  if(mo_idx<0) return index;

  int mo_pdgId = GenPart_pdgId[mo_idx];
  if(mo_pdgId!=GenPart_pdgId[index]) return index;
  
  return findFirstCopy((unsigned int)mo_idx);
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
int HTauTauTreeFromNanoBase::genTauDecayMode(std::vector<unsigned int> &daughterIndexes){

  if(daughterIndexes.empty()) return -99;
 
  int numElectrons      = 0;
  int numMuons          = 0;
  int numChargedPions   = 0;
  int numNeutralPions   = 0;
  int numPhotons        = 0;
  int numNeutrinos      = 0;
  int numOtherParticles = 0;
  
  for(unsigned int idx=0; idx<daughterIndexes.size(); ++idx) {
    int pdg_id = std::abs(GenPart_pdgId[daughterIndexes[idx]]);
    if(pdg_id == 11) numElectrons++;
    else if(pdg_id == 13) numMuons++;
    else if(pdg_id == 211 || pdg_id == 321 ) numChargedPions++; //Count both pi+ and K+
    else if(pdg_id == 111 || pdg_id == 130 || pdg_id == 310 || pdg_id == 311 ) numNeutralPions++; //Count both pi0 and K0 {_L/S}
    else if(pdg_id == 12 || 
	    pdg_id == 14 || 
	    pdg_id == 16) {
      numNeutrinos++;
    }
    else if(pdg_id == 22) numPhotons++;
    else {
      numOtherParticles++;
    }
  }
  if(numElectrons>1){//sometimes there are gamma->ee conversions 
    numPhotons += numElectrons/2;
    numElectrons -= 2*(numElectrons/2);
  }
  //"convert" photons to piZeros
  numNeutralPions += numPhotons/2;
  numPhotons -= 2*(numPhotons/2);
  int tauDecayMode = HTTAnalysis::tauDecayOther;
  if( numOtherParticles == 0 ){
    if( numElectrons == 1 ){
      //--- tau decays into electrons
      tauDecayMode = HTTAnalysis::tauDecaysElectron;
    } else if( numMuons == 1 ){
      //--- tau decays into muons
      tauDecayMode = HTTAnalysis::tauDecayMuon;
    } else {
      //--- hadronic tau decays
      switch ( numChargedPions ){
      case 1 :
	if( numPhotons != 0 ){
	  tauDecayMode =  HTTAnalysis::tauDecayOther;
	  break;
	}
	switch ( numNeutralPions ){
	case 0:
	  tauDecayMode = HTTAnalysis::tauDecay1ChargedPion0PiZero;
	  break;
	case 1:
	  tauDecayMode = HTTAnalysis::tauDecay1ChargedPion1PiZero;
	  break;
	case 2:
	  tauDecayMode = HTTAnalysis::tauDecay1ChargedPion2PiZero;
	  break;
	case 3:
	  tauDecayMode = HTTAnalysis::tauDecay1ChargedPion3PiZero;
	  break;
	case 4:
	  tauDecayMode = HTTAnalysis::tauDecay1ChargedPion4PiZero;
	  break;
	default:
	  tauDecayMode = HTTAnalysis::tauDecayOther;
	  break;
	}
	break;
      case 3 : 
	if( numPhotons != 0 ){
	  tauDecayMode = HTTAnalysis::tauDecayOther;
	  break;
	}
	switch ( numNeutralPions ){
	case 0 : 
	  tauDecayMode = HTTAnalysis::tauDecay3ChargedPion0PiZero;
	  break;
	case 1 : 
	  tauDecayMode = HTTAnalysis::tauDecay3ChargedPion1PiZero;
	  break;
	case 2 : 
	  tauDecayMode = HTTAnalysis::tauDecay3ChargedPion2PiZero;
	  break;
	case 3 : 
	  tauDecayMode = HTTAnalysis::tauDecay3ChargedPion3PiZero;
	  break;
	case 4 : 
	  tauDecayMode = HTTAnalysis::tauDecay3ChargedPion4PiZero;
	  break;
	default:
	  tauDecayMode = HTTAnalysis::tauDecayOther;
	  break;
	}
	break;
      }
    }
  }
  return tauDecayMode;

}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauTauTreeFromNanoBase::eventInJson(){
  //If the jsonVec is empty, then no JSON file was provided so all events should be accepted
  if( jsonVector.empty() )
    return true;

  edm::LuminosityBlockID lumiID(run, luminosityBlock);
  for(std::vector<edm::LuminosityBlockRange>::const_iterator itr = jsonVector.begin();
      itr != jsonVector.end(); ++itr){
    if( edm::contains(*itr,lumiID) ) return true;
  }

  return false;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
