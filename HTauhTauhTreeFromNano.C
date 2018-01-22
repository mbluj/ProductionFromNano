#define HTauhTauhTreeFromNano_cxx
#include "HTauhTauhTreeFromNano.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>

#include <iostream>
#include <fstream>

/////////////////////////////////////////////////
/////////////////////////////////////////////////
bool HTauhTauhTreeFromNano::pairSelection(unsigned int iPair){

  ///Baseline+post sync selection as on
  ///https://twiki.cern.ch/twiki/bin/view/CMS/HiggsToTauTauWorking2015#Baseline_mu_tau_h_AN1
  ///Indexes for multiplexed ID variables taken from   LLRHiggsTauTau/NtupleProducer/plugins/
  ///HTauTauNtuplizer.cc, MuFiller.cc, TauFiller.cc, EleFiller.cc

  if(httPairs_.empty()) return false;

  int pdgIdLeg1 = httPairs_[iPair].getLeg1().getPDGid();
  int pdgIdLeg2 = httPairs_[iPair].getLeg2().getPDGid();
  if( std::abs(pdgIdLeg1)!=15 || std::abs(pdgIdLeg2)!=15 ) return 0;

  int tauIDmask = 0, tauIDmaskMedium = 0 , tauIDmaskLoose = 0;
  for(unsigned int iBit=0;iBit<HTTEvent::ntauIds;iBit++){
    if(HTTEvent::tauIDStrings[iBit]=="byTightIsolationMVArun2v1DBoldDMwLT")
      tauIDmask |= (1<<iBit);
    if(HTTEvent::tauIDStrings[iBit]=="byMediumIsolationMVArun2v1DBoldDMwLT")
      tauIDmaskMedium |= (1<<iBit);
    if(HTTEvent::tauIDStrings[iBit]=="byLooseIsolationMVArun2v1DBoldDMwLT")
      tauIDmaskLoose |= (1<<iBit);
    if(HTTEvent::tauIDStrings[iBit]=="againstMuonLoose3") {
      tauIDmask |= (1<<iBit);
      tauIDmaskMedium |= (1<<iBit);
      tauIDmaskLoose |= (1<<iBit);
    }
    if(HTTEvent::tauIDStrings[iBit]=="againstElectronVLooseMVA6") {
      tauIDmask |= (1<<iBit);
      tauIDmaskMedium |= (1<<iBit);
      tauIDmaskLoose |= (1<<iBit);
    }
  }
  unsigned int indexLeg1 = httPairs_[iPair].getIndexLeg1();
  unsigned int indexLeg2 = httPairs_[iPair].getIndexLeg2();
  //MB sort taus within the pair
  double pt_1 = httLeptonCollection[indexLeg1].getP4().Pt();
  double pt_2 = httLeptonCollection[indexLeg2].getP4().Pt();
  if(pt_2>pt_1){//tau with higher-Pt first
    unsigned int indexLegTmp = indexLeg1;
    indexLeg1 = indexLeg2;
    indexLeg2 = indexLegTmp;
  }
  TLorentzVector tau1P4 = httLeptonCollection[indexLeg1].getP4();
  TLorentzVector tau2P4 = httLeptonCollection[indexLeg2].getP4();

  int tau1ID = (int)httLeptonCollection[indexLeg1].getProperty(PropertyEnum::idAntiMu);
  tau1ID += (int)std::pow(2,HTTEvent::againstEIdOffset)*(int)httLeptonCollection[indexLeg1].getProperty(PropertyEnum::idAntiEle);
  tau1ID += (int)std::pow(2,HTTEvent::mvaIsoIdOffset)*(int)httLeptonCollection[indexLeg1].getProperty(PropertyEnum::idMVAoldDM);
  int tau2ID = (int)httLeptonCollection[indexLeg2].getProperty(PropertyEnum::idAntiMu);
  tau2ID += (int)std::pow(2,HTTEvent::againstEIdOffset)*(int)httLeptonCollection[indexLeg2].getProperty(PropertyEnum::idAntiEle);
  tau2ID += (int)std::pow(2,HTTEvent::mvaIsoIdOffset)*(int)httLeptonCollection[indexLeg2].getProperty(PropertyEnum::idMVAoldDM);

  bool tauBaselineSelection1 = tau1P4.Pt()>50 && std::abs(tau1P4.Eta())<2.1 &&
                               httLeptonCollection[indexLeg1].getProperty(PropertyEnum::idDecayMode)>0 &&
                               std::abs(httLeptonCollection[indexLeg1].getProperty(PropertyEnum::dz))<0.2 &&
                               (int)std::abs(httLeptonCollection[indexLeg1].getProperty(PropertyEnum::charge))==1;
  bool tauBaselineSelection2 = tau2P4.Pt()>40 && std::abs(tau2P4.Eta())<2.1 &&
                               httLeptonCollection[indexLeg2].getProperty(PropertyEnum::idDecayMode)>0 &&
                               std::abs(httLeptonCollection[indexLeg2].getProperty(PropertyEnum::dz))<0.2 &&
                               (int)std::abs(httLeptonCollection[indexLeg2].getProperty(PropertyEnum::charge))==1;

  bool baselinePair = tau1P4.DeltaR(tau2P4) > 0.5;
  bool postSynchTau1 = (tau1ID & tauIDmask) == tauIDmask;
  bool postSynchTau2 = (tau2ID & tauIDmask) == tauIDmask;
  ///
  bool postSynchLooseTau1 = (tau1ID & tauIDmaskLoose) == tauIDmaskLoose;
  bool postSynchLooseTau2 = (tau2ID & tauIDmaskLoose) == tauIDmaskLoose;
  bool postSynchMediumTau1 = (tau1ID & tauIDmaskMedium) == tauIDmaskMedium;
  bool postSynchMediumTau2 = (tau2ID & tauIDmaskMedium) == tauIDmaskMedium;

  httEvent->clearSelectionWord();
  httEvent->setSelectionBit(SelectionBitsEnum::muonBaselineSelection,tauBaselineSelection1);
  httEvent->setSelectionBit(SelectionBitsEnum::tauBaselineSelection,tauBaselineSelection2);
  httEvent->setSelectionBit(SelectionBitsEnum::baselinePair,baselinePair);
  httEvent->setSelectionBit(SelectionBitsEnum::postSynchMuon,postSynchTau1);
  httEvent->setSelectionBit(SelectionBitsEnum::postSynchTau,postSynchTau2);
  httEvent->setSelectionBit(SelectionBitsEnum::extraMuonVeto,thirdLeptonVeto(indexLeg1,indexLeg2,13));
  httEvent->setSelectionBit(SelectionBitsEnum::extraElectronVeto,thirdLeptonVeto(indexLeg1,indexLeg2,11));

  return tauBaselineSelection1 && tauBaselineSelection2 && baselinePair
    //&& ( (postSynchLooseTau1 && postSynchMediumTau2) || (postSynchLooseTau2 && postSynchMediumTau1) )
    //&& !thirdLeptonVeto(indexLeg1,indexLeg2,13)
    //&& !thirdLeptonVeto(indexLeg1,indexLeg2,11)
    && true;
}
/////////////////////////////////////////////////
/////////////////////////////////////////////////
unsigned int HTauhTauhTreeFromNano::bestPair(std::vector<unsigned int> &pairIndexes){

  unsigned int bestIndex = 9999;
  ///Pairs are already sorted during the ntuple creation?
  double iso_1=std::numeric_limits<double>::infinity(), iso_2=std::numeric_limits<double>::infinity(), pt_1=-1, pt_2=-1;
  if(pairIndexes.size()) {
    //return pairIndexes[0];//MB
    for(unsigned int ii=0;ii<2*pairIndexes.size();++ii){
      unsigned int i=(ii<pairIndexes.size()?ii:ii-pairIndexes.size());
      unsigned int iPair = pairIndexes[i];
      unsigned int indexLeg1 = httPairs_[iPair].getIndexLeg1();
      unsigned int indexLeg2 = httPairs_[iPair].getIndexLeg2();
      if(ii>=pairIndexes.size()){//invert legs
	indexLeg1 = httPairs_[iPair].getIndexLeg2();
	indexLeg2 = httPairs_[iPair].getIndexLeg1();
      }
      double pt_1_i = httLeptonCollection[indexLeg1].getP4().Pt();
      double pt_2_i = httLeptonCollection[indexLeg2].getP4().Pt();
      //MB: More isolated for MVAIso means higher value so inverted here to keep standard convention in comparison
      double iso_1_i = -httLeptonCollection[indexLeg1].getProperty(PropertyEnum::rawMVAoldDM);
      double iso_2_i = -httLeptonCollection[indexLeg2].getProperty(PropertyEnum::rawMVAoldDM);

      if(iso_1_i>iso_1) continue;
      if(iso_1_i==iso_1 && pt_1_i<pt_1) continue;
      if(iso_2_i>iso_2) continue;
      if(iso_2_i==iso_2 && pt_2_i<pt_2) continue;
      bestIndex = iPair;
      iso_1 = iso_1_i;
      iso_2 = iso_2_i;
      pt_1 = pt_1_i;
      pt_2 = pt_2_i;
    }
  }
  /*
  if(pairIndexes.size() && bestIndex!=pairIndexes[0]){
    unsigned int iPair = pairIndexes[0];
    unsigned int indexLeg1 = httPairs_[iPair].getIndexLeg1();
    unsigned int indexLeg2 = httPairs_[iPair].getIndexLeg2();
    double pt_1_i = httLeptonCollection[indexLeg1].getP4().Pt();
    double pt_2_i = httLeptonCollection[indexLeg2].getP4().Pt();
    std::cout<<"Pair sorting: "<<std::endl
	     <<"best index = "<<bestIndex<<", index[0] = "<<pairIndexes[0]<<std::endl
	     <<"\tiso1[best]="<<-iso_1<<", iso1[0]="<<httLeptonCollection[indexLeg1].getProperty(PropertyEnum::rawMVAoldDM)<<std::endl
	     <<"\tpt1[best]="<<pt_1<<", pt1[0]="<<pt_1_i<<std::endl
	     <<"\tiso2[best]="<<-iso_2<<", iso1[0]="<<httLeptonCollection[indexLeg2].getProperty(PropertyEnum::rawMVAoldDM)<<std::endl
	     <<"\tpt2[best]="<<pt_2<<", pt1[0]="<<pt_2_i<<std::endl;
  }
  */

  return bestIndex;
};
/////////////////////////////////////////////////
/////////////////////////////////////////////////
