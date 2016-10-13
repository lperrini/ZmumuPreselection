#ifndef OfflineHelper_h
#define OfflineHelper_h

/** \class OfflineProducerHelper
 *
 *  Collection of functions to analyze HTauTau primary trees, without CMSSW
 *  dependencies.
 *
 *  $Date: 2012/06/06 00:27:43 $
 *  $Revision: 1.3 $
 *  \author G. Ortona - LLR
 */

#include <TString.h>
#include "TMath.h"
#include "bigTree.h"
#include "HTauTauTree.h"
#include "TLorentzVector.h"
#include <iostream>
#include <vector>
#include <utility>

using namespace std;

using namespace std;

class OfflineProducerHelper {
 public:
  enum particleType {
  MUON = 0,
  ELECTRON = 1,
  TAU =2
  };
  float m_MVAEleIDCuts[2][2][3] ;

  enum pairType {
    MuHad  = 0,
    EHad   = 1,
    HadHad = 2,
    MuMu   = 3,
    EE     = 4,
    EMu    = 5,
    EEPrompt = 6,
    MuMuPrompt = 7,
    Other  = 8
  };

  OfflineProducerHelper();
  OfflineProducerHelper(TH1F* hCounter, TH1F *htauids);
  OfflineProducerHelper(TH1F* hCounter);

  int FindTriggerNumber(TString triggername);
  bool IsTriggerFired(int triggerbit, int triggerNumber);
  bool IsTriggerFired(int triggerbit, TString triggerName){return IsTriggerFired(triggerbit, FindTriggerNumber(triggerName));}
  int printFiredPaths(int);
  bool isMuon(int type){if(type == MUON)return true; else return false;}
  bool isElectron(int type){if(type == ELECTRON)return true; else return false;}
  bool isTau(int type){if(type == TAU)return true; else return false;}
  int getPairType (int type1, int type2);
  bool checkBit (int word, int bitpos);
  int getTAUidNumber(TString tauIDname);

  bool pairPassBaseline (bigTree* tree, int iPair, TString whatApply = "All");
  bool eleBaseline (bigTree* tree, int iDau, float ptMin, float relIso,  int MVAIDflag = 0, TString whatApply = "All"); // return true if leptons passes the baseline selections
  bool muBaseline (bigTree* tree, int iDau, float ptMin, float etaMax, float relIso, TString whatApply = "All");
  bool tauBaseline (bigTree* tree, int iDau, float ptMin, float etaMax, int againstEleWP, int againstMuWP, float isoRaw3Hits, TString whatApply = "All");
  bool tightEleMVAID (float BDT, float fSCeta);

  int getMothPairType (bigTree* tree, int iMoth);


  bool EleMVAID (float BDT, float eta, float pT, int strength) ;

  TLorentzVector buildDauP4 (bigTree* tree, int iDau); // build daughter 4vector
  TLorentzVector buildMothP4 (bigTree* tree, int iMoth); // build pair 4vector
  TLorentzVector buildGenP4 (bigTree* tree, int iGen); // build pair 4 vector
  bool getBestJets (bigTree* tree, int& jet1, int& jet2, int strategy); //select jets, possibly two b jets, returns true if found, else false
  int getBestPair (bigTree* tree, std::vector<int>& pairIdxs, TString strategy= "OSMaxPt"); // from a vector of indexes to pairs in the evemt retutn indexto the one chose by a stategy
  int getBestPair (bigTree* tree, TString strategy = "OSMaxPt"); // calls theprevious on the whole pair collection of the event
  int getPairByIndexes (bigTree* tree, int dau1, int dau2); // knowing thesons, get the pair formed


  int MCHiggsTauTauDecayMode (bigTree* tree); // find the MC decay of theHiggs to tau in the eventa
  bool getHardTauFinalVisGenProducts (bigTree* tree, int& ind1, int& ind2); //find hard scatter tau decay products and store their indices; return falsa^e ifproblems, else true
  bool drMatchGenReco (bigTree* tree, int iGen, int iReco, float dRcone =0.5);
  int getRecoMatchedToGen (bigTree* tree, int iGen, bool checkId = true, bool checkCharge = false, float dRcone = 0.5);


  ~OfflineProducerHelper(){}

 private:
  vector<TString> triggerlist;
  vector<TString> tauidlist;

};


OfflineProducerHelper::OfflineProducerHelper(){
  const int nTriggers=20;
  TString tmptrigger[nTriggers]={
	  "IsoMu24_eta2p1",
	  "IsoMu27",
	  "IsoMu17_eta2p1",
	  "IsoMu18",
	  "IsoMu22",
	  "Ele32_eta2p1_WP75_Gsf",
	  "Ele22_eta2p1_WP75_Gsf",
	  "Ele32_eta2p1_WPTight_Gsf",
	  "Ele23_WPLoose_Gsf",
	  "IsoMu17_eta2p1_LooseIsoPFTau20",
	  "Ele22_eta2p1_WP75_Gsf_LooseIsoPFTau20",
	  "Ele22_eta2p1_WPLoose_Gsf_LooseIsoPFTau20",
	  "DoubleMediumIsoPFTau35_Trk1_eta2p1_Reg",
	  "DoubleMediumIsoPFTau40_Trk1_eta2p1_Reg",
	  "Mu17_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
	  "Mu8_TrkIsoVVL_Ele17_CaloIdL_TrackIdL_IsoVL",
	  "Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL",
	  "Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL",
  };
  for(int i=0;i<nTriggers;i++){
    //triggerlist[i]=tmptrigger[i];
    tmptrigger[i].Prepend("HLT_");
    tmptrigger[i].Append("_v");
    triggerlist.push_back(tmptrigger[i]);
  }
  tauidlist.push_back("");

  // MVA ele ID from here:
  //  https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Non_triggering_electron_MVA
  // 80%
  m_MVAEleIDCuts[0][0][0] = -0.253 ; // barrel (eta<0.8) pt 5-10 GeV      
  m_MVAEleIDCuts[0][0][1] =  0.081 ; // barrel (eta>0.8) pt 5-10 GeV      
  m_MVAEleIDCuts[0][0][2] = -0.081 ; // endcap pt 5-10 GeV                
  m_MVAEleIDCuts[0][1][0] =  0.965 ; // barrel (eta<0.8) pt above 10 GeV  
  m_MVAEleIDCuts[0][1][1] =  0.917 ; // barrel (eta>0.8) pt above 10 GeV  
  m_MVAEleIDCuts[0][1][2] =  0.683 ; // endcap pt above 10 GeV            

  // 90%
  m_MVAEleIDCuts[1][0][0] = -0.483 ; // barrel (eta<0.8) pt 5-10 GeV     
  m_MVAEleIDCuts[1][0][1] = -0.267 ; // barrel (eta>0.8) pt 5-10 GeV     
  m_MVAEleIDCuts[1][0][2] = -0.323 ; // endcap pt 5-10 GeV               
  m_MVAEleIDCuts[1][1][0] = 0.933  ; // barrel (eta<0.8) pt above 10 GeV 
  m_MVAEleIDCuts[1][1][1] = 0.825  ; // barrel (eta>0.8) pt above 10 GeV 
  m_MVAEleIDCuts[1][1][2] = 0.337  ; // endcap pt above 10 GeV           
}

OfflineProducerHelper::OfflineProducerHelper(TH1F* hCounter, TH1F* hTauIDs){

  for(int itr=1;itr<=hCounter->GetNbinsX();itr++){
    TString binlabel = hCounter->GetXaxis()->GetBinLabel(itr);
    if(binlabel.BeginsWith("HLT"))triggerlist.push_back(hCounter->GetXaxis()->GetBinLabel(itr));
  }

  for(int itr=1;itr<=hTauIDs->GetNbinsX();itr++){
    TString binlabel = hTauIDs->GetXaxis()->GetBinLabel(itr);
    tauidlist.push_back(hTauIDs->GetXaxis()->GetBinLabel(itr));
  }

  // MVA ele ID from here:
  //   //  https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Non_triggering_electron_MVA
  //     // 80%

  m_MVAEleIDCuts[0][0][0] = -0.253 ; // barrel (eta<0.8) pt 5-10 GeV      
  m_MVAEleIDCuts[0][0][1] =  0.081 ; // barrel (eta>0.8) pt 5-10 GeV      
  m_MVAEleIDCuts[0][0][2] = -0.081 ; // endcap pt 5-10 GeV                
  m_MVAEleIDCuts[0][1][0] =  0.965 ; // barrel (eta<0.8) pt above 10 GeV  
  m_MVAEleIDCuts[0][1][1] =  0.917 ; // barrel (eta>0.8) pt above 10 GeV  
  m_MVAEleIDCuts[0][1][2] =  0.683 ; // endcap pt above 10 GeV            

  // 90%
  m_MVAEleIDCuts[1][0][0] = -0.483 ; // barrel (eta<0.8) pt 5-10 GeV     
  m_MVAEleIDCuts[1][0][1] = -0.267 ; // barrel (eta>0.8) pt 5-10 GeV     
  m_MVAEleIDCuts[1][0][2] = -0.323 ; // endcap pt 5-10 GeV               
  m_MVAEleIDCuts[1][1][0] = 0.933  ; // barrel (eta<0.8) pt above 10 GeV 
  m_MVAEleIDCuts[1][1][1] = 0.825  ; // barrel (eta>0.8) pt above 10 GeV 
  m_MVAEleIDCuts[1][1][2] = 0.337  ; // endcap pt above 10 GeV           
}



OfflineProducerHelper::OfflineProducerHelper(TH1F* hCounter){

  for(int itr=1;itr<=hCounter->GetNbinsX();itr++){
    TString binlabel = hCounter->GetXaxis()->GetBinLabel(itr);
    if(binlabel.BeginsWith("HLT"))triggerlist.push_back(hCounter->GetXaxis()->GetBinLabel(itr));
  }
  tauidlist.push_back("");

  // MVA ele ID from here:
  //   //  https://twiki.cern.ch/twiki/bin/view/CMS/MultivariateElectronIdentificationRun2#Non_triggering_electron_MVA
  //     // 80%

  m_MVAEleIDCuts[0][0][0] = -0.253 ; // barrel (eta<0.8) pt 5-10 GeV      
  m_MVAEleIDCuts[0][0][1] =  0.081 ; // barrel (eta>0.8) pt 5-10 GeV      
  m_MVAEleIDCuts[0][0][2] = -0.081 ; // endcap pt 5-10 GeV                
  m_MVAEleIDCuts[0][1][0] =  0.965 ; // barrel (eta<0.8) pt above 10 GeV  
  m_MVAEleIDCuts[0][1][1] =  0.917 ; // barrel (eta>0.8) pt above 10 GeV  
  m_MVAEleIDCuts[0][1][2] =  0.683 ; // endcap pt above 10 GeV            

  // 90%
  m_MVAEleIDCuts[1][0][0] = -0.483 ; // barrel (eta<0.8) pt 5-10 GeV     
  m_MVAEleIDCuts[1][0][1] = -0.267 ; // barrel (eta>0.8) pt 5-10 GeV     
  m_MVAEleIDCuts[1][0][2] = -0.323 ; // endcap pt 5-10 GeV               
  m_MVAEleIDCuts[1][1][0] = 0.933  ; // barrel (eta<0.8) pt above 10 GeV 
  m_MVAEleIDCuts[1][1][1] = 0.825  ; // barrel (eta>0.8) pt above 10 GeV 
  m_MVAEleIDCuts[1][1][2] = 0.337  ; // endcap pt above 10 GeV           
}


int OfflineProducerHelper::FindTriggerNumber(TString triggername){
  for(int it=0;it<triggerlist.size();it++){ 	
  	if(triggerlist.at(it).CompareTo(triggername.Data())==0)return it;
  	else {
  	    TString newName=triggername.Data();
  	    newName.Prepend("HLT_");
  	    newName.Append("_v");
  	    if(triggerlist.at(it).CompareTo(newName.Data())==0)return it;
  	}
  }
  return -1;
}

bool OfflineProducerHelper::IsTriggerFired(int triggerbit, int triggernumber){ 
  if(triggernumber>=0 && triggernumber<triggerlist.size()) return triggerbit & (1 << triggernumber);
  return false;
}

int OfflineProducerHelper::printFiredPaths(int triggerbit){
  int nFired = 0;
  for(int it=0;it<triggerlist.size();it++){ 	
  	if(IsTriggerFired(triggerbit,it)) {
  	  printf("%s\n",triggerlist.at(it).Data());
  	  nFired++;
  	  }
  }
  return nFired;
}

bool OfflineProducerHelper::checkBit (int word, int bitpos)
{
    bool res = word & (1 << bitpos);
    return res;
}

int OfflineProducerHelper::getTAUidNumber(TString tauIDname){
  int ntau = (int)tauidlist.size();
  for(int i=0;i<ntau;i++)
    if(tauidlist.at(i)==tauIDname.Data()) return i;
  return -1;
}


int OfflineProducerHelper::getPairType (int type1, int type2)
{
    int nmu = 0;
    int nele = 0;
    int ntau = 0;
    
    if (isMuon (type1) )     nmu++;
    if (isElectron (type1) ) nele++;
    if (isTau (type1) )      ntau++;

    if (isMuon (type2) )     nmu++;
    if (isElectron (type2) ) nele++;
    if (isTau (type2) )      ntau++;

    if (nmu == 1 && nele == 0 && ntau == 1) return (int) pairType::MuHad;
    if (nmu == 0 && nele == 1 && ntau == 1) return (int) pairType::EHad;
    if (nmu == 0 && nele == 0 && ntau == 2) return (int) pairType::HadHad;
    if (nmu == 2 && nele == 0 && ntau == 0) return (int) pairType::MuMu;
    if (nmu == 0 && nele == 2 && ntau == 0) return (int) pairType::EE;
    if (nmu == 1 && nele == 1 && ntau == 0) return (int) pairType::EMu;
    
    return -1;

}

bool OfflineProducerHelper::pairPassBaseline (bigTree* tree, int iPair, TString whatApply)
{
    int dau1index = tree->indexDau1->at(iPair);
    int dau2index = tree->indexDau2->at(iPair);
    int type1 = tree->particleType->at(dau1index);
    int type2 = tree->particleType->at(dau2index);
    int pairType = getPairType (type1, type2);

    bool isOS = tree->isOSCand->at(iPair);
    if (whatApply.Contains("OScharge") && !isOS) 
        return false; // do not even check the rest if requiring the charge
    if (whatApply.Contains("SScharge") && isOS) 
        return false; // for the same sign selection at the moment full selection over SS pairs

    // pairs are always ordered as: e mu | e tau | mu tau  (e < mu < tau)
    // if same type of particle, highest pt one is the first
    if (pairType == MuHad)
    {
        bool leg1 = muBaseline (tree, dau1index, 18., 2.1, 0.1, whatApply);
        bool leg2 = tauBaseline (tree, dau2index, 20., 2.3, 0, 1, 3.0, whatApply);
        return (leg1 && leg2);
    }

    if (pairType == EHad)
    {
        bool leg1 = eleBaseline (tree, dau1index, 23., 0.1, 0, whatApply);
        bool leg2 = tauBaseline (tree, dau2index, 20., 2.3, 3, 0, 3.0, whatApply);
        return (leg1 && leg2);
    }

    // ordered by pT and not by most isolated, but baseline asked in sync is the same...
    if (pairType == HadHad)
    {
        bool leg1 = tauBaseline (tree, dau1index, 45., 2.1, 0, 0, 2.0, whatApply);
        bool leg2 = tauBaseline (tree, dau2index, 45., 2.1, 0, 0, 2.0, whatApply);
        return (leg1 && leg2);
    }

    if (pairType == EMu)
    {
        bool leg1 = eleBaseline (tree, dau1index, 13., 0.15, 0, whatApply);
        bool leg2 = muBaseline (tree, dau2index, 9., 2.4, 0.15, whatApply);
        return (leg1 && leg2);
    }
    
    // e e, mu mu missing for the moment...
    if (pairType == EE) return false;
    if (pairType == MuMu) return false;
    
    return false;
        
}


bool 
OfflineProducerHelper::eleBaseline (bigTree* tree, int iDau, 
                                    float ptMin, float relIso, int MVAIDflag, 
                                    TString whatApply)
{ 
    float px = tree->daughters_px->at(iDau);
    float py = tree->daughters_py->at(iDau);
    float pz = tree->daughters_pz->at(iDau);
    float e =  tree->daughters_e->at(iDau);

    TLorentzVector p4 (px, py, pz, e);
 
    // bypasser(s) and taker according to the string content
    bool byp_vertexS = false;
    bool byp_idS  = false;
    bool byp_isoS = false;
    bool byp_ptS  = false;
    bool byp_etaS = false;

    // whatApply: use "All", "Iso", "LepID", pTMin", "etaMax", "againstEle", "againstMu", "Vertex", "SScharge"; separate various arguments with a semicolon
    if (!whatApply.Contains("All") && 
        !whatApply.Contains("SScharge") && 
        !whatApply.Contains("OScharge"))
    {
      byp_vertexS = byp_idS = byp_isoS = byp_ptS = byp_etaS = true;
      // set selections
      if (whatApply.Contains("Vertex")) byp_vertexS = false; 
      if (whatApply.Contains("Iso"))    byp_isoS = false; 
      if (whatApply.Contains("LepID"))  byp_idS = false; 
      if (whatApply.Contains("pTMin"))  byp_ptS = false; 
      if (whatApply.Contains("etaMax")) byp_etaS = false;
    }

    bool vertexS = (tree->dxy->at(iDau) < 0.045 && tree->dz->at(iDau) < 0.2) || byp_vertexS;
    bool ptS = (p4.Pt() > ptMin) || byp_ptS;
    bool etaS = (fabs(p4.Eta()) < 2.5) || byp_etaS;
    //bool idS = checkBit (tree->daughters_iseleCUT->at(iDau), 3) || byp_idS; // 3 is TIGHT ele id CUT BASED
    bool idS = EleMVAID (tree->discriminator->at (iDau), tree->daughters_SCeta->at (iDau), p4.Pt (), MVAIDflag) || byp_idS ; // 2015/07/09 PG
    //bool idS = tree->daughters_iseleBDT->at(iDau) || byp_idS; // use it in ntuples produced after 11 June 2015, contains tight WP bool  
    //bool idS = tightEleMVAID (tree->discriminator->at(iDau), TMath::Abs(p4.Eta())) || byp_idS; // APPROX! Using lepton eta and not super cluster eta, discriminator contains ele BDT  
    bool isoS = (tree->combreliso->at(iDau) < relIso) || byp_isoS;
    if (whatApply.Contains ("InvertIzo")) isoS = !isoS ;
    
    bool totalS = (vertexS && idS && isoS && ptS && etaS);
    return totalS;
    
}

bool OfflineProducerHelper::muBaseline (
     bigTree* tree, int iDau, float ptMin, 
     float etaMax, float relIso, TString whatApply)
{
    float px = tree->daughters_px->at(iDau);
    float py = tree->daughters_py->at(iDau);
    float pz = tree->daughters_pz->at(iDau);
    float e =  tree->daughters_e->at(iDau);

    TLorentzVector p4 (px, py, pz, e);
    int discr = tree->daughters_muonID->at(iDau);

    // bypasser(s) according to the string content
    bool byp_vertexS = false;
    bool byp_idS  = false;
    bool byp_isoS = false;
    bool byp_ptS  = false;
    bool byp_etaS = false;

    // whatApply: use "All", "Iso", "LepID", pTMin", "etaMax", "againstEle", 
    // "againstMu", "Vertex"; separate various arguments with a semicolon
    if (!whatApply.Contains("All") && 
        !whatApply.Contains("SScharge") && 
        !whatApply.Contains("OScharge"))
    {
      byp_vertexS = byp_idS = byp_isoS = byp_ptS = byp_etaS = true;
      // set selections
      if (whatApply.Contains("Vertex")) byp_vertexS = false; 
      if (whatApply.Contains("Iso"))    byp_isoS = false; 
      if (whatApply.Contains("LepID"))  byp_idS = false; 
      if (whatApply.Contains("pTMin"))  byp_ptS = false; 
      if (whatApply.Contains("etaMax")) byp_etaS = false;
    }
        
    bool vertexS = (tree->dxy->at(iDau) < 0.045 && tree->dz->at(iDau) < 0.2) || byp_vertexS;
    bool idS = checkBit (discr, 2) || byp_idS; // bit 2 is MEDIUM mu id
    bool isoS = (tree->combreliso->at(iDau) < relIso) || byp_isoS;
    if (whatApply.Contains ("InvertIzo")) isoS = !isoS ;
    bool ptS = (p4.Pt() > ptMin) || byp_ptS;
    bool etaS = (fabs(p4.Eta()) < etaMax) || byp_etaS;
    
    bool totalS = (vertexS && idS && isoS && ptS && etaS);
    return totalS;
}

// againstEleWP: 0 = VLoose, 1 = Loose, 2 = Medium, 3 = Tight, 4 = VTight [all are MVA discr]
// againstMuWP: 0 = Loose, 1 = Tight
bool OfflineProducerHelper::tauBaseline (bigTree* tree, int iDau, float ptMin, 
         float etaMax, int againstEleWP, int againstMuWP, float isoRaw3Hits, 
         TString whatApply)
{
    float px = tree->daughters_px->at(iDau);
    float py = tree->daughters_py->at(iDau);
    float pz = tree->daughters_pz->at(iDau);
    float e =  tree->daughters_e->at(iDau);

    TLorentzVector p4 (px, py, pz, e);

    // bypasser(s) according to the string content
    bool byp_vertexS = false;
    bool byp_dmfS  = false;
    bool byp_agEleS = false;
    bool byp_agMuS = false;
    bool byp_isoS = false;
    bool byp_ptS  = false;
    bool byp_etaS = false;

    // whatApply: use "All", "Iso", "LepID", pTMin", "etaMax", "againstEle", 
    // "againstMu", "Vertex"; separate various arguments with a semicolon
    if (!whatApply.Contains("All") && 
        !whatApply.Contains("SScharge") && 
        !whatApply.Contains("OScharge"))
    {
      byp_vertexS = byp_dmfS = byp_agEleS = byp_agMuS = byp_isoS = byp_ptS = byp_etaS = true;
      // set selections
      if (whatApply.Contains("Vertex")) byp_vertexS = false; 
      if (whatApply.Contains("Iso"))    byp_isoS = false; 
      if (whatApply.Contains("LepID"))  byp_dmfS = false; 
      if (whatApply.Contains("againstEle"))  byp_agEleS = false; 
      if (whatApply.Contains("againstMu"))   byp_agMuS = false;       
      if (whatApply.Contains("pTMin"))  byp_ptS = false; 
      if (whatApply.Contains("etaMax")) byp_etaS = false;
    }


    if (againstEleWP < 0 || againstEleWP > 4) {
        cout << " ** OfflineProducerHelper::tauBaseline: againstEleWP must be between 0 and 4 --> using 0" << endl;
        againstEleWP = 0;
    } 

    if (againstMuWP < 0 || againstMuWP > 1) {
        cout << " ** OfflineProducerHelper::tauBaseline: againstMuWP must be between 0 and 1 --> using 0" << endl;
        againstMuWP = 0;
    }
    
    int agEleVal = 0;
    int agMuVal = 0;
    
    // ag ele:
    if (againstEleWP == 0)      agEleVal = checkBit(tree->tauID->at(iDau),getTAUidNumber("againstElectronVLooseMVA5"));
    else if (againstEleWP == 1) agEleVal = checkBit(tree->tauID->at(iDau),getTAUidNumber("againstElectronLooseMVA5"));
    else if (againstEleWP == 2) agEleVal = checkBit(tree->tauID->at(iDau),getTAUidNumber("againstElectronMediumMVA5"));
    else if (againstEleWP == 3) agEleVal = checkBit(tree->tauID->at(iDau),getTAUidNumber("againstElectronTightMVA5"));
    else if (againstEleWP == 4) agEleVal = checkBit(tree->tauID->at(iDau),getTAUidNumber("againstElectronVTightMVA5"));

    // ag mu:
    if (againstMuWP == 0)      agMuVal = checkBit(tree->tauID->at(iDau),getTAUidNumber("againstMuonLoose3"));
    else if (againstMuWP == 1) agMuVal = checkBit(tree->tauID->at(iDau),getTAUidNumber("againstMuonTight3"));

    bool dmfS = (tree->daughters_decayModeFindingOldDMs->at(iDau) == 1 || tree->daughters_decayModeFindingNewDMs->at(iDau) == 1) || byp_dmfS;
    bool vertexS = (tree->dxy->at(iDau) < 0.045 && tree->dz->at(iDau) < 0.2) || byp_vertexS;
    bool agEleS = (agEleVal == 1) || byp_agEleS; 
    bool agMuS  = (agMuVal == 1) || byp_agMuS; 
    bool isoS = (tree->daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(iDau) < isoRaw3Hits) || byp_isoS;
    if (whatApply.Contains ("InvertIzo")) isoS = !isoS ;

    bool ptS = (p4.Pt() > ptMin) || byp_ptS;
    bool etaS = (fabs(p4.Eta()) < etaMax) || byp_etaS;

    bool totalS = (dmfS && vertexS && agEleS && agMuS && isoS && ptS && etaS);
    return totalS;    
}

// branch isBDT has been updated, this will soon be unnecessary (and eventually obsolete)
bool OfflineProducerHelper::tightEleMVAID (float BDT, float fSCeta)
{
    // POG_MVA_ID_Run2_NonTrig_Tight PHYS14
    if (fSCeta < 0) fSCeta *= -1;
    bool isBDT = false;
    if (fSCeta <0.8) isBDT = (BDT>0.73);
    else if (fSCeta < 1.479) isBDT = (BDT>0.57);
    else isBDT = (BDT >0.05);
    
    return isBDT;
}

int OfflineProducerHelper::getMothPairType (bigTree* tree, int iMoth)
{
    if (iMoth < 0 || iMoth >= tree->indexDau1->size())
    {
      cout << "warning ** getMothPairType: iMoth out of range" << endl;
      return -1;
    }

    int iDau1 = tree->indexDau1->at(iMoth);
    int iDau2 = tree->indexDau2->at(iMoth);
    int type1 = tree->particleType->at(iDau1);
    int type2 = tree->particleType->at(iDau2);

    return (getPairType(type1, type2));
}


bool 
OfflineProducerHelper::EleMVAID (float BDT, float eta, float pT, int strength)
  {
    // 0 = tight (80%), 1 = loose (90%)
    if (pT < 5) return false ;
    int etaBin = 0 ;
    if (fabs (eta > 0.8)) etaBin = 1 ;
    if (fabs (eta > 1.479)) etaBin = 2 ;
    return BDT > m_MVAEleIDCuts[strength][pT>10][etaBin] ;
  }


TLorentzVector OfflineProducerHelper::buildDauP4 (bigTree* tree, int iDau)
{
    float px = tree->daughters_px->at(iDau);
    float py = tree->daughters_py->at(iDau);
    float pz = tree->daughters_pz->at(iDau);
    float e =  tree->daughters_e->at(iDau);

    TLorentzVector p4 (px, py, pz, e);
    return p4;
}

TLorentzVector OfflineProducerHelper::buildMothP4 (bigTree* tree, int iMoth)
{
    float px = tree->mothers_px->at(iMoth);
    float py = tree->mothers_py->at(iMoth);
    float pz = tree->mothers_pz->at(iMoth);
    float e =  tree->mothers_e->at(iMoth);

    TLorentzVector p4 (px, py, pz, e);
    return p4;
}

TLorentzVector OfflineProducerHelper::buildGenP4 (bigTree* tree, int iGen)
{
    float px = tree->genpart_px->at(iGen);
    float py = tree->genpart_py->at(iGen);
    float pz = tree->genpart_pz->at(iGen);
    float e =  tree->genpart_e->at(iGen);

    TLorentzVector p4 (px, py, pz, e);
    return p4;
}

int OfflineProducerHelper::MCHiggsTauTauDecayMode (bigTree* tree)
{
    int decay = -1; // good decays go from 0 to 7, see enum
    for (int i = 0; i < tree->genpart_HZDecayMode->size(); i++)
    {
        int val = tree->genpart_HZDecayMode->at(i);
        if (val >= 0 && val <= 7)
        {
            decay = val;
            break; // I don;t expect more than 1 H to tau tau per event
        }
    }
    return decay; 
}


bool OfflineProducerHelper::getBestJets (bigTree* tree, int& jet1, int& jet2, int strategy)
{
    jet1 = jet2 = -1;
    switch (strategy)
    {
        case(0): // two with highest b score
        {
            int njets = tree->bCSVscore->size();
            if (njets < 2) return false;
            std::vector<std::pair<float, int>> scores;
            for (int i = 0; i < njets; i++) scores.push_back (std::make_pair(tree->bCSVscore->at(i), i));
            std::sort (scores.begin(), scores.end()); // are sorted according to the first index, i.e. the CSV score
            jet1 = (scores.at(njets-1)).second; //leading
            jet2 = (scores.at(njets-2)).second; // subleading
            return true;                            
            break;
        }
        
        default:
        {
            std::cout << "B jet selection strategy" << strategy << " not implemented" << std::endl;
            return false;
        }
    }
}


int OfflineProducerHelper::getBestPair (bigTree* tree, std::vector<int>& pairIdxs, TString strategy)
{
  int strat = -1;
  if (strategy == "OSMaxPt") strat = 0;
  if (strat == -1)
  {
    cout << "Best pair strategy not defined: " << strategy << endl;
    return -1;
  }

  // empty vector --> no pair
  if (pairIdxs.size() == 0) return -1;

  // ======= STRATEGIES ======
  // 0: favor OS, then max scalar sum PT
  if (strat == 0)
  {
    // pairs are already ordered by OS , SS, and inside each category by decreasing scalar sum pT
    // so it is enough to order indexes
    return ( *std::min_element (pairIdxs.begin(), pairIdxs.end()) );
  }

  return -1; // redudant, for safety
}


int OfflineProducerHelper::getBestPair (bigTree* tree, TString strategy)
{
  // FIXME: TO IMPLEMENT YET
  cout << "GET BEST PAIR FUNCTION TO BE IMPLEMENTED" << endl;
  return -1;
}

int OfflineProducerHelper::getPairByIndexes (bigTree* tree, int dau1, int dau2)
{
    int pair = -1;
    for (int iPair = 0; iPair < tree->indexDau1->size(); iPair++)
    {
        int ind1 = tree->indexDau1->at(iPair);
        int ind2 = tree->indexDau2->at(iPair);
        if (ind1 == dau1 && ind2 == dau2) pair = iPair;
        else if (ind2 == dau1 && ind1 == dau2) pair = iPair;
        
        if (pair != -1) break; // don't continue search, pairs are unique
    }
    return pair;
}



bool OfflineProducerHelper::getHardTauFinalVisGenProducts (bigTree* tree, int& ind1, int& ind2)
{
        int finalProds = 0;
        for (int iPart = 0; iPart < tree->genpart_pdg->size(); iPart++)
        {   
            int HInd = tree->genpart_HMothInd->at(iPart);
            int TauInd = tree->genpart_TauMothInd->at(iPart);
            int Pdg = tree->genpart_pdg->at(iPart);
            int aPdg = abs(Pdg);
            bool fromH = ( HInd != -1 );
            bool fromTau = ( TauInd != -1 );
            bool isPromptTauDecayProduct = checkBit (tree->genpart_flags->at(iPart), 5);
            //cout << iPart << " " << Pdg << " | fromH: " << fromH << " | fromTau: " << fromTau << endl;
            
            // e, mu
            if (aPdg == 11 || aPdg == 13)
            {
                if (fromH && fromTau && isPromptTauDecayProduct)
                {
                    if (finalProds == 0)      ind1 = iPart;
                    else if (finalProds == 1) ind2 = iPart;
                    finalProds++;
                }
            }
            
            // tau h
            if (aPdg == 66615)
            {
                bool promptTau = checkBit (tree->genpart_flags->at(TauInd), 0); // check if mother tau is prompt
                if (fromH && fromTau && promptTau)
                {
                    if (finalProds == 0)      ind1 = iPart;
                    else if (finalProds == 1) ind2 = iPart;
                    finalProds++;            
                }
            }
        }
                
        if (finalProds != 2) return false; // I expect 2 and only 2 visible products (e, mu, tauh 66615)
        else return true;
        
}


bool OfflineProducerHelper::drMatchGenReco (bigTree* tree, int iGen, int iReco, float dRcone)
{
    TLorentzVector genP4  (tree->genpart_px->at(iGen), tree->genpart_py->at(iGen), tree->genpart_pz->at(iGen), tree->genpart_e->at(iGen));
    TLorentzVector recoP4 (tree->daughters_px->at(iReco), tree->daughters_py->at(iReco), tree->daughters_pz->at(iReco), tree->daughters_e->at(iReco));
    
    if (genP4.DeltaR(recoP4) < dRcone) return true;
    else return false;
}

int OfflineProducerHelper::getRecoMatchedToGen (bigTree* tree, int iGen, bool checkId, bool checkCharge, float dRcone)
{
    TLorentzVector genP4  (tree->genpart_px->at(iGen), tree->genpart_py->at(iGen), tree->genpart_pz->at(iGen), tree->genpart_e->at(iGen));
    int genID = tree->genpart_pdg->at(iGen);
    // change to tau Id for tauh
    if (abs(genID) == 66615) genID = 15*genID/abs(genID);
    int AgenID = abs(genID);
    
    std::vector < std::pair<float, int> > matchedReco;
    for (int iReco = 0; iReco < tree->daughters_px->size(); iReco++)
    {
        int recoID = tree->PDGIdDaughters->at(iReco);
        bool IDCheck;
        if (checkCharge) IDCheck = ( genID == recoID );
        else             IDCheck = ( AgenID == abs(recoID) );
        
        if (!checkId) IDCheck = true; // bypass this requirement if I don't want ID to be checked        
        if (IDCheck)
        {
            TLorentzVector recoP4 (tree->daughters_px->at(iReco), tree->daughters_py->at(iReco), tree->daughters_pz->at(iReco), tree->daughters_e->at(iReco));
            float dR = genP4.DeltaR(recoP4);
            if (dR < dRcone)   matchedReco.push_back (std::make_pair(dR, iReco));
        }
    }
    
    // find closest matched
    if (matchedReco.size() == 0) return -1;
    std::sort (matchedReco.begin(), matchedReco.end()); // are sorted according to the first index, i.e. the dR
    return ((matchedReco.at(0)).second );
}


#endif

