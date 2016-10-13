#ifndef THTH_TREE_H
#define	THTH_TREE_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TSystem.h"
#include "HTauTauTree.h"
#include "myHelper.h"
#include "OfflineProducerHelper.h"

using namespace std;

unsigned int run, lumi, evt, NUP;
bool secondMuon   = false;
HTauTauTree* tree;
float lheHt; int lheNOutPartons;
float npu, rho, npv=-1, puweight, weight;
float trigweight_1, e_1, px_1, py_1, pz_1, pt_1, phi_1, eta_1, m_1, q_1, d0_1, dZ_1, mt_1, iso_1, l1_decayMode;
bool id_m_highpt_1, id_m_loose_1, id_m_medium_1, id_m_tight_1, id_m_tightnovtx_1, id_m_hihghpt_1, id_e_mva_nt_loose_1, id_e_mva_nt_tight_1, id_e_cut_veto_1, id_e_cut_loose_1, id_e_cut_medium_1, id_e_cut_tight_1;
float e_2, px_2, py_2, pz_2, pt_2, phi_2, eta_2, m_2, q_2, d0_2, dZ_2, mt_2, iso_2, l2_decayMode;
bool id_m_highpt_2, id_m_loose_2, id_m_medium_2, id_m_tight_2, id_m_tightnovtx_2, id_m_hihghpt_2, id_e_mva_nt_loose_2, id_e_mva_nt_tight_2, id_e_cut_veto_2, id_e_cut_loose_2, id_e_cut_medium_2, id_e_cut_tight_2, againstElectronLooseMVA6_2, againstElectronMediumMVA6_2, againstElectronTightMVA6_2, againstElectronVLooseMVA6_2, againstElectronVTightMVA6_2, againstMuonLoose3_2, againstMuonTight3_2, byLooseCombinedIsolationDeltaBetaCorr3Hits_2, byMediumCombinedIsolationDeltaBetaCorr3Hits_2, byTightCombinedIsolationDeltaBetaCorr3Hits_2, byVLooseIsolationMVArun2v1DBoldDMwLT_2, byLooseIsolationMVArun2v1DBoldDMwLT_2, byMediumIsolationMVArun2v1DBoldDMwLT_2, byTightIsolationMVArun2v1DBoldDMwLT_2, byVTightIsolationMVArun2v1DBoldDMwLT_2;
int decayModeFindingOldDMs_2, decayModeFindingNewDMs_2;
float trigweight_2, byCombinedIsolationDeltaBetaCorrRaw3Hits_2, byIsolationMVA3oldDMwoLTraw_2, byIsolationMVA3oldDMwLTraw_2, byIsolationMVA3newDMwoLTraw_2, byIsolationMVA3newDMwLTraw_2, chargedIsoPtSum_2, neutralIsoPtSum_2, puCorrPtSum_2;

float sigmaIetaIeta_1, deltaPhiSuperCluster_1, deltaEtaSuperCluster_1, depositR03tracker_1, depositR03ecal_1, depositR03hcal_1, trackIso_1, ecalIso_1, hcalIso_1; 
float sigmaIetaIeta_2, deltaPhiSuperCluster_2, deltaEtaSuperCluster_2, depositR03tracker_2, depositR03ecal_2, depositR03hcal_2, trackIso_2, ecalIso_2, hcalIso_2;
float pt_tt, m_vis, ptvis;
float m_sv_tr, m_sv_tr_tauUP,m_sv_tr_tauDOWN, m_sv, m_sv_tauUP,m_sv_tauDOWN,pt_sv,pt_sv_tauUP,pt_sv_tauDOWN,ptu_sv,ptu_sv_tauUP,ptu_sv_tauDOWN,eta_sv,eta_sv_tauUP,eta_sv_tauDOWN,etau_sv,etau_sv_tauUP,etau_sv_tauDOWN,phi_sv,phi_sv_tauUP, phi_sv_tauDOWN,phiu_sv,phiu_sv_tauUP,phiu_sv_tauDOWN,fitMETRho_sv,fitMETRho_sv_tauUP,fitMETRho_sv_tauDOWN,fitMETPhi_sv,fitMETPhi_sv_tauUP,fitMETPhi_sv_tauDOWN;
float met_pt, met_phi, met_px, met_py, umet_px, umet_py;
float mvamet, mvametphi, pzetavis, pzetamiss, mvacov00, mvacov01, mvacov11, mvacov10;
int nb_extra_electrons, nb_extra_muons;
float top_reweighting, gen_Z_pt, gen_Z_mass, gen_Z_eta;
bool extraelec_veto, extramuon_veto, dilepton_veto;
int ptype_1, ptype_2;
int   gen_match_1, gen_match_2;
float gen_match_px1, gen_match_py1, gen_match_pz1, gen_match_e1;   
float gen_match_px2, gen_match_py2, gen_match_pz2, gen_match_e2;   
float gen_match_tau_px1, gen_match_tau_py1, gen_match_tau_pz1, gen_match_tau_e1;   
float gen_match_tau_px2, gen_match_tau_py2, gen_match_tau_pz2, gen_match_tau_e2;   
vector<float>  jets_px;
vector<float>  jets_py;
vector<float>  jets_pz;
vector<float>  jets_e;
vector<float>  jets_rawPt;
vector<float>  jets_area;
vector<float>  jets_mT;
vector<int>    jets_Flavour;
vector<int>    jets_HadronFlavour;
vector<int>    jets_genjetIndex;
vector<float>  jets_PUJetID;
vector<float>  jets_PUJetIDupdated;
vector<float>  jets_jecUnc;
vector<float>  bCSVscore;
vector<int>    PFjetID;
vector<float>  jetRawf;
vector<float>  genjet_px;
vector<float>  genjet_py;
vector<float>  genjet_pz;
vector<float>  genjet_e;
vector<int>    genjet_partonFlavour;
vector<int>    genjet_hadronFlavour;

void fillTree(TTree *Run_Tree, int entry_tree, int indice_paire, int indice_tau1, int indice_tau2, TLorentzVector MEt, bool isData, TH1F* hTauID, bool ismc, float pu) {

	run=-999; lumi=-999; evt=-999; NUP=-999;
	lheHt=-999;  lheNOutPartons=-999;
	npu-999, rho=-999, npv-999, puweight=-999, weight=-999;
	e_1=-999, px_1=-999, py_1=-999, pz_1=-999, pt_1=-999, phi_1=-999, eta_1=-999, m_1=-999, q_1=-999, d0_1=-999, dZ_1=-999, mt_1=-999, iso_1=-999, l1_decayMode=-999;
	e_2=-999, px_2=-999, py_2=-999, pz_2=-999, pt_2=-999, phi_2=-999, eta_2=-999, m_2=-999, q_2=-999, d0_2=-999, dZ_2=-999, mt_2=-999, iso_2=-999, l2_decayMode=-999;
	secondMuon  = false;
	trigweight_1=-999;
	byLooseCombinedIsolationDeltaBetaCorr3Hits_2=false, byMediumCombinedIsolationDeltaBetaCorr3Hits_2=false, byTightCombinedIsolationDeltaBetaCorr3Hits_2=false;
	trigweight_2-999, byCombinedIsolationDeltaBetaCorrRaw3Hits_2=-999, byIsolationMVA3oldDMwoLTraw_2=-999, byIsolationMVA3oldDMwLTraw_2=-999, byIsolationMVA3newDMwoLTraw_2=-999, byIsolationMVA3newDMwLTraw_2=-999, chargedIsoPtSum_2=-999, neutralIsoPtSum_2=-999, puCorrPtSum_2=-999;
	byVLooseIsolationMVArun2v1DBoldDMwLT_2=false, byLooseIsolationMVArun2v1DBoldDMwLT_2=false, byMediumIsolationMVArun2v1DBoldDMwLT_2=false, byTightIsolationMVArun2v1DBoldDMwLT_2=false, byVTightIsolationMVArun2v1DBoldDMwLT_2=false;
        sigmaIetaIeta_1=-999, deltaPhiSuperCluster_1=-999, deltaEtaSuperCluster_1=-999, depositR03tracker_1=-999, depositR03ecal_1=-999, depositR03hcal_1=-999, trackIso_1=-999, ecalIso_1=-999, hcalIso_1=-999; 
	sigmaIetaIeta_2=-999, deltaPhiSuperCluster_2=-999, deltaEtaSuperCluster_2=-999, depositR03tracker_2=-999, depositR03ecal_2=-999, depositR03hcal_2=-999, trackIso_2=-999, ecalIso_2=-999, hcalIso_2=-999;
	pt_tt=-999, m_vis=-999, ptvis=-999;
	m_sv_tr=-999, m_sv_tr_tauUP=-999,m_sv_tr_tauDOWN=-999, m_sv=-999, m_sv_tauUP=-999,m_sv_tauDOWN=-999,pt_sv=-999,pt_sv_tauUP=-999,pt_sv_tauDOWN=-999,ptu_sv=-999,ptu_sv_tauUP=-999,ptu_sv_tauDOWN=-999,eta_sv=-999,eta_sv_tauUP=-999,eta_sv_tauDOWN=-999,etau_sv=-999,etau_sv_tauUP=-999,etau_sv_tauDOWN=-999,phi_sv=-999,phi_sv_tauUP=-999, phi_sv_tauDOWN=-999,phiu_sv=-999,phiu_sv_tauUP=-999,phiu_sv_tauDOWN=-999,fitMETRho_sv=-999,fitMETRho_sv_tauUP=-999,fitMETRho_sv_tauDOWN=-999,fitMETPhi_sv=-999,fitMETPhi_sv_tauUP=-999,fitMETPhi_sv_tauDOWN=-999;
	met_pt=-999, met_phi=-999, met_px=-999, met_py=-999, umet_px=-999, umet_py=-999;
	mvamet=-999, mvametphi=-999, pzetavis=-999, pzetamiss=-999, mvacov00=-999, mvacov01=-999, mvacov11=-999, mvacov10=-999;
	nb_extra_electrons=-999, nb_extra_muons=-999;
	top_reweighting=-999, gen_Z_pt=-999, gen_Z_mass=-999, gen_Z_eta=-999;
	extraelec_veto=false, extramuon_veto=false, dilepton_veto=false;
	ptype_1=-999, ptype_2=-999;

	id_m_highpt_1=false, id_m_loose_1=false, id_m_medium_1=false, id_m_tight_1=false, id_m_tightnovtx_1=false, id_m_hihghpt_1=false, id_e_mva_nt_loose_1=false, id_e_mva_nt_tight_1=false, id_e_cut_veto_1=false, id_e_cut_loose_1=false, id_e_cut_medium_1=false, id_e_cut_tight_1=false;
	id_m_highpt_2=false, id_m_loose_2=false, id_m_medium_2=false, id_m_tight_2=false, id_m_tightnovtx_2=false, id_m_hihghpt_2=false, id_e_mva_nt_loose_2=false, id_e_mva_nt_tight_2=false, id_e_cut_veto_2=false, id_e_cut_loose_2=false, id_e_cut_medium_2=false, id_e_cut_tight_2=false, againstElectronLooseMVA6_2=false, againstElectronMediumMVA6_2=false, againstElectronTightMVA6_2=false, againstElectronVLooseMVA6_2=false, againstElectronVTightMVA6_2=false, againstMuonLoose3_2=false, againstMuonTight3_2=false, decayModeFindingOldDMs_2=-999, decayModeFindingNewDMs_2=-999;

	gen_match_1    = -999;
	gen_match_2    = -999;
	gen_match_px1  = -999;
	gen_match_py1  = -999;
	gen_match_pz1  = -999;
	gen_match_e1   = -999;
	gen_match_px2  = -999;
	gen_match_py2  = -999;
	gen_match_pz2  = -999;
	gen_match_e2   = -999;
	gen_match_tau_px1  = -999;
	gen_match_tau_py1  = -999;
	gen_match_tau_pz1  = -999;
	gen_match_tau_e1   = -999;
	gen_match_tau_px2  = -999;
	gen_match_tau_py2  = -999;
	gen_match_tau_pz2  = -999;
	gen_match_tau_e2   = -999;

	jets_px.clear();
	jets_py.clear();
	jets_pz.clear();
	jets_e.clear();
	jets_rawPt.clear();
	jets_area.clear();
	jets_mT.clear();
	jets_Flavour.clear();
	jets_HadronFlavour.clear();
	jets_genjetIndex.clear();
	jets_PUJetID.clear();
	jets_PUJetIDupdated.clear();
	jets_jecUnc.clear();
	bCSVscore.clear();
	PFjetID.clear();
	jetRawf.clear();
	genjet_px.clear();
	genjet_py.clear();
	genjet_pz.clear();
	genjet_e.clear();
	genjet_partonFlavour.clear();
	genjet_hadronFlavour.clear();
	
        run  = tree->RunNumber;
	lumi = tree->lumi;
	evt  = tree->EventNumber;

	for (int iJet = 0; iJet < tree->bDiscriminator->size(); iJet++) {
		jets_px.push_back(tree->jets_px->at(iJet));
		jets_py.push_back(tree->jets_py->at(iJet));
		jets_pz.push_back(tree->jets_pz->at(iJet));
		jets_e.push_back(tree->jets_e->at(iJet));
		jets_rawPt.push_back(tree->jets_rawPt->at(iJet));
		jets_area.push_back(tree->jets_area->at(iJet));
		jets_mT.push_back(tree->jets_mT->at(iJet));
		jets_Flavour.push_back(tree->jets_Flavour->at(iJet));
		jets_HadronFlavour.push_back(tree->jets_HadronFlavour->at(iJet));
		if(ismc) jets_genjetIndex.push_back(tree->jets_genjetIndex->at(iJet));
		jets_PUJetID.push_back(tree->jets_PUJetID->at(iJet));
		jets_PUJetIDupdated.push_back(tree->jets_PUJetIDupdated->at(iJet));
		jets_jecUnc.push_back(tree->jets_jecUnc->at(iJet));
		bCSVscore.push_back(tree->bCSVscore->at(iJet));
		PFjetID.push_back(tree->PFjetID->at(iJet));
		jetRawf.push_back(tree->jetRawf->at(iJet));
	}

	if(ismc){
		if(jets_genjetIndex.size()<tree->genjet_px->size()){
			int diff = tree->genjet_px->size()-jets_genjetIndex.size(); 
			for (int iRJet = 0; iRJet < diff; iRJet++) {
				jets_genjetIndex.push_back(-1);
			}
		}

		for (int iGenJet = 0; iGenJet < tree->genjet_px->size(); iGenJet++) {
			if(jets_genjetIndex.at(iGenJet)>-1){
				genjet_px.push_back(tree->genjet_px->at(jets_genjetIndex.at(iGenJet)));
				genjet_py.push_back(tree->genjet_py->at(jets_genjetIndex.at(iGenJet)));
				genjet_pz.push_back(tree->genjet_pz->at(jets_genjetIndex.at(iGenJet)));
				genjet_e.push_back(tree->genjet_e->at(jets_genjetIndex.at(iGenJet)));
				genjet_partonFlavour.push_back(tree->genjet_partonFlavour->at(jets_genjetIndex.at(iGenJet)));
				genjet_hadronFlavour.push_back(tree->genjet_hadronFlavour->at(jets_genjetIndex.at(iGenJet)));
			}else{
				genjet_px.push_back(-999);
				genjet_py.push_back(-999);
				genjet_pz.push_back(-999);
				genjet_e.push_back(-999);
				genjet_partonFlavour.push_back(-999);
				genjet_hadronFlavour.push_back(-999);

			}
		}
	}

	nb_extra_electrons=0; nb_extra_muons=0;

	mvacov00 = tree->MET_cov00->at(indice_paire);
	mvacov01 = tree->MET_cov01->at(indice_paire);
	mvacov10 = tree->MET_cov10->at(indice_paire);
	mvacov11 = tree->MET_cov11->at(indice_paire);

	m_sv_tr=tree->SVfitTransverseMass->at(indice_paire);
	m_sv_tr_tauUP=tree->SVfitTransverseMassTauUp->at(indice_paire);
	m_sv_tr_tauDOWN=tree->SVfitTransverseMassTauDown->at(indice_paire);
	m_sv=tree->SVfitMass->at(indice_paire);
	m_sv_tauUP=tree->SVfitMassTauUp->at(indice_paire);
	m_sv_tauDOWN=tree->SVfitMassTauDown->at(indice_paire);
	pt_sv=tree->SVfit_pt->at(indice_paire);
	pt_sv_tauUP=tree->SVfit_ptTauUp->at(indice_paire);
	pt_sv_tauDOWN=tree->SVfit_ptTauDown->at(indice_paire);
	ptu_sv=tree->SVfit_ptUnc->at(indice_paire);
	ptu_sv_tauUP=tree->SVfit_ptUncTauUp->at(indice_paire);
	ptu_sv_tauDOWN=tree->SVfit_ptUncTauDown->at(indice_paire);
	eta_sv=tree->SVfit_eta->at(indice_paire);
	eta_sv_tauUP=tree->SVfit_etaTauUp->at(indice_paire);
	eta_sv_tauDOWN=tree->SVfit_etaTauDown->at(indice_paire);
	etau_sv=tree->SVfit_etaUnc->at(indice_paire);
	etau_sv_tauUP=tree->SVfit_etaUncTauUp->at(indice_paire);
	etau_sv_tauDOWN=tree->SVfit_etaUncTauDown->at(indice_paire);
	phi_sv=tree->SVfit_phi->at(indice_paire);
	phi_sv_tauUP=tree->SVfit_phiTauUp->at(indice_paire);
	phi_sv_tauDOWN=tree->SVfit_phiTauDown->at(indice_paire);
	phiu_sv=tree->SVfit_phiUnc->at(indice_paire);
	phiu_sv_tauUP=tree->SVfit_phiUncTauUp->at(indice_paire);
	phiu_sv_tauDOWN=tree->SVfit_phiUncTauDown->at(indice_paire);
	fitMETRho_sv=tree->SVfit_fitMETRho->at(indice_paire);
	fitMETRho_sv_tauUP=tree->SVfit_fitMETRhoTauUp->at(indice_paire);
	fitMETRho_sv_tauDOWN=tree->SVfit_fitMETRhoTauDown->at(indice_paire);
	fitMETPhi_sv=tree->SVfit_fitMETPhi->at(indice_paire);
	fitMETPhi_sv_tauUP=tree->SVfit_fitMETPhiTauUp->at(indice_paire);
	fitMETPhi_sv_tauDOWN=tree->SVfit_fitMETPhiTauDown->at(indice_paire);

	TLorentzVector ETM (tree->METx->at(indice_paire),tree->METy->at(indice_paire),0,0);
	met_pt=ETM.Pt();
	met_phi=ETM.Phi();
	met_px=ETM.Px();
	met_py=ETM.Py();
	umet_px=tree->uncorrMETx->at(indice_paire);
	umet_py=tree->uncorrMETy->at(indice_paire);


	if(ismc) weight=tree->aMCatNLOweight;
	if(ismc) lheHt=tree->lheHt;
	if(ismc) lheNOutPartons=tree->lheNOutPartons;
	puweight=pu;

	TLorentzVector tau1;
	TLorentzVector tau2;
        ptype_1 = tree->particleType->at(indice_tau1);
        ptype_2 = tree->particleType->at(indice_tau2);

        tau2.SetPxPyPzE(tree->daughters_px->at(indice_tau2),tree->daughters_py->at(indice_tau2),tree->daughters_pz->at(indice_tau2),tree->daughters_e->at(indice_tau2));
	tau1.SetPxPyPzE(tree->daughters_px->at(indice_tau1),tree->daughters_py->at(indice_tau1),tree->daughters_pz->at(indice_tau1),tree->daughters_e->at(indice_tau1));

	TLorentzVector lep;
	for (int iDau = 0; iDau < tree->daughters_px->size(); iDau++){

		
                lep.SetPxPyPzE(tree->daughters_px->at(iDau),tree->daughters_py->at(iDau),tree->daughters_pz->at(iDau),tree->daughters_e->at(iDau));
		
                if (!Overlap_2(tau1,lep) && tree->particleType->at(iDau)==0 && tree->combreliso->at(iDau)<0.3 && fabs(tree->dxy->at(iDau))<0.045 && fabs(tree->dz->at(iDau))<0.2 && fabs(lep.Eta())<2.4 && lep.Pt()>15)
			secondMuon=true;
		if (tree->particleType->at(iDau)==1 && tree->combreliso->at(iDau)<0.3 && tree->daughters_iseleWP90->at(iDau) && fabs(tree->dxy->at(iDau))<0.045 && fabs(tree->dz->at(iDau))<0.2 && fabs(lep.Eta())<2.5 && lep.Pt()>10)
			extraelec_veto=true;
		if (!Overlap_015(tau1,lep) && tree->particleType->at(iDau)==0 && tree->combreliso->at(iDau)<0.3 && fabs(tree->dxy->at(iDau))<0.045 && fabs(tree->dz->at(iDau))<0.2 && fabs(lep.Eta())<2.4 && lep.Pt()>10 && ((tree->daughters_muonID->at(iDau) >> 2) & 1))
			extramuon_veto=true;
		//GOOD FOR MMMMMMM !!! if (!Overlap_015(tau1,lep) && !Overlap_015(tau2,lep) && tree->particleType->at(iDau)==0 && tree->combreliso->at(iDau)<0.3 && fabs(tree->dxy->at(iDau))<0.045 && fabs(tree->dz->at(iDau))<0.2 && fabs(lep.Eta())<2.4 && lep.Pt()>10 && ((tree->daughters_muonID->at(iDau) >> 2) & 1))
		//	extramuon_veto=true;
		if (!Overlap_015(tau1,lep) && tree->daughters_charge->at(indice_tau1)*tree->daughters_charge->at(iDau)<0 && tree->particleType->at(iDau)==0 && tree->combreliso->at(iDau)<0.3 && fabs(tree->dxy->at(iDau))<0.045 && fabs(tree->dz->at(iDau))<0.2 && fabs(lep.Eta())<2.4 && lep.Pt()>15 && ((tree->daughters_typeOfMuon->at(iDau) >> 0) & 1) && ((tree->daughters_typeOfMuon->at(iDau) >> 1) & 1) && ((tree->daughters_typeOfMuon->at(iDau) >> 2) & 1))
			dilepton_veto=true;

	}

	int indice_gen_1=0; int indice_gen_2=0; float newdr1=0.5; float newdr2=0.5;
	if (ismc){

                bool isFoundZ = false;
		for (int iGen = 0; iGen < tree->genpart_px->size() && !isFoundZ; iGen++){
			TLorentzVector genZ;
			//cout << " DEBUGGGG " << tree->genpart_pdg->at(iGen) << endl;
			if (tree->genpart_pdg->at(iGen)==23){
				//cout << "DEBUGGGG I have found the Z "<<endl;
				genZ.SetPxPyPzE(tree->genpart_px->at(iGen),tree->genpart_py->at(iGen),tree->genpart_pz->at(iGen),tree->genpart_e->at(iGen));
				gen_Z_pt = genZ.Pt(); gen_Z_mass = genZ.M(); gen_Z_eta = genZ.Eta();		
                                isFoundZ = true;
                                break;
			}
		}
 
		if(!isFoundZ){
                        bool isFound=false;
			TLorentzVector genZ,g1,g2;
			for (int iGen1 = 0; iGen1 < tree->genpart_px->size() && !isFound; iGen1++){
				for (int iGen2 = iGen1+1; iGen2 < tree->genpart_px->size() && !isFound; iGen2++){
					if(fabs(tree->genpart_pdg->at(iGen1))==fabs(tree->genpart_pdg->at(iGen2)) && 
					  (fabs(tree->genpart_pdg->at(iGen1))==11 || fabs(tree->genpart_pdg->at(iGen1))==13 || fabs(tree->genpart_pdg->at(iGen1))==15)){
						//cout << " tree->genpart_pdg->at(iGen1) " << tree->genpart_pdg->at(iGen1) << " tree->genpart_pdg->at(iGen2) " << tree->genpart_pdg->at(iGen2) << endl; 	
						g1.SetPxPyPzE(tree->genpart_px->at(iGen1),tree->genpart_py->at(iGen1),tree->genpart_pz->at(iGen1),tree->genpart_e->at(iGen1));
						g2.SetPxPyPzE(tree->genpart_px->at(iGen2),tree->genpart_py->at(iGen2),tree->genpart_pz->at(iGen2),tree->genpart_e->at(iGen2));
						isFound=true;
						break;
					}
				}
			}
			genZ = g1+g2;
			gen_Z_pt = genZ.Pt(); gen_Z_mass = genZ.M(); gen_Z_eta = genZ.Eta();		
		}


		for (int iGen = 2; iGen < tree->genpart_px->size(); iGen++){
			TLorentzVector g;
			g.SetPxPyPzE(tree->genpart_px->at(iGen),tree->genpart_py->at(iGen),tree->genpart_pz->at(iGen),tree->genpart_e->at(iGen));
			if ((fabs(tree->genpart_pdg->at(iGen))==11 or fabs(tree->genpart_pdg->at(iGen))==13 or fabs(tree->genpart_pdg->at(iGen))==66615)){
				if (g.DeltaR(tau2)<0.5 && g.DeltaR(tau2)<newdr2){
					newdr2=g.DeltaR(tau2); indice_gen_2=iGen;
					gen_match_px2=g.Px() ; gen_match_py2=g.Py(); gen_match_pz2=g.Pz() ; gen_match_e2=g.E(); 
				}
				if (g.DeltaR(tau1)<0.5 && g.DeltaR(tau1)<newdr1){
					newdr1=g.DeltaR(tau1); indice_gen_1=iGen;
					gen_match_px1=g.Px() ; gen_match_py1=g.Py(); gen_match_pz1=g.Pz() ; gen_match_e1=g.E(); 
				}
			}
			else if (fabs(tree->genpart_pdg->at(iGen))==15){
				if (g.DeltaR(tau2)<0.5 && g.DeltaR(tau2)<newdr2){
					newdr2=g.DeltaR(tau2); indice_gen_2=iGen; 
					gen_match_tau_px2=g.Px() ; gen_match_tau_py2=g.Py(); gen_match_tau_pz2=g.Pz() ; gen_match_tau_e2=g.E(); 
				}
				if (g.DeltaR(tau1)<0.5 && g.DeltaR(tau1)<newdr1){
					newdr1=g.DeltaR(tau1); indice_gen_1=iGen;
					gen_match_tau_px1=g.Px() ; gen_match_tau_py1=g.Py(); gen_match_tau_pz1=g.Pz() ; gen_match_tau_e1=g.E(); 
				}
			}
		}
		//cout << " genZ_pt/genZ_mass/genZ_eta "<< gen_Z_pt << "/" << gen_Z_mass << "/" << gen_Z_eta<< endl;

		int indice_gen1=tree->daughters_genindex->at(indice_tau1);
		indice_gen1=indice_gen_1;
		if (indice_gen1>-1){
			if (fabs(tree->genpart_pdg->at(indice_gen1))==11 && sqrt(pow(tree->genpart_px->at(indice_gen1),2)+pow(tree->genpart_py->at(indice_gen1),2))>8 && ((tree->genpart_flags->at(indice_gen1) >> 0) & 1)) gen_match_1=1; 
			else if (fabs(tree->genpart_pdg->at(indice_gen1))==13 && sqrt(pow(tree->genpart_px->at(indice_gen1),2)+pow(tree->genpart_py->at(indice_gen1),2))>8 && ((tree->genpart_flags->at(indice_gen1) >> 0) & 1)) gen_match_1=2;
			else if (fabs(tree->genpart_pdg->at(indice_gen1))==11 && sqrt(pow(tree->genpart_px->at(indice_gen1),2)+pow(tree->genpart_py->at(indice_gen1),2))>8 && ((tree->genpart_flags->at(indice_gen1) >> 0) & 5)) gen_match_1=3;
			else if (fabs(tree->genpart_pdg->at(indice_gen1))==13 && sqrt(pow(tree->genpart_px->at(indice_gen1),2)+pow(tree->genpart_py->at(indice_gen1),2))>8 && ((tree->genpart_flags->at(indice_gen1) >> 0) & 5)) gen_match_1=4;
			else if ((fabs(tree->genpart_pdg->at(indice_gen1))==66615 or fabs(tree->genpart_pdg->at(indice_gen1))==15) && sqrt(pow(tree->genpart_px->at(indice_gen1),2)+pow(tree->genpart_py->at(indice_gen1),2))>15) gen_match_1=5;
			else gen_match_1=6;
		}
		else gen_match_1=6;

		int indice_gen2=tree->daughters_genindex->at(indice_tau2);
		indice_gen2=indice_gen_2;
		if (indice_gen2>-1){
			if (fabs(tree->genpart_pdg->at(indice_gen2))==11 && sqrt(pow(tree->genpart_px->at(indice_gen2),2)+pow(tree->genpart_py->at(indice_gen2),2))>8 && ((tree->genpart_flags->at(indice_gen2) >> 0) & 1)) gen_match_2=1;
			else if (fabs(tree->genpart_pdg->at(indice_gen2))==13 && sqrt(pow(tree->genpart_px->at(indice_gen2),2)+pow(tree->genpart_py->at(indice_gen2),2))>8 && ((tree->genpart_flags->at(indice_gen2) >> 0) & 1)) gen_match_2=2;
			else if (fabs(tree->genpart_pdg->at(indice_gen2))==11 && sqrt(pow(tree->genpart_px->at(indice_gen2),2)+pow(tree->genpart_py->at(indice_gen2),2))>8 && ((tree->genpart_flags->at(indice_gen2) >> 0) & 5)) gen_match_2=3;
			else if (fabs(tree->genpart_pdg->at(indice_gen2))==13 && sqrt(pow(tree->genpart_px->at(indice_gen2),2)+pow(tree->genpart_py->at(indice_gen2),2))>8 && ((tree->genpart_flags->at(indice_gen2) >> 0) & 5)) gen_match_2=4;
			else if ((fabs(tree->genpart_pdg->at(indice_gen2))==66615 or fabs(tree->genpart_pdg->at(indice_gen2))==15) && sqrt(pow(tree->genpart_px->at(indice_gen2),2)+pow(tree->genpart_py->at(indice_gen2),2))>15) gen_match_2=5;
			else gen_match_2=6;
		}
		else gen_match_2=6;
	}

	m_vis=(tau1+tau2).M();

	l1_decayMode=tree->decayMode->at(indice_tau1);
	l2_decayMode=tree->decayMode->at(indice_tau2);

	m_1                 = tau1.M();
	px_1                = tau1.Px();
	py_1                = tau1.Py();
	pz_1                = tau1.Pz();
	e_1                 = tau1.E();
	pt_1                = tau1.Pt();
	phi_1               = tau1.Phi();
	eta_1               = tau1.Eta();
	mt_1                = tree->mT_Dau1->at(indice_paire);
	dZ_1                = tree->dz->at(indice_tau1);
	d0_1                = tree->dxy->at(indice_tau1);
	depositR03tracker_1 = tree->daughters_depositR03_tracker->at(indice_tau1);
	depositR03ecal_1    = tree->daughters_depositR03_ecal->at(indice_tau1);
	depositR03hcal_1    = tree->daughters_depositR03_hcal->at(indice_tau1);
	id_m_loose_1        = tree->discriminator->at(indice_tau1)==1;
	id_m_medium_1       = tree->discriminator->at(indice_tau1)==2;
	id_m_tight_1        = tree->discriminator->at(indice_tau1)==3;
	iso_1               = tree->combreliso->at(indice_tau1);
	q_1                 = tree->daughters_charge->at(indice_tau1);

	m_2                 = tau2.M();
	px_2                = tau2.Px();
	py_2                = tau2.Py();
	pz_2                = tau2.Pz();
	e_2                 = tau2.E();
	pt_2                = tau2.Pt();
	phi_2               = tau2.Phi();
	eta_2               = tau2.Eta();
	mt_2                = tree->mT_Dau2->at(indice_paire);
	dZ_2                = tree->dz->at(indice_tau2);
	d0_2                = tree->dxy->at(indice_tau2);
	depositR03tracker_2 = tree->daughters_depositR03_tracker->at(indice_tau2);
	depositR03ecal_2    = tree->daughters_depositR03_ecal->at(indice_tau2);
	depositR03hcal_2    = tree->daughters_depositR03_hcal->at(indice_tau2);
	id_m_loose_2        = tree->discriminator->at(indice_tau2)==1;
	id_m_medium_2       = tree->discriminator->at(indice_tau2)==2;
	id_m_tight_2        = tree->discriminator->at(indice_tau2)==3;
	iso_2               = tree->combreliso->at(indice_tau2);
	q_2                 = tree->daughters_charge->at(indice_tau2);

	OfflineProducerHelper* HelperTauID = new OfflineProducerHelper(hTauID,hTauID);
	againstMuonTight3_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("againstMuonTight3")) & 1);
	againstMuonLoose3_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("againstMuonLoose3")) & 1);
	againstElectronVLooseMVA6_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("againstElectronVLooseMVA6")) & 1);
	againstElectronLooseMVA6_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("againstElectronLooseMVA6")) & 1);
	againstElectronMediumMVA6_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("againstElectronMediumMVA6")) & 1);
	againstElectronTightMVA6_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("againstElectronTightMVA6")) & 1);
	againstElectronVTightMVA6_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("againstElectronVTightMVA6")) & 1);
	byLooseCombinedIsolationDeltaBetaCorr3Hits_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("byLooseCombinedIsolationDeltaBetaCorr3Hits")) & 1);
	byMediumCombinedIsolationDeltaBetaCorr3Hits_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("byMediumCombinedIsolationDeltaBetaCorr3Hits")) & 1);
 	byTightCombinedIsolationDeltaBetaCorr3Hits_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("byTightCombinedIsolationDeltaBetaCorr3Hits")) & 1);
	byVLooseIsolationMVArun2v1DBoldDMwLT_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("byVLooseIsolationMVArun2v1DBoldDMwLT")) & 1);
	byLooseIsolationMVArun2v1DBoldDMwLT_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("byLooseIsolationMVArun2v1DBoldDMwLT")) & 1);
	byMediumIsolationMVArun2v1DBoldDMwLT_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("byMediumIsolationMVArun2v1DBoldDMwLT")) & 1);
	byTightIsolationMVArun2v1DBoldDMwLT_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("byTightIsolationMVArun2v1DBoldDMwLT")) & 1);
	byVTightIsolationMVArun2v1DBoldDMwLT_2 = ((tree->tauID->at(indice_tau2) >>  HelperTauID->getTAUidNumber("byVTightIsolationMVArun2v1DBoldDMwLT")) & 1);

        byCombinedIsolationDeltaBetaCorrRaw3Hits_2=tree->daughters_byCombinedIsolationDeltaBetaCorrRaw3Hits->at(indice_tau2);
	byIsolationMVA3oldDMwoLTraw_2=tree->daughters_byIsolationMVA3oldDMwoLTraw->at(indice_tau2);
	byIsolationMVA3newDMwoLTraw_2=tree->daughters_byIsolationMVA3newDMwoLTraw->at(indice_tau2);
	byIsolationMVA3oldDMwLTraw_2=tree->daughters_byIsolationMVA3oldDMwLTraw->at(indice_tau2);
	byIsolationMVA3newDMwLTraw_2=tree->daughters_byIsolationMVA3newDMwLTraw->at(indice_tau2);
	neutralIsoPtSum_2=tree->daughters_neutralIsoPtSum->at(indice_tau2);
	chargedIsoPtSum_2=tree->daughters_chargedIsoPtSum->at(indice_tau2);
	puCorrPtSum_2=tree->daughters_puCorrPtSum->at(indice_tau2);
	decayModeFindingOldDMs_2=tree->daughters_decayModeFindingOldDMs->at(indice_tau2);
	decayModeFindingNewDMs_2=tree->daughters_decayModeFindingNewDMs->at(indice_tau2);

	TLorentzVector h=tau1+tau2;

	ptvis=h.Pt();

	float norm_zeta=norm_F(tau1.Px()/tau1.Pt()+tau2.Px()/tau2.Pt(),tau1.Py()/tau1.Pt()+tau2.Py()/tau2.Pt());
	float x_zeta= (tau1.Px()/tau1.Pt()+tau2.Px()/tau2.Pt())/norm_zeta;
	float y_zeta= (tau1.Py()/tau1.Pt()+tau2.Py()/tau2.Pt())/norm_zeta;
	pzetamiss=ETM.Px()*x_zeta+ETM.Py()*y_zeta;
	pzetavis=(tau1.Px()+tau2.Px())*x_zeta+(tau1.Py()+tau2.Py())*y_zeta;

	NUP=tree->NUP;
	pt_tt=(tau1+tau2+ETM).Pt();
	npu=tree->npu;
	npv=tree->npv;
	rho=tree->rho;

	trigweight_1=1.0;
	trigweight_2=1.0;

	Run_Tree->Fill();

}

#endif


