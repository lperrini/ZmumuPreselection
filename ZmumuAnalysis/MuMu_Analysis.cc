#include <iostream>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <utility>
#include <map>
#include <string>
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"
#include "TMath.h"
#include "TSystem.h"
#include "TRandom.h"
#include "TLorentzVector.h"
#include "TRandom3.h"
#include "makeHisto.h"
#include "mutau_Tree.h"
#include "OfflineProducerHelper.h"
#include "LumiReweightingStandAlone.h"

int main(int argc, char** argv) {

	using namespace std;
	myMap1 = new map<string, TH1F*>();
	myMap2 = new map<string, TH2F*>();

	string status_sample  = *(argv + 1);

	bool isMC        = false;
	bool isData_RunC = false;
	bool isData_RunD = false;
	bool isData      = false;

	if (status_sample.compare("mc") == 0) isMC = true;
	if (status_sample.compare("dataC") == 0) isData_RunC = true;
	if (status_sample.compare("dataD") == 0) isData_RunD = true;
        if (isData_RunC || isData_RunD) isData = true;

	string out = *(argv + 2);
	string outname= out;

	string in = *(argv + 3);
	string inname = in;

	TChain *fIn = new TChain("HTauTauTree/HTauTauTree");
	TH1F *evCounter = NULL, *htauid = NULL;

	/* read file list */
	std::ifstream ifs( (inname).c_str() );
	//std::ifstream ifs( (string(std::getenv("CMSSW_BASE"))+"/src/ZmumuAnalysis/"+inname).c_str() );
	std::string tmpName = "";
	int i=0;
	while( ifs.good() ){
		tmpName = "";
		getline(ifs,tmpName);
		if( !tmpName.empty() ){ 
			fIn->Add( tmpName.c_str() );
			TFile *t = new TFile( tmpName.c_str() );
			if(i==0){
				TH1F *m_evCounter   = (TH1F*) t->Get("HTauTauTree/Counters");
				TH1F *m_htauid      = (TH1F*) t->Get("HTauTauTree/TauIDs");
				evCounter = (TH1F*)m_evCounter->Clone();
				htauid = (TH1F*)m_htauid->Clone();
			}
			else{
				TH1F *m_evCounter   = (TH1F*) t->Get("HTauTauTree/Counters");
				TH1F *m_htauid      = (TH1F*) t->Get("HTauTauTree/TauIDs");
				evCounter->Add(m_evCounter);
				htauid->Add(m_htauid);
			}
		}
		i++;
	}

	TTree* treePtr    = (TTree*) fIn;//->Get("HTauTauTree/HTauTauTree");
	//HTauTauTree* tree = new HTauTauTree (treePtr);
	tree = new HTauTauTree(treePtr,isMC);
	OfflineProducerHelper helper;

	//#################################################################################################
	TFile *fout = TFile::Open(outname.c_str(), "RECREATE");
	TTree *Run_Tree = new TTree("RLE_tree", "RLE_tree");
	Run_Tree->Branch("run", &run, "run/I");
	Run_Tree->Branch("lumi", &lumi, "lumi/I");
	Run_Tree->Branch("evt", &evt, "evt/I");

	Run_Tree->Branch("weight", &weight, "weight/F");
	Run_Tree->Branch("lheNOutPartons", &lheNOutPartons, "lheNOutPartons/I");
	Run_Tree->Branch("NUP", &NUP, "NUP/I");
	Run_Tree->Branch("lheHt", &lheHt, "lheHt/F");
	Run_Tree->Branch("genjet_px", &genjet_px);//, "genjet_px/F");
	Run_Tree->Branch("genjet_py", &genjet_py);//, "genjet_py/F");
	Run_Tree->Branch("genjet_pz", &genjet_pz);//, "genjet_pz/F");
	Run_Tree->Branch("genjet_e", &genjet_e);//, "genjet_e/F"); 
	Run_Tree->Branch("genjet_partonFlavour", &genjet_partonFlavour);//, "genjet_partonFlavour/I");
	Run_Tree->Branch("genjet_hadronFlavour", &genjet_hadronFlavour);//, "genjet_hadronFlavour/I");
	Run_Tree->Branch("gen_match_1",  &gen_match_1,    "gen_match_1/I");
	Run_Tree->Branch("gen_match_2",   &gen_match_2,   "gen_match_2/I");
	Run_Tree->Branch("gen_match_px1", &gen_match_px1, "gen_match_px1/F");
	Run_Tree->Branch("gen_match_py1", &gen_match_py1, "gen_match_py1/F");
	Run_Tree->Branch("gen_match_pz1", &gen_match_pz1, "gen_match_pz1/F");
	Run_Tree->Branch("gen_match_e1",  &gen_match_e1,  "gen_match_e1/F");
	Run_Tree->Branch("gen_match_px2", &gen_match_px2, "gen_match_px2/F");
	Run_Tree->Branch("gen_match_py2", &gen_match_py2, "gen_match_py2/F");
	Run_Tree->Branch("gen_match_pz2", &gen_match_pz2, "gen_match_pz2/F");
	Run_Tree->Branch("gen_match_e2",  &gen_match_e2,  "gen_match_e2/F");
	Run_Tree->Branch("gen_match_tau_px1", &gen_match_tau_px1, "gen_match_tau_px1/F");
	Run_Tree->Branch("gen_match_tau_py1", &gen_match_tau_py1, "gen_match_tau_py1/F");
	Run_Tree->Branch("gen_match_tau_pz1", &gen_match_tau_pz1, "gen_match_tau_pz1/F");
	Run_Tree->Branch("gen_match_tau_e1",  &gen_match_tau_e1,  "gen_match_tau_e1/F");
	Run_Tree->Branch("gen_match_tau_px2", &gen_match_tau_px2, "gen_match_tau_px2/F");
	Run_Tree->Branch("gen_match_tau_py2", &gen_match_tau_py2, "gen_match_tau_py2/F");
	Run_Tree->Branch("gen_match_tau_pz2", &gen_match_tau_pz2, "gen_match_tau_pz2/F");
	Run_Tree->Branch("gen_match_tau_e2",  &gen_match_tau_e2,  "gen_match_tau_e2/F");
	Run_Tree->Branch("gen_Z_pt", &gen_Z_pt, "gen_Z_pt/F");
	Run_Tree->Branch("gen_Z_mass", &gen_Z_mass, "gen_Z_mass/F");
	Run_Tree->Branch("gen_Z_eta", &gen_Z_eta, "gen_Z_eta/F");
	
        Run_Tree->Branch("secondMuon", &secondMuon, "secondMuon/O");
	Run_Tree->Branch("dilepton_veto", &dilepton_veto, "dilepton_veto/O");
	Run_Tree->Branch("extraelec_veto", &extraelec_veto, "extraelec_veto/O");
	Run_Tree->Branch("extramuon_veto", &extramuon_veto, "extramuon_veto/O");
	Run_Tree->Branch("puweight", &puweight, "puweight/F");
	Run_Tree->Branch("npv", &npv, "npv/F");
	Run_Tree->Branch("npu", &npu, "npu/F");
	Run_Tree->Branch("rho", &rho, "rho/F");

	Run_Tree->Branch("ptype_1", &ptype_1, "ptype_1/I");
	Run_Tree->Branch("pt_1", &pt_1, "pt_1/F");
	Run_Tree->Branch("px_1", &px_1, "px_1/F");
	Run_Tree->Branch("py_1", &py_1, "py_1/F");
	Run_Tree->Branch("pz_1", &pz_1, "pz_1/F");
	Run_Tree->Branch("e_1", &e_1, "e_1/F");
	Run_Tree->Branch("phi_1", &phi_1, "phi_1/F");
	Run_Tree->Branch("eta_1", &eta_1, "eta_1/F");
	Run_Tree->Branch("m_1", &m_1, "m_1/F");
	Run_Tree->Branch("q_1", &q_1, "q_1/F");
	Run_Tree->Branch("d0_1", &d0_1, "d0_1/F");
	Run_Tree->Branch("dZ_1", &dZ_1, "dZ_1/F");
	Run_Tree->Branch("mt_1", &mt_1, "mt_1/F");
	Run_Tree->Branch("iso_1", &iso_1, "iso_1/F");
	Run_Tree->Branch("l1_decayMode", &l1_decayMode, "l1_decayMode/F");
	Run_Tree->Branch("id_m_loose_1", &id_m_loose_1, "id_m_loose_1/O");
	Run_Tree->Branch("id_m_medium_1", &id_m_medium_1, "id_m_medium_1/O");
	Run_Tree->Branch("id_m_tight_1", &id_m_tight_1, "id_m_tight_1/O");
	Run_Tree->Branch("id_m_tightnovtx_1", &id_m_tightnovtx_1, "id_m_tightnovtx_1/O");
	Run_Tree->Branch("id_m_highpt_1", &id_m_highpt_1, "id_m_highpt_1/O");
	Run_Tree->Branch("id_e_mva_nt_loose_1", &id_e_mva_nt_loose_1, "id_e_mva_nt_loose_1/O");
	Run_Tree->Branch("id_e_mva_nt_tight_1", &id_e_mva_nt_tight_1, "id_e_mva_nt_tight_1/O");
	Run_Tree->Branch("id_e_cut_veto_1", &id_e_cut_veto_1, "id_e_cut_veto_1/O");
	Run_Tree->Branch("id_e_cut_loose_1", &id_e_cut_loose_1, "id_e_cut_loose_1/O");
	Run_Tree->Branch("id_e_cut_medium_1", &id_e_cut_medium_1, "id_e_cut_medium_1/O");
	Run_Tree->Branch("id_e_cut_tight_1", &id_e_cut_tight_1, "id_e_cut_tight_1/O");
	Run_Tree->Branch("trigweight_1", &trigweight_1, "trigweight_1/F");

	Run_Tree->Branch("ptype_2", &ptype_2, "ptype_2/I");
	Run_Tree->Branch("pt_2", &pt_2, "pt_2/F");
	Run_Tree->Branch("phi_2", &phi_2, "phi_2/F");
	Run_Tree->Branch("eta_2", &eta_2, "eta_2/F");
	Run_Tree->Branch("px_2", &px_2, "px_2/F");
	Run_Tree->Branch("py_2", &py_2, "py_2/F");
	Run_Tree->Branch("pz_2", &pz_2, "pz_2/F");
	Run_Tree->Branch("e_2", &e_2, "e_2/F");
	Run_Tree->Branch("m_2", &m_2, "m_2/F");
	Run_Tree->Branch("q_2", &q_2, "q_2/F");
	Run_Tree->Branch("d0_2", &d0_2, "d0_2/F");
	Run_Tree->Branch("dZ_2", &dZ_2, "dZ_2/F");
	Run_Tree->Branch("mt_2", &mt_2, "mt_2/F");
	Run_Tree->Branch("iso_2", &iso_2, "iso_2/F");
	Run_Tree->Branch("l2_decayMode", &l2_decayMode, "l2_decayMode/F");
	Run_Tree->Branch("id_m_loose_2", &id_m_loose_2, "id_m_loose_2/O");
	Run_Tree->Branch("id_m_medium_2", &id_m_medium_2, "id_m_medium_2/O");
	Run_Tree->Branch("id_m_tight_2", &id_m_tight_2, "id_m_tight_2/O");
	Run_Tree->Branch("id_m_tightnovtx_2", &id_m_tightnovtx_2, "id_m_tightnovtx_2/O");
	Run_Tree->Branch("id_m_hihghpt_2", &id_m_highpt_2, "id_m_highpt_2/O");
	Run_Tree->Branch("id_e_mva_nt_loose_2", &id_e_mva_nt_loose_2, "id_e_mva_nt_loose_2/O");
	Run_Tree->Branch("id_e_mva_nt_tight_2", &id_e_mva_nt_tight_2, "id_e_mva_nt_tight_2/O");
	Run_Tree->Branch("id_e_cut_veto_2", &id_e_cut_veto_2, "id_e_cut_veto_2/O");
	Run_Tree->Branch("id_e_cut_loose_2", &id_e_cut_loose_2, "id_e_cut_loose_2/O");
	Run_Tree->Branch("id_e_cut_medium_2", &id_e_cut_medium_2, "id_e_cut_medium_2/O");
	Run_Tree->Branch("id_e_cut_tight_2", &id_e_cut_tight_2, "id_e_cut_tight_2/O");
	Run_Tree->Branch("trigweight_2", &trigweight_2, "trigweight_2/F");
	Run_Tree->Branch("againstElectronLooseMVA6_2", &againstElectronLooseMVA6_2, "againstElectronLooseMVA6_2/O");
	Run_Tree->Branch("againstElectronMediumMVA6_2", &againstElectronMediumMVA6_2, "againstElectronMediumMVA6_2/O");
	Run_Tree->Branch("againstElectronTightMVA6_2", &againstElectronTightMVA6_2, "againstElectronTightMVA6_2/O");
	Run_Tree->Branch("againstElectronVLooseMVA6_2", &againstElectronVLooseMVA6_2, "againstElectronVLooseMVA6_2/O");
	Run_Tree->Branch("againstElectronVTightMVA6_2", &againstElectronVTightMVA6_2, "againstElectronVTightMVA6_2/O");
	Run_Tree->Branch("againstMuonLoose3_2", &againstMuonLoose3_2, "againstMuonLoose3_2/O");
	Run_Tree->Branch("againstMuonTight3_2", &againstMuonTight3_2, "againstMuonTight3_2/O");
        Run_Tree->Branch("byVLooseIsolationMVArun2v1DBoldDMwLT_2", &byVLooseIsolationMVArun2v1DBoldDMwLT_2, "byVLooseIsolationMVArun2v1DBoldDMwLT_2/O");
        Run_Tree->Branch("byLooseIsolationMVArun2v1DBoldDMwLT_2", &byLooseIsolationMVArun2v1DBoldDMwLT_2, "byLooseIsolationMVArun2v1DBoldDMwLT_2/O");
        Run_Tree->Branch("byMediumIsolationMVArun2v1DBoldDMwLT_2", &byMediumIsolationMVArun2v1DBoldDMwLT_2, "byMediumIsolationMVArun2v1DBoldDMwLT_2/O");
        Run_Tree->Branch("byTightIsolationMVArun2v1DBoldDMwLT_2", &byTightIsolationMVArun2v1DBoldDMwLT_2, "byTightIsolationMVArun2v1DBoldDMwLT_2/O");
        Run_Tree->Branch("byVTightIsolationMVArun2v1DBoldDMwLT_2", &byVTightIsolationMVArun2v1DBoldDMwLT_2, "byVTightIsolationMVArun2v1DBoldDMwLT_2/O");
        Run_Tree->Branch("byLooseCombinedIsolationDeltaBetaCorr3Hits_2", &byLooseCombinedIsolationDeltaBetaCorr3Hits_2, "byLooseCombinedIsolationDeltaBetaCorr3Hits_2/O");
	Run_Tree->Branch("byMediumCombinedIsolationDeltaBetaCorr3Hits_2", &byMediumCombinedIsolationDeltaBetaCorr3Hits_2, "byMediumCombinedIsolationDeltaBetaCorr3Hits_2/O");
	Run_Tree->Branch("byTightCombinedIsolationDeltaBetaCorr3Hits_2", &byTightCombinedIsolationDeltaBetaCorr3Hits_2, "byTightCombinedIsolationDeltaBetaCorr3Hits_2/O");
	Run_Tree->Branch("byCombinedIsolationDeltaBetaCorrRaw3Hits_2", &byCombinedIsolationDeltaBetaCorrRaw3Hits_2, "byCombinedIsolationDeltaBetaCorrRaw3Hits_2/F");
	Run_Tree->Branch("byIsolationMVA3oldDMwoLTraw_2", &byIsolationMVA3oldDMwoLTraw_2, "byIsolationMVA3oldDMwoLTraw_2/F");
	Run_Tree->Branch("byIsolationMVA3oldDMwLTraw_2", &byIsolationMVA3oldDMwLTraw_2, "byIsolationMVA3oldDMwLTraw_2/F");
	Run_Tree->Branch("byIsolationMVA3newDMwoLTraw_2", &byIsolationMVA3newDMwoLTraw_2, "byIsolationMVA3newDMwoLTraw_2/F");
	Run_Tree->Branch("byIsolationMVA3newDMwLTraw_2", &byIsolationMVA3newDMwLTraw_2, "byIsolationMVA3newDMwLTraw_2/F");
	Run_Tree->Branch("chargedIsoPtSum_2", &chargedIsoPtSum_2, "chargedIsoPtSum_2/F");
	Run_Tree->Branch("decayModeFindingOldDMs_2", &decayModeFindingOldDMs_2, "decayModeFindingOldDMs_2/I");
	Run_Tree->Branch("decayModeFindingNewDMs_2", &decayModeFindingNewDMs_2, "decayModeFindingNewDMs_2/I");
	Run_Tree->Branch("neutralIsoPtSum_2", &neutralIsoPtSum_2, "neutralIsoPtSum_2/F");
	Run_Tree->Branch("puCorrPtSum_2", &puCorrPtSum_2, "puCorrPtSum_2/F");

	Run_Tree->Branch("pt_tt", &pt_tt, "pt_tt/F");
	Run_Tree->Branch("m_vis", &m_vis, "m_vis/F");
	Run_Tree->Branch("m_sv", &m_sv, "m_sv/F");
	Run_Tree->Branch("m_sv_tauUP", &m_sv_tauUP, "m_sv_tauUP/F");
	Run_Tree->Branch("m_sv_tauDOWN", &m_sv_tauDOWN, "m_sv_tauDOWN/F");
	Run_Tree->Branch("pt_sv", &pt_sv, "pt_sv/F");
	Run_Tree->Branch("ptu_sv", &ptu_sv, "ptu_sv/F");
	Run_Tree->Branch("pt_sv_tauUP", &pt_sv_tauUP, "pt_sv_tauUP/F");
	Run_Tree->Branch("pt_sv_tauDOWN", &pt_sv_tauDOWN, "pt_sv_tauDOWN/F");
	Run_Tree->Branch("ptu_sv_tauUP", &ptu_sv_tauUP, "ptu_sv_tauUP/F");
	Run_Tree->Branch("ptu_sv_tauDOWN", &ptu_sv_tauDOWN, "ptu_sv_tauDOWN/F");
	Run_Tree->Branch("eta_sv", &eta_sv, "eta_sv/F");
	Run_Tree->Branch("etau_sv", &etau_sv, "etau_sv/F");
	Run_Tree->Branch("eta_sv_tauUP", &eta_sv_tauUP, "eta_sv_tauUP/F");
	Run_Tree->Branch("eta_sv_tauDOWN", &eta_sv_tauDOWN, "eta_sv_tauDOWN/F");
	Run_Tree->Branch("etau_sv_tauUP", &etau_sv_tauUP, "etau_sv_tauUP/F");
	Run_Tree->Branch("etau_sv_tauDOWN", &etau_sv_tauDOWN, "etau_sv_tauDOWN/F");
	Run_Tree->Branch("phi_sv", &phi_sv, "phi_sv/F");
	Run_Tree->Branch("phi_sv_tauUP", &phi_sv_tauUP, "phi_sv_tauUP/F");
	Run_Tree->Branch("phi_sv_tauDOWN", &phi_sv_tauDOWN, "phi_sv_tauDOWN/F");
	Run_Tree->Branch("phiu_sv", &phiu_sv, "phiu_sv/F");
	Run_Tree->Branch("phiu_sv_tauUP", &phiu_sv_tauUP, "phiu_sv_tauUP/F");
	Run_Tree->Branch("phiu_sv_tauDOWN", &phiu_sv_tauDOWN, "phiu_sv_tauDOWN/F");
	Run_Tree->Branch("fitMETRho_sv", &fitMETRho_sv, "fitMETRho_sv/F");
	Run_Tree->Branch("fitMETRho_sv_tauUP", &fitMETRho_sv_tauUP, "fitMETRho_sv_tauUP/F");
	Run_Tree->Branch("fitMETRho_sv_tauDOWN", &fitMETRho_sv_tauDOWN, "fitMETRho_sv_tauDOWN/F");
	Run_Tree->Branch("fitMETPhi_sv", &fitMETPhi_sv, "fitMETPhi_sv/F");
	Run_Tree->Branch("fitMETPhi_sv_tauUP", &fitMETPhi_sv_tauUP, "fitMETPhi_sv_tauUP/F");
	Run_Tree->Branch("fitMETPhi_sv_tauDOWN", &fitMETPhi_sv_tauDOWN, "fitMETPhi_sv_tauDOWN/F");
	Run_Tree->Branch("m_sv_tr", &m_sv_tr, "m_sv_tr/F");
	Run_Tree->Branch("m_sv_tr_tauUP", &m_sv_tr_tauUP, "m_sv_tr_tauUP/F");
	Run_Tree->Branch("m_sv_tr_tauDOWN", &m_sv_tr_tauDOWN, "m_sv_tr_tauDOWN/F");

	Run_Tree->Branch("met_pt", &met_pt, "met_pt/F");
	Run_Tree->Branch("met_px", &met_px, "met_px/F");
	Run_Tree->Branch("met_py", &met_py, "met_py/F");
	Run_Tree->Branch("umet_px", &umet_px, "umet_px/F");
	Run_Tree->Branch("umet_py", &umet_py, "umet_py/F");
	Run_Tree->Branch("met_phi", &met_phi, "met_phi/F");
	Run_Tree->Branch("mvamet", &mvamet, "mvamet/F");
	Run_Tree->Branch("mvametphi", &mvametphi, "mvametphi/F");
	Run_Tree->Branch("pzetavis", &pzetavis, "pzetavis/F");
	Run_Tree->Branch("pzetamiss", &pzetamiss, "pzetamiss/F");
	Run_Tree->Branch("mvacov00", &mvacov00, "mvacov00/F");
	Run_Tree->Branch("mvacov10", &mvacov10, "mvacov10/F");
	Run_Tree->Branch("mvacov01", &mvacov01, "mvacov01/F");
	Run_Tree->Branch("mvacov11", &mvacov11, "mvacov11/F");
	Run_Tree->Branch("ptvis", &ptvis, "ptvis/F");
	Run_Tree->Branch("nb_extra_electrons", &nb_extra_electrons, "nb_extra_electrons/I");
	Run_Tree->Branch("nb_extra_muons", &nb_extra_muons, "nb_extra_muons/I");
	
        Run_Tree->Branch("jets_px", &jets_px);//, "jets_px/F");
	Run_Tree->Branch("jets_py", &jets_py);//, "jets_py/F");
	Run_Tree->Branch("jets_pz", &jets_pz);//, "jets_pz/F");
	Run_Tree->Branch("jets_e", &jets_e);//, "jets_e/F");
	Run_Tree->Branch("jets_rawPt", &jets_rawPt);//, "jets_rawPt/F");
	Run_Tree->Branch("jets_area", &jets_area);//, "jets_area/F");
	Run_Tree->Branch("jets_mT", &jets_mT);//, "jets_mT/F");
	Run_Tree->Branch("jets_Flavour", &jets_Flavour);//, "jets_Flavour/I");
	Run_Tree->Branch("jets_HadronFlavour", &jets_HadronFlavour);//, "jets_HadronFlavour/I");
	Run_Tree->Branch("jets_genjetIndex", &jets_genjetIndex);//, "jets_genjetIndex/I");
	Run_Tree->Branch("jets_PUJetID", &jets_PUJetID);//, "jets_PUJetID/F");
	Run_Tree->Branch("jets_PUJetIDupdated", &jets_PUJetIDupdated);//, "jets_PUJetIDupdated/F");
	Run_Tree->Branch("jets_jecUnc", &jets_jecUnc);//, "jets_jecUnc/F");
	Run_Tree->Branch("bCSVscore", &bCSVscore);//, "bCSVscore/F");
	Run_Tree->Branch("PFjetID", &PFjetID);//, "PFjetID/I");
	Run_Tree->Branch("jetRawf", &jetRawf);//, "jetRawf/F");

	//reweight::LumiReWeighting* LumiWeights_12;
	//LumiWeights_12 = new reweight::LumiReWeighting("pileup-hists/MC_Spring15_PU25_Startup.root", "pileup-hists/Data_Pileup_2015D_Nov17.root", "pileup", "pileup");
	//LumiWeights_12 = new reweight::LumiReWeighting(string(std::getenv("CMSSW_BASE"))+"/src/ZmumuAnalysis/pileup-hists/MC_Spring15_PU25_Startup.root", string(std::getenv("CMSSW_BASE"))+"/src/ZmumuAnalysis/pileup-hists/Data_Pileup_2015D_Nov17.root", "pileup", "pileup");

	double weight_plus=0;
	double weight_minus=0;
	//int ntrigger=0;
	//#################################################################################################
	for (int iEntry = 0; iEntry < tree->GetEntries() ; iEntry++)
	{
                //if(iEntry!=19283) continue;
                //cout << " Entry " << iEntry << endl;

		if ( tree->GetEntry(iEntry) < 0 ){
			fprintf(stderr,"Problem at %d\n",i); 
                        continue; 
		}		
                bool print=false;
		float pu=1.0;
		float aMCatNLO=1.0;
		if (isMC) {
			aMCatNLO=tree->aMCatNLOweight;
			//pu = LumiWeights_12->weight(tree->npu);
		}  

		if (iEntry % 1000 == 0) fprintf(stdout, "\r  Processed events: %8d ", iEntry);
                fflush(stdout);
		plotFill("nbevt",0,9,0,9);
		plotFill("nbevt",1,9,0,9,aMCatNLO);
		if(aMCatNLO>0) weight_plus++;
		else weight_minus++;

                if(print) cout << "isMC/isData_RunC/isData_RunD " <<isMC<<"/"<<isData_RunC<<"/"<<isData_RunD << endl;

		bool PassMuTrigger=false;
		if(print) cout<<"Event number "<<tree->EventNumber<<endl;
		if (isMC && helper.IsTriggerFired(tree->triggerbit,helper.FindTriggerNumber("HLT_IsoMu18_v"))) PassMuTrigger=true;
		if (isData_RunC && (helper.IsTriggerFired(tree->triggerbit,helper.FindTriggerNumber("HLT_IsoMu17_eta2p1_LooseIsoPFTau20_v")) || helper.IsTriggerFired(tree->triggerbit,helper.FindTriggerNumber("HLT_IsoMu24_eta2p1_v"))) ) PassMuTrigger=true;
		if (isData_RunD && helper.IsTriggerFired(tree->triggerbit,helper.FindTriggerNumber("HLT_IsoMu18_v"))) PassMuTrigger=true;
		if (!PassMuTrigger) continue; 

                if(print) cout << " trigger passed "<< endl;

		bool isPair=false;
                
		int tmp1index=0; int dau1index_=0; int dau1index=0; 
		int tmp2index=0; int dau2index_=0; int dau2index=0; 
		int pairindex=-1;
		float metx=0; float mety=0;

		//cout << " # particles " << tree->mothers_px->size() << endl;
		for (int iMoth = 0; iMoth < tree->mothers_px->size() && !isPair; iMoth++)
		{

			bool isMM = false;
			bool isMT = false;
			
			tmp1index = tree->indexDau1->at(iMoth);
			tmp2index = tree->indexDau2->at(iMoth);
 
                        //cout << "index 1/2 " << tmp1index << "/" << tmp2index << endl; 
                        if(print) cout << "particle type 1/2 " << tree->particleType->at(tmp1index) << "/" << tree->particleType->at(tmp2index) << endl;

                        if(tree->particleType->at(tmp1index) == 0 && tree->particleType->at(tmp2index) == 2) isMT = true;  
                        if(tree->particleType->at(tmp1index) == 0 && tree->particleType->at(tmp2index) == 0) isMM = true;

                        if(!isMT && !isMM) continue;

                        if(print) cout << " passed the channel selection " << endl;
                 
                        TLorentzVector tmp1,dau1;
			TLorentzVector tmp2,dau2;
			tmp1.SetPxPyPzE(tree->daughters_px->at(tmp1index),tree->daughters_py->at(tmp1index),tree->daughters_pz->at(tmp1index),tree->daughters_e->at(tmp1index));
			tmp2.SetPxPyPzE(tree->daughters_px->at(tmp2index),tree->daughters_py->at(tmp2index),tree->daughters_pz->at(tmp2index),tree->daughters_e->at(tmp2index));
		
			if(isMM){	
				if(tmp1.Pt()>tmp2.Pt()){
					dau1=tmp1;
					dau2=tmp2;
					dau1index_=tmp1index;
					dau2index_=tmp2index;
				}else{
					dau1=tmp2;
					dau2=tmp1;
					dau1index_=tmp2index;
					dau2index_=tmp1index;
				}
			}else if(isMT){
				dau1=tmp1;
				dau2=tmp2;
				dau1index_=tmp1index;
				dau2index_=tmp2index;
			}

			
                        if (Overlap_2(dau1,dau2)) continue;
                        if (print) cout << " overlap passed"<< endl;
			if (isMC && !((tree->daughters_FilterFired->at(dau1index_)  >> helper.FindTriggerNumber("HLT_IsoMu18_v")) & 1)) continue;
			if (isData_RunC && !((tree->daughters_FilterFired->at(dau1index_)  >> helper.FindTriggerNumber("HLT_IsoMu24_eta2p1_v")) & 1)) continue;
			if (isData_RunD && !((tree->daughters_FilterFired->at(dau1index_)  >> helper.FindTriggerNumber("HLT_IsoMu18_v")) & 1)) continue;
			if (isMC && tree->daughters_HLTpt->at(dau1index_)<18) continue;

                        /////// if is MUMU selection
			if(isMM){
				if (fabs(tree->dxy->at(dau1index_))>0.045 or fabs(tree->dxy->at(dau2index_))>0.045) continue;//Dxy
				if (fabs(tree->dz->at(dau1index_))>0.2 or fabs(tree->dz->at(dau2index_))>0.2) continue;//Dz
				if (!((tree->daughters_muonID->at(dau1index_) >> 2) & 1)) continue;
				if (!((tree->daughters_muonID->at(dau2index_) >> 2) & 1)) continue;

			}
			else if(isMT){
				if (fabs(tree->dxy->at(dau1index_))>0.045) continue;//Dxy
				if (fabs(tree->dz->at(dau1index_))>0.2 or fabs(tree->dz->at(dau2index_))>0.2) continue;//Dz
				if (!((tree->daughters_muonID->at(dau1index_) >> 2) & 1)) continue;
				if (tree->daughters_decayModeFindingOldDMs->at(dau2index_)==0) continue;
                                
                        }
			if (fabs(dau1.Eta())>2.5) continue;
			if (fabs(dau2.Eta())>2.5) continue;

                        if(print) cout << "passed eta selection " << endl;

			pairindex=iMoth; 
			dau1index=dau1index_;
			dau2index=dau2index_;
			isPair=true;
		}

		if (!isPair) continue;

                
                metx=tree->METx->at(pairindex);
		mety=tree->METy->at(pairindex);
		TLorentzVector tau1;
		TLorentzVector tau2;
		tau1.SetPxPyPzE(tree->daughters_px->at(dau1index),tree->daughters_py->at(dau1index),tree->daughters_pz->at(dau1index),tree->daughters_e->at(dau1index));
		tau2.SetPxPyPzE(tree->daughters_px->at(dau2index),tree->daughters_py->at(dau2index),tree->daughters_pz->at(dau2index),tree->daughters_e->at(dau2index));
		TLorentzVector met(metx,mety,0,0);

                cout <<"cand1 = " <<tau1.Px()<<"/"<<tau1.Py()<<"/"<<tau1.Pz()<<"/"<<tau1.E()<<endl;
                cout <<"cand2 = " <<tau2.Px()<<"/"<<tau2.Py()<<"/"<<tau2.Pz()<<"/"<<tau2.E()<<endl;

		fillTree(Run_Tree,iEntry,pairindex,dau1index,dau2index,met,isData,htauid,isMC,pu);
	}//loop over events

	plotFill("nbevt",2,9,0,9,(weight_plus-weight_minus));
	fout->cd();
	Run_Tree->Write();
	map<string, TH1F*>::const_iterator iMap1 = myMap1->begin();
	map<string, TH1F*>::const_iterator jMap1 = myMap1->end();
	for (; iMap1 != jMap1; ++iMap1)
		nplot1(iMap1->first)->Write();
	map<string, TH2F*>::const_iterator iMap2 = myMap2->begin();
	map<string, TH2F*>::const_iterator jMap2 = myMap2->end();
	for (; iMap2 != jMap2; ++iMap2)
		nplot2(iMap2->first)->Write();

	//fout->Write();
	fout->Close();
	return 0;
}

