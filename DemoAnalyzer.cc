// -*- C++ -*-
//
// Package:    Demo/DemoAnalyzer
// Class:      DemoAnalyzer
// 
/**\class DemoAnalyzer DemoAnalyzer.cc Demo/DemoAnalyzer/plugins/DemoAnalyzer.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Yongho Jeong
//         Created:  Thu, 04 Feb 2016 16:32:47 GMT
//
//


// system include files
#include <memory>
#include <iostream>
#include <vector>
#include <map>

// user include files
// -------------------------FWCore
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ModuleFactory.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
//-----------------------CommonTools
#include "CommonTools/UtilAlgos/interface/TFileService.h"
//----------------------DataFormats
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
//#include "DataFormats/TrackingRecHit/interface/TrackingRecoHitFwd.h"

#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"
#include "DataFormats/MuonReco/interface/MuonEnergy.h"
#include "DataFormats/MuonReco/interface/MuonTime.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/Lepton.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "SimMuon/MCTruth/interface/MuonAssociatorByHits.h"
#include "SimMuon/MCTruth/plugins/MuonAssociatorEDProducer.h"
#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include <string>
#include <sstream>
#include <TH1.h>
#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <TMath.h>
#include <TLorentzVector.h>

using namespace edm;
using namespace std;
using namespace reco;
using namespace pat;

// class declaration

const int nMax=24000;//over than cmsRun events
int tnm,nm,ng,ngp;
//int met;
float Pmu_pt[nMax], Pmu_eta[nMax], Tgmu_pt[nMax], Tgmu_eta[nMax], Tmet_pt[nMax], Ttr_pt[nMax], Ttr_eta[nMax];
float mllpm, mllpmT, mllpmP, tr_tagPat_mllpm, tr_allproPat_mllpm, tr_proPat_mllpm;
float Tgen_pt[nMax], Tgen_eta[nMax];

class DemoAnalyzer : public edm::EDAnalyzer {
	public:
		//Constructor
		explicit DemoAnalyzer(const edm::ParameterSet&);
		//Destructor
		~DemoAnalyzer();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		TFile *demo_tree;
		TH1D *demohisto;
		TH1D *hpmu_pt, *hpmu_eta, *htr_pt, *htot_tr_pt, *htot_pmu_pt, *htr_eta, *hmllpm, *hmllpmT, *hDeltaT, *hDeltaP, *hmllpmP;
		TH1D *hmllpmT_D, *hmllpmP_D;
		TH1D *hgmu_pt, *hgmu_eta, *htot_gmu_pt; 
		//TH1D *hmet_pt, *hmet_eta;

		TH1D *hpat_g_pt, *hpat_g_eta;
		TH1D *hpat_g_nochar_pt, *hpat_g_nochar_eta;
		TH1D *htr_pro_P, *htr_pro_M, *htr_pro_P_tot, *htr_pro_M_tot;//Track Tag&Probe
		TH1D *htr_div_ProP, *htr_div_ProM;

		TH1D *hpmu_pro_P, *hpmu_pro_M, *hpmu_pro_P_tot, *hpmu_pro_M_tot;
		TH1D *hgmu_pro_P, *hgmu_pro_M, *hgmu_pro_P_tot, *hgmu_pro_M_tot;
		TH1D *hpmu_pro, *hpmu_pro_tot;
		TH1D *hgmu_pro, *hgmu_pro_tot;
		TH1D *hgmu_div_Pro;

		TH1D *hpmu_div_ProP, *hpmu_div_ProM, *hpmu_div_Pro;//PatMuon Tag&Probe
		TH1D *hgmu_div_ProP, *hgmu_div_ProM;
		TH1D *hgpmu_div_pt, *hgpmu_div_eta;
		//------------------------------Track Isolation-----------------------
		TH1D *htrIso_P, *htrIso_M, *htrIso_tagPat_P, *htrIso_tagPat_M, *htrIso_proPat_P, *htrIso_proPat_M;
		TH1D *htrIso_proPat_pt_iso, *htrIso_proPat_eta_iso, *htrIso_proPat_pt_P_iso, *htrIso_proPat_eta_P_iso, *htrIso_proPat_pt_M_iso, *htrIso_proPat_eta_M_iso;
		TH1D *htrIso_allproPat_P, *htrIso_allproPat_M;

		TH1D *htrIso_allproPat_pt_P_nocut, *htrIso_allproPat_eta_P_nocut, *htrIso_allproPat_pt_M_nocut, *htrIso_allproPat_eta_M_nocut, *htrIso_allproPat_pt_nocut, *htrIso_allproPat_eta_nocut;

		TH1D *htrIso_allproPat_pt_P_iso, *htrIso_allproPat_eta_P_iso, *htrIso_allproPat_pt_M_iso, *htrIso_allproPat_eta_M_iso, *htrIso_allproPat_pt_iso, *htrIso_allproPat_eta_iso;

		TH1D *htrIso_proPat_pt_nocut, *htrIso_proPat_eta_nocut, *htrIso_proPat_pt_P_nocut, *htrIso_proPat_eta_P_nocut, *htrIso_proPat_pt_M_nocut, *htrIso_proPat_eta_M_nocut;

		TH1D *htrIso_tag_allproPat_mllpm, *htrIso_tag_proPat_mllpm, *htrIso_tag_failproPat_mllpm, *htrIso_div_mllpm;

		TH1D *htrIso_tag_allproPat_pt_iso, *htrIso_tag_allproPat_eta_iso, *htrIso_tag_allproPat_pt_mi, *htrIso_tag_allproPat_eta_mi;
		TH1D *htrIso_tag_allproPat_pt_nocut, *htrIso_tag_allproPat_eta_nocut, *htrIso_tag_allproPat_pt_mass, *htrIso_tag_allproPat_eta_mass;
		TH1D *htrIso_tag_proPat_pt_iso, *htrIso_tag_proPat_eta_iso, *htrIso_tag_proPat_pt_mi, *htrIso_tag_proPat_eta_mi;
		TH1D *htrIso_tag_proPat_pt_nocut, *htrIso_tag_proPat_eta_nocut, *htrIso_tag_proPat_pt_mass, *htrIso_tag_proPat_eta_mass;

		TH1D *htrIso_div_tag_pt_nocut, *htrIso_div_tag_pt_mass, *htrIso_div_tag_pt_mi, *htrIso_div_tag_pt_iso;
		TH1D *htrIso_div_tag_eta_nocut, *htrIso_div_tag_eta_mass, *htrIso_div_tag_eta_mi, *htrIso_div_tag_eta_iso;

		//------------------------------PAt Muon Isolation------------------------
		TH1D *hpatIso_P, *hpatIso_M;
		//-------------------------------Global Muon Isolation---------------------
		TH1D *hgmuIso_P, *hgmuIso_M, *hgmuIso_pt_P, *hgmuIso_pt_M, *hgmuIso_eta_P, *hgmuIso_eta_M, *hgmuIso_sum_P, *hgmuIso_sum_M;
		TH1D *hgen_pt, *hgen_eta;

		TTree *tree;
		unsigned int minTracks_;
		//bool pro_muM;
		int pro_muP, pro_muM;
		bool pat_g;
		bool SaveHisto;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
DemoAnalyzer::DemoAnalyzer(const edm::ParameterSet& iConfig):
	minTracks_(iConfig.getUntrackedParameter<unsigned int>("minTracks",0))
	//	minTracks_(iConfig.getUntrackedParameter<edm::InputTag>("minTracks"))
{
	//now do what ever initialization is needed
	SaveHisto=iConfig.getParameter<bool>("SaveHisto");

	edm::Service<TFileService> fs;
	if(SaveHisto)demo_tree = new TFile("PAT_tree.root","recreate");
	//	demohisto = fs->make<TH1D>("tracks","Tracks",100,0,5000)
	demohisto = new TH1D("tracks","Tracks",100,0,5000);

	hpmu_pt = new TH1D("hpmu_pt","Muon_P_{t}",25,0,250);
	hpmu_eta = new TH1D("hpmu_eta","Muon_{#eta}",25,-2.4,2.4);
	htot_pmu_pt = new TH1D("htot_pmu_pt","ToTPatMuonCharge P_{t}",20,0,200);

	htr_pt = new TH1D("htr_pt","Track_P_{t}",100,0,200);
	htr_eta = new TH1D("htr_eta","Track_{#eta}",100,-2.4,2.4);
	htot_tr_pt = new TH1D("htot_tr_pt","TOTTrackCharge_P_{t}",100,0,200);

	hgmu_pt = new TH1D("hgmu_pt","GlobalMu_P_{t}",25,0,250);
	hgmu_eta = new TH1D("hgmu_eta","GlobalMu_{#eta}",25,-2.4,2.4);
	htot_gmu_pt = new TH1D("htot_gmu_pt","TOTGlobalMuonCharge P_{t}",100,0,180);

	//----------------------------------------------------
	hmllpm = new TH1D("hmllpm","Di_Muon_{m}",100,0,140);
	hmllpmT = new TH1D("hmllpmT","Di_Track_{mass}",100,0,140);
	hmllpmP = new TH1D("hmllpmP","Di_Muon_{mass}",100,0,140);
	hmllpmT_D = new TH1D("hmllpmT_D","Track_Mass_{dR<0.01}",100,0,140);
	hmllpmP_D = new TH1D("hmllpmP_D","Muon_Mass_{dR<0.01}",100,0,140);

	hDeltaP = new TH1D("hDeltaP","PatMuon_Delta",100,0,5);
	hDeltaT = new TH1D("hDeltaT","Track_Delta",100,0,5);

	htr_pro_P = new TH1D("htr_pro_P","ProbeTrack_{Pat}_{passed +dR<0.01}",30,0,300);
	htr_pro_M = new TH1D("htr_pro_M","ProbeTrack_{Pat}_{passed -dR<0.01}",30,0,300);
	htr_pro_P_tot = new TH1D("htr_pro_P_tot","ProbeTrack_{Pat}_{+total}",30,0,300);
	htr_pro_M_tot = new TH1D("htr_pro_M_tot","ProbeTrack_{Pat}_{-total}",30,0,300);

	htr_div_ProP = new TH1D("htr_div_ProP","Track_{Pat}_{+efficiency}",30,0,300);
	htr_div_ProM = new TH1D("htr_div_ProM","Track_{Pat}_{-efficiency}",30,0,300);

	hpmu_pro_P = new TH1D("hpmu_pro_P","PatProbeMuon_{Global}_{+pass}",30,0,300);
	hpmu_pro_M = new TH1D("hpmu_pro_M","PatProbeMuon_{Global}_{-pass}",30,0,300);
	hpmu_pro_P_tot = new TH1D("hpmu_pro_P_tot","PatProbeMuon_{Global}_{+total}",30,0,300);
	hpmu_pro_M_tot = new TH1D("hpmu_pro_M_tot","PatProbeMuon_{Global}_{-total}",30,0,300);
	hpmu_div_ProP = new TH1D("hpmu_div_ProP","PatProbeMuon_{Global}_{+efficiency}",30,0,300);
	hpmu_div_ProM = new TH1D("hpmu_div_ProM","PatProbeMuon_{Global}_{-efficiency}",30,0,300);

	hgmu_pro_P = new TH1D("hgmu_pro_P","TrackProbeMuon_{Global}_{+pass}",30,0,300);
	hgmu_pro_M = new TH1D("hgmu_pro_M","TrackProbeMuon_{Global}_{-pass}",30,0,300);
	hgmu_pro_P_tot = new TH1D("hgmu_pro_P_tot","TrackProbeMuon_{Global}_{+total}",30,0,300);
	hgmu_pro_M_tot = new TH1D("hgmu_pro_M_tot","TrackProbeMuon_{Global}_{-total}",30,0,300);
	hgmu_div_ProP = new TH1D("hgmu_div_ProP","TrackProbeMuon_{Global}_{+efficiency}",30,0,300);
	hgmu_div_ProM = new TH1D("hgmu_div_ProM","TrackProbeMuon_{Global}_{-eddiciency}",30,0,300);

	hpmu_pro = new TH1D("hpmu_pro","PatProbeMuon_{Global}_{passed_no_charge}",15,0,300);
	hpmu_pro_tot = new TH1D("hpmu_pro_tot","PatProbeMuon_{Global}_{total_no_charge}",15,0,300);
	hpmu_div_Pro = new TH1D("hpmu_div_Pro","PatNoChargeMuon_{Global}_{efficiency}",15,0,300);

	hgmu_pro = new TH1D("hgmu_pro","TrackProbeMuon_{Global}_{passed_no_charge}",15,0,300);
	hgmu_pro_tot = new TH1D("hgmu_pro_tot","TrackProbeMuon_{Global}_{total_no_charge}",15,0,300);
	hgmu_div_Pro = new TH1D("hgmu_div_Pro","TrackNoChargeMuon_{Global}_{efficiency}",15,0,300);

	hpat_g_pt = new TH1D("hpat_g_pt","GlobalMuonInPatMuon_{Pt}",25,0,250);
	hpat_g_eta = new TH1D("hpat_g_eta","GlobalMuonInPatMuon_{#eta}",25,-2.4,2.4);

	hgpmu_div_pt = new TH1D("hgpmu_div_pt","Muon_{Pt}_{Gmu/Pmu}",25,0,250);
	hgpmu_div_eta = new TH1D("hgpmu_div_eta","Muon_{#eta}_{Gmu/Pmu}",25,-2.4,2.4);
	//------------------------------Track Isolation----------------------------
	htrIso_P = new TH1D("htrIso_P","Track_Isolation_{+}",20,0,1);
	htrIso_M = new TH1D("htrIso_M","Track_Isolation_{-}",20,0,1);
	htrIso_tagPat_P = new TH1D("htrIso_tagPat_P","Track_tag_Isolation_{+}",20,0,1);
	htrIso_tagPat_M = new TH1D("htrIso_tagPat_M","Track_tag_Isolation_{-}",20,0,1);
	htrIso_proPat_P = new TH1D("htrIso_proPat_P","Track_probe_Isolation_{+}",20,0,1);
	htrIso_proPat_M = new TH1D("htrIso_proPat_M","Track_probe_Isolation_{-}",20,0,1);
	htrIso_allproPat_P = new TH1D("htrIso_allproPat_P","Track_allprobe_Isolation_{+}",20,0,1);
	htrIso_allproPat_M = new TH1D("htrIso_allproPat_M","Track_allprobe_Isolation_{-}",20,0,1);

	htrIso_proPat_pt_P_iso = new TH1D("htrIso_proPat_pt_P_iso","Track_probe_Pat_{+p}_{t}_{isopt<0.14}",20,0,300);
	htrIso_proPat_pt_M_iso = new TH1D("htrIso_proPat_pt_M_iso","Track_probe_Pat_{-p}_{t}_{isopt<0.14}",20,0,300);
	htrIso_proPat_eta_P_iso = new TH1D("htrIso_proPat_eta_P_iso","Track_probe_Pat_{+#eta}_{isopt<0.14}",20,-3,3);
	htrIso_proPat_eta_M_iso = new TH1D("htrIso_proPat_eta_M_iso","Track_probe_Pat_{-#eta}_{isopt<0.14}",20,-3,3);
	htrIso_proPat_pt_iso = new TH1D("htrIso_proPat_pt_iso","Track_probe_{P}_{t}_{isopt<0.14}",20,0,300);
	htrIso_proPat_eta_iso = new TH1D("htrIso_proPat_eta_iso","Track_probe_{#eta}_{isopt<0.14}",20,-3,3);

	//----------------------Track nocut Isolation----------------------------
	htrIso_proPat_pt_P_nocut = new TH1D("htrIso_proPat_pt_P_nocut","Track_probe_Pat_{+p}_{t}",20,0,300);
	htrIso_proPat_pt_M_nocut = new TH1D("htrIso_proPat_pt_M_nocut","Track_probe_Pat_{-p}_{t}",20,0,300);
	htrIso_proPat_eta_P_nocut = new TH1D("htrIso_proPat_eta_P_nocut","Track_probe_Pat_{+#eta}",20,-3,3);
	htrIso_proPat_eta_M_nocut = new TH1D("htrIso_proPat_eta_M_nocut","Track_probe_Pat_{-#eta}",20,-3,3);
	htrIso_proPat_pt_nocut = new TH1D("htrIso_proPat_pt_nocut","Track_probe_{P}_{t}",20,0,300);
	htrIso_proPat_eta_nocut = new TH1D("htrIso_proPat_eta_nocut","Track_probe_{#eta}",20,-3,3);

	htrIso_allproPat_pt_P_nocut = new TH1D("htrIso_allproPat_pt_P_nocut","Track_probe_Pat_{+P}_{t}",20,0,300);
	htrIso_allproPat_eta_P_nocut = new TH1D("htrIso_allproPat_eta_P_nocut","Track_probe_Pat_{+#eta}",20,-3,3);
	htrIso_allproPat_pt_M_nocut = new TH1D("htrIso_allproPat_pt_M_nocut","Track_probe_Pat_{-P}_{t}",20,0,300);
	htrIso_allproPat_eta_M_nocut = new TH1D("htrIso_allproPat_eta_M_nocut","Track_probe_Pat_{-#eta}",20,-3,3);
	htrIso_allproPat_pt_nocut = new TH1D("htrIso_allproPat_pt_nocut","Track_probe_Pat_{Pt}",20,0,300);
	htrIso_allproPat_eta_nocut = new TH1D("htrIso_allproPat_eta_nocut","Track_probe_Pat_{#eta}",20,-3,3);

	htrIso_allproPat_pt_P_iso = new TH1D("htrIso_allproPat_pt_P_iso","Track_probe_Pat_{+P}_{t} [iso<0.14]",20,0,300);
	htrIso_allproPat_eta_P_iso = new TH1D("htrIso_allproPat_eta_P_iso","Track_probe_Pat_{+#eta} [iso<0.14]",20,-3,3);
	htrIso_allproPat_pt_M_iso = new TH1D("htrIso_allproPat_pt_M_iso","Track_probe_Pat_{-P}_{t} [iso<0.14]",20,0,300);
	htrIso_allproPat_eta_M_iso = new TH1D("htrIso_allproPat_eta_M_iso","Track_probe_Pat_{-#eta} [iso<0.14] []",20,-3,3);
	htrIso_allproPat_pt_iso = new TH1D("htrIso_allproPat_pt_iso","Track_probe_Pat_{Pt} [iso<0.14]}",20,0,300);
	htrIso_allproPat_eta_iso = new TH1D("htrIso_allproPat_eta_iso","Track_probe_Pat_{#eta} [iso<0.14]",20,-3,3);
	//-----------------------Track Isolation efficiency---------------------------------------------

	//--------------------------htrIso_tag_mass distribution-----------------------------------
	htrIso_tag_allproPat_mllpm = new TH1D("htrIso_tag_allproPat_mllpm","Iso_Tag+AllPro_{mass}",50,5,140);
	htrIso_tag_allproPat_pt_iso = new TH1D("htrIso_tag_allproPat_pt_iso","Iso_Tag+Allpro_{pt}_{iso}",20,0,300);
	htrIso_tag_allproPat_eta_iso = new TH1D("htrIso_tag_allproPat_eta_iso","Iso_Tag+Allpro_{#eta}_{iso}",20,-3,3);
	htrIso_tag_allproPat_pt_mi = new TH1D("htrIso_tag_allproPat_pt_mi","Iso_Tag+Allpro_{pt}_{iso,mass}",20,0,300);
	htrIso_tag_allproPat_eta_mi = new TH1D("htrIso_tag_allproPat_eta_mi","Iso_Tag+Allpro_{#eta}_{iso,mass}",20,-3,3);

	htrIso_tag_allproPat_pt_nocut = new TH1D("htrIso_tag_allproPat_pt_nocut","Iso_Tag+Allpro_{pt}",20,0,300);
	htrIso_tag_allproPat_eta_nocut = new TH1D("htrIso_tag_allproPat_eta_nocut","Iso_Tag+Allpro_{#eta}",20,-3,3);
	htrIso_tag_allproPat_pt_mass = new TH1D("htrIso_tag_allproPat_pt_mass","Iso_Tag+Allpro_{pt}_{mass}",20,0,300);
	htrIso_tag_allproPat_eta_mass = new TH1D("htrIso_tag_allproPat_eta_mass","Iso_Tag+Allpro_{#eta}_{mass}",20,-3,3);

	htrIso_tag_proPat_mllpm = new TH1D("htrIso_tag_proPat_mllpm","Iso_Tag_Pro_{mass}",50,5,140);
	htrIso_tag_proPat_pt_iso = new TH1D("htrIso_tag_proPat_pt_iso","Iso_Tag_Pro_{pt}_{iso}",20,0,300);
	htrIso_tag_proPat_eta_iso = new TH1D("htrIso_tag_proPat_eta_iso","Iso_Tag+pro_{#eta}_{iso}",20,-3,3);
	htrIso_tag_proPat_pt_mi = new TH1D("htrIso_tag_proPat_pt_mi","Iso_Tag+pro_{pt}_{iso,mass}",20,0,300);
	htrIso_tag_proPat_eta_mi = new TH1D("htrIso_tag_proPat_eta_mi","Iso_Tag+pro_{#eta}_{iso,mass}",20,-3,3);

	htrIso_tag_proPat_pt_nocut = new TH1D("htrIso_tag_proPat_pt_nocut","Iso_Tag_Pro_{pt}",20,0,300);
	htrIso_tag_proPat_eta_nocut = new TH1D("htrIso_tag_proPat_eta_nocut","Iso_Tag+pro_{#eta}",20,-3,3);
	htrIso_tag_proPat_pt_mass = new TH1D("htrIso_tag_proPat_pt_mass","Iso_Tag+pro_{pt}_{mass}",20,0,300);
	htrIso_tag_proPat_eta_mass = new TH1D("htrIso_tag_proPat_eta_mass","Iso_Tag+pro_{#eta}_{mass}",20,-3,3);

	htrIso_tag_failproPat_mllpm = new TH1D("htrIso_tag_failproPat_mllpm","Iso_Tag_failPro_distribution_{mass}",50,5,140);
	htrIso_div_mllpm = new TH1D("htrIso_div_mllpm","Isolation_mass_{efficiency_all/passed probe}",50,5,140);

	htrIso_div_tag_pt_nocut = new TH1D("htrIso_div_tag_pt_nocut","eff_pro/allpro_Pt_{nocut}",20,0,300);
	htrIso_div_tag_pt_iso = new TH1D("htrIso_div_tag_pt_iso","eff_pro/allpro_Pt_{iso}",20,0,300);
	htrIso_div_tag_pt_mi = new TH1D("htrIso_div_tag_pt_mi","eff_pro/allpro_Pt_{mass,iso}",20,0,300);
	htrIso_div_tag_pt_mass = new TH1D("htrIso_div_tag_pt_mass","eff_pro/allpro_Pt_{mass}",20,0,300);

	htrIso_div_tag_eta_nocut = new TH1D("htrIso_div_tag_eta_nocut","eff_pro/allpro_#eta_{nocut}",20,-3,3);
	htrIso_div_tag_eta_iso = new TH1D("htrIso_div_tag_eta_iso","eff_pro/allpro_#eta_{iso}",20,-3,3);
	htrIso_div_tag_eta_mi = new TH1D("htrIso_div_tag_eta_mi","eff_pro/allpro_#eta_{mass,iso}",20,-3,3);
	htrIso_div_tag_eta_mass = new TH1D("htrIso_div_tag_eta_mass","eff_pro/allpro_#eta_{mass}",20,-3,3);

	//----------------------------------------PAtMuon Isolation----------------------------
	hpatIso_P = new TH1D("hpatIso_P","PatMuon_Isolation_{+}",20,0,1);
	hpatIso_M = new TH1D("hpatIso_M","PatMuon_Isolation_{-}",20,0,1);
	//-----------------------------------------GlobalMuon ISolation -------------------------------
	hgmuIso_P = new TH1D("hgmuIso_P","GlobalMuon_Isolation_{+}",20,0,1);
	hgmuIso_M = new TH1D("hgmuIso_M","GlobalMuon_Isolation_{-}",20,0,1);
	hgmuIso_pt_P = new TH1D("hgmuIso_pt_P","GlobalMuon_Isolation_{+P}_{t}/{Isopt<0.14}",20,0,300);
	hgmuIso_eta_P = new TH1D("hgmuIso_eta_P","GlobalMuon_Isolation_{+#eta}/{Isopt<0.14}",20,-3,3);
	hgmuIso_pt_M = new TH1D("hgmuIso_pt_M","GlobalMuon_Isolation_{-P}_{t}/{Isopt<0.14}",20,0,300);
	hgmuIso_eta_M = new TH1D("hgmuIso_eta_M","GlobalMuon_Isolation_{-#eta}/{Isopt<0.14}",20,-3,3);
	hgmuIso_sum_P = new TH1D("hgmuIso_sum_P","GlobalMuon_Iso_Sum_{+Pt}",20,0,50);
	hgmuIso_sum_M = new TH1D("hgmuIso_sum_M","GlobalMuon_Iso_sum_{-Pt}",20,0,50);
	//---------------------------------------------GenParticle------------------------------------
	hgen_pt = new TH1D("hgen_pt","GenMuon_P_{t}",20,0,300);
	hgen_eta = new TH1D("hgen_eta","GenMuon_#eta",20,-3,3);

	//-------------------------------------------tree Branch-----------------------------------	
	tree = new TTree("tree","example_tree");	
	tree->Branch("nm",&nm,"nm/I");
	tree->Branch("Pmu_pt",Pmu_pt,"Pmu_pt[nm]/F");
	tree->Branch("Pmu_eta",Pmu_eta,"Pmu_eta[nm]/F");

	tree->Branch("tnm",&tnm,"tnm/I");
	tree->Branch("Ttr_pt",Ttr_pt,"Ttr_pt[tnm]/F");
	tree->Branch("Ttr_eta",Ttr_eta,"Ttr_pt[tnm]/F");

	tree->Branch("ng",&ng,"ng/I");
	tree->Branch("Tgmu_pt",Tgmu_pt,"Tgmu_pt[ng]/F");
	tree->Branch("Tgmu_eta",Tgmu_eta,"Tgmu_eta[ng]/F");

	tree->Branch("ngp",&ngp,"ngp/I");
	tree->Branch("Tgen_pt",Tgen_pt,"Tgen_pt[ngp]/F");
	tree->Branch("Tgen_pt",Tgen_pt,"Tgen_eta[ngp]/F");

	tree->Branch("mllpmT",&mllpmT,"mllpmT/F");
	tree->Branch("mllpm",&mllpm,"mllpm/F");
	tree->Branch("mllpmP",&mllpmP,"mllpmP/F");
	tree->Branch("tr_tagPat_mllpm",&tr_tagPat_mllpm,"tr_tagPat_mllpm/F");
	tree->Branch("tr_allproPat_mllpm",&tr_allproPat_mllpm,"tr_allproPat_mllpm/F");
	tree->Branch("tr_proPat_mllpm",&tr_proPat_mllpm,"tr_proPat_mllpm/F");
}

DemoAnalyzer::~DemoAnalyzer()
{
	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
	void
DemoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	vector<float> vpt;

/*	Handle<reco::VertexCollection> Ver;
	iEvent.getByLabel("offlinePrimaryVertices",Ver);
	if(Ver->empty()){cout<<"no PV"<<endl;return;}
	auto pv0 = Ver->front();
*/
	Handle<reco::GenParticleCollection>Gen;
	iEvent.getByLabel("genParticles",Gen);
	double minDRP=100;//set variable of minimum "DeltaR" between two particle.
	double minDRM=100;//set variable of minimum "DeltaR" between two particle.
	unsigned int ngen=0;
	int nGenP=-1, nGenM=-1;
	double geptP=0,geptM=0;
	for(auto gen=Gen->begin();gen!=Gen->end();gen++,ngen++)
	{
		int ndau = gen->numberOfDaughters();
		int id = gen->pdgId();
		int gstat = gen->status();
		if(fabs(gen->pt())>10 && fabs(gen->eta())<2.4 && abs(id)==13 && abs(ndau)!=0 && abs(gstat)!=1)//or status>1
		{
			vpt.push_back(gen->pt());
			if(abs(gen->mother(0)->pdgId())==23)
			{
				if(gen->charge()==1) nGenP=ngen;
				if(gen->charge()==-1) nGenM=ngen;
			}
			if(gen->charge()==1 && geptP<gen->pt()) geptP=gen->pt();
			if(gen->charge()==-1 && geptM<gen->pt()) geptM=gen->pt();
			hgen_pt->Fill(gen->pt());
			hgen_eta->Fill(gen->eta());
		}
		Tgen_pt[ngen] = gen->pt();
		Tgen_eta[ngen] = gen->eta();
	}
	ngp=ngen;

	TLorentzVector GenP, GenM;
	if(nGenP!=-1 && nGenM!=-1)
	{
		//cout<<nGenP<<endl;
		GenP.SetPxPyPzE((*Gen)[nGenP].px(),(*Gen)[nGenP].py(),(*Gen)[nGenP].pz(),(*Gen)[nGenP].energy());
		GenM.SetPxPyPzE((*Gen)[nGenM].px(),(*Gen)[nGenM].py(),(*Gen)[nGenM].pz(),(*Gen)[nGenM].energy());
	}

	reco::CandidateCollection tempLeps;
	//--------------------------------------------------Track section----------------------------------
	Handle<reco::TrackCollection> Tracks;
	iEvent.getByLabel("generalTracks",Tracks);

	if (minTracks_ <=Tracks->size()) {
		LogInfo("Demo")<<"number of Tracks"<<Tracks->size();
	}
	demohisto->Fill(Tracks->size());

	//--------Muon section

	unsigned int ntr=0;//only variable for Collection for roop
	int ntrP=-1, ntrM=-1;
	double trptP=0, trptM=0;

	for(reco::TrackCollection::const_iterator tr=Tracks->begin();tr!=Tracks->end()-1;tr++,ntr++)
	{
		if(fabs(tr->pt())>10&&fabs(tr->eta())<2.4)
		{
			vpt.push_back(tr->pt());
			//cout<<"track_test2"<<endl;
			if(tr->charge()==1&&trptP<tr->pt())//select only one 'best track'
			{
				ntrP=ntr;
				trptP=tr->pt();
			}
			if(tr->charge()==-1&&trptM<tr->pt())//select only one 'best track'
			{
				ntrM=ntr;
				trptM=tr->pt();
				//			cout<<(trptM=tr->pt())<<endl;
			}
			htr_pt->Fill(tr->pt());
			htr_eta->Fill(tr->eta());
			Ttr_pt[ntr] = tr->pt();
			Ttr_eta[ntr] = tr->eta();
		}
	}
	tnm=ntr;

	unsigned int tot_ntr=0;
	if(trptP>10&&trptM>10)
	{
		for(reco::TrackCollection::const_iterator tot_tr=Tracks->begin();tot_tr!=Tracks->end();tot_tr++,tot_ntr++)
		{
			if(fabs(tot_tr->pt())>10)
			{	
				vpt.push_back(tot_tr->pt());
				if(tot_tr->charge()==1)
				{
					trptP=tot_tr->pt();
					ntrP=tot_ntr;
				}
				if(tot_tr->charge()==-1)
				{	
					trptM=tot_tr->pt();
					ntrM=tot_ntr;
				}
				htot_tr_pt->Fill(tot_tr->pt());
			}
		}
	}

	TLorentzVector trP, trM;//for track "tag" LorentzVector
	if(ntrP!=-1&&ntrM!=-1)
	{
		trP.SetPxPyPzE((*Tracks)[ntrP].px(),(*Tracks)[ntrP].py(),(*Tracks)[ntrP].pz(),(*Tracks)[ntrP].p());
		trM.SetPxPyPzE((*Tracks)[ntrM].px(),(*Tracks)[ntrM].py(),(*Tracks)[ntrM].pz(),(*Tracks)[ntrM].p());
		if(fabs((trP+trM).M())>5)
		{
			//cout<<((trP+trM).M())<<endl;		
			hmllpmT->Fill((trP+trM).M());
		}
		mllpmT=(trP+trM).M();
	}

	TLorentzVector trP_P, trM_P;//for track "probe" LorentzVector
	if(ntrP!=-1&&ntrM!=-1)//there is no 'a=b' form in 'if' sentence
	{
		trP_P.SetPxPyPzE((*Tracks)[ntrP].px(),(*Tracks)[ntrP].py(),(*Tracks)[ntrP].pz(),(*Tracks)[ntrP].p());
		trM_P.SetPxPyPzE((*Tracks)[ntrM].px(),(*Tracks)[ntrM].py(),(*Tracks)[ntrM].pz(),(*Tracks)[ntrM].p());
	}
	//if(fabs((trP_P+trM_P).M())>5)htrIso_tagPat_mllpm->Fill((trP_P+trM_P).M());

	Handle<pat::MuonCollection> Pmu;
	iEvent.getByLabel("selectedPatMuons",Pmu);

	reco::CandidateCollection tempLepsP;
	unsigned int pnmu=0;
	int PnmuP=-1,PnmuM=-1;
	double PmuptP=0,PmuptM=0;

	for(pat::MuonCollection::const_iterator pmu=Pmu->begin();pmu!=Pmu->end();pmu++,pnmu++)
	{
		if(fabs(pmu->eta())<2.4&&fabs(pmu->pt())>10)
		{
			vpt.push_back(pmu->pt());
			tempLepsP.push_back(*pmu);
			if(pmu->charge()==1&&PmuptP<pmu->pt())//select only one Muon
			{
				PnmuP=pnmu;//events numbering is global variable
				PmuptP=pmu->pt();
			}
			if(pmu->charge()==-1&&PmuptM<pmu->pt())
			{
				PnmuM=pnmu;//events numbering is global variable
				PmuptM=pmu->pt();
			}
			hpmu_pt->Fill(pmu->pt());
			hpmu_eta->Fill(pmu->eta());
		}
		Pmu_pt[pnmu] = pmu->pt();
		Pmu_eta[pnmu] = pmu->eta();

		if(pmu->isGlobalMuon()) pat_g=true;
		else pat_g=false;
		if(pat_g==true) hpat_g_pt->Fill(pmu->pt());
		if(pat_g==true) hpat_g_eta->Fill(pmu->eta());
	}
	nm=pnmu;

	unsigned int tot_pnmu=0;
	if(PmuptP>10&&PmuptM>10)
	{
		for(pat::MuonCollection::const_iterator tot_pmu=Pmu->begin();tot_pmu!=Pmu->end();tot_pmu++,tot_pnmu++)
		{
			if(fabs(tot_pmu->pt())>10)
			{
				vpt.push_back(tot_pmu->pt());
				if(tot_pmu->charge()==1)
				{
					PmuptP=tot_pmu->pt();
					PnmuP=tot_pnmu;
				}
				if(tot_pmu->charge()==-1)
				{
					PmuptM=tot_pmu->pt();
					PnmuM=tot_pnmu;
				}
				htot_pmu_pt->Fill(tot_pmu->pt());
			}
		}
	}

	TLorentzVector PmuP,PmuM;//for "tag" LorentzVector
	if(PnmuP!=-1&&PnmuM!=-1)
	{
		PmuP.SetPxPyPzE((*Pmu)[PnmuP].px(),(*Pmu)[PnmuP].py(),(*Pmu)[PnmuP].pz(),(*Pmu)[PnmuP].p());
		PmuM.SetPxPyPzE((*Pmu)[PnmuM].px(),(*Pmu)[PnmuM].py(),(*Pmu)[PnmuM].pz(),(*Pmu)[PnmuM].p());
		if(fabs((PmuP+PmuM).M())>5)
		{
			//cout<<((PmuP+PmuM).M())<<endl;
			hmllpmP->Fill((PmuP+PmuM).M());
		}
		mllpmP=(PmuP+PmuM).M();
	}

	if(trP.DeltaR(PmuP)<0.01&&trM.DeltaR(PmuM)<0.01)
	{
		if((trP+trM).M()>5&&(PmuP+PmuM).M()>5)
		{
			//cout<<((trP+trM).M())<<endl;
			hmllpmT_D->Fill((trP+trM).M());
			//hDeltaT->Fill(trP.DeltaR(PmuP));

			//cout<<((PmuP+PmuM).M())<<endl;
			hmllpmP_D->Fill((PmuP+PmuM).M());
			//hDeltaP->Fill(trM.DeltaR(PmuM));
		}
	}

	if(trP.DeltaR(trM)&&PmuP.DeltaR(PmuM))//for searching PP->Z->m+m-
	{
		//cout<<(trP.DeltaR(trM))<<endl;
		hDeltaT->Fill(trP.DeltaR(trM));
		//cout<<(PmuP.DeltaR(PmuM))<<endl;
		hDeltaP->Fill(PmuP.DeltaR(PmuM));
	}

	//---------------------------------select 'isopt' of Track--------------------------------------
	//---------------------------------------plus charge----------------------------------------------
	for(reco::TrackCollection::const_iterator tr=Tracks->begin();tr!=Tracks->end();tr++)
	{
		TLorentzVector t1_P, t2_P;
		double iso_trptP=0;
		if(fabs(tr->pt())>10)
		{
			t1_P.SetPxPyPzE(tr->px(),tr->py(),tr->pz(),tr->p());
			vpt.push_back(tr->pt());
			double isopt=0,sumpt=0;
			if(tr->charge()==1)
			{
				for(reco::TrackCollection::const_iterator tr2=Tracks->begin();tr2!=Tracks->end();tr2++)
				{
					t2_P.SetPxPyPzE(tr2->px(),tr2->py(),tr2->pz(),tr2->p());
					if(t1_P.DeltaR(t2_P)<0.3 && t1_P.DeltaR(t2_P)>0.01)
					{
						sumpt+=tr2->pt();
					}
				}
				isopt=sumpt/tr->pt();
				if(iso_trptP<tr->pt());
				{
					htrIso_P->Fill(isopt);
					//cout<<(isopt)<<endl;
					iso_trptP=tr->pt();
				}
			}
		}
	}
	//--------------------------------------minus charge--------------------------------------------------
	for(reco::TrackCollection::const_iterator tr=Tracks->begin();tr!=Tracks->end();tr++)
	{
		TLorentzVector t1_M, t2_M;
		double iso_trptM=0;
		if(fabs(tr->pt())>10)
		{
			t1_M.SetPxPyPzE(tr->px(),tr->py(),tr->pz(),tr->p());
			vpt.push_back(tr->pt());
			double isopt=0,sumpt=0;
			if(tr->charge()==-1)
			{
				for(reco::TrackCollection::const_iterator tr2=Tracks->begin();tr2!=Tracks->end();tr2++)
				{
					t2_M.SetPxPyPzE(tr2->px(),tr2->py(),tr2->pz(),tr2->p());
					if(t1_M.DeltaR(t2_M)<0.3 && t1_M.DeltaR(t2_M)>0.01)
					{
						sumpt+=tr2->pt();
					}
				}
				isopt=sumpt/tr->pt();
				if(iso_trptM<tr->pt())
				{
					htrIso_M->Fill(isopt);
					//cout<<(isopt)<<endl;
					iso_trptM=tr->pt();
				}
			}
		}
	}

	//---------------------------------Tag & Prove Track with PatMuon-------------------------------
	TLorentzVector tag_trP,iso_trP,tag_trM,iso_trM;//-------the LorentzVector for 'tag'
	double trptP_P=0, trptM_P=0;
	double isopt_tr_tag_P=0, isopt_tr_tag_M=0;
	for(reco::TrackCollection::const_iterator tr=Tracks->begin();tr!=Tracks->end();tr++)
	{
		if(fabs(tr->pt())>10)
		{
			vpt.push_back(tr->pt());
			double isopt=0,sumpt=0;//select iso pt of 'tag'
			if(tr->charge()==1)
			{
				tag_trP.SetPxPyPzE(tr->px(),tr->py(),tr->pz(),tr->p());
				double dR = tag_trP.DeltaR(PmuP);

				for(reco::TrackCollection::const_iterator tr2=Tracks->begin();tr2!=Tracks->end();tr2++)
				{
					iso_trP.SetPxPyPzE(tr2->px(),tr2->py(),tr2->pz(),tr2->p());
					if(tag_trP.DeltaR(iso_trP)<0.3 && tag_trP.DeltaR(iso_trP)>0.01)
					{
						//sumpt-=tr->pt();
						sumpt+=tr2->pt();
					}
				}
				isopt=sumpt/tr->pt();//include muon isolation near target 'tag' and 'probe' and except own 'isopt'
				isopt=isopt_tr_tag_P;
				if(trptP_P<tr->pt() && minDRP>dR && dR<0.01)//tag+ is completely Muon
				{
					htrIso_tagPat_P->Fill(isopt);
					//cout<<(isopt)<<endl;
					trptP_P=tr->pt();
					minDRP=dR;
					trP=tag_trP;//tag_trP is "tag" candidate and save "trP"
					//htrIso_tagPat_mllpm->Fill((trP).M());
				}
			}
			if(tr->charge()==-1)
			{
				tag_trM.SetPxPyPzE(tr->px(),tr->py(),tr->pz(),tr->p());
				double dR = tag_trM.DeltaR(PmuM);
				for(reco::TrackCollection::const_iterator tr2=Tracks->begin();tr2!=Tracks->end();tr2++)
				{
					iso_trM.SetPxPyPzE(tr2->px(),tr2->py(),tr2->pz(),tr2->p());
					if(tag_trM.DeltaR(iso_trM)<0.3 && tag_trM.DeltaR(iso_trM)>0.01)
					{
						sumpt+=tr2->pt();
					}
				}
				isopt=(sumpt/tr->pt());//include muon isolation near target 'tag' and except own 'isopt'
				isopt=isopt_tr_tag_M;
				if(trptM_P<tr->pt() && minDRM>dR && dR<0.01)//tag- is completely muon
				{
					htrIso_tagPat_M->Fill(isopt);
					//cout<<(isopt)<<endl;
					trptM_P=tr->pt();
					minDRM=dR;
					trM=tag_trM;//tag_trM 0s saved "tag" candidate
					//htrIso_tagPat_mllpm->Fill((trM).M());
				}
			}
		}
	}
	if(fabs((trM+trP).M())>5)tr_tagPat_mllpm = (trM+trP).M();//tag for tree

	//cout<<(trP.Pt())<<endl;
	//cout<<(tag_trP.Pt())<<endl;
	//----------------------------Prove-----------------------------
	TLorentzVector probe_trP,probe_trM;//LorantzVector of --probe+, probe+
	double isopt_tr_pro_P=0,isopt_tr_pro_M=0;//------------------------the Isolation for probe--------------------
	for(reco::TrackCollection::const_iterator tr=Tracks->begin();tr!=Tracks->end();tr++)
	{
		if(fabs(tr->pt())>10)
		{
			vpt.push_back(tr->pt());
			if(tr->charge()==-1)
			{
				double isopt=0, sumpt=0;
				probe_trM.SetPxPyPzE(tr->px(),tr->py(),tr->pz(),tr->p());
				for(reco::TrackCollection::const_iterator tr2=Tracks->begin();tr2!=Tracks->end();tr2++)
				{
					iso_trM.SetPxPyPzE(tr2->px(),tr2->py(),tr2->pz(),tr2->p());
					if(probe_trM.DeltaR(iso_trM)<0.3 && probe_trM.DeltaR(iso_trM)>0.01)
					{
						sumpt+=tr2->pt();
					}
				}
				isopt=(sumpt/tr->pt());
				isopt_tr_pro_M=isopt;

				double dR=trP.DeltaR(probe_trM);//probe candidate
				if(trptM_P<tr->pt() && dR>1)
				{
					minDRM=dR;
					trptM_P=tr->pt();
					trM_P=probe_trM;//saved "probe" LorentzVector->probe candidate
					htrIso_allproPat_M->Fill(isopt);

					htrIso_allproPat_pt_M_nocut->Fill(tr->pt());
					htrIso_allproPat_eta_M_nocut->Fill(tr->eta());
					htrIso_allproPat_pt_nocut->Fill(tr->pt());
					htrIso_allproPat_eta_nocut->Fill(tr->eta());

					//htrIso_allproPat_mllpm->Fill(probe_trM.M());

					if(isopt<0.14)
					{
						htrIso_allproPat_pt_M_iso->Fill(tr->pt());
						htrIso_allproPat_eta_M_iso->Fill(tr->eta());
						htrIso_allproPat_pt_iso->Fill(tr->pt());
						htrIso_allproPat_eta_iso->Fill(tr->eta());
					}
				}
			}
			if(tr->charge()==1)
			{
				double isopt=0, sumpt=0;
				probe_trP.SetPxPyPzE(tr->px(),tr->py(),tr->pz(),tr->p());
				for(reco::TrackCollection::const_iterator tr2=Tracks->begin();tr2!=Tracks->end();tr2++)
				{
					iso_trP.SetPxPyPzE(tr2->px(),tr2->py(),tr2->pz(),tr2->p());
					if(probe_trP.DeltaR(iso_trP)<0.3 && probe_trP.DeltaR(iso_trP)>0.01)
					{
						sumpt+=tr2->pt();
					}
				}
				isopt=(sumpt/tr->pt());
				isopt_tr_pro_P=isopt;

				double dR=trM.DeltaR(probe_trP);//probe candidate
				if(trptP_P<tr->pt()&&dR>1)
				{
					minDRP=dR;
					trptP_P=tr->pt();
					trP_P=probe_trP;//saved "probe" LorentzVector->probe candidate
					htrIso_allproPat_P->Fill(isopt);

					htrIso_allproPat_pt_P_nocut->Fill(tr->pt());
					htrIso_allproPat_eta_P_nocut->Fill(tr->eta());
					htrIso_allproPat_pt_nocut->Fill(tr->pt());
					htrIso_allproPat_eta_nocut->Fill(tr->eta());

					//htrIso_allproPat_mllpm->Fill(trP_P.M());

					if(isopt<0.14)
					{
						htrIso_allproPat_pt_P_iso->Fill(tr->pt());
						htrIso_allproPat_eta_P_iso->Fill(tr->eta());
						htrIso_allproPat_pt_iso->Fill(tr->pt());
						htrIso_allproPat_eta_iso->Fill(tr->eta());
					}
				}
			}
		}
	}
	tr_allproPat_mllpm = (trP_P+trM_P).M();//all probe for tree
	if(trP_P.DeltaR(PmuP)<0.01 && trM_P.DeltaR(PmuM)<0.01)tr_proPat_mllpm = (trP_P+trM_P).M();//passed probe for tree

	if(isopt_tr_pro_M<0.14 && isopt_tr_pro_P<0.14)
	{
		//------------------------------------------------all probe------------------------------------
		if(fabs((trM+trP_P).M())>5 && fabs((trM+trP_P).M())<140)htrIso_tag_allproPat_mllpm->Fill((trM+trP_P).M());//all probe Z mass peak plus and minus charge. no dot individually each other
		if(fabs((trP+trM_P).M())>5 && fabs((trP+trM_P).M())<140)htrIso_tag_allproPat_mllpm->Fill((trP+trM_P).M());//all probe Z mass peak plus and minus charge. no dot individually each other
		if(fabs((trM+trP_P).M())>5 && fabs((trM+trP_P).M())<140 && trP_P.DeltaR(PmuP)<0.01)htrIso_tag_proPat_mllpm->Fill((trM+trP_P).M());
		if(fabs((trP+trM_P).M())>5 && fabs((trP+trM_P).M())<140 && trM_P.DeltaR(PmuM)<0.01)htrIso_tag_proPat_mllpm->Fill((trP+trM_P).M());
		if(fabs((trM+trP_P).M())>5 && fabs((trM+trP_P).M())<140 && trP_P.DeltaR(PmuP)>=0.01)htrIso_tag_failproPat_mllpm->Fill((trM+trP_P).M());
		if(fabs((trP+trM_P).M())>5 && fabs((trP+trM_P).M())<140 && trM_P.DeltaR(PmuM)>=0.01)htrIso_tag_failproPat_mllpm->Fill((trP+trM_P).M());

		htrIso_tag_allproPat_pt_iso->Fill((trM+trP_P).Pt());
		htrIso_tag_allproPat_pt_iso->Fill((trP+trM_P).Pt());
		htrIso_tag_allproPat_eta_iso->Fill((trM+trP_P).Eta());
		htrIso_tag_allproPat_eta_iso->Fill((trP+trM_P).Eta());
		if(fabs((trP+trM_P).M())>80 && fabs((trP+trM_P).M())<140)htrIso_tag_allproPat_pt_mi->Fill((trP+trM_P).Pt());
		if(fabs((trM+trP_P).M())>80 && fabs((trM+trP_P).M())<140)htrIso_tag_allproPat_pt_mi->Fill((trM+trP_P).Pt());
		if(fabs((trP+trM_P).M())>80 && fabs((trP+trM_P).M())<140)htrIso_tag_allproPat_eta_mi->Fill((trP+trM_P).Eta());
		if(fabs((trM+trP_P).M())>80 && fabs((trM+trP_P).M())<140)htrIso_tag_allproPat_eta_mi->Fill((trM+trP_P).Eta());
		//-------------------------------------------------------probe------------------------------

		if(trP_P.DeltaR(PmuP)<0.01)htrIso_tag_proPat_pt_iso->Fill((trM+trP_P).Pt());
		if(trM_P.DeltaR(PmuM)<0.01)htrIso_tag_proPat_pt_iso->Fill((trP+trM_P).Pt());
		if(trP_P.DeltaR(PmuP)<0.01)htrIso_tag_proPat_eta_iso->Fill((trM+trP_P).Eta());
		if(trM_P.DeltaR(PmuM)<0.01)htrIso_tag_proPat_eta_iso->Fill((trP+trM_P).Eta());
		if(fabs((trP+trM_P).M())>80 && fabs((trP+trM_P).M())<140 && trM_P.DeltaR(PmuM)<0.01)htrIso_tag_proPat_pt_mi->Fill((trP+trM_P).Pt());
		if(fabs((trM+trP_P).M())>80 && fabs((trM+trP_P).M())<140 && trP_P.DeltaR(PmuP)<0.01)htrIso_tag_proPat_pt_mi->Fill((trM+trP_P).Pt());
		if(fabs((trP+trM_P).M())>80 && fabs((trP+trM_P).M())<140 && trM_P.DeltaR(PmuM)<0.01)htrIso_tag_proPat_eta_mi->Fill((trP+trM_P).Eta());
		if(fabs((trM+trP_P).M())>80 && fabs((trM+trP_P).M())<140 && trP_P.DeltaR(PmuP)<0.01)htrIso_tag_proPat_eta_mi->Fill((trM+trP_P).Eta());
	}

	//---------------------------------------all probe---------------------------------------------
	htrIso_tag_allproPat_pt_nocut->Fill((trM+trP_P).Pt());
	htrIso_tag_allproPat_pt_nocut->Fill((trP+trM_P).Pt());
	htrIso_tag_allproPat_eta_nocut->Fill((trM+trP_P).Eta());
	htrIso_tag_allproPat_eta_nocut->Fill((trP+trM_P).Eta());
	if(fabs((trP+trM_P).M())>80 && fabs((trP+trM_P).M())<140)htrIso_tag_allproPat_pt_mass->Fill((trP+trM_P).Pt());
	if(fabs((trM+trP_P).M())>80 && fabs((trM+trP_P).M())<140)htrIso_tag_allproPat_pt_mass->Fill((trM+trP_P).Pt());
	if(fabs((trP+trM_P).M())>80 && fabs((trP+trM_P).M())<140)htrIso_tag_allproPat_eta_mass->Fill((trP+trM_P).Eta());
	if(fabs((trM+trP_P).M())>80 && fabs((trM+trP_P).M())<140)htrIso_tag_allproPat_eta_mass->Fill((trM+trP_P).Eta());

	//------------------------------------passed probe--------------------------------------------
	if(trP_P.DeltaR(PmuP)<0.01)htrIso_tag_proPat_pt_nocut->Fill((trM+trP_P).Pt());
	if(trM_P.DeltaR(PmuM)<0.01)htrIso_tag_proPat_pt_nocut->Fill((trP+trM_P).Pt());
	if(trP_P.DeltaR(PmuP)<0.01)htrIso_tag_proPat_eta_nocut->Fill((trM+trP_P).Eta());
	if(trM_P.DeltaR(PmuM)<0.01)htrIso_tag_proPat_eta_nocut->Fill((trP+trM_P).Eta());
	if(fabs((trP+trM_P).M())>80 && fabs((trP+trM_P).M())<140 && trM_P.DeltaR(PmuM)<0.01)htrIso_tag_proPat_pt_mass->Fill((trP+trM_P).Pt());
	if(fabs((trM+trP_P).M())>80 && fabs((trM+trP_P).M())<140 && trP_P.DeltaR(PmuP)<0.01)htrIso_tag_proPat_pt_mass->Fill((trM+trP_P).Pt());
	if(fabs((trP+trM_P).M())>80 && fabs((trP+trM_P).M())<140 && trM_P.DeltaR(PmuM)<0.01)htrIso_tag_proPat_eta_mass->Fill((trP+trM_P).Eta());
	if(fabs((trM+trP_P).M())>80 && fabs((trM+trP_P).M())<140 && trP_P.DeltaR(PmuP)<0.01)htrIso_tag_proPat_eta_mass->Fill((trM+trP_P).Eta());

	if(trP_P.DeltaR(PmuP)<0.01&&trP_P.Pt()>10)
	{
		htrIso_proPat_P->Fill(isopt_tr_pro_P);
		htr_pro_P->Fill(trP_P.Pt());//passed probe->Plus charge
		//cout<<(isopt_tr_pro_P)<<endl;
		htrIso_proPat_pt_P_nocut->Fill(trP_P.Pt());
		htrIso_proPat_eta_P_nocut->Fill(trP_P.Eta());
		htrIso_proPat_pt_nocut->Fill(trP_P.Pt());//---------------------total isolation probe no charge-----------
		htrIso_proPat_eta_nocut->Fill(trP_P.Eta());

		if(isopt_tr_pro_P<0.14)
		{
			htrIso_proPat_pt_P_iso->Fill(trP_P.Pt());
			htrIso_proPat_eta_P_iso->Fill(trP_P.Eta());
			htrIso_proPat_pt_iso->Fill(trP_P.Pt());//---------------------total isolation probe no charge-----------
			htrIso_proPat_eta_iso->Fill(trP_P.Eta());
		}
	}

	if(trM_P.DeltaR(PmuM)<0.01&&trM_P.Pt()>10)
	{
		htrIso_proPat_M->Fill(isopt_tr_pro_M);
		htr_pro_M->Fill(trM_P.Pt());//passed probe->Minus charge
		htrIso_proPat_pt_M_nocut->Fill(trM_P.Pt());
		htrIso_proPat_eta_M_nocut->Fill(trM_P.Eta());
		htrIso_proPat_pt_nocut->Fill(trP_P.Pt());//---------------------total isolation probe no charge-----------
		htrIso_proPat_eta_nocut->Fill(trP_P.Eta());

		if(isopt_tr_pro_M<0.14)
		{
			htrIso_proPat_pt_M_iso->Fill(trM_P.Pt());
			htrIso_proPat_eta_M_iso->Fill(trM_P.Eta());
			htrIso_proPat_pt_iso->Fill(trM_P.Pt());
			htrIso_proPat_eta_iso->Fill(trM_P.Eta());
		}
	}

	if(trP_P.Pt()>10)
	{
		htr_pro_P_tot->Fill(trP_P.Pt());
	}
	if(trM_P.Pt()>10)
	{
		htr_pro_M_tot->Fill(trM_P.Pt());
	}

	//-----------------------------Tag & Prove PatMuon with GlobalMuon------------------------------

	TLorentzVector PmuP_GT,PmuM_GT;//The 'tag' of PatMuon for GlobalMuon
	if(PnmuP!=-1&&PnmuM!=-1)
	{
		PmuP_GT.SetPxPyPzE((*Pmu)[PnmuP].px(),(*Pmu)[PnmuP].py(),(*Pmu)[PnmuP].pz(),(*Pmu)[PnmuP].p());
		PmuM_GT.SetPxPyPzE((*Pmu)[PnmuM].px(),(*Pmu)[PnmuM].py(),(*Pmu)[PnmuM].pz(),(*Pmu)[PnmuM].p());
	}

	double PmuptP_P=0,PmuptM_P=0;
	for(pat::MuonCollection::const_iterator pmu=Pmu->begin();pmu!=Pmu->end();pmu++)
	{
		if(fabs(pmu->pt())>10)
		{
			vpt.push_back(pmu->pt());
			double isopt=0,sumpt=0;
			if(pmu->charge()==1)
			{
				double iso_PmuptP=0;
				TLorentzVector pmu1_P, tr2_P;
				pmu1_P.SetPxPyPzE(pmu->px(),pmu->py(),pmu->pz(),pmu->p());
				for(reco::TrackCollection::const_iterator tr2=Tracks->begin();tr2!=Tracks->end();tr2++)
				{
					tr2_P.SetPxPyPzE(tr2->px(),tr2->py(),tr2->pz(),tr2->p());
					if(pmu1_P.DeltaR(tr2_P)<0.3 && pmu1_P.DeltaR(tr2_P)>0.01)
					{
						sumpt+=tr2->pt();
					}
				}
				isopt=sumpt/pmu->pt();
				if(iso_PmuptP<pmu->pt())
				{
					iso_PmuptP=pmu->pt();
					hpatIso_P->Fill(isopt);
				}
			}
			if(pmu->charge()==-1)
			{
				double iso_PmuptM=0;
				TLorentzVector pmu1_M, tr2_M;
				pmu1_M.SetPxPyPzE(pmu->px(),pmu->py(),pmu->pz(),pmu->p());
				for(reco::TrackCollection::const_iterator tr2=Tracks->begin();tr2!=Tracks->end();tr2++)
				{
					tr2_M.SetPxPyPzE(tr2->px(),tr2->py(),tr2->pz(),tr2->p());
					if(pmu1_M.DeltaR(tr2_M)<0.3 && pmu1_M.DeltaR(tr2_M)>0.01)sumpt+=tr2->pt();
				}
				isopt=sumpt/pmu->pt();
				if(iso_PmuptM<pmu->pt())
				{
					iso_PmuptM=pmu->pt();
					hpatIso_M->Fill(isopt);
				}
			}
		}

		if(fabs(pmu->pt())>10&&pmu->isGlobalMuon())
		{
			TLorentzVector tag_PmuP, tag_PmuM;//tag candidate
			vpt.push_back(pmu->pt());
			if(pmu->charge()==1)
			{
				tag_PmuP.SetPxPyPzE(pmu->px(),pmu->py(),pmu->pz(),pmu->p());
				if(PmuptP_P<pmu->pt())// && minDRP>dR && dR<0.01)
				{
					PmuptP_P=pmu->pt();//selected highest pt muon plus charge
					PmuP_GT=tag_PmuP;//save "tag" highest pt muon plus charge at LorentzVector
					//cout<<(PmuP_GT.Pt())<<endl;
				}
			}
			if(pmu->charge()==-1)
			{
				tag_PmuM.SetPxPyPzE(pmu->px(),pmu->py(),pmu->pz(),pmu->p());
				if(PmuptM_P<pmu->pt())// && minDRM>dR && dR<0.01)
				{
					PmuptM_P=pmu->pt();//selected highest pt muon minus charge
					PmuM_GT=tag_PmuM;//save "tag" highest pt muon minus charge at LorentzVector
					//cout<<(PmuM_GT.Pt())<<endl;
				}
			}
		}
	}

	TLorentzVector PmuP_GP,PmuM_GP;//for save "probe" of PatMuon LorentzVector
	if(PnmuP!=-1&&PnmuM!=-1)
	{
		PmuP_GP.SetPxPyPzE((*Pmu)[PnmuP].px(),(*Pmu)[PnmuP].py(),(*Pmu)[PnmuP].pz(),(*Pmu)[PnmuP].p());
		PmuM_GP.SetPxPyPzE((*Pmu)[PnmuM].px(),(*Pmu)[PnmuM].py(),(*Pmu)[PnmuM].pz(),(*Pmu)[PnmuM].p());
	}

	//-----------------------------Prove
	TLorentzVector probe_PmuP, probe_PmuM;
	for(pat::MuonCollection::const_iterator pmu=Pmu->begin(); pmu!=Pmu->end();pmu++)
	{//Normal PatMuon-->'pat'
		if(fabs(pmu->pt())>10)
		{
			vpt.push_back(pmu->pt());
			if(pmu->charge()==-1)
			{
				probe_PmuM.SetPxPyPzE(pmu->px(),pmu->py(),pmu->pz(),pmu->p());
				double dR=PmuP_GT.DeltaR(probe_PmuM);//probe candidate
				if(PmuptM_P<pmu->pt()&&dR>1)
				{
					minDRM=dR;
					PmuptM_P=pmu->pt();
					PmuM_GP=probe_PmuM;//save "probe" highest pt muon plus charge and candidate probe
					//if(pmu->isGlobalMuon())pro_muM=true;//selected "probe" is GlobalMuon
					if(pmu->isGlobalMuon())pro_muM=1;
					//else pro_muM=false;
					else pro_muM=0;
				}
			}
			if(pmu->charge()==1)
			{
				probe_PmuP.SetPxPyPzE(pmu->px(),pmu->py(),pmu->pz(),pmu->p());
				double dR=PmuM_GT.DeltaR(probe_PmuP);//probe candidate
				if(PmuptP_P<pmu->pt()&&dR>1)
				{
					minDRP=dR;
					PmuptP_P=pmu->pt();
					PmuP_GP=probe_PmuP;//save "probe" highest pt muon minus charge and candidate probe
					//					cout<<(PmuP_GP.Pt())<<endl;
					if(pmu->isGlobalMuon())pro_muP=1;//selectec "probe" is GlobalMuon
					else pro_muP=0;
				}
			}
		}
	}
	//--------------------------------passed probe section---------------------
	if(PmuP_GP.Pt()>10)
	{
		hpmu_pro_P_tot->Fill(PmuP_GP.Pt());//total probe +charge
		hpmu_pro_tot->Fill(PmuP_GP.Pt());//no charge total probe
	}
	if(PmuM_GP.Pt()>10)
	{
		hpmu_pro_M_tot->Fill(PmuM_GP.Pt());//total probe -charge
		hpmu_pro_tot->Fill(PmuM_GP.Pt());
	}
	//---------------------------------------total probe section--------------------
	if(PmuP_GP.Pt()>10)
	{
		if(pro_muP==1)hpmu_pro_P->Fill(PmuP_GP.Pt());//passed probe +charge
		if(pro_muP==1)hpmu_pro->Fill(PmuP_GP.Pt());//no charge passed probe
	}
	if(PmuM_GP.Pt()>10)
	{
		if(pro_muM==1)hpmu_pro_M->Fill(PmuM_GP.Pt());//passed probe -charge
		if(pro_muM==1)hpmu_pro->Fill(PmuM_GP.Pt());
	}

	if(tempLepsP.size()==2)
	{
		for(reco::CandidateCollection::const_iterator l1=tempLepsP.begin();l1!=tempLepsP.end()-1;l1++)
			for(reco::CandidateCollection::const_iterator l2=l1+1;l2!=tempLepsP.end(); l2++)
			{
				hmllpm->Fill((l1->p4()+l2->p4()).M());
				mllpm=(l1->p4()+l2->p4()).M();
			}
		//mllpm=(tempLepsP[0].p4()+tempLepsP[1].p4()).M();
	}

	//---------globalmuon section

	Handle <reco::TrackCollection> Gmu;
	iEvent.getByLabel("globalMuons",Gmu);
	//--------------------------------------------isopt of global muon--------------------
	unsigned int ngmu=0;
	int GnmuP=-1,GnmuM=-1;
	double GmuptP=0,GmuptM=0;
	for(reco::TrackCollection::const_iterator gmu=Gmu->begin();gmu!=Gmu->end();gmu++,ngmu++)
	{
		if(fabs(gmu->eta())<2.4&&fabs(gmu->pt())>10)
		{
			double isopt=0,sumpt=0;
			vpt.push_back(gmu->pt());
			if(gmu->charge()==1)
			{
				double iso_gmuptP=0;
				TLorentzVector gmu1_P, tr2_P;
				gmu1_P.SetPxPyPzE(gmu->px(),gmu->py(),gmu->pz(),gmu->p());
				for(reco::TrackCollection::const_iterator tr2=Tracks->begin();tr2!=Tracks->end();tr2++)
				{
					tr2_P.SetPxPyPzE(tr2->px(),tr2->py(),tr2->pz(),tr2->p());
					if(gmu1_P.DeltaR(tr2_P)<0.3 && gmu1_P.DeltaR(tr2_P)>0.01)
					{
						sumpt+=tr2->pt();
						hgmuIso_sum_P->Fill(sumpt);
					}
				}
				isopt=sumpt/gmu->pt();
				if(iso_gmuptP<gmu->pt())
				{
					iso_gmuptP=gmu->pt();
					hgmuIso_P->Fill(isopt);
					if(isopt<0.14)
					{
						hgmuIso_pt_P->Fill(gmu->pt());
						hgmuIso_eta_P->Fill(gmu->eta());
					}
				}
				if(GmuptP<gmu->pt())
				{
					GnmuP=ngmu;
					GmuptP=gmu->pt();
				}
			}
			if(gmu->charge()==-1)
			{
				double iso_gmuptM=0;
				TLorentzVector gmu1_M,tr2_M;
				gmu1_M.SetPxPyPzE(gmu->px(),gmu->py(),gmu->pz(),gmu->p());
				for(reco::TrackCollection::const_iterator tr2=Tracks->begin();tr2!=Tracks->end();tr2++)
				{
					tr2_M.SetPxPyPzE(tr2->px(),tr2->py(),tr2->pz(),tr2->p());
					if(gmu1_M.DeltaR(tr2_M)<0.3 && gmu1_M.DeltaR(tr2_M)>0.01)
					{
						sumpt+=tr2->pt();
						hgmuIso_sum_M->Fill(sumpt);
					}
				}
				isopt=sumpt/gmu->pt();
				if(iso_gmuptM<gmu->pt())
				{
					iso_gmuptM=gmu->pt();
					hgmuIso_M->Fill(isopt);

					if(isopt<0.14)
					{
						hgmuIso_pt_M->Fill(gmu->pt());
						hgmuIso_eta_M->Fill(gmu->eta());
					}
				}
				if(GmuptM<gmu->pt())
				{
					GnmuM=ngmu;
					GmuptM=gmu->pt();
				}
			}
			hgmu_pt->Fill(gmu->pt());
			hgmu_eta->Fill(gmu->eta());
		}
		Tgmu_pt[ngmu] = gmu->pt();
		Tgmu_eta[ngmu] = gmu->eta();
	}
	ng=ngmu;

	unsigned int tot_gnmu=0;
	if(GmuptP>10&&GmuptM>10)
	{
		for(reco::TrackCollection::const_iterator tot_gmu=Gmu->begin();tot_gmu!=Gmu->end();tot_gmu++,tot_gnmu++)
		{
			if(fabs(tot_gmu->pt())>10)
			{
				vpt.push_back(tot_gmu->pt());
				if(tot_gmu->charge()==1)
				{
					GmuptP=tot_gmu->pt();
					GnmuP=tot_gnmu;
				}
				if(tot_gmu->charge()==-1)
				{
					GmuptM=tot_gmu->pt();
					GnmuP=tot_gnmu;
				}	
				htot_gmu_pt->Fill(tot_gmu->pt());
			}
		}
	}

	TLorentzVector GmuP, GmuM;//LorentzVector of GlobalMuon in TrackCollection
	if(GnmuP!=-1&&GnmuM!=-1)
	{
		GmuP.SetPxPyPzE((*Gmu)[GnmuP].px(),(*Gmu)[GnmuP].py(),(*Gmu)[GnmuP].pz(),(*Gmu)[GnmuP].p());
		GmuM.SetPxPyPzE((*Gmu)[GnmuM].px(),(*Gmu)[GnmuM].py(),(*Gmu)[GnmuM].pz(),(*Gmu)[GnmuM].p());
	}

	TLorentzVector trP_G, trM_G;//LorentzVector 'tag' with GlobalMuon
	if(ntrP!=-1&&ntrM!=-1)
	{
		trP_G.SetPxPyPzE((*Tracks)[ntrP].px(),(*Tracks)[ntrP].py(),(*Tracks)[ntrP].pz(),(*Tracks)[ntrP].p());
		trM_G.SetPxPyPzE((*Tracks)[ntrM].px(),(*Tracks)[ntrM].py(),(*Tracks)[ntrM].pz(),(*Tracks)[ntrM].p());
	}

	//------------------------------------------TrackCollection of tag&probe with GlobalMuon

	TLorentzVector tag_trP_G,tag_trM_G;//
	double trptP_G=0,trptM_G=0;
	for(reco::TrackCollection::const_iterator tr=Tracks->begin();tr!=Tracks->end();tr++)
	{
		if(fabs(tr->pt())>10)
		{
			vpt.push_back(tr->pt());
			if(tr->charge()==1)
			{
				tag_trP_G.SetPxPyPzE(tr->px(),tr->py(),tr->pz(),tr->p());
				double dR=tag_trP_G.DeltaR(GmuP);
				if(trptP_G<tr->pt()&&dR<0.01)//'tag+' is completely muon
				{
					trptP_G=tr->pt();
					minDRP=dR;
					trP_G=tag_trP_G;//save a tr+ the 'tag_trP_G' in trP_G
				}
			}
			if(tr->charge()==-1)
			{
				tag_trM_G.SetPxPyPzE(tr->px(),tr->py(),tr->pz(),tr->p());
				double dR=tag_trM_G.DeltaR(GmuM);
				if(trptM_G<tr->pt()&&dR<0.01)//tag- is completely muon
				{
					trptM_G=tr->pt();
					minDRM=dR;
					trM_G=tag_trM_G;//save a tr+ the 'tag_trM_G' in trM_G
				}
			}
		}
	}

	TLorentzVector trP_GP, trM_GP;//for "probe" LorentzVector with GlobaMuon
	if(ntrP!=-1&&ntrM!=-1)
	{
		trP_GP.SetPxPyPzE((*Tracks)[ntrP].px(),(*Tracks)[ntrP].py(),(*Tracks)[ntrP].pz(),(*Tracks)[ntrP].p());
		trM_GP.SetPxPyPzE((*Tracks)[ntrM].px(),(*Tracks)[ntrM].py(),(*Tracks)[ntrM].pz(),(*Tracks)[ntrM].p());
	}

	TLorentzVector probe_trP_G, probe_trM_G;
	for(reco::TrackCollection::const_iterator tr=Tracks->begin();tr!=Tracks->end();tr++)
	{
		if(fabs(tr->pt())>10)
		{
			vpt.push_back(tr->pt());
			if(tr->charge()==-1)
			{
				probe_trM_G.SetPxPyPzE(tr->px(),tr->py(),tr->pz(),tr->p());
				double dR = trP_G.DeltaR(probe_trM_G);//probe candidate minus charge
				if(trptM_G<tr->pt()&&dR>1)
				{
					minDRM=dR;
					trptM_G=tr->pt();
					trM_GP=probe_trM_G;//save 'probe'Lorentzvector
				}
			}
			if(tr->charge()==1)
			{
				probe_trP_G.SetPxPyPzE(tr->px(),tr->py(),tr->pz(),tr->p());
				double dR = trM_G.DeltaR(probe_trP_G);//probe candidate plus charge
				if(trptP_G<tr->pt()&&dR>1)
				{
					minDRP=dR;
					trptP_G=tr->pt();
					trP_GP=probe_trP_G;//save 'probe'LorentzVector
				}
			}
		}
	}
	//---------------------------------passed probe section----------------
	if(trP_GP.DeltaR(GmuP)<0.01&&trP_GP.Pt()>10)//passed probe->Plus charge
	{
		hgmu_pro_P->Fill(trP_GP.Pt());
		hgmu_pro->Fill(trP_GP.Pt());//passed probe plus and minus charge__no charge
		//cout<<(trP_GP.Pt())<<endl;
	}
	if(trM_GP.DeltaR(GmuM)<0.01&&trM_GP.Pt()>10)//passed probe->Minus charge
	{
		hgmu_pro_M->Fill(trM_GP.Pt());
		hgmu_pro->Fill(trM_GP.Pt());
		//cout<<(trP_GP.Pt())<<endl;
	}
	//---------------------------------total probe section-----------------
	if(trP_GP.Pt()>10)
	{
		hgmu_pro_P_tot->Fill(trP_GP.Pt());
		hgmu_pro_tot->Fill(trP_GP.Pt());//total probe plus and minus charge_no charge
		//cout<<(trP_GP.Pt())<<endl;//total probe plus charge
	}
	if(trM_GP.Pt()>10)
	{
		hgmu_pro_M_tot->Fill(trM_GP.Pt());
		hgmu_pro_tot->Fill(trM_GP.Pt());
		//cout<<(trM_GP.Pt())<<endl;//total probe minus charge
	}

	//-----------MET section
	//	Handle<pat::METCollection> pMET;
	//	iEvent.getByLabel("patMETs",pMET);
	//	unsigned int nmet=0;
	//	for(pat::METCollection::const_iterator pMet = pMET->begin();pMet!= pMET->end();pMet++,nmet++)
	//	{
	//		hmet_pt->Fill(pMet->pt());
	//		Tmet_pt[nmet] = pMet->pt();
	//	}
	//	met=nmet;

	tree->Fill();
}
//-----------------------------------------------------

// ------------ method called once each job just before starting event loop  ------------
	void 
DemoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
DemoAnalyzer::endJob() 
{
	htr_div_ProP->Divide(htr_pro_P,htr_pro_P_tot);
	htr_div_ProM->Divide(htr_pro_M,htr_pro_M_tot);
	hpmu_div_ProP->Divide(hpmu_pro_P,hpmu_pro_P_tot);
	hpmu_div_ProM->Divide(hpmu_pro_M,hpmu_pro_M_tot);

	hpmu_pro->Sumw2();
	hpmu_pro_tot->Sumw2();
	hpmu_div_Pro->Divide(hpmu_pro,hpmu_pro_tot,1,1,"b");//Error is larger as value is far from 1, because number of events is not many.

	hgpmu_div_pt->Divide(hpat_g_pt,hpmu_pt);
	hgpmu_div_eta->Divide(hpat_g_eta,hpmu_eta);

	hgmu_div_ProP->Divide(hgmu_pro_P,hgmu_pro_P_tot);
	hgmu_div_ProM->Divide(hgmu_pro_M,hgmu_pro_M_tot);

	hgmu_pro->Sumw2();
	hgmu_pro_tot->Sumw2();
	hgmu_div_Pro->Divide(hgmu_pro,hgmu_pro_tot,1,1,"b");

	htrIso_tag_proPat_mllpm->Sumw2();
	htrIso_tag_allproPat_mllpm->Sumw2();
	htrIso_div_mllpm->Divide(htrIso_tag_proPat_mllpm,htrIso_tag_allproPat_mllpm,1,1,"b");


	htrIso_tag_proPat_pt_nocut->Sumw2();
	htrIso_tag_proPat_pt_iso->Sumw2();
	htrIso_tag_proPat_pt_mi->Sumw2();
	htrIso_tag_proPat_pt_mass->Sumw2();

	htrIso_tag_allproPat_pt_nocut->Sumw2();
	htrIso_tag_allproPat_pt_iso->Sumw2();
	htrIso_tag_allproPat_pt_mi->Sumw2();
	htrIso_tag_allproPat_pt_mass->Sumw2();

	htrIso_div_tag_pt_nocut->Divide(htrIso_tag_proPat_pt_nocut,htrIso_tag_allproPat_pt_nocut,1,1,"b");
	htrIso_div_tag_pt_iso->Divide(htrIso_tag_proPat_pt_iso,htrIso_tag_allproPat_pt_iso,1,1,"b");
	htrIso_div_tag_pt_mi->Divide(htrIso_tag_proPat_pt_mi,htrIso_tag_allproPat_pt_mi,1,1,"b");
	htrIso_div_tag_pt_mass->Divide(htrIso_tag_proPat_pt_mass,htrIso_tag_allproPat_pt_mass,1,1,"b");

	htrIso_tag_proPat_eta_nocut->Sumw2();
	htrIso_tag_proPat_eta_iso->Sumw2();
	htrIso_tag_proPat_eta_mi->Sumw2();
	htrIso_tag_proPat_eta_mass->Sumw2();

	htrIso_tag_allproPat_eta_nocut->Sumw2();
	htrIso_tag_allproPat_eta_iso->Sumw2();
	htrIso_tag_allproPat_eta_mi->Sumw2();
	htrIso_tag_allproPat_eta_mass->Sumw2();

	htrIso_div_tag_eta_nocut->Divide(htrIso_tag_proPat_eta_nocut,htrIso_tag_allproPat_eta_nocut,1,1,"b");
	htrIso_div_tag_eta_iso->Divide(htrIso_tag_proPat_eta_iso,htrIso_tag_allproPat_eta_iso,1,1,"b");
	htrIso_div_tag_eta_mi->Divide(htrIso_tag_proPat_eta_mi,htrIso_tag_allproPat_eta_mi,1,1,"b");
	htrIso_div_tag_eta_mass->Divide(htrIso_tag_proPat_eta_mass,htrIso_tag_allproPat_eta_mass,1,1,"b");

	//htrIso_tag_failproPat_mllpm->Sumw2();

	if(SaveHisto)
	{
		demo_tree->cd();
		demohisto->Write();
		tree->Write();

		hgen_pt->Write();
		hgen_eta->Write();

		hpmu_pt->Write();
		hpmu_eta->Write();
		htot_pmu_pt->Write();

		htr_pt->Write();
		htr_eta->Write();
		htot_tr_pt->Write();

		hgmu_pt->Write();
		hgmu_eta->Write();
		htot_gmu_pt->Write();

		hDeltaT->Write();
		hmllpmT->Write();
		hDeltaP->Write();

		hmllpmP->Write();
		hmllpm->Write();
		hmllpmT_D->Write();
		hmllpmP_D->Write();

		htr_pro_P->Write();htr_pro_M->Write();
		htr_pro_P_tot->Write();htr_pro_M_tot->Write();

		hpmu_pro_P->Write();hpmu_pro_M->Write();
		hpmu_pro_P_tot->Write();hpmu_pro_M_tot->Write();

		hgmu_pro_P->Write();hgmu_pro_M->Write();
		hgmu_pro_P_tot->Write();hgmu_pro_M_tot->Write();

		hpmu_pro->Write();//nocharge probe muon
		hpmu_pro_tot->Write();
		hpmu_div_Pro->Write();// efficiency passed probe of candidate probe.

		hgmu_pro->Write();
		hgmu_pro_tot->Write();
		hgmu_div_Pro->Write();

		htr_div_ProP->Write();htr_div_ProM->Write();
		hpmu_div_ProP->Write();hpmu_div_ProM->Write();
		hgmu_div_ProP->Write();hgmu_div_ProM->Write();

		hpat_g_pt->Write();
		hpat_g_eta->Write();
		hgpmu_div_pt->Write();
		hgpmu_div_eta->Write();
		//--------------------Track Isolation----------------------
		htrIso_P->Write();htrIso_M->Write();
		htrIso_tagPat_P->Write();htrIso_tagPat_M->Write();
		htrIso_allproPat_P->Write();htrIso_allproPat_M->Write();
		htrIso_proPat_P->Write();htrIso_proPat_M->Write();
		//--------------Track Isolation Pt<0.14---------------------
		htrIso_proPat_pt_P_iso->Write();htrIso_proPat_pt_M_iso->Write();
		htrIso_proPat_eta_P_iso->Write();htrIso_proPat_eta_M_iso->Write();
		htrIso_proPat_pt_iso->Write();
		htrIso_proPat_eta_iso->Write();

		//------------------Track nocut Isolation------------------------
		htrIso_proPat_pt_P_nocut->Write();htrIso_proPat_pt_M_nocut->Write();
		htrIso_proPat_eta_P_nocut->Write();htrIso_proPat_eta_M_nocut->Write();
		htrIso_proPat_pt_nocut->Write();
		htrIso_proPat_eta_nocut->Write();

		htrIso_allproPat_pt_nocut->Write();
		htrIso_allproPat_eta_nocut->Write();

		htrIso_tag_allproPat_mllpm->Write();
		htrIso_tag_proPat_mllpm->Write();
		htrIso_tag_failproPat_mllpm->Write();

		htrIso_tag_allproPat_pt_iso->Write();htrIso_tag_allproPat_eta_iso->Write();
		htrIso_tag_allproPat_pt_mi->Write();htrIso_tag_allproPat_eta_mi->Write();
		htrIso_tag_allproPat_pt_nocut->Write();htrIso_tag_allproPat_eta_nocut->Write();
		htrIso_tag_allproPat_pt_mass->Write();htrIso_tag_allproPat_eta_mass->Write();

		htrIso_tag_proPat_pt_iso->Write();htrIso_tag_proPat_eta_iso->Write();
		htrIso_tag_proPat_pt_mi->Write();htrIso_tag_proPat_eta_mi->Write();
		htrIso_tag_proPat_pt_nocut->Write();htrIso_tag_proPat_eta_nocut->Write();
		htrIso_tag_proPat_pt_mass->Write();htrIso_tag_proPat_eta_mass->Write();

		htrIso_div_mllpm->Write();
		htrIso_div_tag_pt_nocut->Write();htrIso_div_tag_pt_iso->Write();
		htrIso_div_tag_pt_mass->Write();htrIso_div_tag_pt_mi->Write();
		htrIso_div_tag_eta_nocut->Write();htrIso_div_tag_eta_iso->Write();
		htrIso_div_tag_eta_mass->Write();htrIso_div_tag_eta_mi->Write();

		//--------------------PatMuon Isolation-----------------------
		hpatIso_P->Write();hpatIso_M->Write();
		//--------------------GlobalMuon Isolation----------------------
		hgmuIso_P->Write();hgmuIso_M->Write();
		hgmuIso_pt_P->Write();hgmuIso_eta_P->Write();
		hgmuIso_pt_M->Write();hgmuIso_eta_M->Write();
		hgmuIso_sum_P->Write();hgmuIso_sum_M->Write();
		demo_tree->Close();
	}
}

// ------------ method called when starting to processes a run  ------------
/*
   void 
   DemoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a run  ------------
/*
   void 
   DemoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when starting to processes a luminosity block  ------------
/*
   void 
   DemoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method called when ending the processing of a luminosity block  ------------
/*
   void 
   DemoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
   {
   }
   */

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
DemoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(DemoAnalyzer);
