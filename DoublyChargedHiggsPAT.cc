// -*- C++ -*-
//
// Package:    DoublyChargedHiggsPAT
// Class:      DoublyChargedHiggsPAT
// 
/**\class DoublyChargedHiggsPAT DoublyChargedHiggsPAT.cc Demo/DoublyChargedHiggsPAT/src/DoublyChargedHiggsPAT.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Jongseok Lee
//         Created:  Mon Oct  5 13:56:37 CEST 2009
// $Id$
//
//


// system include files
#include <memory>
#include <map>
#include <vector>
#include <iostream>
#include <fstream>
 
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
 
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
//#include "FWCore/ParameterSet/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
 
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonIsolation.h"
#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"
#include "DataFormats/Common/interface/Ref.h"
//#include "DataFormats/HepMCCandidate/interface/GenParticleCandidate.h"
//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Candidate/interface/CompositeRefCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

/*
//Jet
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/JetReco/interface/JetCollection.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenericJetCollection.h"
#include "DataFormats/JetReco/interface/GenericJet.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
 
//Jet Corrections
#include "JetMETCorrections/Objects/interface/JetCorrector.h"
#include "JetMETCorrections/Objects/interface/JetCorrectionsRecord.h"
#include "JetMETCorrections/Objects/interface/ChainedJetCorrector.h"

//MET
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/METCollection.h"
#include "DataFormats/METReco/interface/CorrMETData.h"
 
//Vertex
#include "RecoVertex/VertexPrimitives/interface/VertexReconstructor.h"
#include "RecoVertex/VertexPrimitives/interface/VertexFitter.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/AdaptiveVertexReconstructor.h"
#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableAdaptiveReconstructor.h"
#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableAdaptiveFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/TransientTrack/interface/TrackTransientTrack.h"
#include "DataFormats/Common/interface/Handle.h" 
#include "TrackingTools/TransientTrack/interface/GsfTransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackFromFTS.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
*/ 
#include "DataFormats/VertexReco/interface/Vertex.h"
//#include "DataFormats/MuonReco/interface/MuonIsolation.h"
//#include "DataFormats/MuonReco/interface/MuIsoDeposit.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/JetCorrFactors.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Conversion.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

//#include "DataFormats/PatCandidates/interface/Isolation.h"
//#include "DataFormats/PatCandidates/interface/Lepton.h"

#include "DataFormats/PatCandidates/interface/Lepton.h"
//#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
//#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
 
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/HLTResult.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include "HLTrigger/HLTanalyzers/interface/HLTInfo.h"
#include "HLTrigger/HLTfilters/interface/HLTHighLevel.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3D.h>
#include <TMath.h>
#include <TLorentzVector.h>
#include <TTree.h>

#define kElStart	1
#define kElBasic	kElStart+1
#define kElLowInv	kElBasic+1
#define kElZpole	kElLowInv+1
#define kElAllInv	kElZpole+1
#define kElMCTrue	kElAllInv+1
#define kMuStart	kElMCTrue+1
#define kMuBasic	kMuStart+1
#define kMuLowInv	kMuBasic+1
#define kMuZpole	kMuLowInv+1
#define kMuAllInv	kMuZpole+1
#define kMuMCTrue	kMuAllInv+1
#define kTauStart	kMuMCTrue+1
#define kTauBasic	kTauStart+1
#define kTauMCTrue	kTauBasic+1
//
// class decleration
//
using namespace edm;
using namespace std;
using namespace reco;
//using namespace pat;

int events=0;
int ntest=0;
bool intau=true , debug0=false, debug=false, debug2=false;
TString NAME_ch[22] = {"MC bkg","eeee","eeem","eeet","eemm","eemt","eett","emem","emet","emmm","emmt","emtt","etet","etmm","etmt","ettt","mmmm","mmmt","mmtt","mtmt","mttt","tttt"};
int run,lumi,ev,nlep,nlep_acc,nlgen, mP,mM,eP,eM,tP,tM, Mp,Mm,Ep,Em,Tp,Tm, MP,MM,EP,EM,TP,TM; // lp-# of MC gen leptons, LP-# of MC gen leptons in acceptance
int Mp1,Mm1,Ep1,Em1,Tp1,Tm1;
int nlp, nlm, pid1, pid2, pid3, pid4;
float firstpt, secondpt, pt1, pt2, pt3, pt4, iso1, iso2, iso3, iso4, mp1, mn1, mp2, mn2, mp, mn, mp4, mn4, dz, dZ, dilpt, m4, mpp, mpm, mnp, mnm, m4p, m4m;
float eta1, eta2, eta3, eta4;
int nvtx, nvtxg, npum1, npu, npup1, Tnpu;
const int nMAX=100;
float vtxndof[nMAX], vtxnchi2[nMAX], vtxz[nMAX], vtxrho[nMAX], vtxleprho[nMAX], vtxmcrho[nMAX];
float vtxndof_, vtxnchi2_, vtxz_, vtxrho_, vtxleprho_, vtxmcrho_;
float vtxndofmc_, vtxnchi2mc_, vtxzmc_, vtxrhomc_, vtxleprhomc_;
float mu_pt[nMAX], mu_eta[nMAX], mu_dz[nMAX], mu_dxy[nMAX], mu_rho[nMAX], mu_dB[nMAX], mu_sip[nMAX], mu_nchi2[nMAX];
float mu_trackIso[nMAX], mu_ecalIso[nMAX], mu_hcalIso[nMAX], mu_totIso[nMAX], mu_totIso2[nMAX];
float mu_dimass1[nMAX], mu_dimass2[nMAX], mu_dimass3[nMAX], mu_dimass4[nMAX];
int nm, mu_glb[nMAX], mu_trk[nMAX], mu_pf[nMAX], mu_good[nMAX], mu_tight[nMAX], mu_hits_val[nMAX], mu_hits_trk[nMAX], mu_hits_pix[nMAX], mu_hits_lay[nMAX], mu_nstation[nMAX], mu_genmat[nMAX];
float el_pt[nMAX], el_eta[nMAX], el_dz[nMAX], el_dxy[nMAX], el_rho[nMAX], el_dB[nMAX], el_sip[nMAX], el_nchi2[nMAX];
float el_trackIso[nMAX], el_ecalIso[nMAX], el_hcalIso[nMAX], el_totIso[nMAX], el_totIso2[nMAX];
float el_dimass1[nMAX], el_dimass2[nMAX], el_dimass3[nMAX], el_dimass4[nMAX];
int  ne, el_chargeok[nMAX], el_losthit[nMAX], el_eidVeryLoose[nMAX], el_eidLoose[nMAX], el_eidMedium[nMAX], el_eidTight[nMAX], el_eidSuperTight[nMAX], el_eidHyperTight1[nMAX], el_genmat[nMAX];
float el_mvaTrigV0[nMAX];
int el_passConvVeto[nMAX], el_misshits[nMAX];
int nt, tau_signalTracks[nMAX], tau_charge[nMAX];
int tau_againstElectronLoose[nMAX], tau_againstElectronMedium[nMAX], tau_againstElectronTight[nMAX], tau_againstElectronMVA[nMAX];
int tau_againstMuonLoose[nMAX], tau_againstMuonMedium[nMAX], tau_againstMuonTight[nMAX];
int tau_byLooseCombinedIsolationDeltaBetaCorr[nMAX], tau_byMediumCombinedIsolationDeltaBetaCorr[nMAX], tau_byTightCombinedIsolationDeltaBetaCorr[nMAX];
int tau_decayModeFinding[nMAX];
int tau_overlap[nMAX], tau_genmat[nMAX];
float tau_pt[nMAX], tau_eta[nMAX];
int HLT;
float lep_eff_pp, lep_eff_pp_high, lep_eff_pp_low;
float lep_eff_mm, lep_eff_mm_high, lep_eff_mm_low;
float METpt;
float mllpm, mllpmp, mllpmm, mllss, mZ1, mZ2, mW;
// mllpm = mZ, mllpmp = mZ+pt_uncertianty, mllpmm = mZ-pt_uncertianty
		
TString outroot;
TString emt[3]={"e","m","t"};


typedef std::pair<float, int> ptvsindex;

double deltaphi(double ph1, double ph2)
{
	// in ORCA phi = [0,2pi], in TLorentzVector phi = [-pi,pi].
	// With the conversion below deltaPhi works ok despite the
	// 2*pi difference in phi definitions.
	double PI = 3.1415;
	if(ph1 < 0) ph1 += 2*PI;
	if(ph2 < 0) ph2 += 2*PI;
	
	double dphi = fabs(ph1-ph2);
	
	if(dphi > PI) dphi = 2*PI - dphi;
	return dphi;
}
/*
double calcIso(const reco::Candidate &cand) {
        double iso = -1;
        if (abs(cand.pdgId()) == 11) {
                const pat::Electron *el = dynamic_cast<const pat::Electron*>(&cand);
//                iso=(el->hcalIso()+el->ecalIso()+el->trackIso())/el->pt();
                iso=(el->chargedHadronIso()+el->neutralHadronIso()+el->photonIso())/el->pt();
//                iso=el->trackIso()/el->pt();
        }
        if (abs(cand.pdgId()) == 13) {
                const pat::Muon *mu = dynamic_cast<const pat::Muon*>(&cand);
//                iso=mu->isolationR03().sumPt/mu->pt();
//                iso=(mu->hcalIso()+mu->ecalIso()+mu->trackIso())/mu->pt();
                iso=(mu->chargedHadronIso()+mu->neutralHadronIso()+mu->photonIso())/mu->pt();
        }
        if (!cand.pdgId() || abs(cand.pdgId()) == 15)
                iso=0;
        return iso;
}               
*/
class DoublyChargedHiggsPAT : public edm::EDAnalyzer {
   public:
      explicit DoublyChargedHiggsPAT(const edm::ParameterSet&);
      ~DoublyChargedHiggsPAT();


   private:
//      virtual void beginJob(const edm::EventSetup&) ;
      virtual void beginJob();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;

//      virtual void beginRun(edm::Run&, edm::EventSetup const&);
//      virtual void endRun(edm::Run&, edm::EventSetup const&);
      virtual void beginRun(const edm::Run&, const edm::EventSetup&);
      virtual void endRun(const edm::Run&, const edm::EventSetup&);

      // ----------member data ---------------------------
	double matmcl(TLorentzVector t, TLorentzVector *mcmu, int nmcl);
	double calcIso_reco(const reco::Candidate &cand, bool, bool);
	double calcIso_el(const pat::Electron &cand, bool);
	double calcIso_mu(const pat::Muon &cand, bool);
//	double calcIso_pat(const pat::Candidate &cand, bool);
	double eaIso(const pat::Electron  &el, double rho);
	double Zmass(double pt1, double pt2, double phi1, double phi2, double pz1, double pz2, double e1, double e2);
//	double lepton_eff(const reco::Candidate &cand);
//	double lepton_eff(int lep_pid, double lep_pt, double lep_eta);
	double lepton_eff(int lep_pid, double lep_pt, double lep_eta, int lepeff_uncertainty);

	double rho;
        int gench, genmatch1, genmatch2, genmatch, genmatch4;
	int recoch1, recoch;
        int channel[3][3][3][3][3][3];// = {0,}; //channel[eP][mP][tP][eM][mM][tM], lP = # of plus leptons, lM = # of minus leptons
	float dR_matching;

        vector<string> HLT_path;
        int HLT_path_find[100];
        bool isValidHltConfig_;
        HLTConfigProvider  hltConfigProvider_;

        string output;
	TTree *tree;
	TFile *HPP_TTree;
	bool dogen;
	edm::InputTag elLabel, muLabel, tauLabel, vertLabel, metLabel;
	std::string elId;
	bool debug, lowMass, zPole, docoll;//, zVeto;
	double muMinPt, elMinPt, tauMinPt, lowInvCut, zInvCut;
	int mc, id, run, lumi, event, elIdComp;
	double eId,px,py,pz,E,iso,pvndof,pvmaxz,pvmaxd0;//,dZ;
	std::map<std::string,TH1D*> hc_;
	std::map<std::string,TH2D*> hc2_;

	TH1D *hmupt, *hmueta, *helpt, *heleta;

	TH1D *hprobe_best_pt, *hprobe_best_eta, *hprobe_best_m;
	TH1D *hprobe_match_pt, *hprobe_match_eta, *hprobe_match_m;
	TH1D *hprobe_matmc_pt, *hprobe_matmc_eta, *hprobe_matmc_m;

	double DeltaPt; // pt uncertianty
//	TFile *f_3DPU;
//	TH3D *h_3DPU;
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
DoublyChargedHiggsPAT::DoublyChargedHiggsPAT(const edm::ParameterSet& iConfig)
:output(iConfig.getParameter<std::string>("output"))
{
	elLabel 	= iConfig.getParameter<edm::InputTag>("elLabel");
	muLabel 	= iConfig.getParameter<edm::InputTag>("muLabel");
	tauLabel	= iConfig.getParameter<edm::InputTag>("tauLabel");
	elId		= iConfig.getParameter<std::string>("elId");
	vertLabel	= iConfig.getParameter<edm::InputTag>("vertLabel");
	
	debug		= iConfig.getParameter<bool>("debug");
	pvndof		= iConfig.getParameter<double>("pvndof");
	pvmaxz		= iConfig.getParameter<double>("pvmaxz");
	pvmaxd0		= iConfig.getParameter<double>("pvmaxd0");
	muMinPt		= iConfig.getParameter<double>("muMinPt");
	elMinPt		= iConfig.getParameter<double>("elMinPt");
	tauMinPt	= iConfig.getParameter<double>("tauMinPt");
	lowInvCut	= iConfig.getParameter<double>("lowInvCut");
	zInvCut		= iConfig.getParameter<double>("zInvCut");
	elIdComp	= iConfig.getParameter<int>("elIdComp");

	metLabel	= iConfig.getParameter<edm::InputTag>("metLabel");
	docoll		= iConfig.getParameter<bool>("doCollinear");

	dR_matching     = 0.01;

	dogen = iConfig.getParameter<bool>("dogen");
	outroot=output;
	HPP_TTree = new TFile(outroot,"RECREATE");
	tree = new TTree("tree","My tree");
	tree ->Branch("run",&run,"run/I");
	tree ->Branch("lumi",&lumi,"lumi/I");
	tree ->Branch("ev",&ev,"ev/I");
        if(dogen) tree ->Branch("gench",&gench,"gench/I");
//        tree ->Branch("genmatch1",&genmatch1,"genmatch1/I");
//        tree ->Branch("genmatch2",&genmatch2,"genmatch2/I");
        if(dogen) tree ->Branch("genmatch",&genmatch,"genmatch/I");
//        tree ->Branch("genmatch4",&genmatch4,"genmatch4/I");
//        tree ->Branch("recoch1",&recoch1,"recoch1/I");
        tree ->Branch("recoch",&recoch,"recoch/I");
//        tree ->Branch("nlep",&nlep,"nlep/I");
//        tree ->Branch("nlgen",&nlgen,"nlgen/I");
        tree ->Branch("Ep",&Ep,"Ep/I");
        tree ->Branch("Em",&Em,"Em/I");
        tree ->Branch("Mp",&Mp,"Mp/I");
        tree ->Branch("Mm",&Mm,"Mm/I");
        tree ->Branch("Tp",&Tp,"Tp/I");
        tree ->Branch("Tm",&Tm,"Tm/I");

//	tree ->Branch("presel",&presel,"presel/I");
//	tree ->Branch("firstpt",&firstpt,"firstpt/F");
//	tree ->Branch("secondpt",&secondpt,"secondpt/F");
	if(dogen) tree ->Branch("npum1",&npum1,"npum1/I");
	if(dogen) tree ->Branch("npu",&npu,"npu/I");
	if(dogen) tree ->Branch("npup1",&npup1,"npup1/I");
	if(dogen) tree ->Branch("Tnpu",&Tnpu,"Tnpu/I");

	tree ->Branch("nvtx",&nvtx,"nvtx/I");
	tree ->Branch("nvtxg",&nvtxg,"nvtxg/I");
//	tree ->Branch("vtxndof",vtxndof,"vtxndof[nvtx]/F");
//	tree ->Branch("vtxnchi2",vtxnchi2,"vtxnchi2[nvtx]/F");
//	tree ->Branch("vtxz",vtxz,"vtxz[nvtx]/F");
//	tree ->Branch("vtxrho",vtxrho,"vtxrho[nvtx]/F");
//	tree ->Branch("vtxleprho",vtxleprho,"vtxleprho[nvtx]/F");

//	if(dogen) tree ->Branch("vtxmcrho",vtxmcrho,"vtxmcrho[nvtx]/F");
//	tree ->Branch("vtxndof_",&vtxndof_,"vtxndof_/F");
//	tree ->Branch("vtxnchi2_",&vtxnchi2_,"vtxnchi2_/F");
//	tree ->Branch("vtxz_",&vtxz_,"vtxz_/F");
//	tree ->Branch("vtxrho_",&vtxrho_,"vtxrho_/F");
//	tree ->Branch("vtxleprho_",&vtxleprho_,"vtxleprho_/F");
//	if(dogen) tree ->Branch("vtxmcrho_",&vtxmcrho_,"vtxmcrho_/F");

//	tree ->Branch("vtxndof2_",&vtxndof2_,"vtxndof2_/F");
//	tree ->Branch("vtxnchi22_",&vtxnchi22_,"vtxnchi22_/F");
//	tree ->Branch("vtxz2_",&vtxz2_,"vtxz2_/F");
//	tree ->Branch("vtxrho2_",&vtxrho2_,"vtxrho2_/F");
//	tree ->Branch("vtxleprho2_",&vtxleprho2_,"vtxleprho2_/F");
//	if(dogen)
//	{
//		tree ->Branch("vtxndofmc_",&vtxndofmc_,"vtxndofmc_/F");
//		tree ->Branch("vtxnchi2mc_",&vtxnchi2mc_,"vtxnchi2mc_/F");
//		tree ->Branch("vtxzmc_",&vtxzmc_,"vtxzmc_/F");
//		tree ->Branch("vtxrhomc_",&vtxrhomc_,"vtxrhomc_/F");
//		tree ->Branch("vtxleprhomc_",&vtxleprhomc_,"vtxleprhomc_/F");
//	}

	tree ->Branch("nm",&nm,"nm/I");
	tree ->Branch("mu_pt",mu_pt,"mu_pt[nm]/F");
	tree ->Branch("mu_eta",mu_eta,"mu_eta[nm]/F");
	tree ->Branch("mu_dz",mu_dz,"mu_dz[nm]/F");
	tree ->Branch("mu_dxy",mu_dxy,"mu_dxy[nm]/F");
	tree ->Branch("mu_rho",mu_rho,"mu_rho[nm]/F");
	tree ->Branch("mu_dB",mu_dB,"mu_dB[nm]/F");
	tree ->Branch("mu_sip",mu_sip,"mu_sip[nm]/F");
	tree ->Branch("mu_nchi2",mu_nchi2,"mu_nchi2[nm]/F");
//	tree ->Branch("mu_trackIso",mu_trackIso,"mu_trackIso[nm]/F");
//	tree ->Branch("mu_ecalIso",mu_ecalIso,"mu_ecalIso[nm]/F");
//	tree ->Branch("mu_hcalIso",mu_hcalIso,"mu_hcalIso[nm]/F");
	tree ->Branch("mu_totIso",mu_totIso,"mu_totIso[nm]/F");
	tree ->Branch("mu_totIso2",mu_totIso2,"mu_totIso2[nm]/F");
//	tree ->Branch("mu_dimass1",mu_dimass1,"mu_dimass1[nm]/F");
//	tree ->Branch("mu_dimass2",mu_dimass2,"mu_dimass2[nm]/F");
//	tree ->Branch("mu_dimass3",mu_dimass3,"mu_dimass3[nm]/F");
//	tree ->Branch("mu_dimass4",mu_dimass4,"mu_dimass4[nm]/F");
	tree ->Branch("mu_glb",mu_glb,"mu_glb[nm]/I");
	tree ->Branch("mu_trk",mu_trk,"mu_trk[nm]/I");
	tree ->Branch("mu_pf", mu_pf, "mu_pf[nm]/I");
	tree ->Branch("mu_good",mu_good,"mu_good[nm]/I");
	tree ->Branch("mu_tight",mu_tight,"mu_tight[nm]/I");
	tree ->Branch("mu_hits_val",mu_hits_val,"mu_hits_val[nm]/I");
//	tree ->Branch("mu_hits_trk",mu_hits_trk,"mu_hits_trk[nm]/I");
//	tree ->Branch("mu_hits_pix",mu_hits_pix,"mu_hits_pix[nm]/I");
//	tree ->Branch("mu_hits_lay",mu_hits_lay,"mu_hits_lay[nm]/I");
	tree ->Branch("mu_nstation",mu_nstation,"mu_nstation[nm]/I");
	if(dogen) tree ->Branch("mu_genmat",mu_genmat,"mu_genmat[nm]/I");

	tree ->Branch("ne",&ne,"ne/I");
	tree ->Branch("el_pt",el_pt,"el_pt[ne]/F");
	tree ->Branch("el_eta",el_eta,"el_eta[ne]/F");
	tree ->Branch("el_dz",el_dz,"el_dz[ne]/F");
	tree ->Branch("el_dxy",el_dxy,"el_dxy[ne]/F");
	tree ->Branch("el_rho",el_rho,"el_rho[ne]/F");
	tree ->Branch("el_dB",el_dB,"el_dB[ne]/F");
	tree ->Branch("el_sip",el_sip,"el_sip[ne]/F");
	tree ->Branch("el_nchi2",el_nchi2,"el_nchi2[ne]/F");
//	tree ->Branch("el_trackIso",el_trackIso,"el_trackIso[ne]/F");
//	tree ->Branch("el_ecalIso",el_ecalIso,"el_ecalIso[ne]/F");
//	tree ->Branch("el_hcalIso",el_hcalIso,"el_hcalIso[ne]/F");
	tree ->Branch("el_totIso",el_totIso,"el_totIso[ne]/F");
	tree ->Branch("el_totIso2",el_totIso2,"el_totIso2[ne]/F");
//	tree ->Branch("el_dimass1",el_dimass1,"el_dimass1[ne]/F");
//	tree ->Branch("el_dimass2",el_dimass2,"el_dimass2[ne]/F");
//	tree ->Branch("el_dimass3",el_dimass3,"el_dimass3[ne]/F");
//	tree ->Branch("el_dimass4",el_dimass4,"el_dimass4[ne]/F");
	tree ->Branch("el_chargeok",el_chargeok,"el_chargeok[ne]/I");
//	tree ->Branch("el_losthit",el_losthit,"el_losthit[ne]/I");
//	tree ->Branch("el_eidVeryLoose",el_eidVeryLoose,"el_eidVeryLoose[ne]/I");
//	tree ->Branch("el_eidLoose",el_eidLoose,"el_eidLoose[ne]/I");
//	tree ->Branch("el_eidMedium",el_eidMedium,"el_eidMedium[ne]/I");
//	tree ->Branch("el_eidTight",el_eidTight,"el_eidTight[ne]/I");
//	tree ->Branch("el_eidSuperTight",el_eidSuperTight,"el_eidSuperTight[ne]/I");
//	tree ->Branch("el_eidHyperTight1",el_eidHyperTight1,"el_eidHyperTight1[ne]/I");
	tree ->Branch("el_mvaTrigV0",el_mvaTrigV0,"el_mvaTrigV0[ne]/F");
	tree ->Branch("el_misshits",el_misshits,"el_misshits[ne]/I");
	tree ->Branch("el_passConvVeto",el_passConvVeto,"el_passConvVeto[ne]/I");
	if(dogen) tree ->Branch("el_genmat",el_genmat,"el_genmat[ne]/I");

	tree ->Branch("nt",&nt,"nt/I");
	tree ->Branch("tau_pt",tau_pt,"tau_pt[nt]/F");
	tree ->Branch("tau_eta",tau_eta,"tau_eta[nt]/F");
	tree ->Branch("tau_signalTracks",tau_signalTracks,"tau_signalTracks[nt]/I");
	tree ->Branch("tau_charge",tau_charge,"tau_charge[nt]/I");
//	tree ->Branch("tau_againstElectronMVA",tau_againstElectronMVA,"tau_againstElectronMVA[nt]/I");
	tree ->Branch("tau_againstElectronLoose",tau_againstElectronLoose,"tau_againstElectronLoose[nt]/I");
//	tree ->Branch("tau_againstElectronMedium",tau_againstElectronMedium,"tau_againstElectronMedium[nt]/I");
//	tree ->Branch("tau_againstElectronTight",tau_againstElectronTight,"tau_againstElectronTight[nt]/I");
	tree ->Branch("tau_againstMuonLoose",tau_againstMuonLoose,"tau_againstMuonLoose[nt]/I");
//	tree ->Branch("tau_againstMuonMedium",tau_againstMuonMedium,"tau_againstMuonMedium[nt]/I");
//	tree ->Branch("tau_againstMuonTight",tau_againstMuonTight,"tau_againstMuonTight[nt]/I");
	tree ->Branch("tau_byLooseCombinedIsolationDeltaBetaCorr",tau_byLooseCombinedIsolationDeltaBetaCorr,"tau_byLooseCombinedIsolationDeltaBetaCorr[nt]/I");
//	tree ->Branch("tau_byMediumCombinedIsolationDeltaBetaCorr",tau_byMediumCombinedIsolationDeltaBetaCorr,"tau_byMediumCombinedIsolationDeltaBetaCorr[nt]/I");
//	tree ->Branch("tau_byTightCombinedIsolationDeltaBetaCorr",tau_byTightCombinedIsolationDeltaBetaCorr,"tau_byTightCombinedIsolationDeltaBetaCorr[nt]/I");
	tree ->Branch("tau_decayModeFinding",tau_decayModeFinding,"tau_decayModeFinding[nt]/I");
	tree ->Branch("tau_overlap",tau_overlap,"tau_overlap[nt]/I");
	if(dogen) tree ->Branch("tau_genmat",tau_genmat,"tau_genmat[nt]/I");
//	tree ->Branch("tau_",tau_,"tau_[nt]/F");


	tree ->Branch("nlp",&nlp,"nlp/I");
	tree ->Branch("nlm",&nlm,"nlm/I");
	tree ->Branch("pid1",&pid1,"pid1/I");
	tree ->Branch("pid2",&pid2,"pid2/I");
	tree ->Branch("pid3",&pid3,"pid3/I");
	tree ->Branch("pid4",&pid4,"pid4/I");
	tree ->Branch("pt1",&pt1,"pt1/F");
	tree ->Branch("pt2",&pt2,"pt2/F");
	tree ->Branch("pt3",&pt3,"pt3/F");
	tree ->Branch("pt4",&pt4,"pt4/F");
	tree ->Branch("eta1",&eta1,"eta1/F");
	tree ->Branch("eta2",&eta2,"eta2/F");
	tree ->Branch("eta3",&eta3,"eta3/F");
	tree ->Branch("eta4",&eta4,"eta4/F");
	tree ->Branch("iso1",&iso1,"iso1/F");
	tree ->Branch("iso2",&iso2,"iso2/F");
	tree ->Branch("iso3",&iso3,"iso3/F");
	tree ->Branch("iso4",&iso4,"iso4/F");
//	tree ->Branch("mp1",&mp1,"mp1/F");
//	tree ->Branch("mn1",&mn1,"mn1/F");
//	tree ->Branch("mp2",&mp2,"mp2/F");
//	tree ->Branch("mn2",&mn2,"mn2/F");
	tree ->Branch("mp",&mp,"mp/F");
	tree ->Branch("mn",&mn,"mn/F");
	tree ->Branch("mpp",&mpp,"mpp/F");
	tree ->Branch("mpm",&mpm,"mpm/F");
	tree ->Branch("mnp",&mnp,"mnp/F");
	tree ->Branch("mnm",&mnm,"mnm/F");
//	tree ->Branch("mp4",&mp4,"mp4/F");
//	tree ->Branch("mn4",&mn4,"mn4/F");
	tree ->Branch("dz",&dz,"dz/F");
	tree ->Branch("dZ",&dZ,"dZ/F");
        tree ->Branch("dilpt",&dilpt,"dilpt/F");
        tree ->Branch("m4",&m4,"m4/F");
        tree ->Branch("m4p",&m4p,"m4p/F");
        tree ->Branch("m4m",&m4m,"m4m/F");
        tree ->Branch("lep_eff_pp",&lep_eff_pp,"lep_eff_pp/F");
        tree ->Branch("lep_eff_pp_high",&lep_eff_pp_high,"lep_eff_pp_high/F");
        tree ->Branch("lep_eff_pp_low",&lep_eff_pp_low,"lep_eff_pp_low/F");
        tree ->Branch("lep_eff_mm",&lep_eff_mm,"lep_eff_mm/F");
        tree ->Branch("lep_eff_mm_high",&lep_eff_mm_high,"lep_eff_mm_high/F");
        tree ->Branch("lep_eff_mm_low",&lep_eff_mm_low,"lep_eff_mm_low/F");
        tree ->Branch("METpt",&METpt,"METpt/F");
        tree ->Branch("mllpm",&mllpm,"mllpm/F");
        tree ->Branch("mllpmp",&mllpmp,"mllpmp/F");
        tree ->Branch("mllpmm",&mllpmm,"mllpmm/F");
        tree ->Branch("mllss",&mllss,"mllss/F");
        tree ->Branch("mZ1",&mZ1,"mZ1/F");
        tree ->Branch("mZ2",&mZ2,"mZ2/F");
        tree ->Branch("mW",&mW,"mW/F");

//	tree ->Branch("HLT",&HLT,"HLT/I");

//	edm::Service<TFileService> fs;
	hmupt = new TH1D("hmupt","P_{T} of muons",50,0,250);
	hmueta = new TH1D("hmueta","#eta of muons",50,-2.4,2.4);
	helpt = new TH1D("helpt","P_{T} of electrons",50,0,250);
	heleta = new TH1D("heleta","#eta of electrons",50,-2.5,2.5);

	hprobe_best_pt = new TH1D("hprobe_best_pt","P_{T} of probe muons",20,0,100);
	hprobe_best_eta = new TH1D("hprobe_best_eta","#eta of probe muons",20,-2.4,2.4);
	hprobe_best_m = new TH1D("hprobe_best_m","Mass of tag+probe",20,0,200);
	hprobe_match_pt = new TH1D("hprobe_match_pt","P_{T} of probe muons",20,0,100);
	hprobe_match_eta = new TH1D("hprobe_match_eta","#eta of probe muons",20,-2.4,2.4);
	hprobe_match_m = new TH1D("hprobe_match_m","Mass of tag+probe",20,0,200);
	hprobe_matmc_pt = new TH1D("hprobe_matmc_pt","P_{T} of probe muons",20,0,100);
	hprobe_matmc_eta = new TH1D("hprobe_matmc_eta","#eta of probe muons",20,-2.4,2.4);
	hprobe_matmc_m = new TH1D("hprobe_matmc_m","Mass of tag+probe",20,0,200);

	for(int i1=0;i1<3;i1++) for(int i2=0;i2<3;i2++) for(int i3=0;i3<3;i3++) for(int i4=0;i4<3;i4++) for(int i5=0;i5<3;i5++) for(int i6=0;i6<3;i6++) channel[i1][i2][i3][i4][i5][i6] = 0;
        channel[0][0][0][0][0][0]=0;                                //MC bkg
        channel[2][0][0][2][0][0]=1;                                //eeee
        channel[2][0][0][1][1][0]=2;  channel[1][1][0][2][0][0]=2;  //eeem
        channel[2][0][0][1][0][1]=3;  channel[1][0][1][2][0][0]=3;  //eeet
        channel[2][0][0][0][2][0]=4;  channel[0][2][0][2][0][0]=4;  //eemm
        channel[2][0][0][0][1][1]=5;  channel[0][1][1][2][0][0]=5;  //eemt
        channel[2][0][0][0][0][2]=6;  channel[0][0][2][2][0][0]=6;  //eett
        channel[1][1][0][1][1][0]=7;                                //emem
        channel[1][1][0][1][0][1]=8;  channel[1][0][1][1][1][0]=8;  //emet
        channel[1][1][0][0][2][0]=9;  channel[0][2][0][1][1][0]=9;  //emmm
        channel[1][1][0][0][1][1]=10; channel[0][1][1][1][1][0]=10; //emmt
        channel[1][1][0][0][0][2]=11; channel[0][0][2][1][1][0]=11; //emtt
        channel[1][0][1][1][0][1]=12;                               //etet
        channel[1][0][1][0][2][0]=13; channel[0][2][0][1][0][1]=13; //etmm
        channel[1][0][1][0][1][1]=14; channel[0][1][1][1][0][1]=14; //etmt
        channel[1][0][1][0][0][2]=15; channel[0][0][2][1][0][1]=15; //etmm
        channel[0][2][0][0][2][0]=16;                               //mmmm
        channel[0][2][0][0][1][1]=17; channel[0][1][1][0][2][0]=17; //mmmt
        channel[0][2][0][0][0][2]=18; channel[0][0][2][0][2][0]=18; //mmtt
        channel[0][1][1][0][1][1]=19;                               //mtmt
        channel[0][1][1][0][0][2]=20; channel[0][0][2][0][1][1]=20; //mttt
        channel[0][0][2][0][0][2]=21;                               //mmmm

	isValidHltConfig_ = false;
	HLT_path = iConfig.getParameter<vector<string> >("HLT_path");

//	f_3DPU = new TFile("/afs/cern.ch/user/j/jslee/CMSSW_4_4_2/src/Weight3D.root");
//	h_3DPU = (TH3D*)f_3DPU->Get("WHist");

	DeltaPt = 0.01;
}


DoublyChargedHiggsPAT::~DoublyChargedHiggsPAT()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

//---------------matching object with MC muon---------------//
double DoublyChargedHiggsPAT::matmcl(TLorentzVector t, TLorentzVector *mcl, int nmcl)
{
        double trmcl=100;
//        for(int i=0;i<4;i++) if(mcl[i].Pt()>10 && fabs(mcl[i].Eta())<2.4 && mcl[i].DeltaR(ttr)<0.01 && ((i<=1&&tr->charge()==1)||(i>=2&&tr->charge()==-1))) trmcl=true;
//        for(int i=0;i<4;i++) if(trmcl>mcl[i].DeltaR(t)) trmcl=mcl[i].DeltaR(t);
        for(int i=0;i<nmcl;i++)
//        for(int i=0;i<int(sizeof(mcl)/sizeof(mcl[0]));i++)
//	int nloop = sizeof(mcl)/sizeof(mcl[];
//        for(int i=0;i<nloop;i++)
	if(mcl[i].Pt()>0)
	{
//		double dR = sqrt(pow(mcl[i].Phi()-t.Phi(),2)+pow(mcl[i].Eta()-t.Eta(),2));
		double dR = t.DeltaR(mcl[i]);
		if(trmcl>dR) trmcl = dR;
	}
        return trmcl;
}

double DoublyChargedHiggsPAT::calcIso_reco(const reco::Candidate &cand, bool useNew = false, bool RealData = false) {
	double iso = -1;
	double efAreaEl[4] = { 0.101, 0.046, 0.021, 0.040 };
	double efAreaMu[4] = { 0.074, 0.045, 0.022, 0.030 };
	double nefAreaEl[4] = { 0.078, 0.046, 0.026, 0.072 };
	double nefAreaMu[4] = { 0.087, 0.049, 0.042, 0.059 };
	if (useNew) 
	  for (unsigned int i=0; i<4; i++) {
		efAreaEl[i]=nefAreaEl[i];
		efAreaMu[i]=nefAreaMu[i];
	  }
	double ecalRho=0, hcalRho=0;
	if (abs(cand.pdgId()) == 11) {
		const pat::Electron *el = dynamic_cast<const pat::Electron*>(&cand);
		if (el->isEB()) {
			ecalRho = efAreaEl[0] * rho;
			hcalRho = efAreaEl[2] * rho;
		} else {
			ecalRho = efAreaEl[1] * rho;
			hcalRho = efAreaEl[3] * rho;
		}
//		iso=(el->hcalIso()+el->ecalIso()+el->trackIso()-ecalRho-hcalRho)/el->pt();
//                iso=(el->chargedHadronIso() + TMath::Max(0.,el->neutralHadronIso()+el->photonIso()-0.5*el->puChargedHadronIso()))/el->pt();
		double scEta = fabs(el->superCluster()->eta());
		double Aeff = 0;
		if(!RealData)
		{
			if(scEta<1.0) Aeff=0.21;
			if(scEta>1.0 && scEta<1.479) Aeff=0.21;
			if(scEta>1.479 && scEta<2.0) Aeff=0.11;
			if(scEta>2.0 && scEta<2.2) Aeff=0.14;
			if(scEta>2.2 && scEta<2.3) Aeff=0.18;
			if(scEta>2.3 && scEta<2.4) Aeff=0.19;
			if(scEta>2.4) Aeff=0.26;
		}
		if(RealData)
		{
			if(scEta<1.0) Aeff=0.13;
			if(scEta>1.0 && scEta<1.479) Aeff=0.14;
			if(scEta>1.479 && scEta<2.0) Aeff=0.07;
			if(scEta>2.0 && scEta<2.2) Aeff=0.09;
			if(scEta>2.2 && scEta<2.3) Aeff=0.11;
			if(scEta>2.3 && scEta<2.4) Aeff=0.11;
			if(scEta>2.4) Aeff=0.14;
		}
//		if(iEvent.isRealData())  el_totIso[nel]= (it->chargedHadronIso()+it->neutralHadronIso()+it->photonIso())/it->pt();
//		if(!iEvent.isRealData()) el_totIso[nel] = (it->chargedHadronIso() + TMath::Max(0.,it->neutralHadronIso()+it->photonIso()-rho*Aeff))/it->pt();
		iso = (el->chargedHadronIso() + TMath::Max(0.,el->neutralHadronIso()+el->photonIso()-rho*Aeff))/el->pt();
	}
	if (abs(cand.pdgId()) == 13) {
		const pat::Muon *mu = dynamic_cast<const pat::Muon*>(&cand);
		if (fabs(mu->eta()) < 1.479) {
			ecalRho = efAreaMu[0] * rho;
			hcalRho = efAreaMu[2] * rho;
		} else {
			ecalRho = efAreaMu[1] * rho;
			hcalRho = efAreaMu[3] * rho;
		}
//		iso=(mu->hcalIso()+mu->ecalIso()+mu->trackIso() - ecalRho - hcalRho)/mu->pt();
                iso=(mu->chargedHadronIso() + TMath::Max(0.,mu->neutralHadronIso()+mu->photonIso()-0.5*mu->puChargedHadronIso()))/mu->pt();
	}
	if (!cand.pdgId() || abs(cand.pdgId()) == 15) {
		//const pat::Tau *tau = dynamic_cast<const pat::Tau*>(&cand);
                //iso=tau->isolationTracksPtSum()/tau->pt();
		iso=0;
	}
	return iso;
}

double DoublyChargedHiggsPAT::calcIso_el(const pat::Electron &el, bool useNew = false) {
	double iso = -1;
	double efAreaEl[4] = { 0.101, 0.046, 0.021, 0.040 };
	double nefAreaEl[4] = { 0.078, 0.046, 0.026, 0.072 };
	if (useNew) 
	  for (unsigned int i=0; i<4; i++) {
		efAreaEl[i]=nefAreaEl[i];
	  }
	double ecalRho=0, hcalRho=0;
	if (el.isEB()) {
		ecalRho = efAreaEl[0] * rho;
		hcalRho = efAreaEl[2] * rho;
	} else {
		ecalRho = efAreaEl[1] * rho;
		hcalRho = efAreaEl[3] * rho;
	}
//	iso=(el.hcalIso()+el.ecalIso()+el.trackIso()-ecalRho-hcalRho)/el.pt();
        iso=(el.chargedHadronIso() + TMath::Max(0.,el.neutralHadronIso()+el.photonIso()-0.5*el.puChargedHadronIso()))/el.pt();
	return iso;
}

double DoublyChargedHiggsPAT::calcIso_mu(const pat::Muon &mu, bool useNew = false) {
	double iso = -1;
	double efAreaMu[4] = { 0.074, 0.045, 0.022, 0.030 };
	double nefAreaMu[4] = { 0.087, 0.049, 0.042, 0.059 };
	if (useNew) 
	  for (unsigned int i=0; i<4; i++) {
		efAreaMu[i]=nefAreaMu[i];
	  }
	double ecalRho=0, hcalRho=0;
	if (fabs(mu.eta()) < 1.479) {
		ecalRho = efAreaMu[0] * rho;
		hcalRho = efAreaMu[2] * rho;
	} else {
		ecalRho = efAreaMu[1] * rho;
		hcalRho = efAreaMu[3] * rho;
	}
//	iso=(mu.hcalIso()+mu.ecalIso()+mu.trackIso() - ecalRho - hcalRho)/mu.pt();
        iso=(mu.chargedHadronIso() + TMath::Max(0.,mu.neutralHadronIso()+mu.photonIso()-0.5*mu.puChargedHadronIso()))/mu.pt();
	return iso;
}

double DoublyChargedHiggsPAT::eaIso(const pat::Electron &el, double rho) {
        double iso = 99;
	double Aeff=0;
//	if(abs(el.superCluster().eta())<1.0) Aeff=0.21;
//	if(abs(el.superCluster().eta())>1.0 && abs(el.superCluster().eta())<1.479) Aeff=0.21;
//	if(abs(el.superCluster().eta())>1.479 && abs(el.superCluster().eta())<2.0) Aeff=0.11;
//	if(abs(el.superCluster().eta())>2.0 && abs(el.superCluster().eta())<2.2) Aeff=0.14;
//	if(abs(el.superCluster().eta())>2.2 && abs(el.superCluster().eta())<2.3) Aeff=0.18;
//	if(abs(el.superCluster().eta())>2.3 && abs(el.superCluster().eta())<2.4) Aeff=0.19;
//	if(abs(el.superCluster().eta())>2.4) Aeff=0.26;
	iso = (el.chargedHadronIso() + TMath::Max(0.,el.neutralHadronIso()+el.photonIso()-rho*Aeff))/el.pt();
	return iso;
}

/*
double DoublyChargedHiggsPAT::calcIso_pat(const pat::Candidate &cand, bool useNew = false) {
	double iso = -1;
	double efAreaEl[4] = { 0.101, 0.046, 0.021, 0.040 };
	double efAreaMu[4] = { 0.074, 0.045, 0.022, 0.030 };
	double nefAreaEl[4] = { 0.078, 0.046, 0.026, 0.072 };
	double nefAreaMu[4] = { 0.087, 0.049, 0.042, 0.059 };
	if (useNew) 
	  for (unsigned int i=0; i<4; i++) {
		efAreaEl[i]=nefAreaEl[i];
		efAreaMu[i]=nefAreaMu[i];
	  }
	double ecalRho=0, hcalRho=0;
	if (abs(cand.pdgId()) == 11) {
		const pat::Electron *el = dynamic_cast<const pat::Electron*>(&cand);
		if (el->isEB()) {
			ecalRho = efAreaEl[0] * rho;
			hcalRho = efAreaEl[2] * rho;
		} else {
			ecalRho = efAreaEl[1] * rho;
			hcalRho = efAreaEl[3] * rho;
		}
		iso=(el->hcalIso()+el->ecalIso()+el->trackIso()-ecalRho-hcalRho)/el->pt();
	}
	if (abs(cand.pdgId()) == 13) {
		const pat::Muon *mu = dynamic_cast<const pat::Muon*>(&cand);
		if (fabs(mu->eta()) < 1.479) {
			ecalRho = efAreaMu[0] * rho;
			hcalRho = efAreaMu[2] * rho;
		} else {
			ecalRho = efAreaMu[1] * rho;
			hcalRho = efAreaMu[3] * rho;
		}
		iso=(mu->isolationR03().sumPt - ecalRho - hcalRho)/mu->pt();
	}
	if (!cand.pdgId() || abs(cand.pdgId()) == 15) {
		//const pat::Tau *tau = dynamic_cast<const pat::Tau*>(&cand);
                //iso=tau->isolationTracksPtSum()/tau->pt();
		iso=0;
	}
	return iso;
}
*/

double DoublyChargedHiggsPAT::Zmass(double pt1, double pt2, double phi1, double phi2, double pz1, double pz2, double e1, double e2)
{
//	double pt1=0, pt2=0, phi1=0, phi2=0, sumpt=0, sumpz=0, sump=0, sume=0;
	double sumpt=0, sumpz=0, sump=0, sume=0;
//	pt1=l1->pt(), pt2=l2->pt(), phi1=l1->phi(), phi2=l2->phi();
	sumpt = sqrt( pow(pt1*TMath::Cos(phi1)+pt2*TMath::Cos(phi2),2) + pow(pt1*TMath::Sin(phi1)+pt2*TMath::Sin(phi2),2) );
	sumpz = pz1+pz2;
	sump = sqrt(pow(sumpt,2)+pow(sumpz,2));
	sume = e1+e2;
	return sqrt( pow(sume,2) - pow(sump,2) );
}

//double DoublyChargedHiggsPAT::lepton_eff(const reco::Candidate &cand)
//double DoublyChargedHiggsPAT::lepton_eff(int lep_pid, double lep_pt, double lep_eta)
double DoublyChargedHiggsPAT::lepton_eff(int lep_pid, double lep_pt, double lep_eta, int lepeff_uncertainty)
{
	double eff = 1;
//	const int nbinx=7, nbiny=10;
//	double range_pt[nbinx+1] = {7, 10, 15, 20, 30, 40, 50, 200};
//	double range_eta[nbiny+1] = {-2.5, -2, -1.566, -1.4442, -0.8, 0, 0.8, 1.4442, 1.566, 2, 2.5};
//	double range_pt[] = {7, 10, 15, 20, 30, 40, 50, 200};
//	double range_eta[] = {-2.5, -2, -1.566, -1.4442, -0.8, 0, 0.8, 1.4442, 1.566, 2, 2.5};
//	double el_eff[nbinx][nbiny] = {
//		{0.86718, 0.86718, 0.869255, 0.908405, 0.908405, 0.908405, 0.908405, 0.869255, 0.86718, 0.86718}, 
//		{0.946243, 0.946243, 0.967838, 0.970138, 0.970138, 0.970138, 0.970138, 0.967838, 0.946243, 0.946243}, 
//		{0.983049, 0.983049, 0.988719, 0.99047, 0.99047, 0.99047, 0.99047, 0.988719, 0.983049, 0.983049}, 
//		{0.976708, 0.98106, 0.992359, 0.984466, 0.985481, 0.985481, 0.984466, 0.992359, 0.98106, 0.976708}, 
//		{0.969604, 0.981627, 0.987532, 0.988428, 0.991399, 0.991399, 0.988428, 0.987532, 0.981627, 0.969604}, 
//		{0.974259, 0.986409, 0.989732, 0.992509, 0.994047, 0.994047, 0.992509, 0.989732, 0.986409, 0.974259}, 
//		{0.971602, 0.98342, 0.987562, 0.990988, 0.990402, 0.990402, 0.990988, 0.987562, 0.98342, 0.971602}};
//	if (abs(cand.pdgId()) == 11)
	if (abs(lep_pid) == 11)
	{
		int el_npt=7, el_neta=10;
		int el_npti=0, el_netai=0;
		double el_range_pt[] = {7, 10, 15, 20, 30, 40, 50, 200};
		double el_range_eta[] = {-2.5, -2, -1.566, -1.4442, -0.8, 0, 0.8, 1.4442, 1.566, 2, 2.5};
		double electron_sf[7][10] = {
		        {0.86718, 0.86718, 0.869255, 0.908405, 0.908405, 0.908405, 0.908405, 0.869255, 0.86718, 0.86718},
		        {0.946243, 0.946243, 0.967838, 0.970138, 0.970138, 0.970138, 0.970138, 0.967838, 0.946243, 0.946243},
		        {0.983049, 0.983049, 0.988719, 0.99047, 0.99047, 0.99047, 0.99047, 0.988719, 0.983049, 0.983049},
		        {0.976708, 0.98106, 0.992359, 0.984466, 0.985481, 0.985481, 0.984466, 0.992359, 0.98106, 0.976708},
		        {0.969604, 0.981627, 0.987532, 0.988428, 0.991399, 0.991399, 0.988428, 0.987532, 0.981627, 0.969604},
		        {0.974259, 0.986409, 0.989732, 0.992509, 0.994047, 0.994047, 0.992509, 0.989732, 0.986409, 0.974259},
		        {0.971602, 0.98342, 0.987562, 0.990988, 0.990402, 0.990402, 0.990988, 0.987562, 0.98342, 0.971602}
		};
		double electron_sf_error[7][10] = {
		        {0.174107, 0.174107, 0.1781, 0.0999565, 0.0999565, 0.0999565, 0.0999565, 0.1781, 0.174107, 0.174107},
		        {0.0459263, 0.0459263, 0.0159751, 0.0165119, 0.0165119, 0.0165119, 0.0165119, 0.0159751, 0.0459263, 0.0459263},
		        {0.0270321, 0.0270321, 0.0346234, 0.0305257, 0.0305257, 0.0305257, 0.0305257, 0.0346234, 0.0270321, 0.0270321},
		        {0.0254421, 0.0171784, 0.0118175, 0.00553195, 0.00619889, 0.00619889, 0.00553195, 0.0118175, 0.0171784, 0.0254421},
		        {0.00377458, 0.00122883, 0.00393825, 0.00190892, 0.000890952, 0.000890952, 0.00190892, 0.00393825, 0.00122883, 0.00377458},
		        {0.00158529, 0.0020278, 0.0101571, 0.00439135, 0.000851254, 0.000851254, 0.00439135, 0.0101571, 0.0020278, 0.00158529},
		        {0.00910625, 0.0020813, 0.0120684, 0.00475305, 0.000874583, 0.000874583, 0.00475305, 0.0120684, 0.0020813, 0.00910625}
		};

		for(int i=0;i<el_npt;i++)
			if((el_range_pt[i]<=lep_pt && el_range_pt[i+1]>lep_pt) || (i==(el_npt-1) && el_range_pt[i]<=lep_pt)) el_npti = i;
		for(int i=0;i<el_neta;i++)
			if(el_range_eta[i]<=lep_eta && el_range_eta[i+1]>lep_eta) el_netai = i;
		eff = electron_sf[el_npti][el_netai];
		if(lepeff_uncertainty==1)  eff = eff*(1+electron_sf_error[el_npti][el_netai]);
		if(lepeff_uncertainty==-1) eff = eff*(1-electron_sf_error[el_npti][el_netai]);
//		double pti = lep_pt, etai = lep_eta;
//		int nbinxi=-1, nbinyi=-1;
//		for(int i=0;i<nbinx;i++)
//		{
//		        if((range_pt[i]<=pti && range_pt[i+1]>pti) || (i==(nbinx-1) && range_pt[i]<=pti)) nbinxi = i;
//		        if(range_pt[i]<=pti && range_pt[i+1]>pti) nbinxi = i;
//		}
//		for(int i=0;i<nbiny;i++)
//		{
//		        if(range_eta[i]<=etai && range_eta[i+1]>etai) nbinyi = i;
//		}
//		if(nbinxi>=0&&nbinyi>=0) eff = el_eff[nbinxi][nbinyi];
	}
//	if (abs(cand.pdgId()) == 13)
	if (abs(lep_pid) == 13)
	{
/*
		int mu_npt=10, mu_neta=4;
		double mu_range_pt[] = {10,20,25,30,35,40,50,60,90,140,300};
		double mu_range_eta[] = {0,0.9,1.2,2.1,2.4};
		double muID_sf[10][4] = {
		        {0.984868,0.986855,1.01235,0.994},
		        {0.988681,0.987375,1.00155,0.994},
		        {0.993889,0.994212,0.999149,0.994},
		        {0.994164,0.990593,0.997573,0.994},
		        {0.994084,0.990353,0.996585,0.994},
		        {0.99247,0.989641,0.997431,0.994},
		        {0.990978,0.991311,0.997521,0.994},
		        {0.990444,0.98631,0.993942,0.994},
		        {1.00385,1.01191,1.01922,0.994},
			{1.02798,0.955563,1.01648,0.994}
		};
		double muID_sf_error_high[10][4] = {
		        {0.00499703,0.00684363,0.00379483,0.00701515},
		        {0.00163837,0.00249344,0.00141401,0.0029204},
		        {0.000769964,0.0013776,0.000856059,0.00178972},
		        {0.000535011,0.00103566,0.000703718,0.00150357},
		        {0.000412461,0.000764144,0.000584787,0.00132915},
		        {0.000265719,0.000477195,0.000238895,0.000982007},
		        {0.000652323,0.00122485,0.000921858,0.00264014},
		        {0.00103232,0.0019409,0.00154479,0.0049719},
		        {0.00313606,0.00619454,0.00535135,0.0103516},
		        {0.0173816,0.0320736,0.0300988,0.124062}
		};
		double muID_sf_error_low[10][4] = {
		        {0.00498232,0.00681157,0.00378731,0.00698443},
		        {0.00164344,0.00250396,0.00141755,0.00292784},
		        {0.000773876,0.00138759,0.000859489,0.00179876},
		        {0.00053711,0.00104393,0.000706583,0.00151252},
		        {0.000414817,0.000769155,0.000587615,0.00133643},
		        {0.000266026,0.000479106,0.000239149,0.00098553},
		        {0.000655459,0.00123488,0.000926239,0.00265544},
		        {0.00103791,0.00196002,0.00155262,0.00498929},
		        {0.00316356,0.00628191,0.00538083,0.0156803},
		        {0.0172591,0.0345151,0.0291727,0.160031}
		};
		double muISO_sf[10][4] = {
		        {0.94705,0.951836,0.980045,1.025},
		        {0.974978,0.988368,0.997342,1.025},
		        {0.997129,1.00083,1.00784,1.025},
		        {0.993863,0.998546,1.00685,1.025},
		        {0.993442,0.99914,1.0037,1.025},
		        {0.994101,0.998176,1.00209,1.025},
		        {0.995544,0.998696,1.00125,1.025},
		        {0.999036,0.999132,1.00065,1.025},
		        {1.00104,0.999559,0.999878,1.025},
			{1.0003,0.996767,0.99989,1.025}
		};
		double muISO_sf_error_high[10][4] = {
		        {0.00431858,0.00479878,0.00231474,0.0055448},
		        {0.00217679,0.00306216,0.00152893,0.00406506},
		        {0.00122609,0.00201228,0.00103338,0.00283551},
		        {0.000797998,0.00143375,0.000790693,0.00215544},
		        {0.000542762,0.000924067,0.000562107,0.00166087},
		        {0.000273166,0.00043507,0.000260123,0.000918489},
		        {0.000497448,0.000794029,0.000479668,0.00166931},
		        {0.000577833,0.000970408,0.000595745,0.00209412},
		        {0.00108927,0.00170772,0.00118927,0.00435686},
		        {0.00247052,0.00432471,0.00235911,0.0149747}
		};
		double muISO_sf_error_low[10][4] = {
		        {0.00432435,0.00482091,0.00232347,0.00556479},
		        {0.00218211,0.00308007,0.00153378,0.00407779},
		        {0.0012289,0.0020214,0.00103712,0.00284741},
		        {0.000800002,0.00144036,0.000794208,0.00216398},
		        {0.000544127,0.000929838,0.000564385,0.00166844},
		        {0.000273836,0.000437463,0.000260438,0.000923371},
		        {0.000500851,0.000805465,0.000484969,0.00169185},
		        {0.000584723,0.000994457,0.000606404,0.00214643},
		        {0.00113296,0.00185362,0.00126147,0.00470127},
		        {0.00280672,0.00558541,0.00320899,0.0208796}
		};
		int mu_npti=0, mu_netai=0;
		for(int i=0;i<mu_npt;i++)
			if((mu_range_pt[i]<=lep_pt && mu_range_pt[i+1]>lep_pt) || (i==(mu_npt-1) && mu_range_pt[i]<=lep_pt)) mu_npti = i;
		for(int i=0;i<mu_neta;i++)
			if(mu_range_eta[i]<=fabs(lep_eta) && mu_range_eta[i+1]>fabs(lep_eta)) mu_netai = i;
		eff = muID_sf[mu_npti][mu_netai]*muISO_sf[mu_npti][mu_netai];
//		if(lepeff_uncertainty==1)  eff = eff*(1+muID_sf_error_high[mu_npti][mu_netai])*(1+muISO_sf_error_high[mu_npti][mu_netai]);
//		if(lepeff_uncertainty==-1) eff = eff*(1-muID_sf_error_low[mu_npti][mu_netai])*(1-muISO_sf_error_low[mu_npti][mu_netai]);
		if(lepeff_uncertainty==1)  eff = eff*(1+sqrt(pow(muID_sf_error_high[mu_npti][mu_netai]),2)+pow(muISO_sf_error_high[mu_npti][mu_netai],2));
		if(lepeff_uncertainty==-1) eff = eff*(1-sqrt(pow(muID_sf_error_high[mu_npti][mu_netai]),2)+pow(muISO_sf_error_high[mu_npti][mu_netai],2));
*/
		double mu_eta = lep_eta, muID_sf_error=1, muISO_sf_error=1;
		if(fabs(mu_eta)<0.9) 
		{
			eff = 0.9925*0.9959;
			muID_sf_error = sqrt(pow(0.005,2)+pow(0.0002,2));
			muISO_sf_error = sqrt(pow(0.002,2)+pow(0.00002,2));
		}
		if(fabs(mu_eta)>0.9 && fabs(mu_eta)<1.2)
		{
			eff = 0.9928*0.9878;
			muID_sf_error = sqrt(pow(0.005,2)+pow(0.0003,2));
			muISO_sf_error = sqrt(pow(0.002,2)+pow(0.0003,2));
		}
		if(fabs(mu_eta)>1.2 && fabs(mu_eta)<2.1)
		{
			eff = 0.9960*1.0027;
			muID_sf_error = sqrt(pow(0.005,2)+pow(0.0003,2));
			muISO_sf_error = sqrt(pow(0.002,2)+pow(0.0002,2));
		}
		if(fabs(mu_eta)>2.1 && fabs(mu_eta)<2.4)
		{
			eff = 0.9952*1.0633;
			muID_sf_error = sqrt(pow(0.005,2)+pow(0.0006,2));
			muISO_sf_error = sqrt(pow(0.002,2)+pow(0.0007,2));
		}
		if(fabs(mu_eta)>=2.4) eff = 1;
		if(lepeff_uncertainty==1)  eff = eff*(1+sqrt(pow(muID_sf_error,2)+pow(muISO_sf_error,2)));
		if(lepeff_uncertainty==-1) eff = eff*(1-sqrt(pow(muID_sf_error,2)+pow(muISO_sf_error,2)));
//		if(fabs(mu_eta)<0.9) eff = 0.9959;
//		if(fabs(mu_eta)>0.9 && fabs(mu_eta)<1.2) eff = 0.9878;
//		if(fabs(mu_eta)>1.2 && fabs(mu_eta)<2.1) eff = 1.0027;
//		if(fabs(mu_eta)>2.1 && fabs(mu_eta)<2.4) eff = 1.0633;
//		if(fabs(mu_eta)>=2.4) eff = 1;
	}
//	if (abs(lep_pid) == 13)
	if (abs(lep_pid) == 15)
	{
		eff = 1;
	}
	return eff;
}

//
// member functions
//

// ------------ method called to for each event  ------------
void
DoublyChargedHiggsPAT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
//	events++; hevents->Fill(1); if(int(events)%1000==0) cout<<"Event "<<events<<endl;
	events++; if(int(events)%1000==0) cout<<"Event "<<events<<endl;
/*
	nlep=0, nlep_acc=0, mP=0,mM=0,eP=0,eM=0,tP=0,tM=0, Mp=0,Mm=0,Ep=0,Em=0,Tp=0,Tm=0, MP=0,MM=0,EP=0,EM=0,TP=0,TM=0;
	Mp1=0,Mm1=0,Ep1=0,Em1=0,Tp1=0,Tm1=0;
	pid1=0, pid2=0, pid3=0, pid4=0;
//	HLT_Mu3=0, HLT_Mu5=0, HLT_Mu7=0, HLT_Mu9=0, HLT_Mu11=0, HLT_IsoMu9=0, HLT_DoubleMu0=0, HLT_DoubleMu3=0, HLT_Ele10_SW_L1R=0, HLT_Ele12_SW_TightEleId_L1R=0, HLT_Ele12_SW_TightEleIdIsol_L1R=0, HLT_Ele12_SW_TightEleIdIsol_NoDEtaInEE_L1R=0, HLT_Ele17_SW_L1R=0, HLT_Ele17_SW_CaloEleId_L1R=0, HLT_Ele17_SW_LooseEleId_L1R=0, HLT_Ele17_SW_EleId_L1R=0, HLT_Ele22_SW_CaloEleId_L1R=0, HLT_Ele40_SW_L1R=0, HLT_DoubleEle4_SW_eeRes_L1R=0, HLT_DoubleEle10_SW_L1R=0, HLT_SingleIsoTau20_Trk5_MET20=0, HLT_SingleIsoTau20_Trk15_MET20=0, HLT_SingleIsoTau30_Trk5_MET20=0, HLT_SingleIsoTau30_Trk5_L120or30=0, HLT_DoubleLooseIsoTau15_OneLeg_Trk5=0, HLT_DoubleLooseIsoTau15_Trk5=0;
	nlp=0, nlm=0, firstpt=0, secondpt=0, pt1=0, pt2=0, pt3=0, pt4=0, iso1=99, iso2=99, iso3=99, iso4=99;
	eta1=99, eta2=99, eta3=99, eta4=99;
	mp1=-1, mn1=-1, mp2=-1, mn2=-1, mp=-1, mn=-1, mp4=-1, mn4=-1;
	mpp=-1, mpm=-1, mnp=-1, mnm=-1, m4p=-1, m4m=-1;
	npum1=0, npu=0, npup1=0, Tnpu=0;
	nvtx=0, nvtxg=0;
	nm=0;
	ne=0;
	nt=0;
	HLT=0;
	lep_eff_pp=1, lep_eff_pp_high=1, lep_eff_pp_low=1;
	lep_eff_mm=1, lep_eff_mm_high=1, lep_eff_mm_low=1;
	METpt=0;
	mllpm=-1; mllpmp=-1; mllpmm=-1; mllss=-1; mZ1=-1; mZ2=-1, mW=-1;
	for(int i=0;i<nMAX;i++) 
	{
		vtxndof[i]=0; vtxnchi2[i]=99; vtxz[i]=99; vtxrho[i]=99; vtxleprho[i]=0; vtxmcrho[i]=0;
		mu_pt[i]=-1, mu_eta[i]=10, mu_dz[i]=99, mu_dxy[i]=99, mu_rho[i]=99, mu_dB[i]=99, mu_sip[i]=99, mu_nchi2[i]=99;
		mu_trackIso[i]=99, mu_ecalIso[i]=99, mu_hcalIso[i]=99, mu_totIso[i]=99, mu_totIso2[i]=99;
		mu_dimass1[i]=999, mu_dimass2[i]=999, mu_dimass3[i]=999, mu_dimass4[i]=999;
		mu_glb[i]=0, mu_trk[i]=0, mu_pf[i]=0,  mu_good[i]=0, mu_tight[i]=0;
		mu_hits_val[i]=0, mu_hits_trk[i]=0, mu_hits_pix[i]=0, mu_hits_lay[i]=0, mu_nstation[i]=0;
		mu_genmat[i]=0;
		el_pt[i]=-1; el_eta[i]=10; el_dz[i]=99; el_dxy[i]=99; el_rho[i]=99; el_dB[i]=99; el_sip[i]=99; el_nchi2[i]=99;
		el_trackIso[i]=99, el_ecalIso[i]=99, el_hcalIso[i]=99, el_totIso[i]=99, el_totIso2[i]=99;
		el_dimass1[i]=999, el_dimass2[i]=999, el_dimass3[i]=999, el_dimass4[i]=999;
		el_chargeok[i]=0; el_losthit[i]=99; 
		el_eidVeryLoose[i]=0; el_eidLoose[i]=0; el_eidMedium[i]=0; el_eidTight[i]=0; el_eidSuperTight[i]=0; el_eidHyperTight1[i]=0; el_mvaTrigV0[i]=0;
		el_passConvVeto[i]=0; el_misshits[i]=0;
		el_genmat[i]=0;
		tau_signalTracks[i]=0, tau_charge[i]=0;
		tau_againstElectronLoose[i]=0, tau_againstElectronMedium[i]=0, tau_againstElectronTight[i]=0, tau_againstElectronMVA[i]=0;
		tau_againstMuonLoose[i]=0, tau_againstMuonMedium[i]=0, tau_againstMuonTight[i]=0;
		tau_decayModeFinding[i]=0;
		tau_overlap[i]=0, tau_genmat[i]=0;
		tau_byLooseCombinedIsolationDeltaBetaCorr[i]=0, tau_byMediumCombinedIsolationDeltaBetaCorr[i]=0, tau_byTightCombinedIsolationDeltaBetaCorr[i]=0;
		tau_pt[i]=0, tau_eta[i]=0;
	}
	vtxndof_=0, vtxnchi2_=99, vtxz_=99, vtxrho_=99, vtxleprho_=99, vtxmcrho_=99;
//	vtxndof2_=0, vtxnchi22_=99, vtxz2_=99, vtxrho2_=99, vtxleprho2_=99;
	vtxndofmc_=0, vtxnchi2mc_=99, vtxzmc_=99, vtxrhomc_=99, vtxleprhomc_=99;
//	presel=0;
	dz=99999, dZ=99999, dilpt=-1, m4=-1;
	gench=0, genmatch1=0, genmatch2=0, genmatch=0, genmatch4=0;
	recoch1=0, recoch=0;
	using namespace edm;
	using namespace std;
//	using namespace reco;
//	using namespace pat;

	run = iEvent.id().run();
	lumi = iEvent.id().luminosityBlock();
	ev = iEvent.id().event();

	if(debug0) cout<<"01"<<endl;

	rho=0;
//	double PUeventReweight = 0;
	if (!iEvent.isRealData()) 
	{
		Handle<vector< PileupSummaryInfo > >  PupInfo;
		iEvent.getByLabel("addPileupInfo", PupInfo);
      	
		vector<PileupSummaryInfo>::const_iterator PVI;
		for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
		{
			if (PVI->getBunchCrossing()==-1) npum1 += PVI->getPU_NumInteractions();
			if (PVI->getBunchCrossing()==0)  npu   += PVI->getPU_NumInteractions();
			if (PVI->getBunchCrossing()==1)  npup1 += PVI->getPU_NumInteractions();
			if (PVI->getBunchCrossing()==0)  Tnpu  += PVI->getTrueNumInteractions();
			if (abs(PVI->getBunchCrossing())>1) cout<<"BunchCrossing too big"<<endl;
		}
	}
	Handle<double> rhoHandle;
//	iEvent.getByLabel(InputTag("kt6corPFJets","rho",""),rhoHandle);
	iEvent.getByLabel(InputTag("kt6PFJets","rho",""),rhoHandle);
	rho = *rhoHandle;
//	cout<<"rho : "<<rho<<endl;

	TLorentzVector LVmet;
	Handle<pat::METCollection> MetC;
	iEvent.getByLabel("patMETs",MetC);
//	cout<<"MET "<<MetC->front.et()<<endl;
//	cout<<"nMetC "<<MetC->size()<<endl;
	for (pat::METCollection::const_iterator met0=MetC->begin(); met0!=MetC->end(); met0++)
	{
		METpt = met0->pt();
		LVmet.SetPtEtaPhiM(met0->pt(),met0->eta(),met0->phi(),0);
//		cout<<"MET : pt "<<met0->pt()<<", et "<<met0->et()<<", sumEt "<<met0->sumEt()<<", sig "<<met0->mEtSig()<<", significance "<<met0->significance()<<", eta "<<met0->eta()<<", phi "<<met0->phi()<<endl;
	}
//{
//--------------------MC gen particles---------------------// 
//--------------------MC gen particles---------------------// 
//--------------------MC gen particles---------------------// 
	TLorentzVector mcmu[4], mcH[2];//, mcH[2];
	TLorentzVector mcl[100], mcme[4], mcm[4], mce[4], mct[4];
	int Muidx=0, MuidxP=0, MuidxM=0, Hidx=0;
	int nmcl=0, nmcme=0, nmcm=0, nmce=0, nmct=0, nmcl_acceptance=0;
	int nmc_pteta=0;
	const Candidate *genlep;
	const Candidate *mclep[100];
	int count_genlep = 0;

//	if(dogen)
	if(!iEvent.isRealData())
	{
		Handle<GenParticleCollection> GenCollection;
		iEvent.getByLabel("genParticles", GenCollection);

		const Candidate * mcmu1=0;
                for (GenParticleCollection::const_iterator gen=GenCollection->begin(); gen!=GenCollection->end(); gen++)
                {
                        const Candidate * mom  = gen->mother();
                        const Candidate *mom2=0;
                        const Candidate *mom3=0;
                        const Candidate *mom4=0;
                        const Candidate *mom5=0;
                        int momid=0, mom2id=0, mom3id=0, mom4id=0, mom5id=0;
                        if(mom!=0) {momid = mom->pdgId(); mom2 = (*mom).mother();} if(mom!=0 && mom2!=0) mom2id = mom2->pdgId();
                        if(mom2!=0) {mom2id = mom2->pdgId(); mom3 = (*mom2).mother();} if(mom2!=0 && mom3!=0) mom3id = mom3->pdgId();
                        if(mom3!=0) {mom3id = mom3->pdgId(); mom4 = (*mom3).mother();} if(mom3!=0 && mom4!=0) mom4id = mom4->pdgId();
                        if(mom4!=0) {mom4id = mom4->pdgId(); mom5 = (*mom4).mother();} if(mom4!=0 && mom5!=0) mom5id = mom5->pdgId();

                        int id = gen->pdgId();  //status=1 : MC final particle
                        int ndau = gen->numberOfDaughters();
			int status = gen->status();
//			if(events<=10) {cout<<mom2id<<", "<<momid<<", "<<id<<" : "; for(int i=0;i<ndau;i++) cout<<gen->daughter(i)->pdgId()<<", "; cout<<endl; }
//			if(events<=100) {cout<<mom5id<<", "<<mom4id<<", "<<mom3id<<", "<<mom2id<<", "<<momid<<", "<<id<<" ("<<status<<")"<<" : "; for(int i=0;i<ndau;i++) cout<<gen->daughter(i)->pdgId()<<", "; cout<<endl; }
//                      if(gen->status()==1 && ndau==0 && abs(id)==13 && Muidx<4)
//			if(events==1&&momid==0&&mom2id==0&&id==2212) cout<<"proton momentum : "<<gen->p()<<endl;
//                      if(ndau==0 && abs(id)==13 && Muidx<4 && abs(momid)!=15 && abs(mom2id)!=15)
			if(ndau==0 && abs(id)==13 && Muidx<4)
//			if( abs(id)==13 && (momid==9900041||momid==-9900041||momid==9900042||momid==-9900042||momid==23||(abs(momid)==24&&abs(mom2id)==6)||abs(momid)==5) && Muidx<4) 
                        {
                                if(Muidx==0) mcmu1 = &(*gen);
                                if(id<0 && MuidxP<2) {mcmu[MuidxP].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->energy()); MuidxP++; Muidx++;}
                                if(id>0 && MuidxM<2) {mcmu[MuidxM+2].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->energy()); MuidxM++; Muidx++;}
                        }
			if( (id==9900041||id==-9900041||id==9900042||id==-9900042) && ndau==3 && Hidx<2)
			{
				if(id>0) mcH[0].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->energy());
				if(id<0) mcH[1].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->energy());
				Hidx++;
			}
                        if(ndau==0 && (abs(momid)==9900041||abs(momid)==9900042||abs(mom2id)==9900041||abs(mom2id)==9900042||abs(mom3id)==9900041||abs(mom3id)==9900042))
                        {
                                if(abs(id)==11) {mce[nmce].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->energy()); nmce++;}
                                if(abs(id)==13) {mcm[nmcm].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->energy()); nmcm++;}
                                if(abs(id)==15) {mct[nmct].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->energy()); nmct++;}
                                if(abs(id)==11||abs(id)==13)
                                {
                                        mcme[nmcme].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->energy());
                                        nmcme++;
//                                        mcl[nmcl].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->p());
//                                        nmcl++;
//                                        if(fabs(gen->eta())<2.4) nmcl_acceptance++;
                                }
//                                if(abs(id)==15)
//                                {
//                                        mcl[nmcl].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->p());
//                                        nmcl++;
//                                        if(fabs(gen->eta())<5.0) nmcl_acceptance++;
//                                }
			}
			if(abs(momid)==9900041||abs(momid)==9900042||abs(momid)==37)
			{
				if(id==-11) {eP++; if(fabs(gen->eta())<2.4) EP++;}
				if(id==11)  {eM++; if(fabs(gen->eta())<2.4) EM++;}
				if(id==-13) {mP++; if(fabs(gen->eta())<2.5) MP++;}
				if(id==13)  {mM++; if(fabs(gen->eta())<2.5) MM++;}
				if(id==-15) {tP++; if(fabs(gen->eta())<5.0) TP++;}
				if(id==15)  {tM++; if(fabs(gen->eta())<5.0) TM++;}
				if(count_genlep==0) genlep = &(*gen);
				count_genlep++;
				//cout<<"ev "<<ev<<", id "<<id<<", "<<gen->vertex().x()<<" "<<gen->vertex().y()<<" "<<gen->vertex().z()<<endl;
                        }
//	                if(ndau==0&&(abs(id)==11||abs(id)==13)&&nmcl<100)
	                if(( (ndau==0&&(abs(id)==11||abs(id)==13)) || (abs(id)==15&&status==2) )&&nmcl<100)
                        {
                                mcl[nmcl].SetPxPyPzE(gen->px(),gen->py(),gen->pz(),gen->energy());
				mclep[nmcl] = &(*gen);
                                nmcl++;
                                if(fabs(gen->eta())<2.4) nmcl_acceptance++;
                                if(abs(id)==11&&gen->pt()>10&&fabs(gen->eta())<2.4) nmc_pteta++;
                                if(abs(id)==13&&gen->pt()>10&&fabs(gen->eta())<2.5) nmc_pteta++;
			}
//if(ndau==0 && (abs(id)==13 || abs(id)==11)) cout<<"Event "<<events<<", id "<<id<<", pt "<<gen->pt()<<", vertex : x "<<gen->vx()<<", y "<<gen->vy()<<", z "<<gen->vz()<<endl;
                }
		nlep_acc = nmcl_acceptance;
		nlgen = nmc_pteta;

		gench = channel[eP][mP][tP][eM][mM][tM];
		if(gench>21) cout<<"Event "<<events<<", gench "<<gench<<" / "<<eP<<mP<<tP<<eM<<mM<<tM<<endl;
	} // dogen
//	cout<<"Event "<<events<<", mcmu pt : "; {for(int i=0;i<4;i++) cout<<mcmu[i].Pt()<<", "; cout<<endl;}

//--------------------global muons---------------------// 
//--------------------global muons---------------------// 
//--------------------global muons---------------------// 
//	iEvent.getByLabel("selectedLayer1Muons",Pgmu);
//	iEvent.getByLabel("selectedPatMuons",Pgmu);
//	iEvent.getByLabel("cleanLayer1Muons",Pgmu);
//	iEvent.getByLabel("selectedLayer1Jets",Pjet);
//	iEvent.getByLabel("selectedPatJets",Pjet);
//	iEvent.getByLabel("cleanLayer1Jets",Pjet);
//	int NGMU=0;
//	for (pat::MuonCollection::const_iterator gmu =Pgmu->begin();  gmu != Pgmu->end(); gmu++)
//	{
//		if(gmu->isGlobalMuon()) NGMU++;
//	}
//	Handle<reco::TrackCollection> TrCollection;
//	iEvent.getByLabel("generalTracks",TrCollection);
	if(debug0) cout<<"02"<<endl;

	if(debug0) cout<<"18"<<endl;
*/
/*
//--------------------Tag and Probe---------------------// 
//--------------------Tag and Probe---------------------// 
//--------------------Tag and Probe---------------------// 
	Handle<reco::TrackCollection> TrCollection;
	iEvent.getByLabel("generalTracks",TrCollection);

//        if(muonsptsortedi2.size()>=1)
//        {       
//		const Muon *gmui2 = & (*Pgmu)[muonsptsortedi2[i].second];
//	}
	if(Zmumu)
	for (pat::MuonCollection::const_iterator gmu =Pgmu->begin();  gmu != Pgmu->end(); gmu++)
	{
		if( gmu->isGlobalMuon() && gmu->pt()>10 && fabs(gmu->eta())<2.4 && (gmu->isolationR03().sumPt)/gmu->pt()<0.1 )
		{
			TLorentzVector ti2(gmu->px(),gmu->py(),gmu->pz(),gmu->p()); // tag
			TLorentzVector pi2; // probe
			bool probe=false;
			double maxpt=0, maxeta=0, maxm=0;
			for (TrackCollection::const_iterator tr=TrCollection->begin(); tr!=TrCollection->end(); tr++)
			{
				TLorentzVector ci2(tr->px(),tr->py(),tr->pz(),tr->p()); // probe candidate
//				cout<<"bool : "<<tr->innerOk();
				if( tr->pt()>10 && fabs(tr->eta())<2.4 && isolPt(ci2,TrCollection)<0.1 && ti2.DeltaPhi(ci2)>2.0 && fabs((ti2+ci2).M()-91.1876)<20 && tr->pt()>maxpt)
				{
					maxpt=tr->pt();
					maxeta=tr->eta();
					maxm=(ti2+ci2).M();
					pi2.SetPxPyPzE(tr->px(),tr->py(),tr->pz(),tr->p()); // probe
					probe=true;
				}
			}
			double probe_match=100;
			if(probe)
			{
				hprobe_best_pt->Fill(maxpt);
				hprobe_best_eta->Fill(maxeta);
				hprobe_best_m->Fill(maxm);
				for (pat::MuonCollection::const_iterator gmu =Pgmu->begin();  gmu != Pgmu->end(); gmu++)
				{
					if( gmu->isGlobalMuon() && gmu->pt()>10 && fabs(gmu->eta())<2.4 && (gmu->isolationR03().sumPt)/gmu->pt()<0.1 )
					{
						TLorentzVector gi2(gmu->px(),gmu->py(),gmu->pz(),gmu->p());
						if(probe_match>gi2.DeltaR(pi2)) probe_match=gi2.DeltaR(pi2);
					}
				}
//cout<<probe_match<<endl;
				if(probe_match<0.01)
				{
					hprobe_match_pt->Fill(maxpt);
					hprobe_match_eta->Fill(maxeta);
					hprobe_match_m->Fill(maxm);
				}
//cout<<matmcl(pi2,&(*mcmu))<<endl;
				if(MCon && matmcl(pi2,&(*mcmu))<0.01)
				{
					hprobe_matmc_pt->Fill(maxpt);
					hprobe_matmc_eta->Fill(maxeta);
					hprobe_matmc_m->Fill(maxm);
				}
			}
		}
	}
*/
/*
	// Collection that will keep the final leptons that are saved out
	auto_ptr<reco::CandidateCollection> selLeps(new reco::CandidateCollection());	
	auto_ptr<reco::CandidateCollection> pLeps(new reco::CandidateCollection());	// plus leptons
	auto_ptr<reco::CandidateCollection> mLeps(new reco::CandidateCollection());	// minus leptons
//	auto_ptr<reco::CandidateCollection> taus(new reco::CandidateCollection());

	// Collection that will keep the temporary pre-selection of leptons for multi lepton cleaning phase
	reco::CandidateCollection tempLeps, tempLepsP, tempLepsM;
	reco::CandidateCollection taus;
	reco::CandidateCollection Lep3;
	vector<unsigned int> ita, itaup, itaum;

	// Read in electrons
	Handle<pat::ElectronCollection> elIn;
	iEvent.getByLabel(elLabel,elIn);
//	iEvent.getByLabel("selectedPatElectrons",elIn);
	// Read in muons
	Handle<pat::MuonCollection> muIn;
	iEvent.getByLabel(muLabel,muIn);
//	iEvent.getByLabel("selectedPatMuons",muIn);
	// Cut on primary vertex
	Handle<vector<reco::Vertex> > pvs;
	iEvent.getByLabel(vertLabel,pvs);
//	if (!pvs->size()) return;
	if(debug0) cout<<"19"<<endl;

	Handle<reco::BeamSpot> bsHandle;
	iEvent.getByLabel("offlineBeamSpot", bsHandle);
	const reco::BeamSpot &beamspot = *bsHandle.product();

//	Handle<pat::ConversionCollection> hConversions;
//	iEvent.getByLabel("patConversions", hConversions);
	Handle<reco::ConversionCollection> hConversions;
	iEvent.getByLabel("allConversions", hConversions);

//	if(events==1)
//	cout<<"# of vertices : "<<pvs->size()<<endl;

	bool foundPV = false;
	int npv=0;
	int pvswap = 0;
	int pvswapok = 0;
	double sptSel = 0;
	double spt = 0;
	double pv_sumpt = 0, pv_sumpt_max = 0;
	double min_rho = 999, min_rho_mc = 999;
	reco::Vertex pv, pv2, pvmc;
	for (vector<reco::Vertex>::const_iterator it = pvs->begin(); it != pvs->end(); it++,nvtx++) {
		if (it->isValid() && it->ndof() > pvndof && fabs(it->z()) < pvmaxz && it->position().rho() < pvmaxd0) {
			spt=0;
//			for (vector<reco::TrackBaseRef>::const_iterator tr = it->tracks_begin(); tr != it->tracks_end(); tr++)
//				spt+=(*tr)->pt();
//				spt+=tr->pt();
//			npv++;
//			if (foundPV && pv.chi2() > it->chi2()) {
//			  pvswap++;
//			  if (spt > sptSel) pvswapok++;
//			}
//			if (!foundPV) {
//			  pv=*it;
//			  sptSel=spt;
//			  foundPV = true;
//			}
		}
		vtxndof[nvtx] = it->ndof();
//		vtxnchi2[nvtx] = it->chi2()/it->ndof();
		vtxnchi2[nvtx] = it->normalizedChi2();
		vtxz[nvtx] = it->z();
		vtxrho[nvtx] = it->position().rho();
		pv_sumpt = 0;
//		if (it->isValid() && (it->chi2()/it->ndof())<6 && fabs(it->z())<20 && it->position().rho()>0.42 && it->position().rho()<0.51)
		if (it->isValid() && (it->ndof())>4 && fabs(it->z())<24 && it->position().rho()<2)
		{
			nvtxg++;
			for (vector<reco::TrackBaseRef>::const_iterator tr = it->tracks_begin(); tr != it->tracks_end(); tr++) 
			{
				pv_sumpt += (*tr)->pt();
			}
			if(pv_sumpt_max < pv_sumpt)
			{
				pv=*it;
				pv_sumpt_max = pv_sumpt;
			}
		}
//		double vtx_lep_rho = 0, count_lep = 0;
//		for (pat::ElectronCollection::const_iterator el = elIn->begin(); el != elIn->end(); el++)
//		{
//			if(el->pt()>10 && fabs(el->eta())<2.5) vtx_lep_rho += (it->position() - el->vertex()).rho();
//			if(el->pt()>10 && fabs(el->eta())<2.5) count_lep++;
//		}
//		for (pat::MuonCollection::const_iterator mu = muIn->begin(); mu != muIn->end(); mu++)
//		{
//			if(mu->isGlobalMuon() && mu->pt()>10 && fabs(mu->eta())<2.4) vtx_lep_rho += (it->position() - mu->vertex()).rho();
//			if(mu->isGlobalMuon() && mu->pt()>10 && fabs(mu->eta())<2.4) count_lep++;
//		}
//		if(count_lep>=1) vtx_lep_rho = vtx_lep_rho/count_lep;
//		if(count_lep>=1 && vtx_lep_rho==0) vtx_lep_rho = 998;
//		vtxleprho[nvtx] = vtx_lep_rho;
//		if(min_rho>vtx_lep_rho) 
//		{
//			min_rho = vtx_lep_rho;
//			vtxleprho_ = min_rho;
//			pv=*it;
//		}
//		for (reco::Vertex::trackRef_iterator tr = it->tracks_begin(); tr != it->tracks_end(); it++) pv_sumpt += (*tr)->pt();
//		if(pv_sumpt_max<pv_sumpt)
//		{
//			pv2=*it;
//			vtxleprho2_ = vtx_lep_rho;
//		}
////		if(dogen)
//		if(!iEvent.isRealData())
//		{
//			double vtx_lep_rho_mc = 0, count_lep_mc = 0;
////			for(int i=0;i<nmcl;i++)
////			{
////				if(abs(mclep[i]->pdgId())==11 && mclep[i]->pt()>10 && fabs(mclep[i]->eta())<2.5) vtx_lep_rho_mc += (it->position() - mclep[i]->vertex()).rho();
////				if(abs(mclep[i]->pdgId())==11 && mclep[i]->pt()>10 && fabs(mclep[i]->eta())<2.5) count_lep_mc++;
////				if(abs(mclep[i]->pdgId())==13 && mclep[i]->pt()>10 && fabs(mclep[i]->eta())<2.4) vtx_lep_rho_mc += (it->position() - mclep[i]->vertex()).rho();
////				if(abs(mclep[i]->pdgId())==13 && mclep[i]->pt()>10 && fabs(mclep[i]->eta())<2.4) count_lep_mc++;
////			}
////			if(count_lep_mc>=1 && vtx_lep_rho_mc==0) vtx_lep_rho_mc = 999;
//			if(count_genlep>0) vtx_lep_rho_mc = (it->position() - genlep->vertex()).rho();
//			vtxmcrho[nvtx] = vtx_lep_rho_mc;
//			if(min_rho_mc>vtx_lep_rho_mc) 
//			{
//				min_rho_mc = vtx_lep_rho_mc;
//				vtxmcrho_ = min_rho_mc;
//				pvmc=*it;
//				vtxleprhomc_ = vtx_lep_rho;
//			}
//		}
	}
	if(nvtx>=1)
	{
		vtxndof_ = pv.ndof();
//		vtxnchi2_ = pv.chi2()/pv.ndof();
		vtxnchi2_ = pv.normalizedChi2();
		vtxz_ = pv.z();
		vtxrho_ = pv.position().rho();
//		vtxndof2_ = pv2.ndof();
//		vtxnchi22_ = pv2.normalizedChi2();
//		vtxz2_ = pv2.z();
//		vtxrho2_ = pv2.position().rho();
	}
//	if(nvtx>=1 && dogen)
	if(nvtx>=1 && !iEvent.isRealData())
	{
		vtxndofmc_ = pvmc.ndof();
		vtxnchi2mc_ = pvmc.normalizedChi2();
		vtxzmc_ = pvmc.z();
		vtxrhomc_ = pvmc.position().rho();
	}
	if(debug0) cout<<"20"<<endl;

	int lposition[3][1000]={{0,0,0},}, pposition[3][1000]={{0,0,0},}, mposition[3][1000]={{0,0,0},};
	unsigned int nel=0, nmu=0, ntau=0, nelp=0, nmup=0, ntaup=0;
	unsigned int nelP=0, nmuP=0, ntauP=0, nelM=0, nmuM=0, ntauM=0;
	// Read in electrons
//	Handle<pat::ElectronCollection> elIn;
//	iEvent.getByLabel(elLabel,elIn);

//	hc_["filtEff"]->Fill(kElStart,elIn->size()); 

	for (pat::ElectronCollection::const_iterator it = elIn->begin(); it != elIn->end(); it++,nel++)	 {
//cout<<(it->chargedHadronIso() + TMath::Max(0.,it->neutralHadronIso()+it->photonIso()-0.5*it->puChargedHadronIso()))<<endl;
		if(firstpt<it->pt()) {secondpt=firstpt; firstpt=it->pt();}
		else if(secondpt<it->pt()) secondpt=it->pt();
		el_pt[nel] = it->pt();
		el_eta[nel] = it->eta();
		el_dz[nel] = fabs( it->vertex().z() - pv.z() );
		el_dxy[nel] = it->gsfTrack()->dxy(pv.position());
		el_rho[nel] = (it->vertex() - pv.position()).rho();
		el_dB[nel] = fabs(it->dB());
		el_sip[nel] = fabs(it->dB(pat::Electron::PV3D))/fabs(it->edB(pat::Electron::PV3D));
		el_nchi2[nel] = it->gsfTrack()->chi2()/it->gsfTrack()->ndof();
//cout<<it->gsfTrack()->chi2()<<", "<<it->gsfTrack()->ndof()<<", "<<el_nchi2[nel]<<endl;
		el_chargeok[nel] = it->isGsfCtfScPixChargeConsistent();
		el_losthit[nel] = it->gsfTrack()->lost();
//		el_eidVeryLoose[nel] = ((int)it->electronID("eidVeryLoose")&1);
//		el_eidLoose[nel] = ((int)it->electronID("eidLoose")&1);
		el_mvaTrigV0[nel] = it->electronID("mvaTrigV0");
//		el_eidMedium[nel] = ((int)it->electronID("eidMedium")&1);
//		el_eidTight[nel] = ((int)it->electronID("eidTight")&1);
//		el_eidSuperTight[nel] = ((int)it->electronID("eidSuperTight")&1);
//		el_eidHyperTight1[nel] = ((int)it->electronID("eidHyperTight1")&1);
		el_misshits[nel] = it->gsfTrack()->trackerExpectedHitsInner().numberOfHits();
//		bool passconversionveto = !ConversionTools::hasMatchedConversion(dynamic_cast<reco::GsfElectron const&>(*(it->originalObjectRef())),hConversions,pv.position(),true, 2.0, 1e-6, 0);
		bool passconversionveto = !ConversionTools::hasMatchedConversion(dynamic_cast<reco::GsfElectron const&>(*(it->originalObjectRef())),hConversions,beamspot.position(),true, 2.0, 1e-6, 0);
		if(passconversionveto) el_passConvVeto[nel] = 1;
		else el_passConvVeto[nel] = 0;
//cout<<"ev "<<events<<", nel "<<nel+1<<", miss-hits "<<el_misshits[nel]<<", passconversionveto "<<passconversionveto<<endl;
//cout<<"	eidVeryLoose "<<it->electronID("eidLoose")<<endl;
//		if(it->isElectronIDAvailable("eidLoose")) cout<<"eidLoose "<<it->electronID("eidLoose")<<endl;
//		if(it->isElectronIDAvailable("eidTight")) cout<<"eidTight "<<it->electronID("eidTight")<<endl;
//cout<<"	eidSuperTight "<<it->electronID("eidSuperTight")<<endl;
//cout<<"	eidHyperTight1 "<<it->electronID("eidHyperTight1")<<endl;
//		if(it->isElectronIDAvailable("mvaTrigV0")) cout<<"mvaTrigV0 "<<it->electronID("mvaTrigV0")<<endl;
//cout<<"	mvaID "<<it->electronID("mvaID")<<endl;
//cout<<"	mvaTrigV0 "<<it->electronID("mvaTrigV0")<<endl;
//cout<<"	mvaNonTrigV0 "<<it->electronID("mvaNonTrigV0")<<endl;
//cout<<"	eidVBTFRel95 "<<it->electronID("eidVBTFRel95")<<endl;
//cout<<"	eidCiCHZZMedium "<<it->electronID("eidCiCHZZMedium")<<endl;
//cout<<"	eidCiCMedium "<<it->electronID("eidCiCMedium")<<endl;
//cout<<"	eid "<<it->electronID("eid")<<endl;
//		el_trackIso[nel] = it->trackIso()/it->pt();
//		el_ecalIso[nel] = it->ecalIso()/it->pt();
//		el_hcalIso[nel] = it->hcalIso()/it->pt();
//		if(iEvent.isRealData())  el_totIso[nel] = (it->trackIso()+it->ecalIso()+it->hcalIso())/it->pt();
//		if(!iEvent.isRealData()) el_totIso[nel] = calcIso_el(*it);
//		el_totIso[nel] =  (it->chargedHadronIso() + TMath::Max(0.,it->neutralHadronIso()+it->photonIso()-0.5*it->puChargedHadronIso()))/it->pt();
		el_totIso2[nel] = (it->chargedHadronIso() + TMath::Max(0.,it->neutralHadronIso()+it->photonIso()-0.5*it->puChargedHadronIso()))/it->pt();
		double scEta = fabs(it->superCluster()->eta());
		double Aeff = 0;
		if(!iEvent.isRealData())
		{
			if(scEta<1.0) Aeff=0.21;
			if(scEta>1.0 && scEta<1.479) Aeff=0.21;
			if(scEta>1.479 && scEta<2.0) Aeff=0.11;
			if(scEta>2.0 && scEta<2.2) Aeff=0.14;
			if(scEta>2.2 && scEta<2.3) Aeff=0.18;
			if(scEta>2.3 && scEta<2.4) Aeff=0.19;
			if(scEta>2.4) Aeff=0.26;
		}
		if(iEvent.isRealData())
		{
			if(scEta<1.0) Aeff=0.13;
			if(scEta>1.0 && scEta<1.479) Aeff=0.14;
			if(scEta>1.479 && scEta<2.0) Aeff=0.07;
			if(scEta>2.0 && scEta<2.2) Aeff=0.09;
			if(scEta>2.2 && scEta<2.3) Aeff=0.11;
			if(scEta>2.3 && scEta<2.4) Aeff=0.11;
			if(scEta>2.4) Aeff=0.14;
		}
//		if(iEvent.isRealData())  el_totIso[nel]= (it->chargedHadronIso()+it->neutralHadronIso()+it->photonIso())/it->pt();
//		if(!iEvent.isRealData()) el_totIso[nel] = (it->chargedHadronIso() + TMath::Max(0.,it->neutralHadronIso()+it->photonIso()-rho*Aeff))/it->pt();
		el_totIso[nel] = (it->chargedHadronIso() + TMath::Max(0.,it->neutralHadronIso()+it->photonIso()-rho*Aeff))/it->pt();
//cout<<el_totIso[nel]<<", "<<el_totIso2[nel]<<", "<<(it->trackIso()+it->ecalIso()+it->hcalIso())/it->pt()<<", pu "<<it->puChargedHadronIso()<<endl;
		bool genmatch = false;
//		if(dogen)
		if(!iEvent.isRealData())
			for(int i=0;i<nmcl;i++)
				if( abs(mclep[i]->pdgId())==11 && sqrt(pow(it->eta()-mclep[i]->eta(),2) + pow(deltaphi(it->phi(),mclep[i]->phi()),2))<0.01 && it->charge()==mclep[i]->charge() ) 
					genmatch = true;
		if(genmatch) el_genmat[nel] = 1;
//		if (it->pt() > elMinPt && ((int)it->electronID(elId) & elIdComp) == elIdComp && fabs(it->eta()) < 2.5 && it->gsfTrack()->dxy(pv.position()) < 0.2 && it->isGsfCtfScPixChargeConsistent()) 
//		if (abs(it->pdgId()) == 11 && it->pt() > elMinPt && ((int)it->electronID(elId) & elIdComp) == elIdComp && fabs(it->eta()) < 2.5 && it->isGsfCtfScPixChargeConsistent()) 
//		if (abs(it->pdgId()) == 11 && it->pt() > elMinPt && fabs(it->eta()) < 2.5) 
//		if (abs(it->pdgId()) == 11 && it->pt() > elMinPt && fabs(it->eta()) < 2.5 && el_chargeok[nel] && el_eidMedium[nel]) 
//		if (abs(it->pdgId()) == 11 && it->pt() > elMinPt && fabs(it->eta()) < 2.5 && fabs(el_dxy[nel])<0.02 && el_sip[nel]<4 && el_nchi2[nel]<10 && el_chargeok[nel] && el_eidLoose[nel] ) 
//		if (abs(it->pdgId()) == 11 && it->pt() > 5 && fabs(it->eta()) < 2.5 && el_chargeok[nel] && el_eidLoose[nel] ) 
//		if (abs(it->pdgId()) == 11 && it->pt() > 10 && fabs(it->eta()) < 2.5 && el_chargeok[nel] && el_mvaTrigV0[nel] > 0.5) 
//		if (abs(it->pdgId()) == 11 && it->pt() > 20 && ( (fabs(it->eta())<0.8&&el_mvaTrigV0[nel]>-0.34) || (fabs(it->eta())>0.8&&fabs(it->eta())<1.479&&el_mvaTrigV0[nel]>-0.65) || (fabs(it->eta())>1.479&&fabs(it->eta())<2.5&&el_mvaTrigV0[nel]>0.60) )&&el_sip[nel]<4&&el_misshits[nel]<=1&&el_passConvVeto[nel]) 
		if (abs(it->pdgId()) == 11 && it->pt() > 20 && ( (fabs(it->eta())<0.8&&el_mvaTrigV0[nel]>-0.34) || (fabs(it->eta())>0.8&&fabs(it->eta())<1.479&&el_mvaTrigV0[nel]>-0.65) || (fabs(it->eta())>1.479&&fabs(it->eta())<2.5&&el_mvaTrigV0[nel]>0.60) )&&el_sip[nel]<4&&el_misshits[nel]<=1&&el_passConvVeto[nel]&&el_totIso[nel]<0.4) 
		{
//			hc2_["etoPV"]->Fill((pv.position() - it->vertex()).rho(), (pv.position() - it->vertex()).R());
			tempLeps.push_back(*it);
//			hc_["filtEff"]->Fill(kElBasic);
//			lposition[0][nel]=1;
			lposition[0][nelp]=nel; nelp++;
			if(it->charge()== 1) {pposition[0][nelP]=nel; nelP++; tempLepsP.push_back(*it);}
			if(it->charge()==-1) {mposition[0][nelM]=nel; nelM++; tempLepsM.push_back(*it);}
			if(debug2&&events<1000) cout<<"e"<<nel<<", ";
			helpt->Fill(it->pt());
			heleta->Fill(it->eta());
//			if(debug2&&events<1000) cout<<"e"<<nel*it->charge()<<", ";
			unsigned int nel2=0;
//			for (pat::ElectronCollection::const_iterator it2 = elIn->begin(); it2 != elIn->end(); it2++,nel2++)
//				if (abs(it2->pdgId()) == 11 && it2->pt() > elMinPt && fabs(it2->eta()) < 2.5 && fabs(it2->gsfTrack()->dxy(pv.position()))<0.02 && fabs(it2->dB(pat::Electron::PV3D))/fabs(it2->edB(pat::Electron::PV3D))<4 && (it2->gsfTrack()->chi2()/it2->gsfTrack()->ndof())<10 && it2->isGsfCtfScPixChargeConsistent() && ((int)it2->electronID("eidLoose")&1) && nel!=nel2)
////				if (abs(it2->pdgId()) == 11 && it2->pt() > elMinPt && fabs(it2->eta()) < 2.5 && fabs(it2->gsfTrack()->dxy(pv.position()))<0.02 && ((int)it2->electronID("eidMedium")&1) && nel!=nel2)
//				{
//					if(el_dimass1[nel] > (it->p4()+it2->p4()).M()) el_dimass1[nel] = (it->p4()+it2->p4()).M();
//					if(el_dimass2[nel] > (it->p4()+it2->p4()).M() && (it->charge() == it2->charge())) el_dimass2[nel] = (it->p4()+it2->p4()).M();
//					if(el_dimass3[nel] > (it->p4()+it2->p4()).M() && (it->charge() != it2->charge())) el_dimass3[nel] = (it->p4()+it2->p4()).M();
//				}
//			el_dimass4[nel] = el_dimass1[nel];
//			for (pat::MuonCollection::const_iterator it3 = muIn->begin(); it3 != muIn->end(); it3++)
//				if(abs(it3->pdgId()) == 13 && it3->pt() > muMinPt && fabs(it3->eta()) < 2.4 && it3->isGlobalMuon() && it3->globalTrack().isNonnull() && it3->normChi2()<10 && it3->globalTrack()->dxy(pv.position())<0.02 && fabs(it3->dB(pat::Muon::PV3D))/fabs(it3->edB(pat::Muon::PV3D))<4 && muon::isGoodMuon(*it3, muon::GlobalMuonPromptTight) && it3->globalTrack()->numberOfValidHits() > 20 && it3->numberOfMatchedStations() >= 1)
//				{
//					if(el_dimass4[nel] > (it->p4()+it3->p4()).M()) el_dimass4[nel] = (it->p4()+it3->p4()).M();
//				}
		}
	}
	ne = nel;
//	if(elIn->size()>=1) {const pat::Electron *el = & (*elIn)[0]; cout<<el->hcalIso()<<", "<<el->ecalIso()<<", "<<el->trackIso()<<endl;}
	if(debug0) cout<<"21"<<endl;
	// Read in muons
//	Handle<pat::MuonCollection> muIn;
//	iEvent.getByLabel(muLabel,muIn);

//	hc_["filtEff"]->Fill(kMuStart,muIn->size()); 
	bool good = false;
	for (pat::MuonCollection::const_iterator it = muIn->begin(); it != muIn->end(); it++,nmu++) {
//	if(debug0) cout<<"21.1"<<endl;
//if((it->pfIsolationR04().sumChargedHadronPt + it->pfIsolationR04().sumNeutralHadronEt + it->pfIsolationR04().sumPhotonEt + it->pfIsolationR04().sumPUPt)>0) cout<<it->pfIsolationR04().sumChargedHadronPt<<", "<<it->pfIsolationR04().sumNeutralHadronEt<<", "<<it->pfIsolationR04().sumPhotonEt<<", "<<it->pfIsolationR04().sumPUPt<<endl;
//cout<<it->chargedHadronIso()<<", "<<it->neutralHadronIso()<<", "<<it->photonIso()<<", "<<it->puChargedHadronIso()<<endl;
//cout<<(it->chargedHadronIso() + TMath::Max(0.,it->neutralHadronIso()+it->photonIso()-0.5*it->puChargedHadronIso()))<<endl;
//	if(debug0) cout<<"21.5"<<endl;
//cout<<"isisGlobalMuon "<<it->isGlobalMuon()<<", isTrackerMuon "<<it->isTrackerMuon()<<", isPFMuon "<<it->isPFMuon()<<endl;
		if(firstpt<it->pt()) {secondpt=firstpt; firstpt=it->pt();}
		else if(secondpt<it->pt()) secondpt=it->pt();
		mu_pt[nmu] = it->pt();
		mu_eta[nmu] = it->eta();
		if(it->isGlobalMuon()) mu_glb[nmu] = 1;
		if(it->isTrackerMuon()) mu_trk[nmu] = 1;
		if(it->isPFMuon()) mu_pf[nmu] = 1;
		if(muon::isGoodMuon(*it, muon::GlobalMuonPromptTight)) mu_good[nmu] = 1;
		if(muon::isTightMuon(*it, pv)) mu_tight[nmu] = 1;
		mu_dz[nmu] = fabs( it->vertex().z() - pv.z() );
		if(it->globalTrack().isNonnull()) 
		{
			mu_dxy[nmu] = it->globalTrack()->dxy(pv.position());
			mu_hits_val[nmu] = it->globalTrack()->numberOfValidHits();
			mu_nchi2[nmu] = it->normChi2();
			mu_hits_trk[nmu] = it->globalTrack()->hitPattern().numberOfValidTrackerHits();
			mu_hits_pix[nmu] = it->globalTrack()->hitPattern().numberOfValidPixelHits();
			mu_hits_lay[nmu] = it->globalTrack()->hitPattern().trackerLayersWithMeasurement();
//cout<<mu_hits_val[nmu]<<", "<<mu_hits_trk[nmu]<<", "<<mu_hits_pix[nmu]<<", "<<mu_hits_trk[nmu]+mu_hits_pix[nmu]<<", "<<it->globalTrack()->hitPattern().numberOfValidMuonHits()<<endl;
		}
		mu_rho[nmu] = (it->vertex() - pv.position()).rho();
		double mu_rho_trk = 99;
		if(it->track().isNonnull())
//		if(muon::isGoodMuon(*it, muon::GlobalMuonPromptTight) && it->track().isNonnull())
		{
//			cout<<it->track()->numberOfValidHits()<<endl;
//			mu_rho_trk = it->track()->dxy(pv.position());
//			mu_hits_trk[nmu] = it->track()->hitPattern().numberOfValidTrackerHits();
//			mu_hits_pix[nmu] = it->track()->hitPattern().numberOfValidPixelHits();
//			mu_hits_lay[nmu] = it->track()->hitPattern().trackerLayersWithMeasurement();
		}
//cout<<mu_dxy[nmu]<<", "<<mu_rho_trk<<", "<<mu_rho[nmu]<<" // "<<endl;
		mu_nstation[nmu] = it->numberOfMatchedStations();
		mu_dB[nmu] = fabs(it->dB());
		mu_sip[nmu] = fabs(it->dB(pat::Muon::PV3D))/fabs(it->edB(pat::Muon::PV3D));
//		mu_trackIso[nmu] = it->trackIso()/it->pt();
//		mu_ecalIso[nmu] = it->ecalIso()/it->pt();
//		mu_hcalIso[nmu] = it->hcalIso()/it->pt();
//		if(iEvent.isRealData())  mu_totIso[nmu] = (it->trackIso()+it->ecalIso()+it->hcalIso())/it->pt();
//		if(!iEvent.isRealData()) mu_totIso[nmu] = calcIso_mu(*it);
//		if(iEvent.isRealData()) mu_totIso[nmu]= (it->chargedHadronIso()+it->neutralHadronIso()+it->photonIso())/it->pt();
//		if(!iEvent.isRealData()) mu_totIso[nmu] = (it->chargedHadronIso() + TMath::Max(0.,it->neutralHadronIso()+it->photonIso()-0.5*it->puChargedHadronIso()))/it->pt();
		mu_totIso[nmu] = (it->chargedHadronIso() + TMath::Max(0.,it->neutralHadronIso()+it->photonIso()-0.5*it->puChargedHadronIso()))/it->pt();
		mu_totIso2[nmu]= (it->chargedHadronIso()+it->neutralHadronIso()+it->photonIso())/it->pt();
//cout<<mu_dxy[nmu]<<", "<<mu_rho[nmu]<<", "<<it->dB(pat::Muon::PV2D)<<endl;
		bool genmatch = false;
//		if(dogen)
		if(!iEvent.isRealData())
			for(int i=0;i<nmcl;i++)
				if( abs(mclep[i]->pdgId())==13 && sqrt(pow(it->eta()-mclep[i]->eta(),2) + pow(deltaphi(it->phi(),mclep[i]->phi()),2))<0.01 && it->charge()==mclep[i]->charge() ) 
					genmatch = true;
		if(genmatch) mu_genmat[nmu] = 1;
//		if (it->pt() > muMinPt && fabs(it->eta()) < 2.4 && it->isGlobalMuon() && muon::isGoodMuon(*it, muon::GlobalMuonPromptTight) && it->globalTrack().isNonnull() && it->globalTrack()->numberOfValidHits() > 10 && fabs(it->globalTrack()->dxy(pv.position())) < 0.2 ) {
//		if (abs(it->pdgId()) == 13 && it->pt() > muMinPt && fabs(it->eta()) < 2.4 && it->isGlobalMuon() && muon::isGoodMuon(*it, muon::GlobalMuonPromptTight) && it->globalTrack().isNonnull() && it->globalTrack()->numberOfValidHits() > 10) {
//		if (abs(it->pdgId()) == 13 && it->pt() > muMinPt && fabs(it->eta()) < 2.4 && it->isGlobalMuon()) {
//		if (abs(it->pdgId())==13 && it->pt() > muMinPt && fabs(it->eta())<2.4 && it->isGlobalMuon() && it->globalTrack().isNonnull() && mu_nchi2[nmu]<10 && mu_dxy[nmu]<0.02 && mu_sip[nmu]<4 && mu_good[nmu] && mu_hits_val[nmu]>20 && mu_nstation[nmu]>0)
//		if (abs(it->pdgId())==13 && it->pt() > muMinPt && fabs(it->eta())<2.4 && it->isGlobalMuon() && it->globalTrack().isNonnull() && mu_good[nmu] && mu_hits_val[nmu]>20 && mu_nstation[nmu]>0)
//		if (abs(it->pdgId())==13 && it->pt() > 5 && fabs(it->eta())<2.4 && it->isGlobalMuon() && it->globalTrack().isNonnull() && mu_good[nmu] && mu_hits_val[nmu]>20 && mu_nstation[nmu]>0)
//		if (abs(it->pdgId())==13 && it->pt() > 10 && fabs(it->eta())<2.4 && it->isPFMuon() && it->globalTrack().isNonnull() && mu_good[nmu] && mu_hits_val[nmu]>20 && mu_nstation[nmu]>0)
//		if (abs(it->pdgId())==13 && it->pt() > 10 && fabs(it->eta())<2.4 && mu_glb[nmu] && it->isPFMuon() && it->globalTrack().isNonnull() && mu_good[nmu] && mu_nstation[nmu]>0)
//		if (abs(it->pdgId())==13 && it->pt() > 10 && fabs(it->eta())<2.4 && mu_glb[nmu] && it->isPFMuon() && it->globalTrack().isNonnull() && mu_tight[nmu] && mu_nstation[nmu]>0)
//		if (abs(it->pdgId())==13 && it->pt() > 10 && fabs(it->eta())<2.4 && mu_tight[nmu] && mu_sip[nmu]<4)
//		if (abs(it->pdgId())==13 && it->pt() > 20 && fabs(it->eta())<2.4 && mu_tight[nmu])
		if (abs(it->pdgId())==13 && it->pt() > 20 && fabs(it->eta())<2.4 && mu_tight[nmu] && mu_totIso[nmu]<0.12)
		{
//			hc2_["mtoPV"]->Fill((pv.position() - it->vertex()).rho(), (pv.position() - it->vertex()).R());
			tempLeps.push_back(*it);
//			hc_["filtEff"]->Fill(kMuBasic);
//			lposition[1][nmu]=1;
			lposition[1][nmup]=nmu;	nmup++;
			if(it->charge()== 1) {pposition[1][nmuP]=nmu; nmuP++; tempLepsP.push_back(*it);}
			if(it->charge()==-1) {mposition[1][nmuM]=nmu; nmuM++; tempLepsM.push_back(*it);}
			hmupt->Fill(it->pt());
			hmueta->Fill(it->eta());
//			if(debug2) cout<<"m"<<nmu*it->charge()<<", ";
			if(debug2&&events<1000) cout<<"m"<<nmu<<", ";
			unsigned int nmu2=0;
//			for (pat::MuonCollection::const_iterator it2 = muIn->begin(); it2 != muIn->end(); it2++,nmu2++)
////				if(abs(it2->pdgId()) == 13 && it2->pt() > muMinPt && fabs(it2->eta()) < 2.4 && it2->isGlobalMuon() && it2->globalTrack().isNonnull() && it2->normChi2()<10 && it2->globalTrack()->dxy(pv.position())<0.02 && fabs(it2->dB(pat::Muon::PV3D))/fabs(it2->edB(pat::Muon::PV3D))<4 && muon::isGoodMuon(*it2, muon::GlobalMuonPromptTight) && it2->globalTrack()->numberOfValidHits() > 20 && it2->numberOfMatchedStations() >= 1 && nmu!=nmu2)
//				if(abs(it2->pdgId()) == 13 && it2->pt() > muMinPt && fabs(it2->eta()) < 2.4 && it2->isPFMuon() && it2->globalTrack().isNonnull() && muon::isGoodMuon(*it2, muon::GlobalMuonPromptTight) && it2->globalTrack()->numberOfValidHits() > 20 && it2->numberOfMatchedStations() > 0 && nmu!=nmu2)
////				if(abs(it2->pdgId()) == 13 && it2->pt() > muMinPt && fabs(it2->eta()) < 2.4 && it2->isGlobalMuon() && it2->globalTrack().isNonnull() && muon::isGoodMuon(*it2, muon::GlobalMuonPromptTight) && it2->globalTrack()->numberOfValidHits() > 20 && it2->numberOfMatchedStations() >= 1 && nmu!=nmu2)
//				{
//					if(mu_dimass1[nmu] > (it->p4()+it2->p4()).M()) mu_dimass1[nmu] = (it->p4()+it2->p4()).M();
//					if(mu_dimass2[nmu] > (it->p4()+it2->p4()).M() && (it->charge() == it2->charge())) mu_dimass2[nmu] = (it->p4()+it2->p4()).M();
//					if(mu_dimass3[nmu] > (it->p4()+it2->p4()).M() && (it->charge() != it2->charge())) mu_dimass3[nmu] = (it->p4()+it2->p4()).M();
//				}
//			mu_dimass4[nmu] = mu_dimass1[nmu];
//			for (pat::ElectronCollection::const_iterator it3 = elIn->begin(); it3 != elIn->end(); it3++)
//				if (abs(it3->pdgId()) == 11 && it3->pt() > elMinPt && fabs(it3->eta()) < 2.5 && fabs(it3->gsfTrack()->dxy(pv.position()))<0.02 && fabs(it3->dB(pat::Electron::PV3D))/fabs(it3->edB(pat::Electron::PV3D))<4 && (it3->gsfTrack()->chi2()/it3->gsfTrack()->ndof())<10 && it3->isGsfCtfScPixChargeConsistent() && ((int)it3->electronID("eidLoose")&1))
////				if (abs(it3->pdgId()) == 11 && it3->pt() > elMinPt && fabs(it3->eta()) < 2.5 && fabs(it3->gsfTrack()->dxy(pv.position()))<0.02 && ((int)it3->electronID("eidMedium")&1))
//				{
//					if(mu_dimass4[nmu] > (it->p4()+it3->p4()).M()) mu_dimass4[nmu] = (it->p4()+it3->p4()).M();
//				}
		}
	}
	nm = nmu;
	if(debug0) cout<<"22"<<endl;
//	if(nmuP>=2&&nmuM>=2) cout<<pposition[1][0]<<" "<<pposition[1][1]<<" "<<mposition[1][0]<<" "<<mposition[1][1]<<endl;
//	if(muIn->size()>=1) {const pat::Muon *gmuP = & (*muIn)[0]; cout<<gmuP->isolationR03().sumPt<<endl;}
	// Read in taus
	Handle<pat::TauCollection> tauIn;
	iEvent.getByLabel(tauLabel,tauIn);
      
//	hc_["filtEff"]->Fill(kTauStart,tauIn->size()); 
	bool overlap = false;
	if(intau)
	for (pat::TauCollection::const_iterator it = tauIn->begin(); it != tauIn->end(); it++,ntau++) {
//		if (debug) {
//			cout << "Tau: " << it->pdgId() << ": " << it->isJet() << ": " << it->p4() << " q: " << it->charge() << " trk: " << it->signalTracks().size();
//			cout << " against: (" << it->tauID("againstElectronMVA") << "," << it->tauID("againstMuonTight") << "," << it->tauID("leadingTrackFinding") << "," << it->tauID("leadingTrackPtCut");
//			cout << "," << it->tauID("trackIsolation") << "," << it->tauID("ecalIsolation") << ")" << " MC: " << it->genJet() << endl;
//		}
		if(it->pt() > 20 && fabs(it->eta()) < 2.4)
		{
			tau_pt[ntau] = it->pt();
			tau_eta[ntau] = it->eta();
			tau_signalTracks[ntau] = it->signalTracks().size();
			tau_charge[ntau] = it->charge();
			if(it->tauID("againstElectronLoose")) tau_againstElectronLoose[ntau] = 1;
//			if(it->tauID("againstElectronMedium")) tau_againstElectronMedium[ntau] = 1;
//			if(it->tauID("againstElectronTight")) tau_againstElectronTight[ntau] = 1;
//			if(it->tauID("againstElectronMVA")) tau_againstElectronMVA[ntau] = 1;
			if(it->tauID("againstMuonLoose")) tau_againstMuonLoose[ntau] = 1;
//			if(it->tauID("againstMuonMedium")) tau_againstMuonMedium[ntau] = 1;
//			if(it->tauID("againstMuonTight")) tau_againstMuonTight[ntau] = 1;
//			if(it->tauID("againstElectronMVA") || it->tauID("againstMuonTight")) tau_overlap[ntau] = 1;
			if(it->tauID("byLooseCombinedIsolationDeltaBetaCorr")) tau_byLooseCombinedIsolationDeltaBetaCorr[ntau] = 1;
//			if(it->tauID("byMediumCombinedIsolationDeltaBetaCorr")) tau_byMediumCombinedIsolationDeltaBetaCorr[ntau] = 1;
//			if(it->tauID("byTightCombinedIsolationDeltaBetaCorr")) tau_byTightCombinedIsolationDeltaBetaCorr[ntau] = 1;
			if(it->tauID("decayModeFinding")) tau_decayModeFinding[ntau] = 1;
		}
		bool genmatch = false;
//		if(dogen)
		if(!iEvent.isRealData())
			for(int i=0;i<nmcl;i++)
//				if( abs(mclep[i]->pdgId())==15 && sqrt(pow(it->eta()-mclep[i]->eta(),2) + pow(deltaphi(it->phi(),mclep[i]->phi()),2))<0.1 && it->charge()==mclep[i]->charge() ) 
				if( it->pt() > 20 && fabs(it->eta()) < 2.4 && abs(mclep[i]->pdgId())==15 && sqrt(pow(it->eta()-mclep[i]->eta(),2) + pow(deltaphi(it->phi(),mclep[i]->phi()),2))<0.1 && it->charge()==mclep[i]->charge() ) 
					genmatch = true;
		if(genmatch) tau_genmat[ntau] = 1;
//		if (it->pt() > tauMinPt && ( it->signalTracks().size() == 1 || it->signalTracks().size() == 3 ) && 
//			fabs(it->eta()) < 2.1 && ( fabs(it->eta()) < 1.46 || fabs(it->eta()) > 1.558 ) && abs(it->charge()) == 1 &&
//			it->tauID("againstElectronMVA") && it->tauID("againstMuonTight") && it->tauID("byTaNCfrOnePercent") ) {
//		if (it->pt() > tauMinPt && fabs(it->eta()) < 2.1 && ( fabs(it->eta()) < 1.46 || fabs(it->eta()) > 1.558 ) && abs(it->charge())==1 && tau_overlap[ntau])
		overlap = false;
		if (it->pt() > 20 && fabs(it->eta()) < 2.4 && abs(it->charge())==1 && (tau_againstElectronLoose[ntau] || tau_againstMuonLoose[ntau]) && tau_byLooseCombinedIsolationDeltaBetaCorr[ntau] && tau_decayModeFinding[ntau])
		{
			for (unsigned int l = 0; l < tempLeps.size(); l++) 
			  if (reco::deltaR(tempLeps[l].eta(),tempLeps[l].phi(),it->eta(),it->phi()) < 0.1) 
				overlap = true;
//			if (!overlap) 
			{
			  tempLeps.push_back(*it);
//			  hc_["filtEff"]->Fill(kTauBasic);
			  lposition[2][ntaup]=ntau; ntaup++;
			  if(it->charge()== 1) {pposition[2][ntauP]=ntau; ntauP++; tempLepsP.push_back(*it);}
			  if(it->charge()==-1) {mposition[2][ntauM]=ntau; ntauM++; tempLepsM.push_back(*it);}
			  if(debug2&&events<1000) cout<<"t"<<ntau<<", ";
			  if(debug2&&events<1000) cout<<"t"<<ntau*it->charge()<<", ";
			}
		}
		if(!overlap) tau_overlap[ntau] = 1;
	}
	nt = ntau;
	if(debug0) cout<<"23"<<endl;
	if(debug2&&events<1000)
	{
		if((nelP+nelM+nmuP+nmuM+ntauP+ntauM)>0) cout<<" ==> passed lepton id 1"<<endl;
		for(unsigned int i=0;i<nelP;i++) cout<<"e"<<pposition[0][i]<<", ";
		for(unsigned int i=0;i<nelM;i++) cout<<"e"<<mposition[0][i]<<", ";
		for(unsigned int i=0;i<nmuP;i++) cout<<"m"<<pposition[1][i]<<", ";
		for(unsigned int i=0;i<nmuM;i++) cout<<"m"<<mposition[1][i]<<", ";
		for(unsigned int i=0;i<ntauP;i++) cout<<"t"<<pposition[2][i]<<", ";
		for(unsigned int i=0;i<ntauM;i++) cout<<"t"<<mposition[2][i]<<", ";
		if((nelP+nelM+nmuP+nmuM+ntauP+ntauM)>0) cout<<" ==> passed lepton id 2"<<endl;
	}
	if((nelp+nmup+ntaup)!=tempLeps.size()) cout<<"check for # of temp leptons 1"<<endl;
	if((nelP+nelM+nmuP+nmuM+ntauP+ntauM)!=(nelp+nmup+ntaup)) cout<<"check for # of temp leptons 2"<<endl;
	int Lposition[3][1000]={{0,0,0},}, Pposition[3][1000]={{0,0,0},}, Mposition[3][1000]={{0,0,0},};
	int Nel=0, Nmu=0, Ntau=0, Pel=0, Pmu=0, Ptau=0, Mel=0, Mmu=0, Mtau=0;
	unsigned int npl=0, nml=0;
//	double mZ=14000;
	int selnum[4] = {-1,-1,-1,-1};
	double selpt[4] = {0,};
	for(int i=0;i<4;i++) {
		for (unsigned int l1 = 0; l1 < tempLeps.size(); l1++) {
			bool matchnum = false;
			for(int j=0;j<4;j++) if((int)l1 == selnum[j]) matchnum = true;
			if(!matchnum && selpt[i] < tempLeps[l1].pt()) {
				selnum[i] = l1;
				selpt[i] = tempLeps[l1].pt();
			}
		}
	}
	int selnumP[2] = {-1,-1};
	double selptP[2] = {0,};
	for(int i=0;i<2;i++) {
		for (unsigned int l1 = 0; l1 < tempLepsP.size(); l1++) {
			bool matchnum = false;
			for(int j=0;j<2;j++) if((int)l1 == selnumP[j]) matchnum = true;
			if(!matchnum && selptP[i] < tempLepsP[l1].pt()) {
				selnumP[i] = l1;
				selptP[i] = tempLepsP[l1].pt();
			}
		}
	}
	int selnumM[2] = {-1,-1};
	double selptM[2] = {0,};
	for(int i=0;i<2;i++) {
		for (unsigned int l1 = 0; l1 < tempLepsM.size(); l1++) {
			bool matchnum = false;
			for(int j=0;j<2;j++) if((int)l1 == selnumM[j]) matchnum = true;
			if(!matchnum && selptM[i] < tempLepsM[l1].pt()) {
				selnumM[i] = l1;
				selptM[i] = tempLepsM[l1].pt();
			}
		}
	}
//	if(tempLepsP.size()>=2 && tempLepsM.size()>=2) {
//		for (unsigned int l1 = 0; l1 < tempLeps.size(); l1++) cout<<tempLeps[l1].pt()*tempLeps[l1].charge()<<" "; cout<<" // ";
//		for(int i=0;i<2;i++) cout<<tempLepsP[selnumP[i]].pt()*tempLepsP[selnumP[i]].charge()<<" ";
//		for(int i=0;i<2;i++) cout<<tempLepsM[selnumM[i]].pt()*tempLepsM[selnumM[i]].charge()<<" "; 
//		cout<<endl;
//	}
//	for (unsigned int l1 = 0; l1 < tempLeps.size(); l1++) cout<<lepsptsorted[l1].first<<" "<<lepsptsorted[l1].second<<", "; cout<<end;
	// Perform cleaning from b-decays
	for (unsigned int l1 = 0; l1 < tempLeps.size(); l1++) {
		// New lepton so set checks to false
//		if(abs(tempLeps[l1].pdgId())!=11&&abs(tempLeps[l1].pdgId()!=13)&&ntest<10) {cout<<"pdgId "<<tempLeps[l1].pdgId()<<", charge "<<tempLeps[l1].charge()<<endl; ntest++;}
		lowMass = false;
		zPole = false;
//		if (tempLeps[l1].pdgId()) {
			// Iterate over all leptons (except the one we're running on) and check both low mass cutoff as well as zPole closeness
			// Also, this is only run in case it's e or mu, no point for taus...
//			for (unsigned int l2 = 0; l2 < tempLeps.size(); l2++) //, l2 != l1; l2++)
//				if ( l2 != l1 ) {
//					if ( (tempLeps[l1].p4() + tempLeps[l2].p4()).M() < lowInvCut ) lowMass = true; 
//					if ( fabs((tempLeps[l1].p4() + tempLeps[l2].p4()).M() - 91.188) < zInvCut && tempLeps[l1].pdgId() == -tempLeps[l2].pdgId() ) zPole = true;
//					if ( tempLeps[l1].pdgId() == -tempLeps[l2].pdgId() && fabs(mZ - 91.188) > fabs((tempLeps[l1].p4() + tempLeps[l2].p4()).M() - 91.188)) {
//						pt1 = tempLeps[l1].pt();
//						pt2 = tempLeps[l2].pt();
//						iso1 = calcIso_reco(tempLeps[l1]);
//						iso2 = calcIso_reco(tempLeps[l2]);
//						mZ=(tempLeps[l1].p4() + tempLeps[l2].p4()).M();
//						mp=mZ;
//						int pid=abs(tempLeps[l1].pdgId());
//						if(pid==11) {Ep=1; Em=1; Mp=0; Mm=0; Tp=0; Tm=0;}
//						if(pid==13) {Ep=0; Em=0; Mp=1; Mm=1; Tp=0; Tm=0;}
//						if(pid==15) {Ep=0; Em=0; Mp=0; Mm=0; Tp=1; Tm=1;}
//					}
//				}
//		}
		if (!lowMass && !zPole) {
			selLeps->push_back(tempLeps[l1]);
			if(tempLeps[l1].charge()== 1) pLeps->push_back(tempLeps[l1]);
			if(tempLeps[l1].charge()==-1) mLeps->push_back(tempLeps[l1]);
			if(l1<nelp) {Lposition[0][Nel]=lposition[0][l1]; Nel++;}
			if(l1>=nelp&&l1<(nelp+nmup)) {Lposition[1][Nmu]=lposition[1][l1-nelp]; Nmu++;}
			if(l1>=(nelp+nmup)&&l1<(nelp+nmup+ntaup))
			{
				Lposition[2][Ntau]=lposition[2][l1-nelp-nmup];
//				taus.push_back((*tauIn)[Lposition[2][Ntau]]);
				if(tempLeps[l1].charge()== 1) itaup.push_back(Pel+Pmu+Ptau);
				if(tempLeps[l1].charge()==-1) itaum.push_back(Mel+Mmu+Mtau);
				Ntau++;
			}
//			if(abs(tempLeps[l1].pdgId())==15) {taus->push_back(tempLeps[l1]); itau.push_back(PMtau); PMtau++;}
			if(tempLeps[l1].charge()== 1)
			{
				if(npl<nelP) {Pposition[0][Pel]=pposition[0][npl]; Pel++;}// cout<<pposition[0][npl]<<" ";}
				if(npl>=nelP&&npl<(nelP+nmuP)) {Pposition[1][Pmu]=pposition[1][npl-nelP]; Pmu++;}// cout<<"(npl-nelP "<<npl-nelP<<") "<<pposition[1][npl-nelP]<<" ";}
				if(npl>=(nelP+nmuP)) {Pposition[2][Ptau]=pposition[2][npl-nelP-nmuP]; Ptau++;}// cout<<pposition[2][npl-nelP-nmuP]<<" ";}
			}
			if(tempLeps[l1].charge()==-1)
			{
				if(nml<nelM) {Mposition[0][Mel]=mposition[0][nml]; Mel++;}// cout<<mposition[0][nml]<<" ";}
				if(nml>=nelM&&nml<(nelM+nmuM)) {Mposition[1][Mmu]=mposition[1][nml-nelM]; Mmu++;}// cout<<"(nml-nelM "<<nml-nelM<<") "<<mposition[1][nml-nelM]<<" ";}
				if(nml>=(nelM+nmuM)) {Mposition[2][Mtau]=mposition[2][nml-nelM-nmuM]; Mtau++;}// cout<<mposition[2][nml-nelM-nmuM]<<" ";}
			}
		}
//		if (tempLeps[l1].pdgId()<0) npl++;
//		if (tempLeps[l1].pdgId()>0) nml++;
		if(tempLeps[l1].charge()== 1) npl++;
		if(tempLeps[l1].charge()==-1) nml++;
//		if(abs(tempLeps[l1].pdgId())==13) cout<<tempLeps[l1].pfIsolationR04().sumChargedHadronPt<<", "<<tempLeps[l1].pfIsolationR04().sumNeutralHadronEt<<", "<<tempLeps[l1].pfIsolationR04().sumPhotonEt<<", "<<tempLeps[l1].pfIsolationR04().sumPUPt<<endl;
	}
	if(debug0) cout<<"24"<<endl;
	if (debug) cout << "Final selection content: " << selLeps->size() << endl;	
	if(debug2&&events<1000)
	{
		for(int i=0;i<Pel;i++) cout<<"e"<<Pposition[0][i]<<", ";
		for(int i=0;i<Mel;i++) cout<<"e"<<Mposition[0][i]<<", ";
		for(int i=0;i<Pmu;i++) cout<<"m"<<Pposition[1][i]<<", ";
		for(int i=0;i<Mmu;i++) cout<<"m"<<Mposition[1][i]<<", ";
		for(int i=0;i<Ptau;i++) cout<<"t"<<Pposition[2][i]<<", ";
		for(int i=0;i<Mtau;i++) cout<<"t"<<Mposition[2][i]<<", ";
		if((Pel+Mel+Pmu+Mmu+Ptau+Mtau)>0) cout<<" ==> preselection"<<endl;
	}
	if(tempLeps.size()==2)
		for (reco::CandidateCollection::const_iterator l1 = tempLeps.begin(); l1 != tempLeps.end()-1; l1++) 
			for (reco::CandidateCollection::const_iterator l2 = l1+1; l2 != tempLeps.end(); l2++) 
			{
				if ( l1->charge() != l2->charge() && l1->pdgId() == -l2->pdgId() ) 
				{
					pt1 = l1->pt();
					pt2 = l2->pt();
					eta1 = l1->eta();
					eta2 = l2->eta();
//					if(iEvent.isRealData())  iso1 = calcIso(*l1);
//					if(iEvent.isRealData())  iso2 = calcIso(*l2);
//					if(!iEvent.isRealData()) iso1 = calcIso_reco(*l1);
//					if(!iEvent.isRealData()) iso2 = calcIso_reco(*l2);
					iso1 = calcIso_reco(*l1, false, iEvent.isRealData());
					iso2 = calcIso_reco(*l2, false, iEvent.isRealData());
//					mp=(l1->p4() + l2->p4()).M();
					mllpm=(l1->p4() + l2->p4()).M();
					int pid=abs(l1->pdgId());
					if(pid==11) {Ep=1; Em=1;}
					if(pid==13) {Mp=1; Mm=1;}
					if(pid==15) {Tp=1; Tm=1;}
					lep_eff_pp = lepton_eff(pid,l1->pt(),l1->eta(),0)*lepton_eff(pid,l2->pt(),l2->eta(),0);
					lep_eff_pp_high = lepton_eff(pid,l1->pt(),l1->eta(),1)*lepton_eff(pid,l2->pt(),l2->eta(),1);
					lep_eff_pp_low = lepton_eff(pid,l1->pt(),l1->eta(),-1)*lepton_eff(pid,l2->pt(),l2->eta(),-1);
//cout<<"ev "<<events<<", pid "<<pid<<", pt1 "<<pt1<<", eta1 "<<l1->eta()<<", eff "<<lepton_eff(pid,l1->pt(),l1->eta())<<", pt2 "<<pt2<<", eta2 "<<l2->eta()<<", eff "<<lepton_eff(pid,l2->pt(),l2->eta())<<", sf "<<lep_eff_pp<<endl;
//					double pt1=0, pt2=0, phi1=0, phi2=0, sumpt=0, sumpz=0, sump=0, sume=0;
//					pt1=l1->pt(), pt2=l2->pt(), phi1=l1->phi(), phi2=l2->phi();
//					sumpt = sqrt( pow(pt1*TMath::Cos(phi1)+pt2*TMath::Cos(phi2),2) + pow(pt1*TMath::Sin(phi1)+pt2*TMath::Sin(phi2),2) );
//					sumpz = l1->pz()+l2->pz();
//					sump = sqrt(pow(sumpt,2)+pow(sumpz,2));
//					sume = l1->energy()+l2->energy();
//					double mpp=0, mpm=0;
//					mpp = Zmass(l1->pt()*(1+DeltaPt),l2->pt()*(1+DeltaPt),l1->phi(),l2->phi(),l1->pz(),l2->pz(),l1->energy(),l2->energy());
//					mpm = Zmass(l1->pt()*(1-DeltaPt),l2->pt()*(1-DeltaPt),l1->phi(),l2->phi(),l1->pz(),l2->pz(),l1->energy(),l2->energy());

//					reco::Candidate::LorentzVector  L1(l1->px()*DeltaPt,l1->py()*DeltaPt,0,0), L2(l2->px()*DeltaPt,l2->py()*DeltaPt,0,0);
//					reco::Candidate::LorentzVector mp_vp = l1->p4()+l2->p4()+L1+L2, mp_vm = l1->p4()+l2->p4()-L1-L2;
//					mpp = mp_vp.mass();
//					mpm = mp_vm.mass();
					TLorentzVector L1p, L1m, L2p, L2m;
					L1p.SetPtEtaPhiM(l1->pt()*(1+DeltaPt),l1->eta(),l1->phi(),0);
					L1m.SetPtEtaPhiM(l1->pt()*(1-DeltaPt),l1->eta(),l1->phi(),0);
					L2p.SetPtEtaPhiM(l2->pt()*(1+DeltaPt),l2->eta(),l2->phi(),0);
					L2m.SetPtEtaPhiM(l2->pt()*(1-DeltaPt),l2->eta(),l2->phi(),0);
//					mpp = (L1p+L2p).M();
//					mpm = (L1m+L2m).M();
					mllpmp = (L1p+L2p).M();
					mllpmm = (L1m+L2m).M();
//cout<<"ev "<<events<<", mp "<<mp<<", mpp "<<mpp<<", mpp2 "<<(L1p+L2p).M()<<", mpm "<<mpm<<", mpm2 "<<(L1m+L2m).M()<<endl;
//cout<<"ev "<<events<<", mp "<<mp<<", mpp "<<mpp<<", mpp2 "<<v1p.mass()<<" mpm "<<mpm<<", mpm2 "<<v1m.mass()<<endl;

//cout<<"ev "<<events<<", mp "<<mp<<", mp2 "<<Zmass(l1->pt(),l2->pt(),l1->phi(),l2->phi(),l1->pz(),l2->pz(),l1->energy(),l2->energy())<<", mpp "<<mpp<<", mpm "<<mpm<<endl;
//cout<<"ev "<<events<<", l1 p4 "<<l1->p4()<<", l2 p4 "<<l2->p4()<<", mp "<<mp<<", pt "<<(l1->p4()+l2->p4()).pt()<<", pt2 "<<sumpt<<endl;
//cout<<"ev "<<events<<", l1 p4 "<<l1->p4()<<", l2 p4 "<<l2->p4()<<", mp "<<mp<<", mp2 "<<sqrt( pow(sume,2) - pow(sump,2) )<<endl;
//cout<<"ev "<<events<<", l1 p4 "<<l1->p4()<<", energy "<<l1->energy()<<", mass "<<l1->p4().M()<<", mass2 "<<sqrt(pow(l1->energy(),2)-pow(l1->p(),2))<<endl;
//cout<<"ev "<<events<<", l1 p4 "<<l1->p4()<<", pt "<<l1->pt()<<", pz "<<l1->pz()<<", p "<<l1->p()<<", p2 "<<sqrt(pow(l1->pt(),2)+pow(l1->pz(),2))<<endl;
				}
				if ( l1->charge() == l2->charge() && l1->pdgId() == l2->pdgId() ) 
				{
					pt1 = l1->pt();
					pt2 = l2->pt();
					eta1 = l1->eta();
					eta2 = l2->eta();
					mllss=(l1->p4() + l2->p4()).M();
					int pid=abs(l1->pdgId());
					if(pid==11&&l1->charge()==1)  {Ep=2; mp=mllss;}
					if(pid==11&&l1->charge()==-1) {Em=2; mn=mllss;}
					if(pid==13&&l1->charge()==1)  {Mp=2; mp=mllss;}
					if(pid==13&&l1->charge()==-1) {Mm=2; mn=mllss;}
					if(pid==15&&l1->charge()==1)  {Tp=2; mp=mllss;}
					if(pid==15&&l1->charge()==-1) {Tm=2; mn=mllss;}
					lep_eff_pp = lepton_eff(pid,l1->pt(),l1->eta(),0)*lepton_eff(pid,l2->pt(),l2->eta(),0);
					lep_eff_pp_high = lepton_eff(pid,l1->pt(),l1->eta(),1)*lepton_eff(pid,l2->pt(),l2->eta(),1);
					lep_eff_pp_low = lepton_eff(pid,l1->pt(),l1->eta(),-1)*lepton_eff(pid,l2->pt(),l2->eta(),-1);
				}
			}
//	if(debug2) cout<<"e size "<<elIn->size()<<", m size "<<muIn->size()<<", t size "<<tauIn->size()<<" ==> no cut"<<endl;
//	if(debug2) cout<<"e size "<<nelp<<", m size "<<nmup<<", t size "<<ntaup<<" ==> lepton id"<<endl;
//	if(debug2) cout<<"e size "<<Nel<<", m size "<<Nmu<<", t size "<<Ntau<<" ==> preselection"<<endl;
//	hc_["nLeps"]->Fill(selLeps->size());
//	if (debug) cout << "So far events that contribtued to nLeps: " << hc_["nLeps"]->Integral() << endl;
//	iEvent.put(selLeps);
	nlep=selLeps->size();
//	if(selLeps->size()>=3&&pvs->size()) nlep=selLeps->size();
//	if(pvs->size()) nlep=selLeps->size();
//	if(selLeps->size()<3) return;
//cout<<"nlp "<<pLeps->size()<<", nlm "<<mLeps->size()<<" : ";
	// Check for Z veto
	double dztmp = 9999;
	if(selLeps->size()>=3)
	for (reco::CandidateCollection::const_iterator l1 = selLeps->begin(); l1 != selLeps->end()-1; l1++) 
		for (reco::CandidateCollection::const_iterator l2 = l1+1; l2 != selLeps->end(); l2++) {
			if ( l1->pdgId() && l1->pdgId() == -l2->pdgId() ) {
				dztmp = fabs((l1->p4()+l2->p4()).M() - 91.2);
				if (dztmp < dz) dz=dztmp;
			}
		}
	if(debug0) cout<<"25"<<endl;

	// If taus are present reconstruc tau leptons
//	Handle<pat::METCollection> mets;
//	iEvent.getByLabel(metLabel,mets);

//	double k1=1, k2=1;
//	if (docoll) {
//		if (taus.size() == 2) {
//			double ang = deltaphi(taus[0].p4().phi(),taus[1].p4().phi());
//			if (debug) cout << "Angle between two taus: " << ang << endl;
//		}
//		if (taus.size() == 1 && deltaphi(taus[0].p4().phi(),mets->at(0).phi()) < 3.1415/2) {
//			k1 = 1+( mets->at(0).px() * taus[0].px() + mets->at(0).py() * taus[0].py() ) / ( taus[0].pt()*taus[0].pt() );
//		} else if (taus.size() == 2) {
//			double s1,s2;
//			s1 = taus[0].px()/taus[0].py();
//			s2 = taus[1].px()/taus[1].py();
//			k1 = 1 + ( mets->at(0).px()*s1 - mets->at(0).py()*s1*s2 ) / (taus[0].px() * (s1-s2) );
//			k2 = 1 + ( mets->at(0).py()*s1 - mets->at(0).px() ) / (taus[1].py() * (s1-s2) );
//		}
//	}
	if(debug0) cout<<"26"<<endl;

//	unsigned int np1=0, np2=1, nm1=0, nm2=1, nl1=0, nl2=0, nl3=0, nl4=0, nl[4]={0,};
	unsigned int np1=0, np2=1, nm1=0, nm2=1, nl[4]={0,}, NL[4]={0,};
	double mR=1e10;
	double pt1stP=0, pt2ndP=0, pt1stM=0, pt2ndM=0;
	reco::CompositeRefCandidate lep1[4], lep[4];
	reco::CandidateCollection::const_iterator it, it1, it2, it3, it4;
	std::vector<double> ptV;//, etaV, isoV, pidV;
	std::vector<double> isoV;
	std::vector<double> pidV;
	std::vector<double> etaV;
	if(pLeps->size()>=2&&mLeps->size()>=2)
	{
		for(it1=pLeps->begin(),np1=0; it1!=pLeps->end(); it1++,np1++) for(it2=(it1+1),np2=(np1+1); it2!=pLeps->end(); it2++,np2++)
		{
//			cout<<"m(H++) "<<(it1->p4()+it2->p4()).M()<<", ";
			for(it3=mLeps->begin(),nm1=0; it3!=mLeps->end(); it3++,nm1++) for(it4=(it3+1),nm2=(nm1+1); it4!=mLeps->end(); it4++,nm2++)
			{
				if(debug0) cout<<"26.01"<<endl;
//				double kk1=1, kk2=1, kk3=1, kk4=1;
				// Sanity check is that k1, k2 > 1 (neutrinos contribute some non-zero energy ...)
//				if (taus.size() == 1 && k1 > 1) {
//					if(Ptau==1) {if(np1==itaup[0]) kk1=k1; if(np2==itaup[0]) kk2=k1;}
//					if(Mtau==1) {if(nm1==itaum[0]) kk3=k1; if(nm2==itaum[0]) kk4=k1;}
//				}
//				if (taus.size() == 2 && k1 > 1 && k2 > 1) {
//					if(Ptau==2) {if(np1==itaup[0]) kk1=k1; if(np2==itaup[1]) kk2=k2;}
//					if(Mtau==2) {if(nm1==itaum[0]) kk3=k1; if(nm2==itaum[1]) kk4=k2;}
//				}
//				tmp.setP4(kk1*p1->p4()+kk2*p2->p4()); // By default 4 vectors are not calculated for composites
				double mHPP=(it1->p4()+it2->p4()).M();
				double mHMM=(it3->p4()+it4->p4()).M();
//				double mHPP=(kk1*it1->p4()+kk2*it2->p4()).M();
//				double mHMM=(kk3*it3->p4()+kk4*it4->p4()).M();
//				if(fabs(mHPP/mHMM-1)<mR) {mR=fabs(mHPP/mHMM-1); nl[0]=np1; nl[1]=np2; nl[2]=nm1; nl[3]=nm2; mp=mHPP; mn=mHMM;}
				if(fabs(mHPP-mHMM)/(mHPP+mHMM)<mR) {
				if(debug0) cout<<"26.02"<<endl;
					mR=fabs(mHPP-mHMM)/(mHPP+mHMM); 
					nl[0]=np1; nl[1]=np2; nl[2]=nm1; nl[3]=nm2; 
					mp1=mHPP; 
					mn1=mHMM;
//					m4 = (kk1*it1->p4()+kk2*it2->p4()+kk3*it3->p4()+kk4*it4->p4()).M();
//					m4 = (it1->p4()+it2->p4()+it3->p4()+it4->p4()).M();
//					if(dilpt < (it1->p4()+it2->p4()).Pt()) dilpt = (it1->p4()+it2->p4()).Pt();
//					if(dilpt < (it3->p4()+it4->p4()).Pt()) dilpt = (it3->p4()+it4->p4()).Pt();
//					double dZ13 = fabs((it1->p4()+it3->p4()).M()-91.2);
//					double dZ24 = fabs((it2->p4()+it4->p4()).M()-91.2);
//					double dZ14 = fabs((it1->p4()+it4->p4()).M()-91.2);
//					double dZ23 = fabs((it2->p4()+it3->p4()).M()-91.2);
//					if( (fabs((it1->p4()+it3->p4()).M()-91.2)+fabs((it2->p4()+it4->p4()).M()-91.2)) < (fabs((it1->p4()+it4->p4()).M()-91.2)+fabs((it2->p4()+it3->p4()).M()-91.2)) )
//					if((dZ13+dZ24) < (dZ14+dZ23)) {
//						if(dZ13 < dZ24) dZ = dZ13;  
//						else dZ = dZ24;
//					}
//					if((dZ14+dZ23) < (dZ13+dZ24)) {
//					else {
//						if(dZ14 < dZ23) dZ = dZ14;
//						else dZ = dZ23;
//					}
				}
				lep1[0].setP4(it1->p4());
				lep1[1].setP4(it2->p4());
				lep1[2].setP4(it3->p4());
				lep1[3].setP4(it4->p4());
//				cout<<"m(H--) "<<(it1->p4()+it2->p4()).M()<<", "<<endl;
//				cout<<np1<<" "<<np2<<" "<<nm1<<" "<<nm2<<endl;
				if(debug0) cout<<"26.03"<<endl;
				if(pt1stP<it1->pt()) {pt2ndP=pt1stP; pt1stP=it1->pt(); NL[1]=NL[0]; NL[0]=np1; lep[1].setP4(lep[0].p4()); lep[0].setP4(it1->p4());}
				else if(pt2ndP<it1->pt()) {pt2ndP=it1->pt(); NL[1]=np1; lep[1].setP4(it1->p4());}
				if(pt1stP<it2->pt()) {pt2ndP=pt1stP; pt1stP=it2->pt(); NL[1]=NL[0]; NL[0]=np2; lep[1].setP4(lep[0].p4()); lep[0].setP4(it2->p4());}
				else if(pt2ndP<it2->pt()) {pt2ndP=it2->pt(); NL[1]=np2; lep[1].setP4(it2->p4());}
				if(pt1stM<it3->pt()) {pt2ndM=pt1stM; pt1stM=it3->pt(); NL[3]=NL[2]; NL[2]=nm1; lep[3].setP4(lep[2].p4()); lep[2].setP4(it3->p4());}
				else if(pt2ndM<it3->pt()) {pt2ndM=it3->pt(); NL[3]=nm1; lep[3].setP4(it3->p4());}
				if(pt1stM<it4->pt()) {pt2ndM=pt1stM; pt1stM=it4->pt(); NL[3]=NL[2]; NL[2]=nm2; lep[3].setP4(lep[2].p4()); lep[2].setP4(it4->p4());}
				else if(pt2ndM<it4->pt()) {pt2ndM=it4->pt(); NL[3]=nm2; lep[3].setP4(it4->p4());}
				if(debug0) cout<<"26.04"<<endl;
			}
		}
		if(debug0) cout<<"26.05"<<endl;
//		TLorentzVector tlep1[4];
		unsigned int nLep=0;
		for(it=pLeps->begin(),nLep=0; it!=pLeps->end(); it++,nLep++)
		{
			if(nLep==nl[0]||nLep==nl[1])
			{
//				ptV.push_back(it->pt());
				pidV.push_back(it->pdgId());
				etaV.push_back(it->eta());
				if((int)nLep<Pel)
				{
//					const pat::Electron *el = & (*elIn)[Pposition[0][nLep]]; 
//					isoV.push_back((el->hcalIso()+el->ecalIso()+el->trackIso())/el->pt());
//					isoV.push_back((el->trackIso())/el->pt());
//					tlep1[Ep1+Mp1+Tp1].SetPxPyPzE(el->px(), el->py(), el->pz(), el->energy());
					Ep1++;
				}
				if((int)nLep>=Pel&&(int)nLep<(Pel+Pmu))
				{
//					const pat::Muon *mu = & (*muIn)[Pposition[1][nLep-Pel]];
//					isoV.push_back(mu->isolationR03().sumPt/mu->pt());
//					tlep1[Ep1+Mp1+Tp1].SetPxPyPzE(mu->px(), mu->py(), mu->pz(), mu->energy());
					Mp1++;
				}
				if((int)nLep>=(Pel+Pmu))
				{
//					isoV.push_back(0);
					Tp1++;
				}
			}
			if(nLep==NL[0]||nLep==NL[1])
			{
				ptV.push_back(it->pt());
				etaV.push_back(it->eta());
				if((int)nLep<Pel)
				{
					const pat::Electron *el = & (*elIn)[Pposition[0][nLep]]; 
//					if(iEvent.isRealData())  isoV.push_back((el->trackIso()+el->ecalIso()+el->hcalIso())/el->pt());
//					if(iEvent.isRealData())  isoV.push_back((el->chargedHadronIso()+el->neutralHadronIso()+el->photonIso())/el->pt());
//					if(!iEvent.isRealData()) isoV.push_back(calcIso_el(*el));
//					isoV.push_back(calcIso_el(*el));
					double scEta = fabs(el->superCluster()->eta());
					double Aeff = 0;
					if(!iEvent.isRealData())
					{
						if(scEta<1.0) Aeff=0.21;
						if(scEta>1.0 && scEta<1.479) Aeff=0.21;
						if(scEta>1.479 && scEta<2.0) Aeff=0.11;
						if(scEta>2.0 && scEta<2.2) Aeff=0.14;
						if(scEta>2.2 && scEta<2.3) Aeff=0.18;
						if(scEta>2.3 && scEta<2.4) Aeff=0.19;
						if(scEta>2.4) Aeff=0.26;
					}
					if(iEvent.isRealData())
					{
						if(scEta<1.0) Aeff=0.13;
						if(scEta>1.0 && scEta<1.479) Aeff=0.14;
						if(scEta>1.479 && scEta<2.0) Aeff=0.07;
						if(scEta>2.0 && scEta<2.2) Aeff=0.09;
						if(scEta>2.2 && scEta<2.3) Aeff=0.11;
						if(scEta>2.3 && scEta<2.4) Aeff=0.11;
						if(scEta>2.4) Aeff=0.14;
					}
					isoV.push_back((el->chargedHadronIso() + TMath::Max(0.,el->neutralHadronIso()+el->photonIso()-rho*Aeff))/el->pt());
					lep_eff_pp = lep_eff_pp*lepton_eff(11,it->pt(),it->eta(),0);
					lep_eff_pp_high = lep_eff_pp_high*lepton_eff(11,it->pt(),it->eta(),1);
					lep_eff_pp_low = lep_eff_pp_low*lepton_eff(11,it->pt(),it->eta(),-1);
					Ep++;
				}
				if((int)nLep>=Pel&&(int)nLep<(Pel+Pmu))
				{
					const pat::Muon *mu = & (*muIn)[Pposition[1][nLep-Pel]];
//					if(iEvent.isRealData())  isoV.push_back((mu->trackIso()+mu->ecalIso()+mu->hcalIso())/mu->pt());
//					if(iEvent.isRealData())  isoV.push_back((mu->chargedHadronIso()+mu->neutralHadronIso()+mu->photonIso())/mu->pt());
//					if(!iEvent.isRealData()) isoV.push_back(calcIso_mu(*mu));
//					isoV.push_back(calcIso_mu(*mu));
//					if(!iEvent.isRealData()) isoV.push_back((mu->chargedHadronIso() + TMath::Max(0.,mu->neutralHadronIso()+mu->photonIso()-0.5*mu->puChargedHadronIso()))/mu->pt());
					isoV.push_back((mu->chargedHadronIso() + TMath::Max(0.,mu->neutralHadronIso()+mu->photonIso()-0.5*mu->puChargedHadronIso()))/mu->pt());
					lep_eff_pp = lep_eff_pp*lepton_eff(13,it->pt(),it->eta(),0);
					lep_eff_pp_high = lep_eff_pp_high*lepton_eff(13,it->pt(),it->eta(),1);
					lep_eff_pp_low = lep_eff_pp_low*lepton_eff(13,it->pt(),it->eta(),-1);
					Mp++;
				}
				if((int)nLep>=(Pel+Pmu))
				{
	      				const pat::Tau *tau = & (*tauIn)[Pposition[2][nLep-Pel-Pmu]];
					isoV.push_back(0);
					Tp++;
				}
			}
		}
		if(debug0) cout<<"26.06"<<endl;
		int nl2[4]={0,};
		for(int i=0;i<2;i++)
		{
			if((int)nl[i]<Pel) nl2[i]=0;
			if((int)nl[i]>=Pel&&(int)nl[i]<(Pel+Pmu)) nl2[i]=1;
			if((int)nl[i]>=(Pel+Pmu)) nl2[i]=2;
			if((int)nl[i+2]<Mel) nl2[i+2]=0;
			if((int)nl[i+2]>=Mel&&(int)nl[i+2]<(Mel+Mmu)) nl2[i+2]=1;
			if((int)nl[i+2]>=(Mel+Mmu)) nl2[i+2]=2;
		}
		if(debug2&&events<1000)
		{
			for(int i=0;i<2;i++) cout<<emt[nl2[i]]<<Pposition[nl2[i]][nl[i]]<<", ";
			for(int i=2;i<4;i++) cout<<emt[nl2[i]]<<Mposition[nl2[i]][nl[i]]<<", ";
			cout<<" ==> best lepton candidates"<<endl<<endl;
		}
		for(it=mLeps->begin(),nLep=0; it!=mLeps->end(); it++,nLep++)
		{
			if(nLep==nl[2]||nLep==nl[3])
			{
				if(debug0) cout<<"26.07"<<endl;
//				ptV.push_back(it->pt());
				pidV.push_back(it->pdgId());
				etaV.push_back(it->eta());
				if((int)nLep<Mel)
				{
//					const pat::Electron *el = & (*elIn)[Mposition[0][nLep]]; 
//					isoV.push_back((el->hcalIso()+el->ecalIso()+el->trackIso())/el->pt());
//					isoV.push_back((el->trackIso())/el->pt());
//					tlep[Ep1+Mp1+Tp1+Em1+Mm1+Tm1].SetPxPyPzE(el->px(), el->py(), el->pz(), el->energy());
					Em1++;
				}
				if((int)nLep>=Mel&&(int)nLep<(Mel+Mmu))
				{
//					const pat::Muon *mu = & (*muIn)[Mposition[1][nLep-Mel]];
//					isoV.push_back(mu->isolationR03().sumPt/mu->pt());
//					tlep[Ep1+Mp1+Tp1+Em1+Mm1+Tm1].SetPxPyPzE(mu->px(), mu->py(), mu->pz(), mu->energy());
					Mm1++;
				}
				if((int)nLep>=(Mel+Mmu))
				{
//					isoV.push_back(0);
					Tm1++;
				}
			}
			if(nLep==NL[2]||nLep==NL[3])
			{
				ptV.push_back(it->pt());
				etaV.push_back(it->eta());
				if((int)nLep<Mel)
				{
					const pat::Electron *el = & (*elIn)[Mposition[0][nLep]]; 
//					if(iEvent.isRealData())  isoV.push_back((el->trackIso()+el->ecalIso()+el->hcalIso())/el->pt());
//					if(iEvent.isRealData())  isoV.push_back((el->chargedHadronIso()+el->neutralHadronIso()+el->photonIso())/el->pt());
//					if(!iEvent.isRealData()) isoV.push_back(calcIso_el(*el));
//					isoV.push_back(calcIso_el(*el));
					double scEta = fabs(el->superCluster()->eta());
					double Aeff = 0;
					if(!iEvent.isRealData())
					{
						if(scEta<1.0) Aeff=0.21;
						if(scEta>1.0 && scEta<1.479) Aeff=0.21;
						if(scEta>1.479 && scEta<2.0) Aeff=0.11;
						if(scEta>2.0 && scEta<2.2) Aeff=0.14;
						if(scEta>2.2 && scEta<2.3) Aeff=0.18;
						if(scEta>2.3 && scEta<2.4) Aeff=0.19;
						if(scEta>2.4) Aeff=0.26;
					}
					if(iEvent.isRealData())
					{
						if(scEta<1.0) Aeff=0.13;
						if(scEta>1.0 && scEta<1.479) Aeff=0.14;
						if(scEta>1.479 && scEta<2.0) Aeff=0.07;
						if(scEta>2.0 && scEta<2.2) Aeff=0.09;
						if(scEta>2.2 && scEta<2.3) Aeff=0.11;
						if(scEta>2.3 && scEta<2.4) Aeff=0.11;
						if(scEta>2.4) Aeff=0.14;
					}
					isoV.push_back((el->chargedHadronIso() + TMath::Max(0.,el->neutralHadronIso()+el->photonIso()-rho*Aeff))/el->pt());
					lep_eff_mm = lep_eff_mm*lepton_eff(11,it->pt(),it->eta(),0);
					lep_eff_mm_high = lep_eff_mm_high*lepton_eff(11,it->pt(),it->eta(),1);
					lep_eff_mm_low = lep_eff_mm_low*lepton_eff(11,it->pt(),it->eta(),-1);
					Em++;
				}
				if((int)nLep>=Mel&&(int)nLep<(Mel+Mmu))
				{
					const pat::Muon *mu = & (*muIn)[Mposition[1][nLep-Mel]];
//					if(iEvent.isRealData())  isoV.push_back((mu->trackIso()+mu->ecalIso()+mu->hcalIso())/mu->pt());
//					if(iEvent.isRealData())  isoV.push_back((mu->chargedHadronIso()+mu->neutralHadronIso()+mu->photonIso())/mu->pt());
//					if(!iEvent.isRealData()) isoV.push_back(calcIso_mu(*mu));
//					isoV.push_back(calcIso_mu(*mu));
//					if(!iEvent.isRealData()) isoV.push_back((mu->chargedHadronIso() + TMath::Max(0.,mu->neutralHadronIso()+mu->photonIso()-0.5*mu->puChargedHadronIso()))/mu->pt());
					isoV.push_back((mu->chargedHadronIso() + TMath::Max(0.,mu->neutralHadronIso()+mu->photonIso()-0.5*mu->puChargedHadronIso()))/mu->pt());
					lep_eff_mm = lep_eff_mm*lepton_eff(13,it->pt(),it->eta(),0);
					lep_eff_mm_high = lep_eff_mm_high*lepton_eff(13,it->pt(),it->eta(),1);
					lep_eff_mm_low = lep_eff_mm_low*lepton_eff(13,it->pt(),it->eta(),-1);
					Mm++;
				}
				if((int)nLep>=(Mel+Mmu))
				{
					const pat::Tau *tau = & (*tauIn)[Mposition[2][nLep-Mel-Mmu]];
					isoV.push_back(0);
					Tm++;
				}
			}
		}
		if(debug0) cout<<"26.08"<<endl;
//		sort(ptV.begin(), ptV.end());
		sort(isoV.begin(), isoV.end());
//		sort(pidV.begin(), pidV.end());
		pt1 = ptV[0];
		pt2 = ptV[1];
		pt3 = ptV[2];
		pt4 = ptV[3];
		iso1 = isoV[0];
		iso2 = isoV[1];
		iso3 = isoV[2];
		iso4 = isoV[3];
		pid1 = pidV[0];
		pid2 = pidV[1];
		pid3 = pidV[2];
		pid4 = pidV[3];
		eta1 = etaV[0];
		eta2 = etaV[1];
		eta3 = etaV[2];
		eta4 = etaV[3];
//		cout<<pt1<<" "<<pt2<<" "<<pt3<<" "<<pt4<<" || "<<iso1<<" "<<iso2<<" "<<iso3<<" "<<iso4<<endl;
//		cout<<Pposition[1][0]<<" "<<Pposition[1][1]<<" "<<Mposition[1][0]<<" "<<Mposition[1][1]<<endl;
		recoch1 = channel[Ep1][Mp1][Tp1][Em1][Mm1][Tm1];
		recoch = channel[Ep][Mp][Tp][Em][Mm][Tm];
		if(recoch>21) cout<<"Event "<<events<<", recoch "<<recoch<<" / "<<Ep<<Mp<<Tp<<Em<<Mm<<Tm<<endl;

		TLorentzVector tlep1[4];
		int nmatch1=0;
		if(debug0) cout<<"26.09"<<endl;
		for(int i=0;i<4;i++)
		{
			tlep1[i].SetPxPyPzE(lep1[i].px(), lep1[i].py(), lep1[i].pz(), lep1[i].energy());
//			if(debug0) cout<<"26.09.01, nmatch1 "<<nmatch1<<endl;
//			if(debug0) cout<<"26.04, tlep["<<i<<"] pt "<<tlep[i].Pt()<<endl;
			if(matmcl(tlep1[i],&(*mcl),nmcl)<dR_matching) nmatch1++;
//			if(debug0) cout<<"26.09.02, nmatch1 "<<nmatch1<<endl;
//        		double dRmcmu=100, dRmin=0;
//		        for(int l=0;l<nmcm;l++)
//			{
//				if(debug0) cout<<"26.07.01, mcmu pt"<<mcm[l].Pt()<<endl;
//				dRmin=sqrt(pow(mcm[l].Phi()-tlep[i].Phi(),2)+pow(mcm[l].Eta()-tlep[i].Eta(),2));
//				if(debug0) cout<<"26.07.02, dRmin "<<dRmin<<endl;
//				if(dRmcmu>dRmin) dRmcmu=dRmin;
//				if(debug0) cout<<"26.07.03"<<endl;
//			}
//			if(debug0) cout<<"26.07.04, dRmcmu "<<dRmcmu<<endl;
//			if(dRmcmu<dR_matching) nmatch1++;
		}
		if(nmatch1==4) genmatch1=1;
		if(debug0) cout<<"26.10"<<endl;

		mp = (lep[0].p4()+lep[1].p4()).M();
		mn = (lep[2].p4()+lep[3].p4()).M();
		m4 = (lep[0].p4()+lep[1].p4()+lep[2].p4()+lep[3].p4()).M();
		if(dilpt < (lep[0].p4()+lep[1].p4()).Pt()) dilpt = (lep[0].p4()+lep[1].p4()).Pt();
		if(dilpt < (lep[2].p4()+lep[3].p4()).Pt()) dilpt = (lep[2].p4()+lep[3].p4()).Pt();
		double mZ13 = (lep[0].p4()+lep[2].p4()).M();
		double mZ24 = (lep[1].p4()+lep[3].p4()).M();
		double mZ14 = (lep[0].p4()+lep[3].p4()).M();
		double mZ23 = (lep[1].p4()+lep[2].p4()).M();
		double dZ13 = fabs(mZ13-91.2);
		double dZ24 = fabs(mZ24-91.2);
		double dZ14 = fabs(mZ14-91.2);
		double dZ23 = fabs(mZ23-91.2);
		if((dZ13+dZ24) < (dZ14+dZ23)) {
			if(dZ13 < dZ24) dZ = dZ13;  
			else dZ = dZ24;
			mZ1 = mZ13;
			mZ2 = mZ24;
		}
		else {
			if(dZ14 < dZ23) dZ = dZ14;
			else dZ = dZ23;
			mZ1 = mZ14;
			mZ2 = mZ23;
		}
		if(debug0) cout<<"26.11"<<endl;
		TLorentzVector tlep[4];
		int nmatch=0;
		for(int i=0;i<4;i++) tlep[i].SetPxPyPzE(lep[i].px(), lep[i].py(), lep[i].pz(), lep[i].energy());
		for(int i=0;i<4;i++) if(matmcl(tlep[i],&(*mcl),nmcl)<dR_matching) nmatch++;
//		cout<<"Event "<<events<<", "; for(int i=0;i<4;i++) cout<<matmcl(tlep[i],&(*mcl),nmcl)<<", "; cout<<endl;
		if(nmatch==4) genmatch=1;

//		reco::Candidate::LorentzVector  L0(lep[0].px()*DeltaPt,lep[0].py()*DeltaPt,0,0);
//		reco::Candidate::LorentzVector  L1(lep[1].px()*DeltaPt,lep[1].py()*DeltaPt,0,0);
//		reco::Candidate::LorentzVector  L2(lep[2].px()*DeltaPt,lep[2].py()*DeltaPt,0,0);
//		reco::Candidate::LorentzVector  L3(lep[3].px()*DeltaPt,lep[3].py()*DeltaPt,0,0);
//		reco::Candidate::LorentzVector  mp_vp=lep[0].p4()+lep[1].p4()+L0+L1;
//		reco::Candidate::LorentzVector  mp_vm=lep[0].p4()+lep[1].p4()-L0-L1;
//		reco::Candidate::LorentzVector  mn_vp=lep[2].p4()+lep[3].p4()+L2+L3;
//		reco::Candidate::LorentzVector  mn_vm=lep[2].p4()+lep[3].p4()-L2-L3;
//		mpp = mp_vp.mass();
//		mpm = mp_vm.mass();
//		mnp = mn_vp.mass();
//		mnm = mn_vm.mass();

//		reco::Candidate::LorentzVector  m4_vp=lep[0].p4()+lep[1].p4()+lep[2].p4()+lep[3].p4()+L0+L1+L2+L3;
//		reco::Candidate::LorentzVector  m4_vm=lep[0].p4()+lep[1].p4()+lep[2].p4()+lep[3].p4()-L0-L1-L2-L3;
//		m4p = m4_vp.mass();
//		m4m = m4_vm.mass();

		TLorentzVector L0p, L0m, L1p, L1m, L2p, L2m, L3p, L3m;
		L0p.SetPtEtaPhiM(lep[0].pt()*(1+DeltaPt),lep[0].eta(),lep[0].phi(),0);
		L0m.SetPtEtaPhiM(lep[0].pt()*(1-DeltaPt),lep[0].eta(),lep[0].phi(),0);
		L1p.SetPtEtaPhiM(lep[1].pt()*(1+DeltaPt),lep[1].eta(),lep[1].phi(),0);
		L1m.SetPtEtaPhiM(lep[1].pt()*(1-DeltaPt),lep[1].eta(),lep[1].phi(),0);
		L2p.SetPtEtaPhiM(lep[2].pt()*(1+DeltaPt),lep[2].eta(),lep[2].phi(),0);
		L2m.SetPtEtaPhiM(lep[2].pt()*(1-DeltaPt),lep[2].eta(),lep[2].phi(),0);
		L3p.SetPtEtaPhiM(lep[3].pt()*(1+DeltaPt),lep[3].eta(),lep[3].phi(),0);
		L3m.SetPtEtaPhiM(lep[3].pt()*(1-DeltaPt),lep[3].eta(),lep[3].phi(),0);
		mpp = (L0p+L1p).M();
		mpm = (L0m+L1m).M();
		mnp = (L2p+L3p).M();
		mnm = (L2m+L3m).M();
		m4p = (L0p+L1p+L2p+L3p).M();
		m4m = (L0m+L1m+L2m+L3m).M();
//cout<<"ev "<<events<<", mp "<<mp<<", mpp "<<mpp<<", mpp2 "<<(L1p+L2p).M()<<", mpm "<<mpm<<", mpm2 "<<(L1m+L2m).M()<<endl;
	}
	if(debug0) cout<<"27"<<endl;
	if((pLeps->size()==2&&mLeps->size()==1)||(pLeps->size()==1&&mLeps->size()==2))
	{
		if(selLeps->size()!=3) cout<<"check for size of pLeps, mLeps and selLeps"<<endl;
		if(selLeps->size()==3) m4 = ((*selLeps)[0].p4()+(*selLeps)[1].p4()+(*selLeps)[2].p4()).M();
//		cout<<"m3 : "<<m4<<endl;
	}
//	if(pLeps->size()==2&&mLeps->size()==1) cout<<"m(H++) : "<<((*pLeps)[0].p4()+(*pLeps)[1].p4()).M()<<endl; //mp
//	if(pLeps->size()==1&&mLeps->size()==2) cout<<"m(H--) : "<<((*mLeps)[0].p4()+(*mLeps)[1].p4()).M()<<endl; //mn
	if(pLeps->size()==2&&mLeps->size()==1) {
		pt1 = (*pLeps)[0].pt();
		pt2 = (*pLeps)[1].pt();
		pt3 = (*mLeps)[0].pt();
		mp = ((*pLeps)[0].p4()+(*pLeps)[1].p4()).M();
		mllss = mp;
		TLorentzVector LVlep;
//		LVlep.SetPtEtaPhiM((*mLeps)[0].pt(),(*mLeps)[0].eta(),(*mLeps)[0].phi(),0);
		double mz1 = ((*pLeps)[0].p4()+(*mLeps)[0].p4()).M();
		double mz2 = ((*pLeps)[1].p4()+(*mLeps)[0].p4()).M();
		double dZ00 = fabs(mz1-91.2);
		double dZ10 = fabs(mz2-91.2);
		if(dZ00 < dZ10) {dZ = dZ00; mllpm = mz1; LVlep.SetPtEtaPhiM((*pLeps)[1].pt(),(*pLeps)[1].eta(),(*pLeps)[1].phi(),0);}
		else {dZ = dZ10; mllpm = mz2; LVlep.SetPtEtaPhiM((*pLeps)[0].pt(),(*pLeps)[0].eta(),(*pLeps)[0].phi(),0);}
		double pt00 = ((*pLeps)[0].p4()+(*mLeps)[0].p4()).Pt();
		double pt10 = ((*pLeps)[1].p4()+(*mLeps)[0].p4()).Pt();
		if(pt00 < pt10) dilpt = pt10;
		else dilpt = pt00;
		m4 = ((*pLeps)[0].p4()+(*pLeps)[1].p4()+(*mLeps)[0].p4()).M();
//		eP=nelP, eM=nelM;
//		mP=nmuP, mM=nmuM;
//		tP=0, tM=0;
		Ep=nelP, Em=nelM;
		Mp=nmuP, Mm=nmuM;
		Tp=0, Tm=0;

		mW = sqrt( 2*LVmet.Pt()*LVlep.Pt()*(1-TMath::Cos(LVmet.DeltaPhi(LVlep))) );
	}
	if(pLeps->size()==1&&mLeps->size()==2) {
		pt1 = (*pLeps)[0].pt();
		pt2 = (*mLeps)[0].pt();
		pt3 = (*mLeps)[1].pt();
		mn=((*mLeps)[0].p4()+(*mLeps)[1].p4()).M();
		mllss = mn;
		TLorentzVector LVlep;
//		LVlep.SetPtEtaPhiM((*mLeps)[0].pt(),(*mLeps)[0].eta(),(*mLeps)[0].phi(),0);
		double mz1 = ((*pLeps)[0].p4()+(*mLeps)[0].p4()).M();
		double mz2 = ((*pLeps)[0].p4()+(*mLeps)[1].p4()).M();
		double dZ00 = fabs(mz1-91.2);
		double dZ01 = fabs(mz2-91.2);
		if(dZ00 < dZ01) {dZ = dZ00; mllpm = mz1; LVlep.SetPtEtaPhiM((*mLeps)[1].pt(),(*mLeps)[1].eta(),(*mLeps)[1].phi(),0);}
		else {dZ = dZ01; mllpm = mz2; LVlep.SetPtEtaPhiM((*mLeps)[0].pt(),(*mLeps)[0].eta(),(*mLeps)[0].phi(),0);}
		double pt00 = ((*pLeps)[0].p4()+(*mLeps)[0].p4()).Pt();
		double pt01 = ((*pLeps)[0].p4()+(*mLeps)[1].p4()).Pt();
		if(pt00 < pt01) dilpt = pt01;
		else dilpt = pt00;
		m4 = ((*pLeps)[0].p4()+(*mLeps)[0].p4()+(*mLeps)[1].p4()).M();
//		eP=nelP, eM=nelM;
//		mP=nmuP, mM=nmuM;
//		tP=0, tM=0;
		Ep=nelP, Em=nelM;
		Mp=nmuP, Mm=nmuM;
		Tp=0, Tm=0;

		mW = sqrt( 2*LVmet.Pt()*LVlep.Pt()*(1-TMath::Cos(LVmet.DeltaPhi(LVlep))) );
	}
	for(reco::CandidateCollection::const_iterator it=selLeps->begin(); it!=selLeps->end(); it++)
	{
		if (debug) cout<<"pdgId "<<it->pdgId()<<", pt "<<it->pt()<<", charge "<<it->charge()<<", eta "<<it->eta()<<endl;
//		const reco::Candidate *c = it;
//		const Electron *el = dynamic_cast<const Electron*>(c);
//		cout<<el->hcalIso()<<endl;
//		cout<<dynamic_cast<const Electron*>(c)->hcalIso()<<endl;
//		if(abs(it->pdgId())==11) cout<<"hcalIso "<<it->hcalIso()<<", ecalIs "<<it->ecalIso()<<", trackIso "<<it->trackIso()<<endl;
//		if(abs(it->pdgId())==13) cout<<"sumPt "<<it->isolationR03().sumPt<<endl;
		if(it->charge()== 1) nlp++;
		if(it->charge()==-1) nlm++;
//		ptV.push_back(it->pt());
	}
	if((unsigned int)nlp!=pLeps->size()) cout<<"ev "<<events<<" : check for plus selected leptons"<<endl;
	if((unsigned int)nlm!=mLeps->size()) cout<<"ev "<<events<<" : check for minus selected leptons"<<endl;

	if(debug0) cout<<"28"<<endl;
	if (debug) cout<<"ev "<<events<<", ";
//}
        Handle<TriggerResults> trig;
        iEvent.getByLabel(InputTag("TriggerResults","","HLT"),trig);
//      iEvent.getByLabel("TriggerResults",trig);
        const TriggerResults* trigResult = trig.product();
//      if(events_hlt0==1) cout<<"hlttrigger size : "<<trigResult->size()<<endl;
        TriggerNames triggerNames_ = iEvent.triggerNames(*trig);
        TriggerResultsByName triggerNames2_(&(*trig),&triggerNames_);
        bool HLT_pass = false;
        if(trigResult->wasrun() && trigResult->accept())
        {               
                const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigResult);
                TriggerResultsByName trigNames2(&(*trigResult),&trigNames);
                for(unsigned int i=0;i<HLT_path.size();i++)
                {
//			if(HLT_path_find[i]) if(trigNames2.accept(HLT_path[i]) && hltConfigProvider_.prescaleValue(iEvent,iSetup,HLT_path[i])==1) HLT_pass = true;
			if(HLT_path_find[i]) if(trigNames2.accept(HLT_path[i])) HLT_pass = true;
//                      if(trigNames2.find(HLT_path[i]) != std::string::npos) if(hltConfigProvider_.prescaleValue(iEvent,iSetup,trigNames2)==1) HLT_pass = true;
//			if(trigname.find(HLT_path[i]) != std::string::npos)
//			for(unsigned int j=0;j<trigResult->size();j++)
//			{
//				string triggerName = trigNames.triggerName(j);
//				if(triggerName.find(HLT_path[i]) != std::string::npos) if(hltConfigProvider_.prescaleValue(iEvent,iSetup,triggerName)==1) HLT_pass = true;
//				cout<<triggerName<<", "<<HLT_path[i]<<", "<<hltConfigProvider_.prescaleValue(iEvent,iSetup,triggerName)<<endl;
//			}
                }
	}
	if(HLT_pass) HLT = 1;

	tree->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
   Handle<ExampleData> pIn;
   iEvent.getByLabel("example",pIn);
#endif
   
#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
*/
}


// ------------ method called once each job just before starting event loop  ------------
void 
DoublyChargedHiggsPAT::beginJob()
//DoublyChargedHiggsPAT::beginJob(const edm::EventSetup&)
{
	cout<<"Event start"<<endl;
}

void
//DoublyChargedHiggsPAT::beginRun(edm::Run &run, edm::EventSetup const &iSetup)
DoublyChargedHiggsPAT::beginRun(const edm::Run &run, const edm::EventSetup &iSetup)
{
//      bool isConfigChanged = false;
        bool isConfigChanged = true;
 
        // isValidHltConfig_ used to short-circuit analyze() in case of problems
        isValidHltConfig_ = hltConfigProvider_.init( run, iSetup, "HLT", isConfigChanged );
 
//      isValidHltConfig_ = hltConfigProvider_.init( run, iSetup, "TriggerResults", isConfigChanged );
//      isValidHltConfig_ = hltConfigProvider_.init( run, iSetup, "TEST", isConfigChanged );
        //std::cout << "hlt config trigger is valid??" << isValidHltConfig_ << std::endl; 
        if(isValidHltConfig_) cout<<"initialisation has succeeded! - run "<<run.id().run()<<endl;
 
//        for(int i=0;i<ntrig;i++) HLT_path_find[i] = 0;
        for(unsigned int i=0;i<100;i++)
        {
                HLT_path_find[i] = 0;
        }
 
        const std::vector<std::string>& triggerNames3 = hltConfigProvider_.triggerNames();
        cout<<"Used paths : "<<endl;
        for (size_t ts = 0; ts< triggerNames3.size() ; ts++)
        {       
                string trigname = triggerNames3[ts];
//		cout<<trigname<<endl;
                for(unsigned int i=0;i<HLT_path.size();i++)
                        if (trigname.find(HLT_path[i]) != std::string::npos)
                        {       
				cout<<trigname<<endl;
                                HLT_path_find[i] = 1; 
                                HLT_path[i] = trigname;
                        }
        }
}

void
//DoublyChargedHiggsPAT::endRun(edm::Run &, edm::EventSetup const &)
DoublyChargedHiggsPAT::endRun(const edm::Run &, const edm::EventSetup &)
{
}
// ------------ method called once each job just after ending the event loop  ------------
void 
DoublyChargedHiggsPAT::endJob() {
//	HPP_TTree->Write();
        HPP_TTree->cd();
	hmupt->Write();
	hmueta->Write();
	helpt->Write();
	heleta->Write();
	hprobe_best_pt->Write();
	hprobe_best_eta->Write();
	hprobe_best_m->Write();
	hprobe_match_pt->Write();
	hprobe_match_eta->Write();
	hprobe_match_m->Write();
	hprobe_matmc_pt->Write();
	hprobe_matmc_eta->Write();
	hprobe_matmc_m->Write();
        tree->Write();
        HPP_TTree->Close();
	cout<<"end job"<<endl<<"Event "<<events<<endl;
}

//define this as a plug-in
DEFINE_FWK_MODULE(DoublyChargedHiggsPAT);
