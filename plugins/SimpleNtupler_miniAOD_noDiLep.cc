#include <boost/foreach.hpp>
#include "TTree.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/MuonReco/interface/MuonCocktails.h"
#include "DataFormats/PatCandidates/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/AsymFunctions.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/TrackUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/DileptonUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/GeneralUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/HardInteraction.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/PATUtilities.h"
#include "SUSYBSMAnalysis/Zprime2muAnalysis/src/ToConcrete.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TLorentzVector.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

class SimpleNtupler_miniAOD_noDiLep : public edm::EDAnalyzer {
public:
  explicit SimpleNtupler_miniAOD_noDiLep(const edm::ParameterSet&);
  ~SimpleNtupler_miniAOD_noDiLep() { delete hardInteraction; }
  void analyze(const edm::Event&, const edm::EventSetup&);
  void endJob();
  TString replace_all(const TString& a, const TString& b, const TString& c);
private:
  struct tree_t {

    unsigned run;
    unsigned lumi;
    unsigned event;
    float genWeight;
    int nQuark_initialState;
    int nGluon_initialState;
    int nGisr;
    float beamspot_x;
    float beamspot_x_err;
    float beamspot_y;
    float beamspot_y_err;
    float beamspot_z;
    float beamspot_z_err;
    int nvertices;
    int dil_chosen;
    float dil_mass;
    float dil_pt;
    float dil_rap;
    float dil_eta;
    float dil_phi;
    float dil_dR;
    float dil_dPhi;
    float dil_lep_pt[2];
    float cos_angle;
    float vertex_chi2;
    float cos_cs;
    float chi_dilepton;
    float phi_cs;
    float vertex_m;
    float vertex_m_err;
    float vertex_x;
    float vertex_x_err;
    float vertex_y;
    float vertex_y_err;
    float vertex_z;
    float vertex_z_err;
    int lep_id[2];
    float lep_p[2];
    float lep_pt[2];
    float lep_pt_err[2];
    float lep_px[2];
    float lep_py[2];
    float lep_pz[2];
    float lep_E[2];
    float lep_eta[2];
    float lep_phi[2];
    float lep_dxy[2];
    float lep_dz[2];
    float lep_qOverPt[2];
    float lep_tk_p[2];
    float lep_tk_pt[2];
    float lep_tk_pt_err[2];
    float lep_tk_px[2];
    float lep_tk_py[2];
    float lep_tk_pz[2];
    float lep_tk_eta[2];
    float lep_tk_phi[2];
    float lep_tk_dz[2];
    float lep_tk_vz[2];
    float lep_tk_chi2[2];
    float lep_tk_ndf[2];
    float lep_tk_qOverPt[2];
    float lep_glb_p[2];
    float lep_glb_pt[2];
    float lep_glb_pt_err[2];
    float lep_glb_px[2];
    float lep_glb_py[2];
    float lep_glb_pz[2];
    float lep_glb_eta[2];
    float lep_glb_phi[2];
    float lep_glb_chi2[2];
    float lep_glb_ndf[2];
    float lep_glb_qOverPt[2];
    float lep_tpfms_p[2];
    float lep_tpfms_pt[2];
    float lep_tpfms_pt_err[2];
    float lep_tpfms_px[2];
    float lep_tpfms_py[2];
    float lep_tpfms_pz[2];
    float lep_tpfms_eta[2];
    float lep_tpfms_phi[2];
    float lep_tpfms_chi2[2];
    float lep_tpfms_ndf[2];
    float lep_tpfms_qOverPt[2];
    float lep_picky_p[2];
    float lep_picky_pt[2];
    float lep_picky_pt_err[2];
    float lep_picky_px[2];
    float lep_picky_py[2];
    float lep_picky_pz[2];
    float lep_picky_eta[2];
    float lep_picky_phi[2];
    float lep_picky_chi2[2];
    float lep_picky_ndf[2];
    float lep_picky_qOverPt[2];
    float lep_cocktail_p[2];
    float lep_cocktail_pt[2];
    float lep_cocktail_pt_err[2];
    float lep_cocktail_px[2];
    float lep_cocktail_py[2];
    float lep_cocktail_pz[2];
    float lep_cocktail_eta[2];
    float lep_cocktail_phi[2];
    float lep_cocktail_chi2[2];
    float lep_cocktail_ndf[2];
    float lep_cocktail_qOverPt[2];
    short lep_cocktail_choice[2];
    float lep_tuneP_p[2];
    float lep_tuneP_pt[2];
    float lep_tuneP_pt_err[2];
    float lep_tuneP_px[2];
    float lep_tuneP_py[2];
    float lep_tuneP_pz[2];
    float lep_tuneP_eta[2];
    float lep_tuneP_phi[2];
    float lep_tuneP_dz[2];
    float lep_tuneP_vz[2];
    float lep_tuneP_chi2[2];
    float lep_tuneP_ndf[2];
    float lep_tuneP_qOverPt[2];
    float lep_triggerMatchPt[2];
    float lep_triggerMatchEta[2];
    float lep_chi2dof[2];
    float lep_dB[2];
    float lep_sumPt[2];
    float lep_emEt[2];
    float lep_hadEt[2];
    float lep_hoEt[2];
    float lep_pfIso[2];
    float lep_pfIsoDB[2];
    int lep_timeNdof[2];
    float lep_timeInOut[2];
    float lep_timeOutIn[2];
    float lep_timeInOutErr[2];
    float lep_timeOutInErr[2];
    int lep_heep_id[2];
    float lep_min_muon_dR[2];
    short lep_tk_numberOfValidTrackerHits[2];
    short lep_tk_numberOfValidTrackerLayers[2];
    short lep_tk_numberOfValidPixelHits[2];
    short lep_glb_numberOfValidTrackerHits[2]; 
    short lep_glb_numberOfValidTrackerLayers[2]; 
    short lep_glb_numberOfValidPixelHits[2];
    short lep_glb_numberOfValidMuonHits[2];
    short lep_glb_numberOfValidMuonDTHits[2];
    short lep_glb_numberOfValidMuonCSCHits[2];
    short lep_glb_numberOfValidMuonRPCHits[2];
    short lep_glb_muonStationsWithValidHits[2];
    short lep_glb_dtStationsWithValidHits[2];
    short lep_glb_cscStationsWithValidHits[2];
    short lep_glb_rpcStationsWithValidHits[2];
    short lep_glb_innermostMuonStationWithValidHits[2];
    short lep_glb_outermostMuonStationWithValidHits[2];
    short lep_numberOfMatches[2];
    short lep_numberOfMatchedStations[2];
    short lep_numberOfMatchedRPCLayers[2];
    unsigned int lep_stationMask[2];
    int lep_numberOfChambers[2];
    int lep_numberOfChambersNoRPC[2];
    unsigned int lep_stationGapMaskDistance[2];
    unsigned int lep_stationGapMaskPull[2];
    bool lep_isGlobalMuon[2];
    bool lep_isTrackerMuon[2];
    bool GoodDataRan;
    bool GoodVtx;
    bool METFilter;
    float gen_res_mass;
    float gen_res_pt;
    float gen_res_rap;
    float gen_res_eta;
    float gen_res_phi;
    float gen_dil_mass;
    float gen_dil_pt;
    float gen_dil_rap;
    float gen_dil_eta;
    float gen_dil_phi;
    float gen_dil_dR;
    float gen_dil_dPhi;
    float gen_lep_p[2];
    float gen_lep_pt[2];
    float gen_lep_px[2];
    float gen_lep_py[2];
    float gen_lep_pz[2];
    float gen_lep_E[2];
    float gen_lep_eta[2];
    float gen_lep_phi[2];
    float gen_lep_qOverPt[2];
    float gen_lep_noib_p[2];
    float gen_lep_noib_pt[2];
    float gen_lep_noib_px[2];
    float gen_lep_noib_py[2];
    float gen_lep_noib_pz[2];
    float gen_lep_noib_E[2];
    float gen_lep_noib_eta[2];
    float gen_lep_noib_phi[2];
    float gen_lep_noib_qOverPt[2];
    float met_pt;
    float met_phi;
    int nJets;
    float jet_pt[10];
    float jet_eta[10];
    float jet_phi[10];
    float jet_btag[10];
    //Added jet information
    float jet_E[10];

  };

  tree_t t;
  TTree* tree;

  const static int nCuts = 10;
  int cutflow_study_dilep[nCuts];

  TH1D* NPV_no_cuts;
  TH1D* NPV_dilepton_loop;
  TH1D* NPV_cuts;
  TH1F* TH1F_cutlfow;

  TH1F* TH1F_muonPt_all;
  TH1F* TH1F_muonPt_trigg;
  TH1F* TH1F_muonPt_one;
  TH1F* TH1F_muonPt_two;
  TH1F* TH1F_muonPt_oppSign;
  TH1F* TH1F_muonPt_deltaR;
  TH1F* TH1F_muonPt_mass;

  TH1F* TH1F_muonEta_all;
  TH1F* TH1F_muonEta_trigg;
  TH1F* TH1F_muonEta_one;
  TH1F* TH1F_muonEta_two;
  TH1F* TH1F_muonEta_oppSign;
  TH1F* TH1F_muonEta_deltaR;
  TH1F* TH1F_muonEta_mass;

  const edm::InputTag mu_src;
  const edm::InputTag beamspot_src;
  const edm::InputTag met_src;
  const edm::InputTag jet_src;
  const edm::InputTag vertices_src;
  const bool fill_gen_info;
  const edm::InputTag TriggerResults_src;
  const edm::InputTag genEventInfo_;
  std::vector<edm::InputTag> filterTags;
  HardInteraction* hardInteraction;
  //const edm::InputTag TriggerResults_src; 
};

TString SimpleNtupler_miniAOD_noDiLep::replace_all(const TString& a, const TString& b, const TString& c) {
  TString ret = a;
  ret.ReplaceAll(b, c);
  return ret;
}

SimpleNtupler_miniAOD_noDiLep::SimpleNtupler_miniAOD_noDiLep(const edm::ParameterSet& cfg)
  : mu_src(cfg.getParameter<edm::InputTag>("mu_src")),
    beamspot_src(cfg.getParameter<edm::InputTag>("beamspot_src")),
    met_src(cfg.getParameter<edm::InputTag>("met_src")),
    jet_src(cfg.getParameter<edm::InputTag>("jet_src")),
    vertices_src(cfg.getParameter<edm::InputTag>("vertices_src")),
    fill_gen_info(cfg.existsAs<edm::ParameterSet>("hardInteraction")),
    TriggerResults_src(cfg.getParameter<edm::InputTag>("TriggerResults_src")),
    genEventInfo_(cfg.getUntrackedParameter<edm::InputTag>("genEventInfo")),
    //filterTags(cfg.getParameter<std::vector<edm::InputTag> > ("metFilter")),  
    hardInteraction(fill_gen_info ? new HardInteraction(cfg.getParameter<edm::ParameterSet>("hardInteraction")) : 0)
{

  consumes<std::vector<pat::Muon>>(mu_src);
  consumes<std::vector<pat::MET>>(met_src);
  consumes<std::vector<pat::Jet>>(jet_src);
  consumes<reco::BeamSpot>(beamspot_src);
  consumes<reco::VertexCollection>(vertices_src);
  consumes<edm::TriggerResults>(TriggerResults_src);
  consumes<GenEventInfoProduct>(genEventInfo_);
  if (fill_gen_info) consumes<std::vector<reco::GenParticle>>(hardInteraction->src);
 
  edm::Service<TFileService> fs;
  tree = fs->make<TTree>("t", "");
  tree->Branch("run", &t.run, "run/i");
  tree->Branch("lumi", &t.lumi, "lumi/i");
  tree->Branch("event", &t.event, "event/i");
  tree->Branch("beamspot_x", &t.beamspot_x, "beamspot_x/F");
  tree->Branch("beamspot_x_err", &t.beamspot_x_err, "beamspot_x_err/F");
  tree->Branch("beamspot_y", &t.beamspot_y, "beamspot_y/F");
  tree->Branch("beamspot_y_err", &t.beamspot_y_err, "beamspot_y_err/F");
  tree->Branch("beamspot_z", &t.beamspot_z, "beamspot_z/F");
  tree->Branch("beamspot_z_err", &t.beamspot_z_err, "beamspot_z_err/F");
  tree->Branch("nvertices", &t.nvertices, "nvertices/I");
  tree->Branch("dil_chosen", &t.dil_chosen, "dil_chosen/I");
  tree->Branch("dil_mass", &t.dil_mass, "dil_mass/F");
  tree->Branch("dil_pt", &t.dil_pt, "dil_pt/F");
  tree->Branch("dil_rap", &t.dil_rap, "dil_rap/F");
  tree->Branch("dil_eta", &t.dil_eta, "dil_eta/F");
  tree->Branch("dil_phi", &t.dil_phi, "dil_phi/F");
  tree->Branch("dil_dR", &t.dil_dR, "dil_dR/F");
  tree->Branch("dil_dPhi", &t.dil_dPhi, "dil_dPhi/F");
  tree->Branch("dil_lep_pt", t.dil_lep_pt, "dil_lep_pt[2]/F");
  tree->Branch("cos_angle", &t.cos_angle, "cos_angle/F");
  tree->Branch("vertex_chi2", &t.vertex_chi2, "vertex_chi2/F");
  tree->Branch("cos_cs", &t.cos_cs, "cos_cs/F");
  tree->Branch("chi_dilepton", &t.chi_dilepton, "chi_dilepton/F");
  tree->Branch("phi_cs", &t.phi_cs, "phi_cs/F");
  tree->Branch("vertex_m", &t.vertex_m, "vertex_m/F");
  tree->Branch("vertex_m_err", &t.vertex_m_err, "vertex_m_err/F");
  tree->Branch("vertex_x", &t.vertex_x, "vertex_x/F");
  tree->Branch("vertex_x_err", &t.vertex_x_err, "vertex_x_err/F");
  tree->Branch("vertex_y", &t.vertex_y, "vertex_y/F");
  tree->Branch("vertex_y_err", &t.vertex_y_err, "vertex_y_err/F");
  tree->Branch("vertex_z", &t.vertex_z, "vertex_z/F");
  tree->Branch("vertex_z_err", &t.vertex_z_err, "vertex_z_err/F");
  tree->Branch("lep_id", t.lep_id, "lep_id[2]/I");
  tree->Branch("lep_heep_id", t.lep_heep_id, "lep_heep_id[2]/I");
  tree->Branch("lep_p", t.lep_p, "lep_p[2]/F");
  tree->Branch("lep_pt", t.lep_pt, "lep_pt[2]/F");
  tree->Branch("lep_pt_err", t.lep_pt_err, "lep_pt_err[2]/F");
  tree->Branch("lep_px", t.lep_px, "lep_px[2]/F");
  tree->Branch("lep_py", t.lep_py, "lep_py[2]/F");
  tree->Branch("lep_pz", t.lep_pz, "lep_pz[2]/F");
  tree->Branch("lep_E", t.lep_E, "lep_E[2]/F");
  tree->Branch("lep_eta", t.lep_eta, "lep_eta[2]/F");
  tree->Branch("lep_phi", t.lep_phi, "lep_phi[2]/F");
  tree->Branch("lep_dxy", t.lep_dxy, "lep_dxy[2]/F");
  tree->Branch("lep_dz", t.lep_dz, "lep_dz[2]/F");
  tree->Branch("lep_qOverPt", t.lep_qOverPt, "lep_qOverPt[2]/F");
  tree->Branch("lep_tk_p", t.lep_tk_p, "lep_tk_p[2]/F");
  tree->Branch("lep_tk_pt", t.lep_tk_pt, "lep_tk_pt[2]/F");
  tree->Branch("lep_tk_pt_err", t.lep_tk_pt_err, "lep_tk_pt_err[2]/F");
  tree->Branch("lep_tk_px", t.lep_tk_px, "lep_tk_px[2]/F");
  tree->Branch("lep_tk_py", t.lep_tk_py, "lep_tk_py[2]/F");
  tree->Branch("lep_tk_pz", t.lep_tk_pz, "lep_tk_pz[2]/F");
  tree->Branch("lep_tk_eta", t.lep_tk_eta, "lep_tk_eta[2]/F");
  tree->Branch("lep_tk_phi", t.lep_tk_phi, "lep_tk_phi[2]/F");
  tree->Branch("lep_tk_dz", t.lep_tk_dz, "lep_tk_dz[2]/F");
  tree->Branch("lep_tk_vz", t.lep_tk_vz, "lep_tk_vz[2]/F");
  tree->Branch("lep_tk_chi2", t.lep_tk_chi2, "lep_tk_chi2[2]/F");
  tree->Branch("lep_tk_ndf", t.lep_tk_ndf, "lep_tk_ndf[2]/F");
  tree->Branch("lep_tk_qOverPt", t.lep_tk_qOverPt, "lep_tk_qOverPt[2]/F");
  tree->Branch("lep_glb_p", t.lep_glb_p, "lep_glb_p[2]/F");
  tree->Branch("lep_glb_pt", t.lep_glb_pt, "lep_glb_pt[2]/F");
  tree->Branch("lep_glb_pt_err", t.lep_glb_pt_err, "lep_glb_pt_err[2]/F");
  tree->Branch("lep_glb_px", t.lep_glb_px, "lep_glb_px[2]/F");
  tree->Branch("lep_glb_py", t.lep_glb_py, "lep_glb_py[2]/F");
  tree->Branch("lep_glb_pz", t.lep_glb_pz, "lep_glb_pz[2]/F");
  tree->Branch("lep_glb_eta", t.lep_glb_eta, "lep_glb_eta[2]/F");
  tree->Branch("lep_glb_phi", t.lep_glb_phi, "lep_glb_phi[2]/F");
  tree->Branch("lep_glb_chi2", t.lep_glb_chi2, "lep_glb_chi2[2]/F");
  tree->Branch("lep_glb_ndf", t.lep_glb_ndf, "lep_glb_ndf[2]/F");
  tree->Branch("lep_glb_qOverPt", t.lep_glb_qOverPt, "lep_glb_qOverPt[2]/F");
  tree->Branch("lep_tpfms_p", t.lep_tpfms_p, "lep_tpfms_p[2]/F");
  tree->Branch("lep_tpfms_pt", t.lep_tpfms_pt, "lep_tpfms_pt[2]/F");
  tree->Branch("lep_tpfms_pt_err", t.lep_tpfms_pt_err, "lep_tpfms_pt_err[2]/F");
  tree->Branch("lep_tpfms_px", t.lep_tpfms_px, "lep_tpfms_px[2]/F");
  tree->Branch("lep_tpfms_py", t.lep_tpfms_py, "lep_tpfms_py[2]/F");
  tree->Branch("lep_tpfms_pz", t.lep_tpfms_pz, "lep_tpfms_pz[2]/F");
  tree->Branch("lep_tpfms_eta", t.lep_tpfms_eta, "lep_tpfms_eta[2]/F");
  tree->Branch("lep_tpfms_phi", t.lep_tpfms_phi, "lep_tpfms_phi[2]/F");
  tree->Branch("lep_tpfms_chi2", t.lep_tpfms_chi2, "lep_tpfms_chi2[2]/F");
  tree->Branch("lep_tpfms_ndf", t.lep_tpfms_ndf, "lep_tpfms_ndf[2]/F");
  tree->Branch("lep_tpfms_qOverPt", t.lep_tpfms_qOverPt, "lep_tpfms_qOverPt[2]/F");
  tree->Branch("lep_picky_p", t.lep_picky_p, "lep_picky_p[2]/F");
  tree->Branch("lep_picky_pt", t.lep_picky_pt, "lep_picky_pt[2]/F");
  tree->Branch("lep_picky_pt_err", t.lep_picky_pt_err, "lep_picky_pt_err[2]/F");
  tree->Branch("lep_picky_px", t.lep_picky_px, "lep_picky_px[2]/F");
  tree->Branch("lep_picky_py", t.lep_picky_py, "lep_picky_py[2]/F");
  tree->Branch("lep_picky_pz", t.lep_picky_pz, "lep_picky_pz[2]/F");
  tree->Branch("lep_picky_eta", t.lep_picky_eta, "lep_picky_eta[2]/F");
  tree->Branch("lep_picky_phi", t.lep_picky_phi, "lep_picky_phi[2]/F");
  tree->Branch("lep_picky_chi2", t.lep_picky_chi2, "lep_picky_chi2[2]/F");
  tree->Branch("lep_picky_ndf", t.lep_picky_ndf, "lep_picky_ndf[2]/F");
  tree->Branch("lep_picky_qOverPt", t.lep_picky_qOverPt, "lep_picky_qOverPt[2]/F");
  tree->Branch("lep_cocktail_p", t.lep_cocktail_p, "lep_cocktail_p[2]/F");
  tree->Branch("lep_cocktail_pt", t.lep_cocktail_pt, "lep_cocktail_pt[2]/F");
  tree->Branch("lep_cocktail_pt_err", t.lep_cocktail_pt_err, "lep_cocktail_pt_err[2]/F");
  tree->Branch("lep_cocktail_px", t.lep_cocktail_px, "lep_cocktail_px[2]/F");
  tree->Branch("lep_cocktail_py", t.lep_cocktail_py, "lep_cocktail_py[2]/F");
  tree->Branch("lep_cocktail_pz", t.lep_cocktail_pz, "lep_cocktail_pz[2]/F");
  tree->Branch("lep_cocktail_eta", t.lep_cocktail_eta, "lep_cocktail_eta[2]/F");
  tree->Branch("lep_cocktail_phi", t.lep_cocktail_phi, "lep_cocktail_phi[2]/F");
  tree->Branch("lep_cocktail_chi2", t.lep_cocktail_chi2, "lep_cocktail_chi2[2]/F");
  tree->Branch("lep_cocktail_ndf", t.lep_cocktail_ndf, "lep_cocktail_ndf[2]/F");
  tree->Branch("lep_cocktail_qOverPt", t.lep_cocktail_qOverPt, "lep_cocktail_qOverPt[2]/F");
  tree->Branch("lep_cocktail_choice", t.lep_cocktail_choice, "lep_cocktail_choice[2]/S");
  tree->Branch("lep_tuneP_p", t.lep_tuneP_p, "lep_tuneP_p[2]/F");
  tree->Branch("lep_tuneP_pt", t.lep_tuneP_pt, "lep_tuneP_pt[2]/F");
  tree->Branch("lep_tuneP_pt_err", t.lep_tuneP_pt_err, "lep_tuneP_pt_err[2]/F");
  tree->Branch("lep_tuneP_px", t.lep_tuneP_px, "lep_tuneP_px[2]/F");
  tree->Branch("lep_tuneP_py", t.lep_tuneP_py, "lep_tuneP_py[2]/F");
  tree->Branch("lep_tuneP_pz", t.lep_tuneP_pz, "lep_tuneP_pz[2]/F");
  tree->Branch("lep_tuneP_eta", t.lep_tuneP_eta, "lep_tuneP_eta[2]/F");
  tree->Branch("lep_tuneP_phi", t.lep_tuneP_phi, "lep_tuneP_phi[2]/F");
  tree->Branch("lep_tuneP_chi2", t.lep_tuneP_chi2, "lep_tuneP_chi2[2]/F");
  tree->Branch("lep_tuneP_ndf", t.lep_tuneP_ndf, "lep_tuneP_ndf[2]/F");
  tree->Branch("lep_tuneP_qOverPt", t.lep_tuneP_qOverPt, "lep_tuneP_qOverPt[2]/F");
  tree->Branch("lep_triggerMatchPt", t.lep_triggerMatchPt, "lep_triggerMatchPt[2]/F");
  tree->Branch("lep_triggerMatchEta", t.lep_triggerMatchEta, "lep_triggerMatchEta[2]/F");
  tree->Branch("lep_chi2dof", t.lep_chi2dof, "lep_chi2dof[2]/F");
  tree->Branch("lep_dB", t.lep_dB, "lep_dB[2]/F");
  tree->Branch("lep_sumPt", t.lep_sumPt, "lep_sumPt[2]/F");
  tree->Branch("lep_emEt", t.lep_emEt, "lep_emEt[2]/F");
  tree->Branch("lep_hadEt", t.lep_hadEt, "lep_hadEt[2]/F");
  tree->Branch("lep_hoEt", t.lep_hoEt, "lep_hoEt[2]/F");
  tree->Branch("lep_pfIso", t.lep_pfIso, "lep_pfIso[2]/F");
  tree->Branch("lep_pfIsoDB", t.lep_pfIsoDB, "lep_pfIsoDB[2]/F");
  tree->Branch("lep_timeNdof", t.lep_timeNdof, "lep_timeNdof[2]/I");
  tree->Branch("lep_timeInOut", t.lep_timeInOut, "lep_timeInOut[2]/F");
  tree->Branch("lep_timeOutIn", t.lep_timeOutIn, "lep_timeOutIn[2]/F");
  tree->Branch("lep_timeInOutErr", t.lep_timeInOutErr, "lep_timeInOutErr[2]/F");
  tree->Branch("lep_timeOutInErr", t.lep_timeOutInErr, "lep_timeOutInErr[2]/F");
  tree->Branch("lep_min_muon_dR", t.lep_min_muon_dR, "lep_min_muon_dR[2]/F");
  tree->Branch("lep_tk_numberOfValidTrackerHits", t.lep_tk_numberOfValidTrackerHits, "lep_tk_numberOfValidTrackerHits[2]/S");
  tree->Branch("lep_tk_numberOfValidTrackerLayers", t.lep_tk_numberOfValidTrackerLayers, "lep_tk_numberOfValidTrackerLayers[2]/S");
  tree->Branch("lep_tk_numberOfValidPixelHits", t.lep_tk_numberOfValidPixelHits, "lep_tk_numberOfValidPixelHits[2]/S");
  tree->Branch("lep_glb_numberOfValidTrackerHits", t.lep_glb_numberOfValidTrackerHits, "lep_glb_numberOfValidTrackerHits[2]/S");
  tree->Branch("lep_glb_numberOfValidTrackerLayers", t.lep_glb_numberOfValidTrackerLayers, "lep_glb_numberOfValidTrackerLayers[2]/S");
  tree->Branch("lep_glb_numberOfValidPixelHits", t.lep_glb_numberOfValidPixelHits, "lep_glb_numberOfValidPixelHits[2]/S");
  tree->Branch("lep_glb_numberOfValidMuonHits", t.lep_glb_numberOfValidMuonHits, "lep_glb_numberOfValidMuonHits[2]/S");
  tree->Branch("lep_glb_numberOfValidMuonDTHits", t.lep_glb_numberOfValidMuonDTHits, "lep_glb_numberOfValidMuonDTHits[2]/S");
  tree->Branch("lep_glb_numberOfValidMuonCSCHits", t.lep_glb_numberOfValidMuonCSCHits, "lep_glb_numberOfValidMuonCSCHits[2]/S");
  tree->Branch("lep_glb_numberOfValidMuonRPCHits", t.lep_glb_numberOfValidMuonRPCHits, "lep_glb_numberOfValidMuonRPCHits[2]/S");
  tree->Branch("lep_glb_muonStationsWithValidHits", t.lep_glb_muonStationsWithValidHits, "lep_glb_muonStationsWithValidHits[2]/S");
  tree->Branch("lep_glb_dtStationsWithValidHits", t.lep_glb_dtStationsWithValidHits, "lep_glb_dtStationsWithValidHits[2]/S");
  tree->Branch("lep_glb_cscStationsWithValidHits", t.lep_glb_cscStationsWithValidHits, "lep_glb_cscStationsWithValidHits[2]/S");
  tree->Branch("lep_glb_rpcStationsWithValidHits", t.lep_glb_rpcStationsWithValidHits, "lep_glb_rpcStationsWithValidHits[2]/S");
  tree->Branch("lep_glb_innermostMuonStationWithValidHits", t.lep_glb_innermostMuonStationWithValidHits, "lep_glb_innermostMuonStationWithValidHits[2]/S");
  tree->Branch("lep_glb_outermostMuonStationWithValidHits", t.lep_glb_outermostMuonStationWithValidHits, "lep_glb_outermostMuonStationWithValidHits[2]/S");
  tree->Branch("lep_numberOfMatches", t.lep_numberOfMatches, "lep_numberOfMatches[2]/S");
  tree->Branch("lep_numberOfMatchedStations", t.lep_numberOfMatchedStations, "lep_numberOfMatchedStations[2]/S");
  tree->Branch("lep_numberOfMatchedRPCLayers",t.lep_numberOfMatchedRPCLayers, "lep_numberOfMatchedRPCLayers[2]/S");
  tree->Branch("lep_stationMask", t.lep_stationMask, "lep_stationMask[2]/I");
  tree->Branch("lep_numberOfChambers", t.lep_numberOfChambers, "lep_numberOfChambers[2]/I");
  tree->Branch("lep_numberOfChambersNoRPC", t.lep_numberOfChambersNoRPC, "lep_numberOfChambersNoRPC[2]/I");
  tree->Branch("lep_stationGapMaskDistance", t.lep_stationGapMaskDistance, "lep_stationGapMaskDistance[2]/I");
  tree->Branch("lep_stationGapMaskPull", t.lep_stationGapMaskPull, "lep_stationGapMaskPull[2]/I");
  tree->Branch("lep_isGlobalMuon", t.lep_isGlobalMuon, "lep_isGlobalMuon[2]/O");
  tree->Branch("lep_isTrackerMuon", t.lep_isTrackerMuon, "lep_isTrackerMuon[2]/O");
  tree->Branch("GoodDataRan", &t.GoodDataRan, "GoodDataRan/O");
  tree->Branch("GoodVtx", &t.GoodVtx, "GoodVtx/O");
  tree->Branch("METFilter", &t.METFilter, "METFilter/O");
  tree->Branch("met_pt", &t.met_pt, "met_pt/F");
  tree->Branch("met_phi", &t.met_phi, "met_phi/F");
  tree->Branch("nJets", &t.nJets, "nJets/I");
  tree->Branch("jet_pt", t.jet_pt, "jet_pt[10]/F");
  tree->Branch("jet_eta", t.jet_eta, "jet_eta[10]/F");
  tree->Branch("jet_phi", t.jet_phi, "jet_phi[10]/F");
  tree->Branch("jet_btag", t.jet_btag, "jet_btag[10]/F");
  //Add jet information
  tree->Branch("jet_E", t.jet_E, "jet_E[10]/F");
  //////////////////////////////////////////////////////
  if (fill_gen_info) {
    tree->Branch("genWeight", &t.genWeight, "genWeight/F");
    tree->Branch("nQuark_initialState", &t.nQuark_initialState, "nQuark_initialState/I");
    tree->Branch("nGluon_initialState", &t.nGluon_initialState, "nGluon_initialState/I");
    tree->Branch("nGisr", &t.nGisr, "nGisr/I");
    tree->Branch("gen_res_mass", &t.gen_res_mass, "gen_res_mass/F");
    tree->Branch("gen_res_pt", &t.gen_res_pt, "gen_res_pt/F");
    tree->Branch("gen_res_rap", &t.gen_res_rap, "gen_res_rap/F");
    tree->Branch("gen_res_eta", &t.gen_res_eta, "gen_res_eta/F");
    tree->Branch("gen_res_phi", &t.gen_res_phi, "gen_res_phi/F");
    tree->Branch("gen_dil_mass", &t.gen_dil_mass, "gen_dil_mass/F");
    tree->Branch("gen_dil_pt", &t.gen_dil_pt, "gen_dil_pt/F");
    tree->Branch("gen_dil_rap", &t.gen_dil_rap, "gen_dil_rap/F");
    tree->Branch("gen_dil_eta", &t.gen_dil_eta, "gen_dil_eta/F");
    tree->Branch("gen_dil_phi", &t.gen_dil_phi, "gen_dil_phi/F");
    tree->Branch("gen_dil_dR", &t.gen_dil_dR, "gen_dil_dR/F");
    tree->Branch("gen_dil_dPhi", &t.gen_dil_dPhi, "gen_dil_dPhi/F");
    tree->Branch("gen_lep_p", t.gen_lep_p, "gen_lep_p[2]/F");
    tree->Branch("gen_lep_pt", t.gen_lep_pt, "gen_lep_pt[2]/F");
    tree->Branch("gen_lep_px", t.gen_lep_px, "gen_lep_px[2]/F");
    tree->Branch("gen_lep_py", t.gen_lep_py, "gen_lep_py[2]/F");
    tree->Branch("gen_lep_pz", t.gen_lep_pz, "gen_lep_pz[2]/F");
    tree->Branch("gen_lep_E", t.gen_lep_E, "gen_lep_E[2]/F");
    tree->Branch("gen_lep_eta", t.gen_lep_eta, "gen_lep_eta[2]/F");
    tree->Branch("gen_lep_phi", t.gen_lep_phi, "gen_lep_phi[2]/F");
    tree->Branch("gen_lep_qOverPt", t.gen_lep_qOverPt, "gen_lep_qOverPt[2]/F");
    tree->Branch("gen_lep_noib_pt", t.gen_lep_noib_pt, "gen_lep_noib_pt[2]/F");
    tree->Branch("gen_lep_noib_px", t.gen_lep_noib_px, "gen_lep_noib_px[2]/F");
    tree->Branch("gen_lep_noib_py", t.gen_lep_noib_py, "gen_lep_noib_py[2]/F");
    tree->Branch("gen_lep_noib_pz", t.gen_lep_noib_pz, "gen_lep_noib_pz[2]/F");
    tree->Branch("gen_lep_noib_E", t.gen_lep_noib_E, "gen_lep_noib_E[2]/F");
    tree->Branch("gen_lep_noib_eta", t.gen_lep_noib_eta, "gen_lep_noib_eta[2]/F");
    tree->Branch("gen_lep_noib_phi", t.gen_lep_noib_phi, "gen_lep_noib_phi[2]/F");
    tree->Branch("gen_lep_noib_qOverPt", t.gen_lep_noib_qOverPt, "gen_lep_noib_qOverPt[2]/F");
  }

  tree->SetAlias("OppSign",  "lep_id[0]*lep_id[1] < 0");
  tree->SetAlias("SameSign", "lep_id[0]*lep_id[1] > 0");
  tree->SetAlias("Dimu",     "abs(lep_id[0]*lep_id[1]) == 169");
  tree->SetAlias("Emu",      "abs(lep_id[0]*lep_id[1]) == 143");

#define offlineMinPt "53"
#define triggerMatchMinPt "50"

//#define offlineMinPt "25"
//#define triggerMatchMinPt "25"
#define triggerMatchMaxEta "2.1"

//  tree->SetAlias("trigger_match_0", "lep_triggerMatchPt[0] > " triggerMatchMinPt " && abs(lep_triggerMatchEta[0]) < " triggerMatchMaxEta);
//  tree->SetAlias("trigger_match_1", "lep_triggerMatchPt[1] > " triggerMatchMinPt " && abs(lep_triggerMatchEta[1]) < " triggerMatchMaxEta);
  tree->SetAlias("trigger_match_0", "lep_triggerMatchPt[0] > " triggerMatchMinPt );
  tree->SetAlias("trigger_match_1", "lep_triggerMatchPt[1] > " triggerMatchMinPt );
  tree->SetAlias("triggerMatched", "trigger_match_0 || trigger_match_1");

  // tree->SetAlias("GoodData", "GoodDataRan && HLTPhysicsDeclared && NoScraping && GoodVtx");
  tree->SetAlias("GoodData", "GoodDataRan && HLTPhysicsDeclared && GoodVtx");

  tree->SetAlias("extraDimuonCuts", "cos_angle > -0.9998 && vertex_chi2 < 20");

  TString loose_2010 =
    "lep_isGlobalMuon[X] && "           \
    "lep_pt[X] > " offlineMinPt " && "          \
    "lep_tk_numberOfValidTrackerHits[X] >= 10 && "      \
    "lep_sumPt[X] / lep_tk_pt[X] < 0.1";

  TString tight_2010 =
    "abs(lep_dB[X]) < 0.2 && "            \
    "lep_chi2dof[X] < 10 && "           \
    "lep_tk_numberOfValidPixelHits[X] >= 1 && "       \
    "lep_glb_muonStationsWithValidHits[X] >= 2 && "     \
    "lep_isTrackerMuon[X] && "            \
    "lep_triggerMatchPt[X] > " triggerMatchMinPt;

  TString vbtf =
    "lep_isGlobalMuon[X] && "           \
    "lep_isTrackerMuon[X] && "            \
    "lep_tk_pt[X] > " offlineMinPt " && "       \
    "abs(lep_tk_eta[X]) < 2.1 && "          \
    "abs(lep_dB[X]) < 0.2 && "            \
    "lep_sumPt[X] < 3 && "            \
    "lep_glb_numberOfValidTrackerHits[X] > 10 && "      \
    "lep_glb_numberOfValidPixelHits[X] > 0 && "       \
    "lep_glb_numberOfValidMuonHits[X] > 0 && "        \
    "lep_numberOfMatches[X] > 1";

  TString loose_no_iso =
    "lep_isGlobalMuon[X] && "           \
    "lep_isTrackerMuon[X] && "            \
    "lep_pt[X] > " offlineMinPt " && "          \
    "abs(lep_dB[X]) < 0.2 && "            \
    "lep_glb_numberOfValidTrackerLayers[X] > 5 && "     \
    "lep_glb_numberOfValidPixelHits[X] >= 1 && "      \
    "lep_glb_numberOfValidMuonHits[X] > 0 && "        \
    "lep_numberOfMatchedStations[X] > 1 && "                            \
    "lep_pt_err[X] / lep_pt[X] < 0.3";

  TString tight_2015 =
    "lep_isGlobalMuon[X] && "           \
    "lep_isTrackerMuon[X] && "            \
    "lep_tuneP_pt[X] > " offlineMinPt " && "        \
    "abs(lep_dB[X]) < 0.2 && "            \
    "lep_glb_numberOfValidTrackerLayers[X] > 5 && "     \
    "lep_glb_numberOfValidPixelHits[X] >= 1 && "      \
    "lep_glb_numberOfValidMuonHits[X] > 0 && "        \
    "lep_numberOfMatchedStations[X] > 1 && "                            \
    "lep_tuneP_pt_err[X] / lep_tuneP_pt[X] < 0.3 && "       \
    "lep_sumPt[X] / lep_tk_pt[X] < 0.1";

  TString loose_no_iso_pt_ptErr = 
    "lep_isGlobalMuon[X] && "           \
    "lep_isTrackerMuon[X] && "            \
    "lep_glb_numberOfValidTrackerLayers[X] > 5 && "     \
    "lep_glb_numberOfValidPixelHits[X] >= 1 && "      \
    "lep_glb_numberOfValidMuonHits[X] > 0 && "        \
    "lep_numberOfMatchedStations[X] > 1 && "        \
    "abs(lep_dB[X]) < 0.2";

  TString loose_2012 = loose_no_iso + " && lep_sumPt[X] / lep_tk_pt[X] < 0.1";

  TString loose_new(loose_2012);
  loose_new.ReplaceAll("lep_glb_numberOfValidTrackerLayers[X] > 5", 
           "lep_glb_numberOfValidTrackerLayers[X] > 8");
  loose_new.ReplaceAll(" && lep_pt_err[X] / lep_pt[X] < 0.3", "");

  TString loose_2011eps(loose_new);
  loose_2011eps.ReplaceAll("lep_glb_numberOfValidTrackerLayers", 
         "lep_glb_numberOfValidTrackerHits");

  tree->SetAlias("loose_2010_0",    replace_all(loose_2010,    "[X]", "[0]"));
  tree->SetAlias("loose_2010_1",    replace_all(loose_2010,    "[X]", "[1]"));
  tree->SetAlias("tight_2010_0",    replace_all(tight_2010,    "[X]", "[0]"));
  tree->SetAlias("tight_2010_1",    replace_all(tight_2010,    "[X]", "[1]"));
  tree->SetAlias("vbtf_0",          replace_all(vbtf,          "[X]", "[0]"));
  tree->SetAlias("vbtf_1",          replace_all(vbtf,          "[X]", "[1]"));
  tree->SetAlias("loose_2011eps_0", replace_all(loose_2011eps, "[X]", "[0]"));
  tree->SetAlias("loose_2011eps_1", replace_all(loose_2011eps, "[X]", "[1]"));
  tree->SetAlias("loose_new_0",     replace_all(loose_new,     "[X]", "[0]"));
  tree->SetAlias("loose_new_1",     replace_all(loose_new,     "[X]", "[1]"));
  tree->SetAlias("loose_no_iso_0",  replace_all(loose_no_iso,  "[X]", "[0]"));
  tree->SetAlias("loose_no_iso_1",  replace_all(loose_no_iso,  "[X]", "[1]"));
  tree->SetAlias("loose_2012_0",    replace_all(loose_2012,    "[X]", "[0]"));
  tree->SetAlias("loose_2012_1",    replace_all(loose_2012,    "[X]", "[1]"));
  tree->SetAlias("tight_2015_0",    replace_all(tight_2015,    "[X]", "[0]"));
  tree->SetAlias("tight_2015_1",    replace_all(tight_2015,    "[X]", "[1]"));
  tree->SetAlias("loose_no_iso_pt_ptErr_0", replace_all(loose_no_iso_pt_ptErr, "[X]", "[0]"));
  tree->SetAlias("loose_no_iso_pt_ptErr_1", replace_all(loose_no_iso_pt_ptErr, "[X]", "[1]"));

  tree->SetAlias("OurSel2010",
     "loose_2010_0 && loose_2010_1 && "     \
     "(tight_2010_0 || tight_2010_1) && "     \
     "OppSign && "            \
     "extraDimuonCuts && "          \
     "GoodData");

  tree->SetAlias("VBTFSel",
     "vbtf_0 && vbtf_1 && "         \
     "triggerMatched && "         \
     "OppSign");

  tree->SetAlias("OurSel2011EPS",
     "loose_2011eps_0 && loose_2011eps_1 && "   \
     "triggerMatched && "         \
     "OppSign && "            \
     "extraDimuonCuts && "          \
     "GoodData");

  tree->SetAlias("OurSelNewNoSign",
     "loose_new_0 && loose_new_1 && "     \
     "triggerMatched && "         \
     "extraDimuonCuts && "          \
     "GoodData");

  tree->SetAlias("OurSelNew",   "OurSelNewNoSign && OppSign");
  tree->SetAlias("OurSelNewSS", "OurSelNewNoSign && SameSign");

  tree->SetAlias("OurSel2012NoSign",
     "loose_2012_0 && loose_2012_1 && "     \
     "triggerMatched && "         \
     "extraDimuonCuts && "          \
     "GoodData");

  tree->SetAlias("OurSel2012",   "OurSel2012NoSign && OppSign");
  tree->SetAlias("OurSel2012SS", "OurSel2012NoSign && SameSign");

  tree->SetAlias("OurSel2015",
                 "OppSign && "            \
                 "tight_2015_0 && tight_2015_1 && "     \
                 "triggerMatched && "           \
                 "extraDimuonCuts && "          \
                 "GoodData");

  // For e-mu dileptons, below we always put the muon in [0] and the
  // electron in [1], so don't have to check the other combination.
  tree->SetAlias("EmuSelNoSign",
     "abs(lep_id[1]) == 11 && "       \
     "lep_heep_id[1] == 0 && "        \
     "loose_2012_0 && "         \
     "trigger_match_0 && "          \
     "GoodData");

  tree->SetAlias("EmuSel", "EmuSelNoSign && OppSign");

  //NPV_no_cuts = new TH1D("NPV_no_cuts", "NPV_no_cuts", 80, 0,80);


  for (int i=0; i<nCuts;i++){
    cutflow_study_dilep[i] = 0;
  }

  NPV_no_cuts = fs->make<TH1D>("NPV_no_cuts", "NPV_no_cuts", 80, 0,80);
  NPV_dilepton_loop = fs->make<TH1D>("NPV_dilepton_loop", "NPV_dilepton_loop", 80, 0,80);
  NPV_cuts = fs->make<TH1D>("NPV_cuts", "NPV_cuts", 80, 0,80);

  TH1F_cutlfow = fs->make<TH1F>("TH1F_cutlfow", "TH1F_cutlfow", 10, -.5,8.5);

  TH1F_muonPt_all = fs->make<TH1F>("muonPt_all", "muonPt_all", 100, 0, 150);
  TH1F_muonPt_trigg = fs->make<TH1F>("muonPt_trigg", "muonPt_trigg", 100, 0, 150);
  TH1F_muonPt_one = fs->make<TH1F>("muonPt_one", "muonPt_one", 100, 0, 150);
  TH1F_muonPt_two = fs->make<TH1F>("muonPt_two", "muonPt_two", 100, 0, 150);
  TH1F_muonPt_oppSign = fs->make<TH1F>("muonPt_oppSign", "muonPt_oppSign", 100, 0, 150);
  TH1F_muonPt_deltaR = fs->make<TH1F>("muonPt_deltaR", "muonPt_deltaR", 100, 0, 150);
  TH1F_muonPt_mass = fs->make<TH1F>("muonPt_mass", "muonPt_mass", 100, 0, 150);

  TH1F_muonEta_all = fs->make<TH1F>("muonEta_all", "muonEta_all", 100, -2.5,2.5);
  TH1F_muonEta_trigg = fs->make<TH1F>("muonEta_trigg", "muonEta_trigg", 100, -2.5,2.5);
  TH1F_muonEta_one = fs->make<TH1F>("muonEta_one", "muonEta_one", 100, -2.5,2.5);
  TH1F_muonEta_two = fs->make<TH1F>("muonEta_two", "muonEta_two", 100, -2.5,2.5);
  TH1F_muonEta_oppSign = fs->make<TH1F>("muonEta_oppSign", "muonEta_oppSign", 100, -2.5,2.5);
  TH1F_muonEta_deltaR = fs->make<TH1F>("muonEta_deltaR", "muonEta_deltaR", 100, -2.5,2.5);
  TH1F_muonEta_mass = fs->make<TH1F>("muonEta_mass", "muonEta_mass", 100, -2.5,2.5);


}

template <typename T>
float userFloat(const T& patobj, const char* name, float def=-999.) {
  return patobj.hasUserFloat(name) ? patobj.userFloat(name) : def;
}

template <typename T>
int userInt(const T& patobj, const char* name, int def=-999) {
  return patobj.hasUserInt(name) ? patobj.userInt(name) : def;
}

void SimpleNtupler_miniAOD_noDiLep::analyze(const edm::Event& event, const edm::EventSetup&) {
  memset(&t, 0, sizeof(tree_t));

  //
  // Event Information
  //
  t.run = event.id().run();
  t.lumi = event.luminosityBlock();
  t.event = event.id().event();

  //std::cout << t.run << " " << t.lumi << " " << t.event << std::endl;

  // Get Trigger information



/*
  edm::Handle<edm::TriggerResults> respat;
  //event.getByLabel(edm::InputTag("TriggerResults", "", "PAT"), respat);
  event.getByLabel(TriggerResults_src, respat);
  
  const edm::TriggerNames& namespat = event.triggerNames(*respat);
 

  if (namespat.triggerIndex("Flag_goodVertices") < respat->size()) {
    t.GoodDataRan = 1;
    t.GoodVtx = respat->accept(namespat.triggerIndex("Flag_goodVertices"));
    bool metFilterAccept = true;
    for ( std::vector<edm::InputTag>::iterator filterTag_i = filterTags.begin(); filterTag_i != filterTags.end(); ++filterTag_i ) {
      std::string filterTag = (*filterTag_i).label(); 
      metFilterAccept  *= respat->accept(namespat.triggerIndex(filterTag)); 
    }
    t.METFilter = metFilterAccept;
  }*/

  // Get Beamspot information
  edm::Handle<reco::BeamSpot> bs;
  event.getByLabel(beamspot_src, bs);
  t.beamspot_x     = bs->x0();
  t.beamspot_x_err = bs->x0Error();
  t.beamspot_y     = bs->y0();
  t.beamspot_y_err = bs->y0Error();
  t.beamspot_z     = bs->z0();
  t.beamspot_z_err = bs->z0Error();

  // Get Vertex information
  edm::Handle<reco::VertexCollection> pvs;
  event.getByLabel(vertices_src, pvs);
  t.nvertices = 0;
  BOOST_FOREACH(const reco::Vertex& vtx, *pvs)
    if (vtx.ndof() > 4 && fabs(vtx.z()) <= 24 && fabs(vtx.position().rho()) <= 2)
      t.nvertices += 1;

  NPV_no_cuts->Fill(t.nvertices);
  

  if (fill_gen_info) {


    //access initial state
    edm::Handle<GenEventInfoProduct> genEvtInfo;
    event.getByLabel( "generator", genEvtInfo) ;
    double theWeight = genEvtInfo->weight();

    const gen::PdfInfo *pdfInfo = (genEvtInfo->hasPDF()) ? genEvtInfo->pdf() : 0;

    int nQuark = 0;
    int nGluon = 0;
    ////std::cout << pdfInfo->id.first << std::endl;
    ////std::cout << pdfInfo->id.second << std::endl;

    if (abs(pdfInfo->id.first) < 9) nQuark = nQuark+1;
    if (abs(pdfInfo->id.first) == 21) nGluon = nGluon+1;

    if (abs(pdfInfo->id.second) < 9) nQuark = nQuark+1;
    if (abs(pdfInfo->id.second) == 21) nGluon = nGluon+1;

    t.nQuark_initialState =   nQuark;

    t.nGluon_initialState = nGluon;


    //edm::Handle<std::vector<reco::GenParticle>> genPart;
    //event.getByLabel( "prunedGenParticles", genPart) ;
    //const reco::GenParticle genParticles = *(genPart.product());  
//
//
//
//
    //////std::cout << theWeight << std::endl;
    //for(size_t i=0; i<2;i++){
    //  //std::cout << i << std::endl; //(*genParticles)[i].pdgId() << std::endl;
    //}



    
    // This only works for DY/Z'/RSG events, and really just for PYTHIA!
    hardInteraction->Fill(event);
    int EventWeight = 1.;
    edm::Handle<GenEventInfoProduct> gen_ev_info;
    event.getByLabel(genEventInfo_, gen_ev_info);
    EventWeight = gen_ev_info->weight();
    t.genWeight = ( EventWeight > 0 ) ? 1 : -1;

    t.genWeight = theWeight;
    //t.genWeight = EventWeight;
    
    
    //
    // Store Generator Level information
    //
//     if(hardInteraction->IsValid()){
  if(hardInteraction->IsValidForRes()){

      //t.nQuark_initialState = hardInteraction->nQinitialState;
      //t.nGluon_initialState = hardInteraction->nGinitialState;
      //t.nGisr = hardInteraction->gISR;

      ////std::cout << "tetst" << std::endl;
      ////std::cout << "nquark " << t.nQuark_initialState << std::endl;
      ////std::cout << "nglu " <<  t.nGluon_initialState << std::endl;
      ////std::cout << "gisr " <<  t.nGisr << std::endl;

      t.gen_res_mass = hardInteraction->resonance->mass();
      t.gen_res_pt   = hardInteraction->resonance->pt();
      t.gen_res_rap  = hardInteraction->resonance->rapidity();
      t.gen_res_eta  = hardInteraction->resonance->eta();
      t.gen_res_phi  = hardInteraction->resonance->phi();
      
      t.gen_dil_mass = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).mass();//hardInteraction->dilepton().mass();
      t.gen_dil_pt   = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).pt();
      t.gen_dil_rap  = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).Rapidity();
      t.gen_dil_eta  = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).eta();
      t.gen_dil_phi  = (hardInteraction->lepPlusNoIB->p4() + hardInteraction->lepMinusNoIB->p4()).phi();
      t.gen_dil_dR   = deltaR(*hardInteraction->lepMinusNoIB, *hardInteraction->lepPlusNoIB);
      t.gen_dil_dPhi = deltaPhi(*hardInteraction->lepMinusNoIB, *hardInteraction->lepPlusNoIB);
//       
      t.gen_lep_p[0]  = hardInteraction->lepMinusNoIB->p();
      t.gen_lep_pt[0]  = hardInteraction->lepMinusNoIB->pt();
      t.gen_lep_px[0]  = hardInteraction->lepMinusNoIB->px();
      t.gen_lep_py[0]  = hardInteraction->lepMinusNoIB->py();
      t.gen_lep_pz[0]  = hardInteraction->lepMinusNoIB->pz();
      t.gen_lep_E[0]  = hardInteraction->lepMinusNoIB->energy();
      t.gen_lep_eta[0] = hardInteraction->lepMinusNoIB->eta();
      t.gen_lep_phi[0] = hardInteraction->lepMinusNoIB->phi();
      t.gen_lep_qOverPt[0] = hardInteraction->lepMinusNoIB->charge() / hardInteraction->lepMinusNoIB->pt();
//       
      t.gen_lep_p[1]  = hardInteraction->lepPlusNoIB->p();
      t.gen_lep_pt[1]  = hardInteraction->lepPlusNoIB->pt();
      t.gen_lep_px[1]  = hardInteraction->lepPlusNoIB->px();
      t.gen_lep_py[1]  = hardInteraction->lepPlusNoIB->py();
      t.gen_lep_pz[1]  = hardInteraction->lepPlusNoIB->pz();
      t.gen_lep_E[1]  = hardInteraction->lepPlusNoIB->energy();
      t.gen_lep_eta[1] = hardInteraction->lepPlusNoIB->eta();
      t.gen_lep_phi[1] = hardInteraction->lepPlusNoIB->phi();
      t.gen_lep_qOverPt[1] = hardInteraction->lepPlusNoIB->charge() / hardInteraction->lepPlusNoIB->pt();
      
      /*
       t.gen_lep_noib_pt[0]  = hardInteraction->lepMinusNoIB->pt();
       t.gen_lep_noib_px[0]  = hardInteraction->lepMinusNoIB->px();
       t.gen_lep_noib_py[0]  = hardInteraction->lepMinusNoIB->py();
       t.gen_lep_noib_pz[0]  = hardInteraction->lepMinusNoIB->pz();
       t.gen_lep_noib_e[0]  = hardInteraction->lepMinusNoIB->energy();
       t.gen_lep_noib_eta[0] = hardInteraction->lepMinusNoIB->eta();
       t.gen_lep_noib_phi[0] = hardInteraction->lepMinusNoIB->phi();
       
       t.gen_lep_noib_pt[1]  = hardInteraction->lepPlusNoIB->pt();
       t.gen_lep_noib_px[1]  = hardInteraction->lepMinusNoIB->px();
       t.gen_lep_noib_py[1]  = hardInteraction->lepMinusNoIB->py();
       t.gen_lep_noib_pz[1]  = hardInteraction->lepMinusNoIB->pz();
       t.gen_lep_noib_e[1]  = hardInteraction->lepMinusNoIB->energy();
       t.gen_lep_noib_eta[1] = hardInteraction->lepPlusNoIB->eta();
       t.gen_lep_noib_phi[1] = hardInteraction->lepPlusNoIB->phi();
      */
    } // end if hardInteraction->IsValid()
    
  } // end if fill_gen_info
  


  //
  //trigger
  //
  edm::Handle<edm::TriggerResults> respat;
  event.getByLabel(TriggerResults_src, respat);

  const edm::TriggerNames &names = event.triggerNames(*respat);

  bool passTrigger = false;

  bool passTriggerMu = false;
  bool passTriggerTkMu = false;
  std::string name = "";
  for (unsigned int i = 0, n = respat->size(); i < n; ++i) {
     name = names.triggerName(i);
  

    bool pass = false;
    int accept = respat->accept(i) ;
    int prescale = 0;
    if (accept ==1 ) pass = true;
    ////std::cout << name << std::endl;
    if(pass && (name.find("HLT_Mu50_v") !=std::string::npos )  ) {
      passTriggerMu = true;
    }
    if(pass && (name.find("HLT_TkMu50_v") !=std::string::npos )  ) {
      passTriggerTkMu = true;
    }
  }
  if (passTriggerMu || passTriggerTkMu){
    passTrigger = true;
  }

  //
  // Get dilepton collection
  //

  bool atLeastOneHighPtMuon = false;
  bool atLeastTwoHighPtMuon = false;
  bool atLeastTwoOppositeSignHighPtMuon = false;
  bool passDeltaR = false;
  bool passMassCut = false;

  edm::Handle<std::vector<pat::Muon>> muons;
  event.getByLabel(mu_src, muons);

  std::vector<TLorentzVector> muonTLVector = {};
  std::vector<int> muonCharge = {};
  TLorentzVector muonTL;
  double leadingPt = -1;
  int count = 0;

  for (std::vector<pat::Muon>::const_iterator itMuon = muons->begin(); itMuon != muons->end(); itMuon++) {
    ////std::cout << "muon pt " << itMuon->pt() << std::endl;
    if(count == 0) leadingPt = itMuon->pt();
    count = count + 1;
  
    //std::cout << itMuon->pt() << " " << itMuon->charge() << " " << itMuon->CutBasedIdGlobalHighPt << " " << itMuon->isolationR03().sumPt/ itMuon->pt() << " " << abs(itMuon->eta())  << std::endl;
    if(not passTrigger) continue;

    //muCut = "Muon_pt > 53 and Muon_highPtId > 1 and Muon_pfRelIso03_chg < 0.1 and -2.4 < Muon_eta and Muon_eta < 2.4"
    if (not (itMuon->pt() > 53)) continue;
    if (not (itMuon->CutBasedIdGlobalHighPt > 1) ) continue;


    //old high pt id:
   // if (not (itMuon->isGlobalMuon())) continue;
    //if (not (itMuon->globalTrack()->hitPattern().numberOfValidMuonHits() > 0)) continue;
    //if (not (itMuon->numberOfMatchedStations() > 1)) continue;
    //if (not (itMuon->tunePMuonBestTrack()->ptError()/itMuon->tunePMuonBestTrack()->pt() < .3)) continue;
    //if (not (abs(itMuon->dB()) < 2)) continue;
    //if (not (abs(itMuon->dZ()) < 2)) continue;
    //if (not (abs(itMuon->innerTrack()->dz(vertex->position()) < .5)) ) continue;   
    //if (not (itMuon->globalTrack()->hitPattern().numberOfValidPixelHits() > 0)) continue;
    //if (not (itMuon->globalTrack()->hitPattern().trackerLayersWithMeasurement() > 5)) continue;


    if (not (itMuon->pfIsolationR03().sumChargedHadronPt/ itMuon->pt() < 0.10)) continue; //probably not real rel iso, i think it needs to be tracker pt
 //   if (not (itMuon->isolationR03().sumPt/ itMuon->pt() < 0.10)) continue; //probably not real rel iso, i think it needs to be tracker pt
    if (not (abs(itMuon->eta()) < 2.4)) continue;

    muonTL.SetPtEtaPhiM(itMuon->pt(),itMuon->eta(),itMuon->phi(),.105658);
    muonTLVector.push_back( muonTL );
    muonCharge.push_back(itMuon->charge());


  }

  //at least one

  if(muonTLVector.size() > 0){
    atLeastOneHighPtMuon =  true;
  }


//at least 2
  if(muonTLVector.size() > 1){
     atLeastTwoHighPtMuon =  true;
  }


  TLorentzVector positiveMuonTL;
  TLorentzVector negativeMuonTL;
  int charge= 0;
  if(muonTLVector.size() > 1){
    int index = 0;
     for (std::vector<TLorentzVector>::const_iterator itMuon = muonTLVector.begin(); itMuon != muonTLVector.end(); itMuon++) {
      //std::cout << itMuon->Pt() << std::endl;
      charge = muonCharge.at(index);
      index = index+1;

      if (charge > 0 and positiveMuonTL.Pt() < itMuon->Pt()){
        positiveMuonTL.SetPtEtaPhiM(itMuon->Pt(),itMuon->Phi(),itMuon->Eta(),.105658);
      }

      if (charge < 0 and negativeMuonTL.Pt() < itMuon->Pt()){
        negativeMuonTL.SetPtEtaPhiM(itMuon->Pt(),itMuon->Phi(),itMuon->Eta(),.105658);
      }
     }
  }

  //opposite sign
  if (atLeastTwoHighPtMuon && positiveMuonTL.Pt() > 0 && negativeMuonTL.Pt() > 0){
    atLeastTwoOppositeSignHighPtMuon =  true;
  }

  //deltaR
  if (atLeastTwoOppositeSignHighPtMuon && negativeMuonTL.DeltaR(positiveMuonTL) > .1){
    passDeltaR =  true;
  }

//mass
  if (passDeltaR && (negativeMuonTL + positiveMuonTL).M() > 55){
    passMassCut = true;
  }



  int cut_index = 0;
  cutflow_study_dilep[cut_index] = cutflow_study_dilep[cut_index] + 1;
  TH1F_cutlfow->Fill(cut_index);
  cut_index = cut_index+1;
  //std::cout << "Event" << std::endl;


  for (std::vector<pat::Muon>::const_iterator itMuon = muons->begin(); itMuon != muons->end(); itMuon++) {
    TH1F_muonPt_all->Fill(itMuon->pt());
    TH1F_muonEta_all->Fill(itMuon->eta());
  }

  if (passTrigger){
    cutflow_study_dilep[cut_index] = cutflow_study_dilep[cut_index] + 1;
    TH1F_cutlfow->Fill(cut_index);
    //std::cout << "passTrigger" << std::endl;
    for (std::vector<pat::Muon>::const_iterator itMuon = muons->begin(); itMuon != muons->end(); itMuon++) {
      TH1F_muonPt_trigg->Fill(itMuon->pt());
      TH1F_muonEta_trigg->Fill(itMuon->eta());
    }
  }
  cut_index = cut_index+1;


  if (atLeastOneHighPtMuon){
    cutflow_study_dilep[cut_index] = cutflow_study_dilep[cut_index] + 1;
    TH1F_cutlfow->Fill(cut_index);
    //std::cout << "oneMuon" << std::endl;
    for (std::vector<pat::Muon>::const_iterator itMuon = muons->begin(); itMuon != muons->end(); itMuon++) {
      TH1F_muonPt_one->Fill(itMuon->pt());
      TH1F_muonEta_one->Fill(itMuon->eta());
    }
  }
  cut_index = cut_index+1;


  if (atLeastTwoHighPtMuon){
    cutflow_study_dilep[cut_index] = cutflow_study_dilep[cut_index] + 1;
    TH1F_cutlfow->Fill(cut_index);
    //std::cout << "twoMuon" << std::endl;
    for (std::vector<pat::Muon>::const_iterator itMuon = muons->begin(); itMuon != muons->end(); itMuon++) {
      TH1F_muonPt_two->Fill(itMuon->pt());
      TH1F_muonEta_two->Fill(itMuon->eta());
    }
  }
  cut_index = cut_index+1;

  if (atLeastTwoOppositeSignHighPtMuon){
    cutflow_study_dilep[cut_index] = cutflow_study_dilep[cut_index] + 1;
    TH1F_cutlfow->Fill(cut_index);
    //std::cout << "oppSign" << std::endl;
    for (std::vector<pat::Muon>::const_iterator itMuon = muons->begin(); itMuon != muons->end(); itMuon++) {
         TH1F_muonPt_oppSign->Fill(itMuon->pt());
         TH1F_muonEta_oppSign->Fill(itMuon->eta());
    }
  }
  cut_index = cut_index+1;

  if (passDeltaR){
    cutflow_study_dilep[cut_index] = cutflow_study_dilep[cut_index] + 1;
    TH1F_cutlfow->Fill(cut_index);
    //std::cout << "deltaR" << std::endl;
    for (std::vector<pat::Muon>::const_iterator itMuon = muons->begin(); itMuon != muons->end(); itMuon++) {
         TH1F_muonPt_deltaR->Fill(itMuon->pt());
         TH1F_muonEta_deltaR->Fill(itMuon->eta());
    }
  }
  cut_index = cut_index+1;

  if (passMassCut){
    cutflow_study_dilep[cut_index] = cutflow_study_dilep[cut_index] + 1;
    TH1F_cutlfow->Fill(cut_index);
    //std::cout << "mass" << std::endl;
    for (std::vector<pat::Muon>::const_iterator itMuon = muons->begin(); itMuon != muons->end(); itMuon++) {
      TH1F_muonPt_mass->Fill(itMuon->pt());
      TH1F_muonEta_mass->Fill(itMuon->eta());
    }
  }
  cut_index = cut_index+1;






  if (1==2){
  std::cout << t.run << "," << t.lumi << "," << t.event <<  ","  << passTrigger << "," << atLeastOneHighPtMuon << "," << atLeastTwoHighPtMuon << "," << atLeastTwoOppositeSignHighPtMuon << "," << passDeltaR << "," << passMassCut << "," << positiveMuonTL.Pt() << "," << negativeMuonTL.Pt() << "," << leadingPt <<  std::endl;
 }

 tree->Fill();
  // 
  // Put additional event info here
  // MET, Jets, etc.
  //

} // end SimpleNtupler_miniAOD_noDiLep::analyze

void SimpleNtupler_miniAOD_noDiLep::endJob() 
{
  for(int cutIndex = 0; cutIndex < nCuts; cutIndex++){
    std::cout << "nCuts " << cutIndex << " nEvents " << cutflow_study_dilep[cutIndex] << std::endl;
  }
}

DEFINE_FWK_MODULE(SimpleNtupler_miniAOD_noDiLep);


