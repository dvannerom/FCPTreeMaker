// -*- C++ -*-
//
// Package:    AnalyzerArea/FCPTreeMaker
// Class:      FCPTreeMaker
// 
/**\class FCPTreeMaker FCPTreeMaker.cc AnalyzerArea/FCPTreeMaker/plugins/FCPTreeMaker.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  David Vannerom
//         Created:  Fri, 27 Apr 2018 12:49:22 GMT
//
//

// system include files
#include <memory>
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/METReco/interface/MET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/angle.h"
#include "DataFormats/Common/interface/Association.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

// Geometry
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include "Geometry/Records/interface/RPCRecoGeometryRcd.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/Records/interface/DTRecoGeometryRcd.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"

// ROOT
#include "TTree.h"
#include "TLorentzVector.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class FCPTreeMaker : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
   public:
      explicit FCPTreeMaker(const edm::ParameterSet&);
      ~FCPTreeMaker();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

      // ----------member data ---------------------------

      bool isMC;

      // HLT paths
      const edm::InputTag triggerResultsTag;
      edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken;

      uint8_t HLT_Mu50, HLT_TkMu50, HLT_IsoMu24, HLT_IsoTkMu24;
      uint8_t HLT_MET75_IsoTrk50;
      uint8_t HLT_PFMETNoMu90_PFMHTNoMu90_IDTight, HLT_PFMETNoMu100_PFMHTNoMu100_IDTight, HLT_PFMETNoMu110_PFMHTNoMu110_IDTight, HLT_PFMETNoMu120_PFMHTNoMu120_IDTight;
      uint8_t HLT_PFMET170_NotCleaned, HLT_PFMET170_HBHECleaned, HLT_PFMET170_BeamHaloCleaned, HLT_PFMET170_HBHE_BeamHaloCleaned;

      // MET filters
      const edm::InputTag metFiltersTag;
      edm::EDGetTokenT<edm::TriggerResults> metFiltersToken;

      uint8_t filter_hbhenoise, filter_hbheiso;
      uint8_t filter_goodVertices;
      uint8_t filter_globaltighthalo2016;
      uint8_t filter_ecaltp, filter_eebadsc;

      // Jets
      const edm::InputTag jetsTag;
      edm::EDGetTokenT<std::vector<reco::PFJet> > jetsToken;

      uint32_t nJets;
      std::vector<double> jet_pt, jet_eta, jet_phi;
      std::vector<double> jet_chargedEmEF, jet_chargedHadEF, jet_neutralEmEF, jet_neutralHadEF;
      std::vector<double> jet_muonEF, jet_electronEF;
      std::vector<int> jet_nConstituents;

      // Muons
      const edm::InputTag muonsTag;
      edm::EDGetTokenT<std::vector<reco::Muon> > muonsToken;

      uint32_t nMuons;
      std::vector<double> muon_pt, muon_eta, muon_theta, muon_phi;
      std::vector<bool> muon_isTracker, muon_isStandAlone, muon_isGlobal, muon_isPF;
      std::vector<bool> muon_isLoose, muon_isMedium, muon_isTight;
      std::vector<double> muon_emEt, muon_hadEt, muon_sumPt, muon_PFIso;
      std::vector<double> muon_timeAtIpInOut;

      // Tracks
      const edm::InputTag tracksTag;
      edm::EDGetTokenT<std::vector<reco::Track> > tracksToken;

      uint32_t nTracks;
      std::vector<int> track_index;
      std::vector<double> track_pt, track_eta, track_theta, track_phi;
      std::vector<int> track_nValidHits, track_nLostHits;
      std::vector<double> track_normalizedChi2, track_dxy, track_dz;
      std::vector<int> track_nDeDx, track_nDeDxLower2, track_nPixelHits;
      std::vector<double> track_meanDeDx;
      std::vector<bool> track_isGlobalMuonMatched;
      std::vector<double> track_alphaMax;

      // Primary vertices
      const edm::InputTag primaryVerticesTag;
      edm::EDGetTokenT<std::vector<reco::Vertex> > primaryVerticesToken;

      uint32_t nVertices;

      // PFMET
      const edm::InputTag metTag;
      edm::EDGetTokenT<edm::View<reco::PFMET> > metToken;

      double met, metPhi, metSumEt;
      double HT;

      // dEdx
      const edm::InputTag dedxTag;
      edm::EDGetTokenT<edm::Association<std::vector<reco::DeDxHitInfo> > > dedxToken;

      std::vector<double> dEdx, dEdx_eta, dEdx_phi;
      std::vector<int> dEdx_trackIndex;
      std::vector<bool> dEdx_isPixelBarrel, dEdx_isPixelEndcaps;
      std::vector<bool> dEdx_isTIB, dEdx_isTID, dEdx_isTOB, dEdx_isTEC;
      std::vector<int> dEdx_layerIndex;

      // GEN particles
      const edm::InputTag GenTag;
      edm::EDGetTokenT<std::vector<reco::GenParticle> > GenToken;

      uint32_t nTauprime;
      std::vector<double> Tauprime_pt, Tauprime_eta, Tauprime_phi;
      double MTauTauprime;

      // Tree
      TTree* tree;

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
FCPTreeMaker::FCPTreeMaker(const edm::ParameterSet& iConfig):
   isMC(iConfig.getParameter<bool>("isMC")),
   triggerResultsTag(iConfig.getParameter<edm::InputTag>("triggerResults")),
   metFiltersTag(iConfig.getParameter<edm::InputTag>("metFilters")),
   jetsTag(iConfig.getParameter<edm::InputTag>("jets")),
   muonsTag(iConfig.getParameter<edm::InputTag>("muons")),
   tracksTag(iConfig.getParameter<edm::InputTag>("tracks")),
   primaryVerticesTag(iConfig.getParameter<edm::InputTag>("primaryVertices")),
   metTag(iConfig.getParameter<edm::InputTag>("met")),
   dedxTag(iConfig.getParameter<edm::InputTag>("dedx")),
   GenTag(iConfig.getParameter<edm::InputTag>("genparticles"))
{
   //now do what ever initialization is needed
   usesResource("TFileService");

   triggerResultsToken = consumes<edm::TriggerResults> (triggerResultsTag);
   metFiltersToken = consumes<edm::TriggerResults> (metFiltersTag);
   muonsToken = consumes<std::vector<reco::Muon> > (muonsTag);
   jetsToken = consumes<std::vector<reco::PFJet> > (jetsTag);
   tracksToken = consumes<std::vector<reco::Track> > (tracksTag);
   primaryVerticesToken = consumes<std::vector<reco::Vertex> > (primaryVerticesTag);
   metToken = consumes<edm::View<reco::PFMET> > (metTag);
   dedxToken = consumes<edm::Association<std::vector<reco::DeDxHitInfo> > > (dedxTag);
   if(isMC) GenToken = consumes<std::vector<reco::GenParticle> > (GenTag);

}


FCPTreeMaker::~FCPTreeMaker()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
FCPTreeMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   // HLT
   Handle<TriggerResults> triggerResultsH;
   iEvent.getByToken(triggerResultsToken, triggerResultsH);

   HLT_Mu50 = 0; HLT_TkMu50 = 0; HLT_IsoMu24 = 0; HLT_IsoTkMu24 = 0;
   HLT_MET75_IsoTrk50 = 0;
   HLT_PFMETNoMu90_PFMHTNoMu90_IDTight = 0; HLT_PFMETNoMu100_PFMHTNoMu100_IDTight = 0; HLT_PFMETNoMu110_PFMHTNoMu110_IDTight = 0; HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = 0;
   HLT_PFMET170_NotCleaned = 0; HLT_PFMET170_HBHECleaned = 0; HLT_PFMET170_BeamHaloCleaned = 0; HLT_PFMET170_HBHE_BeamHaloCleaned = 0;

   const edm::TriggerNames &triggernames = iEvent.triggerNames(*triggerResultsH);
   for(unsigned int i=0; i<triggerResultsH->size(); i++){
     if(triggernames.triggerName(i).find("HLT_Mu50")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_Mu50 = 1;
     if(triggernames.triggerName(i).find("HLT_TkMu50")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_TkMu50 = 1;
     if(triggernames.triggerName(i).find("HLT_IsoMu24")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_IsoMu24 = 1;
     if(triggernames.triggerName(i).find("HLT_IsoTkMu24")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_IsoTkMu24 = 1;
     if(triggernames.triggerName(i).find("HLT_MET75_IsoTrk50")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_MET75_IsoTrk50 = 1;
     if(triggernames.triggerName(i).find("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_PFMETNoMu90_PFMHTNoMu90_IDTight = 1;
     if(triggernames.triggerName(i).find("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_PFMETNoMu100_PFMHTNoMu100_IDTight = 1;
     if(triggernames.triggerName(i).find("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_PFMETNoMu110_PFMHTNoMu110_IDTight = 1;
     if(triggernames.triggerName(i).find("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_PFMETNoMu120_PFMHTNoMu120_IDTight = 1;
     if(triggernames.triggerName(i).find("HLT_PFMET170_NotCleaned")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_PFMET170_NotCleaned = 1;
     if(triggernames.triggerName(i).find("HLT_PFMET170_HBHECleaned")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_PFMET170_HBHECleaned = 1;
     if(triggernames.triggerName(i).find("HLT_PFMET170_BeamHaloCleaned")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_PFMET170_BeamHaloCleaned = 1;
     if(triggernames.triggerName(i).find("HLT_PFMET170_HBHE_BeamHaloCleaned")!=std::string::npos && triggerResultsH->accept(i)==1) HLT_PFMET170_HBHE_BeamHaloCleaned = 1;
   }

   // MET filters
   Handle<TriggerResults> metFiltersH;
   iEvent.getByToken(metFiltersToken, metFiltersH);

   filter_hbhenoise = 0; filter_hbheiso = 0;
   filter_goodVertices = 0;
   filter_globaltighthalo2016 = 0;
   filter_ecaltp = 0; filter_eebadsc = 0;

   const edm::TriggerNames &filternames = iEvent.triggerNames(*metFiltersH);
   for(unsigned int i=0; i<metFiltersH->size(); i++){
     if(filternames.triggerName(i).find("Flag_HBHENoiseFilter")!=std::string::npos && metFiltersH->accept(i)==1) filter_hbhenoise = 1;
     if(filternames.triggerName(i).find("Flag_HBHENoiseIsoFilter")!=std::string::npos && metFiltersH->accept(i)==1) filter_hbheiso = 1;
     if(filternames.triggerName(i).find("Flag_goodVertices")!=std::string::npos && metFiltersH->accept(i)==1) filter_goodVertices = 1;
     if(filternames.triggerName(i).find("Flag_globalTightHalo2016Filter")!=std::string::npos && metFiltersH->accept(i)==1) filter_globaltighthalo2016 = 1;
     if(filternames.triggerName(i).find("Flag_EcalDeadCellTriggerPrimitiveFilter")!=std::string::npos && metFiltersH->accept(i)==1) filter_ecaltp = 1;
     if(filternames.triggerName(i).find("Flag_eeBadScFilter")!=std::string::npos && metFiltersH->accept(i)==1) filter_eebadsc = 1;
   }

   // Primary vertices
   Handle<std::vector<reco::Vertex> > primaryVerticesH;
   iEvent.getByToken(primaryVerticesToken, primaryVerticesH);

   nVertices = 0;

   for(const reco::Vertex &vertex : *primaryVerticesH){
     double R = sqrt(pow(vertex.x(),2)+pow(vertex.y(),2));
     if(!vertex.isFake() && vertex.ndof()>4 && fabs(vertex.z())<24 && R<2) nVertices++;
   }

   // Muons
   Handle<std::vector<reco::Muon> > muonsH;
   iEvent.getByToken(muonsToken, muonsH);

   muon_pt.clear(); muon_eta.clear(); muon_phi.clear(); muon_theta.clear();
   muon_isTracker.clear(); muon_isStandAlone.clear(); muon_isGlobal.clear(); muon_isPF.clear();
   muon_isLoose.clear(); muon_isMedium.clear(); muon_isTight.clear();
   muon_emEt.clear(); muon_hadEt.clear(); muon_sumPt.clear(); muon_PFIso.clear();
   muon_timeAtIpInOut.clear();

   for(const reco::Muon &muon : *muonsH){
     if(muon.pt()<20) continue;
     muon_pt.push_back(muon.pt());
     muon_eta.push_back(muon.eta());
     muon_phi.push_back(muon.phi());
     muon_theta.push_back(muon.theta());
     muon_isTracker.push_back(muon.isTrackerMuon());
     muon_isStandAlone.push_back(muon.isStandAloneMuon());
     muon_isGlobal.push_back(muon.isGlobalMuon());
     muon_isPF.push_back(muon.isPFMuon());
     muon_isLoose.push_back(muon::isLooseMuon(muon));
     muon_isMedium.push_back(muon::isMediumMuon(muon));
     // Find muon closest PV - needed for the tight muon identification
     double d3Dmin = 1;
     bool isTightMuon_tmp = false;
     for(const reco::Vertex &vertex : *primaryVerticesH){
       double R = sqrt(pow(vertex.x(),2)+pow(vertex.y(),2));
       if(vertex.isFake() || vertex.ndof()<=4 || fabs(vertex.z())>24 || R>2 || !muon.globalTrack()) continue;
       double d3D = sqrt(pow(muon.globalTrack()->vx()-vertex.x(),2)+pow(muon.globalTrack()->vy()-vertex.y(),2)+pow(muon.globalTrack()->vz()-vertex.z(),2));
       if(d3D<d3Dmin){
         d3Dmin = d3D;
         const reco::Vertex &muonVertex = vertex;
         isTightMuon_tmp = muon::isTightMuon(muon,muonVertex);
       }
     }
     if(isTightMuon_tmp) muon_isTight.push_back(true);
     else muon_isTight.push_back(false);
     muon_emEt.push_back(muon.isolationR03().emEt);
     muon_hadEt.push_back(muon.isolationR03().hadEt);
     muon_sumPt.push_back(muon.isolationR03().sumPt);
     muon_PFIso.push_back((muon.pfIsolationR04().sumChargedHadronPt + std::max(0., muon.pfIsolationR04().sumNeutralHadronEt + muon.pfIsolationR04().sumPhotonEt - 0.5*muon.pfIsolationR04().sumPUPt))/muon.pt());
     muon_timeAtIpInOut.push_back(muon.time().timeAtIpInOut);
   }
   nMuons = muon_pt.size();

   // Jets
   Handle<std::vector<reco::PFJet> > jetsH;
   iEvent.getByToken(jetsToken, jetsH);

   jet_pt.clear(); jet_eta.clear(); jet_phi.clear();
   jet_chargedEmEF.clear(); jet_chargedHadEF.clear(); jet_neutralEmEF.clear(); jet_neutralHadEF.clear();
   jet_electronEF.clear(); jet_muonEF.clear();
   jet_nConstituents.clear();
   HT = 0;

   for(const reco::PFJet &jet : *jetsH){
     bool isLooseJet = false;
     if(fabs(jet.eta())<=2.7){
       if((jet.neutralHadronEnergyFraction()<0.99 && jet.neutralEmEnergyFraction()<0.99 && jet.nConstituents()>1) &&
          ((fabs(jet.eta())<=2.4 && jet.chargedHadronEnergyFraction()>0 && jet.chargedEmEnergyFraction()<0.99 && jet.chargedMultiplicity()>0) || fabs(jet.eta())>2.4)) isLooseJet = true;
     }
     else if(fabs(jet.eta())>2.7 && fabs(jet.eta())<=3){
       if(jet.neutralHadronEnergyFraction()<0.98 && jet.neutralEmEnergyFraction()>0.01 && jet.neutralMultiplicity()>2) isLooseJet = true;
     }
     else{
       if(jet.neutralEmEnergyFraction()<0.9 && jet.neutralMultiplicity()>10) isLooseJet = true;
     }
     bool isTightMuonMatched = false;
     for(size_t imuon=0; imuon<muon_eta.size(); imuon++){
       if(!muon_isTight.at(imuon)) continue;
       if(deltaR(muon_eta.at(imuon),muon_phi.at(imuon),jet.eta(),jet.phi())<0.4) isTightMuonMatched = true;
     }
     if(!isLooseJet || isTightMuonMatched || jet.pt()<20) continue;
     HT += jet.pt();
     jet_pt.push_back(jet.pt());
     jet_eta.push_back(jet.eta());
     jet_phi.push_back(jet.phi());
     jet_chargedEmEF.push_back(jet.chargedEmEnergyFraction());
     jet_chargedHadEF.push_back(jet.chargedHadronEnergyFraction());
     jet_neutralEmEF.push_back(jet.neutralEmEnergyFraction());
     jet_neutralHadEF.push_back(jet.neutralHadronEnergyFraction());
     jet_electronEF.push_back(jet.electronEnergyFraction());
     jet_muonEF.push_back(jet.muonEnergyFraction());
     jet_nConstituents.push_back(jet.nConstituents());
   }
   nJets = jet_pt.size();

   // Tracks and hits
   Handle<std::vector<reco::Track> > tracksH;
   iEvent.getByToken(tracksToken, tracksH);

   ESHandle<TrackerGeometry> theTrackerGeometry;
   iSetup.get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);

   //Retrieve tracker topology from geometry
   edm::ESHandle<TrackerTopology> tTopoHandle;
   iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
   const TrackerTopology* const tTopo = tTopoHandle.product();

   track_index.clear();
   track_pt.clear(); track_eta.clear(); track_phi.clear(); track_theta.clear();
   track_nValidHits.clear(); track_nLostHits.clear(); track_normalizedChi2.clear(); track_dxy.clear(); track_dz.clear();
   track_nDeDx.clear(); track_nDeDxLower2.clear(); track_meanDeDx.clear(); track_nPixelHits.clear();
   track_isGlobalMuonMatched.clear();
   track_alphaMax.clear();

   // dedx
   Handle<edm::Association<std::vector<reco::DeDxHitInfo> > > dedxH;
   iEvent.getByToken(dedxToken, dedxH);

   dEdx.clear(); dEdx_trackIndex.clear(); dEdx_eta.clear(); dEdx_phi.clear();
   dEdx_isPixelBarrel.clear(); dEdx_isPixelEndcaps.clear();
   dEdx_isTIB.clear(); dEdx_isTID.clear(); dEdx_isTOB.clear(); dEdx_isTEC.clear();
   dEdx_layerIndex.clear();

   unsigned int iref = 0, itrack = 0;
   for(const reco::Track &track : *tracksH){
     reco::TrackRef trackref = reco::TrackRef(tracksH.product(),iref);
     iref++;

     // Skip track if it doesn't exist or if it has low pT
     if(trackref.isNull() || track.pt()<20) continue;
     // Define associated dE/dx hits collection
     const reco::DeDxHitInfo *dedx = NULL;
     reco::DeDxHitInfoRef dedxref = dedxH->get(trackref.key());
     if(!dedxref.isNull()) dedx = &(*dedxref);
     // Skip track if no hits are associated
     if(!dedx) continue;
     // Compute dxy and dz w.r.t the closest primary vertex
     double d3Dmin = 1;
     double dxy = 1, dz = 1;
     for(const reco::Vertex &vertex : *primaryVerticesH){
       double R = sqrt(pow(vertex.x(),2)+pow(vertex.y(),2));
       if(vertex.isFake() || vertex.ndof()<=4 || fabs(vertex.z())>24 || R>2) continue;
       double d3D = sqrt(pow(track.vx()-vertex.x(),2)+pow(track.vy()-vertex.y(),2)+pow(track.vz()-vertex.z(),2));
       if(d3D<d3Dmin){
         d3Dmin = d3D;
         dxy = sqrt(pow(track.vx()-vertex.x(),2)+pow(track.vy()-vertex.y(),2));
         dz = fabs(track.vz()-vertex.z());
       }
     }
     // Skip if no close good PV was found
     if(d3Dmin==1) continue;
     // Fill with standard variables
     track_index.push_back(itrack);
     track_pt.push_back(track.pt());
     track_eta.push_back(track.eta());
     track_phi.push_back(track.phi());
     track_theta.push_back(track.theta());
     track_nValidHits.push_back(track.numberOfValidHits());
     track_nLostHits.push_back(track.numberOfLostHits());
     track_normalizedChi2.push_back(track.normalizedChi2());
     track_dxy.push_back(dxy);
     track_dz.push_back(dz);

     // Compute max 3D angle between current track and any other pT>35 GeV track
     double track_alphaMax_tmp = 0;
     for(const reco::Track &track2 : *tracksH){
       if(track2.pt()<35) continue;
       if(angle(track.px(),track.py(),track.pz(),track2.px(),track2.py(),track2.pz())>track_alphaMax_tmp) track_alphaMax_tmp = angle(track.px(),track.py(),track.pz(),track2.px(),track2.py(),track2.pz());
     }

     double deltaRmin = 1;
     double deltaR_tmp = 0;
     for(const reco::Muon &muon : *muonsH){
       if(!muon.isGlobalMuon()) continue;

       // Compute deltaR(track, global muons) to find the best matching global muon
       deltaR_tmp = deltaR(track.eta(),track.phi(),muon.eta(),muon.phi());
       if(deltaR_tmp<0.02 && deltaR_tmp<deltaRmin) deltaRmin = deltaR_tmp;

       // Compute max 3D angle between current track and any pT>35 GeV standalone muon
       if(muon.pt()>35 && muon.isStandAloneMuon() && angle(track.px(),track.py(),track.pz(),muon.px(),muon.py(),muon.pz())>track_alphaMax_tmp) track_alphaMax_tmp = angle(track.px(),track.py(),track.pz(),muon.px(),muon.py(),muon.pz());
     }
     if(deltaRmin<1) track_isGlobalMuonMatched.push_back(true);
     else track_isGlobalMuonMatched.push_back(false);

     track_alphaMax.push_back(track_alphaMax_tmp);

     // Fill dE/dx related variables for hits associated to the current track
     int nDeDxLower2 = 0;
     int nPixelHits = 0;
     double dE_tot = 0;
     double dx_tot = 0;
     for(size_t i=0; i<dedx->size(); i++){
       dEdx_trackIndex.push_back(itrack);
       const DetId dedxId = dedx->detId(i);
       const auto dedxPos = theTrackerGeometry->idToDet(dedxId)->position();
       dEdx_eta.push_back(dedxPos.eta());
       dEdx_phi.push_back(dedxPos.phi());
       double charge;
       if(dedxId.subdetId()<3){
         dEdx_isTIB.push_back(false); dEdx_isTID.push_back(false); dEdx_isTOB.push_back(false); dEdx_isTEC.push_back(false);
         if(dedxId.subdetId()==1){
           dEdx_isPixelBarrel.push_back(true); dEdx_isPixelEndcaps.push_back(false);
           PXBDetId PixelBarrelDetId(dedxId);
           dEdx_layerIndex.push_back(PixelBarrelDetId.layer());
         }
         else if(dedxId.subdetId()==2){
           dEdx_isPixelBarrel.push_back(false); dEdx_isPixelEndcaps.push_back(true);
           PXFDetId PixelEndcapsDetId(dedxId);
           dEdx_layerIndex.push_back(PixelEndcapsDetId.disk());
         }
         nPixelHits++;
         charge = 3.61e-06*dedx->charge(i);
       }
       else{
         dEdx_isPixelBarrel.push_back(false); dEdx_isPixelEndcaps.push_back(false);
         SiStripDetId TrackerDetId(dedxId);
         if(TrackerDetId.subDetector()==3){ 
           dEdx_isTIB.push_back(true); dEdx_isTID.push_back(false); dEdx_isTOB.push_back(false); dEdx_isTEC.push_back(false);
           dEdx_layerIndex.push_back(tTopo->tibLayer(TrackerDetId));
         }
         else if(TrackerDetId.subDetector()==4){
           dEdx_isTIB.push_back(false); dEdx_isTID.push_back(true); dEdx_isTOB.push_back(false); dEdx_isTEC.push_back(false);
           dEdx_layerIndex.push_back(tTopo->tidWheel(TrackerDetId));
         }
         else if(TrackerDetId.subDetector()==5){
           dEdx_isTIB.push_back(false); dEdx_isTID.push_back(false); dEdx_isTOB.push_back(true); dEdx_isTEC.push_back(false);
           dEdx_layerIndex.push_back(tTopo->tobLayer(TrackerDetId));
         }
         else if(TrackerDetId.subDetector()==6){
           dEdx_isTIB.push_back(false); dEdx_isTID.push_back(false); dEdx_isTOB.push_back(false); dEdx_isTEC.push_back(true);
           dEdx_layerIndex.push_back(tTopo->tecWheel(TrackerDetId));
         }
         charge = 3.61e-06*265*dedx->charge(i);
       }
       dEdx.push_back(charge/(dedx->pathlength(i)));
       if((charge/(dedx->pathlength(i)))<2) nDeDxLower2++;
       dE_tot += charge;
       dx_tot += dedx->pathlength(i);
     }
     track_nDeDx.push_back(dedx->size());
     track_nDeDxLower2.push_back(nDeDxLower2);
     track_nPixelHits.push_back(nPixelHits);
     track_meanDeDx.push_back(dE_tot/dx_tot);
     itrack++;
   }
   nTracks = track_pt.size();

   Handle<View<reco::PFMET> > metH;
   iEvent.getByToken(metToken, metH);

   const reco::PFMET &recoMET = metH->front();
   met = recoMET.pt();
   metPhi = recoMET.phi();
   metSumEt = recoMET.sumEt();

   // Gen particles
   if(isMC){
     Handle<std::vector<reco::GenParticle> > GenH;
     iEvent.getByToken(GenToken, GenH);

     TLorentzVector v_Tau, v_Tauprime;

     Tauprime_pt.clear(); Tauprime_eta.clear(); Tauprime_phi.clear();

     for(const reco::GenParticle &gen : *GenH){
       if(gen.pdgId()==17){
         v_Tau.SetPtEtaPhiE(gen.pt(),gen.eta(),gen.phi(),gen.energy());
         Tauprime_pt.push_back(gen.pt());
         Tauprime_eta.push_back(gen.eta());
         Tauprime_phi.push_back(gen.phi());
       }
       else if(gen.pdgId()==-17){
         v_Tauprime.SetPtEtaPhiE(gen.pt(),gen.eta(),gen.phi(),gen.energy());
         Tauprime_pt.push_back(gen.pt());
         Tauprime_eta.push_back(gen.eta());
         Tauprime_phi.push_back(gen.phi());
       }
     }
     nTauprime = Tauprime_pt.size();
     MTauTauprime = (v_Tau+v_Tauprime).M();
   }

   // Filling the tree
   if(nTracks>0 && HLT_Mu50>0) tree->Fill();

}


// ------------ method called once each job just before starting event loop  ------------
void 
FCPTreeMaker::beginJob()
{
   std::cout << "isMC = " << isMC << std::endl;

   edm::Service<TFileService> fs;
   tree = fs->make<TTree>("tree", "tree");

   tree->Branch("HLT_Mu50", &HLT_Mu50, "HLT_Mu50/b");
   tree->Branch("HLT_TkMu50", &HLT_TkMu50, "HLT_TkMu50/b");
   tree->Branch("HLT_IsoMu24", &HLT_IsoMu24, "HLT_IsoMu24/b");
   tree->Branch("HLT_IsoTkMu24", &HLT_IsoTkMu24, "HLT_IsoTkMu24/b");
   tree->Branch("HLT_MET75_IsoTrk50", &HLT_MET75_IsoTrk50, "HLT_MET75_IsoTrk50/b");
   tree->Branch("HLT_PFMETNoMu90_PFMHTNoMu90_IDTight", &HLT_PFMETNoMu90_PFMHTNoMu90_IDTight, "HLT_PFMETNoMu90_PFMHTNoMu90_IDTight/b");
   tree->Branch("HLT_PFMETNoMu100_PFMHTNoMu100_IDTight", &HLT_PFMETNoMu100_PFMHTNoMu100_IDTight, "HLT_PFMETNoMu100_PFMHTNoMu100_IDTight/b");
   tree->Branch("HLT_PFMETNoMu110_PFMHTNoMu110_IDTight", &HLT_PFMETNoMu110_PFMHTNoMu110_IDTight, "HLT_PFMETNoMu110_PFMHTNoMu110_IDTight/b");
   tree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight", &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight, "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight/b");
   tree->Branch("HLT_PFMET170_NotCleaned", &HLT_PFMET170_NotCleaned, "HLT_PFMET170_NotCleaned/b");
   tree->Branch("HLT_PFMET170_HBHECleaned", &HLT_PFMET170_HBHECleaned, "HLT_PFMET170_HBHECleaned/b");
   tree->Branch("HLT_PFMET170_BeamHaloCleaned", &HLT_PFMET170_BeamHaloCleaned, "HLT_PFMET170_BeamHaloCleaned/b");
   tree->Branch("HLT_PFMET170_HBHE_BeamHaloCleaned", &HLT_PFMET170_HBHE_BeamHaloCleaned, "HLT_PFMET170_HBHE_BeamHaloCleaned/b");

   tree->Branch("nJets", &nJets, "nJets/i");
   tree->Branch("jet_pt","std::vector<double>", &jet_pt);
   tree->Branch("jet_eta","std::vector<double>", &jet_eta);
   tree->Branch("jet_phi","std::vector<double>", &jet_phi);
   tree->Branch("jet_chargedEmEF","std::vector<double>", &jet_chargedEmEF);
   tree->Branch("jet_chargedHadEF","std::vector<double>", &jet_chargedHadEF);
   tree->Branch("jet_neutralEmEF","std::vector<double>", &jet_neutralEmEF);
   tree->Branch("jet_neutralHadEF","std::vector<double>", &jet_neutralHadEF);
   tree->Branch("jet_muonEF","std::vector<double>", &jet_muonEF);
   tree->Branch("jet_electronEF","std::vector<double>", &jet_electronEF);
   tree->Branch("jet_nConstituents","std::vector<int>", &jet_nConstituents);

   tree->Branch("nMuons", &nMuons, "nMuons/i");
   tree->Branch("muon_pt","std::vector<double>", &muon_pt);
   tree->Branch("muon_eta","std::vector<double>", &muon_eta);
   tree->Branch("muon_phi","std::vector<double>", &muon_phi);
   tree->Branch("muon_theta","std::vector<double>", &muon_theta);
   tree->Branch("muon_isTracker", "std::vector<bool>", &muon_isTracker);
   tree->Branch("muon_isStandAlone", "std::vector<bool>", &muon_isStandAlone);
   tree->Branch("muon_isGlobal", "std::vector<bool>", &muon_isGlobal);
   tree->Branch("muon_isPF", "std::vector<bool>", &muon_isPF);
   tree->Branch("muon_isLoose", "std::vector<bool>", &muon_isLoose);
   tree->Branch("muon_isMedium", "std::vector<bool>", &muon_isMedium);
   tree->Branch("muon_isTight", "std::vector<bool>", &muon_isTight);
   tree->Branch("muon_emEt", "std::vector<double>", &muon_emEt);
   tree->Branch("muon_hadEt", "std::vector<double>", &muon_hadEt);
   tree->Branch("muon_sumPt", "std::vector<double>", &muon_sumPt);
   tree->Branch("muon_PFIso", "std::vector<double>", &muon_PFIso);
   tree->Branch("muon_timeAtIpInOut", "std::vector<double>", &muon_timeAtIpInOut);

   tree->Branch("nTracks", &nTracks, "nTracks/i");
   tree->Branch("track_index","std::vector<int>", &track_index);
   tree->Branch("track_pt","std::vector<double>", &track_pt);
   tree->Branch("track_eta","std::vector<double>", &track_eta);
   tree->Branch("track_phi","std::vector<double>", &track_phi);
   tree->Branch("track_nValidHits","std::vector<int>", &track_nValidHits);
   tree->Branch("track_nLostHits","std::vector<int>", &track_nLostHits);
   tree->Branch("track_normalizedChi2","std::vector<double>", &track_normalizedChi2);
   tree->Branch("track_dxy","std::vector<double>", &track_dxy);
   tree->Branch("track_dz","std::vector<double>", &track_dz);
   tree->Branch("track_nDeDx","std::vector<int>", &track_nDeDx);
   tree->Branch("track_nDeDxLower2","std::vector<int>", &track_nDeDxLower2);
   tree->Branch("track_nPixelHits","std::vector<int>",&track_nPixelHits);
   tree->Branch("track_meanDeDx","std::vector<double>", &track_meanDeDx);
   tree->Branch("track_isGlobalMuonMatched", "std::vector<bool>", &track_isGlobalMuonMatched);
   tree->Branch("track_alphaMax", "std::vector<double>", &track_alphaMax);

   tree->Branch("nVertices", &nVertices, "nVertices/i");

   tree->Branch("dEdx", "std::vector<double>", &dEdx);
   tree->Branch("dEdx_eta", "std::vector<double>", &dEdx_eta);
   tree->Branch("dEdx_phi", "std::vector<double>", &dEdx_phi);
   tree->Branch("dEdx_trackIndex", "std::vector<int>", &dEdx_trackIndex);
   tree->Branch("dEdx_isPixelBarrel", "std::vector<bool>", &dEdx_isPixelBarrel);
   tree->Branch("dEdx_isPixelEndcaps", "std::vector<bool>", &dEdx_isPixelEndcaps);
   tree->Branch("dEdx_isTIB", "std::vector<bool>", &dEdx_isTIB);
   tree->Branch("dEdx_isTID", "std::vector<bool>", &dEdx_isTID);
   tree->Branch("dEdx_isTOB", "std::vector<bool>", &dEdx_isTOB);
   tree->Branch("dEdx_isTEC", "std::vector<bool>", &dEdx_isTEC);
   tree->Branch("dEdx_layerIndex", "std::vector<int>", &dEdx_layerIndex);

   tree->Branch("met", &met, "met/D");
   tree->Branch("metPhi", &metPhi, "metPhi/D");
   tree->Branch("metSumEt", &metSumEt, "metSumEt/D");
   tree->Branch("HT", &HT, "HT/D");

   tree->Branch("filter_hbhenoise", &filter_hbhenoise, "filter_hbhenoise/b");
   tree->Branch("filter_hbheiso", &filter_hbheiso, "filter_hbheiso/b");
   tree->Branch("filter_goodVertices", &filter_goodVertices, "filter_goodVertices/b");
   tree->Branch("filter_globaltighthalo2016", &filter_globaltighthalo2016, "filter_globaltighthalo2016/b");
   tree->Branch("filter_ecaltp", &filter_ecaltp, "filter_ecaltp/b");
   tree->Branch("filter_eebadsc", &filter_eebadsc, "filter_eebadsc/b");

   if(isMC){
     tree->Branch("nTauprime", &nTauprime, "nTauprime/i");
     tree->Branch("Tauprime_pt","std::vector<double>", &Tauprime_pt);
     tree->Branch("Tauprime_eta","std::vector<double>", &Tauprime_eta);
     tree->Branch("Tauprime_phi","std::vector<double>", &Tauprime_phi);
     tree->Branch("MTauTauprime", &MTauTauprime, "MTauTauprime/D");
   }
}

// ------------ method called once each job just after ending the event loop  ------------
void 
FCPTreeMaker::endJob() 
{
}



// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
FCPTreeMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(FCPTreeMaker);
