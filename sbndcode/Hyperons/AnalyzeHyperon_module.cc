#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/FindOneP.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

// LArSoft Includes
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/PFParticleMetadata.h"
#include "lardataobj/RecoBase/Slice.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/Shower.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbndcode/CosmicId/Algs/PandoraNuScoreCosmicIdAlg.h"
#include "larsim/Utils/TruthMatchUtils.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

// Root Includes
#include "TLorentzVector.h"
#include <iostream>
#include <vector>
#include <set>
#include <map>
#include <TTree.h>
#include <TH1.h>
#include <fstream>
#include <string>
#include <exception>
#include <limits> 

namespace hyperon {

    class AnalyzeHyperon : public art::EDAnalyzer {

    public:
        explicit AnalyzeHyperon(fhicl::ParameterSet const& p);
        AnalyzeHyperon(AnalyzeHyperon const&) = delete;
        AnalyzeHyperon(AnalyzeHyperon&&) = delete;

        void beginJob() override;
        void analyze(art::Event const& e) override;
        void endSubRun(art::SubRun const& sr) override;

    private:
        // The ROOT Tree(s)
        TTree* fTree;
        TTree* fSubRunTree;

        // SubRun Info
        int fRun_sr;
        int fSubRun_sr;
        double fPOT;

        // Event info
        int fRun;
        int fEvent;

        // Labels
        std::string fSliceLabel = "pandora";
        std::string fPFParticleLabel = "pandora";
        std::string fTrackLabel = "pandoraTrack";
        std::string fShowerLabel = "pandoraShower";
        std::string fCrumbsLabel = "crumbs";
        std::string fHitLabel = "gaushit";
        std::string fPIDLabel = "pandoraPid";
        std::string fOpT0Label = "opt0finder";

        std::vector<int> fNu_PDG;
        std::vector<int> fNu_Mode;
        std::vector<int> fNu_NuMode;
        std::vector<double> fNu_E;
        std::vector<int> fNu_TrackID;

        std::vector<int> particle_interaction_index;

        std::vector<int> neutrino_daughters;
        std::vector<int> particlePDG;
        std::vector<int> status_code;
        std::vector<int> particleMotherID;
        std::vector<int> particleTrackID;
        std::vector<int> particleMotherPDG;
        std::vector<double> particleKE;
        std::vector<double> particleP;
        std::vector<double> particleE;
        std::vector<double> particlePx;
        std::vector<double> particlePy;
        std::vector<double> particlePz;

        // True Geant4 Particle Blocks
        std::vector<int> geant_PDG;
        std::vector<int> geant_TrackID;
        std::vector<int> geant_MotherID;
        std::vector<double> geant_startX;
        std::vector<double> geant_startY;
        std::vector<double> geant_startZ;
        std::vector<double> geant_endX;
        std::vector<double> geant_endY;
        std::vector<double> geant_endZ;
        std::vector<double> geant_P;
        std::vector<double> geant_Px;
        std::vector<double> geant_Py;
        std::vector<double> geant_Pz;
        std::vector<int> geant_TrueHitCount;
        std::vector<int> geant_HitsU;
        std::vector<int> geant_HitsV;
        std::vector<int> geant_HitsW;

        // True Interaction Vertex
        std::vector<double> fTrueVtxX;
        std::vector<double> fTrueVtxY;
        std::vector<double> fTrueVtxZ;
        std::vector<double> fTrueVtxT;
        

        // Reconstructed Slices
        std::vector<int> fSlice_ID;
        std::vector<float> fSlice_NuScore;
        std::vector<int>   fSlice_TotalHits;
        std::vector<int>   fSlice_TrueNuHits;
        int                fEvent_TotalTrueNuHits; 


        std::vector<float> fSlice_VtxX;
        std::vector<float> fSlice_VtxY;
        std::vector<float> fSlice_VtxZ;

        std::vector<float> fSlice_Opt0Score;

        // Reconstructed Tracks
        std::vector<int> fTrack_ID;
        std::vector<int> fTrack_SliceID;
        std::vector<int> fTrack_TrueGeantID;
        std::vector<int> fTrack_IsPrimary;

        // Reconstructed Track Geometry s
        std::vector<float> fTrack_Length;
        std::vector<float> fTrack_StartX;
        std::vector<float> fTrack_StartY;
        std::vector<float> fTrack_StartZ;
        std::vector<float> fTrack_EndX;
        std::vector<float> fTrack_EndY;
        std::vector<float> fTrack_EndZ;

        // Reconstructed Showers
        std::vector<int> fShower_ID;
        std::vector<int> fShower_SliceID;
        std::vector<int> fShower_TrueGeantID;

        std::vector<float> fTrack_TrackScore;
        std::vector<float> fShower_TrackScore;
        std::vector<int> fShower_IsPrimary;
        std::vector<float> fTrack_Chi2Proton;
        std::vector<float> fTrack_Chi2Muon;



    };

    AnalyzeHyperon::AnalyzeHyperon(fhicl::ParameterSet const& p)
        : EDAnalyzer(p) {}

    void AnalyzeHyperon::beginJob() {
        art::ServiceHandle<art::TFileService> tfs;
        fTree = tfs->make<TTree>("tree", "Hyperon Slicing Performance Tree");

        fTree->Branch("run", &fRun, "run/I");
        fTree->Branch("event", &fEvent, "event/I");

        fTree->Branch("nu_pdg", &fNu_PDG);
        fTree->Branch("nu_mode", &fNu_Mode);
        fTree->Branch("nu_numode", &fNu_NuMode);
        fTree->Branch("nu_energy", &fNu_E);
        fTree->Branch("nu_track_id", &fNu_TrackID);
        fTree->Branch("NeutrinoDaughters", &neutrino_daughters);

        fTree->Branch("particle_interaction_index", &particle_interaction_index);

        fTree->Branch("TrueVtxX", &fTrueVtxX);
        fTree->Branch("TrueVtxY", &fTrueVtxY);
        fTree->Branch("TrueVtxZ", &fTrueVtxZ);
        fTree->Branch("TrueVtxT", &fTrueVtxT);

        fTree->Branch("particlePDG", &particlePDG);
        fTree->Branch("status_code", &status_code);
        fTree->Branch("particleMotherID", &particleMotherID);
        fTree->Branch("particleTrackID", &particleTrackID);
        fTree->Branch("particleMotherPDG", &particleMotherPDG);
        fTree->Branch("particleKE", &particleKE);
        fTree->Branch("particleP", &particleP);
        fTree->Branch("particleE", &particleE);
        fTree->Branch("particlePx", &particlePx);
        fTree->Branch("particlePy", &particlePy);
        fTree->Branch("particlePz", &particlePz);

        // Geant4 Identity & Coordinates
        fTree->Branch("geant_PDG", &geant_PDG);
        fTree->Branch("geant_TrackID", &geant_TrackID);
        fTree->Branch("geant_MotherID", &geant_MotherID);
        fTree->Branch("geant_startX", &geant_startX);
        fTree->Branch("geant_startY", &geant_startY);
        fTree->Branch("geant_startZ", &geant_startZ);
        fTree->Branch("geant_endX", &geant_endX);
        fTree->Branch("geant_endY", &geant_endY);
        fTree->Branch("geant_endZ", &geant_endZ);
        fTree->Branch("geant_P", &geant_P);
        fTree->Branch("geant_Px", &geant_Px);
        fTree->Branch("geant_Py", &geant_Py);
        fTree->Branch("geant_Pz", &geant_Pz);
        fTree->Branch("geant_TrueHitCount", &geant_TrueHitCount);
        fTree->Branch("geant_HitsU", &geant_HitsU);
        fTree->Branch("geant_HitsV", &geant_HitsV);
        fTree->Branch("geant_HitsW", &geant_HitsW);

        // Flat Slices
        fTree->Branch("slice_ID", &fSlice_ID);
        fTree->Branch("slice_NuScore", &fSlice_NuScore);
        fTree->Branch("slice_TotalHits", &fSlice_TotalHits);
        fTree->Branch("slice_TrueNuHits", &fSlice_TrueNuHits);
        fTree->Branch("event_TotalTrueNuHits", &fEvent_TotalTrueNuHits, "event_TotalTrueNuHits/I");

        fTree->Branch("slice_VtxX", &fSlice_VtxX);
        fTree->Branch("slice_VtxY", &fSlice_VtxY);
        fTree->Branch("slice_VtxZ", &fSlice_VtxZ);

        fTree->Branch("slice_Opt0Score", &fSlice_Opt0Score);

        // Flat Tracks
        fTree->Branch("track_ID", &fTrack_ID);
        fTree->Branch("track_SliceID", &fTrack_SliceID);
        fTree->Branch("track_TrueGeantID", &fTrack_TrueGeantID);
        fTree->Branch("track_IsPrimary", &fTrack_IsPrimary);

        fTree->Branch("track_Length", &fTrack_Length);
        fTree->Branch("track_StartX", &fTrack_StartX);
        fTree->Branch("track_StartY", &fTrack_StartY);
        fTree->Branch("track_StartZ", &fTrack_StartZ);
        fTree->Branch("track_EndX", &fTrack_EndX);
        fTree->Branch("track_EndY", &fTrack_EndY);
        fTree->Branch("track_EndZ", &fTrack_EndZ);

        // Flat Showers
        fTree->Branch("shower_ID", &fShower_ID);
        fTree->Branch("shower_SliceID", &fShower_SliceID);
        fTree->Branch("shower_TrueGeantID", &fShower_TrueGeantID);
        fTree->Branch("shower_IsPrimary", &fShower_IsPrimary);

        fTree->Branch("track_TrackScore", &fTrack_TrackScore);
        fTree->Branch("shower_TrackScore", &fShower_TrackScore);

        fTree->Branch("track_Chi2Proton", &fTrack_Chi2Proton);
        fTree->Branch("track_Chi2Muon", &fTrack_Chi2Muon);

        // Subrun Metadata Tree
        fSubRunTree = tfs->make<TTree>("subrunTree", "SubRun Level Info");
        fSubRunTree->Branch("run", &fRun_sr, "run/I");
        fSubRunTree->Branch("subRun", &fSubRun_sr, "subRun/I");
        fSubRunTree->Branch("pot", &fPOT, "pot/D");
    }

    void AnalyzeHyperon::analyze(art::Event const& e) {

        art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
        art::ServiceHandle<cheat::BackTrackerService> bt_serv;

        // Reset Generator & Event Info
        fNu_PDG.clear();
        fNu_TrackID.clear();
        fNu_Mode.clear();
        fNu_NuMode.clear();
        fNu_E.clear();
        fTrueVtxX.clear();
        fTrueVtxY.clear();
        fTrueVtxZ.clear();
        fTrueVtxT.clear();

        particle_interaction_index.clear();

        neutrino_daughters.clear();
        particlePDG.clear();
        status_code.clear();
        particleMotherID.clear();
        particleTrackID.clear();
        particleMotherPDG.clear();
        particleKE.clear();
        particleP.clear();
        particleE.clear();
        particlePx.clear();
        particlePy.clear();
        particlePz.clear();

        // Reset Geant4 Data
        geant_PDG.clear();
        geant_TrackID.clear();
        geant_MotherID.clear();
        geant_startX.clear();
        geant_startY.clear();
        geant_startZ.clear();
        geant_endX.clear();
        geant_endY.clear();
        geant_endZ.clear();
        geant_P.clear();
        geant_Px.clear();
        geant_Py.clear();
        geant_Pz.clear();
        geant_TrueHitCount.clear();
        geant_HitsU.clear();
        geant_HitsV.clear();
        geant_HitsW.clear();

        // Reset Slices & Tracks & Showers
        fSlice_ID.clear();
        fSlice_NuScore.clear();
        fSlice_TotalHits.clear();
        fSlice_TrueNuHits.clear();
        fEvent_TotalTrueNuHits = 0;

        fTrack_ID.clear();
        fTrack_SliceID.clear();
        fTrack_TrueGeantID.clear();
        fTrack_IsPrimary.clear();

        fShower_ID.clear();
        fShower_SliceID.clear();
        fShower_TrueGeantID.clear();
        fShower_IsPrimary.clear();

        fTrack_TrackScore.clear();
        fShower_TrackScore.clear();

        fSlice_Opt0Score.clear();

        fTrack_Chi2Proton.clear();
        fTrack_Chi2Muon.clear();

        fRun = e.run();
        fEvent = e.event();



        // Set to keep track of Geant4 particles originating from our neutrino
        std::set<int> nu_g4_track_ids;
        // target beam neutrino interaction.
        art::Ptr<simb::MCTruth> beamMCTruth;

        /* ********************************************************
        Truth Block
        ******************************************************** */ 
        art::Handle<std::vector<simb::MCTruth>> mctruthListHandle;
        if (e.getByLabel("generator", mctruthListHandle) && !mctruthListHandle->empty()) {
            art::FindManyP<simb::MCParticle, sim::GeneratedParticleInfo> geantAssns(
                mctruthListHandle, e, "largeant"
            );
            int truth_iter = 0;
            int interaction_idx = 0;
            for(size_t i = 0; i < mctruthListHandle->size(); ++i){

                art::Ptr<simb::MCTruth> truth(mctruthListHandle, i);
                if (truth->Origin() != simb::kBeamNeutrino) continue;

                beamMCTruth = truth;

                auto const& neutrino = truth->GetNeutrino();
                fNu_TrackID.push_back(neutrino.Nu().TrackId());
                fNu_PDG.push_back(neutrino.Nu().PdgCode());
                fNu_E.push_back(neutrino.Nu().E()); 

                auto const& nu_pos = neutrino.Nu().Position();
                fTrueVtxX.push_back(nu_pos.X());
                fTrueVtxY.push_back(nu_pos.Y());
                fTrueVtxZ.push_back(nu_pos.Z());
                fTrueVtxT.push_back(nu_pos.T());

                int num_neutrino_daughters = neutrino.Nu().NumberDaughters();
                for (int j = 0; j < num_neutrino_daughters; j++) {
                    neutrino_daughters.push_back(neutrino.Nu().Daughter(j));
                }

                fNu_Mode.push_back(neutrino.InteractionType());
                fNu_NuMode.push_back(neutrino.Mode());

                // Generator particles
                std::map<int, int> gen_track_to_pdg;
                for(int k = 0; k < truth->NParticles(); k++){
                    gen_track_to_pdg[truth->GetParticle(i).TrackId()] = truth->GetParticle(k).PdgCode();
                }
                for(int ip = 0; ip < truth->NParticles(); ip++){
                    simb::MCParticle particle = truth->GetParticle(ip);

                    if (particle.StatusCode() == 1){
                        particlePDG.push_back(particle.PdgCode());
                        status_code.push_back(particle.StatusCode());
                        particleMotherID.push_back(particle.Mother());
                        particleTrackID.push_back(particle.TrackId());
                        particleKE.push_back(particle.E() - particle.Mass());
                        particleP.push_back(particle.P());
                        particleE.push_back(particle.E());
                        particlePx.push_back(particle.Px());
                        particlePy.push_back(particle.Py());
                        particlePz.push_back(particle.Pz());
                        particle_interaction_index.push_back(interaction_idx);


                        int m_pdg = -999;
                        int m_id = particle.Mother();

                        if (gen_track_to_pdg.count(m_id)) {
                            m_pdg = gen_track_to_pdg[m_id];
                        }
                        particleMotherPDG.push_back(m_pdg);

                    }

                    
                }
                
                // Geant4 Stage particles
                std::vector<art::Ptr<simb::MCParticle>> associatedGeantParticles = geantAssns.at(truth_iter);
                for (auto const& geant_part : associatedGeantParticles) {
                    geant_PDG.push_back(geant_part->PdgCode());
                    geant_TrackID.push_back(geant_part->TrackId());
                    geant_MotherID.push_back(geant_part->Mother());

                    geant_startX.push_back(geant_part->Vx());
                    geant_startY.push_back(geant_part->Vy());
                    geant_startZ.push_back(geant_part->Vz());

                    geant_endX.push_back(geant_part->EndX());
                    geant_endY.push_back(geant_part->EndY());
                    geant_endZ.push_back(geant_part->EndZ());

                    geant_P.push_back(geant_part->P());
                    geant_Px.push_back(geant_part->Px());
                    geant_Py.push_back(geant_part->Py());
                    geant_Pz.push_back(geant_part->Pz());

                    nu_g4_track_ids.insert(geant_part->TrackId());
                }
                interaction_idx++;
                truth_iter++;
            }
        }

        geant_HitsU.resize(geant_TrackID.size(), 0);
        geant_HitsV.resize(geant_TrackID.size(), 0);
        geant_HitsW.resize(geant_TrackID.size(), 0);

        /* ********************************************************
        Slicing, NuScores, and Hit Profiling
        ******************************************************** */
        auto sliceHandle = e.getValidHandle<std::vector<recob::Slice>>(fSliceLabel);
        auto pfpHandle   = e.getValidHandle<std::vector<recob::PFParticle>>(fPFParticleLabel);
        auto trackHandle = e.getValidHandle<std::vector<recob::Track>>(fTrackLabel);

        art::FindManyP<recob::PFParticle> pfpsFromSlices(sliceHandle, e, fSliceLabel);
        art::FindManyP<larpandoraobj::PFParticleMetadata> pfpMetadataAssoc(pfpHandle, e, fPFParticleLabel);
        art::FindManyP<recob::PFParticle> pfpsFromTracks(trackHandle, e, fTrackLabel);
        art::FindManyP<recob::Slice> slicesFromPFPs(pfpHandle, e, fPFParticleLabel);
        art::FindManyP<recob::Hit> trackHitAssoc(trackHandle, e, fTrackLabel);
        art::FindManyP<recob::Vertex> pfpVertexAssoc(pfpHandle, e, fPFParticleLabel);
        
        std::unique_ptr<art::FindManyP<sbn::OpT0Finder>> opt0Assns;
        try { opt0Assns = std::make_unique<art::FindManyP<sbn::OpT0Finder>>(sliceHandle, e, fOpT0Label); } catch(...) {}

        // Needed for purity/completeness calculations
        art::FindManyP<recob::Hit> hitsFromSlices(sliceHandle, e, fSliceLabel);


        art::Handle<std::vector<recob::Hit>> globalHitHandle;
        
        // Total true neutrino hits in the whole event. This is used to calculate the purity/completeness of each slice.
        fEvent_TotalTrueNuHits = 0;
        e.getByLabel(fHitLabel, globalHitHandle);

        std::map<int, size_t> sliceToNuPFP;

        art::FindManyP<simb::MCParticle, anab::BackTrackerHitMatchingData> hitTruthAssns(globalHitHandle, e, "gaushitTruthMatch");

        if (e.getByLabel(fHitLabel, globalHitHandle)) {

            if (!hitTruthAssns.isValid()) {
                mf::LogVerbatim("AnalyzeHyperon") << "WARNING: Could not find gaushitTruthMatch associations!";
            } 
            else {
                // Loop over every reconstructed hit inside the entire event data block
                for (size_t j = 0; j < globalHitHandle->size(); ++j) {
                    art::Ptr<recob::Hit> hit(globalHitHandle, j);
                    bool hit_belongs_to_nu = false;
                    
                    auto const& particles = hitTruthAssns.at(hit.key());
                    std::set<int> counted_geant_indices;


                    for (const auto& true_part : particles) {
                        int current_id = std::abs(true_part->TrackId());
                        const art::Ptr<simb::MCTruth> hitMCTruth = pi_serv->TrackIdToMCTruth_P(std::abs(true_part->TrackId()));
                        
                        // Check if its origin is our Beam Neutrino - if not, record this!
                        if (hitMCTruth.isNonnull() && hitMCTruth->Origin() == simb::kBeamNeutrino) {
                            hit_belongs_to_nu = true;

                            // Climb the Geant4 family tree
                            while (current_id > 0) {
                                auto it = std::find(geant_TrackID.begin(), geant_TrackID.end(), current_id);
                                if (it != geant_TrackID.end()) {
                                    counted_geant_indices.insert(std::distance(geant_TrackID.begin(), it));
                                    break; 
                                }
                                const simb::MCParticle* part = pi_serv->TrackIdToParticle_P(current_id);
                                if (part && part->Mother() != 0 && part->Mother() != current_id) {
                                    current_id = std::abs(part->Mother());
                                } else {
                                    break; 
                                }
                            }   
                        }
                    }
                    
                    if (hit_belongs_to_nu) {
                        fEvent_TotalTrueNuHits++;
                        geo::View_t view = hit->View();
                        for (int idx : counted_geant_indices) {
                            if (view == geo::kU) geant_HitsU[idx]++;
                            else if (view == geo::kV) geant_HitsV[idx]++;
                            else if (view == geo::kW) geant_HitsW[idx]++;
                        }
                    }
                }
            }
        }

        for (size_t i = 0; i < sliceHandle->size(); ++i) {
            art::Ptr<recob::Slice> slice(sliceHandle, i);
            fSlice_ID.push_back(slice->ID());

            float current_nu_score = -1.0;
            float current_vtx_x = -999.0;
            float current_vtx_y = -999.0;
            float current_vtx_z = -999.0;
            auto slicePFPs = pfpsFromSlices.at(i);

            for (const auto& pfp : slicePFPs) {

                if (pfp->IsPrimary() && std::abs(pfp->PdgCode()) == 14) {
                    sliceToNuPFP[slice->ID()] = pfp->Self();
                    auto metadataVec = pfpMetadataAssoc.at(pfp.key());
                    for (const auto& metadata : metadataVec) {
                        const auto& propertiesMap = metadata->GetPropertiesMap();
                        if (propertiesMap.find("NuScore") != propertiesMap.end()) {
                            current_nu_score = propertiesMap.at("NuScore");
                        }
                        // Some check for clear comsic needed - TJ
                    }

                    auto vertices = pfpVertexAssoc.at(pfp.key());
                    if (!vertices.empty()) {
                        double vtx[3];
                        vertices.front()->XYZ(vtx);
                        current_vtx_x = vtx[0];
                        current_vtx_y = vtx[1];
                        current_vtx_z = vtx[2];
                    }

                }
            } // loop over pfps in slice


            float current_opt0_score = -1.0;
            if (opt0Assns && opt0Assns->isValid()) {
                auto opt0s = opt0Assns->at(i);
                if (!opt0s.empty()) {
                    current_opt0_score = opt0s.front()->score;
                }
            }
            fSlice_Opt0Score.push_back(current_opt0_score);

            fSlice_NuScore.push_back(current_nu_score);
            fSlice_VtxX.push_back(current_vtx_x);
            fSlice_VtxY.push_back(current_vtx_y);
            fSlice_VtxZ.push_back(current_vtx_z);
            
            // Profile hits within this slice
            auto const& sliceHits = hitsFromSlices.at(i);
            int total_slice_hits = sliceHits.size();
            int true_nu_slice_hits = 0;
            // Loop over every hit in the slice
            for (const auto& hitPtr : sliceHits) {
                bool hit_belongs_to_nu = false;
                
                if (hitTruthAssns.isValid()) {
                    auto const& particles = hitTruthAssns.at(hitPtr.key());
                    
                    for (const auto& true_part : particles) {
                        const art::Ptr<simb::MCTruth> hitMCTruth = pi_serv->TrackIdToMCTruth_P(std::abs(true_part->TrackId()));
                        
                        if (hitMCTruth.isNonnull() && hitMCTruth->Origin() == simb::kBeamNeutrino) {
                            hit_belongs_to_nu = true;
                            break;
                        }
                    }
                }
                
                if (hit_belongs_to_nu) {
                    true_nu_slice_hits++;
                }
            }

            fSlice_TotalHits.push_back(total_slice_hits);
            fSlice_TrueNuHits.push_back(true_nu_slice_hits);




        } // loop over slices

        /* ********************************************************
        Track & Shower Reconstruction Diagnostics
        ******************************************************** */
        art::Handle<std::vector<recob::Shower>> showerHandle;
        bool has_showers = e.getByLabel(fShowerLabel, showerHandle);

        std::unique_ptr<art::FindManyP<anab::ParticleID>> trackPIDAssns;
        trackPIDAssns = std::make_unique<art::FindManyP<anab::ParticleID>>(trackHandle, e, fPIDLabel); 

        std::unique_ptr<art::FindManyP<recob::Hit>> showerHitAssoc;
        std::unique_ptr<art::FindManyP<recob::PFParticle>> pfpsFromShowers;
        
        if (has_showers) {
            showerHitAssoc = std::make_unique<art::FindManyP<recob::Hit>>(showerHandle, e, fShowerLabel);
            pfpsFromShowers = std::make_unique<art::FindManyP<recob::PFParticle>>(showerHandle, e, fShowerLabel);
        }

        // Start track loop 
        for (size_t i = 0; i < trackHandle->size(); ++i) {
            art::Ptr<recob::Track> track(trackHandle, i);
            fTrack_ID.push_back(track->ID());

            // Export Track Geometry
            fTrack_Length.push_back(track->Length());
            fTrack_StartX.push_back(track->Vertex().X());
            fTrack_StartY.push_back(track->Vertex().Y());
            fTrack_StartZ.push_back(track->Vertex().Z());
            fTrack_EndX.push_back(track->End().X());
            fTrack_EndY.push_back(track->End().Y());
            fTrack_EndZ.push_back(track->End().Z());

            int best_geant_id = -999;
            int max_hits = -1;
            std::map<int, int> trackTruthMap;
            auto trackHits = trackHitAssoc.at(i);

            if (hitTruthAssns.isValid()) {
                for (const auto& hitPtr : trackHits) {
                    auto const& particles = hitTruthAssns.at(hitPtr.key());
                    for (const auto& true_part : particles) {
                        int current_id = std::abs(true_part->TrackId());
                        // Climb the Geant4 family tree
                        while (current_id > 0) {
                            if (std::find(geant_TrackID.begin(), geant_TrackID.end(), current_id) != geant_TrackID.end()) {
                                trackTruthMap[current_id]++;
                                break; 
                            }
                            const simb::MCParticle* part = pi_serv->TrackIdToParticle_P(current_id);
                            if (part && part->Mother() != 0 && part->Mother() != current_id) {
                                current_id = std::abs(part->Mother());
                            } else {
                                break; 
                            }
                        }
                    }
                }
            }
            
            for (auto const& [tid, count] : trackTruthMap) {
                if (count > max_hits) { max_hits = count; best_geant_id = tid; }
            }
            fTrack_TrueGeantID.push_back(best_geant_id);

            // PID extraction
            float chi2_proton = -1.0;
            float chi2_muon = -1.0;

            if (trackPIDAssns && trackPIDAssns->isValid()) {
                auto pids = trackPIDAssns->at(i);
                float p_chi2_sum = 0; int p_chi2_n = 0;
                float m_chi2_sum = 0; int m_chi2_n = 0;

                std::cout<<"We are calculating pids :) "<<std::endl;
                
                for (auto const& pid : pids) {
                    for (auto const& score : pid->ParticleIDAlgScores()) {
                        if (score.fAlgName.find("Chi2") != std::string::npos && score.fValue > 0) {
                            if (score.fAssumedPdg == 2212) {
                                p_chi2_sum += score.fValue; p_chi2_n++;
                            } else if (score.fAssumedPdg == 13) {
                                m_chi2_sum += score.fValue; m_chi2_n++;
                            }
                        }
                    }
                }
                if (p_chi2_n > 0) chi2_proton = p_chi2_sum / p_chi2_n;
                if (m_chi2_n > 0) chi2_muon = m_chi2_sum / m_chi2_n;
            }
            fTrack_Chi2Proton.push_back(chi2_proton);
            fTrack_Chi2Muon.push_back(chi2_muon);

            float t_score = -1.0;
            int slice_id = -999;
            int is_primary = 0;
            auto pfps = pfpsFromTracks.at(i);
            
            if (!pfps.empty()) {
                auto metas = pfpMetadataAssoc.at(pfps.front().key());
                if (!metas.empty() && metas.front()->GetPropertiesMap().count("TrackScore")) {
                    t_score = metas.front()->GetPropertiesMap().at("TrackScore");
                }
                auto slices = slicesFromPFPs.at(pfps.front().key());
                if (!slices.empty()){
                    slice_id = slices.front()->ID();
                    if (sliceToNuPFP.count(slice_id) && pfps.front()->Parent() == sliceToNuPFP[slice_id]) {
                        is_primary = 1;
                    }
                
                }
            }
            fTrack_TrackScore.push_back(t_score);
            fTrack_SliceID.push_back(slice_id);
            fTrack_IsPrimary.push_back(is_primary);
        } // End track loop

        // Shower loop
        if (has_showers) {
            for (size_t i = 0; i < showerHandle->size(); ++i) {
                art::Ptr<recob::Shower> shower(showerHandle, i);
                fShower_ID.push_back(shower->ID());

                int best_geant_id = -999;
                int max_hits = -1;
                std::map<int, int> showerTruthMap;
                auto showerHits = showerHitAssoc->at(i);

                // FIX: Using hitTruthAssns pointer instead of broken HitToTrackIDEs!
                if (hitTruthAssns.isValid()) {
                    for (const auto& hitPtr : showerHits) {
                        auto const& particles = hitTruthAssns.at(hitPtr.key());
                        for (const auto& true_part : particles) {
                            int current_id = std::abs(true_part->TrackId());
                            while (current_id > 0) {
                                if (std::find(geant_TrackID.begin(), geant_TrackID.end(), current_id) != geant_TrackID.end()) {
                                    showerTruthMap[current_id]++;
                                    break; 
                                }
                                const simb::MCParticle* part = pi_serv->TrackIdToParticle_P(current_id);
                                if (part && part->Mother() != 0 && part->Mother() != current_id) {
                                    current_id = std::abs(part->Mother());
                                } else {
                                    break; 
                                }
                            }
                        }
                    }
                }
                for (auto const& [tid, count] : showerTruthMap) {
                    if (count > max_hits) { max_hits = count; best_geant_id = tid; }
                }
                fShower_TrueGeantID.push_back(best_geant_id);

                float t_score = -1.0;
                int slice_id = -999;
                int is_primary = 0;
                auto pfps = pfpsFromShowers->at(i);
                
                if (!pfps.empty()) {
                    auto metas = pfpMetadataAssoc.at(pfps.front().key());
                    if (!metas.empty() && metas.front()->GetPropertiesMap().count("TrackScore")) {
                        t_score = metas.front()->GetPropertiesMap().at("TrackScore");
                    }
                    auto slices = slicesFromPFPs.at(pfps.front().key());
                    if (!slices.empty()){
                        slice_id = slices.front()->ID();
                        if (sliceToNuPFP.count(slice_id) && pfps.front()->Parent() == sliceToNuPFP[slice_id]) {
                            is_primary = 1;
                        }
                    }
                }
                fShower_TrackScore.push_back(t_score);
                fShower_SliceID.push_back(slice_id);
                fShower_IsPrimary.push_back(is_primary);
            }
        } // End shower loop

        
        fTree->Fill();
    }
    
    void AnalyzeHyperon::endSubRun(art::SubRun const& sr) {
        fRun_sr = sr.run();
        fSubRun_sr = sr.subRun();
        fPOT = 0.0;

        art::Handle<sumdata::POTSummary> potHandle;
        if (sr.getByLabel("generator", potHandle)) {
            fPOT = potHandle->totpot;
        } else {
            mf::LogVerbatim("AnalyzeHyperon") << "Warning: No POTSummary found in this SubRun!";
        }

        fSubRunTree->Fill();
    } 

    DEFINE_ART_MODULE(AnalyzeHyperon)
}