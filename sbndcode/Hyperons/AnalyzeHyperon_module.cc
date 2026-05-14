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


#include "art_root_io/TFileService.h"
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
#include "lardataobj/RecoBase/Vertex.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "sbnobj/Common/Reco/OpT0FinderResult.h"
#include "sbnobj/Common/Reco/CRUMBSResult.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "sbnobj/Common/Reco/SimpleFlashMatchVars.h"
#include "sbndcode/CosmicId/Algs/PandoraNuScoreCosmicIdAlg.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "larsim/Utils/TruthMatchUtils.h"

// // Root Includes
#include "TLorentzVector.h"
#include <iostream>
#include <vector>
#include <TTree.h>
#include <TH1.h>
#include <fstream>
#include <string>
#include <exception>

namespace hyperon {

    class AnalyzeHyperon : public art::EDAnalyzer {

    public:
        explicit AnalyzeHyperon(fhicl::ParameterSet const& p);
        AnalyzeHyperon(AnalyzeHyperon const&) = delete;
        AnalyzeHyperon(AnalyzeHyperon&&) = delete;

        void beginJob() override;
        void analyze(art::Event const& e) override;

    private:
        // The ROOT Tree
        TTree* fTree;

        // Event info
        int fRun;
        int fEvent;

        // Primary interaction vertex
        float fNuVertexX;
        float fNuVertexY;
        float fNuVertexZ;

        // All particle vertices
        std::vector<int>   fTruePartPDG;
        std::vector<int>   fTruePartMother;
        std::vector<int>   fTruePartTrackID;
        std::vector<float> fTruePartStartX;
        std::vector<float> fTruePartStartY;
        std::vector<float> fTruePartStartZ;

        // Labels
        std::string fSliceLabel = "pandora";
        std::string fPFParticleLabel = "pandora";
        std::string fTrackLabel = "pandoraTrack";
        std::string fCrumbsLabel = "crumbs";

        // Reco Tree Variables
        int fNSlices;
        std::vector<int> fSliceParentPDG; // What is pandora's PDG guess?
        std::vector<float> fSliceNuScore; // The BDT score (0.0 -> 1.0)

        // Reco Tree Variables (Track Level - Inside Neutrino Slice)
        int fNTracks;
        std::vector<int>   fNuTrackID;
        std::vector<float> fNuTrackStartX;
        std::vector<float> fNuTrackStartY;
        std::vector<float> fNuTrackStartZ;
        std::vector<float> fNuTrackLength;

        // True track stuff
        std::vector<int> fNuTrackTruePDG;
        std::vector<int> fNuTrackTrueMotherTrackID;
        std::vector<int> fNuTrackTrueID;

        std::vector<float> fSliceVertexX;
        std::vector<float> fSliceVertexY;
        std::vector<float> fSliceVertexZ;

        // Simple tag for true neutrino slice
        bool fHasTrueNuSlice;

        std::vector<int> fSliceIsTrueNeutrino;


        std::vector<float> SliceCRUMBSScore;

        
    };

    AnalyzeHyperon::AnalyzeHyperon(fhicl::ParameterSet const& p)
        : EDAnalyzer(p) {}

    
    void AnalyzeHyperon::beginJob() {
        art::ServiceHandle<art::TFileService> tfs;
        fTree = tfs->make<TTree>("tree", "Minimal Hyperon Truth Tree");

        fTree->Branch("run", &fRun, "run/I");
        fTree->Branch("event", &fEvent, "event/I");

        // Primary Vertex Branches - we need these
        fTree->Branch("nuVertexX", &fNuVertexX, "nuVertexX/F");
        fTree->Branch("nuVertexY", &fNuVertexY, "nuVertexY/F");
        fTree->Branch("nuVertexZ", &fNuVertexZ, "nuVertexZ/F");

        // All particle
        fTree->Branch("TruePartPDG", &fTruePartPDG);
        fTree->Branch("TruePartMother", &fTruePartMother);
        fTree->Branch("TruePartTrackID", &fTruePartTrackID);
        fTree->Branch("TruePartStartX", &fTruePartStartX);
        fTree->Branch("TruePartStartY", &fTruePartStartY);
        fTree->Branch("TruePartStartZ", &fTruePartStartZ);

        // Reco level slice stuff
        fTree->Branch("NSlices", &fNSlices, "NSlices/I");
        fTree->Branch("SliceParentPDG", &fSliceParentPDG);
        fTree->Branch("SliceNuScore", &fSliceNuScore);

        // Reco level track stuff (Neutrino children only)
        fTree->Branch("NTracks", &fNTracks, "NTracks/I");
        fTree->Branch("NuTrackID", &fNuTrackID);
        fTree->Branch("NuTrackStartX", &fNuTrackStartX);
        fTree->Branch("NuTrackStartY", &fNuTrackStartY);
        fTree->Branch("NuTrackStartZ", &fNuTrackStartZ);
        fTree->Branch("NuTrackLength", &fNuTrackLength);

        // True track stuff
        fTree->Branch("NuTrackTruePDG", &fNuTrackTruePDG);
        fTree->Branch("NuTrackTrueMotherTrackID", &fNuTrackTrueMotherTrackID);
        fTree->Branch("NuTrackTrueID", &fNuTrackTrueID);

        fTree->Branch("SliceVertexX", &fSliceVertexX);
        fTree->Branch("SliceVertexY", &fSliceVertexY);
        fTree->Branch("SliceVertexZ", &fSliceVertexZ);


        // Simple truth tag for slice
        fTree->Branch("HasTrueNuSlice", &fHasTrueNuSlice, "HasTrueNuSlice/O");
        fTree->Branch("SliceIsTrueNeutrino",&fSliceIsTrueNeutrino);

        // CRUMBS
        fTree->Branch("SliceCRUMBSScore", &SliceCRUMBSScore);
        
    }

    // Analyze, this runs on every event
    void AnalyzeHyperon::analyze(art::Event const& e) {

        // Cleanup
        fTruePartPDG.clear();
        fTruePartMother.clear();
        fTruePartTrackID.clear();
        fTruePartStartX.clear();
        fTruePartStartY.clear();
        fTruePartStartZ.clear();

        fSliceParentPDG.clear();
        fSliceNuScore.clear();
        fNSlices = 0;

        fNuTrackID.clear();
        fNuTrackStartX.clear();
        fNuTrackStartY.clear();
        fNuTrackStartZ.clear();
        fNuTrackLength.clear();
        fNTracks = 0;

        fNuTrackTruePDG.clear();
        fNuTrackTrueID.clear();
        fNuTrackTrueMotherTrackID.clear();

        fHasTrueNuSlice = false;
        fSliceIsTrueNeutrino.clear();

        fSliceVertexX.clear();
        fSliceVertexY.clear();
        fSliceVertexZ.clear();

        SliceCRUMBSScore.clear();

        fRun = e.run();
        fEvent = e.event();


        /* 
        ********************************************************
        1. Get the Primary Neutrino Vertex (Truth stuff)
        ********************************************************
        */ 
        art::Handle<std::vector<simb::MCTruth>> mctruthListHandle;
        if (e.getByLabel("generator", mctruthListHandle) && !mctruthListHandle->empty()) {
            
            // Grab the first neutrino interaction in the event
            const simb::MCTruth& truth = mctruthListHandle->at(0);
            
            

            const simb::MCNeutrino& nu = truth.GetNeutrino();



            
            // Save the exact coordinate the neutrino interacted
            fNuVertexX = nu.Nu().Vx();
            fNuVertexY = nu.Nu().Vy();
            fNuVertexZ = nu.Nu().Vz();
        }

        /*
        ********************************************************
        2. Get All Particle Vertices (MCParticle stuff)
        ********************************************************
        */
        art::Handle<std::vector<simb::MCParticle>> mcpartListHandle;
        if (e.getByLabel("largeant", mcpartListHandle)) {
            
            // Loop over every single particle simulated by Geant4
            for (auto const& particle : *mcpartListHandle) {
                
                fTruePartPDG.push_back(particle.PdgCode());
                fTruePartMother.push_back(particle.Mother());
                fTruePartTrackID.push_back(particle.TrackId());
                
                // particle.Vx() gives the starting position of this specific track
                fTruePartStartX.push_back(particle.Vx());
                fTruePartStartY.push_back(particle.Vy());
                fTruePartStartZ.push_back(particle.Vz());
            }
        }


        /*
        ********************************************************
        3. Slices & Tracks (Reco stuff)
        ********************************************************
        */

        // Get the Slices, PFParticles and tracks from the event
        art::Handle<std::vector<recob::Slice>> sliceHandle; 
        art::Handle<std::vector<recob::PFParticle>> pfpHandle;
        art::Handle<std::vector<recob::Track>> trackHandle;

        auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(e);
        art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;

        if (e.getByLabel(fSliceLabel, sliceHandle) && e.getByLabel(fPFParticleLabel, pfpHandle) && e.getByLabel(fTrackLabel, trackHandle)) {
        
            fNSlices = sliceHandle->size();

            // Find associations
            // "Give me the PFParticles for each Slice"
            art::FindManyP<recob::PFParticle> sliceToPFP(sliceHandle, e, fSliceLabel);
            // "Give me the Metadata nuscore for each PFParticle"
            art::FindManyP<larpandoraobj::PFParticleMetadata> pfpToMetadata(pfpHandle, e, fPFParticleLabel);
            // "Give me the Track associated with this PFParticle"
            art::FindManyP<recob::Track> pfpToTrack(pfpHandle, e, fTrackLabel);
            // "Give me all the Hits that make up this Track"
            art::FindManyP<recob::Hit> trackToHit(trackHandle, e, fTrackLabel);
            // Give me all the vertices associated with this pfp
            art::FindManyP<recob::Vertex> pfpToVertex(pfpHandle, e, fPFParticleLabel);
            // Give me all the hits associated with this slice.
            art::FindManyP<recob::Hit> sliceToHit(sliceHandle, e, fSliceLabel);

            art::FindManyP<sbn::CRUMBSResult> fmCrumbs(sliceHandle, e, fCrumbsLabel);

            // slice loop
            for (size_t i_slice = 0; i_slice < sliceHandle->size(); ++i_slice) {
                
                // Get all articles in this slice
                std::vector<art::Ptr<recob::PFParticle>> sliceParticles = sliceToPFP.at(i_slice);
                
                bool isNeutrinoSlice = false;

                // Find parent and check if it is a neutrino
                for (const art::Ptr<recob::PFParticle>& pfp : sliceParticles) {
                    
                    float crumbs_score = -999.0;
                    if (fmCrumbs.isValid()){

                        std::cout<<"CRUMBS Association is valid!!"<<std::endl;
                    }
                    if (fmCrumbs.at(i_slice).size() > 0){

                        std::cout<<"There are crumbs scores for these slices!!!"<<std::endl;
                    }
                    if (fmCrumbs.isValid() && fmCrumbs.at(i_slice).size() > 0) {
                        crumbs_score = fmCrumbs.at(i_slice).front()->score;
                    }
                    SliceCRUMBSScore.push_back(crumbs_score);
                    if (pfp->IsPrimary()) {
                        
                        // Pandora's guess
                        fSliceParentPDG.push_back(pfp->PdgCode());
                        
                        // Fill nuScore
                        float nuScore = -1.0; 
                        std::vector<art::Ptr<larpandoraobj::PFParticleMetadata>> metadataVec = pfpToMetadata.at(pfp.key());
                        
                        if (!metadataVec.empty()) {
                            const auto& propertiesMap = metadataVec.front()->GetPropertiesMap();
                            if (propertiesMap.find("NuScore") != propertiesMap.end()) {
                                nuScore = propertiesMap.at("NuScore");
                            }
                        }
                        
                        fSliceNuScore.push_back(nuScore);

                        std::vector<art::Ptr<recob::Vertex>> vertices = pfpToVertex.at(pfp.key());

                        if (!vertices.empty()) {
                            fSliceVertexX.push_back(vertices.front()->position().X());
                            fSliceVertexY.push_back(vertices.front()->position().Y());
                            fSliceVertexZ.push_back(vertices.front()->position().Z());
                        } else {
                            // If Pandora failed to make a vertex, use a dummy value
                            fSliceVertexX.push_back(-9999.0);
                            fSliceVertexY.push_back(-9999.0);
                            fSliceVertexZ.push_back(-9999.0);
                        }

                        int isTrueNu = 0;
                        
                        // Get all hits in this specific slice
                        std::vector<art::Ptr<recob::Hit>> sliceHits = sliceToHit.at(i_slice);
                        
                        // Ask BackTracker for the Geant4 ID that contributed the most hits
                        int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, sliceHits, 1);

                        // If a valid particle won the election
                        if (TruthMatchUtils::Valid(g4id)) {
                            // Trace it back to the event generator
                            const art::Ptr<simb::MCTruth> truth = pi_serv->TrackIdToMCTruth_P(g4id);
                            
                            if (truth.isNonnull()) {
                                // Was it generated as a beam neutrino?
                                if (truth->Origin() == simb::kBeamNeutrino) {
                                    isTrueNu = 1;
                                }
                            }
                        }
                        
                        fSliceIsTrueNeutrino.push_back(isTrueNu);


                        // Flag slice
                        if (pfp->PdgCode()== 14 || pfp->PdgCode() == 12) {

                            isNeutrinoSlice = true;

                        }
                        
                        // We found the primary particle for this slice, so we can stop looking
                        break; 
                    }
                }
                // Extract tracks if they are a neutrino slice
                if (isNeutrinoSlice) {
                    
                    for (const art::Ptr<recob::PFParticle>& pfp : sliceParticles) {

                        // Neutrinos don't have tracks
                        if (pfp->IsPrimary()) continue;
                        
                        std::vector<art::Ptr<recob::Track>> tracks = pfpToTrack.at(pfp.key());

                        if (!tracks.empty()) {
                            const art::Ptr<recob::Track>& track = tracks.front();
                            
                            fNuTrackID.push_back(track->ID());
                            fNuTrackStartX.push_back(track->Vertex().X());
                            fNuTrackStartY.push_back(track->Vertex().Y());
                            fNuTrackStartZ.push_back(track->Vertex().Z());
                            fNuTrackLength.push_back(track->Length());
                            // Only fill NTracks when slice is neutrino
                            fNTracks++;

                            // Back tracker
                            int truePDG = 0;
                            int trueMotherTrackID = 0;
                            int trueTrackID = 0;
                            std::vector<art::Ptr<recob::Hit>> trackHits = trackToHit.at(track.key());

                            // Dom's utility
                            int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockData, trackHits, 1);

                            // If a valid true particle won the election
                            if (TruthMatchUtils::Valid(g4id)) {
                                
                                const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);
                                
                                if (matched_mcparticle) {
                                    truePDG = matched_mcparticle->PdgCode();
                                    trueMotherTrackID = matched_mcparticle->Mother();
                                    trueTrackID = matched_mcparticle->TrackId();

                                }
                            }
            
                            fNuTrackTruePDG.push_back(truePDG);
                            fNuTrackTrueMotherTrackID.push_back(trueMotherTrackID);
                            fNuTrackTrueID.push_back(trueTrackID);

                        }

                    }
                    // Found neutrino slice
                    break;
                }
        }
    }
        
        fTree->Fill();
    }

    DEFINE_ART_MODULE(AnalyzeHyperon)
}