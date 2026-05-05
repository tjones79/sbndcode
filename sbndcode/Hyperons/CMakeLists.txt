set(
  MODULE_LIBRARIES
  larcorealg::Geometry
  larcore::Geometry_Geometry_service
  larsim::Simulation
  lardataobj::Simulation
  larsim::MCCheater_BackTrackerService_service
  larsim::Utils
  lardata::Utilities
  lardata::ArtDataHelper
  larevt::Filters
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  lardataobj::MCBase
  sbnobj::Common_Reco
  larreco::RecoAlg
  larsim::MCCheater_BackTrackerService_service
  lardata::Utilities
  lardata::ArtDataHelper
  larevt::Filters
  lardataobj::RawData
  lardataobj::RecoBase
  lardataobj::AnalysisBase
  lardataobj::MCBase
  sbnobj::Common_Reco
  larreco::RecoAlg
  lardata::RecoObjects
  larpandora::LArPandoraInterface
  nusimdata::SimulationBase
  art::Framework_Core
  art::Framework_Principal
  art::Framework_Services_Registry
  art::Persistency_Common
  art::Persistency_Provenance
  art::Utilities
  messagefacility::MF_MessageLogger
  fhiclcpp::fhiclcpp
  ROOT::Geom
  ROOT::XMLIO
  ROOT::Gdml
  ROOT::Core
  sbndcode_RecoUtils
  )

find_package(cetmodules REQUIRED)
cet_build_plugin(AnalyzeEvents art::Module IMPL_SOURCE AnalyzeEvents_module.cc LIBRARIES ${MODULE_LIBRARIES} )
cet_build_plugin(AnalyzeHyperon art::Module IMPL_SOURCE AnalyzeHyperon_module.cc LIBRARIES ${MODULE_LIBRARIES})

install_headers()
install_fhicl()
install_source()