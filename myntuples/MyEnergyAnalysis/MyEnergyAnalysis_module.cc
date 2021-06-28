/**
 * @file   MyEnergyAnalysis_module.cc
 * @brief  A file to read and analyze art::Event records from a DUNE FD MC file,
 * @author Wei Shi (wei.shi.1@stonybrook.edu)
 *
 * Adapted from https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/_AnalysisExample_
 */

// Include headers: starting from LArSoft and going up the software
// layers (nusimdata, art, etc.), ending with C++ is standard.

// LArSoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcorealg/Geometry/GeometryCore.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "larsim/Simulation/LArG4Parameters.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

// Framework includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/Exception.h"

// Utility libraries
#include "cetlib/pow.h" // cet::sum_of_squares()
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// ROOT includes
#include "TH1.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"

// C++ includes
#include <cmath>
#include <map>

namespace {

  // This is a local namespace. Stuff declared here will not be
  // visible beyond this file. We will define functions at the end,
  // but we declare them here so that the module can freely use them.

  // Utility function to get the diagonal of the detector
  double DetectorDiagonal(geo::GeometryCore const& geom);

  // Sort MC particles based on its start momentum P(0)
  bool MomentumOrderMCParticle(const simb::MCParticle*, const simb::MCParticle*);

  // Comparison routine for using std::lower/upper_bound to search TDCIDE vectors.
  bool TDCIDETimeCompare(const sim::TDCIDE&, const sim::TDCIDE&);

} // local namespace

// An outside package call this module like lar::example::MyEnergyAnalysis

namespace lar {
  namespace example {

    // BEGIN MyEnergyAnalysis group
    // -----------------------------------------------
    // class definition
    /**
     * This class produces a ROOT tree that contains information
     * from the generated/simulated and reconstructed particles.
     *
     * Configuration parameters
     * =========================
     *
     * - *SimulationLabel* (string, default: "largeant"): tag of the input data
     *   product with the detector simulation information (typically an instance
     *   of the LArG4 module)
     *
     * - *HitLabel* (string, mandatory): tag of the input data product with
     *   reconstructed hits
     *
     * - *ClusterLabel* (string, mandatory): tag of the input data product with
     *   reconstructed clusters
     */
    class MyEnergyAnalysis : public art::EDAnalyzer {
    public:

      // This structure describes the configuration parameters of the module.
      // Any missing or unknown parameters will generate a configuration error.

      struct Config {

        // Save some typing:
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // One Atom for each parameter
        fhicl::Atom<art::InputTag> SimulationLabel{
          Name("SimulationLabel"),
          Comment("tag of the input data product with the detector simulation "
                  "information")};

        fhicl::Atom<art::InputTag> HitLabel{
          Name("HitLabel"),
          Comment("tag of the input data product with reconstructed hits")};

        fhicl::Atom<art::InputTag> ClusterLabel{
          Name("ClusterLabel"),
          Comment("tag of the input data product with reconstructed clusters")};

      }; // Config

      using Parameters = art::EDAnalyzer::Table<Config>;

      /// Constructor: configures the module (see the Config structure above)
      explicit MyEnergyAnalysis(Parameters const& config);

      // This method is called once, at the start of the job. In this
      // example, it will define the histograms and n-tuples we'll
      // write.
      virtual void beginJob() override;

      // This method is called once, at the start of each run. It's a
      // good place to read databases or files that may have
      // run-dependent information.
      virtual void beginRun(const art::Run& run) override;

      // The analysis routine, called once per event.
      virtual void analyze(const art::Event& event) override;

    private:

      // The parameters we'll read from the .fcl file.
      art::InputTag fSimulationProducerLabel; // The name of the producer that tracked simulated particles through the detector
      art::InputTag fHitProducerLabel;        // The name of the producer that created hits
      art::InputTag fClusterProducerLabel;    // The name of the producer that created clusters

      // The n-tuple we'll create.
      TTree* fNtuple;

      // Event info
      int fEvent;  // number of the event being processed
      int fRun;    // number of the run being processed
      int fSubRun; // number of the sub-run being processed

      // Variables related to simulation
      int fSimPDG;                       // MCParticle PDG ID
      int fSimTrackID;                   // GEANT ID of the particle being processed
      int fSim_nEle;                     // Number of Sim electrons (e+/e-) in the event
      int fSim_nMu;                      // Number of Sim muons (mu+/mu-) in the event
      int fSim_nTau;                     // Number of Sim tau leptons (+/-) in the event
      int fSim_nPhoton;                  // Number of Sim photons in the event
      int fSim_nPionNeutral;             // Number of Sim pi+/pi- in the event
      int fSim_nPionCharged;             // Number of Sim pi0 in the event
      int fSim_nNeutron;                 // Number of Sim neutrons in the event
      int fSim_nProton;                  // Number of Sim protons in the event
      double fSim_mu_start_vx;           // x position of the muon trajectory start
      double fSim_mu_start_vy;           // y .....................................
      double fSim_mu_start_vz;           // z .....................................
      double fSim_mu_end_vx;             // x position of the muon trajectory end
      double fSim_mu_end_vy;             // y ...................................
      double fSim_mu_end_vz;             // z ...................................
      double fSim_mu_start_px;           // x momentum of the muon trajectory start
      double fSim_mu_start_py;           // y .....................................
      double fSim_mu_start_pz;           // z .....................................
      double fSim_mu_end_px;             // x momentum of the muon trajectory end
      double fSim_mu_end_py;             // y ...................................
      double fSim_mu_end_pz;             // z ...................................
      double fSim_mu_start_4position[4]; // (x,y,z,t) of the muon trajectory start
      double fSim_mu_end_4position[4];   // ................................ end
      double fSim_mu_start_4mommenta[4]; // (Px,Py,Pz,E) of the muon trajectory start
      double fSim_mu_end_4mommenta[4];   // ................................... end
      double fSim_mu_track_length;       // muon track length
      double fSim_hadronic_Edep;         // hadron deposited energy on collection plane

      // Variables related to reconstruction

      // Other variables that will be shared between different methods.
      geo::GeometryCore const* fGeometryService; // pointer to Geometry provider
      double fElectronsToGeV;                    // conversion factor for no. of ionization electrons to energy deposited in GeV
      int fTriggerOffset;                        // (units of ticks) time of expected neutrino event

    }; // class MyEnergyAnalysis

    // END MyEnergyAnalysis group
    // -------------------------------------------------

    //-----------------------------------------------------------------------
    // class implementation

    //-----------------------------------------------------------------------
    // Constructor

    MyEnergyAnalysis::MyEnergyAnalysis(Parameters const& config)
      : EDAnalyzer(config)
      , fSimulationProducerLabel(config().SimulationLabel())
      , fHitProducerLabel(config().HitLabel())
      , fClusterProducerLabel(config().ClusterLabel())
    {
      // Get a pointer to the geometry service provider.
      fGeometryService = lar::providerFrom<geo::Geometry>();

      auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
      fTriggerOffset = trigger_offset(clock_data);

      // Tell beforehand all the data the module is going to read ("consumes") or
      // might read ("may_consume").
      consumes<std::vector<simb::MCParticle>>(fSimulationProducerLabel);
      consumes<std::vector<sim::SimChannel>>(fSimulationProducerLabel);
      consumes<art::Assns<simb::MCTruth, simb::MCParticle>>(fSimulationProducerLabel);
      consumes<std::vector<recob::Hit>>(fHitProducerLabel);
      consumes<std::vector<recob::Cluster>>(fClusterProducerLabel);
      consumes<art::Assns<recob::Cluster, recob::Hit>>(fHitProducerLabel);
    }

    //-----------------------------------------------------------------------
    void MyEnergyAnalysis::beginJob()
    {
      // Get the detector length
      const double detectorLength = DetectorDiagonal(*fGeometryService);
      std::cout << "Detector length=" << detectorLength << " cm" << std::endl;

      // Access art's TFileService, which will handle creating and writing
      // histograms and n-tuples for us.
      art::ServiceHandle<art::TFileService const> tfs;

      // Define n-tuples
      fNtuple = tfs->make<TTree>("MyTree", "MyTree");

      fNtuple->Branch("Event",                 &fEvent,                  "Event/I");
      fNtuple->Branch("SubRun",                &fSubRun,                 "SubRun/I");
      fNtuple->Branch("Run",                   &fRun,                    "Run/I");
      // Simulation branches Sim*
      fNtuple->Branch("Sim_nMu",               &fSim_nMu,                "Sim_nMu/I");
      // muon position
      fNtuple->Branch("Sim_mu_start_vx",       &fSim_mu_start_vx,        "Sim_mu_start_vx/D");
      fNtuple->Branch("Sim_mu_start_vy",       &fSim_mu_start_vy,        "Sim_mu_start_vy/D");
      fNtuple->Branch("Sim_mu_start_vz",       &fSim_mu_start_vz,        "Sim_mu_start_vz/D");
      fNtuple->Branch("Sim_mu_end_vx",         &fSim_mu_end_vx,          "Sim_mu_end_vx/D");
      fNtuple->Branch("Sim_mu_end_vy",         &fSim_mu_end_vy,          "Sim_mu_end_vy/D");
      fNtuple->Branch("Sim_mu_end_vz",         &fSim_mu_end_vz,          "Sim_mu_end_vz/D");
      // muon momentum
      fNtuple->Branch("Sim_mu_start_px",       &fSim_mu_start_px,        "Sim_mu_start_px/D");
      fNtuple->Branch("Sim_mu_start_py",       &fSim_mu_start_py,        "Sim_mu_start_py/D");
      fNtuple->Branch("Sim_mu_start_pz",       &fSim_mu_start_pz,        "Sim_mu_start_pz/D");
      fNtuple->Branch("Sim_mu_end_px",         &fSim_mu_end_px,          "Sim_mu_end_px/D");
      fNtuple->Branch("Sim_mu_end_py",         &fSim_mu_end_py,          "Sim_mu_end_py/D");
      fNtuple->Branch("Sim_mu_end_pz",         &fSim_mu_end_pz,          "Sim_mu_end_pz/D");
      // Write arrays giving the address of the array: simply the array name.
      fNtuple->Branch("Sim_mu_start_4position", fSim_mu_start_4position, "Sim_mu_start_4position[4]/D");
      fNtuple->Branch("Sim_mu_end_4position",   fSim_mu_end_4position,   "Sim_mu_end_4position[4]/D");
      fNtuple->Branch("Sim_mu_start_4mommenta", fSim_mu_start_4mommenta, "Sim_mu_start_4mommenta[4]/D");
      fNtuple->Branch("Sim_mu_end_4mommenta",   fSim_mu_end_4mommenta,   "Sim_mu_end_4mommenta[4]/D");
      fNtuple->Branch("Sim_mu_track_length",   &fSim_mu_track_length,    "Sim_mu_track_length/D");
      fNtuple->Branch("Sim_hadronic_Edep",     &fSim_hadronic_Edep,      "Sim_hadronic_Edep/D");

      fNtuple->Branch("Sim_nEle",              &fSim_nEle,                "Sim_nEle/I");
      fNtuple->Branch("Sim_nTau",              &fSim_nTau,                "Sim_nTau/I");
      fNtuple->Branch("Sim_nPhoton",           &fSim_nPhoton,             "Sim_nPhoton/I");
      fNtuple->Branch("Sim_nPionNeutral",      &fSim_nPionNeutral,        "Sim_nPionNeutral/I");
      fNtuple->Branch("Sim_nPionCharged",      &fSim_nPionCharged,        "Sim_nPionCharged/I");
      fNtuple->Branch("Sim_nNeutron",          &fSim_nNeutron,            "Sim_nNeutron/I");
      fNtuple->Branch("Sim_nProton",           &fSim_nProton,             "Sim_nProton/I");

      // Reconstruction branches

    }

    //-----------------------------------------------------------------------
    void MyEnergyAnalysis::beginRun(const art::Run& /*run*/)
    {
      // Conversion factor for no. of ionization electrons to energy deposited in GeV
      // The ultimate source of this conversion factor is
      // ${LARCOREOBJ_INC}/larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h.
      art::ServiceHandle<sim::LArG4Parameters const> larParameters;
      fElectronsToGeV = 1. / larParameters->GeVToElectrons();
    }

    //-----------------------------------------------------------------------
    void MyEnergyAnalysis::analyze(const art::Event& event)
    {
      // Fetching basic event information.
      fEvent = event.id().event();
      fRun = event.run();
      fSubRun = event.subRun();

      // Read multiple objects associated with one event:
      // define a "handle" to point to a vector of the objects.
      art::Handle<std::vector<simb::MCParticle>> particleHandle;

      // Then fill the vector with all the objects
      if (!event.getByLabel(fSimulationProducerLabel, particleHandle)) {
        // If no MCParticles in an event, throw an exception to force this module to stop.
        throw cet::exception("MyEnergyAnalysis") << " No simb::MCParticle objects in this event - " << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

      // Get all the simulated channels for the event. These channels
      // include the energy deposited for each simulated track.
      auto simChannelHandle = event.getValidHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);

      // Create a map pf MCParticle to its track ID, to be used for hadronic part later
      std::map<int, const simb::MCParticle*> particleMap;

      //
      // Process Sim MC particles info
      //

      // Store specific particles
      std::vector<const simb::MCParticle*> SimElectrons;
      std::vector<const simb::MCParticle*> SimMuons;
      std::vector<const simb::MCParticle*> SimTaus;
      std::vector<const simb::MCParticle*> SimPhotons;
      std::vector<const simb::MCParticle*> SimNeutralPions;
      std::vector<const simb::MCParticle*> SimChargedPions;
      std::vector<const simb::MCParticle*> SimNeutrons;
      std::vector<const simb::MCParticle*> SimProtons;

      // Loop over the list of particles in the event
      for (auto const& particle : (*particleHandle)) {

        // For the methods you can call for MCParticle, see ${NUSIMDATA_INC}/nusimdata/SimulationBase/MCParticle.h.
        fSimTrackID = particle.TrackId();

        // Add the address of the MCParticle to the map, with the track ID as the key.
        particleMap[fSimTrackID] = &particle;

        // Only for primary particles in the event
        fSimPDG = particle.PdgCode();
        if ( particle.Process() == "primary" ) {
          if ( abs(fSimPDG) == 11 )   SimElectrons.push_back(&particle);
          if ( abs(fSimPDG) == 13 )   SimMuons.push_back(&particle);
          if ( abs(fSimPDG) == 15 )   SimTaus.push_back(&particle);
          if ( abs(fSimPDG) == 22 )   SimPhotons.push_back(&particle);
          if ( abs(fSimPDG) == 111 )  SimNeutralPions.push_back(&particle);
          if ( abs(fSimPDG) == 211 )  SimChargedPions.push_back(&particle);
          if ( abs(fSimPDG) == 2112 ) SimNeutrons.push_back(&particle);
          if ( abs(fSimPDG) == 2212 ) SimProtons.push_back(&particle);
        }

      } // end loop over all particles in the event.

      fSim_nEle         = SimElectrons.size();
      fSim_nMu          = SimMuons.size();
      fSim_nTau         = SimTaus.size();
      fSim_nPhoton      = SimPhotons.size();
      fSim_nPionNeutral = SimNeutralPions.size();
      fSim_nPionCharged = SimChargedPions.size();
      fSim_nNeutron     = SimNeutrons.size();
      fSim_nProton      = SimProtons.size();

      // If multiple sim muons present in event, sort momentum from high to low
      if ( fSim_nMu > 1 ) std::sort(SimMuons.begin(), SimMuons.end(), MomentumOrderMCParticle);

      // Initialize 4-vectors
      for (int i = 0; i < 4; i++) {
        fSim_mu_start_4position[i] = -99.;
        fSim_mu_end_4position[i]   = -99.;
        fSim_mu_start_4mommenta[i] = -99.;
        fSim_mu_end_4mommenta[i]   = -99.;
      }

      // Store info for leading momentum sim muon
      if ( fSim_nMu > 0 ) {

        const simb::MCParticle& leadingmu = *(SimMuons[0]);

        // A muon track consists of a set of 4-positions and 4-mommenta.
        const size_t numberTrajectoryPoints = leadingmu.NumberTrajectoryPoints();

        // For trajectories, as for vectors and arrays, the first point is #0, not #1.
        const int last = numberTrajectoryPoints - 1;
        const TLorentzVector& positionStart = leadingmu.Position(0);
        const TLorentzVector& positionEnd = leadingmu.Position(last);
        const TLorentzVector& momentumStart = leadingmu.Momentum(0);
        const TLorentzVector& momentumEnd = leadingmu.Momentum(last);

        // Fill position and momentum components
        fSim_mu_start_vx = leadingmu.Vx(0);
        fSim_mu_start_vy = leadingmu.Vy(0);
        fSim_mu_start_vz = leadingmu.Vz(0);
        fSim_mu_end_vx = leadingmu.Vx(last);
        fSim_mu_end_vy = leadingmu.Vy(last);
        fSim_mu_end_vz = leadingmu.Vz(last);
        fSim_mu_start_px = leadingmu.Px(0);
        fSim_mu_start_py = leadingmu.Py(0);
        fSim_mu_start_pz = leadingmu.Pz(0);
        fSim_mu_end_px = leadingmu.Px(last);
        fSim_mu_end_py = leadingmu.Py(last);
        fSim_mu_end_pz = leadingmu.Pz(last);

        // Fill arrays with the 4-values.
        positionStart.GetXYZT(fSim_mu_start_4position);
        positionEnd.GetXYZT(fSim_mu_end_4position);
        momentumStart.GetXYZT(fSim_mu_start_4mommenta);
        momentumEnd.GetXYZT(fSim_mu_end_4mommenta);

        // Calculate length using spherical cooridnate system: assume straight track? time negligible?
        const double trackLength = (positionEnd - positionStart).Rho();
        fSim_mu_track_length = trackLength;
      }
      else {
        fSim_mu_start_vx           = -99.;
        fSim_mu_start_vy           = -99.;
        fSim_mu_start_vz           = -99.;
        fSim_mu_end_vx             = -99.;
        fSim_mu_end_vy             = -99.;
        fSim_mu_end_vz             = -99.;
        fSim_mu_start_px           = -99.;
        fSim_mu_start_py           = -99.;
        fSim_mu_start_pz           = -99.;
        fSim_mu_end_px             = -99.;
        fSim_mu_end_py             = -99.;
        fSim_mu_end_pz             = -99.;
        fSim_mu_track_length       = -99.;
      }

      //
      // Calculate sim hadronic deposit energy
      //

      // Loop over the SimChannel objects in the event to look at the energy deposited by particle's track.
      for (auto const& channel : (*simChannelHandle)) {

        // Get the numeric ID associated with this channel.
        // See methods at https://internal.dunescience.org/doxygen/SimChannel_8h_source.html
        auto const channelNumber = channel.Channel();

        // Note: There is more than one plane that reacts to charge in the TPC. We only want to include the
        // energy from the collection plane: geo::kCollection defined in
        // ${LARCOREOBJ_INC}/larcoreobj/SimpleTypesAndConstants/geo_types.h
        // Another way to makre sure the plane is correct is go through channel -> wire -> plane ID, and require planeID to be 0. (Sanil's method)
        // May need to check if two methods give the same energy deposits.
        if (fGeometryService->SignalType(channelNumber) != geo::kCollection) continue;

        // Each channel has a map inside it that connects a time slice to energy deposits in the detector.
        // The full type of this map is std::map<unsigned short, std::vector<sim::IDE>>; we'll use "auto" here
        auto const& timeSlices = channel.TDCIDEMap();

        for (auto const& timeSlice : timeSlices) {

          // For the timeSlices map, the 'first' is a time slice number; The 'second' is a vector of IDE objects.
          auto const& energyDeposits = timeSlice.second;

          // An "energy deposit" object stores how much charge/energy was deposited in a small volume, by which particle, and where.
          // The type of 'energyDeposit' will be sim::IDE, here use auto.
          for (auto const& energyDeposit : energyDeposits) {

            auto search = particleMap.find( energyDeposit.trackID );
            if ( search == particleMap.end() ) continue;
            // "search" points to a pair in the map: <track ID, MCParticle*>
            const simb::MCParticle& particle = *((*search).second);

            // If it's not leptons, we think it's hadronic
            if ( particle.Process() == "primary" && abs(particle.PdgCode()) != 11 && abs(particle.PdgCode()) != 13 && abs(particle.PdgCode()) != 15 ){

              // Note another way to get E deposit is simply energyDeposit.energy
              fSim_hadronic_Edep += energyDeposit.numElectrons * fElectronsToGeV;

              // Get the (x,y,z) of the energy deposit.
              // TVector3 location(energyDeposit.x, energyDeposit.y, energyDeposit.z);

              // Do we need hadroic E depisit position? weighted?

            } // end if hadronic

          }   // end For each energy deposit
        }     // end For each time slice
      }       // end For each SimChannel

      //
      // Access reco Hits.
      //
      // See ${LARDATAOBJ_INC}/lardataobj/RecoBase/Hit.h
      art::Handle<std::vector<recob::Hit>> hitHandle;
      if (!event.getByLabel(fHitProducerLabel, hitHandle)) return;

      auto const clock_data = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);

      for (auto const& hit : (*hitHandle)) {

        // The channel associated with this hit.
        auto hitChannelNumber = hit.Channel();

        // We have a hit. For this example let's just focus on the hits in the collection plane.
        if (fGeometryService->SignalType(hitChannelNumber) != geo::kCollection) continue;

        // std::cout << "Hit in collection plane" << std::endl;

        // In reconstruction, the channel waveforms are truncated. So
        // we have to adjust the Hit TDC ticks to match those of the
        // SimChannels, which were created during simulation.
        typedef sim::SimChannel::StoredTDC_t TDC_t;
        TDC_t start_tdc = clock_data.TPCTick2TDC(hit.StartTick());
        TDC_t end_tdc = clock_data.TPCTick2TDC(hit.EndTick());
        TDC_t hitStart_tdc = clock_data.TPCTick2TDC(hit.PeakTime() - 3. * hit.SigmaPeakTime());
        TDC_t hitEnd_tdc = clock_data.TPCTick2TDC(hit.PeakTime() + 3. * hit.SigmaPeakTime());

        start_tdc = std::max(start_tdc, hitStart_tdc);
        end_tdc = std::min(end_tdc, hitEnd_tdc);

        // Search the SimChannels for matching channel number, then look
        // at the tracks inside the channel.
        for (auto const& channel : (*simChannelHandle)) {
          auto simChannelNumber = channel.Channel();
          if (simChannelNumber != hitChannelNumber) continue;

          // std::cout << "Matched channel no.: " << simChannelNumber << std::endl;

          // The time slices in this sim channel.
          auto const& timeSlices = channel.TDCIDEMap();

          // Find the range of time slices in the sim channel that correspond to the range of hit times.
          sim::TDCIDE startTime;
          sim::TDCIDE endTime;
          startTime.first = start_tdc;
          endTime.first = end_tdc;
          // Find a pointer to the first time slice with time >= start_tdc.
          auto const startPointer = std::lower_bound(timeSlices.begin(), timeSlices.end(), startTime, TDCIDETimeCompare);
          // Find the last time slice with time < end_tdc.
          auto const endPointer = std::upper_bound(startPointer, timeSlices.end(), endTime, TDCIDETimeCompare);

          // Did we find anything? If not, skip.
          if (startPointer == timeSlices.end() || startPointer == endPointer) continue;
          // std::cout << "Time slice start = " << (*startPointer).first << ", end = " << (*endPointer).first  << std::endl;

          // Loop over the time slices we found that match the hit times.
          for (auto slicePointer = startPointer; slicePointer != endPointer; slicePointer++) {
            auto const timeSlice = *slicePointer;

            // Loop over the energy deposits from sim channel.
            auto const& energyDeposits = timeSlice.second;
            for (auto const& energyDeposit : energyDeposits) {

              // Search the map for the track ID associated with this energy deposit.
              auto search = particleMap.find(energyDeposit.trackID);

              // Did we find this track ID in the particle map?
              if (search == particleMap.end()) continue;

              // "search" points to a pair in the map: <track ID, MCParticle*>
              const simb::MCParticle& particle = *((*search).second);

              // muon
              if ( particle.Process() != "primary" || abs(particle.PdgCode()) != 13 ) continue;

            } // loop over energy deposits
          }   // loop over time slices
        }     // for each SimChannel
      }       // for each Hit

      // In general, objects in the LArSoft reconstruction chain are linked using the art::Assns class:
      // <https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/Using_art_in_LArSoft#artAssns>

      // The following statement will find the simb::MCTruth associated with the simb::MCParticle
      const art::FindManyP<simb::MCTruth> findManyTruth(particleHandle, event, fSimulationProducerLabel);

      if (!findManyTruth.isValid()) {
        std::cout << "findManyTruth simb::MCTruth for simb::MCParticle failed!" << std::endl;
      }

      size_t particle_index = 0; // only look at first particle in particleHandle's vector.
      auto const& truth = findManyTruth.at(particle_index);

      // Make sure there's no problem.
      if (truth.empty()) {
        std::cout << "Particle ID=" << particleHandle->at(particle_index).TrackId() << " has no primary!" << std::endl;
      }

      // std::cout << "Sec. particle ID=" << particleHandle->at(particle_index).TrackId() << " <-> primary gen PDG=" << truth[0]->GetParticle(0).PdgCode() << std::endl;

      // Example to find what hits are associated with clusters.
      art::Handle<std::vector<recob::Cluster>> clusterHandle;
      if (!event.getByLabel(fClusterProducerLabel, clusterHandle)) return;

      // Note that this is not as trivial a query as the one with MCTruth, since it's
      // possible for a hit to be assigned to more than one cluster.
      const art::FindManyP<recob::Hit> findManyHits(clusterHandle, event, fClusterProducerLabel);

      if (!findManyHits.isValid()) {
        std::cout << "findManyHits recob::Hit for recob::Cluster failed;" << " cluster label='" << fClusterProducerLabel << "'"<< std::endl;
      }

      // Loop over the clusters to see the hits associated with each one.
      /*
      for (size_t cluster_index = 0; cluster_index != clusterHandle->size(); cluster_index++) {

        auto const& hits = findManyHits.at(cluster_index);

        // A vector of pointers to the hits associated with the cluster,
        std::cout << "Cluster ID=" << clusterHandle->at(cluster_index).ID() << " has " << hits.size() << " hits" << std::endl;
      }*/

      fNtuple->Fill();

    } // MyEnergyAnalysis::analyze()

    // This macro has to be defined for this module to be invoked from a
    // .fcl file; see MyEnergyAnalysis.fcl for more information.
    DEFINE_ART_MODULE(MyEnergyAnalysis)

  } // namespace example
} // namespace lar

// Back to our local namespace.
namespace {

  double DetectorDiagonal(geo::GeometryCore const& geom)
  {
    const double length = geom.DetLength();
    const double width = 2. * geom.DetHalfWidth();
    const double height = 2. * geom.DetHalfHeight();

    return std::sqrt(cet::sum_of_squares(length, width, height));
  }

  bool MomentumOrderMCParticle(const simb::MCParticle* p1, const simb::MCParticle* p2) {
    return ( p1->P(0) > p2->P(0) );
  }

  bool TDCIDETimeCompare(const sim::TDCIDE& lhs, const sim::TDCIDE& rhs)
  {
    return lhs.first < rhs.first;
  }

} // local namespace
