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

  // Ancestor Mother is primary lepton
  bool IsAncestorMotherPrimaryLep(const simb::MCParticle&, int, std::map<int, const simb::MCParticle*>);

} // local namespace

// An outside package call this module like lar::example::MyEnergyAnalysis

namespace lar {
  namespace example {

    // BEGIN MyEnergyAnalysis group
    // -----------------------------------------------
    // class definition
    //
    // This class produces a ROOT tree that contains information
    // from the generated/simulated and reconstructed particles.
    //
    // Configuration parameters
    // =========================
    //
    // - GenieGenModuleLabel (string, default: "generator"): tag of the input data
    //   product with the event generator information
    //
    // - SimulationLabel (string, default: "largeant"): tag of the input data
    //   product with the detector simulation information (typically an instance
    //   of the LArG4 module)
    //
    class MyEnergyAnalysis : public art::EDAnalyzer {
    public:

      // This structure describes the configuration parameters of the module.
      // Any missing or unknown parameters will generate a configuration error.

      struct Config {

        // Save some typing:
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        // One Atom for each parameter
        fhicl::Atom<art::InputTag> GenieGenModuleLabel{
          Name("GenieGenModuleLabel"),
          Comment("tag of the input data product with the event generator "
                  "information")};

        fhicl::Atom<art::InputTag> SimulationLabel{
          Name("SimulationLabel"),
          Comment("tag of the input data product with the detector simulation "
                  "information")};

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

      // The parameters we will read from the .fcl file.
      art::InputTag fGenieGenModuleLabel;     // The name of the producer that generated particles e.g. GENIE
      art::InputTag fSimulationProducerLabel; // The name of the producer that tracked simulated particles through the detector

      // The n-tuple to create
      TTree* fNtuple;

      // Event info
      int fEvent;  // number of the event being processed
      int fRun;    // number of the run being processed
      int fSubRun; // number of the sub-run being processed

      // Add true nu information
      double eP, eN, ePip, ePim, ePi0, eOther;    // Energy of particles
      int nP, nN, nPip, nPim, nPi0, nOther;                            // number of particles
      double E_vis_true;                 // True vis energy [GeV]

      //
      // Variables related to geneator/simulation
      //
      int fSimPDG;                       // MCParticle PDG ID
      std::vector<int> fSimP_TrackID_vec;
      std::vector<int> EDep_TrackID_vec;
      std::vector<int> fSimP_PDG_vec;
      std::vector<int> fSimP_Mom_vec;
      std::vector<int> fSimP_SC_vec;
      std::vector<float> fSimP_vtx_x_vec;
      std::vector<float> fSimP_vtx_y_vec;
      std::vector<float> fSimP_vtx_z_vec;
      std::vector<float> fSimP_ptot_vec;
      std::vector<float> fSimP_px_vec;
      std::vector<float> fSimP_py_vec;
      std::vector<float> fSimP_pz_vec;
      std::vector<float> fSimP_E_vec;
      std::vector<float> fSimP_M_vec;
      std::vector<float> fSimP_Ek_vec;

      int fSimTrackID;                   // GEANT ID of the particle being processed
      int EDepTrackID;
      int primarylep_trkID;
      int fSim_nEle;                     // No. of Sim electrons (e+/e-)
      int fSim_nNue;                     // No. of Sim electron neutrinos (nue and nuebar)
      int fSim_nMu;                      // No. of Sim muons (mu+/mu-)
      int fSim_nNumu;                    // No. of Sim muon neutrinos (numu and numubar)
      int fSim_nTau;                     // No. of Sim tau leptons (+/-)
      int fSim_nNutau;                   // No. of Sim tau neutrinos (nutau and nutaubar)
      int fSim_nPhoton;                  // No. of Sim photons
      int fSim_nPionNeutral;             // No. of Sim pi+/pi-
      int fSim_nPionCharged;             // No. of Sim pi0
      int fSim_nNeutron;                 // No. of Sim neutrons
      int fSim_nProton;                  // No. of Sim protons
      double fSim_LepE, fSim_HadE;       // Energy of Sim lep and had

      int fCCNC_truth;     		 //0=CC 1=NC
      int fMode_truth;     		 //0=QE/El, 1=RES, 2=DIS, 3=Coherent production
      int fInteractionType; 		 // Interaction type
      double fNuvtxx_truth; 		 //Genie true neutrino interaction vertex x
      double fNuvtxy_truth; 		 //Genie true neutrino interaction vertex y
      double fNuvtxz_truth;  		 //Genie true neutrino interaction vertex z
      int fNuPDG;                        // Generator level neutrino PDG code
      int fLepPDG;                       // Generator level outgoing lepton PDG code
      double fLepMomX, fLepMomY, fLepMomZ;      // Generator level outgoing lepton momentum
      double fLepvtx_x, fLepvtx_y, fLepvtx_z;               // Generator level outgoing lepton vtx
      double fVis_LepE;                      // Generator level neutrino lepton energy [GeV]
      double fLepMass;                      // Generator level neutrino lepton mass [GeV]
      int fStatusCode;                    // Generator level neutrino lepton statuscode
      double fLepNuAngle;                // Angle b/w nu and lepton



      double fGen_numu_E;                // Energy of generator level neutrino [GeV]
      double fSim_numu_E;                // Energy of leading muon (anti) neutrino
      double fSim_mu_start_vx;           // x position of the muon trajectory start
      double fSim_mu_start_vy;           // y .....................................
      double fSim_mu_start_vz;           // z .....................................
      double fSim_mu_start_4position[4]; // (x,y,z,t) of the muon trajectory start
      double fSim_mu_end_vx;             // x position of the muon trajectory end
      double fSim_mu_end_vy;             // y ...................................
      double fSim_mu_end_vz;             // z ...................................
      double fSim_mu_end_4position[4];   // ................................ end
      double fSim_mu_start_px;           // x momentum of the muon trajectory start
      double fSim_mu_start_py;           // y .....................................
      double fSim_mu_start_pz;           // z .....................................
      double fSim_mu_start_E;            // Energy of leading mu
      double fSim_mu_start_4mommenta[4]; // (Px,Py,Pz,E) of the muon trajectory start
      double fSim_mu_end_px;             // x momentum of the muon trajectory end
      double fSim_mu_end_py;             // y ...................................
      double fSim_mu_end_pz;             // z ...................................
      double fSim_mu_end_E;              // Energy of leading mu
      double fSim_mu_end_4mommenta[4];   // ................................... end
      double fSim_mu_track_length;       // leading mu track length

      /*
      double fSim_mu_Edep_a1;                // muon energy deposit [GeV]: total amount of electrons reaching the readout channel
      double fSim_mu_Edep_a2;                // muon energy deposit [MeV]: total amount of energy released by ionizations in the event (from Geant4 simulation)
      double fSim_mu_Edep_NonCollectionPlane_a2; // [MeV]
      double fSim_mu_Edep_a2_debug;          // [MeV]
      double fSim_mu_Edep_b2_debug;          // [MeV]
      double fSim_mu_Edep_b1;                // [GeV]
      double fSim_mu_Edep_NonCollectionPlane_b2; // [MeV]
      */
      double fSim_mu_Edep_b2;                // [MeV]

      // Two ways (a, b) to access collection plane +
      // Two ways (1, 2) of get E deposit for sim::IDE
      /*
      // Method a
      double fSim_hadronic_Edep_a1;      // total amount of electrons reaching the readout channel [GeV]
      double fSim_hadronic_Edep_a2;      // total amount of energy released by ionizations in the event (from Geant4 simulation) [MeV]
      double fSim_hadronic_Edep_NonCollectionPlane_a2;  // [MeV]
      double fSim_hadronic_Edep_a2_debug;          // [MeV]
      int fSim_n_hadronic_Edep_a;        // Number of hadronic energy deposits
      std::vector<float> fSim_hadronic_hit_x_a;      // Store position x for each energy deposit
      std::vector<float> fSim_hadronic_hit_y_a;
      std::vector<float> fSim_hadronic_hit_z_a;
      std::vector<float> fSim_hadronic_hit_Edep_a1;
      std::vector<float> fSim_hadronic_hit_Edep_a2;
      */
      // Method b
      //double fSim_hadronic_Edep_b1;
      double fSim_hadronic_Edep_b2;
      //double fSim_hadronic_Edep_NonCollectionPlane_b2;  // [MeV]
      //double fSim_hadronic_Edep_b2_debug;          // [MeV]
      int fSim_n_hadronic_Edep_b;        // Number of hadronic energy deposits
      std::vector<float> fSim_hadronic_hit_x_b;
      std::vector<float> fSim_hadronic_hit_y_b;
      std::vector<float> fSim_hadronic_hit_z_b;
      //std::vector<float> fSim_hadronic_hit_Edep_b1;
      std::vector<float> fSim_hadronic_hit_Edep_b2;

      //
      // Other variables that will be shared between different methods.
      //
      geo::GeometryCore const* fGeometryService; // pointer to Geometry provider
      double fElectronsToGeV;                    // conversion factor for no. of ionization electrons to energy deposited in GeV

      // True info for each particle generated
      std::vector<int> fP_PDG;                         // PDG code for each particle
      std::vector<int> fP_TrackID;                     // TrackID for each particle
      int fP_num;                         // Number of types of particle
      std::vector<int> fP_StatusCode;                  // Status code for each particle, https://internal.dunescience.org/doxygen/GENIEGen__module_8cc_source.html
      std::vector<float> fP_vtx_x;                    // Position: x component for each particle
      std::vector<float> fP_vtx_y;                    // Position: y component for each particle
      std::vector<float> fP_vtx_z;                    // Position: z component for each particle
      std::vector<float> fP_ptot;                     // Total momentum for each particle
      std::vector<float> fP_px;                       // Momentum: x component for each particle
      std::vector<float> fP_py;                       // Momentum: y component for each particle
      std::vector<float> fP_pz;                       // Momentum: z component for each particle
      std::vector<float> fP_E;                        // Energy for each particle [GeV]
      std::vector<float> fP_mass;                     // Mass for each particle [GeV/c^2]
      std::vector<float> fP_Ek;                       // Kinetic Energy for each particle [GeV]
      std::vector<int>   fP_mother;                    // Find the parent of the produced particle. -1 means this particle has no mother

      // True info for energy
      double fTrue_HadE;                              // True had E by adding all fP_E (!=lepton)
      double fTrue_LepE;                              // True Lep E by adding all fP_E (==lepton)
      double fVis_HadE;                               // Visible had E

    }; // class MyEnergyAnalysis

    // END MyEnergyAnalysis group
    // -------------------------------------------------

    //-----------------------------------------------------------------------
    // class implementation

    //-----------------------------------------------------------------------
    // Constructor

    MyEnergyAnalysis::MyEnergyAnalysis(Parameters const& config)
      : EDAnalyzer(config)
      , fGenieGenModuleLabel(config().GenieGenModuleLabel())
      , fSimulationProducerLabel(config().SimulationLabel())
    {
      // Get a pointer to the geometry service provider.
      fGeometryService = lar::providerFrom<geo::Geometry>();

      // Tell beforehand all the data the module is going to read ("consumes") or
      // might read ("may_consume").
      consumes<std::vector<simb::MCTruth>>(fGenieGenModuleLabel);
      consumes<std::vector<simb::MCParticle>>(fSimulationProducerLabel);
      consumes<std::vector<sim::SimChannel>>(fSimulationProducerLabel);
      consumes<art::Assns<simb::MCTruth, simb::MCParticle>>(fSimulationProducerLabel);
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

      fNtuple->Branch("Event",                    &fEvent,                  "Event/I");
      fNtuple->Branch("SubRun",                   &fSubRun,                 "SubRun/I");
      fNtuple->Branch("Run",                      &fRun,                    "Run/I");
      // Add true nu information
      fNtuple->Branch("Vis_LepE",          &fVis_LepE,         "Vis_LepE/D");
      fNtuple->Branch("LepMass",           &fLepMass,          "LepMass/D");

      fNtuple->Branch("eP",            &eP,            "eP/D");
      fNtuple->Branch("eN",            &eN,            "eN/D");
      fNtuple->Branch("ePip",          &ePip,          "ePip/D");
      fNtuple->Branch("ePim",          &ePim,          "ePim/D");
      fNtuple->Branch("ePi0",          &ePi0,          "ePi0/D");
      fNtuple->Branch("eOther",        &eOther,        "eOther/D");
      fNtuple->Branch("nP",            &nP,            "nP/I");
      fNtuple->Branch("nN",            &nN,            "nN/I");
      fNtuple->Branch("nPip",          &nPip,          "nPip/I");
      fNtuple->Branch("nPim",          &nPim,          "nPim/I");
      fNtuple->Branch("nPi0",          &nPi0,          "nPi0/I");
      fNtuple->Branch("nOther",        &nOther,        "nOther/D");
      fNtuple->Branch("E_vis_true",    &E_vis_true,    "E_vis_true/D");

      // GEN neutrino E
      fNtuple->Branch("Gen_numu_E",               &fGen_numu_E,             "Gen_numu_E/D");
      fNtuple->Branch("CCNC_truth",               &fCCNC_truth,             "CCNC_truth/I");
      fNtuple->Branch("Mode_truth",               &fMode_truth,             "Mode_truth/I");
      fNtuple->Branch("InteractionType",          &fInteractionType,        "InteractionType/I");
      fNtuple->Branch("Nuvtxx_truth",             &fNuvtxx_truth,           "Nuvtxx_truth/D");
      fNtuple->Branch("Nuvtxy_truth",             &fNuvtxy_truth,           "Nuvtxy_truth/D");
      fNtuple->Branch("Nuvtxz_truth",             &fNuvtxz_truth,           "Nuvtxz_truth/D");
      // Generator level PDG code
      fNtuple->Branch("LepPDG",        	          &fLepPDG,     	          "LepPDG/I");
      fNtuple->Branch("neuPDG",         	        &fNuPDG,             	    "neuPDG/I");
      fNtuple->Branch("LepNuAngle",               &fLepNuAngle,             "LepNuAngle/D");
      fNtuple->Branch("LepMomX",                  &fLepMomX,                "LepMomX/D");
      fNtuple->Branch("LepMomY",                  &fLepMomY,                "LepMomY/D");
      fNtuple->Branch("LepMomZ",                  &fLepMomZ,                "LepMomZ/D");
      fNtuple->Branch("Lepvtx_x",                 &fLepvtx_x,               "Lepvtx_x/D");
      fNtuple->Branch("Lepvtx_y",                 &fLepvtx_y,               "Lepvtx_y/D");
      fNtuple->Branch("Lepvtx_z",                 &fLepvtx_z,               "Lepvtx_z/D");
      fNtuple->Branch("StatusCode",         	    &fStatusCode,             "StatusCode/I");

      // Simulation branches Sim*
      fNtuple->Branch("SimP_TrackID_vec",              &fSimP_TrackID_vec);
      fNtuple->Branch("SimP_PDG_vec",                  &fSimP_PDG_vec);
      fNtuple->Branch("SimP_Mom_vec",                  &fSimP_Mom_vec);
      fNtuple->Branch("SimP_SC_vec",                   &fSimP_SC_vec);
      fNtuple->Branch("SimP_vtx_x_vec",                &fSimP_vtx_x_vec);
      fNtuple->Branch("SimP_vtx_y_vec",                &fSimP_vtx_y_vec);
      fNtuple->Branch("SimP_vtx_z_vec",                &fSimP_vtx_z_vec);
      fNtuple->Branch("SimP_ptot_vec",                 &fSimP_ptot_vec);
      fNtuple->Branch("SimP_px_vec",                   &fSimP_px_vec);
      fNtuple->Branch("SimP_py_vec",                   &fSimP_py_vec);
      fNtuple->Branch("SimP_pz_vec",                   &fSimP_pz_vec);
      fNtuple->Branch("SimP_E_vec",                    &fSimP_E_vec);
      fNtuple->Branch("SimP_M_vec",                    &fSimP_M_vec);
      fNtuple->Branch("SimP_Ek_vec",                   &fSimP_Ek_vec);

      fNtuple->Branch("Sim_nEle",                 &fSim_nEle,               "Sim_nEle/I");
      fNtuple->Branch("Sim_nNue",                 &fSim_nNue,               "Sim_nNue/I");
      fNtuple->Branch("Sim_nMu",                  &fSim_nMu,                "Sim_nMu/I");
      fNtuple->Branch("Sim_nNumu",                &fSim_nNumu,              "Sim_nNumu/I");
      fNtuple->Branch("Sim_nTau",                 &fSim_nTau,               "Sim_nTau/I");
      fNtuple->Branch("Sim_nNutau",               &fSim_nNutau,             "Sim_nNutau/I");
      fNtuple->Branch("Sim_nPhoton",              &fSim_nPhoton,            "Sim_nPhoton/I");
      fNtuple->Branch("Sim_nPionNeutral",         &fSim_nPionNeutral,       "Sim_nPionNeutral/I");
      fNtuple->Branch("Sim_nPionCharged",         &fSim_nPionCharged,       "Sim_nPionCharged/I");
      fNtuple->Branch("Sim_nNeutron",             &fSim_nNeutron,           "Sim_nNeutron/I");
      fNtuple->Branch("Sim_nProton",              &fSim_nProton,            "Sim_nProton/I");
      fNtuple->Branch("Sim_LepE",                 &fSim_LepE,               "Sim_LepE/D");
      fNtuple->Branch("Sim_HadE",                 &fSim_HadE,               "Sim_HadE/D");

      // GEANT level neutrino E
      fNtuple->Branch("Sim_numu_E",               &fSim_numu_E,             "Sim_numu_E/D");
      // muon position
      fNtuple->Branch("Sim_mu_start_vx",          &fSim_mu_start_vx,        "Sim_mu_start_vx/D");
      fNtuple->Branch("Sim_mu_start_vy",          &fSim_mu_start_vy,        "Sim_mu_start_vy/D");
      fNtuple->Branch("Sim_mu_start_vz",          &fSim_mu_start_vz,        "Sim_mu_start_vz/D");
      fNtuple->Branch("Sim_mu_end_vx",            &fSim_mu_end_vx,          "Sim_mu_end_vx/D");
      fNtuple->Branch("Sim_mu_end_vy",            &fSim_mu_end_vy,          "Sim_mu_end_vy/D");
      fNtuple->Branch("Sim_mu_end_vz",            &fSim_mu_end_vz,          "Sim_mu_end_vz/D");
      // muon momentum
      fNtuple->Branch("Sim_mu_start_px",          &fSim_mu_start_px,        "Sim_mu_start_px/D");
      fNtuple->Branch("Sim_mu_start_py",          &fSim_mu_start_py,        "Sim_mu_start_py/D");
      fNtuple->Branch("Sim_mu_start_pz",          &fSim_mu_start_pz,        "Sim_mu_start_pz/D");
      fNtuple->Branch("Sim_mu_start_E",           &fSim_mu_start_E,         "Sim_mu_start_E/D");
      fNtuple->Branch("Sim_mu_end_px",            &fSim_mu_end_px,          "Sim_mu_end_px/D");
      fNtuple->Branch("Sim_mu_end_py",            &fSim_mu_end_py,          "Sim_mu_end_py/D");
      fNtuple->Branch("Sim_mu_end_pz",            &fSim_mu_end_pz,          "Sim_mu_end_pz/D");
      fNtuple->Branch("Sim_mu_end_E",             &fSim_mu_end_E,           "Sim_mu_end_E/D");
      fNtuple->Branch("Sim_mu_start_4position",    fSim_mu_start_4position, "Sim_mu_start_4position[4]/D");
      fNtuple->Branch("Sim_mu_end_4position",      fSim_mu_end_4position,   "Sim_mu_end_4position[4]/D");
      fNtuple->Branch("Sim_mu_start_4mommenta",    fSim_mu_start_4mommenta, "Sim_mu_start_4mommenta[4]/D");
      fNtuple->Branch("Sim_mu_end_4mommenta",      fSim_mu_end_4mommenta,   "Sim_mu_end_4mommenta[4]/D");
      fNtuple->Branch("Sim_mu_track_length",      &fSim_mu_track_length,    "Sim_mu_track_length/D");
      fNtuple->Branch("Sim_mu_Edep_b2",           &fSim_mu_Edep_b2,         "Sim_mu_Edep_b2/D");
      /*
      fNtuple->Branch("Sim_mu_Edep_a1",           &fSim_mu_Edep_a1,            "Sim_mu_Edep_a1/D");
      fNtuple->Branch("Sim_mu_Edep_a2",           &fSim_mu_Edep_a2,            "Sim_mu_Edep_a2/D");
      fNtuple->Branch("Sim_mu_Edep_NonCollectionPlane_a2",           &fSim_mu_Edep_NonCollectionPlane_a2,            "fSim_mu_Edep_NonCollectionPlane_a2/D");
      fNtuple->Branch("Sim_mu_Edep_b1",           &fSim_mu_Edep_b1,            "Sim_mu_Edep_b1/D");
      fNtuple->Branch("Sim_mu_Edep_NonCollectionPlane_b2",           &fSim_mu_Edep_NonCollectionPlane_b2,            "fSim_mu_Edep_NonCollectionPlane_b2/D");
      fNtuple->Branch("Sim_mu_Edep_a2_debug",           &fSim_mu_Edep_a2_debug,            "fSim_mu_Edep_a2_debug/D");
      fNtuple->Branch("Sim_mu_Edep_b2_debug",           &fSim_mu_Edep_b2_debug,            "fSim_mu_Edep_b2_debug/D");

      fNtuple->Branch("Sim_hadronic_Edep_a1",     &fSim_hadronic_Edep_a1,   "Sim_hadronic_Edep_a1/D");
      fNtuple->Branch("Sim_hadronic_Edep_a2",     &fSim_hadronic_Edep_a2,   "Sim_hadronic_Edep_a2/D");
      fNtuple->Branch("Sim_hadronic_Edep_NonCollectionPlane_a2",     &fSim_hadronic_Edep_NonCollectionPlane_a2,   "Sim_hadronic_Edep_NonCollectionPlane_a2/D");
      fNtuple->Branch("Sim_n_hadronic_Edep_a",    &fSim_n_hadronic_Edep_a,  "Sim_n_hadronic_Edep_a/I");
      fNtuple->Branch("Sim_hadronic_hit_x_a",     &fSim_hadronic_hit_x_a);
      fNtuple->Branch("Sim_hadronic_hit_y_a",     &fSim_hadronic_hit_y_a);
      fNtuple->Branch("Sim_hadronic_hit_z_a",     &fSim_hadronic_hit_z_a);
      fNtuple->Branch("Sim_hadronic_hit_Edep_a1", &fSim_hadronic_hit_Edep_a1);
      fNtuple->Branch("Sim_hadronic_hit_Edep_a2", &fSim_hadronic_hit_Edep_a2);

      fNtuple->Branch("Sim_hadronic_Edep_NonCollectionPlane_b2",     &fSim_hadronic_Edep_NonCollectionPlane_b2,   "Sim_hadronic_Edep_NonCollectionPlane_b2/D");
      fNtuple->Branch("Sim_hadronic_Edep_a2_debug",     &fSim_hadronic_Edep_a2_debug,   "Sim_hadronic_Edep_a2_debug/D");
      fNtuple->Branch("Sim_hadronic_Edep_b2_debug",     &fSim_hadronic_Edep_b2_debug,   "Sim_hadronic_Edep_b2_debug/D");
      fNtuple->Branch("Sim_hadronic_hit_Edep_b1", &fSim_hadronic_hit_Edep_b1);
      fNtuple->Branch("Sim_hadronic_Edep_b1",     &fSim_hadronic_Edep_b1,   "Sim_hadronic_Edep_b1/D");
      */

      fNtuple->Branch("Sim_hadronic_Edep_b2",     &fSim_hadronic_Edep_b2,   "Sim_hadronic_Edep_b2/D");
      fNtuple->Branch("Sim_n_hadronic_Edep_b",    &fSim_n_hadronic_Edep_b,  "Sim_n_hadronic_Edep_b/I");
      fNtuple->Branch("Sim_hadronic_hit_x_b",     &fSim_hadronic_hit_x_b);
      fNtuple->Branch("Sim_hadronic_hit_y_b",     &fSim_hadronic_hit_y_b);
      fNtuple->Branch("Sim_hadronic_hit_z_b",     &fSim_hadronic_hit_z_b);
      fNtuple->Branch("Sim_hadronic_hit_Edep_b2", &fSim_hadronic_hit_Edep_b2);

      // True info for each particle
      fNtuple->Branch("P_num",                    &fP_num,                "P_num/I");
      fNtuple->Branch("P_mother",                 &fP_mother);
      fNtuple->Branch("P_TrackID",                &fP_TrackID);
      fNtuple->Branch("P_PDG",        	          &fP_PDG);
      fNtuple->Branch("P_StatusCode",             &fP_StatusCode);
      fNtuple->Branch("P_vtx_x",                  &fP_vtx_x);
      fNtuple->Branch("P_vtx_y",                  &fP_vtx_y);
      fNtuple->Branch("P_vtx_z",                  &fP_vtx_z);
      fNtuple->Branch("P_ptot",                   &fP_ptot);
      fNtuple->Branch("P_px",                     &fP_px);
      fNtuple->Branch("P_py",                     &fP_py);
      fNtuple->Branch("P_pz",                     &fP_pz);
      fNtuple->Branch("P_E",                      &fP_E);
      fNtuple->Branch("P_mass",                   &fP_mass);
      fNtuple->Branch("P_Ek",                     &fP_Ek);

      // Reconstruction branches
      fNtuple->Branch("True_HadE",                &fTrue_HadE,            "True_HadE/D");
      fNtuple->Branch("True_LepE",                &fTrue_LepE,            "True_LepE/D");
      fNtuple->Branch("Vis_HadE",                 &fVis_HadE,             "Vis_HadE/D");

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

      // Initialize
      fGen_numu_E                = 0.;
      fCCNC_truth     		       = -9999.;
      fMode_truth     		       = -9999.;
      fInteractionType 		       = -9999.;
      fNuvtxx_truth 		         = -9999.;
      fNuvtxy_truth         		 = -9999.;
      fNuvtxz_truth  	        	 = -9999.;
      fSim_numu_E                = 0.;
      fSim_mu_start_vx           = -9999.;
      fSim_mu_start_vy           = -9999.;
      fSim_mu_start_vz           = -9999.;
      fSim_mu_end_vx             = -9999.;
      fSim_mu_end_vy             = -9999.;
      fSim_mu_end_vz             = -9999.;
      fSim_mu_start_px           = -9999.;
      fSim_mu_start_py           = -9999.;
      fSim_mu_start_pz           = -9999.;
      fSim_mu_start_E            = -9999.;
      fSim_mu_end_px             = -9999.;
      fSim_mu_end_py             = -9999.;
      fSim_mu_end_pz             = -9999.;
      fSim_mu_end_E              = -9999.;
      fSim_mu_track_length       = -9999.;
      fSim_LepE                  = 0.;
      fSim_HadE                  = 0.;

      // Initialize true info
      fLepNuAngle = -9999.;
      fLepMomX    = -9999.;
      fLepMomY    = -9999.;
      fLepMomZ    = -9999.;
      fLepvtx_x    = -9999.;
      fLepvtx_y    = -9999.;
      fLepvtx_z    = -9999.;
      fVis_LepE    = -9999.;
      fLepMass     = -9999.;

      fP_num        = 0;
      fP_PDG.clear();
      fP_mother.clear();
      fP_TrackID.clear();
      fP_StatusCode.clear();
      fP_vtx_x.clear();
      fP_vtx_y.clear();
      fP_vtx_z.clear();
      fP_ptot.clear();
      fP_px.clear();
      fP_py.clear();
      fP_pz.clear();
      fP_E.clear();
      fP_mass.clear();
      fP_Ek.clear();

      fSimP_TrackID_vec.clear();
      EDep_TrackID_vec.clear();
      fSimP_PDG_vec.clear();
      fSimP_Mom_vec.clear();
      fSimP_SC_vec.clear();
      fSimP_vtx_x_vec.clear();
      fSimP_vtx_y_vec.clear();
      fSimP_vtx_z_vec.clear();
      fSimP_ptot_vec.clear();
      fSimP_px_vec.clear();
      fSimP_py_vec.clear();
      fSimP_pz_vec.clear();
      fSimP_E_vec.clear();
      fSimP_M_vec.clear();
      fSimP_Ek_vec.clear();


      // Initialize mu deposit energy
      /*
      fSim_mu_Edep_a1 = 0.;                // muon energy deposit [GeV]: total amount of electrons reaching the readout channel
      fSim_mu_Edep_a2 = 0.;                // muon energy deposit [MeV]: total amount of energy released by ionizations in the event (from Geant4 simulation)
      fSim_mu_Edep_NonCollectionPlane_a2 = 0.;
      fSim_mu_Edep_b1 = 0.;
      fSim_mu_Edep_NonCollectionPlane_b2 = 0.;
      fSim_mu_Edep_a2_debug = 0.;
      fSim_mu_Edep_b2_debug = 0.;
      fSim_hadronic_Edep_a1      = 0.; // This initilization is necessary
      fSim_hadronic_Edep_a2      = 0.;
      fSim_hadronic_Edep_NonCollectionPlane_a2 = 0.;
      fSim_hadronic_Edep_b1      = 0.;
      fSim_hadronic_Edep_NonCollectionPlane_b2 = 0.;
      fSim_hadronic_Edep_a2_debug = 0.;
      fSim_hadronic_Edep_b2_debug = 0.;
      */
      fSim_mu_Edep_b2        = 0.;
      fSim_hadronic_Edep_b2  = 0.;

      for (int i = 0; i < 4; i++) {
        fSim_mu_start_4position[i] = -9999.;
        fSim_mu_end_4position[i]   = -9999.;
        fSim_mu_start_4mommenta[i] = -9999.;
        fSim_mu_end_4mommenta[i]   = -9999.;
      }
      /*
      fSim_hadronic_hit_x_a.clear();
      fSim_hadronic_hit_y_a.clear();
      fSim_hadronic_hit_z_a.clear();
      fSim_hadronic_hit_Edep_a1.clear();
      fSim_hadronic_hit_Edep_a2.clear();

      fSim_hadronic_hit_Edep_b1.clear();
      */
      fSim_hadronic_hit_x_b.clear();
      fSim_hadronic_hit_y_b.clear();
      fSim_hadronic_hit_z_b.clear();

      fSim_hadronic_hit_Edep_b2.clear();

      // LArSoft data products: https://larsoft.org/important-concepts-in-larsoft/data-products/

      //
      // Process generator level info
      //

      // c.f. https://github.com/DUNE/dunetpc/blob/master/dune/FDSensOpt/CAFMaker_module.cc#L720
      //      https://github.com/DUNE/dunetpc/blob/master/dune/FDSensOpt/NueAna_module.cc#L639
      art::Handle<std::vector<simb::MCTruth>> mctruthListHandle; // Generator level truth
      std::vector<art::Ptr<simb::MCTruth>> mclist;
      if ( event.getByLabel(fGenieGenModuleLabel, mctruthListHandle) ) art::fill_ptr_vector(mclist, mctruthListHandle);

      // There could be more than one MCTruth, e.g., you might have multiple neutrino interactions per spill,
      // in which case you’d run GENIE multiple times and have one MCTruth per interaction.
      // Or you might want one MCTruth information for the GENIE event and another that overlays cosmic simulation or data onto the same event
      if ( mclist.size() )
      {
        fGen_numu_E = mclist[0]->GetNeutrino().Nu().E(); // true neutrino energy
        fCCNC_truth   = mclist[0]->GetNeutrino().CCNC(); // CC or NC interaction
      	fMode_truth   = mclist[0]->GetNeutrino().Mode(); // Interaction mode (QE/1-pi/DIS...)
        fInteractionType = mclist[0]->GetNeutrino().InteractionType(); // Interaction type
        fNuvtxx_truth = mclist[0]->GetNeutrino().Nu().Vx(); //Genie true neutrino interaction vertex x
	      fNuvtxy_truth = mclist[0]->GetNeutrino().Nu().Vy(); //Genie true neutrino interaction vertex y
      	fNuvtxz_truth = mclist[0]->GetNeutrino().Nu().Vz(); //Genie true neutrino interaction vertex z
      	fNuPDG    = mclist[0]->GetNeutrino().Nu().PdgCode(); // Generator level neutrino PDG code
	      fLepPDG     = mclist[0]->GetNeutrino().Lepton().PdgCode(); // Generator level lepton PDG code
        fLepMomX    = mclist[0]->GetNeutrino().Lepton().Momentum().X(); // Generator level lepton momentum x
        fLepMomY    = mclist[0]->GetNeutrino().Lepton().Momentum().Y();  // Generator level lepton momentum y
        fLepMomZ    = mclist[0]->GetNeutrino().Lepton().Momentum().Z();  // Generator level lepton momentum z
        fLepvtx_x       = mclist[0]->GetNeutrino().Lepton().Vx();   // Generator level lepton vtx x
        fLepvtx_y       = mclist[0]->GetNeutrino().Lepton().Vy();   // Generator level lepton vtx y
        fLepvtx_z       = mclist[0]->GetNeutrino().Lepton().Vz();   // Generator level lepton vtx z
        fLepMass        = mclist[0]->GetNeutrino().Lepton().Mass();
        fVis_LepE       = mclist[0]->GetNeutrino().Lepton().Momentum().T() - fLepMass; // Generator level neutrino lepton kinetic energy
        fStatusCode     = mclist[0]->GetNeutrino().Lepton().StatusCode();  // Generator level neutrino lepton statuscode
        fLepNuAngle = mclist[0]->GetNeutrino().Nu().Momentum().Vect().Angle(mclist[0]->GetNeutrino().Lepton().Momentum().Vect()); // Angle b/w nu and lepton
      }
      // Is evt vtx GetNeutrino().Nu().Vx()?


      // Add true particle counts

      eP = 0.;
      eN = 0.;
      ePip = 0.;
      ePim = 0.;
      ePi0 = 0.;
      eOther = 0.;

      nP = 0;
      nN = 0;
      nPip = 0;
      nPim = 0;
      nPi0 = 0;
      nOther = 0;

      fP_num = mclist[0]->NParticles();
      // std::cout << "fP_num: " << fP_num << "\n\n";

      // Initialize
      fTrue_HadE = 0.;
      fTrue_LepE = 0.;
      fVis_HadE = 0.;
      double proton_mass = 0.93827; // GeV

      // Choose CC event only
      if ( fCCNC_truth ==0 )
      {
        for ( int p = 0; p < mclist[0]->NParticles(); p++ )
        {
          fP_TrackID.push_back(mclist[0]->GetParticle(p).TrackId());
          fP_PDG.push_back(mclist[0]->GetParticle(p).PdgCode());
          fP_mother.push_back(mclist[0]->GetParticle(p).Mother());
          fP_StatusCode.push_back(mclist[0]->GetParticle(p).StatusCode());
          fP_vtx_x.push_back(mclist[0]->GetParticle(p).Vx());
          fP_vtx_y.push_back(mclist[0]->GetParticle(p).Vy());
          fP_vtx_z.push_back(mclist[0]->GetParticle(p).Vz());
          fP_ptot.push_back(mclist[0]->GetParticle(p).P());
          fP_px.push_back(mclist[0]->GetParticle(p).Px());
          fP_py.push_back(mclist[0]->GetParticle(p).Py());
          fP_pz.push_back(mclist[0]->GetParticle(p).Pz());
          fP_E.push_back(mclist[0]->GetParticle(p).E());
          fP_mass.push_back(mclist[0]->GetParticle(p).Mass());
          fP_Ek.push_back(fP_E.at(p) - fP_mass.at(p));

          // Stable Final State
          // The sum of true energy of hadrons and leptons should be true nu energy minus binding energy
          // Paper related to the binding energy: https://link.springer.com/article/10.1140/epjc/s10052-019-6750-3
          // Calculate true vis had E
          if ( fP_StatusCode.at(p) == 1 ) // Stable Final State
          {

            // Claculate true Lep E
            if ( abs(fP_PDG.at(p)) == 13 )
            {
              fTrue_LepE += fP_E.at(p);
            }

            if ( abs(fP_PDG.at(p)) <= 999 && abs(fP_PDG.at(p)) >=100 ) // kPdgMeson
            {
              fTrue_HadE += fP_E.at(p);
            }
            else if (fP_PDG.at(p) == 2212 || fP_PDG.at(p) == 2112) // kPdgProton or kPdgNeutron
            {
              fTrue_HadE += fP_Ek.at(p);
            }
            else if ( fP_PDG.at(p) <= 9999 && fP_PDG.at(p) >= 1000  ) // kPdgBaryon except proton and neutron
            {
              fTrue_HadE += fP_Ek.at(p) + ( fP_mass.at(p) - proton_mass );
            }
            else if ( fP_PDG.at(p) >= -9999 && fP_PDG.at(p) <= -1000  ) // kPdgAntiBaryon except proton and neutron, antihyperon
            {
              fTrue_HadE += fP_Ek.at(p) + 2*fP_mass.at(p) + ( fP_mass.at(p) - proton_mass );
            }
            else if (fP_PDG.at(p) == 22) // kPdgGamma
            {
              fTrue_HadE += fP_E.at(p);
            }



            if ( fP_PDG.at(p) == 2212 ) // kPdgProton
            {
              eP += fP_Ek.at(p);
              nP++;
            }
            else if ( fP_PDG.at(p) == 2112 ) // kPdgNeutron
            {
              eN += fP_Ek.at(p);
              nN++;
            }
            else if ( fP_PDG.at(p) == 211 ) // kPdgPiP
            {
              ePip += fP_Ek.at(p);
              nPip++;
            }
            else if ( fP_PDG.at(p) == -211 ) // kPdgPiM
            {
              ePim += fP_Ek.at(p);
              nPim++;
            }
            else if ( fP_PDG.at(p) == 111 ) // kPdgPi0
            {
              ePi0 += fP_Ek.at(p);
              nPi0++;
            }
            else if ( fP_PDG.at(p) == 321 || fP_PDG.at(p) == -321 || fP_PDG.at(p) == 311 || fP_PDG.at(p) == -311 || fP_PDG.at(p) == 130 || fP_PDG.at(p) == 310 || fP_PDG.at(p) == 22 || (fP_PDG.at(p)>=100 && fP_PDG.at(p)<=9999) || (fP_PDG.at(p)>=-9999 && fP_PDG.at(p)<=-100)) // kPdgKP, kPdgKM, kPdgK0, kPdgAntiK0, kPdgK0L, kPdgK0S, kPdgGamma, IsHadron(pdg)
            {
              eOther += fP_Ek.at(p);
              nOther++;
            }
          } //end kIStHadronInTheNucleus
        } // end mclist[0]->NParticles() loop

        // True visible energy:
        double pi0_mass = 0.134977; // GeV
        fVis_HadE = eP + ePip + ePim + ePi0 + eOther + nPi0 * pi0_mass;
        E_vis_true = fVis_LepE + fVis_HadE; // KE of leptons and hadrons
        // neutron will not deposit, so it cannot be counted in the E_vis_true
        // VisTrue_NDFD = LepE + HadE,
        // HadE = eP + ePip + ePim + ePi0 + (0.135 * nipi0) + eother

      } // end CC events selection

      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      //------------------------------------------------------------------------
      // Get all the simulated channels for the event. These channels
      // include the energy deposited for each simulated track.
      auto simChannelHandle = event.getValidHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);

      // Create a map pf MCParticle to its track ID, to be used for hadronic part later
      std::map<int, const simb::MCParticle*> particleMap;
      // Create a map of energy deposits to its track ID
      std::map<int, double> EDepMap;

      //
      // Process Sim MCparticles info
      //

      art::Handle<std::vector<simb::MCParticle>> particleHandle; // GEANT 4 level truth

      // Then fill the vector with all the objects
      if ( !event.getByLabel(fSimulationProducerLabel, particleHandle) ) {
        // If no MCParticles in an event, throw an exception to force this module to stop.
        throw cet::exception("MyEnergyAnalysis") << " No simb::MCParticle objects in this event - " << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

      // Store specific particles
      std::vector<const simb::MCParticle*> SimElectrons;
      std::vector<const simb::MCParticle*> SimNues;
      std::vector<const simb::MCParticle*> SimMuons;
      std::vector<const simb::MCParticle*> SimNumus;
      std::vector<const simb::MCParticle*> SimTaus;
      std::vector<const simb::MCParticle*> SimNutaus;
      std::vector<const simb::MCParticle*> SimPhotons;
      std::vector<const simb::MCParticle*> SimNeutralPions;
      std::vector<const simb::MCParticle*> SimChargedPions;
      std::vector<const simb::MCParticle*> SimNeutrons;
      std::vector<const simb::MCParticle*> SimProtons;

      // Loop over the list of particles in the event
      // GENIE: primary process; GEANT4: primary+secondary
      for ( auto const& particle : (*particleHandle) ) {

        // For the methods you can call for MCParticle, see ${NUSIMDATA_INC}/nusimdata/SimulationBase/MCParticle.h.
        fSimTrackID = particle.TrackId();
        fSimP_TrackID_vec.push_back(fSimTrackID);
        // Add the address of the MCParticle to the map, with the track ID as the key.
        particleMap[fSimTrackID] = &particle;

        // Only for primary particles in the event
        fSimPDG = particle.PdgCode();
        fSimP_PDG_vec.push_back(fSimPDG);
        fSimP_Mom_vec.push_back(particle.Mother());
        fSimP_SC_vec.push_back(particle.StatusCode());
        fSimP_vtx_x_vec.push_back(particle.Vx());
        fSimP_vtx_y_vec.push_back(particle.Vy());
        fSimP_vtx_z_vec.push_back(particle.Vz());
        fSimP_ptot_vec.push_back(particle.P());
        fSimP_px_vec.push_back(particle.Px());
        fSimP_py_vec.push_back(particle.Py());
        fSimP_pz_vec.push_back(particle.Pz());
        fSimP_E_vec.push_back(particle.E());
        fSimP_M_vec.push_back(particle.Mass());
        fSimP_Ek_vec.push_back(particle.E()-particle.Mass());

        // Take note of primary lepton track id, to be used later
        if ( particle.Process() == "primary" && ( abs(fSimPDG) == 13 || abs(fSimPDG) == 11 || abs(fSimPDG) == 15 ) ) {
          primarylep_trkID = fSimTrackID;
          if (false) std::cout << "primarylep_trkID: " << primarylep_trkID << std::endl; // the primary lep should always have trk id = 1
        }

        // Calculate sim_lepE and sim_hadE
        if ( particle.StatusCode() == 1 )
        {
         // Sim_LepE
         if ( abs(fSimPDG) == 13 ) fSim_LepE += particle.E();
         // Sim_HadE
         if ( abs(fSimPDG) <= 999 && abs(fSimPDG) >=100 ) // kPdgMeson
            {
              fSim_HadE += particle.E();
            }
            else if (fSimPDG == 2212 || fSimPDG == 2112) // kPdgProton or kPdgNeutron
            {
              fSim_HadE += particle.E()-particle.Mass();
            }
            else if ( fSimPDG <= 9999 && fSimPDG >= 1000  ) // kPdgBaryon except proton and neutron
            {
              fSim_HadE += particle.E()-particle.Mass() + ( particle.Mass() - proton_mass );
            }
            else if ( fSimPDG >= -9999 && fSimPDG <= -1000  ) // kPdgAntiBaryon except proton and neutron, antihyperon
            {
              fSim_HadE += particle.E()-particle.Mass() + 2*particle.Mass() + ( particle.Mass() - proton_mass );
            }
            else if (fSimPDG == 22) // kPdgGamma
            {
              fSim_HadE += particle.E();
            }
        }

        if ( particle.Process() == "primary" ) {
          if ( abs(fSimPDG) == 11 )   SimElectrons.push_back(&particle);
          if ( abs(fSimPDG) == 12 )   SimNues.push_back(&particle);
          if ( abs(fSimPDG) == 13 )   SimMuons.push_back(&particle);
          if ( abs(fSimPDG) == 14 )   SimNumus.push_back(&particle);
          if ( abs(fSimPDG) == 15 )   SimTaus.push_back(&particle);
          if ( abs(fSimPDG) == 16 )   SimNutaus.push_back(&particle);
          if ( abs(fSimPDG) == 22 )   SimPhotons.push_back(&particle);
          if ( abs(fSimPDG) == 111 )  SimNeutralPions.push_back(&particle);
          if ( abs(fSimPDG) == 211 )  SimChargedPions.push_back(&particle);
          if ( abs(fSimPDG) == 2112 ) SimNeutrons.push_back(&particle);
          if ( abs(fSimPDG) == 2212 ) SimProtons.push_back(&particle);
        }

      } // end loop over all particles in the event.

      fSim_nEle         = SimElectrons.size();
      fSim_nNue         = SimNues.size();
      fSim_nMu          = SimMuons.size();
      fSim_nNumu        = SimNumus.size();
      fSim_nTau         = SimTaus.size();
      fSim_nNutau       = SimNutaus.size();
      fSim_nPhoton      = SimPhotons.size();
      fSim_nPionNeutral = SimNeutralPions.size();
      fSim_nPionCharged = SimChargedPions.size();
      fSim_nNeutron     = SimNeutrons.size();
      fSim_nProton      = SimProtons.size();

      // If multiple sim particles present in event, sort momentum from high to low
      if ( fSim_nNumu > 1 ) std::sort(SimNumus.begin(), SimNumus.end(), MomentumOrderMCParticle);
      if ( fSim_nMu > 1 )   std::sort(SimMuons.begin(), SimMuons.end(), MomentumOrderMCParticle);

      // Store info for leading E sim numu GEANT 4 level
      if ( fSim_nNumu > 0 ) {
        const simb::MCParticle& leadingnumu = *(SimNumus[0]);
        fSim_numu_E = leadingnumu.E();
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
        fSim_mu_start_vx = leadingmu.Vx(0); // unit?
        fSim_mu_start_vy = leadingmu.Vy(0);
        fSim_mu_start_vz = leadingmu.Vz(0);
        fSim_mu_end_vx   = leadingmu.Vx(last);
        fSim_mu_end_vy   = leadingmu.Vy(last);
        fSim_mu_end_vz   = leadingmu.Vz(last);
        fSim_mu_start_px = leadingmu.Px(0);
        fSim_mu_start_py = leadingmu.Py(0);
        fSim_mu_start_pz = leadingmu.Pz(0);
        fSim_mu_start_E  = leadingmu.E(0);
        fSim_mu_end_px   = leadingmu.Px(last);
        fSim_mu_end_py   = leadingmu.Py(last);
        fSim_mu_end_pz   = leadingmu.Pz(last);
        fSim_mu_end_E    = leadingmu.E(last);

        // Fill arrays with the 4-values.
        positionStart.GetXYZT(fSim_mu_start_4position);
        positionEnd.GetXYZT(fSim_mu_end_4position);
        momentumStart.GetXYZT(fSim_mu_start_4mommenta);
        momentumEnd.GetXYZT(fSim_mu_end_4mommenta);

        // Calculate length using spherical cooridnate system: assume straight track? time negligible?
        const double trackLength = (positionEnd - positionStart).Rho();
        fSim_mu_track_length = trackLength;
      } // End if muon exists

      //
      // Calculate sim hadronic deposit energy
      //

      // Loop over the SimChannel objects in the event to look at the energy deposited by particle's track.
      for ( auto const& channel : (*simChannelHandle) )
      {

        // Get the numeric ID associated with this channel.
        // See methods at https://internal.dunescience.org/doxygen/SimChannel_8h_source.html
        auto const channelNumber = channel.Channel();

        // Each channel has a map inside it that connects a time slice to energy deposits in the detector.
        // The full type of this map is std::map<unsigned short, std::vector<sim::IDE>>; we'll use "auto" here
        auto const& timeSlices = channel.TDCIDEMap();
        for ( auto const& timeSlice : timeSlices ) {

          // For the timeSlices map, the 'first' is a time slice number; The 'second' is a vector of IDE objects.
          auto const& energyDeposits = timeSlice.second;

          // An "energy deposit" object stores how much charge/energy was deposited in a small volume, by which particle, and where.
          // The type of 'energyDeposit' will be sim::IDE, here use auto.
          /*
          for ( auto const& energyDeposit : energyDeposits )
          {

            auto search = particleMap.find( energyDeposit.trackID );
            if ( search == particleMap.end() ) continue;

            // "search" points to a pair in the map: <track ID, MCParticle*>
            const simb::MCParticle& particle = *((*search).second);

            // Deposit energy for lepton
            if ( particle.Process() == "primary" && ( abs(particle.PdgCode()) == 13 || abs(particle.PdgCode()) == 11 || abs(particle.PdgCode()) == 15 ))
             // if ( abs(particle.PdgCode()) == 13)
            {
              // Method a
              if ( fGeometryService->SignalType(channelNumber) == geo::kCollection )
              {
                fSim_mu_Edep_a1 += energyDeposit.numElectrons * fElectronsToGeV;
                fSim_mu_Edep_a2 += energyDeposit.energy;
              }
              else{
                fSim_mu_Edep_NonCollectionPlane_a2 += energyDeposit.energy;
              }
              // Method b
              std::vector<geo::WireID> const Wires = fGeometryService->ChannelToWire(channelNumber);
              if ( Wires[0].planeID().Plane == 0 )
              {
                fSim_mu_Edep_b1 += energyDeposit.numElectrons * fElectronsToGeV;
                fSim_mu_Edep_b2 += energyDeposit.energy;
              } // end if access plane info via channel -> wire -> plane ID is 0
              else{
                fSim_mu_Edep_NonCollectionPlane_b2 += energyDeposit.energy;
              }

            }// end muon energy deposit
            else
            {
              // If it's not from primary leptons, count it as hadronic
              // Method a: only include the energy from the collection plane: geo::kCollection defined in
              // ${LARCOREOBJ_INC}/larcoreobj/SimpleTypesAndConstants/geo_types.h
              //
              if ( fGeometryService->SignalType(channelNumber) == geo::kCollection )
              {
                //
                // Note: there are also two ways to get E deposit for sim::IDE
                //
                fSim_hadronic_Edep_a1 += energyDeposit.numElectrons * fElectronsToGeV;
                fSim_hadronic_Edep_a2 += energyDeposit.energy;

                // Store position and E for each deposit
                fSim_hadronic_hit_x_a.push_back(energyDeposit.x);
                fSim_hadronic_hit_y_a.push_back(energyDeposit.y);
                fSim_hadronic_hit_z_a.push_back(energyDeposit.z);
                fSim_hadronic_hit_Edep_a1.push_back(energyDeposit.numElectrons * fElectronsToGeV);
                fSim_hadronic_hit_Edep_a2.push_back(energyDeposit.energy);
              } // end if access plane info via channel signal type
              else{
                fSim_hadronic_Edep_NonCollectionPlane_a2 += energyDeposit.energy;
              }

              //
              // Method b: navigate via channel -> wire -> plane ID, and require planeID to be 0.
              //
              std::vector<geo::WireID> const Wires = fGeometryService->ChannelToWire(channelNumber);
              if ( Wires[0].planeID().Plane == 0 )
              {

                fSim_hadronic_Edep_b1 += energyDeposit.numElectrons * fElectronsToGeV;
                fSim_hadronic_Edep_b2 += energyDeposit.energy;

                // Store position and E for each deposit
                fSim_hadronic_hit_x_b.push_back(energyDeposit.x);
                fSim_hadronic_hit_y_b.push_back(energyDeposit.y);
                fSim_hadronic_hit_z_b.push_back(energyDeposit.z);
                fSim_hadronic_hit_Edep_b1.push_back(energyDeposit.numElectrons * fElectronsToGeV);
                fSim_hadronic_hit_Edep_b2.push_back(energyDeposit.energy);
              } // end if access plane info via channel -> wire -> plane ID is 0
              else{
                fSim_hadronic_Edep_NonCollectionPlane_b2 += energyDeposit.energy;
              }

            } // end if hadronic
          }   // end For each energy deposit
          */
          for ( auto const& energyDeposit : energyDeposits ) {
            /*
            // Method a
            if ( fGeometryService->SignalType(channelNumber) == geo::kCollection )
            {
              // Still do the search, but now only for primary lepton
              auto search = particleMap.find( energyDeposit.trackID );

              if ( search != particleMap.end() ) { // found match in map

                const simb::MCParticle& particle = *((*search).second);

                // if it's primary lepton
                if ( particle.Process() == "primary" && ( abs(particle.PdgCode()) == 13 || abs(particle.PdgCode()) == 11 || abs(particle.PdgCode()) == 15 ))
                {
                  fSim_mu_Edep_a2_debug += energyDeposit.energy;
                  // now continue to the next energy deposit
                  continue;
                } // end if it's primary lepton

                // if it's not a primary lepton, do nothing
              }// end found match

              // If the energyDeposit made this far, count it as hadronic deposits (primary+secondary), do not involve particleMap
              fSim_hadronic_Edep_a2_debug += energyDeposit.energy;
            } // end method a
            */
            // Method b: First check if it's on collection plane
            std::vector<geo::WireID> const Wires = fGeometryService->ChannelToWire(channelNumber);
            if ( Wires[0].planeID().Plane == 0 ) {

              // All EM shower are treated as secondary interactions, and their particles are not saved in the MC particle list
              // Still do the search, but now only for primary lepton (particleMap trkID is always positive)
              // Also search for EM shower particles from primary lepton, these deposits has trkID that's negative of the primary lepton trkID
              auto search = particleMap.find( abs(energyDeposit.trackID) );


              if ( search != particleMap.end() ) { // found match in map

                const simb::MCParticle& particle = *((*search).second);
                // if the energy deposit is from primary lepton,
                // or its ancestor mother particle is the primary lepton (e.g., from muon decays)
                if ( ( particle.Process() == "primary" && ( abs(particle.PdgCode()) == 13 || abs(particle.PdgCode()) == 11 || abs(particle.PdgCode()) == 15 ) ) || IsAncestorMotherPrimaryLep(particle, primarylep_trkID, particleMap) )
                {
                  fSim_mu_Edep_b2 += energyDeposit.energy;
                  // now continue to the next energy deposit
                  continue;
                } // end lepton dep e
                // if it's not, do nothing
              } // end found match

              // If the energyDeposit made this far, it's counted as hadronic deposits (primary+secondary), do not involve particleMap
              fSim_hadronic_Edep_b2 += energyDeposit.energy;
              fSim_hadronic_hit_x_b.push_back(energyDeposit.x);
              fSim_hadronic_hit_y_b.push_back(energyDeposit.y);
              fSim_hadronic_hit_z_b.push_back(energyDeposit.z);
              fSim_hadronic_hit_Edep_b2.push_back(energyDeposit.energy);

              // Store trackID for debug
              EDepTrackID = energyDeposit.trackID;
              auto exist = EDepMap.find( EDepTrackID );
              // if can't find, store it and add the edep
              if ( exist == EDepMap.end() ) {
                EDep_TrackID_vec.push_back( EDepTrackID ); // negative trackID exists
                EDepMap[EDepTrackID] = energyDeposit.energy;
              } else { // find the same track id
                EDepMap[EDepTrackID] += energyDeposit.energy;
              }

            } // end plane == 0
          } // end energy deposit loop
        } // end For each time slice
      } // end For each SimChannel

      fSim_n_hadronic_Edep_b = fSim_hadronic_hit_x_b.size();

      // Print out EDepMap
      //if ( false ) std::cout << "fGen_numu_E: "<< fGen_numu_E << ", fSim_mu_Edep_b2_debug: " << fSim_mu_Edep_b2_debug<< ", fSim_hadronic_Edep_b2_debug: " << fSim_hadronic_Edep_b2_debug << ", Tot had E track id: " << EDep_TrackID_vec.size() << std::endl;
      if ( false ) {
        for (long unsigned int i=0; i< EDep_TrackID_vec.size(); i++) {
          std::cout << "Evt track id: " << EDep_TrackID_vec.at(i) << std::endl;
        }
        std::map<int, double>::iterator it;
        std::cout << "TrackID" << " | " << "Tot EDep" << std::endl;
        for (it = EDepMap.begin(); it != EDepMap.end(); it++) std::cout << "    " << it->first << " | " << it->second << std::endl;
      }

      // In general, objects in the LArSoft reconstruction chain are linked using the art::Assns class:
      // <https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/Using_art_in_LArSoft#artAssns>
      // The following statement will find the simb::MCTruth associated with the simb::MCParticle
      const art::FindManyP<simb::MCTruth> findManyTruth(particleHandle, event, fSimulationProducerLabel);

      if ( !findManyTruth.isValid() ) {
        std::cout << "findManyTruth simb::MCTruth for simb::MCParticle failed!" << std::endl;
      }

      size_t particle_index = 0; // only look at first particle in particleHandle's vector.
      auto const& truth = findManyTruth.at(particle_index);

      // Make sure there's no problem.
      if ( truth.empty() ) {
        std::cout << "Particle ID=" << particleHandle->at(particle_index).TrackId() << " has no primary!" << std::endl;
      }

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

  // If this returns true, then the energy deposit is associated with primary lepton
  bool IsAncestorMotherPrimaryLep(const simb::MCParticle& p1, int primarylep_trkID, std::map<int, const simb::MCParticle*> particleMap) {
    int MothertrkID = p1.Mother();
    // Immediate mother is the primary lep
    if ( MothertrkID == primarylep_trkID )  return true;
    // Immediate mother is not primary lep, but other primary particles from genie
    else if ( MothertrkID == 0 ) return false;
    // Keep looking upstream, find it in particleMap
    else {
      auto tmp_search = particleMap.find( MothertrkID ); // this search must be found, can't be null
      const simb::MCParticle& tmp_mother = *((*tmp_search).second);
      return IsAncestorMotherPrimaryLep(tmp_mother, primarylep_trkID, particleMap);
    }
  } // end GetAncestorMotherTrkID

} // local namespace
