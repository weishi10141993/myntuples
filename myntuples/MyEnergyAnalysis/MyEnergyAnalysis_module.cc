/**
 * @file   MyEnergyAnalysis_module.cc
 * @brief  A basic "skeleton" to read in art::Event records from a file,
 *         access their information, and do something with them.
 * @ingroup MyEnergyAnalysis
 * @author William Seligman (seligman@nevis.columbia.edu)
 *
 * See
 * <https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/Using_art_in_LArSoft>
 * for a description of the ART classes used here.
 *
 * Almost everything you see in the code below may have to be changed
 * by you to suit your task. The example task is to make histograms
 * and n-tuples related to @f$ dE/dx @f$ of particle tracks in the detector.
 *
 * As you try to understand why things are done a certain way in this
 * example ("What's all this stuff about 'auto const&'?"), it will help
 * to read ADDITIONAL_NOTES.md in the same directory as this file.
 *
 * Also note that, despite our efforts, the documentation and the practices in
 * this code may fall out of date. In doubt, ask!
 *
 * The last revision of this code was done in August 2017 with LArSoft
 * v06_44_00.
 *
 * This is in-source documentation in Doxygen format. Doxygen is a
 * utility that creates web pages from specially-formatted comments in
 * program code. If your package is ever added to an official LArSoft
 * release, the doxygen-style comments will be added to the LArSoft
 * doxygen pages at <http://nusoft.fnal.gov/larsoft/doxsvn/html/>.
 *
 * You can see the doxygenated version of the following lines at
 * <http://nusoft.fnal.gov/larsoft/doxsvn/html/classlar_1_1example_1_1MyEnergyAnalysis.html>
 *
 * Doxygen tips:
 *
 * In general, comments that start with "///" or in C++ comment blocks
 * (like this one) will appear in the web pages created by Doxygen.
 *
 * A comment starting with "///" will be associated wth the next code
 * line that follows it, while one starting with "///<" will be
 * associated with the one above it.
 *
 * For more on doxygen, see the documentation at
 * <http://www.stack.nl/~dimitri/doxygen/manual/index.html>
 *
 */

// Always include headers defining everything you use. Starting from
// LArSoft and going up the software layers (nusimdata, art, etc.)
// ending with C++ is standard.

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

// ROOT includes. Note: To look up the properties of the ROOT classes,
// use the ROOT web site; e.g.,
// <https://root.cern.ch/doc/master/annotated.html>
#include "TH1.h"
#include "TLorentzVector.h"
#include "TTree.h"
#include "TVector3.h"

// C++ includes
#include <cmath>
#include <map>

namespace {

  // This is a local namespace (as opposed to the one below, which has
  // the nested name lar::example::).
  //
  // The stuff you declare in an local namespace will not be
  // visible beyond this file (more technically, this "translation
  // unit", which is everything that gets included in the .cc file from
  // now on). In this way, any functions you define in this namespace
  // won't affect the environment of other modules.

  // We will define functions at the end, but we declare them here so
  // that the module can freely use them.

  /// Utility function to get the diagonal of the detector
  /// @ingroup MyEnergyAnalysis
  double DetectorDiagonal(geo::GeometryCore const& geom);

  /// Comparison routine for using std::lower/upper_bound to search
  /// TDCIDE vectors.
  /// @ingroup MyEnergyAnalysis
  bool TDCIDETimeCompare(const sim::TDCIDE&, const sim::TDCIDE&);

} // local namespace

// It is good programming practice to put your code into a namespace,
// to make sure that no method or function you define will conflict
// with anything in any other module with the same name. Here we
// follow the LArSoft standard of declaring a main namespace of "lar"
// with a nested namespace of "example" because we're in the
// larexamples product. If an outside package wanted to call this
// module, it would have to use something like
// lar::example::MyEnergyAnalysis.

namespace lar {
  namespace example {

    // BEGIN MyEnergyAnalysis group
    // -----------------------------------------------
    /// @ingroup MyEnergyAnalysis
    /// @{
    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    // class definition
    /**
     * @brief Example analyzer module.
     * @see @ref MyEnergyAnalysis "analysis module example overview"
     *
     * This class extracts information from the generated and reconstructed
     * particles.
     *
     * It produces histograms for the simulated particles in the input file:
     * - PDG ID (flavor) of all particles
     * - momentum of the primary particles selected to have a specific PDG ID
     * - length of the selected particle trajectory
     *
     * It also produces two ROOT trees.
     *
     * The first ROOT tree contains information on the simulated
     * particles, including "dEdx", a binned histogram of collected
     * charge as function of track range.
     *
     * The second ROOT tree contains information on the reconstructed hits.
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
     *
     * - *PDGcode* (integer, mandatory): particle type (PDG ID) to be selected;
     *   only primary particles of this type will be considered
     *
     * - *BinSize* (real, mandatory): @f$ dx @f$ [cm] used for the @f$ dE/dx @f$
     *   calculation
     *
     */
    class MyEnergyAnalysis : public art::EDAnalyzer {
    public:
      // This structure describes the configuration parameters of the
      // module.  Any missing or unknown parameters will generate a
      // configuration error.
      //
      // With an additional trick (see below) it allows configuration
      // documentation to be displayed as a command-line option (see
      // below).
      //
      // Note that, in this example, the Name() string (that is, the
      // name you call a parameter in the FHiCL configuration file) and
      // the data member name in the Config struct (that is, the name
      // you access that parameter in your C++ code) match. This is not
      // required, but it makes it easier to remember them.
      //
      // More details at:
      // https://cdcvs.fnal.gov/redmine/projects/fhicl-cpp/wiki/Configuration_validation_and_fhiclcpp_types
      //
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

        fhicl::Atom<int> PDGcode{
          Name("PDGcode"),
          Comment("particle type (PDG ID) of the primary particle to be selected")};

        fhicl::Atom<double> BinSize{Name("BinSize"),
                                    Comment("dx [cm] used for the dE/dx calculation")};

      }; // Config

      // If we define "Parameters" in this way, art will know to use the
      // "Comment" items above for the module description. See what
      // this command does:
      //
      // lar --print-description MyEnergyAnalysis
      //
      // The details of the art side of this trick are at:
      //
      // https://cdcvs.fnal.gov/redmine/projects/art/wiki/Configuration_validation_and_description
      //
      using Parameters = art::EDAnalyzer::Table<Config>;

      // -------------------------------------------------------------------
      // -------------------------------------------------------------------
      // Standard constructor for an ART module with configuration validation;
      // we don't need a special destructor here.

      /// Constructor: configures the module (see the Config structure above)
      explicit MyEnergyAnalysis(Parameters const& config);

      // The following methods have a standard meaning and a standard signature
      // inherited from the framework (art::EDAnalyzer class).

      // - The "virtual" keyword is a reminder that the function we are
      //   dealing with is, in fact, virtual. You don't need to
      //   understand it now, but it's very important when you write a
      //   new algorithm.

      // - The "override" keyword, new in C++ 2011, is an important
      //   safety measure that ensures that the method we are going to
      //   write will override a matching method in the base class. For
      //   example, if we mistyped it as

      //   virtual void beginJob() const;

      //   (adding "const"), the compiler will be very happy with it,
      //   but art will not touch it, because art needs a "void
      //   beginJob()" (non-const) and it will find one in the base
      //   class (void art::EDAnalyzer::beginJob()) and will silently
      //   use that one instead. If you accidentally use:

      //   virtual void beginJob() const override;

      //   the compiler will immediately complain with us that this
      //   method is overriding nothing, hinting to some mistake (the
      //   spurious "const" in this case).

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
      // The stuff below is the part you'll most likely have to change to
      // go from this custom example to your own task.

      // The parameters we'll read from the .fcl file.
      art::InputTag fSimulationProducerLabel; ///< The name of the producer that tracked
                                              ///< simulated particles through the detector
      art::InputTag fHitProducerLabel;        ///< The name of the producer that created hits
      art::InputTag fClusterProducerLabel;    ///< The name of the producer that
                                              ///< created clusters
      int fSelectedPDG;                       ///< PDG code of particle we'll focus on
      double fBinSize;                        ///< For dE/dx work: the value of dx.

      // Pointers to the histograms we'll create.
      TH1D* fPDGCodeHist;     ///< PDG code of all particles
      TH1D* fMomentumHist;    ///< momentum [GeV] of all selected particles
      TH1D* fTrackLengthHist; ///< true length [cm] of all selected particles

      // The n-tuples we'll create.
      TTree* fSimulationNtuple;     ///< tuple with simulated data
      TTree* fReconstructionNtuple; ///< tuple with reconstructed data

      // The comment lines with the @ symbols define groups in doxygen.
      /// @name The variables that will go into both n-tuples.
      /// @{
      int fEvent;  ///< number of the event being processed
      int fRun;    ///< number of the run being processed
      int fSubRun; ///< number of the sub-run being processed
      /// @}

      /// @name The variables that will go into the simulation n-tuple.
      /// @{
      int fSimPDG;     ///< PDG ID of the particle being processed
      int fSimTrackID; ///< GEANT ID of the particle being processed

      // Arrays for 4-vectors: (x,y,z,t) and (Px,Py,Pz,E).
      // Note: old-style C++ arrays are considered obsolete. However,
      // to create simple n-tuples, we still need to use them.
      double fStartXYZT[4]; ///< (x,y,z,t) of the true start of the particle
      double fEndXYZT[4];   ///< (x,y,z,t) of the true end of the particle
      double fStartPE[4];   ///< (Px,Py,Pz,E) at the true start of the particle
      double fEndPE[4];     ///< (Px,Py,Pz,E) at the true end of the particle

      /// Number of dE/dx bins in a given track.
      int fSimNdEdxBins;

      /// The vector that will be used to accumulate dE/dx values as a function
      /// of range.
      std::vector<double> fSimdEdxBins;
      /// @}

      /// @name Variables used in the reconstruction n-tuple
      /// @{
      int fRecoPDG;     ///< PDG ID of the particle being processed
      int fRecoTrackID; ///< GEANT ID of the particle being processed

      /// Number of dE/dx bins in a given track.
      int fRecoNdEdxBins;

      /// The vector that will be used to accumulate dE/dx values as a function
      /// of range.
      std::vector<double> fRecodEdxBins;

      /// @}

      // Other variables that will be shared between different methods.
      geo::GeometryCore const* fGeometryService; ///< pointer to Geometry provider
      double fElectronsToGeV;                    ///< conversion factor
      int fTriggerOffset; ///< (units of ticks) time of expected neutrino event

    }; // class MyEnergyAnalysis

    /// @}
    // END MyEnergyAnalysis group
    // -------------------------------------------------

    //-----------------------------------------------------------------------
    //-----------------------------------------------------------------------
    // class implementation

    //-----------------------------------------------------------------------
    // Constructor
    //
    // Note that config is a Table<Config>, and to access the Config
    // value we need to use an operator: "config()". In the same way,
    // each element in Config is an Atom<Type>, so to access the type we
    // again use the call operator, e.g. "SimulationLabel()".
    //
    MyEnergyAnalysis::MyEnergyAnalysis(Parameters const& config)
      : EDAnalyzer(config)
      , fSimulationProducerLabel(config().SimulationLabel())
      , fHitProducerLabel(config().HitLabel())
      , fClusterProducerLabel(config().ClusterLabel())
      , fSelectedPDG(config().PDGcode())
      , fBinSize(config().BinSize())
    {
      // Get a pointer to the geometry service provider.
      fGeometryService = lar::providerFrom<geo::Geometry>();

      auto const clock_data =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataForJob();
      fTriggerOffset = trigger_offset(clock_data);

      // Since art 2.8, you can and should tell beforehand, here in the
      // constructor, all the data the module is going to read ("consumes") or
      // might read
      // ("may_consume"). Diligence here will in the future help the framework
      // execute modules in parallel, making sure the order is correct.
      consumes<std::vector<simb::MCParticle>>(fSimulationProducerLabel);
      consumes<std::vector<sim::SimChannel>>(fSimulationProducerLabel);
      consumes<art::Assns<simb::MCTruth, simb::MCParticle>>(fSimulationProducerLabel);
      consumes<std::vector<recob::Hit>>(fHitProducerLabel);
      consumes<std::vector<recob::Cluster>>(fClusterProducerLabel);
      consumes<art::Assns<recob::Cluster, recob::Hit>>(fHitProducerLabel);
    }

    //-----------------------------------------------------------------------
    void
    MyEnergyAnalysis::beginJob()
    {
      // Get the detector length, to determine the maximum bin edge of one
      // of the histograms.
      const double detectorLength = DetectorDiagonal(*fGeometryService);

      // Access ART's TFileService, which will handle creating and writing
      // histograms and n-tuples for us.
      art::ServiceHandle<art::TFileService const> tfs;

      // For TFileService, the arguments to 'make<whatever>' are the
      // same as those passed to the 'whatever' constructor, provided
      // 'whatever' is a ROOT class that TFileService recognizes.

      // Define the histograms. Putting semi-colons around the title
      // causes it to be displayed as the x-axis label if the histogram
      // is drawn (the format is "title;label on abscissae;label on ordinates").
      fPDGCodeHist = tfs->make<TH1D>("pdgcodes", ";PDG Code;", 5000, -2500, 2500);
      fMomentumHist = tfs->make<TH1D>("mom", ";particle Momentum (GeV);", 100, 0., 10.);
      fTrackLengthHist =
        tfs->make<TH1D>("length", ";particle track length (cm);", 200, 0, detectorLength);

      // Define our n-tuples, which are limited forms of ROOT
      // TTrees. Start with the TTree itself.
      fSimulationNtuple =
        tfs->make<TTree>("MyEnergyAnalysisSimulation", "MyEnergyAnalysisSimulation");
      fReconstructionNtuple =
        tfs->make<TTree>("MyEnergyAnalysisReconstruction", "MyEnergyAnalysisReconstruction");

      // Define the branches (columns) of our simulation n-tuple. To
      // write a variable, we give the address of the variable to
      // TTree::Branch.
      fSimulationNtuple->Branch("Event", &fEvent, "Event/I");
      fSimulationNtuple->Branch("SubRun", &fSubRun, "SubRun/I");
      fSimulationNtuple->Branch("Run", &fRun, "Run/I");
      fSimulationNtuple->Branch("TrackID", &fSimTrackID, "TrackID/I");
      fSimulationNtuple->Branch("PDG", &fSimPDG, "PDG/I");
      // When we write arrays, we give the address of the array to
      // TTree::Branch; in C++ this is simply the array name.
      fSimulationNtuple->Branch("StartXYZT", fStartXYZT, "StartXYZT[4]/D");
      fSimulationNtuple->Branch("EndXYZT", fEndXYZT, "EndXYZT[4]/D");
      fSimulationNtuple->Branch("StartPE", fStartPE, "StartPE[4]/D");
      fSimulationNtuple->Branch("EndPE", fEndPE, "EndPE[4]/D");
      // For a variable-length array: include the number of bins.
      fSimulationNtuple->Branch("NdEdx", &fSimNdEdxBins, "NdEdx/I");
      // ROOT branches can contain std::vector objects.
      fSimulationNtuple->Branch("dEdx", &fSimdEdxBins);

      // A similar definition for the reconstruction n-tuple. Note that we
      // use some of the same variables in both n-tuples.
      fReconstructionNtuple->Branch("Event", &fEvent, "Event/I");
      fReconstructionNtuple->Branch("SubRun", &fSubRun, "SubRun/I");
      fReconstructionNtuple->Branch("Run", &fRun, "Run/I");
      fReconstructionNtuple->Branch("TrackID", &fRecoTrackID, "TrackID/I");
      fReconstructionNtuple->Branch("PDG", &fRecoPDG, "PDG/I");
      fReconstructionNtuple->Branch("NdEdx", &fRecoNdEdxBins, "NdEdx/I");
      fReconstructionNtuple->Branch("dEdx", &fRecodEdxBins);
    }

    //-----------------------------------------------------------------------
    // art expects this function to have a art::Run argument; C++
    // expects us to use all the arguments we are given, or it will
    // generate an "unused variable" warning. But we don't actually need
    // nor use the art::Run object in this example. The trick to prevent
    // that warning is to omit (or comment out) the name of the
    // parameter.

    void
    MyEnergyAnalysis::beginRun(const art::Run& /*run*/)
    {
      // How to convert from number of electrons to GeV. The ultimate
      // source of this conversion factor is
      // ${LARCOREOBJ_INC}/larcoreobj/SimpleTypesAndConstants/PhysicalConstants.h.
      // But sim::LArG4Parameters might in principle ask a database for it.
      art::ServiceHandle<sim::LArG4Parameters const> larParameters;
      fElectronsToGeV = 1. / larParameters->GeVToElectrons();
    }

    //-----------------------------------------------------------------------
    void
    MyEnergyAnalysis::analyze(const art::Event& event)
    {
      // Start by fetching some basic event information for our n-tuple.
      fEvent = event.id().event();
      fRun = event.run();
      fSubRun = event.subRun();

      // This is the standard method of reading multiple objects
      // associated with the same event; see
      // <https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/Using_art_in_LArSoft>
      // for more information.
      //
      // Define a "handle" to point to a vector of the objects.
      art::Handle<std::vector<simb::MCParticle>> particleHandle;

      // Then tell the event to fill the vector with all the objects of
      // that type produced by a particular producer.
      //
      // Note that if I don't test whether this is successful, and there
      // aren't any simb::MCParticle objects, then the first time we
      // access particleHandle, art will display a "ProductNotFound"
      // exception message and, depending on art settings, it may skip
      // all processing for the rest of this event (including any
      // subsequent analysis steps) or stop the execution.
      if (!event.getByLabel(fSimulationProducerLabel, particleHandle)) {
        // If we have no MCParticles at all in an event, then we're in
        // big trouble. Throw an exception to force this module to
        // stop. Try to include enough information for the user to
        // figure out what's going on. Note that when we throw a
        // cet::exception, the run and event number will automatically
        // be displayed.
        //
        // __LINE__ and __FILE__ are values computed by the compiler.

        throw cet::exception("MyEnergyAnalysis")
          << " No simb::MCParticle objects in this event - "
          << " Line " << __LINE__ << " in file " << __FILE__ << std::endl;
      }

      // Get all the simulated channels for the event. These channels
      // include the energy deposited for each simulated track.
      //
      // Here we use a different method to access objects:
      // art::ValidHandle. Using this method, if there aren't any
      // objects of the given type (sim::SimChannel in this case) in the
      // input file for this art::Event, art will throw a
      // ProductNotFound exception.
      //
      // The "auto" type means that the C++ compiler will determine the
      // appropriate type for us, based on the return type of
      // art::Event::getValidHandle<T>(). The "auto" keyword is a great
      // timesaver, especially with frameworks like LArSoft which often
      // have complicated data product structures.

      auto simChannelHandle =
        event.getValidHandle<std::vector<sim::SimChannel>>(fSimulationProducerLabel);

      //
      // Let's compute the variables for the simulation n-tuple first.
      //

      // The MCParticle objects are not necessarily in any particular
      // order. Since we may have to search the list of particles, let's
      // put them into a map, a sorted container that will make
      // searching fast and easy. To save both space and time, the map
      // will not contain a copy of the MCParticle, but a pointer to it.
      std::map<int, const simb::MCParticle*> particleMap;

      // This is a "range-based for loop" in the 2011 version of C++; do
      // a web search on "c++ range based for loop" for more
      // information. Here's how it breaks down:

      // - A range-based for loop operates on a container.
      //   "particleHandle" is not a container; it's a pointer to a
      //   container. If we want C++ to "see" a container, we have to
      //   dereference the pointer, like this: *particleHandle.

      // - The loop variable that is set to each object in the container
      //   is named "particle". As for the loop variable's type:

      //   - To save a little bit of typing, and to make the code easier
      //     to maintain, we're going to let the C++ compiler deduce the
      //     type of what's in the container (simb::MCParticle objects
      //     in this case), so we use "auto".

      //   - We do _not_ want to change the contents of the container,
      //     so we use the "const" keyword to make sure.

      //   - We don't need to copy each object from the container into
      //     the variable "particle". It's sufficient to refer to the
      //     object by its address. So we use the reference operator "&"
      //     to tell the compiler to just copy the address, not the
      //     entire object.

      // It sounds complicated, but the end result is that we loop over
      // the list of particles in the art::Event in the most efficient
      // way possible.

      for (auto const& particle : (*particleHandle)) {
        // For the methods you can call to get particle information,
        // see ${NUSIMDATA_INC}/nusimdata/SimulationBase/MCParticle.h.
        fSimTrackID = particle.TrackId();

        // Add the address of the MCParticle to the map, with the
        // track ID as the key.
        particleMap[fSimTrackID] = &particle;

        // Histogram the PDG code of every particle in the event.
        fSimPDG = particle.PdgCode();
        fPDGCodeHist->Fill(fSimPDG);

        // For this example, we want to fill the n-tuples and histograms
        // only with information from the primary particles in the
        // event, whose PDG codes match a value supplied in the .fcl file.
        if (particle.Process() != "primary" || fSimPDG != fSelectedPDG) continue;

        // A particle has a trajectory, consisting of a set of
        // 4-positions and 4-mommenta.
        const size_t numberTrajectoryPoints = particle.NumberTrajectoryPoints();

        // For trajectories, as for vectors and arrays, the first
        // point is #0, not #1.
        const int last = numberTrajectoryPoints - 1;
        const TLorentzVector& positionStart = particle.Position(0);
        const TLorentzVector& positionEnd = particle.Position(last);
        const TLorentzVector& momentumStart = particle.Momentum(0);
        const TLorentzVector& momentumEnd = particle.Momentum(last);

        // Make a histogram of the starting momentum.
        fMomentumHist->Fill(momentumStart.P());

        // Fill arrays with the 4-values. (Don't be fooled by
        // the name of the method; it just puts the numbers from
        // the 4-vector into the array.)
        positionStart.GetXYZT(fStartXYZT);
        positionEnd.GetXYZT(fEndXYZT);
        momentumStart.GetXYZT(fStartPE);
        momentumEnd.GetXYZT(fEndPE);

        // Use a polar-coordinate view of the 4-vectors to
        // get the track length.
        const double trackLength = (positionEnd - positionStart).Rho();

        // Let's print some debug information in the job output to see
        // that everything is fine. MF_LOG_DEBUG() is a messagefacility
        // macro that prints stuff when the message level is set to
        // standard_debug in the .fcl file.
        MF_LOG_DEBUG("MyEnergyAnalysis") << "Track length: " << trackLength << " cm";

        // Fill a histogram of the track length.
        fTrackLengthHist->Fill(trackLength);

        MF_LOG_DEBUG("MyEnergyAnalysis")
          << "track ID=" << fSimTrackID << " (PDG ID: " << fSimPDG << ") " << trackLength
          << " cm long, momentum " << momentumStart.P() << " GeV/c, has " << numberTrajectoryPoints
          << " trajectory points";

        // Determine the number of dE/dx bins for the n-tuple.
        fSimNdEdxBins = int(trackLength / fBinSize) + 1;
        // Initialize the vector of dE/dx bins to be empty.
        fSimdEdxBins.clear();

        // To look at the energy deposited by this particle's track,
        // we loop over the SimChannel objects in the event.
        for (auto const& channel : (*simChannelHandle)) {
          // Get the numeric ID associated with this channel. (The
          // channel number is a 32-bit unsigned int, which normally
          // requires a special data type. Let's use "auto" so we
          // don't have to remember "raw::ChannelID_t".
          auto const channelNumber = channel.Channel();

          // A little care: There is more than one plane that reacts
          // to charge in the TPC. We only want to include the
          // energy from the collection plane. Note:
          // geo::kCollection is defined in
          // ${LARCOREOBJ_INC}/larcoreobj/SimpleTypesAndConstants/geo_types.h
          if (fGeometryService->SignalType(channelNumber) != geo::kCollection) continue;

          // Each channel has a map inside it that connects a time
          // slice to energy deposits in the detector. We'll use
          // "auto", but it's worth noting that the full type of
          // this map is
          // std::map<unsigned short, std::vector<sim::IDE>>
          auto const& timeSlices = channel.TDCIDEMap();

          // For every time slice in this channel:
          for (auto const& timeSlice : timeSlices) {
            // Each entry in a map is a pair<first,second>. For
            // the timeSlices map, the 'first' is a time slice
            // number. The 'second' is a vector of IDE objects.
            auto const& energyDeposits = timeSlice.second;

            // Loop over the energy deposits. An "energy deposit"
            // object is something that knows how much
            // charge/energy was deposited in a small volume, by
            // which particle, and where. The type of
            // 'energyDeposit' will be sim::IDE, which is defined
            // in
            // ${LARDATAOBJ_INC}/lardataobj/Simulation/SimChannel.h.
            for (auto const& energyDeposit : energyDeposits) {
              // Check if the track that deposited the
              // energy matches the track of the particle.
              if (energyDeposit.trackID != fSimTrackID) continue;
              // Get the (x,y,z) of the energy deposit.
              TVector3 location(energyDeposit.x, energyDeposit.y, energyDeposit.z);

              // Distance from the start of the track.
              const double distance = (location - positionStart.Vect()).Mag();

              // Into which bin of the dE/dx array do we add the energy?
              const unsigned int bin = (unsigned int)(distance / fBinSize);

              // Is the dE/dx array big enough to include this bin?
              if (fSimdEdxBins.size() < bin + 1) {
                //  Increase the array size to accomodate
                //  the new bin, padding it with zeros.
                fSimdEdxBins.resize(bin + 1, 0.);
              }

              // Add the energy deposit to that bin. (If you look at the
              // definition of sim::IDE, you'll see that there's another
              // way to get the energy. Are the two methods equivalent?
              // Compare the results and see!)
              fSimdEdxBins[bin] += energyDeposit.numElectrons * fElectronsToGeV;

            } // For each energy deposit
          }   // For each time slice
        }     // For each SimChannel

        // At this point we've filled in the values of all the
        // variables and arrays we want to write to the n-tuple
        // for this particle. The following command actually
        // writes those values.
        fSimulationNtuple->Fill();

      } // loop over all particles in the event.

      //
      // Reconstruction n-tuple
      //

      // All of the above is based on objects entirely produced by the
      // simulation. Let's try to do something based on reconstructed
      // objects. A Hit (see ${LARDATAOBJ_INC}/lardataobj/RecoBase/Hit.h)
      // is a 2D object in a plane.

      // This code duplicates much of the code in the previous big
      // simulation loop, and will produce the similar results. (It
      // won't be identical, because of shaping and other factors; not
      // every bit of charge in a channel ends up contributing to a
      // hit.) The point is to show different methods of accessing
      // information, not to produce profound physics results -- that
      // part is up to you!

      // For the rest of this method, I'm going to assume you've read
      // the comments in previous section; I won't repeat all the C++
      // coding tricks and whatnot.

      // Start by reading the Hits. We don't use art::ValidHandle here
      // because there might be no hits in the input; e.g., we ran the
      // simulation but did not run the reconstruction steps yet. If
      // there are no hits we may as well skip the rest of this module.

      art::Handle<std::vector<recob::Hit>> hitHandle;
      if (!event.getByLabel(fHitProducerLabel, hitHandle)) return;

      // Our goal is to accumulate the dE/dx of any particles associated
      // with the hits that match our criteria: primary particles with
      // the PDG code from the .fcl file. I don't know how many such
      // particles there will be in a given event. We'll use a map, with
      // track ID as the key, to hold the vectors of dE/dx information.
      std::map<int, std::vector<double>> dEdxMap;

      auto const clock_data =
        art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(event);

      // For every Hit:
      for (auto const& hit : (*hitHandle)) {
        // The channel associated with this hit.
        auto hitChannelNumber = hit.Channel();

        // We have a hit. For this example let's just focus on the
        // hits in the collection plane.
        if (fGeometryService->SignalType(hitChannelNumber) != geo::kCollection) continue;

        MF_LOG_DEBUG("MyEnergyAnalysis") << "Hit in collection plane" << std::endl;

        // In a few lines we're going to look for possible energy
        // deposits that correspond to that hit. Determine a
        // reasonable range of times that might correspond to those
        // energy deposits.

        // In reconstruction, the channel waveforms are truncated. So
        // we have to adjust the Hit TDC ticks to match those of the
        // SimChannels, which were created during simulation.

        // Save a bit of typing, while still allowing for potential
        // changes in the definitions of types in
        // $LARDATAOBJ_DIR/source/lardataobj/Simulation/SimChannel.h

        typedef sim::SimChannel::StoredTDC_t TDC_t;
        TDC_t start_tdc = clock_data.TPCTick2TDC(hit.StartTick());
        TDC_t end_tdc = clock_data.TPCTick2TDC(hit.EndTick());
        TDC_t hitStart_tdc = clock_data.TPCTick2TDC(hit.PeakTime() - 3. * hit.SigmaPeakTime());
        TDC_t hitEnd_tdc = clock_data.TPCTick2TDC(hit.PeakTime() + 3. * hit.SigmaPeakTime());

        start_tdc = std::max(start_tdc, hitStart_tdc);
        end_tdc = std::min(end_tdc, hitEnd_tdc);

        // In the simulation section, we started with particles to find
        // channels with a matching track ID. Now we search in reverse:
        // search the SimChannels for matching channel number, then look
        // at the tracks inside the channel.

        for (auto const& channel : (*simChannelHandle)) {
          auto simChannelNumber = channel.Channel();
          if (simChannelNumber != hitChannelNumber) continue;

          MF_LOG_DEBUG("MyEnergyAnalysis")
            << "SimChannel number = " << simChannelNumber << std::endl;

          // The time slices in this channel.
          auto const& timeSlices = channel.TDCIDEMap();

          // We want to look over the range of time slices in this
          // channel that correspond to the range of hit times.

          // To do this, we're going to use some fast STL search
          // methods; STL algorithms are usually faster than the
          // ones you might write yourself. The price for this speed
          // is a bit of code complexity: in particular, we need to
          // a custom comparison function, TDCIDETimeCompare, to
          // define a "less-than" function for the searches.

          // For an introduction to STL algorithms, see
          // <http://www.learncpp.com/cpp-tutorial/16-4-stl-algorithms-overview/>.
          // For a list of available STL algorithms, see
          // <http://en.cppreference.com/w/cpp/algorithm>

          // We have to create "dummy" time slices for the search.
          sim::TDCIDE startTime;
          sim::TDCIDE endTime;
          startTime.first = start_tdc;
          endTime.first = end_tdc;

          // Here are the fast searches:

          // Find a pointer to the first channel with time >= start_tdc.
          auto const startPointer =
            std::lower_bound(timeSlices.begin(), timeSlices.end(), startTime, TDCIDETimeCompare);

          // From that time slice, find the last channel with time < end_tdc.
          auto const endPointer =
            std::upper_bound(startPointer, timeSlices.end(), endTime, TDCIDETimeCompare);

          // Did we find anything? If not, skip.
          if (startPointer == timeSlices.end() || startPointer == endPointer) continue;
          MF_LOG_DEBUG("MyEnergyAnalysis")
            << "Time slice start = " << (*startPointer).first << std::endl;

          // Loop over the channel times we found that match the hit
          // times.
          for (auto slicePointer = startPointer; slicePointer != endPointer; slicePointer++) {
            auto const timeSlice = *slicePointer;
            auto time = timeSlice.first;

            // How to debug a problem: Lots of print statements. There are
            // debuggers such as gdb, but they can be tricky to use with
            // shared libraries and don't work if you're using software
            // that was compiled somewhere else (e.g., you're accessing
            // LArSoft libraries via CVMFS).

            // The MF_LOG_DEBUG statements below are left over from when I
            // was trying to solve a problem about hit timing. I could
            // have deleted them, but decided to leave them to demonsrate
            // what a typical debugging process looks like.

            MF_LOG_DEBUG("MyEnergyAnalysis")
              << "Hit index = " << hit.LocalIndex() << " channel number = " << hitChannelNumber
              << " start TDC tick = " << hit.StartTick() << " end TDC tick = " << hit.EndTick()
              << " peak TDC tick = " << hit.PeakTime()
              << " sigma peak time = " << hit.SigmaPeakTime()
              << " adjusted start TDC tick = " << clock_data.TPCTick2TDC(hit.StartTick())
              << " adjusted end TDC tick = " << clock_data.TPCTick2TDC(hit.EndTick())
              << " adjusted peak TDC tick = " << clock_data.TPCTick2TDC(hit.PeakTime())
              << " adjusted start_tdc = " << start_tdc << " adjusted end_tdc = " << end_tdc
              << " time = " << time << std::endl;

            // Loop over the energy deposits.
            auto const& energyDeposits = timeSlice.second;
            for (auto const& energyDeposit : energyDeposits) {
              // Remember that map of MCParticles we created
              // near the top of this method? Now we can use
              // it. Search the map for the track ID associated
              // with this energy deposit. Since a map is
              // automatically sorted, we can use a fast binary
              // search method, 'find()'.

              // By the way, the type of "search" is an iterator
              // (to be specific, it's an
              // std::map<int,simb::MCParticle*>::const_iterator,
              // which makes you thankful for the "auto"
              // keyword). If you're going to work with C++
              // containers, you'll have to learn about
              // iterators eventually; do a web search on "STL
              // iterator" to get started.
              auto search = particleMap.find(energyDeposit.trackID);

              // Did we find this track ID in the particle map?
              // It's possible for the answer to be "no"; some
              // particles are too low in kinetic energy to be
              // written by the simulation (see
              // ${LARSIM_DIR}/job/simulationservices.fcl,
              // parameter ParticleKineticEnergyCut).
              if (search == particleMap.end()) continue;

              // "search" points to a pair in the map: <track ID, MCParticle*>
              int trackID = (*search).first;
              const simb::MCParticle& particle = *((*search).second);

              // Is this a primary particle, with a PDG code that
              // matches the user input?
              if (particle.Process() != "primary" || particle.PdgCode() != fSelectedPDG) continue;

              // Determine the dE/dx of this particle.
              const TLorentzVector& positionStart = particle.Position(0);
              TVector3 location(energyDeposit.x, energyDeposit.y, energyDeposit.z);
              double distance = (location - positionStart.Vect()).Mag();
              unsigned int bin = int(distance / fBinSize);
              double energy = energyDeposit.numElectrons * fElectronsToGeV;

              // A feature of maps: if we refer to
              // dEdxMap[trackID], and there's no such entry in
              // the map yet, it will be automatically created
              // with a zero-size vector. Test to see if the
              // vector for this track ID is big enough.
              //
              // dEdxMap is a map, which is a slow container
              // compared to a vector. If we are going to access
              // the same element over and over, it is a good
              // idea to find that element once, and then refer
              // to that item directly. Since we don't care
              // about the specific type of dEdxMap[trackID] (a
              // vector, by the way), we can use "auto" to save
              // some time. This must be a reference, since we
              // want to change the original value in the map,
              // and can't be constant.
              auto& track_dEdX = dEdxMap[trackID];
              if (track_dEdX.size() < bin + 1) {
                // Increase the vector size, padding it with
                // zeroes.
                track_dEdX.resize(bin + 1, 0);
              }

              // Add the energy to the dE/dx bins for this track.
              track_dEdX[bin] += energy;

            } // loop over energy deposits
          }   // loop over time slices
        }     // for each SimChannel
      }       // for each Hit

      // We have a map of dE/dx vectors. Write each one of them to the
      // reconstruction n-tuple.
      for (const auto& dEdxEntry : dEdxMap) {
        // Here, the map entries are <first,second>=<track ID, dE/dx vector>
        fRecoTrackID = dEdxEntry.first;

        // This is an example of how we'd pick out the PDG code if
        // there are multiple particle types or tracks in a single
        // event allowed in the n-tuple.
        fRecoPDG = particleMap[fRecoTrackID]->PdgCode();

        // Get the number of bins for this track.
        const std::vector<double>& dEdx = dEdxEntry.second;
        fRecoNdEdxBins = dEdx.size();

        // Copy this track's dE/dx information.
        fRecodEdxBins = dEdx;

        // At this point, we've filled in all the reconstruction
        // n-tuple's variables. Write it.
        fReconstructionNtuple->Fill();
      }

      // Think about the two big loops above, One starts from the
      // particles then looks at the channels; the other starts with the
      // hits and backtracks to the particles. What links the
      // information in simb::MCParticle and sim::SimChannel is the
      // track ID number assigned by the LArG4 simulation; what links
      // sim::SimChannel and recob::Hit is the channel ID.

      // In general, that is not how objects in the LArSoft
      // reconstruction chain are linked together. Most of them are
      // linked using associations and the art::Assns class:
      // <https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/Using_art_in_LArSoft#artAssns>

      // The web page may a bit difficult to understand (at least, it is
      // for me), so let's try to simplify things:

      // - If you're doing an analysis task, you don't have to worry
      //   about creating art::Assns objects.

      // - You don't have to read the art:Assns objects on your
      //   own. There are helper classes (FindXXX) which will do that for
      //   you.

      // - There's only one helper you need: art::FindManyP. It will do
      //   what you want with a minimum of fuss.

      // - Associations are symmetric. If you see an
      //   art::Assns<ClassA,ClassB>, the helper classes will find all of
      //   ClassA associated with ClassB or all of ClassB associated with
      //   ClassA.

      // - To know what associations exist, you have to be a 'code
      //   detective'. The important clue is to look for a 'produces' line
      //   in the code that writes objects to an art::Event. For example,
      //   in ${LARSIM_DIR}/source/larsim/LArG4/LArG4_module.cc, you'll see this
      //   line:

      //   produces< art::Assns<simb::MCTruth, simb::MCParticle> >();

      //   That means a simulated event will have an association between
      //   simb::MCTruth (the primary particle produced by the event
      //   generator) and the simb::MCParticle objects (the secondary
      //   particles produced in the detector simulation).

      // Let's try it. The following statement will find the
      // simb::MCTruth objects associated with the simb::MCParticle
      // objects in the event (referenced by particleHandle):

      const art::FindManyP<simb::MCTruth> findManyTruth(
        particleHandle, event, fSimulationProducerLabel);

      // Note that we still have to supply the module label of the step
      // that created the association. Also note that we did not have to
      // explicitly read in the simb::MCTruth objects from the
      // art::Event object 'event'; FindManyP did that for us.

      // Also note that at this point art::FindManyP has already found
      // all the simb::MCTruth associated with each of the particles in
      // particleHandle. This is a slow process, so in general you want
      // to do it only once. If we had a loop over the particles, we
      // would still do this outside that loop.

      // Now we can query the 'findManyTruth' object to access the
      // information. First, check that there wasn't some kind of error:

      if (!findManyTruth.isValid()) {
        mf::LogError("MyEnergyAnalysis")
          << "findManyTruth simb::MCTruth for simb::MCParticle failed!";
      }

      // I'm going to be lazy, and just look at the simb::MCTruth object
      // associated with the first simb::MCParticle we read in. (The
      // main reason I'm being lazy is that if I used the
      // single-particle generator in prodsingle.fcl, every particle in
      // the event is going to be associated with just the one primary
      // particle from the event generator.)

      size_t particle_index = 0; // look at first particle in
                                 // particleHandle's vector.

      // I'm using "auto" to save on typing. The result of
      // FindManyP::at() is a vector of pointers, in this case
      // simb::MCTruth*. In this case it will be a vector with just one
      // entry; I could have used art::FindOneP instead. (This will be a
      // vector of art::Ptr, which is a type of smart pointer; see
      // <https://cdcvs.fnal.gov/redmine/projects/larsoft/wiki/Using_art_in_LArSoft#artPtrltTgt-and-artPtrVectorltTgt>
      // To avoid unnecessary copying, and since art::FindManyP returns
      // a constant reference, use "auto const&".

      auto const& truth = findManyTruth.at(particle_index);

      // Make sure there's no problem.
      if (truth.empty()) {
        mf::LogError("MyEnergyAnalysis")
          << "Particle ID=" << particleHandle->at(particle_index).TrackId() << " has no primary!";
      }

      // Use the message facility to write output. I don't have to write
      // the event, run, or subrun numbers; the message facility takes
      // care of that automatically. I'm "going at warp speed" with the
      // vectors, pointers, and methods; see
      // ${NUSIMDATA_INC}/nusimdata/SimulationBase/MCTruth.h and
      // ${NUSIMDATA_INC}/nusimdata/SimulationBase/MCParticle.h for the
      // nested calls I'm using.
      mf::LogInfo("MyEnergyAnalysis")
        << "Particle ID=" << particleHandle->at(particle_index).TrackId()
        << " primary PDG code=" << truth[0]->GetParticle(0).PdgCode();

      // Let's try a slightly more realistic example. Suppose I want to
      // read in the clusters, and learn what hits are associated with
      // them. Then I could backtrack from those hits to determine the
      // dE/dx of the particles in the clusters. (Don't worry; I won't
      // put you through the n-tuple creation procedure for a third
      // time.)

      // First, read in the clusters (if there are any).
      art::Handle<std::vector<recob::Cluster>> clusterHandle;
      if (!event.getByLabel(fClusterProducerLabel, clusterHandle)) return;

      // Now use the associations to find which hits are associated with
      // which clusters. Note that this is not as trivial a query as the
      // one with MCTruth, since it's possible for a hit to be assigned
      // to more than one cluster.

      // We have to include fClusterProducerLabel, since that's the step
      // that created the art::Assns<recob::Hit,recob::Cluster> object;
      // look at the modules in ${LARRECO_DIR}/source/larreco/ClusterFinder/
      // and search for the 'produces' lines. (I did not know this before I
      // wrote these comments. I had to be a code detective and use UNIX
      // tools like 'grep' and 'find' to locate those routines.)
      const art::FindManyP<recob::Hit> findManyHits(clusterHandle, event, fClusterProducerLabel);

      if (!findManyHits.isValid()) {
        mf::LogError("MyEnergyAnalysis") << "findManyHits recob::Hit for recob::Cluster failed;"
                                        << " cluster label='" << fClusterProducerLabel << "'";
      }

      // Now we'll loop over the clusters to see the hits associated
      // with each one. Note that we're not using a range-based for
      // loop. That's because FindManyP::at() expects a number as an
      // argument, so we might as well go through the cluster objects
      // via numerical index instead.
      for (size_t cluster_index = 0; cluster_index != clusterHandle->size(); cluster_index++) {
        // In this case, FindManyP::at() will return a vector of
        // pointers to recob::Hit that corresponds to the
        // "cluster_index"-th cluster.
        auto const& hits = findManyHits.at(cluster_index);

        // We have a vector of pointers to the hits associated
        // with the cluster, but for this example I'm not going to
        // do anything fancy with them. I'll just print out how
        // many there are.
        mf::LogInfo("MyEnergyAnalysis") << "Cluster ID=" << clusterHandle->at(cluster_index).ID()
                                       << " has " << hits.size() << " hits";
      }

    } // MyEnergyAnalysis::analyze()

    // This macro has to be defined for this module to be invoked from a
    // .fcl file; see MyEnergyAnalysis.fcl for more information.
    DEFINE_ART_MODULE(MyEnergyAnalysis)

  } // namespace example
} // namespace lar

// Back to our local namespace.
namespace {

  // Define a local function to calculate the detector diagonal.
  double
  DetectorDiagonal(geo::GeometryCore const& geom)
  {
    const double length = geom.DetLength();
    const double width = 2. * geom.DetHalfWidth();
    const double height = 2. * geom.DetHalfHeight();

    return std::sqrt(cet::sum_of_squares(length, width, height));
  } // DetectorDiagonal()

  // Define a comparison function to use in std::upper_bound and
  // std::lower_bound searches above.
  bool
  TDCIDETimeCompare(const sim::TDCIDE& lhs, const sim::TDCIDE& rhs)
  {
    return lhs.first < rhs.first;
  }
} // local namespace
