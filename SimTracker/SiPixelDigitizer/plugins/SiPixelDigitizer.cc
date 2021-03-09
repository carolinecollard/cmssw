// -*- C++ -*-
//
// Package:    SiPixelDigitizer
// Class:      SiPixelDigitizer
//
/**\class SiPixelDigitizer SiPixelDigitizer.cc SimTracker/SiPixelDigitizer/src/SiPixelDigitizer.cc

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Michele Pioppi-INFN perugia
//   Modifications: Freya Blekman - Cornell University
//         Created:  Mon Sep 26 11:08:32 CEST 2005
//
//

// system include files
#include <memory>
#include <set>

// user include files
#include "SiPixelDigitizer.h"
#include "PixelDigiAddTempInfo.h"
#include "SiPixelDigitizerAlgorithm.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigiCollection.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"
#include "DataFormats/GeometryVector/interface/LocalVector.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "Geometry/CommonTopologies/interface/PixelTopology.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetType.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "SimGeneral/MixingModule/interface/PileUpEventPrincipal.h"
#include "DataFormats/SiPixelDetId/interface/PixelFEDChannel.h"
//Random Number
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Utilities/interface/StreamID.h"
#include "FWCore/Utilities/interface/Exception.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
//using namespace std;

namespace cms {
  SiPixelDigitizer::SiPixelDigitizer(const edm::ParameterSet& iConfig,
                                     edm::ProducesCollector producesCollector,
                                     edm::ConsumesCollector& iC)
      : firstInitializeEvent_(true),
        firstFinalizeEvent_(true),
        _pixeldigialgo(),
        hitsProducer(iConfig.getParameter<std::string>("hitsProducer")),
        trackerContainers(iConfig.getParameter<std::vector<std::string> >("RoutList")),
        geometryType(iConfig.getParameter<std::string>("PixGeometryType")),
        pilotBlades(iConfig.exists("enablePilotBlades") ? iConfig.getParameter<bool>("enablePilotBlades") : false),
        NumberOfEndcapDisks(iConfig.exists("NumPixelEndcap") ? iConfig.getParameter<int>("NumPixelEndcap") : 2) {
    edm::LogInfo("PixelDigitizer ") << "Enter the Pixel Digitizer";

    const std::string alias("simSiPixelDigis");

    producesCollector.produces<edm::DetSetVector<PixelDigi> >().setBranchAlias(alias);
    producesCollector.produces<edm::DetSetVector<PixelDigiSimLink> >().setBranchAlias(alias + "siPixelDigiSimLink");

    for (auto const& trackerContainer : trackerContainers) {
      edm::InputTag tag(hitsProducer, trackerContainer);
      iC.consumes<std::vector<PSimHit> >(edm::InputTag(hitsProducer, trackerContainer));
    }
    edm::Service<edm::RandomNumberGenerator> rng;
    if (!rng.isAvailable()) {
      throw cms::Exception("Configuration")
          << "SiPixelDigitizer requires the RandomNumberGeneratorService\n"
             "which is not present in the configuration file.  You must add the service\n"
             "in the configuration file or remove the modules that require it.";
    }

    _pixeldigialgo = std::make_unique<SiPixelDigitizerAlgorithm>(iConfig);
    if (NumberOfEndcapDisks != 2)
      producesCollector.produces<PixelFEDChannelCollection>();


    // init of counters
    Val1Digi1Simhit=0;
    Val1Digi2Simhit=0;
    Val1Digi3Simhit=0;
    Val1Digi4Simhit=0;

    NduplHit=0;

    NPrim1SimHit=0;
    NPrim2SimHit=0;
    NPrim2SimHit10=0;
    NPrim1if2SimHit=0;
    NStd2SimHit=0;
    NStd3SimHit=0;
    NStd4SimHit=0;
    NnoStd2SimHit=0;
    NnoStd3SimHit=0;
    NnoStd4SimHit=0;


  }

  SiPixelDigitizer::~SiPixelDigitizer() { 

     edm::LogInfo("PixelDigitizer ") << "Destruct the Pixel Digitizer"; 
          std::cout << "  ------------------------------------------ " << std::endl;
          std::cout << " Summary of Investigations in PixelDigitizer " << std::endl;

          std::cout << "  Val1Digi1Simhit ... "  <<  Val1Digi1Simhit << "  " << Val1Digi2Simhit << " " << Val1Digi3Simhit << " " << Val1Digi4Simhit << std::endl;


          // print Duplication results 
          // Val1Digi1Simhit = number of different digis
          // Val1Digi2Simhit = number of digis appearing more than 1
          // Val1Digi3Simhit = number of digis appearing more than 2
          // Val1Digi4Simhit = number of digis appearing more than 3
          int ValExactly1Simhit = Val1Digi1Simhit - Val1Digi2Simhit;
          int ValExactly2Simhit = Val1Digi2Simhit - Val1Digi3Simhit;
          int ValExactly3Simhit = Val1Digi3Simhit - Val1Digi4Simhit;
          //int check_the_sum = ValExactly1Simhit + 2* ValExactly2Simhit + 3* ValExactly3Simhit;
          //if (check_the_sum!=allDigi) std::cout << " !!!!!!!!! problem in decoding the number of digi in multiple apparition " 
          //                                          << check_the_sum << "  != " << allDigi << std::endl;
          std::cout <<  "   RESULTS : Fraction of unique digi-simhit association " << 1.*ValExactly1Simhit/(1.*Val1Digi1Simhit) 
                    << " over " << Val1Digi1Simhit << " digis "<< std::endl; 
          std::cout <<  "   RESULTS : Fraction of double digi-simhit association " << 1.*ValExactly2Simhit/(1.*Val1Digi1Simhit)  << std::endl;
          if (ValExactly3Simhit>0) std::cout <<  "   RESULTS : Fraction of triple  digi-simhit association " << 1.*ValExactly3Simhit/(1.*Val1Digi1Simhit)  << std::endl;
          if (Val1Digi4Simhit>0)   std::cout <<  "   RESULTS : Fraction of >3 digi-simhit association " << 1.*Val1Digi4Simhit/(1.*Val1Digi1Simhit)  << std::endl;


          std::cout <<  "   For information : Duplication of hit " << NduplHit << std::endl;

           // decode extra checks
           //
           std::cout << "   RESULTS : Fraction of Primary for 1st hit " << 1.*NPrim1SimHit/(1.*Val1Digi1Simhit) << std::endl;
           std::cout << "   RESULTS : Fraction of Primary for 1st hit if 2 SimHits " << 1.*NPrim1if2SimHit/(1.*NStd2SimHit+1.*NnoStd2SimHit) << std::endl;
           std::cout << "   RESULTS : Fraction of Primary for 2nd hit " << 1.*NPrim2SimHit/(1.*NStd2SimHit+1.*NnoStd2SimHit) << std::endl;
           std::cout << "   RESULTS : Fraction of Primary for 2nd hit with 1st Primary too " << 1.*NPrim2SimHit10/(1.*NStd2SimHit+1.*NnoStd2SimHit) << std::endl;
           std::cout << "   RESULTS : Fraction of same trackID for 2 SimHits " << 1.*NStd2SimHit/(1.*NStd2SimHit+1.*NnoStd2SimHit) << std::endl;
           if (ValExactly3Simhit>0) std::cout <<  "   RESULTS : Fraction of same trackID for 3 SimHits " << 1.*NStd3SimHit/(1.*NStd3SimHit+1.*NnoStd3SimHit) << std::endl;
           if (Val1Digi4Simhit>0) std::cout <<  "   RESULTS : Fraction of same trackID for >=4 SimHits " << 1.*NStd4SimHit/(1.*NStd4SimHit+1.*NnoStd4SimHit) << std::endl;
    
          std::cout << "  ------------------------------------------ " << std::endl;
  }

  //
  // member functions
  //

  void SiPixelDigitizer::accumulatePixelHits(edm::Handle<std::vector<PSimHit> > hSimHits,
                                             size_t globalSimHitIndex,
                                             const unsigned int tofBin,
                                             edm::EventSetup const& iSetup) {
    if (hSimHits.isValid()) {
      std::set<unsigned int> detIds;
      std::vector<PSimHit> const& simHits = *hSimHits.product();
      edm::ESHandle<TrackerTopology> tTopoHand;
      iSetup.get<TrackerTopologyRcd>().get(tTopoHand);
      const TrackerTopology* tTopo = tTopoHand.product();


      _pixeldigialgo->fillSimHitMaps(simHits, tofBin);

      for (std::vector<PSimHit>::const_iterator it = simHits.begin(), itEnd = simHits.end(); it != itEnd;
           ++it, ++globalSimHitIndex) {
        unsigned int detId = (*it).detUnitId();
        if (detIds.insert(detId).second) {
          // The insert succeeded, so this detector element has not yet been processed.
          assert(detectorUnits[detId]);
          if (detectorUnits[detId] &&
              detectorUnits[detId]
                  ->type()
                  .isTrackerPixel()) {  // this test could be avoided and changed into a check of pixdet!=0
            std::map<unsigned int, PixelGeomDetUnit const*>::iterator itDet = detectorUnits.find(detId);
            if (itDet == detectorUnits.end())
              continue;
            auto pixdet = itDet->second;
            assert(pixdet != nullptr);
            //access to magnetic field in global coordinates
            GlobalVector bfield = pSetup->inTesla(pixdet->surface().position());
            LogDebug("PixelDigitizer ") << "B-field(T) at " << pixdet->surface().position()
                                        << "(cm): " << pSetup->inTesla(pixdet->surface().position());

            std::cout << " in SiPixelDigitizer::accumulatePixelHits  globalSimHitIndex = "  << globalSimHitIndex << std::endl;
            _pixeldigialgo->accumulateSimHits(
                it, itEnd, globalSimHitIndex, tofBin, pixdet, bfield, tTopo, randomEngine_);
          }
        }
      }
    }
  }

  void SiPixelDigitizer::initializeEvent(edm::Event const& e, edm::EventSetup const& iSetup) {
    if (firstInitializeEvent_) {
      _pixeldigialgo->init(iSetup);
      firstInitializeEvent_ = false;
    }

    std::cout << " SiPixelDigitizer initializeEvent " << std::endl;
/*
    Val1Digi1Simhit=0;
    Val1Digi2Simhit=0;
    Val1Digi3Simhit=0;
    Val1Digi4Simhit=0;

    NPrim1SimHit=0;
    NStd2SimHit=0;
    NStd3SimHit=0;
    NStd4SimHit=0;
    NnoStd2SimHit=0;
    NnoStd3SimHit=0;
    NnoStd4SimHit=0;
*/

    // Make sure that the first crossing processed starts indexing the sim hits from zero.
    // This variable is used so that the sim hits from all crossing frames have sequential
    // indices used to create the digi-sim link (if configured to do so) rather than starting
    // from zero for each crossing.
    crossingSimHitIndexOffset_.clear();

    // Cache random number engine
    edm::Service<edm::RandomNumberGenerator> rng;
    randomEngine_ = &rng->getEngine(e.streamID());

    _pixeldigialgo->initializeEvent();
    iSetup.get<TrackerDigiGeometryRecord>().get(geometryType, pDD);
    iSetup.get<IdealMagneticFieldRecord>().get(pSetup);
    edm::ESHandle<TrackerTopology> tTopoHand;
    iSetup.get<TrackerTopologyRcd>().get(tTopoHand);
    const TrackerTopology* tTopo = tTopoHand.product();

    // FIX THIS! We only need to clear and (re)fill this map when the geometry type IOV changes.  Use ESWatcher to determine this.
    if (true) {  // Replace with ESWatcher
      detectorUnits.clear();
      for (const auto& iu : pDD->detUnits()) {
        unsigned int detId = iu->geographicalId().rawId();
        if (iu->type().isTrackerPixel()) {
          auto pixdet = dynamic_cast<const PixelGeomDetUnit*>(iu);
          assert(pixdet != nullptr);
          if (iu->subDetector() ==
              GeomDetEnumerators::SubDetector::PixelEndcap) {  // true ONLY for the phase 0 pixel deetctor
            unsigned int disk = tTopo->layer(detId);           // using the generic layer method
            //if using pilot blades, then allowing it for current detector only
            if ((disk == 3) && ((!pilotBlades) && (NumberOfEndcapDisks == 2)))
              continue;
          }
          detectorUnits.insert(std::make_pair(detId, pixdet));
        }
      }
    }
  }

  void SiPixelDigitizer::accumulate(edm::Event const& iEvent, edm::EventSetup const& iSetup) {
    // Step A: Get Inputs
    for (vstring::const_iterator i = trackerContainers.begin(), iEnd = trackerContainers.end(); i != iEnd; ++i) {
      edm::Handle<std::vector<PSimHit> > simHits;
      edm::InputTag tag(hitsProducer, *i);

      iEvent.getByLabel(tag, simHits);
      unsigned int tofBin = PixelDigiSimLink::LowTof;
      if ((*i).find(std::string("HighTof")) != std::string::npos)
        tofBin = PixelDigiSimLink::HighTof;
             std::cout << "in SiPixelDigitizer::accumulate Event " << tag.encode() << " crossingSimHitIndexOffset_ " << crossingSimHitIndexOffset_[tag.encode()] << std::endl;
      accumulatePixelHits(simHits, crossingSimHitIndexOffset_[tag.encode()], tofBin, iSetup);
      // Now that the hits have been processed, I'll add the amount of hits in this crossing on to
      // the global counter. Next time accumulateStripHits() is called it will count the sim hits
      // as though they were on the end of this collection.
      // Note that this is only used for creating digi-sim links (if configured to do so).
             std::cout << "in SiPixelDigitizer::accumulate Event " << tag.encode() << std::endl;
             std::cout << "index offset, current hit count = " << crossingSimHitIndexOffset_[tag.encode()] << ", " << simHits->size() 
                       << "-> total " << crossingSimHitIndexOffset_[tag.encode()]+simHits->size()  << std::endl;
      if (simHits.isValid())
        crossingSimHitIndexOffset_[tag.encode()] += simHits->size();
    }
  }

  void SiPixelDigitizer::accumulate(PileUpEventPrincipal const& iEvent,
                                    edm::EventSetup const& iSetup,
                                    edm::StreamID const& streamID) {
    // Step A: Get Inputs
    for (vstring::const_iterator i = trackerContainers.begin(), iEnd = trackerContainers.end(); i != iEnd; ++i) {
      edm::Handle<std::vector<PSimHit> > simHits;
      edm::InputTag tag(hitsProducer, *i);

      iEvent.getByLabel(tag, simHits);
      unsigned int tofBin = PixelDigiSimLink::LowTof;
      if ((*i).find(std::string("HighTof")) != std::string::npos)
        tofBin = PixelDigiSimLink::HighTof;
             std::cout << "in SiPixelDigitizer::accumulate PileUpEventPrincipal " << tag.encode() << " crossingSimHitIndexOffset_ " << crossingSimHitIndexOffset_[tag.encode()] << std::endl;

      // test Caro
      if (simHits.isValid()) {
        std::vector<PSimHit> const& ZesimHits = *simHits.product();
        int printnum=0;
        for (std::vector<PSimHit>::const_iterator it = ZesimHits.begin(), itEnd = ZesimHits.end(); it != itEnd;
           ++it, ++printnum) {
                unsigned int detId = (*it).detUnitId();
                std::cout <<  " loop hit " << printnum << " "  << detId << " Eloss  " << (*it).energyLoss() << " Points " << (*it).entryPoint() << " " << (*it).exitPoint() << std::endl;
        }
      }
      // end test




      accumulatePixelHits(simHits, crossingSimHitIndexOffset_[tag.encode()], tofBin, iSetup);
      // Now that the hits have been processed, I'll add the amount of hits in this crossing on to
      // the global counter. Next time accumulateStripHits() is called it will count the sim hits
      // as though they were on the end of this collection.
      // Note that this is only used for creating digi-sim links (if configured to do so).
             std::cout << "in SiPixelDigitizer::accumulate PileUpEventPrincipal "  << iEvent.bunchCrossing() << "  " << tag.encode() << std::endl;
             std::cout << "index offset, current hit count = " << crossingSimHitIndexOffset_[tag.encode()] << ", " << simHits->size() 
                       << "-> total " << crossingSimHitIndexOffset_[tag.encode()]+simHits->size()  << std::endl;
      if (simHits.isValid())
        crossingSimHitIndexOffset_[tag.encode()] += simHits->size();
    }
  }

  // ------------ method called to produce the data  ------------
  void SiPixelDigitizer::finalizeEvent(edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::ESHandle<TrackerTopology> tTopoHand;
    iSetup.get<TrackerTopologyRcd>().get(tTopoHand);
    const TrackerTopology* tTopo = tTopoHand.product();

    std::vector<edm::DetSet<PixelDigi> > theDigiVector;
    std::vector<edm::DetSet<PixelDigiSimLink> > theDigiLinkVector;

    if (firstFinalizeEvent_) {
      const unsigned int bunchspace = PileupInfo_->getMix_bunchSpacing();
      _pixeldigialgo->init_DynIneffDB(iSetup, bunchspace);
      firstFinalizeEvent_ = false;
    }
    _pixeldigialgo->calculateInstlumiFactor(PileupInfo_.get());




    if (_pixeldigialgo->killBadFEDChannels()) {
      std::unique_ptr<PixelFEDChannelCollection> PixelFEDChannelCollection_ =
          _pixeldigialgo->chooseScenario(PileupInfo_.get(), randomEngine_);
      if (PixelFEDChannelCollection_ == nullptr) {
        throw cms::Exception("NullPointerError") << "PixelFEDChannelCollection not set in chooseScenario function.\n";
      }
      iEvent.put(std::move(PixelFEDChannelCollection_));
    }

    for (const auto& iu : pDD->detUnits()) {
      if (iu->type().isTrackerPixel()) {
        //

        edm::DetSet<PixelDigi> collector(iu->geographicalId().rawId());
        edm::DetSet<PixelDigiSimLink> linkcollector(iu->geographicalId().rawId());
//        edm::DetSet<PixelDigiAddTempInfo> tempcollector(iu->geographicalId().rawId());
        std::vector<PixelDigiAddTempInfo> tempcollector;

        _pixeldigialgo->digitize(
            dynamic_cast<const PixelGeomDetUnit*>(iu), collector.data, linkcollector.data,
//                        tempcollector.data, tTopo, randomEngine_);
                        tempcollector, tTopo, randomEngine_);
        if (!collector.data.empty()) {
          theDigiVector.push_back(std::move(collector));
        }
        if (!linkcollector.data.empty()) {
          theDigiLinkVector.push_back(std::move(linkcollector));
        }
//        if (!tempcollector.data.empty()) {
        if (tempcollector.size()>0) {
//          std::cout << " Caro : New class is not empty " << tempcollector.size() << std::endl;
//          std::cout << "  ------------------------------------------ " << std::endl;
          std::vector<PixelDigiAddTempInfo>::const_iterator loopNewClass;
          unsigned int channelPrevious=-1;
          unsigned int channel2Previous=-1;
          unsigned int channel3Previous=-1;
          size_t hitFirstOne=-1;
          int trackIdFirstOne=-1;
          int processTFirstOne=-1;
          int nLoopChan=0;
          int allDigi=0;
          for (loopNewClass = tempcollector.begin(); loopNewClass != tempcollector.end(); ++loopNewClass)  {  // ITERATOR OVER DET IDs
              if (channelPrevious==loopNewClass->channel() && hitFirstOne!=loopNewClass->hitIndex() ) {
                Val1Digi2Simhit++;               
                if (channel2Previous==loopNewClass->channel()) {
                    Val1Digi3Simhit++;
                    if (channel3Previous==loopNewClass->channel())  {
                        Val1Digi4Simhit++;
                    }
                }
                if (trackIdFirstOne==loopNewClass->trackID()) { 
                     if (nLoopChan==1) NStd2SimHit++; 
                     if (nLoopChan==2) NStd3SimHit++; 
                     else NStd4SimHit++; 
                }
                else {
                     if (nLoopChan==1) NnoStd2SimHit++; 
                     if (nLoopChan==2) NnoStd3SimHit++; 
                     else NnoStd4SimHit++; 
                }
                if (nLoopChan==1 && processTFirstOne==0) NPrim1if2SimHit++;
                if (nLoopChan==1 && loopNewClass->processType()==0) NPrim2SimHit++;
                if (nLoopChan==1 && loopNewClass->processType()==0 && processTFirstOne==0) NPrim2SimHit10++;
                nLoopChan++;
                std::cout << " digi " << loopNewClass->channel()  << " appears more than once in the list (" << nLoopChan << "x)" << std::endl;
              }
              else if (channelPrevious==loopNewClass->channel() && hitFirstOne==loopNewClass->hitIndex() ) {
                   NduplHit++;
              }
              else {
                 Val1Digi1Simhit++;
                 // extra info used for checks
                 if (loopNewClass->processType() == 0) NPrim1SimHit++;
                 hitFirstOne= loopNewClass->hitIndex();
                 trackIdFirstOne=loopNewClass->trackID();
                 processTFirstOne=loopNewClass->processType();
                 nLoopChan=1;
              }

              channel3Previous=channel2Previous;
              channel2Previous=channelPrevious;
              channelPrevious=loopNewClass->channel();
              allDigi++;
          }
          if (allDigi!= (int) tempcollector.size()) std::cout << " !!!!!!!  problem : not looping on all the digi of the new class" << std::endl;
  
/*
          std::cout << "  Val1Digi1Simhit ... "  <<  Val1Digi1Simhit << "  " << Val1Digi2Simhit << " " << Val1Digi3Simhit << " " << Val1Digi4Simhit << std::endl;


          // print Duplication results 
          // Val1Digi1Simhit = number of different digis
          // Val1Digi2Simhit = number of digis appearing more than 1
          // Val1Digi3Simhit = number of digis appearing more than 2
          // Val1Digi4Simhit = number of digis appearing more than 3
          int ValExactly1Simhit = Val1Digi1Simhit - Val1Digi2Simhit;
          int ValExactly2Simhit = Val1Digi2Simhit - Val1Digi3Simhit;
          int ValExactly3Simhit = Val1Digi3Simhit - Val1Digi4Simhit;
          //int check_the_sum = ValExactly1Simhit + 2* ValExactly2Simhit + 3* ValExactly3Simhit;
          //if (check_the_sum!=allDigi) std::cout << " !!!!!!!!! problem in decoding the number of digi in multiple apparition " 
          //                                          << check_the_sum << "  != " << allDigi << std::endl;
          std::cout <<  "   RESULTS : Fraction of unique digi-simhit association " << 1.*ValExactly1Simhit/(1.*Val1Digi1Simhit) 
                    << " over " << Val1Digi1Simhit << " digis "<< std::endl; 
          std::cout <<  "   RESULTS : Fraction of double digi-simhit association " << 1.*ValExactly2Simhit/(1.*Val1Digi1Simhit)  << std::endl;
          if (ValExactly3Simhit>0) std::cout <<  "   RESULTS : Fraction of triple  digi-simhit association " << 1.*ValExactly3Simhit/(1.*Val1Digi1Simhit)  << std::endl;
          if (Val1Digi4Simhit>0)   std::cout <<  "   RESULTS : Fraction of >3 digi-simhit association " << 1.*Val1Digi4Simhit/(1.*Val1Digi1Simhit)  << std::endl;


           // decode extra checks
           //
           std::cout << "   RESULTS : Fraction of Primary for 1st hit " << 1.*NPrim1SimHit/(1.*Val1Digi1Simhit) << std::endl;
           std::cout << "   RESULTS : Fraction of same trackID for 2 SimHits " << 1.*NStd2SimHit/(1.*NStd2SimHit+1.*NnoStd2SimHit) << std::endl;
           if (ValExactly3Simhit>0) std::cout <<  "   RESULTS : Fraction of same trackID for 3 SimHits " << 1.*NStd3SimHit/(1.*NStd3SimHit+1.*NnoStd3SimHit) << std::endl;
           if (Val1Digi4Simhit>0) std::cout <<  "   RESULTS : Fraction of same trackID for >=4 SimHits " << 1.*NStd4SimHit/(1.*NStd4SimHit+1.*NnoStd4SimHit) << std::endl;
    
          std::cout << "  ------------------------------------------ " << std::endl;
*/
        }
      }
    }
    _pixeldigialgo->ResetSimHitMaps();


// Test for duplication of digi 

/*
    // std::cout << "  ------------------------------------------ " << std::endl;
    int listDigiSize = theDigiVector.size();
    // std::cout << "  Digi : " << theDigiVector.size()  << std::endl;

    int doubleDigi=0;
    int allDigi=0;
    edm::DetSetVector<PixelDigi>::const_iterator digiIter;
    for (digiIter = theDigiVector.begin(); digiIter != theDigiVector.end(); ++digiIter)  {  // ITERATOR OVER DET IDs
       uint32_t detid = digiIter->id;
       int channelPrevious=-1;
       // std::cout << "  detid  "  << detid  << std::endl;
       edm::DetSet<PixelDigi>::const_iterator ipix;  // ITERATOR OVER DIGI DATA
       for (ipix = digiIter->data.begin(); ipix != digiIter->end(); ++ipix) {
          //std::cout << ipix->row()  << "    "  << ipix->column()  << "    "  <<  ipix->adc() << std::endl;
          if (channelPrevious==ipix->channel()) {
              std::cout << " detid " << detid << "    " <<  ipix->row()  << "    "  << ipix->column()  << " appears twice in the list " << std::endl;
              doubleDigi++;
          }
          channelPrevious=ipix->channel();
          allDigi++;
       }
    }
    std::cout <<  doubleDigi << " double apparitions in the digi list over " << allDigi << " digis (among " << listDigiSize << " different detIDs)" << std::endl;
    // std::cout << "  ------------------------------------------ " << std::endl;
*/

// end of Test

    // Step C: create collection with the cache vector of DetSet
    std::unique_ptr<edm::DetSetVector<PixelDigi> > output(new edm::DetSetVector<PixelDigi>(theDigiVector));
    std::unique_ptr<edm::DetSetVector<PixelDigiSimLink> > outputlink(
        new edm::DetSetVector<PixelDigiSimLink>(theDigiLinkVector));

    // Step D: write output to file
    iEvent.put(std::move(output));
    iEvent.put(std::move(outputlink));

    randomEngine_ = nullptr;  // to prevent access outside event
  }
}  // namespace cms
