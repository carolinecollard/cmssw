// -*- C++ -*-
//
// Package:    checkAddInfo
// Class:      checkAddInfo
//
/**\class checkAddInfo checkAddInfo.cc 

 Description: Test pixel digis. 
*/
//
// Original Author:  Caroline Collard
//         Created:  March 2021
// inspired by PixelDigiTest.c
//
// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/Common/interface/Handle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/InputTag.h"

// my includes
//#include "Geometry/CommonDetUnit/interface/GeomDet.h"

#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonDetUnit/interface/PixelGeomDetType.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackerCommon/interface/PixelBarrelName.h"

// data formats
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/SiPixelDigi/interface/PixelDigi.h"
#include "DataFormats/DetId/interface/DetId.h"

#include "SimDataFormats/TrackerDigiSimLink/interface/PixelDigiSimLink.h"

// For the big pixel recongnition
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"

// for simulated Tracker hits
//#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

// To use root histos
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

// For ROOT
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TF1.h>
#include <TH2F.h>
#include <TH1F.h>

//#define HISTOS
#define USE_SIM_LINKS

using namespace std;

//
// class decleration
//

class checkAddInfo : public edm::EDAnalyzer {
public:
  explicit checkAddInfo(const edm::ParameterSet &);
  ~checkAddInfo() override;
  void beginJob() override;
  void analyze(const edm::Event &, const edm::EventSetup &) override;
  void endJob() override;

private:
  // ----------member data ---------------------------
  bool PRINT;
  edm::EDGetTokenT<edm::DetSetVector<PixelDigi>> tPixelDigi;

  edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> trackerTopoToken;
  edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> trackerGeomToken;

#ifdef USE_SIM_LINKS
  edm::EDGetTokenT<edm::DetSetVector<PixelDigiSimLink>> tPixelDigiSimLink;
#endif

#ifdef HISTOS

  TH1F *hdetunit;
  TH2F *hpixMap1;

#endif

  edm::InputTag src_;
  int NtotalNumOfDigis = 0;
  int Nn_bingo=0;
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
checkAddInfo::checkAddInfo(const edm::ParameterSet &iConfig) {

  PRINT = iConfig.getUntrackedParameter<bool>("Verbosity", false);
  src_ = iConfig.getParameter<edm::InputTag>("src");
  tPixelDigi = consumes<edm::DetSetVector<PixelDigi>>(src_);

  trackerTopoToken = esConsumes<TrackerTopology, TrackerTopologyRcd>();
  trackerGeomToken = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();

#ifdef USE_SIM_LINKS
  tPixelDigiSimLink = consumes<edm::DetSetVector<PixelDigiSimLink>>(src_);
#endif

  cout << " Construct checkAddInfo " << endl;
}

checkAddInfo::~checkAddInfo() {
  cout << " Destroy checkAddInfo " << endl;
}

//
// member functions
//
// ------------ method called at the begining   ------------
void checkAddInfo::beginJob() {
  using namespace edm;
  cout << "Initialize checkAddInfo " << endl;

#ifdef HISTOS
  edm::Service<TFileService> fs;

  // example of histograms
  hdetunit = fs->make<TH1F>("hdetunit", "Det unit", 1000, 302000000., 302300000.);
  hpixMap1 = fs->make<TH2F>("hpixMap1", " ", 416, 0., 416., 160, 0., 160.);
  hpixMap1->SetOption("colz");

#endif
}

// ------------ method called to produce the data  ------------
void checkAddInfo::analyze(const edm::Event &iEvent, const edm::EventSetup &iSetup) {
  //Retrieve tracker topology from geometry
  edm::ESHandle<TrackerTopology> tTopo = iSetup.getHandle(trackerTopoToken);

  using namespace edm;
  if (PRINT)
    cout << " Analyze checkAddInfo " << endl;

  int run       = iEvent.id().run();
  int event = iEvent.id().event();

  // Get digis
  edm::Handle<edm::DetSetVector<PixelDigi>> pixelDigis;
  iEvent.getByToken(tPixelDigi, pixelDigis);

#ifdef USE_SIM_LINKS
  // Get simlink data
  edm::Handle<edm::DetSetVector<PixelDigiSimLink>> pixelSimLinks;
  iEvent.getByToken(tPixelDigiSimLink, pixelSimLinks);
#endif

  // Get event setup (to get global transformation)
  edm::ESHandle<TrackerGeometry> geom = iSetup.getHandle(trackerGeomToken);
  const TrackerGeometry &theTracker(*geom);

  int numberOfDetUnits = 0;
  int totalNumOfDigis = 0;
  int totalNumOfSimhits = 0;
  int n_bingo=0;

  // Iterate on detector units
  edm::DetSetVector<PixelDigi>::const_iterator DSViter;
  for (DSViter = pixelDigis->begin(); DSViter != pixelDigis->end(); DSViter++) {
    bool valid = false;
    unsigned int detid = DSViter->id;  // = rawid
    DetId detId(detid);
    //const GeomDetUnit      * geoUnit = geom->idToDetUnit( detId );
    //const PixelGeomDetUnit * pixDet  = dynamic_cast<const PixelGeomDetUnit*>(geoUnit);
    unsigned int detType = detId.det();     // det type, tracker=1
    unsigned int subid = detId.subdetId();  //subdetector type, barrel=1

//    if (PRINT)
//      cout << "Det: " << detId.rawId() << " " << detId.null() << " " << detType << " " << subid << endl;

#ifdef HISTOS
    hdetunit->Fill(float(detid));
#endif

    if (detType != 1)
      continue;  // look only at tracker
    ++numberOfDetUnits;


#ifdef USE_SIM_LINKS
    // Look at simlink information (simulated data only)
    //
        int numberOfSimLinks = 0;
    edm::DetSetVector<PixelDigiSimLink>::const_iterator isearch = pixelSimLinks->find(detid);

    if (isearch != pixelSimLinks->end()) {  //if it is not empty
      edm::DetSet<PixelDigiSimLink> link_detset = (*pixelSimLinks)[detid];
      edm::DetSet<PixelDigiSimLink>::const_iterator it;
      // Loop over DigisSimLink in this det unit
      for (it = link_detset.data.begin(); it != link_detset.data.end(); it++) {
        numberOfSimLinks++;
        totalNumOfSimhits++;
        unsigned int chan = it->channel();
        unsigned int simTrack = it->SimTrackId();
        float frac = it->fraction();
        unsigned int rawID = it->eventId().rawId();

//        if (PRINT) cout << " Sim link " << numberOfSimLinks << " " << chan << " " << simTrack << " " << frac << "  rawID " << rawID <<  endl;
        
        // test if there are duplicated channel 
//        int numbloop=0;
        edm::DetSet<PixelDigiSimLink>::const_iterator it2;
        for (it2=link_detset.data.begin(); it2 != it ; it2++) {
//           numbloop++;
//           if (PRINT) cout << " debug loop interne " <<  numbloop << " " << it2->channel() << " " << it2->SimTrackId() << " " << it2->fraction() << "  rawID " << it2->eventId().rawId() <<  endl;
           if (it2->channel() == chan && it2->eventId().rawId() !=rawID) {
//              if (PRINT) cout << " ah ah ah -----  BINGO ------ "  << it2->channel() << " " << chan << " " << it2->eventId().rawId() << " " << rawID 
//                              << " for detID " << detid << " in event " << event <<  endl; 
              n_bingo++;
              Nn_bingo++;
           }
        }
       
      }  // end simlink det loop
    }  // end simlink if

#endif  // USE_SIM_LINKS

    unsigned int numberOfDigis = 0;

    // Look at digis now
    edm::DetSet<PixelDigi>::const_iterator di;
    for (di = DSViter->data.begin(); di != DSViter->data.end(); di++) {
      numberOfDigis++;
      totalNumOfDigis++;
      NtotalNumOfDigis++;
      int adc = di->adc();     // charge, modifued to unsiged short
      int col = di->column();  // column
      int row = di->row();     // row
      int layer=0;
      int ladder=0;
      int zindex=0;

      if (subid == 1) { 
        // Convert to online
        PixelBarrelName pbn(detid);
        layer = pbn.layerName();
        ladder = pbn.ladderName();
        zindex = tTopo->pxfModule(detid); 
      }

      // channel index needed to look for the simlink to simtracks
      int channel = PixelChannelIdentifier::pixelToChannel(row, col);
//      if (PRINT)
//        cout << numberOfDigis << " Col: " << col << " Row: " << row << " ADC: " << adc << " channel = " << channel
//             << endl;

      if (col > 415)
        cout << " Error: column index too large " << col << " Barrel layer, ladder, module " << layer << " " << ladder
             << " " << zindex << endl;
      if (row > 159)
        cout << " Error: row index too large " << row << endl;

#ifdef HISTOS
      if (layer == 1) {
          hpixMap1->Fill(float(col), float(row));
      } 
#endif
    }

  }  // end for det-units

  if (PRINT) {
    cout << " Number of full det-units = " << numberOfDetUnits << " total digis = " << totalNumOfDigis ;
#ifdef USE_SIM_LINKS
    cout << " total simhits = " << totalNumOfSimhits << endl;
    cout << " test frequency of 2 PU giving signal in the same digi " << (n_bingo*1.)/(totalNumOfDigis*1.); 
    
#endif
    cout <<  endl;

  }


}
// ------------ method called to at the end of the job  ------------
void checkAddInfo::endJob() {
#ifdef USE_SIM_LINKS
    cout << " Average frequency of 2 PU giving signal in the same digi " << (Nn_bingo*1.)/(NtotalNumOfDigis*1.) << " over " << NtotalNumOfDigis << " digis "<< endl;
    
#endif
  cout << " End checkAddInfo " << endl;
  //hFile->Write();
  //hFile->Close();
}

//define this as a plug-in
DEFINE_FWK_MODULE(checkAddInfo);
