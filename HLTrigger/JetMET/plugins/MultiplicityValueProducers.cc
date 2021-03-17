#include "MultiplicityValueProducers.h"
#include "DataFormats/Phase2TrackerCluster/interface/Phase2TrackerCluster1D.h"
#include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h"
#include "FWCore/Framework/interface/MakerMacros.h"

typedef MultiplicityValueProducerFromNestedCollection<SiPixelClusterCollectionNew, double>
    SiPixelClusterMultiplicityValueProducer;
DEFINE_FWK_MODULE(SiPixelClusterMultiplicityValueProducer);

typedef MultiplicityValueProducerFromNestedCollection<Phase2TrackerCluster1DCollectionNew, double>
    SiPhase2TrackerClusterMultiplicityValueProducer;
DEFINE_FWK_MODULE(SiPhase2TrackerClusterMultiplicityValueProducer);
