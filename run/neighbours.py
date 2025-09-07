from Gaudi.Configuration import *

# DD4hep geometry service
from Configurables import GeoSvc
geoservice = GeoSvc("GeoSvc")
path_to_detector = os.environ.get("K4GEO", "")
detectors_to_use = [
    'FCCee/ALLEGRO/compact/ALLEGRO_o1_v04/ALLEGRO_o1_v04.xml'
]
geoservice.detectors = [
    os.path.join(path_to_detector, _det) for _det in detectors_to_use
]
geoservice.OutputLevel = INFO


# Geant4 service
# Configures the Geant simulation: geometry, physics list and user actions
from Configurables import CreateFCCeeCaloNeighbours
neighbours = CreateFCCeeCaloNeighbours("neighbours", 
#                                       outputFileName="neighbours_map_ecalB_thetamodulemerged_hcalB_hcalEndcap.root",
                                       outputFileName="map.root",
                                       readoutNames=["ECalBarrelModuleThetaMerged","HCalBarrelReadout","HCalEndcapReadout"],
                                       systemNames=["system","system","system"],
                                       systemValues=[4,8,9],
                                       activeFieldNames=["layer","layer","layer"],
                                       activeVolumesNumbers=[11,13,37],
                                       connectBarrels=True, 
                                       connectHCal=True,
                                       OutputLevel=INFO)

# ApplicationMgr
from Configurables import ApplicationMgr
ApplicationMgr( TopAlg = [],
                EvtSel = 'NONE',
                EvtMax   = 1,
                # order is important, as GeoSvc is needed by G4SimSvc
                ExtSvc = [geoservice, neighbours],
                OutputLevel=INFO
)
