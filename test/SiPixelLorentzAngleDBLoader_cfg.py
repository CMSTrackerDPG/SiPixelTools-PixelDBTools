#
import FWCore.ParameterSet.Config as cms
process = cms.Process("SiPixelLorentzAngleLoader")

process.load("Configuration.Geometry.GeometryRecoDB_cff")

process.load("Configuration.StandardSequences.MagneticField_cff")
#process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cfi")
process.load("CalibTracker.Configuration.TrackerAlignment.TrackerAlignment_Fake_cff")
#process.load("Geometry.TrackerGeometryBuilder.trackerGeometry_cfi")
#process.load("Geometry.TrackerNumberingBuilder.trackerNumberingGeometry_cfi")
process.load("CondTools.SiPixel.SiPixelGainCalibrationService_cfi")
process.load("CondCore.DBCommon.CondDBCommon_cfi")



process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
#process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#from Configuration.AlCa.GlobalTag import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run1_data', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_design', '')

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.destinations = cms.untracked.vstring("cout")
process.MessageLogger.cout = cms.untracked.PSet(threshold = cms.untracked.string("ERROR"))

process.source = cms.Source("EmptyIOVSource",
    firstValue = cms.uint64(1),
    lastValue = cms.uint64(1),
    timetype = cms.string('runnumber'),
    interval = cms.uint64(1)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
)

# has to be deleted if it exist  
file = "la.db"
sqlfile = "sqlite_file:" + file
print '\n-> Uploading into file %s, i.e. %s\n' % (file, sqlfile)

##### DATABASE CONNNECTION AND INPUT TAGS ######
process.PoolDBOutputService = cms.Service("PoolDBOutputService",
    BlobStreamerName = cms.untracked.string('TBufferBlobStreamingService'),
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('.'),
        connectionRetrialPeriod = cms.untracked.int32(10),
        idleConnectionCleanupPeriod = cms.untracked.int32(10),
        messageLevel = cms.untracked.int32(1),
        enablePoolAutomaticCleanUp = cms.untracked.bool(False),
        enableConnectionSharing = cms.untracked.bool(True),
        connectionRetrialTimeOut = cms.untracked.int32(60),
        connectionTimeOut = cms.untracked.int32(0),
        enableReadOnlySessionOnUpdateConnection = cms.untracked.bool(False)
    ),
    timetype = cms.untracked.string('runnumber'),
    connect = cms.string(sqlfile),
    toPut = cms.VPSet(
        cms.PSet(
            record = cms.string('SiPixelLorentzAngleRcd'),
            tag = cms.string('SiPixelLorentzAngle_test')
#	     tag = cms.string("SiPixelLorentzAngle_fromAlignment_v01_mc")	
#	     tag = cms.string("SiPixelLorentzAngle_fromAlignment_v01")	
#	     tag = cms.string("SiPixelLorentzAngle_forWidth_v01_mc")
#	     tag = cms.string("SiPixelLorentzAngle_forWidth_v01")
        ),
        cms.PSet(
            record = cms.string('SiPixelLorentzAngleSimRcd'),
            tag = cms.string('SiPixelLorentzAngleSim_test')
        ),
       )
)


###### LORENTZ ANGLE OBJECT ######
process.SiPixelLorentzAngle = cms.EDAnalyzer("SiPixelLorentzAngleDBLoader",
# common input for all rings 
    bPixLorentzAnglePerTesla = cms.double(0.10),
    fPixLorentzAnglePerTesla = cms.double(0.06),
#    bPixLorentzAnglePerTesla = cms.double(0.05),
#    fPixLorentzAnglePerTesla = cms.double(0.03),
# enter -9999 if individual input for rings 
#    bPixLorentzAnglePerTesla = cms.double(-9999.),
#    fPixLorentzAnglePerTesla = cms.double(-9999.),

    #in case of PSet
    BPixParameters = cms.untracked.VPSet(
        cms.PSet(
            layer = cms.uint32(1),
            module = cms.uint32(1),
            angle = cms.double(0.0948)
        ),
        cms.PSet(
            layer = cms.uint32(1),
            module = cms.uint32(2),
            angle = cms.double(0.0948)
        ),
        cms.PSet(
            layer = cms.uint32(1),
            module = cms.uint32(3),
            angle = cms.double(0.0948)
        ),
        cms.PSet(
            layer = cms.uint32(1),
            module = cms.uint32(4),
            angle = cms.double(0.0948)
        ),
        cms.PSet(
            layer = cms.uint32(1),
            module = cms.uint32(5),
            angle = cms.double(0.0964)
        ),
        cms.PSet(
            layer = cms.uint32(1),
            module = cms.uint32(6),
            angle = cms.double(0.0964)
        ),
        cms.PSet(
            layer = cms.uint32(1),
            module = cms.uint32(7),
            angle = cms.double(0.0964)
        ),
        cms.PSet(
            layer = cms.uint32(1),
            module = cms.uint32(8),
            angle = cms.double(0.0964)
        ),
        cms.PSet(
            layer = cms.uint32(2),
            module = cms.uint32(1),
            angle = cms.double(0.0916)
        ),
        cms.PSet(
            layer = cms.uint32(2),
            module = cms.uint32(2),
            angle = cms.double(0.0916)
        ),
        cms.PSet(
            layer = cms.uint32(2),
            module = cms.uint32(3),
            angle = cms.double(0.0916)
        ),
        cms.PSet(
            layer = cms.uint32(2),
            module = cms.uint32(4),
            angle = cms.double(0.0916)
        ),
        cms.PSet(
            layer = cms.uint32(2),
            module = cms.uint32(5),
            angle = cms.double(0.0931)
        ),
        cms.PSet(
            layer = cms.uint32(2),
            module = cms.uint32(6),
            angle = cms.double(0.0931)
        ),
        cms.PSet(
            layer = cms.uint32(2),
            module = cms.uint32(7),
            angle = cms.double(0.0931)
        ),
        cms.PSet(
            layer = cms.uint32(2),
            module = cms.uint32(8),
            angle = cms.double(0.0931)
        ),
        cms.PSet(
            layer = cms.uint32(3),
            module = cms.uint32(1),
            angle = cms.double(0.0920)
        ),
        cms.PSet(
            layer = cms.uint32(3),
            module = cms.uint32(2),
            angle = cms.double(0.0920)
        ),
        cms.PSet(
            layer = cms.uint32(3),
            module = cms.uint32(3),
            angle = cms.double(0.0920)
        ),
        cms.PSet(
            layer = cms.uint32(3),
            module = cms.uint32(4),
            angle = cms.double(0.0920)
        ),
        cms.PSet(
            layer = cms.uint32(3),
            module = cms.uint32(5),
            angle = cms.double(0.0935)
        ),
        cms.PSet(
            layer = cms.uint32(3),
            module = cms.uint32(6),
            angle = cms.double(0.0935)
        ),
        cms.PSet(
            layer = cms.uint32(3),
            module = cms.uint32(7),
            angle = cms.double(0.0935)
        ),
        cms.PSet(
            layer = cms.uint32(3),
            module = cms.uint32(8),
            angle = cms.double(0.0935)
        ),
    ),
    FPixParameters = cms.untracked.VPSet(
        cms.PSet(
            side = cms.uint32(1),
            disk = cms.uint32(1),
            HVgroup = cms.uint32(1),
            angle = cms.double(0.081)
        ),
        cms.PSet(
            side = cms.uint32(1),
            disk = cms.uint32(2),
            HVgroup = cms.uint32(1),
            angle = cms.double(0.081)
        ),
        cms.PSet(
            side = cms.uint32(2),
            disk = cms.uint32(1),
            HVgroup = cms.uint32(1),
            angle = cms.double(0.081)
        ),
        cms.PSet(
            side = cms.uint32(2),
            disk = cms.uint32(2),
            HVgroup = cms.uint32(1),
            angle = cms.double(0.081)
        ),
        cms.PSet(
            side = cms.uint32(1),
            disk = cms.uint32(1),
            HVgroup = cms.uint32(2),
            angle = cms.double(0.081)
        ),
        cms.PSet(
            side = cms.uint32(1),
            disk = cms.uint32(2),
            HVgroup = cms.uint32(2),
            angle = cms.double(0.081)
        ),
        cms.PSet(
            side = cms.uint32(2),
            disk = cms.uint32(1),
            HVgroup = cms.uint32(2),
            angle = cms.double(0.081)
        ),
        cms.PSet(
            side = cms.uint32(2),
            disk = cms.uint32(2),
            HVgroup = cms.uint32(2),
            angle = cms.double(0.081)
        ),
    ),
    #in case lorentz angle values for bpix should be read from file -> not implemented yet
    useFile = cms.bool(False),
    record = cms.untracked.string('SiPixelLorentzAngleRcd'),  
    fileName = cms.string('lorentzFit.txt')	
)

process.SiPixelLorentzAngleSim = cms.EDAnalyzer("SiPixelLorentzAngleDBLoader",
   # magneticField = cms.double(3.8),
    bPixLorentzAnglePerTesla = cms.double(0.10),
    fPixLorentzAnglePerTesla = cms.double(0.06),
    #in case lorentz angle values for bpix should be read from file -> not implemented yet
    useFile = cms.bool(False),
    record = cms.untracked.string('SiPixelLorentzAngleSimRcd'),
    fileName = cms.string('lorentzFit.txt'),	
    #in case of PSet
    BPixParameters = cms.untracked.VPSet(
        cms.PSet(
            layer = cms.uint32(0),
            module = cms.uint32(0),
            angle = cms.double(0.0)
        ),
    ),
    FPixParameters = cms.untracked.VPSet(
        cms.PSet(
            side = cms.uint32(0),
            disk = cms.uint32(0),
            HVgroup = cms.uint32(0),
            angle = cms.double(0.0)
        ),
    ),

)


process.p = cms.Path(
#    process.SiPixelLorentzAngleSim
    process.SiPixelLorentzAngle
    )

