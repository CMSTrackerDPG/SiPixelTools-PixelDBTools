#ifndef SiPixelLorentzAngleDBLoader_h
#define SiPixelLorentzAngleDBLoader_h

#include <map>

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/DetId/interface/DetId.h"

// Magnetic field
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"


class SiPixelLorentzAngleDBLoader : public edm::one::EDAnalyzer<>
{
 public:
  
  explicit SiPixelLorentzAngleDBLoader(const edm::ParameterSet& conf);
  
  ~SiPixelLorentzAngleDBLoader() override;  
  void analyze(const edm::Event& e, const edm::EventSetup& c) override;
  
  
 private:
  const edm::ESGetToken<TrackerGeometry, TrackerDigiGeometryRecord> tkGeomToken_;
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tkTopoToken_;
  int HVgroup(int panel, int module);

  std::vector< std::pair<uint32_t, float> > detid_la;
  edm::ParameterSet conf_;
  //double magneticField_;
  std::string recordName_;

  typedef std::vector< edm::ParameterSet > Parameters;
  Parameters BPixParameters_;
  Parameters FPixParameters_;
  Parameters ModuleParameters_;

  float bPixLorentzAnglePerTesla_;
  float fPixLorentzAnglePerTesla_;

  std::string fileName_;
  bool useFile_;
};


#endif
