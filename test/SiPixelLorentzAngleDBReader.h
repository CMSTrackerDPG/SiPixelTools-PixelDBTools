#ifndef SiPixelLorentzAngleDBReader_H
#define SiPixelLorentzAngleDBReader_H

#include <string>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TROOT.h"
#include "TFile.h"
#include "TH2F.h"


//
//
// class decleration
//
class SiPixelLorentzAngleDBReader : public edm::one::EDAnalyzer<edm::one::SharedResources> {

  public:
  explicit SiPixelLorentzAngleDBReader(const edm::ParameterSet& );
  ~SiPixelLorentzAngleDBReader() override;
  
  void analyze(const edm::Event&, const edm::EventSetup&) override;

  private:
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tkTopoToken_;
  const edm::ESGetToken<SiPixelLorentzAngle, SiPixelLorentzAngleRcd> siPixelLAToken_;
  const edm::ESGetToken<SiPixelLorentzAngle, SiPixelLorentzAngleSimRcd> siPixelSimLAToken_;
  bool printdebug_;
  TH1F* LorentzAngleBarrel_;
  TH1F* LorentzAngleForward_;
  TH1F *LABPixL1_[8], *LABPixL2_[8], *LABPixL3_[8], *LABPixL4_[8];

  bool useSimRcd_;
  std::string tagLabel_;
  std::vector<unsigned int> l1New, l2New, l3New;
};

#endif
