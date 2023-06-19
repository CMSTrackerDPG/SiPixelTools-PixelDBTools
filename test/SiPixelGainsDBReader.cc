#include <memory>

#include "SiPixelGainsDBReader.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
//#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/CommonTopologies/interface/PixelTopology.h"

#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationForHLTService.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationForHLTSimService.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationOfflineService.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationOfflineSimService.h"
#include "CalibTracker/SiPixelESProducers/interface/SiPixelGainCalibrationService.h"

#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "DataFormats/SiPixelDetId/interface/PixelSubdetector.h"
#include "DataFormats/TrackerCommon/interface/PixelBarrelName.h"
#include "DataFormats/TrackerCommon/interface/PixelEndcapName.h"

#include "DataFormats/DetId/interface/DetId.h"


using namespace std;

SiPixelGainsDBReader::SiPixelGainsDBReader(const edm::ParameterSet& conf): conf_(conf) {

  std::string payloadType = conf_.getParameter<std::string>("payloadType");
  PRINT       = conf_.getParameter<bool>("verbose");
  bool simRcd = conf_.getParameter<bool>("useSimRcd");
  vcalIncluded = conf_.getUntrackedParameter<bool>("vcalIncluded",true);

  cout<<" type = "<<payloadType<<" sim/reco = "<< simRcd<<" PRINT = "<<PRINT <<endl;

  trackerTopoToken_ = esConsumes<TrackerTopology, TrackerTopologyRcd>();
  trackerGeomToken_ = esConsumes<TrackerGeometry, TrackerDigiGeometryRecord>();

  if (strcmp(payloadType.c_str(), "HLT") == 0) { // HLT
    if( simRcd )  //Sim
      SiPixelGainCalibrationService_ = new  SiPixelGainCalibrationForHLTSimService(conf_,consumesCollector());
    else
      SiPixelGainCalibrationService_ = new SiPixelGainCalibrationForHLTService(conf_,consumesCollector());

  } else if (strcmp(payloadType.c_str(), "Offline") == 0) {  // Offline 
    if(simRcd)  //Sim
      SiPixelGainCalibrationService_ = new SiPixelGainCalibrationOfflineSimService(conf_,consumesCollector());
    else
      SiPixelGainCalibrationService_ = new SiPixelGainCalibrationOfflineService(conf_,consumesCollector());

  } else if (strcmp(payloadType.c_str(), "Full") == 0) { //Full, not used at the moment
    SiPixelGainCalibrationService_ = new SiPixelGainCalibrationService(conf,consumesCollector());
  }
  
}

void
SiPixelGainsDBReader::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  
  //const bool PRINT = false;
  const bool phase1=true;
  cout<<PRINT<<endl;

  //Create Subdirectories
  edm::Service<TFileService> fs;
  TFileDirectory subDirPed = fs->mkdir("Pedestals");
  TFileDirectory subDirGain = fs->mkdir("Gains");
  char name[128];

  unsigned int nmodules = 0, nmodulesB=0, nmodulesF=0;
  uint32_t nchannels = 0, nchannelsB=0, nchannelsF=0;
  uint32_t ndead=0, ndeadInL1=0,ndeadInL2=0,ndeadInL3=0,ndeadInL4=0,ndeadInF=0;
  uint32_t nnoisy=0;
  
  // Get the calibration data
  SiPixelGainCalibrationService_->setESObjects(iSetup);
  edm::LogInfo("SiPixelGainsDBReader") << "[SiPixelGainsDBReader::beginJob] End Reading Gainsects" << std::endl;


  // Get the Geometry
  //iSetup.get<TrackerDigiGeometryRecord>().get( tkgeom );     
  edm::ESHandle<TrackerGeometry> tkgeom = iSetup.getHandle(trackerGeomToken_);
  //const TrackerGeometry &tkgeom(*geom);

  std::cout <<" There are "<<tkgeom->dets().size() <<" detectors"<<std::endl;

  //Retrieve tracker topology from geometry
  //edm::ESHandle<TrackerTopology> tTopoH;
  //iSetup.get<TrackerTopologyRcd>().get(tTopoH);
  edm::ESHandle<TrackerTopology> tTopoH = iSetup.getHandle(trackerTopoToken_);
  const TrackerTopology *tTopo=tTopoH.product();

  // Get the list of DetId's
  std::vector<uint32_t> vdetId_ = SiPixelGainCalibrationService_->getDetIds();

  // range of gain histos 
  float maxGain = 10.;
  if(vcalIncluded) maxGain=500.;
  //Create histograms
  _TH1F_Dead_sum =  fs->make<TH1F>("Summary_dead","Dead pixel fraction (0=dead, 1=alive)",vdetId_.size()+1,0,vdetId_.size()+1);
  _TH1F_Dead_all =  fs->make<TH1F>("DeadAll","Dead pixel fraction (0=dead, 1=alive)",50,0.,conf_.getUntrackedParameter<double>("maxRangeDeadPixHist",0.001));
  _TH1F_Noisy_sum = fs->make<TH1F>("Summary_noisy","Noisy pixel fraction (0=noisy, 1=alive)",vdetId_.size()+1,0,vdetId_.size()+1);
  _TH1F_Noisy_all = fs->make<TH1F>("NoisyAll","Noisy pixel fraction (0=noisy, 1=alive)",50,0.,conf_.getUntrackedParameter<double>("maxRangeDeadPixHist",0.001));
  _TH1F_Gains_sum =  fs->make<TH1F>("Summary_Gain","Gain Summary", vdetId_.size()+1,0,vdetId_.size()+1);
  _TH1F_Pedestals_sum =  fs->make<TH1F>("Summary_Pedestal","Pedestal Summary", vdetId_.size()+1,0,vdetId_.size()+1);
  _TH1F_Pedestals_all = fs->make<TH1F>("PedestalsAll","all Pedestals",500,-250,250);
  _TH1F_Pedestals_bpix = fs->make<TH1F>("PedestalsBpix","bpix Pedestals",500,-250,250);
  _TH1F_Pedestals_fpix = fs->make<TH1F>("PedestalsFpix","fpix Pedestals",500,-250,250);
  _TH1F_Gains_all = fs->make<TH1F>("GainsAll","all Gains",500,0,maxGain);
  _TH1F_Gains_bpix = fs->make<TH1F>("GainsBpix", "bpix Gains", 500, 0, maxGain);
  _TH1F_Gains_fpix = fs->make<TH1F>("GainsFpix", "fpix Gains", 500, 0, maxGain);

  GainsL1 = fs->make<TH1D>("GainsL1", "Gains L1", 500, 0, maxGain);
  GainsL2 = fs->make<TH1D>("GainsL2", "Gains L2", 500, 0, maxGain);
  GainsL3 = fs->make<TH1D>("GainsL3", "Gains L3", 500, 0, maxGain);
  GainsL4 = fs->make<TH1D>("GainsL4", "Gains L4", 500, 0, maxGain);
  PedsL1 = fs->make<TH1D>("PedsL1", "Pedestals L1", 500, -250, 250);
  PedsL2 = fs->make<TH1D>("PedsL2", "Pedestals L2", 500, -250, 250);
  PedsL3 = fs->make<TH1D>("PedsL3", "Pedestals L3", 500, -250, 250);
  PedsL4 = fs->make<TH1D>("PedsL4", "Pedestals L4", 500, -250, 250);
  ADCL1 = fs->make<TH1D>("ADCL1", "ADC@elec=2000 L1", 500,0, 100);
  ADCL2 = fs->make<TH1D>("ADCL2", "ADC@elec=2000 L2", 500, 0, 100);
  ADCL3 = fs->make<TH1D>("ADCL3", "ADC@elec=2000 L3", 500, 0, 100);
  ADCL4 = fs->make<TH1D>("ADCL4", "ADC@elec=2000 L4", 500, 0, 100);
  ADCF  = fs->make<TH1D>("ADCF", "ADC@elec=2000 FPix", 500, 0, 100);
  VcalL1 = fs->make<TH1D>("VcalL1","elec@adc=0 L1", 500, -7500, 2500);
  VcalL2 = fs->make<TH1D>("VcalL2","elec@adc=0 L2", 500, -7500, 2500);
  VcalL3 = fs->make<TH1D>("VcalL3","elec@adc=0 L3", 500, -7500, 2500);
  VcalL4 = fs->make<TH1D>("VcalL4","elec@adc=0 L4", 500, -7500, 2500);
  VcalF  = fs->make<TH1D>("VcalF", "elec@adc=0 FPix",500,-7500, 2500);

  GainsSumL1 = fs->make<TH1D>("GainsSumL1", "Gains av L1", 500, 0,maxGain);
  GainsSumL2 = fs->make<TH1D>("GainsSumL2", "Gains av L2", 500, 0,maxGain);
  GainsSumL3 = fs->make<TH1D>("GainsSumL3", "Gains av L3", 500, 0,maxGain);
  GainsSumL4 = fs->make<TH1D>("GainsSumL4", "Gains av L4", 500, 0,maxGain);
  GainsSumF = fs->make<TH1D>("GainsSumF", "Gains av FPix", 500, 0,maxGain);
  PedsSumL1 = fs->make<TH1D>("PedsSumL1", "Pedestals av L1", 500, -250.,250.);
  PedsSumL2 = fs->make<TH1D>("PedsSumL2", "Pedestals av L2", 500, -250.,250.);
  PedsSumL3 = fs->make<TH1D>("PedsSumL3", "Pedestals av L3", 500, -250.,250.);
  PedsSumL4 = fs->make<TH1D>("PedsSumL4", "Pedestals av L4", 500, -250.,250.);
  PedsSumF  = fs->make<TH1D>("PedsSumF", "Pedestals av FPix",500, -250.,250.);
 
  ADCSumL1 = fs->make<TH1D>("ADCSumL1", "ADC@elec=2000 av L1", 500, 0, 100);
  ADCSumL2 = fs->make<TH1D>("ADCSumL2", "ADC@elec=2000 av L2", 500, 0, 100);
  ADCSumL3 = fs->make<TH1D>("ADCSumL3", "ADC@elec=2000 av L3", 500, 0, 100);
  ADCSumL4 = fs->make<TH1D>("ADCSumL4", "ADC@elec=2000 av L4", 500, 0, 100);
  ADCSumF = fs->make<TH1D>("ADCSumF", "ADC@elec=2000 av FPix", 500, 0, 100);
  VcalSumL1 = fs->make<TH1D>("VcalSumL1", "VCAL@adc=0 av L1", 500, -7500, 2500);
  VcalSumL2 = fs->make<TH1D>("VcalSumL2", "VCAL@adc=0 av L2", 500, -7500, 2500);
  VcalSumL3 = fs->make<TH1D>("VcalSumL3", "VCAL@adc=0 av L3", 500, -7500, 2500);
  VcalSumL4 = fs->make<TH1D>("VcalSumL4", "VCAL@adc=0 av L4", 500, -7500, 2500);
  VcalSumF = fs->make<TH1D>("VcalSumF", "VCAL@adc=0 av FPix", 500, -7500, 2500);

  const bool fillTree = false;
  uint32_t detid;
  double gainmeanfortree, gainrmsfortree, pedmeanfortree, pedrmsfortree;
  TTree *tree = new TTree("tree","tree");
  if(fillTree) {
    tree->Branch("detid",&detid,"detid/I");
    tree->Branch("ped_mean",&pedmeanfortree,"ped_mean/D");
    tree->Branch("ped_rms",&pedrmsfortree,"ped_rms/D");
    tree->Branch("gain_mean",&gainmeanfortree,"gain_mean/D");
    tree->Branch("gain_rms",&gainrmsfortree,"gain_rms/D");
  }
  // Loop over DetId's
  int ibin = 1;
  for (std::vector<uint32_t>::const_iterator detid_iter=vdetId_.begin();detid_iter!=vdetId_.end();detid_iter++){
    bool select = false;
    detid = *detid_iter;
    //cout<<detid<<endl;

    sprintf(name,"Pedestals_%d",detid);
    _TH1F_Pedestals_m[detid] = subDirPed.make<TH1F>(name,name,350,-100.,250.);    
    sprintf(name,"Gains_%d",detid);
    _TH1F_Gains_m[detid] = subDirGain.make<TH1F>(name,name,100,0.,maxGain); 

    DetId detIdObject(detid);
    const PixelGeomDetUnit* _PixelGeomDetUnit = dynamic_cast<const PixelGeomDetUnit*>(tkgeom->idToDetUnit(DetId(detid)));
    if (_PixelGeomDetUnit==0){
      edm::LogError("SiPixelCondObjHLTDisplay")<<"[SiPixelCondObjHLTReader::beginJob] the detID "<<detid<<" doesn't seem to belong to Tracker"<<std::endl; 
      continue;
    }     
         
    _deadfrac_m[detid]=0.;
    _noisyfrac_m[detid]=0.;

    nmodules++;

    const GeomDetUnit      * geoUnit = tkgeom->idToDetUnit( detIdObject );
    const PixelGeomDetUnit * pixDet  = dynamic_cast<const PixelGeomDetUnit*>(geoUnit);
    const PixelTopology & topol = pixDet->specificTopology();       

    int layer = -1;
    string name;
    if ( detIdObject.subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel) ) { // BPIX

      nmodulesB++;
      // Barell layer = 1,2,3
      unsigned int layerC1 = tTopo->pxbLayer(detIdObject);
      // Barrel ladder id 1-20,32,44.
      unsigned int ladderC1 = tTopo->pxbLadder(detIdObject);
      // Barrel Z-index=1,8
      unsigned int zindex1 = tTopo->pxbModule(detIdObject);

      // example printouts
      if(layerC1==1 && ladderC1==1 && zindex1==1) select=true;
      if(layerC1==2 && ladderC1==1 && zindex1==1) select=true;
      if(layerC1==3 && ladderC1==1 && zindex1==1) select=true;
      if(layerC1==4 && ladderC1==1 && zindex1==1) select=true;

      // Online 
      PixelBarrelName pbn(detIdObject,tTopo,phase1); // use det-id
      PixelBarrelName::Shell sh = pbn.shell(); //enum
      int sector = pbn.sectorName();
      int ladder = pbn.ladderName();
      layer  = pbn.layerName();
      int module = pbn.moduleName();
      //bool half  = pbn.isHalfModule();
      name= pbn.name();
      //PixelModuleName::ModuleType moduleType = pbn.moduleType();

      // Shell { mO = 1, mI = 2 , pO =3 , pI =4 };
      int shell = int(sh);
      // change the module sign for z<0
      if(shell==1 || shell==2) module = -module;
      // change ladeer sign for Outer )x<0)
      if(shell==1 || shell==3) ladder = -ladder;

      //if(layer==1 && ladder==1 && module==1) select=true;
      if(PRINT&&select) cout<<" shell "<<sh<<"("<<shell<<") "<<sector<<" "<<layer<<" "<<ladder<<" "
			    <<module<<" "<<name<<endl;

    } else {  // endcaps

      nmodulesF++;
      // ids 
      int disk=tTopo->pxfDisk(detIdObject); //1,2,3
      int blade=tTopo->pxfBlade(detIdObject); //1-24
      int plaq=tTopo->pxfModule(detIdObject); //
      int side=tTopo->pxfSide(detIdObject); //size=1 for -z, 2 for +z
      int panel=tTopo->pxfPanel(detIdObject); //panel=1

      //if(disk==1 && blade==1 && side==1 && panel==1 && plaq==1 ) select=true;
      //if(disk==2 && blade==1 && side==1 && panel==1 && plaq==1 ) select=true;
      //if(disk==3 && blade==1 && side==1 && panel==1 && plaq==1 ) select=true;
      if(PRINT&&select) cout<<"endcap, size "<<side<<" disk "<<disk<<", blade "
		    <<blade<<", panel "
		    <<panel<<", plaq "<<plaq<<endl;
 

      // Convert to online 
      //PixelEndcapName pen(detIdObject,tTopo,phase1); // use det-id
      //string nameF = pen.name();
      //int diskName = pen.diskName();
      //int bladeName = pen.bladeName();
      //int pannelName = pen.pannelName();
      //int plaquetteName = pen.plaquetteName();
      //PixelEndcapName::HalfCylinder part = pen.halfCylinder();
      //PixelModuleName::ModuleType moduleType = pen.moduleType();
      //if(PRINT&&select) cout<<nameF<<" "<<diskName<<" "<<bladeName<<" "<<pannelName<<" "
      //	    <<plaquetteName<<" "<<part<<" "<<moduleType<<endl;

    }

    // Get the module sizes.
    int nrows = topol.nrows();      // rows in x
    int ncols = topol.ncolumns();   // cols in y
    float nchannelspermod=0, nchannelsperroc=0;
    float gains=0, peds=0;
    int ndeadInModule=0;
    int rowMin=0, rowMax=0;
    //cout<<detid<<endl;

    int rocIdOld=-1, layerOld=-1;
    for(int n=0; n<2; ++n) {
      if(n==0) {rowMin=0;  rowMax=80;}
      else     {rowMin=80; rowMax=160;}
      for(int col_iter=0; col_iter<ncols; col_iter++) {
	for(int row_iter=rowMin; row_iter<rowMax; row_iter++) {
	  int rocId = int(col_iter/52) + int(row_iter/80)*8;
	  //cout<<detid<<" "<<col_iter<<" "<<row_iter<<" "<<rocId<<" "<<rocIdOld<<" "<<nchannelsperroc<<endl;
	  if( (rocId!=rocIdOld) && nchannelsperroc>0) {  // new roc
	    //cout<<"new roc "<<endl; int dum; cin>>dum;
	    //if(dum==-1) return;

	    gains /=float(nchannelsperroc);
	    peds /=float(nchannelsperroc);
	    float elecs = -1. * peds * gains;
	    // adc for 2000e (at threshold)
	    float adcs = (2000/gains) + peds;
	    
	    // ladder is wrong here, review
	    if(layerOld>0) { // BPIX
	      if(layerOld==1)      { GainsSumL1->Fill(gains); PedsSumL1->Fill(peds);ADCSumL1->Fill(adcs);VcalSumL1->Fill(elecs);
	      }else if(layerOld==2) { GainsSumL2->Fill(gains); PedsSumL2->Fill(peds);ADCSumL2->Fill(adcs);VcalSumL2->Fill(elecs);}
	      else if(layerOld==3) { GainsSumL3->Fill(gains); PedsSumL3->Fill(peds);ADCSumL3->Fill(adcs);VcalSumL3->Fill(elecs);}
	      else if(layerOld==4) { GainsSumL4->Fill(gains); PedsSumL4->Fill(peds);ADCSumL4->Fill(adcs);VcalSumL4->Fill(elecs);}
	    } else { GainsSumF->Fill(gains);  PedsSumF->Fill(peds); ADCSumF->Fill(adcs); VcalSumF->Fill(elecs);}
	    nchannelsperroc=0;
	    gains=0.; peds=0.;
	  } // new roc 

	  rocIdOld=rocId; layerOld=layer;
	nchannelsperroc++;
	nchannelspermod++;
	nchannels++;
	 
	 if(SiPixelGainCalibrationService_->isDead(detid,col_iter,row_iter)) {
	    //	    std::cout << "found dead pixel " << detid << " " <<col_iter << "," << row_iter << std::endl;
	   ndead++;
	   ndeadInModule++;
	   _deadfrac_m[detid]++;
	   continue;

	 } else if(SiPixelGainCalibrationService_->isNoisy(detid,col_iter,row_iter)) {
	    //	    std::cout << "found noisy pixel " << detid << " " <<col_iter << "," << row_iter << std::endl;
	   nnoisy++;
	   _noisyfrac_m[detid]++;
	   continue;
	 }

	 float gain  = SiPixelGainCalibrationService_->getGain(detid, col_iter, row_iter);
	 _TH1F_Gains_m[detid]->Fill( gain );
	 _TH1F_Gains_all->Fill(gain);
	 gains += gain;

	 //if (detIdObject.subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel))
	 // _TH1F_Gains_bpix->Fill(gain);
	 //if (detIdObject.subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap))
	 //_TH1F_Gains_fpix->Fill(gain);

	 float ped  = SiPixelGainCalibrationService_->getPedestal(detid, col_iter, row_iter);
	 _TH1F_Pedestals_m[detid]->Fill( ped );
       	 _TH1F_Pedestals_all->Fill(ped);
	 peds += ped;

         //if (detIdObject.subdetId() == static_cast<int>(PixelSubdetector::PixelBarrel))
	 //_TH1F_Pedestals_bpix->Fill(ped);
	 //if (detIdObject.subdetId() == static_cast<int>(PixelSubdetector::PixelEndcap))
	 //_TH1F_Pedestals_fpix->Fill(ped);

	 // electrons = gain * (adc - pedestal)
	 // adc =  electrons/gain + pedestal
	 // electrons for adc=0
	 float elec = -1. * ped * gain;
	 // adc for 2000e (at threshold)
	 float adc = (2000/gain) + ped;

	 if(layer>0) { // BPIX
	   _TH1F_Gains_bpix->Fill(gain);
	   _TH1F_Pedestals_bpix->Fill(ped);
	   nchannelsB++;
	   if(layer==1)      { GainsL1->Fill(gain); PedsL1->Fill(ped);ADCL1->Fill(adc);VcalL1->Fill(elec);}
	   else if(layer==2) { GainsL2->Fill(gain); PedsL2->Fill(ped);ADCL2->Fill(adc);VcalL2->Fill(elec);}
	   else if(layer==3) { GainsL3->Fill(gain); PedsL3->Fill(ped);ADCL3->Fill(adc);VcalL3->Fill(elec);}
	   else if(layer==4) { GainsL4->Fill(gain); PedsL4->Fill(ped);ADCL4->Fill(adc);VcalL4->Fill(elec);}
	   else {cout<<" what else? "<<layer<<endl;}

	 } else {  // FPix
	   nchannelsF++;
	   _TH1F_Gains_fpix->Fill(gain); _TH1F_Pedestals_fpix->Fill(ped); ADCF->Fill(adc); VcalF->Fill(elec);
	 }
	 
	 if(PRINT && select && (col_iter%10==0) && (row_iter%20==0) ) 
	   std::cout <<" DetId "<<detid<<" "<<name<<" Col "<<col_iter<<" Row "<<row_iter
		     <<" Ped "<<ped<<" Gain "<<gain<<std::endl;	 
	 //if( (detid==303054856 && col_iter==109 && row_iter==49 ) || 
	 //  (detid==303046684 && col_iter==150 && row_iter==65 ) ||
	 //  (detid==303042580 && col_iter==94 && row_iter==137 ) ) 
	 //std::cout <<" DetId "<<detid<<" "<<name<<" Col "<<col_iter<<" Row "<<row_iter
	 // <<" Ped "<<ped<<" Gain "<<gain<<std::endl;	 
	} // row
      } // up/down 
      } // col

    if( ndeadInModule>0 && ndeadInModule<66460) { // model not completely missing/empty
      cout<<" Dead pixels in module "<<detid<<" "<<ndeadInModule<<endl;
      if     (layer==1) ndeadInL1 += ndeadInModule;
      else if(layer==2) ndeadInL2 += ndeadInModule;
      else if(layer==3) ndeadInL3 += ndeadInModule;
      else if(layer==4) ndeadInL4 += ndeadInModule;
      else              ndeadInF += ndeadInModule;
    }

    _deadfrac_m[detid]/=nchannelspermod;
    _noisyfrac_m[detid]/=nchannelspermod;
    _TH1F_Dead_sum->SetBinContent(ibin,_deadfrac_m[detid]);
    _TH1F_Dead_all->Fill(_deadfrac_m[detid]);
    _TH1F_Noisy_sum->SetBinContent(ibin,_noisyfrac_m[detid]);
    _TH1F_Noisy_all->Fill(_noisyfrac_m[detid]);
    _TH1F_Gains_sum->SetBinContent(ibin,_TH1F_Gains_m[detid]->GetMean());
    _TH1F_Gains_sum->SetBinError(ibin,_TH1F_Gains_m[detid]->GetRMS());
    _TH1F_Pedestals_sum->SetBinContent(ibin,_TH1F_Pedestals_m[detid]->GetMean());
    _TH1F_Pedestals_sum->SetBinError(ibin,_TH1F_Pedestals_m[detid]->GetRMS());

    gainmeanfortree = _TH1F_Gains_m[detid]->GetMean();
    gainrmsfortree  = _TH1F_Gains_m[detid]->GetRMS();
    pedmeanfortree  = _TH1F_Pedestals_m[detid]->GetMean();
    pedrmsfortree   = _TH1F_Pedestals_m[detid]->GetRMS();

    //std::cout<<"DetId "<<detid<<"       GainMean "<<gainmeanfortree<<" RMS "<<gainrmsfortree<<"      PedMean "<<pedmeanfortree<<" RMS "<<pedrmsfortree<<std::endl;
    if(fillTree) tree->Fill();
    ibin++;
  } // modules
    
  std::cout<<" ---> SUMMARY :"<<std::endl;
  std::cout<<"Modules "<<nmodules<<" BPix " <<nmodulesB<<" FPix "<<nmodulesF<< std::endl;
  std::cout<<"Pixels "<<nchannels<<" BPix " <<nchannelsB<<" FPix "<<nchannelsF<< std::endl;
  std::cout<<"Encounted "<<ndead<<" dead pixels "<<nnoisy<<" noisy pixels"<<std::endl;

  std::cout<<"The Gain Mean is "<<_TH1F_Gains_all->GetMean()<<" with rms "<<_TH1F_Gains_all->GetRMS()<<" entries "<<_TH1F_Gains_all->GetEntries()<<std::endl;
  std::cout<<"         in BPIX "<<_TH1F_Gains_bpix->GetMean()<<" with rms "<<_TH1F_Gains_bpix->GetRMS()<<" entries "<<_TH1F_Gains_bpix->GetEntries()<<std::endl;
  std::cout<<"         in FPIX "<<_TH1F_Gains_fpix->GetMean()<<" with rms "<<_TH1F_Gains_fpix->GetRMS()<<" entries "<<_TH1F_Gains_fpix->GetEntries()<<std::endl;
  std::cout<<"The Ped Mean is "<<_TH1F_Pedestals_all->GetMean()<<" with rms "<<_TH1F_Pedestals_all->GetRMS()<<" entries "<<_TH1F_Pedestals_all->GetEntries()<<std::endl;
  std::cout<<"         in BPIX "<<_TH1F_Pedestals_bpix->GetMean()<<" with rms "<<_TH1F_Pedestals_bpix->GetRMS()<<" entries "<<_TH1F_Pedestals_bpix->GetEntries()<<std::endl;
  std::cout<<"         in FPIX "<<_TH1F_Pedestals_fpix->GetMean()<<" with rms "<<_TH1F_Pedestals_fpix->GetRMS()<<" entries "<<_TH1F_Pedestals_fpix->GetEntries()<<std::endl;

  std::cout<<"BPIX1 Gain= "<<GainsL1->GetMean()<<" rms "<<GainsL1->GetRMS()
	   <<" Ped ="<<PedsL1->GetMean()<<" rms "<<PedsL1->GetRMS()<<" entries "<<GainsL1->GetEntries()<<std::endl;
  std::cout<<"BPIX2 Gain= "<<GainsL2->GetMean()<<" rms "<<GainsL2->GetRMS()
	   <<" Ped ="<<PedsL2->GetMean()<<" rms "<<PedsL2->GetRMS()<<" entries "<<GainsL2->GetEntries()<<std::endl;
  std::cout<<"BPIX3 Gain= "<<GainsL3->GetMean()<<" rms "<<GainsL3->GetRMS()
	   <<" Ped ="<<PedsL3->GetMean()<<" rms "<<PedsL3->GetRMS()<<" entries "<<GainsL3->GetEntries()<<std::endl;
  std::cout<<"BPIX4 Gain= "<<GainsL4->GetMean()<<" rms "<<GainsL4->GetRMS()
	   <<" Ped ="<<PedsL4->GetMean()<<" rms "<<PedsL4->GetRMS()<<" entries "<<GainsL4->GetEntries()<<std::endl;
  int tmp = ndeadInL1+ndeadInL2+ndeadInL3+ndeadInL4+ndeadInF;
 cout<<" gains missing in L1 "<<ndeadInL1<<" L2 "<<ndeadInL2<<" L3 "<<ndeadInL3
     <<" L4 "<<ndeadInL4<<" FPix "<<ndeadInF<<" total "<<tmp<<endl;
}

// ------------ method called once each job just before starting event loop  ------------

void SiPixelGainsDBReader::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void SiPixelGainsDBReader::endJob() {}

void SiPixelGainsDBReader::beginRun(const edm::Run&, const edm::EventSetup& iSetup) {}
void SiPixelGainsDBReader::endRun(const edm::Run&, const edm::EventSetup& iSetup) {}

SiPixelGainsDBReader::~SiPixelGainsDBReader() {}



//define this as a plug-in
DEFINE_FWK_MODULE(SiPixelGainsDBReader);
