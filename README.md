# PixelDBTools

Various test programs to monitor the content of pixel DB payloads.


## Software setup

Prepare your working directory with CMSSW

```
export SCRAM_ARCH=slc7_amd64_gcc900
cmsrel CMSSW_11_3_0_pre4
cd CMSSW_11_3_0_pre4/src
cmsenv
git clone https://github.com/CMSTrackerDPG/SiPixelTools-PixelDBTools.git SiPixelTools/PixelDBTools
scram b -j 8
cd SiPixelTools/PixelDBTools/test/
```

Move here all DB test files  
1) LorentzAngle Reader and Loader - done

## Instructions on how to create a SiPixelQuality DB

The updated `.cc` and `.h` files have not been pushed to the master CMSSW repository, they are just kept here:  

https://github.com/CMSTrackerDPG/SiPixelTools-PixelDBTools/blob/master/test/SiPixelBadModuleByHandBuilder.cc  
and  
https://github.com/CMSTrackerDPG/SiPixelTools-PixelDBTools/blob/master/test/SiPixelBadModuleByHandBuilder.h

so the first step is to overwrite these files in the CMSSW `CondTools/SiPixel/test/` and then modify

https://github.com/CMSTrackerDPG/SiPixelTools-PixelDBTools/blob/master/test/0SiPixelBadModuleByHandBuilder_cfg_Run305064.py  

line 30 (sqlite filename, make sure that this didnt exists in the folder where you're working)  
line 33 (tagname)  
line 163 (input file list)

and then you have to have file list similar to SiPixelQuality_phase1_Run3Beginning_forDigitizer_fixed.txt which is a list according to the DQM naming convention where there is a whitespace between the ROC and its number.

