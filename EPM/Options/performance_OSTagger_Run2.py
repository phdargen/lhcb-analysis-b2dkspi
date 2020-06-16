#####################################
# FLAVOURTAGGINGTOOLS CONFIGURATION #
#####################################

########
# DATA #
########

#  This is the file/directory that you want to run (LAST ONE IS TAKEN):
#  if a directory is given all root files in it will be read:
  
datafile = "b2dkspi_sw.root"
TupleName = "DecayTree"

###########
# GENERIC #
###########

CalibrationMode = "Bd"
DoCalibrations = 0
CalibrationLink = "MISTAG"
CalibrationDegree = 1
CalibrationModel = "POLY"
UseNewtonRaphson = 0

Selection = "run == 2 && KsCat == 0 && B_PV_TAUERR < 0.15 && B_PV_TAU > 0.2"
#DBGLEVEL    =  5
Nmax        = -1  # Events to run, -1 means all

##################
# SET LHCB STYLE #
##################
PlotLabel = "LHCb"
PlotTitle = 0
PlotExtension = ".eps"
PlotStatBox = 0

###################
# SIMPLEEVALUATOR #
###################

BranchID             = "B_ID"
UseWeight            = 1
BranchWeight         = "N_B_sw_noKs"
UseTau  = 1
TypeTau = "Double_t"
TauUnits = "ps"
BranchTau = "B_PV_TAU"
UseTauErr = 1
TypeTauErr = "Double_t"
BranchTauErr = "B_PV_TAUERR"

#ResolutionGaussian1_A = 0.0076
#ResolutionGaussian1_B = 0.958
ResolutionGaussian1_A = 0.0097
ResolutionGaussian1_B = 0.915
DrawOscillationPlots = 1
#OscillationPlotsMaximum = 1.1

###################
# SPECIFY TAGGERS #
###################

# TaggerName_NumBins = 20

### TAGGERS USED IN Bs->DsK
OS_Muon_Use = 1
OS_Muon_TypeDec          = "Int_t"
OS_Muon_BranchDec        = "B_OS_Muon_TAGDEC"
OS_Muon_TypeProb        = "Double_t"
OS_Muon_BranchProb      = "B_OS_Muon_TAGETA"

OS_Kaon_Use = 1
OS_Kaon_TypeDec        = "Int_t"
OS_Kaon_BranchDec      = "B_OS_Kaon_TAGDEC"
OS_Kaon_TypeProb      = "Double_t"
OS_Kaon_BranchProb    = "B_OS_Kaon_TAGETA"

OS_Electron_Use = 1
OS_Electron_TypeDec      = "Int_t"
OS_Electron_BranchDec    = "B_OS_Electron_TAGDEC"
OS_Electron_TypeProb    = "Double_t"
OS_Electron_BranchProb  = "B_OS_Electron_TAGETA"

#OS_Kaon_Use = 1
#OS_Kaon_TypeDec        = "Short_t"
#OS_Kaon_BranchDec      = "Bs_OS_Kaon_DEC"
#OS_Kaon_TypeProb      = "Float_t"
#OS_Kaon_BranchProb    = "Bs_OS_Kaon_PROB"

VtxCharge_Use = 1
VtxCharge_TypeDec     = "Int_t"
VtxCharge_BranchDec   = "B_VtxCharge_TAGDEC"
VtxCharge_TypeProb   = "Double_t"
VtxCharge_BranchProb = "B_VtxCharge_TAGETA"

#SS_nnetKaon_Use = 1
#SS_nnetKaon_TypeDec      = "Int_t"
#SS_nnetKaon_BranchDec    = "SS_Kaon_DEC"
#SS_nnetKaon_TypeProb    = "Double_t"
#SS_nnetKaon_BranchProb  = "SS_Kaon_PROB"

#SS_Kaon_Use = 1
#SS_Kaon_TypeDec      = "Short_t"
#SS_Kaon_BranchDec    = "Bs_SS_Kaon_DEC"
#SS_Kaon_TypeProb    = "Float_t"
#SS_Kaon_BranchProb  = "Bs_SS_Kaon_PROB"

## BUGGY ?
OS_Charm_Use = 1
OS_Charm_TypeDec = "Int_t"
OS_Charm_BranchDec = "B_OS_Charm_TAGDEC"
OS_Charm_TypeProb = "Double_t"
OS_Charm_BranchProb = "B_OS_Charm_TAGETA"

#OS_Combination_Use  = 1
#OS_Combination_TypeDec	= "Int_t"
#OS_Combination_BranchDec  = "Bs_TAGDECISION_OS"
#OS_Combination_TypeProb	= "Double_t"
#OS_Combination_BranchProb = "Bs_TAGOMEGA_OS"


### OS AND OS+SS COMBINATION

#PerformOfflineCombination_OS = 1
#OS_Muon_InOSComb = 1
#OS_Electron_InOSComb = 1
#OS_nnetKaon_InOSComb = 1
#VtxCharge_InOSComb = 1

#PerformOfflineCombination_OSplusSS = 1
#OS_Combination_InComb = 1
#SS_nnetKaon_InComb = 1

############################
## SAVE CALIBRATION OUTPUT #
############################

#WriteCalibratedMistagBranches = 1
#OS_Combination_Write = 1
#OS_Muon_Write = 1
#OS_Electron_Write = 1
#OS_nnetKaon_Write = 1
#VtxCharge_Write = 1
##SS_nnetKaon_Write = 1
##Combination_Write = 1

#CalibratedOutputFile = "OS_combo_Run2.root"


############################
# CALIBRATION INPUT VALUES #
############################

#import EspressoCalibrations_OS_Run1.py
#SaveCalibrationsToXML = 1

OS_Muon_CalibrationArchive = "out_OS_Run2/OS_Muon_Calibration.xml"
OS_Electron_CalibrationArchive = "out_OS_Run2/OS_Electron_Calibration.xml"
OS_Kaon_CalibrationArchive = "out_OS_Run2/OS_Kaon_Calibration.xml"
VtxCharge_CalibrationArchive = "out_OS_Run2/VtxCharge_Calibration.xml"
OS_Charm_CalibrationArchive = "out_OS_Run2/OS_Charm_Calibration.xml"

