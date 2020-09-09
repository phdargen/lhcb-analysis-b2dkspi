# compile code:
(lhcb-proxy-init)
source make.sh

# Run over stripped data: 
Removes not needed branches and applies loose preselection (choose one of the options):

./MiniMaker "B2DKspi_LL\DD" "Data\MC" 11\12\15\16\17\18

or submit all jobs to batch:

source mini_mc.sh

source mini_data.sh


# train BDT   

Requires TMVA 4.2.0:

source setup.sh /work/dargent/TMVA-v4.2.0/

root -l 

.L TMVAClassification.cpp 

TMVAClassification("BDTG","Run1")

# apply BDT

root -l 

.L TMVAClassificationApplication.cpp

TMVAClassificationApplication("Signal/Norm","Data/MC","BDTG")


