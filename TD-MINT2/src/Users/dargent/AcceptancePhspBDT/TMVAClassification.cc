#include <cstdlib>
#include <iostream>
#include <map>
#include <math.h>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TMVAGui.C"
#include "../../../../../TMVA-v4.2.0/test/variables.C"
#include "../../../../../TMVA-v4.2.0/test/efficiencies.C"
#include "../../../../../TMVA-v4.2.0/test/mvas.C"
#include "../../../../../TMVA-v4.2.0/test/correlations.C"
#include "../../../../../TMVA-v4.2.0/test/mvaeffs.C"
#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#endif

void TMVAClassification( TString myMethodList = "BDTG", TString run = "all", TString Ks = "all")
{
   TChain* background = new TChain("MCDecayTree");
   background->Add("GenMC.root");
  
   TChain* signal = new TChain("DecayTree");
   signal->Add("SelMC.root");

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   // --- Boosted Decision Trees
   Use["BDT"]             = 1; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting 
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.

   //TString outDir = "plots";
   TString outDir = "figs/";
   outDir +=  myMethodList + "_" + run + "_" + Ks;
 
   TString outfileName = "TMVA_Bs2DsKpipi_" +  myMethodList + "_" + run + "_" + Ks + ".root";
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. 
   TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification_" + run + "_" + Ks, outputFile,
                                               "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );

   signal->SetBranchStatus("*",0);  // disable all branches   signal->SetBranchStatus("s_*",1); 
   signal->SetBranchStatus("*m_*",1); 
   signal->SetBranchStatus("*mp_*",1);     
   signal->SetBranchStatus("*cos*",1); 
   signal->SetBranchStatus("KsCat",1); 

   background->SetBranchStatus("*",0);  // disable all branches
   background->SetBranchStatus("*m_*",1);
   background->SetBranchStatus("*mp_*",1);
   background->SetBranchStatus("*cos*",1);
   background->SetBranchStatus("KsCat",1); 

   // Define the input variables that shall be used for the MVA training
   //factory->AddVariable( "TRUE_m_DKs", "m(DKs)", "MeV", 'F' );
   //factory->AddVariable( "TRUE_m_Dpi", "m(D#pi)", "MeV", 'F' );
   //factory->AddVariable( "TRUE_m_Kspi", "m(Ks#pi)", "MeV", 'F' );

   factory->AddVariable( "TRUE_mp_DKs", "m'(DKs)", "MeV", 'F' );
   factory->AddVariable( "TRUE_mp_Dpi", "m'(D#pi)", "MeV", 'F' );
   factory->AddVariable( "TRUE_mp_Kspi", "m'(Ks#pi)", "MeV", 'F' );

   factory->AddVariable( "TRUE_cos_DKs", "#theta'(DKs)", "MeV", 'F' );
   factory->AddVariable( "TRUE_cos_Dpi", "#theta'(D#pi)", "MeV", 'F' );
   factory->AddVariable( "TRUE_cos_Kspi", "#theta'(Ks#pi)", "MeV", 'F' );
    
   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

   // Apply additional cuts on the signal and background samples (can be different)
   TCut mycuts = ""; 
   if(Ks == "LL") mycuts += "KsCat == 0";
   if(Ks == "DD") mycuts += "KsCat == 1";
   //mycuts = "run == " + run.ReplaceAll("run","")  ;
   //mycuts += "TriggerCat == " + trigger.ReplaceAll("t","")  ;
   TCut mycutb = "";
   
   factory->AddSignalTree    ( signal,     signalWeight     );
   factory->AddBackgroundTree( background, backgroundWeight );
   if(Ks == "all")factory->PrepareTrainingAndTestTree( mycuts, mycutb, "nTrain_Signal=1800:nTrain_Background=150000:nTest_Background=50000:SplitMode=Random:NormMode=NumEvents:!V" );
   if(Ks == "LL")factory->PrepareTrainingAndTestTree( mycuts, mycutb, "nTrain_Signal=500:nTrain_Background=150000:nTest_Background=50000:SplitMode=Random:NormMode=NumEvents:!V" );
   if(Ks == "DD")factory->PrepareTrainingAndTestTree( mycuts, mycutb, "nTrain_Signal=1300:nTrain_Background=150000:nTest_Background=50000:SplitMode=Random:NormMode=NumEvents:!V" );

   factory->SetSignalWeightExpression("weight");
   // ---- Book MVA methods

   // Boosted Decision Trees
   if (Use["BDTG"]) // Gradient Boost
	//       factory->BookMethod( TMVA::Types::kBDT, "BDTG",
	//                            "!H:!V:NTrees=400:MinNodeSize=0%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=40:MaxDepth=2:NegWeightTreatment=Pray:nEventsMin=500" );
      factory->BookMethod( TMVA::Types::kBDT, "BDTG",
                           "!H:!V:NTrees=500:MinNodeSize=0.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=40:MaxDepth=2:NegWeightTreatment=Pray" );

   if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=500:MinNodeSize=3.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=80" );

   if (Use["BDTB"]) // Bagging
      factory->BookMethod( TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

   if (Use["BDTD"]) // Decorrelation + Adaptive Boost
      factory->BookMethod( TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

   if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory->BookMethod( TMVA::Types::kBDT, "BDTMitFisher",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

   // Train MVAs using the set of training events
   factory->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory;

   variables(outfileName,"InputVariables_Id", "TMVA Input Variables",kFALSE, kTRUE, outDir);
   correlations( outfileName,  kFALSE, kFALSE, kTRUE ,outDir);
   efficiencies( outfileName,  2, kTRUE ,outDir);
   mvas( outfileName, CompareType,  kTRUE , outDir, true);

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVAGui( outfileName );
}

void trainAll( TString myMethodList = "BDTG", TString trainOn = "MC") {
	gROOT->SetBatch(true);
	TMVAClassification( myMethodList, "all",  "all" );
    TMVAClassification( myMethodList, "all",  "LL" );
    TMVAClassification( myMethodList, "all",  "DD" );
}
