/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <math.h>

#include "TMath.h"

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include <TChain.h>
#include "TStopwatch.h"

//#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace std;
using namespace TMVA;

void TMVAClassificationApplication(TString decay = "Signal", TString dataType = "Data", TString myMethod = "BDTG", TString trainedOn = "MC" ) 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------
   TChain* theTree = new TChain("DecayTree");


   TString outFileName = "BDT/";

   if(decay == "Signal" && dataType == "Data"){ 	  
	theTree->Add("Preselected/Data_b2dkspi_LL_11.root");
    theTree->Add("Preselected/Data_b2dkspi_LL_12.root");
    theTree->Add("Preselected/Data_b2dkspi_LL_15.root");
    theTree->Add("Preselected/Data_b2dkspi_LL_16.root");
    theTree->Add("Preselected/Data_b2dkspi_LL_17.root");
    theTree->Add("Preselected/Data_b2dkspi_LL_18.root");

    theTree->Add("Preselected/Data_b2dkspi_DD_11.root");
    theTree->Add("Preselected/Data_b2dkspi_DD_12.root");
    theTree->Add("Preselected/Data_b2dkspi_DD_15.root");
    theTree->Add("Preselected/Data_b2dkspi_DD_16.root");
    theTree->Add("Preselected/Data_b2dkspi_DD_17.root");
    theTree->Add("Preselected/Data_b2dkspi_DD_18.root");
       
    outFileName += "signal_data.root";
   }

   else if(decay == "Signal" && dataType == "MC"){ 	  
	theTree->Add("Preselected/MC_b2dkspi_DD_12.root");
	outFileName += "signal_mc.root";
   }
   else {
	cout << "Unknown options, I'll crash now." << endl;
	throw "ERROR";
   }

   // Ouput tree
   TFile *hFile = new TFile(outFileName,"RECREATE");
   TTree* tree = theTree->CloneTree(0);

   // This loads the library
   TMVA::Tools::Instance();

   // --- Create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   Float_t r_log_B_FDCHI2_OWNPV;
   Float_t r_log_B_IPCHI2_OWNPV;
   Float_t r_log_B_DIRA_OWNPV;
   
   Float_t r_log_Ks_PT;
   Float_t r_log_Ks_RFD;
   Float_t r_log_D_RFD; 
    
   Float_t r_PV_CHI2NDOF;
   Float_t r_log_pi_ProbNNpi;
   Float_t r_log_K_D_ProbNNk;
    
   Float_t r_log_min_IPCHI2; 
   Float_t r_maxCos;
   
   reader->AddVariable( "log_B_FDCHI2_OWNPV := log(B_FDCHI2_OWNPV)",&r_log_B_FDCHI2_OWNPV );
   reader->AddVariable( "log_B_IPCHI2_OWNPV := log(B_IPCHI2_OWNPV)",&r_log_B_IPCHI2_OWNPV );
   reader->AddVariable( "log_B_DIRA := log(1-B_DIRA_OWNPV)",&r_log_B_DIRA_OWNPV );

   reader->AddVariable( "log_Ks_PT:=log(Ks_PT)",&r_log_Ks_PT);
   reader->AddVariable( "log_Ks_RFD:=log(Ks_RFD)",&r_log_Ks_RFD);
   reader->AddVariable( "log_D_RFD:=log(D_RFD)",&r_log_D_RFD);
    
   reader->AddVariable( "PV_CHI2NDOF", &r_PV_CHI2NDOF );
   reader->AddVariable( "log_pi_ProbNNpi := log(1-pi_ProbNNpi)",&r_log_pi_ProbNNpi);
   reader->AddVariable( "log_K_D_ProbNNk := log(1-K_D_ProbNNk)",&r_log_K_D_ProbNNk);

   reader->AddVariable( "log_min_IPCHI2 := log(track_min_IPCHI2)",&r_log_min_IPCHI2);
   reader->AddVariable( "maxCos2", &r_maxCos );

   //reader->AddVariable( "Bs_ptasy_1.00",&r_Bs_ptasy );
   //reader->AddVariable("max_ghostProb",&r_max_ghostProb);

   // --- Book the MVA methods
   TString prefix = "weights/TMVAClassification_"+trainedOn+ "_";

   std::vector<TString> weightFiles;
   //weightFiles.push_back("all_all_all");
   weightFiles.push_back("run1_LL_all");
   weightFiles.push_back("run2_LL_all");
   weightFiles.push_back("run1_DD_all");
   weightFiles.push_back("run2_DD_all");

   for(int i= 0 ; i < weightFiles.size(); i++) 
        reader->BookMVA( myMethod + weightFiles[i], prefix + weightFiles[i] + "_" + myMethod + ".weights.xml" ); 

    Double_t B_FDCHI2_OWNPV;
    Double_t B_IPCHI2_OWNPV;
    Double_t B_DIRA_OWNPV;
    
    Double_t Ks_PT;
    Double_t Ks_RFD;
    Double_t D_RFD; 
    
    Double_t PV_CHI2NDOF;
    Double_t pi_ProbNNpi;
    Double_t K_D_ProbNNk;
    
    Double_t min_IPCHI2; 
    Double_t maxCos;

    theTree->SetBranchAddress( "B_FDCHI2_OWNPV", &B_FDCHI2_OWNPV );
    theTree->SetBranchAddress( "B_IPCHI2_OWNPV", &B_IPCHI2_OWNPV );
    theTree->SetBranchAddress( "B_DIRA_OWNPV", &B_DIRA_OWNPV );

    theTree->SetBranchAddress( "Ks_PT", &Ks_PT );
    theTree->SetBranchAddress( "Ks_RFD", &Ks_RFD );
    theTree->SetBranchAddress( "D_RFD", &D_RFD );
    
    theTree->SetBranchAddress( "PV_CHI2NDOF", &PV_CHI2NDOF );
    theTree->SetBranchAddress( "pi_ProbNNpi", &pi_ProbNNpi );
    theTree->SetBranchAddress( "K_D_ProbNNk", &K_D_ProbNNk );

    theTree->SetBranchAddress( "track_min_IPCHI2", &min_IPCHI2 );
    theTree->SetBranchAddress( "maxCos2", &maxCos );
   
   Int_t year, run, Ds_finalState, TriggerCat, KsCat; 
   ULong64_t eventNumber;
   theTree->SetBranchAddress( "year", &year );
   theTree->SetBranchAddress( "run", &run );
   theTree->SetBranchAddress( "TriggerCat", &TriggerCat );
   theTree->SetBranchAddress( "eventNumber", &eventNumber );
   theTree->SetBranchAddress( "KsCat", &KsCat );
    
   Double_t FullDTF_status;
   Double_t DTF_status;
   Double_t PV_status;
   Double_t B_DTF_MM;
   theTree->SetBranchAddress( "FullDTF_status", &FullDTF_status );
   theTree->SetBranchAddress( "DTF_status", &DTF_status );
   theTree->SetBranchAddress( "PV_status", &PV_status );
   theTree->SetBranchAddress( "B_DTF_MM", &B_DTF_MM );

   //output file---------------------------------------------------------------------------------------------------------------------------------------
   Float_t BDTG_response;
   double BDTG;
   tree->Branch("BDTG_response",&BDTG_response, "BDTG_response/F");
   tree->Branch("BDTG",&BDTG, "BDTG/D");

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;

   TStopwatch sw;
   sw.Start();
   int N = theTree->GetEntries();
   for (Long64_t ievt=0; ievt< N ;ievt++) {

      if (ievt%5000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

        theTree->GetEntry(ievt);
        if(FullDTF_status > 1)continue;
        if(DTF_status > 1)continue;
        if(PV_status > 1)continue;
        if(TMath::IsNaN(B_DTF_MM))continue;

        r_log_B_FDCHI2_OWNPV = float(log(B_FDCHI2_OWNPV));
        r_log_B_IPCHI2_OWNPV = float(log(B_IPCHI2_OWNPV));
        r_log_B_DIRA_OWNPV = float(log(1.-B_DIRA_OWNPV));

        r_log_Ks_PT = float(log(Ks_PT));
        r_log_Ks_RFD = float(log(Ks_RFD));
        r_log_D_RFD = float(log(D_RFD));
       
        r_PV_CHI2NDOF = float(PV_CHI2NDOF);
        r_log_pi_ProbNNpi = float(log(1-pi_ProbNNpi));
        r_log_K_D_ProbNNk = float(log(1-K_D_ProbNNk));
       
        r_log_min_IPCHI2 = float(log(min_IPCHI2)); 
        r_maxCos = float(maxCos);

        TString methodName = myMethod + "run";
        methodName += run;
        //methodName += "_t"; 
        //methodName += TriggerCat;
        //if(eventNumber % 2 == 0) methodName += "_odd" ;
        //else methodName += "_even";
        if(KsCat == 0) methodName += "_LL" ;
        else methodName += "_DD";
        methodName += "_all" ;
       
        //TString methodName = myMethod + "all_all_all";
        BDTG_response=reader->EvaluateMVA(methodName);
        BDTG = double(BDTG_response);
        tree->Fill();    
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   tree->Write();
   delete reader;
   hFile->Close();
    
   cout << "Wrote to file: " << outFileName << endl;
   cout << "==> TMVAClassificationApplication is done!" << endl << endl;
} 

void applyToAll(TString myMethod = "BDTG", TString trainedOn = "MC" ){

    TMVAClassificationApplication("Signal", "Data", myMethod, trainedOn );
    TMVAClassificationApplication("Signal", "MC", myMethod, trainedOn );    
}
