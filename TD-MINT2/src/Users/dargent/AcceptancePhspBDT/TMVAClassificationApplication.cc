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

void TMVAClassificationApplication(TString myMethod = "BDTG", TString applyTo = "Gen") 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------
   TChain* theTree;
   TString outFileName;

   if(applyTo == "Gen"){
   	theTree = new TChain("MCDecayTree");
   	outFileName = "GenMC_BDT.root";
   	theTree->Add("GenMC.root");
   }
   else if(applyTo == "Sel"){
   	theTree = new TChain("DecayTree");
   	outFileName = "SelMC_BDT.root";
   	theTree->Add("SelMC.root");
   }
   else if(applyTo == "MINT"){
   	theTree = new TChain("DalitzEventList");
   	outFileName = "MintMC_BDT.root";
   	theTree->Add("MintMC.root");
   }

   // Disable not needed branches
   // Ouput tree
   theTree->SetBranchStatus("*",0);
   theTree->SetBranchStatus("*m*",1);
   if(applyTo == "Gen")theTree->SetBranchStatus("*TRUEP*",1);
   if(applyTo == "MINT")theTree->SetBranchStatus("*",1);

   TFile *hFile = new TFile(outFileName,"RECREATE");
   TTree* tree = theTree->CloneTree(0);

   // This loads the library
   TMVA::Tools::Instance();

   // --- Create the Reader object
   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   Float_t r_m_Dpi;
   Float_t r_m_Kspi;
   Float_t r_m_DKs;

   //reader->AddVariable( "TRUE_m_DKs", &r_m_DKs );
   reader->AddVariable( "TRUE_m_Dpi", &r_m_Dpi );
   reader->AddVariable( "TRUE_m_Kspi", &r_m_Kspi );
   
   // --- Book the MVA methods
   TString prefix = "weights/TMVAClassification_all_all";
   reader->BookMVA( myMethod, prefix + "_" + myMethod + ".weights.xml" ); 

   Double_t m_Dpi;
   Double_t m_Kspi;
   Double_t m_DKs;

   theTree->SetBranchAddress( "TRUE_m_Dpi", &m_Dpi );
   theTree->SetBranchAddress( "TRUE_m_Kspi", &m_Kspi );
   theTree->SetBranchAddress( "TRUE_m_DKs", &m_DKs );

   //output file---------------------------------------------------------------------------------------------------------------------------------------
   double BDTG;
   tree->Branch("eff_BDTG",&BDTG, "eff_BDTG/D");

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   int N = theTree->GetEntries();
   for (Long64_t ievt=0; ievt< N ;ievt++) {
       if (ievt%5000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
       theTree->GetEntry(ievt);

       //r_m_DKs= float(m_DKs);
       r_m_Dpi= float(m_Dpi);
       r_m_Kspi= float(m_Kspi);
       	
       TString methodName = myMethod;
       BDTG = double(reader->EvaluateMVA(methodName));
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

void applyToAll(TString myMethod = "BDTG"){
	TMVAClassificationApplication(myMethod, "Gen" );
	TMVAClassificationApplication(myMethod, "Sel" );
	TMVAClassificationApplication(myMethod, "MINT" );
}
