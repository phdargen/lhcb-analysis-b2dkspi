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
   theTree->SetBranchStatus("*cos*",1);
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

   Float_t r_mp_Dpi;
   Float_t r_mp_Kspi;
   Float_t r_mp_DKs;
   Float_t r_cos_Dpi;
   Float_t r_cos_Kspi;
   Float_t r_cos_DKs;

   //reader->AddVariable( "TRUE_m_DKs", &r_m_DKs );
   //reader->AddVariable( "TRUE_m_Dpi", &r_m_Dpi );
   //reader->AddVariable( "TRUE_m_Kspi", &r_m_Kspi );

   reader->AddVariable( "TRUE_mp_DKs", &r_mp_DKs );
   reader->AddVariable( "TRUE_mp_Dpi", &r_mp_Dpi );
   reader->AddVariable( "TRUE_mp_Kspi", &r_mp_Kspi );

   reader->AddVariable( "TRUE_cos_DKs", &r_cos_DKs );
   reader->AddVariable( "TRUE_cos_Dpi", &r_cos_Dpi );
   reader->AddVariable( "TRUE_cos_Kspi", &r_cos_Kspi );
    
   // --- Book the MVA methods
   //TString prefix = "weights/TMVAClassification_all_all";
   TString prefix = "weights/TMVAClassification_";
   std::vector<TString> weightFiles;
   weightFiles.push_back("all_all");
   weightFiles.push_back("all_LL");
   weightFiles.push_back("all_DD");
   for(int i= 0 ; i < weightFiles.size(); i++) 
        reader->BookMVA( myMethod + weightFiles[i], prefix + weightFiles[i] + "_" + myMethod + ".weights.xml" ); 

   //reader->BookMVA( myMethod, prefix + "_" + myMethod + ".weights.xml" ); 

   Double_t m_Dpi;
   Double_t m_Kspi;
   Double_t m_DKs;
   Double_t mp_Dpi;
   Double_t mp_Kspi;
   Double_t mp_DKs;
   Double_t cos_Dpi;
   Double_t cos_Kspi;
   Double_t cos_DKs;
    
   theTree->SetBranchAddress( "TRUE_m_Dpi", &m_Dpi );
   theTree->SetBranchAddress( "TRUE_m_Kspi", &m_Kspi );
   theTree->SetBranchAddress( "TRUE_m_DKs", &m_DKs );

   theTree->SetBranchAddress( "TRUE_mp_Dpi", &mp_Dpi );
   theTree->SetBranchAddress( "TRUE_mp_Kspi", &mp_Kspi );
   theTree->SetBranchAddress( "TRUE_mp_DKs", &mp_DKs );
 
   theTree->SetBranchAddress( "TRUE_cos_Dpi", &cos_Dpi );
   theTree->SetBranchAddress( "TRUE_cos_Kspi", &cos_Kspi );
   theTree->SetBranchAddress( "TRUE_cos_DKs", &cos_DKs );
    
   //output file---------------------------------------------------------------------------------------------------------------------------------------
   double BDTG,BDTG_DD,BDTG_LL;
   tree->Branch("eff_BDTG",&BDTG, "eff_BDTG/D");
   tree->Branch("eff_BDTG_LL",&BDTG_LL, "eff_BDTG_LL/D");
   tree->Branch("eff_BDTG_DD",&BDTG_DD, "eff_BDTG_DD/D");

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   int N = theTree->GetEntries();
   for (Long64_t ievt=0; ievt< N ;ievt++) {
       if (ievt%5000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;
       theTree->GetEntry(ievt);

       r_m_DKs= float(m_DKs);
       r_m_Dpi= float(m_Dpi);
       r_m_Kspi= float(m_Kspi);

       r_mp_DKs= float(mp_DKs);
       r_mp_Dpi= float(mp_Dpi);
       r_mp_Kspi= float(mp_Kspi);

       r_cos_DKs= float(cos_DKs);
       r_cos_Dpi= float(cos_Dpi);
       r_cos_Kspi= float(cos_Kspi);

       //TString methodName = myMethod;
       BDTG = double(reader->EvaluateMVA(myMethod + weightFiles[0]));
       BDTG_LL = double(reader->EvaluateMVA(myMethod + weightFiles[1]));
       BDTG_DD = double(reader->EvaluateMVA(myMethod + weightFiles[2]));
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
