//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Jul 24 11:58:28 2018 by ROOT version 5.34/10
// from TTree MinuitParameterSetNtp/MinuitParameterSetNtp
// found on file: signal_toy/pull_par0_1.root
//////////////////////////////////////////////////////////

#ifndef pull_h
#define pull_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class pull {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Double_t        C_mean;
   Double_t        C_init;
   Double_t        C_err;
   Double_t        C_pull;
   Double_t        D_mean;
   Double_t        D_init;
   Double_t        D_err;
   Double_t        D_pull;
   Double_t        D_bar_mean;
   Double_t        D_bar_init;
   Double_t        D_bar_err;
   Double_t        D_bar_pull;
   Double_t        S_mean;
   Double_t        S_init;
   Double_t        S_err;
   Double_t        S_pull;
   Double_t        S_bar_mean;
   Double_t        S_bar_init;
   Double_t        S_bar_err;
   Double_t        S_bar_pull;

   Double_t	   delta_pull_C;
   Double_t        delta_pull_D;
   Double_t        delta_pull_S;
   Double_t        delta_pull_D_bar;
   Double_t        delta_pull_S_bar;

   Double_t	   delta_pull_C0;
   Double_t        delta_pull_D0;
   Double_t        delta_pull_S0;
   Double_t        delta_pull_D_bar0;
   Double_t        delta_pull_S_bar0;
   Double_t        delta_pull_C1;
   Double_t        delta_pull_D1;
   Double_t        delta_pull_S1;
   Double_t        delta_pull_D_bar1;
   Double_t        delta_pull_S_bar1;
   Double_t        delta_pull_C2;
   Double_t        delta_pull_D2;
   Double_t        delta_pull_S2;
   Double_t        delta_pull_D_bar2;
   Double_t        delta_pull_S_bar2;
   Double_t        delta_pull_C3;
   Double_t        delta_pull_D3;
   Double_t        delta_pull_S3;
   Double_t        delta_pull_D_bar3;
   Double_t        delta_pull_S_bar3;

   // List of branches
   TBranch        *b_C_mean;   //!
   TBranch        *b_C_init;   //!
   TBranch        *b_C_err;   //!
   TBranch        *b_C_pull;   //!
   TBranch        *b_D_mean;   //!
   TBranch        *b_D_init;   //!
   TBranch        *b_D_err;   //!
   TBranch        *b_D_pull;   //!
   TBranch        *b_D_bar_mean;   //!
   TBranch        *b_D_bar_init;   //!
   TBranch        *b_D_bar_err;   //!
   TBranch        *b_D_bar_pull;   //!
   TBranch        *b_S_mean;   //!
   TBranch        *b_S_init;   //!
   TBranch        *b_S_err;   //!
   TBranch        *b_S_pull;   //!
   TBranch        *b_S_bar_mean;   //!
   TBranch        *b_S_bar_init;   //!
   TBranch        *b_S_bar_err;   //!
   TBranch        *b_S_bar_pull;   //!
   TBranch        *b___noname0;   //!

   TBranch	  *b_delta_pull_C;
   TBranch        *b_delta_pull_D;
   TBranch        *b_delta_pull_S;
   TBranch        *b_delta_pull_D_bar;
   TBranch        *b_delta_pull_S_bar;

   TBranch	  *b_delta_pull_C0;
   TBranch        *b_delta_pull_D0;
   TBranch        *b_delta_pull_S0;
   TBranch        *b_delta_pull_D_bar0;
   TBranch        *b_delta_pull_S_bar0;
   TBranch        *b_delta_pull_C1;
   TBranch        *b_delta_pull_D1;
   TBranch        *b_delta_pull_S1;
   TBranch        *b_delta_pull_D_bar1;
   TBranch        *b_delta_pull_S_bar1;
   TBranch        *b_delta_pull_C2;
   TBranch        *b_delta_pull_D2;
   TBranch        *b_delta_pull_S2;
   TBranch        *b_delta_pull_D_bar2;
   TBranch        *b_delta_pull_S_bar2;
   TBranch        *b_delta_pull_C3;
   TBranch        *b_delta_pull_D3;
   TBranch        *b_delta_pull_S3;
   TBranch        *b_delta_pull_D_bar3;
   TBranch        *b_delta_pull_S_bar3;

   pull();
   virtual ~pull();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(string parName);
   virtual void     LoopSyst_noChol(string parName);
   virtual void     getShift(string parName);
   virtual void     makeDeltaPulls_Col();
   virtual void     makeDeltaPulls_noCol();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef pull_cxx

pull::pull() : fChain(0)
{


}	


pull::pull_DsK_noChol() : fChain(0) 
{

   TChain* chain =  new TChain("MinuitParameterSetNtp");

   chain->Add("signal_toy_noChol/pull_*.root"); //1

   tree= (TTree*)  chain;
   Init(tree);

   Loop("noChol");
}

pull::pull_Chol() : fChain(0) 
{
    makeDeltaPulls_Col();

   TChain* chain =  new TChain("MinuitParameterSetNtp");

   chain->Add("tdfit_deltaPulls_Chol.root"); //1

   tree= (TTree*)  chain;
   Init(tree);

   LoopSyst("systChol");
  // getShift("all");
}


pull::pull_noChol() : fChain(0) 
{
    makeDeltaPulls_noCol();

   TChain* chain =  new TChain("MinuitParameterSetNtp");

   chain->Add("tdfit_deltaPulls_noChol.root"); //1

   tree= (TTree*)  chain;
   Init(tree);

   LoopSyst_noChol("systNoChol");
  // getShift("all");
}

pull::pull_noSyst() : fChain(0) 
{
   TChain* chain =  new TChain("MinuitParameterSetNtp");

   chain->Add("/auto/data/kecke/BsDsKpipi/toys/signal_noSyst/pull_*.root"); //1

   tree= (TTree*)  chain;
   Init(tree);

   Loop("all");
  // getShift("all");
}

pull::pull_c0() : fChain(0) 
{
   TChain* chain =  new TChain("MinuitParameterSetNtp");
   chain->Add("signal_toy_fitterValid/pull_par0_*.root"); //1

   tree= (TTree*)  chain;
   Init(tree);

   Loop("c0");
//   getShift("c0");
}

pull::pull_c1() : fChain(0) 
{
   TChain* chain =  new TChain("MinuitParameterSetNtp");
   chain->Add("signal_toy_fitterValid/pull_par1_*.root"); //1

   tree= (TTree*)  chain;
   Init(tree);

   Loop("c1");
 //  getShift("c1");
}

pull::pull_c2() : fChain(0) 
{
   TChain* chain =  new TChain("MinuitParameterSetNtp");
   chain->Add("signal_toy_fitterValid/pull_par2_*.root"); //1

   tree= (TTree*)  chain;
   Init(tree);

   Loop("c2");
 //  getShift("c2");
}

pull::pull_c3() : fChain(0) 
{
   TChain* chain =  new TChain("MinuitParameterSetNtp");
   chain->Add("signal_toy_fitterValid/pull_par3_*.root"); //1

   tree= (TTree*)  chain;
   Init(tree);

   Loop("c3");
 //  getShift("c3");
}

pull::~pull()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t pull::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t pull::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void pull::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("C_mean", &C_mean, &b_C_mean);
   fChain->SetBranchAddress("C_init", &C_init, &b_C_init);
   fChain->SetBranchAddress("C_err", &C_err, &b_C_err);
   fChain->SetBranchAddress("C_pull", &C_pull, &b_C_pull);
   fChain->SetBranchAddress("D_mean", &D_mean, &b_D_mean);
   fChain->SetBranchAddress("D_init", &D_init, &b_D_init);
   fChain->SetBranchAddress("D_err", &D_err, &b_D_err);
   fChain->SetBranchAddress("D_pull", &D_pull, &b_D_pull);
   fChain->SetBranchAddress("D_bar_mean", &D_bar_mean, &b_D_bar_mean);
   fChain->SetBranchAddress("D_bar_init", &D_bar_init, &b_D_bar_init);
   fChain->SetBranchAddress("D_bar_err", &D_bar_err, &b_D_bar_err);
   fChain->SetBranchAddress("D_bar_pull", &D_bar_pull, &b_D_bar_pull);
   fChain->SetBranchAddress("S_mean", &S_mean, &b_S_mean);
   fChain->SetBranchAddress("S_init", &S_init, &b_S_init);
   fChain->SetBranchAddress("S_err", &S_err, &b_S_err);
   fChain->SetBranchAddress("S_pull", &S_pull, &b_S_pull);
   fChain->SetBranchAddress("S_bar_mean", &S_bar_mean, &b_S_bar_mean);
   fChain->SetBranchAddress("S_bar_init", &S_bar_init, &b_S_bar_init);
   fChain->SetBranchAddress("S_bar_err", &S_bar_err, &b_S_bar_err);
   fChain->SetBranchAddress("S_bar_pull", &S_bar_pull, &b_S_bar_pull);

   fChain->SetBranchAddress("delta_pull_C", &delta_pull_C, &b_delta_pull_C);
   fChain->SetBranchAddress("delta_pull_D", &delta_pull_D, &b_delta_pull_D);
   fChain->SetBranchAddress("delta_pull_S", &delta_pull_S, &b_delta_pull_S);
   fChain->SetBranchAddress("delta_pull_D_bar", &delta_pull_D_bar, &b_delta_pull_D_bar);
   fChain->SetBranchAddress("delta_pull_S_bar", &delta_pull_S_bar, &b_delta_pull_S_bar);

   fChain->SetBranchAddress("delta_pull_C0", &delta_pull_C0, &b_delta_pull_C0);
   fChain->SetBranchAddress("delta_pull_D0", &delta_pull_D0, &b_delta_pull_D0);
   fChain->SetBranchAddress("delta_pull_S0", &delta_pull_S0, &b_delta_pull_S0);
   fChain->SetBranchAddress("delta_pull_D_bar0", &delta_pull_D_bar0, &b_delta_pull_D_bar0);
   fChain->SetBranchAddress("delta_pull_S_bar0", &delta_pull_S_bar0, &b_delta_pull_S_bar0);
   fChain->SetBranchAddress("delta_pull_C1", &delta_pull_C1, &b_delta_pull_C1);
   fChain->SetBranchAddress("delta_pull_D1", &delta_pull_D1, &b_delta_pull_D1);
   fChain->SetBranchAddress("delta_pull_S1", &delta_pull_S1, &b_delta_pull_S1);
   fChain->SetBranchAddress("delta_pull_D_bar1", &delta_pull_D_bar1, &b_delta_pull_D_bar1);
   fChain->SetBranchAddress("delta_pull_S_bar1", &delta_pull_S_bar1, &b_delta_pull_S_bar1);
   fChain->SetBranchAddress("delta_pull_C2", &delta_pull_C2, &b_delta_pull_C2);
   fChain->SetBranchAddress("delta_pull_D2", &delta_pull_D2, &b_delta_pull_D2);
   fChain->SetBranchAddress("delta_pull_S2", &delta_pull_S2, &b_delta_pull_S2);
   fChain->SetBranchAddress("delta_pull_D_bar2", &delta_pull_D_bar2, &b_delta_pull_D_bar2);
   fChain->SetBranchAddress("delta_pull_S_bar2", &delta_pull_S_bar2, &b_delta_pull_S_bar2);
   fChain->SetBranchAddress("delta_pull_C3", &delta_pull_C3, &b_delta_pull_C3);
   fChain->SetBranchAddress("delta_pull_D3", &delta_pull_D3, &b_delta_pull_D3);
   fChain->SetBranchAddress("delta_pull_S3", &delta_pull_S3, &b_delta_pull_S3);
   fChain->SetBranchAddress("delta_pull_D_bar3", &delta_pull_D_bar3, &b_delta_pull_D_bar3);
   fChain->SetBranchAddress("delta_pull_S_bar3", &delta_pull_S_bar3, &b_delta_pull_S_bar3);


   Notify();
}

Bool_t pull::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void pull::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t pull::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef pull_cxx
