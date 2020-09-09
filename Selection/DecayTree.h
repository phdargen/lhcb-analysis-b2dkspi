//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Feb 28 10:42:57 2020 by ROOT version 5.34/36
// from TTree DecayTree/DecayTree
// found on file: b2dkspi.root
//////////////////////////////////////////////////////////

#ifndef DecayTree_h
#define DecayTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TLorentzVector.h"

using namespace std;

struct Decay{
    enum  Type { B2DKspi_LL, B2DKspi_DD, B2DKsK_LL, B2DKsK_DD };
};

struct Year{
    enum  Type { y11 = 11, y12 = 12, y15 = 15, y16 = 16, y17 = 17, y18 = 18 };
};

struct DataType{
    enum  Type { mc, data };
};

struct McEventType{
    enum  Type { BdDKspi, BsDKsK, BsDstKsK, BdDstKspi, BsDstKspi };
};

static const double massKs = 497.611;
static const double massKaon = 493.68;
static const double massPion = 139.57;
static const double massProton = 938.27;
static const double massMuon = 195.65838;
static const double massPhi = 1019.46;
static const double massKstar = 895.81;
static const double massDs = 1968.30;
static const double massDminus = 1869.61;
static const double massLambda_c = 2286.46;
static const double c_light = 0.299792458;

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
   const Int_t kMaxB_ENDVERTEX_COV = 1;
   const Int_t kMaxB_OWNPV_COV = 1;
   const Int_t kMaxB_TOPPV_COV = 1;
   const Int_t kMaxD_ENDVERTEX_COV = 1;
   const Int_t kMaxD_OWNPV_COV = 1;
   const Int_t kMaxD_TOPPV_COV = 1;
   const Int_t kMaxD_ORIVX_COV = 1;
   const Int_t kMaxK_D_OWNPV_COV = 1;
   const Int_t kMaxK_D_TOPPV_COV = 1;
   const Int_t kMaxK_D_ORIVX_COV = 1;
   const Int_t kMaxpi1_D_OWNPV_COV = 1;
   const Int_t kMaxpi1_D_TOPPV_COV = 1;
   const Int_t kMaxpi1_D_ORIVX_COV = 1;
   const Int_t kMaxpi2_D_OWNPV_COV = 1;
   const Int_t kMaxpi2_D_TOPPV_COV = 1;
   const Int_t kMaxpi2_D_ORIVX_COV = 1;
   const Int_t kMaxKs_ENDVERTEX_COV = 1;
   const Int_t kMaxKs_OWNPV_COV = 1;
   const Int_t kMaxKs_TOPPV_COV = 1;
   const Int_t kMaxKs_ORIVX_COV = 1;
   const Int_t kMaxpip_Ks_OWNPV_COV = 1;
   const Int_t kMaxpip_Ks_TOPPV_COV = 1;
   const Int_t kMaxpip_Ks_ORIVX_COV = 1;
   const Int_t kMaxpim_Ks_OWNPV_COV = 1;
   const Int_t kMaxpim_Ks_TOPPV_COV = 1;
   const Int_t kMaxpim_Ks_ORIVX_COV = 1;
   const Int_t kMaxpi_OWNPV_COV = 1;
   const Int_t kMaxpi_TOPPV_COV = 1;
   const Int_t kMaxpi_ORIVX_COV = 1;

class DecayTree {
public :
    TString _inFileLoc;
    TString _outFileLoc;
    TString _outFileName;
    TTree          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain
    
    Decay::Type _decay;
    Year::Type _year;
    DataType::Type _data;
    McEventType::Type _mcEventType;
    TString _polarity;
    
    DecayTree(Decay::Type decay, Year::Type year, DataType::Type dataType, TString polarity = "Both", TString inFileLoc = "/auto/data/dargent/B2psiKpipi/", TString outFileLoc = "/auto/data/dargent/B2psiKpipi/", McEventType::Type mcEventType = McEventType::BdDKspi );
    
    virtual ~DecayTree();
    virtual Int_t    GetEntry(Long64_t entry);
    virtual Long64_t LoadTree(Long64_t entry);
    virtual void     Init();
    virtual void     Loop();
    virtual Bool_t   Notify();
    virtual void     Show(Long64_t entry = -1);
    virtual TTree* GetInputTree();
    
    inline  Bool_t TriggerCuts(Long64_t entry);
    inline  Bool_t LooseCuts(Long64_t entry);
    inline  Bool_t LooseCutsLTU(Long64_t entry);
    inline void set_LorentzVectors();
    
   TLorentzVector pi,pip_Ks,pim_Ks,pi1_D,pi2_D,K_D,Ks,D;
   TLorentzVector FullDTF_pi,FullDTF_pip_Ks,FullDTF_pim_Ks,FullDTF_pi1_D,FullDTF_pi2_D,FullDTF_K_D,FullDTF_Ks,FullDTF_D;
   TLorentzVector FullBsDTF_pi,FullBsDTF_pip_Ks,FullBsDTF_pim_Ks,FullBsDTF_pi1_D,FullBsDTF_pi2_D,FullBsDTF_K_D,FullBsDTF_Ks,FullBsDTF_D;
    
   TLorentzVector K_fromD_asP_MissID,K_fromD_asPi_MissID;
   TLorentzVector pi_asP_MissID,pi_asK_MissID;
   TLorentzVector pi1_fromD_asP_MissID,pi1_fromD_asK_MissID;
   TLorentzVector pi2_fromD_asP_MissID,pi2_fromD_asK_MissID;
    
   // Declaration of leaf types
   Double_t        B_ETA;
   Double_t        B_MINIP;
   Double_t        B_MINIPCHI2;
   Double_t        B_MINIPNEXTBEST;
   Double_t        B_MINIPCHI2NEXTBEST;
   Int_t           nPV;
   Float_t         B_AllIP[100];   //[nPV]
   Float_t         B_AllIPchi2[100];   //[nPV]
   Float_t         B_AllDIRA[100];   //[nPV]
   Double_t        B_ENDVERTEX_X;
   Double_t        B_ENDVERTEX_Y;
   Double_t        B_ENDVERTEX_Z;
   Double_t        B_ENDVERTEX_XERR;
   Double_t        B_ENDVERTEX_YERR;
   Double_t        B_ENDVERTEX_ZERR;
   Double_t        B_ENDVERTEX_CHI2;
   Int_t           B_ENDVERTEX_NDOF;
   Float_t         B_ENDVERTEX_COV_[3][3];
   Double_t        B_OWNPV_X;
   Double_t        B_OWNPV_Y;
   Double_t        B_OWNPV_Z;
   Double_t        B_OWNPV_XERR;
   Double_t        B_OWNPV_YERR;
   Double_t        B_OWNPV_ZERR;
   Double_t        B_OWNPV_CHI2;
   Int_t           B_OWNPV_NDOF;
   Float_t         B_OWNPV_COV_[3][3];
   Double_t        B_IP_OWNPV;
   Double_t        B_IPCHI2_OWNPV;
   Double_t        B_FD_OWNPV;
   Double_t        B_FDCHI2_OWNPV;
   Double_t        B_DIRA_OWNPV;
   Double_t        B_TOPPV_X;
   Double_t        B_TOPPV_Y;
   Double_t        B_TOPPV_Z;
   Double_t        B_TOPPV_XERR;
   Double_t        B_TOPPV_YERR;
   Double_t        B_TOPPV_ZERR;
   Double_t        B_TOPPV_CHI2;
   Int_t           B_TOPPV_NDOF;
   Float_t         B_TOPPV_COV_[3][3];
   Double_t        B_IP_TOPPV;
   Double_t        B_IPCHI2_TOPPV;
   Double_t        B_FD_TOPPV;
   Double_t        B_FDCHI2_TOPPV;
   Double_t        B_DIRA_TOPPV;
   Double_t        B_P;
   Double_t        B_PT;
   Double_t        B_PE;
   Double_t        B_PX;
   Double_t        B_PY;
   Double_t        B_PZ;
   Double_t        B_MM;
   Double_t        B_MMERR;
   Double_t        B_M;
   Int_t           B_ID;
   Int_t           B_TAGDECISION;
   Double_t        B_TAGOMEGA;
   Int_t           B_TAGDECISION_OS;
   Double_t        B_TAGOMEGA_OS;
   Int_t           B_TAGGER;
   Short_t         B_OS_Muon_DEC;
   Float_t         B_OS_Muon_PROB;
   Short_t         B_OS_Electron_DEC;
   Float_t         B_OS_Electron_PROB;
   Short_t         B_OS_Kaon_DEC;
   Float_t         B_OS_Kaon_PROB;
   Short_t         B_SS_Kaon_DEC;
   Float_t         B_SS_Kaon_PROB;
   Short_t         B_SS_Pion_DEC;
   Float_t         B_SS_Pion_PROB;
   Short_t         B_SS_PionBDT_DEC;
   Float_t         B_SS_PionBDT_PROB;
   Short_t         B_VtxCharge_DEC;
   Float_t         B_VtxCharge_PROB;
   Short_t         B_OS_nnetKaon_DEC;
   Float_t         B_OS_nnetKaon_PROB;
   Short_t         B_SS_nnetKaon_DEC;
   Float_t         B_SS_nnetKaon_PROB;
   Short_t         B_SS_Proton_DEC;
   Float_t         B_SS_Proton_PROB;
   Short_t         B_OS_Charm_DEC;
   Float_t         B_OS_Charm_PROB;
   Double_t        B_cpx_1_00;
   Double_t        B_cpy_1_00;
   Double_t        B_cpz_1_00;
   Double_t        B_cpt_1_00;
   Double_t        B_cp_1_00;
   Int_t           B_cmult_1_00;
   Double_t        B_pxasy_1_00;
   Double_t        B_pyasy_1_00;
   Double_t        B_pzasy_1_00;
   Double_t        B_pasy_1_00;
   Double_t        B_ptasy_1_00;
   Double_t        B_DOCA1;
   Double_t        B_TAU;
   Double_t        B_TAUERR;
   Int_t           B_BDTF_nPV;
   Float_t         B_BDTF_Dplus_M[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_Dplus_MERR[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_Dplus_P[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_Dplus_PERR[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_Dplus_ctau[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_Dplus_ctauErr[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_Dplus_decayLength[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_Dplus_decayLengthErr[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_KS0_M[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_KS0_MERR[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_KS0_P[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_KS0_PERR[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_KS0_ctau[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_KS0_ctauErr[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_KS0_decayLength[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_KS0_decayLengthErr[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_M[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_MERR[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_P[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_PERR[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_PV_X[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_PV_Y[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_PV_Z[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_PV_key[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_chi2[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_ctau[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_ctauErr[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_decayLength[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_decayLengthErr[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_nDOF[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_nIter[100];   //[B_BDTF_nPV]
   Float_t         B_BDTF_status[100];   //[B_BDTF_nPV]
   Int_t           B_DTF_nPV;
   Float_t         B_DTF_Dplus_M[100];   //[B_DTF_nPV]
   Float_t         B_DTF_Dplus_MERR[100];   //[B_DTF_nPV]
   Float_t         B_DTF_Dplus_P[100];   //[B_DTF_nPV]
   Float_t         B_DTF_Dplus_PERR[100];   //[B_DTF_nPV]
   Float_t         B_DTF_Dplus_ctau[100];   //[B_DTF_nPV]
   Float_t         B_DTF_Dplus_ctauErr[100];   //[B_DTF_nPV]
   Float_t         B_DTF_Dplus_decayLength[100];   //[B_DTF_nPV]
   Float_t         B_DTF_Dplus_decayLengthErr[100];   //[B_DTF_nPV]
   Float_t         B_DTF_KS0_M[100];   //[B_DTF_nPV]
   Float_t         B_DTF_KS0_MERR[100];   //[B_DTF_nPV]
   Float_t         B_DTF_KS0_P[100];   //[B_DTF_nPV]
   Float_t         B_DTF_KS0_PERR[100];   //[B_DTF_nPV]
   Float_t         B_DTF_KS0_ctau[100];   //[B_DTF_nPV]
   Float_t         B_DTF_KS0_ctauErr[100];   //[B_DTF_nPV]
   Float_t         B_DTF_KS0_decayLength[100];   //[B_DTF_nPV]
   Float_t         B_DTF_KS0_decayLengthErr[100];   //[B_DTF_nPV]
   Float_t         B_DTF_M[100];   //[B_DTF_nPV]
   Float_t         B_DTF_MERR[100];   //[B_DTF_nPV]
   Float_t         B_DTF_P[100];   //[B_DTF_nPV]
   Float_t         B_DTF_PERR[100];   //[B_DTF_nPV]
   Float_t         B_DTF_PV_X[100];   //[B_DTF_nPV]
   Float_t         B_DTF_PV_Y[100];   //[B_DTF_nPV]
   Float_t         B_DTF_PV_Z[100];   //[B_DTF_nPV]
   Float_t         B_DTF_PV_key[100];   //[B_DTF_nPV]
   Float_t         B_DTF_chi2[100];   //[B_DTF_nPV]
   Float_t         B_DTF_ctau[100];   //[B_DTF_nPV]
   Float_t         B_DTF_ctauErr[100];   //[B_DTF_nPV]
   Float_t         B_DTF_decayLength[100];   //[B_DTF_nPV]
   Float_t         B_DTF_decayLengthErr[100];   //[B_DTF_nPV]
   Float_t         B_DTF_nDOF[100];   //[B_DTF_nPV]
   Float_t         B_DTF_nIter[100];   //[B_DTF_nPV]
   Float_t         B_DTF_status[100];   //[B_DTF_nPV]

   Int_t           B_FullDTF_nPV;
   Float_t         B_FullDTF_Dplus_Kplus_ID[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_Kplus_PE[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_Kplus_PX[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_Kplus_PY[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_Kplus_PZ[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_M[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_MERR[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_P[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_PERR[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_ctau[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_ctauErr[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_decayLength[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_decayLengthErr[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_piplus_0_ID[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_piplus_0_PE[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_piplus_0_PX[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_piplus_0_PY[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_piplus_0_PZ[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_piplus_ID[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_piplus_PE[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_piplus_PX[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_piplus_PY[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_Dplus_piplus_PZ[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_M[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_MERR[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_P[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_PERR[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_ctau[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_ctauErr[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_decayLength[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_decayLengthErr[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_piplus_0_ID[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_piplus_0_PE[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_piplus_0_PX[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_piplus_0_PY[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_piplus_0_PZ[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_piplus_ID[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_piplus_PE[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_piplus_PX[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_piplus_PY[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_KS0_piplus_PZ[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_M[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_MERR[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_P[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_PERR[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_PV_X[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_PV_Y[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_PV_Z[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_PV_key[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_chi2[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_ctau[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_ctauErr[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_decayLength[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_decayLengthErr[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_nDOF[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_nIter[100];   //[B_FullDTF_nPV]
   Float_t         B_FullDTF_status[100];   //[B_FullDTF_nPV]
    Float_t         B_FullDTF_piplus_ID[100];   //[B_FullDTF_nPV]
    Float_t         B_FullDTF_piplus_PE[100];   //[B_FullDTF_nPV]
    Float_t         B_FullDTF_piplus_PX[100];   //[B_FullDTF_nPV]
    Float_t         B_FullDTF_piplus_PY[100];   //[B_FullDTF_nPV]
    Float_t         B_FullDTF_piplus_PZ[100];   //[B_FullDTF_nPV]
    
    Int_t           B_FullBsDTF_nPV;
    Float_t         B_FullBsDTF_Dplus_Kplus_ID[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_Kplus_PE[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_Kplus_PX[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_Kplus_PY[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_Kplus_PZ[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_M[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_MERR[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_P[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_PERR[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_ctau[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_ctauErr[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_decayLength[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_decayLengthErr[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_piplus_0_ID[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_piplus_0_PE[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_piplus_0_PX[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_piplus_0_PY[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_piplus_0_PZ[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_piplus_ID[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_piplus_PE[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_piplus_PX[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_piplus_PY[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_Dplus_piplus_PZ[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_M[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_MERR[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_P[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_PERR[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_ctau[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_ctauErr[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_decayLength[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_decayLengthErr[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_piplus_0_ID[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_piplus_0_PE[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_piplus_0_PX[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_piplus_0_PY[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_piplus_0_PZ[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_piplus_ID[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_piplus_PE[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_piplus_PX[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_piplus_PY[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_KS0_piplus_PZ[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_M[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_MERR[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_P[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_PERR[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_PV_X[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_PV_Y[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_PV_Z[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_PV_key[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_chi2[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_ctau[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_ctauErr[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_decayLength[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_decayLengthErr[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_nDOF[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_nIter[100];   //[B_FullBsDTF_nPV]
    Float_t         B_FullBsDTF_piplus_ID[100];   //[B_FullDTF_nPV]
    Float_t         B_FullBsDTF_piplus_PE[100];   //[B_FullDTF_nPV]
    Float_t         B_FullBsDTF_piplus_PX[100];   //[B_FullDTF_nPV]
    Float_t         B_FullBsDTF_piplus_PY[100];   //[B_FullDTF_nPV]
    Float_t         B_FullBsDTF_piplus_PZ[100];   //[B_FullDTF_nPV]
    Float_t         B_FullBsDTF_status[100];   //[B_FullBsDTF_nPV]
     
   Int_t           B_PV_nPV;
   Float_t         B_PV_Dplus_Kplus_ID[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_Kplus_PE[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_Kplus_PX[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_Kplus_PY[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_Kplus_PZ[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_M[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_MERR[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_P[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_PERR[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_ctau[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_ctauErr[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_decayLength[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_decayLengthErr[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_piplus_0_ID[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_piplus_0_PE[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_piplus_0_PX[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_piplus_0_PY[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_piplus_0_PZ[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_piplus_ID[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_piplus_PE[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_piplus_PX[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_piplus_PY[100];   //[B_PV_nPV]
   Float_t         B_PV_Dplus_piplus_PZ[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_M[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_MERR[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_P[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_PERR[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_ctau[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_ctauErr[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_decayLength[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_decayLengthErr[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_piplus_0_ID[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_piplus_0_PE[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_piplus_0_PX[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_piplus_0_PY[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_piplus_0_PZ[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_piplus_ID[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_piplus_PE[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_piplus_PX[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_piplus_PY[100];   //[B_PV_nPV]
   Float_t         B_PV_KS0_piplus_PZ[100];   //[B_PV_nPV]
   Float_t         B_PV_M[100];   //[B_PV_nPV]
   Float_t         B_PV_MERR[100];   //[B_PV_nPV]
   Float_t         B_PV_P[100];   //[B_PV_nPV]
   Float_t         B_PV_PERR[100];   //[B_PV_nPV]
   Float_t         B_PV_PV_X[100];   //[B_PV_nPV]
   Float_t         B_PV_PV_Y[100];   //[B_PV_nPV]
   Float_t         B_PV_PV_Z[100];   //[B_PV_nPV]
   Float_t         B_PV_PV_key[100];   //[B_PV_nPV]
   Float_t         B_PV_chi2[100];   //[B_PV_nPV]
   Float_t         B_PV_ctau[100];   //[B_PV_nPV]
   Float_t         B_PV_ctauErr[100];   //[B_PV_nPV]
   Float_t         B_PV_decayLength[100];   //[B_PV_nPV]
   Float_t         B_PV_decayLengthErr[100];   //[B_PV_nPV]
   Float_t         B_PV_nDOF[100];   //[B_PV_nPV]
   Float_t         B_PV_nIter[100];   //[B_PV_nPV]
   Float_t         B_PV_piplus_ID[100];   //[B_PV_nPV]
   Float_t         B_PV_piplus_PE[100];   //[B_PV_nPV]
   Float_t         B_PV_piplus_PX[100];   //[B_PV_nPV]
   Float_t         B_PV_piplus_PY[100];   //[B_PV_nPV]
   Float_t         B_PV_piplus_PZ[100];   //[B_PV_nPV]
   Float_t         B_PV_status[100];   //[B_PV_nPV]
   Bool_t          B_L0Global_Dec;
   Bool_t          B_L0Global_TIS;
   Bool_t          B_L0Global_TOS;
   Bool_t          B_Hlt1Global_Dec;
   Bool_t          B_Hlt1Global_TIS;
   Bool_t          B_Hlt1Global_TOS;
   Bool_t          B_Hlt1Phys_Dec;
   Bool_t          B_Hlt1Phys_TIS;
   Bool_t          B_Hlt1Phys_TOS;
   Bool_t          B_Hlt2Global_Dec;
   Bool_t          B_Hlt2Global_TIS;
   Bool_t          B_Hlt2Global_TOS;
   Bool_t          B_Hlt2Phys_Dec;
   Bool_t          B_Hlt2Phys_TIS;
   Bool_t          B_Hlt2Phys_TOS;
    
    Bool_t          B_Hlt1TrackAllL0Decision_Dec;
    Bool_t          B_Hlt1TrackAllL0Decision_TIS;
    Bool_t          B_Hlt1TrackAllL0Decision_TOS;

    Bool_t          B_Hlt2Topo2BodyBBDTDecision_Dec;
    Bool_t          B_Hlt2Topo2BodyBBDTDecision_TIS;
    Bool_t          B_Hlt2Topo2BodyBBDTDecision_TOS;
    Bool_t          B_Hlt2Topo3BodyBBDTDecision_Dec;
    Bool_t          B_Hlt2Topo3BodyBBDTDecision_TIS;
    Bool_t          B_Hlt2Topo3BodyBBDTDecision_TOS;
    Bool_t          B_Hlt2Topo4BodyBBDTDecision_Dec;
    Bool_t          B_Hlt2Topo4BodyBBDTDecision_TIS;
    Bool_t          B_Hlt2Topo4BodyBBDTDecision_TOS;
    
   Bool_t          B_L0HadronDecision_Dec;
   Bool_t          B_L0HadronDecision_TIS;
   Bool_t          B_L0HadronDecision_TOS;
   Bool_t          B_L0MuonDecision_Dec;
   Bool_t          B_L0MuonDecision_TIS;
   Bool_t          B_L0MuonDecision_TOS;
   Bool_t          B_L0DiMuonDecision_Dec;
   Bool_t          B_L0DiMuonDecision_TIS;
   Bool_t          B_L0DiMuonDecision_TOS;
   Bool_t          B_L0ElectronDecision_Dec;
   Bool_t          B_L0ElectronDecision_TIS;
   Bool_t          B_L0ElectronDecision_TOS;
   Bool_t          B_L0PhotonDecision_Dec;
   Bool_t          B_L0PhotonDecision_TIS;
   Bool_t          B_L0PhotonDecision_TOS;
   Bool_t          B_Hlt1TrackMVADecision_Dec;
   Bool_t          B_Hlt1TrackMVADecision_TIS;
   Bool_t          B_Hlt1TrackMVADecision_TOS;
   Bool_t          B_Hlt1TwoTrackMVADecision_Dec;
   Bool_t          B_Hlt1TwoTrackMVADecision_TIS;
   Bool_t          B_Hlt1TwoTrackMVADecision_TOS;
   Bool_t          B_Hlt2IncPhiDecision_Dec;
   Bool_t          B_Hlt2IncPhiDecision_TIS;
   Bool_t          B_Hlt2IncPhiDecision_TOS;
   Bool_t          B_Hlt2PhiIncPhiDecision_Dec;
   Bool_t          B_Hlt2PhiIncPhiDecision_TIS;
   Bool_t          B_Hlt2PhiIncPhiDecision_TOS;
   Bool_t          B_Hlt2Topo2BodyDecision_Dec;
   Bool_t          B_Hlt2Topo2BodyDecision_TIS;
   Bool_t          B_Hlt2Topo2BodyDecision_TOS;
   Bool_t          B_Hlt2Topo3BodyDecision_Dec;
   Bool_t          B_Hlt2Topo3BodyDecision_TIS;
   Bool_t          B_Hlt2Topo3BodyDecision_TOS;
   Bool_t          B_Hlt2Topo4BodyDecision_Dec;
   Bool_t          B_Hlt2Topo4BodyDecision_TIS;
   Bool_t          B_Hlt2Topo4BodyDecision_TOS;
   Double_t        D_ETA;
   Double_t        D_MINIP;
   Double_t        D_MINIPCHI2;
   Double_t        D_MINIPNEXTBEST;
   Double_t        D_MINIPCHI2NEXTBEST;
   Float_t         D_AllIP[100];   //[nPV]
   Float_t         D_AllIPchi2[100];   //[nPV]
   Float_t         D_AllDIRA[100];   //[nPV]
   Double_t        D_ENDVERTEX_X;
   Double_t        D_ENDVERTEX_Y;
   Double_t        D_ENDVERTEX_Z;
   Double_t        D_ENDVERTEX_XERR;
   Double_t        D_ENDVERTEX_YERR;
   Double_t        D_ENDVERTEX_ZERR;
   Double_t        D_ENDVERTEX_CHI2;
   Int_t           D_ENDVERTEX_NDOF;
   Float_t         D_ENDVERTEX_COV_[3][3];
   Double_t        D_OWNPV_X;
   Double_t        D_OWNPV_Y;
   Double_t        D_OWNPV_Z;
   Double_t        D_OWNPV_XERR;
   Double_t        D_OWNPV_YERR;
   Double_t        D_OWNPV_ZERR;
   Double_t        D_OWNPV_CHI2;
   Int_t           D_OWNPV_NDOF;
   Float_t         D_OWNPV_COV_[3][3];
   Double_t        D_IP_OWNPV;
   Double_t        D_IPCHI2_OWNPV;
   Double_t        D_FD_OWNPV;
   Double_t        D_FDCHI2_OWNPV;
   Double_t        D_DIRA_OWNPV;
   Double_t        D_TOPPV_X;
   Double_t        D_TOPPV_Y;
   Double_t        D_TOPPV_Z;
   Double_t        D_TOPPV_XERR;
   Double_t        D_TOPPV_YERR;
   Double_t        D_TOPPV_ZERR;
   Double_t        D_TOPPV_CHI2;
   Int_t           D_TOPPV_NDOF;
   Float_t         D_TOPPV_COV_[3][3];
   Double_t        D_IP_TOPPV;
   Double_t        D_IPCHI2_TOPPV;
   Double_t        D_FD_TOPPV;
   Double_t        D_FDCHI2_TOPPV;
   Double_t        D_DIRA_TOPPV;
   Double_t        D_ORIVX_X;
   Double_t        D_ORIVX_Y;
   Double_t        D_ORIVX_Z;
   Double_t        D_ORIVX_XERR;
   Double_t        D_ORIVX_YERR;
   Double_t        D_ORIVX_ZERR;
   Double_t        D_ORIVX_CHI2;
   Int_t           D_ORIVX_NDOF;
   Float_t         D_ORIVX_COV_[3][3];
   Double_t        D_IP_ORIVX;
   Double_t        D_IPCHI2_ORIVX;
   Double_t        D_FD_ORIVX;
   Double_t        D_FDCHI2_ORIVX;
   Double_t        D_DIRA_ORIVX;
   Double_t        D_P;
   Double_t        D_PT;
   Double_t        D_PE;
   Double_t        D_PX;
   Double_t        D_PY;
   Double_t        D_PZ;
   Double_t        D_MM;
   Double_t        D_MMERR;
   Double_t        D_M;
   Int_t           D_ID;
   Double_t        D_cpx_1_00;
   Double_t        D_cpy_1_00;
   Double_t        D_cpz_1_00;
   Double_t        D_cpt_1_00;
   Double_t        D_cp_1_00;
   Int_t           D_cmult_1_00;
   Double_t        D_pxasy_1_00;
   Double_t        D_pyasy_1_00;
   Double_t        D_pzasy_1_00;
   Double_t        D_pasy_1_00;
   Double_t        D_ptasy_1_00;
   Double_t        D_DOCA1;
   Double_t        D_TAU;
   Double_t        D_TAUERR;
   Double_t        K_D_ETA;
   Double_t        K_D_MC12TuneV2_ProbNNe;
   Double_t        K_D_MC12TuneV2_ProbNNmu;
   Double_t        K_D_MC12TuneV2_ProbNNpi;
   Double_t        K_D_MC12TuneV2_ProbNNk;
   Double_t        K_D_MC12TuneV2_ProbNNp;
   Double_t        K_D_MC12TuneV2_ProbNNghost;
   Double_t        K_D_MC12TuneV3_ProbNNe;
   Double_t        K_D_MC12TuneV3_ProbNNmu;
   Double_t        K_D_MC12TuneV3_ProbNNpi;
   Double_t        K_D_MC12TuneV3_ProbNNk;
   Double_t        K_D_MC12TuneV3_ProbNNp;
   Double_t        K_D_MC12TuneV3_ProbNNghost;
   Double_t        K_D_MC12TuneV4_ProbNNe;
   Double_t        K_D_MC12TuneV4_ProbNNmu;
   Double_t        K_D_MC12TuneV4_ProbNNpi;
   Double_t        K_D_MC12TuneV4_ProbNNk;
   Double_t        K_D_MC12TuneV4_ProbNNp;
   Double_t        K_D_MC12TuneV4_ProbNNghost;
   Double_t        K_D_MC15TuneV1_ProbNNe;
   Double_t        K_D_MC15TuneV1_ProbNNmu;
   Double_t        K_D_MC15TuneV1_ProbNNpi;
   Double_t        K_D_MC15TuneV1_ProbNNk;
   Double_t        K_D_MC15TuneV1_ProbNNp;
   Double_t        K_D_MC15TuneV1_ProbNNghost;
   Double_t        K_D_MINIP;
   Double_t        K_D_MINIPCHI2;
   Double_t        K_D_MINIPNEXTBEST;
   Double_t        K_D_MINIPCHI2NEXTBEST;
   Float_t         K_D_AllIP[100];   //[nPV]
   Float_t         K_D_AllIPchi2[100];   //[nPV]
   Double_t        K_D_OWNPV_X;
   Double_t        K_D_OWNPV_Y;
   Double_t        K_D_OWNPV_Z;
   Double_t        K_D_OWNPV_XERR;
   Double_t        K_D_OWNPV_YERR;
   Double_t        K_D_OWNPV_ZERR;
   Double_t        K_D_OWNPV_CHI2;
   Int_t           K_D_OWNPV_NDOF;
   Float_t         K_D_OWNPV_COV_[3][3];
   Double_t        K_D_IP_OWNPV;
   Double_t        K_D_IPCHI2_OWNPV;
   Double_t        K_D_TOPPV_X;
   Double_t        K_D_TOPPV_Y;
   Double_t        K_D_TOPPV_Z;
   Double_t        K_D_TOPPV_XERR;
   Double_t        K_D_TOPPV_YERR;
   Double_t        K_D_TOPPV_ZERR;
   Double_t        K_D_TOPPV_CHI2;
   Int_t           K_D_TOPPV_NDOF;
   Float_t         K_D_TOPPV_COV_[3][3];
   Double_t        K_D_IP_TOPPV;
   Double_t        K_D_IPCHI2_TOPPV;
   Double_t        K_D_ORIVX_X;
   Double_t        K_D_ORIVX_Y;
   Double_t        K_D_ORIVX_Z;
   Double_t        K_D_ORIVX_XERR;
   Double_t        K_D_ORIVX_YERR;
   Double_t        K_D_ORIVX_ZERR;
   Double_t        K_D_ORIVX_CHI2;
   Int_t           K_D_ORIVX_NDOF;
   Float_t         K_D_ORIVX_COV_[3][3];
   Double_t        K_D_IP_ORIVX;
   Double_t        K_D_IPCHI2_ORIVX;
   Double_t        K_D_P;
   Double_t        K_D_PT;
   Double_t        K_D_PE;
   Double_t        K_D_PX;
   Double_t        K_D_PY;
   Double_t        K_D_PZ;
   Double_t        K_D_M;
   Int_t           K_D_ID;
   Double_t        K_D_PIDe;
   Double_t        K_D_PIDmu;
   Double_t        K_D_PIDK;
   Double_t        K_D_PIDp;
   Double_t        K_D_ProbNNe;
   Double_t        K_D_ProbNNk;
   Double_t        K_D_ProbNNp;
   Double_t        K_D_ProbNNpi;
   Double_t        K_D_ProbNNmu;
   Double_t        K_D_ProbNNghost;
   Bool_t          K_D_hasMuon;
   Bool_t          K_D_isMuon;
   Bool_t          K_D_hasRich;
   Bool_t          K_D_UsedRichAerogel;
   Bool_t          K_D_UsedRich1Gas;
   Bool_t          K_D_UsedRich2Gas;
   Bool_t          K_D_RichAboveElThres;
   Bool_t          K_D_RichAboveMuThres;
   Bool_t          K_D_RichAbovePiThres;
   Bool_t          K_D_RichAboveKaThres;
   Bool_t          K_D_RichAbovePrThres;
   Bool_t          K_D_hasCalo;
   Int_t           K_D_TRACK_Type;
   Int_t           K_D_TRACK_Key;
   Double_t        K_D_TRACK_CHI2NDOF;
   Double_t        K_D_TRACK_PCHI2;
   Double_t        K_D_TRACK_MatchCHI2;
   Double_t        K_D_TRACK_GhostProb;
   Double_t        K_D_TRACK_CloneDist;
   Double_t        K_D_TRACK_Likelihood;
   Double_t        K_D_cpx_1_00;
   Double_t        K_D_cpy_1_00;
   Double_t        K_D_cpz_1_00;
   Double_t        K_D_cpt_1_00;
   Double_t        K_D_cp_1_00;
   Int_t           K_D_cmult_1_00;
   Double_t        K_D_pxasy_1_00;
   Double_t        K_D_pyasy_1_00;
   Double_t        K_D_pzasy_1_00;
   Double_t        K_D_pasy_1_00;
   Double_t        K_D_ptasy_1_00;
   Double_t        pi1_D_ETA;
   Double_t        pi1_D_MC12TuneV2_ProbNNe;
   Double_t        pi1_D_MC12TuneV2_ProbNNmu;
   Double_t        pi1_D_MC12TuneV2_ProbNNpi;
   Double_t        pi1_D_MC12TuneV2_ProbNNk;
   Double_t        pi1_D_MC12TuneV2_ProbNNp;
   Double_t        pi1_D_MC12TuneV2_ProbNNghost;
   Double_t        pi1_D_MC12TuneV3_ProbNNe;
   Double_t        pi1_D_MC12TuneV3_ProbNNmu;
   Double_t        pi1_D_MC12TuneV3_ProbNNpi;
   Double_t        pi1_D_MC12TuneV3_ProbNNk;
   Double_t        pi1_D_MC12TuneV3_ProbNNp;
   Double_t        pi1_D_MC12TuneV3_ProbNNghost;
   Double_t        pi1_D_MC12TuneV4_ProbNNe;
   Double_t        pi1_D_MC12TuneV4_ProbNNmu;
   Double_t        pi1_D_MC12TuneV4_ProbNNpi;
   Double_t        pi1_D_MC12TuneV4_ProbNNk;
   Double_t        pi1_D_MC12TuneV4_ProbNNp;
   Double_t        pi1_D_MC12TuneV4_ProbNNghost;
   Double_t        pi1_D_MC15TuneV1_ProbNNe;
   Double_t        pi1_D_MC15TuneV1_ProbNNmu;
   Double_t        pi1_D_MC15TuneV1_ProbNNpi;
   Double_t        pi1_D_MC15TuneV1_ProbNNk;
   Double_t        pi1_D_MC15TuneV1_ProbNNp;
   Double_t        pi1_D_MC15TuneV1_ProbNNghost;
   Double_t        pi1_D_MINIP;
   Double_t        pi1_D_MINIPCHI2;
   Double_t        pi1_D_MINIPNEXTBEST;
   Double_t        pi1_D_MINIPCHI2NEXTBEST;
   Float_t         pi1_D_AllIP[100];   //[nPV]
   Float_t         pi1_D_AllIPchi2[100];   //[nPV]
   Double_t        pi1_D_OWNPV_X;
   Double_t        pi1_D_OWNPV_Y;
   Double_t        pi1_D_OWNPV_Z;
   Double_t        pi1_D_OWNPV_XERR;
   Double_t        pi1_D_OWNPV_YERR;
   Double_t        pi1_D_OWNPV_ZERR;
   Double_t        pi1_D_OWNPV_CHI2;
   Int_t           pi1_D_OWNPV_NDOF;
   Float_t         pi1_D_OWNPV_COV_[3][3];
   Double_t        pi1_D_IP_OWNPV;
   Double_t        pi1_D_IPCHI2_OWNPV;
   Double_t        pi1_D_TOPPV_X;
   Double_t        pi1_D_TOPPV_Y;
   Double_t        pi1_D_TOPPV_Z;
   Double_t        pi1_D_TOPPV_XERR;
   Double_t        pi1_D_TOPPV_YERR;
   Double_t        pi1_D_TOPPV_ZERR;
   Double_t        pi1_D_TOPPV_CHI2;
   Int_t           pi1_D_TOPPV_NDOF;
   Float_t         pi1_D_TOPPV_COV_[3][3];
   Double_t        pi1_D_IP_TOPPV;
   Double_t        pi1_D_IPCHI2_TOPPV;
   Double_t        pi1_D_ORIVX_X;
   Double_t        pi1_D_ORIVX_Y;
   Double_t        pi1_D_ORIVX_Z;
   Double_t        pi1_D_ORIVX_XERR;
   Double_t        pi1_D_ORIVX_YERR;
   Double_t        pi1_D_ORIVX_ZERR;
   Double_t        pi1_D_ORIVX_CHI2;
   Int_t           pi1_D_ORIVX_NDOF;
   Float_t         pi1_D_ORIVX_COV_[3][3];
   Double_t        pi1_D_IP_ORIVX;
   Double_t        pi1_D_IPCHI2_ORIVX;
   Double_t        pi1_D_P;
   Double_t        pi1_D_PT;
   Double_t        pi1_D_PE;
   Double_t        pi1_D_PX;
   Double_t        pi1_D_PY;
   Double_t        pi1_D_PZ;
   Double_t        pi1_D_M;
   Int_t           pi1_D_ID;
   Double_t        pi1_D_PIDe;
   Double_t        pi1_D_PIDmu;
   Double_t        pi1_D_PIDK;
   Double_t        pi1_D_PIDp;
   Double_t        pi1_D_ProbNNe;
   Double_t        pi1_D_ProbNNk;
   Double_t        pi1_D_ProbNNp;
   Double_t        pi1_D_ProbNNpi;
   Double_t        pi1_D_ProbNNmu;
   Double_t        pi1_D_ProbNNghost;
   Bool_t          pi1_D_hasMuon;
   Bool_t          pi1_D_isMuon;
   Bool_t          pi1_D_hasRich;
   Bool_t          pi1_D_UsedRichAerogel;
   Bool_t          pi1_D_UsedRich1Gas;
   Bool_t          pi1_D_UsedRich2Gas;
   Bool_t          pi1_D_RichAboveElThres;
   Bool_t          pi1_D_RichAboveMuThres;
   Bool_t          pi1_D_RichAbovePiThres;
   Bool_t          pi1_D_RichAboveKaThres;
   Bool_t          pi1_D_RichAbovePrThres;
   Bool_t          pi1_D_hasCalo;
   Int_t           pi1_D_TRACK_Type;
   Int_t           pi1_D_TRACK_Key;
   Double_t        pi1_D_TRACK_CHI2NDOF;
   Double_t        pi1_D_TRACK_PCHI2;
   Double_t        pi1_D_TRACK_MatchCHI2;
   Double_t        pi1_D_TRACK_GhostProb;
   Double_t        pi1_D_TRACK_CloneDist;
   Double_t        pi1_D_TRACK_Likelihood;
   Double_t        pi1_D_cpx_1_00;
   Double_t        pi1_D_cpy_1_00;
   Double_t        pi1_D_cpz_1_00;
   Double_t        pi1_D_cpt_1_00;
   Double_t        pi1_D_cp_1_00;
   Int_t           pi1_D_cmult_1_00;
   Double_t        pi1_D_pxasy_1_00;
   Double_t        pi1_D_pyasy_1_00;
   Double_t        pi1_D_pzasy_1_00;
   Double_t        pi1_D_pasy_1_00;
   Double_t        pi1_D_ptasy_1_00;
   Double_t        pi2_D_ETA;
   Double_t        pi2_D_MC12TuneV2_ProbNNe;
   Double_t        pi2_D_MC12TuneV2_ProbNNmu;
   Double_t        pi2_D_MC12TuneV2_ProbNNpi;
   Double_t        pi2_D_MC12TuneV2_ProbNNk;
   Double_t        pi2_D_MC12TuneV2_ProbNNp;
   Double_t        pi2_D_MC12TuneV2_ProbNNghost;
   Double_t        pi2_D_MC12TuneV3_ProbNNe;
   Double_t        pi2_D_MC12TuneV3_ProbNNmu;
   Double_t        pi2_D_MC12TuneV3_ProbNNpi;
   Double_t        pi2_D_MC12TuneV3_ProbNNk;
   Double_t        pi2_D_MC12TuneV3_ProbNNp;
   Double_t        pi2_D_MC12TuneV3_ProbNNghost;
   Double_t        pi2_D_MC12TuneV4_ProbNNe;
   Double_t        pi2_D_MC12TuneV4_ProbNNmu;
   Double_t        pi2_D_MC12TuneV4_ProbNNpi;
   Double_t        pi2_D_MC12TuneV4_ProbNNk;
   Double_t        pi2_D_MC12TuneV4_ProbNNp;
   Double_t        pi2_D_MC12TuneV4_ProbNNghost;
   Double_t        pi2_D_MC15TuneV1_ProbNNe;
   Double_t        pi2_D_MC15TuneV1_ProbNNmu;
   Double_t        pi2_D_MC15TuneV1_ProbNNpi;
   Double_t        pi2_D_MC15TuneV1_ProbNNk;
   Double_t        pi2_D_MC15TuneV1_ProbNNp;
   Double_t        pi2_D_MC15TuneV1_ProbNNghost;
   Double_t        pi2_D_MINIP;
   Double_t        pi2_D_MINIPCHI2;
   Double_t        pi2_D_MINIPNEXTBEST;
   Double_t        pi2_D_MINIPCHI2NEXTBEST;
   Float_t         pi2_D_AllIP[100];   //[nPV]
   Float_t         pi2_D_AllIPchi2[100];   //[nPV]
   Double_t        pi2_D_OWNPV_X;
   Double_t        pi2_D_OWNPV_Y;
   Double_t        pi2_D_OWNPV_Z;
   Double_t        pi2_D_OWNPV_XERR;
   Double_t        pi2_D_OWNPV_YERR;
   Double_t        pi2_D_OWNPV_ZERR;
   Double_t        pi2_D_OWNPV_CHI2;
   Int_t           pi2_D_OWNPV_NDOF;
   Float_t         pi2_D_OWNPV_COV_[3][3];
   Double_t        pi2_D_IP_OWNPV;
   Double_t        pi2_D_IPCHI2_OWNPV;
   Double_t        pi2_D_TOPPV_X;
   Double_t        pi2_D_TOPPV_Y;
   Double_t        pi2_D_TOPPV_Z;
   Double_t        pi2_D_TOPPV_XERR;
   Double_t        pi2_D_TOPPV_YERR;
   Double_t        pi2_D_TOPPV_ZERR;
   Double_t        pi2_D_TOPPV_CHI2;
   Int_t           pi2_D_TOPPV_NDOF;
   Float_t         pi2_D_TOPPV_COV_[3][3];
   Double_t        pi2_D_IP_TOPPV;
   Double_t        pi2_D_IPCHI2_TOPPV;
   Double_t        pi2_D_ORIVX_X;
   Double_t        pi2_D_ORIVX_Y;
   Double_t        pi2_D_ORIVX_Z;
   Double_t        pi2_D_ORIVX_XERR;
   Double_t        pi2_D_ORIVX_YERR;
   Double_t        pi2_D_ORIVX_ZERR;
   Double_t        pi2_D_ORIVX_CHI2;
   Int_t           pi2_D_ORIVX_NDOF;
   Float_t         pi2_D_ORIVX_COV_[3][3];
   Double_t        pi2_D_IP_ORIVX;
   Double_t        pi2_D_IPCHI2_ORIVX;
   Double_t        pi2_D_P;
   Double_t        pi2_D_PT;
   Double_t        pi2_D_PE;
   Double_t        pi2_D_PX;
   Double_t        pi2_D_PY;
   Double_t        pi2_D_PZ;
   Double_t        pi2_D_M;
   Int_t           pi2_D_ID;
   Double_t        pi2_D_PIDe;
   Double_t        pi2_D_PIDmu;
   Double_t        pi2_D_PIDK;
   Double_t        pi2_D_PIDp;
   Double_t        pi2_D_ProbNNe;
   Double_t        pi2_D_ProbNNk;
   Double_t        pi2_D_ProbNNp;
   Double_t        pi2_D_ProbNNpi;
   Double_t        pi2_D_ProbNNmu;
   Double_t        pi2_D_ProbNNghost;
   Bool_t          pi2_D_hasMuon;
   Bool_t          pi2_D_isMuon;
   Bool_t          pi2_D_hasRich;
   Bool_t          pi2_D_UsedRichAerogel;
   Bool_t          pi2_D_UsedRich1Gas;
   Bool_t          pi2_D_UsedRich2Gas;
   Bool_t          pi2_D_RichAboveElThres;
   Bool_t          pi2_D_RichAboveMuThres;
   Bool_t          pi2_D_RichAbovePiThres;
   Bool_t          pi2_D_RichAboveKaThres;
   Bool_t          pi2_D_RichAbovePrThres;
   Bool_t          pi2_D_hasCalo;
   Int_t           pi2_D_TRACK_Type;
   Int_t           pi2_D_TRACK_Key;
   Double_t        pi2_D_TRACK_CHI2NDOF;
   Double_t        pi2_D_TRACK_PCHI2;
   Double_t        pi2_D_TRACK_MatchCHI2;
   Double_t        pi2_D_TRACK_GhostProb;
   Double_t        pi2_D_TRACK_CloneDist;
   Double_t        pi2_D_TRACK_Likelihood;
   Double_t        pi2_D_cpx_1_00;
   Double_t        pi2_D_cpy_1_00;
   Double_t        pi2_D_cpz_1_00;
   Double_t        pi2_D_cpt_1_00;
   Double_t        pi2_D_cp_1_00;
   Int_t           pi2_D_cmult_1_00;
   Double_t        pi2_D_pxasy_1_00;
   Double_t        pi2_D_pyasy_1_00;
   Double_t        pi2_D_pzasy_1_00;
   Double_t        pi2_D_pasy_1_00;
   Double_t        pi2_D_ptasy_1_00;
   Double_t        Ks_ETA;
   Double_t        Ks_MINIP;
   Double_t        Ks_MINIPCHI2;
   Double_t        Ks_MINIPNEXTBEST;
   Double_t        Ks_MINIPCHI2NEXTBEST;
   Float_t         Ks_AllIP[100];   //[nPV]
   Float_t         Ks_AllIPchi2[100];   //[nPV]
   Float_t         Ks_AllDIRA[100];   //[nPV]
   Double_t        Ks_ENDVERTEX_X;
   Double_t        Ks_ENDVERTEX_Y;
   Double_t        Ks_ENDVERTEX_Z;
   Double_t        Ks_ENDVERTEX_XERR;
   Double_t        Ks_ENDVERTEX_YERR;
   Double_t        Ks_ENDVERTEX_ZERR;
   Double_t        Ks_ENDVERTEX_CHI2;
   Int_t           Ks_ENDVERTEX_NDOF;
   Float_t         Ks_ENDVERTEX_COV_[3][3];
   Double_t        Ks_OWNPV_X;
   Double_t        Ks_OWNPV_Y;
   Double_t        Ks_OWNPV_Z;
   Double_t        Ks_OWNPV_XERR;
   Double_t        Ks_OWNPV_YERR;
   Double_t        Ks_OWNPV_ZERR;
   Double_t        Ks_OWNPV_CHI2;
   Int_t           Ks_OWNPV_NDOF;
   Float_t         Ks_OWNPV_COV_[3][3];
   Double_t        Ks_IP_OWNPV;
   Double_t        Ks_IPCHI2_OWNPV;
   Double_t        Ks_FD_OWNPV;
   Double_t        Ks_FDCHI2_OWNPV;
   Double_t        Ks_DIRA_OWNPV;
   Double_t        Ks_TOPPV_X;
   Double_t        Ks_TOPPV_Y;
   Double_t        Ks_TOPPV_Z;
   Double_t        Ks_TOPPV_XERR;
   Double_t        Ks_TOPPV_YERR;
   Double_t        Ks_TOPPV_ZERR;
   Double_t        Ks_TOPPV_CHI2;
   Int_t           Ks_TOPPV_NDOF;
   Float_t         Ks_TOPPV_COV_[3][3];
   Double_t        Ks_IP_TOPPV;
   Double_t        Ks_IPCHI2_TOPPV;
   Double_t        Ks_FD_TOPPV;
   Double_t        Ks_FDCHI2_TOPPV;
   Double_t        Ks_DIRA_TOPPV;
   Double_t        Ks_ORIVX_X;
   Double_t        Ks_ORIVX_Y;
   Double_t        Ks_ORIVX_Z;
   Double_t        Ks_ORIVX_XERR;
   Double_t        Ks_ORIVX_YERR;
   Double_t        Ks_ORIVX_ZERR;
   Double_t        Ks_ORIVX_CHI2;
   Int_t           Ks_ORIVX_NDOF;
   Float_t         Ks_ORIVX_COV_[3][3];
   Double_t        Ks_IP_ORIVX;
   Double_t        Ks_IPCHI2_ORIVX;
   Double_t        Ks_FD_ORIVX;
   Double_t        Ks_FDCHI2_ORIVX;
   Double_t        Ks_DIRA_ORIVX;
   Double_t        Ks_P;
   Double_t        Ks_PT;
   Double_t        Ks_PE;
   Double_t        Ks_PX;
   Double_t        Ks_PY;
   Double_t        Ks_PZ;
   Double_t        Ks_MM;
   Double_t        Ks_MMERR;
   Double_t        Ks_M;
   Int_t           Ks_ID;
   Double_t        Ks_cpx_1_00;
   Double_t        Ks_cpy_1_00;
   Double_t        Ks_cpz_1_00;
   Double_t        Ks_cpt_1_00;
   Double_t        Ks_cp_1_00;
   Int_t           Ks_cmult_1_00;
   Double_t        Ks_pxasy_1_00;
   Double_t        Ks_pyasy_1_00;
   Double_t        Ks_pzasy_1_00;
   Double_t        Ks_pasy_1_00;
   Double_t        Ks_ptasy_1_00;
   Double_t        Ks_DOCA1;
   Double_t        Ks_TAU;
   Double_t        Ks_TAUERR;
   Double_t        pip_Ks_ETA;
   Double_t        pip_Ks_MC12TuneV2_ProbNNe;
   Double_t        pip_Ks_MC12TuneV2_ProbNNmu;
   Double_t        pip_Ks_MC12TuneV2_ProbNNpi;
   Double_t        pip_Ks_MC12TuneV2_ProbNNk;
   Double_t        pip_Ks_MC12TuneV2_ProbNNp;
   Double_t        pip_Ks_MC12TuneV2_ProbNNghost;
   Double_t        pip_Ks_MC12TuneV3_ProbNNe;
   Double_t        pip_Ks_MC12TuneV3_ProbNNmu;
   Double_t        pip_Ks_MC12TuneV3_ProbNNpi;
   Double_t        pip_Ks_MC12TuneV3_ProbNNk;
   Double_t        pip_Ks_MC12TuneV3_ProbNNp;
   Double_t        pip_Ks_MC12TuneV3_ProbNNghost;
   Double_t        pip_Ks_MC12TuneV4_ProbNNe;
   Double_t        pip_Ks_MC12TuneV4_ProbNNmu;
   Double_t        pip_Ks_MC12TuneV4_ProbNNpi;
   Double_t        pip_Ks_MC12TuneV4_ProbNNk;
   Double_t        pip_Ks_MC12TuneV4_ProbNNp;
   Double_t        pip_Ks_MC12TuneV4_ProbNNghost;
   Double_t        pip_Ks_MC15TuneV1_ProbNNe;
   Double_t        pip_Ks_MC15TuneV1_ProbNNmu;
   Double_t        pip_Ks_MC15TuneV1_ProbNNpi;
   Double_t        pip_Ks_MC15TuneV1_ProbNNk;
   Double_t        pip_Ks_MC15TuneV1_ProbNNp;
   Double_t        pip_Ks_MC15TuneV1_ProbNNghost;
   Double_t        pip_Ks_MINIP;
   Double_t        pip_Ks_MINIPCHI2;
   Double_t        pip_Ks_MINIPNEXTBEST;
   Double_t        pip_Ks_MINIPCHI2NEXTBEST;
   Float_t         pip_Ks_AllIP[100];   //[nPV]
   Float_t         pip_Ks_AllIPchi2[100];   //[nPV]
   Double_t        pip_Ks_OWNPV_X;
   Double_t        pip_Ks_OWNPV_Y;
   Double_t        pip_Ks_OWNPV_Z;
   Double_t        pip_Ks_OWNPV_XERR;
   Double_t        pip_Ks_OWNPV_YERR;
   Double_t        pip_Ks_OWNPV_ZERR;
   Double_t        pip_Ks_OWNPV_CHI2;
   Int_t           pip_Ks_OWNPV_NDOF;
   Float_t         pip_Ks_OWNPV_COV_[3][3];
   Double_t        pip_Ks_IP_OWNPV;
   Double_t        pip_Ks_IPCHI2_OWNPV;
   Double_t        pip_Ks_TOPPV_X;
   Double_t        pip_Ks_TOPPV_Y;
   Double_t        pip_Ks_TOPPV_Z;
   Double_t        pip_Ks_TOPPV_XERR;
   Double_t        pip_Ks_TOPPV_YERR;
   Double_t        pip_Ks_TOPPV_ZERR;
   Double_t        pip_Ks_TOPPV_CHI2;
   Int_t           pip_Ks_TOPPV_NDOF;
   Float_t         pip_Ks_TOPPV_COV_[3][3];
   Double_t        pip_Ks_IP_TOPPV;
   Double_t        pip_Ks_IPCHI2_TOPPV;
   Double_t        pip_Ks_ORIVX_X;
   Double_t        pip_Ks_ORIVX_Y;
   Double_t        pip_Ks_ORIVX_Z;
   Double_t        pip_Ks_ORIVX_XERR;
   Double_t        pip_Ks_ORIVX_YERR;
   Double_t        pip_Ks_ORIVX_ZERR;
   Double_t        pip_Ks_ORIVX_CHI2;
   Int_t           pip_Ks_ORIVX_NDOF;
   Float_t         pip_Ks_ORIVX_COV_[3][3];
   Double_t        pip_Ks_IP_ORIVX;
   Double_t        pip_Ks_IPCHI2_ORIVX;
   Double_t        pip_Ks_P;
   Double_t        pip_Ks_PT;
   Double_t        pip_Ks_PE;
   Double_t        pip_Ks_PX;
   Double_t        pip_Ks_PY;
   Double_t        pip_Ks_PZ;
   Double_t        pip_Ks_M;
   Int_t           pip_Ks_ID;
   Double_t        pip_Ks_PIDe;
   Double_t        pip_Ks_PIDmu;
   Double_t        pip_Ks_PIDK;
   Double_t        pip_Ks_PIDp;
   Double_t        pip_Ks_ProbNNe;
   Double_t        pip_Ks_ProbNNk;
   Double_t        pip_Ks_ProbNNp;
   Double_t        pip_Ks_ProbNNpi;
   Double_t        pip_Ks_ProbNNmu;
   Double_t        pip_Ks_ProbNNghost;
   Bool_t          pip_Ks_hasMuon;
   Bool_t          pip_Ks_isMuon;
   Bool_t          pip_Ks_hasRich;
   Bool_t          pip_Ks_UsedRichAerogel;
   Bool_t          pip_Ks_UsedRich1Gas;
   Bool_t          pip_Ks_UsedRich2Gas;
   Bool_t          pip_Ks_RichAboveElThres;
   Bool_t          pip_Ks_RichAboveMuThres;
   Bool_t          pip_Ks_RichAbovePiThres;
   Bool_t          pip_Ks_RichAboveKaThres;
   Bool_t          pip_Ks_RichAbovePrThres;
   Bool_t          pip_Ks_hasCalo;
   Int_t           pip_Ks_TRACK_Type;
   Int_t           pip_Ks_TRACK_Key;
   Double_t        pip_Ks_TRACK_CHI2NDOF;
   Double_t        pip_Ks_TRACK_PCHI2;
   Double_t        pip_Ks_TRACK_MatchCHI2;
   Double_t        pip_Ks_TRACK_GhostProb;
   Double_t        pip_Ks_TRACK_CloneDist;
   Double_t        pip_Ks_TRACK_Likelihood;
   Double_t        pip_Ks_cpx_1_00;
   Double_t        pip_Ks_cpy_1_00;
   Double_t        pip_Ks_cpz_1_00;
   Double_t        pip_Ks_cpt_1_00;
   Double_t        pip_Ks_cp_1_00;
   Int_t           pip_Ks_cmult_1_00;
   Double_t        pip_Ks_pxasy_1_00;
   Double_t        pip_Ks_pyasy_1_00;
   Double_t        pip_Ks_pzasy_1_00;
   Double_t        pip_Ks_pasy_1_00;
   Double_t        pip_Ks_ptasy_1_00;
   Double_t        pim_Ks_ETA;
   Double_t        pim_Ks_MC12TuneV2_ProbNNe;
   Double_t        pim_Ks_MC12TuneV2_ProbNNmu;
   Double_t        pim_Ks_MC12TuneV2_ProbNNpi;
   Double_t        pim_Ks_MC12TuneV2_ProbNNk;
   Double_t        pim_Ks_MC12TuneV2_ProbNNp;
   Double_t        pim_Ks_MC12TuneV2_ProbNNghost;
   Double_t        pim_Ks_MC12TuneV3_ProbNNe;
   Double_t        pim_Ks_MC12TuneV3_ProbNNmu;
   Double_t        pim_Ks_MC12TuneV3_ProbNNpi;
   Double_t        pim_Ks_MC12TuneV3_ProbNNk;
   Double_t        pim_Ks_MC12TuneV3_ProbNNp;
   Double_t        pim_Ks_MC12TuneV3_ProbNNghost;
   Double_t        pim_Ks_MC12TuneV4_ProbNNe;
   Double_t        pim_Ks_MC12TuneV4_ProbNNmu;
   Double_t        pim_Ks_MC12TuneV4_ProbNNpi;
   Double_t        pim_Ks_MC12TuneV4_ProbNNk;
   Double_t        pim_Ks_MC12TuneV4_ProbNNp;
   Double_t        pim_Ks_MC12TuneV4_ProbNNghost;
   Double_t        pim_Ks_MC15TuneV1_ProbNNe;
   Double_t        pim_Ks_MC15TuneV1_ProbNNmu;
   Double_t        pim_Ks_MC15TuneV1_ProbNNpi;
   Double_t        pim_Ks_MC15TuneV1_ProbNNk;
   Double_t        pim_Ks_MC15TuneV1_ProbNNp;
   Double_t        pim_Ks_MC15TuneV1_ProbNNghost;
   Double_t        pim_Ks_MINIP;
   Double_t        pim_Ks_MINIPCHI2;
   Double_t        pim_Ks_MINIPNEXTBEST;
   Double_t        pim_Ks_MINIPCHI2NEXTBEST;
   Float_t         pim_Ks_AllIP[100];   //[nPV]
   Float_t         pim_Ks_AllIPchi2[100];   //[nPV]
   Double_t        pim_Ks_OWNPV_X;
   Double_t        pim_Ks_OWNPV_Y;
   Double_t        pim_Ks_OWNPV_Z;
   Double_t        pim_Ks_OWNPV_XERR;
   Double_t        pim_Ks_OWNPV_YERR;
   Double_t        pim_Ks_OWNPV_ZERR;
   Double_t        pim_Ks_OWNPV_CHI2;
   Int_t           pim_Ks_OWNPV_NDOF;
   Float_t         pim_Ks_OWNPV_COV_[3][3];
   Double_t        pim_Ks_IP_OWNPV;
   Double_t        pim_Ks_IPCHI2_OWNPV;
   Double_t        pim_Ks_TOPPV_X;
   Double_t        pim_Ks_TOPPV_Y;
   Double_t        pim_Ks_TOPPV_Z;
   Double_t        pim_Ks_TOPPV_XERR;
   Double_t        pim_Ks_TOPPV_YERR;
   Double_t        pim_Ks_TOPPV_ZERR;
   Double_t        pim_Ks_TOPPV_CHI2;
   Int_t           pim_Ks_TOPPV_NDOF;
   Float_t         pim_Ks_TOPPV_COV_[3][3];
   Double_t        pim_Ks_IP_TOPPV;
   Double_t        pim_Ks_IPCHI2_TOPPV;
   Double_t        pim_Ks_ORIVX_X;
   Double_t        pim_Ks_ORIVX_Y;
   Double_t        pim_Ks_ORIVX_Z;
   Double_t        pim_Ks_ORIVX_XERR;
   Double_t        pim_Ks_ORIVX_YERR;
   Double_t        pim_Ks_ORIVX_ZERR;
   Double_t        pim_Ks_ORIVX_CHI2;
   Int_t           pim_Ks_ORIVX_NDOF;
   Float_t         pim_Ks_ORIVX_COV_[3][3];
   Double_t        pim_Ks_IP_ORIVX;
   Double_t        pim_Ks_IPCHI2_ORIVX;
   Double_t        pim_Ks_P;
   Double_t        pim_Ks_PT;
   Double_t        pim_Ks_PE;
   Double_t        pim_Ks_PX;
   Double_t        pim_Ks_PY;
   Double_t        pim_Ks_PZ;
   Double_t        pim_Ks_M;
   Int_t           pim_Ks_ID;
   Double_t        pim_Ks_PIDe;
   Double_t        pim_Ks_PIDmu;
   Double_t        pim_Ks_PIDK;
   Double_t        pim_Ks_PIDp;
   Double_t        pim_Ks_ProbNNe;
   Double_t        pim_Ks_ProbNNk;
   Double_t        pim_Ks_ProbNNp;
   Double_t        pim_Ks_ProbNNpi;
   Double_t        pim_Ks_ProbNNmu;
   Double_t        pim_Ks_ProbNNghost;
   Bool_t          pim_Ks_hasMuon;
   Bool_t          pim_Ks_isMuon;
   Bool_t          pim_Ks_hasRich;
   Bool_t          pim_Ks_UsedRichAerogel;
   Bool_t          pim_Ks_UsedRich1Gas;
   Bool_t          pim_Ks_UsedRich2Gas;
   Bool_t          pim_Ks_RichAboveElThres;
   Bool_t          pim_Ks_RichAboveMuThres;
   Bool_t          pim_Ks_RichAbovePiThres;
   Bool_t          pim_Ks_RichAboveKaThres;
   Bool_t          pim_Ks_RichAbovePrThres;
   Bool_t          pim_Ks_hasCalo;
   Int_t           pim_Ks_TRACK_Type;
   Int_t           pim_Ks_TRACK_Key;
   Double_t        pim_Ks_TRACK_CHI2NDOF;
   Double_t        pim_Ks_TRACK_PCHI2;
   Double_t        pim_Ks_TRACK_MatchCHI2;
   Double_t        pim_Ks_TRACK_GhostProb;
   Double_t        pim_Ks_TRACK_CloneDist;
   Double_t        pim_Ks_TRACK_Likelihood;
   Double_t        pim_Ks_cpx_1_00;
   Double_t        pim_Ks_cpy_1_00;
   Double_t        pim_Ks_cpz_1_00;
   Double_t        pim_Ks_cpt_1_00;
   Double_t        pim_Ks_cp_1_00;
   Int_t           pim_Ks_cmult_1_00;
   Double_t        pim_Ks_pxasy_1_00;
   Double_t        pim_Ks_pyasy_1_00;
   Double_t        pim_Ks_pzasy_1_00;
   Double_t        pim_Ks_pasy_1_00;
   Double_t        pim_Ks_ptasy_1_00;
   Double_t        pi_ETA;
   Double_t        pi_MC12TuneV2_ProbNNe;
   Double_t        pi_MC12TuneV2_ProbNNmu;
   Double_t        pi_MC12TuneV2_ProbNNpi;
   Double_t        pi_MC12TuneV2_ProbNNk;
   Double_t        pi_MC12TuneV2_ProbNNp;
   Double_t        pi_MC12TuneV2_ProbNNghost;
   Double_t        pi_MC12TuneV3_ProbNNe;
   Double_t        pi_MC12TuneV3_ProbNNmu;
   Double_t        pi_MC12TuneV3_ProbNNpi;
   Double_t        pi_MC12TuneV3_ProbNNk;
   Double_t        pi_MC12TuneV3_ProbNNp;
   Double_t        pi_MC12TuneV3_ProbNNghost;
   Double_t        pi_MC12TuneV4_ProbNNe;
   Double_t        pi_MC12TuneV4_ProbNNmu;
   Double_t        pi_MC12TuneV4_ProbNNpi;
   Double_t        pi_MC12TuneV4_ProbNNk;
   Double_t        pi_MC12TuneV4_ProbNNp;
   Double_t        pi_MC12TuneV4_ProbNNghost;
   Double_t        pi_MC15TuneV1_ProbNNe;
   Double_t        pi_MC15TuneV1_ProbNNmu;
   Double_t        pi_MC15TuneV1_ProbNNpi;
   Double_t        pi_MC15TuneV1_ProbNNk;
   Double_t        pi_MC15TuneV1_ProbNNp;
   Double_t        pi_MC15TuneV1_ProbNNghost;
   Double_t        pi_MINIP;
   Double_t        pi_MINIPCHI2;
   Double_t        pi_MINIPNEXTBEST;
   Double_t        pi_MINIPCHI2NEXTBEST;
   Float_t         pi_AllIP[100];   //[nPV]
   Float_t         pi_AllIPchi2[100];   //[nPV]
   Double_t        pi_OWNPV_X;
   Double_t        pi_OWNPV_Y;
   Double_t        pi_OWNPV_Z;
   Double_t        pi_OWNPV_XERR;
   Double_t        pi_OWNPV_YERR;
   Double_t        pi_OWNPV_ZERR;
   Double_t        pi_OWNPV_CHI2;
   Int_t           pi_OWNPV_NDOF;
   Float_t         pi_OWNPV_COV_[3][3];
   Double_t        pi_IP_OWNPV;
   Double_t        pi_IPCHI2_OWNPV;
   Double_t        pi_TOPPV_X;
   Double_t        pi_TOPPV_Y;
   Double_t        pi_TOPPV_Z;
   Double_t        pi_TOPPV_XERR;
   Double_t        pi_TOPPV_YERR;
   Double_t        pi_TOPPV_ZERR;
   Double_t        pi_TOPPV_CHI2;
   Int_t           pi_TOPPV_NDOF;
   Float_t         pi_TOPPV_COV_[3][3];
   Double_t        pi_IP_TOPPV;
   Double_t        pi_IPCHI2_TOPPV;
   Double_t        pi_ORIVX_X;
   Double_t        pi_ORIVX_Y;
   Double_t        pi_ORIVX_Z;
   Double_t        pi_ORIVX_XERR;
   Double_t        pi_ORIVX_YERR;
   Double_t        pi_ORIVX_ZERR;
   Double_t        pi_ORIVX_CHI2;
   Int_t           pi_ORIVX_NDOF;
   Float_t         pi_ORIVX_COV_[3][3];
   Double_t        pi_IP_ORIVX;
   Double_t        pi_IPCHI2_ORIVX;
   Double_t        pi_P;
   Double_t        pi_PT;
   Double_t        pi_PE;
   Double_t        pi_PX;
   Double_t        pi_PY;
   Double_t        pi_PZ;
   Double_t        pi_M;
   Int_t           pi_ID;
   Double_t        pi_PIDe;
   Double_t        pi_PIDmu;
   Double_t        pi_PIDK;
   Double_t        pi_PIDp;
   Double_t        pi_ProbNNe;
   Double_t        pi_ProbNNk;
   Double_t        pi_ProbNNp;
   Double_t        pi_ProbNNpi;
   Double_t        pi_ProbNNmu;
   Double_t        pi_ProbNNghost;
   Bool_t          pi_hasMuon;
   Bool_t          pi_isMuon;
   Bool_t          pi_hasRich;
   Bool_t          pi_UsedRichAerogel;
   Bool_t          pi_UsedRich1Gas;
   Bool_t          pi_UsedRich2Gas;
   Bool_t          pi_RichAboveElThres;
   Bool_t          pi_RichAboveMuThres;
   Bool_t          pi_RichAbovePiThres;
   Bool_t          pi_RichAboveKaThres;
   Bool_t          pi_RichAbovePrThres;
   Bool_t          pi_hasCalo;
   Int_t           pi_TRACK_Type;
   Int_t           pi_TRACK_Key;
   Double_t        pi_TRACK_CHI2NDOF;
   Double_t        pi_TRACK_PCHI2;
   Double_t        pi_TRACK_MatchCHI2;
   Double_t        pi_TRACK_GhostProb;
   Double_t        pi_TRACK_CloneDist;
   Double_t        pi_TRACK_Likelihood;
   Double_t        pi_cpx_1_00;
   Double_t        pi_cpy_1_00;
   Double_t        pi_cpz_1_00;
   Double_t        pi_cpt_1_00;
   Double_t        pi_cp_1_00;
   Int_t           pi_cmult_1_00;
   Double_t        pi_pxasy_1_00;
   Double_t        pi_pyasy_1_00;
   Double_t        pi_pzasy_1_00;
   Double_t        pi_pasy_1_00;
   Double_t        pi_ptasy_1_00;
   UInt_t          nCandidate;
   ULong64_t       totCandidates;
   ULong64_t       EventInSequence;
   UInt_t          runNumber;
   ULong64_t       eventNumber;
   UInt_t          BCID;
   Int_t           BCType;
   UInt_t          OdinTCK;
   UInt_t          L0DUTCK;
   UInt_t          HLT1TCK;
   UInt_t          HLT2TCK;
   ULong64_t       GpsTime;
   Short_t         Polarity;
   Float_t         PVX[100];   //[nPV]
   Float_t         PVY[100];   //[nPV]
   Float_t         PVZ[100];   //[nPV]
   Float_t         PVXERR[100];   //[nPV]
   Float_t         PVYERR[100];   //[nPV]
   Float_t         PVZERR[100];   //[nPV]
   Float_t         PVCHI2[100];   //[nPV]
   Float_t         PVNDOF[100];   //[nPV]
   Float_t         PVNTRACKS[100];   //[nPV]
   Int_t           nPVs;
   Int_t           nTracks;
   Int_t           nLongTracks;
   Int_t           nDownstreamTracks;
   Int_t           nUpstreamTracks;
   Int_t           nVeloTracks;
   Int_t           nTTracks;
   Int_t           nBackTracks;
   Int_t           nRich1Hits;
   Int_t           nRich2Hits;
   Int_t           nVeloClusters;
   Int_t           nITClusters;
   Int_t           nTTClusters;
   Int_t           nOTClusters;
   Int_t           nSPDHits;
   Int_t           nMuonCoordsS0;
   Int_t           nMuonCoordsS1;
   Int_t           nMuonCoordsS2;
   Int_t           nMuonCoordsS3;
   Int_t           nMuonCoordsS4;
   Int_t           nMuonTracks;

   // List of branches
   TBranch        *b_B_ETA;   //!
   TBranch        *b_B_MINIP;   //!
   TBranch        *b_B_MINIPCHI2;   //!
   TBranch        *b_B_MINIPNEXTBEST;   //!
   TBranch        *b_B_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_B_AllIP;   //!
   TBranch        *b_B_AllIPchi2;   //!
   TBranch        *b_B_AllDIRA;   //!
   TBranch        *b_B_ENDVERTEX_X;   //!
   TBranch        *b_B_ENDVERTEX_Y;   //!
   TBranch        *b_B_ENDVERTEX_Z;   //!
   TBranch        *b_B_ENDVERTEX_XERR;   //!
   TBranch        *b_B_ENDVERTEX_YERR;   //!
   TBranch        *b_B_ENDVERTEX_ZERR;   //!
   TBranch        *b_B_ENDVERTEX_CHI2;   //!
   TBranch        *b_B_ENDVERTEX_NDOF;   //!
   TBranch        *b_B_ENDVERTEX_COV_;   //!
   TBranch        *b_B_OWNPV_X;   //!
   TBranch        *b_B_OWNPV_Y;   //!
   TBranch        *b_B_OWNPV_Z;   //!
   TBranch        *b_B_OWNPV_XERR;   //!
   TBranch        *b_B_OWNPV_YERR;   //!
   TBranch        *b_B_OWNPV_ZERR;   //!
   TBranch        *b_B_OWNPV_CHI2;   //!
   TBranch        *b_B_OWNPV_NDOF;   //!
   TBranch        *b_B_OWNPV_COV_;   //!
   TBranch        *b_B_IP_OWNPV;   //!
   TBranch        *b_B_IPCHI2_OWNPV;   //!
   TBranch        *b_B_FD_OWNPV;   //!
   TBranch        *b_B_FDCHI2_OWNPV;   //!
   TBranch        *b_B_DIRA_OWNPV;   //!
   TBranch        *b_B_TOPPV_X;   //!
   TBranch        *b_B_TOPPV_Y;   //!
   TBranch        *b_B_TOPPV_Z;   //!
   TBranch        *b_B_TOPPV_XERR;   //!
   TBranch        *b_B_TOPPV_YERR;   //!
   TBranch        *b_B_TOPPV_ZERR;   //!
   TBranch        *b_B_TOPPV_CHI2;   //!
   TBranch        *b_B_TOPPV_NDOF;   //!
   TBranch        *b_B_TOPPV_COV_;   //!
   TBranch        *b_B_IP_TOPPV;   //!
   TBranch        *b_B_IPCHI2_TOPPV;   //!
   TBranch        *b_B_FD_TOPPV;   //!
   TBranch        *b_B_FDCHI2_TOPPV;   //!
   TBranch        *b_B_DIRA_TOPPV;   //!
   TBranch        *b_B_P;   //!
   TBranch        *b_B_PT;   //!
   TBranch        *b_B_PE;   //!
   TBranch        *b_B_PX;   //!
   TBranch        *b_B_PY;   //!
   TBranch        *b_B_PZ;   //!
   TBranch        *b_B_MM;   //!
   TBranch        *b_B_MMERR;   //!
   TBranch        *b_B_M;   //!
   TBranch        *b_B_ID;   //!
   TBranch        *b_B_TAGDECISION;   //!
   TBranch        *b_B_TAGOMEGA;   //!
   TBranch        *b_B_TAGDECISION_OS;   //!
   TBranch        *b_B_TAGOMEGA_OS;   //!
   TBranch        *b_B_TAGGER;   //!
   TBranch        *b_B_OS_Muon_DEC;   //!
   TBranch        *b_B_OS_Muon_PROB;   //!
   TBranch        *b_B_OS_Electron_DEC;   //!
   TBranch        *b_B_OS_Electron_PROB;   //!
   TBranch        *b_B_OS_Kaon_DEC;   //!
   TBranch        *b_B_OS_Kaon_PROB;   //!
   TBranch        *b_B_SS_Kaon_DEC;   //!
   TBranch        *b_B_SS_Kaon_PROB;   //!
   TBranch        *b_B_SS_Pion_DEC;   //!
   TBranch        *b_B_SS_Pion_PROB;   //!
   TBranch        *b_B_SS_PionBDT_DEC;   //!
   TBranch        *b_B_SS_PionBDT_PROB;   //!
   TBranch        *b_B_VtxCharge_DEC;   //!
   TBranch        *b_B_VtxCharge_PROB;   //!
   TBranch        *b_B_OS_nnetKaon_DEC;   //!
   TBranch        *b_B_OS_nnetKaon_PROB;   //!
   TBranch        *b_B_SS_nnetKaon_DEC;   //!
   TBranch        *b_B_SS_nnetKaon_PROB;   //!
   TBranch        *b_B_SS_Proton_DEC;   //!
   TBranch        *b_B_SS_Proton_PROB;   //!
   TBranch        *b_B_OS_Charm_DEC;   //!
   TBranch        *b_B_OS_Charm_PROB;   //!
   TBranch        *b_B_cpx_1_00;   //!
   TBranch        *b_B_cpy_1_00;   //!
   TBranch        *b_B_cpz_1_00;   //!
   TBranch        *b_B_cpt_1_00;   //!
   TBranch        *b_B_cp_1_00;   //!
   TBranch        *b_B_cmult_1_00;   //!
   TBranch        *b_B_pxasy_1_00;   //!
   TBranch        *b_B_pyasy_1_00;   //!
   TBranch        *b_B_pzasy_1_00;   //!
   TBranch        *b_B_pasy_1_00;   //!
   TBranch        *b_B_ptasy_1_00;   //!
   TBranch        *b_B_DOCA1;   //!
   TBranch        *b_B_TAU;   //!
   TBranch        *b_B_TAUERR;   //!
   TBranch        *b_B_BDTF_nPV;   //!
   TBranch        *b_B_BDTF_Dplus_M;   //!
   TBranch        *b_B_BDTF_Dplus_MERR;   //!
   TBranch        *b_B_BDTF_Dplus_P;   //!
   TBranch        *b_B_BDTF_Dplus_PERR;   //!
   TBranch        *b_B_BDTF_Dplus_ctau;   //!
   TBranch        *b_B_BDTF_Dplus_ctauErr;   //!
   TBranch        *b_B_BDTF_Dplus_decayLength;   //!
   TBranch        *b_B_BDTF_Dplus_decayLengthErr;   //!
   TBranch        *b_B_BDTF_KS0_M;   //!
   TBranch        *b_B_BDTF_KS0_MERR;   //!
   TBranch        *b_B_BDTF_KS0_P;   //!
   TBranch        *b_B_BDTF_KS0_PERR;   //!
   TBranch        *b_B_BDTF_KS0_ctau;   //!
   TBranch        *b_B_BDTF_KS0_ctauErr;   //!
   TBranch        *b_B_BDTF_KS0_decayLength;   //!
   TBranch        *b_B_BDTF_KS0_decayLengthErr;   //!
   TBranch        *b_B_BDTF_M;   //!
   TBranch        *b_B_BDTF_MERR;   //!
   TBranch        *b_B_BDTF_P;   //!
   TBranch        *b_B_BDTF_PERR;   //!
   TBranch        *b_B_BDTF_PV_X;   //!
   TBranch        *b_B_BDTF_PV_Y;   //!
   TBranch        *b_B_BDTF_PV_Z;   //!
   TBranch        *b_B_BDTF_PV_key;   //!
   TBranch        *b_B_BDTF_chi2;   //!
   TBranch        *b_B_BDTF_ctau;   //!
   TBranch        *b_B_BDTF_ctauErr;   //!
   TBranch        *b_B_BDTF_decayLength;   //!
   TBranch        *b_B_BDTF_decayLengthErr;   //!
   TBranch        *b_B_BDTF_nDOF;   //!
   TBranch        *b_B_BDTF_nIter;   //!
   TBranch        *b_B_BDTF_status;   //!
   TBranch        *b_B_DTF_nPV;   //!
   TBranch        *b_B_DTF_Dplus_M;   //!
   TBranch        *b_B_DTF_Dplus_MERR;   //!
   TBranch        *b_B_DTF_Dplus_P;   //!
   TBranch        *b_B_DTF_Dplus_PERR;   //!
   TBranch        *b_B_DTF_Dplus_ctau;   //!
   TBranch        *b_B_DTF_Dplus_ctauErr;   //!
   TBranch        *b_B_DTF_Dplus_decayLength;   //!
   TBranch        *b_B_DTF_Dplus_decayLengthErr;   //!
   TBranch        *b_B_DTF_KS0_M;   //!
   TBranch        *b_B_DTF_KS0_MERR;   //!
   TBranch        *b_B_DTF_KS0_P;   //!
   TBranch        *b_B_DTF_KS0_PERR;   //!
   TBranch        *b_B_DTF_KS0_ctau;   //!
   TBranch        *b_B_DTF_KS0_ctauErr;   //!
   TBranch        *b_B_DTF_KS0_decayLength;   //!
   TBranch        *b_B_DTF_KS0_decayLengthErr;   //!
   TBranch        *b_B_DTF_M;   //!
   TBranch        *b_B_DTF_MERR;   //!
   TBranch        *b_B_DTF_P;   //!
   TBranch        *b_B_DTF_PERR;   //!
   TBranch        *b_B_DTF_PV_X;   //!
   TBranch        *b_B_DTF_PV_Y;   //!
   TBranch        *b_B_DTF_PV_Z;   //!
   TBranch        *b_B_DTF_PV_key;   //!
   TBranch        *b_B_DTF_chi2;   //!
   TBranch        *b_B_DTF_ctau;   //!
   TBranch        *b_B_DTF_ctauErr;   //!
   TBranch        *b_B_DTF_decayLength;   //!
   TBranch        *b_B_DTF_decayLengthErr;   //!
   TBranch        *b_B_DTF_nDOF;   //!
   TBranch        *b_B_DTF_nIter;   //!
   TBranch        *b_B_DTF_status;   //!
    
   TBranch        *b_B_FullDTF_nPV;   //!
   TBranch        *b_B_FullDTF_Dplus_Kplus_ID;   //!
   TBranch        *b_B_FullDTF_Dplus_Kplus_PE;   //!
   TBranch        *b_B_FullDTF_Dplus_Kplus_PX;   //!
   TBranch        *b_B_FullDTF_Dplus_Kplus_PY;   //!
   TBranch        *b_B_FullDTF_Dplus_Kplus_PZ;   //!
   TBranch        *b_B_FullDTF_Dplus_M;   //!
   TBranch        *b_B_FullDTF_Dplus_MERR;   //!
   TBranch        *b_B_FullDTF_Dplus_P;   //!
   TBranch        *b_B_FullDTF_Dplus_PERR;   //!
   TBranch        *b_B_FullDTF_Dplus_ctau;   //!
   TBranch        *b_B_FullDTF_Dplus_ctauErr;   //!
   TBranch        *b_B_FullDTF_Dplus_decayLength;   //!
   TBranch        *b_B_FullDTF_Dplus_decayLengthErr;   //!
   TBranch        *b_B_FullDTF_Dplus_piplus_0_ID;   //!
   TBranch        *b_B_FullDTF_Dplus_piplus_0_PE;   //!
   TBranch        *b_B_FullDTF_Dplus_piplus_0_PX;   //!
   TBranch        *b_B_FullDTF_Dplus_piplus_0_PY;   //!
   TBranch        *b_B_FullDTF_Dplus_piplus_0_PZ;   //!
   TBranch        *b_B_FullDTF_Dplus_piplus_ID;   //!
   TBranch        *b_B_FullDTF_Dplus_piplus_PE;   //!
   TBranch        *b_B_FullDTF_Dplus_piplus_PX;   //!
   TBranch        *b_B_FullDTF_Dplus_piplus_PY;   //!
   TBranch        *b_B_FullDTF_Dplus_piplus_PZ;   //!
   TBranch        *b_B_FullDTF_KS0_M;   //!
   TBranch        *b_B_FullDTF_KS0_MERR;   //!
   TBranch        *b_B_FullDTF_KS0_P;   //!
   TBranch        *b_B_FullDTF_KS0_PERR;   //!
   TBranch        *b_B_FullDTF_KS0_ctau;   //!
   TBranch        *b_B_FullDTF_KS0_ctauErr;   //!
   TBranch        *b_B_FullDTF_KS0_decayLength;   //!
   TBranch        *b_B_FullDTF_KS0_decayLengthErr;   //!
   TBranch        *b_B_FullDTF_KS0_piplus_0_ID;   //!
   TBranch        *b_B_FullDTF_KS0_piplus_0_PE;   //!
   TBranch        *b_B_FullDTF_KS0_piplus_0_PX;   //!
   TBranch        *b_B_FullDTF_KS0_piplus_0_PY;   //!
   TBranch        *b_B_FullDTF_KS0_piplus_0_PZ;   //!
   TBranch        *b_B_FullDTF_KS0_piplus_ID;   //!
   TBranch        *b_B_FullDTF_KS0_piplus_PE;   //!
   TBranch        *b_B_FullDTF_KS0_piplus_PX;   //!
   TBranch        *b_B_FullDTF_KS0_piplus_PY;   //!
   TBranch        *b_B_FullDTF_KS0_piplus_PZ;   //!
   TBranch        *b_B_FullDTF_M;   //!
   TBranch        *b_B_FullDTF_MERR;   //!
   TBranch        *b_B_FullDTF_P;   //!
   TBranch        *b_B_FullDTF_PERR;   //!
   TBranch        *b_B_FullDTF_PV_X;   //!
   TBranch        *b_B_FullDTF_PV_Y;   //!
   TBranch        *b_B_FullDTF_PV_Z;   //!
   TBranch        *b_B_FullDTF_PV_key;   //!
   TBranch        *b_B_FullDTF_chi2;   //!
   TBranch        *b_B_FullDTF_ctau;   //!
   TBranch        *b_B_FullDTF_ctauErr;   //!
   TBranch        *b_B_FullDTF_decayLength;   //!
   TBranch        *b_B_FullDTF_decayLengthErr;   //!
   TBranch        *b_B_FullDTF_nDOF;   //!
   TBranch        *b_B_FullDTF_nIter;   //!
   TBranch        *b_B_FullDTF_piplus_ID;   //!
   TBranch        *b_B_FullDTF_piplus_PE;   //!
   TBranch        *b_B_FullDTF_piplus_PX;   //!
   TBranch        *b_B_FullDTF_piplus_PY;   //!
   TBranch        *b_B_FullDTF_piplus_PZ;   //!
   TBranch        *b_B_FullDTF_status;   //!

    TBranch        *b_B_FullBsDTF_nPV;   //!
    TBranch        *b_B_FullBsDTF_Dplus_Kplus_ID;   //!
    TBranch        *b_B_FullBsDTF_Dplus_Kplus_PE;   //!
    TBranch        *b_B_FullBsDTF_Dplus_Kplus_PX;   //!
    TBranch        *b_B_FullBsDTF_Dplus_Kplus_PY;   //!
    TBranch        *b_B_FullBsDTF_Dplus_Kplus_PZ;   //!
    TBranch        *b_B_FullBsDTF_Dplus_M;   //!
    TBranch        *b_B_FullBsDTF_Dplus_MERR;   //!
    TBranch        *b_B_FullBsDTF_Dplus_P;   //!
    TBranch        *b_B_FullBsDTF_Dplus_PERR;   //!
    TBranch        *b_B_FullBsDTF_Dplus_ctau;   //!
    TBranch        *b_B_FullBsDTF_Dplus_ctauErr;   //!
    TBranch        *b_B_FullBsDTF_Dplus_decayLength;   //!
    TBranch        *b_B_FullBsDTF_Dplus_decayLengthErr;   //!
    TBranch        *b_B_FullBsDTF_Dplus_piplus_0_ID;   //!
    TBranch        *b_B_FullBsDTF_Dplus_piplus_0_PE;   //!
    TBranch        *b_B_FullBsDTF_Dplus_piplus_0_PX;   //!
    TBranch        *b_B_FullBsDTF_Dplus_piplus_0_PY;   //!
    TBranch        *b_B_FullBsDTF_Dplus_piplus_0_PZ;   //!
    TBranch        *b_B_FullBsDTF_Dplus_piplus_ID;   //!
    TBranch        *b_B_FullBsDTF_Dplus_piplus_PE;   //!
    TBranch        *b_B_FullBsDTF_Dplus_piplus_PX;   //!
    TBranch        *b_B_FullBsDTF_Dplus_piplus_PY;   //!
    TBranch        *b_B_FullBsDTF_Dplus_piplus_PZ;   //!
    TBranch        *b_B_FullBsDTF_KS0_M;   //!
    TBranch        *b_B_FullBsDTF_KS0_MERR;   //!
    TBranch        *b_B_FullBsDTF_KS0_P;   //!
    TBranch        *b_B_FullBsDTF_KS0_PERR;   //!
    TBranch        *b_B_FullBsDTF_KS0_ctau;   //!
    TBranch        *b_B_FullBsDTF_KS0_ctauErr;   //!
    TBranch        *b_B_FullBsDTF_KS0_decayLength;   //!
    TBranch        *b_B_FullBsDTF_KS0_decayLengthErr;   //!
    TBranch        *b_B_FullBsDTF_KS0_piplus_0_ID;   //!
    TBranch        *b_B_FullBsDTF_KS0_piplus_0_PE;   //!
    TBranch        *b_B_FullBsDTF_KS0_piplus_0_PX;   //!
    TBranch        *b_B_FullBsDTF_KS0_piplus_0_PY;   //!
    TBranch        *b_B_FullBsDTF_KS0_piplus_0_PZ;   //!
    TBranch        *b_B_FullBsDTF_KS0_piplus_ID;   //!
    TBranch        *b_B_FullBsDTF_KS0_piplus_PE;   //!
    TBranch        *b_B_FullBsDTF_KS0_piplus_PX;   //!
    TBranch        *b_B_FullBsDTF_KS0_piplus_PY;   //!
    TBranch        *b_B_FullBsDTF_KS0_piplus_PZ;   //!
    TBranch        *b_B_FullBsDTF_M;   //!
    TBranch        *b_B_FullBsDTF_MERR;   //!
    TBranch        *b_B_FullBsDTF_P;   //!
    TBranch        *b_B_FullBsDTF_PERR;   //!
    TBranch        *b_B_FullBsDTF_PV_X;   //!
    TBranch        *b_B_FullBsDTF_PV_Y;   //!
    TBranch        *b_B_FullBsDTF_PV_Z;   //!
    TBranch        *b_B_FullBsDTF_PV_key;   //!
    TBranch        *b_B_FullBsDTF_chi2;   //!
    TBranch        *b_B_FullBsDTF_ctau;   //!
    TBranch        *b_B_FullBsDTF_ctauErr;   //!
    TBranch        *b_B_FullBsDTF_decayLength;   //!
    TBranch        *b_B_FullBsDTF_decayLengthErr;   //!
    TBranch        *b_B_FullBsDTF_nDOF;   //!
    TBranch        *b_B_FullBsDTF_nIter;   //!
    TBranch        *b_B_FullBsDTF_piplus_ID;   //!
    TBranch        *b_B_FullBsDTF_piplus_PE;   //!
    TBranch        *b_B_FullBsDTF_piplus_PX;   //!
    TBranch        *b_B_FullBsDTF_piplus_PY;   //!
    TBranch        *b_B_FullBsDTF_piplus_PZ;   //!
    TBranch        *b_B_FullBsDTF_status;   //!

    
   TBranch        *b_B_PV_nPV;   //!
   TBranch        *b_B_PV_Dplus_Kplus_ID;   //!
   TBranch        *b_B_PV_Dplus_Kplus_PE;   //!
   TBranch        *b_B_PV_Dplus_Kplus_PX;   //!
   TBranch        *b_B_PV_Dplus_Kplus_PY;   //!
   TBranch        *b_B_PV_Dplus_Kplus_PZ;   //!
   TBranch        *b_B_PV_Dplus_M;   //!
   TBranch        *b_B_PV_Dplus_MERR;   //!
   TBranch        *b_B_PV_Dplus_P;   //!
   TBranch        *b_B_PV_Dplus_PERR;   //!
   TBranch        *b_B_PV_Dplus_ctau;   //!
   TBranch        *b_B_PV_Dplus_ctauErr;   //!
   TBranch        *b_B_PV_Dplus_decayLength;   //!
   TBranch        *b_B_PV_Dplus_decayLengthErr;   //!
   TBranch        *b_B_PV_Dplus_piplus_0_ID;   //!
   TBranch        *b_B_PV_Dplus_piplus_0_PE;   //!
   TBranch        *b_B_PV_Dplus_piplus_0_PX;   //!
   TBranch        *b_B_PV_Dplus_piplus_0_PY;   //!
   TBranch        *b_B_PV_Dplus_piplus_0_PZ;   //!
   TBranch        *b_B_PV_Dplus_piplus_ID;   //!
   TBranch        *b_B_PV_Dplus_piplus_PE;   //!
   TBranch        *b_B_PV_Dplus_piplus_PX;   //!
   TBranch        *b_B_PV_Dplus_piplus_PY;   //!
   TBranch        *b_B_PV_Dplus_piplus_PZ;   //!
   TBranch        *b_B_PV_KS0_M;   //!
   TBranch        *b_B_PV_KS0_MERR;   //!
   TBranch        *b_B_PV_KS0_P;   //!
   TBranch        *b_B_PV_KS0_PERR;   //!
   TBranch        *b_B_PV_KS0_ctau;   //!
   TBranch        *b_B_PV_KS0_ctauErr;   //!
   TBranch        *b_B_PV_KS0_decayLength;   //!
   TBranch        *b_B_PV_KS0_decayLengthErr;   //!
   TBranch        *b_B_PV_KS0_piplus_0_ID;   //!
   TBranch        *b_B_PV_KS0_piplus_0_PE;   //!
   TBranch        *b_B_PV_KS0_piplus_0_PX;   //!
   TBranch        *b_B_PV_KS0_piplus_0_PY;   //!
   TBranch        *b_B_PV_KS0_piplus_0_PZ;   //!
   TBranch        *b_B_PV_KS0_piplus_ID;   //!
   TBranch        *b_B_PV_KS0_piplus_PE;   //!
   TBranch        *b_B_PV_KS0_piplus_PX;   //!
   TBranch        *b_B_PV_KS0_piplus_PY;   //!
   TBranch        *b_B_PV_KS0_piplus_PZ;   //!
   TBranch        *b_B_PV_M;   //!
   TBranch        *b_B_PV_MERR;   //!
   TBranch        *b_B_PV_P;   //!
   TBranch        *b_B_PV_PERR;   //!
   TBranch        *b_B_PV_PV_X;   //!
   TBranch        *b_B_PV_PV_Y;   //!
   TBranch        *b_B_PV_PV_Z;   //!
   TBranch        *b_B_PV_PV_key;   //!
   TBranch        *b_B_PV_chi2;   //!
   TBranch        *b_B_PV_ctau;   //!
   TBranch        *b_B_PV_ctauErr;   //!
   TBranch        *b_B_PV_decayLength;   //!
   TBranch        *b_B_PV_decayLengthErr;   //!
   TBranch        *b_B_PV_nDOF;   //!
   TBranch        *b_B_PV_nIter;   //!
   TBranch        *b_B_PV_piplus_ID;   //!
   TBranch        *b_B_PV_piplus_PE;   //!
   TBranch        *b_B_PV_piplus_PX;   //!
   TBranch        *b_B_PV_piplus_PY;   //!
   TBranch        *b_B_PV_piplus_PZ;   //!
   TBranch        *b_B_PV_status;   //!
    
    TBranch        *b_B_Hlt1TrackAllL0Decision_Dec;   //!
    TBranch        *b_B_Hlt1TrackAllL0Decision_TIS;   //!
    TBranch        *b_B_Hlt1TrackAllL0Decision_TOS;   //!
    TBranch        *b_B_Hlt2Topo2BodyBBDTDecision_Dec;   //!
    TBranch        *b_B_Hlt2Topo2BodyBBDTDecision_TIS;   //!
    TBranch        *b_B_Hlt2Topo2BodyBBDTDecision_TOS;   //!
    TBranch        *b_B_Hlt2Topo3BodyBBDTDecision_Dec;   //!
    TBranch        *b_B_Hlt2Topo3BodyBBDTDecision_TIS;   //!
    TBranch        *b_B_Hlt2Topo3BodyBBDTDecision_TOS;   //!
    TBranch        *b_B_Hlt2Topo4BodyBBDTDecision_Dec;   //!
    TBranch        *b_B_Hlt2Topo4BodyBBDTDecision_TIS;   //!
    TBranch        *b_B_Hlt2Topo4BodyBBDTDecision_TOS;   //!
    
   TBranch        *b_B_L0Global_Dec;   //!
   TBranch        *b_B_L0Global_TIS;   //!
   TBranch        *b_B_L0Global_TOS;   //!
   TBranch        *b_B_Hlt1Global_Dec;   //!
   TBranch        *b_B_Hlt1Global_TIS;   //!
   TBranch        *b_B_Hlt1Global_TOS;   //!
   TBranch        *b_B_Hlt1Phys_Dec;   //!
   TBranch        *b_B_Hlt1Phys_TIS;   //!
   TBranch        *b_B_Hlt1Phys_TOS;   //!
   TBranch        *b_B_Hlt2Global_Dec;   //!
   TBranch        *b_B_Hlt2Global_TIS;   //!
   TBranch        *b_B_Hlt2Global_TOS;   //!
   TBranch        *b_B_Hlt2Phys_Dec;   //!
   TBranch        *b_B_Hlt2Phys_TIS;   //!
   TBranch        *b_B_Hlt2Phys_TOS;   //!
   TBranch        *b_B_L0HadronDecision_Dec;   //!
   TBranch        *b_B_L0HadronDecision_TIS;   //!
   TBranch        *b_B_L0HadronDecision_TOS;   //!
   TBranch        *b_B_L0MuonDecision_Dec;   //!
   TBranch        *b_B_L0MuonDecision_TIS;   //!
   TBranch        *b_B_L0MuonDecision_TOS;   //!
   TBranch        *b_B_L0DiMuonDecision_Dec;   //!
   TBranch        *b_B_L0DiMuonDecision_TIS;   //!
   TBranch        *b_B_L0DiMuonDecision_TOS;   //!
   TBranch        *b_B_L0ElectronDecision_Dec;   //!
   TBranch        *b_B_L0ElectronDecision_TIS;   //!
   TBranch        *b_B_L0ElectronDecision_TOS;   //!
   TBranch        *b_B_L0PhotonDecision_Dec;   //!
   TBranch        *b_B_L0PhotonDecision_TIS;   //!
   TBranch        *b_B_L0PhotonDecision_TOS;   //!
   TBranch        *b_B_Hlt1TrackMVADecision_Dec;   //!
   TBranch        *b_B_Hlt1TrackMVADecision_TIS;   //!
   TBranch        *b_B_Hlt1TrackMVADecision_TOS;   //!
   TBranch        *b_B_Hlt1TwoTrackMVADecision_Dec;   //!
   TBranch        *b_B_Hlt1TwoTrackMVADecision_TIS;   //!
   TBranch        *b_B_Hlt1TwoTrackMVADecision_TOS;   //!
   TBranch        *b_B_Hlt2IncPhiDecision_Dec;   //!
   TBranch        *b_B_Hlt2IncPhiDecision_TIS;   //!
   TBranch        *b_B_Hlt2IncPhiDecision_TOS;   //!
   TBranch        *b_B_Hlt2PhiIncPhiDecision_Dec;   //!
   TBranch        *b_B_Hlt2PhiIncPhiDecision_TIS;   //!
   TBranch        *b_B_Hlt2PhiIncPhiDecision_TOS;   //!
   TBranch        *b_B_Hlt2Topo2BodyDecision_Dec;   //!
   TBranch        *b_B_Hlt2Topo2BodyDecision_TIS;   //!
   TBranch        *b_B_Hlt2Topo2BodyDecision_TOS;   //!
   TBranch        *b_B_Hlt2Topo3BodyDecision_Dec;   //!
   TBranch        *b_B_Hlt2Topo3BodyDecision_TIS;   //!
   TBranch        *b_B_Hlt2Topo3BodyDecision_TOS;   //!
   TBranch        *b_B_Hlt2Topo4BodyDecision_Dec;   //!
   TBranch        *b_B_Hlt2Topo4BodyDecision_TIS;   //!
   TBranch        *b_B_Hlt2Topo4BodyDecision_TOS;   //!
   TBranch        *b_D_ETA;   //!
   TBranch        *b_D_MINIP;   //!
   TBranch        *b_D_MINIPCHI2;   //!
   TBranch        *b_D_MINIPNEXTBEST;   //!
   TBranch        *b_D_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_D_AllIP;   //!
   TBranch        *b_D_AllIPchi2;   //!
   TBranch        *b_D_AllDIRA;   //!
   TBranch        *b_D_ENDVERTEX_X;   //!
   TBranch        *b_D_ENDVERTEX_Y;   //!
   TBranch        *b_D_ENDVERTEX_Z;   //!
   TBranch        *b_D_ENDVERTEX_XERR;   //!
   TBranch        *b_D_ENDVERTEX_YERR;   //!
   TBranch        *b_D_ENDVERTEX_ZERR;   //!
   TBranch        *b_D_ENDVERTEX_CHI2;   //!
   TBranch        *b_D_ENDVERTEX_NDOF;   //!
   TBranch        *b_D_ENDVERTEX_COV_;   //!
   TBranch        *b_D_OWNPV_X;   //!
   TBranch        *b_D_OWNPV_Y;   //!
   TBranch        *b_D_OWNPV_Z;   //!
   TBranch        *b_D_OWNPV_XERR;   //!
   TBranch        *b_D_OWNPV_YERR;   //!
   TBranch        *b_D_OWNPV_ZERR;   //!
   TBranch        *b_D_OWNPV_CHI2;   //!
   TBranch        *b_D_OWNPV_NDOF;   //!
   TBranch        *b_D_OWNPV_COV_;   //!
   TBranch        *b_D_IP_OWNPV;   //!
   TBranch        *b_D_IPCHI2_OWNPV;   //!
   TBranch        *b_D_FD_OWNPV;   //!
   TBranch        *b_D_FDCHI2_OWNPV;   //!
   TBranch        *b_D_DIRA_OWNPV;   //!
   TBranch        *b_D_TOPPV_X;   //!
   TBranch        *b_D_TOPPV_Y;   //!
   TBranch        *b_D_TOPPV_Z;   //!
   TBranch        *b_D_TOPPV_XERR;   //!
   TBranch        *b_D_TOPPV_YERR;   //!
   TBranch        *b_D_TOPPV_ZERR;   //!
   TBranch        *b_D_TOPPV_CHI2;   //!
   TBranch        *b_D_TOPPV_NDOF;   //!
   TBranch        *b_D_TOPPV_COV_;   //!
   TBranch        *b_D_IP_TOPPV;   //!
   TBranch        *b_D_IPCHI2_TOPPV;   //!
   TBranch        *b_D_FD_TOPPV;   //!
   TBranch        *b_D_FDCHI2_TOPPV;   //!
   TBranch        *b_D_DIRA_TOPPV;   //!
   TBranch        *b_D_ORIVX_X;   //!
   TBranch        *b_D_ORIVX_Y;   //!
   TBranch        *b_D_ORIVX_Z;   //!
   TBranch        *b_D_ORIVX_XERR;   //!
   TBranch        *b_D_ORIVX_YERR;   //!
   TBranch        *b_D_ORIVX_ZERR;   //!
   TBranch        *b_D_ORIVX_CHI2;   //!
   TBranch        *b_D_ORIVX_NDOF;   //!
   TBranch        *b_D_ORIVX_COV_;   //!
   TBranch        *b_D_IP_ORIVX;   //!
   TBranch        *b_D_IPCHI2_ORIVX;   //!
   TBranch        *b_D_FD_ORIVX;   //!
   TBranch        *b_D_FDCHI2_ORIVX;   //!
   TBranch        *b_D_DIRA_ORIVX;   //!
   TBranch        *b_D_P;   //!
   TBranch        *b_D_PT;   //!
   TBranch        *b_D_PE;   //!
   TBranch        *b_D_PX;   //!
   TBranch        *b_D_PY;   //!
   TBranch        *b_D_PZ;   //!
   TBranch        *b_D_MM;   //!
   TBranch        *b_D_MMERR;   //!
   TBranch        *b_D_M;   //!
   TBranch        *b_D_ID;   //!
   TBranch        *b_D_cpx_1_00;   //!
   TBranch        *b_D_cpy_1_00;   //!
   TBranch        *b_D_cpz_1_00;   //!
   TBranch        *b_D_cpt_1_00;   //!
   TBranch        *b_D_cp_1_00;   //!
   TBranch        *b_D_cmult_1_00;   //!
   TBranch        *b_D_pxasy_1_00;   //!
   TBranch        *b_D_pyasy_1_00;   //!
   TBranch        *b_D_pzasy_1_00;   //!
   TBranch        *b_D_pasy_1_00;   //!
   TBranch        *b_D_ptasy_1_00;   //!
   TBranch        *b_D_DOCA1;   //!
   TBranch        *b_D_TAU;   //!
   TBranch        *b_D_TAUERR;   //!
   TBranch        *b_K_D_ETA;   //!
   TBranch        *b_K_D_MC12TuneV2_ProbNNe;   //!
   TBranch        *b_K_D_MC12TuneV2_ProbNNmu;   //!
   TBranch        *b_K_D_MC12TuneV2_ProbNNpi;   //!
   TBranch        *b_K_D_MC12TuneV2_ProbNNk;   //!
   TBranch        *b_K_D_MC12TuneV2_ProbNNp;   //!
   TBranch        *b_K_D_MC12TuneV2_ProbNNghost;   //!
   TBranch        *b_K_D_MC12TuneV3_ProbNNe;   //!
   TBranch        *b_K_D_MC12TuneV3_ProbNNmu;   //!
   TBranch        *b_K_D_MC12TuneV3_ProbNNpi;   //!
   TBranch        *b_K_D_MC12TuneV3_ProbNNk;   //!
   TBranch        *b_K_D_MC12TuneV3_ProbNNp;   //!
   TBranch        *b_K_D_MC12TuneV3_ProbNNghost;   //!
   TBranch        *b_K_D_MC12TuneV4_ProbNNe;   //!
   TBranch        *b_K_D_MC12TuneV4_ProbNNmu;   //!
   TBranch        *b_K_D_MC12TuneV4_ProbNNpi;   //!
   TBranch        *b_K_D_MC12TuneV4_ProbNNk;   //!
   TBranch        *b_K_D_MC12TuneV4_ProbNNp;   //!
   TBranch        *b_K_D_MC12TuneV4_ProbNNghost;   //!
   TBranch        *b_K_D_MC15TuneV1_ProbNNe;   //!
   TBranch        *b_K_D_MC15TuneV1_ProbNNmu;   //!
   TBranch        *b_K_D_MC15TuneV1_ProbNNpi;   //!
   TBranch        *b_K_D_MC15TuneV1_ProbNNk;   //!
   TBranch        *b_K_D_MC15TuneV1_ProbNNp;   //!
   TBranch        *b_K_D_MC15TuneV1_ProbNNghost;   //!
   TBranch        *b_K_D_MINIP;   //!
   TBranch        *b_K_D_MINIPCHI2;   //!
   TBranch        *b_K_D_MINIPNEXTBEST;   //!
   TBranch        *b_K_D_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_K_D_AllIP;   //!
   TBranch        *b_K_D_AllIPchi2;   //!
   TBranch        *b_K_D_OWNPV_X;   //!
   TBranch        *b_K_D_OWNPV_Y;   //!
   TBranch        *b_K_D_OWNPV_Z;   //!
   TBranch        *b_K_D_OWNPV_XERR;   //!
   TBranch        *b_K_D_OWNPV_YERR;   //!
   TBranch        *b_K_D_OWNPV_ZERR;   //!
   TBranch        *b_K_D_OWNPV_CHI2;   //!
   TBranch        *b_K_D_OWNPV_NDOF;   //!
   TBranch        *b_K_D_OWNPV_COV_;   //!
   TBranch        *b_K_D_IP_OWNPV;   //!
   TBranch        *b_K_D_IPCHI2_OWNPV;   //!
   TBranch        *b_K_D_TOPPV_X;   //!
   TBranch        *b_K_D_TOPPV_Y;   //!
   TBranch        *b_K_D_TOPPV_Z;   //!
   TBranch        *b_K_D_TOPPV_XERR;   //!
   TBranch        *b_K_D_TOPPV_YERR;   //!
   TBranch        *b_K_D_TOPPV_ZERR;   //!
   TBranch        *b_K_D_TOPPV_CHI2;   //!
   TBranch        *b_K_D_TOPPV_NDOF;   //!
   TBranch        *b_K_D_TOPPV_COV_;   //!
   TBranch        *b_K_D_IP_TOPPV;   //!
   TBranch        *b_K_D_IPCHI2_TOPPV;   //!
   TBranch        *b_K_D_ORIVX_X;   //!
   TBranch        *b_K_D_ORIVX_Y;   //!
   TBranch        *b_K_D_ORIVX_Z;   //!
   TBranch        *b_K_D_ORIVX_XERR;   //!
   TBranch        *b_K_D_ORIVX_YERR;   //!
   TBranch        *b_K_D_ORIVX_ZERR;   //!
   TBranch        *b_K_D_ORIVX_CHI2;   //!
   TBranch        *b_K_D_ORIVX_NDOF;   //!
   TBranch        *b_K_D_ORIVX_COV_;   //!
   TBranch        *b_K_D_IP_ORIVX;   //!
   TBranch        *b_K_D_IPCHI2_ORIVX;   //!
   TBranch        *b_K_D_P;   //!
   TBranch        *b_K_D_PT;   //!
   TBranch        *b_K_D_PE;   //!
   TBranch        *b_K_D_PX;   //!
   TBranch        *b_K_D_PY;   //!
   TBranch        *b_K_D_PZ;   //!
   TBranch        *b_K_D_M;   //!
   TBranch        *b_K_D_ID;   //!
   TBranch        *b_K_D_PIDe;   //!
   TBranch        *b_K_D_PIDmu;   //!
   TBranch        *b_K_D_PIDK;   //!
   TBranch        *b_K_D_PIDp;   //!
   TBranch        *b_K_D_ProbNNe;   //!
   TBranch        *b_K_D_ProbNNk;   //!
   TBranch        *b_K_D_ProbNNp;   //!
   TBranch        *b_K_D_ProbNNpi;   //!
   TBranch        *b_K_D_ProbNNmu;   //!
   TBranch        *b_K_D_ProbNNghost;   //!
   TBranch        *b_K_D_hasMuon;   //!
   TBranch        *b_K_D_isMuon;   //!
   TBranch        *b_K_D_hasRich;   //!
   TBranch        *b_K_D_UsedRichAerogel;   //!
   TBranch        *b_K_D_UsedRich1Gas;   //!
   TBranch        *b_K_D_UsedRich2Gas;   //!
   TBranch        *b_K_D_RichAboveElThres;   //!
   TBranch        *b_K_D_RichAboveMuThres;   //!
   TBranch        *b_K_D_RichAbovePiThres;   //!
   TBranch        *b_K_D_RichAboveKaThres;   //!
   TBranch        *b_K_D_RichAbovePrThres;   //!
   TBranch        *b_K_D_hasCalo;   //!
   TBranch        *b_K_D_TRACK_Type;   //!
   TBranch        *b_K_D_TRACK_Key;   //!
   TBranch        *b_K_D_TRACK_CHI2NDOF;   //!
   TBranch        *b_K_D_TRACK_PCHI2;   //!
   TBranch        *b_K_D_TRACK_MatchCHI2;   //!
   TBranch        *b_K_D_TRACK_GhostProb;   //!
   TBranch        *b_K_D_TRACK_CloneDist;   //!
   TBranch        *b_K_D_TRACK_Likelihood;   //!
   TBranch        *b_K_D_cpx_1_00;   //!
   TBranch        *b_K_D_cpy_1_00;   //!
   TBranch        *b_K_D_cpz_1_00;   //!
   TBranch        *b_K_D_cpt_1_00;   //!
   TBranch        *b_K_D_cp_1_00;   //!
   TBranch        *b_K_D_cmult_1_00;   //!
   TBranch        *b_K_D_pxasy_1_00;   //!
   TBranch        *b_K_D_pyasy_1_00;   //!
   TBranch        *b_K_D_pzasy_1_00;   //!
   TBranch        *b_K_D_pasy_1_00;   //!
   TBranch        *b_K_D_ptasy_1_00;   //!
   TBranch        *b_pi1_D_ETA;   //!
   TBranch        *b_pi1_D_MC12TuneV2_ProbNNe;   //!
   TBranch        *b_pi1_D_MC12TuneV2_ProbNNmu;   //!
   TBranch        *b_pi1_D_MC12TuneV2_ProbNNpi;   //!
   TBranch        *b_pi1_D_MC12TuneV2_ProbNNk;   //!
   TBranch        *b_pi1_D_MC12TuneV2_ProbNNp;   //!
   TBranch        *b_pi1_D_MC12TuneV2_ProbNNghost;   //!
   TBranch        *b_pi1_D_MC12TuneV3_ProbNNe;   //!
   TBranch        *b_pi1_D_MC12TuneV3_ProbNNmu;   //!
   TBranch        *b_pi1_D_MC12TuneV3_ProbNNpi;   //!
   TBranch        *b_pi1_D_MC12TuneV3_ProbNNk;   //!
   TBranch        *b_pi1_D_MC12TuneV3_ProbNNp;   //!
   TBranch        *b_pi1_D_MC12TuneV3_ProbNNghost;   //!
   TBranch        *b_pi1_D_MC12TuneV4_ProbNNe;   //!
   TBranch        *b_pi1_D_MC12TuneV4_ProbNNmu;   //!
   TBranch        *b_pi1_D_MC12TuneV4_ProbNNpi;   //!
   TBranch        *b_pi1_D_MC12TuneV4_ProbNNk;   //!
   TBranch        *b_pi1_D_MC12TuneV4_ProbNNp;   //!
   TBranch        *b_pi1_D_MC12TuneV4_ProbNNghost;   //!
   TBranch        *b_pi1_D_MC15TuneV1_ProbNNe;   //!
   TBranch        *b_pi1_D_MC15TuneV1_ProbNNmu;   //!
   TBranch        *b_pi1_D_MC15TuneV1_ProbNNpi;   //!
   TBranch        *b_pi1_D_MC15TuneV1_ProbNNk;   //!
   TBranch        *b_pi1_D_MC15TuneV1_ProbNNp;   //!
   TBranch        *b_pi1_D_MC15TuneV1_ProbNNghost;   //!
   TBranch        *b_pi1_D_MINIP;   //!
   TBranch        *b_pi1_D_MINIPCHI2;   //!
   TBranch        *b_pi1_D_MINIPNEXTBEST;   //!
   TBranch        *b_pi1_D_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_pi1_D_AllIP;   //!
   TBranch        *b_pi1_D_AllIPchi2;   //!
   TBranch        *b_pi1_D_OWNPV_X;   //!
   TBranch        *b_pi1_D_OWNPV_Y;   //!
   TBranch        *b_pi1_D_OWNPV_Z;   //!
   TBranch        *b_pi1_D_OWNPV_XERR;   //!
   TBranch        *b_pi1_D_OWNPV_YERR;   //!
   TBranch        *b_pi1_D_OWNPV_ZERR;   //!
   TBranch        *b_pi1_D_OWNPV_CHI2;   //!
   TBranch        *b_pi1_D_OWNPV_NDOF;   //!
   TBranch        *b_pi1_D_OWNPV_COV_;   //!
   TBranch        *b_pi1_D_IP_OWNPV;   //!
   TBranch        *b_pi1_D_IPCHI2_OWNPV;   //!
   TBranch        *b_pi1_D_TOPPV_X;   //!
   TBranch        *b_pi1_D_TOPPV_Y;   //!
   TBranch        *b_pi1_D_TOPPV_Z;   //!
   TBranch        *b_pi1_D_TOPPV_XERR;   //!
   TBranch        *b_pi1_D_TOPPV_YERR;   //!
   TBranch        *b_pi1_D_TOPPV_ZERR;   //!
   TBranch        *b_pi1_D_TOPPV_CHI2;   //!
   TBranch        *b_pi1_D_TOPPV_NDOF;   //!
   TBranch        *b_pi1_D_TOPPV_COV_;   //!
   TBranch        *b_pi1_D_IP_TOPPV;   //!
   TBranch        *b_pi1_D_IPCHI2_TOPPV;   //!
   TBranch        *b_pi1_D_ORIVX_X;   //!
   TBranch        *b_pi1_D_ORIVX_Y;   //!
   TBranch        *b_pi1_D_ORIVX_Z;   //!
   TBranch        *b_pi1_D_ORIVX_XERR;   //!
   TBranch        *b_pi1_D_ORIVX_YERR;   //!
   TBranch        *b_pi1_D_ORIVX_ZERR;   //!
   TBranch        *b_pi1_D_ORIVX_CHI2;   //!
   TBranch        *b_pi1_D_ORIVX_NDOF;   //!
   TBranch        *b_pi1_D_ORIVX_COV_;   //!
   TBranch        *b_pi1_D_IP_ORIVX;   //!
   TBranch        *b_pi1_D_IPCHI2_ORIVX;   //!
   TBranch        *b_pi1_D_P;   //!
   TBranch        *b_pi1_D_PT;   //!
   TBranch        *b_pi1_D_PE;   //!
   TBranch        *b_pi1_D_PX;   //!
   TBranch        *b_pi1_D_PY;   //!
   TBranch        *b_pi1_D_PZ;   //!
   TBranch        *b_pi1_D_M;   //!
   TBranch        *b_pi1_D_ID;   //!
   TBranch        *b_pi1_D_PIDe;   //!
   TBranch        *b_pi1_D_PIDmu;   //!
   TBranch        *b_pi1_D_PIDK;   //!
   TBranch        *b_pi1_D_PIDp;   //!
   TBranch        *b_pi1_D_ProbNNe;   //!
   TBranch        *b_pi1_D_ProbNNk;   //!
   TBranch        *b_pi1_D_ProbNNp;   //!
   TBranch        *b_pi1_D_ProbNNpi;   //!
   TBranch        *b_pi1_D_ProbNNmu;   //!
   TBranch        *b_pi1_D_ProbNNghost;   //!
   TBranch        *b_pi1_D_hasMuon;   //!
   TBranch        *b_pi1_D_isMuon;   //!
   TBranch        *b_pi1_D_hasRich;   //!
   TBranch        *b_pi1_D_UsedRichAerogel;   //!
   TBranch        *b_pi1_D_UsedRich1Gas;   //!
   TBranch        *b_pi1_D_UsedRich2Gas;   //!
   TBranch        *b_pi1_D_RichAboveElThres;   //!
   TBranch        *b_pi1_D_RichAboveMuThres;   //!
   TBranch        *b_pi1_D_RichAbovePiThres;   //!
   TBranch        *b_pi1_D_RichAboveKaThres;   //!
   TBranch        *b_pi1_D_RichAbovePrThres;   //!
   TBranch        *b_pi1_D_hasCalo;   //!
   TBranch        *b_pi1_D_TRACK_Type;   //!
   TBranch        *b_pi1_D_TRACK_Key;   //!
   TBranch        *b_pi1_D_TRACK_CHI2NDOF;   //!
   TBranch        *b_pi1_D_TRACK_PCHI2;   //!
   TBranch        *b_pi1_D_TRACK_MatchCHI2;   //!
   TBranch        *b_pi1_D_TRACK_GhostProb;   //!
   TBranch        *b_pi1_D_TRACK_CloneDist;   //!
   TBranch        *b_pi1_D_TRACK_Likelihood;   //!
   TBranch        *b_pi1_D_cpx_1_00;   //!
   TBranch        *b_pi1_D_cpy_1_00;   //!
   TBranch        *b_pi1_D_cpz_1_00;   //!
   TBranch        *b_pi1_D_cpt_1_00;   //!
   TBranch        *b_pi1_D_cp_1_00;   //!
   TBranch        *b_pi1_D_cmult_1_00;   //!
   TBranch        *b_pi1_D_pxasy_1_00;   //!
   TBranch        *b_pi1_D_pyasy_1_00;   //!
   TBranch        *b_pi1_D_pzasy_1_00;   //!
   TBranch        *b_pi1_D_pasy_1_00;   //!
   TBranch        *b_pi1_D_ptasy_1_00;   //!
   TBranch        *b_pi2_D_ETA;   //!
   TBranch        *b_pi2_D_MC12TuneV2_ProbNNe;   //!
   TBranch        *b_pi2_D_MC12TuneV2_ProbNNmu;   //!
   TBranch        *b_pi2_D_MC12TuneV2_ProbNNpi;   //!
   TBranch        *b_pi2_D_MC12TuneV2_ProbNNk;   //!
   TBranch        *b_pi2_D_MC12TuneV2_ProbNNp;   //!
   TBranch        *b_pi2_D_MC12TuneV2_ProbNNghost;   //!
   TBranch        *b_pi2_D_MC12TuneV3_ProbNNe;   //!
   TBranch        *b_pi2_D_MC12TuneV3_ProbNNmu;   //!
   TBranch        *b_pi2_D_MC12TuneV3_ProbNNpi;   //!
   TBranch        *b_pi2_D_MC12TuneV3_ProbNNk;   //!
   TBranch        *b_pi2_D_MC12TuneV3_ProbNNp;   //!
   TBranch        *b_pi2_D_MC12TuneV3_ProbNNghost;   //!
   TBranch        *b_pi2_D_MC12TuneV4_ProbNNe;   //!
   TBranch        *b_pi2_D_MC12TuneV4_ProbNNmu;   //!
   TBranch        *b_pi2_D_MC12TuneV4_ProbNNpi;   //!
   TBranch        *b_pi2_D_MC12TuneV4_ProbNNk;   //!
   TBranch        *b_pi2_D_MC12TuneV4_ProbNNp;   //!
   TBranch        *b_pi2_D_MC12TuneV4_ProbNNghost;   //!
   TBranch        *b_pi2_D_MC15TuneV1_ProbNNe;   //!
   TBranch        *b_pi2_D_MC15TuneV1_ProbNNmu;   //!
   TBranch        *b_pi2_D_MC15TuneV1_ProbNNpi;   //!
   TBranch        *b_pi2_D_MC15TuneV1_ProbNNk;   //!
   TBranch        *b_pi2_D_MC15TuneV1_ProbNNp;   //!
   TBranch        *b_pi2_D_MC15TuneV1_ProbNNghost;   //!
   TBranch        *b_pi2_D_MINIP;   //!
   TBranch        *b_pi2_D_MINIPCHI2;   //!
   TBranch        *b_pi2_D_MINIPNEXTBEST;   //!
   TBranch        *b_pi2_D_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_pi2_D_AllIP;   //!
   TBranch        *b_pi2_D_AllIPchi2;   //!
   TBranch        *b_pi2_D_OWNPV_X;   //!
   TBranch        *b_pi2_D_OWNPV_Y;   //!
   TBranch        *b_pi2_D_OWNPV_Z;   //!
   TBranch        *b_pi2_D_OWNPV_XERR;   //!
   TBranch        *b_pi2_D_OWNPV_YERR;   //!
   TBranch        *b_pi2_D_OWNPV_ZERR;   //!
   TBranch        *b_pi2_D_OWNPV_CHI2;   //!
   TBranch        *b_pi2_D_OWNPV_NDOF;   //!
   TBranch        *b_pi2_D_OWNPV_COV_;   //!
   TBranch        *b_pi2_D_IP_OWNPV;   //!
   TBranch        *b_pi2_D_IPCHI2_OWNPV;   //!
   TBranch        *b_pi2_D_TOPPV_X;   //!
   TBranch        *b_pi2_D_TOPPV_Y;   //!
   TBranch        *b_pi2_D_TOPPV_Z;   //!
   TBranch        *b_pi2_D_TOPPV_XERR;   //!
   TBranch        *b_pi2_D_TOPPV_YERR;   //!
   TBranch        *b_pi2_D_TOPPV_ZERR;   //!
   TBranch        *b_pi2_D_TOPPV_CHI2;   //!
   TBranch        *b_pi2_D_TOPPV_NDOF;   //!
   TBranch        *b_pi2_D_TOPPV_COV_;   //!
   TBranch        *b_pi2_D_IP_TOPPV;   //!
   TBranch        *b_pi2_D_IPCHI2_TOPPV;   //!
   TBranch        *b_pi2_D_ORIVX_X;   //!
   TBranch        *b_pi2_D_ORIVX_Y;   //!
   TBranch        *b_pi2_D_ORIVX_Z;   //!
   TBranch        *b_pi2_D_ORIVX_XERR;   //!
   TBranch        *b_pi2_D_ORIVX_YERR;   //!
   TBranch        *b_pi2_D_ORIVX_ZERR;   //!
   TBranch        *b_pi2_D_ORIVX_CHI2;   //!
   TBranch        *b_pi2_D_ORIVX_NDOF;   //!
   TBranch        *b_pi2_D_ORIVX_COV_;   //!
   TBranch        *b_pi2_D_IP_ORIVX;   //!
   TBranch        *b_pi2_D_IPCHI2_ORIVX;   //!
   TBranch        *b_pi2_D_P;   //!
   TBranch        *b_pi2_D_PT;   //!
   TBranch        *b_pi2_D_PE;   //!
   TBranch        *b_pi2_D_PX;   //!
   TBranch        *b_pi2_D_PY;   //!
   TBranch        *b_pi2_D_PZ;   //!
   TBranch        *b_pi2_D_M;   //!
   TBranch        *b_pi2_D_ID;   //!
   TBranch        *b_pi2_D_PIDe;   //!
   TBranch        *b_pi2_D_PIDmu;   //!
   TBranch        *b_pi2_D_PIDK;   //!
   TBranch        *b_pi2_D_PIDp;   //!
   TBranch        *b_pi2_D_ProbNNe;   //!
   TBranch        *b_pi2_D_ProbNNk;   //!
   TBranch        *b_pi2_D_ProbNNp;   //!
   TBranch        *b_pi2_D_ProbNNpi;   //!
   TBranch        *b_pi2_D_ProbNNmu;   //!
   TBranch        *b_pi2_D_ProbNNghost;   //!
   TBranch        *b_pi2_D_hasMuon;   //!
   TBranch        *b_pi2_D_isMuon;   //!
   TBranch        *b_pi2_D_hasRich;   //!
   TBranch        *b_pi2_D_UsedRichAerogel;   //!
   TBranch        *b_pi2_D_UsedRich1Gas;   //!
   TBranch        *b_pi2_D_UsedRich2Gas;   //!
   TBranch        *b_pi2_D_RichAboveElThres;   //!
   TBranch        *b_pi2_D_RichAboveMuThres;   //!
   TBranch        *b_pi2_D_RichAbovePiThres;   //!
   TBranch        *b_pi2_D_RichAboveKaThres;   //!
   TBranch        *b_pi2_D_RichAbovePrThres;   //!
   TBranch        *b_pi2_D_hasCalo;   //!
   TBranch        *b_pi2_D_TRACK_Type;   //!
   TBranch        *b_pi2_D_TRACK_Key;   //!
   TBranch        *b_pi2_D_TRACK_CHI2NDOF;   //!
   TBranch        *b_pi2_D_TRACK_PCHI2;   //!
   TBranch        *b_pi2_D_TRACK_MatchCHI2;   //!
   TBranch        *b_pi2_D_TRACK_GhostProb;   //!
   TBranch        *b_pi2_D_TRACK_CloneDist;   //!
   TBranch        *b_pi2_D_TRACK_Likelihood;   //!
   TBranch        *b_pi2_D_cpx_1_00;   //!
   TBranch        *b_pi2_D_cpy_1_00;   //!
   TBranch        *b_pi2_D_cpz_1_00;   //!
   TBranch        *b_pi2_D_cpt_1_00;   //!
   TBranch        *b_pi2_D_cp_1_00;   //!
   TBranch        *b_pi2_D_cmult_1_00;   //!
   TBranch        *b_pi2_D_pxasy_1_00;   //!
   TBranch        *b_pi2_D_pyasy_1_00;   //!
   TBranch        *b_pi2_D_pzasy_1_00;   //!
   TBranch        *b_pi2_D_pasy_1_00;   //!
   TBranch        *b_pi2_D_ptasy_1_00;   //!
   TBranch        *b_Ks_ETA;   //!
   TBranch        *b_Ks_MINIP;   //!
   TBranch        *b_Ks_MINIPCHI2;   //!
   TBranch        *b_Ks_MINIPNEXTBEST;   //!
   TBranch        *b_Ks_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_Ks_AllIP;   //!
   TBranch        *b_Ks_AllIPchi2;   //!
   TBranch        *b_Ks_AllDIRA;   //!
   TBranch        *b_Ks_ENDVERTEX_X;   //!
   TBranch        *b_Ks_ENDVERTEX_Y;   //!
   TBranch        *b_Ks_ENDVERTEX_Z;   //!
   TBranch        *b_Ks_ENDVERTEX_XERR;   //!
   TBranch        *b_Ks_ENDVERTEX_YERR;   //!
   TBranch        *b_Ks_ENDVERTEX_ZERR;   //!
   TBranch        *b_Ks_ENDVERTEX_CHI2;   //!
   TBranch        *b_Ks_ENDVERTEX_NDOF;   //!
   TBranch        *b_Ks_ENDVERTEX_COV_;   //!
   TBranch        *b_Ks_OWNPV_X;   //!
   TBranch        *b_Ks_OWNPV_Y;   //!
   TBranch        *b_Ks_OWNPV_Z;   //!
   TBranch        *b_Ks_OWNPV_XERR;   //!
   TBranch        *b_Ks_OWNPV_YERR;   //!
   TBranch        *b_Ks_OWNPV_ZERR;   //!
   TBranch        *b_Ks_OWNPV_CHI2;   //!
   TBranch        *b_Ks_OWNPV_NDOF;   //!
   TBranch        *b_Ks_OWNPV_COV_;   //!
   TBranch        *b_Ks_IP_OWNPV;   //!
   TBranch        *b_Ks_IPCHI2_OWNPV;   //!
   TBranch        *b_Ks_FD_OWNPV;   //!
   TBranch        *b_Ks_FDCHI2_OWNPV;   //!
   TBranch        *b_Ks_DIRA_OWNPV;   //!
   TBranch        *b_Ks_TOPPV_X;   //!
   TBranch        *b_Ks_TOPPV_Y;   //!
   TBranch        *b_Ks_TOPPV_Z;   //!
   TBranch        *b_Ks_TOPPV_XERR;   //!
   TBranch        *b_Ks_TOPPV_YERR;   //!
   TBranch        *b_Ks_TOPPV_ZERR;   //!
   TBranch        *b_Ks_TOPPV_CHI2;   //!
   TBranch        *b_Ks_TOPPV_NDOF;   //!
   TBranch        *b_Ks_TOPPV_COV_;   //!
   TBranch        *b_Ks_IP_TOPPV;   //!
   TBranch        *b_Ks_IPCHI2_TOPPV;   //!
   TBranch        *b_Ks_FD_TOPPV;   //!
   TBranch        *b_Ks_FDCHI2_TOPPV;   //!
   TBranch        *b_Ks_DIRA_TOPPV;   //!
   TBranch        *b_Ks_ORIVX_X;   //!
   TBranch        *b_Ks_ORIVX_Y;   //!
   TBranch        *b_Ks_ORIVX_Z;   //!
   TBranch        *b_Ks_ORIVX_XERR;   //!
   TBranch        *b_Ks_ORIVX_YERR;   //!
   TBranch        *b_Ks_ORIVX_ZERR;   //!
   TBranch        *b_Ks_ORIVX_CHI2;   //!
   TBranch        *b_Ks_ORIVX_NDOF;   //!
   TBranch        *b_Ks_ORIVX_COV_;   //!
   TBranch        *b_Ks_IP_ORIVX;   //!
   TBranch        *b_Ks_IPCHI2_ORIVX;   //!
   TBranch        *b_Ks_FD_ORIVX;   //!
   TBranch        *b_Ks_FDCHI2_ORIVX;   //!
   TBranch        *b_Ks_DIRA_ORIVX;   //!
   TBranch        *b_Ks_P;   //!
   TBranch        *b_Ks_PT;   //!
   TBranch        *b_Ks_PE;   //!
   TBranch        *b_Ks_PX;   //!
   TBranch        *b_Ks_PY;   //!
   TBranch        *b_Ks_PZ;   //!
   TBranch        *b_Ks_MM;   //!
   TBranch        *b_Ks_MMERR;   //!
   TBranch        *b_Ks_M;   //!
   TBranch        *b_Ks_ID;   //!
   TBranch        *b_Ks_cpx_1_00;   //!
   TBranch        *b_Ks_cpy_1_00;   //!
   TBranch        *b_Ks_cpz_1_00;   //!
   TBranch        *b_Ks_cpt_1_00;   //!
   TBranch        *b_Ks_cp_1_00;   //!
   TBranch        *b_Ks_cmult_1_00;   //!
   TBranch        *b_Ks_pxasy_1_00;   //!
   TBranch        *b_Ks_pyasy_1_00;   //!
   TBranch        *b_Ks_pzasy_1_00;   //!
   TBranch        *b_Ks_pasy_1_00;   //!
   TBranch        *b_Ks_ptasy_1_00;   //!
   TBranch        *b_Ks_DOCA1;   //!
   TBranch        *b_Ks_TAU;   //!
   TBranch        *b_Ks_TAUERR;   //!
   TBranch        *b_pip_Ks_ETA;   //!
   TBranch        *b_pip_Ks_MC12TuneV2_ProbNNe;   //!
   TBranch        *b_pip_Ks_MC12TuneV2_ProbNNmu;   //!
   TBranch        *b_pip_Ks_MC12TuneV2_ProbNNpi;   //!
   TBranch        *b_pip_Ks_MC12TuneV2_ProbNNk;   //!
   TBranch        *b_pip_Ks_MC12TuneV2_ProbNNp;   //!
   TBranch        *b_pip_Ks_MC12TuneV2_ProbNNghost;   //!
   TBranch        *b_pip_Ks_MC12TuneV3_ProbNNe;   //!
   TBranch        *b_pip_Ks_MC12TuneV3_ProbNNmu;   //!
   TBranch        *b_pip_Ks_MC12TuneV3_ProbNNpi;   //!
   TBranch        *b_pip_Ks_MC12TuneV3_ProbNNk;   //!
   TBranch        *b_pip_Ks_MC12TuneV3_ProbNNp;   //!
   TBranch        *b_pip_Ks_MC12TuneV3_ProbNNghost;   //!
   TBranch        *b_pip_Ks_MC12TuneV4_ProbNNe;   //!
   TBranch        *b_pip_Ks_MC12TuneV4_ProbNNmu;   //!
   TBranch        *b_pip_Ks_MC12TuneV4_ProbNNpi;   //!
   TBranch        *b_pip_Ks_MC12TuneV4_ProbNNk;   //!
   TBranch        *b_pip_Ks_MC12TuneV4_ProbNNp;   //!
   TBranch        *b_pip_Ks_MC12TuneV4_ProbNNghost;   //!
   TBranch        *b_pip_Ks_MC15TuneV1_ProbNNe;   //!
   TBranch        *b_pip_Ks_MC15TuneV1_ProbNNmu;   //!
   TBranch        *b_pip_Ks_MC15TuneV1_ProbNNpi;   //!
   TBranch        *b_pip_Ks_MC15TuneV1_ProbNNk;   //!
   TBranch        *b_pip_Ks_MC15TuneV1_ProbNNp;   //!
   TBranch        *b_pip_Ks_MC15TuneV1_ProbNNghost;   //!
   TBranch        *b_pip_Ks_MINIP;   //!
   TBranch        *b_pip_Ks_MINIPCHI2;   //!
   TBranch        *b_pip_Ks_MINIPNEXTBEST;   //!
   TBranch        *b_pip_Ks_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_pip_Ks_AllIP;   //!
   TBranch        *b_pip_Ks_AllIPchi2;   //!
   TBranch        *b_pip_Ks_OWNPV_X;   //!
   TBranch        *b_pip_Ks_OWNPV_Y;   //!
   TBranch        *b_pip_Ks_OWNPV_Z;   //!
   TBranch        *b_pip_Ks_OWNPV_XERR;   //!
   TBranch        *b_pip_Ks_OWNPV_YERR;   //!
   TBranch        *b_pip_Ks_OWNPV_ZERR;   //!
   TBranch        *b_pip_Ks_OWNPV_CHI2;   //!
   TBranch        *b_pip_Ks_OWNPV_NDOF;   //!
   TBranch        *b_pip_Ks_OWNPV_COV_;   //!
   TBranch        *b_pip_Ks_IP_OWNPV;   //!
   TBranch        *b_pip_Ks_IPCHI2_OWNPV;   //!
   TBranch        *b_pip_Ks_TOPPV_X;   //!
   TBranch        *b_pip_Ks_TOPPV_Y;   //!
   TBranch        *b_pip_Ks_TOPPV_Z;   //!
   TBranch        *b_pip_Ks_TOPPV_XERR;   //!
   TBranch        *b_pip_Ks_TOPPV_YERR;   //!
   TBranch        *b_pip_Ks_TOPPV_ZERR;   //!
   TBranch        *b_pip_Ks_TOPPV_CHI2;   //!
   TBranch        *b_pip_Ks_TOPPV_NDOF;   //!
   TBranch        *b_pip_Ks_TOPPV_COV_;   //!
   TBranch        *b_pip_Ks_IP_TOPPV;   //!
   TBranch        *b_pip_Ks_IPCHI2_TOPPV;   //!
   TBranch        *b_pip_Ks_ORIVX_X;   //!
   TBranch        *b_pip_Ks_ORIVX_Y;   //!
   TBranch        *b_pip_Ks_ORIVX_Z;   //!
   TBranch        *b_pip_Ks_ORIVX_XERR;   //!
   TBranch        *b_pip_Ks_ORIVX_YERR;   //!
   TBranch        *b_pip_Ks_ORIVX_ZERR;   //!
   TBranch        *b_pip_Ks_ORIVX_CHI2;   //!
   TBranch        *b_pip_Ks_ORIVX_NDOF;   //!
   TBranch        *b_pip_Ks_ORIVX_COV_;   //!
   TBranch        *b_pip_Ks_IP_ORIVX;   //!
   TBranch        *b_pip_Ks_IPCHI2_ORIVX;   //!
   TBranch        *b_pip_Ks_P;   //!
   TBranch        *b_pip_Ks_PT;   //!
   TBranch        *b_pip_Ks_PE;   //!
   TBranch        *b_pip_Ks_PX;   //!
   TBranch        *b_pip_Ks_PY;   //!
   TBranch        *b_pip_Ks_PZ;   //!
   TBranch        *b_pip_Ks_M;   //!
   TBranch        *b_pip_Ks_ID;   //!
   TBranch        *b_pip_Ks_PIDe;   //!
   TBranch        *b_pip_Ks_PIDmu;   //!
   TBranch        *b_pip_Ks_PIDK;   //!
   TBranch        *b_pip_Ks_PIDp;   //!
   TBranch        *b_pip_Ks_ProbNNe;   //!
   TBranch        *b_pip_Ks_ProbNNk;   //!
   TBranch        *b_pip_Ks_ProbNNp;   //!
   TBranch        *b_pip_Ks_ProbNNpi;   //!
   TBranch        *b_pip_Ks_ProbNNmu;   //!
   TBranch        *b_pip_Ks_ProbNNghost;   //!
   TBranch        *b_pip_Ks_hasMuon;   //!
   TBranch        *b_pip_Ks_isMuon;   //!
   TBranch        *b_pip_Ks_hasRich;   //!
   TBranch        *b_pip_Ks_UsedRichAerogel;   //!
   TBranch        *b_pip_Ks_UsedRich1Gas;   //!
   TBranch        *b_pip_Ks_UsedRich2Gas;   //!
   TBranch        *b_pip_Ks_RichAboveElThres;   //!
   TBranch        *b_pip_Ks_RichAboveMuThres;   //!
   TBranch        *b_pip_Ks_RichAbovePiThres;   //!
   TBranch        *b_pip_Ks_RichAboveKaThres;   //!
   TBranch        *b_pip_Ks_RichAbovePrThres;   //!
   TBranch        *b_pip_Ks_hasCalo;   //!
   TBranch        *b_pip_Ks_TRACK_Type;   //!
   TBranch        *b_pip_Ks_TRACK_Key;   //!
   TBranch        *b_pip_Ks_TRACK_CHI2NDOF;   //!
   TBranch        *b_pip_Ks_TRACK_PCHI2;   //!
   TBranch        *b_pip_Ks_TRACK_MatchCHI2;   //!
   TBranch        *b_pip_Ks_TRACK_GhostProb;   //!
   TBranch        *b_pip_Ks_TRACK_CloneDist;   //!
   TBranch        *b_pip_Ks_TRACK_Likelihood;   //!
   TBranch        *b_pip_Ks_cpx_1_00;   //!
   TBranch        *b_pip_Ks_cpy_1_00;   //!
   TBranch        *b_pip_Ks_cpz_1_00;   //!
   TBranch        *b_pip_Ks_cpt_1_00;   //!
   TBranch        *b_pip_Ks_cp_1_00;   //!
   TBranch        *b_pip_Ks_cmult_1_00;   //!
   TBranch        *b_pip_Ks_pxasy_1_00;   //!
   TBranch        *b_pip_Ks_pyasy_1_00;   //!
   TBranch        *b_pip_Ks_pzasy_1_00;   //!
   TBranch        *b_pip_Ks_pasy_1_00;   //!
   TBranch        *b_pip_Ks_ptasy_1_00;   //!
   TBranch        *b_pim_Ks_ETA;   //!
   TBranch        *b_pim_Ks_MC12TuneV2_ProbNNe;   //!
   TBranch        *b_pim_Ks_MC12TuneV2_ProbNNmu;   //!
   TBranch        *b_pim_Ks_MC12TuneV2_ProbNNpi;   //!
   TBranch        *b_pim_Ks_MC12TuneV2_ProbNNk;   //!
   TBranch        *b_pim_Ks_MC12TuneV2_ProbNNp;   //!
   TBranch        *b_pim_Ks_MC12TuneV2_ProbNNghost;   //!
   TBranch        *b_pim_Ks_MC12TuneV3_ProbNNe;   //!
   TBranch        *b_pim_Ks_MC12TuneV3_ProbNNmu;   //!
   TBranch        *b_pim_Ks_MC12TuneV3_ProbNNpi;   //!
   TBranch        *b_pim_Ks_MC12TuneV3_ProbNNk;   //!
   TBranch        *b_pim_Ks_MC12TuneV3_ProbNNp;   //!
   TBranch        *b_pim_Ks_MC12TuneV3_ProbNNghost;   //!
   TBranch        *b_pim_Ks_MC12TuneV4_ProbNNe;   //!
   TBranch        *b_pim_Ks_MC12TuneV4_ProbNNmu;   //!
   TBranch        *b_pim_Ks_MC12TuneV4_ProbNNpi;   //!
   TBranch        *b_pim_Ks_MC12TuneV4_ProbNNk;   //!
   TBranch        *b_pim_Ks_MC12TuneV4_ProbNNp;   //!
   TBranch        *b_pim_Ks_MC12TuneV4_ProbNNghost;   //!
   TBranch        *b_pim_Ks_MC15TuneV1_ProbNNe;   //!
   TBranch        *b_pim_Ks_MC15TuneV1_ProbNNmu;   //!
   TBranch        *b_pim_Ks_MC15TuneV1_ProbNNpi;   //!
   TBranch        *b_pim_Ks_MC15TuneV1_ProbNNk;   //!
   TBranch        *b_pim_Ks_MC15TuneV1_ProbNNp;   //!
   TBranch        *b_pim_Ks_MC15TuneV1_ProbNNghost;   //!
   TBranch        *b_pim_Ks_MINIP;   //!
   TBranch        *b_pim_Ks_MINIPCHI2;   //!
   TBranch        *b_pim_Ks_MINIPNEXTBEST;   //!
   TBranch        *b_pim_Ks_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_pim_Ks_AllIP;   //!
   TBranch        *b_pim_Ks_AllIPchi2;   //!
   TBranch        *b_pim_Ks_OWNPV_X;   //!
   TBranch        *b_pim_Ks_OWNPV_Y;   //!
   TBranch        *b_pim_Ks_OWNPV_Z;   //!
   TBranch        *b_pim_Ks_OWNPV_XERR;   //!
   TBranch        *b_pim_Ks_OWNPV_YERR;   //!
   TBranch        *b_pim_Ks_OWNPV_ZERR;   //!
   TBranch        *b_pim_Ks_OWNPV_CHI2;   //!
   TBranch        *b_pim_Ks_OWNPV_NDOF;   //!
   TBranch        *b_pim_Ks_OWNPV_COV_;   //!
   TBranch        *b_pim_Ks_IP_OWNPV;   //!
   TBranch        *b_pim_Ks_IPCHI2_OWNPV;   //!
   TBranch        *b_pim_Ks_TOPPV_X;   //!
   TBranch        *b_pim_Ks_TOPPV_Y;   //!
   TBranch        *b_pim_Ks_TOPPV_Z;   //!
   TBranch        *b_pim_Ks_TOPPV_XERR;   //!
   TBranch        *b_pim_Ks_TOPPV_YERR;   //!
   TBranch        *b_pim_Ks_TOPPV_ZERR;   //!
   TBranch        *b_pim_Ks_TOPPV_CHI2;   //!
   TBranch        *b_pim_Ks_TOPPV_NDOF;   //!
   TBranch        *b_pim_Ks_TOPPV_COV_;   //!
   TBranch        *b_pim_Ks_IP_TOPPV;   //!
   TBranch        *b_pim_Ks_IPCHI2_TOPPV;   //!
   TBranch        *b_pim_Ks_ORIVX_X;   //!
   TBranch        *b_pim_Ks_ORIVX_Y;   //!
   TBranch        *b_pim_Ks_ORIVX_Z;   //!
   TBranch        *b_pim_Ks_ORIVX_XERR;   //!
   TBranch        *b_pim_Ks_ORIVX_YERR;   //!
   TBranch        *b_pim_Ks_ORIVX_ZERR;   //!
   TBranch        *b_pim_Ks_ORIVX_CHI2;   //!
   TBranch        *b_pim_Ks_ORIVX_NDOF;   //!
   TBranch        *b_pim_Ks_ORIVX_COV_;   //!
   TBranch        *b_pim_Ks_IP_ORIVX;   //!
   TBranch        *b_pim_Ks_IPCHI2_ORIVX;   //!
   TBranch        *b_pim_Ks_P;   //!
   TBranch        *b_pim_Ks_PT;   //!
   TBranch        *b_pim_Ks_PE;   //!
   TBranch        *b_pim_Ks_PX;   //!
   TBranch        *b_pim_Ks_PY;   //!
   TBranch        *b_pim_Ks_PZ;   //!
   TBranch        *b_pim_Ks_M;   //!
   TBranch        *b_pim_Ks_ID;   //!
   TBranch        *b_pim_Ks_PIDe;   //!
   TBranch        *b_pim_Ks_PIDmu;   //!
   TBranch        *b_pim_Ks_PIDK;   //!
   TBranch        *b_pim_Ks_PIDp;   //!
   TBranch        *b_pim_Ks_ProbNNe;   //!
   TBranch        *b_pim_Ks_ProbNNk;   //!
   TBranch        *b_pim_Ks_ProbNNp;   //!
   TBranch        *b_pim_Ks_ProbNNpi;   //!
   TBranch        *b_pim_Ks_ProbNNmu;   //!
   TBranch        *b_pim_Ks_ProbNNghost;   //!
   TBranch        *b_pim_Ks_hasMuon;   //!
   TBranch        *b_pim_Ks_isMuon;   //!
   TBranch        *b_pim_Ks_hasRich;   //!
   TBranch        *b_pim_Ks_UsedRichAerogel;   //!
   TBranch        *b_pim_Ks_UsedRich1Gas;   //!
   TBranch        *b_pim_Ks_UsedRich2Gas;   //!
   TBranch        *b_pim_Ks_RichAboveElThres;   //!
   TBranch        *b_pim_Ks_RichAboveMuThres;   //!
   TBranch        *b_pim_Ks_RichAbovePiThres;   //!
   TBranch        *b_pim_Ks_RichAboveKaThres;   //!
   TBranch        *b_pim_Ks_RichAbovePrThres;   //!
   TBranch        *b_pim_Ks_hasCalo;   //!
   TBranch        *b_pim_Ks_TRACK_Type;   //!
   TBranch        *b_pim_Ks_TRACK_Key;   //!
   TBranch        *b_pim_Ks_TRACK_CHI2NDOF;   //!
   TBranch        *b_pim_Ks_TRACK_PCHI2;   //!
   TBranch        *b_pim_Ks_TRACK_MatchCHI2;   //!
   TBranch        *b_pim_Ks_TRACK_GhostProb;   //!
   TBranch        *b_pim_Ks_TRACK_CloneDist;   //!
   TBranch        *b_pim_Ks_TRACK_Likelihood;   //!
   TBranch        *b_pim_Ks_cpx_1_00;   //!
   TBranch        *b_pim_Ks_cpy_1_00;   //!
   TBranch        *b_pim_Ks_cpz_1_00;   //!
   TBranch        *b_pim_Ks_cpt_1_00;   //!
   TBranch        *b_pim_Ks_cp_1_00;   //!
   TBranch        *b_pim_Ks_cmult_1_00;   //!
   TBranch        *b_pim_Ks_pxasy_1_00;   //!
   TBranch        *b_pim_Ks_pyasy_1_00;   //!
   TBranch        *b_pim_Ks_pzasy_1_00;   //!
   TBranch        *b_pim_Ks_pasy_1_00;   //!
   TBranch        *b_pim_Ks_ptasy_1_00;   //!
   TBranch        *b_pi_ETA;   //!
   TBranch        *b_pi_MC12TuneV2_ProbNNe;   //!
   TBranch        *b_pi_MC12TuneV2_ProbNNmu;   //!
   TBranch        *b_pi_MC12TuneV2_ProbNNpi;   //!
   TBranch        *b_pi_MC12TuneV2_ProbNNk;   //!
   TBranch        *b_pi_MC12TuneV2_ProbNNp;   //!
   TBranch        *b_pi_MC12TuneV2_ProbNNghost;   //!
   TBranch        *b_pi_MC12TuneV3_ProbNNe;   //!
   TBranch        *b_pi_MC12TuneV3_ProbNNmu;   //!
   TBranch        *b_pi_MC12TuneV3_ProbNNpi;   //!
   TBranch        *b_pi_MC12TuneV3_ProbNNk;   //!
   TBranch        *b_pi_MC12TuneV3_ProbNNp;   //!
   TBranch        *b_pi_MC12TuneV3_ProbNNghost;   //!
   TBranch        *b_pi_MC12TuneV4_ProbNNe;   //!
   TBranch        *b_pi_MC12TuneV4_ProbNNmu;   //!
   TBranch        *b_pi_MC12TuneV4_ProbNNpi;   //!
   TBranch        *b_pi_MC12TuneV4_ProbNNk;   //!
   TBranch        *b_pi_MC12TuneV4_ProbNNp;   //!
   TBranch        *b_pi_MC12TuneV4_ProbNNghost;   //!
   TBranch        *b_pi_MC15TuneV1_ProbNNe;   //!
   TBranch        *b_pi_MC15TuneV1_ProbNNmu;   //!
   TBranch        *b_pi_MC15TuneV1_ProbNNpi;   //!
   TBranch        *b_pi_MC15TuneV1_ProbNNk;   //!
   TBranch        *b_pi_MC15TuneV1_ProbNNp;   //!
   TBranch        *b_pi_MC15TuneV1_ProbNNghost;   //!
   TBranch        *b_pi_MINIP;   //!
   TBranch        *b_pi_MINIPCHI2;   //!
   TBranch        *b_pi_MINIPNEXTBEST;   //!
   TBranch        *b_pi_MINIPCHI2NEXTBEST;   //!
   TBranch        *b_pi_AllIP;   //!
   TBranch        *b_pi_AllIPchi2;   //!
   TBranch        *b_pi_OWNPV_X;   //!
   TBranch        *b_pi_OWNPV_Y;   //!
   TBranch        *b_pi_OWNPV_Z;   //!
   TBranch        *b_pi_OWNPV_XERR;   //!
   TBranch        *b_pi_OWNPV_YERR;   //!
   TBranch        *b_pi_OWNPV_ZERR;   //!
   TBranch        *b_pi_OWNPV_CHI2;   //!
   TBranch        *b_pi_OWNPV_NDOF;   //!
   TBranch        *b_pi_OWNPV_COV_;   //!
   TBranch        *b_pi_IP_OWNPV;   //!
   TBranch        *b_pi_IPCHI2_OWNPV;   //!
   TBranch        *b_pi_TOPPV_X;   //!
   TBranch        *b_pi_TOPPV_Y;   //!
   TBranch        *b_pi_TOPPV_Z;   //!
   TBranch        *b_pi_TOPPV_XERR;   //!
   TBranch        *b_pi_TOPPV_YERR;   //!
   TBranch        *b_pi_TOPPV_ZERR;   //!
   TBranch        *b_pi_TOPPV_CHI2;   //!
   TBranch        *b_pi_TOPPV_NDOF;   //!
   TBranch        *b_pi_TOPPV_COV_;   //!
   TBranch        *b_pi_IP_TOPPV;   //!
   TBranch        *b_pi_IPCHI2_TOPPV;   //!
   TBranch        *b_pi_ORIVX_X;   //!
   TBranch        *b_pi_ORIVX_Y;   //!
   TBranch        *b_pi_ORIVX_Z;   //!
   TBranch        *b_pi_ORIVX_XERR;   //!
   TBranch        *b_pi_ORIVX_YERR;   //!
   TBranch        *b_pi_ORIVX_ZERR;   //!
   TBranch        *b_pi_ORIVX_CHI2;   //!
   TBranch        *b_pi_ORIVX_NDOF;   //!
   TBranch        *b_pi_ORIVX_COV_;   //!
   TBranch        *b_pi_IP_ORIVX;   //!
   TBranch        *b_pi_IPCHI2_ORIVX;   //!
   TBranch        *b_pi_P;   //!
   TBranch        *b_pi_PT;   //!
   TBranch        *b_pi_PE;   //!
   TBranch        *b_pi_PX;   //!
   TBranch        *b_pi_PY;   //!
   TBranch        *b_pi_PZ;   //!
   TBranch        *b_pi_M;   //!
   TBranch        *b_pi_ID;   //!
   TBranch        *b_pi_PIDe;   //!
   TBranch        *b_pi_PIDmu;   //!
   TBranch        *b_pi_PIDK;   //!
   TBranch        *b_pi_PIDp;   //!
   TBranch        *b_pi_ProbNNe;   //!
   TBranch        *b_pi_ProbNNk;   //!
   TBranch        *b_pi_ProbNNp;   //!
   TBranch        *b_pi_ProbNNpi;   //!
   TBranch        *b_pi_ProbNNmu;   //!
   TBranch        *b_pi_ProbNNghost;   //!
   TBranch        *b_pi_hasMuon;   //!
   TBranch        *b_pi_isMuon;   //!
   TBranch        *b_pi_hasRich;   //!
   TBranch        *b_pi_UsedRichAerogel;   //!
   TBranch        *b_pi_UsedRich1Gas;   //!
   TBranch        *b_pi_UsedRich2Gas;   //!
   TBranch        *b_pi_RichAboveElThres;   //!
   TBranch        *b_pi_RichAboveMuThres;   //!
   TBranch        *b_pi_RichAbovePiThres;   //!
   TBranch        *b_pi_RichAboveKaThres;   //!
   TBranch        *b_pi_RichAbovePrThres;   //!
   TBranch        *b_pi_hasCalo;   //!
   TBranch        *b_pi_TRACK_Type;   //!
   TBranch        *b_pi_TRACK_Key;   //!
   TBranch        *b_pi_TRACK_CHI2NDOF;   //!
   TBranch        *b_pi_TRACK_PCHI2;   //!
   TBranch        *b_pi_TRACK_MatchCHI2;   //!
   TBranch        *b_pi_TRACK_GhostProb;   //!
   TBranch        *b_pi_TRACK_CloneDist;   //!
   TBranch        *b_pi_TRACK_Likelihood;   //!
   TBranch        *b_pi_cpx_1_00;   //!
   TBranch        *b_pi_cpy_1_00;   //!
   TBranch        *b_pi_cpz_1_00;   //!
   TBranch        *b_pi_cpt_1_00;   //!
   TBranch        *b_pi_cp_1_00;   //!
   TBranch        *b_pi_cmult_1_00;   //!
   TBranch        *b_pi_pxasy_1_00;   //!
   TBranch        *b_pi_pyasy_1_00;   //!
   TBranch        *b_pi_pzasy_1_00;   //!
   TBranch        *b_pi_pasy_1_00;   //!
   TBranch        *b_pi_ptasy_1_00;   //!
   TBranch        *b_nCandidate;   //!
   TBranch        *b_totCandidates;   //!
   TBranch        *b_EventInSequence;   //!
   TBranch        *b_runNumber;   //!
   TBranch        *b_eventNumber;   //!
   TBranch        *b_BCID;   //!
   TBranch        *b_BCType;   //!
   TBranch        *b_OdinTCK;   //!
   TBranch        *b_L0DUTCK;   //!
   TBranch        *b_HLT1TCK;   //!
   TBranch        *b_HLT2TCK;   //!
   TBranch        *b_GpsTime;   //!
   TBranch        *b_Polarity;   //!
   TBranch        *b_PVX;   //!
   TBranch        *b_PVY;   //!
   TBranch        *b_PVZ;   //!
   TBranch        *b_PVXERR;   //!
   TBranch        *b_PVYERR;   //!
   TBranch        *b_PVZERR;   //!
   TBranch        *b_PVCHI2;   //!
   TBranch        *b_PVNDOF;   //!
   TBranch        *b_PVNTRACKS;   //!
   TBranch        *b_nPVs;   //!
   TBranch        *b_nTracks;   //!
   TBranch        *b_nLongTracks;   //!
   TBranch        *b_nDownstreamTracks;   //!
   TBranch        *b_nUpstreamTracks;   //!
   TBranch        *b_nVeloTracks;   //!
   TBranch        *b_nTTracks;   //!
   TBranch        *b_nBackTracks;   //!
   TBranch        *b_nRich1Hits;   //!
   TBranch        *b_nRich2Hits;   //!
   TBranch        *b_nVeloClusters;   //!
   TBranch        *b_nITClusters;   //!
   TBranch        *b_nTTClusters;   //!
   TBranch        *b_nOTClusters;   //!
   TBranch        *b_nSPDHits;   //!
   TBranch        *b_nMuonCoordsS0;   //!
   TBranch        *b_nMuonCoordsS1;   //!
   TBranch        *b_nMuonCoordsS2;   //!
   TBranch        *b_nMuonCoordsS3;   //!
   TBranch        *b_nMuonCoordsS4;   //!
   TBranch        *b_nMuonTracks;   //!

};

#endif

#ifdef DecayTree_cxx
DecayTree::~DecayTree()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t DecayTree::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t DecayTree::LoadTree(Long64_t entry)
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

void DecayTree::Init()
{
   TTree* tree = this->GetInputTree();
   cout << "Found files, now init" << endl;
    
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("B_ETA", &B_ETA, &b_B_ETA);
   fChain->SetBranchAddress("B_MINIP", &B_MINIP, &b_B_MINIP);
   fChain->SetBranchAddress("B_MINIPCHI2", &B_MINIPCHI2, &b_B_MINIPCHI2);
   fChain->SetBranchAddress("B_MINIPNEXTBEST", &B_MINIPNEXTBEST, &b_B_MINIPNEXTBEST);
   fChain->SetBranchAddress("B_MINIPCHI2NEXTBEST", &B_MINIPCHI2NEXTBEST, &b_B_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   //fChain->SetBranchAddress("B_AllIP", B_AllIP, &b_B_AllIP);
   //fChain->SetBranchAddress("B_AllIPchi2", B_AllIPchi2, &b_B_AllIPchi2);
   //fChain->SetBranchAddress("B_AllDIRA", B_AllDIRA, &b_B_AllDIRA);
   fChain->SetBranchAddress("B_ENDVERTEX_X", &B_ENDVERTEX_X, &b_B_ENDVERTEX_X);
   fChain->SetBranchAddress("B_ENDVERTEX_Y", &B_ENDVERTEX_Y, &b_B_ENDVERTEX_Y);
   fChain->SetBranchAddress("B_ENDVERTEX_Z", &B_ENDVERTEX_Z, &b_B_ENDVERTEX_Z);
   fChain->SetBranchAddress("B_ENDVERTEX_XERR", &B_ENDVERTEX_XERR, &b_B_ENDVERTEX_XERR);
   fChain->SetBranchAddress("B_ENDVERTEX_YERR", &B_ENDVERTEX_YERR, &b_B_ENDVERTEX_YERR);
   fChain->SetBranchAddress("B_ENDVERTEX_ZERR", &B_ENDVERTEX_ZERR, &b_B_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("B_ENDVERTEX_CHI2", &B_ENDVERTEX_CHI2, &b_B_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("B_ENDVERTEX_NDOF", &B_ENDVERTEX_NDOF, &b_B_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("B_ENDVERTEX_COV_", B_ENDVERTEX_COV_, &b_B_ENDVERTEX_COV_);
   fChain->SetBranchAddress("B_OWNPV_X", &B_OWNPV_X, &b_B_OWNPV_X);
   fChain->SetBranchAddress("B_OWNPV_Y", &B_OWNPV_Y, &b_B_OWNPV_Y);
   fChain->SetBranchAddress("B_OWNPV_Z", &B_OWNPV_Z, &b_B_OWNPV_Z);
   fChain->SetBranchAddress("B_OWNPV_XERR", &B_OWNPV_XERR, &b_B_OWNPV_XERR);
   fChain->SetBranchAddress("B_OWNPV_YERR", &B_OWNPV_YERR, &b_B_OWNPV_YERR);
   fChain->SetBranchAddress("B_OWNPV_ZERR", &B_OWNPV_ZERR, &b_B_OWNPV_ZERR);
   fChain->SetBranchAddress("B_OWNPV_CHI2", &B_OWNPV_CHI2, &b_B_OWNPV_CHI2);
   fChain->SetBranchAddress("B_OWNPV_NDOF", &B_OWNPV_NDOF, &b_B_OWNPV_NDOF);
   fChain->SetBranchAddress("B_OWNPV_COV_", B_OWNPV_COV_, &b_B_OWNPV_COV_);
   fChain->SetBranchAddress("B_IP_OWNPV", &B_IP_OWNPV, &b_B_IP_OWNPV);
   fChain->SetBranchAddress("B_IPCHI2_OWNPV", &B_IPCHI2_OWNPV, &b_B_IPCHI2_OWNPV);
   fChain->SetBranchAddress("B_FD_OWNPV", &B_FD_OWNPV, &b_B_FD_OWNPV);
   fChain->SetBranchAddress("B_FDCHI2_OWNPV", &B_FDCHI2_OWNPV, &b_B_FDCHI2_OWNPV);
   fChain->SetBranchAddress("B_DIRA_OWNPV", &B_DIRA_OWNPV, &b_B_DIRA_OWNPV);
   fChain->SetBranchAddress("B_TOPPV_X", &B_TOPPV_X, &b_B_TOPPV_X);
   fChain->SetBranchAddress("B_TOPPV_Y", &B_TOPPV_Y, &b_B_TOPPV_Y);
   fChain->SetBranchAddress("B_TOPPV_Z", &B_TOPPV_Z, &b_B_TOPPV_Z);
   fChain->SetBranchAddress("B_TOPPV_XERR", &B_TOPPV_XERR, &b_B_TOPPV_XERR);
   fChain->SetBranchAddress("B_TOPPV_YERR", &B_TOPPV_YERR, &b_B_TOPPV_YERR);
   fChain->SetBranchAddress("B_TOPPV_ZERR", &B_TOPPV_ZERR, &b_B_TOPPV_ZERR);
   fChain->SetBranchAddress("B_TOPPV_CHI2", &B_TOPPV_CHI2, &b_B_TOPPV_CHI2);
   fChain->SetBranchAddress("B_TOPPV_NDOF", &B_TOPPV_NDOF, &b_B_TOPPV_NDOF);
   fChain->SetBranchAddress("B_TOPPV_COV_", B_TOPPV_COV_, &b_B_TOPPV_COV_);
   fChain->SetBranchAddress("B_IP_TOPPV", &B_IP_TOPPV, &b_B_IP_TOPPV);
   fChain->SetBranchAddress("B_IPCHI2_TOPPV", &B_IPCHI2_TOPPV, &b_B_IPCHI2_TOPPV);
   fChain->SetBranchAddress("B_FD_TOPPV", &B_FD_TOPPV, &b_B_FD_TOPPV);
   fChain->SetBranchAddress("B_FDCHI2_TOPPV", &B_FDCHI2_TOPPV, &b_B_FDCHI2_TOPPV);
   fChain->SetBranchAddress("B_DIRA_TOPPV", &B_DIRA_TOPPV, &b_B_DIRA_TOPPV);
   fChain->SetBranchAddress("B_P", &B_P, &b_B_P);
   fChain->SetBranchAddress("B_PT", &B_PT, &b_B_PT);
   fChain->SetBranchAddress("B_PE", &B_PE, &b_B_PE);
   fChain->SetBranchAddress("B_PX", &B_PX, &b_B_PX);
   fChain->SetBranchAddress("B_PY", &B_PY, &b_B_PY);
   fChain->SetBranchAddress("B_PZ", &B_PZ, &b_B_PZ);
   fChain->SetBranchAddress("B_MM", &B_MM, &b_B_MM);
   fChain->SetBranchAddress("B_MMERR", &B_MMERR, &b_B_MMERR);
   fChain->SetBranchAddress("B_M", &B_M, &b_B_M);
   fChain->SetBranchAddress("B_ID", &B_ID, &b_B_ID);
   /*
   fChain->SetBranchAddress("B_TAGDECISION", &B_TAGDECISION, &b_B_TAGDECISION);
   fChain->SetBranchAddress("B_TAGOMEGA", &B_TAGOMEGA, &b_B_TAGOMEGA);
   fChain->SetBranchAddress("B_TAGDECISION_OS", &B_TAGDECISION_OS, &b_B_TAGDECISION_OS);
   fChain->SetBranchAddress("B_TAGOMEGA_OS", &B_TAGOMEGA_OS, &b_B_TAGOMEGA_OS);
   fChain->SetBranchAddress("B_TAGGER", &B_TAGGER, &b_B_TAGGER);
   fChain->SetBranchAddress("B_OS_Muon_DEC", &B_OS_Muon_DEC, &b_B_OS_Muon_DEC);
   fChain->SetBranchAddress("B_OS_Muon_PROB", &B_OS_Muon_PROB, &b_B_OS_Muon_PROB);
   fChain->SetBranchAddress("B_OS_Electron_DEC", &B_OS_Electron_DEC, &b_B_OS_Electron_DEC);
   fChain->SetBranchAddress("B_OS_Electron_PROB", &B_OS_Electron_PROB, &b_B_OS_Electron_PROB);
   fChain->SetBranchAddress("B_OS_Kaon_DEC", &B_OS_Kaon_DEC, &b_B_OS_Kaon_DEC);
   fChain->SetBranchAddress("B_OS_Kaon_PROB", &B_OS_Kaon_PROB, &b_B_OS_Kaon_PROB);
   fChain->SetBranchAddress("B_SS_Kaon_DEC", &B_SS_Kaon_DEC, &b_B_SS_Kaon_DEC);
   fChain->SetBranchAddress("B_SS_Kaon_PROB", &B_SS_Kaon_PROB, &b_B_SS_Kaon_PROB);
   fChain->SetBranchAddress("B_SS_Pion_DEC", &B_SS_Pion_DEC, &b_B_SS_Pion_DEC);
   fChain->SetBranchAddress("B_SS_Pion_PROB", &B_SS_Pion_PROB, &b_B_SS_Pion_PROB);
   fChain->SetBranchAddress("B_SS_PionBDT_DEC", &B_SS_PionBDT_DEC, &b_B_SS_PionBDT_DEC);
   fChain->SetBranchAddress("B_SS_PionBDT_PROB", &B_SS_PionBDT_PROB, &b_B_SS_PionBDT_PROB);
   fChain->SetBranchAddress("B_VtxCharge_DEC", &B_VtxCharge_DEC, &b_B_VtxCharge_DEC);
   fChain->SetBranchAddress("B_VtxCharge_PROB", &B_VtxCharge_PROB, &b_B_VtxCharge_PROB);
   fChain->SetBranchAddress("B_OS_nnetKaon_DEC", &B_OS_nnetKaon_DEC, &b_B_OS_nnetKaon_DEC);
   fChain->SetBranchAddress("B_OS_nnetKaon_PROB", &B_OS_nnetKaon_PROB, &b_B_OS_nnetKaon_PROB);
   fChain->SetBranchAddress("B_SS_nnetKaon_DEC", &B_SS_nnetKaon_DEC, &b_B_SS_nnetKaon_DEC);
   fChain->SetBranchAddress("B_SS_nnetKaon_PROB", &B_SS_nnetKaon_PROB, &b_B_SS_nnetKaon_PROB);
   fChain->SetBranchAddress("B_SS_Proton_DEC", &B_SS_Proton_DEC, &b_B_SS_Proton_DEC);
   fChain->SetBranchAddress("B_SS_Proton_PROB", &B_SS_Proton_PROB, &b_B_SS_Proton_PROB);
   fChain->SetBranchAddress("B_OS_Charm_DEC", &B_OS_Charm_DEC, &b_B_OS_Charm_DEC);
   fChain->SetBranchAddress("B_OS_Charm_PROB", &B_OS_Charm_PROB, &b_B_OS_Charm_PROB);
   fChain->SetBranchAddress("B_cpx_1.00", &B_cpx_1_00, &b_B_cpx_1_00);
   fChain->SetBranchAddress("B_cpy_1.00", &B_cpy_1_00, &b_B_cpy_1_00);
   fChain->SetBranchAddress("B_cpz_1.00", &B_cpz_1_00, &b_B_cpz_1_00);
   fChain->SetBranchAddress("B_cpt_1.00", &B_cpt_1_00, &b_B_cpt_1_00);
   fChain->SetBranchAddress("B_cp_1.00", &B_cp_1_00, &b_B_cp_1_00);
   fChain->SetBranchAddress("B_cmult_1.00", &B_cmult_1_00, &b_B_cmult_1_00);
   fChain->SetBranchAddress("B_pxasy_1.00", &B_pxasy_1_00, &b_B_pxasy_1_00);
   fChain->SetBranchAddress("B_pyasy_1.00", &B_pyasy_1_00, &b_B_pyasy_1_00);
   fChain->SetBranchAddress("B_pzasy_1.00", &B_pzasy_1_00, &b_B_pzasy_1_00);
   fChain->SetBranchAddress("B_pasy_1.00", &B_pasy_1_00, &b_B_pasy_1_00);
   fChain->SetBranchAddress("B_ptasy_1.00", &B_ptasy_1_00, &b_B_ptasy_1_00); */
   fChain->SetBranchAddress("B_DOCA1", &B_DOCA1, &b_B_DOCA1);
   fChain->SetBranchAddress("B_TAU", &B_TAU, &b_B_TAU);
   fChain->SetBranchAddress("B_TAUERR", &B_TAUERR, &b_B_TAUERR);
   fChain->SetBranchAddress("B_BDTF_nPV", &B_BDTF_nPV, &b_B_BDTF_nPV);
   fChain->SetBranchAddress("B_BDTF_Dplus_M", B_BDTF_Dplus_M, &b_B_BDTF_Dplus_M);
   fChain->SetBranchAddress("B_BDTF_Dplus_MERR", B_BDTF_Dplus_MERR, &b_B_BDTF_Dplus_MERR);
   fChain->SetBranchAddress("B_BDTF_Dplus_P", B_BDTF_Dplus_P, &b_B_BDTF_Dplus_P);
   fChain->SetBranchAddress("B_BDTF_Dplus_PERR", B_BDTF_Dplus_PERR, &b_B_BDTF_Dplus_PERR);
   fChain->SetBranchAddress("B_BDTF_Dplus_ctau", B_BDTF_Dplus_ctau, &b_B_BDTF_Dplus_ctau);
   fChain->SetBranchAddress("B_BDTF_Dplus_ctauErr", B_BDTF_Dplus_ctauErr, &b_B_BDTF_Dplus_ctauErr);
   fChain->SetBranchAddress("B_BDTF_Dplus_decayLength", B_BDTF_Dplus_decayLength, &b_B_BDTF_Dplus_decayLength);
   fChain->SetBranchAddress("B_BDTF_Dplus_decayLengthErr", B_BDTF_Dplus_decayLengthErr, &b_B_BDTF_Dplus_decayLengthErr);
   fChain->SetBranchAddress("B_BDTF_KS0_M", B_BDTF_KS0_M, &b_B_BDTF_KS0_M);
   fChain->SetBranchAddress("B_BDTF_KS0_MERR", B_BDTF_KS0_MERR, &b_B_BDTF_KS0_MERR);
   fChain->SetBranchAddress("B_BDTF_KS0_P", B_BDTF_KS0_P, &b_B_BDTF_KS0_P);
   fChain->SetBranchAddress("B_BDTF_KS0_PERR", B_BDTF_KS0_PERR, &b_B_BDTF_KS0_PERR);
   fChain->SetBranchAddress("B_BDTF_KS0_ctau", B_BDTF_KS0_ctau, &b_B_BDTF_KS0_ctau);
   fChain->SetBranchAddress("B_BDTF_KS0_ctauErr", B_BDTF_KS0_ctauErr, &b_B_BDTF_KS0_ctauErr);
   fChain->SetBranchAddress("B_BDTF_KS0_decayLength", B_BDTF_KS0_decayLength, &b_B_BDTF_KS0_decayLength);
   fChain->SetBranchAddress("B_BDTF_KS0_decayLengthErr", B_BDTF_KS0_decayLengthErr, &b_B_BDTF_KS0_decayLengthErr);
   fChain->SetBranchAddress("B_BDTF_M", B_BDTF_M, &b_B_BDTF_M);
   fChain->SetBranchAddress("B_BDTF_MERR", B_BDTF_MERR, &b_B_BDTF_MERR);
   fChain->SetBranchAddress("B_BDTF_P", B_BDTF_P, &b_B_BDTF_P);
   fChain->SetBranchAddress("B_BDTF_PERR", B_BDTF_PERR, &b_B_BDTF_PERR);
   fChain->SetBranchAddress("B_BDTF_PV_X", B_BDTF_PV_X, &b_B_BDTF_PV_X);
   fChain->SetBranchAddress("B_BDTF_PV_Y", B_BDTF_PV_Y, &b_B_BDTF_PV_Y);
   fChain->SetBranchAddress("B_BDTF_PV_Z", B_BDTF_PV_Z, &b_B_BDTF_PV_Z);
   fChain->SetBranchAddress("B_BDTF_PV_key", B_BDTF_PV_key, &b_B_BDTF_PV_key);
   fChain->SetBranchAddress("B_BDTF_chi2", B_BDTF_chi2, &b_B_BDTF_chi2);
   fChain->SetBranchAddress("B_BDTF_ctau", B_BDTF_ctau, &b_B_BDTF_ctau);
   fChain->SetBranchAddress("B_BDTF_ctauErr", B_BDTF_ctauErr, &b_B_BDTF_ctauErr);
   fChain->SetBranchAddress("B_BDTF_decayLength", B_BDTF_decayLength, &b_B_BDTF_decayLength);
   fChain->SetBranchAddress("B_BDTF_decayLengthErr", B_BDTF_decayLengthErr, &b_B_BDTF_decayLengthErr);
   fChain->SetBranchAddress("B_BDTF_nDOF", B_BDTF_nDOF, &b_B_BDTF_nDOF);
   fChain->SetBranchAddress("B_BDTF_nIter", B_BDTF_nIter, &b_B_BDTF_nIter);
   fChain->SetBranchAddress("B_BDTF_status", B_BDTF_status, &b_B_BDTF_status);
   fChain->SetBranchAddress("B_DTF_nPV", &B_DTF_nPV, &b_B_DTF_nPV);
   fChain->SetBranchAddress("B_DTF_Dplus_M", B_DTF_Dplus_M, &b_B_DTF_Dplus_M);
   fChain->SetBranchAddress("B_DTF_Dplus_MERR", B_DTF_Dplus_MERR, &b_B_DTF_Dplus_MERR);
   fChain->SetBranchAddress("B_DTF_Dplus_P", B_DTF_Dplus_P, &b_B_DTF_Dplus_P);
   fChain->SetBranchAddress("B_DTF_Dplus_PERR", B_DTF_Dplus_PERR, &b_B_DTF_Dplus_PERR);
   fChain->SetBranchAddress("B_DTF_Dplus_ctau", B_DTF_Dplus_ctau, &b_B_DTF_Dplus_ctau);
   fChain->SetBranchAddress("B_DTF_Dplus_ctauErr", B_DTF_Dplus_ctauErr, &b_B_DTF_Dplus_ctauErr);
   fChain->SetBranchAddress("B_DTF_Dplus_decayLength", B_DTF_Dplus_decayLength, &b_B_DTF_Dplus_decayLength);
   fChain->SetBranchAddress("B_DTF_Dplus_decayLengthErr", B_DTF_Dplus_decayLengthErr, &b_B_DTF_Dplus_decayLengthErr);
   fChain->SetBranchAddress("B_DTF_KS0_M", B_DTF_KS0_M, &b_B_DTF_KS0_M);
   fChain->SetBranchAddress("B_DTF_KS0_MERR", B_DTF_KS0_MERR, &b_B_DTF_KS0_MERR);
   fChain->SetBranchAddress("B_DTF_KS0_P", B_DTF_KS0_P, &b_B_DTF_KS0_P);
   fChain->SetBranchAddress("B_DTF_KS0_PERR", B_DTF_KS0_PERR, &b_B_DTF_KS0_PERR);
   fChain->SetBranchAddress("B_DTF_KS0_ctau", B_DTF_KS0_ctau, &b_B_DTF_KS0_ctau);
   fChain->SetBranchAddress("B_DTF_KS0_ctauErr", B_DTF_KS0_ctauErr, &b_B_DTF_KS0_ctauErr);
   fChain->SetBranchAddress("B_DTF_KS0_decayLength", B_DTF_KS0_decayLength, &b_B_DTF_KS0_decayLength);
   fChain->SetBranchAddress("B_DTF_KS0_decayLengthErr", B_DTF_KS0_decayLengthErr, &b_B_DTF_KS0_decayLengthErr);
   fChain->SetBranchAddress("B_DTF_M", B_DTF_M, &b_B_DTF_M);
   fChain->SetBranchAddress("B_DTF_MERR", B_DTF_MERR, &b_B_DTF_MERR);
   fChain->SetBranchAddress("B_DTF_P", B_DTF_P, &b_B_DTF_P);
   fChain->SetBranchAddress("B_DTF_PERR", B_DTF_PERR, &b_B_DTF_PERR);
   fChain->SetBranchAddress("B_DTF_PV_X", B_DTF_PV_X, &b_B_DTF_PV_X);
   fChain->SetBranchAddress("B_DTF_PV_Y", B_DTF_PV_Y, &b_B_DTF_PV_Y);
   fChain->SetBranchAddress("B_DTF_PV_Z", B_DTF_PV_Z, &b_B_DTF_PV_Z);
   fChain->SetBranchAddress("B_DTF_PV_key", B_DTF_PV_key, &b_B_DTF_PV_key);
   fChain->SetBranchAddress("B_DTF_chi2", B_DTF_chi2, &b_B_DTF_chi2);
   fChain->SetBranchAddress("B_DTF_ctau", B_DTF_ctau, &b_B_DTF_ctau);
   fChain->SetBranchAddress("B_DTF_ctauErr", B_DTF_ctauErr, &b_B_DTF_ctauErr);
   fChain->SetBranchAddress("B_DTF_decayLength", B_DTF_decayLength, &b_B_DTF_decayLength);
   fChain->SetBranchAddress("B_DTF_decayLengthErr", B_DTF_decayLengthErr, &b_B_DTF_decayLengthErr);
   fChain->SetBranchAddress("B_DTF_nDOF", B_DTF_nDOF, &b_B_DTF_nDOF);
   fChain->SetBranchAddress("B_DTF_nIter", B_DTF_nIter, &b_B_DTF_nIter);
   fChain->SetBranchAddress("B_DTF_status", B_DTF_status, &b_B_DTF_status);
   fChain->SetBranchAddress("B_FullDTF_nPV", &B_FullDTF_nPV, &b_B_FullDTF_nPV);
   fChain->SetBranchAddress("B_FullDTF_Dplus_Kplus_ID", B_FullDTF_Dplus_Kplus_ID, &b_B_FullDTF_Dplus_Kplus_ID);
   fChain->SetBranchAddress("B_FullDTF_Dplus_Kplus_PE", B_FullDTF_Dplus_Kplus_PE, &b_B_FullDTF_Dplus_Kplus_PE);
   fChain->SetBranchAddress("B_FullDTF_Dplus_Kplus_PX", B_FullDTF_Dplus_Kplus_PX, &b_B_FullDTF_Dplus_Kplus_PX);
   fChain->SetBranchAddress("B_FullDTF_Dplus_Kplus_PY", B_FullDTF_Dplus_Kplus_PY, &b_B_FullDTF_Dplus_Kplus_PY);
   fChain->SetBranchAddress("B_FullDTF_Dplus_Kplus_PZ", B_FullDTF_Dplus_Kplus_PZ, &b_B_FullDTF_Dplus_Kplus_PZ);
   fChain->SetBranchAddress("B_FullDTF_Dplus_M", B_FullDTF_Dplus_M, &b_B_FullDTF_Dplus_M);
   fChain->SetBranchAddress("B_FullDTF_Dplus_MERR", B_FullDTF_Dplus_MERR, &b_B_FullDTF_Dplus_MERR);
   fChain->SetBranchAddress("B_FullDTF_Dplus_P", B_FullDTF_Dplus_P, &b_B_FullDTF_Dplus_P);
   fChain->SetBranchAddress("B_FullDTF_Dplus_PERR", B_FullDTF_Dplus_PERR, &b_B_FullDTF_Dplus_PERR);
   fChain->SetBranchAddress("B_FullDTF_Dplus_ctau", B_FullDTF_Dplus_ctau, &b_B_FullDTF_Dplus_ctau);
   fChain->SetBranchAddress("B_FullDTF_Dplus_ctauErr", B_FullDTF_Dplus_ctauErr, &b_B_FullDTF_Dplus_ctauErr);
   fChain->SetBranchAddress("B_FullDTF_Dplus_decayLength", B_FullDTF_Dplus_decayLength, &b_B_FullDTF_Dplus_decayLength);
   fChain->SetBranchAddress("B_FullDTF_Dplus_decayLengthErr", B_FullDTF_Dplus_decayLengthErr, &b_B_FullDTF_Dplus_decayLengthErr);
   fChain->SetBranchAddress("B_FullDTF_Dplus_piplus_0_ID", B_FullDTF_Dplus_piplus_0_ID, &b_B_FullDTF_Dplus_piplus_0_ID);
   fChain->SetBranchAddress("B_FullDTF_Dplus_piplus_0_PE", B_FullDTF_Dplus_piplus_0_PE, &b_B_FullDTF_Dplus_piplus_0_PE);
   fChain->SetBranchAddress("B_FullDTF_Dplus_piplus_0_PX", B_FullDTF_Dplus_piplus_0_PX, &b_B_FullDTF_Dplus_piplus_0_PX);
   fChain->SetBranchAddress("B_FullDTF_Dplus_piplus_0_PY", B_FullDTF_Dplus_piplus_0_PY, &b_B_FullDTF_Dplus_piplus_0_PY);
   fChain->SetBranchAddress("B_FullDTF_Dplus_piplus_0_PZ", B_FullDTF_Dplus_piplus_0_PZ, &b_B_FullDTF_Dplus_piplus_0_PZ);
   fChain->SetBranchAddress("B_FullDTF_Dplus_piplus_ID", B_FullDTF_Dplus_piplus_ID, &b_B_FullDTF_Dplus_piplus_ID);
   fChain->SetBranchAddress("B_FullDTF_Dplus_piplus_PE", B_FullDTF_Dplus_piplus_PE, &b_B_FullDTF_Dplus_piplus_PE);
   fChain->SetBranchAddress("B_FullDTF_Dplus_piplus_PX", B_FullDTF_Dplus_piplus_PX, &b_B_FullDTF_Dplus_piplus_PX);
   fChain->SetBranchAddress("B_FullDTF_Dplus_piplus_PY", B_FullDTF_Dplus_piplus_PY, &b_B_FullDTF_Dplus_piplus_PY);
   fChain->SetBranchAddress("B_FullDTF_Dplus_piplus_PZ", B_FullDTF_Dplus_piplus_PZ, &b_B_FullDTF_Dplus_piplus_PZ);
   fChain->SetBranchAddress("B_FullDTF_KS0_M", B_FullDTF_KS0_M, &b_B_FullDTF_KS0_M);
   fChain->SetBranchAddress("B_FullDTF_KS0_MERR", B_FullDTF_KS0_MERR, &b_B_FullDTF_KS0_MERR);
   fChain->SetBranchAddress("B_FullDTF_KS0_P", B_FullDTF_KS0_P, &b_B_FullDTF_KS0_P);
   fChain->SetBranchAddress("B_FullDTF_KS0_PERR", B_FullDTF_KS0_PERR, &b_B_FullDTF_KS0_PERR);
   fChain->SetBranchAddress("B_FullDTF_KS0_ctau", B_FullDTF_KS0_ctau, &b_B_FullDTF_KS0_ctau);
   fChain->SetBranchAddress("B_FullDTF_KS0_ctauErr", B_FullDTF_KS0_ctauErr, &b_B_FullDTF_KS0_ctauErr);
   fChain->SetBranchAddress("B_FullDTF_KS0_decayLength", B_FullDTF_KS0_decayLength, &b_B_FullDTF_KS0_decayLength);
   fChain->SetBranchAddress("B_FullDTF_KS0_decayLengthErr", B_FullDTF_KS0_decayLengthErr, &b_B_FullDTF_KS0_decayLengthErr);
   fChain->SetBranchAddress("B_FullDTF_KS0_piplus_0_ID", B_FullDTF_KS0_piplus_0_ID, &b_B_FullDTF_KS0_piplus_0_ID);
   fChain->SetBranchAddress("B_FullDTF_KS0_piplus_0_PE", B_FullDTF_KS0_piplus_0_PE, &b_B_FullDTF_KS0_piplus_0_PE);
   fChain->SetBranchAddress("B_FullDTF_KS0_piplus_0_PX", B_FullDTF_KS0_piplus_0_PX, &b_B_FullDTF_KS0_piplus_0_PX);
   fChain->SetBranchAddress("B_FullDTF_KS0_piplus_0_PY", B_FullDTF_KS0_piplus_0_PY, &b_B_FullDTF_KS0_piplus_0_PY);
   fChain->SetBranchAddress("B_FullDTF_KS0_piplus_0_PZ", B_FullDTF_KS0_piplus_0_PZ, &b_B_FullDTF_KS0_piplus_0_PZ);
   fChain->SetBranchAddress("B_FullDTF_KS0_piplus_ID", B_FullDTF_KS0_piplus_ID, &b_B_FullDTF_KS0_piplus_ID);
   fChain->SetBranchAddress("B_FullDTF_KS0_piplus_PE", B_FullDTF_KS0_piplus_PE, &b_B_FullDTF_KS0_piplus_PE);
   fChain->SetBranchAddress("B_FullDTF_KS0_piplus_PX", B_FullDTF_KS0_piplus_PX, &b_B_FullDTF_KS0_piplus_PX);
   fChain->SetBranchAddress("B_FullDTF_KS0_piplus_PY", B_FullDTF_KS0_piplus_PY, &b_B_FullDTF_KS0_piplus_PY);
   fChain->SetBranchAddress("B_FullDTF_KS0_piplus_PZ", B_FullDTF_KS0_piplus_PZ, &b_B_FullDTF_KS0_piplus_PZ);
   fChain->SetBranchAddress("B_FullDTF_M", B_FullDTF_M, &b_B_FullDTF_M);
   fChain->SetBranchAddress("B_FullDTF_MERR", B_FullDTF_MERR, &b_B_FullDTF_MERR);
   fChain->SetBranchAddress("B_FullDTF_P", B_FullDTF_P, &b_B_FullDTF_P);
   fChain->SetBranchAddress("B_FullDTF_PERR", B_FullDTF_PERR, &b_B_FullDTF_PERR);
   fChain->SetBranchAddress("B_FullDTF_PV_X", B_FullDTF_PV_X, &b_B_FullDTF_PV_X);
   fChain->SetBranchAddress("B_FullDTF_PV_Y", B_FullDTF_PV_Y, &b_B_FullDTF_PV_Y);
   fChain->SetBranchAddress("B_FullDTF_PV_Z", B_FullDTF_PV_Z, &b_B_FullDTF_PV_Z);
   fChain->SetBranchAddress("B_FullDTF_PV_key", B_FullDTF_PV_key, &b_B_FullDTF_PV_key);
   fChain->SetBranchAddress("B_FullDTF_chi2", B_FullDTF_chi2, &b_B_FullDTF_chi2);
   fChain->SetBranchAddress("B_FullDTF_ctau", B_FullDTF_ctau, &b_B_FullDTF_ctau);
   fChain->SetBranchAddress("B_FullDTF_ctauErr", B_FullDTF_ctauErr, &b_B_FullDTF_ctauErr);
   fChain->SetBranchAddress("B_FullDTF_decayLength", B_FullDTF_decayLength, &b_B_FullDTF_decayLength);
   fChain->SetBranchAddress("B_FullDTF_decayLengthErr", B_FullDTF_decayLengthErr, &b_B_FullDTF_decayLengthErr);
   fChain->SetBranchAddress("B_FullDTF_nDOF", B_FullDTF_nDOF, &b_B_FullDTF_nDOF);
   fChain->SetBranchAddress("B_FullDTF_nIter", B_FullDTF_nIter, &b_B_FullDTF_nIter);
    fChain->SetBranchAddress("B_FullDTF_status", B_FullDTF_status, &b_B_FullDTF_status);

    fChain->SetBranchAddress("B_FullBsDTF_nPV", &B_FullBsDTF_nPV, &b_B_FullBsDTF_nPV);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_Kplus_ID", B_FullBsDTF_Dplus_Kplus_ID, &b_B_FullBsDTF_Dplus_Kplus_ID);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_Kplus_PE", B_FullBsDTF_Dplus_Kplus_PE, &b_B_FullBsDTF_Dplus_Kplus_PE);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_Kplus_PX", B_FullBsDTF_Dplus_Kplus_PX, &b_B_FullBsDTF_Dplus_Kplus_PX);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_Kplus_PY", B_FullBsDTF_Dplus_Kplus_PY, &b_B_FullBsDTF_Dplus_Kplus_PY);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_Kplus_PZ", B_FullBsDTF_Dplus_Kplus_PZ, &b_B_FullBsDTF_Dplus_Kplus_PZ);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_M", B_FullBsDTF_Dplus_M, &b_B_FullBsDTF_Dplus_M);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_MERR", B_FullBsDTF_Dplus_MERR, &b_B_FullBsDTF_Dplus_MERR);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_P", B_FullBsDTF_Dplus_P, &b_B_FullBsDTF_Dplus_P);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_PERR", B_FullBsDTF_Dplus_PERR, &b_B_FullBsDTF_Dplus_PERR);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_ctau", B_FullBsDTF_Dplus_ctau, &b_B_FullBsDTF_Dplus_ctau);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_ctauErr", B_FullBsDTF_Dplus_ctauErr, &b_B_FullBsDTF_Dplus_ctauErr);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_decayLength", B_FullBsDTF_Dplus_decayLength, &b_B_FullBsDTF_Dplus_decayLength);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_decayLengthErr", B_FullBsDTF_Dplus_decayLengthErr, &b_B_FullBsDTF_Dplus_decayLengthErr);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_piplus_0_ID", B_FullBsDTF_Dplus_piplus_0_ID, &b_B_FullBsDTF_Dplus_piplus_0_ID);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_piplus_0_PE", B_FullBsDTF_Dplus_piplus_0_PE, &b_B_FullBsDTF_Dplus_piplus_0_PE);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_piplus_0_PX", B_FullBsDTF_Dplus_piplus_0_PX, &b_B_FullBsDTF_Dplus_piplus_0_PX);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_piplus_0_PY", B_FullBsDTF_Dplus_piplus_0_PY, &b_B_FullBsDTF_Dplus_piplus_0_PY);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_piplus_0_PZ", B_FullBsDTF_Dplus_piplus_0_PZ, &b_B_FullBsDTF_Dplus_piplus_0_PZ);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_piplus_ID", B_FullBsDTF_Dplus_piplus_ID, &b_B_FullBsDTF_Dplus_piplus_ID);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_piplus_PE", B_FullBsDTF_Dplus_piplus_PE, &b_B_FullBsDTF_Dplus_piplus_PE);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_piplus_PX", B_FullBsDTF_Dplus_piplus_PX, &b_B_FullBsDTF_Dplus_piplus_PX);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_piplus_PY", B_FullBsDTF_Dplus_piplus_PY, &b_B_FullBsDTF_Dplus_piplus_PY);
    fChain->SetBranchAddress("B_FullBsDTF_Dplus_piplus_PZ", B_FullBsDTF_Dplus_piplus_PZ, &b_B_FullBsDTF_Dplus_piplus_PZ);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_M", B_FullBsDTF_KS0_M, &b_B_FullBsDTF_KS0_M);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_MERR", B_FullBsDTF_KS0_MERR, &b_B_FullBsDTF_KS0_MERR);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_P", B_FullBsDTF_KS0_P, &b_B_FullBsDTF_KS0_P);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_PERR", B_FullBsDTF_KS0_PERR, &b_B_FullBsDTF_KS0_PERR);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_ctau", B_FullBsDTF_KS0_ctau, &b_B_FullBsDTF_KS0_ctau);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_ctauErr", B_FullBsDTF_KS0_ctauErr, &b_B_FullBsDTF_KS0_ctauErr);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_decayLength", B_FullBsDTF_KS0_decayLength, &b_B_FullBsDTF_KS0_decayLength);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_decayLengthErr", B_FullBsDTF_KS0_decayLengthErr, &b_B_FullBsDTF_KS0_decayLengthErr);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_piplus_0_ID", B_FullBsDTF_KS0_piplus_0_ID, &b_B_FullBsDTF_KS0_piplus_0_ID);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_piplus_0_PE", B_FullBsDTF_KS0_piplus_0_PE, &b_B_FullBsDTF_KS0_piplus_0_PE);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_piplus_0_PX", B_FullBsDTF_KS0_piplus_0_PX, &b_B_FullBsDTF_KS0_piplus_0_PX);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_piplus_0_PY", B_FullBsDTF_KS0_piplus_0_PY, &b_B_FullBsDTF_KS0_piplus_0_PY);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_piplus_0_PZ", B_FullBsDTF_KS0_piplus_0_PZ, &b_B_FullBsDTF_KS0_piplus_0_PZ);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_piplus_ID", B_FullBsDTF_KS0_piplus_ID, &b_B_FullBsDTF_KS0_piplus_ID);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_piplus_PE", B_FullBsDTF_KS0_piplus_PE, &b_B_FullBsDTF_KS0_piplus_PE);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_piplus_PX", B_FullBsDTF_KS0_piplus_PX, &b_B_FullBsDTF_KS0_piplus_PX);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_piplus_PY", B_FullBsDTF_KS0_piplus_PY, &b_B_FullBsDTF_KS0_piplus_PY);
    fChain->SetBranchAddress("B_FullBsDTF_KS0_piplus_PZ", B_FullBsDTF_KS0_piplus_PZ, &b_B_FullBsDTF_KS0_piplus_PZ);
    fChain->SetBranchAddress("B_FullBsDTF_M", B_FullBsDTF_M, &b_B_FullBsDTF_M);
    fChain->SetBranchAddress("B_FullBsDTF_MERR", B_FullBsDTF_MERR, &b_B_FullBsDTF_MERR);
    fChain->SetBranchAddress("B_FullBsDTF_P", B_FullBsDTF_P, &b_B_FullBsDTF_P);
    fChain->SetBranchAddress("B_FullBsDTF_PERR", B_FullBsDTF_PERR, &b_B_FullBsDTF_PERR);
    fChain->SetBranchAddress("B_FullBsDTF_PV_X", B_FullBsDTF_PV_X, &b_B_FullBsDTF_PV_X);
    fChain->SetBranchAddress("B_FullBsDTF_PV_Y", B_FullBsDTF_PV_Y, &b_B_FullBsDTF_PV_Y);
    fChain->SetBranchAddress("B_FullBsDTF_PV_Z", B_FullBsDTF_PV_Z, &b_B_FullBsDTF_PV_Z);
    fChain->SetBranchAddress("B_FullBsDTF_PV_key", B_FullBsDTF_PV_key, &b_B_FullBsDTF_PV_key);
    fChain->SetBranchAddress("B_FullBsDTF_chi2", B_FullBsDTF_chi2, &b_B_FullBsDTF_chi2);
    fChain->SetBranchAddress("B_FullBsDTF_ctau", B_FullBsDTF_ctau, &b_B_FullBsDTF_ctau);
    fChain->SetBranchAddress("B_FullBsDTF_ctauErr", B_FullBsDTF_ctauErr, &b_B_FullBsDTF_ctauErr);
    fChain->SetBranchAddress("B_FullBsDTF_decayLength", B_FullBsDTF_decayLength, &b_B_FullBsDTF_decayLength);
    fChain->SetBranchAddress("B_FullBsDTF_decayLengthErr", B_FullBsDTF_decayLengthErr, &b_B_FullBsDTF_decayLengthErr);
    fChain->SetBranchAddress("B_FullBsDTF_nDOF", B_FullBsDTF_nDOF, &b_B_FullBsDTF_nDOF);
    fChain->SetBranchAddress("B_FullBsDTF_nIter", B_FullBsDTF_nIter, &b_B_FullBsDTF_nIter);
    fChain->SetBranchAddress("B_FullBsDTF_status", B_FullBsDTF_status, &b_B_FullBsDTF_status);        
        
    if(_decay==Decay::B2DKspi_LL || _decay==Decay::B2DKspi_DD){
        fChain->SetBranchAddress("B_FullDTF_piplus_ID", B_FullDTF_piplus_ID, &b_B_FullDTF_piplus_ID);
        fChain->SetBranchAddress("B_FullDTF_piplus_PE", B_FullDTF_piplus_PE, &b_B_FullDTF_piplus_PE);
        fChain->SetBranchAddress("B_FullDTF_piplus_PX", B_FullDTF_piplus_PX, &b_B_FullDTF_piplus_PX);
        fChain->SetBranchAddress("B_FullDTF_piplus_PY", B_FullDTF_piplus_PY, &b_B_FullDTF_piplus_PY);
        fChain->SetBranchAddress("B_FullDTF_piplus_PZ", B_FullDTF_piplus_PZ, &b_B_FullDTF_piplus_PZ);

        fChain->SetBranchAddress("B_FullBsDTF_piplus_ID", B_FullBsDTF_piplus_ID, &b_B_FullBsDTF_piplus_ID);
        fChain->SetBranchAddress("B_FullBsDTF_piplus_PE", B_FullBsDTF_piplus_PE, &b_B_FullBsDTF_piplus_PE);
        fChain->SetBranchAddress("B_FullBsDTF_piplus_PX", B_FullBsDTF_piplus_PX, &b_B_FullBsDTF_piplus_PX);
        fChain->SetBranchAddress("B_FullBsDTF_piplus_PY", B_FullBsDTF_piplus_PY, &b_B_FullBsDTF_piplus_PY);
        fChain->SetBranchAddress("B_FullBsDTF_piplus_PZ", B_FullBsDTF_piplus_PZ, &b_B_FullBsDTF_piplus_PZ);        
    }
    else {
        fChain->SetBranchAddress("B_FullDTF_Kplus_ID", B_FullDTF_piplus_ID, &b_B_FullDTF_piplus_ID);
        fChain->SetBranchAddress("B_FullDTF_Kplus_PE", B_FullDTF_piplus_PE, &b_B_FullDTF_piplus_PE);
        fChain->SetBranchAddress("B_FullDTF_Kplus_PX", B_FullDTF_piplus_PX, &b_B_FullDTF_piplus_PX);
        fChain->SetBranchAddress("B_FullDTF_Kplus_PY", B_FullDTF_piplus_PY, &b_B_FullDTF_piplus_PY);
        fChain->SetBranchAddress("B_FullDTF_Kplus_PZ", B_FullDTF_piplus_PZ, &b_B_FullDTF_piplus_PZ);
        
        fChain->SetBranchAddress("B_FullBsDTF_Kplus_ID", B_FullBsDTF_piplus_ID, &b_B_FullBsDTF_piplus_ID);
        fChain->SetBranchAddress("B_FullBsDTF_Kplus_PE", B_FullBsDTF_piplus_PE, &b_B_FullBsDTF_piplus_PE);
        fChain->SetBranchAddress("B_FullBsDTF_Kplus_PX", B_FullBsDTF_piplus_PX, &b_B_FullBsDTF_piplus_PX);
        fChain->SetBranchAddress("B_FullBsDTF_Kplus_PY", B_FullBsDTF_piplus_PY, &b_B_FullBsDTF_piplus_PY);
        fChain->SetBranchAddress("B_FullBsDTF_Kplus_PZ", B_FullBsDTF_piplus_PZ, &b_B_FullBsDTF_piplus_PZ);
    }

    //fChain->SetBranchAddress("B_PV_piplus_ID", B_PV_piplus_ID, &b_B_PV_piplus_ID);
    //fChain->SetBranchAddress("B_PV_piplus_PE", B_PV_piplus_PE, &b_B_PV_piplus_PE);
    //fChain->SetBranchAddress("B_PV_piplus_PX", B_PV_piplus_PX, &b_B_PV_piplus_PX);
    //fChain->SetBranchAddress("B_PV_piplus_PY", B_PV_piplus_PY, &b_B_PV_piplus_PY);
    //fChain->SetBranchAddress("B_PV_piplus_PZ", B_PV_piplus_PZ, &b_B_PV_piplus_PZ);
    //fChain->SetBranchAddress("B_PV_Kplus_ID", B_PV_Kplus_ID, &b_B_PV_Kplus_ID);
    //fChain->SetBranchAddress("B_PV_Kplus_PE", B_PV_Kplus_PE, &b_B_PV_Kplus_PE);
    //fChain->SetBranchAddress("B_PV_Kplus_PX", B_PV_Kplus_PX, &b_B_PV_Kplus_PX);
    //fChain->SetBranchAddress("B_PV_Kplus_PY", B_PV_Kplus_PY, &b_B_PV_Kplus_PY);
    //fChain->SetBranchAddress("B_PV_Kplus_PZ", B_PV_Kplus_PZ, &b_B_PV_Kplus_PZ);
    fChain->SetBranchAddress("B_PV_nPV", &B_PV_nPV, &b_B_PV_nPV);
   fChain->SetBranchAddress("B_PV_Dplus_Kplus_ID", B_PV_Dplus_Kplus_ID, &b_B_PV_Dplus_Kplus_ID);
   fChain->SetBranchAddress("B_PV_Dplus_Kplus_PE", B_PV_Dplus_Kplus_PE, &b_B_PV_Dplus_Kplus_PE);
   fChain->SetBranchAddress("B_PV_Dplus_Kplus_PX", B_PV_Dplus_Kplus_PX, &b_B_PV_Dplus_Kplus_PX);
   fChain->SetBranchAddress("B_PV_Dplus_Kplus_PY", B_PV_Dplus_Kplus_PY, &b_B_PV_Dplus_Kplus_PY);
   fChain->SetBranchAddress("B_PV_Dplus_Kplus_PZ", B_PV_Dplus_Kplus_PZ, &b_B_PV_Dplus_Kplus_PZ);
   fChain->SetBranchAddress("B_PV_Dplus_M", B_PV_Dplus_M, &b_B_PV_Dplus_M);
   fChain->SetBranchAddress("B_PV_Dplus_MERR", B_PV_Dplus_MERR, &b_B_PV_Dplus_MERR);
   fChain->SetBranchAddress("B_PV_Dplus_P", B_PV_Dplus_P, &b_B_PV_Dplus_P);
   fChain->SetBranchAddress("B_PV_Dplus_PERR", B_PV_Dplus_PERR, &b_B_PV_Dplus_PERR);
   fChain->SetBranchAddress("B_PV_Dplus_ctau", B_PV_Dplus_ctau, &b_B_PV_Dplus_ctau);
   fChain->SetBranchAddress("B_PV_Dplus_ctauErr", B_PV_Dplus_ctauErr, &b_B_PV_Dplus_ctauErr);
   fChain->SetBranchAddress("B_PV_Dplus_decayLength", B_PV_Dplus_decayLength, &b_B_PV_Dplus_decayLength);
   fChain->SetBranchAddress("B_PV_Dplus_decayLengthErr", B_PV_Dplus_decayLengthErr, &b_B_PV_Dplus_decayLengthErr);
   fChain->SetBranchAddress("B_PV_Dplus_piplus_0_ID", B_PV_Dplus_piplus_0_ID, &b_B_PV_Dplus_piplus_0_ID);
   fChain->SetBranchAddress("B_PV_Dplus_piplus_0_PE", B_PV_Dplus_piplus_0_PE, &b_B_PV_Dplus_piplus_0_PE);
   fChain->SetBranchAddress("B_PV_Dplus_piplus_0_PX", B_PV_Dplus_piplus_0_PX, &b_B_PV_Dplus_piplus_0_PX);
   fChain->SetBranchAddress("B_PV_Dplus_piplus_0_PY", B_PV_Dplus_piplus_0_PY, &b_B_PV_Dplus_piplus_0_PY);
   fChain->SetBranchAddress("B_PV_Dplus_piplus_0_PZ", B_PV_Dplus_piplus_0_PZ, &b_B_PV_Dplus_piplus_0_PZ);
   fChain->SetBranchAddress("B_PV_Dplus_piplus_ID", B_PV_Dplus_piplus_ID, &b_B_PV_Dplus_piplus_ID);
   fChain->SetBranchAddress("B_PV_Dplus_piplus_PE", B_PV_Dplus_piplus_PE, &b_B_PV_Dplus_piplus_PE);
   fChain->SetBranchAddress("B_PV_Dplus_piplus_PX", B_PV_Dplus_piplus_PX, &b_B_PV_Dplus_piplus_PX);
   fChain->SetBranchAddress("B_PV_Dplus_piplus_PY", B_PV_Dplus_piplus_PY, &b_B_PV_Dplus_piplus_PY);
   fChain->SetBranchAddress("B_PV_Dplus_piplus_PZ", B_PV_Dplus_piplus_PZ, &b_B_PV_Dplus_piplus_PZ);
   fChain->SetBranchAddress("B_PV_KS0_M", B_PV_KS0_M, &b_B_PV_KS0_M);
   fChain->SetBranchAddress("B_PV_KS0_MERR", B_PV_KS0_MERR, &b_B_PV_KS0_MERR);
   fChain->SetBranchAddress("B_PV_KS0_P", B_PV_KS0_P, &b_B_PV_KS0_P);
   fChain->SetBranchAddress("B_PV_KS0_PERR", B_PV_KS0_PERR, &b_B_PV_KS0_PERR);
   fChain->SetBranchAddress("B_PV_KS0_ctau", B_PV_KS0_ctau, &b_B_PV_KS0_ctau);
   fChain->SetBranchAddress("B_PV_KS0_ctauErr", B_PV_KS0_ctauErr, &b_B_PV_KS0_ctauErr);
   fChain->SetBranchAddress("B_PV_KS0_decayLength", B_PV_KS0_decayLength, &b_B_PV_KS0_decayLength);
   fChain->SetBranchAddress("B_PV_KS0_decayLengthErr", B_PV_KS0_decayLengthErr, &b_B_PV_KS0_decayLengthErr);
   fChain->SetBranchAddress("B_PV_KS0_piplus_0_ID", B_PV_KS0_piplus_0_ID, &b_B_PV_KS0_piplus_0_ID);
   fChain->SetBranchAddress("B_PV_KS0_piplus_0_PE", B_PV_KS0_piplus_0_PE, &b_B_PV_KS0_piplus_0_PE);
   fChain->SetBranchAddress("B_PV_KS0_piplus_0_PX", B_PV_KS0_piplus_0_PX, &b_B_PV_KS0_piplus_0_PX);
   fChain->SetBranchAddress("B_PV_KS0_piplus_0_PY", B_PV_KS0_piplus_0_PY, &b_B_PV_KS0_piplus_0_PY);
   fChain->SetBranchAddress("B_PV_KS0_piplus_0_PZ", B_PV_KS0_piplus_0_PZ, &b_B_PV_KS0_piplus_0_PZ);
   fChain->SetBranchAddress("B_PV_KS0_piplus_ID", B_PV_KS0_piplus_ID, &b_B_PV_KS0_piplus_ID);
   fChain->SetBranchAddress("B_PV_KS0_piplus_PE", B_PV_KS0_piplus_PE, &b_B_PV_KS0_piplus_PE);
   fChain->SetBranchAddress("B_PV_KS0_piplus_PX", B_PV_KS0_piplus_PX, &b_B_PV_KS0_piplus_PX);
   fChain->SetBranchAddress("B_PV_KS0_piplus_PY", B_PV_KS0_piplus_PY, &b_B_PV_KS0_piplus_PY);
   fChain->SetBranchAddress("B_PV_KS0_piplus_PZ", B_PV_KS0_piplus_PZ, &b_B_PV_KS0_piplus_PZ);
   fChain->SetBranchAddress("B_PV_M", B_PV_M, &b_B_PV_M);
   fChain->SetBranchAddress("B_PV_MERR", B_PV_MERR, &b_B_PV_MERR);
   fChain->SetBranchAddress("B_PV_P", B_PV_P, &b_B_PV_P);
   fChain->SetBranchAddress("B_PV_PERR", B_PV_PERR, &b_B_PV_PERR);
   fChain->SetBranchAddress("B_PV_PV_X", B_PV_PV_X, &b_B_PV_PV_X);
   fChain->SetBranchAddress("B_PV_PV_Y", B_PV_PV_Y, &b_B_PV_PV_Y);
   fChain->SetBranchAddress("B_PV_PV_Z", B_PV_PV_Z, &b_B_PV_PV_Z);
   fChain->SetBranchAddress("B_PV_PV_key", B_PV_PV_key, &b_B_PV_PV_key);
   fChain->SetBranchAddress("B_PV_chi2", B_PV_chi2, &b_B_PV_chi2);
   fChain->SetBranchAddress("B_PV_ctau", B_PV_ctau, &b_B_PV_ctau);
   fChain->SetBranchAddress("B_PV_ctauErr", B_PV_ctauErr, &b_B_PV_ctauErr);
   fChain->SetBranchAddress("B_PV_decayLength", B_PV_decayLength, &b_B_PV_decayLength);
   fChain->SetBranchAddress("B_PV_decayLengthErr", B_PV_decayLengthErr, &b_B_PV_decayLengthErr);
   fChain->SetBranchAddress("B_PV_nDOF", B_PV_nDOF, &b_B_PV_nDOF);
   fChain->SetBranchAddress("B_PV_nIter", B_PV_nIter, &b_B_PV_nIter);
   fChain->SetBranchAddress("B_PV_status", B_PV_status, &b_B_PV_status);
   fChain->SetBranchAddress("B_L0Global_Dec", &B_L0Global_Dec, &b_B_L0Global_Dec);
   fChain->SetBranchAddress("B_L0Global_TIS", &B_L0Global_TIS, &b_B_L0Global_TIS);
   fChain->SetBranchAddress("B_L0Global_TOS", &B_L0Global_TOS, &b_B_L0Global_TOS);
   fChain->SetBranchAddress("B_Hlt1Global_Dec", &B_Hlt1Global_Dec, &b_B_Hlt1Global_Dec);
   fChain->SetBranchAddress("B_Hlt1Global_TIS", &B_Hlt1Global_TIS, &b_B_Hlt1Global_TIS);
   fChain->SetBranchAddress("B_Hlt1Global_TOS", &B_Hlt1Global_TOS, &b_B_Hlt1Global_TOS);
   fChain->SetBranchAddress("B_Hlt1Phys_Dec", &B_Hlt1Phys_Dec, &b_B_Hlt1Phys_Dec);
   fChain->SetBranchAddress("B_Hlt1Phys_TIS", &B_Hlt1Phys_TIS, &b_B_Hlt1Phys_TIS);
   fChain->SetBranchAddress("B_Hlt1Phys_TOS", &B_Hlt1Phys_TOS, &b_B_Hlt1Phys_TOS);
   fChain->SetBranchAddress("B_Hlt2Global_Dec", &B_Hlt2Global_Dec, &b_B_Hlt2Global_Dec);
   fChain->SetBranchAddress("B_Hlt2Global_TIS", &B_Hlt2Global_TIS, &b_B_Hlt2Global_TIS);
   fChain->SetBranchAddress("B_Hlt2Global_TOS", &B_Hlt2Global_TOS, &b_B_Hlt2Global_TOS);
   fChain->SetBranchAddress("B_Hlt2Phys_Dec", &B_Hlt2Phys_Dec, &b_B_Hlt2Phys_Dec);
   fChain->SetBranchAddress("B_Hlt2Phys_TIS", &B_Hlt2Phys_TIS, &b_B_Hlt2Phys_TIS);
   fChain->SetBranchAddress("B_Hlt2Phys_TOS", &B_Hlt2Phys_TOS, &b_B_Hlt2Phys_TOS);
   fChain->SetBranchAddress("B_L0HadronDecision_Dec", &B_L0HadronDecision_Dec, &b_B_L0HadronDecision_Dec);
   fChain->SetBranchAddress("B_L0HadronDecision_TIS", &B_L0HadronDecision_TIS, &b_B_L0HadronDecision_TIS);
   fChain->SetBranchAddress("B_L0HadronDecision_TOS", &B_L0HadronDecision_TOS, &b_B_L0HadronDecision_TOS);
   fChain->SetBranchAddress("B_L0MuonDecision_Dec", &B_L0MuonDecision_Dec, &b_B_L0MuonDecision_Dec);
   fChain->SetBranchAddress("B_L0MuonDecision_TIS", &B_L0MuonDecision_TIS, &b_B_L0MuonDecision_TIS);
   fChain->SetBranchAddress("B_L0MuonDecision_TOS", &B_L0MuonDecision_TOS, &b_B_L0MuonDecision_TOS);
   fChain->SetBranchAddress("B_L0DiMuonDecision_Dec", &B_L0DiMuonDecision_Dec, &b_B_L0DiMuonDecision_Dec);
   fChain->SetBranchAddress("B_L0DiMuonDecision_TIS", &B_L0DiMuonDecision_TIS, &b_B_L0DiMuonDecision_TIS);
   fChain->SetBranchAddress("B_L0DiMuonDecision_TOS", &B_L0DiMuonDecision_TOS, &b_B_L0DiMuonDecision_TOS);
   fChain->SetBranchAddress("B_L0ElectronDecision_Dec", &B_L0ElectronDecision_Dec, &b_B_L0ElectronDecision_Dec);
   fChain->SetBranchAddress("B_L0ElectronDecision_TIS", &B_L0ElectronDecision_TIS, &b_B_L0ElectronDecision_TIS);
   fChain->SetBranchAddress("B_L0ElectronDecision_TOS", &B_L0ElectronDecision_TOS, &b_B_L0ElectronDecision_TOS);
   fChain->SetBranchAddress("B_L0PhotonDecision_Dec", &B_L0PhotonDecision_Dec, &b_B_L0PhotonDecision_Dec);
   fChain->SetBranchAddress("B_L0PhotonDecision_TIS", &B_L0PhotonDecision_TIS, &b_B_L0PhotonDecision_TIS);
   fChain->SetBranchAddress("B_L0PhotonDecision_TOS", &B_L0PhotonDecision_TOS, &b_B_L0PhotonDecision_TOS);
   if(_year < 15){
       fChain->SetBranchAddress("B_Hlt1TrackAllL0Decision_TOS", &B_Hlt1TrackAllL0Decision_TOS, &b_B_Hlt1TrackAllL0Decision_TOS); 
       fChain->SetBranchAddress("B_Hlt2Topo2BodyBBDTDecision_TOS", &B_Hlt2Topo2BodyBBDTDecision_TOS, &b_B_Hlt2Topo2BodyBBDTDecision_TOS);
       fChain->SetBranchAddress("B_Hlt2Topo3BodyBBDTDecision_TOS", &B_Hlt2Topo3BodyBBDTDecision_TOS, &b_B_Hlt2Topo3BodyBBDTDecision_TOS);
       fChain->SetBranchAddress("B_Hlt2Topo4BodyBBDTDecision_TOS", &B_Hlt2Topo4BodyBBDTDecision_TOS, &b_B_Hlt2Topo4BodyBBDTDecision_TOS);
   }
   else{
       fChain->SetBranchAddress("B_Hlt1TrackMVADecision_Dec", &B_Hlt1TrackMVADecision_Dec, &b_B_Hlt1TrackMVADecision_Dec);
       fChain->SetBranchAddress("B_Hlt1TrackMVADecision_TIS", &B_Hlt1TrackMVADecision_TIS, &b_B_Hlt1TrackMVADecision_TIS);
       fChain->SetBranchAddress("B_Hlt1TrackMVADecision_TOS", &B_Hlt1TrackMVADecision_TOS, &b_B_Hlt1TrackMVADecision_TOS);
       fChain->SetBranchAddress("B_Hlt1TwoTrackMVADecision_Dec", &B_Hlt1TwoTrackMVADecision_Dec, &b_B_Hlt1TwoTrackMVADecision_Dec);
       fChain->SetBranchAddress("B_Hlt1TwoTrackMVADecision_TIS", &B_Hlt1TwoTrackMVADecision_TIS, &b_B_Hlt1TwoTrackMVADecision_TIS);
       fChain->SetBranchAddress("B_Hlt1TwoTrackMVADecision_TOS", &B_Hlt1TwoTrackMVADecision_TOS, &b_B_Hlt1TwoTrackMVADecision_TOS);
       fChain->SetBranchAddress("B_Hlt2Topo2BodyDecision_Dec", &B_Hlt2Topo2BodyDecision_Dec, &b_B_Hlt2Topo2BodyDecision_Dec);
       fChain->SetBranchAddress("B_Hlt2Topo2BodyDecision_TIS", &B_Hlt2Topo2BodyDecision_TIS, &b_B_Hlt2Topo2BodyDecision_TIS);
       fChain->SetBranchAddress("B_Hlt2Topo2BodyDecision_TOS", &B_Hlt2Topo2BodyDecision_TOS, &b_B_Hlt2Topo2BodyDecision_TOS);
       fChain->SetBranchAddress("B_Hlt2Topo3BodyDecision_Dec", &B_Hlt2Topo3BodyDecision_Dec, &b_B_Hlt2Topo3BodyDecision_Dec);
       fChain->SetBranchAddress("B_Hlt2Topo3BodyDecision_TIS", &B_Hlt2Topo3BodyDecision_TIS, &b_B_Hlt2Topo3BodyDecision_TIS);
       fChain->SetBranchAddress("B_Hlt2Topo3BodyDecision_TOS", &B_Hlt2Topo3BodyDecision_TOS, &b_B_Hlt2Topo3BodyDecision_TOS);
       fChain->SetBranchAddress("B_Hlt2Topo4BodyDecision_Dec", &B_Hlt2Topo4BodyDecision_Dec, &b_B_Hlt2Topo4BodyDecision_Dec);
       fChain->SetBranchAddress("B_Hlt2Topo4BodyDecision_TIS", &B_Hlt2Topo4BodyDecision_TIS, &b_B_Hlt2Topo4BodyDecision_TIS);
       fChain->SetBranchAddress("B_Hlt2Topo4BodyDecision_TOS", &B_Hlt2Topo4BodyDecision_TOS, &b_B_Hlt2Topo4BodyDecision_TOS);
       fChain->SetBranchAddress("B_Hlt2PhiIncPhiDecision_Dec", &B_Hlt2PhiIncPhiDecision_Dec, &b_B_Hlt2PhiIncPhiDecision_Dec);
       fChain->SetBranchAddress("B_Hlt2PhiIncPhiDecision_TIS", &B_Hlt2PhiIncPhiDecision_TIS, &b_B_Hlt2PhiIncPhiDecision_TIS);
       fChain->SetBranchAddress("B_Hlt2PhiIncPhiDecision_TOS", &B_Hlt2PhiIncPhiDecision_TOS, &b_B_Hlt2PhiIncPhiDecision_TOS);
   }    
   fChain->SetBranchAddress("B_Hlt2IncPhiDecision_Dec", &B_Hlt2IncPhiDecision_Dec, &b_B_Hlt2IncPhiDecision_Dec);
   fChain->SetBranchAddress("B_Hlt2IncPhiDecision_TIS", &B_Hlt2IncPhiDecision_TIS, &b_B_Hlt2IncPhiDecision_TIS);
   fChain->SetBranchAddress("B_Hlt2IncPhiDecision_TOS", &B_Hlt2IncPhiDecision_TOS, &b_B_Hlt2IncPhiDecision_TOS);
   fChain->SetBranchAddress("D_ETA", &D_ETA, &b_D_ETA);
   fChain->SetBranchAddress("D_MINIP", &D_MINIP, &b_D_MINIP);
   fChain->SetBranchAddress("D_MINIPCHI2", &D_MINIPCHI2, &b_D_MINIPCHI2);
   fChain->SetBranchAddress("D_MINIPNEXTBEST", &D_MINIPNEXTBEST, &b_D_MINIPNEXTBEST);
   fChain->SetBranchAddress("D_MINIPCHI2NEXTBEST", &D_MINIPCHI2NEXTBEST, &b_D_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("D_ENDVERTEX_X", &D_ENDVERTEX_X, &b_D_ENDVERTEX_X);
   fChain->SetBranchAddress("D_ENDVERTEX_Y", &D_ENDVERTEX_Y, &b_D_ENDVERTEX_Y);
   fChain->SetBranchAddress("D_ENDVERTEX_Z", &D_ENDVERTEX_Z, &b_D_ENDVERTEX_Z);
   fChain->SetBranchAddress("D_ENDVERTEX_XERR", &D_ENDVERTEX_XERR, &b_D_ENDVERTEX_XERR);
   fChain->SetBranchAddress("D_ENDVERTEX_YERR", &D_ENDVERTEX_YERR, &b_D_ENDVERTEX_YERR);
   fChain->SetBranchAddress("D_ENDVERTEX_ZERR", &D_ENDVERTEX_ZERR, &b_D_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("D_ENDVERTEX_CHI2", &D_ENDVERTEX_CHI2, &b_D_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("D_ENDVERTEX_NDOF", &D_ENDVERTEX_NDOF, &b_D_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("D_ENDVERTEX_COV_", D_ENDVERTEX_COV_, &b_D_ENDVERTEX_COV_);
   fChain->SetBranchAddress("D_OWNPV_X", &D_OWNPV_X, &b_D_OWNPV_X);
   fChain->SetBranchAddress("D_OWNPV_Y", &D_OWNPV_Y, &b_D_OWNPV_Y);
   fChain->SetBranchAddress("D_OWNPV_Z", &D_OWNPV_Z, &b_D_OWNPV_Z);
   fChain->SetBranchAddress("D_OWNPV_XERR", &D_OWNPV_XERR, &b_D_OWNPV_XERR);
   fChain->SetBranchAddress("D_OWNPV_YERR", &D_OWNPV_YERR, &b_D_OWNPV_YERR);
   fChain->SetBranchAddress("D_OWNPV_ZERR", &D_OWNPV_ZERR, &b_D_OWNPV_ZERR);
   fChain->SetBranchAddress("D_OWNPV_CHI2", &D_OWNPV_CHI2, &b_D_OWNPV_CHI2);
   fChain->SetBranchAddress("D_OWNPV_NDOF", &D_OWNPV_NDOF, &b_D_OWNPV_NDOF);
   fChain->SetBranchAddress("D_OWNPV_COV_", D_OWNPV_COV_, &b_D_OWNPV_COV_);
   fChain->SetBranchAddress("D_IP_OWNPV", &D_IP_OWNPV, &b_D_IP_OWNPV);
   fChain->SetBranchAddress("D_IPCHI2_OWNPV", &D_IPCHI2_OWNPV, &b_D_IPCHI2_OWNPV);
   fChain->SetBranchAddress("D_FD_OWNPV", &D_FD_OWNPV, &b_D_FD_OWNPV);
   fChain->SetBranchAddress("D_FDCHI2_OWNPV", &D_FDCHI2_OWNPV, &b_D_FDCHI2_OWNPV);
   fChain->SetBranchAddress("D_DIRA_OWNPV", &D_DIRA_OWNPV, &b_D_DIRA_OWNPV);
   fChain->SetBranchAddress("D_TOPPV_X", &D_TOPPV_X, &b_D_TOPPV_X);
   fChain->SetBranchAddress("D_TOPPV_Y", &D_TOPPV_Y, &b_D_TOPPV_Y);
   fChain->SetBranchAddress("D_TOPPV_Z", &D_TOPPV_Z, &b_D_TOPPV_Z);
   fChain->SetBranchAddress("D_TOPPV_XERR", &D_TOPPV_XERR, &b_D_TOPPV_XERR);
   fChain->SetBranchAddress("D_TOPPV_YERR", &D_TOPPV_YERR, &b_D_TOPPV_YERR);
   fChain->SetBranchAddress("D_TOPPV_ZERR", &D_TOPPV_ZERR, &b_D_TOPPV_ZERR);
   fChain->SetBranchAddress("D_TOPPV_CHI2", &D_TOPPV_CHI2, &b_D_TOPPV_CHI2);
   fChain->SetBranchAddress("D_TOPPV_NDOF", &D_TOPPV_NDOF, &b_D_TOPPV_NDOF);
   fChain->SetBranchAddress("D_TOPPV_COV_", D_TOPPV_COV_, &b_D_TOPPV_COV_);
   fChain->SetBranchAddress("D_IP_TOPPV", &D_IP_TOPPV, &b_D_IP_TOPPV);
   fChain->SetBranchAddress("D_IPCHI2_TOPPV", &D_IPCHI2_TOPPV, &b_D_IPCHI2_TOPPV);
   fChain->SetBranchAddress("D_FD_TOPPV", &D_FD_TOPPV, &b_D_FD_TOPPV);
   fChain->SetBranchAddress("D_FDCHI2_TOPPV", &D_FDCHI2_TOPPV, &b_D_FDCHI2_TOPPV);
   fChain->SetBranchAddress("D_DIRA_TOPPV", &D_DIRA_TOPPV, &b_D_DIRA_TOPPV);
   fChain->SetBranchAddress("D_ORIVX_X", &D_ORIVX_X, &b_D_ORIVX_X);
   fChain->SetBranchAddress("D_ORIVX_Y", &D_ORIVX_Y, &b_D_ORIVX_Y);
   fChain->SetBranchAddress("D_ORIVX_Z", &D_ORIVX_Z, &b_D_ORIVX_Z);
   fChain->SetBranchAddress("D_ORIVX_XERR", &D_ORIVX_XERR, &b_D_ORIVX_XERR);
   fChain->SetBranchAddress("D_ORIVX_YERR", &D_ORIVX_YERR, &b_D_ORIVX_YERR);
   fChain->SetBranchAddress("D_ORIVX_ZERR", &D_ORIVX_ZERR, &b_D_ORIVX_ZERR);
   fChain->SetBranchAddress("D_ORIVX_CHI2", &D_ORIVX_CHI2, &b_D_ORIVX_CHI2);
   fChain->SetBranchAddress("D_ORIVX_NDOF", &D_ORIVX_NDOF, &b_D_ORIVX_NDOF);
   fChain->SetBranchAddress("D_ORIVX_COV_", D_ORIVX_COV_, &b_D_ORIVX_COV_);
   fChain->SetBranchAddress("D_IP_ORIVX", &D_IP_ORIVX, &b_D_IP_ORIVX);
   fChain->SetBranchAddress("D_IPCHI2_ORIVX", &D_IPCHI2_ORIVX, &b_D_IPCHI2_ORIVX);
   fChain->SetBranchAddress("D_FD_ORIVX", &D_FD_ORIVX, &b_D_FD_ORIVX);
   fChain->SetBranchAddress("D_FDCHI2_ORIVX", &D_FDCHI2_ORIVX, &b_D_FDCHI2_ORIVX);
   fChain->SetBranchAddress("D_DIRA_ORIVX", &D_DIRA_ORIVX, &b_D_DIRA_ORIVX);
   fChain->SetBranchAddress("D_P", &D_P, &b_D_P);
   fChain->SetBranchAddress("D_PT", &D_PT, &b_D_PT);
   fChain->SetBranchAddress("D_PE", &D_PE, &b_D_PE);
   fChain->SetBranchAddress("D_PX", &D_PX, &b_D_PX);
   fChain->SetBranchAddress("D_PY", &D_PY, &b_D_PY);
   fChain->SetBranchAddress("D_PZ", &D_PZ, &b_D_PZ);
   fChain->SetBranchAddress("D_MM", &D_MM, &b_D_MM);
   fChain->SetBranchAddress("D_MMERR", &D_MMERR, &b_D_MMERR);
   fChain->SetBranchAddress("D_M", &D_M, &b_D_M);
   fChain->SetBranchAddress("D_ID", &D_ID, &b_D_ID);
   fChain->SetBranchAddress("D_DOCA1", &D_DOCA1, &b_D_DOCA1);
   fChain->SetBranchAddress("D_TAU", &D_TAU, &b_D_TAU);
   fChain->SetBranchAddress("D_TAUERR", &D_TAUERR, &b_D_TAUERR);
   fChain->SetBranchAddress("K_D_ETA", &K_D_ETA, &b_K_D_ETA);
   fChain->SetBranchAddress("K_D_MC12TuneV2_ProbNNe", &K_D_MC12TuneV2_ProbNNe, &b_K_D_MC12TuneV2_ProbNNe);
   fChain->SetBranchAddress("K_D_MC12TuneV2_ProbNNmu", &K_D_MC12TuneV2_ProbNNmu, &b_K_D_MC12TuneV2_ProbNNmu);
   fChain->SetBranchAddress("K_D_MC12TuneV2_ProbNNpi", &K_D_MC12TuneV2_ProbNNpi, &b_K_D_MC12TuneV2_ProbNNpi);
   fChain->SetBranchAddress("K_D_MC12TuneV2_ProbNNk", &K_D_MC12TuneV2_ProbNNk, &b_K_D_MC12TuneV2_ProbNNk);
   fChain->SetBranchAddress("K_D_MC12TuneV2_ProbNNp", &K_D_MC12TuneV2_ProbNNp, &b_K_D_MC12TuneV2_ProbNNp);
   fChain->SetBranchAddress("K_D_MC12TuneV2_ProbNNghost", &K_D_MC12TuneV2_ProbNNghost, &b_K_D_MC12TuneV2_ProbNNghost);
   fChain->SetBranchAddress("K_D_MC12TuneV3_ProbNNe", &K_D_MC12TuneV3_ProbNNe, &b_K_D_MC12TuneV3_ProbNNe);
   fChain->SetBranchAddress("K_D_MC12TuneV3_ProbNNmu", &K_D_MC12TuneV3_ProbNNmu, &b_K_D_MC12TuneV3_ProbNNmu);
   fChain->SetBranchAddress("K_D_MC12TuneV3_ProbNNpi", &K_D_MC12TuneV3_ProbNNpi, &b_K_D_MC12TuneV3_ProbNNpi);
   fChain->SetBranchAddress("K_D_MC12TuneV3_ProbNNk", &K_D_MC12TuneV3_ProbNNk, &b_K_D_MC12TuneV3_ProbNNk);
   fChain->SetBranchAddress("K_D_MC12TuneV3_ProbNNp", &K_D_MC12TuneV3_ProbNNp, &b_K_D_MC12TuneV3_ProbNNp);
   fChain->SetBranchAddress("K_D_MC12TuneV3_ProbNNghost", &K_D_MC12TuneV3_ProbNNghost, &b_K_D_MC12TuneV3_ProbNNghost);
   fChain->SetBranchAddress("K_D_MC12TuneV4_ProbNNe", &K_D_MC12TuneV4_ProbNNe, &b_K_D_MC12TuneV4_ProbNNe);
   fChain->SetBranchAddress("K_D_MC12TuneV4_ProbNNmu", &K_D_MC12TuneV4_ProbNNmu, &b_K_D_MC12TuneV4_ProbNNmu);
   fChain->SetBranchAddress("K_D_MC12TuneV4_ProbNNpi", &K_D_MC12TuneV4_ProbNNpi, &b_K_D_MC12TuneV4_ProbNNpi);
   fChain->SetBranchAddress("K_D_MC12TuneV4_ProbNNk", &K_D_MC12TuneV4_ProbNNk, &b_K_D_MC12TuneV4_ProbNNk);
   fChain->SetBranchAddress("K_D_MC12TuneV4_ProbNNp", &K_D_MC12TuneV4_ProbNNp, &b_K_D_MC12TuneV4_ProbNNp);
   fChain->SetBranchAddress("K_D_MC12TuneV4_ProbNNghost", &K_D_MC12TuneV4_ProbNNghost, &b_K_D_MC12TuneV4_ProbNNghost);
   fChain->SetBranchAddress("K_D_MC15TuneV1_ProbNNe", &K_D_MC15TuneV1_ProbNNe, &b_K_D_MC15TuneV1_ProbNNe);
   fChain->SetBranchAddress("K_D_MC15TuneV1_ProbNNmu", &K_D_MC15TuneV1_ProbNNmu, &b_K_D_MC15TuneV1_ProbNNmu);
   fChain->SetBranchAddress("K_D_MC15TuneV1_ProbNNpi", &K_D_MC15TuneV1_ProbNNpi, &b_K_D_MC15TuneV1_ProbNNpi);
   fChain->SetBranchAddress("K_D_MC15TuneV1_ProbNNk", &K_D_MC15TuneV1_ProbNNk, &b_K_D_MC15TuneV1_ProbNNk);
   fChain->SetBranchAddress("K_D_MC15TuneV1_ProbNNp", &K_D_MC15TuneV1_ProbNNp, &b_K_D_MC15TuneV1_ProbNNp);
   fChain->SetBranchAddress("K_D_MC15TuneV1_ProbNNghost", &K_D_MC15TuneV1_ProbNNghost, &b_K_D_MC15TuneV1_ProbNNghost);
   fChain->SetBranchAddress("K_D_MINIP", &K_D_MINIP, &b_K_D_MINIP);
   fChain->SetBranchAddress("K_D_MINIPCHI2", &K_D_MINIPCHI2, &b_K_D_MINIPCHI2);
   fChain->SetBranchAddress("K_D_MINIPNEXTBEST", &K_D_MINIPNEXTBEST, &b_K_D_MINIPNEXTBEST);
   fChain->SetBranchAddress("K_D_MINIPCHI2NEXTBEST", &K_D_MINIPCHI2NEXTBEST, &b_K_D_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("K_D_OWNPV_X", &K_D_OWNPV_X, &b_K_D_OWNPV_X);
   fChain->SetBranchAddress("K_D_OWNPV_Y", &K_D_OWNPV_Y, &b_K_D_OWNPV_Y);
   fChain->SetBranchAddress("K_D_OWNPV_Z", &K_D_OWNPV_Z, &b_K_D_OWNPV_Z);
   fChain->SetBranchAddress("K_D_OWNPV_XERR", &K_D_OWNPV_XERR, &b_K_D_OWNPV_XERR);
   fChain->SetBranchAddress("K_D_OWNPV_YERR", &K_D_OWNPV_YERR, &b_K_D_OWNPV_YERR);
   fChain->SetBranchAddress("K_D_OWNPV_ZERR", &K_D_OWNPV_ZERR, &b_K_D_OWNPV_ZERR);
   fChain->SetBranchAddress("K_D_OWNPV_CHI2", &K_D_OWNPV_CHI2, &b_K_D_OWNPV_CHI2);
   fChain->SetBranchAddress("K_D_OWNPV_NDOF", &K_D_OWNPV_NDOF, &b_K_D_OWNPV_NDOF);
   fChain->SetBranchAddress("K_D_OWNPV_COV_", K_D_OWNPV_COV_, &b_K_D_OWNPV_COV_);
   fChain->SetBranchAddress("K_D_IP_OWNPV", &K_D_IP_OWNPV, &b_K_D_IP_OWNPV);
   fChain->SetBranchAddress("K_D_IPCHI2_OWNPV", &K_D_IPCHI2_OWNPV, &b_K_D_IPCHI2_OWNPV);
   fChain->SetBranchAddress("K_D_TOPPV_X", &K_D_TOPPV_X, &b_K_D_TOPPV_X);
   fChain->SetBranchAddress("K_D_TOPPV_Y", &K_D_TOPPV_Y, &b_K_D_TOPPV_Y);
   fChain->SetBranchAddress("K_D_TOPPV_Z", &K_D_TOPPV_Z, &b_K_D_TOPPV_Z);
   fChain->SetBranchAddress("K_D_TOPPV_XERR", &K_D_TOPPV_XERR, &b_K_D_TOPPV_XERR);
   fChain->SetBranchAddress("K_D_TOPPV_YERR", &K_D_TOPPV_YERR, &b_K_D_TOPPV_YERR);
   fChain->SetBranchAddress("K_D_TOPPV_ZERR", &K_D_TOPPV_ZERR, &b_K_D_TOPPV_ZERR);
   fChain->SetBranchAddress("K_D_TOPPV_CHI2", &K_D_TOPPV_CHI2, &b_K_D_TOPPV_CHI2);
   fChain->SetBranchAddress("K_D_TOPPV_NDOF", &K_D_TOPPV_NDOF, &b_K_D_TOPPV_NDOF);
   fChain->SetBranchAddress("K_D_TOPPV_COV_", K_D_TOPPV_COV_, &b_K_D_TOPPV_COV_);
   fChain->SetBranchAddress("K_D_IP_TOPPV", &K_D_IP_TOPPV, &b_K_D_IP_TOPPV);
   fChain->SetBranchAddress("K_D_IPCHI2_TOPPV", &K_D_IPCHI2_TOPPV, &b_K_D_IPCHI2_TOPPV);
   fChain->SetBranchAddress("K_D_ORIVX_X", &K_D_ORIVX_X, &b_K_D_ORIVX_X);
   fChain->SetBranchAddress("K_D_ORIVX_Y", &K_D_ORIVX_Y, &b_K_D_ORIVX_Y);
   fChain->SetBranchAddress("K_D_ORIVX_Z", &K_D_ORIVX_Z, &b_K_D_ORIVX_Z);
   fChain->SetBranchAddress("K_D_ORIVX_XERR", &K_D_ORIVX_XERR, &b_K_D_ORIVX_XERR);
   fChain->SetBranchAddress("K_D_ORIVX_YERR", &K_D_ORIVX_YERR, &b_K_D_ORIVX_YERR);
   fChain->SetBranchAddress("K_D_ORIVX_ZERR", &K_D_ORIVX_ZERR, &b_K_D_ORIVX_ZERR);
   fChain->SetBranchAddress("K_D_ORIVX_CHI2", &K_D_ORIVX_CHI2, &b_K_D_ORIVX_CHI2);
   fChain->SetBranchAddress("K_D_ORIVX_NDOF", &K_D_ORIVX_NDOF, &b_K_D_ORIVX_NDOF);
   fChain->SetBranchAddress("K_D_ORIVX_COV_", K_D_ORIVX_COV_, &b_K_D_ORIVX_COV_);
   fChain->SetBranchAddress("K_D_IP_ORIVX", &K_D_IP_ORIVX, &b_K_D_IP_ORIVX);
   fChain->SetBranchAddress("K_D_IPCHI2_ORIVX", &K_D_IPCHI2_ORIVX, &b_K_D_IPCHI2_ORIVX);
   fChain->SetBranchAddress("K_D_P", &K_D_P, &b_K_D_P);
   fChain->SetBranchAddress("K_D_PT", &K_D_PT, &b_K_D_PT);
   fChain->SetBranchAddress("K_D_PE", &K_D_PE, &b_K_D_PE);
   fChain->SetBranchAddress("K_D_PX", &K_D_PX, &b_K_D_PX);
   fChain->SetBranchAddress("K_D_PY", &K_D_PY, &b_K_D_PY);
   fChain->SetBranchAddress("K_D_PZ", &K_D_PZ, &b_K_D_PZ);
   fChain->SetBranchAddress("K_D_M", &K_D_M, &b_K_D_M);
   fChain->SetBranchAddress("K_D_ID", &K_D_ID, &b_K_D_ID);
   fChain->SetBranchAddress("K_D_PIDe", &K_D_PIDe, &b_K_D_PIDe);
   fChain->SetBranchAddress("K_D_PIDmu", &K_D_PIDmu, &b_K_D_PIDmu);
   fChain->SetBranchAddress("K_D_PIDK", &K_D_PIDK, &b_K_D_PIDK);
   fChain->SetBranchAddress("K_D_PIDp", &K_D_PIDp, &b_K_D_PIDp);
   fChain->SetBranchAddress("K_D_ProbNNe", &K_D_ProbNNe, &b_K_D_ProbNNe);
   fChain->SetBranchAddress("K_D_ProbNNk", &K_D_ProbNNk, &b_K_D_ProbNNk);
   fChain->SetBranchAddress("K_D_ProbNNp", &K_D_ProbNNp, &b_K_D_ProbNNp);
   fChain->SetBranchAddress("K_D_ProbNNpi", &K_D_ProbNNpi, &b_K_D_ProbNNpi);
   fChain->SetBranchAddress("K_D_ProbNNmu", &K_D_ProbNNmu, &b_K_D_ProbNNmu);
   fChain->SetBranchAddress("K_D_ProbNNghost", &K_D_ProbNNghost, &b_K_D_ProbNNghost);
   fChain->SetBranchAddress("K_D_hasMuon", &K_D_hasMuon, &b_K_D_hasMuon);
   fChain->SetBranchAddress("K_D_isMuon", &K_D_isMuon, &b_K_D_isMuon);
   fChain->SetBranchAddress("K_D_hasRich", &K_D_hasRich, &b_K_D_hasRich);
   fChain->SetBranchAddress("K_D_UsedRichAerogel", &K_D_UsedRichAerogel, &b_K_D_UsedRichAerogel);
   fChain->SetBranchAddress("K_D_UsedRich1Gas", &K_D_UsedRich1Gas, &b_K_D_UsedRich1Gas);
   fChain->SetBranchAddress("K_D_UsedRich2Gas", &K_D_UsedRich2Gas, &b_K_D_UsedRich2Gas);
   fChain->SetBranchAddress("K_D_RichAboveElThres", &K_D_RichAboveElThres, &b_K_D_RichAboveElThres);
   fChain->SetBranchAddress("K_D_RichAboveMuThres", &K_D_RichAboveMuThres, &b_K_D_RichAboveMuThres);
   fChain->SetBranchAddress("K_D_RichAbovePiThres", &K_D_RichAbovePiThres, &b_K_D_RichAbovePiThres);
   fChain->SetBranchAddress("K_D_RichAboveKaThres", &K_D_RichAboveKaThres, &b_K_D_RichAboveKaThres);
   fChain->SetBranchAddress("K_D_RichAbovePrThres", &K_D_RichAbovePrThres, &b_K_D_RichAbovePrThres);
   fChain->SetBranchAddress("K_D_hasCalo", &K_D_hasCalo, &b_K_D_hasCalo);
   fChain->SetBranchAddress("K_D_TRACK_Type", &K_D_TRACK_Type, &b_K_D_TRACK_Type);
   fChain->SetBranchAddress("K_D_TRACK_Key", &K_D_TRACK_Key, &b_K_D_TRACK_Key);
   fChain->SetBranchAddress("K_D_TRACK_CHI2NDOF", &K_D_TRACK_CHI2NDOF, &b_K_D_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("K_D_TRACK_PCHI2", &K_D_TRACK_PCHI2, &b_K_D_TRACK_PCHI2);
   fChain->SetBranchAddress("K_D_TRACK_MatchCHI2", &K_D_TRACK_MatchCHI2, &b_K_D_TRACK_MatchCHI2);
   fChain->SetBranchAddress("K_D_TRACK_GhostProb", &K_D_TRACK_GhostProb, &b_K_D_TRACK_GhostProb);
   fChain->SetBranchAddress("K_D_TRACK_CloneDist", &K_D_TRACK_CloneDist, &b_K_D_TRACK_CloneDist);
   fChain->SetBranchAddress("K_D_TRACK_Likelihood", &K_D_TRACK_Likelihood, &b_K_D_TRACK_Likelihood);
   fChain->SetBranchAddress("pi1_D_ETA", &pi1_D_ETA, &b_pi1_D_ETA);
   fChain->SetBranchAddress("pi1_D_MC12TuneV2_ProbNNe", &pi1_D_MC12TuneV2_ProbNNe, &b_pi1_D_MC12TuneV2_ProbNNe);
   fChain->SetBranchAddress("pi1_D_MC12TuneV2_ProbNNmu", &pi1_D_MC12TuneV2_ProbNNmu, &b_pi1_D_MC12TuneV2_ProbNNmu);
   fChain->SetBranchAddress("pi1_D_MC12TuneV2_ProbNNpi", &pi1_D_MC12TuneV2_ProbNNpi, &b_pi1_D_MC12TuneV2_ProbNNpi);
   fChain->SetBranchAddress("pi1_D_MC12TuneV2_ProbNNk", &pi1_D_MC12TuneV2_ProbNNk, &b_pi1_D_MC12TuneV2_ProbNNk);
   fChain->SetBranchAddress("pi1_D_MC12TuneV2_ProbNNp", &pi1_D_MC12TuneV2_ProbNNp, &b_pi1_D_MC12TuneV2_ProbNNp);
   fChain->SetBranchAddress("pi1_D_MC12TuneV2_ProbNNghost", &pi1_D_MC12TuneV2_ProbNNghost, &b_pi1_D_MC12TuneV2_ProbNNghost);
   fChain->SetBranchAddress("pi1_D_MC12TuneV3_ProbNNe", &pi1_D_MC12TuneV3_ProbNNe, &b_pi1_D_MC12TuneV3_ProbNNe);
   fChain->SetBranchAddress("pi1_D_MC12TuneV3_ProbNNmu", &pi1_D_MC12TuneV3_ProbNNmu, &b_pi1_D_MC12TuneV3_ProbNNmu);
   fChain->SetBranchAddress("pi1_D_MC12TuneV3_ProbNNpi", &pi1_D_MC12TuneV3_ProbNNpi, &b_pi1_D_MC12TuneV3_ProbNNpi);
   fChain->SetBranchAddress("pi1_D_MC12TuneV3_ProbNNk", &pi1_D_MC12TuneV3_ProbNNk, &b_pi1_D_MC12TuneV3_ProbNNk);
   fChain->SetBranchAddress("pi1_D_MC12TuneV3_ProbNNp", &pi1_D_MC12TuneV3_ProbNNp, &b_pi1_D_MC12TuneV3_ProbNNp);
   fChain->SetBranchAddress("pi1_D_MC12TuneV3_ProbNNghost", &pi1_D_MC12TuneV3_ProbNNghost, &b_pi1_D_MC12TuneV3_ProbNNghost);
   fChain->SetBranchAddress("pi1_D_MC12TuneV4_ProbNNe", &pi1_D_MC12TuneV4_ProbNNe, &b_pi1_D_MC12TuneV4_ProbNNe);
   fChain->SetBranchAddress("pi1_D_MC12TuneV4_ProbNNmu", &pi1_D_MC12TuneV4_ProbNNmu, &b_pi1_D_MC12TuneV4_ProbNNmu);
   fChain->SetBranchAddress("pi1_D_MC12TuneV4_ProbNNpi", &pi1_D_MC12TuneV4_ProbNNpi, &b_pi1_D_MC12TuneV4_ProbNNpi);
   fChain->SetBranchAddress("pi1_D_MC12TuneV4_ProbNNk", &pi1_D_MC12TuneV4_ProbNNk, &b_pi1_D_MC12TuneV4_ProbNNk);
   fChain->SetBranchAddress("pi1_D_MC12TuneV4_ProbNNp", &pi1_D_MC12TuneV4_ProbNNp, &b_pi1_D_MC12TuneV4_ProbNNp);
   fChain->SetBranchAddress("pi1_D_MC12TuneV4_ProbNNghost", &pi1_D_MC12TuneV4_ProbNNghost, &b_pi1_D_MC12TuneV4_ProbNNghost);
   fChain->SetBranchAddress("pi1_D_MC15TuneV1_ProbNNe", &pi1_D_MC15TuneV1_ProbNNe, &b_pi1_D_MC15TuneV1_ProbNNe);
   fChain->SetBranchAddress("pi1_D_MC15TuneV1_ProbNNmu", &pi1_D_MC15TuneV1_ProbNNmu, &b_pi1_D_MC15TuneV1_ProbNNmu);
   fChain->SetBranchAddress("pi1_D_MC15TuneV1_ProbNNpi", &pi1_D_MC15TuneV1_ProbNNpi, &b_pi1_D_MC15TuneV1_ProbNNpi);
   fChain->SetBranchAddress("pi1_D_MC15TuneV1_ProbNNk", &pi1_D_MC15TuneV1_ProbNNk, &b_pi1_D_MC15TuneV1_ProbNNk);
   fChain->SetBranchAddress("pi1_D_MC15TuneV1_ProbNNp", &pi1_D_MC15TuneV1_ProbNNp, &b_pi1_D_MC15TuneV1_ProbNNp);
   fChain->SetBranchAddress("pi1_D_MC15TuneV1_ProbNNghost", &pi1_D_MC15TuneV1_ProbNNghost, &b_pi1_D_MC15TuneV1_ProbNNghost);
   fChain->SetBranchAddress("pi1_D_MINIP", &pi1_D_MINIP, &b_pi1_D_MINIP);
   fChain->SetBranchAddress("pi1_D_MINIPCHI2", &pi1_D_MINIPCHI2, &b_pi1_D_MINIPCHI2);
   fChain->SetBranchAddress("pi1_D_MINIPNEXTBEST", &pi1_D_MINIPNEXTBEST, &b_pi1_D_MINIPNEXTBEST);
   fChain->SetBranchAddress("pi1_D_MINIPCHI2NEXTBEST", &pi1_D_MINIPCHI2NEXTBEST, &b_pi1_D_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("pi1_D_OWNPV_X", &pi1_D_OWNPV_X, &b_pi1_D_OWNPV_X);
   fChain->SetBranchAddress("pi1_D_OWNPV_Y", &pi1_D_OWNPV_Y, &b_pi1_D_OWNPV_Y);
   fChain->SetBranchAddress("pi1_D_OWNPV_Z", &pi1_D_OWNPV_Z, &b_pi1_D_OWNPV_Z);
   fChain->SetBranchAddress("pi1_D_OWNPV_XERR", &pi1_D_OWNPV_XERR, &b_pi1_D_OWNPV_XERR);
   fChain->SetBranchAddress("pi1_D_OWNPV_YERR", &pi1_D_OWNPV_YERR, &b_pi1_D_OWNPV_YERR);
   fChain->SetBranchAddress("pi1_D_OWNPV_ZERR", &pi1_D_OWNPV_ZERR, &b_pi1_D_OWNPV_ZERR);
   fChain->SetBranchAddress("pi1_D_OWNPV_CHI2", &pi1_D_OWNPV_CHI2, &b_pi1_D_OWNPV_CHI2);
   fChain->SetBranchAddress("pi1_D_OWNPV_NDOF", &pi1_D_OWNPV_NDOF, &b_pi1_D_OWNPV_NDOF);
   fChain->SetBranchAddress("pi1_D_OWNPV_COV_", pi1_D_OWNPV_COV_, &b_pi1_D_OWNPV_COV_);
   fChain->SetBranchAddress("pi1_D_IP_OWNPV", &pi1_D_IP_OWNPV, &b_pi1_D_IP_OWNPV);
   fChain->SetBranchAddress("pi1_D_IPCHI2_OWNPV", &pi1_D_IPCHI2_OWNPV, &b_pi1_D_IPCHI2_OWNPV);
   fChain->SetBranchAddress("pi1_D_TOPPV_X", &pi1_D_TOPPV_X, &b_pi1_D_TOPPV_X);
   fChain->SetBranchAddress("pi1_D_TOPPV_Y", &pi1_D_TOPPV_Y, &b_pi1_D_TOPPV_Y);
   fChain->SetBranchAddress("pi1_D_TOPPV_Z", &pi1_D_TOPPV_Z, &b_pi1_D_TOPPV_Z);
   fChain->SetBranchAddress("pi1_D_TOPPV_XERR", &pi1_D_TOPPV_XERR, &b_pi1_D_TOPPV_XERR);
   fChain->SetBranchAddress("pi1_D_TOPPV_YERR", &pi1_D_TOPPV_YERR, &b_pi1_D_TOPPV_YERR);
   fChain->SetBranchAddress("pi1_D_TOPPV_ZERR", &pi1_D_TOPPV_ZERR, &b_pi1_D_TOPPV_ZERR);
   fChain->SetBranchAddress("pi1_D_TOPPV_CHI2", &pi1_D_TOPPV_CHI2, &b_pi1_D_TOPPV_CHI2);
   fChain->SetBranchAddress("pi1_D_TOPPV_NDOF", &pi1_D_TOPPV_NDOF, &b_pi1_D_TOPPV_NDOF);
   fChain->SetBranchAddress("pi1_D_TOPPV_COV_", pi1_D_TOPPV_COV_, &b_pi1_D_TOPPV_COV_);
   fChain->SetBranchAddress("pi1_D_IP_TOPPV", &pi1_D_IP_TOPPV, &b_pi1_D_IP_TOPPV);
   fChain->SetBranchAddress("pi1_D_IPCHI2_TOPPV", &pi1_D_IPCHI2_TOPPV, &b_pi1_D_IPCHI2_TOPPV);
   fChain->SetBranchAddress("pi1_D_ORIVX_X", &pi1_D_ORIVX_X, &b_pi1_D_ORIVX_X);
   fChain->SetBranchAddress("pi1_D_ORIVX_Y", &pi1_D_ORIVX_Y, &b_pi1_D_ORIVX_Y);
   fChain->SetBranchAddress("pi1_D_ORIVX_Z", &pi1_D_ORIVX_Z, &b_pi1_D_ORIVX_Z);
   fChain->SetBranchAddress("pi1_D_ORIVX_XERR", &pi1_D_ORIVX_XERR, &b_pi1_D_ORIVX_XERR);
   fChain->SetBranchAddress("pi1_D_ORIVX_YERR", &pi1_D_ORIVX_YERR, &b_pi1_D_ORIVX_YERR);
   fChain->SetBranchAddress("pi1_D_ORIVX_ZERR", &pi1_D_ORIVX_ZERR, &b_pi1_D_ORIVX_ZERR);
   fChain->SetBranchAddress("pi1_D_ORIVX_CHI2", &pi1_D_ORIVX_CHI2, &b_pi1_D_ORIVX_CHI2);
   fChain->SetBranchAddress("pi1_D_ORIVX_NDOF", &pi1_D_ORIVX_NDOF, &b_pi1_D_ORIVX_NDOF);
   fChain->SetBranchAddress("pi1_D_ORIVX_COV_", pi1_D_ORIVX_COV_, &b_pi1_D_ORIVX_COV_);
   fChain->SetBranchAddress("pi1_D_IP_ORIVX", &pi1_D_IP_ORIVX, &b_pi1_D_IP_ORIVX);
   fChain->SetBranchAddress("pi1_D_IPCHI2_ORIVX", &pi1_D_IPCHI2_ORIVX, &b_pi1_D_IPCHI2_ORIVX);
   fChain->SetBranchAddress("pi1_D_P", &pi1_D_P, &b_pi1_D_P);
   fChain->SetBranchAddress("pi1_D_PT", &pi1_D_PT, &b_pi1_D_PT);
   fChain->SetBranchAddress("pi1_D_PE", &pi1_D_PE, &b_pi1_D_PE);
   fChain->SetBranchAddress("pi1_D_PX", &pi1_D_PX, &b_pi1_D_PX);
   fChain->SetBranchAddress("pi1_D_PY", &pi1_D_PY, &b_pi1_D_PY);
   fChain->SetBranchAddress("pi1_D_PZ", &pi1_D_PZ, &b_pi1_D_PZ);
   fChain->SetBranchAddress("pi1_D_M", &pi1_D_M, &b_pi1_D_M);
   fChain->SetBranchAddress("pi1_D_ID", &pi1_D_ID, &b_pi1_D_ID);
   fChain->SetBranchAddress("pi1_D_PIDe", &pi1_D_PIDe, &b_pi1_D_PIDe);
   fChain->SetBranchAddress("pi1_D_PIDmu", &pi1_D_PIDmu, &b_pi1_D_PIDmu);
   fChain->SetBranchAddress("pi1_D_PIDK", &pi1_D_PIDK, &b_pi1_D_PIDK);
   fChain->SetBranchAddress("pi1_D_PIDp", &pi1_D_PIDp, &b_pi1_D_PIDp);
   fChain->SetBranchAddress("pi1_D_ProbNNe", &pi1_D_ProbNNe, &b_pi1_D_ProbNNe);
   fChain->SetBranchAddress("pi1_D_ProbNNk", &pi1_D_ProbNNk, &b_pi1_D_ProbNNk);
   fChain->SetBranchAddress("pi1_D_ProbNNp", &pi1_D_ProbNNp, &b_pi1_D_ProbNNp);
   fChain->SetBranchAddress("pi1_D_ProbNNpi", &pi1_D_ProbNNpi, &b_pi1_D_ProbNNpi);
   fChain->SetBranchAddress("pi1_D_ProbNNmu", &pi1_D_ProbNNmu, &b_pi1_D_ProbNNmu);
   fChain->SetBranchAddress("pi1_D_ProbNNghost", &pi1_D_ProbNNghost, &b_pi1_D_ProbNNghost);
   fChain->SetBranchAddress("pi1_D_hasMuon", &pi1_D_hasMuon, &b_pi1_D_hasMuon);
   fChain->SetBranchAddress("pi1_D_isMuon", &pi1_D_isMuon, &b_pi1_D_isMuon);
   fChain->SetBranchAddress("pi1_D_hasRich", &pi1_D_hasRich, &b_pi1_D_hasRich);
   fChain->SetBranchAddress("pi1_D_UsedRichAerogel", &pi1_D_UsedRichAerogel, &b_pi1_D_UsedRichAerogel);
   fChain->SetBranchAddress("pi1_D_UsedRich1Gas", &pi1_D_UsedRich1Gas, &b_pi1_D_UsedRich1Gas);
   fChain->SetBranchAddress("pi1_D_UsedRich2Gas", &pi1_D_UsedRich2Gas, &b_pi1_D_UsedRich2Gas);
   fChain->SetBranchAddress("pi1_D_RichAboveElThres", &pi1_D_RichAboveElThres, &b_pi1_D_RichAboveElThres);
   fChain->SetBranchAddress("pi1_D_RichAboveMuThres", &pi1_D_RichAboveMuThres, &b_pi1_D_RichAboveMuThres);
   fChain->SetBranchAddress("pi1_D_RichAbovePiThres", &pi1_D_RichAbovePiThres, &b_pi1_D_RichAbovePiThres);
   fChain->SetBranchAddress("pi1_D_RichAboveKaThres", &pi1_D_RichAboveKaThres, &b_pi1_D_RichAboveKaThres);
   fChain->SetBranchAddress("pi1_D_RichAbovePrThres", &pi1_D_RichAbovePrThres, &b_pi1_D_RichAbovePrThres);
   fChain->SetBranchAddress("pi1_D_hasCalo", &pi1_D_hasCalo, &b_pi1_D_hasCalo);
   fChain->SetBranchAddress("pi1_D_TRACK_Type", &pi1_D_TRACK_Type, &b_pi1_D_TRACK_Type);
   fChain->SetBranchAddress("pi1_D_TRACK_Key", &pi1_D_TRACK_Key, &b_pi1_D_TRACK_Key);
   fChain->SetBranchAddress("pi1_D_TRACK_CHI2NDOF", &pi1_D_TRACK_CHI2NDOF, &b_pi1_D_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("pi1_D_TRACK_PCHI2", &pi1_D_TRACK_PCHI2, &b_pi1_D_TRACK_PCHI2);
   fChain->SetBranchAddress("pi1_D_TRACK_MatchCHI2", &pi1_D_TRACK_MatchCHI2, &b_pi1_D_TRACK_MatchCHI2);
   fChain->SetBranchAddress("pi1_D_TRACK_GhostProb", &pi1_D_TRACK_GhostProb, &b_pi1_D_TRACK_GhostProb);
   fChain->SetBranchAddress("pi1_D_TRACK_CloneDist", &pi1_D_TRACK_CloneDist, &b_pi1_D_TRACK_CloneDist);
   fChain->SetBranchAddress("pi1_D_TRACK_Likelihood", &pi1_D_TRACK_Likelihood, &b_pi1_D_TRACK_Likelihood);
   fChain->SetBranchAddress("pi2_D_ETA", &pi2_D_ETA, &b_pi2_D_ETA);
   fChain->SetBranchAddress("pi2_D_MC12TuneV2_ProbNNe", &pi2_D_MC12TuneV2_ProbNNe, &b_pi2_D_MC12TuneV2_ProbNNe);
   fChain->SetBranchAddress("pi2_D_MC12TuneV2_ProbNNmu", &pi2_D_MC12TuneV2_ProbNNmu, &b_pi2_D_MC12TuneV2_ProbNNmu);
   fChain->SetBranchAddress("pi2_D_MC12TuneV2_ProbNNpi", &pi2_D_MC12TuneV2_ProbNNpi, &b_pi2_D_MC12TuneV2_ProbNNpi);
   fChain->SetBranchAddress("pi2_D_MC12TuneV2_ProbNNk", &pi2_D_MC12TuneV2_ProbNNk, &b_pi2_D_MC12TuneV2_ProbNNk);
   fChain->SetBranchAddress("pi2_D_MC12TuneV2_ProbNNp", &pi2_D_MC12TuneV2_ProbNNp, &b_pi2_D_MC12TuneV2_ProbNNp);
   fChain->SetBranchAddress("pi2_D_MC12TuneV2_ProbNNghost", &pi2_D_MC12TuneV2_ProbNNghost, &b_pi2_D_MC12TuneV2_ProbNNghost);
   fChain->SetBranchAddress("pi2_D_MC12TuneV3_ProbNNe", &pi2_D_MC12TuneV3_ProbNNe, &b_pi2_D_MC12TuneV3_ProbNNe);
   fChain->SetBranchAddress("pi2_D_MC12TuneV3_ProbNNmu", &pi2_D_MC12TuneV3_ProbNNmu, &b_pi2_D_MC12TuneV3_ProbNNmu);
   fChain->SetBranchAddress("pi2_D_MC12TuneV3_ProbNNpi", &pi2_D_MC12TuneV3_ProbNNpi, &b_pi2_D_MC12TuneV3_ProbNNpi);
   fChain->SetBranchAddress("pi2_D_MC12TuneV3_ProbNNk", &pi2_D_MC12TuneV3_ProbNNk, &b_pi2_D_MC12TuneV3_ProbNNk);
   fChain->SetBranchAddress("pi2_D_MC12TuneV3_ProbNNp", &pi2_D_MC12TuneV3_ProbNNp, &b_pi2_D_MC12TuneV3_ProbNNp);
   fChain->SetBranchAddress("pi2_D_MC12TuneV3_ProbNNghost", &pi2_D_MC12TuneV3_ProbNNghost, &b_pi2_D_MC12TuneV3_ProbNNghost);
   fChain->SetBranchAddress("pi2_D_MC12TuneV4_ProbNNe", &pi2_D_MC12TuneV4_ProbNNe, &b_pi2_D_MC12TuneV4_ProbNNe);
   fChain->SetBranchAddress("pi2_D_MC12TuneV4_ProbNNmu", &pi2_D_MC12TuneV4_ProbNNmu, &b_pi2_D_MC12TuneV4_ProbNNmu);
   fChain->SetBranchAddress("pi2_D_MC12TuneV4_ProbNNpi", &pi2_D_MC12TuneV4_ProbNNpi, &b_pi2_D_MC12TuneV4_ProbNNpi);
   fChain->SetBranchAddress("pi2_D_MC12TuneV4_ProbNNk", &pi2_D_MC12TuneV4_ProbNNk, &b_pi2_D_MC12TuneV4_ProbNNk);
   fChain->SetBranchAddress("pi2_D_MC12TuneV4_ProbNNp", &pi2_D_MC12TuneV4_ProbNNp, &b_pi2_D_MC12TuneV4_ProbNNp);
   fChain->SetBranchAddress("pi2_D_MC12TuneV4_ProbNNghost", &pi2_D_MC12TuneV4_ProbNNghost, &b_pi2_D_MC12TuneV4_ProbNNghost);
   fChain->SetBranchAddress("pi2_D_MC15TuneV1_ProbNNe", &pi2_D_MC15TuneV1_ProbNNe, &b_pi2_D_MC15TuneV1_ProbNNe);
   fChain->SetBranchAddress("pi2_D_MC15TuneV1_ProbNNmu", &pi2_D_MC15TuneV1_ProbNNmu, &b_pi2_D_MC15TuneV1_ProbNNmu);
   fChain->SetBranchAddress("pi2_D_MC15TuneV1_ProbNNpi", &pi2_D_MC15TuneV1_ProbNNpi, &b_pi2_D_MC15TuneV1_ProbNNpi);
   fChain->SetBranchAddress("pi2_D_MC15TuneV1_ProbNNk", &pi2_D_MC15TuneV1_ProbNNk, &b_pi2_D_MC15TuneV1_ProbNNk);
   fChain->SetBranchAddress("pi2_D_MC15TuneV1_ProbNNp", &pi2_D_MC15TuneV1_ProbNNp, &b_pi2_D_MC15TuneV1_ProbNNp);
   fChain->SetBranchAddress("pi2_D_MC15TuneV1_ProbNNghost", &pi2_D_MC15TuneV1_ProbNNghost, &b_pi2_D_MC15TuneV1_ProbNNghost);
   fChain->SetBranchAddress("pi2_D_MINIP", &pi2_D_MINIP, &b_pi2_D_MINIP);
   fChain->SetBranchAddress("pi2_D_MINIPCHI2", &pi2_D_MINIPCHI2, &b_pi2_D_MINIPCHI2);
   fChain->SetBranchAddress("pi2_D_MINIPNEXTBEST", &pi2_D_MINIPNEXTBEST, &b_pi2_D_MINIPNEXTBEST);
   fChain->SetBranchAddress("pi2_D_MINIPCHI2NEXTBEST", &pi2_D_MINIPCHI2NEXTBEST, &b_pi2_D_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("pi2_D_OWNPV_X", &pi2_D_OWNPV_X, &b_pi2_D_OWNPV_X);
   fChain->SetBranchAddress("pi2_D_OWNPV_Y", &pi2_D_OWNPV_Y, &b_pi2_D_OWNPV_Y);
   fChain->SetBranchAddress("pi2_D_OWNPV_Z", &pi2_D_OWNPV_Z, &b_pi2_D_OWNPV_Z);
   fChain->SetBranchAddress("pi2_D_OWNPV_XERR", &pi2_D_OWNPV_XERR, &b_pi2_D_OWNPV_XERR);
   fChain->SetBranchAddress("pi2_D_OWNPV_YERR", &pi2_D_OWNPV_YERR, &b_pi2_D_OWNPV_YERR);
   fChain->SetBranchAddress("pi2_D_OWNPV_ZERR", &pi2_D_OWNPV_ZERR, &b_pi2_D_OWNPV_ZERR);
   fChain->SetBranchAddress("pi2_D_OWNPV_CHI2", &pi2_D_OWNPV_CHI2, &b_pi2_D_OWNPV_CHI2);
   fChain->SetBranchAddress("pi2_D_OWNPV_NDOF", &pi2_D_OWNPV_NDOF, &b_pi2_D_OWNPV_NDOF);
   fChain->SetBranchAddress("pi2_D_OWNPV_COV_", pi2_D_OWNPV_COV_, &b_pi2_D_OWNPV_COV_);
   fChain->SetBranchAddress("pi2_D_IP_OWNPV", &pi2_D_IP_OWNPV, &b_pi2_D_IP_OWNPV);
   fChain->SetBranchAddress("pi2_D_IPCHI2_OWNPV", &pi2_D_IPCHI2_OWNPV, &b_pi2_D_IPCHI2_OWNPV);
   fChain->SetBranchAddress("pi2_D_TOPPV_X", &pi2_D_TOPPV_X, &b_pi2_D_TOPPV_X);
   fChain->SetBranchAddress("pi2_D_TOPPV_Y", &pi2_D_TOPPV_Y, &b_pi2_D_TOPPV_Y);
   fChain->SetBranchAddress("pi2_D_TOPPV_Z", &pi2_D_TOPPV_Z, &b_pi2_D_TOPPV_Z);
   fChain->SetBranchAddress("pi2_D_TOPPV_XERR", &pi2_D_TOPPV_XERR, &b_pi2_D_TOPPV_XERR);
   fChain->SetBranchAddress("pi2_D_TOPPV_YERR", &pi2_D_TOPPV_YERR, &b_pi2_D_TOPPV_YERR);
   fChain->SetBranchAddress("pi2_D_TOPPV_ZERR", &pi2_D_TOPPV_ZERR, &b_pi2_D_TOPPV_ZERR);
   fChain->SetBranchAddress("pi2_D_TOPPV_CHI2", &pi2_D_TOPPV_CHI2, &b_pi2_D_TOPPV_CHI2);
   fChain->SetBranchAddress("pi2_D_TOPPV_NDOF", &pi2_D_TOPPV_NDOF, &b_pi2_D_TOPPV_NDOF);
   fChain->SetBranchAddress("pi2_D_TOPPV_COV_", pi2_D_TOPPV_COV_, &b_pi2_D_TOPPV_COV_);
   fChain->SetBranchAddress("pi2_D_IP_TOPPV", &pi2_D_IP_TOPPV, &b_pi2_D_IP_TOPPV);
   fChain->SetBranchAddress("pi2_D_IPCHI2_TOPPV", &pi2_D_IPCHI2_TOPPV, &b_pi2_D_IPCHI2_TOPPV);
   fChain->SetBranchAddress("pi2_D_ORIVX_X", &pi2_D_ORIVX_X, &b_pi2_D_ORIVX_X);
   fChain->SetBranchAddress("pi2_D_ORIVX_Y", &pi2_D_ORIVX_Y, &b_pi2_D_ORIVX_Y);
   fChain->SetBranchAddress("pi2_D_ORIVX_Z", &pi2_D_ORIVX_Z, &b_pi2_D_ORIVX_Z);
   fChain->SetBranchAddress("pi2_D_ORIVX_XERR", &pi2_D_ORIVX_XERR, &b_pi2_D_ORIVX_XERR);
   fChain->SetBranchAddress("pi2_D_ORIVX_YERR", &pi2_D_ORIVX_YERR, &b_pi2_D_ORIVX_YERR);
   fChain->SetBranchAddress("pi2_D_ORIVX_ZERR", &pi2_D_ORIVX_ZERR, &b_pi2_D_ORIVX_ZERR);
   fChain->SetBranchAddress("pi2_D_ORIVX_CHI2", &pi2_D_ORIVX_CHI2, &b_pi2_D_ORIVX_CHI2);
   fChain->SetBranchAddress("pi2_D_ORIVX_NDOF", &pi2_D_ORIVX_NDOF, &b_pi2_D_ORIVX_NDOF);
   fChain->SetBranchAddress("pi2_D_ORIVX_COV_", pi2_D_ORIVX_COV_, &b_pi2_D_ORIVX_COV_);
   fChain->SetBranchAddress("pi2_D_IP_ORIVX", &pi2_D_IP_ORIVX, &b_pi2_D_IP_ORIVX);
   fChain->SetBranchAddress("pi2_D_IPCHI2_ORIVX", &pi2_D_IPCHI2_ORIVX, &b_pi2_D_IPCHI2_ORIVX);
   fChain->SetBranchAddress("pi2_D_P", &pi2_D_P, &b_pi2_D_P);
   fChain->SetBranchAddress("pi2_D_PT", &pi2_D_PT, &b_pi2_D_PT);
   fChain->SetBranchAddress("pi2_D_PE", &pi2_D_PE, &b_pi2_D_PE);
   fChain->SetBranchAddress("pi2_D_PX", &pi2_D_PX, &b_pi2_D_PX);
   fChain->SetBranchAddress("pi2_D_PY", &pi2_D_PY, &b_pi2_D_PY);
   fChain->SetBranchAddress("pi2_D_PZ", &pi2_D_PZ, &b_pi2_D_PZ);
   fChain->SetBranchAddress("pi2_D_M", &pi2_D_M, &b_pi2_D_M);
   fChain->SetBranchAddress("pi2_D_ID", &pi2_D_ID, &b_pi2_D_ID);
   fChain->SetBranchAddress("pi2_D_PIDe", &pi2_D_PIDe, &b_pi2_D_PIDe);
   fChain->SetBranchAddress("pi2_D_PIDmu", &pi2_D_PIDmu, &b_pi2_D_PIDmu);
   fChain->SetBranchAddress("pi2_D_PIDK", &pi2_D_PIDK, &b_pi2_D_PIDK);
   fChain->SetBranchAddress("pi2_D_PIDp", &pi2_D_PIDp, &b_pi2_D_PIDp);
   fChain->SetBranchAddress("pi2_D_ProbNNe", &pi2_D_ProbNNe, &b_pi2_D_ProbNNe);
   fChain->SetBranchAddress("pi2_D_ProbNNk", &pi2_D_ProbNNk, &b_pi2_D_ProbNNk);
   fChain->SetBranchAddress("pi2_D_ProbNNp", &pi2_D_ProbNNp, &b_pi2_D_ProbNNp);
   fChain->SetBranchAddress("pi2_D_ProbNNpi", &pi2_D_ProbNNpi, &b_pi2_D_ProbNNpi);
   fChain->SetBranchAddress("pi2_D_ProbNNmu", &pi2_D_ProbNNmu, &b_pi2_D_ProbNNmu);
   fChain->SetBranchAddress("pi2_D_ProbNNghost", &pi2_D_ProbNNghost, &b_pi2_D_ProbNNghost);
   fChain->SetBranchAddress("pi2_D_hasMuon", &pi2_D_hasMuon, &b_pi2_D_hasMuon);
   fChain->SetBranchAddress("pi2_D_isMuon", &pi2_D_isMuon, &b_pi2_D_isMuon);
   fChain->SetBranchAddress("pi2_D_hasRich", &pi2_D_hasRich, &b_pi2_D_hasRich);
   fChain->SetBranchAddress("pi2_D_UsedRichAerogel", &pi2_D_UsedRichAerogel, &b_pi2_D_UsedRichAerogel);
   fChain->SetBranchAddress("pi2_D_UsedRich1Gas", &pi2_D_UsedRich1Gas, &b_pi2_D_UsedRich1Gas);
   fChain->SetBranchAddress("pi2_D_UsedRich2Gas", &pi2_D_UsedRich2Gas, &b_pi2_D_UsedRich2Gas);
   fChain->SetBranchAddress("pi2_D_RichAboveElThres", &pi2_D_RichAboveElThres, &b_pi2_D_RichAboveElThres);
   fChain->SetBranchAddress("pi2_D_RichAboveMuThres", &pi2_D_RichAboveMuThres, &b_pi2_D_RichAboveMuThres);
   fChain->SetBranchAddress("pi2_D_RichAbovePiThres", &pi2_D_RichAbovePiThres, &b_pi2_D_RichAbovePiThres);
   fChain->SetBranchAddress("pi2_D_RichAboveKaThres", &pi2_D_RichAboveKaThres, &b_pi2_D_RichAboveKaThres);
   fChain->SetBranchAddress("pi2_D_RichAbovePrThres", &pi2_D_RichAbovePrThres, &b_pi2_D_RichAbovePrThres);
   fChain->SetBranchAddress("pi2_D_hasCalo", &pi2_D_hasCalo, &b_pi2_D_hasCalo);
   fChain->SetBranchAddress("pi2_D_TRACK_Type", &pi2_D_TRACK_Type, &b_pi2_D_TRACK_Type);
   fChain->SetBranchAddress("pi2_D_TRACK_Key", &pi2_D_TRACK_Key, &b_pi2_D_TRACK_Key);
   fChain->SetBranchAddress("pi2_D_TRACK_CHI2NDOF", &pi2_D_TRACK_CHI2NDOF, &b_pi2_D_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("pi2_D_TRACK_PCHI2", &pi2_D_TRACK_PCHI2, &b_pi2_D_TRACK_PCHI2);
   fChain->SetBranchAddress("pi2_D_TRACK_MatchCHI2", &pi2_D_TRACK_MatchCHI2, &b_pi2_D_TRACK_MatchCHI2);
   fChain->SetBranchAddress("pi2_D_TRACK_GhostProb", &pi2_D_TRACK_GhostProb, &b_pi2_D_TRACK_GhostProb);
   fChain->SetBranchAddress("pi2_D_TRACK_CloneDist", &pi2_D_TRACK_CloneDist, &b_pi2_D_TRACK_CloneDist);
   fChain->SetBranchAddress("pi2_D_TRACK_Likelihood", &pi2_D_TRACK_Likelihood, &b_pi2_D_TRACK_Likelihood);
   fChain->SetBranchAddress("Ks_ETA", &Ks_ETA, &b_Ks_ETA);
   fChain->SetBranchAddress("Ks_MINIP", &Ks_MINIP, &b_Ks_MINIP);
   fChain->SetBranchAddress("Ks_MINIPCHI2", &Ks_MINIPCHI2, &b_Ks_MINIPCHI2);
   fChain->SetBranchAddress("Ks_MINIPNEXTBEST", &Ks_MINIPNEXTBEST, &b_Ks_MINIPNEXTBEST);
   fChain->SetBranchAddress("Ks_MINIPCHI2NEXTBEST", &Ks_MINIPCHI2NEXTBEST, &b_Ks_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("Ks_ENDVERTEX_X", &Ks_ENDVERTEX_X, &b_Ks_ENDVERTEX_X);
   fChain->SetBranchAddress("Ks_ENDVERTEX_Y", &Ks_ENDVERTEX_Y, &b_Ks_ENDVERTEX_Y);
   fChain->SetBranchAddress("Ks_ENDVERTEX_Z", &Ks_ENDVERTEX_Z, &b_Ks_ENDVERTEX_Z);
   fChain->SetBranchAddress("Ks_ENDVERTEX_XERR", &Ks_ENDVERTEX_XERR, &b_Ks_ENDVERTEX_XERR);
   fChain->SetBranchAddress("Ks_ENDVERTEX_YERR", &Ks_ENDVERTEX_YERR, &b_Ks_ENDVERTEX_YERR);
   fChain->SetBranchAddress("Ks_ENDVERTEX_ZERR", &Ks_ENDVERTEX_ZERR, &b_Ks_ENDVERTEX_ZERR);
   fChain->SetBranchAddress("Ks_ENDVERTEX_CHI2", &Ks_ENDVERTEX_CHI2, &b_Ks_ENDVERTEX_CHI2);
   fChain->SetBranchAddress("Ks_ENDVERTEX_NDOF", &Ks_ENDVERTEX_NDOF, &b_Ks_ENDVERTEX_NDOF);
   fChain->SetBranchAddress("Ks_ENDVERTEX_COV_", Ks_ENDVERTEX_COV_, &b_Ks_ENDVERTEX_COV_);
   fChain->SetBranchAddress("Ks_OWNPV_X", &Ks_OWNPV_X, &b_Ks_OWNPV_X);
   fChain->SetBranchAddress("Ks_OWNPV_Y", &Ks_OWNPV_Y, &b_Ks_OWNPV_Y);
   fChain->SetBranchAddress("Ks_OWNPV_Z", &Ks_OWNPV_Z, &b_Ks_OWNPV_Z);
   fChain->SetBranchAddress("Ks_OWNPV_XERR", &Ks_OWNPV_XERR, &b_Ks_OWNPV_XERR);
   fChain->SetBranchAddress("Ks_OWNPV_YERR", &Ks_OWNPV_YERR, &b_Ks_OWNPV_YERR);
   fChain->SetBranchAddress("Ks_OWNPV_ZERR", &Ks_OWNPV_ZERR, &b_Ks_OWNPV_ZERR);
   fChain->SetBranchAddress("Ks_OWNPV_CHI2", &Ks_OWNPV_CHI2, &b_Ks_OWNPV_CHI2);
   fChain->SetBranchAddress("Ks_OWNPV_NDOF", &Ks_OWNPV_NDOF, &b_Ks_OWNPV_NDOF);
   fChain->SetBranchAddress("Ks_OWNPV_COV_", Ks_OWNPV_COV_, &b_Ks_OWNPV_COV_);
   fChain->SetBranchAddress("Ks_IP_OWNPV", &Ks_IP_OWNPV, &b_Ks_IP_OWNPV);
   fChain->SetBranchAddress("Ks_IPCHI2_OWNPV", &Ks_IPCHI2_OWNPV, &b_Ks_IPCHI2_OWNPV);
   fChain->SetBranchAddress("Ks_FD_OWNPV", &Ks_FD_OWNPV, &b_Ks_FD_OWNPV);
   fChain->SetBranchAddress("Ks_FDCHI2_OWNPV", &Ks_FDCHI2_OWNPV, &b_Ks_FDCHI2_OWNPV);
   fChain->SetBranchAddress("Ks_DIRA_OWNPV", &Ks_DIRA_OWNPV, &b_Ks_DIRA_OWNPV);
   fChain->SetBranchAddress("Ks_TOPPV_X", &Ks_TOPPV_X, &b_Ks_TOPPV_X);
   fChain->SetBranchAddress("Ks_TOPPV_Y", &Ks_TOPPV_Y, &b_Ks_TOPPV_Y);
   fChain->SetBranchAddress("Ks_TOPPV_Z", &Ks_TOPPV_Z, &b_Ks_TOPPV_Z);
   fChain->SetBranchAddress("Ks_TOPPV_XERR", &Ks_TOPPV_XERR, &b_Ks_TOPPV_XERR);
   fChain->SetBranchAddress("Ks_TOPPV_YERR", &Ks_TOPPV_YERR, &b_Ks_TOPPV_YERR);
   fChain->SetBranchAddress("Ks_TOPPV_ZERR", &Ks_TOPPV_ZERR, &b_Ks_TOPPV_ZERR);
   fChain->SetBranchAddress("Ks_TOPPV_CHI2", &Ks_TOPPV_CHI2, &b_Ks_TOPPV_CHI2);
   fChain->SetBranchAddress("Ks_TOPPV_NDOF", &Ks_TOPPV_NDOF, &b_Ks_TOPPV_NDOF);
   fChain->SetBranchAddress("Ks_TOPPV_COV_", Ks_TOPPV_COV_, &b_Ks_TOPPV_COV_);
   fChain->SetBranchAddress("Ks_IP_TOPPV", &Ks_IP_TOPPV, &b_Ks_IP_TOPPV);
   fChain->SetBranchAddress("Ks_IPCHI2_TOPPV", &Ks_IPCHI2_TOPPV, &b_Ks_IPCHI2_TOPPV);
   fChain->SetBranchAddress("Ks_FD_TOPPV", &Ks_FD_TOPPV, &b_Ks_FD_TOPPV);
   fChain->SetBranchAddress("Ks_FDCHI2_TOPPV", &Ks_FDCHI2_TOPPV, &b_Ks_FDCHI2_TOPPV);
   fChain->SetBranchAddress("Ks_DIRA_TOPPV", &Ks_DIRA_TOPPV, &b_Ks_DIRA_TOPPV);
   fChain->SetBranchAddress("Ks_ORIVX_X", &Ks_ORIVX_X, &b_Ks_ORIVX_X);
   fChain->SetBranchAddress("Ks_ORIVX_Y", &Ks_ORIVX_Y, &b_Ks_ORIVX_Y);
   fChain->SetBranchAddress("Ks_ORIVX_Z", &Ks_ORIVX_Z, &b_Ks_ORIVX_Z);
   fChain->SetBranchAddress("Ks_ORIVX_XERR", &Ks_ORIVX_XERR, &b_Ks_ORIVX_XERR);
   fChain->SetBranchAddress("Ks_ORIVX_YERR", &Ks_ORIVX_YERR, &b_Ks_ORIVX_YERR);
   fChain->SetBranchAddress("Ks_ORIVX_ZERR", &Ks_ORIVX_ZERR, &b_Ks_ORIVX_ZERR);
   fChain->SetBranchAddress("Ks_ORIVX_CHI2", &Ks_ORIVX_CHI2, &b_Ks_ORIVX_CHI2);
   fChain->SetBranchAddress("Ks_ORIVX_NDOF", &Ks_ORIVX_NDOF, &b_Ks_ORIVX_NDOF);
   fChain->SetBranchAddress("Ks_ORIVX_COV_", Ks_ORIVX_COV_, &b_Ks_ORIVX_COV_);
   fChain->SetBranchAddress("Ks_IP_ORIVX", &Ks_IP_ORIVX, &b_Ks_IP_ORIVX);
   fChain->SetBranchAddress("Ks_IPCHI2_ORIVX", &Ks_IPCHI2_ORIVX, &b_Ks_IPCHI2_ORIVX);
   fChain->SetBranchAddress("Ks_FD_ORIVX", &Ks_FD_ORIVX, &b_Ks_FD_ORIVX);
   fChain->SetBranchAddress("Ks_FDCHI2_ORIVX", &Ks_FDCHI2_ORIVX, &b_Ks_FDCHI2_ORIVX);
   fChain->SetBranchAddress("Ks_DIRA_ORIVX", &Ks_DIRA_ORIVX, &b_Ks_DIRA_ORIVX);
   fChain->SetBranchAddress("Ks_P", &Ks_P, &b_Ks_P);
   fChain->SetBranchAddress("Ks_PT", &Ks_PT, &b_Ks_PT);
   fChain->SetBranchAddress("Ks_PE", &Ks_PE, &b_Ks_PE);
   fChain->SetBranchAddress("Ks_PX", &Ks_PX, &b_Ks_PX);
   fChain->SetBranchAddress("Ks_PY", &Ks_PY, &b_Ks_PY);
   fChain->SetBranchAddress("Ks_PZ", &Ks_PZ, &b_Ks_PZ);
   fChain->SetBranchAddress("Ks_MM", &Ks_MM, &b_Ks_MM);
   fChain->SetBranchAddress("Ks_MMERR", &Ks_MMERR, &b_Ks_MMERR);
   fChain->SetBranchAddress("Ks_M", &Ks_M, &b_Ks_M);
   fChain->SetBranchAddress("Ks_ID", &Ks_ID, &b_Ks_ID);
   fChain->SetBranchAddress("Ks_DOCA1", &Ks_DOCA1, &b_Ks_DOCA1);
   fChain->SetBranchAddress("Ks_TAU", &Ks_TAU, &b_Ks_TAU);
   fChain->SetBranchAddress("Ks_TAUERR", &Ks_TAUERR, &b_Ks_TAUERR);
   fChain->SetBranchAddress("pip_Ks_ETA", &pip_Ks_ETA, &b_pip_Ks_ETA);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV2_ProbNNe", &pip_Ks_MC12TuneV2_ProbNNe, &b_pip_Ks_MC12TuneV2_ProbNNe);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV2_ProbNNmu", &pip_Ks_MC12TuneV2_ProbNNmu, &b_pip_Ks_MC12TuneV2_ProbNNmu);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV2_ProbNNpi", &pip_Ks_MC12TuneV2_ProbNNpi, &b_pip_Ks_MC12TuneV2_ProbNNpi);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV2_ProbNNk", &pip_Ks_MC12TuneV2_ProbNNk, &b_pip_Ks_MC12TuneV2_ProbNNk);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV2_ProbNNp", &pip_Ks_MC12TuneV2_ProbNNp, &b_pip_Ks_MC12TuneV2_ProbNNp);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV2_ProbNNghost", &pip_Ks_MC12TuneV2_ProbNNghost, &b_pip_Ks_MC12TuneV2_ProbNNghost);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV3_ProbNNe", &pip_Ks_MC12TuneV3_ProbNNe, &b_pip_Ks_MC12TuneV3_ProbNNe);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV3_ProbNNmu", &pip_Ks_MC12TuneV3_ProbNNmu, &b_pip_Ks_MC12TuneV3_ProbNNmu);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV3_ProbNNpi", &pip_Ks_MC12TuneV3_ProbNNpi, &b_pip_Ks_MC12TuneV3_ProbNNpi);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV3_ProbNNk", &pip_Ks_MC12TuneV3_ProbNNk, &b_pip_Ks_MC12TuneV3_ProbNNk);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV3_ProbNNp", &pip_Ks_MC12TuneV3_ProbNNp, &b_pip_Ks_MC12TuneV3_ProbNNp);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV3_ProbNNghost", &pip_Ks_MC12TuneV3_ProbNNghost, &b_pip_Ks_MC12TuneV3_ProbNNghost);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV4_ProbNNe", &pip_Ks_MC12TuneV4_ProbNNe, &b_pip_Ks_MC12TuneV4_ProbNNe);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV4_ProbNNmu", &pip_Ks_MC12TuneV4_ProbNNmu, &b_pip_Ks_MC12TuneV4_ProbNNmu);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV4_ProbNNpi", &pip_Ks_MC12TuneV4_ProbNNpi, &b_pip_Ks_MC12TuneV4_ProbNNpi);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV4_ProbNNk", &pip_Ks_MC12TuneV4_ProbNNk, &b_pip_Ks_MC12TuneV4_ProbNNk);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV4_ProbNNp", &pip_Ks_MC12TuneV4_ProbNNp, &b_pip_Ks_MC12TuneV4_ProbNNp);
   fChain->SetBranchAddress("pip_Ks_MC12TuneV4_ProbNNghost", &pip_Ks_MC12TuneV4_ProbNNghost, &b_pip_Ks_MC12TuneV4_ProbNNghost);
   fChain->SetBranchAddress("pip_Ks_MC15TuneV1_ProbNNe", &pip_Ks_MC15TuneV1_ProbNNe, &b_pip_Ks_MC15TuneV1_ProbNNe);
   fChain->SetBranchAddress("pip_Ks_MC15TuneV1_ProbNNmu", &pip_Ks_MC15TuneV1_ProbNNmu, &b_pip_Ks_MC15TuneV1_ProbNNmu);
   fChain->SetBranchAddress("pip_Ks_MC15TuneV1_ProbNNpi", &pip_Ks_MC15TuneV1_ProbNNpi, &b_pip_Ks_MC15TuneV1_ProbNNpi);
   fChain->SetBranchAddress("pip_Ks_MC15TuneV1_ProbNNk", &pip_Ks_MC15TuneV1_ProbNNk, &b_pip_Ks_MC15TuneV1_ProbNNk);
   fChain->SetBranchAddress("pip_Ks_MC15TuneV1_ProbNNp", &pip_Ks_MC15TuneV1_ProbNNp, &b_pip_Ks_MC15TuneV1_ProbNNp);
   fChain->SetBranchAddress("pip_Ks_MC15TuneV1_ProbNNghost", &pip_Ks_MC15TuneV1_ProbNNghost, &b_pip_Ks_MC15TuneV1_ProbNNghost);
   fChain->SetBranchAddress("pip_Ks_MINIP", &pip_Ks_MINIP, &b_pip_Ks_MINIP);
   fChain->SetBranchAddress("pip_Ks_MINIPCHI2", &pip_Ks_MINIPCHI2, &b_pip_Ks_MINIPCHI2);
   fChain->SetBranchAddress("pip_Ks_MINIPNEXTBEST", &pip_Ks_MINIPNEXTBEST, &b_pip_Ks_MINIPNEXTBEST);
   fChain->SetBranchAddress("pip_Ks_MINIPCHI2NEXTBEST", &pip_Ks_MINIPCHI2NEXTBEST, &b_pip_Ks_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("pip_Ks_OWNPV_X", &pip_Ks_OWNPV_X, &b_pip_Ks_OWNPV_X);
   fChain->SetBranchAddress("pip_Ks_OWNPV_Y", &pip_Ks_OWNPV_Y, &b_pip_Ks_OWNPV_Y);
   fChain->SetBranchAddress("pip_Ks_OWNPV_Z", &pip_Ks_OWNPV_Z, &b_pip_Ks_OWNPV_Z);
   fChain->SetBranchAddress("pip_Ks_OWNPV_XERR", &pip_Ks_OWNPV_XERR, &b_pip_Ks_OWNPV_XERR);
   fChain->SetBranchAddress("pip_Ks_OWNPV_YERR", &pip_Ks_OWNPV_YERR, &b_pip_Ks_OWNPV_YERR);
   fChain->SetBranchAddress("pip_Ks_OWNPV_ZERR", &pip_Ks_OWNPV_ZERR, &b_pip_Ks_OWNPV_ZERR);
   fChain->SetBranchAddress("pip_Ks_OWNPV_CHI2", &pip_Ks_OWNPV_CHI2, &b_pip_Ks_OWNPV_CHI2);
   fChain->SetBranchAddress("pip_Ks_OWNPV_NDOF", &pip_Ks_OWNPV_NDOF, &b_pip_Ks_OWNPV_NDOF);
   fChain->SetBranchAddress("pip_Ks_OWNPV_COV_", pip_Ks_OWNPV_COV_, &b_pip_Ks_OWNPV_COV_);
   fChain->SetBranchAddress("pip_Ks_IP_OWNPV", &pip_Ks_IP_OWNPV, &b_pip_Ks_IP_OWNPV);
   fChain->SetBranchAddress("pip_Ks_IPCHI2_OWNPV", &pip_Ks_IPCHI2_OWNPV, &b_pip_Ks_IPCHI2_OWNPV);
   fChain->SetBranchAddress("pip_Ks_TOPPV_X", &pip_Ks_TOPPV_X, &b_pip_Ks_TOPPV_X);
   fChain->SetBranchAddress("pip_Ks_TOPPV_Y", &pip_Ks_TOPPV_Y, &b_pip_Ks_TOPPV_Y);
   fChain->SetBranchAddress("pip_Ks_TOPPV_Z", &pip_Ks_TOPPV_Z, &b_pip_Ks_TOPPV_Z);
   fChain->SetBranchAddress("pip_Ks_TOPPV_XERR", &pip_Ks_TOPPV_XERR, &b_pip_Ks_TOPPV_XERR);
   fChain->SetBranchAddress("pip_Ks_TOPPV_YERR", &pip_Ks_TOPPV_YERR, &b_pip_Ks_TOPPV_YERR);
   fChain->SetBranchAddress("pip_Ks_TOPPV_ZERR", &pip_Ks_TOPPV_ZERR, &b_pip_Ks_TOPPV_ZERR);
   fChain->SetBranchAddress("pip_Ks_TOPPV_CHI2", &pip_Ks_TOPPV_CHI2, &b_pip_Ks_TOPPV_CHI2);
   fChain->SetBranchAddress("pip_Ks_TOPPV_NDOF", &pip_Ks_TOPPV_NDOF, &b_pip_Ks_TOPPV_NDOF);
   fChain->SetBranchAddress("pip_Ks_TOPPV_COV_", pip_Ks_TOPPV_COV_, &b_pip_Ks_TOPPV_COV_);
   fChain->SetBranchAddress("pip_Ks_IP_TOPPV", &pip_Ks_IP_TOPPV, &b_pip_Ks_IP_TOPPV);
   fChain->SetBranchAddress("pip_Ks_IPCHI2_TOPPV", &pip_Ks_IPCHI2_TOPPV, &b_pip_Ks_IPCHI2_TOPPV);
   fChain->SetBranchAddress("pip_Ks_ORIVX_X", &pip_Ks_ORIVX_X, &b_pip_Ks_ORIVX_X);
   fChain->SetBranchAddress("pip_Ks_ORIVX_Y", &pip_Ks_ORIVX_Y, &b_pip_Ks_ORIVX_Y);
   fChain->SetBranchAddress("pip_Ks_ORIVX_Z", &pip_Ks_ORIVX_Z, &b_pip_Ks_ORIVX_Z);
   fChain->SetBranchAddress("pip_Ks_ORIVX_XERR", &pip_Ks_ORIVX_XERR, &b_pip_Ks_ORIVX_XERR);
   fChain->SetBranchAddress("pip_Ks_ORIVX_YERR", &pip_Ks_ORIVX_YERR, &b_pip_Ks_ORIVX_YERR);
   fChain->SetBranchAddress("pip_Ks_ORIVX_ZERR", &pip_Ks_ORIVX_ZERR, &b_pip_Ks_ORIVX_ZERR);
   fChain->SetBranchAddress("pip_Ks_ORIVX_CHI2", &pip_Ks_ORIVX_CHI2, &b_pip_Ks_ORIVX_CHI2);
   fChain->SetBranchAddress("pip_Ks_ORIVX_NDOF", &pip_Ks_ORIVX_NDOF, &b_pip_Ks_ORIVX_NDOF);
   fChain->SetBranchAddress("pip_Ks_ORIVX_COV_", pip_Ks_ORIVX_COV_, &b_pip_Ks_ORIVX_COV_);
   fChain->SetBranchAddress("pip_Ks_IP_ORIVX", &pip_Ks_IP_ORIVX, &b_pip_Ks_IP_ORIVX);
   fChain->SetBranchAddress("pip_Ks_IPCHI2_ORIVX", &pip_Ks_IPCHI2_ORIVX, &b_pip_Ks_IPCHI2_ORIVX);
   fChain->SetBranchAddress("pip_Ks_P", &pip_Ks_P, &b_pip_Ks_P);
   fChain->SetBranchAddress("pip_Ks_PT", &pip_Ks_PT, &b_pip_Ks_PT);
   fChain->SetBranchAddress("pip_Ks_PE", &pip_Ks_PE, &b_pip_Ks_PE);
   fChain->SetBranchAddress("pip_Ks_PX", &pip_Ks_PX, &b_pip_Ks_PX);
   fChain->SetBranchAddress("pip_Ks_PY", &pip_Ks_PY, &b_pip_Ks_PY);
   fChain->SetBranchAddress("pip_Ks_PZ", &pip_Ks_PZ, &b_pip_Ks_PZ);
   fChain->SetBranchAddress("pip_Ks_M", &pip_Ks_M, &b_pip_Ks_M);
   fChain->SetBranchAddress("pip_Ks_ID", &pip_Ks_ID, &b_pip_Ks_ID);
   fChain->SetBranchAddress("pip_Ks_PIDe", &pip_Ks_PIDe, &b_pip_Ks_PIDe);
   fChain->SetBranchAddress("pip_Ks_PIDmu", &pip_Ks_PIDmu, &b_pip_Ks_PIDmu);
   fChain->SetBranchAddress("pip_Ks_PIDK", &pip_Ks_PIDK, &b_pip_Ks_PIDK);
   fChain->SetBranchAddress("pip_Ks_PIDp", &pip_Ks_PIDp, &b_pip_Ks_PIDp);
   fChain->SetBranchAddress("pip_Ks_ProbNNe", &pip_Ks_ProbNNe, &b_pip_Ks_ProbNNe);
   fChain->SetBranchAddress("pip_Ks_ProbNNk", &pip_Ks_ProbNNk, &b_pip_Ks_ProbNNk);
   fChain->SetBranchAddress("pip_Ks_ProbNNp", &pip_Ks_ProbNNp, &b_pip_Ks_ProbNNp);
   fChain->SetBranchAddress("pip_Ks_ProbNNpi", &pip_Ks_ProbNNpi, &b_pip_Ks_ProbNNpi);
   fChain->SetBranchAddress("pip_Ks_ProbNNmu", &pip_Ks_ProbNNmu, &b_pip_Ks_ProbNNmu);
   fChain->SetBranchAddress("pip_Ks_ProbNNghost", &pip_Ks_ProbNNghost, &b_pip_Ks_ProbNNghost);
   fChain->SetBranchAddress("pip_Ks_hasMuon", &pip_Ks_hasMuon, &b_pip_Ks_hasMuon);
   fChain->SetBranchAddress("pip_Ks_isMuon", &pip_Ks_isMuon, &b_pip_Ks_isMuon);
   fChain->SetBranchAddress("pip_Ks_hasRich", &pip_Ks_hasRich, &b_pip_Ks_hasRich);
   fChain->SetBranchAddress("pip_Ks_UsedRichAerogel", &pip_Ks_UsedRichAerogel, &b_pip_Ks_UsedRichAerogel);
   fChain->SetBranchAddress("pip_Ks_UsedRich1Gas", &pip_Ks_UsedRich1Gas, &b_pip_Ks_UsedRich1Gas);
   fChain->SetBranchAddress("pip_Ks_UsedRich2Gas", &pip_Ks_UsedRich2Gas, &b_pip_Ks_UsedRich2Gas);
   fChain->SetBranchAddress("pip_Ks_RichAboveElThres", &pip_Ks_RichAboveElThres, &b_pip_Ks_RichAboveElThres);
   fChain->SetBranchAddress("pip_Ks_RichAboveMuThres", &pip_Ks_RichAboveMuThres, &b_pip_Ks_RichAboveMuThres);
   fChain->SetBranchAddress("pip_Ks_RichAbovePiThres", &pip_Ks_RichAbovePiThres, &b_pip_Ks_RichAbovePiThres);
   fChain->SetBranchAddress("pip_Ks_RichAboveKaThres", &pip_Ks_RichAboveKaThres, &b_pip_Ks_RichAboveKaThres);
   fChain->SetBranchAddress("pip_Ks_RichAbovePrThres", &pip_Ks_RichAbovePrThres, &b_pip_Ks_RichAbovePrThres);
   fChain->SetBranchAddress("pip_Ks_hasCalo", &pip_Ks_hasCalo, &b_pip_Ks_hasCalo);
   fChain->SetBranchAddress("pip_Ks_TRACK_Type", &pip_Ks_TRACK_Type, &b_pip_Ks_TRACK_Type);
   fChain->SetBranchAddress("pip_Ks_TRACK_Key", &pip_Ks_TRACK_Key, &b_pip_Ks_TRACK_Key);
   fChain->SetBranchAddress("pip_Ks_TRACK_CHI2NDOF", &pip_Ks_TRACK_CHI2NDOF, &b_pip_Ks_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("pip_Ks_TRACK_PCHI2", &pip_Ks_TRACK_PCHI2, &b_pip_Ks_TRACK_PCHI2);
   fChain->SetBranchAddress("pip_Ks_TRACK_MatchCHI2", &pip_Ks_TRACK_MatchCHI2, &b_pip_Ks_TRACK_MatchCHI2);
   fChain->SetBranchAddress("pip_Ks_TRACK_GhostProb", &pip_Ks_TRACK_GhostProb, &b_pip_Ks_TRACK_GhostProb);
   fChain->SetBranchAddress("pip_Ks_TRACK_CloneDist", &pip_Ks_TRACK_CloneDist, &b_pip_Ks_TRACK_CloneDist);
   fChain->SetBranchAddress("pip_Ks_TRACK_Likelihood", &pip_Ks_TRACK_Likelihood, &b_pip_Ks_TRACK_Likelihood);
   fChain->SetBranchAddress("pim_Ks_ETA", &pim_Ks_ETA, &b_pim_Ks_ETA);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV2_ProbNNe", &pim_Ks_MC12TuneV2_ProbNNe, &b_pim_Ks_MC12TuneV2_ProbNNe);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV2_ProbNNmu", &pim_Ks_MC12TuneV2_ProbNNmu, &b_pim_Ks_MC12TuneV2_ProbNNmu);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV2_ProbNNpi", &pim_Ks_MC12TuneV2_ProbNNpi, &b_pim_Ks_MC12TuneV2_ProbNNpi);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV2_ProbNNk", &pim_Ks_MC12TuneV2_ProbNNk, &b_pim_Ks_MC12TuneV2_ProbNNk);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV2_ProbNNp", &pim_Ks_MC12TuneV2_ProbNNp, &b_pim_Ks_MC12TuneV2_ProbNNp);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV2_ProbNNghost", &pim_Ks_MC12TuneV2_ProbNNghost, &b_pim_Ks_MC12TuneV2_ProbNNghost);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV3_ProbNNe", &pim_Ks_MC12TuneV3_ProbNNe, &b_pim_Ks_MC12TuneV3_ProbNNe);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV3_ProbNNmu", &pim_Ks_MC12TuneV3_ProbNNmu, &b_pim_Ks_MC12TuneV3_ProbNNmu);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV3_ProbNNpi", &pim_Ks_MC12TuneV3_ProbNNpi, &b_pim_Ks_MC12TuneV3_ProbNNpi);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV3_ProbNNk", &pim_Ks_MC12TuneV3_ProbNNk, &b_pim_Ks_MC12TuneV3_ProbNNk);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV3_ProbNNp", &pim_Ks_MC12TuneV3_ProbNNp, &b_pim_Ks_MC12TuneV3_ProbNNp);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV3_ProbNNghost", &pim_Ks_MC12TuneV3_ProbNNghost, &b_pim_Ks_MC12TuneV3_ProbNNghost);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV4_ProbNNe", &pim_Ks_MC12TuneV4_ProbNNe, &b_pim_Ks_MC12TuneV4_ProbNNe);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV4_ProbNNmu", &pim_Ks_MC12TuneV4_ProbNNmu, &b_pim_Ks_MC12TuneV4_ProbNNmu);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV4_ProbNNpi", &pim_Ks_MC12TuneV4_ProbNNpi, &b_pim_Ks_MC12TuneV4_ProbNNpi);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV4_ProbNNk", &pim_Ks_MC12TuneV4_ProbNNk, &b_pim_Ks_MC12TuneV4_ProbNNk);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV4_ProbNNp", &pim_Ks_MC12TuneV4_ProbNNp, &b_pim_Ks_MC12TuneV4_ProbNNp);
   fChain->SetBranchAddress("pim_Ks_MC12TuneV4_ProbNNghost", &pim_Ks_MC12TuneV4_ProbNNghost, &b_pim_Ks_MC12TuneV4_ProbNNghost);
   fChain->SetBranchAddress("pim_Ks_MC15TuneV1_ProbNNe", &pim_Ks_MC15TuneV1_ProbNNe, &b_pim_Ks_MC15TuneV1_ProbNNe);
   fChain->SetBranchAddress("pim_Ks_MC15TuneV1_ProbNNmu", &pim_Ks_MC15TuneV1_ProbNNmu, &b_pim_Ks_MC15TuneV1_ProbNNmu);
   fChain->SetBranchAddress("pim_Ks_MC15TuneV1_ProbNNpi", &pim_Ks_MC15TuneV1_ProbNNpi, &b_pim_Ks_MC15TuneV1_ProbNNpi);
   fChain->SetBranchAddress("pim_Ks_MC15TuneV1_ProbNNk", &pim_Ks_MC15TuneV1_ProbNNk, &b_pim_Ks_MC15TuneV1_ProbNNk);
   fChain->SetBranchAddress("pim_Ks_MC15TuneV1_ProbNNp", &pim_Ks_MC15TuneV1_ProbNNp, &b_pim_Ks_MC15TuneV1_ProbNNp);
   fChain->SetBranchAddress("pim_Ks_MC15TuneV1_ProbNNghost", &pim_Ks_MC15TuneV1_ProbNNghost, &b_pim_Ks_MC15TuneV1_ProbNNghost);
   fChain->SetBranchAddress("pim_Ks_MINIP", &pim_Ks_MINIP, &b_pim_Ks_MINIP);
   fChain->SetBranchAddress("pim_Ks_MINIPCHI2", &pim_Ks_MINIPCHI2, &b_pim_Ks_MINIPCHI2);
   fChain->SetBranchAddress("pim_Ks_MINIPNEXTBEST", &pim_Ks_MINIPNEXTBEST, &b_pim_Ks_MINIPNEXTBEST);
   fChain->SetBranchAddress("pim_Ks_MINIPCHI2NEXTBEST", &pim_Ks_MINIPCHI2NEXTBEST, &b_pim_Ks_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("pim_Ks_OWNPV_X", &pim_Ks_OWNPV_X, &b_pim_Ks_OWNPV_X);
   fChain->SetBranchAddress("pim_Ks_OWNPV_Y", &pim_Ks_OWNPV_Y, &b_pim_Ks_OWNPV_Y);
   fChain->SetBranchAddress("pim_Ks_OWNPV_Z", &pim_Ks_OWNPV_Z, &b_pim_Ks_OWNPV_Z);
   fChain->SetBranchAddress("pim_Ks_OWNPV_XERR", &pim_Ks_OWNPV_XERR, &b_pim_Ks_OWNPV_XERR);
   fChain->SetBranchAddress("pim_Ks_OWNPV_YERR", &pim_Ks_OWNPV_YERR, &b_pim_Ks_OWNPV_YERR);
   fChain->SetBranchAddress("pim_Ks_OWNPV_ZERR", &pim_Ks_OWNPV_ZERR, &b_pim_Ks_OWNPV_ZERR);
   fChain->SetBranchAddress("pim_Ks_OWNPV_CHI2", &pim_Ks_OWNPV_CHI2, &b_pim_Ks_OWNPV_CHI2);
   fChain->SetBranchAddress("pim_Ks_OWNPV_NDOF", &pim_Ks_OWNPV_NDOF, &b_pim_Ks_OWNPV_NDOF);
   fChain->SetBranchAddress("pim_Ks_OWNPV_COV_", pim_Ks_OWNPV_COV_, &b_pim_Ks_OWNPV_COV_);
   fChain->SetBranchAddress("pim_Ks_IP_OWNPV", &pim_Ks_IP_OWNPV, &b_pim_Ks_IP_OWNPV);
   fChain->SetBranchAddress("pim_Ks_IPCHI2_OWNPV", &pim_Ks_IPCHI2_OWNPV, &b_pim_Ks_IPCHI2_OWNPV);
   fChain->SetBranchAddress("pim_Ks_TOPPV_X", &pim_Ks_TOPPV_X, &b_pim_Ks_TOPPV_X);
   fChain->SetBranchAddress("pim_Ks_TOPPV_Y", &pim_Ks_TOPPV_Y, &b_pim_Ks_TOPPV_Y);
   fChain->SetBranchAddress("pim_Ks_TOPPV_Z", &pim_Ks_TOPPV_Z, &b_pim_Ks_TOPPV_Z);
   fChain->SetBranchAddress("pim_Ks_TOPPV_XERR", &pim_Ks_TOPPV_XERR, &b_pim_Ks_TOPPV_XERR);
   fChain->SetBranchAddress("pim_Ks_TOPPV_YERR", &pim_Ks_TOPPV_YERR, &b_pim_Ks_TOPPV_YERR);
   fChain->SetBranchAddress("pim_Ks_TOPPV_ZERR", &pim_Ks_TOPPV_ZERR, &b_pim_Ks_TOPPV_ZERR);
   fChain->SetBranchAddress("pim_Ks_TOPPV_CHI2", &pim_Ks_TOPPV_CHI2, &b_pim_Ks_TOPPV_CHI2);
   fChain->SetBranchAddress("pim_Ks_TOPPV_NDOF", &pim_Ks_TOPPV_NDOF, &b_pim_Ks_TOPPV_NDOF);
   fChain->SetBranchAddress("pim_Ks_TOPPV_COV_", pim_Ks_TOPPV_COV_, &b_pim_Ks_TOPPV_COV_);
   fChain->SetBranchAddress("pim_Ks_IP_TOPPV", &pim_Ks_IP_TOPPV, &b_pim_Ks_IP_TOPPV);
   fChain->SetBranchAddress("pim_Ks_IPCHI2_TOPPV", &pim_Ks_IPCHI2_TOPPV, &b_pim_Ks_IPCHI2_TOPPV);
   fChain->SetBranchAddress("pim_Ks_ORIVX_X", &pim_Ks_ORIVX_X, &b_pim_Ks_ORIVX_X);
   fChain->SetBranchAddress("pim_Ks_ORIVX_Y", &pim_Ks_ORIVX_Y, &b_pim_Ks_ORIVX_Y);
   fChain->SetBranchAddress("pim_Ks_ORIVX_Z", &pim_Ks_ORIVX_Z, &b_pim_Ks_ORIVX_Z);
   fChain->SetBranchAddress("pim_Ks_ORIVX_XERR", &pim_Ks_ORIVX_XERR, &b_pim_Ks_ORIVX_XERR);
   fChain->SetBranchAddress("pim_Ks_ORIVX_YERR", &pim_Ks_ORIVX_YERR, &b_pim_Ks_ORIVX_YERR);
   fChain->SetBranchAddress("pim_Ks_ORIVX_ZERR", &pim_Ks_ORIVX_ZERR, &b_pim_Ks_ORIVX_ZERR);
   fChain->SetBranchAddress("pim_Ks_ORIVX_CHI2", &pim_Ks_ORIVX_CHI2, &b_pim_Ks_ORIVX_CHI2);
   fChain->SetBranchAddress("pim_Ks_ORIVX_NDOF", &pim_Ks_ORIVX_NDOF, &b_pim_Ks_ORIVX_NDOF);
   fChain->SetBranchAddress("pim_Ks_ORIVX_COV_", pim_Ks_ORIVX_COV_, &b_pim_Ks_ORIVX_COV_);
   fChain->SetBranchAddress("pim_Ks_IP_ORIVX", &pim_Ks_IP_ORIVX, &b_pim_Ks_IP_ORIVX);
   fChain->SetBranchAddress("pim_Ks_IPCHI2_ORIVX", &pim_Ks_IPCHI2_ORIVX, &b_pim_Ks_IPCHI2_ORIVX);
   fChain->SetBranchAddress("pim_Ks_P", &pim_Ks_P, &b_pim_Ks_P);
   fChain->SetBranchAddress("pim_Ks_PT", &pim_Ks_PT, &b_pim_Ks_PT);
   fChain->SetBranchAddress("pim_Ks_PE", &pim_Ks_PE, &b_pim_Ks_PE);
   fChain->SetBranchAddress("pim_Ks_PX", &pim_Ks_PX, &b_pim_Ks_PX);
   fChain->SetBranchAddress("pim_Ks_PY", &pim_Ks_PY, &b_pim_Ks_PY);
   fChain->SetBranchAddress("pim_Ks_PZ", &pim_Ks_PZ, &b_pim_Ks_PZ);
   fChain->SetBranchAddress("pim_Ks_M", &pim_Ks_M, &b_pim_Ks_M);
   fChain->SetBranchAddress("pim_Ks_ID", &pim_Ks_ID, &b_pim_Ks_ID);
   fChain->SetBranchAddress("pim_Ks_PIDe", &pim_Ks_PIDe, &b_pim_Ks_PIDe);
   fChain->SetBranchAddress("pim_Ks_PIDmu", &pim_Ks_PIDmu, &b_pim_Ks_PIDmu);
   fChain->SetBranchAddress("pim_Ks_PIDK", &pim_Ks_PIDK, &b_pim_Ks_PIDK);
   fChain->SetBranchAddress("pim_Ks_PIDp", &pim_Ks_PIDp, &b_pim_Ks_PIDp);
   fChain->SetBranchAddress("pim_Ks_ProbNNe", &pim_Ks_ProbNNe, &b_pim_Ks_ProbNNe);
   fChain->SetBranchAddress("pim_Ks_ProbNNk", &pim_Ks_ProbNNk, &b_pim_Ks_ProbNNk);
   fChain->SetBranchAddress("pim_Ks_ProbNNp", &pim_Ks_ProbNNp, &b_pim_Ks_ProbNNp);
   fChain->SetBranchAddress("pim_Ks_ProbNNpi", &pim_Ks_ProbNNpi, &b_pim_Ks_ProbNNpi);
   fChain->SetBranchAddress("pim_Ks_ProbNNmu", &pim_Ks_ProbNNmu, &b_pim_Ks_ProbNNmu);
   fChain->SetBranchAddress("pim_Ks_ProbNNghost", &pim_Ks_ProbNNghost, &b_pim_Ks_ProbNNghost);
   fChain->SetBranchAddress("pim_Ks_hasMuon", &pim_Ks_hasMuon, &b_pim_Ks_hasMuon);
   fChain->SetBranchAddress("pim_Ks_isMuon", &pim_Ks_isMuon, &b_pim_Ks_isMuon);
   fChain->SetBranchAddress("pim_Ks_hasRich", &pim_Ks_hasRich, &b_pim_Ks_hasRich);
   fChain->SetBranchAddress("pim_Ks_UsedRichAerogel", &pim_Ks_UsedRichAerogel, &b_pim_Ks_UsedRichAerogel);
   fChain->SetBranchAddress("pim_Ks_UsedRich1Gas", &pim_Ks_UsedRich1Gas, &b_pim_Ks_UsedRich1Gas);
   fChain->SetBranchAddress("pim_Ks_UsedRich2Gas", &pim_Ks_UsedRich2Gas, &b_pim_Ks_UsedRich2Gas);
   fChain->SetBranchAddress("pim_Ks_RichAboveElThres", &pim_Ks_RichAboveElThres, &b_pim_Ks_RichAboveElThres);
   fChain->SetBranchAddress("pim_Ks_RichAboveMuThres", &pim_Ks_RichAboveMuThres, &b_pim_Ks_RichAboveMuThres);
   fChain->SetBranchAddress("pim_Ks_RichAbovePiThres", &pim_Ks_RichAbovePiThres, &b_pim_Ks_RichAbovePiThres);
   fChain->SetBranchAddress("pim_Ks_RichAboveKaThres", &pim_Ks_RichAboveKaThres, &b_pim_Ks_RichAboveKaThres);
   fChain->SetBranchAddress("pim_Ks_RichAbovePrThres", &pim_Ks_RichAbovePrThres, &b_pim_Ks_RichAbovePrThres);
   fChain->SetBranchAddress("pim_Ks_hasCalo", &pim_Ks_hasCalo, &b_pim_Ks_hasCalo);
   fChain->SetBranchAddress("pim_Ks_TRACK_Type", &pim_Ks_TRACK_Type, &b_pim_Ks_TRACK_Type);
   fChain->SetBranchAddress("pim_Ks_TRACK_Key", &pim_Ks_TRACK_Key, &b_pim_Ks_TRACK_Key);
   fChain->SetBranchAddress("pim_Ks_TRACK_CHI2NDOF", &pim_Ks_TRACK_CHI2NDOF, &b_pim_Ks_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("pim_Ks_TRACK_PCHI2", &pim_Ks_TRACK_PCHI2, &b_pim_Ks_TRACK_PCHI2);
   fChain->SetBranchAddress("pim_Ks_TRACK_MatchCHI2", &pim_Ks_TRACK_MatchCHI2, &b_pim_Ks_TRACK_MatchCHI2);
   fChain->SetBranchAddress("pim_Ks_TRACK_GhostProb", &pim_Ks_TRACK_GhostProb, &b_pim_Ks_TRACK_GhostProb);
   fChain->SetBranchAddress("pim_Ks_TRACK_CloneDist", &pim_Ks_TRACK_CloneDist, &b_pim_Ks_TRACK_CloneDist);
   fChain->SetBranchAddress("pim_Ks_TRACK_Likelihood", &pim_Ks_TRACK_Likelihood, &b_pim_Ks_TRACK_Likelihood);
   fChain->SetBranchAddress("pi_ETA", &pi_ETA, &b_pi_ETA);
   fChain->SetBranchAddress("pi_MC12TuneV2_ProbNNe", &pi_MC12TuneV2_ProbNNe, &b_pi_MC12TuneV2_ProbNNe);
   fChain->SetBranchAddress("pi_MC12TuneV2_ProbNNmu", &pi_MC12TuneV2_ProbNNmu, &b_pi_MC12TuneV2_ProbNNmu);
   fChain->SetBranchAddress("pi_MC12TuneV2_ProbNNpi", &pi_MC12TuneV2_ProbNNpi, &b_pi_MC12TuneV2_ProbNNpi);
   fChain->SetBranchAddress("pi_MC12TuneV2_ProbNNk", &pi_MC12TuneV2_ProbNNk, &b_pi_MC12TuneV2_ProbNNk);
   fChain->SetBranchAddress("pi_MC12TuneV2_ProbNNp", &pi_MC12TuneV2_ProbNNp, &b_pi_MC12TuneV2_ProbNNp);
   fChain->SetBranchAddress("pi_MC12TuneV2_ProbNNghost", &pi_MC12TuneV2_ProbNNghost, &b_pi_MC12TuneV2_ProbNNghost);
   fChain->SetBranchAddress("pi_MC12TuneV3_ProbNNe", &pi_MC12TuneV3_ProbNNe, &b_pi_MC12TuneV3_ProbNNe);
   fChain->SetBranchAddress("pi_MC12TuneV3_ProbNNmu", &pi_MC12TuneV3_ProbNNmu, &b_pi_MC12TuneV3_ProbNNmu);
   fChain->SetBranchAddress("pi_MC12TuneV3_ProbNNpi", &pi_MC12TuneV3_ProbNNpi, &b_pi_MC12TuneV3_ProbNNpi);
   fChain->SetBranchAddress("pi_MC12TuneV3_ProbNNk", &pi_MC12TuneV3_ProbNNk, &b_pi_MC12TuneV3_ProbNNk);
   fChain->SetBranchAddress("pi_MC12TuneV3_ProbNNp", &pi_MC12TuneV3_ProbNNp, &b_pi_MC12TuneV3_ProbNNp);
   fChain->SetBranchAddress("pi_MC12TuneV3_ProbNNghost", &pi_MC12TuneV3_ProbNNghost, &b_pi_MC12TuneV3_ProbNNghost);
   fChain->SetBranchAddress("pi_MC12TuneV4_ProbNNe", &pi_MC12TuneV4_ProbNNe, &b_pi_MC12TuneV4_ProbNNe);
   fChain->SetBranchAddress("pi_MC12TuneV4_ProbNNmu", &pi_MC12TuneV4_ProbNNmu, &b_pi_MC12TuneV4_ProbNNmu);
   fChain->SetBranchAddress("pi_MC12TuneV4_ProbNNpi", &pi_MC12TuneV4_ProbNNpi, &b_pi_MC12TuneV4_ProbNNpi);
   fChain->SetBranchAddress("pi_MC12TuneV4_ProbNNk", &pi_MC12TuneV4_ProbNNk, &b_pi_MC12TuneV4_ProbNNk);
   fChain->SetBranchAddress("pi_MC12TuneV4_ProbNNp", &pi_MC12TuneV4_ProbNNp, &b_pi_MC12TuneV4_ProbNNp);
   fChain->SetBranchAddress("pi_MC12TuneV4_ProbNNghost", &pi_MC12TuneV4_ProbNNghost, &b_pi_MC12TuneV4_ProbNNghost);
   fChain->SetBranchAddress("pi_MC15TuneV1_ProbNNe", &pi_MC15TuneV1_ProbNNe, &b_pi_MC15TuneV1_ProbNNe);
   fChain->SetBranchAddress("pi_MC15TuneV1_ProbNNmu", &pi_MC15TuneV1_ProbNNmu, &b_pi_MC15TuneV1_ProbNNmu);
   fChain->SetBranchAddress("pi_MC15TuneV1_ProbNNpi", &pi_MC15TuneV1_ProbNNpi, &b_pi_MC15TuneV1_ProbNNpi);
   fChain->SetBranchAddress("pi_MC15TuneV1_ProbNNk", &pi_MC15TuneV1_ProbNNk, &b_pi_MC15TuneV1_ProbNNk);
   fChain->SetBranchAddress("pi_MC15TuneV1_ProbNNp", &pi_MC15TuneV1_ProbNNp, &b_pi_MC15TuneV1_ProbNNp);
   fChain->SetBranchAddress("pi_MC15TuneV1_ProbNNghost", &pi_MC15TuneV1_ProbNNghost, &b_pi_MC15TuneV1_ProbNNghost);
   fChain->SetBranchAddress("pi_MINIP", &pi_MINIP, &b_pi_MINIP);
   fChain->SetBranchAddress("pi_MINIPCHI2", &pi_MINIPCHI2, &b_pi_MINIPCHI2);
   fChain->SetBranchAddress("pi_MINIPNEXTBEST", &pi_MINIPNEXTBEST, &b_pi_MINIPNEXTBEST);
   fChain->SetBranchAddress("pi_MINIPCHI2NEXTBEST", &pi_MINIPCHI2NEXTBEST, &b_pi_MINIPCHI2NEXTBEST);
   fChain->SetBranchAddress("pi_OWNPV_X", &pi_OWNPV_X, &b_pi_OWNPV_X);
   fChain->SetBranchAddress("pi_OWNPV_Y", &pi_OWNPV_Y, &b_pi_OWNPV_Y);
   fChain->SetBranchAddress("pi_OWNPV_Z", &pi_OWNPV_Z, &b_pi_OWNPV_Z);
   fChain->SetBranchAddress("pi_OWNPV_XERR", &pi_OWNPV_XERR, &b_pi_OWNPV_XERR);
   fChain->SetBranchAddress("pi_OWNPV_YERR", &pi_OWNPV_YERR, &b_pi_OWNPV_YERR);
   fChain->SetBranchAddress("pi_OWNPV_ZERR", &pi_OWNPV_ZERR, &b_pi_OWNPV_ZERR);
   fChain->SetBranchAddress("pi_OWNPV_CHI2", &pi_OWNPV_CHI2, &b_pi_OWNPV_CHI2);
   fChain->SetBranchAddress("pi_OWNPV_NDOF", &pi_OWNPV_NDOF, &b_pi_OWNPV_NDOF);
   fChain->SetBranchAddress("pi_OWNPV_COV_", pi_OWNPV_COV_, &b_pi_OWNPV_COV_);
   fChain->SetBranchAddress("pi_IP_OWNPV", &pi_IP_OWNPV, &b_pi_IP_OWNPV);
   fChain->SetBranchAddress("pi_IPCHI2_OWNPV", &pi_IPCHI2_OWNPV, &b_pi_IPCHI2_OWNPV);
   fChain->SetBranchAddress("pi_TOPPV_X", &pi_TOPPV_X, &b_pi_TOPPV_X);
   fChain->SetBranchAddress("pi_TOPPV_Y", &pi_TOPPV_Y, &b_pi_TOPPV_Y);
   fChain->SetBranchAddress("pi_TOPPV_Z", &pi_TOPPV_Z, &b_pi_TOPPV_Z);
   fChain->SetBranchAddress("pi_TOPPV_XERR", &pi_TOPPV_XERR, &b_pi_TOPPV_XERR);
   fChain->SetBranchAddress("pi_TOPPV_YERR", &pi_TOPPV_YERR, &b_pi_TOPPV_YERR);
   fChain->SetBranchAddress("pi_TOPPV_ZERR", &pi_TOPPV_ZERR, &b_pi_TOPPV_ZERR);
   fChain->SetBranchAddress("pi_TOPPV_CHI2", &pi_TOPPV_CHI2, &b_pi_TOPPV_CHI2);
   fChain->SetBranchAddress("pi_TOPPV_NDOF", &pi_TOPPV_NDOF, &b_pi_TOPPV_NDOF);
   fChain->SetBranchAddress("pi_TOPPV_COV_", pi_TOPPV_COV_, &b_pi_TOPPV_COV_);
   fChain->SetBranchAddress("pi_IP_TOPPV", &pi_IP_TOPPV, &b_pi_IP_TOPPV);
   fChain->SetBranchAddress("pi_IPCHI2_TOPPV", &pi_IPCHI2_TOPPV, &b_pi_IPCHI2_TOPPV);
   fChain->SetBranchAddress("pi_ORIVX_X", &pi_ORIVX_X, &b_pi_ORIVX_X);
   fChain->SetBranchAddress("pi_ORIVX_Y", &pi_ORIVX_Y, &b_pi_ORIVX_Y);
   fChain->SetBranchAddress("pi_ORIVX_Z", &pi_ORIVX_Z, &b_pi_ORIVX_Z);
   fChain->SetBranchAddress("pi_ORIVX_XERR", &pi_ORIVX_XERR, &b_pi_ORIVX_XERR);
   fChain->SetBranchAddress("pi_ORIVX_YERR", &pi_ORIVX_YERR, &b_pi_ORIVX_YERR);
   fChain->SetBranchAddress("pi_ORIVX_ZERR", &pi_ORIVX_ZERR, &b_pi_ORIVX_ZERR);
   fChain->SetBranchAddress("pi_ORIVX_CHI2", &pi_ORIVX_CHI2, &b_pi_ORIVX_CHI2);
   fChain->SetBranchAddress("pi_ORIVX_NDOF", &pi_ORIVX_NDOF, &b_pi_ORIVX_NDOF);
   fChain->SetBranchAddress("pi_ORIVX_COV_", pi_ORIVX_COV_, &b_pi_ORIVX_COV_);
   fChain->SetBranchAddress("pi_IP_ORIVX", &pi_IP_ORIVX, &b_pi_IP_ORIVX);
   fChain->SetBranchAddress("pi_IPCHI2_ORIVX", &pi_IPCHI2_ORIVX, &b_pi_IPCHI2_ORIVX);
   fChain->SetBranchAddress("pi_P", &pi_P, &b_pi_P);
   fChain->SetBranchAddress("pi_PT", &pi_PT, &b_pi_PT);
   fChain->SetBranchAddress("pi_PE", &pi_PE, &b_pi_PE);
   fChain->SetBranchAddress("pi_PX", &pi_PX, &b_pi_PX);
   fChain->SetBranchAddress("pi_PY", &pi_PY, &b_pi_PY);
   fChain->SetBranchAddress("pi_PZ", &pi_PZ, &b_pi_PZ);
   fChain->SetBranchAddress("pi_M", &pi_M, &b_pi_M);
   fChain->SetBranchAddress("pi_ID", &pi_ID, &b_pi_ID);
   fChain->SetBranchAddress("pi_PIDe", &pi_PIDe, &b_pi_PIDe);
   fChain->SetBranchAddress("pi_PIDmu", &pi_PIDmu, &b_pi_PIDmu);
   fChain->SetBranchAddress("pi_PIDK", &pi_PIDK, &b_pi_PIDK);
   fChain->SetBranchAddress("pi_PIDp", &pi_PIDp, &b_pi_PIDp);
   fChain->SetBranchAddress("pi_ProbNNe", &pi_ProbNNe, &b_pi_ProbNNe);
   fChain->SetBranchAddress("pi_ProbNNk", &pi_ProbNNk, &b_pi_ProbNNk);
   fChain->SetBranchAddress("pi_ProbNNp", &pi_ProbNNp, &b_pi_ProbNNp);
   fChain->SetBranchAddress("pi_ProbNNpi", &pi_ProbNNpi, &b_pi_ProbNNpi);
   fChain->SetBranchAddress("pi_ProbNNmu", &pi_ProbNNmu, &b_pi_ProbNNmu);
   fChain->SetBranchAddress("pi_ProbNNghost", &pi_ProbNNghost, &b_pi_ProbNNghost);
   fChain->SetBranchAddress("pi_hasMuon", &pi_hasMuon, &b_pi_hasMuon);
   fChain->SetBranchAddress("pi_isMuon", &pi_isMuon, &b_pi_isMuon);
   fChain->SetBranchAddress("pi_hasRich", &pi_hasRich, &b_pi_hasRich);
   fChain->SetBranchAddress("pi_UsedRichAerogel", &pi_UsedRichAerogel, &b_pi_UsedRichAerogel);
   fChain->SetBranchAddress("pi_UsedRich1Gas", &pi_UsedRich1Gas, &b_pi_UsedRich1Gas);
   fChain->SetBranchAddress("pi_UsedRich2Gas", &pi_UsedRich2Gas, &b_pi_UsedRich2Gas);
   fChain->SetBranchAddress("pi_RichAboveElThres", &pi_RichAboveElThres, &b_pi_RichAboveElThres);
   fChain->SetBranchAddress("pi_RichAboveMuThres", &pi_RichAboveMuThres, &b_pi_RichAboveMuThres);
   fChain->SetBranchAddress("pi_RichAbovePiThres", &pi_RichAbovePiThres, &b_pi_RichAbovePiThres);
   fChain->SetBranchAddress("pi_RichAboveKaThres", &pi_RichAboveKaThres, &b_pi_RichAboveKaThres);
   fChain->SetBranchAddress("pi_RichAbovePrThres", &pi_RichAbovePrThres, &b_pi_RichAbovePrThres);
   fChain->SetBranchAddress("pi_hasCalo", &pi_hasCalo, &b_pi_hasCalo);
   fChain->SetBranchAddress("pi_TRACK_Type", &pi_TRACK_Type, &b_pi_TRACK_Type);
   fChain->SetBranchAddress("pi_TRACK_Key", &pi_TRACK_Key, &b_pi_TRACK_Key);
   fChain->SetBranchAddress("pi_TRACK_CHI2NDOF", &pi_TRACK_CHI2NDOF, &b_pi_TRACK_CHI2NDOF);
   fChain->SetBranchAddress("pi_TRACK_PCHI2", &pi_TRACK_PCHI2, &b_pi_TRACK_PCHI2);
   fChain->SetBranchAddress("pi_TRACK_MatchCHI2", &pi_TRACK_MatchCHI2, &b_pi_TRACK_MatchCHI2);
   fChain->SetBranchAddress("pi_TRACK_GhostProb", &pi_TRACK_GhostProb, &b_pi_TRACK_GhostProb);
   fChain->SetBranchAddress("pi_TRACK_CloneDist", &pi_TRACK_CloneDist, &b_pi_TRACK_CloneDist);
   fChain->SetBranchAddress("pi_TRACK_Likelihood", &pi_TRACK_Likelihood, &b_pi_TRACK_Likelihood);
   fChain->SetBranchAddress("nCandidate", &nCandidate, &b_nCandidate);
   fChain->SetBranchAddress("totCandidates", &totCandidates, &b_totCandidates);
   fChain->SetBranchAddress("EventInSequence", &EventInSequence, &b_EventInSequence);
   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("eventNumber", &eventNumber, &b_eventNumber);
   fChain->SetBranchAddress("BCID", &BCID, &b_BCID);
   fChain->SetBranchAddress("BCType", &BCType, &b_BCType);
   fChain->SetBranchAddress("OdinTCK", &OdinTCK, &b_OdinTCK);
   fChain->SetBranchAddress("L0DUTCK", &L0DUTCK, &b_L0DUTCK);
   fChain->SetBranchAddress("HLT1TCK", &HLT1TCK, &b_HLT1TCK);
   fChain->SetBranchAddress("HLT2TCK", &HLT2TCK, &b_HLT2TCK);
   fChain->SetBranchAddress("GpsTime", &GpsTime, &b_GpsTime);
   fChain->SetBranchAddress("Polarity", &Polarity, &b_Polarity);
   fChain->SetBranchAddress("PVX", PVX, &b_PVX);
   fChain->SetBranchAddress("PVY", PVY, &b_PVY);
   fChain->SetBranchAddress("PVZ", PVZ, &b_PVZ);
   fChain->SetBranchAddress("PVXERR", PVXERR, &b_PVXERR);
   fChain->SetBranchAddress("PVYERR", PVYERR, &b_PVYERR);
   fChain->SetBranchAddress("PVZERR", PVZERR, &b_PVZERR);
   fChain->SetBranchAddress("PVCHI2", PVCHI2, &b_PVCHI2);
   fChain->SetBranchAddress("PVNDOF", PVNDOF, &b_PVNDOF);
   fChain->SetBranchAddress("PVNTRACKS", PVNTRACKS, &b_PVNTRACKS);
   fChain->SetBranchAddress("nPVs", &nPVs, &b_nPVs);
   fChain->SetBranchAddress("nTracks", &nTracks, &b_nTracks);
   fChain->SetBranchAddress("nLongTracks", &nLongTracks, &b_nLongTracks);
   fChain->SetBranchAddress("nDownstreamTracks", &nDownstreamTracks, &b_nDownstreamTracks);
   fChain->SetBranchAddress("nUpstreamTracks", &nUpstreamTracks, &b_nUpstreamTracks);
   fChain->SetBranchAddress("nVeloTracks", &nVeloTracks, &b_nVeloTracks);
   fChain->SetBranchAddress("nTTracks", &nTTracks, &b_nTTracks);
   fChain->SetBranchAddress("nBackTracks", &nBackTracks, &b_nBackTracks);
   fChain->SetBranchAddress("nRich1Hits", &nRich1Hits, &b_nRich1Hits);
   fChain->SetBranchAddress("nRich2Hits", &nRich2Hits, &b_nRich2Hits);
   fChain->SetBranchAddress("nVeloClusters", &nVeloClusters, &b_nVeloClusters);
   fChain->SetBranchAddress("nITClusters", &nITClusters, &b_nITClusters);
   fChain->SetBranchAddress("nTTClusters", &nTTClusters, &b_nTTClusters);
   fChain->SetBranchAddress("nOTClusters", &nOTClusters, &b_nOTClusters);
   fChain->SetBranchAddress("nSPDHits", &nSPDHits, &b_nSPDHits);
   fChain->SetBranchAddress("nMuonCoordsS0", &nMuonCoordsS0, &b_nMuonCoordsS0);
   fChain->SetBranchAddress("nMuonCoordsS1", &nMuonCoordsS1, &b_nMuonCoordsS1);
   fChain->SetBranchAddress("nMuonCoordsS2", &nMuonCoordsS2, &b_nMuonCoordsS2);
   fChain->SetBranchAddress("nMuonCoordsS3", &nMuonCoordsS3, &b_nMuonCoordsS3);
   fChain->SetBranchAddress("nMuonCoordsS4", &nMuonCoordsS4, &b_nMuonCoordsS4);
   fChain->SetBranchAddress("nMuonTracks", &nMuonTracks, &b_nMuonTracks);
   Notify();
}

Bool_t DecayTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void DecayTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}

#endif // #ifdef DecayTree_cxx
