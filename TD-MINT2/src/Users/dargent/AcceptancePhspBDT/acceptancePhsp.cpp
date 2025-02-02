// Phaspace acceptance studies
// author: Philippe d'Argent, Matthieu Kecke
#include "Mint/DalitzEvent.h"
#include "Mint/DalitzEventList.h"
#include "Mint/DiskResidentEventList.h"
#include "TTree.h"
#include <iostream>
#include "TVector3.h"
#include "TLorentzVector.h"
#include "Mint/DalitzEventPattern.h"
#include "Mint/CLHEPSystemOfUnits.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "Mint/PlotSet.h"
#include "TChain.h"
#include "Mint/NamedParameter.h"
#include <TROOT.h>
#include <sstream>
#include "TProfile.h"
#include "Mint/NamedParameter.h"
#include "Mint/HyperHistogram.h"
#include "Mint/Utils.h"

using namespace std;
using namespace MINT;

double dotProduct( const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& pX )
{
    return -p1.Dot( p2 ) + p1.Dot( pX ) * p2.Dot( pX ) / pX.Dot( pX );
}

double helicityCosine( const TLorentzVector& p1, const TLorentzVector& p2, const TLorentzVector& pX)
{
    return dotProduct(p1, p2, pX) / sqrt( dotProduct( p1, p1, pX ) * dotProduct( p2, p2, pX ) );
}

double cosThetaAngle(TLorentzVector p0,TLorentzVector p1,TLorentzVector p2,TLorentzVector p3){
 	TLorentzVector pD = p0 + p1 + p2 + p3 ;
 	p0.Boost( - pD.BoostVector() );
 	p1.Boost( - pD.BoostVector() );
 	p2.Boost( - pD.BoostVector() );
 	p3.Boost( - pD.BoostVector() );

	TVector3 mother = (p0+p1).Vect().Unit();
	p0.Boost( - (p0+p1).BoostVector());
	TVector3 daughter = p0.Vect().Unit();
	
	return mother.Dot(daughter);
}

double acoplanarityAngle(TLorentzVector p0,TLorentzVector p1,TLorentzVector p2,TLorentzVector p3){
 	TLorentzVector pD = p0 + p1 + p2 + p3 ;
 	p0.Boost( - pD.BoostVector() );
 	p1.Boost( - pD.BoostVector() );
 	p2.Boost( - pD.BoostVector() );
 	p3.Boost( - pD.BoostVector() );
 	TVector3 e1 = (p0.Vect().Cross( p1.Vect() )).Unit();
 	TVector3 e2 = (p2.Vect().Cross( p3.Vect() )).Unit();
 	//return t1.Angle( t2 ); 	
	TVector3 ez=  (p3+p2).Vect().Unit();

        double cosPhi= e1.Dot(e2);
	double sinPhi = (e1.Cross(e2)).Dot(ez);
	double phi= acos(cosPhi);
	return (sinPhi > 0.0 ? phi : -phi);
}

double cosThetaAngle(const DalitzEvent& evt, int a, int b, int c, int d){
	TLorentzVector p0 = evt.p(a);
  	TLorentzVector p1 = evt.p(b) ;
  	TLorentzVector p2 = evt.p(c) ;
 	TLorentzVector p3 = evt.p(d) ;
 	TLorentzVector pD = p0 + p1 + p2 + p3 ;
 	p0.Boost( - pD.BoostVector() );
 	p1.Boost( - pD.BoostVector() );
 	p2.Boost( - pD.BoostVector() );
 	p3.Boost( - pD.BoostVector() );

	TVector3 mother = (p0+p1).Vect().Unit();
	p0.Boost( - (p0+p1).BoostVector());
	TVector3 daughter = p0.Vect().Unit();
	
	return mother.Dot(daughter);
}

double acoplanarityAngle(const DalitzEvent& evt, int a, int b, int c, int d){
	TLorentzVector p0 = evt.p(a);
  	TLorentzVector p1 = evt.p(b) ;
  	TLorentzVector p2 = evt.p(c) ;
 	TLorentzVector p3 = evt.p(d) ;
 	TLorentzVector pD = p0 + p1 + p2 + p3 ;
 	p0.Boost( - pD.BoostVector() );
 	p1.Boost( - pD.BoostVector() );
 	p2.Boost( - pD.BoostVector() );
 	p3.Boost( - pD.BoostVector() );
 	TVector3 e1 = (p0.Vect().Cross( p1.Vect() )).Unit();
 	TVector3 e2 = (p2.Vect().Cross( p3.Vect() )).Unit();
 	//return t1.Angle( t2 ); 	
	TVector3 ez=  (p3+p2).Vect().Unit();

        double cosPhi= e1.Dot(e2);
	double sinPhi = (e1.Cross(e2)).Dot(ez);
	double phi= acos(cosPhi);
	return (sinPhi > 0.0 ? phi : -phi);
}

void prepareFilesForBDT(){

    DalitzEventPattern pdg(511 ,-411 ,310 ,211);
    NamedParameter<string> InputGenMC("InputGenMC", (std::string) "Gen_11166161.root");
    TChain* tree_gen=new TChain("MCDecayTreeTuple/MCDecayTree");
    tree_gen->Add(((string)InputGenMC).c_str());
    double Ks_gen[5]; 
    double pi_gen[5]; 
    double D_gen[5]; 
    
    tree_gen->SetBranchAddress("KS0_TRUEP_X",&Ks_gen[0]);
    tree_gen->SetBranchAddress("KS0_TRUEP_Y",&Ks_gen[1]);
    tree_gen->SetBranchAddress("KS0_TRUEP_Z",&Ks_gen[2]); 
    tree_gen->SetBranchAddress("KS0_TRUEP_E",&Ks_gen[3]); 
	
    tree_gen->SetBranchAddress("piplus_TRUEP_X",&pi_gen[0]);
    tree_gen->SetBranchAddress("piplus_TRUEP_Y",&pi_gen[1]);
    tree_gen->SetBranchAddress("piplus_TRUEP_Z",&pi_gen[2]); 
    tree_gen->SetBranchAddress("piplus_TRUEP_E",&pi_gen[3]); 
    
    tree_gen->SetBranchAddress("Dminus_TRUEP_X",&D_gen[0]);
    tree_gen->SetBranchAddress("Dminus_TRUEP_Y",&D_gen[1]);
    tree_gen->SetBranchAddress("Dminus_TRUEP_Z",&D_gen[2]); 
    tree_gen->SetBranchAddress("Dminus_TRUEP_E",&D_gen[3]); 

    TFile* output_gen = new TFile("GenMC.root","RECREATE");
    TTree* summary_tree_gen = tree_gen->CloneTree(0);
    double m_DKs,m_Dpi,m_Kspi;
    summary_tree_gen->Branch("TRUE_m_DKs", &m_DKs, "TRUE_m_DKs/D");
    summary_tree_gen->Branch("TRUE_m_Dpi", &m_Dpi, "TRUE_m_Dpi/D");
    summary_tree_gen->Branch("TRUE_m_Kspi", &m_Kspi, "TRUE_m_Kspi/D");

    double mp_DKs,mp_Dpi,mp_Kspi;
    summary_tree_gen->Branch("TRUE_mp_DKs", &mp_DKs, "TRUE_mp_DKs/D");
    summary_tree_gen->Branch("TRUE_mp_Dpi", &mp_Dpi, "TRUE_mp_Dpi/D");
    summary_tree_gen->Branch("TRUE_mp_Kspi", &mp_Kspi, "TRUE_mp_Kspi/D");

    double cos_DKs,cos_Dpi,cos_Kspi;
    summary_tree_gen->Branch("TRUE_cos_DKs", &cos_DKs, "TRUE_cos_DKs/D");
    summary_tree_gen->Branch("TRUE_cos_Dpi", &cos_Dpi, "TRUE_cos_Dpi/D");
    summary_tree_gen->Branch("TRUE_cos_Kspi", &cos_Kspi, "TRUE_cos_Kspi/D");

    int KsCat = -1;
    summary_tree_gen->Branch("KsCat", &KsCat, "KsCat/D");
    
    for(int i=0; i< tree_gen->GetEntries(); i++)
    {	
        tree_gen->GetEntry(i);

        TLorentzVector Ks_p(Ks_gen[0],Ks_gen[1],Ks_gen[2],Ks_gen[3]);
        TLorentzVector pi_p(pi_gen[0],pi_gen[1],pi_gen[2],pi_gen[3]);
        TLorentzVector D_p(D_gen[0],D_gen[1],D_gen[2],D_gen[3]);
        
        m_DKs = (D_p+Ks_p).M();
        m_Dpi = (D_p+pi_p).M();
        m_Kspi = (pi_p+Ks_p).M();
        
        // 511 ,-411 ,310 ,211 
        // B0, D-, Ks, pip
        mp_DKs = 1./TMath::Pi() * TMath::ACos(2. * ( m_DKs - sqrt(pdg.sijMin(1,2))*MeV ) / (sqrt(pdg.sijMax(1,2))*MeV - sqrt(pdg.sijMin(1,2))*MeV)  - 1.);  
        mp_Dpi = 1./TMath::Pi() * TMath::ACos(2. * ( m_Dpi - sqrt(pdg.sijMin(1,3))*MeV ) / (sqrt(pdg.sijMax(1,3))*MeV - sqrt(pdg.sijMin(1,3))*MeV)  - 1.);  
        mp_Kspi = 1./TMath::Pi()* TMath::ACos(2. * ( m_Kspi - sqrt(pdg.sijMin(3,2))*MeV ) / (sqrt(pdg.sijMax(3,2))*MeV - sqrt(pdg.sijMin(3,2))*MeV)  - 1.);  

        cos_DKs = 1./TMath::Pi() * TMath::ACos(helicityCosine(D_p,Ks_p,D_p+pi_p));
        cos_Dpi = 1./TMath::Pi() * TMath::ACos(helicityCosine(D_p,pi_p,D_p+Ks_p));
        cos_Kspi = 1./TMath::Pi() * TMath::ACos(helicityCosine(Ks_p,D_p,Ks_p+pi_p));

        summary_tree_gen->Fill();
    }
    summary_tree_gen->Write();
    output_gen->Close();

    NamedParameter<string> InputSelectedMC("InputSelectedMC", (std::string) "../../../../../Selection/Final/signal_mc.root");
    TChain* tree=new TChain("DecayTree");
    tree->Add(((string)InputSelectedMC).c_str());
    double Ks[5]; 
    double pi[5]; 
    double D[5];
    int cat;
    double BDTG;
    
    tree->SetBranchAddress("B_BKGCAT",&cat);
    tree->SetBranchAddress("BDTG",&BDTG);
    
    tree->SetBranchAddress("Ks_TRUEP_X",&Ks[0]);
    tree->SetBranchAddress("Ks_TRUEP_Y",&Ks[1]);
    tree->SetBranchAddress("Ks_TRUEP_Z",&Ks[2]); 
    tree->SetBranchAddress("Ks_TRUEP_E",&Ks[3]); 
	
    tree->SetBranchAddress("pi_TRUEP_X",&pi[0]);
    tree->SetBranchAddress("pi_TRUEP_Y",&pi[1]);
    tree->SetBranchAddress("pi_TRUEP_Z",&pi[2]); 
    tree->SetBranchAddress("pi_TRUEP_E",&pi[3]); 
	
    tree->SetBranchAddress("D_TRUEP_X",&D[0]);
    tree->SetBranchAddress("D_TRUEP_Y",&D[1]);
    tree->SetBranchAddress("D_TRUEP_Z",&D[2]); 
    tree->SetBranchAddress("D_TRUEP_E",&D[3]); 
    
    TFile* output = new TFile("SelMC.root","RECREATE");
    TTree* summary_tree = tree->CloneTree(0);
    summary_tree->Branch("TRUE_m_DKs", &m_DKs, "TRUE_m_DKs/D");
    summary_tree->Branch("TRUE_m_Dpi", &m_Dpi, "TRUE_m_Dpi/D");
    summary_tree->Branch("TRUE_m_Kspi", &m_Kspi, "TRUE_m_Kspi/D");

    summary_tree->Branch("TRUE_mp_DKs", &mp_DKs, "TRUE_mp_DKs/D");
    summary_tree->Branch("TRUE_mp_Dpi", &mp_Dpi, "TRUE_mp_Dpi/D");
    summary_tree->Branch("TRUE_mp_Kspi", &mp_Kspi, "TRUE_mp_Kspi/D");
    
    summary_tree->Branch("TRUE_cos_DKs", &cos_DKs, "TRUE_cos_DKs/D");
    summary_tree->Branch("TRUE_cos_Dpi", &cos_Dpi, "TRUE_cos_Dpi/D");
    summary_tree->Branch("TRUE_cos_Kspi", &cos_Kspi, "TRUE_cos_Kspi/D");
    
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++)
    {
        tree->GetEntry(i);
        if(cat>0)continue;        
        TLorentzVector Ks_p(Ks[0],Ks[1],Ks[2],Ks[3]);
        TLorentzVector pi_p(pi[0],pi[1],pi[2],pi[3]);
        TLorentzVector D_p(D[0],D[1],D[2],D[3]);
        
        m_DKs = (D_p+Ks_p).M();
        m_Dpi = (D_p+pi_p).M();
        m_Kspi = (pi_p+Ks_p).M();
        
        mp_DKs = 1./TMath::Pi() * TMath::ACos(2. * ( m_DKs - sqrt(pdg.sijMin(1,2))*MeV ) / (sqrt(pdg.sijMax(1,2))*MeV - sqrt(pdg.sijMin(1,2))*MeV)  - 1.);  
        mp_Dpi = 1./TMath::Pi() * TMath::ACos(2. * ( m_Dpi - sqrt(pdg.sijMin(1,3))*MeV ) / (sqrt(pdg.sijMax(1,3))*MeV - sqrt(pdg.sijMin(1,3))*MeV)  - 1.);  
        mp_Kspi = 1./TMath::Pi()* TMath::ACos(2. * ( m_Kspi - sqrt(pdg.sijMin(3,2))*MeV ) / (sqrt(pdg.sijMax(3,2))*MeV - sqrt(pdg.sijMin(3,2))*MeV)  - 1.);  
        
        cos_DKs = 1./TMath::Pi() * TMath::ACos(helicityCosine(D_p,Ks_p,D_p+pi_p));
        cos_Dpi = 1./TMath::Pi() * TMath::ACos(helicityCosine(D_p,pi_p,D_p+Ks_p));
        cos_Kspi = 1./TMath::Pi() * TMath::ACos(helicityCosine(Ks_p,D_p,Ks_p+pi_p));
        
        summary_tree->Fill();
    }
    summary_tree->Write();
    output->Close();

    DalitzEventList eventList;
    TFile *file_MINT =  TFile::Open("SignalIntegrationEvents.root");
    TTree* tree_MINT=dynamic_cast<TTree*>(file_MINT->Get("DalitzEventList"));

    TFile* output_MINT = new TFile("MintMC.root","RECREATE");
    TTree* summary_tree_MINT = tree_MINT->CloneTree();
    eventList.fromNtuple(tree_MINT,1);
    cout << " I've got " << eventList.size() << " events." << endl;
    file_MINT->Close();

    TBranch* b_m_DKs = summary_tree_MINT->Branch("TRUE_m_DKs", &m_DKs, "TRUE_m_DKs/D");
    TBranch* b_m_Dpi = summary_tree_MINT->Branch("TRUE_m_Dpi", &m_Dpi, "TRUE_m_Dpi/D");
    TBranch* b_m_Kspi = summary_tree_MINT->Branch("TRUE_m_Kspi", &m_Kspi, "TRUE_m_Kspi/D");

    TBranch* b_mp_DKs = summary_tree_MINT->Branch("TRUE_mp_DKs", &mp_DKs, "TRUE_mp_DKs/D");
    TBranch* b_mp_Dpi = summary_tree_MINT->Branch("TRUE_mp_Dpi", &mp_Dpi, "TRUE_mp_Dpi/D");
    TBranch* b_mp_Kspi = summary_tree_MINT->Branch("TRUE_mp_Kspi", &mp_Kspi, "TRUE_mp_Kspi/D");
    
    TBranch* b_cos_DKs = summary_tree_MINT->Branch("TRUE_cos_DKs", &cos_DKs, "TRUE_cos_DKs/D");
    TBranch* b_cos_Dpi = summary_tree_MINT->Branch("TRUE_cos_Dpi", &cos_Dpi, "TRUE_cos_Dpi/D");
    TBranch* b_cos_Kspi = summary_tree_MINT->Branch("TRUE_cos_Kspi", &cos_Kspi, "TRUE_cos_Kspi/D");
    
    for(int i=0; i<eventList.size();i++){
        m_DKs = sqrt(eventList[i].s(1,2))*MeV;
        m_Dpi = sqrt(eventList[i].s(1,3))*MeV;
        m_Kspi = sqrt(eventList[i].s(2,3))*MeV;
        
        mp_DKs = 1./TMath::Pi() * TMath::ACos(2. * ( m_DKs - sqrt(pdg.sijMin(1,2))*MeV ) / (sqrt(pdg.sijMax(1,2))*MeV - sqrt(pdg.sijMin(1,2))*MeV)  - 1.);  
        mp_Dpi = 1./TMath::Pi() * TMath::ACos(2. * ( m_Dpi - sqrt(pdg.sijMin(1,3))*MeV ) / (sqrt(pdg.sijMax(1,3))*MeV - sqrt(pdg.sijMin(1,3))*MeV)  - 1.);  
        mp_Kspi = 1./TMath::Pi()* TMath::ACos(2. * ( m_Kspi - sqrt(pdg.sijMin(3,2))*MeV ) / (sqrt(pdg.sijMax(3,2))*MeV - sqrt(pdg.sijMin(3,2))*MeV)  - 1.);  

        TLorentzVector Ks_p = eventList[i].p(2);
        TLorentzVector pi_p = eventList[i].p(3);
        TLorentzVector D_p = eventList[i].p(1);
        
        cos_DKs = 1./TMath::Pi() * TMath::ACos(helicityCosine(D_p,Ks_p,D_p+pi_p));
        cos_Dpi = 1./TMath::Pi() * TMath::ACos(helicityCosine(D_p,pi_p,D_p+Ks_p));
        cos_Kspi = 1./TMath::Pi() * TMath::ACos(helicityCosine(Ks_p,D_p,Ks_p+pi_p));

        summary_tree_MINT->GetEntry(i);
        b_m_DKs->Fill();
        b_m_Dpi->Fill();
        b_m_Kspi->Fill();

        b_mp_DKs->Fill();
        b_mp_Dpi->Fill();
        b_mp_Kspi->Fill();

        b_cos_DKs->Fill();
        b_cos_Dpi->Fill();
        b_cos_Kspi->Fill();        
    }
    summary_tree_MINT->Write();
    output_MINT->Close();
}

void reweightGen(){
	
    Double_t m_Dpi;
    Double_t m_Kspi;
    Double_t m_DKs;
    Double_t mp_Dpi;
    Double_t mp_Kspi;
    Double_t mp_DKs;
    Double_t cos_Dpi;
    Double_t cos_Kspi;
    Double_t cos_DKs;
    
    double BDTG,BDTG_LL,BDTG_DD;
    double weight = 1;
    double genPdf = 1;
        
    TChain* tree_gen=new TChain("MCDecayTree");
    tree_gen->Add("GenMC_BDT.root");
    tree_gen->SetBranchAddress( "TRUE_m_Dpi", &m_Dpi );
    tree_gen->SetBranchAddress( "TRUE_m_Kspi", &m_Kspi );
    tree_gen->SetBranchAddress( "TRUE_m_DKs", &m_DKs );
    tree_gen->SetBranchAddress( "TRUE_mp_Dpi", &mp_Dpi );
    tree_gen->SetBranchAddress( "TRUE_mp_Kspi", &mp_Kspi );
    tree_gen->SetBranchAddress( "TRUE_mp_DKs", &mp_DKs );
    tree_gen->SetBranchAddress( "TRUE_cos_Dpi", &cos_Dpi );
    tree_gen->SetBranchAddress( "TRUE_cos_Kspi", &cos_Kspi );
    tree_gen->SetBranchAddress( "TRUE_cos_DKs", &cos_DKs );
    tree_gen->SetBranchAddress("eff_BDTG",&BDTG);
    tree_gen->SetBranchAddress("eff_BDTG_LL",&BDTG_LL);
    tree_gen->SetBranchAddress("eff_BDTG_DD",&BDTG_DD);

    TChain* tree=new TChain("DecayTree");
    tree->Add("SelMC_BDT.root");
    tree->SetBranchAddress("TRUE_m_Dpi", &m_Dpi );
    tree->SetBranchAddress("TRUE_m_Kspi", &m_Kspi );
    tree->SetBranchAddress("TRUE_m_DKs", &m_DKs );
    tree->SetBranchAddress("TRUE_mp_Dpi", &mp_Dpi );
    tree->SetBranchAddress("TRUE_mp_Kspi", &mp_Kspi );
    tree->SetBranchAddress("TRUE_mp_DKs", &mp_DKs );
    tree->SetBranchAddress("TRUE_cos_Dpi", &cos_Dpi );
    tree->SetBranchAddress("TRUE_cos_Kspi", &cos_Kspi );
    tree->SetBranchAddress("TRUE_cos_DKs", &cos_DKs );
    tree->SetBranchAddress("eff_BDTG",&BDTG);
    tree->SetBranchAddress("eff_BDTG_LL",&BDTG_LL);
    tree->SetBranchAddress("eff_BDTG_DD",&BDTG_DD);

    TChain* tree_MINT=new TChain("DalitzEventList");
    tree_MINT->Add("MintMC_BDT.root");
    tree_MINT->SetBranchAddress( "TRUE_m_Dpi", &m_Dpi );
    tree_MINT->SetBranchAddress( "TRUE_m_Kspi", &m_Kspi );
    tree_MINT->SetBranchAddress( "TRUE_m_DKs", &m_DKs );    
    tree_MINT->SetBranchAddress( "TRUE_mp_Dpi", &mp_Dpi );
    tree_MINT->SetBranchAddress( "TRUE_mp_Kspi", &mp_Kspi );
    tree_MINT->SetBranchAddress( "TRUE_mp_DKs", &mp_DKs );
    tree_MINT->SetBranchAddress( "TRUE_cos_Dpi", &cos_Dpi );
    tree_MINT->SetBranchAddress( "TRUE_cos_Kspi", &cos_Kspi );
    tree_MINT->SetBranchAddress( "TRUE_cos_DKs", &cos_DKs );
    tree_MINT->SetBranchAddress("eff_BDTG",&BDTG);
    tree_MINT->SetBranchAddress("eff_BDTG_LL",&BDTG_LL);
    tree_MINT->SetBranchAddress("eff_BDTG_DD",&BDTG_DD);
    tree_MINT->SetBranchAddress("weight",&weight);
    tree_MINT->SetBranchAddress("genPdf",&genPdf);

    double nBins = 50;
    TH1D* h_DKs = new TH1D("m_DKs","; m(DK_{s}) [GeV]; Yield",nBins,2,5.5);
    TH1D* h_Dpi = new TH1D("m_Dpi","; m(D#pi) [GeV]; Yield",nBins,1,5.5);
    TH1D* h_Kspi = new TH1D("m_Kspi","; m(K_{s}#pi) [GeV]; Yield",nBins,0,4);
    
    TH1D* h_DKs_gen = (TH1D*) h_DKs->Clone("m_DKs_gen");
    TH1D* h_Dpi_gen = (TH1D*) h_Dpi->Clone("m_Dpi_gen");
    TH1D* h_Kspi_gen = (TH1D*) h_Kspi->Clone("m_Kspi_gen");
    TH1D* h_DKs_gen_rw = (TH1D*) h_DKs->Clone("m_DKs_gen_rw");
    TH1D* h_Dpi_gen_rw = (TH1D*) h_Dpi->Clone("m_Dpi_gen_rw");
    TH1D* h_Kspi_gen_rw = (TH1D*) h_Kspi->Clone("m_Kspi_gen_rw");

    TH1D* h_DKs_MINT = (TH1D*) h_DKs->Clone("m_DKs_MINT");
    TH1D* h_Dpi_MINT = (TH1D*) h_Dpi->Clone("m_Dpi_MINT");
    TH1D* h_Kspi_MINT = (TH1D*) h_Kspi->Clone("m_Kspi_MINT");
    TH1D* h_DKs_MINT_rw = (TH1D*) h_DKs->Clone("m_DKs_MINT_rw");
    TH1D* h_Dpi_MINT_rw = (TH1D*) h_Dpi->Clone("m_Dpi_MINT_rw");
    TH1D* h_Kspi_MINT_rw = (TH1D*) h_Kspi->Clone("m_Kspi_MINT_rw");

    TH1D* bdt_gen= new TH1D("","; BTDG; Efficiency (norm.)",nBins,-0.65,0.25);
    TH1D* bdt_sel= new TH1D("","; BTDG; Events (norm.)",nBins,-0.65,0.25);
    TH1D* bdt_MINT= new TH1D("","",nBins,-0.65,0.25);

    for(int i=0; i< tree_gen->GetEntries(); i++)
    {	
        tree_gen->GetEntry(i);
        bdt_gen->Fill(BDTG);
    }

    for(int i=0; i< tree->GetEntries(); i++)
    {	
        tree->GetEntry(i);
        bdt_sel->Fill(BDTG);
        h_DKs->Fill(m_DKs/1000.);
        h_Dpi->Fill(m_Dpi/1000.);
        h_Kspi->Fill(m_Kspi/1000.);
    }
    	
    TH1D *h_weight = (TH1D*)bdt_gen->Clone();
    h_weight->Divide(bdt_sel,bdt_gen);
    h_weight->Scale(1./h_weight->Integral());
    TCanvas* c = new TCanvas();
    h_weight->Draw();
    c->Print("eff.eps");

    for(int i=0; i< tree_gen->GetEntries(); i++)
    {	
        tree_gen->GetEntry(i);
        h_DKs_gen->Fill(m_DKs/1000.);
        h_Dpi_gen->Fill(m_Dpi/1000.);
        h_Kspi_gen->Fill(m_Kspi/1000.);

        double weight = h_weight->GetBinContent(h_weight->FindBin(BDTG));	
        h_DKs_gen_rw->Fill(m_DKs/1000.,weight);
        h_Dpi_gen_rw->Fill(m_Dpi/1000.,weight);
        h_Kspi_gen_rw->Fill(m_Kspi/1000.,weight);
    }

    DalitzEventPattern pdg(511 ,-411 ,310 ,211);
    TFile *file_integ =  TFile::Open("SignalIntegrationEvents.root");
    TTree* tree_integ=dynamic_cast<TTree*>(file_integ->Get("DalitzEventList"));
    DalitzEventList eventList,eventList_rw,eventList_rw_CP;
    eventList.fromNtuple(tree_integ,1);
    eventList[0].print();
    if(tree_MINT->GetEntries() != tree_integ->GetEntries()) throw "ERROR";

    TFile* out = new TFile("SignalIntegrationEvents_AccBDT.root","RECREATE");
    TTree* out_tree = new TTree("DalitzEventList","DalitzEventList");
    double D[4]; 
    double Ks[4]; 
    double pi[4]; 
    
    out_tree->Branch("weight",&weight);
    out_tree->Branch("genPdf",&genPdf);
    out_tree->Branch("Ks_PX",&Ks[0]);
    out_tree->Branch("Ks_PY",&Ks[1]);
    out_tree->Branch("Ks_PZ",&Ks[2]); 
    out_tree->Branch("Ks_PE",&Ks[3]); 
    out_tree->Branch("pi_PX",&pi[0]);
    out_tree->Branch("pi_PY",&pi[1]);
    out_tree->Branch("pi_PZ",&pi[2]); 
    out_tree->Branch("pi_PE",&pi[3]); 
    out_tree->Branch("D_PX",&D[0]);
    out_tree->Branch("D_PY",&D[1]);
    out_tree->Branch("D_PZ",&D[2]); 
    out_tree->Branch("D_PE",&D[3]); 
    
    for(int i=0; i< tree_MINT->GetEntries(); i++)
    {	
        tree_MINT->GetEntry(i);
        bdt_MINT->Fill(BDTG);

        h_DKs_MINT->Fill(m_DKs/1000.);
        h_Dpi_MINT->Fill(m_Dpi/1000.);
        h_Kspi_MINT->Fill(m_Kspi/1000.);

        double w = h_weight->GetBinContent(h_weight->FindBin(BDTG));	
        h_DKs_MINT_rw->Fill(m_DKs/1000.,w);
        h_Dpi_MINT_rw->Fill(m_Dpi/1000.,w);
        h_Kspi_MINT_rw->Fill(m_Kspi/1000.,w);

        DalitzEvent evt(eventList[i]);
        evt.setWeight(evt.getWeight()*w);

        //if(evt.getGeneratorPdfRelativeToPhaseSpace()>1000.)continue;

        eventList_rw.Add(evt);
        //evt.CP_conjugateYourself();
        //eventList_rw_CP.Add(evt);
        
        weight *= w;
        D[0] = evt.p(1).X();
        D[1] = evt.p(1).Y();
        D[2] = evt.p(1).Z();
        D[3] = evt.p(1).E();
        
        Ks[0] = evt.p(2).X();
        Ks[1] = evt.p(2).Y();
        Ks[2] = evt.p(2).Z();
        Ks[3] = evt.p(2).E();
        
        pi[0] = evt.p(3).X();
        pi[1] = evt.p(3).Y();
        pi[2] = evt.p(3).Z();
        pi[3] = evt.p(3).E();
                
        for(int j=0; j<4; j++){
            D[j] /= 1000.;
            Ks[j] /= 1000.;
            pi[j] /= 1000.;
        }
        out_tree->Fill();
    }
    
    eventList_rw.saveAsNtuple("MINT_SignalIntegrationEvents_AccBDT.root");
    //eventList_rw_CP.saveAsNtuple("SignalIntegrationEvents_AccBDT_CP2.root");

    bdt_sel->SetLineColor(kBlue);
    bdt_sel->DrawNormalized("hist",1);
    bdt_gen->SetLineColor(kRed);
    bdt_gen->SetMarkerColor(kRed);
    bdt_gen->DrawNormalized("same",1);
    //bdt_MINT->DrawNormalized("histsame",1);
    c->Print("BDTG.eps");

    h_DKs->DrawNormalized("",1);
    h_DKs_gen->SetLineColor(kRed);
    h_DKs_gen->DrawNormalized("histsame",1);
    h_DKs_gen_rw->SetLineColor(kBlue);
    h_DKs_gen_rw->DrawNormalized("histsame",1);
    c->Print("m_DKs.eps");
    TH1D* h_DKs_eff = (TH1D*)h_DKs->Clone();
    h_DKs_eff->Divide(h_DKs,h_DKs_gen);
    h_DKs_eff->GetYaxis()->SetTitle("Efficiency (norm.)");
    h_DKs_eff->DrawNormalized("e",1);
    TH1D* h_DKs_MINT_eff = (TH1D*)h_DKs_MINT->Clone();
    h_DKs_MINT_eff->Divide(h_DKs_MINT_rw,h_DKs_MINT);
    h_DKs_MINT_eff->SetLineColor(kBlue);
    h_DKs_MINT_eff->DrawNormalized("histsame",1);
    c->Print("eff_m_DKs.eps");

    h_Dpi->DrawNormalized("",1);
    h_Dpi_gen->SetLineColor(kRed);
    h_Dpi_gen->DrawNormalized("histsame",1);
    h_Dpi_gen_rw->SetLineColor(kBlue);
    h_Dpi_gen_rw->DrawNormalized("histsame",1);
    c->Print("m_Dpi.eps");
    TH1D* h_Dpi_eff = (TH1D*)h_Dpi->Clone();
    h_Dpi_eff->Divide(h_Dpi,h_Dpi_gen);
    h_Dpi_eff->GetYaxis()->SetTitle("Efficiency (norm.)");
    h_Dpi_eff->DrawNormalized("e",1);
    TH1D* h_Dpi_MINT_eff = (TH1D*)h_Dpi_MINT->Clone();
    h_Dpi_MINT_eff->Divide(h_Dpi_MINT_rw,h_Dpi_MINT);
    h_Dpi_MINT_eff->SetLineColor(kBlue);
    h_Dpi_MINT_eff->DrawNormalized("histsame",1);
    c->Print("eff_m_Dpi.eps");

    h_Kspi->DrawNormalized("",1);
    h_Kspi_gen->SetLineColor(kRed);
    h_Kspi_gen->DrawNormalized("histsame",1);
    h_Kspi_gen_rw->SetLineColor(kBlue);
    h_Kspi_gen_rw->DrawNormalized("histsame",1);
    c->Print("m_Kspi.eps");
    TH1D* h_Kspi_eff = (TH1D*)h_Kspi->Clone();
    h_Kspi_eff->Divide(h_Kspi,h_Kspi_gen);
    h_Kspi_eff->GetYaxis()->SetTitle("Efficiency (norm.)");
    h_Kspi_eff->DrawNormalized("e",1);
    TH1D* h_Kspi_MINT_eff = (TH1D*)h_Kspi_MINT->Clone();
    h_Kspi_MINT_eff->Divide(h_Kspi_MINT_rw,h_Kspi_MINT);
    h_Kspi_MINT_eff->SetLineColor(kBlue);
    h_Kspi_MINT_eff->DrawNormalized("histsame",1);
    c->Print("eff_m_Kspi.eps");
}

int main(int argc, char** argv){

    time_t startTime = time(0);
    TH1::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    gStyle->SetPalette(1);

    //prepareFilesForBDT();
    reweightGen();

    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
