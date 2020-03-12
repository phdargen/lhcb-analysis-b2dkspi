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
    
    for(int i=0; i< tree_gen->GetEntries(); i++)
    {	
        tree_gen->GetEntry(i);

        TLorentzVector Ks_p(Ks_gen[0],Ks_gen[1],Ks_gen[2],Ks_gen[3]);
        TLorentzVector pi_p(pi_gen[0],pi_gen[1],pi_gen[2],pi_gen[3]);
        TLorentzVector D_p(D_gen[0],D_gen[1],D_gen[2],D_gen[3]);
        
        m_DKs = (D_p+Ks_p).M();
        m_Dpi = (D_p+pi_p).M();
        m_Kspi = (pi_p+Ks_p).M();
        
        summary_tree_gen->Fill();
    }
    summary_tree_gen->Write();
    output_gen->Close();

    NamedParameter<string> InputSelectedMC("InputSelectedMC", (std::string) "../../../../../Selection/BDT/signal_mc.root");
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
        
        summary_tree->Fill();
    }
    summary_tree->Write();
    output->Close();

    DalitzEventPattern pdg(511 ,-411 ,310 ,211);
    DalitzEventList eventList;
    TFile *file_MINT =  TFile::Open("../LASSO/SignalIntegrationEvents.root");
    TTree* tree_MINT=dynamic_cast<TTree*>(file_MINT->Get("DalitzEventList"));

    TFile* output_MINT = new TFile("MintMC.root","RECREATE");
    TTree* summary_tree_MINT = tree_MINT->CloneTree();
    eventList.fromNtuple(tree_MINT,1);
    cout << " I've got " << eventList.size() << " events." << endl;
    file_MINT->Close();

    TBranch* b_m_DKs = summary_tree_MINT->Branch("TRUE_m_DKs", &m_DKs, "TRUE_m_DKs/D");
    TBranch* b_m_Dpi = summary_tree_MINT->Branch("TRUE_m_Dpi", &m_Dpi, "TRUE_m_Dpi/D");
    TBranch* b_m_Kspi = summary_tree_MINT->Branch("TRUE_m_Kspi", &m_Kspi, "TRUE_m_Kspi/D");

    for(int i=0; i<eventList.size();i++){
        m_DKs = sqrt(eventList[i].s(1,2))*MeV;
        m_Dpi = sqrt(eventList[i].s(1,3))*MeV;
        m_Kspi = sqrt(eventList[i].s(2,3))*MeV;

        summary_tree_MINT->GetEntry(i);
        b_m_DKs->Fill();
        b_m_Dpi->Fill();
        b_m_Kspi->Fill();
    }
    summary_tree_MINT->Write();
    output_MINT->Close();
}

void reweightGen(){
	
    Double_t m_Dpi;
    Double_t m_Kspi;
    Double_t m_DKs;
    double BDTG;
        
    TChain* tree_gen=new TChain("MCDecayTree");
    tree_gen->Add("GenMC_BDT.root");
    tree_gen->SetBranchAddress( "TRUE_m_Dpi", &m_Dpi );
    tree_gen->SetBranchAddress( "TRUE_m_Kspi", &m_Kspi );
    tree_gen->SetBranchAddress( "TRUE_m_DKs", &m_DKs );
    tree_gen->SetBranchAddress("eff_BDTG",&BDTG);

    TChain* tree=new TChain("DecayTree");
    tree->Add("SelMC_BDT.root");
    tree->SetBranchAddress( "TRUE_m_Dpi", &m_Dpi );
    tree->SetBranchAddress( "TRUE_m_Kspi", &m_Kspi );
    tree->SetBranchAddress( "TRUE_m_DKs", &m_DKs );
    tree->SetBranchAddress("eff_BDTG",&BDTG);

    TChain* tree_MINT=new TChain("DalitzEventList");
    tree_MINT->Add("MintMC_BDT.root");
    tree_MINT->SetBranchAddress( "TRUE_m_Dpi", &m_Dpi );
    tree_MINT->SetBranchAddress( "TRUE_m_Kspi", &m_Kspi );
    tree_MINT->SetBranchAddress( "TRUE_m_DKs", &m_DKs );    
    tree_MINT->SetBranchAddress("eff_BDTG",&BDTG);

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
    TFile *file_integ =  TFile::Open("../LASSO/SignalIntegrationEvents.root");
    TTree* tree_integ=dynamic_cast<TTree*>(file_integ->Get("DalitzEventList"));
    DalitzEventList eventList,eventList_rw,eventList_rw_CP;
    eventList.fromNtuple(tree_integ,1);
    eventList[0].print();
    if(tree_MINT->GetEntries() != tree_integ->GetEntries()) throw "ERROR";

    for(int i=0; i< tree_MINT->GetEntries(); i++)
    {	
        tree_MINT->GetEntry(i);
        bdt_MINT->Fill(BDTG);

        h_DKs_MINT->Fill(m_DKs/1000.);
        h_Dpi_MINT->Fill(m_Dpi/1000.);
        h_Kspi_MINT->Fill(m_Kspi/1000.);

        double weight = h_weight->GetBinContent(h_weight->FindBin(BDTG));	
        h_DKs_MINT_rw->Fill(m_DKs/1000.,weight);
        h_Dpi_MINT_rw->Fill(m_Dpi/1000.,weight);
        h_Kspi_MINT_rw->Fill(m_Kspi/1000.,weight);

        DalitzEvent evt(eventList[i]);
        evt.setWeight(evt.getWeight()*weight);
        eventList_rw.Add(evt);
        evt.CP_conjugateYourself();
        eventList_rw_CP.Add(evt);
    }
    
    eventList_rw.saveAsNtuple("SignalIntegrationEvents_AccBDT.root");
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
