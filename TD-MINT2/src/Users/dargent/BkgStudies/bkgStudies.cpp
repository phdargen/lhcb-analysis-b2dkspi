#include <cmath>
#include <algorithm>
#include <iostream>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMath.h>
#include <TFitResult.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TMultiGraph.h>
#include <TNtuple.h>
#include "TRandom3.h"
#include <sstream>
#include <ctime>
#include "Mint/NamedParameter.h"
#include "Mint/Utils.h"

using namespace std;
using namespace MINT;

void BkgStudies(){
    
    TChain* tree = 0;
    tree =new TChain("DecayTree");    
    tree->Add("/eos/lhcb/user/p/phdargen/B2DKspi/BDT/B2DKspi_data.root");
    //tree->Add("../../../../../Selection/Preselected/Data_B2DKspi_LL_16.root");
    //tree->Add("../MassFits/b2dkspi_sw.root");

    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("*M",1) ;
    tree->SetBranchStatus("*MM",1) ;
    tree->SetBranchStatus( "*P", 1 );
    tree->SetBranchStatus( "*PX", 1 );
    tree->SetBranchStatus( "*PY", 1);
    tree->SetBranchStatus( "*PZ", 1);
    tree->SetBranchStatus( "*PE", 1);
    tree->SetBranchStatus( "*PT", 1 );
    tree->SetBranchStatus( "*PIDK*", 1);
    tree->SetBranchStatus( "*PIDp*", 1);
    tree->SetBranchStatus( "*ProbNN*", 1);
    
    //define Lorentz vectors
    TLorentzVector pi,pip_Ks,pim_Ks,pi1_D,pi2_D,K_D,Ks,D;
    TLorentzVector FullDTF_pi,FullDTF_pip_Ks,FullDTF_pim_Ks,FullDTF_pi1_D,FullDTF_pi2_D,FullDTF_K_D,FullDTF_Ks,FullDTF_D;
    TLorentzVector FullBsDTF_pi,FullBsDTF_pip_Ks,FullBsDTF_pim_Ks,FullBsDTF_pi1_D,FullBsDTF_pi2_D,FullBsDTF_K_D,FullBsDTF_Ks,FullBsDTF_D;
    
    //misIDs
    TLorentzVector K_fromD_asP_MissID,K_fromD_asPi_MissID;
    TLorentzVector pi_asP_MissID,pi_asK_MissID;
    TLorentzVector pi1_fromD_asP_MissID,pi1_fromD_asK_MissID;
    TLorentzVector pi2_fromD_asP_MissID,pi2_fromD_asK_MissID;

    //masses
    double massKaon = 493.68;
    double massPion = 139.57;
    double massProton = 938.27;
    double massPhi = 1019.46;
    double massKstar = 895.81;
    double massDs = 1968.30;
    double massDminus = 1869.61;
    double massLambda_c = 2286.46;
    double massD0 = 1864.83;
    double massKs = 497.611;
    double massB0 = 5279.65;
    double massBs = 5366.88;
    double massBp = 5279.34;

    
    //define observables
    Double_t pi_PX;
    Double_t pi_PY;
    Double_t pi_PZ;
    Double_t pip_Ks_PX;
    Double_t pip_Ks_PY;
    Double_t pip_Ks_PZ;
    Double_t pim_Ks_PX;
    Double_t pim_Ks_PY;
    Double_t pim_Ks_PZ;
    
    Double_t pi1_D_PX;
    Double_t pi1_D_PY;
    Double_t pi1_D_PZ;
    Double_t pi2_D_PX;
    Double_t pi2_D_PY;
    Double_t pi2_D_PZ;
    Double_t K_D_PX;
    Double_t K_D_PY;
    Double_t K_D_PZ;
    
    Double_t pi_minus_fromDs_PIDp;
    Double_t K_minus_fromDs_PIDK;
    Double_t K_minus_fromDs_PIDp;
    Double_t B_MM;
    
    //link to tree
    tree->SetBranchAddress( "pi_PX" , &pi_PX );
    tree->SetBranchAddress( "pi_PY" , &pi_PY );
    tree->SetBranchAddress( "pi_PZ" , &pi_PZ );
    tree->SetBranchAddress( "pip_Ks_PX" , &pip_Ks_PX );
    tree->SetBranchAddress( "pip_Ks_PY" , &pip_Ks_PY );
    tree->SetBranchAddress( "pip_Ks_PZ" , &pip_Ks_PZ );
    tree->SetBranchAddress( "pim_Ks_PX" , &pim_Ks_PX );
    tree->SetBranchAddress( "pim_Ks_PY" , &pim_Ks_PY );
    tree->SetBranchAddress( "pim_Ks_PZ" , &pim_Ks_PZ );   
    tree->SetBranchAddress( "pi1_D_PX" , &pi1_D_PX );
    tree->SetBranchAddress( "pi1_D_PY" , &pi1_D_PY );
    tree->SetBranchAddress( "pi1_D_PZ" , &pi1_D_PZ );    
    tree->SetBranchAddress( "pi2_D_PX" , &pi2_D_PX );
    tree->SetBranchAddress( "pi2_D_PZ" , &pi2_D_PZ );
    tree->SetBranchAddress( "pi2_D_PY" , &pi2_D_PY );
    tree->SetBranchAddress( "K_D_PX" , &K_D_PX );
    tree->SetBranchAddress( "K_D_PY" , &K_D_PY );
    tree->SetBranchAddress( "K_D_PZ" , &K_D_PZ );
    
    Double_t pi_ProbNNpi,pi_ProbNNk,pi_ProbNNp;
    Double_t pi1_D_ProbNNpi,pi1_D_ProbNNk,pi1_D_ProbNNp;
    Double_t pi2_D_ProbNNpi,pi2_D_ProbNNk,pi2_D_ProbNNp;
    tree->SetBranchAddress( "pi_ProbNNpi" , &pi_ProbNNpi );
    tree->SetBranchAddress( "pi_ProbNNk" , &pi_ProbNNk );
    tree->SetBranchAddress( "pi_ProbNNp" , &pi_ProbNNp );
    tree->SetBranchAddress( "pi1_D_ProbNNpi" , &pi1_D_ProbNNpi );
    tree->SetBranchAddress( "pi1_D_ProbNNk" , &pi1_D_ProbNNk );
    tree->SetBranchAddress( "pi1_D_ProbNNp" , &pi1_D_ProbNNp );
    tree->SetBranchAddress( "pi2_D_ProbNNpi" , &pi2_D_ProbNNpi );
    tree->SetBranchAddress( "pi2_D_ProbNNk" , &pi2_D_ProbNNk );
    tree->SetBranchAddress( "pi2_D_ProbNNp" , &pi2_D_ProbNNp );

    tree->SetBranchAddress( "B_MM" , &B_MM );
    
    //misID histos
    
    // D-
    TH1D* mass_D2KKpi_as_D2Kpipi_beforeVeto = new TH1D(" ", ";m(K_{#pi}K^{-}#pi) - m_{D_{s}} [MeV/c^{2}];Events", 50, -100,100);
    TH1D* mass_D2KKpi_as_D2Kpipi_afterVeto = new TH1D(" ", ";m(K_{#pi}K^{-}#pi) - m_{D} [MeV/c^{2}];Events", 50, -100,100);

    TH1D* mass_D2pipipi_as_D2Kpipi_beforeVeto = new TH1D(" ", ";m(#pi^{-}_{K}#pi^{+}#pi^{+}) - m_{D_{s}} [MeV/c^{2}];Events", 50, -150,150);
    TH1D* mass_D2pipipi_as_D2Kpipi_afterVeto = new TH1D(" ", ";m(#pi^{-}_{K}#pi^{+}#pi^{+}) - m_{D_{s}} [MeV/c^{2}];Events", 50, -150,150);

    TH1D* mass_Lc2ppipi_as_D2Kpipi_beforeVeto = new TH1D(" ", ";m(p_{K}#pi^{+}#pi^{+}) - m_{#Lambda_{c}} [MeV/c^{2}];Events", 50, -150,150);
    TH1D* mass_Lc2ppipi_as_D2Kpipi_afterVeto = new TH1D(" ", ";m(p_{K}#pi^{+}#pi^{+}) - m_{#Lambda_{c}} [MeV/c^{2}];Events", 50, -150,150);

    TH1D* mass_Lc2Kppi_as_D2Kpipi_beforeVeto = new TH1D(" ", ";m(p_{#pi}K^{-}#pi^{+}) - m_{#Lambda_{c}} [MeV/c^{2}];Events", 50, -150,150);
    TH1D* mass_Lc2Kppi_as_D2Kpipi_afterVeto = new TH1D(" ", ";m(p_{#pi}K^{-}#pi^{+}) - m_{#Lambda_{c}} [MeV/c^{2}];Events", 50, -150,150);

    TH1D* mass_D2Kpi_as_D2Kpipi_beforeVeto = new TH1D(" ", ";m(K^{-}#pi^{+}) - m_{D^{0}} [MeV/c^{2}];Events", 50, -100,100);
    TH1D* mass_D2Kpi_as_D2Kpipi_afterVeto = new TH1D(" ", ";m(K^{-}#pi^{+}) - m_{D^{0}} [MeV/c^{2}];Events", 50, -100,100);

    //
    TH1D* mass_phi_as_Kpi_beforeVeto = new TH1D(" ", ";m(K^{-}#pi_{K}^{+}) - m_{#phi(1020)} [MeV/c^{2}];Events", 50, -100,100);
    TH1D* mass_phi_as_Kpi_afterVeto = new TH1D(" ", ";m(K^{-}#pi_{K}^{+}) - m_{#phi(1020)} [MeV/c^{2}];Events", 50, -100,100);

    TH1D* mass_phi_as_Kpi_beforeVeto2 = new TH1D(" ", ";m(K^{-}#pi_{K}^{+}) - m_{#phi(1020)} [MeV/c^{2}];Events", 50, -100,100);
    TH1D* mass_phi_as_Kpi_afterVeto2 = new TH1D(" ", ";m(K^{-}#pi_{K}^{+}) - m_{#phi(1020)} [MeV/c^{2}];Events", 50, -100,100);

    // B   
    TH1D* mass_B2DKsK_as_DKspi_beforeVeto = new TH1D(" ", ";m(D K_{s} K_{#pi}) - m_{B_{s}^{0}} [MeV/c^{2}];Events", 50, -100,200);
    TH1D* mass_B2DKsK_as_DKspi_afterVeto = new TH1D(" ", ";m(D K_{s} K_{#pi}) - m_{B^{0}} [MeV/c^{2}];Events", 50, -100,200);

    TH1D* mass_B2Dpi_as_DKspi_beforeVeto = new TH1D(" ", ";m(D^{-}#pi^{+}) - m_{B^{0}} [MeV/c^{2}];Events", 50, -600,100);
    TH1D* mass_B2Dpi_as_DKspi_afterVeto = new TH1D(" ", ";m(D^{-}#pi^{+}) - m_{B^{0}} [MeV/c^{2}];Events", 50, -600,100);

    TH1D* mass_B2DK_as_DKspi_beforeVeto = new TH1D(" ", ";m(D^{-}K_{S}) - m_{B^{+}} [MeV/c^{2}];Events", 50, -100,100);
    TH1D* mass_B2DK_as_DKspi_afterVeto = new TH1D(" ", ";m(D^{-}K_{S}) - m_{B^{+}} [MeV/c^{2}];Events", 50, -100,100);

    TH1D* mass_D2Kspi_as_Kspi_beforeVeto = new TH1D(" ", ";m(K_{S}#pi^{+}) - m_{D^{+}} [MeV/c^{2}];Events", 50, -200,200);
    TH1D* mass_D2Kspi_as_Kspi_afterVeto = new TH1D(" ", ";m(K_{S}#pi^{+}) - m_{D^{+}} [MeV/c^{2}];Events", 50, -200,200);

    TH1D* mass_Lc2pKs_as_Kspi_beforeVeto = new TH1D(" ", ";m(K_{S}p_{#pi}^{+}) - m_{#Lambda_{c}} [MeV/c^{2}];Events", 50, -100,100);
    TH1D* mass_Lc2pKs_as_Kspi_afterVeto = new TH1D(" ", ";m(K_{S}p_{#pi}^{+}) - m_{#Lambda_{c}} [MeV/c^{2}];Events", 50, -100,100);

        
    //loop over events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++)
    {
        if (0ul == (i % 10000ul)) std::cout << "Read event " << i << "/" << numEvents << std::endl;
        tree->GetEntry(i);
        //if(B_MM < 5200 || B_MM > 5700) continue;
        
        //fill the Lorentz vectors
        pi.SetXYZM(pi_PX,pi_PY,pi_PZ,massPion);        
        pip_Ks.SetXYZM(pip_Ks_PX,pip_Ks_PY,pip_Ks_PZ,massPion);
        pim_Ks.SetXYZM(pim_Ks_PX,pim_Ks_PY,pim_Ks_PZ,massPion);
        
        pi1_D.SetXYZM(pi1_D_PX,pi1_D_PY,pi1_D_PZ,massPion);
        pi2_D.SetXYZM(pi2_D_PX,pi2_D_PY,pi2_D_PZ,massPion);
        K_D.SetXYZM(K_D_PX,K_D_PY,K_D_PZ,massKaon);
        
        Ks = pip_Ks + pim_Ks;
        D = pi1_D + pi2_D + K_D;
         
        //misID Lorentz vectors
        K_fromD_asP_MissID.SetXYZM(K_D_PX,K_D_PY,K_D_PZ, massProton);
        K_fromD_asPi_MissID.SetXYZM(K_D_PX,K_D_PY,K_D_PZ,massPion);
        
        pi_asP_MissID.SetXYZM(pi_PX,pi_PY,pi_PZ,massProton);        
        pi_asK_MissID.SetXYZM(pi_PX,pi_PY,pi_PZ,massKaon);        

        pi1_fromD_asP_MissID.SetXYZM(pi1_D_PX,pi1_D_PY,pi1_D_PZ,massProton);        
        pi1_fromD_asK_MissID.SetXYZM(pi1_D_PX,pi1_D_PY,pi1_D_PZ,massKaon);        
        
        pi2_fromD_asP_MissID.SetXYZM(pi2_D_PX,pi2_D_PY,pi2_D_PZ,massProton);        
        pi2_fromD_asK_MissID.SetXYZM(pi2_D_PX,pi2_D_PY,pi2_D_PZ,massKaon);        
        
        
        //fill histos without veto cuts
        mass_D2KKpi_as_D2Kpipi_beforeVeto->Fill((K_D + pi1_fromD_asK_MissID + pi2_D).M()-massDs);
        mass_D2KKpi_as_D2Kpipi_beforeVeto->Fill((K_D + pi2_fromD_asK_MissID + pi1_D).M()-massDs);         
        
        mass_D2pipipi_as_D2Kpipi_beforeVeto->Fill( (K_fromD_asPi_MissID + pi1_D + pi2_D).M() - massDs);
        
        mass_Lc2ppipi_as_D2Kpipi_beforeVeto->Fill( (K_fromD_asP_MissID + pi1_D + pi2_D).M() - massLambda_c );

        mass_Lc2Kppi_as_D2Kpipi_beforeVeto->Fill( (K_D + pi1_fromD_asP_MissID + pi2_D).M() - massLambda_c );
        mass_Lc2Kppi_as_D2Kpipi_beforeVeto->Fill( (K_D + pi2_fromD_asP_MissID + pi1_D).M() - massLambda_c );
        
        mass_D2Kpi_as_D2Kpipi_beforeVeto->Fill( (K_D + pi1_D).M() - massD0 );
        mass_D2Kpi_as_D2Kpipi_beforeVeto->Fill( (K_D + pi2_D).M() - massD0 );
        
        
        mass_phi_as_Kpi_beforeVeto->Fill( (K_D+pi1_fromD_asK_MissID).M() - massPhi );
        mass_phi_as_Kpi_beforeVeto->Fill( (K_D+pi2_fromD_asK_MissID).M() - massPhi );

        mass_phi_as_Kpi_beforeVeto2->Fill( (K_D+pi_asK_MissID).M() - massPhi );

        
        mass_B2DKsK_as_DKspi_beforeVeto->Fill( (Ks+D+pi_asK_MissID).M() - massBs );
        
        mass_B2Dpi_as_DKspi_beforeVeto->Fill( (D+pi).M() - massB0 );
        
        mass_B2DK_as_DKspi_beforeVeto->Fill( (D+Ks).M() - massB0 );
        
        mass_D2Kspi_as_Kspi_beforeVeto->Fill( (Ks+pi).M() - massDminus );
        
        mass_Lc2pKs_as_Kspi_beforeVeto->Fill( (Ks+pi_asP_MissID).M() - massLambda_c );
        
        //fill histos with veto cuts
        if(pi1_D_ProbNNpi/(pi1_D_ProbNNpi+pi1_D_ProbNNk) > 0.6 || abs((K_D + pi1_fromD_asK_MissID + pi2_D).M()-massDs) > 25 )mass_D2KKpi_as_D2Kpipi_afterVeto->Fill((K_D + pi1_fromD_asK_MissID + pi2_D).M()-massDs);
        if(pi2_D_ProbNNpi/(pi2_D_ProbNNpi+pi2_D_ProbNNk) > 0.6 || abs((K_D + pi2_fromD_asK_MissID + pi1_D).M()-massDs) > 25 )mass_D2KKpi_as_D2Kpipi_afterVeto->Fill((K_D + pi2_fromD_asK_MissID + pi1_D).M()-massDs);         
        
        if( pi1_D_ProbNNpi/(pi1_D_ProbNNpi+pi1_D_ProbNNp) > 0.75 || abs((K_D + pi1_fromD_asP_MissID + pi2_D).M() - massLambda_c) > 25 )mass_Lc2Kppi_as_D2Kpipi_afterVeto->Fill( (K_D + pi1_fromD_asP_MissID + pi2_D).M() - massLambda_c );
        if( pi2_D_ProbNNpi/(pi2_D_ProbNNpi+pi2_D_ProbNNp) > 0.75 || abs((K_D + pi2_fromD_asP_MissID + pi1_D).M() - massLambda_c) > 25 )mass_Lc2Kppi_as_D2Kpipi_afterVeto->Fill( (K_D + pi2_fromD_asP_MissID + pi1_D).M() - massLambda_c );
        
        if(pi_ProbNNpi>0.8)mass_B2DKsK_as_DKspi_afterVeto->Fill( (Ks+D+pi_asK_MissID).M() - massBs );
        
        if((D+pi).M() < 4900 )mass_B2Dpi_as_DKspi_afterVeto->Fill( (D+pi).M() - massB0 );
        
        if(abs((K_D+pi1_fromD_asK_MissID).M() - massPhi) > 10 )mass_phi_as_Kpi_afterVeto->Fill( (K_D+pi1_fromD_asK_MissID).M() - massPhi );
        if( abs((K_D+pi2_fromD_asK_MissID).M() - massPhi ) > 10 )mass_phi_as_Kpi_afterVeto->Fill( (K_D+pi2_fromD_asK_MissID).M() - massPhi );
    }
    
    TCanvas* c = new TCanvas();
    
    mass_D2KKpi_as_D2Kpipi_beforeVeto->SetMinimum(0);
    mass_D2KKpi_as_D2Kpipi_beforeVeto->SetLineColor(kBlack);
    mass_D2KKpi_as_D2Kpipi_beforeVeto->SetMarkerColor(kBlack);
    mass_D2KKpi_as_D2Kpipi_beforeVeto->Draw("");
    mass_D2KKpi_as_D2Kpipi_afterVeto->SetLineColor(kBlue);
    mass_D2KKpi_as_D2Kpipi_afterVeto->SetMarkerColor(kBlue);
    mass_D2KKpi_as_D2Kpipi_afterVeto->Draw("histsame");
    c->Print("plots/D2KKpi_as_D2Kpipi.eps");
    
    mass_D2pipipi_as_D2Kpipi_beforeVeto->SetMinimum(0);
    mass_D2pipipi_as_D2Kpipi_beforeVeto->SetLineColor(kBlack);
    mass_D2pipipi_as_D2Kpipi_beforeVeto->SetMarkerColor(kBlack);
    mass_D2pipipi_as_D2Kpipi_beforeVeto->Draw("");
    mass_D2pipipi_as_D2Kpipi_afterVeto->SetLineColor(kBlue);
    mass_D2pipipi_as_D2Kpipi_afterVeto->SetMarkerColor(kBlue);
    mass_D2pipipi_as_D2Kpipi_afterVeto->Draw("histsame");
    c->Print("plots/D2pipipi_as_D2Kpipi.eps");

    mass_Lc2ppipi_as_D2Kpipi_beforeVeto->SetMinimum(0);
    mass_Lc2ppipi_as_D2Kpipi_beforeVeto->SetLineColor(kBlack);
    mass_Lc2ppipi_as_D2Kpipi_beforeVeto->SetMarkerColor(kBlack);
    mass_Lc2ppipi_as_D2Kpipi_beforeVeto->Draw("");
    mass_Lc2ppipi_as_D2Kpipi_afterVeto->SetLineColor(kBlue);
    mass_Lc2ppipi_as_D2Kpipi_afterVeto->SetMarkerColor(kBlue);
    mass_Lc2ppipi_as_D2Kpipi_afterVeto->Draw("histsame");
    c->Print("plots/Lc2ppipi_as_D2Kpipi_beforeVeto.eps");

    mass_Lc2Kppi_as_D2Kpipi_beforeVeto->SetMinimum(0);
    mass_Lc2Kppi_as_D2Kpipi_beforeVeto->SetLineColor(kBlack);
    mass_Lc2Kppi_as_D2Kpipi_beforeVeto->SetMarkerColor(kBlack);
    mass_Lc2Kppi_as_D2Kpipi_beforeVeto->Draw("");
    mass_Lc2Kppi_as_D2Kpipi_afterVeto->SetLineColor(kBlue);
    mass_Lc2Kppi_as_D2Kpipi_afterVeto->SetMarkerColor(kBlue);
    mass_Lc2Kppi_as_D2Kpipi_afterVeto->Draw("histsame");
    c->Print("plots/Lc2Kppi_as_D2Kpipi_beforeVeto.eps");
    
    mass_D2Kpi_as_D2Kpipi_beforeVeto->SetMinimum(0);
    mass_D2Kpi_as_D2Kpipi_beforeVeto->SetLineColor(kBlack);
    mass_D2Kpi_as_D2Kpipi_beforeVeto->SetMarkerColor(kBlack);
    mass_D2Kpi_as_D2Kpipi_beforeVeto->Draw("");
    mass_D2Kpi_as_D2Kpipi_afterVeto->SetLineColor(kBlue);
    mass_D2Kpi_as_D2Kpipi_afterVeto->SetMarkerColor(kBlue);
    mass_D2Kpi_as_D2Kpipi_afterVeto->Draw("histsame");
    c->Print("plots/D2Kpi_as_D2Kpipi_beforeVeto.eps");
    

    mass_phi_as_Kpi_beforeVeto->SetMinimum(0);
    mass_phi_as_Kpi_beforeVeto->SetLineColor(kBlack);
    mass_phi_as_Kpi_beforeVeto->SetMarkerColor(kBlack);
    mass_phi_as_Kpi_beforeVeto->Draw("");
    mass_phi_as_Kpi_afterVeto->SetLineColor(kBlue);
    mass_phi_as_Kpi_afterVeto->SetMarkerColor(kBlue);
    mass_phi_as_Kpi_afterVeto->Draw("histsame");
    c->Print("plots/phi_as_Kpi_beforeVeto.eps");

    mass_phi_as_Kpi_beforeVeto2->SetMinimum(0);
    mass_phi_as_Kpi_beforeVeto2->SetLineColor(kBlack);
    mass_phi_as_Kpi_beforeVeto2->SetMarkerColor(kBlack);
    mass_phi_as_Kpi_beforeVeto2->Draw("");
    mass_phi_as_Kpi_afterVeto2->SetLineColor(kBlue);
    mass_phi_as_Kpi_afterVeto2->SetMarkerColor(kBlue);
    mass_phi_as_Kpi_afterVeto2->Draw("histsame");
    c->Print("plots/phi_as_Kpi_beforeVeto2.eps");


    mass_B2DKsK_as_DKspi_beforeVeto->SetMinimum(0);
    mass_B2DKsK_as_DKspi_beforeVeto->SetLineColor(kBlack);
    mass_B2DKsK_as_DKspi_beforeVeto->SetMarkerColor(kBlack);
    mass_B2DKsK_as_DKspi_beforeVeto->Draw("");
    mass_B2DKsK_as_DKspi_afterVeto->SetLineColor(kBlue);
    mass_B2DKsK_as_DKspi_afterVeto->SetMarkerColor(kBlue);
    mass_B2DKsK_as_DKspi_afterVeto->Draw("histsame");
    c->Print("plots/B2DKsK_as_DKspi_beforeVeto.eps");

    mass_B2Dpi_as_DKspi_beforeVeto->SetMinimum(0);
    mass_B2Dpi_as_DKspi_beforeVeto->SetLineColor(kBlack);
    mass_B2Dpi_as_DKspi_beforeVeto->SetMarkerColor(kBlack);
    mass_B2Dpi_as_DKspi_beforeVeto->Draw("");
    mass_B2Dpi_as_DKspi_afterVeto->SetLineColor(kBlue);
    mass_B2Dpi_as_DKspi_afterVeto->SetMarkerColor(kBlue);
    mass_B2Dpi_as_DKspi_afterVeto->Draw("histsame");
    c->Print("plots/B2Dpi_as_DKspi_beforeVeto.eps");

    mass_B2DK_as_DKspi_beforeVeto->SetMinimum(0);
    mass_B2DK_as_DKspi_beforeVeto->SetLineColor(kBlack);
    mass_B2DK_as_DKspi_beforeVeto->SetMarkerColor(kBlack);
    mass_B2DK_as_DKspi_beforeVeto->Draw("");
    mass_B2DK_as_DKspi_afterVeto->SetLineColor(kBlue);
    mass_B2DK_as_DKspi_afterVeto->SetMarkerColor(kBlue);
    mass_B2DK_as_DKspi_afterVeto->Draw("histsame");
    c->Print("plots/B2DK_as_DKspi_beforeVeto.eps");

    mass_D2Kspi_as_Kspi_beforeVeto->SetMinimum(0);
    mass_D2Kspi_as_Kspi_beforeVeto->SetLineColor(kBlack);
    mass_D2Kspi_as_Kspi_beforeVeto->SetMarkerColor(kBlack);
    mass_D2Kspi_as_Kspi_beforeVeto->Draw("");
    mass_D2Kspi_as_Kspi_afterVeto->SetLineColor(kBlue);
    mass_D2Kspi_as_Kspi_afterVeto->SetMarkerColor(kBlue);
    mass_D2Kspi_as_Kspi_afterVeto->Draw("histsame");
    c->Print("plots/D2Kspi_as_Kspi_beforeVeto.eps");

    mass_Lc2pKs_as_Kspi_beforeVeto->SetMinimum(0);
    mass_Lc2pKs_as_Kspi_beforeVeto->SetLineColor(kBlack);
    mass_Lc2pKs_as_Kspi_beforeVeto->SetMarkerColor(kBlack);
    mass_Lc2pKs_as_Kspi_beforeVeto->Draw("");
    mass_Lc2pKs_as_Kspi_afterVeto->SetLineColor(kBlue);
    mass_Lc2pKs_as_Kspi_afterVeto->SetMarkerColor(kBlue);
    mass_Lc2pKs_as_Kspi_afterVeto->Draw("histsame");
    c->Print("plots/Lc2pKs_as_Kspi_beforeVeto.eps");
}


int main(int argc, char** argv){

	time_t startTime = time(0);
	TH1::SetDefaultSumw2();
	TH2::SetDefaultSumw2();
	gROOT->ProcessLine(".x ../lhcbStyle.C");

    BkgStudies();

	std::cout << "==============================================" << std::endl;
	std::cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << std::endl;
	std::cout << "==============================================" << std::endl;

	return 0;
}
