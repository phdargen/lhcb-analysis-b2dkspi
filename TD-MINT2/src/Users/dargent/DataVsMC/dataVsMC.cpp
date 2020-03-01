// MC studies
// author: Philippe d'Argent, Matthieu Kecke
#include <cmath>
#include <algorithm>
#include <iostream>
#include <TChain.h>
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TLorentzVector.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TNtuple.h>
#include "TRandom3.h"
#include <sstream>
#include "TProfile.h"
#include "Mint/NamedParameter.h"
#include "Mint/HyperHistogram.h"
#include "Mint/Utils.h"
#include <fstream>

using namespace std;
using namespace MINT;

/// HFLAV summer 17 values
static const double tau = 1.509;
static const double dgamma = 0.09; 
static const double deltaMs = 17.757;

static const double massKaon = 493.68;
static const double massPion = 139.57;
static const double massBs = 5366.89;
static const double massDs = 1968.30;

void plot(TTree* tree, TTree* treeMC, TString Branch,TString TitleX, int bins, double min, double max, TString weightA, TString weightB, TString newWeightB, TString label, bool log = false, bool legendLeft = false){
        
    /// options
    NamedParameter<string> legTitle("legTitle", (std::string) "");
    NamedParameter<string> nameA("nameA", (std::string) "MC");
    NamedParameter<string> nameB("nameB", (std::string) "Data");
    NamedParameter<string> OutputDir("OutputDir", (std::string) "final/", (char*) 0);
    NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 0);

    cout << "Plotting " << Branch << endl;

    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus(Branch,1);
    if(weightA != "noweight")tree->SetBranchStatus(weightA,1);
    double var;
    float varF[100];
    Short_t varS;
    int varI;
    double sw = 1;
    if(Branch == "BDTG_response" || Branch == "Bs_SS_nnetKaon_PROB" || Branch == "Bs_DTF_MERR")tree->SetBranchAddress(Branch,&varF);
    else if(Branch == "Bs_TAGDECISION_OS" || Branch == "Ds_finalState" || Branch == "TriggerCat" || Branch == "OS_Combination_DEC" || Branch == "SS_Kaon_DEC" || Branch == "run")tree->SetBranchAddress(Branch,&varI);
    else if( A_is_in_B("_DEC",(string)Branch))tree->SetBranchAddress(Branch,&varS);
    else tree->SetBranchAddress(Branch,&var);
    tree->SetBranchAddress(weightA,&sw);
    
    treeMC->SetBranchStatus("*",0);
    treeMC->SetBranchStatus(Branch,1);
    if(weightB != "noweight")treeMC->SetBranchStatus(weightB,1);
    if(newWeightB != "noweight")treeMC->SetBranchStatus(newWeightB,1);
    double varMC;
    float varMCF[100];
    Short_t varMCS;
    int varMCI;
    double w = 1;
    double new_w = 1;
    if(Branch == "BDTG_response" || Branch == "Bs_SS_nnetKaon_PROB" || Branch == "Bs_DTF_MERR")treeMC->SetBranchAddress(Branch,&varMCF);
    else if(Branch == "Bs_TAGDECISION_OS" || Branch == "Ds_finalState" || Branch == "TriggerCat" || Branch == "OS_Combination_DEC" || Branch == "SS_Kaon_DEC" || Branch == "run")treeMC->SetBranchAddress(Branch,&varMCI);
//     else if( A_is_in_B("_DEC",(string)Branch))treeMC->SetBranchAddress(Branch,&varMCS);
    else treeMC->SetBranchAddress(Branch,&varMC);         
    if(weightB != "noweight")treeMC->SetBranchAddress(weightB,&w);
    if(newWeightB != "noweight")treeMC->SetBranchAddress(newWeightB,&new_w);

    ///Make histograms
    TString title= ";"+TitleX+";Yield (a.u.)";
    TH1D h(Branch,title,bins,min,max);
    TH1D h_MC(Branch+"_MC",title,bins,min,max);
    TH1D h_MC_rw(Branch+"_MC_rw",title,bins,min,max);
    
    ///loop over data events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++)
    {	
        if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);
        if(Branch == "BDTG_response" || Branch == "Bs_SS_nnetKaon_PROB" || Branch == "Bs_DTF_MERR")var = (double)varF[0];
        else if(Branch == "Bs_TAGDECISION_OS" || Branch == "Ds_finalState" || Branch == "TriggerCat" || Branch == "run")var = (double)varI;
        else if( A_is_in_B("_DEC",(string)Branch))var = (double)varS;
        h.Fill(var,sw);
    }
    
    ///loop over MC events
    int numEventsMC = treeMC->GetEntries();
    for(int i=0; i< numEventsMC; i++)
    {	
        if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
        treeMC->GetEntry(i);
        if(Branch == "BDTG_response" || Branch == "Bs_SS_nnetKaon_PROB"|| Branch == "Bs_DTF_MERR")varMC = (double)varMCF[0];
        else if(Branch == "Bs_TAGDECISION_OS" || Branch == "Ds_finalState" || Branch == "TriggerCat" || Branch == "run")varMC = (double)varMCI;
        else if( A_is_in_B("_DEC",(string)Branch))varMC = (double)varMCS;
        h_MC.Fill(varMC,w);
        h_MC_rw.Fill(varMC,new_w);
    }
    
    ///Plot it
    TCanvas c;
    
    h.Scale(1./h.Integral());
    h_MC.Scale(1./h_MC.Integral());
    h_MC_rw.Scale(1./h_MC_rw.Integral());
    double maxY= h.GetMaximum();
    if(h_MC.GetMaximum()>maxY)maxY=h_MC.GetMaximum();
    h.SetMinimum(0.);
    if(log){
        h.SetMinimum(0.0001);
        gPad->SetLogy(1);
    }
    else gPad->SetLogy(0);
    h.SetMaximum(maxY*1.4);
    h.SetLineColor(kBlack);
    h.Draw("");
    h_MC.SetMarkerColor(kRed);
    h_MC.SetLineColor(kRed);
    h_MC.Draw("esame");
    h_MC_rw.SetLineColor(kBlue);
    h_MC_rw.SetMarkerColor(kBlue);
    if(newWeightB != weightB && newWeightB != "noweight")h_MC_rw.Draw("esame");
    
    double KolmoTest = h.KolmogorovTest(&h_MC);
    double KolmoTest_rw = h.KolmogorovTest(&h_MC_rw);

    TLegend* leg;
    if(legendLeft)leg = new TLegend(0.15,0.6,0.45,0.9,"");
    else leg = new TLegend(0.55,0.6,0.85,0.9,"");
    leg->SetLineStyle(0);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetTextFont(22);
    leg->SetTextColor(1);
    leg->SetTextSize(0.04);
    leg->SetTextAlign(12);

    if((string)legTitle != "")leg->AddEntry((TObject*)0,((string)legTitle).c_str(), "");
    leg->AddEntry(&h,((string)nameB).c_str(),"LEP");
    
    TLegendEntry* le = leg->AddEntry(&h_MC,((string)nameA).c_str(),"LEP");
    le->SetTextColor(kRed);    

    stringstream ss ;
    TString leg_kol = "KS-Test : ";
    ss << std::fixed << std::setprecision(4) << KolmoTest ;
    leg_kol += ss.str();    
//     le = leg->AddEntry((TObject*)0, leg_kol, "");
    le->SetTextColor(kRed);    

    if(newWeightB != weightB && newWeightB != "noweight"){
        le = leg->AddEntry(&h_MC_rw, ((string)nameA + " (reweighted)").c_str(),"LEP");
        le->SetTextColor(kBlue);  
    }
    ss.str("");
    leg_kol = "KS-Test : ";
    ss << std::fixed << std::setprecision(4) << KolmoTest_rw ;
    leg_kol += ss.str();    
    if(newWeightB != weightB && newWeightB != "noweight"){
        TLegendEntry* le = leg->AddEntry((TObject*)0, leg_kol, "");
        le->SetTextColor(kBlue);    
    }
    leg->Draw(); 
    
    cout << endl;
    if(A_is_in_B("_1.00",(string)Branch))Branch.ReplaceAll("_1.00","");
    c.Print((string)OutputDir + label + "_"+Branch+".eps");
    if(updateAnaNotePlots)c.Print("../../../../../TD-AnaNote/latex/figs/dataVsMC/" + (string)OutputDir + label + "_"+Branch+".pdf" );

}

void compare(TString fileA, TString fileB, TString weightA, TString weightB, TString newWeightB, TString CutA = "", TString CutB = "", int Year = -1, TString finalState = "all", int Trigger = -1, TString label = ""){
    
    // Cuts
    TString Cut;
    if(Year>10)Cut += " year == " + anythingToString(Year);
    else if(Year == -1)Cut += " year > " + anythingToString(Year);
    else Cut += " run == " + anythingToString(Year);
    if(finalState == "KKpi")Cut += " && Ds_finalState < 3 ";
    else if(finalState == "pipipi")Cut += " && Ds_finalState == 3 ";
    else if(finalState == "Kpipi")Cut += " && Ds_finalState == 4 ";
    if(Trigger != -1) Cut += " && TriggerCat == " + anythingToString(Trigger);
    
    if(CutA != "")CutA += " && ";
    CutA += Cut;
    if(CutB != "")CutB += " && ";
    CutB += Cut;
    
    cout << endl << "Comparing file " << endl << fileA << " ( " << CutA << " ) " << endl; 
    cout << " to " << endl << fileB << " ( " << CutB << " ) " << endl << endl;
    
    ///Load files
    TFile* fA = new TFile(fileA,"OPEN");
    TTree* treeA =dynamic_cast<TTree*>(fA->Get("DecayTree"));        

    TFile* fB = new TFile(fileB,"OPEN");
    TTree* treeB =dynamic_cast<TTree*>(fB->Get("DecayTree"));        
    
    treeA->SetBranchStatus("*ENDVERTEX*",0);
    treeA->SetBranchStatus("*OWNPV*",0);
    treeA->SetBranchStatus("*ORIVX*",0);
    treeA->SetBranchStatus("*TOPPV*",0);
    treeA->SetBranchStatus("*TRUE*VERTEX*",0);
    treeA->SetBranchStatus("Bs_*_DEC",0);
    treeA->SetBranchStatus("Bs_*_PROB",0);
    treeA->SetBranchStatus("Bs_B0DTF_*",0);
    treeA->SetBranchStatus("Bs_DTF_*",0);
    treeA->SetBranchStatus("Bs_BsDTF_*",0);
    treeA->SetBranchStatus("Bs_PV_*",0);
    treeA->SetBranchStatus("*BsTaggingTool*",0);
    treeA->SetBranchStatus("*_PP_*",0);
    treeA->SetBranchStatus("*ProtoParticles*",0);
    treeA->SetBranchStatus("*gen*",0);
    treeA->SetBranchStatus("*corr*",0);
    treeA->SetBranchStatus("*CHI2*",1);
    treeA->SetBranchStatus("*TAU*",1);
    treeA->SetBranchStatus("Bs_DTF_MM*",1);
    treeA->SetBranchStatus("*DIRA*",1);

    treeB->SetBranchStatus("*ENDVERTEX*",0);
    treeB->SetBranchStatus("*OWNPV*",0);
    treeB->SetBranchStatus("*ORIVX*",0);
    treeB->SetBranchStatus("*TOPPV*",0);
    treeB->SetBranchStatus("*TRUE*VERTEX*",0);
    treeB->SetBranchStatus("Bs_*_DEC",0);
    treeB->SetBranchStatus("Bs_*_PROB",0);
    treeB->SetBranchStatus("Bs_B0DTF_*",0);
    treeB->SetBranchStatus("Bs_DTF_*",0);
    treeB->SetBranchStatus("Bs_BsDTF_*",0);
    treeB->SetBranchStatus("Bs_PV_*",0);
    treeB->SetBranchStatus("*BsTaggingTool*",0);
    treeB->SetBranchStatus("*_PP_*",0);
    treeB->SetBranchStatus("*ProtoParticles*",0);
    treeB->SetBranchStatus("*gen*",0);
    treeB->SetBranchStatus("*corr*",0);
    treeB->SetBranchStatus("*CHI2*",1);
    treeB->SetBranchStatus("*TAU*",1);
    treeB->SetBranchStatus("Bs_DTF_MM*",1);
    treeB->SetBranchStatus("*DIRA*",1);

    TFile* output = new TFile("dummy.root","RECREATE");
    TTree* new_treeA = treeA->CopyTree(CutA);
    TTree* new_treeB = treeB->CopyTree(CutB);

    /// Options
    NamedParameter<int> nBins("nBins", 40); 
    
    label += "Ds2";
    if(finalState != "")label +=  finalState;
    else label += "all";
    if(Year>-1) label += "_" + anythingToString(Year);
    if(Trigger>-1) label += "_t" + anythingToString(Trigger);
    
    TString Decay, selection;
    if(A_is_in_B("signal",(string)fileA) && A_is_in_B("signal",(string)fileB)) Decay = "signal";
    if(A_is_in_B("norm",(string)fileA) && A_is_in_B("norm",(string)fileB)) Decay = "norm";
    if(A_is_in_B("Final",(string)fileA) && A_is_in_B("Final",(string)fileB)) selection = "Final";

    /// Bs
    plot(new_treeA,new_treeB,"Bs_PT","p_{T}(B) [MeV]",nBins,0,40000,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_P","p(B) [MeV]",nBins,0,900000,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_ETA","#eta(B)",nBins,1,6,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_FDCHI2_OWNPV","#chi^{2}_{FD}(B)",nBins,0,100000,weightA, weightB, newWeightB, label,true);     plot(new_treeA,new_treeB,"Bs_ENDVERTEX_CHI2","#chi^{2}_{vtx}(B)",nBins,0,35,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_BsDTF_TAU","t(B) [ps]",nBins,0.4,10.,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_BsDTF_TAUERR","#sigma_{t}(B) [ps]",nBins,0,0.1,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_DTF_TAU","t(B) [ps]",nBins,0.4,10.,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_DTF_TAUERR","#sigma_{t}(B) [ps]",nBins,0,0.1,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_ptasy_1.00","B_ptasy_1.00",nBins, -1, 1 ,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"Bs_DTF_MMERR","#sigma_{m} [MeV]",nBins,4.,25.,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Bs_DTF_MM","m(B) [MeV]",100.,5200.,5700.,"noweight", "noweight", "noweight", label);

    /// BDT
//     plot(new_treeA,new_treeB,"PV_CHI2NDOF","DTF #chi^{2}",nBins,0.,7,weightA, weightB, newWeightB, label);    
    plot(new_treeA,new_treeB,"Bs_IPCHI2_OWNPV","#chi^{2}_{IP}(B)",nBins,0,20,weightA, weightB, newWeightB, label,true);
    plot(new_treeA,new_treeB,"Bs_DIRA_OWNPV","DIRA(B)",nBins,0.99997,1,weightA, weightB, newWeightB, label,true,true);

//     plot(new_treeA,new_treeB,"XsDaughters_min_IPCHI2","X_{s} min(#chi^{2}_{IP})",nBins, 0, 10000 ,weightA, weightB, newWeightB, label,true);
//     if(Decay == "norm")plot(new_treeA,new_treeB,"a_1_1260_plus_ptasy_1.00","Xs_ptasy_1.00",nBins, -1, 1. ,weightA, weightB, newWeightB, label,false,true);
//     else plot(new_treeA,new_treeB,"K_1_1270_plus_ptasy_1.00","Xs_ptasy_1.00",nBins, -1, 1 ,weightA, weightB, newWeightB, label,false,true);
//     plot(new_treeA,new_treeB,"Xs_max_DOCA","X_{s} max DOCA [mm]",nBins, 0, 0.4 ,weightA, weightB, newWeightB, label);

     plot(new_treeA,new_treeB,"track_min_IPCHI2","min(#chi^{2}_{IP})",nBins, 0, 10000 ,weightA, weightB, newWeightB, label,true);
//     plot(new_treeA,new_treeB,"DsDaughters_min_IPCHI2","D_{s} min(#chi^{2}_{IP})",nBins, 0, 10000 ,weightA, weightB, newWeightB, label,true);
//     plot(new_treeA,new_treeB,"Ds_ptasy_1.00","Ds_ptasy_1.00",nBins, -1, 1 ,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"Ds_FDCHI2_ORIVX","#chi^{2}_{FD}(D_{s})",nBins,0,40000,weightA, weightB, newWeightB, label,true);
//     plot(new_treeA,new_treeB,"Ds_RFD","Ds RFD",nBins,0,10,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Ds_PT","p_{T}(Ds) [MeV]",nBins,0,40000,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Ds_ETA","#eta(Ds)",nBins,1,6,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Ds_PV_TAU","t(Ds) [ps]",nBins,0.,5.,weightA, weightB, newWeightB, label);

    plot(new_treeA,new_treeB,"pi_PT","p_{T}(pi) [MeV]",nBins,0,40000,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"pi_ETA","#eta(pi)",nBins,1,6,weightA, weightB, newWeightB, label);


//     plot(new_treeA,new_treeB,"maxCos","maxCos",nBins,-1,1,weightA, weightB, newWeightB, label);    
//     plot(new_treeA,new_treeB,"max_ghostProb","max(Track_ghostProb)",nBins,0,0.4,weightA, weightB, newWeightB, label);
//     plot(new_treeA,new_treeB,"max_ProbNNghost","max(Track_ghostProb)",nBins,0,0.4,weightA, weightB, newWeightB, label);
//     plot(new_treeA,new_treeB,"track_min_PT","min(p_{T})  [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);

    /// Tagging
    plot(new_treeA,new_treeB,"OS_Combination_DEC","q_{OS}",8,-1.5,6.5,weightA, weightB, newWeightB, label);    
    plot(new_treeA,new_treeB,"OS_Combination_PROB","#eta_{OS}",nBins,0.,0.499999,weightA, weightB, newWeightB, label, false, true);    
    plot(new_treeA,new_treeB,"SS_Kaon_DEC","q_{SS}",8,-1.5,6.5,weightA, weightB, newWeightB, label);    
    plot(new_treeA,new_treeB,"SS_Kaon_PROB","#eta_{SS}",nBins,0.,0.499999,weightA, weightB, newWeightB, label, false, true);  
    
    /// Ds
//     plot(new_treeA,new_treeB,"Ds_PV_MM","Ds_PV_MM",nBins,1950,1990,weightA, weightB, newWeightB, label);
//     plot(new_treeA,new_treeB,"Ds_m12","Ds_m12",nBins,900,1900,weightA, weightB, newWeightB, label);
//     plot(new_treeA,new_treeB,"Ds_m13","Ds_m13",nBins,600,1600,weightA, weightB, newWeightB, label);
//     plot(new_treeA,new_treeB,"Ds_PT","p_{T}(D_{s}) [MeV]",nBins,0,40000,weightA, weightB, newWeightB, label);
//     plot(new_treeA,new_treeB,"Ds_ETA","#eta(D_{s})",nBins,1,6,weightA, weightB, newWeightB, label);

    /// Trigger
    plot(new_treeA,new_treeB,"TriggerCat","Trigger category",8,-0.5,7.5,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"run","Run",8,-0.5,7.5,weightA, weightB, newWeightB, label);
    plot(new_treeA,new_treeB,"Ds_finalState","D_{s} final state",8,-0.5,7.5,weightA, weightB, newWeightB, label);
    
    /// Dalitz
    plot(new_treeA,new_treeB,"m_Kpipi","m(K^{+}#pi^{+}#pi^{-})[MeV]",nBins,900,1900,weightA, weightB, newWeightB, label,false);
    plot(new_treeA,new_treeB,"m_Kpi","m(K^{+}#pi^{-})[MeV]",nBins,600,1200,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"m_pipi","m(#pi^{+}#pi^{-})[MeV]",nBins,200,1200,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"m_Dspipi","m(D_{s}^{-}#pi^{+}#pi^{-})[MeV]",nBins,1900,5550,weightA, weightB, newWeightB, label,false,true);
    plot(new_treeA,new_treeB,"m_Dspi","m(D_{s}^{-}#pi^{+})[MeV]",nBins,0,5500,weightA, weightB, newWeightB, label,false,true);
//     plot(new_treeA,new_treeB,"m_DsK","m(D_{s}^{-}K^{+})[MeV]",nBins,0,5500,weightA, weightB, newWeightB, label,false,true);
//     plot(new_treeA,new_treeB,"m_DsKpi","m(D_{s}^{-}K^{+}#pi^{-})[MeV]",nBins,1900,5500,weightA, weightB, newWeightB, label,false,true);
//     plot(new_treeA,new_treeB,"m_DsKpip","m(D_{s}^{-}K^{+}#pi^{+})[MeV]",nBins,1900,5500,weightA, weightB, newWeightB, label,false,true);

    if(Decay == "signal"){
//         plot(new_treeA,new_treeB,"K_plus_PT","p_{T}(K^{+}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"K_plus_ETA","#eta(K^{+})",nBins,1,6,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"K_plus_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
//         plot(new_treeA,new_treeB,"K_plus_PIDK","DLL_{K#pi}(K^{+}) ",nBins,-20,100,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"K_plus_TRACK_GhostProb","ghost prob (K^{+})",nBins,0,0.4,weightA, weightB, newWeightB, label);
        //plot(new_treeA,new_treeB,"K_plus_ProbNNghost","ghost prob (K^{+})",nBins,0,0.4,weightA, weightB, newWeightB, label);

//         plot(new_treeA,new_treeB,"pi_plus_PT","p_{T}(#pi^{+}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"pi_plus_ETA","#eta(#pi^{+})",nBins,1,6,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"pi_plus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{+})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
//         plot(new_treeA,new_treeB,"pi_plus_PIDK","DLL_{K#pi}(#pi^{+}) ",nBins,-100,100,weightA, weightB, newWeightB, label);
        //plot(new_treeA,new_treeB,"pi_plus_TRACK_GhostProb","ghost prob (#pi^{+})",nBins,0,0.4,weightA, weightB, newWeightB, label);

//         plot(new_treeA,new_treeB,"pi_minus_PT","p_{T}(#pi^{-}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"pi_minus_ETA","#eta(#pi^{-})",nBins,1,6,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"pi_minus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
//         plot(new_treeA,new_treeB,"pi_minus_PIDK","DLL_{K#pi}(#pi^{-}) ",nBins,-100,100,weightA, weightB, newWeightB, label);
        //plot(new_treeA,new_treeB,"pi_minus_TRACK_GhostProb","ghost prob (#pi^{-})",nBins,0,0.4,weightA, weightB, newWeightB, label);
    }   
    else if(Decay == "norm") {
//         plot(new_treeA,new_treeB,"pi_plus1_PT","p_{T}(K^{+}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"pi_plus1_ETA","#eta(K^{+})",nBins,1,6,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"pi_plus1_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
//         plot(new_treeA,new_treeB,"pi_plus1_PIDK","DLL_{K#pi}(#pi^{+}) ",nBins,-100,100,weightA, weightB, newWeightB, label);
       // plot(new_treeA,new_treeB,"pi_plus1_TRACK_GhostProb","ghost prob (K^{+})",nBins,0,0.4,weightA, weightB, newWeightB, label);

//         plot(new_treeA,new_treeB,"pi_plus2_PT","p_{T}(#pi^{+}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"pi_plus2_ETA","#eta(#pi^{+})",nBins,1,6,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"pi_plus2_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{+})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
//         plot(new_treeA,new_treeB,"pi_plus2_PIDK","DLL_{K#pi}(#pi^{+}) ",nBins,-100,100,weightA, weightB, newWeightB, label);
        //plot(new_treeA,new_treeB,"pi_plus2_TRACK_GhostProb","ghost prob (#pi^{+})",nBins,0,0.4,weightA, weightB, newWeightB, label);

//         plot(new_treeA,new_treeB,"pi_minus_PT","p_{T}(#pi^{-}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"pi_minus_ETA","#eta(#pi^{-})",nBins,1,6,weightA, weightB, newWeightB, label);
//         plot(new_treeA,new_treeB,"pi_minus_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
//         plot(new_treeA,new_treeB,"pi_minus_PIDK","DLL_{K#pi}(#pi^{-}) ",nBins,-100,100,weightA, weightB, newWeightB, label);
        //plot(new_treeA,new_treeB,"pi_minus_TRACK_GhostProb","ghost prob (#pi^{-})",nBins,0,0.4,weightA, weightB, newWeightB, label);
    }
    if(finalState == "KKpi") {
        plot(new_treeA,new_treeB,"K_plus_fromDs_PT","p_{T}(K^{+} from D_{s}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_plus_fromDs_ETA","#eta(K^{+} from D_{s})",nBins,1,6,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_plus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{+} from D_{s})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
        //plot(new_treeA,new_treeB,"K_plus_fromDs_TRACK_GhostProb","ghost prob (K^{+} from D_{s})",nBins,0,0.4,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_plus_fromDs_PIDK","DLL_{K#pi}(K^{+}) ",nBins,-20,100,weightA, weightB, newWeightB, label);

        plot(new_treeA,new_treeB,"pi_minus_fromDs_PT","p_{T}(#pi^{-} from D_{s}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_minus_fromDs_ETA","#eta(#pi^{-} from D_{s})",nBins,1,6,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_minus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(#pi^{-} from D_{s})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
        //plot(new_treeA,new_treeB,"pi_minus_fromDs_TRACK_GhostProb","ghost prob (#pi^{-} from D_{s})",nBins,0,0.4,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"pi_minus_fromDs_PIDK","DLL_{K#pi}(#pi^{-} from D_{s}) ",nBins,-100,100,weightA, weightB, newWeightB, label);

        plot(new_treeA,new_treeB,"K_minus_fromDs_PT","p_{T}(K^{-} from D_{s}) [MeV]",nBins,0,10000,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_minus_fromDs_ETA","#eta(K^{-} from D_{s})",nBins,1,6,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_minus_fromDs_IPCHI2_OWNPV","#chi^{2}_{IP}(K^{-} from D_{s})",nBins,0,10000,weightA, weightB, newWeightB, label,true);
        //plot(new_treeA,new_treeB,"K_minus_fromDs_TRACK_GhostProb","ghost prob (K^{-} from D_{s})",nBins,0,0.4,weightA, weightB, newWeightB, label);
        plot(new_treeA,new_treeB,"K_minus_fromDs_PIDK","DLL_{K#pi}(K^{-} from D_{s}) ",nBins,-20,100,weightA, weightB, newWeightB, label);
    }    
    plot(new_treeA,new_treeB,"NTracks","# of tracks",nBins,0,550,weightA, weightB, newWeightB, label);
    //if(selection == "Final") 
    plot(new_treeA,new_treeB,"BDTG_response","BDTG",nBins,0,1.,weightA, weightB, newWeightB, label,false,true);
}

void plotPID(TTree* tree, TTree* treeMC, TTree* treeMC_PIDGen, TTree* treeMC_PIDCorr, TTree* treeMC_noPID, 
 TString Branch,TString TitleX, int bins, double min, double max, TString weightData, TString weightMC, TString label, bool log = false, bool legendLeft = false, bool eff = false){
        
    /// options
    NamedParameter<string> legTitle("legTitle", (std::string) "");
    NamedParameter<string> OutputDir("OutputDir", (std::string) "final/", (char*) 0);
    NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 0);

    cout << "Plotting " << Branch << endl;

    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus(Branch,1);
    tree->SetBranchStatus(weightData,1);
    double var;
    double sw;
    tree->SetBranchAddress(Branch,&var);
    tree->SetBranchAddress(weightData,&sw);
    
    treeMC->SetBranchStatus("*",0);
    treeMC->SetBranchStatus(Branch,1);
    if(weightMC != "noweight")treeMC->SetBranchStatus(weightMC,1);
    double varMC;
    double w = 1;
    treeMC->SetBranchAddress(Branch,&varMC);         
    if(weightMC != "noweight")treeMC->SetBranchAddress(weightMC,&w);

    treeMC_PIDGen->SetBranchStatus("*",0);
    treeMC_PIDGen->SetBranchStatus(Branch,1);
    if(weightMC != "noweight")treeMC_PIDGen->SetBranchStatus(weightMC,1);
    treeMC_PIDGen->SetBranchAddress(Branch,&varMC);         
    if(weightMC != "noweight")treeMC_PIDGen->SetBranchAddress(weightMC,&w);

    treeMC_PIDCorr->SetBranchStatus("*",0);
    treeMC_PIDCorr->SetBranchStatus(Branch,1);
    if(weightMC != "noweight")treeMC_PIDCorr->SetBranchStatus(weightMC,1);
    treeMC_PIDCorr->SetBranchAddress(Branch,&varMC);         
    if(weightMC != "noweight")treeMC_PIDCorr->SetBranchAddress(weightMC,&w);

    treeMC_noPID->SetBranchStatus("*",0);
    treeMC_noPID->SetBranchStatus(Branch,1);
    if(weightMC != "noweight")treeMC_noPID->SetBranchStatus(weightMC,1);
    treeMC_noPID->SetBranchAddress(Branch,&varMC);         
    if(weightMC != "noweight")treeMC_noPID->SetBranchAddress(weightMC,&w);

    ///Make histograms
    TString title;
    if(eff)title = ";"+TitleX+";#epsilon_{PID} [norm.]";
    else title = ";"+TitleX+";Yield [norm.]";
    TH1D h(Branch,title,bins,min,max);
    TH1D h_MC(Branch+"_MC",title,bins,min,max);
    TH1D h_MC_PIDGen(Branch+"_MC_PIDGen",title,bins,min,max);
    TH1D h_MC_PIDCorr(Branch+"_MC_PIDCorr",title,bins,min,max);
    TH1D h_MC_noPID(Branch+"_MC_noPID",title,bins,min,max);
    
    ///loop over data events
    int numEvents = tree->GetEntries();
    for(int i=0; i< numEvents; i++)
    {	
        if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEvents << endl;
        tree->GetEntry(i);
        h.Fill(var,sw);
    }
    
    ///loop over MC events
    int numEventsMC = treeMC->GetEntries();
    for(int i=0; i< numEventsMC; i++)
    {	
        if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC << endl;
        treeMC->GetEntry(i);
        h_MC.Fill(varMC,w);
    }

    int numEventsMC_PIDGen = treeMC_PIDGen->GetEntries();
    for(int i=0; i< numEventsMC_PIDGen; i++)
    {	
        if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC_PIDGen << endl;
        treeMC_PIDGen->GetEntry(i);
        h_MC_PIDGen.Fill(varMC,w);
    }

    int numEventsMC_PIDCorr = treeMC_PIDCorr->GetEntries();
    for(int i=0; i< numEventsMC_PIDCorr; i++)
    {	
        if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC_PIDCorr << endl;
        treeMC_PIDCorr->GetEntry(i);
        h_MC_PIDCorr.Fill(varMC,w);
    }

    int numEventsMC_noPID = treeMC_noPID->GetEntries();
    for(int i=0; i< numEventsMC_noPID; i++)
    {	
        if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEventsMC_noPID << endl;
        treeMC_noPID->GetEntry(i);
        h_MC_noPID.Fill(varMC,w);
    }
        
    /// Print efficiencies
    cout << endl << "PID efficiencies: " << endl;
    cout << " MC = " << h_MC.Integral()/h_MC_noPID.Integral() * 100. << " ( % ) " <<  endl;
    cout << " MC (PIDGen) = " << h_MC_PIDGen.Integral()/h_MC_noPID.Integral() * 100. << " ( % ) " <<  endl;
    cout << " MC (PIDCorr) = " << h_MC_PIDCorr.Integral()/h_MC_noPID.Integral() * 100. << " ( % ) " <<  endl << endl;

    ///Plot it
    TCanvas c;
    
    h.Scale(1./h.Integral());
    h_MC.Scale(1./h_MC.Integral());
    h_MC_PIDGen.Scale(1./h_MC_PIDGen.Integral());
    h_MC_PIDCorr.Scale(1./h_MC_PIDCorr.Integral());
    h_MC_noPID.Scale(1./h_MC_noPID.Integral());

    h_MC.SetMarkerColor(kRed);
    h_MC.SetLineColor(kRed);
    h_MC_PIDGen.SetMarkerColor(kBlue);
    h_MC_PIDGen.SetLineColor(kBlue);
    h_MC_PIDCorr.SetMarkerColor(kMagenta+1);
    h_MC_PIDCorr.SetLineColor(kMagenta+1);
	
    if(eff){

	h_MC.Divide(&h_MC,&h_MC_noPID);
	h_MC.SetMinimum(0.);
	h_MC.SetMaximum(2.5);
	h_MC.Draw("e");
	
	h_MC_PIDGen.Divide(&h_MC_PIDGen,&h_MC_noPID);	
	h_MC_PIDGen.SetMinimum(0.);
	h_MC_PIDGen.SetMaximum(2.5);
	h_MC_PIDGen.Draw("esame");
	
	h_MC_PIDCorr.Divide(&h_MC_PIDCorr,&h_MC_noPID);	
	h_MC_PIDCorr.SetMinimum(0.);
	h_MC_PIDCorr.SetMaximum(2.5);
	h_MC_PIDCorr.Draw("esame");
    }
    else {
	
	double maxY= h.GetMaximum();
	if(h_MC.GetMaximum()>maxY)maxY=h_MC.GetMaximum();
	h.SetMinimum(0.);
	if(log){
		h.SetMinimum(0.0001);
		gPad->SetLogy(1);
	}
	else gPad->SetLogy(0);
	
	h.SetMaximum(maxY*1.4);
	h.SetLineColor(kBlack);
	h.Draw("");
	h_MC.Draw("esame");
	h_MC_PIDGen.Draw("esame");
	h_MC_PIDCorr.Draw("esame");
    }
    TLegend* leg;
    if(legendLeft)leg = new TLegend(0.15,0.6,0.45,0.9,"");
    else leg = new TLegend(0.55,0.6,0.85,0.9,"");
    leg->SetLineStyle(0);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetTextFont(22);
    leg->SetTextColor(1);
    leg->SetTextSize(0.04);
    leg->SetTextAlign(12);

    if((string)legTitle != "")leg->AddEntry((TObject*)0,((string)legTitle).c_str(), "");
    if(!eff)leg->AddEntry(&h,"Data","LEP");
	
    TLegendEntry* le = leg->AddEntry(&h_MC,"MC","LEP");
    le->SetTextColor(kRed);    
	
    le = leg->AddEntry(&h_MC_PIDGen,"MC (PIDGen)","LEP");
    le->SetTextColor(kBlue);    
	
    le = leg->AddEntry(&h_MC_PIDCorr,"MC (PIDCorr)","LEP");
    le->SetTextColor(kMagenta+1);    
	
    leg->Draw(); 
	
    if(eff)label = "eff_" + label;
    if(updateAnaNotePlots)c.Print("../../../../../TD-AnaNote/latex/figs/dataVsMC/" + (string)OutputDir + label + "_"+Branch+".pdf" );
    c.Print((string)OutputDir + label + "_"+Branch+".eps"); 
  
}

void comparePID(TString fileA, TString fileB, TString weightA, TString weightB, TString CutA = "", TString CutB = "", int Year = -1, TString finalState = "all", int Trigger = -1, TString label = ""){
    
    // Cuts
    TString Cut;
    if(Year>10)Cut += " year == " + anythingToString(Year);
    else if(Year == -1)Cut += " year > " + anythingToString(Year);
    else Cut += " run == " + anythingToString(Year);
    if(finalState == "KKpi")Cut += " && Ds_finalState < 3 ";
    else if(finalState == "pipipi")Cut += " && Ds_finalState == 3 ";
    else if(finalState == "Kpipi")Cut += " && Ds_finalState == 4 ";
    if(Trigger != -1) Cut += " && TriggerCat == " + anythingToString(Trigger);
    
    if(CutA != "")CutA += " && ";
    CutA += Cut;
    if(CutB != "")CutB += " && ";
    CutB += Cut;
   
    TString fileC = fileB;
    fileC.ReplaceAll("PIDMC","PIDGen"); 
    TString fileD = fileB;
    fileD.ReplaceAll("_PIDMC",""); 
    TString fileE = fileB;
    fileE.ReplaceAll("PIDMC","noPID"); 

    cout << endl << "Comparing file " << endl << fileA << " ( " << CutA << " ) " << endl; 
    cout << " to " << endl << fileB << " ( " << CutB << " ) " << endl ;
    cout << " to " << endl << fileC << " ( " << CutB << " ) " << endl ;
    cout << " to " << endl << fileD << " ( " << CutB << " ) " << endl << endl;
    cout << " normalize to " << endl << fileE << " ( " << CutB << " ) " << endl << endl;

    ///Load files
    TChain* treeA = new TChain("DecayTree");
    treeA->Add(fileA);
    
    TChain* treeB= new TChain("DecayTree");
    treeB->Add(fileB);

    TChain* treeC= new TChain("DecayTree");
    treeC->Add(fileC);

    TChain* treeD= new TChain("DecayTree");
    treeD->Add(fileD);

    TChain* treeE= new TChain("DecayTree");
    treeE->Add(fileE);
    
    TFile* output = new TFile("dummy.root","RECREATE");
    TTree* new_treeA = treeA->CopyTree(CutA);
    TTree* new_treeB = treeB->CopyTree(CutB);
    TTree* new_treeC = treeC->CopyTree(CutB);
    TTree* new_treeD = treeD->CopyTree(CutB);
    TTree* new_treeE = treeE->CopyTree(CutB);

    /// Options
    NamedParameter<int> nBins("nBins", 40); 
    
    label += "PID_Ds2";
    if(finalState != "")label +=  finalState;
    else label += "all";
    if(Year>-1) label += "_" + anythingToString(Year);
    if(Trigger>-1) label += "_t" + anythingToString(Trigger);
    
    TString Decay, selection;
    if(A_is_in_B("signal",(string)fileA) && A_is_in_B("signal",(string)fileB)) Decay = "signal";
    if(A_is_in_B("norm",(string)fileA) && A_is_in_B("norm",(string)fileB)) Decay = "norm";
    if(A_is_in_B("Final",(string)fileA) && A_is_in_B("Final",(string)fileB)) selection = "Final";

    /// Bs
    plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"Bs_PT","p_{T}(B) [MeV]",nBins,0,40000,weightA, weightB, label);
    plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"Bs_ETA","#eta(B)",nBins,1,6,weightA, weightB, label);
    plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"Bs_DTF_TAU","t(B) [ps]",nBins,0.,10.,weightA, weightB, label,false,true,true);
    plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"Bs_DTF_TAUERR","#sigma_{t}(B) [ps]",nBins,0,0.15,weightA, weightB, label);

    /// Dalitz
    plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"m_Kpipi","m(K^{+}#pi^{+}#pi^{-})[MeV]",nBins,1000,1950,weightA, weightB, label,false,true,true);
    plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"m_Kpi","m(K^{+}#pi^{-})[MeV]",nBins,massKaon+massPion,1200,weightA, weightB, label,false,true,true);
    plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"m_pipi","m(#pi^{+}#pi^{-})[MeV]",nBins,2.*massPion,1200,weightA, weightB, label,false,true,true);
    plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"m_Dspipi","m(D_{s}^{-}#pi^{+}#pi^{-})[MeV]",nBins,2400,massBs-massKaon,weightA, weightB, label,false,true,true);
    plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"m_Dspi","m(D_{s}^{-}#pi^{+})[MeV]",nBins,massDs+massPion,massBs-massKaon-massPion,weightA, weightB, label,false,true,true);
    plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"m_DsK","m(D_{s}^{-}K^{+})[MeV]",nBins,0,5500,weightA, weightB, label,false,true,true);
    plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"m_DsKpi","m(D_{s}^{-}K^{+}#pi^{-})[MeV]",nBins,1900,5500,weightA, weightB, label,false,true,true);
    plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"m_DsKpip","m(D_{s}^{-}K^{+}#pi^{+})[MeV]",nBins,1900,5500,weightA, weightB, label,false,true,true);
    
    if(Decay == "signal"){
        plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"K_plus_PIDK","DLL_{K#pi}(K^{+}) ",nBins,10,100,weightA, weightB, label,false);
        plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"pi_plus_PIDK","DLL_{K#pi}(#pi^{+}) ",nBins,-100,20,weightA, weightB, label,false,true);
        plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"pi_minus_PIDK","DLL_{K#pi}(#pi^{-}) ",nBins,-100,20,weightA, weightB, label,false,true);
    }   
    else if(Decay == "norm") {
        plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"pi_plus1_PIDK","DLL_{K#pi}(#pi^{+}) ",nBins,-100,20,weightA, weightB, label,false,true);
        plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"pi_plus2_PIDK","DLL_{K#pi}(#pi^{+}) ",nBins,-100,20,weightA, weightB, label,false,true);
        plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"pi_minus_PIDK","DLL_{K#pi}(#pi^{-}) ",nBins,-100,20,weightA, weightB, label,false,true);
    }
    if(finalState == "KKpi") {
        plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"K_plus_fromDs_PIDK","DLL_{K#pi}(K^{+}) ",nBins,-10,100,weightA, weightB, label,false);
        plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"pi_minus_fromDs_PIDK","DLL_{K#pi}(#pi^{-} from D_{s}) ",nBins,-100,20,weightA, weightB, label,false,true);
	plotPID(new_treeA,new_treeB,new_treeC,new_treeD,new_treeE,"K_minus_fromDs_PIDK","DLL_{K#pi}(K^{-} from D_{s}) ",nBins,-10,100,weightA, weightB, label,false);
    }    
}


void applyCorrectionHisto(vector<TString> vars, int Year, TString FinalState, int Trigger, TString ApplyTo, TString& weightVar, TString NewWeightVar){

    NamedParameter<double> maxWeight("maxWeight",100.);
    NamedParameter<string> OutputDir("OutputDir", (std::string) "final/", (char*) 0);
    const int dim = vars.size();

    double sumw = 0;
    double sumw2 = 0;
    double sumw_old = 0;
    double sumw2_old = 0;
    
    TString label;
    for(int j = 0; j < dim; j++) label += "_" + vars[j] ;
    label += "_Ds2" + FinalState;
    if(Year != -1)label += "_" + anythingToString(Year);
    if(Trigger != -1)label+= "_t" + anythingToString(Trigger);
    
    HyperHistogram hist_weights((string)OutputDir + "weights/weights" + label + ".root");

    TFile* f = new TFile(ApplyTo,"UPDATE");
    TTree* treeMC =dynamic_cast<TTree*>(f->Get("DecayTree"));

    cout << "Apply correction histo " << (string)OutputDir + "weights/weights" + label + ".root" << endl;
    cout << "to file " << ApplyTo << endl;
    
    vector<double> var_MC(dim,0.);
    double weight = 1;
    int Ds_finalState_MC, year_MC, trigger_MC;
    for(int j = 0; j < dim; j++)treeMC->SetBranchAddress(vars[j],&var_MC[j]);
    if(weightVar != "noweight")treeMC->SetBranchAddress(weightVar,&weight);
    treeMC->SetBranchAddress("Ds_finalState",&Ds_finalState_MC);
    if(Year>10)treeMC->SetBranchAddress("year",&year_MC);
    else treeMC->SetBranchAddress("run",&year_MC);
    treeMC->SetBranchAddress("TriggerCat",&trigger_MC);

    vector<double> old_weights;
    for (int i = 0; i < treeMC->GetEntries(); i++){
        treeMC->GetEntry(i);
    	old_weights.push_back(weight);
        sumw_old += weight;
        sumw2_old += weight*weight;
    }

    TBranch* br = (TBranch*)treeMC->GetListOfBranches()->FindObject(NewWeightVar);
    if(br != 0)treeMC->SetBranchStatus(NewWeightVar,0);
    if(NewWeightVar != weightVar)weightVar = NewWeightVar;
    /// Removed for performance reasons: Now existing branch cannot be overwritten ! 
    //TTree* summary_tree = treeMC->CloneTree();
    TBranch* b_w = treeMC->Branch(NewWeightVar,&weight,NewWeightVar+"/D"); 

    HyperPointSet points_MC( dim );

    for (int i = 0; i < treeMC->GetEntries(); i++){
        treeMC->GetEntry(i);

        if( (Year != year_MC && Year != -1) 
           || (FinalState == "KKpi" && Ds_finalState_MC > 2) || (FinalState == "pipipi" && Ds_finalState_MC > 3) 
           || (Trigger != trigger_MC && Trigger != -1) ) {
            weight = old_weights[i];
            b_w->Fill();
            sumw += weight;
            sumw2 += weight*weight;
            continue;
        }
        else if(FinalState == "Kpipi" && Ds_finalState_MC > 4) throw "undefined final state";
        
        HyperPoint point( dim );
        for(int j = 0; j < dim; j++)point.at(j)= var_MC[j]; 

        double w = 1.;
        int bin = hist_weights.getBinning().getBinNum(point);
        if(hist_weights.checkBinNumber(bin)!= bin){
            w = 1; //? should't happen
            cout << "ERROR:: Event outside limits" << endl;
        }else w = std::min((double)maxWeight,hist_weights.getBinContent(bin));
	
        if(w < 0) {
            w = 0.;
            cout << "ERROR:: Negative weight" << endl;
        }
        
        weight = w * old_weights[i];
        b_w->Fill();
        sumw += weight;
        sumw2 += weight*weight;
    }

   cout << "Effective weight before reweighting = " << sumw_old/sumw2_old << endl; 
   cout << "Effective weight after reweighting = " << sumw/sumw2 << endl; 

   treeMC->Write();
   f->Close();

   return;
}

void produceCorrectionHisto(vector<TString> vars, vector<double> min, vector<double> max, int Year, TString FinalState, int Trigger, TString ReweightFromA, TString weightVarA, TString ReweightToB, TString weightVarB){

    /// Options
    NamedParameter<string> OutputDir("OutputDir", (std::string) "final/", (char*) 0);
    NamedParameter<int> minEventsPerBin("minEventsPerBin", 60); 
    NamedParameter<int> maxBinsPerDim("maxBinsPerDim", 200); 
    NamedParameter<int> updateAnaNoteHistos("updateAnaNoteHistos", 0);

    /// Check vector sizes
    if(vars.size() * min.size() * max.size() != pow(vars.size(),3) ){
        cout << "ERROR:: Different number of vars and limits";
        throw "ERROR";
    }

    /// Get dimension and minimum bin width
    const int dim = vars.size();
    vector<double> vec_minBinWidths(dim,0.);
    for(int i = 0; i < dim; i++)vec_minBinWidths[i]= (max[i]-min[i])/(double)maxBinsPerDim;
    HyperPoint minBinWidths(vec_minBinWidths);    

    HyperPoint Min(min);
    HyperPoint Max(max);
    HyperCuboid limits(Min, Max);

    /// Get data           
    TFile* f = new TFile(ReweightToB,"OPEN");
    TTree* tree =dynamic_cast<TTree*>(f->Get("DecayTree"));    

    tree->SetBranchStatus("*",0);
    for(int i = 0; i < dim; i++)tree->SetBranchStatus(vars[i],1);
    if((string)weightVarB != "noweight")tree->SetBranchStatus(weightVarB,1);
    tree->SetBranchStatus("Ds_finalState",1);
    tree->SetBranchStatus("year",1);
    tree->SetBranchStatus("run",1);
    tree->SetBranchStatus("TriggerCat",1);

    vector<double> var_data(dim,0.);
    double sw = 1;
    int Ds_finalState, year, trigger;
    for(int i = 0; i < dim; i++)tree->SetBranchAddress(vars[i],&var_data[i]);
    if((string)weightVarB != "noweight")tree->SetBranchAddress(weightVarB,&sw);
    tree->SetBranchAddress("Ds_finalState",&Ds_finalState);
    if(Year>10)tree->SetBranchAddress("year",&year);
    else tree->SetBranchAddress("run",&year);
    tree->SetBranchAddress("TriggerCat",&trigger);

    HyperPointSet points( dim );
    for (int i = 0; i < tree->GetEntries(); i++){
    
        tree->GetEntry(i);
	
        if(Year != year && Year != -1) continue;
    	if(FinalState == "KKpi" && Ds_finalState > 2) continue;
        else if(FinalState == "pipipi" && Ds_finalState > 3) continue;
        else if(FinalState == "Kpipi" && Ds_finalState > 4) throw "undefined final state";
        if(Trigger != trigger && Trigger != -1) continue;

        HyperPoint point( dim );
        for(int j = 0; j < dim; j++)point.at(j)= var_data[j];
        point.addWeight(sw);
        points.push_back(point);
    }
    
    /// Get MC
    TFile* fMC = new TFile(ReweightFromA,"OPEN");
    TTree* treeMC =dynamic_cast<TTree*>(fMC->Get("DecayTree"));    

    treeMC->SetBranchStatus("*",0);
    for(int j = 0; j < dim; j++)treeMC->SetBranchStatus(vars[j],1);
    if(weightVarA != "noweight")treeMC->SetBranchStatus(weightVarA,1);
    treeMC->SetBranchStatus("Ds_finalState",1);
    treeMC->SetBranchStatus("year",1);
    treeMC->SetBranchStatus("run",1);
    treeMC->SetBranchStatus("TriggerCat",1);

    vector<double> var_MC(dim,0.);
    double weight = 1;
    int Ds_finalState_MC, year_MC, trigger_MC;
    for(int j = 0; j < dim; j++)treeMC->SetBranchAddress(vars[j],&var_MC[j]);
    if(weightVarA != "noweight")treeMC->SetBranchAddress(weightVarA,&weight);
    treeMC->SetBranchAddress("Ds_finalState",&Ds_finalState_MC);
    if(Year>10)treeMC->SetBranchAddress("year",&year_MC);
    else treeMC->SetBranchAddress("run",&year_MC);
    treeMC->SetBranchAddress("TriggerCat",&trigger_MC);

    HyperPointSet points_MC( dim );
    for (int i = 0; i < treeMC->GetEntries(); i++){
    
        treeMC->GetEntry(i);

        if(Year != year_MC && Year != -1 ) continue;
        if(FinalState == "KKpi" && Ds_finalState_MC > 2) continue;
        else if(FinalState == "pipipi" && Ds_finalState_MC > 3) continue;
        else if(FinalState == "Kpipi" && Ds_finalState_MC > 4) throw "undefined final state";
        if(Trigger != trigger_MC && Trigger != -1) continue;

        HyperPoint point( dim );
        for(int j = 0; j < dim; j++)point.at(j)= var_MC[j]; 
        point.addWeight(weight);
        points_MC.push_back(point);
    }

    /// Define binning based on sample with smaller statistic
    HyperHistogram* histMC;
    HyperHistogram* hist;
    
    TString label;
    for(int j = 0; j < dim; j++) label += "_" + vars[j] ;
    label += "_Ds2" + FinalState;
    if(Year != -1)label += "_" + anythingToString(Year);
    if(Trigger != -1)label+= "_t" + anythingToString(Trigger);

    cout << "minEventsPerBin " << minEventsPerBin << endl;
    
    if(points.getSumW()>points_MC.getSumW()){
         histMC= new HyperHistogram(limits, points_MC, 
                             /*** Name of the binning algorithm you want to use     */
                             HyperBinningAlgorithms::LIKELIHOOD, 
                             /***  The minimum number of events allowed in each bin */
                             /***  from the HyperPointSet provided (points1)        */
                             AlgOption::MinBinContent      (minEventsPerBin),    
                             /*** This minimum bin width allowed. Can also pass a   */
                             /*** HyperPoint if you would like different min bin    */
                             /*** widths for each dimension                         */
                             AlgOption::MinBinWidth        (minBinWidths),
                             /*** If you want to use the sum of weights rather than */
                             /*** the number of events, set this to true.           */    
                             AlgOption::UseWeights         (true),
                             /*** Some algorithms use a random number generator. Set*/
                             /*** the seed here                                     */
                             AlgOption::RandomSeed         (1),
                             /*** What dimesnion would you like to split first? Only*/
                             /*** applies to certain algortihms                     */
                             AlgOption::StartDimension     (0)
                             /*** What dimesnions would you like to bin in?         */
                             //AlgOption::BinningDimensions  (binningDims),
                             /*** Setting this option will make the agorithm draw   */
                             /*** the binning scheme at each iteration              */
                             //AlgOption::DrawAlgorithm("Algorithm")
                             );
            histMC->setNames(HyperName(vars));
            hist= new HyperHistogram( histMC->getBinning() );
            hist->fill(points); 
            /// Draw binning
            if(dim < 3)histMC->draw((string)OutputDir + "weights/binning" + label);
    }
    else {
        hist= new HyperHistogram(limits, points, 
                                   /*** Name of the binning algorithm you want to use     */
                                   HyperBinningAlgorithms::LIKELIHOOD, 
                                   /***  The minimum number of events allowed in each bin */
                                   /***  from the HyperPointSet provided (points1)        */
                                   AlgOption::MinBinContent      (minEventsPerBin),    
                                   /*** This minimum bin width allowed. Can also pass a   */
                                   /*** HyperPoint if you would like different min bin    */
                                   /*** widths for each dimension                         */
                                   AlgOption::MinBinWidth        (minBinWidths),
                                   /*** If you want to use the sum of weights rather than */
                                   /*** the number of events, set this to true.           */    
                                   AlgOption::UseWeights         (true),
                                   /*** Some algorithms use a random number generator. Set*/
                                   /*** the seed here                                     */
                                   AlgOption::RandomSeed         (1),
                                   /*** What dimesnion would you like to split first? Only*/
                                   /*** applies to certain algortihms                     */
                                   AlgOption::StartDimension     (0)
                                   /*** What dimesnions would you like to bin in?         */
                                   //AlgOption::BinningDimensions  (binningDims),
                                   /*** Setting this option will make the agorithm draw   */
                                   /*** the binning scheme at each iteration              */
                                   //AlgOption::DrawAlgorithm("Algorithm")
                                   );
        hist->setNames(HyperName(vars));
        histMC= new HyperHistogram( hist->getBinning() );
        histMC->fill(points_MC); 
        /// Draw binning
        if(dim < 3)hist->draw((string)OutputDir + "weights/binning" + label );
    }
        
    /// Draw density
    hist->normalise(1);
    histMC->normalise(1);
    if(dim < 3)histMC->drawDensity((string)OutputDir + "weights/density_mc" + label );
    if(dim < 3)hist->drawDensity((string)OutputDir + "weights/density_data" + label );

    /// Produce MC correction histo 
    hist->divide(*histMC);
    if(dim < 3){
        hist->draw((string)OutputDir + "weights/weights" + label );
        if(updateAnaNoteHistos)hist->draw("../../../../../TD-AnaNote/latex/figs/dataVsMC/" +(string)OutputDir + "weights/weights" + label );
    }
    hist->save((string)OutputDir + "weights/weights" + label + ".root" );
    
    return;
}

void createSubset(TString file, TString newfile, TString Cut){
    TChain* tree = new TChain("DecayTree");
    tree->Add(file);
    TFile* output = new TFile(newfile,"RECREATE");
    TTree* new_tree = tree->CopyTree(Cut);
    new_tree->Write();
    output->Close();
}



void plotEff(TString Branch,TString TitleX, int bins, double min, double max, vector<TTree*> trees, vector<TString> weights, vector<TString> titles, vector<int> colors, TString label, bool log = false, bool legendLeft = false, bool eff = false){
        
    /// options
    NamedParameter<string> effSubscript("effSubscript", (std::string) "");
    NamedParameter<string> legTitle("legTitle", (std::string) "");
    NamedParameter<string> OutputDir("OutputDir", (std::string) "final/", (char*) 0);
    NamedParameter<int> useLTweight("useLTweight", 0);
    NamedParameter<int> updateAnaNotePlots("updateAnaNotePlots", 0);

    cout << "Plotting " << Branch << endl;

    double var;
    double w=1;
    for(int n = 0; n < trees.size(); n++){
	trees[n]->SetBranchStatus("*",0);
	trees[n]->SetBranchStatus(Branch,1);
	if(weights[n]!="noweight")trees[n]->SetBranchStatus(weights[n],1);
	trees[n]->SetBranchAddress(Branch,&var);
	if(weights[n]!="noweight")trees[n]->SetBranchAddress(weights[n],&w);
    }
	
    ///Make histograms
    TString title;
    if(eff)title = ";"+TitleX+";#varepsilon [a.u.]";
    else title = ";"+TitleX+";Yield [norm.]";
    if(eff & ((string)effSubscript != ""))title.ReplaceAll("epsilon",("epsilon_{"+(string)effSubscript+"}").c_str());
    if(useLTweight)eff = false;

    vector<TH1D*> histos;
    for(int n = 0; n < trees.size(); n++){
	histos.push_back( new TH1D(Branch+anythingToString(n),title,bins,min,max));
	
	int numEvents = trees[n]->GetEntries();
	for(int i=0; i< numEvents; i++)
	{	
		if (0ul == (i % 100000ul)) cout << "Read event " << i << "/" << numEvents << endl;
		trees[n]->GetEntry(i);
		double lt_w = 1.;
		if(useLTweight && Branch == "Bs_DTF_TAU")lt_w = exp(-var/tau)*cosh(dgamma/2.*var);
		if(weights[n]!="noweight")histos[n]->Fill(var,w/lt_w);
		else histos[n]->Fill(var,1./lt_w);
	}

    }    

/*        
    /// Print efficiencies
    cout << endl << "PID efficiencies: " << endl;
    cout << " MC = " << h_MC.Integral()/h_MC_noPID.Integral() * 100. << " ( % ) " <<  endl;
    cout << " MC (PIDGen) = " << h_MC_PIDGen.Integral()/h_MC_noPID.Integral() * 100. << " ( % ) " <<  endl;
    cout << " MC (PIDCorr) = " << h_MC_PIDCorr.Integral()/h_MC_noPID.Integral() * 100. << " ( % ) " <<  endl << endl;
*/
    ///Plot it
    TCanvas c;
    TLegend* leg;
    if(legendLeft)leg = new TLegend(0.15,0.7,0.45,0.9,"");
    else leg = new TLegend(0.7,0.7,0.85,0.9,"");
    leg->SetLineStyle(0);
    leg->SetLineColor(0);
    leg->SetFillColor(0);
    leg->SetTextFont(22);
    leg->SetTextColor(1);
    leg->SetTextSize(0.04);
    leg->SetTextAlign(12);
    if((string)legTitle != "")leg->AddEntry((TObject*)0,((string)legTitle).c_str(), "");

    if(eff){
	gPad->SetLogy(0);
	for(int n = 0; n < trees.size(); n++){
		histos[n]->Scale(1./histos[n]->Integral());
		if(n==0)continue;
		histos[n]->Divide(histos[n],histos[0]);
		histos[n]->SetMinimum(0.);
		histos[n]->SetMaximum(2.5);

		histos[n]->SetMarkerColor(colors[n]);
		histos[n]->SetLineColor(colors[n]);
		if(n==1)histos[n]->Draw("e");
		else histos[n]->Draw("esame");
		TLegendEntry* le = leg->AddEntry(histos[n],titles[n],"LEP");
		le->SetTextColor(colors[n]); 
	}
    }
    else{
	for(int n = 0; n < trees.size(); n++){
	
		histos[n]->Scale(1./histos[n]->Integral());    
		histos[n]->SetMinimum(0.);
		if(log){
			histos[n]->SetMinimum(0.001);
			gPad->SetLogy(1);
		}
		else gPad->SetLogy(0);
		histos[n]->SetMaximum(histos[n]->GetMaximum()*2.0);
		histos[n]->SetMarkerColor(colors[n]);
		histos[n]->SetLineColor(colors[n]);
		if(n==0)histos[n]->Draw("e");
		else histos[n]->Draw("esame");
		TLegendEntry* le = leg->AddEntry(histos[n],titles[n],"LEP");
		le->SetTextColor(colors[n]); 
	}
    }
    leg->Draw(); 
	
    if(eff)label = "eff_" + label;
    if(updateAnaNotePlots)c.Print("../../../../../TD-AnaNote/latex/figs/dataVsMC/" + (string)OutputDir + label + "_"+Branch+".pdf" );
    c.Print((string)OutputDir + label + "_"+Branch+".eps"); 

    if(useLTweight){
	ofstream datafile;
// 	if(updateAnaNote) datafile.open(("../../../../../TD-AnaNote/latex/tables/Acceptance/"+(string)BinningName+"/splineCoeffs_"+ label + ".tex").c_str(),std::ofstream::trunc);
	//else 
	datafile.open("table.tex",std::ofstream::trunc);
	datafile << "\\begin{table}[h]" << "\n";
	datafile << "\\centering" << "\n";
	datafile << "\\caption{} " << "\n";
	datafile << "\\begin{tabular}{ l" ;
	for(int i = 0; i<histos.size(); i++)datafile << " l ";
	datafile << "}" << "\n";
	datafile << "\\hline" << "\n";
	datafile << "\\hline" << "\n";
	datafile << "KS-Test &" ;
	titles[0].ReplaceAll("#","\\") ;
	for(int i = 1; i<histos.size(); i++){
		datafile << "$" << titles[i].ReplaceAll("#","\\")  << "$" ;
		if(i< histos.size()-1 ) datafile << " & " ;
	}
	datafile << "\\\\ \\hline" << "\n";

	for(int i = 0; i<histos.size()-1; i++){
		datafile << "$" << titles[i] << "$ & " ;
		for(int j=1; j < histos.size(); j++){
			//if(j>i)datafile << std::fixed << std::setprecision(3) << histos[i]->Chi2Test(histos[j],"WWP")/(bins-1.);
			if(j>i)datafile << std::fixed << std::setprecision(3) << histos[i]->KolmogorovTest(histos[j]);
			if(j< histos.size()-1)datafile << " & " ;
			else datafile << "\\\\" << "\n";
		}
	}
	datafile << "\\hline" << "\n";
	datafile << "\\hline" << "\n";
	datafile << "\\end{tabular}" << "\n";
	datafile << "\\label{table:}" << "\n";
	datafile << "\\end{table}" << "\n";
    }
  
}

void compareEff(vector<TString> files, vector<TString> weights, vector<TString> Cuts, vector<TString> titles, vector<int> colors, int Year = -1, TString finalState = "all", int Trigger = -1, TString label = ""){
    
    // Cuts
    TString Cut;
    if(Year>10)Cut += " year == " + anythingToString(Year);
    else if(Year == -1)Cut += " year > " + anythingToString(Year);
    else Cut += " run == " + anythingToString(Year);
    if(finalState == "KKpi")Cut += " && Ds_finalState < 3 ";
    else if(finalState == "pipipi")Cut += " && Ds_finalState == 3 ";
    else if(finalState == "Kpipi")Cut += " && Ds_finalState == 4 ";
    if(Trigger != -1) Cut += " && TriggerCat == " + anythingToString(Trigger);

    TFile* output = new TFile("dummy.root","RECREATE");
    vector<TTree*> new_trees;

    cout << endl << "Comparing files " << endl;
    for(int i=0; i < files.size(); i++){
    
	TString cut = Cuts[i];
	if(cut != "")cut += " && ";
	cut += Cut;	

	///Load files
	TFile *file = new TFile(files[i]);
	TTree* tree = (TTree*) file->Get("DecayTree");

	output->cd();
	TTree* new_tree = tree->CopyTree(cut);
	new_trees.push_back(new_tree);

	cout << files[i] << " ( " << cut << " ) " << endl; 
    }

    /// Options
    NamedParameter<int> nBins("nBins", 40); 
    NamedParameter<int> plotOnlyTau("plotOnlyTau", 0); 

    label += "Ds2";
    if(finalState != "")label +=  finalState;
    else label += "all";
    if(Year>-1) label += "_" + anythingToString(Year);
    if(Trigger>-1) label += "_t" + anythingToString(Trigger);
    
    TString Decay, selection;
    if(A_is_in_B("signal",(string)files[0]) && A_is_in_B("signal",(string)files[1])) Decay = "signal";
    if(A_is_in_B("norm",(string)files[0])&& A_is_in_B("norm",(string)files[1])) Decay = "norm";
    if(A_is_in_B("Final",(string)files[0]) && A_is_in_B("Final",(string)files[1])) selection = "Final";

    /// TAU
    plotEff("Bs_DTF_TAU","t(B) [ps]",nBins,0.,10.,new_trees, weights, titles, colors, label,false,true,true);
    if(plotOnlyTau)return;

    //if(selection == "Final")plotEff("BDTG_response","BDTG",nBins,0,1.,new_trees, weights, titles, colors, label,false,true);

    /// Bs
    plotEff("Bs_PT","p_{T}(B) [MeV]",nBins,0,40000,new_trees, weights, titles, colors,  label);
    plotEff("Bs_ETA","#eta(B)",nBins,1,6,new_trees, weights, titles, colors, label);
    plotEff("Bs_DTF_TAUERR","#sigma_{t}(B) [ps]",nBins,0,0.15,new_trees, weights, titles, colors, label);

    /// Dalitz
    plotEff("m_Kpipi","m(K^{+}#pi^{+}#pi^{-})[MeV]",nBins,1000,1950,new_trees, weights, titles, colors, label,false,false,true);
    plotEff("m_Kpi","m(K^{+}#pi^{-})[MeV]",nBins,massKaon+massPion,1200,new_trees, weights, titles, colors, label,false,false,true);
    plotEff("m_pipi","m(#pi^{+}#pi^{-})[MeV]",nBins,2.*massPion,1200,new_trees, weights, titles, colors, label,false,true,true);
    plotEff("m_Dspipi","m(D_{s}^{-}#pi^{+}#pi^{-})[MeV]",nBins,2400,massBs-massKaon,new_trees, weights, titles, colors, label,false,false,true);
    plotEff("m_Dspi","m(D_{s}^{-}#pi^{+})[MeV]",nBins,massDs+massPion,massBs-massKaon-massPion,new_trees, weights, titles, colors, label,false,true,true);
    plotEff("m_DsK","m(D_{s}^{-}K^{+})[MeV]",nBins,0,5500,new_trees, weights, titles, colors, label,false,true,true);
    plotEff("m_DsKpi","m(D_{s}^{-}K^{+}#pi^{-})[MeV]",nBins,1900,5500,new_trees, weights, titles, colors, label,false,true,true);
    plotEff("m_DsKpip","m(D_{s}^{-}K^{+}#pi^{+})[MeV]",nBins,1900,5500,new_trees, weights, titles, colors, label,false,true,true);
}



void rescaleMC(string fileA, string fileB){

    TFile* f = new TFile(fileA.c_str(),"OPEN");
    TTree* tree =dynamic_cast<TTree*>(f->Get("DecayTree"));
	
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("Ds_finalState",1);
    tree->SetBranchStatus("year",1);
    tree->SetBranchStatus("run",1);
    tree->SetBranchStatus("TriggerCat",1);
    tree->SetBranchStatus("weight*",1);

    double sw;
    int Ds_finalState, year, trigger, run;
    tree->SetBranchAddress("weight_rw",&sw);
    tree->SetBranchAddress("Ds_finalState",&Ds_finalState);
    tree->SetBranchAddress("year",&year);
    tree->SetBranchAddress("run",&run);
    tree->SetBranchAddress("TriggerCat",&trigger);

    double sumw = 0;
    double sumw2 = 0;
   
    double sumDs0 = 0;
    double sumDs1 = 0;
    double sumDs2 = 0;
    double sumDs3 = 0;
    double sumDs4 = 0;

    double sumRun1 = 0;
    double sumRun2 = 0;

    double sumTrigger0 = 0;
    double sumTrigger1 = 0;

    TH3D* h = new TH3D("","",2,-0.5,1.5,2,0.5,2.5,5,-0.5,4.5);
    TH2D* h2 = new TH2D("","",2,-0.5,1.5,2,0.5,2.5);

    for (int i = 0; i < tree->GetEntries(); i++){
    
        tree->GetEntry(i);

	sumw += sw;
	sumw2 += sw*sw;

	if(run==1){ 
		sumRun1 += sw;
	}
	else{
		 sumRun2 += sw;
	}
	if(Ds_finalState==0)sumDs0 += sw;
	if(Ds_finalState==1)sumDs1 += sw;
	if(Ds_finalState==2)sumDs2 += sw;
	if(Ds_finalState==3)sumDs3 += sw;
	if(Ds_finalState==4)sumDs4 += sw;

	if(trigger==0)sumTrigger0 += sw;
	if(trigger==1)sumTrigger1 += sw;

	h->Fill((double)trigger,(double)run,(double)Ds_finalState,sw);
	h2->Fill((double)trigger,(double)run,sw);
    }

    TFile* f2 = new TFile(fileB.c_str(),"OPEN");
    TTree* tree2 =dynamic_cast<TTree*>(f2->Get("DecayTree"));
	
    tree2->SetBranchStatus("*",0);
    tree2->SetBranchStatus("Ds_finalState",1);
    tree2->SetBranchStatus("year",1);
    tree2->SetBranchStatus("run",1);
    tree2->SetBranchStatus("TriggerCat",1);
    tree2->SetBranchStatus("weight",1);
    tree2->SetBranchStatus("N_Bs_sw",1);

    tree2->SetBranchAddress("N_Bs_sw",&sw);
    tree2->SetBranchAddress("Ds_finalState",&Ds_finalState);
    tree2->SetBranchAddress("year",&year);
    tree2->SetBranchAddress("run",&run);
    tree2->SetBranchAddress("TriggerCat",&trigger);

    double sumw_data = 0;
    double sumw2_data = 0;
   
    double sumDs0_data = 0;
    double sumDs1_data = 0;
    double sumDs2_data = 0;
    double sumDs3_data = 0;
    double sumDs4_data = 0;

    double sumRun1_data = 0;
    double sumRun2_data = 0;

    double sumTrigger0_data = 0;
    double sumTrigger1_data = 0;

    TH3D* h_data = new TH3D("","",2,-0.5,1.5,2,0.5,2.5,5,-0.5,4.5);
    TH2D* h2_data = new TH2D("","",2,-0.5,1.5,2,0.5,2.5);

    for (int i = 0; i < tree2->GetEntries(); i++){
    
        tree2->GetEntry(i);

	sumw_data += sw;
	sumw2_data += sw*sw;

	if(run==1){ 
		sumRun1_data += sw;
	}
	else{
		sumRun2_data += sw;
	}
	
	if(Ds_finalState==0)sumDs0_data += sw;
	if(Ds_finalState==1)sumDs1_data += sw;
	if(Ds_finalState==2)sumDs2_data += sw;
	if(Ds_finalState==3)sumDs3_data += sw;
	if(Ds_finalState==4)sumDs4_data += sw;

	if(trigger==0)sumTrigger0_data += sw;
	if(trigger==1)sumTrigger1_data += sw;

	h_data->Fill((double)trigger,(double)run,(double)Ds_finalState,sw);
	h2_data->Fill((double)trigger,(double)run,sw);
    }

    double scaleRun1 = sumRun1_data/sumw_data / (sumRun1/sumw);
    double scaleRun2 = sumRun2_data/sumw_data / (sumRun2/sumw);

    double scaleDs0 = sumDs0_data/sumw_data / (sumDs0/sumw);
    double scaleDs1 = sumDs1_data/sumw_data / (sumDs1/sumw);
    double scaleDs2 = sumDs2_data/sumw_data / (sumDs2/sumw);
    double scaleDs3 = sumDs3_data/sumw_data / (sumDs3/sumw);
    double scaleDs4 = sumDs4_data/sumw_data / (sumDs4/sumw);

    double scaleTrigger0 = sumTrigger0_data/sumw_data / (sumTrigger0/sumw);
    double scaleTrigger1 = sumTrigger1_data/sumw_data / (sumTrigger1/sumw);

    h->Scale(1./h->Integral());
    h_data->Scale(1./h_data->Integral());
    TH3D *h_weight =(TH3D *)h->Clone();
    h_weight->Divide(h_data,h);

    h2->Scale(1./h2->Integral());
    h2_data->Scale(1./h2_data->Integral());
    TH2D *h2_weight =(TH2D *)h2->Clone();
    h2_weight->Divide(h2_data,h2);

    TFile* output=new TFile( ((TString)fileA).ReplaceAll(".root","_scaled.root"),"RECREATE");
    double weight_scaled;

    tree->SetBranchStatus("*",1);
    tree->SetBranchStatus("*ENDVERTEX*",0);
    tree->SetBranchStatus("*OWNPV*",0);
    tree->SetBranchStatus("*ORIVX*",0);
    tree->SetBranchStatus("*TOPPV*",0);
    tree->SetBranchStatus("*TRUE*VERTEX*",0);
    tree->SetBranchStatus("Bs_*_DEC",0);
    tree->SetBranchStatus("Bs_*_PROB",0);
    tree->SetBranchStatus("Bs_B0DTF_*",0);
    tree->SetBranchStatus("Bs_DTF_*",0);
    tree->SetBranchStatus("Bs_BsDTF_*",0);
    tree->SetBranchStatus("Bs_PV_*",0);
    tree->SetBranchStatus("*BsTaggingTool*",0);
    tree->SetBranchStatus("*_PP_*",0);
    tree->SetBranchStatus("*ProtoParticles*",0);
    tree->SetBranchStatus("*gen*",0);
    tree->SetBranchStatus("*corr*",0);
    tree->SetBranchStatus("*CHI2*",1);
    tree->SetBranchStatus("*TAU*",1);
    tree->SetBranchStatus("Bs_DTF_MM*",1);
    tree->SetBranchStatus("*DIRA*",1);

    tree->SetBranchStatus("weight",0);
    TTree* new_tree = tree->CopyTree("");
    TBranch* b_weight_scaled = new_tree->Branch("weight", &weight_scaled, "weight/D");
    tree->SetBranchStatus("weight",1);

    double sumw_new = 0;
    double sumw2_new = 0;

    for (int i = 0; i < tree->GetEntries(); i++){
    
        tree->GetEntry(i);
	
	weight_scaled = sw;
// 	weight_scaled *= h_weight->GetBinContent(h_weight->FindBin((double)trigger,(double)run,(double)Ds_finalState));
	weight_scaled *= h2_weight->GetBinContent(h2_weight->FindBin((double)trigger,(double)run));

	sumw_new += weight_scaled;
	sumw2_new += weight_scaled*weight_scaled;

	b_weight_scaled->Fill();
    }

    double sumDs0_new = 0;
    double sumDs1_new = 0;
    double sumDs2_new = 0;
    double sumDs3_new = 0;
    double sumDs4_new = 0;
    double sumRun1_new = 0;
    double sumRun2_new = 0;
    double sumTrigger0_new = 0;
    double sumTrigger1_new = 0;
    
    for (int i = 0; i < new_tree->GetEntries(); i++){
    
        new_tree->GetEntry(i);

	if(run==1)sumRun1_new += weight_scaled;
	else sumRun2_new += weight_scaled;

	if(Ds_finalState==0)sumDs0_new += weight_scaled;
	if(Ds_finalState==1)sumDs1_new += weight_scaled;
	if(Ds_finalState==2)sumDs2_new += weight_scaled;
	if(Ds_finalState==3)sumDs3_new += weight_scaled;
	if(Ds_finalState==4)sumDs4_new += weight_scaled;

	if(trigger==0)sumTrigger0_new += weight_scaled;
	if(trigger==1)sumTrigger1_new += weight_scaled;
    }

    cout << "scales before rw " << endl;
    cout << scaleRun1 << endl;
    cout << scaleRun2 << endl << endl;

    cout << scaleDs0 << endl;
    cout << scaleDs1 << endl;
    cout << scaleDs2 << endl;
    cout << scaleDs3 << endl;
    cout << scaleDs4 << endl << endl;

    cout << scaleTrigger0 << endl;
    cout << scaleTrigger1 << endl << endl;

    cout << "scales after rw " << endl;
    scaleRun1 = sumRun1_data/sumw_data / (sumRun1_new/sumw_new);
    scaleRun2 = sumRun2_data/sumw_data / (sumRun2_new/sumw_new);

    scaleDs0 = sumDs0_data/sumw_data / (sumDs0_new/sumw_new);
    scaleDs1 = sumDs1_data/sumw_data / (sumDs1_new/sumw_new);
    scaleDs2 = sumDs2_data/sumw_data / (sumDs2_new/sumw_new);
    scaleDs3 = sumDs3_data/sumw_data / (sumDs3_new/sumw_new);
    scaleDs4 = sumDs4_data/sumw_data / (sumDs4_new/sumw_new);

    scaleTrigger0 = sumTrigger0_data/sumw_data / (sumTrigger0_new/sumw_new);
    scaleTrigger1 = sumTrigger1_data/sumw_data / (sumTrigger1_new/sumw_new);


    cout << scaleRun1 << endl;
    cout << scaleRun2 << endl << endl;

    cout << scaleDs0 << endl;
    cout << scaleDs1 << endl;
    cout << scaleDs2 << endl;
    cout << scaleDs3 << endl;
    cout << scaleDs4 << endl << endl;

    cout << scaleTrigger0 << endl;
    cout << scaleTrigger1 << endl;


    cout << "effective weights" << endl;
    cout << sumw/sumw2 << endl;
    cout << sumw_new/sumw2_new << endl;

    new_tree->Write();
    output->Close();
}

int main(int argc, char** argv){
    
    time_t startTime = time(0);
    
//     rescaleMC("/auto/data/dargent/BsDsKpipi/Final/MC/signal.root", "/auto/data/dargent/BsDsKpipi/Final/Data/signal.root");
//     rescaleMC("/auto/data/dargent/BsDsKpipi/Final/MC/norm.root", "/auto/data/dargent/BsDsKpipi/Final/Data/norm.root");

//     rescaleMC("/auto/data/dargent/BsDsKpipi/Final/MC/signal_PIDGen.root", "/auto/data/dargent/BsDsKpipi/Final/Data/signal.root");
//     rescaleMC("/auto/data/dargent/BsDsKpipi/Final/MC/norm_PIDGen.root", "/auto/data/dargent/BsDsKpipi/Final/Data/norm.root");
// 
//     rescaleMC("/auto/data/dargent/BsDsKpipi/Final/MC/signal_PIDMC.root", "/auto/data/dargent/BsDsKpipi/Final/Data/signal.root");
//     rescaleMC("/auto/data/dargent/BsDsKpipi/Final/MC/norm_PIDMC.root", "/auto/data/dargent/BsDsKpipi/Final/Data/norm.root");
// 
//     rescaleMC("/auto/data/dargent/BsDsKpipi/Final/MC/signal_noBDT.root", "/auto/data/dargent/BsDsKpipi/Final/Data/signal.root");
//     rescaleMC("/auto/data/dargent/BsDsKpipi/Final/MC/norm_noBDT.root", "/auto/data/dargent/BsDsKpipi/Final/Data/norm.root");
// 
//    return 0;

    //createSubset("../Files/Final/Data/norm.root","../Files/Final/Data/norm_t0.root","TriggerCat == 0");
    //createSubset("../Files/Final/Data/norm.root","../Files/Final/Data/norm_t1.root","TriggerCat == 1");
    //createSubset("../Files/Final/Data/norm.root","../Files/Final/Data/norm_r1.root","run == 1");
    //createSubset("../Files/Final/Data/norm.root","../Files/Final/Data/norm_r2.root","run == 2");

    /// Options for reweighting
    NamedParameter<int> checkPID("checkPID", 0); 
    NamedParameter<int> checkEff("checkEff", 0); 

    NamedParameter<string> ReweightFromA("ReweightFromA", (std::string) "/auto/data/dargent/BsDsKpipi/Preselected/MC/norm.root");
    NamedParameter<string> ReweightToB("ReweightToB", (std::string) "/auto/data/dargent/BsDsKpipi/Preselected/Data/norm.root");
    NamedParameter<string> ApplyWeightToC("ApplyWeightToC", (std::string) "");

    NamedParameter<string> weightVarA("weightVarA", (std::string) "weight");
    NamedParameter<string> newWeightVarA("newWeightVarA", (std::string) "weight");
    NamedParameter<string> weightVarB("weightVarB", (std::string) "N_Bs_sw");
    NamedParameter<string> weightVarC("weightVarC", (std::string) "noweight");
    NamedParameter<string> newWeightVarC("newWeightVarC", (std::string) "noweight");

    NamedParameter<string> cutA("cutA", (std::string) "");
    NamedParameter<string> cutB("cutB", (std::string) "");

    NamedParameter<int> nIterations("nIterations", 1); 
    NamedParameter<int> reweight("reweight", 1); 
    NamedParameter<int> reweightInBinsOfRun("reweightInBinsOfRun", 1); 
    NamedParameter<int> reweightInBinsOfFinalState("reweightInBinsOfFinalState", 1); 
    NamedParameter<int> reweightInBinsOfTrigger("reweightInBinsOfTrigger", 1); 

    NamedParameter<int> reweightVarSet1("reweightVarSet1", 1); 
    NamedParameter<int> reweightVarSet2("reweightVarSet2", 1); 
    NamedParameter<int> reweightVarSet3("reweightVarSet3", 1); 
    NamedParameter<int> reweightVarSet4("reweightVarSet4", 1); 
    NamedParameter<int> reweightVarSet5("reweightVarSet5", 1); 

    /// Options for efficiency comparisons
    NamedParameter<string> file1("file1", (std::string) "");
    NamedParameter<string> file2("file2", (std::string) "");
    NamedParameter<string> file3("file3", (std::string) "");
    NamedParameter<string> file4("file4", (std::string) "");
    NamedParameter<string> file5("file5", (std::string) "");

    NamedParameter<string> weight1("weight1", (std::string) "noweight");
    NamedParameter<string> weight2("weight2", (std::string) "noweight");
    NamedParameter<string> weight3("weight3", (std::string) "noweight");
    NamedParameter<string> weight4("weight4", (std::string) "noweight");
    NamedParameter<string> weight5("weight5", (std::string) "noweight");

    NamedParameter<string> cut1("cut1", (std::string) "");
    NamedParameter<string> cut2("cut2", (std::string) "");
    NamedParameter<string> cut3("cut3", (std::string) "");
    NamedParameter<string> cut4("cut4", (std::string) "");
    NamedParameter<string> cut5("cut5", (std::string) "");

    NamedParameter<string> title1("title1", (std::string) "");
    NamedParameter<string> title2("title2", (std::string) "");
    NamedParameter<string> title3("title3", (std::string) "");
    NamedParameter<string> title4("title4", (std::string) "");
    NamedParameter<string> title5("title5", (std::string) "");


    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    gStyle->SetPalette(1);
    
    vector<int> years;
    if(reweightInBinsOfRun==1){
        years.push_back(1); // Means run 1
        years.push_back(2); // Means run 2
    }
    else if(reweightInBinsOfRun==0){
        //years.push_back(11);
        //years.push_back(12);
        //years.push_back(15);
        years.push_back(16);
        //years.push_back(17);
    }
    else years.push_back(-1); // Means all 
    
    vector<TString> Ds_finalStates;
    if(reweightInBinsOfFinalState){
        Ds_finalStates.push_back("KKpi");
        Ds_finalStates.push_back("pipipi");
        Ds_finalStates.push_back("Kpipi");
    }
    else Ds_finalStates.push_back("all");
    
    vector<int> trigger;
    if(reweightInBinsOfTrigger){
        trigger.push_back(0);
        trigger.push_back(1);
    }
    else trigger.push_back(-1);

    /// Define reweighting vars
    vector<TString>  vars_1, vars_2, vars_3, vars_4, vars_5;
    vector<double>   min_1, min_2, min_3, min_4, min_5;
    vector<double>   max_1, max_2, max_3, max_4, max_5;
     
    vars_1.push_back("Bs_PT");
    min_1.push_back(0.);
    max_1.push_back(100000.);
    vars_1.push_back("Bs_ETA");
    min_1.push_back(1.5);
    max_1.push_back(5.5);
    vars_1.push_back("NTracks");
    min_1.push_back(0.);
    max_1.push_back(1000.);
    
    vars_2.push_back("NTracks");
    min_2.push_back(0.);
    max_2.push_back(1000.);
//     vars_2.push_back("max_ghostProb");
//     min_2.push_back(0.);
//     max_2.push_back(0.4);

    //vars_3.push_back("BDTG");
    //min_3.push_back(0.);
    //max_3.push_back(1.0);
    vars_3.push_back("Bs_PT");
    min_3.push_back(0.);
    max_3.push_back(200000.);
    //vars_3.push_back("Bs_ETA");
    //min_3.push_back(1.5);
    //max_3.push_back(6.5);
    vars_3.push_back("Bs_BsDTF_TAUERR");
    min_3.push_back(0.);
    max_3.push_back(0.15);


    vars_4.push_back("Ds_FDCHI2_ORIVX");
    min_4.push_back(0.);
    max_4.push_back(40000.);

//     vars_4.push_back("Bs_PT");
//     min_4.push_back(0.);
//     max_4.push_back(200000.);
     vars_4.push_back("Bs_BsDTF_TAUERR");
     min_4.push_back(0.);
     max_4.push_back(0.15);
    
//     vars_5.push_back("Bs_PT");
//     min_5.push_back(0.);
//     max_5.push_back(200000.);
//     vars_5.push_back("Ds_PT");
//     min_5.push_back(0.);
//     max_5.push_back(200000.);
    vars_5.push_back("pi_PT");
    min_5.push_back(0.);
    max_5.push_back(200000.);

    vector< vector<TString> > vars_set;
    vector< vector<double> > min_set;
    vector< vector<double> > max_set;

    if(reweightVarSet1)vars_set.push_back(vars_1);
    if(reweightVarSet2)vars_set.push_back(vars_2);
    if(reweightVarSet3)vars_set.push_back(vars_3);
    if(reweightVarSet4)vars_set.push_back(vars_4);
    if(reweightVarSet5)vars_set.push_back(vars_5);
    
    if(reweightVarSet1)min_set.push_back(min_1);
    if(reweightVarSet2)min_set.push_back(min_2);
    if(reweightVarSet3)min_set.push_back(min_3);
    if(reweightVarSet4)min_set.push_back(min_4);
    if(reweightVarSet5)min_set.push_back(min_5);

    if(reweightVarSet1)max_set.push_back(max_1);
    if(reweightVarSet2)max_set.push_back(max_2);
    if(reweightVarSet3)max_set.push_back(max_3);
    if(reweightVarSet4)max_set.push_back(max_4);
    if(reweightVarSet5)max_set.push_back(max_5);

    /// Compare PID/ Eff

    if(checkEff){ 
	vector<TString> files;
	files.push_back((string) file1);
	files.push_back((string) file2);
	if((string)file3 != "")files.push_back((string) file3);
	if((string)file4 != "")files.push_back((string) file4);
	if((string)file5 != "")files.push_back((string) file5);
	
	vector<TString> weights;
	weights.push_back((string) weight1);
	weights.push_back((string) weight2);
	if((string)file3 != "")weights.push_back((string) weight3);
	if((string)file4 != "")weights.push_back((string) weight4);
	if((string)file5 != "")weights.push_back((string) weight5);
	
	vector<TString> cuts;
	cuts.push_back((string) cut1);
	cuts.push_back((string) cut2);
	if((string)file3 != "")cuts.push_back((string) cut3);
	if((string)file4 != "")cuts.push_back((string) cut4);
	if((string)file5 != "")cuts.push_back((string) cut5);
	
	vector<TString> titles;
	titles.push_back((string)title1);
	titles.push_back((string)title2);
	if((string)file3 != "")titles.push_back((string)title3);
	if((string)file4 != "")titles.push_back((string)title4);
	if((string)file5 != "")titles.push_back((string)title5);
	
	vector<int> colors;
	colors.push_back(kBlack);
	colors.push_back(kRed);
	colors.push_back(kBlue);
	colors.push_back(kMagenta+1);
	colors.push_back(kGreen+3);
		
	for(int i= 0; i < years.size(); i++) for(int j= 0; j < Ds_finalStates.size(); j++)for(int k= 0; k < trigger.size(); k++){ 
			compareEff(files,weights,cuts, titles, colors, years[i], Ds_finalStates[j],trigger[k]);
	}
	return 0;

    }

    if(checkPID){  
	for(int i= 0; i < years.size(); i++) for(int j= 0; j < Ds_finalStates.size(); j++)for(int k= 0; k < trigger.size(); k++){ 
		comparePID((string) ReweightToB, (string) ReweightFromA, (string) weightVarB, (string) weightVarA, (string) cutB, (string) cutA, years[i], Ds_finalStates[j],trigger[k]);
	}
	return 0;
    }

    /// Produce MC correction histos and apply weights
    /// Weights are applied on top of each other with the previous weighting applied
    TString weightA = (string) weightVarA;
    TString weightC = (string) weightVarC;
    
    if(reweight)for(int n= 0; n < nIterations; n++)
			for(int i= 0; i < years.size(); i++) 
				for(int j= 0; j < Ds_finalStates.size(); j++)	
					for(int k= 0; k < trigger.size(); k++)    
						for(int l =0; l < vars_set.size(); l++)	{	
							produceCorrectionHisto(vars_set[l],min_set[l],max_set[l],years[i],Ds_finalStates[j],trigger[k], (string)ReweightFromA, weightA, (string) ReweightToB, (string) weightVarB);
							applyCorrectionHisto(vars_set[l],years[i],Ds_finalStates[j],trigger[k],(string)ReweightFromA, weightA, (string)newWeightVarA);    
							if((string)ApplyWeightToC != "")applyCorrectionHisto(vars_set[l],years[i],Ds_finalStates[j],trigger[k],(string)ApplyWeightToC, weightC, (string)newWeightVarC);
						}
						
    /// Draw comparison plots 
    for(int i= 0; i < years.size(); i++) for(int j= 0; j < Ds_finalStates.size(); j++)for(int k= 0; k < trigger.size(); k++){ 
        compare((string) ReweightToB, (string) ReweightFromA, (string) weightVarB, (string) weightVarA, (string) newWeightVarA, (string) cutB, (string) cutA, years[i], Ds_finalStates[j],trigger[k]);
    }

    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
