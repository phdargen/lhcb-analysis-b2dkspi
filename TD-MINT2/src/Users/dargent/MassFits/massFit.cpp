// Fits B mass distribution and calculates sweights
// author: Philippe d'Argent, Matthieu Kecke
#include <cmath>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>
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
#include <TFitResult.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TNtuple.h>
#include "TRandom3.h"
#include <sstream>
#include <RooDataSet.h>
#include <RooMCStudy.h>
#include "RooGaussModel.h"
#include "RooExponential.h"
#include "RooGenericPdf.h"
#include "RooChebychev.h"
#include "RooAddModel.h"
#include "RooPolynomial.h"
#include "RooTruthModel.h"
#include "RooFitResult.h"
#include "RooDecay.h"
#include "RooPlot.h"
#include "RooGaussian.h"
#include "RooProdPdf.h"
#include "RooDstD0BG.h"
#include "RooAddPdf.h"
#include "RooExtendPdf.h"
#include "RooDataHist.h"
#include "RooCBShape.h"
#include "RooCategory.h"
#include "RooMappedCategory.h"
#include "RooFormulaVar.h"
#include "RooSimultaneous.h"
#include "RooHist.h"
#include "RooStats/SPlot.h"
#include "RooTreeDataStore.h"
#include "RooBinning.h"
#include "RooBifurGauss.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooConstVar.h"
#include "RooWorkspace.h"
#include "RooNDKeysPdf.h"
#include "RooKeysPdf.h"
#include "RooJohnsonSU.h"
#include "RooSimPdfBuilder.h"
#include "RooNumIntConfig.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <ctime>
#include "Mint/NamedParameter.h"
#include "Mint/Utils.h"
#include "Mint/RooHILLdini.h"
#include "Mint/RooHORNSdini.h"
//#include "Mint/CLHEPPhysicalConstants.h"
//#include "Mint/CLHEPSystemOfUnits.h"
using namespace std;
using namespace RooFit ;
using namespace RooStats;
using namespace MINT;

vector<TString> str_year;
vector<TString> str_run;
vector<TString> str_KsCat;
vector<TString> str_trigger;
static const double massKaon = 493.68;
static const double massPion = 139.57;

vector<double> fitPartRecoBkgShape(){
     
    double min_MM = 5050;
    double max_MM = 5280;
    RooRealVar Bs_MM("B_MM", "m(D*K#pi)", min_MM, max_MM,"MeV");
    //mean of gaussians
    RooRealVar mean1("mean1","mu", 5059.,5040.,5070.);
    RooRealVar mean2("mean2","mu", 5182.,5140.,5205.);
    RooRealVar mean3("mean3","mu", 5285.,5270.,5300.);
    //width of gaussians
    RooRealVar sigmaL1("sigma_{1L}", "sigmaL1", 25.9,15.,40.);
    RooRealVar sigmaR1("sigma_{1R}", "sigmaR1", 99.4,50.,115.);
    RooRealVar sigmaL2("sigma_{2L}", "sigmaL2", 13.1,5.,100.);
    RooRealVar sigmaR2("sigma_{2R}", "sigmaR2", 49.5,25.,70.);
    RooRealVar sigmaL3("sigma_{3L}", "sigmaL3", 107.,10.,125.);
    RooRealVar sigmaR3("sigma_{3R}", "sigmaR3", 21.1,5.,33.);
    //bifurcated gaussians
    RooBifurGauss BifGauss1("BifGauss1","BifGauss1", Bs_MM, mean1, sigmaL1,sigmaR1);
    RooBifurGauss BifGauss2("BifGauss2","BifGauss2", Bs_MM, mean2, sigmaL2,sigmaR2);
    RooBifurGauss BifGauss3("BifGauss3","BifGauss3", Bs_MM, mean3, sigmaL3,sigmaR3);
    //fractions of gauss functions
    RooRealVar f_1("f_{1}", "fraction1", 0.405, 0., 1.);
    RooRealVar f_2("f_{2}", "fraction2", 0.1, 0., 1.);
    //add all gaussians
    RooAbsPdf* pdf=new RooAddPdf("BkgShape", "BkgShape", RooArgList(BifGauss1, BifGauss2, BifGauss3), RooArgList(f_1,f_2),kTRUE);
    
    /// Load file
    NamedParameter<string> inFileNamePartReco("inFileNamePartReco",(string)"../../../../../Selection/bkg/bkg_bs2dstarkspi.root");
    TFile* f = new TFile(((string)inFileNamePartReco).c_str());
    TTree* tree = (TTree*) f->Get("DecayTree");
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("B*M",1);
    tree->SetBranchStatus("*BKGCAT",1);
    tree->SetBranchStatus("*TRUEID",1);

    TFile* output = new TFile("dummy.root","RECREATE");    
    TTree* out_tree = tree->CopyTree(("B_MM >= " + anythingToString((double)min_MM) + " && B_MM <= " + anythingToString((double)max_MM) + "&& B_BKGCAT == 50" ).c_str() );
    RooDataSet* data = new RooDataSet("data","data",RooArgSet(Bs_MM),Import(*out_tree));
    
    /// Fit
    RooFitResult *result = pdf->fitTo(*data,Save(kTRUE),NumCPU(2));
    cout << "result is --------------- "<<endl;
    result->Print(); 
    
    //plot mass distribution and fit results
    TCanvas* c1= new TCanvas("");
    RooPlot* frame_m= Bs_MM.frame();
    frame_m->SetTitle("");
    frame_m->GetYaxis()->SetTitle("Yield (norm.)");
    frame_m->GetXaxis()->SetTitle("m(D*K#pi) [MeV]");
    data->plotOn(frame_m,Name("data"),MarkerSize(1),Binning(50));
    pdf->plotOn(frame_m,Name("pdf"),LineColor(kBlue),LineWidth(3));
    pdf->plotOn(frame_m,Components(BifGauss1),LineColor(kBlue),LineStyle(kDashed),LineWidth(1));
    pdf->plotOn(frame_m,Components(BifGauss2),LineColor(kRed),LineStyle(kDashed),LineWidth(1));
    pdf->plotOn(frame_m,Components(BifGauss3),LineColor(kGreen),LineStyle(kDashed),LineWidth(1));
    frame_m->Draw();
    c1->Print("plots/bkg_bs2dstarkspi.eps");
    
    /// Return fit parameters
    vector<double> params;
    params.push_back(mean1.getVal());
    params.push_back(mean2.getVal());
    params.push_back(mean3.getVal());
    params.push_back(sigmaL1.getVal());
    params.push_back(sigmaR1.getVal());
    params.push_back(sigmaL2.getVal());
    params.push_back(sigmaR2.getVal());
    params.push_back(sigmaL3.getVal());
    params.push_back(sigmaR3.getVal());
    params.push_back(f_1.getVal());
    params.push_back(f_2.getVal());
    
    f->Close();
    
    return params;
}

vector<double> fitSignalShape(TString channel = "signal"){
    
    /// Options
    NamedParameter<string> cut_BDT("cutMC_BDT",(string)"BDTG > 0.");
    NamedParameter<int> fitPreselected("fitPreselected", 0);
    NamedParameter<int> sWeight("sWeightMC", 0);
    NamedParameter<int> useDoubleRJ("useDoubleRJ", 0);
    NamedParameter<string> initialState("initialState",(string)"Bd");
    NamedParameter<string> finalState("finalState",(string)"DKspi");
    TString outDir = (string)initialState + "2" + (string)finalState + "/";
    
    double min_MM = 5220. ;
    double max_MM = 5340. ;
    if(channel == "Ks"){
        min_MM = 450 ;
        max_MM = 550. ;        
    }
    
    /// Load file
    NamedParameter<string> inFileNameMC("inFileNameMC",(string)"../../../../../Selection/BDT/signal_mc.root");
    TFile* f = new TFile(((string)inFileNameMC).c_str());
    TTree* tree = (TTree*) f->Get("DecayTree");
    
    if(!sWeight){
        tree->SetBranchStatus("*",0);
        tree->SetBranchStatus("B_BKGCAT",1);
        tree->SetBranchStatus("B_TRUEID",1);
        if(!fitPreselected)tree->SetBranchStatus("BDTG",1);
        tree->SetBranchStatus("run",1);
        tree->SetBranchStatus("year",1);
        tree->SetBranchStatus("TriggerCat",1);
    }
    else {
        tree->SetBranchStatus("*ENDVERTEX*",0);
        tree->SetBranchStatus("*OWNPV*",0);
        tree->SetBranchStatus("*ORIVX*",0);
        tree->SetBranchStatus("*TOPPV*",0);
        tree->SetBranchStatus("*TRUE*VERTEX*",0);
        tree->SetBranchStatus("B_*_DEC",0);
        tree->SetBranchStatus("B_*_PROB",0);
        tree->SetBranchStatus("B_B0DTF_*",0);
        tree->SetBranchStatus("B_DTF_*",0);
        tree->SetBranchStatus("B_BsDTF_*",0);
        tree->SetBranchStatus("B_PV_*",0);
        tree->SetBranchStatus("*BTaggingTool*",0);
        tree->SetBranchStatus("*_PP_*",0);
        tree->SetBranchStatus("*ProtoParticles*",0);
        tree->SetBranchStatus("*gen*",0);
        tree->SetBranchStatus("*corr*",0);
        tree->SetBranchStatus("*CHI2*",1);
        tree->SetBranchStatus("*TAU*",1);
        tree->SetBranchStatus("*DIRA*",1);
    }
    tree->SetBranchStatus("weight",0);
    tree->SetBranchStatus("*m_*",1);
    tree->SetBranchStatus("*MM*",1);
    tree->SetBranchStatus("*ProbNN*",1);

    TFile* output;
    if(!sWeight || fitPreselected) output = new TFile("dummy.root","RECREATE");
    else if(channel == "signal") output = new TFile(TString(inFileNameMC).ReplaceAll("/BDT/","/Final/"),"RECREATE");
    
    static const double massKs = 497.611;
    static const double massPi = ((string)finalState == "DKspi") ? 139.57 : 493.68;
    static const double massD = 1869.65;
    static const double massB = ((string)initialState == "Bd") ? 5279.65 : 5366.88;
    TString prefix_DTF = ((string)initialState == "Bd") ? "FullDTF_" : "FullBsDTF_";

    TString cutPhsp = "( ";
    cutPhsp +=  prefix_DTF+"m_Dpi >=" +  anythingToString(massD+massPi) + " && " + prefix_DTF+"m_Dpi <=" + anythingToString(massB - massKs);
    cutPhsp +=  " && " + prefix_DTF+"m_Kspi >=" +  anythingToString(massKs+massPi) + " && " + prefix_DTF+"m_Kspi <=" +  anythingToString(massB - massD);
    cutPhsp +=  " && " + prefix_DTF+"m_DKs >=" +  anythingToString(massD+massKs) + " && " + prefix_DTF+"m_DKs <=" + anythingToString(massB - massPi);    
    cutPhsp += " ) ";
    
    TTree* out_tree;
    if(fitPreselected)out_tree = tree->CopyTree(("B_DTF_MM >= " + anythingToString((double)min_MM) + " && B_DTF_MM <= " + anythingToString((double)max_MM) ).c_str() );
    else if(channel == "signal") out_tree = tree->CopyTree("B_DTF_MM >= " + anythingToString((double)min_MM) + " && B_DTF_MM <= " + anythingToString((double)max_MM) + " && " + (string)cut_BDT + " && " + cutPhsp  );
    else out_tree = tree->CopyTree("Ks_PV_MM >= " + anythingToString((double)min_MM) + " && Ks_PV_MM <= " + anythingToString((double)max_MM) + " && " + (string)cut_BDT + " && " + cutPhsp );
    
    int bkgCAT,BKGCAT,Bs_TRUEID;
    double sw;
    TBranch* b_bkgCAT = out_tree->Branch("bkgCAT",&bkgCAT,"bkgCAT/I");
    out_tree->SetBranchAddress("B_BKGCAT",&BKGCAT);
    out_tree->SetBranchAddress("B_TRUEID",&Bs_TRUEID);
    TBranch* b_w = out_tree->Branch("weight", &sw, "weight/D");
    
    out_tree->SetBranchStatus("*",0);
    out_tree->SetBranchStatus("B_BKGCAT",1);
    out_tree->SetBranchStatus("bkgCAT",1);
    out_tree->SetBranchStatus("weight",1);
    out_tree->SetBranchStatus("B_DTF_MM",1);
    out_tree->SetBranchStatus("Ks_PV_MM",1);
    
    for(int i= 0; i< out_tree->GetEntries();i++){
        out_tree->GetEntry(i);
        if(BKGCAT == 0 )bkgCAT= 0;
        else if(BKGCAT == 60 || ( BKGCAT == 50 && abs(Bs_TRUEID) == 521) )bkgCAT= 1;
        //else bkgCAT= 1;
        b_bkgCAT->Fill();
    }
    
    TString channelString;
    if(channel == "norm") channelString = "" ;
    if(channel == "signal" && (string)finalState == "DKspi") channelString = "m(DK_{s}#pi)" ;
    else if(channel == "signal") channelString = "m(DK_{s}K)" ;    
    if(channel == "Ks") channelString = "m(K_{s})" ;
    
    TString massVar = "B_DTF_MM";
    if(channel == "Ks") massVar = "Ks_PV_MM";
    RooRealVar DTF_Bs_M(massVar, channelString, min_MM, max_MM,"MeV");
    RooCategory Bs_BKGCAT("bkgCAT","bkgCAT");
    Bs_BKGCAT.defineType("signal",0);
    Bs_BKGCAT.defineType("ghost",1);
    
    RooArgList list =  RooArgList(DTF_Bs_M,Bs_BKGCAT);
    RooDataSet* data = new RooDataSet("data","data",list,Import(*out_tree));
    
    /// Signal pdf
    RooRealVar mean("mean", "mean", 5279, min_MM, max_MM); 
    if(channel == "Ks")mean.setVal(massKaon);
    RooRealVar sigma("sigma", "sigma", 20.,0.,80.); 
    RooRealVar sigma2("sigma2", "sigma2", 50.,0.,80.); 
    RooRealVar gamma("gamma", "gamma", -0.5,-5,5.); 
    RooRealVar delta("delta", "delta", 0.5,-5,5.); 
    RooRealVar gamma2("gamma2", "gamma2", -0.5,-5,5.); 
    RooRealVar delta2("delta2", "delta2", 0.5,-5,5.); 
    RooJohnsonSU* signal1= new RooJohnsonSU("signal1","signal1",DTF_Bs_M, mean,sigma,gamma,delta);
    RooJohnsonSU* signal2= new RooJohnsonSU("signal2","signal2",DTF_Bs_M, mean,sigma2,gamma2,delta2);
    RooRealVar f1("f1", "f1", 0.9);
    RooAddPdf* signal = new RooAddPdf("signal", "signal", RooArgList(*signal1,*signal2), RooArgList(f1));
    if(!useDoubleRJ){
        f1.setVal(1);
        f1.setConstant();
        sigma2.setConstant();
        gamma2.setConstant();
        delta2.setConstant();
    }
    
    RooRealVar mean_ghost("mean_ghost", "mean_ghost", 5279, 5250 , 5300); 
    RooRealVar sigma_ghost("sigma_ghost", "sigma_ghost", 20.,0.,80.); 
    RooRealVar gamma_ghost("gamma_ghost", "gamma_ghost", -0.5,-5,5.); 
    RooRealVar delta_ghost("delta_ghost", "delta_ghost", 0.5,-5,5.); 
    RooJohnsonSU* signal_ghost= new RooJohnsonSU("signal","signal_ghost",DTF_Bs_M, mean,sigma_ghost,gamma,delta);
    
    /// Bkg pdf
    RooRealVar c0_ghost("c0_ghost", "c0_ghost", .0,-10,10); 
    RooRealVar c1_ghost("c1_ghost", "c1_ghost", .0,-10,10); 
    RooRealVar c2_ghost("c2_ghost", "c2_ghost", .0,-10,10); 
    RooChebychev* bkg_ghost= new RooChebychev("bkg_ghost","bkg_ghost",DTF_Bs_M, RooArgList(c0_ghost));
    
    /// Total pdf
    RooRealVar n_sig("n_sig", "n_sig", data->numEntries()*0.9, 0., data->numEntries());
    RooRealVar n_sig_ghost("n_sig_ghost", "n_sig_ghost", data->numEntries()*0.5, 0., data->numEntries());
    RooRealVar n_bkg_ghost("n_bkg_ghost", "n_bkg_ghost", data->numEntries()/10., 0., data->numEntries());
    
    RooAddPdf* pdf_signal = new RooAddPdf("pdf_signal", "pdf_signal", RooArgList(*signal), RooArgList(n_sig));
    RooAddPdf* pdf_ghost = new RooAddPdf("pdf_ghost", "pdf_ghost", RooArgList(*signal_ghost, *bkg_ghost), RooArgList(n_sig_ghost, n_bkg_ghost));
    
    RooSimultaneous* simPdf = new RooSimultaneous("simPdf","simPdf",Bs_BKGCAT);
    simPdf->addPdf(*pdf_signal,"signal");
    simPdf->addPdf(*pdf_ghost,"ghost");
    
    /// Fit
    RooFitResult* result = simPdf->fitTo(*data,Save(kTRUE),NumCPU(2),Extended(kTRUE));
    result->Print();
    
    if(sWeight && channel == "signal"){
        /// Calculate ghost weights
        RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"bkgCAT==bkgCAT::ghost");
        SPlot sPlot("sPlot","sPlot",*data_slice,pdf_ghost,RooArgList(n_sig_ghost,n_bkg_ghost)); 
        
        int n_ij = 0;  /// labels entry number of data slice
        for(int n = 0; n < out_tree->GetEntries(); n++){
            b_bkgCAT->GetEntry(n);
            if(bkgCAT == 1){
                sw = sPlot.GetSWeight(n_ij,"n_sig_ghost_sw");
                n_ij++;
            }
            else sw = 1.;
            b_w->Fill();
        }
    }
    /// Plotting
    TCanvas* c = new TCanvas();
    
    RooPlot* frame= DTF_Bs_M.frame();
    frame->SetTitle("");
    data->plotOn(frame,Name("data"),Binning(50),Cut("bkgCAT==bkgCAT::signal"));
    simPdf->plotOn(frame,Name("signal"),ProjWData(Bs_BKGCAT,*data),Slice(Bs_BKGCAT,"signal"));
    frame->Draw();
    c->Print(outDir+channel+"_MC.eps");
    
    RooPlot* frame2= DTF_Bs_M.frame();
    data->plotOn(frame2,Name("data"),Binning(50),Cut("bkgCAT==bkgCAT::ghost"));
    simPdf->plotOn(frame2,Name("signal"),Slice(Bs_BKGCAT,"ghost"),ProjWData(Bs_BKGCAT,*data));
    simPdf->plotOn(frame2,Name("bkg"),ProjWData(Bs_BKGCAT,*data),Slice(Bs_BKGCAT,"ghost"),LineColor(kRed),LineStyle(kDashed),Components("bkg_ghost"));
    frame2->Draw();
    c->Print(outDir+channel+"_MC_ghost.eps");
    
    /// Return fit params
    vector<double> params;
    params.push_back(mean.getVal());
    params.push_back(sigma.getVal());
    params.push_back(gamma.getVal());
    params.push_back(delta.getVal());
    params.push_back(sigma2.getVal());
    params.push_back(gamma2.getVal());
    params.push_back(delta2.getVal());
    params.push_back(f1.getVal());
    
    cout << endl << "Fraction of signal classified as ghosts = " << n_sig_ghost.getVal()/(n_sig.getVal()+n_sig_ghost.getVal()) << endl;
    
    out_tree->SetBranchStatus("*",1);
    out_tree->Write();
    output->Close();
    return params;
}

void fitSignal(){

	///Options
	NamedParameter<int> numCPU("numCPU", 1);
	NamedParameter<int> sWeight("sWeightSignal", 0);
	NamedParameter<int> nBins("nBins", 80);
	NamedParameter<double> min_MM("min_MM",5100.);
	NamedParameter<double> max_MM("max_MM",5700.);
	NamedParameter<string> cut_BDT("cut_BDT",(string)"");
	NamedParameter<string> inFileName("inFileNameSignal",(string)"../../../../../Selection/BDT/signal_data.root");
	NamedParameter<string> outFileName("outFileNameSignal",(string)"./signal_data.root");
    NamedParameter<int> fixSignalShapeFromMC("fixSignalShapeFromMC", 1);
    NamedParameter<int> fixExpBkgFromSidebands("fixExpBkgFromSidebands", 0);

	///Load file
	TFile *file= new TFile(((string)inFileName).c_str());
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("B_DTF_MM",1);
    tree->SetBranchStatus("BDTG",1);
    tree->SetBranchStatus("m_Dpi",1);
    tree->SetBranchStatus("D_PV_MM",1);
    tree->SetBranchStatus("Ks_PV_MM",1);
    tree->SetBranchStatus("pi_ProbNNpi",1);
    tree->SetBranchStatus("pi_ProbNNk",1);
    tree->SetBranchStatus("KsCat",1);
    tree->SetBranchStatus("run",1);
    tree->SetBranchStatus("TriggerCat",1);

    RooRealVar DTF_B_M("B_DTF_MM", "m(D^{-}K_{S}#pi^{+})", min_MM, max_MM,"MeV/c^{2}");	
    RooRealVar BDTG("BDTG", "BDTG", 0.);    
    RooRealVar m_Dpi("m_Dpi", "m_Dpi", 0.);    
    RooRealVar D_PV_MM("D_PV_MM", "D_PV_MM", 0.);    
    RooRealVar Ks_PV_MM("Ks_PV_MM", "Ks_PV_MM", 0.);    
    RooRealVar pi_ProbNNpi("pi_ProbNNpi", "pi_ProbNNpi", 0.);    
    RooRealVar pi_ProbNNk("pi_ProbNNk", "pi_ProbNNk", 0.);    
    
    ///define categories
    RooCategory run("run","run") ;
    run.defineType("Run1",1) ;
    run.defineType("Run2",2) ;
    RooCategory TriggerCat("TriggerCat","TriggerCat") ;
    TriggerCat.defineType("t0",0);
    TriggerCat.defineType("t1",1);
    RooCategory KsCat("KsCat","KsCat") ;
    KsCat.defineType("LL",0);
    KsCat.defineType("DD",1);

    RooArgList list =  RooArgList(DTF_B_M,BDTG,m_Dpi,D_PV_MM,Ks_PV_MM,pi_ProbNNk,pi_ProbNNpi);
    RooArgList list2 =  RooArgList(KsCat,run,TriggerCat);
    list.add(list2);
    RooDataSet*  data = new RooDataSet("data","data",tree,list,((string)cut_BDT).c_str());	
    
	/// Signal Pdf
    vector<double> sig_params = fitSignalShape("signal");
    RooRealVar mean_MC("mean_MC", "#mu MC", sig_params[0]); 
    RooRealVar sigma_MC("sigma_MC", "#sigma MC", sig_params[1]);
    RooRealVar sigma2_MC("sigma2_MC", "#sigma2 MC", sig_params[4]);
    RooRealVar scale_mean("scale_mean", "scale #mu",1.); 
    RooRealVar scale_sigma("scale_sigma", "scale #sigma", 1.,0.5,2);
    
    RooFormulaVar mean("mean","@0 * @1", RooArgSet(scale_mean,mean_MC)); 
    RooFormulaVar sigma("sigma","@0 * @1", RooArgSet(scale_sigma,sigma_MC)); 
    RooFormulaVar sigma2("sigma2","@0 * @1", RooArgSet(scale_sigma,sigma2_MC)); 
    
    RooRealVar alpha("alpha", "#alpha", sig_params[2],-5.,5.); 
    RooRealVar beta("beta", "#beta", sig_params[3],-5.,5.);
    RooRealVar alpha2("alpha2", "#alpha2", sig_params[5],-5.,5.); 
    RooRealVar beta2("beta2", "#beta2", sig_params[6],-5.,5.);
    RooRealVar f_RJ("f_RJ", "f_RJ", sig_params[7]);
    
    if(fixSignalShapeFromMC){
        alpha.setConstant();
        alpha2.setConstant();
        beta.setConstant();
        beta2.setConstant();
    } 
	RooJohnsonSU signal("signal","signal",DTF_B_M, mean,sigma,alpha,beta);

	/// Bs pdf
	RooFormulaVar mean_Bs("mean_B0","@0 + @1", RooArgSet(mean,RooConst(87.33))); 
	RooJohnsonSU signal_Bs("signal_Bs","signal_Bs",DTF_B_M, mean_Bs,sigma,alpha,beta);
	
	/// Combinatorial bkg pdf
	RooRealVar exp_par("exp_par","exp_par",-1.6508e-03,-10.,10.);
    RooExponential bkg_exp("bkg_exp","bkg_exp",DTF_B_M,exp_par);
    bkg_exp.fitTo(*data,Save(kTRUE),Range(5450.,max_MM));
    if(fixExpBkgFromSidebands)exp_par.setConstant();

    /// Part. reco bkg
    vector<double> bkg_partReco_params = fitPartRecoBkgShape();
    RooRealVar mean1("mean1","mean1", bkg_partReco_params[0]);
    RooRealVar mean2("mean2","mean2", bkg_partReco_params[1]);
    RooRealVar mean3("mean3","mean3", bkg_partReco_params[2]);
    RooRealVar mean1Shifted("mean1Shifted","mean1Shifted", mean1.getVal() - 87.33 );
    RooRealVar mean2Shifted("mean2Shifted","mean2Shifted", mean2.getVal() - 87.33 );
    RooRealVar mean3Shifted("mean3Shifted","mean3Shifted", mean3.getVal() - 87.33 );
    RooRealVar sigmaL1("sigmaL1", "sigmaL1",  bkg_partReco_params[3]);
    RooRealVar sigmaR1("sigmaR1", "sigmaR1",  bkg_partReco_params[4]);
    RooRealVar sigmaL2("sigmaR1", "sigmaL2",  bkg_partReco_params[5]);
    RooRealVar sigmaR2("sigmaR2", "sigmaR2",  bkg_partReco_params[6]);
    RooRealVar sigmaL3("sigmaL3", "sigmaL3",  bkg_partReco_params[7]);
    RooRealVar sigmaR3("sigmaR3", "sigmaR3",  bkg_partReco_params[8]);
    RooRealVar f_1("f_1", "f_1", bkg_partReco_params[9]);
    RooRealVar f_2("f_2", "f_2", bkg_partReco_params[10]);
    
    RooBifurGauss BifGauss1_Bs("BifGauss1_Bs","BifGauss1_Bs", DTF_B_M, mean1, sigmaL1,sigmaR1);
    RooBifurGauss BifGauss2_Bs("BifGauss2_Bs","BifGauss2_Bs", DTF_B_M, mean2, sigmaL2,sigmaR2);
    RooBifurGauss BifGauss3_Bs("BifGauss3_Bs","BifGauss3_Bs", DTF_B_M, mean3, sigmaL3,sigmaR3);
    RooAddPdf bkg_partReco_Bs("bkg_partReco_Bs", "bkg_partReco_Bs", RooArgList(BifGauss1_Bs, BifGauss2_Bs, BifGauss3_Bs), RooArgList(f_1,f_2),kTRUE);

    RooBifurGauss BifGauss1_B0("BifGauss1_B0","BifGauss1_B0", DTF_B_M, mean1Shifted, sigmaL1,sigmaR1);
    RooBifurGauss BifGauss2_B0("BifGauss2_B0","BifGauss2_B0", DTF_B_M, mean2Shifted, sigmaL2,sigmaR2);
    RooBifurGauss BifGauss3_B0("BifGauss3_B0","BifGauss3_B0", DTF_B_M, mean3Shifted, sigmaL3,sigmaR3);
    RooAddPdf bkg_partReco_B0("bkg_partReco_B0", "bkg_partReco_B0", RooArgList(BifGauss1_B0, BifGauss2_B0, BifGauss3_B0), RooArgList(f_1,f_2),kTRUE);

    RooRealVar partReco_f("partReco_f", "partReco_f", 0.5,0,1);
    RooAddPdf* bkg_partReco = new RooAddPdf("bkg_partReco", "bkg_partReco", RooArgList(bkg_partReco_Bs, bkg_partReco_B0), RooArgList(partReco_f));
    
	/// Total pdf
    RooRealVar n_sig("n_sig", "n_sig", data->numEntries()*0.15, 0., data->numEntries());
    RooRealVar n_sig_Bs("n_sig_Bs", "n_sig_B0", data->numEntries()*0.05, 0., data->numEntries());
    RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()*0.7, 0., data->numEntries());
    RooRealVar n_partReco_bkg("n_partReco_bkg", "n_partReco_bkg", data->numEntries()*0.05, 0., data->numEntries() );

    RooAddPdf* pdf = new RooAddPdf("pdf", "pdf", RooArgList(signal,signal_Bs, bkg_exp, *bkg_partReco), RooArgList(n_sig,n_sig_Bs, n_bkg, n_partReco_bkg));
    
    /// Generate simultaneous pdf out of prototype pdf 
    RooSimPdfBuilder* mgr = new RooSimPdfBuilder(*pdf) ;
    RooArgSet* config = mgr->createProtoBuildConfig() ;
    config->setStringValue("physModels","pdf") ;
    config->setStringValue("splitCats" ,"run KsCat") ;
    config->setStringValue("pdf", "run            : scale_mean, scale_sigma "
                               "KsCat :  exp_par "  
                               "run,KsCat : n_sig, n_sig_Bs, n_bkg, n_partReco_bkg ") ; 
    
    RooSimultaneous* simPdf  = mgr->buildPdf(*config,data) ;
    simPdf->Print("v") ;
    RooArgSet* fitParams = simPdf->getParameters(data);
    fitParams->Print("v") ;

    for(int i=0; i<str_run.size(); i++) for(int j=0; j<str_KsCat.size(); j++){
        /// Set start values for yields
        RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"run==run::" + str_run[i] + " && KsCat == KsCat::" + str_KsCat[j]);
        ((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_KsCat[j] + "}"))->setVal(data_slice->numEntries()*0.15);
        ((RooRealVar*) fitParams->find("n_sig_Bs_{"+ str_run[i] + ";" + str_KsCat[j] + "}"))->setVal(data_slice->numEntries()*0.1);
        ((RooRealVar*) fitParams->find("n_bkg_{"+ str_run[i] + ";" + str_KsCat[j] + "}"))->setVal(data_slice->numEntries()*0.7);
        ((RooRealVar*) fitParams->find("n_partReco_bkg_{"+ str_run[i] + ";" + str_KsCat[j] + "}"))->setVal(data_slice->numEntries()*0.05);
    }

    /// Perform fit
	RooFitResult* result = simPdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(numCPU));
	cout << "result is --------------- "<<endl;
	result->Print("v");

	/// Plot combined data and fit
	TCanvas* c = new TCanvas();
	RooPlot* frame= DTF_B_M.frame();
	frame->SetTitle("");
	data->plotOn(frame,Name("data"),MarkerSize(1),Binning(nBins));
	//simPdf->plotOn(frame,Name("pdf"),ProjWData(run,*data),LineColor(kBlue+1),LineWidth(3));
    
    TLegend leg(0.6,0.4,0.9,0.9,"");
    leg.SetLineStyle(0);
    leg.SetLineColor(0);
    leg.SetFillColor(0);
    leg.SetTextFont(22);
    leg.SetTextColor(1);
    leg.SetTextSize(0.05);
    leg.SetTextAlign(12);
    leg.AddEntry(frame->findObject("data"),"LHCb data","ep");
    leg.AddEntry(frame->findObject("pdf"),"Fit","l");
    
    /// Plot components 
    /// argh fuck RooFit, is this seriously so complicated ? 
    TString last_name_signal,last_name_signal_Bs,last_name_exp_bkg,last_name_partReco_bkg;
    
    for(int i=0; i<str_run.size(); i++) for(int j=0; j<str_KsCat.size(); j++){
            RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_run[i] + ";" + str_KsCat[j] + "}");
            RooArgList pdf_slice_comp( pdf_slice->pdfList());
            
            TString name_signal("signal_"+ anythingToString(i)+ "_" + anythingToString(j));
            TString name_signal_Bs("signal_Bs_"+ anythingToString(i)+ "_" + anythingToString(j));
            TString name_exp_bkg("bkg_"+ anythingToString(i)+ "_" + anythingToString(j));
            TString name_partReco_bkg("partReco_bkg_"+ anythingToString(i)+ "_" + anythingToString(j));
            
            if(i==0 && j == 0 ){
                pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
                
                pdf_slice->plotOn(frame,Name(name_signal_Bs),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
                
                pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
                
                pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
                
            }
            else{
                pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal));
                
                pdf_slice->plotOn(frame,Name(name_signal_Bs),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal_Bs));
                
                pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_exp_bkg));
                
                pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_partReco_bkg));
                
            }
            
            if(i== str_run.size()-1 && j == str_KsCat.size() -1 ) {
                
                pdf_slice->plotOn(frame,Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),FillColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg),DrawOption("F"),FillStyle(1001),LineColor(kGray+3));            
                pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),LineColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg));
                
                pdf_slice->plotOn(frame,Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),FillColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal),DrawOption("F"),FillStyle(3353),LineColor(kRed+1));
                pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),LineColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal));
                
                pdf_slice->plotOn(frame,Name(name_signal_Bs),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),FillColor(kGreen+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_Bs),DrawOption("F"),FillStyle(3335),LineColor(kGreen+3));
                pdf_slice->plotOn(frame,Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),LineColor(kGreen+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_Bs));
                
                pdf_slice->plotOn(frame,Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),LineColor(kBlack),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_exp_bkg),LineStyle(kDashed));
                
                
                leg.AddEntry(frame->findObject(name_signal),"#font[132]{B^{0}#rightarrowD^{-}K_{s}^{0}#pi^{+}}","f");
                leg.AddEntry(frame->findObject(name_signal_Bs),"#font[132]{B_{s}#rightarrowD^{-}K_{s}^{0}#pi^{+}}","f");
                leg.AddEntry(frame->findObject(name_exp_bkg),"Comb. bkg.","l");
                leg.AddEntry(frame->findObject(name_partReco_bkg),"Part. reco. bkg.","f");
            }
            else {
                last_name_signal = name_signal;
                last_name_signal_Bs = name_signal_Bs;
                last_name_exp_bkg = name_exp_bkg;
                last_name_partReco_bkg = name_partReco_bkg;
            }
    }
    
    simPdf->plotOn(frame,Name("pdf"),ProjWData(run,*data),LineColor(kBlue+1),LineWidth(3));  //,VisualizeError(*result,1,kTRUE)
    frame->Draw();
    leg.Draw();
    c->Print("plots/signal.eps");

    double chi2 = 0.;
    double covmatr = result->covQual();
    double edm = result->edm();
    RooHist* hpull  = frame->pullHist("data","pdf");
    frame= DTF_B_M.frame();
    frame->addPlotable(hpull,"P") ;
    frame->Draw();
    c->Print("plots/signal_pull.eps");
    
	/// Output file
	TFile *output;
	TTree* out_tree;
	double sw,sw_Bs;
	TBranch *b_sw, *b_w;
    Int_t t_run, t_KsCat, t_TriggerCat;
    TBranch *b_run, *b_KsCat, *b_TriggerCat;

	if(sWeight){
        //SPlot sPlot("sPlot","sPlot",*data,pdf,RooArgList(n_sig,n_sig_Bs, n_bkg, n_partReco_bkg)); 
		output = new TFile(((string)outFileName).c_str(),"RECREATE");
		tree->SetBranchStatus("*",1);
		tree->SetBranchStatus("weight",0);
		tree->SetBranchStatus("N_B_sw",0);

        out_tree = tree->CopyTree(("B_DTF_MM >= " + anythingToString((double)min_MM) + " && B_DTF_MM <= " + anythingToString((double)max_MM) + " && " + (string)cut_BDT).c_str());
        b_sw = out_tree->Branch("N_B_sw", &sw, "N_B_sw/D");
        b_w = out_tree->Branch("weight", &sw_Bs, "weight/D");
        
        out_tree->SetBranchAddress("run", &t_run, &b_run);
        out_tree->SetBranchAddress("KsCat", &t_KsCat, &b_KsCat);
        out_tree->SetBranchAddress("TriggerCat", &t_TriggerCat, &b_TriggerCat);

        if(out_tree->GetEntries() != data->numEntries()) {
            cout << "ERROR:: Different number of events in input and outputfile ! " << endl;
            cout << out_tree->GetEntries() << endl;
            cout << data->numEntries() << endl;
            throw "ERROR";
        }
	}
    
    double weights[(int)data->numEntries()];
    double weights_Bs[(int)data->numEntries()];    
    /// Calculate total signal yield
    double signal_yield = 0.;
    double comb_bkg_yield = 0.;
    int n_sig_perFit = 0.;
    int n_sig_perFit_err = 0.;
    DTF_B_M.setRange("signal_range",mean.getVal()-45.,mean.getVal()+45.);
    
    /// Loop over pdf slices
    for(int i=0; i<str_run.size(); i++)for(int j=0; j<str_KsCat.size(); j++){
            
                TLatex* lhcbtext = new TLatex();
                lhcbtext->SetTextFont(22);
                lhcbtext->SetTextColor(1);
                lhcbtext->SetTextSize(0.07);
                lhcbtext->SetTextAlign(13);
                lhcbtext->SetNDC(1);
                
                /// Get pdf slice
                RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_run[i] + ";" + str_KsCat[j]+ "}");
                RooArgList pdf_slice_comp( pdf_slice->pdfList());
                pdf_slice->Print();
                /// Plot data and pdf slices
                frame=DTF_B_M.frame();
                data->plotOn(frame,Name("data_slice2"),Cut("run==run::" + str_run[i] + " && KsCat==KsCat::" + str_KsCat[j]),MarkerSize(1),Binning(nBins));
                pdf_slice->plotOn(frame,Name("pdf_slice2"),LineColor(kBlue+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected));
                
                pdf_slice->plotOn(frame,LineColor(kGray+3),Components(RooArgSet(bkg_partReco_Bs,bkg_partReco_B0)),Normalization(1.,RooAbsReal::RelativeExpected));
                pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(1001),FillColor(kGray+3),Components(RooArgSet(bkg_partReco_Bs,bkg_partReco_B0)),Normalization(1.,RooAbsReal::RelativeExpected));
           
                pdf_slice->plotOn(frame,LineColor(kRed+1),Components("signal_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
                pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3353),FillColor(kRed+1),Components("signal_{" + str_run[i] + ";" + str_KsCat[j]+ "}"),Normalization(1.,RooAbsReal::RelativeExpected));
                
                pdf_slice->plotOn(frame,LineColor(kGreen+3),Components("signal_Bs_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
                pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3335),FillColor(kGreen+3),Components("signal_Bs_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
                
                pdf_slice->plotOn(frame,LineStyle(kDashed),LineColor(kBlack),Components("bkg_exp_{" + str_run[i] + ";" + str_KsCat[j]+ "}"),Normalization(1.,RooAbsReal::RelativeExpected));
                //data->plotOn(frame,Name("data_slice2"),Cut("year==year::" + str_year[i] + " && Ds_finalState == Ds_finalState::" + str_Ds[j]),MarkerSize(1),Binning(nBins));
                frame->Draw();
                chi2 += frame->chiSquare("pdf_slice2","data_slice2");
                cout << endl << "chi2/nbin = " << frame->chiSquare("pdf_slice2","data_slice2") << endl << endl;    
                
                n_sig_perFit = ((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_KsCat[j] + "}"))->getVal();
                n_sig_perFit_err = ((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_KsCat[j]+ "}"))->getError();
                TString name_n_signal("N_{sig} =" + anythingToString(n_sig_perFit) + "#pm" + anythingToString(n_sig_perFit_err));
                
                /// Print label
                TString label = "LHCb " + str_run[i];
                label.ReplaceAll("1","I");
                label.ReplaceAll("2","II");
                lhcbtext->SetTextFont(22);
                lhcbtext->DrawLatex(0.6,0.85,label.ReplaceAll("Run","Run-"));
                if(str_KsCat[j]=="LL")label = "LL";
                else label = "DD";
                lhcbtext->SetTextFont(132);
                lhcbtext->DrawLatex(0.6,0.78,label);
                //if(str_trigger[k]=="t0")lhcbtext->DrawLatex(0.6,0.68,"L0-TOS");
                //else lhcbtext->DrawLatex(0.6,0.68,"L0-TIS");
                lhcbtext->DrawLatex(0.6,0.68, name_n_signal);
                
                c->Print("plots/signal_" + str_run[i] + "_" + str_KsCat[j] + ".eps");
                hpull = frame->pullHist("data_slice2","pdf_slice2") ;
                frame= DTF_B_M.frame();
                frame->addPlotable(hpull,"P") ;
                frame->Draw();
                c->Print("plots/signal_pull_" + str_run[i] + "_" + str_KsCat[j]  + ".eps");
                
                /// Get signal yield
                signal_yield += ((RooRealVar*) fitParams->find("n_sig_{" + str_run[i] + ";" + str_KsCat[j] + "}"))->getVal();
                RooAbsPdf* pdf_slice_comb_bkg = (RooAbsPdf*) pdf_slice_comp.find("bkg_exp_{" + str_run[i] + ";" + str_KsCat[j]+ "}");
                comb_bkg_yield += pdf_slice_comb_bkg->createIntegral(DTF_B_M,NormSet(DTF_B_M),Range("signal_range"))->getVal()*((RooRealVar*) fitParams->find("n_bkg_{" + str_run[i] + ";" + str_KsCat[j]+  "}"))->getVal();
                
                /// Calculate sWeights
                if(sWeight){
                    RooArgList yield_list(
                                          *((RooRealVar*) fitParams->find("n_sig_{" + str_run[i] + ";" + str_KsCat[j]+ "}")),
                                          *((RooRealVar*) fitParams->find("n_sig_Bs_{" + str_run[i] + ";" + str_KsCat[j]+ "}")),
                                          *((RooRealVar*) fitParams->find("n_bkg_{" + str_run[i] + ";" + str_KsCat[j]+ "}")),
                                          *((RooRealVar*) fitParams->find("n_partReco_bkg_{" + str_run[i] + ";" + str_KsCat[j]+  "}"))
                                          );
                    RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"run==run::" + str_run[i] + " && KsCat==KsCat::" + str_KsCat[j]);
                    SPlot sPlot("sPlot","sPlot",*data_slice,pdf_slice,yield_list); 
                    
                    /// Plot the sWeight distributions as a function of mass
                    TH2 * swHist = (TH2*)data_slice->createHistogram("B_DTF_MM,n_sig_{" + str_run[i] + ";" + str_KsCat[j] + "}" + "_sw");
                    swHist->GetYaxis()->SetTitle("Signal sWeights");
                    swHist->Draw();
                    c->Print("plots/signal_sweight_" + str_run[i] + "_" + str_KsCat[j] + ".eps");
                    
                    /// Save sWeights
                    /// Messy and dangerous hack but works for now
                    int n_ij = 0;  /// labels entry number of data slice
                    for(int n = 0; n < out_tree->GetEntries(); n++){
                        b_run->GetEntry(n);
                        b_KsCat->GetEntry(n);
                        b_TriggerCat->GetEntry(n);
                        if(t_run == i+1 && t_KsCat == j){
                            weights[n] = sPlot.GetSWeight(n_ij,"n_sig_{" + str_run[i] + ";" + str_KsCat[j] + "}" + "_sw");
                            weights_Bs[n] = sPlot.GetSWeight(n_ij,"n_sig_Bs_{" + str_run[i] + ";" + str_KsCat[j]+ "}" + "_sw");
                            n_ij++;
                        }
                    }    
                }
    }
     
    if(sWeight){
        for(int n = 0; n < out_tree->GetEntries(); n++){
            sw = weights[n];
            sw_Bs = weights_Bs[n];
            b_sw->Fill();
            b_w->Fill();
        }
        out_tree->Write();
        output->Close();
        cout << endl;
        cout << "Created file " << ((string)outFileName).c_str()  << endl << endl;
    }
    chi2 = chi2*nBins/(str_run.size()*str_KsCat.size()*nBins-fitParams->selectByAttrib("Constant",kFALSE)->getSize());
    cout << endl; cout<<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl<<endl;
    cout << "Total signal yield = " << signal_yield << endl;
}

void fitSignal2D(){
    
    ///Options
    NamedParameter<int> numCPU("numCPU", 1);
    NamedParameter<int> sWeight("sWeightSignal", 0);
    NamedParameter<int> sWeightNew("sWeightSignalNew", 0);

    NamedParameter<int> nBins("nBins", 80);
    NamedParameter<double> min_MM("min_MM",5100.);
    NamedParameter<double> max_MM("max_MM",5700.);

    NamedParameter<double> signalRange_B_min_MM("signalRange_B_min_MM",5100.);
    NamedParameter<double> signalRange_B_max_MM("signalRange_B_max_MM",5700.);
    NamedParameter<double> signalRange_Ks_min_MM("signalRange_Ks_min_MM",475);
    NamedParameter<double> signalRange_Ks_max_MM("signalRange_Ks_max_MM",525);

    NamedParameter<string> cut_BDT("cut_BDT",(string)"");
    NamedParameter<string> inFileName("inFileNameSignal",(string)"../../../../../Selection/BDT/signal_data.root");
    NamedParameter<string> outFileName("outFileNameSignal",(string)"./signal_data.root");
    NamedParameter<int> fixSignalShapeFromMC("fixSignalShapeFromMC", 1);
    NamedParameter<int> fixExpBkgFromSidebands("fixExpBkgFromSidebands", 0);
    NamedParameter<string> initialState("initialState",(string)"Bd");
    NamedParameter<string> finalState("finalState",(string)"DKspi");
    TString outDir = (string)initialState + "2" + (string)finalState + "/";
    ///Load file
    TFile *file= new TFile(((string)inFileName).c_str());
    TTree* tree = (TTree*) file->Get("DecayTree");    
    tree->SetBranchStatus("*",0);
    tree->SetBranchStatus("B_DTF_MM",1);
    tree->SetBranchStatus("BDTG",1);
    tree->SetBranchStatus("D_PV_MM",1);
    tree->SetBranchStatus("Ks_PV_MM",1);
    tree->SetBranchStatus("pi_ProbNNpi",1);
    tree->SetBranchStatus("pi_ProbNNk",1);
    tree->SetBranchStatus("KsCat",1);
    tree->SetBranchStatus("run",1);
    tree->SetBranchStatus("TriggerCat",1);
    tree->SetBranchStatus("*m_*",1);
    
    RooRealVar DTF_B_M("B_DTF_MM", "m(D^{-}K_{S}#pi^{+})", min_MM, max_MM,"MeV/c^{2}");    
    RooRealVar BDTG("BDTG", "BDTG", 0.);    
    RooRealVar D_PV_MM("D_PV_MM", "D_PV_MM", 0.);    
    RooRealVar Ks_PV_MM("Ks_PV_MM", "Ks_PV_MM", 475 ,525);    
    RooRealVar pi_ProbNNpi("pi_ProbNNpi", "pi_ProbNNpi", 0.);    
    RooRealVar pi_ProbNNk("pi_ProbNNk", "pi_ProbNNk", 0.);    

    RooRealVar m_Dpi("m_Dpi", "m_Dpi", 0.);    
    RooRealVar m_Kspi("m_Kspi", "m_Kspi", 0.);    
    RooRealVar m_DKs("m_DKs", "m_DKs", 0.);    

    static const double massKs = 497.611;
    static const double massPi = ((string)finalState == "DKspi") ? 139.57 : 493.68;
    static const double massD = 1869.65;
    static const double massB = ((string)initialState == "Bd") ? 5279.65 : 5366.88;
    TString prefix_DTF = ((string)initialState == "Bd") ? "FullDTF_" : "FullBsDTF_";
    RooRealVar FullDTF_m_Dpi(prefix_DTF+"m_Dpi", prefix_DTF+"m_Dpi", massD+massPi, massB - massKs );    
    RooRealVar FullDTF_m_Kspi(prefix_DTF+"m_Kspi", prefix_DTF+"m_Kspi", massKs+massPi, massB - massD);    
    RooRealVar FullDTF_m_DKs(prefix_DTF+"m_DKs", prefix_DTF+"m_DKs", massD+massKs, massB - massPi);    
    TString cutPhsp = "( ";
    cutPhsp +=  prefix_DTF+"m_Dpi >=" +  anythingToString(massD+massPi) + " && " + prefix_DTF+"m_Dpi <=" + anythingToString(massB - massKs);
    cutPhsp +=  " && " + prefix_DTF+"m_Kspi >=" +  anythingToString(massKs+massPi) + " && " + prefix_DTF+"m_Kspi <=" +  anythingToString(massB - massD);
    cutPhsp +=  " && " + prefix_DTF+"m_DKs >=" +  anythingToString(massD+massKs) + " && " + prefix_DTF+"m_DKs <=" + anythingToString(massB - massPi);    
    cutPhsp += " ) ";
    cout << cutPhsp << endl;
    
    ///define categories
    RooCategory run("run","run") ;
    run.defineType("Run1",1) ;
    run.defineType("Run2",2) ;
    RooCategory TriggerCat("TriggerCat","TriggerCat") ;
    TriggerCat.defineType("t0",0);
    TriggerCat.defineType("t1",1);
    RooCategory KsCat("KsCat","KsCat") ;
    KsCat.defineType("LL",0);
    KsCat.defineType("DD",1);
    
    RooArgList list =  RooArgList(DTF_B_M,BDTG,D_PV_MM,Ks_PV_MM,pi_ProbNNk,pi_ProbNNpi);
    RooArgList list2 =  RooArgList(KsCat,run,TriggerCat);
    RooArgList list3 = RooArgList(m_Kspi,FullDTF_m_Kspi);
    RooArgList list4 = RooArgList(m_Dpi,FullDTF_m_Dpi);
    RooArgList list5 = RooArgList(m_DKs,FullDTF_m_DKs);
    list.add(list2);
    list.add(list3);
    list.add(list4);
    list.add(list5);
    RooDataSet*  data = new RooDataSet("data","data",tree,list,((string)cut_BDT).c_str());    
    
    /// Signal Pdf
    vector<double> sig_params = fitSignalShape("signal");
    RooRealVar mean_MC("mean_MC", "#mu MC", sig_params[0]); 
    RooRealVar sigma_MC("sigma_MC", "#sigma MC", sig_params[1]);
    RooRealVar sigma2_MC("sigma2_MC", "#sigma2 MC", sig_params[4]);
    RooRealVar scale_mean("scale_mean", "scale #mu",1.); 
    RooRealVar scale_sigma("scale_sigma", "scale #sigma", 1.,0.5,2);
    
    RooFormulaVar mean("mean","@0 * @1", RooArgSet(scale_mean,mean_MC)); 
    RooFormulaVar sigma("sigma","@0 * @1", RooArgSet(scale_sigma,sigma_MC)); 
    RooFormulaVar sigma2("sigma2","@0 * @1", RooArgSet(scale_sigma,sigma2_MC)); 
    
    RooRealVar alpha("alpha", "#alpha", sig_params[2],-5.,5.); 
    RooRealVar beta("beta", "#beta", sig_params[3],-5.,5.);
    RooRealVar alpha2("alpha2", "#alpha2", sig_params[5],-5.,5.); 
    RooRealVar beta2("beta2", "#beta2", sig_params[6],-5.,5.);
    RooRealVar f_RJ("f_RJ", "f_RJ", sig_params[7]);
    
    if(fixSignalShapeFromMC){
        alpha.setConstant();
        alpha2.setConstant();
        beta.setConstant();
        beta2.setConstant();
    } 
    RooJohnsonSU signal_B("signal_B","signal_B",DTF_B_M, mean,sigma,alpha,beta);

    //Ks
    vector<double> sig_params_Ks = fitSignalShape("Ks");
    RooRealVar mean_MC_Ks("mean_MC_Ks", "mean_MC_Ks", sig_params_Ks[0]); 
    RooRealVar sigma_MC_Ks("sigma_MC_Ks", "sigma_MC_Ks", sig_params_Ks[1]);
    RooRealVar scale_mean_Ks("scale_mean_Ks", "scale_mean_Ks",1.); 
    RooRealVar scale_sigma_Ks("scale_sigma_Ks", "scale_sigma_Ks", 1.,0.5,2);
    RooFormulaVar mean_Ks("mean_Ks","@0 * @1", RooArgSet(scale_mean_Ks,mean_MC_Ks)); 
    RooFormulaVar sigma_Ks("sigma_Ks","@0 * @1", RooArgSet(scale_sigma_Ks,sigma_MC_Ks));     
    RooRealVar alpha_Ks("alpha_Ks", "alpha_Ks", sig_params_Ks[2],-5.,5.); 
    RooRealVar beta_Ks("beta_Ks", "beta_Ks", sig_params_Ks[3],-5.,5.);    
    if(fixSignalShapeFromMC){
        alpha_Ks.setConstant();
        beta_Ks.setConstant();
    } 
    RooJohnsonSU signal_Ks("signal_Ks","signal_Ks",Ks_PV_MM, mean_Ks,sigma_Ks,alpha_Ks,beta_Ks);
    
    /// Bs pdf
    RooFormulaVar mean_Bs("mean_B0","@0 + @1", RooArgSet(mean,RooConst(87.33))); 
    RooJohnsonSU signal_Bs("signal_Bs","signal_Bs",DTF_B_M, mean_Bs,sigma,alpha,beta);

    /// Combinatorial bkg pdf
    RooRealVar exp_par("exp_par","exp_par",-1.6508e-03,-10.,10.);
    RooExponential bkg_exp("bkg_exp","bkg_exp",DTF_B_M,exp_par);
    bkg_exp.fitTo(*data,Save(kTRUE),Range(5450.,max_MM));
    if(fixExpBkgFromSidebands)exp_par.setConstant();
    
    RooRealVar c0_Ks("c0_Ks", "c0_Ks", .0,-10,10); 
    RooRealVar c1_Ks("c1_Ks", "c1_Ks", .0,-10,10); 
    RooRealVar c2_Ks("c2_Ks", "c2_Ks", .0,-10,10); 
    RooChebychev bkg_Ks("bkg_Ks","bkg_Ks",Ks_PV_MM, RooArgList(c0_Ks));
    
    /// Part. reco bkg
    vector<double> bkg_partReco_params = fitPartRecoBkgShape();
    RooRealVar mean1("mean1","mean1", bkg_partReco_params[0]);
    RooRealVar mean2("mean2","mean2", bkg_partReco_params[1]);
    RooRealVar mean3("mean3","mean3", bkg_partReco_params[2]);
    RooRealVar mean1Shifted("mean1Shifted","mean1Shifted", mean1.getVal() - 87.33 );
    RooRealVar mean2Shifted("mean2Shifted","mean2Shifted", mean2.getVal() - 87.33 );
    RooRealVar mean3Shifted("mean3Shifted","mean3Shifted", mean3.getVal() - 87.33 );
    RooRealVar sigmaL1("sigmaL1", "sigmaL1",  bkg_partReco_params[3]);
    RooRealVar sigmaR1("sigmaR1", "sigmaR1",  bkg_partReco_params[4]);
    RooRealVar sigmaL2("sigmaR1", "sigmaL2",  bkg_partReco_params[5]);
    RooRealVar sigmaR2("sigmaR2", "sigmaR2",  bkg_partReco_params[6]);
    RooRealVar sigmaL3("sigmaL3", "sigmaL3",  bkg_partReco_params[7]);
    RooRealVar sigmaR3("sigmaR3", "sigmaR3",  bkg_partReco_params[8]);
    RooRealVar f_1("f_1", "f_1", bkg_partReco_params[9]);
    RooRealVar f_2("f_2", "f_2", bkg_partReco_params[10]);
    
    RooBifurGauss BifGauss1_Bs("BifGauss1_Bs","BifGauss1_Bs", DTF_B_M, mean1, sigmaL1,sigmaR1);
    RooBifurGauss BifGauss2_Bs("BifGauss2_Bs","BifGauss2_Bs", DTF_B_M, mean2, sigmaL2,sigmaR2);
    RooBifurGauss BifGauss3_Bs("BifGauss3_Bs","BifGauss3_Bs", DTF_B_M, mean3, sigmaL3,sigmaR3);
    RooAddPdf bkg_partReco_Bs("bkg_partReco_Bs", "bkg_partReco_Bs", RooArgList(BifGauss1_Bs, BifGauss2_Bs, BifGauss3_Bs), RooArgList(f_1,f_2),kTRUE);
    
    RooBifurGauss BifGauss1_B0("BifGauss1_B0","BifGauss1_B0", DTF_B_M, mean1Shifted, sigmaL1,sigmaR1);
    RooBifurGauss BifGauss2_B0("BifGauss2_B0","BifGauss2_B0", DTF_B_M, mean2Shifted, sigmaL2,sigmaR2);
    RooBifurGauss BifGauss3_B0("BifGauss3_B0","BifGauss3_B0", DTF_B_M, mean3Shifted, sigmaL3,sigmaR3);
    RooAddPdf bkg_partReco_B0("bkg_partReco_B0", "bkg_partReco_B0", RooArgList(BifGauss1_B0, BifGauss2_B0, BifGauss3_B0), RooArgList(f_1,f_2),kTRUE);
    
    RooRealVar partReco_f("partReco_f", "partReco_f", 0.5,0,1);
    RooAddPdf bkg_partReco_B("bkg_partReco_B", "bkg_partReco_B", RooArgList(bkg_partReco_Bs, bkg_partReco_B0), RooArgList(partReco_f));
    
    // 2D pdfs
    RooProdPdf signal2D("signal2D","signal2D",signal_B,signal_Ks);
    RooProdPdf signal_Bs2D("signal_Bs2D","signal_Bs2D",signal_Bs,signal_Ks);

    RooProdPdf signal_B_noKs("signal_B_noKs","signal_B_noKs",signal_B,bkg_Ks);
    RooProdPdf signal_Bs_noKs("signal_Bs_noKs","signal_Bs_noKs",signal_Bs,bkg_Ks);
    
    RooRealVar f_comb_Ks("f_comb_Ks", "f_comb_Ks", 0.5,0,1);
    RooAddPdf bkg_comb_Ks("bkg_comb_Ks", "bkg_comb_Ks", RooArgList(signal_Ks, bkg_Ks), RooArgList(f_comb_Ks));
    RooProdPdf bkg_comb2D("bkg_comb2D","bkg_comb2D",bkg_exp,bkg_comb_Ks);

    RooRealVar f_partReco_Ks("f_partReco_Ks", "f_partReco_Ks", 0.5,0,1);
    RooAddPdf bkg_partReco_Ks("bkg_partReco_Ks", "bkg_partReco_Ks", RooArgList(signal_Ks, bkg_Ks), RooArgList(f_partReco_Ks));
    RooProdPdf bkg_partReco2D("bkg_partReco2D","bkg_partReco2D",bkg_partReco_B,bkg_partReco_Ks);
    
    /// Total pdf
    RooRealVar n_sig("n_sig", "n_sig", data->numEntries()*0.15, 0., data->numEntries());
    RooRealVar n_sig_Bs("n_sig_Bs", "n_sig_B0", data->numEntries()*0.05, 0., data->numEntries());
    RooRealVar n_sig_noKs("n_sig_noKs", "n_sig_noKs", data->numEntries()*0.15, 0., data->numEntries());
    RooRealVar n_sig_Bs_noKs("n_sig_Bs_noKs", "n_sig_B0_noKs", data->numEntries()*0.05, 0., data->numEntries());
    RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()*0.7, 0., data->numEntries());
    RooRealVar n_partReco_bkg("n_partReco_bkg", "n_partReco_bkg", data->numEntries()*0.05, 0., data->numEntries() );

    RooAddPdf* pdf = new RooAddPdf("pdf", "pdf", RooArgList(signal2D,signal_Bs2D,bkg_comb2D, bkg_partReco2D, signal_B_noKs,signal_Bs_noKs), RooArgList(n_sig,n_sig_Bs, n_bkg, n_partReco_bkg, n_sig_noKs, n_sig_Bs_noKs));

    /// Generate simultaneous pdf out of prototype pdf 
    RooSimPdfBuilder* mgr = new RooSimPdfBuilder(*pdf) ;
    RooArgSet* config = mgr->createProtoBuildConfig() ;
    config->setStringValue("physModels","pdf") ;
    config->setStringValue("splitCats" ,"run KsCat") ;
    config->setStringValue("pdf", 
                           "KsCat :  scale_sigma_Ks, exp_par, c0_Ks, f_comb_Ks, f_partReco_Ks "  
                           "run,KsCat : scale_mean, scale_sigma, n_sig, n_sig_Bs, n_sig_noKs, n_sig_Bs_noKs, n_bkg, n_partReco_bkg ") ; 
    
    RooSimultaneous* simPdf  = mgr->buildPdf(*config,data) ;
    simPdf->Print("v") ;
    RooArgSet* fitParams = simPdf->getParameters(data);
    fitParams->Print("v") ;
    
    for(int i=0; i<str_run.size(); i++) for(int j=0; j<str_KsCat.size(); j++){
        /// Set start values for yields
        RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"run==run::" + str_run[i] + " && KsCat == KsCat::" + str_KsCat[j]);
        ((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_KsCat[j] + "}"))->setVal(data_slice->numEntries()*0.15);
        ((RooRealVar*) fitParams->find("n_sig_Bs_{"+ str_run[i] + ";" + str_KsCat[j] + "}"))->setVal(data_slice->numEntries()*0.1);
        ((RooRealVar*) fitParams->find("n_sig_noKs_{"+ str_run[i] + ";" + str_KsCat[j] + "}"))->setVal(data_slice->numEntries()*0.15);
        ((RooRealVar*) fitParams->find("n_sig_Bs_noKs_{"+ str_run[i] + ";" + str_KsCat[j] + "}"))->setVal(data_slice->numEntries()*0.1);
        ((RooRealVar*) fitParams->find("n_bkg_{"+ str_run[i] + ";" + str_KsCat[j] + "}"))->setVal(data_slice->numEntries()*0.7);
        ((RooRealVar*) fitParams->find("n_partReco_bkg_{"+ str_run[i] + ";" + str_KsCat[j] + "}"))->setVal(data_slice->numEntries()*0.05);
        
        ((RooRealVar*) fitParams->find("c0_Ks_DD"))->setVal(0.);
        ((RooRealVar*) fitParams->find("c0_Ks_DD"))->setConstant();
    }
    
    /// Perform fit
    RooFitResult* result = simPdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(numCPU));
    cout << "result is --------------- "<<endl;
    result->Print("v");
    
    /// Plot combined data and fit
    TCanvas* c = new TCanvas();
    RooPlot* frame= DTF_B_M.frame();
    frame->SetTitle("");
    data->plotOn(frame,Name("data"),MarkerSize(1),Binning(nBins));
    RooPlot* frame_Ks= Ks_PV_MM.frame();
    frame_Ks->SetTitle("");
    data->plotOn(frame_Ks,Name("data"),MarkerSize(1),Binning(nBins));
    
    TLegend leg(0.6,0.4,0.9,0.9,"");
    leg.SetLineStyle(0);
    leg.SetLineColor(0);
    leg.SetFillColor(0);
    leg.SetTextFont(22);
    leg.SetTextColor(1);
    leg.SetTextSize(0.05);
    leg.SetTextAlign(12);
    leg.AddEntry(frame->findObject("data"),"LHCb data","ep");
    leg.AddEntry(frame->findObject("pdf"),"Fit","l");
    
    /// Plot components 
    /// argh fuck RooFit, is this seriously so complicated ? 
    TString last_name_pdf,last_name_signal,last_name_signal_Bs,last_name_exp_bkg,last_name_partReco_bkg,last_name_signal_noKs,last_name_signal_Bs_noKs;
    vector<RooPlot*> frames;
    frames.push_back(frame);
    frames.push_back(frame_Ks);
    
    for(int k=0; k<frames.size(); k++)for(int i=0; i<str_run.size(); i++) for(int j=0; j<str_KsCat.size(); j++){
        RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_run[i] + ";" + str_KsCat[j] + "}");
        RooArgList pdf_slice_comp( pdf_slice->pdfList());
        
        TString name_signal("signal_"+ anythingToString(i)+ "_" + anythingToString(j));
        TString name_signal_Bs("signal_Bs_"+ anythingToString(i)+ "_" + anythingToString(j));
        TString name_exp_bkg("bkg_"+ anythingToString(i)+ "_" + anythingToString(j));
        TString name_partReco_bkg("partReco_bkg_"+ anythingToString(i)+ "_" + anythingToString(j));
        TString name_signal_noKs("signal_noKs_"+ anythingToString(i)+ "_" + anythingToString(j));
        TString name_signal_Bs_noKs("signal_Bs_noKs"+ anythingToString(i)+ "_" + anythingToString(j));
        
        if(i==0 && j == 0 ){
            pdf_slice->plotOn(frames[k],Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
            
            pdf_slice->plotOn(frames[k],Name(name_signal_Bs),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
            
            pdf_slice->plotOn(frames[k],Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
            
            pdf_slice->plotOn(frames[k],Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
            
            pdf_slice->plotOn(frames[k],Name(name_signal_noKs),Components(RooArgSet(pdf_slice_comp[4])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
            
            pdf_slice->plotOn(frames[k],Name(name_signal_Bs_noKs),Components(RooArgSet(pdf_slice_comp[5])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());

        }
        else{
            pdf_slice->plotOn(frames[k],Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal));
            
            pdf_slice->plotOn(frames[k],Name(name_signal_Bs),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal_Bs));
            
            pdf_slice->plotOn(frames[k],Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_exp_bkg));
            
            pdf_slice->plotOn(frames[k],Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_partReco_bkg));
            pdf_slice->plotOn(frames[k],Name(name_signal_noKs),Components(RooArgSet(pdf_slice_comp[4])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal_noKs));
            
            pdf_slice->plotOn(frames[k],Name(name_signal_Bs_noKs),Components(RooArgSet(pdf_slice_comp[5])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible(),AddTo(last_name_signal_Bs_noKs));
            
        }
        if(i== str_run.size()-1 && j == str_KsCat.size() -1 ) {
            
            pdf_slice->plotOn(frames[k],Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),FillColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg),DrawOption("F"),FillStyle(1001),LineColor(kGray+3));            
            pdf_slice->plotOn(frames[k],Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),LineColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg));
            
            pdf_slice->plotOn(frames[k],Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),FillColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal),DrawOption("F"),FillStyle(3353),LineColor(kRed+1));
            pdf_slice->plotOn(frames[k],Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),LineColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal));
            
            pdf_slice->plotOn(frames[k],Name(name_signal_Bs),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),FillColor(kGreen+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_Bs),DrawOption("F"),FillStyle(3335),LineColor(kGreen+1));
            pdf_slice->plotOn(frames[k],Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),LineColor(kGreen+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_Bs));
            
            pdf_slice->plotOn(frames[k],Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),LineColor(kBlack),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_exp_bkg),LineStyle(kDashed));
            
            pdf_slice->plotOn(frames[k],Name(name_signal_noKs),Components(RooArgSet(pdf_slice_comp[4])),ProjWData(run,*data),LineColor(kRed+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_noKs),LineStyle(kDashed));

            pdf_slice->plotOn(frames[k],Name(name_signal_Bs_noKs),Components(RooArgSet(pdf_slice_comp[5])),ProjWData(run,*data),LineColor(kGreen+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_Bs_noKs),LineStyle(kDashed));
            
            if(k==0){
                leg.AddEntry(frames[k]->findObject(name_signal),"#font[132]{B^{0}#rightarrowD^{-}K_{s}^{0}#pi^{+}}","f");
                leg.AddEntry(frames[k]->findObject(name_signal_Bs),"#font[132]{B_{s}#rightarrowD^{-}K_{s}^{0}#pi^{+}}","f");
                leg.AddEntry(frames[k]->findObject(name_signal_noKs),"#font[132]{B^{0}#rightarrowD^{-}#pi^{+}#pi^{-}#pi^{+}}","l");
                leg.AddEntry(frames[k]->findObject(name_signal_Bs_noKs),"#font[132]{B_{s}#rightarrowD^{-}#pi^{+}#pi^{-}#pi^{+}}","l");
                leg.AddEntry(frames[k]->findObject(name_exp_bkg),"Comb. bkg.","l");
                leg.AddEntry(frames[k]->findObject(name_partReco_bkg),"Part. reco. bkg.","f");
            }
        }
        else {
            last_name_signal = name_signal;
            last_name_signal_Bs = name_signal_Bs;
            last_name_exp_bkg = name_exp_bkg;
            last_name_partReco_bkg = name_partReco_bkg;
            last_name_signal_noKs = name_signal_noKs;
            last_name_signal_Bs_noKs = name_signal_Bs_noKs;
        }
    }
    simPdf->plotOn(frame,Name("pdf"),ProjWData(run,*data),LineColor(kBlue+1),LineWidth(3));  //,VisualizeError(*result,1,kTRUE)
    frame->Draw();
    leg.Draw();
    c->Print(outDir+"signal.eps");
    
    double chi2 = 0.;
    double covmatr = result->covQual();
    double edm = result->edm();
    RooHist* hpull  = frame->pullHist("data","pdf");
    frame= DTF_B_M.frame();
    frame->addPlotable(hpull,"P") ;
    frame->Draw();
    //c->Print(outDir+"signal_pull.eps");

    simPdf->plotOn(frame_Ks,Name("pdf"),ProjWData(run,*data),LineColor(kBlue+1),LineWidth(3));  //,VisualizeError(*result,1,kTRUE)
    frame_Ks->Draw();
    //leg.Draw();
    c->Print(outDir+"Ks.eps");
    hpull  = frame_Ks->pullHist("data","pdf");
    frame_Ks= Ks_PV_MM.frame();
    frame_Ks->addPlotable(hpull,"P") ;
    frame_Ks->Draw();
    //c->Print(outDir+"Ks_pull.eps");


    // Plots integrated over Run
    for(int j=0; j<str_KsCat.size(); j++)for(int k=0; k<frames.size(); k++){
            if(k==0)frames[k] =DTF_B_M.frame();
            else frames[k] =Ks_PV_MM.frame();

            data->plotOn(frames[k],Name("data_slice_Ks"),Cut("KsCat==KsCat::" + str_KsCat[j]),MarkerSize(1),Binning(nBins));

            for(int i=0; i<str_run.size(); i++){

                RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_run[i] + ";" + str_KsCat[j] + "}");
                RooArgList pdf_slice_comp( pdf_slice->pdfList());

                TString name_pdf("pdf_"+ anythingToString(i)+ "_" + anythingToString(j));
                TString name_signal("signal_"+ anythingToString(i)+ "_" + anythingToString(j));
                TString name_signal_Bs("signal_Bs_"+ anythingToString(i)+ "_" + anythingToString(j));
                TString name_exp_bkg("bkg_"+ anythingToString(i)+ "_" + anythingToString(j));
                TString name_partReco_bkg("partReco_bkg_"+ anythingToString(i)+ "_" + anythingToString(j));
                TString name_signal_noKs("signal_noKs_"+ anythingToString(i)+ "_" + anythingToString(j));
                TString name_signal_Bs_noKs("signal_Bs_noKs"+ anythingToString(i)+ "_" + anythingToString(j));
                
                if(i==0){
                    pdf_slice->plotOn(frames[k],Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
                    
                    pdf_slice->plotOn(frames[k],Name(name_signal_Bs),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
                    
                    pdf_slice->plotOn(frames[k],Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
                    
                    pdf_slice->plotOn(frames[k],Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
                    
                    pdf_slice->plotOn(frames[k],Name(name_signal_noKs),Components(RooArgSet(pdf_slice_comp[4])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());
                    
                    pdf_slice->plotOn(frames[k],Name(name_signal_Bs_noKs),Components(RooArgSet(pdf_slice_comp[5])),ProjWData(run,*data),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());

                    pdf_slice->plotOn(frames[k],Name(name_pdf),LineColor(kBlue+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),Invisible());

                    last_name_pdf = name_pdf;
                    last_name_signal = name_signal;
                    last_name_signal_Bs = name_signal_Bs;
                    last_name_exp_bkg = name_exp_bkg;
                    last_name_partReco_bkg = name_partReco_bkg;
                    last_name_signal_noKs = name_signal_noKs;
                    last_name_signal_Bs_noKs = name_signal_Bs_noKs;
                }
                else{            
                    pdf_slice->plotOn(frames[k],Name(name_partReco_bkg),Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),FillColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg),DrawOption("F"),FillStyle(1001),LineColor(kGray+3));            
                    pdf_slice->plotOn(frames[k],Components(RooArgSet(pdf_slice_comp[3])),ProjWData(run,*data),LineColor(kGray+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_partReco_bkg));
                    
                    pdf_slice->plotOn(frames[k],Name(name_signal),Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),FillColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal),DrawOption("F"),FillStyle(3353),LineColor(kRed+1));
                    pdf_slice->plotOn(frames[k],Components(RooArgSet(pdf_slice_comp[0])),ProjWData(run,*data),LineColor(kRed+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal));
                    
                    pdf_slice->plotOn(frames[k],Name(name_signal_Bs),Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),FillColor(kGreen+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_Bs),DrawOption("F"),FillStyle(3335),LineColor(kGreen+1));
                    pdf_slice->plotOn(frames[k],Components(RooArgSet(pdf_slice_comp[1])),ProjWData(run,*data),LineColor(kGreen+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_Bs));
                    
                    pdf_slice->plotOn(frames[k],Name(name_exp_bkg),Components(RooArgSet(pdf_slice_comp[2])),ProjWData(run,*data),LineColor(kBlack),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_exp_bkg),LineStyle(kDashed));
                    
                    pdf_slice->plotOn(frames[k],Name(name_signal_noKs),Components(RooArgSet(pdf_slice_comp[4])),ProjWData(run,*data),LineColor(kRed+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_noKs),LineStyle(kDashed));

                    pdf_slice->plotOn(frames[k],Name(name_signal_Bs_noKs),Components(RooArgSet(pdf_slice_comp[5])),ProjWData(run,*data),LineColor(kGreen+3),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_signal_Bs_noKs),LineStyle(kDashed));
                    
                    pdf_slice->plotOn(frames[k],Name(name_pdf),LineColor(kBlue+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected),AddTo(last_name_pdf));
                }
            }
            frames[k]->Draw();
            if(k==0)leg.Draw();
            c->Print(outDir+"signal_" + anythingToString(k) + "_" + str_KsCat[j] + ".eps");
    }
    







    
    /// Output file
    TFile *output;
    TTree* out_tree;
    double sw,sw_Bs,sw_noKs,sw_Bs_noKs;
    TBranch *b_sw, *b_sw_Bs,*b_sw_noKs, *b_sw_Bs_noKs;
    double sw2,sw2_Bs,sw2_noKs,sw2_Bs_noKs;
    TBranch *b_sw2, *b_sw2_Bs,*b_sw2_noKs, *b_sw2_Bs_noKs;

    Int_t t_run, t_KsCat, t_TriggerCat;
    TBranch *b_run, *b_KsCat, *b_TriggerCat;
    
    if(sWeight || sWeightNew){
        output = new TFile(((string)outFileName).c_str(),"RECREATE");
        tree->SetBranchStatus("*",1);
        tree->SetBranchStatus("weight",0);
        tree->SetBranchStatus("N_B_sw",0);
        
        out_tree = tree->CopyTree("Ks_PV_MM >= 475 && Ks_PV_MM <= 525 && B_DTF_MM >= " + anythingToString((double)min_MM) + " && B_DTF_MM <= " + anythingToString((double)max_MM) + " && " + (string)cut_BDT + " && " + cutPhsp);
        b_sw = out_tree->Branch("N_B_sw", &sw, "N_B_sw/D");
        b_sw_Bs = out_tree->Branch("N_Bs_sw", &sw_Bs, "N_Bs_sw/D");
        b_sw_noKs = out_tree->Branch("N_B_noKs_sw", &sw_noKs, "N_B_noKs_sw/D");
        b_sw_Bs_noKs = out_tree->Branch("N_Bs_noKs_sw", &sw_Bs_noKs, "N_Bs_noKs_sw/D");

        b_sw2 = out_tree->Branch("N_B_sw2", &sw2, "N_B_sw2/D");
        b_sw2_Bs = out_tree->Branch("N_Bs_sw2", &sw2_Bs, "N_Bs_sw2/D");
        b_sw2_noKs = out_tree->Branch("N_B_noKs_sw2", &sw2_noKs, "N_B_noKs_sw2/D");
        b_sw2_Bs_noKs = out_tree->Branch("N_Bs_noKs_sw2", &sw2_Bs_noKs, "N_Bs_noKs_sw2/D");

        out_tree->SetBranchAddress("run", &t_run, &b_run);
        out_tree->SetBranchAddress("KsCat", &t_KsCat, &b_KsCat);
        out_tree->SetBranchAddress("TriggerCat", &t_TriggerCat, &b_TriggerCat);
        
        if(out_tree->GetEntries() != data->numEntries()) {
            cout << "ERROR:: Different number of events in input and outputfile ! " << endl;
            cout << out_tree->GetEntries() << endl;
            cout << data->numEntries() << endl;
            throw "ERROR";
        }
    }
    
    double weights[(int)data->numEntries()];
    double weights_Bs[(int)data->numEntries()];    
    double weights_noKs[(int)data->numEntries()];
    double weights_Bs_noKs[(int)data->numEntries()];    

    double weights2[(int)data->numEntries()];
    double weights2_Bs[(int)data->numEntries()];    
    double weights2_noKs[(int)data->numEntries()];
    double weights2_Bs_noKs[(int)data->numEntries()];    

    /// Calculate total signal yield
    double signal_yield = 0.;
    double comb_bkg_yield = 0.;
    int n_sig_perFit = 0.;
    int n_sig_perFit_err = 0.;
    DTF_B_M.setRange("signal_range",mean.getVal()-45,mean.getVal()+45);
    DTF_B_M.setRange("signal_rangeLL",(double)signalRange_B_min_MM,(double)signalRange_B_max_MM);
    Ks_PV_MM.setRange("signal_rangeLL",(double)signalRange_Ks_min_MM,(double)signalRange_Ks_max_MM);
    DTF_B_M.setRange("signal_rangeDD",(double)signalRange_B_min_MM,(double)signalRange_B_max_MM);
    Ks_PV_MM.setRange("signal_rangeDD ",Ks_PV_MM.getMin() ,Ks_PV_MM.getMax());

    /// Loop over pdf slices
    for(int i=0; i<str_run.size(); i++)for(int j=0; j<str_KsCat.size(); j++){
        
        TLatex* lhcbtext = new TLatex();
        lhcbtext->SetTextFont(22);
        lhcbtext->SetTextColor(1);
        lhcbtext->SetTextSize(0.07);
        lhcbtext->SetTextAlign(13);
        lhcbtext->SetNDC(1);
        
        /// Get pdf slice
        RooAddPdf* pdf_slice = (RooAddPdf*)simPdf->getPdf("{" + str_run[i] + ";" + str_KsCat[j]+ "}");
        RooArgList pdf_slice_comp( pdf_slice->pdfList());
        pdf_slice->Print();
        
        /// Plot data and pdf slices
        frame=DTF_B_M.frame();
        data->plotOn(frame,Name("data_slice2"),Cut("run==run::" + str_run[i] + " && KsCat==KsCat::" + str_KsCat[j]),MarkerSize(1),Binning(nBins));
        pdf_slice->plotOn(frame,Name("pdf_slice2"),LineColor(kBlue+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected));
        
        pdf_slice->plotOn(frame,LineColor(kGray+3),Components("bkg_partReco2D_{" + str_run[i] + ";" + str_KsCat[j]+ "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(1001),FillColor(kGray+3),Components("bkg_partReco2D_{" + str_run[i] + ";" + str_KsCat[j]+ "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        
        pdf_slice->plotOn(frame,LineColor(kRed+1),Components("signal2D_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3353),FillColor(kRed+1),Components("signal2D_{" + str_run[i] + ";" + str_KsCat[j]+ "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        
        pdf_slice->plotOn(frame,LineColor(kGreen+1),Components("signal_Bs2D_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3335),FillColor(kGreen+1),Components("signal_Bs2D_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        
        pdf_slice->plotOn(frame,LineColor(kRed+3),Components("signal_B_noKs_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected),LineStyle(kDashed));

        pdf_slice->plotOn(frame,LineColor(kGreen+3),Components("signal_Bs_noKs_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected),LineStyle(kDashed));

        pdf_slice->plotOn(frame,LineStyle(kDashed),LineColor(kBlack),Components("bkg_comb2D_{" + str_run[i] + ";" + str_KsCat[j]+ "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        frame->Draw();
        
        chi2 += frame->chiSquare("pdf_slice2","data_slice2");
        cout << endl << "chi2/nbin = " << frame->chiSquare("pdf_slice2","data_slice2") << endl << endl;    
        n_sig_perFit = ((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_KsCat[j] + "}"))->getVal();
        n_sig_perFit_err = ((RooRealVar*) fitParams->find("n_sig_{"+ str_run[i] + ";" + str_KsCat[j]+ "}"))->getError();
        TString name_n_signal("N_{sig} =" + anythingToString(n_sig_perFit) + "#pm" + anythingToString(n_sig_perFit_err));
        
        /// Print label
        TString label = "LHCb " + str_run[i];
        label.ReplaceAll("1","I");
        label.ReplaceAll("2","II");
        lhcbtext->SetTextFont(22);
        lhcbtext->DrawLatex(0.6,0.85,label.ReplaceAll("Run","Run-"));
        if(str_KsCat[j]=="LL")label = "LL";
        else label = "DD";
        lhcbtext->SetTextFont(132);
        lhcbtext->DrawLatex(0.6,0.78,label);
        //if(str_trigger[k]=="t0")lhcbtext->DrawLatex(0.6,0.68,"L0-TOS");
        //else lhcbtext->DrawLatex(0.6,0.68,"L0-TIS");
        lhcbtext->DrawLatex(0.6,0.68, name_n_signal);
        
        c->Print(outDir+"signal_" + str_run[i] + "_" + str_KsCat[j] + ".eps");
        hpull = frame->pullHist("data_slice2","pdf_slice2") ;
        frame= DTF_B_M.frame();
        frame->addPlotable(hpull,"P") ;
        frame->Draw();
        //c->Print(outDir+"signal_pull_" + str_run[i] + "_" + str_KsCat[j]  + ".eps");
        
        /// Ks frame
        frame=Ks_PV_MM.frame();
        data->plotOn(frame,Name("data_slice2"),Cut("run==run::" + str_run[i] + " && KsCat==KsCat::" + str_KsCat[j]),MarkerSize(1),Binning(nBins));
        pdf_slice->plotOn(frame,Name("pdf_slice2"),LineColor(kBlue+1),LineWidth(3),Normalization(1.,RooAbsReal::RelativeExpected));
        
        pdf_slice->plotOn(frame,LineColor(kGray+3),Components("bkg_partReco2D_{" + str_run[i] + ";" + str_KsCat[j]+ "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(1001),FillColor(kGray+3),Components("bkg_partReco2D_{" + str_run[i] + ";" + str_KsCat[j]+ "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        
        pdf_slice->plotOn(frame,LineColor(kRed+1),Components("signal2D_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3353),FillColor(kRed+1),Components("signal2D_{" + str_run[i] + ";" + str_KsCat[j]+ "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        
        pdf_slice->plotOn(frame,LineColor(kGreen+1),Components("signal_Bs2D_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        pdf_slice->plotOn(frame,DrawOption("F"),FillStyle(3335),FillColor(kGreen+1),Components("signal_Bs2D_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        
        pdf_slice->plotOn(frame,LineColor(kRed+3),Components("signal_B_noKs_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected),LineStyle(kDashed));
        
        pdf_slice->plotOn(frame,LineColor(kGreen+3),Components("signal_Bs_noKs_{" + str_run[i] + ";" + str_KsCat[j] + "}"),Normalization(1.,RooAbsReal::RelativeExpected),LineStyle(kDashed));
        
        pdf_slice->plotOn(frame,LineStyle(kDashed),LineColor(kBlack),Components("bkg_comb2D_{" + str_run[i] + ";" + str_KsCat[j]+ "}"),Normalization(1.,RooAbsReal::RelativeExpected));
        frame->Draw();
        
        cout << endl << "chi2/nbin (Ks) = " << frame->chiSquare("pdf_slice2","data_slice2") << endl << endl;    
                
        c->Print(outDir+"Ks_" + str_run[i] + "_" + str_KsCat[j] + ".eps");
        hpull = frame->pullHist("data_slice2","pdf_slice2") ;
        frame= Ks_PV_MM.frame();
        frame->addPlotable(hpull,"P") ;
        frame->Draw();
        //c->Print(outDir+"Ks_pull_" + str_run[i] + "_" + str_KsCat[j]  + ".eps");
        
        /// Get signal yield
        signal_yield += ((RooRealVar*) fitParams->find("n_sig_{" + str_run[i] + ";" + str_KsCat[j] + "}"))->getVal();
        RooAbsPdf* pdf_slice_comb_bkg = (RooAbsPdf*) pdf_slice_comp.find("bkg_comb2D_{" + str_run[i] + ";" + str_KsCat[j]+ "}");
        comb_bkg_yield += pdf_slice_comb_bkg->createIntegral(DTF_B_M,NormSet(DTF_B_M),Range("signal_range"))->getVal()*((RooRealVar*) fitParams->find("n_bkg_{" + str_run[i] + ";" + str_KsCat[j]+  "}"))->getVal();
        
        /// Calculate sWeights
        if(sWeight){
            RooArgList yield_list(
                                  *((RooRealVar*) fitParams->find("n_sig_{" + str_run[i] + ";" + str_KsCat[j]+ "}")),
                                  *((RooRealVar*) fitParams->find("n_sig_Bs_{" + str_run[i] + ";" + str_KsCat[j]+ "}")),
                                  *((RooRealVar*) fitParams->find("n_sig_noKs_{" + str_run[i] + ";" + str_KsCat[j]+ "}")),
                                  *((RooRealVar*) fitParams->find("n_sig_Bs_noKs_{" + str_run[i] + ";" + str_KsCat[j]+ "}")),
                                  *((RooRealVar*) fitParams->find("n_bkg_{" + str_run[i] + ";" + str_KsCat[j]+ "}")),
                                  *((RooRealVar*) fitParams->find("n_partReco_bkg_{" + str_run[i] + ";" + str_KsCat[j]+  "}"))
                                  );
            RooDataSet* data_slice = new RooDataSet("data_slice","data_slice",data,list,"run==run::" + str_run[i] + " && KsCat==KsCat::" + str_KsCat[j]);
            SPlot sPlot("sPlot","sPlot",*data_slice,pdf_slice,yield_list); 
            
            /// Plot the sWeight distributions as a function of mass
            TH2 * swHist = (TH2*)data_slice->createHistogram("B_DTF_MM,n_sig_{" + str_run[i] + ";" + str_KsCat[j] + "}" + "_sw");
            swHist->GetYaxis()->SetTitle("Signal sWeights");
            swHist->Draw();
            c->Print("plots/signal_sweight_" + str_run[i] + "_" + str_KsCat[j] + ".eps");
            
            /// Save sWeights
            /// Messy and dangerous hack but works for now
            int n_ij = 0;  /// labels entry number of data slice
            for(int n = 0; n < out_tree->GetEntries(); n++){
                b_run->GetEntry(n);
                b_KsCat->GetEntry(n);
                b_TriggerCat->GetEntry(n);
                if(t_run == i+1 && t_KsCat == j){
                    weights[n] = sPlot.GetSWeight(n_ij,"n_sig_{" + str_run[i] + ";" + str_KsCat[j] + "}" + "_sw");
                    weights_Bs[n] = sPlot.GetSWeight(n_ij,"n_sig_Bs_{" + str_run[i] + ";" + str_KsCat[j]+ "}" + "_sw");
                    weights_noKs[n] = sPlot.GetSWeight(n_ij,"n_sig_noKs_{" + str_run[i] + ";" + str_KsCat[j] + "}" + "_sw");
                    weights_Bs_noKs[n] = sPlot.GetSWeight(n_ij,"n_sig_Bs_noKs_{" + str_run[i] + ";" + str_KsCat[j]+ "}" + "_sw");
                    n_ij++;
                }
            }    
        }

        if(sWeightNew){
            RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-10);
            RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-10);
            RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("extrapolation","WynnEpsilon");
            RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("maxSteps","1000");
            RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("minSteps","0");
            RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setCatLabel("method","21Points");
            RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 1000);
            
            const int dim = pdf_slice_comp.getSize();
            TMatrixD W_kl(dim,dim);
            W_kl.Zero();
            RooArgList yield_list( pdf_slice->coefList());
            cout << "calculating sWeights " << endl;
            pdf_slice->Print();
            yield_list.Print();
                
            RooArgSet* pdfObs = pdf_slice->getObservables(data);
            pdfObs->Print();

            string range = "signal_rangeLL";
            if(j==1)range = "signal_rangeDD";
            vector<double> f_k;
            for(int k = 0 ; k< dim; k++) f_k.push_back(((RooAbsPdf*)pdf_slice_comp.at(k))->createIntegral(RooArgSet(DTF_B_M,Ks_PV_MM),NormSet(RooArgSet(DTF_B_M,Ks_PV_MM)),Range(range.c_str()))->getVal());

            for(int n = 0; n < out_tree->GetEntries(); n++){
                b_run->GetEntry(n);
                b_KsCat->GetEntry(n);
                b_TriggerCat->GetEntry(n);
                if(t_run == i+1 && t_KsCat == j){
                    *pdfObs = *data->get(n);
                    double m_B = ((RooRealVar*)((RooArgList*)pdfObs)->at(0))->getVal();
                    double m_Ks = ((RooRealVar*)((RooArgList*)pdfObs)->at(1))->getVal();
                    if(j==0)if( m_B < (double)signalRange_B_min_MM || m_B > (double)signalRange_B_max_MM  || m_Ks < (double)signalRange_Ks_min_MM || m_Ks > (double)signalRange_Ks_max_MM  ) continue;
                    if(j==1)if( m_B < (double)signalRange_B_min_MM || m_B > (double)signalRange_B_max_MM  || m_Ks < Ks_PV_MM.getMin() || m_Ks > Ks_PV_MM.getMax() ) continue;

                    double denom = 0.;
                    for(int k = 0 ; k< dim; k++){ 
                        double pdf_k = ((RooAbsPdf*)pdf_slice_comp.at(k))->getVal(pdfObs);
                        double yield_k = ((RooAbsReal*)yield_list.at(k))->getVal();
                        yield_k *=  f_k[k];
                        denom += yield_k * pdf_k;
                    }
                    
                    for(int k = 0 ; k< dim; k++)for(int l = 0 ; l< dim; l++){ 
                        double pdf_k = ((RooAbsPdf*)pdf_slice_comp.at(k))->getVal(pdfObs);
                        double pdf_l = ((RooAbsPdf*)pdf_slice_comp.at(l))->getVal(pdfObs);
                        W_kl[k][l] += pdf_k * pdf_l / (denom*denom);
                    }
                }
            }
            W_kl.Print();
            W_kl.Invert();
            W_kl.Print();
            
            TMatrixD alpha(1,dim);
            alpha.Zero();
            alpha[0][0] = 1.;
            alpha *= W_kl;
            alpha.Print();

            TMatrixD alpha_Bs(1,dim);
            alpha_Bs.Zero();
            alpha_Bs[0][1] = 1.;
            alpha_Bs *= W_kl;

            TMatrixD alpha_noKs(1,dim);
            alpha_noKs.Zero();
            alpha_noKs[0][2] = 1.;
            alpha_noKs *= W_kl;

            TMatrixD alpha_Bs_noKs(1,dim);
            alpha_Bs_noKs.Zero();
            alpha_Bs_noKs[0][3] = 1.;
            alpha_Bs_noKs *= W_kl;

            double sum_sw = 0;
            double sum_sw2 = 0;
            
            for(int n = 0; n < out_tree->GetEntries(); n++){
                b_run->GetEntry(n);
                b_KsCat->GetEntry(n);
                b_TriggerCat->GetEntry(n);
                if(t_run == i+1 && t_KsCat == j){
                    *pdfObs = *data->get(n);
                    double m_B = ((RooRealVar*)((RooArgList*)pdfObs)->at(0))->getVal();
                    double m_Ks = ((RooRealVar*)((RooArgList*)pdfObs)->at(1))->getVal();
                    if( m_B < (double)signalRange_B_min_MM || m_B > (double)signalRange_B_max_MM  || m_Ks < (double)signalRange_Ks_min_MM || m_Ks > (double)signalRange_Ks_max_MM  ){
                        weights2[n] = 0;
                        weights2_Bs[n] = 0;
                        weights2_noKs[n] = 0;
                        weights2_Bs_noKs[n] = 0;
                        continue;
                    } 

                    double num = 0.;
                    double num_Bs = 0.;
                    double num_noKs = 0.;
                    double num_Bs_noKs = 0.;
                    double den = 0.;
                    for(int k = 0 ; k< dim; k++){
                        double pdf_k = ((RooAbsPdf*)pdf_slice_comp.at(k))->getVal(pdfObs);
                        double yield_k = ((RooAbsReal*)yield_list.at(k))->getVal();
                        yield_k *=  f_k[k];
                        num += alpha[0][k] * pdf_k;
                        num_Bs += alpha_Bs[0][k] * pdf_k;
                        num_noKs += alpha_noKs[0][k] * pdf_k;
                        num_Bs_noKs += alpha_Bs_noKs[0][k] * pdf_k;
                        den += yield_k * pdf_k;
                    }
                    weights2[n] = num/den;
                    sum_sw += weights2[n];                        
                    sum_sw2 += weights2[n]*weights2[n];        
                    weights2_Bs[n] = num_Bs/den;
                    weights2_noKs[n] = num_noKs/den;
                    weights2_Bs_noKs[n] = num_Bs_noKs/den;
                }
            }      
             cout << "sum_sw = " << sum_sw << endl << endl;        } 
    }

    if(sWeight || sWeightNew){
        for(int n = 0; n < out_tree->GetEntries(); n++){
            sw = weights[n];
            sw_Bs = weights_Bs[n];
            sw_noKs = weights_noKs[n];
            sw_Bs_noKs = weights_Bs_noKs[n];

            sw2 = weights2[n];
            sw2_Bs = weights2_Bs[n];
            sw2_noKs = weights2_noKs[n];
            sw2_Bs_noKs = weights2_Bs_noKs[n];
    
            b_sw->Fill();
            b_sw_Bs->Fill();
            b_sw_noKs->Fill();
            b_sw_Bs_noKs->Fill();

            b_sw2->Fill();
            b_sw2_Bs->Fill();
            b_sw2_noKs->Fill();
            b_sw2_Bs_noKs->Fill();
        }
        out_tree->Write();
        output->Close();
        cout << endl;
        cout << "Created file " << ((string)outFileName).c_str()  << endl << endl;
    }
    chi2 = chi2*nBins/(str_run.size()*str_KsCat.size()*nBins-fitParams->selectByAttrib("Constant",kFALSE)->getSize());
    cout << endl; cout<<"Chi2 data= " << chi2 <<" Cov Matrix: "<<covmatr<<" EDM: "<<edm<<endl<<endl;
    cout << "Total signal yield = " << signal_yield << endl;
}

void plotDalitz(){
    
    TFile* f = new TFile("b2dkspi_sw.root");
    TTree* t = (TTree*) f->Get("DecayTree");

    double sw; 
    t->SetBranchAddress("N_B_sw",&sw);
    double m_DKs,m_Dpi,m_Kspi;
    t->SetBranchAddress("FullDTF_m_DKs",&m_DKs);
    t->SetBranchAddress("FullDTF_m_Dpi",&m_Dpi);
    t->SetBranchAddress("FullDTF_m_Kspi",&m_Kspi);
    
    double nBins = 60;
    TH1D* h_DKs = new TH1D("h_DKs","; m(DK_{s}) [GeV]; Yield",nBins,2,5.5);
    TH1D* h_Dpi = new TH1D("h_Dpi","; m(D#pi) [GeV]; Yield",nBins,1,5.5);
    TH1D* h_Kspi = new TH1D("h_Kspi","; m(K_{s}#pi) [GeV]; Yield",nBins,0,4);
    
    TH2D* h_Kspi_DKs = new TH2D("h_Kspi_DKs","; m(K_{s}#pi) [GeV]; m(DK_{s}) [GeV];  Yield",nBins,0,4,nBins,2,5.5);
    TH2D* h_Kspi_Dpi = new TH2D("h_Kspi_Dpi","; m(K_{s}#pi) [GeV]; m(D#pi) [GeV];  Yield",nBins,0,4,nBins,2,5.5);
    TH2D* h_Dpi_DKs = new TH2D("h_Dpi_DKs","; m(D#pi) [GeV]; m(DK_{s}) [GeV];  Yield",nBins,2,5.5,nBins,2,5.5);

    for (int i=0; i<t->GetEntries(); i++) {
        t->GetEntry(i);
        
        //if(m_Kspi<1000)continue;
        
        h_DKs->Fill(m_DKs/1000.,sw);
        h_Dpi->Fill(m_Dpi/1000.,sw);
        h_Kspi->Fill(m_Kspi/1000.,sw);
        h_Kspi_DKs->Fill(m_Kspi/1000.,m_DKs/1000,sw);
        h_Kspi_Dpi->Fill(m_Kspi/1000.,m_Dpi/1000,sw);
        h_Dpi_DKs->Fill(m_Dpi/1000.,m_DKs/1000,sw);
    }
    
    TCanvas* c = new TCanvas();

    h_DKs->SetMinimum(0.);
    h_Dpi->SetMinimum(0.);
    h_Kspi->SetMinimum(0.);
    h_Kspi_DKs->SetMinimum(0.);
    h_Kspi_DKs->SetMarkerSize(0.4);
    h_Kspi_Dpi->SetMinimum(0.);
    h_Kspi_Dpi->SetMarkerSize(0.4);
    h_Dpi_DKs->SetMinimum(0.);
    h_Dpi_DKs->SetMarkerSize(0.4);
    
    h_DKs->Draw();
    c->Print("plots/m_DKs.eps");
    h_Dpi->Draw();
    c->Print("plots/m_Dpi.eps");
    h_Kspi->Draw();
    c->Print("plots/m_Kspi.eps");
   
    gPad->SetLogz(1);
    
    h_Kspi_DKs->Draw("colz");
    c->Print("plots/m_Kspi_DKs.eps");
    h_Kspi_DKs->Draw();
    c->Print("plots/m_Kspi_DKs_scatter.eps");

    h_Kspi_Dpi->Draw("colz");
    c->Print("plots/h_Kspi_Dpi.eps");
    h_Kspi_Dpi->Draw();
    c->Print("plots/m_Kspi_Dpi_scatter.eps");

    h_Dpi_DKs->Draw("colz");
    c->Print("plots/h_Dpi_DKs.eps");
    h_Dpi_DKs->Draw();
    c->Print("plots/m_Dpi_DKs_scatter.eps");
}

void testSweights(int step = 0){
       
    RooAbsReal::defaultIntegratorConfig()->setEpsAbs(1e-9);
    RooAbsReal::defaultIntegratorConfig()->setEpsRel(1e-9);
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("extrapolation","WynnEpsilon");
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("maxSteps","1000");
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooIntegrator1D").setCatLabel("minSteps","0");
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setCatLabel("method","21Points");
    RooAbsReal::defaultIntegratorConfig()->getConfigSection("RooAdaptiveGaussKronrodIntegrator1D").setRealValue("maxSeg", 1000);
    
    /// Signal pdf
    TString channelString = "m(D_{s}K#pi#pi)" ;
    RooRealVar DTF_Bs_M("Bs_DTF_MM", channelString, 5200, 5700,"MeV/c^{2}");
    
    RooRealVar mean("mean", "mean", 5367,5350.,5390.); 
    RooRealVar sigma("sigma", "sigma", 20.,0.,80.); 
    RooGaussian* signal = new RooGaussian("signal","signal",DTF_Bs_M, mean,sigma);
    mean.setConstant();
    sigma.setConstant();
    
    /// Bkg pdf
    RooRealVar c0("c0", "c0", .0,-10,10); 
    RooRealVar c1("c1", "c1", .0,-10,10); 
    RooRealVar c2("c2", "c2", .0,-10,10); 
    RooChebychev* bkg= new RooChebychev("bkg","bkg",DTF_Bs_M, RooArgList(c0,c1));
    
    RooRealVar exp_par("exp_par","exp_par",-1.6508e-03,-10.,10.);
    RooExponential bkg_exp("bkg_exp","bkg_exp",DTF_Bs_M,exp_par);
    
    NamedParameter<double>  f_bkg("readBkgData::f_bkg", 0.);
    RooRealVar mean_bkg("mean_bkg", "mean_bkg", 5241.,5350.,5390.); 
    RooRealVar sigma_bkg("sigma_bkg", "sigma_bkg", 37.,0.,80.); 
    RooGaussian* partReco = new RooGaussian("partReco","partReco",DTF_Bs_M, mean_bkg,sigma_bkg);
    //     if(f_bkg==0.){
    mean_bkg.setConstant();
    sigma_bkg.setConstant();
    //     }
    /// Total pdf
    int N = 10000;
    RooRealVar n_sig("n_sig", "n_sig", N*0.5, 0., N);
    RooRealVar n_bkg("n_bkg", "n_bkg", N*0.5, 0., N);
    RooRealVar n_bkg2("n_bkg2", "n_bkg2", N*0.25, 0., N);
    
    RooAddPdf* pdf;
    if(f_bkg==0.) pdf = new RooAddPdf("pdf", "pdf", RooArgList(*signal, *bkg), RooArgList(n_sig, n_bkg));
    else pdf = new RooAddPdf("pdf", "pdf", RooArgList(*signal, bkg_exp, *partReco), RooArgList(n_sig, n_bkg,n_bkg2));
    
    // generate toys
    RooRealVar t("t","t", 0,10,"");
    RooRealVar tau_sig("tau_sig","tau_sig",-1.);
    RooRealVar tau_bkg("tau_bkg","tau_bkg",0);

    RooExponential signal_t("signal_t","signal_t",t,tau_sig);
    RooProdPdf signal2D("signal2D","signal2D",*signal,signal_t);

    RooExponential bkg_t("bkg_t","bkg_t",t,tau_bkg);
    RooProdPdf bkg2D("bkg2D","bkg2D",*bkg,bkg_t);
    
    RooAddPdf* pdf2D = new RooAddPdf("pdf2D", "pdf2D", RooArgList(signal2D, bkg2D), RooArgList(n_sig, n_bkg));
    RooDataSet* data = pdf2D->generate(RooArgList(DTF_Bs_M,t),N);
    
    /// Fit
    RooFitResult* result = pdf->fitTo(*data,Save(kTRUE),NumCPU(6),Extended(kTRUE));
    result->Print();
    
    /// Calculate sweights
    TH2D *h_sw = new TH2D("h_sw","h_sw",100,5200,5700,100,-5,5);
    TH2D *h_sw2 = new TH2D("h_sw2","h_sw2",100,5200,5700,100,-5,5);
    TH1D *h_t = new TH1D("h_t","h_t",100,0,10);
    TH1D *h_t_sw = new TH1D("h_t_sw","h_t_sw",100,0,10);
    TH1D *h_t_sw2 = new TH1D("h_t_sw2","h_t_sw2",100,0,10);

    SPlot *sPlot;
    if(f_bkg==0.)sPlot= new SPlot("sPlot","sPlot",*data,pdf,RooArgList(n_sig,n_bkg)); 
    else sPlot= new SPlot("sPlot","sPlot",*data,pdf,RooArgList(n_sig,n_bkg,n_bkg2));    
        
    RooArgList pdf_slice_comp(pdf->pdfList());
    const int dim = pdf_slice_comp.getSize();
    TMatrixD W_kl(dim,dim);
    W_kl.Zero();
    RooArgList yield_list( pdf->coefList());
    cout << "calculating sWeights " << endl;
    pdf->Print();
    yield_list.Print();
        
    RooArgSet* pdfObs = pdf2D->getObservables(data) ;
    DTF_Bs_M.setRange("signal",5350,5400);
    vector<double> f_k;
    for(int k = 0 ; k< dim; k++) f_k.push_back(((RooAbsPdf*)pdf_slice_comp.at(k))->createIntegral(DTF_Bs_M,NormSet(DTF_Bs_M),Range("signal_range"))->getVal());
    for(int n = 0; n < N; n++){
            *pdfObs = *data->get(n);
            if(((RooRealVar*)((RooArgList*)pdfObs)->at(0))->getVal() < 5350) continue;
            if(((RooRealVar*)((RooArgList*)pdfObs)->at(0))->getVal() > 5400) continue;         
            
            double denom = 0.;
            for(int k = 0 ; k< dim; k++){ 
                double pdf_k = ((RooAbsPdf*)pdf_slice_comp.at(k))->getVal(pdfObs);
                double yield_k = ((RooAbsReal*)yield_list.at(k))->getVal();
                yield_k *=  f_k[k];
                denom += yield_k * pdf_k;
            }
            
            for(int k = 0 ; k< dim; k++)for(int l = 0 ; l< dim; l++){ 
                double pdf_k = ((RooAbsPdf*)pdf_slice_comp.at(k))->getVal(pdfObs);
                double pdf_l = ((RooAbsPdf*)pdf_slice_comp.at(l))->getVal(pdfObs);
                W_kl[k][l] += pdf_k * pdf_l / (denom*denom);
            }
    }
    W_kl.Print();
    W_kl.Invert();
    W_kl.Print();
    
    TMatrixD alpha(1,dim);
    alpha.Zero();
    alpha[0][0] = 1.;
    alpha *= W_kl;
    alpha.Print();

    vector<double> w(dim);
    double sum_sw = 0;
    double sum_sw2 = 0;
    
    for(int n = 0; n < N; n++){
        *pdfObs = *data->get(n);

        double num = 0.;
        double den = 0.;
        for(int k = 0 ; k< dim; k++){
            double pdf_k = ((RooAbsPdf*)pdf_slice_comp.at(k))->getVal(pdfObs);
            double yield_k = ((RooAbsReal*)yield_list.at(k))->getVal();
            yield_k *=  f_k[k];
            num += alpha[0][k] * pdf_k;
            den += yield_k * pdf_k;
        }

        sum_sw += sPlot->GetSWeight(n,"n_sig_sw");
        h_sw->Fill(((RooRealVar*)((RooArgList*)pdfObs)->at(0))->getVal(), sPlot->GetSWeight(n,"n_sig_sw"));        
        h_t->Fill(((RooRealVar*)((RooArgList*)pdfObs)->at(1))->getVal());
        h_t_sw->Fill(((RooRealVar*)((RooArgList*)pdfObs)->at(1))->getVal(), sPlot->GetSWeight(n,"n_sig_sw"));
        
        if(((RooRealVar*)((RooArgList*)pdfObs)->at(0))->getVal() < 5350) continue;
        if(((RooRealVar*)((RooArgList*)pdfObs)->at(0))->getVal() > 5400) continue;         

        h_t_sw2->Fill(((RooRealVar*)((RooArgList*)pdfObs)->at(1))->getVal(), num/den);        
        h_sw2->Fill(((RooRealVar*)((RooArgList*)pdfObs)->at(0))->getVal(), num/den);
        sum_sw2 += num/den;        
        
    }
    
    cout << sum_sw << endl;
    cout << sum_sw2 << endl;
        
    /// Plotting
    TCanvas* c = new TCanvas();
    
    RooPlot* frame= DTF_Bs_M.frame();
    frame->SetTitle("");
    data->plotOn(frame,Name("data"),Binning(100));
    pdf->plotOn(frame,Name("pdf"),LineColor(kBlue));
    pdf->plotOn(frame,Name("pdf"),LineColor(kRed),LineStyle(kDashed),Components("signal"));
    pdf->plotOn(frame,Name("pdf"),LineColor(kRed),LineStyle(kDashed),Components("bkg_exp"));
    pdf->plotOn(frame,Name("pdf"),LineColor(kRed),LineStyle(kDashed),Components("bkg"));
    frame->Draw();
    c->Print("toyFit.eps");
    
    h_sw->Draw("");
    h_sw2->SetLineColor(kRed);
    h_sw2->SetMarkerColor(kRed);
    h_sw2->Draw("same");
    c->Print("weights.eps");

    h_t->Draw("");
    h_t_sw->SetLineColor(kBlue);
    h_t_sw->Draw("histcsame");
    h_t_sw2->SetLineColor(kRed);
    h_t_sw2->Draw("histcsame");
    c->Print("t.eps");
}

void makeReducedFile(){
        NamedParameter<string> outFileName("outFileNameSignal",(string)"./signal_data.root");
        NamedParameter<string> cut("makeReducedFile::Cut",(string)"N_B_sw2 != 0");

        TFile *file= new TFile(((string)outFileName).c_str());
        TTree* tree = (TTree*) file->Get("DecayTree");    

        TFile *output = new TFile((TString((string)outFileName).ReplaceAll(".root","_reduced.root")),"RECREATE");
        tree->SetBranchStatus("*",0);
        tree->SetBranchStatus("*MM*",1);
        tree->SetBranchStatus("*sw*",1);
        tree->SetBranchStatus("m_*",1);
        tree->SetBranchStatus("Full*",1);
        tree->SetBranchStatus("B_Full*TAU*",1);

        TTree* out_tree = tree->CopyTree(((string)cut).c_str());
        out_tree->Write();
        output->Close();
        cout << endl;
        cout << "Created file " << TString((string)outFileName).ReplaceAll(".root","_reduced.root")  << endl << endl;
}


int main(int argc, char** argv){

    time_t startTime = time(0);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");

    str_year.push_back(TString("y11"));
    str_year.push_back(TString("y12"));
    str_year.push_back(TString("y15"));
    str_year.push_back(TString("y16"));
    str_year.push_back(TString("y17"));
    str_year.push_back(TString("y18"));

    str_run.push_back(TString("Run1"));
    str_run.push_back(TString("Run2"));
    
    str_KsCat.push_back(TString("LL"));
    str_KsCat.push_back(TString("DD"));

    str_trigger.push_back(TString("t0"));
    str_trigger.push_back(TString("t1"));

    
    //fitSignalShape("signal");
    //fitPartRecoBkgShape();
    fitSignal2D();
    //plotDalitz();
    NamedParameter<int> reduceFile("makeReducedFile", 0);
    if(reduceFile)makeReducedFile();
 
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
