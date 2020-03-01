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
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include <ctime>
#include "Mint/NamedParameter.h"
#include "Mint/Utils.h"
#include "Mint/RooHILLdini.h"
#include "Mint/RooHORNSdini.h"

using namespace std;
using namespace RooFit ;
using namespace RooStats;
using namespace MINT;

vector<TString> str_year;
vector<TString> str_run;
vector<TString> str_Ds;
vector<TString> str_trigger;
static const double massKaon = 493.68;
static const double massPion = 139.57;

void fitSignal(){

	///Options
	NamedParameter<int> numCPU("numCPU", 6);
	NamedParameter<int> sWeight("sWeightSignal", 0);
	NamedParameter<int> nBins("nBins", 80);
	NamedParameter<double> min_MM("min_MM",5100.);
	NamedParameter<double> max_MM("max_MM",5700.);
	NamedParameter<string> cut_BDT("cut_BDT",(string)"");
	NamedParameter<string> inFileName("inFileNameSignal",(string)"/auto/data/dargent/BsDsKpipi/BDT/Data/signal.root");
	NamedParameter<string> outFileName("outFileNameSignal",(string)"/auto/data/dargent/BsDsKpipi/Final/Data/signal.root");

	///Load file
	TFile *file= new TFile(((string)inFileName).c_str());
	TTree* tree = (TTree*) file->Get("DecayTree");	
   	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("B_DTF_MM",1);
	
    RooRealVar DTF_B_M("B_DTF_MM", "m(D^{-}K_{S}#pi^{+})", min_MM, max_MM,"MeV/c^{2}");	
	RooArgList list =  RooArgList(DTF_B_M);
	//RooArgList list2 =  RooArgList(Bs_BsDTF_TAU,Bs_BsDTF_TAUERR,m_Kpipi,m_Kpi,m_pipi,Ds_FDCHI2_ORIVX,pi_minus_PIDK);
	//list.add(list2);
	RooDataSet*  data = new RooDataSet("data","data",tree,list,((string)cut_BDT).c_str() );	
    
	/// Signal Pdf
	RooRealVar mean("mean","mean", 5279, 5250 , 5300); 
	RooRealVar sigma("sigma","sigma", 30,0,100);
    RooRealVar alpha("alpha", "alpha", 0., -100, 100); 
	RooRealVar beta("beta", "beta", 0., -100, 100); 
	RooJohnsonSU signal1("signal1","signal1",DTF_B_M, mean,sigma,alpha,beta);

	/// B0 pdf
	RooFormulaVar mean_B0("mean_B0","@0 - @1", RooArgSet(mean,RooConst(87.33))); 
	RooJohnsonSU signal1_B0("signal1_B0","signal1_B0",DTF_B_M, mean_B0,sigma,alpha,beta);
	
	/// Combinatorial bkg pdf
	RooRealVar exp_par("exp_par","exp_par",-1.6508e-03,-10.,10.);
    RooExponential bkg_exp1("bkg_exp1","bkg_exp1",DTF_B_M,exp_par);


	/// Perform fit
    RooRealVar n_sig("n_sig", "n_sig", data->numEntries()*0.15, 0., data->numEntries());
    RooRealVar n_sig_B0("n_sig_B0", "n_sig_B0", data->numEntries()*0.15, 0., data->numEntries());
    RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()*0.7, 0., data->numEntries());
    
    RooAddPdf* pdf = new RooAddPdf("pdf", "pdf", RooArgList(signal1,signal1_B0, bkg_exp1), RooArgList(n_sig,n_sig_B0, n_bkg));

	RooFitResult* result = pdf->fitTo(*data,Save(kTRUE),Extended(kTRUE),NumCPU(numCPU));
	cout << "result is --------------- "<<endl;
	result->Print("v");

	/// Plot combined data and fit
	TCanvas* c = new TCanvas();
	RooPlot* frame= DTF_B_M.frame();
	frame->SetTitle("");
	data->plotOn(frame,Name("data"),MarkerSize(1),Binning(nBins));
	pdf->plotOn(frame,Name("pdf"),LineColor(kBlue+1),LineWidth(3));

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

    frame->Draw();
	leg.Draw();
	c->Print("eps/signal.eps");

	/// Output file
	TFile *output;
	TTree* out_tree;
	double sw,sw_B0;
	Int_t t_year, t_run, t_Ds_finalState, t_TriggerCat;
	TBranch *b_sw, *b_w, *b_year, *b_run, *b_Ds_finalState, *b_TriggerCat;

	if(sWeight){
		output = new TFile(((string)outFileName).c_str(),"RECREATE");
		tree->SetBranchStatus("*",1);
		tree->SetBranchStatus("weight",0);
		tree->SetBranchStatus("N_Bs_sw",0);

        out_tree = tree->CopyTree(("B_DTF_MM >= " + anythingToString((double)min_MM) + " && B_DTF_MM <= " + anythingToString((double)max_MM) + " && " + (string)cut_BDT).c_str());
        b_sw = out_tree->Branch("N_B_sw", &sw, "N_B_sw/D");
        b_w = out_tree->Branch("weight", &sw_B0, "weight/D");
        if(out_tree->GetEntries() != data->numEntries()) {
            cout << "ERROR:: Different number of events in input and outputfile ! " << endl;
            cout << out_tree->GetEntries() << endl;
            cout << data->numEntries() << endl;
            throw "ERROR";
        }
    }    
    
	double weights[(int)data->numEntries()];
	double weights_B0[(int)data->numEntries()];

	/// Calculate total signal yield
	double signal_yield = 0.;
	double comb_bkg_yield = 0.;
	int n_sig_perFit = 0.;
	int n_sig_perFit_err = 0.;
	DTF_B_M.setRange("signal_range",mean.getVal()-45.,mean.getVal()+45.);

	if(sWeight){
		for(int n = 0; n < out_tree->GetEntries(); n++){
			sw = weights[n];
			sw_B0 = weights_B0[n];
			b_sw->Fill();
			b_w->Fill();
		}
	 	out_tree->Write();
   		output->Close();
		cout << endl;
		cout << "Created file " << outFileName.c_str()  << endl << endl;
	}
}

int main(int argc, char** argv){

    time_t startTime = time(0);
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");

    fitSignal();
 
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
