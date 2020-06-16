// Time fits
// author: Philippe d'Argent
#include "Mint/FitParameter.h"
#include "Mint/NamedParameter.h"
#include "Mint/Minimiser.h"
#include "Mint/Neg2LL.h"
#include "Mint/Neg2LLSum.h"
#include "Mint/Neg2LLMultiConstraint.h"
#include "Mint/DalitzEventList.h"
#include "Mint/NamedDecayTreeList.h"
#include "Mint/DecayTree.h"
#include "Mint/DiskResidentEventList.h"
#include "Mint/CLHEPPhysicalConstants.h"
#include "Mint/CLHEPSystemOfUnits.h"
#include "Mint/PdfBase.h"
#include "Mint/DalitzPdfBase.h"
#include "Mint/DalitzPdfBaseFastInteg.h"
#include "Mint/DalitzPdfBaseFlexiFastInteg.h"
#include "Mint/FitAmplitude.h"
#include "Mint/FitAmpSum.h"
#include "Mint/FitAmpIncoherentSum.h"
#include "Mint/DalitzEvent.h"
#include "Mint/AmpRatios.h"
#include "Mint/IEventGenerator.h"
#include "Mint/DalitzBWBoxSet.h"
#include "Mint/DalitzBoxSet.h"
#include "Mint/SignalGenerator.h"
#include "Mint/FromFileGenerator.h"
#include "Mint/DalitzSumPdf.h"
#include "Mint/cexp.h"
#include "Mint/DalitzPdfNormChecker.h"
#include "Mint/IFastAmplitudeIntegrable.h"
#include "Mint/DalitzPdfSaveInteg.h"
#include "Mint/Chi2Binning.h"
#include "Mint/FitAmpIncoherentSum.h"
#include "Mint/FitAmpList.h"
#include "Mint/DalitzPdfBaseMCInteg.h"
#include "Mint/TimePdfMaster.h"
#include "Mint/FullTimePdf.h"

#include "RooRealConstant.h"
#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDecay.h"
#include "RooBDecay.h"
#include "RooPlot.h"
#include "RooEffProd.h"
#include "RooGenericPdf.h"
#include "RooGaussModel.h"
#include "RooProdPdf.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooUniform.h"
#include "RooExponential.h"
#include "RooRandom.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "RooMultiVarGaussian.h"
#include "RooTruthModel.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "Mint/RooCubicSplineFun.h"
#include "Mint/RooCubicSplineKnot.h"
#include "Mint/RooGaussEfficiencyModel.h"
#include "Mint/DecRateCoeff_Bd.h"
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TNtupleD.h"
#include "TTree.h"
#include "TFile.h"
#include <TStyle.h>
#include <TROOT.h>
#include "TRandom2.h"
#include "TRandom3.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>
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
#include <TLegend.h>
#include <TLegendEntry.h>
#include <TPaveText.h>
#include <TNtuple.h>
#include "TRandom3.h"
#include <sstream>
#include "TProfile.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TDecompChol.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit;
using namespace RooStats;
using namespace MINT;

void fullTimeFit(int step=0, string mode = "fit"){

    /// Options
    NamedParameter<int> updateAnaNote("updateAnaNote", 0);
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    int seed = RandomSeed + step;
    ranLux.SetSeed((int)seed);
    gRandom = &ranLux;
    RooRandom::randomGenerator()->SetSeed(seed);

    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    DalitzEventPattern pat_CP = pat.makeCPConjugate();

    NamedParameter<string> InputFileName("InputFileName", (std::string) "/auto/data/dargent/BsDsKpipi/Final/signal_tagged.root", (char*) 0);
    NamedParameter<string> channel("channel", (std::string) "norm", (char*) 0);
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> max_TAU_ForMixingPlot("max_TAU_ForMixingPlot", 12.);
    NamedParameter<double> min_TAUERR("min_TAUERR", 0.);
    NamedParameter<double> max_TAUERR("max_TAUERR", 0.1);

    NamedParameter<int>  useTrueTau("useTrueTau", 0);
    NamedParameter<int>  useTrueTagging("useTrueTagging", 0);

    NamedParameter<double> w_min("w_min", 0.);
    NamedParameter<double> w_max("w_max", 1.);

    NamedParameter<double> min_year("min_year", 10);
    NamedParameter<double> max_year("max_year", 20);
    NamedParameter<double> min_trigger("min_trigger", -10);
    NamedParameter<double> max_trigger("max_trigger", 10);

    NamedParameter<int>  nBinst("nBinst", 40);
    NamedParameter<int>  nBinsAsym("nBinsAsym", 10);
    NamedParameter<int>  doPlots("doPlots", 1);

    NamedParameter<int>  do2DScan("do2DScan", 0);
    NamedParameter<int>  doSimFit("doSimFit", 0);
    NamedParameter<int>  doSimFitInBins("doSimFitInBins", 0);

    NamedParameter<int>  doBootstrap("doBootstrap", 0);
    NamedParameter<int>  N_bootstrap("N_bootstrap", 10000);

    NamedParameter<int>  doToyStudy("doToyStudy", 0);
    NamedParameter<double> N_scale_toys("N_scale_toys", 1);
    NamedParameter<int>  addBkgToToys("addBkgToToys", 0);

    NamedParameter<int>  useGaussConstrainsTagging("useGaussConstrainsTagging", 0);
    NamedParameter<int>  useGaussConstrainsBias("useGaussConstrainsBias", 0);

    NamedParameter<int>  doAccSystematics("doAccSystematics", 0);
    NamedParameter<int>  useCholDec("useCholDec", 0);
    NamedParameter<int>  varPerParChol("varPerParChol", 100);
    int chol_index = (step-1)/varPerParChol ;

    NamedParameter<string> doSystematic("doSystematic", (std::string) "", (char*) 0);
    NamedParameter<int>  randomizeStartVals("randomizeStartVals", 0);

    NamedParameter<string> weightName("weightName", (std::string) "N_Bs_sw", (char*) 0);

    NamedParameter<int>  dilutionWeight("dilutionWeight", 1);
    NamedParameter<int>  N_plot_it("N_plot_it",10);
    NamedParameter<int> scale_asym("scale_asym", 1);

    /// Common fit parameters
    FitParameter  r("r",1,0.,0.1);
    FitParameter  delta("delta",1,100.,1.);
    FitParameter  gamma("gamma",1,70,1.);
    FitParameter  k("k",1,1,1.);
        
    FitParameter  C("C",1,0.,0.1);
    FitParameter  D("D",1,0.,0.1);
    FitParameter  D_bar("D_bar",1,0.,0.1);
    FitParameter  S("S",1,0.,0.1);
    FitParameter  S_bar("S_bar",1,0.,0.1);
    
    FitParameter  Gamma("Gamma",2,0.6629,0.0018);
    FitParameter  dGamma("dGamma",2,-0.088,0.006);
    FitParameter  dm("dm",2,17.757,0.1);
    
    FitParameter  offset_mean_dt("offset_mean_dt",1,0,0.1);
    FitParameter  scale_mean_dt("scale_mean_dt",1,0,0.1);
    FitParameter  scale_mean_2_dt("scale_mean_2_dt",1,0,0.1);
    FitParameter  offset_sigma_dt("offset_sigma_dt",1,0.,0.1);
    FitParameter  scale_sigma_dt("scale_sigma_dt",1,1.,0.1);
    FitParameter  scale_sigma_2_dt("scale_sigma_2_dt",1,0.,0.1);
    FitParameter  offset_sigma2_dt("offset_sigma2_dt",1,0.,0.1);
    FitParameter  scale_sigma2_dt("scale_sigma2_dt",1,1.,0.1);
    FitParameter  scale_sigma2_2_dt("scale_sigma2_2_dt",1,0.,0.1);
    FitParameter  offset_sigma3_dt("offset_sigma3_dt",1,0.,0.1);
    FitParameter  scale_sigma3_dt("scale_sigma3_dt",1,1.,0.1);
    FitParameter  scale_sigma3_2_dt("scale_sigma3_2_dt",1,0.,0.1);
    FitParameter  offset_f_dt("offset_f_dt",1,1,0.1);
    FitParameter  scale_f_dt("scale_f_dt",1,0.,0.1);
    FitParameter  scale_f_2_dt("scale_f_2_dt",1,0.,0.1);
    FitParameter  offset_f2_dt("offset_f2_dt",1,0.,0.1);
    FitParameter  scale_f2_dt("scale_f2_dt",1,0.,0.1);
    FitParameter  scale_f2_2_dt("scale_f2_2_dt",1,0.,0.1);

    FitParameter  p0_os("p0_os",1,0.,0.);
    FitParameter  p1_os("p1_os",1,1.,0.);
    FitParameter  delta_p0_os("delta_p0_os",1,0.,0.);
    FitParameter  delta_p1_os("delta_p1_os",1,0.,0.);
    FitParameter  avg_eta_os("avg_eta_os",1,0.4,0.);
    FitParameter  tageff_os("tageff_os",1,0.4,0.);
    FitParameter  tageff_asym_os("tageff_asym_os",1,0.,0.);
    FitParameter  p0_ss("p0_ss",1,0.,0.);
    FitParameter  p1_ss("p1_ss",1,1.,0.);
    FitParameter  delta_p0_ss("delta_p0_ss",1,0.,0.);
    FitParameter  delta_p1_ss("delta_p1_ss",1,0.,0.);
    FitParameter  avg_eta_ss("avg_eta_ss",1,0.4,0.);
    FitParameter  tageff_ss("tageff_ss",1,0.7,0.);
    FitParameter  tageff_asym_ss("tageff_asym_ss",1,0.,0.);
    FitParameter  production_asym("production_asym",1,0.,0.);
    FitParameter  detection_asym("detection_asym",1,0.,0.);
    
    FitParameter  c0("c0",1,1,0.1);
    FitParameter  c1("c1",1,1,0.1);
    FitParameter  c2("c2",1,1,0.1);
    FitParameter  c3("c3",1,1,0.1);
    FitParameter  c4("c4",1,1,0.1);
    FitParameter  c5("c5",1,1,0.1);
    FitParameter  c6("c6",1,1,0.1);
    FitParameter  c7("c7",1,1,0.1);
    FitParameter  c8("c8",1,1,0.1);
    FitParameter  c9("c9",1,1,0.1);
    
    /// Fit parameters per Run    
    FitParameter  offset_mean_dt_Run1("offset_mean_dt_Run1",1,0,0.1);
    FitParameter  scale_mean_dt_Run1("scale_mean_dt_Run1",1,0,0.1);
    FitParameter  scale_mean_2_dt_Run1("scale_mean_2_dt_Run1",1,0,0.1);
    FitParameter  offset_sigma_dt_Run1("offset_sigma_dt_Run1",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run1("scale_sigma_dt_Run1",1,1.,0.1);
    FitParameter  scale_sigma_2_dt_Run1("scale_sigma_2_dt_Run1",1,0.,0.1);
    FitParameter  offset_sigma2_dt_Run1("offset_sigma2_dt_Run1",1,0.,0.1);
    FitParameter  scale_sigma2_dt_Run1("scale_sigma2_dt_Run1",1,1.,0.1);
    FitParameter  scale_sigma2_2_dt_Run1("scale_sigma2_2_dt_Run1",1,0.,0.1);
    FitParameter  offset_sigma3_dt_Run1("offset_sigma3_dt_Run1",1,0.,0.1);
    FitParameter  scale_sigma3_dt_Run1("scale_sigma3_dt_Run1",1,1.,0.1);
    FitParameter  scale_sigma3_2_dt_Run1("scale_sigma3_2_dt_Run1",1,0.,0.1);
    FitParameter  offset_f_dt_Run1("offset_f_dt_Run1",1,1,0.1);
    FitParameter  scale_f_dt_Run1("scale_f_dt_Run1",1,0.,0.1);
    FitParameter  scale_f_2_dt_Run1("scale_f_2_dt_Run1",1,0.,0.1);
    FitParameter  offset_f2_dt_Run1("offset_f2_dt_Run1",1,0.,0.1);
    FitParameter  scale_f2_dt_Run1("scale_f2_dt_Run1",1,0.,0.1);
    FitParameter  scale_f2_2_dt_Run1("scale_f2_2_dt_Run1",1,0.,0.1);

    FitParameter  p0_os_Run1("p0_os_Run1",1,0.,0.);
    FitParameter  p1_os_Run1("p1_os_Run1",1,0.,0.);
    FitParameter  delta_p0_os_Run1("delta_p0_os_Run1",1,0.,0.);
    FitParameter  delta_p1_os_Run1("delta_p1_os_Run1",1,0.,0.);
    FitParameter  avg_eta_os_Run1("avg_eta_os_Run1",1,0.,0.);
    FitParameter  tageff_os_Run1("tageff_os_Run1",1,0.,0.);
    FitParameter  tageff_asym_os_Run1("tageff_asym_os_Run1",1,0.,0.);
    FitParameter  p0_ss_Run1("p0_ss_Run1",1,0.,0.);
    FitParameter  p1_ss_Run1("p1_ss_Run1",1,0.,0.);
    FitParameter  delta_p0_ss_Run1("delta_p0_ss_Run1",1,0.,0.);
    FitParameter  delta_p1_ss_Run1("delta_p1_ss_Run1",1,0.,0.);
    FitParameter  avg_eta_ss_Run1("avg_eta_ss_Run1",1,0.,0.);
    FitParameter  tageff_ss_Run1("tageff_ss_Run1",1,0.,0.);
    FitParameter  tageff_asym_ss_Run1("tageff_asym_ss_Run1",1,0.,0.);
    FitParameter  production_asym_Run1("production_asym_Run1",1,0.,0.);
    FitParameter  detection_asym_Run1("detection_asym_Run1",1,0.,0.);
    
    FitParameter  offset_mean_dt_Run2("offset_mean_dt_Run2",1,0,0.1);
    FitParameter  scale_mean_dt_Run2("scale_mean_dt_Run2",1,0,0.1);
    FitParameter  scale_mean_2_dt_Run2("scale_mean_2_dt_Run2",1,0,0.1);
    FitParameter  offset_sigma_dt_Run2("offset_sigma_dt_Run2",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run2("scale_sigma_dt_Run2",1,1.,0.1);
    FitParameter  scale_sigma_2_dt_Run2("scale_sigma_2_dt_Run2",1,0.,0.1);
    FitParameter  offset_sigma2_dt_Run2("offset_sigma2_dt_Run2",1,0.,0.1);
    FitParameter  scale_sigma2_dt_Run2("scale_sigma2_dt_Run2",1,1.,0.1);
    FitParameter  scale_sigma2_2_dt_Run2("scale_sigma2_2_dt_Run2",1,0.,0.1);
    FitParameter  offset_sigma3_dt_Run2("offset_sigma3_dt_Run2",1,0.,0.1);
    FitParameter  scale_sigma3_dt_Run2("scale_sigma3_dt_Run2",1,1.,0.1);
    FitParameter  scale_sigma3_2_dt_Run2("scale_sigma3_2_dt_Run2",1,0.,0.1);
    FitParameter  offset_f_dt_Run2("offset_f_dt_Run2",1,1,0.1);
    FitParameter  scale_f_dt_Run2("scale_f_dt_Run2",1,0.,0.1);
    FitParameter  scale_f_2_dt_Run2("scale_f_2_dt_Run2",1,0.,0.1);
    FitParameter  offset_f2_dt_Run2("offset_f2_dt_Run2",1,0.,0.1);
    FitParameter  scale_f2_dt_Run2("scale_f2_dt_Run2",1,0.,0.1);
    FitParameter  scale_f2_2_dt_Run2("scale_f2_2_dt_Run2",1,0.,0.1);

    FitParameter  p0_os_Run2("p0_os_Run2",1,0.,0.);
    FitParameter  p1_os_Run2("p1_os_Run2",1,0.,0.);
    FitParameter  delta_p0_os_Run2("delta_p0_os_Run2",1,0.,0.);
    FitParameter  delta_p1_os_Run2("delta_p1_os_Run2",1,0.,0.);
    FitParameter  avg_eta_os_Run2("avg_eta_os_Run2",1,0.,0.);
    FitParameter  tageff_os_Run2("tageff_os_Run2",1,0.,0.);
    FitParameter  tageff_asym_os_Run2("tageff_asym_os_Run2",1,0.,0.);
    FitParameter  p0_ss_Run2("p0_ss_Run2",1,0.,0.);
    FitParameter  p1_ss_Run2("p1_ss_Run2",1,0.,0.);
    FitParameter  delta_p0_ss_Run2("delta_p0_ss_Run2",1,0.,0.);
    FitParameter  delta_p1_ss_Run2("delta_p1_ss_Run2",1,0.,0.);
    FitParameter  avg_eta_ss_Run2("avg_eta_ss_Run2",1,0.,0.);
    FitParameter  tageff_ss_Run2("tageff_ss_Run2",1,0.,0.);
    FitParameter  tageff_asym_ss_Run2("tageff_asym_ss_Run2",1,0.,0.);
    FitParameter  production_asym_Run2("production_asym_Run2",1,0.,0.);
    FitParameter  detection_asym_Run2("detection_asym_Run2",1,0.,0.);

    /// Fit parameters per run and trigger cat
    FitParameter  c0_Run1("c0_Run1",1,1,0.1);
    FitParameter  c1_Run1("c1_Run1",1,1,0.1);
    FitParameter  c2_Run1("c2_Run1",1,1,0.1);
    FitParameter  c3_Run1("c3_Run1",1,1,0.1);
    FitParameter  c4_Run1("c4_Run1",1,1,0.1);
    FitParameter  c5_Run1("c5_Run1",1,1,0.1);
    FitParameter  c6_Run1("c6_Run1",1,1,0.1);
    FitParameter  c7_Run1("c7_Run1",1,1,0.1);
    FitParameter  c8_Run1("c8_Run1",1,1,0.1);
    FitParameter  c9_Run1("c9_Run1",1,1,0.1);
    
    FitParameter  c0_Run2("c0_Run2",1,1,0.1);
    FitParameter  c1_Run2("c1_Run2",1,1,0.1);
    FitParameter  c2_Run2("c2_Run2",1,1,0.1);
    FitParameter  c3_Run2("c3_Run2",1,1,0.1);
    FitParameter  c4_Run2("c4_Run2",1,1,0.1);
    FitParameter  c5_Run2("c5_Run2",1,1,0.1);
    FitParameter  c6_Run2("c6_Run2",1,1,0.1);
    FitParameter  c7_Run2("c7_Run2",1,1,0.1);
    FitParameter  c8_Run2("c8_Run2",1,1,0.1);
    FitParameter  c9_Run2("c9_Run2",1,1,0.1);

    //FullTimePdf_mod t_pdf(r,delta,gamma,k);
    string marginalPdfsPrefix = "comb";
    FullTimePdf t_pdf(C, D, D_bar, S, S_bar, k,
                      Gamma, dGamma, dm
                      ,offset_mean_dt,scale_mean_dt,scale_mean_2_dt
                      ,offset_sigma_dt, scale_sigma_dt, scale_sigma_2_dt
                      ,offset_sigma2_dt, scale_sigma2_dt, scale_sigma2_2_dt
                      ,offset_sigma3_dt, scale_sigma3_dt, scale_sigma3_2_dt
                      ,offset_f_dt, scale_f_dt, scale_f_2_dt
                      ,offset_f2_dt, scale_f2_dt, scale_f2_2_dt
                      ,c0, c1, c2 ,c3, c4, c5
                      ,c6, c7, c8, c9,
                      p0_os, p1_os, delta_p0_os, delta_p1_os, 
                      avg_eta_os, tageff_os, tageff_asym_os, 
                      p0_ss, p1_ss, delta_p0_ss, delta_p1_ss, 
                      avg_eta_ss, tageff_ss, tageff_asym_ss, 
                      production_asym, detection_asym, marginalPdfsPrefix );

    /// Simultaneous pdfs
    FullTimePdf t_pdf_Run1(C, D, D_bar, S, S_bar, k,
                      Gamma, dGamma, dm
                      ,offset_mean_dt_Run1,scale_mean_dt_Run1,scale_mean_2_dt_Run1
                      ,offset_sigma_dt_Run1, scale_sigma_dt_Run1, scale_sigma_2_dt_Run1
                      ,offset_sigma2_dt_Run1, scale_sigma2_dt_Run1, scale_sigma2_2_dt_Run1
                      ,offset_sigma3_dt_Run1, scale_sigma3_dt_Run1, scale_sigma3_2_dt_Run1
                      ,offset_f_dt_Run1, scale_f_dt_Run1, scale_f_2_dt_Run1
                      ,offset_f2_dt_Run1, scale_f2_dt_Run1, scale_f2_2_dt_Run1
                      ,c0_Run1, c1_Run1, c2_Run1 ,c3_Run1, c4_Run1, c5_Run1
                      ,c6_Run1, c7_Run1, c8_Run1, c9_Run1,
                      p0_os_Run1, p1_os_Run1, delta_p0_os_Run1, delta_p1_os_Run1, 
                      avg_eta_os_Run1, tageff_os_Run1, tageff_asym_os_Run1, 
                      p0_ss_Run1, p1_ss_Run1, delta_p0_ss_Run1, delta_p1_ss_Run1, 
                      avg_eta_ss_Run1, tageff_ss_Run1, tageff_asym_ss_Run1, 
                      production_asym_Run1, detection_asym_Run1, "Run1" );
    
    
    FullTimePdf t_pdf_Run2(C, D, D_bar, S, S_bar, k,
                              Gamma, dGamma, dm
			      ,offset_mean_dt_Run2,scale_mean_dt_Run2,scale_mean_2_dt_Run2
			      ,offset_sigma_dt_Run2, scale_sigma_dt_Run2, scale_sigma_2_dt_Run2
			      ,offset_sigma2_dt_Run2, scale_sigma2_dt_Run2, scale_sigma2_2_dt_Run2
			      ,offset_sigma3_dt_Run2, scale_sigma3_dt_Run2, scale_sigma3_2_dt_Run2
			      ,offset_f_dt_Run2, scale_f_dt_Run2, scale_f_2_dt_Run2
			      ,offset_f2_dt_Run2, scale_f2_dt_Run2, scale_f2_2_dt_Run2
                              ,c0_Run2, c1_Run2, c2_Run2 ,c3_Run2, c4_Run2, c5_Run2
                              ,c6_Run2, c7_Run2, c8_Run2, c9_Run2,
                              p0_os_Run2, p1_os_Run2, delta_p0_os_Run2, delta_p1_os_Run2, 
                              avg_eta_os_Run2, tageff_os_Run2, tageff_asym_os_Run2, 
                              p0_ss_Run2, p1_ss_Run2, delta_p0_ss_Run2, delta_p1_ss_Run2, 
                              avg_eta_ss_Run2, tageff_ss_Run2, tageff_asym_ss_Run2, 
                              production_asym_Run2, detection_asym_Run2, "Run2" );
    
    /// Randomize start vals
    MinuitParameterSet* mps = MinuitParameterSet::getDefaultSet();
    if(randomizeStartVals && (string)channel=="signal"){
	double val = gRandom->Uniform(0,1);
	mps->getParPtr("C")->setCurrentFitVal(val);
	((FitParameter*)mps->getParPtr("C"))->setInit(val);	

	val = gRandom->Uniform(-1,1);
	mps->getParPtr("D")->setCurrentFitVal(val);
	((FitParameter*)mps->getParPtr("D"))->setInit(val);	

	val = gRandom->Uniform(-1,1);
	mps->getParPtr("D_bar")->setCurrentFitVal(val);
	((FitParameter*)mps->getParPtr("D_bar"))->setInit(val);	

	val = gRandom->Uniform(-1,1);
	mps->getParPtr("S")->setCurrentFitVal(val);
	((FitParameter*)mps->getParPtr("S"))->setInit(val);	

	val = gRandom->Uniform(-1,1);
	mps->getParPtr("S_bar")->setCurrentFitVal(val);
	((FitParameter*)mps->getParPtr("S_bar"))->setInit(val);	
    }

    /// Load data
    double t,dt,mB;
    int f;
    int q_OS;
    double Bs_ID,Ds_ID;
    int Bs_TRUEID,bkgCAT;
    Int_t q_SS;
    double eta_OS;
    Double_t eta_SS;
    double sw;
    int year,run,Ds_finalState,trigger;
    double Ks[4];
    double pi[4];
    double Dm[4];
    int KsCat;
    
    TChain* tree_norm=new TChain("DecayTree");
	if(mode == "fit" && doToyStudy == 1){
		if(addBkgToToys)tree_norm->Add(((string)OutputDir+"sw_toys_"+anythingToString((int)step)+".root").c_str());
		else tree_norm->Add(((string)OutputDir+"toys_"+anythingToString((int)step)+".root").c_str());
	}
	else tree_norm->Add(((string)InputFileName).c_str());
	tree_norm->SetBranchStatus("*",0);
	tree_norm->SetBranchStatus("*sw*",1);
	tree_norm->SetBranchStatus("year",1);
	tree_norm->SetBranchStatus("*DEC",1);
	tree_norm->SetBranchStatus("*PROB",1);
	tree_norm->SetBranchStatus("*OS*",1);
	tree_norm->SetBranchStatus("*TAU*",1);
	tree_norm->SetBranchStatus("*ID*",1);
	tree_norm->SetBranchStatus("weight",1);
	tree_norm->SetBranchStatus("TriggerCat",1);
	tree_norm->SetBranchStatus("run",1);
	tree_norm->SetBranchStatus("KsCat",1);
	tree_norm->SetBranchStatus("FullDTF_*P*",1);
	tree_norm->SetBranchStatus("B_DTF_MM",1);
	if(useTrueTagging)tree_norm->SetBranchStatus("B_TRUEID",1);
	if(useTrueTagging)tree_norm->SetBranchStatus("bkgCAT",1);
	
	tree_norm->SetBranchAddress("B_DTF_MM",&mB);
	if(useTrueTau)tree_norm->SetBranchAddress("B_TRUETAU",&t);
	else tree_norm->SetBranchAddress("B_FullDTF_TAU",&t);
	tree_norm->SetBranchAddress("B_FullDTF_TAUERR",&dt);
	if(useTrueTagging)tree_norm->SetBranchAddress("B_TRUEID",&Bs_TRUEID);
	if(useTrueTagging)tree_norm->SetBranchAddress("bkgCAT",&bkgCAT);
	tree_norm->SetBranchAddress("D_ID",&f);
	tree_norm->SetBranchAddress("OS_Combination_DEC",&q_OS);
	tree_norm->SetBranchAddress("OS_Combination_PROB",&eta_OS);
	tree_norm->SetBranchAddress("SS_Combination_DEC",&q_SS);
	tree_norm->SetBranchAddress("SS_Combination_PROB",&eta_SS);
	tree_norm->SetBranchAddress(((string)weightName).c_str(),&sw);
	tree_norm->SetBranchAddress("year",&year);
	tree_norm->SetBranchAddress("run",&run);
	tree_norm->SetBranchAddress("TriggerCat",&trigger);
    tree_norm->SetBranchAddress("KsCat",&KsCat);
    tree_norm->SetBranchAddress("FullDTF_Ks_PX",&Ks[0]);
    tree_norm->SetBranchAddress("FullDTF_Ks_PY",&Ks[1]);
    tree_norm->SetBranchAddress("FullDTF_Ks_PZ",&Ks[2]);
    tree_norm->SetBranchAddress("FullDTF_Ks_PE",&Ks[3]);
    tree_norm->SetBranchAddress("FullDTF_D_PX",&Dm[0]);
    tree_norm->SetBranchAddress("FullDTF_D_PY",&Dm[1]);
    tree_norm->SetBranchAddress("FullDTF_D_PZ",&Dm[2]);
    tree_norm->SetBranchAddress("FullDTF_D_PE",&Dm[3]);
    tree_norm->SetBranchAddress("FullDTF_pi_PX",&pi[0]);
    tree_norm->SetBranchAddress("FullDTF_pi_PY",&pi[1]);
    tree_norm->SetBranchAddress("FullDTF_pi_PZ",&pi[2]);
    tree_norm->SetBranchAddress("FullDTF_pi_PE",&pi[3]);
    
    DalitzEventList eventList,eventList_Run1,eventList_Run2;
    int N_sample = tree_norm->GetEntries();
    if(N_bootstrap == -1)N_bootstrap = N_sample;

    vector<int> b_indices;
    while( b_indices.size() < N_bootstrap )b_indices.push_back(TMath::Nint(ranLux.Uniform(0,N_sample-1)));
    sort(b_indices.begin(), b_indices.end());
    if(doBootstrap)N_sample = b_indices.size();

    TRandom3 rndm;
    TRandom3 randRes(seed);
    int badEvents = 0;
    
    for(int i=0; i< N_sample; i++)
    {	
        if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << N_sample << endl;
        if(doBootstrap) tree_norm->GetEntry(b_indices[i]);
        else tree_norm->GetEntry(i);
        
        if(useTrueTau)t= t*1000.;//+ gRandom->Gaus(0.,offset_sigma_dt_Run2_17+scale_sigma_dt_Run2_17*dt);

        double sign = 1.;
        TLorentzVector Ks_p(sign*Ks[0],sign*Ks[1],sign*Ks[2],Ks[3]);
        TLorentzVector D_p(sign*Dm[0],sign*Dm[1],sign*Dm[2],Dm[3]);
        TLorentzVector pi_p(sign*pi[0],sign*pi[1],sign*pi[2],pi[3]);
        TLorentzVector B_p = Ks_p + pi_p + D_p;
        // array of vectors
        vector<TLorentzVector> vectorOfvectors;        
        vectorOfvectors.push_back(B_p*MeV);
        vectorOfvectors.push_back(D_p*MeV);
        vectorOfvectors.push_back(Ks_p*MeV);
        vectorOfvectors.push_back(pi_p*MeV);
        
        DalitzEvent evt;
        if(f < 0)evt = DalitzEvent(pat, vectorOfvectors);
        else evt = DalitzEvent(pat_CP, vectorOfvectors);
        if( evt.s(1,3) < pat.sijMin(1,3) || evt.s(1,3) > pat.sijMax(1,3) || evt.s(1,2) < pat.sijMin(1,2) || evt.s(1,2) > pat.sijMax(1,2)  || evt.s(2,3) < pat.sijMin(2,3) || evt.s(2,3) > pat.sijMax(2,3) || TMath::IsNaN(dt) || TMath::IsNaN(t) )
        {
            badEvents++;
            continue;
        }
        if(abs(sqrt(evt.s(2,3))-1869.61)<20)continue;
        if(abs(sqrt(evt.s(2,3))-1968.30)<20)continue;        

        if(useTrueTagging){
            if(bkgCAT==1)continue;
            if(Bs_TRUEID == 531){ 
                q_OS=1;
                q_SS=1;
            }
            else if(Bs_TRUEID == -531) {
                q_OS = -1;
                q_SS = -1;
            }
            else continue;
            eta_OS = 0.;
            eta_SS = 0.;
        }
        if(t < min_TAU || t > max_TAU )continue;
        if( dt < min_TAUERR || dt > max_TAUERR )continue;
        if(year < min_year || year > max_year) continue;
        if(trigger < min_trigger || trigger > max_trigger) continue;
        if(eta_SS < w_min || eta_SS > w_max )continue;
        if((string)channel=="norm" && KsCat==1)continue;
        
        evt.setWeight(sw);
        evt.setValueInVector(0, t);
        evt.setValueInVector(1, dt);   
        if(f<0)evt.setValueInVector(2, 1);
        else if(f > 0)evt.setValueInVector(2, -1);
        else {
            cout << "ERROR:: Undefined final state " << f << endl;  
            throw "ERROR";
        }
        evt.setValueInVector(3, q_OS);
        evt.setValueInVector(4, eta_OS);
        evt.setValueInVector(5, q_SS);
        evt.setValueInVector(6, eta_SS);
        evt.setValueInVector(7, run);
        evt.setValueInVector(8, trigger);
        evt.setValueInVector(9, mB);
        eventList.Add(evt);

        if(run == 1) eventList_Run1.Add(evt);
        else if(run == 2) eventList_Run2.Add(evt);
    }
    cout << "Bad events = " << badEvents << " (" << badEvents/(double)eventList.size() * 100. << " %)" << endl;
    
    /// Fit with MINT Pdf
    Neg2LL neg2LL(t_pdf, eventList);    
    Neg2LL neg2LL_Run1(t_pdf_Run1, eventList_Run1);
    Neg2LL neg2LL_Run2(t_pdf_Run2, eventList_Run2);
  
    Neg2LLSum neg2LL_sim;
    if(eventList_Run1.size()>0)neg2LL_sim.add(&neg2LL_Run1);
    if(eventList_Run2.size()>0)neg2LL_sim.add(&neg2LL_Run2);

     Neg2LLMultiConstraint constrains_tagging_Run1(MinuitParameterSet::getDefaultSet(),"_Tagging_Run1");
     Neg2LLMultiConstraint constrains_tagging_Run2(MinuitParameterSet::getDefaultSet(),"_Tagging_Run2");
     if(useGaussConstrainsTagging){
 	    neg2LL_sim.add(&constrains_tagging_Run1);
 	    neg2LL_sim.add(&constrains_tagging_Run2);
     }

    Neg2LLMultiConstraint constrains_Acc(MinuitParameterSet::getDefaultSet(),"_Acc");
    if(doAccSystematics && mode == "fit"){
	if(useCholDec){
		if(chol_index > constrains_Acc.getNumberParams()-1){
			cout << "ERROR:: Invalid cholesky index ! "<< endl;
			throw "ERROR";
		}
		constrains_Acc.smearInputValuesChol(chol_index,(step-1) - chol_index * varPerParChol, 1);
	}
	else
	{ 
		constrains_Acc.smearInputValues();
	}
    }

    Neg2LLMultiConstraint constrains_sys(MinuitParameterSet::getDefaultSet(),("_" + (string)doSystematic).c_str());
    if((string)doSystematic != "" && mode == "fit")constrains_sys.smearInputValues();

    Neg2LLMultiConstraint constrains_bias(MinuitParameterSet::getDefaultSet(),"_bias");
    if(useGaussConstrainsBias){
 	    neg2LL_sim.add(&constrains_bias);
     }

    Minimiser mini;
    if(doSimFit)mini.attachFunction(&neg2LL_sim);
    else mini.attachFunction(&neg2LL);
    if(mode == "fit"){
        mini.doFit();
    	mini.printResultVsInput();
    }

    /// Plot
    TCanvas* c = new TCanvas();
        
    TH1D* h_t = new TH1D("h_t",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_mixed = new TH1D("h_t_mixed",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed = new TH1D("h_t_unmixed",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_untagegged = new TH1D("h_t_untagegged",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_t_OS = new TH1D("h_t_OS",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_mixed_OS = new TH1D("h_t_mixed_OS",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed_OS = new TH1D("h_t_unmixed_OS",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_SS = new TH1D("h_t_SS",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_mixed_SS = new TH1D("h_t_mixed_SS",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed_SS = new TH1D("h_t_unmixed_SS",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_t_mp = new TH1D("h_t_mp",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_0p = new TH1D("h_t_0p",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_pp = new TH1D("h_t_pp",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_mm = new TH1D("h_t_mm",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_0m = new TH1D("h_t_0m",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_pm = new TH1D("h_t_pm",";t (ps);Yield  (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_t_p = new TH1D("h_t_p",";t (ps);Yield  (a.u.) ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_t_m = new TH1D("h_t_m",";t (ps);Yield  (a.u.) ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_t_N = new TH1D("h_t_N",";t (ps);Yield  (a.u.) ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_t_Nbar = new TH1D("h_t_Nbar",";t (ps);Yield  (a.u.) ",nBinsAsym,min_TAU,max_TAU);

    TH1D* h_dt = new TH1D("h_dt",";#sigma_{t} (ps);Yield  (a.u.) ",nBinst,0,0.15);
    TH1D* h_eta_OS = new TH1D("h_eta_OS",";#eta_{OS};Yield  (a.u.) ",nBinst,0,0.5);
    TH1D* h_eta_SS = new TH1D("h_eta_SS",";#eta_{SS};Yield  (a.u.) ",nBinst,0,0.5);
    TH1D* h_q_OS = new TH1D("h_q_OS",";q_{OS};Yield  (a.u.) ",3,-1.5,1.5);
    TH1D* h_q_SS = new TH1D("h_q_SS",";q_{SS};Yield  (a.u.) ",3,-1.5,1.5);
    TH1D* h_f = new TH1D("h_f",";q_{f};Yield  (a.u.) ",2,-2,2);

    TH1D* h_N_mixed = new TH1D("h_N_mixed",";t modulo (2#pi/#Deltam_{s}) (ps); A_{mix} ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_N_unmixed = (TH1D*) h_N_mixed->Clone("h_N_unmixed");
    TH1D* h_N_mixed_OS = new TH1D("h_N_mixed_OS",";t modulo (2#pi/#Deltam_{s}) (ps); A_{mix} ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_N_unmixed_OS = (TH1D*) h_N_mixed->Clone("h_N_unmixed_OS");
    TH1D* h_N_mixed_SS = new TH1D("h_N_mixed_SS",";t modulo (2#pi/#Deltam_{s}) (ps); A_{mix} ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_N_unmixed_SS = (TH1D*) h_N_mixed->Clone("h_N_unmixed_SS");

    TH1D* h_N_mixed_p = (TH1D*) h_N_mixed->Clone("h_N_mixed_p");
    TH1D* h_N_unmixed_p = (TH1D*) h_N_mixed->Clone("h_N_unmixed_p");
    TH1D* h_N_mixed_m = (TH1D*) h_N_mixed->Clone("h_N_mixed_m");
    TH1D* h_N_unmixed_m = (TH1D*) h_N_mixed->Clone("h_N_unmixed_m");

    TH1D* h_N_mixed_p_unfolded = new TH1D("h_N_mixed_p_unfolded",";t (ps);A_{CP}(t) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_N_unmixed_p_unfolded = (TH1D*) h_N_mixed_p_unfolded->Clone("h_N_unmixed_p_fit");
    TH1D* h_N_mixed_m_unfolded = (TH1D*) h_N_mixed_p_unfolded->Clone("h_N_mixed_m_fit");
    TH1D* h_N_unmixed_m_unfolded = (TH1D*) h_N_mixed_p_unfolded->Clone("h_N_unmixed_m_fit");

    double N_OS = 0;
    double N_SS = 0;
    double N_OS_SS = 0;
    double N = 0;

    double w_OS = 0;
    double w_SS = 0;
    double w_OS_SS = 0;
 	
    double D_OS = 0;
    double D_SS = 0;
    double D_OS_SS = 0;
    double D_comb = 0;

    double N_OS_all = 0;
    double N_SS_all = 0;
    double w_OS_all = 0;
    double w_SS_all = 0;
    double D_OS_all = 0;
    double D_SS_all = 0;

    double N_Run1 = 0;
    double N_Run2 = 0;    
    double N_Run1_t0 = 0;
    double N_Run1_t1 = 0;
    double N_Run2_t0 = 0;
    double N_Run2_t1 = 0;

    double sigma_t_eff = 0;

    for (unsigned int i=0; i<eventList.size(); i++) {
 
        N += eventList[i].getWeight();
        if(!doSimFit)sigma_t_eff += t_pdf.getCalibratedResolution(eventList[i].getValueFromVector(1)) * eventList[i].getWeight();

        h_t->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
        h_dt->Fill(eventList[i].getValueFromVector(1),eventList[i].getWeight());
        if(eventList[i].getValueFromVector(3) != 0)h_eta_OS->Fill(eventList[i].getValueFromVector(4),eventList[i].getWeight());
        if(eventList[i].getValueFromVector(5) != 0)h_eta_SS->Fill(eventList[i].getValueFromVector(6),eventList[i].getWeight());

        int f_evt = eventList[i].getValueFromVector(2);
        int q1 = eventList[i].getValueFromVector(3);
        int q2 = eventList[i].getValueFromVector(5);   
        int q_eff = 0;
        double w_eff = 0.5;
        int run_evt = eventList[i].getValueFromVector(7);   
        int trigger_evt = eventList[i].getValueFromVector(8);   
        
        h_q_OS->Fill(q1,eventList[i].getWeight());
        h_q_SS->Fill(q2,eventList[i].getWeight());
        h_f->Fill(f_evt,eventList[i].getWeight());

        if( q1 != 0)h_t_OS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
        if( q2 != 0)h_t_SS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());

        if(run_evt==1 && trigger_evt == 0)N_Run1_t0 += eventList[i].getWeight();
        else if(run_evt==1 && trigger_evt == 1)N_Run1_t1 += eventList[i].getWeight();
        else if(run_evt==2 && trigger_evt == 0)N_Run2_t0 += eventList[i].getWeight();
        else if(run_evt==2 && trigger_evt == 1)N_Run2_t1 += eventList[i].getWeight();

        std::pair<double, double> calibrated_mistag_os;
        std::pair<double, double> calibrated_mistag_ss;
        if(doSimFit){
            if(run_evt==1){
                calibrated_mistag_os = t_pdf_Run1.getCalibratedMistag_OS(eventList[i]);
                calibrated_mistag_ss = t_pdf_Run1.getCalibratedMistag_SS(eventList[i]);
            }
            else{
                calibrated_mistag_os = t_pdf_Run2.getCalibratedMistag_OS(eventList[i]);
                calibrated_mistag_ss = t_pdf_Run2.getCalibratedMistag_SS(eventList[i]);                
            }
        }

        double p = ( (1.-q1)/2. + q1 * (1.- calibrated_mistag_os.first )) * ( (1.-q2)/2. + q2 * (1.- calibrated_mistag_ss.first ));
        double p_bar = ( (1.+q1)/2. - q1 * (1.- calibrated_mistag_os.second )) * ( (1.+q2)/2. - q2 * (1.- calibrated_mistag_ss.second ));
            
        if( p/(p+p_bar) > 0.5 ){ 
            q_eff = 1;
            w_eff = 1-p/(p+p_bar);
        }
        else if( p/(p+p_bar) < 0.5 ){
            q_eff = -1;
            w_eff = p/(p+p_bar);
        }

        if(q1 != 0 && q2 != 0){
            if(q_eff != 0){
                N_OS_SS += eventList[i].getWeight();
                w_OS_SS += w_eff * eventList[i].getWeight();
                D_OS_SS += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
            }
        }
        else if( q1 != 0){
                //q_eff = q1;  flip tag ???
                N_OS += eventList[i].getWeight();
                w_OS += w_eff * eventList[i].getWeight(); 
                D_OS += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
        }
        else if( q2 != 0){
                //q_eff = q2;
                N_SS += eventList[i].getWeight();
	    	w_SS += w_eff * eventList[i].getWeight(); 
                D_SS += pow(1.-2.*w_eff,2)* eventList[i].getWeight(); 
        } 

        D_comb += pow(1.-2.*w_eff,2)* eventList[i].getWeight();

	double D_res = 1.;
        if(doSimFit){
            if(run_evt==1){
                D_res = exp(-pow(t_pdf_Run1.getCalibratedResolution(eventList[i].getValueFromVector(1))*dm,2)/2.);
            }
            else if(run_evt==2){
                D_res = exp(-pow(t_pdf_Run2.getCalibratedResolution(eventList[i].getValueFromVector(1))*dm,2)/2.);
            }
	}

        double D_tag = 0.;
        if(q_eff != 0) D_tag = (1.-2.*abs(w_eff));
        double D_tot = D_tag * D_res;
        if(!dilutionWeight)D_tot = 1.;
        
        if(q1 != 0) N_OS_all += eventList[i].getWeight();
        if(q2 != 0){
                //if(q2 > 0 && calibrated_mistag_ss.first < 0.5) 
		N_SS_all += eventList[i].getWeight();
                //else if(q2 < 0 && calibrated_mistag_ss.second < 0.5) N_SS_all += eventList[i].getWeight();
        }    
            
        if(q1>0){
			w_OS_all +=  calibrated_mistag_os.first * eventList[i].getWeight();
			D_OS_all +=  pow(1.-2.*calibrated_mistag_os.first,2)* eventList[i].getWeight();
        } 
        else if(q1<0){
			w_OS_all +=  calibrated_mistag_os.second * eventList[i].getWeight();	
			D_OS_all +=  pow(1.-2.*calibrated_mistag_os.second,2)* eventList[i].getWeight();
        }

        if(q2>0){
			w_SS_all +=  calibrated_mistag_ss.first * eventList[i].getWeight();
			D_SS_all +=  pow(1.-2.*calibrated_mistag_ss.first,2)* eventList[i].getWeight();
        } 
        else if(q2<0){
			w_SS_all +=  calibrated_mistag_ss.second * eventList[i].getWeight();	
			D_SS_all +=  pow(1.-2.*calibrated_mistag_ss.second,2)* eventList[i].getWeight();
        }
   
        if((string)channel=="signal"){

            if(f_evt == 1){ 
                h_t_p->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
            }
            else h_t_m->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
        
            if(q_eff == 1)h_t_N->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
            else if(q_eff == -1)h_t_Nbar->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);

                if(q_eff==-1 && f_evt == 1){ 
                h_t_mp->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                if(w_eff<w_max){
                    h_N_mixed_p->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                    h_N_mixed_p_unfolded->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                }

                }
                else if(q_eff==0 && f_evt == 1)h_t_0p->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                else if(q_eff==1 && f_evt == 1){
                            h_t_pp->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                                if(w_eff<w_max){
                                        h_N_unmixed_p->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                                        h_N_unmixed_p_unfolded->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                                }
                }
                else if(q_eff==-1 && f_evt == -1){
                        h_t_mm->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                        if(w_eff<w_max){
                                h_N_unmixed_m->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                                h_N_unmixed_m_unfolded->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                        }
               }
               else if(q_eff==0 && f_evt == -1)h_t_0m->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
               else if(q_eff==1 && f_evt == -1){
                        h_t_pm->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                        if(w_eff<w_max){
                            h_N_mixed_m->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                            h_N_mixed_m_unfolded->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                        }
                    }
        }
//         else { 	
            if(q_eff == 0)h_t_untagegged->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight());
                else if(q_eff*f_evt > 0  ){
                        h_t_mixed->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                        //if(w_eff<w_max)
			h_N_mixed->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);

                        if(q1 != 0)h_t_mixed_OS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
			if(q1 != 0)h_N_mixed_OS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);

                        if(q2 != 0)h_t_mixed_SS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
			if(q2 != 0)h_N_mixed_SS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                    }
                else {
                    h_t_unmixed->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                    //if(w_eff<w_max)
		    h_N_unmixed->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);

                    if(q1 != 0)h_t_unmixed_OS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
		    if(q1 != 0)h_N_unmixed_OS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);

                    if(q2 != 0)h_t_unmixed_SS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
	            if(q2 != 0)h_N_unmixed_SS->Fill(eventList[i].getValueFromVector(0),eventList[i].getWeight()*D_tot);
                }
//         }
  
    }     

    cout << "tree size = " << eventList.size() << endl;
    cout << "sumw = " << N << endl << endl;
    cout << "N_Run1_t0 =" << N_Run1_t0/N <<  endl;
    cout << "N_Run1_t1 =" << N_Run1_t1/N <<  endl;
    cout << "N_Run2_t0 =" << N_Run2_t0/N <<  endl;
    cout << "N_Run2_t1 =" << N_Run2_t1/N <<  endl;
    cout << "sigma_t_eff = " << sigma_t_eff/N << endl << endl;
    cout << "N_OS = " << N_OS_all << endl;
    cout << "N_SS = " << N_SS_all << endl << endl;
    cout << "eff_OS = " <<(N_OS_all)/N << endl;
    cout << "eff_SS = " <<(N_SS_all)/N << endl;
        
    cout << "Tagging perfromance " << endl << endl;        
    cout << "Tagger | eff_tag | <w> | e_eff " <<  endl;
    cout << "OS  | " << (N_OS+N_OS_SS)/N << " | " <<  (w_OS_all)/(N_OS+N_OS_SS) << " | " << D_OS_all/N << endl;
    cout << "SS  | " << (N_SS+N_OS_SS)/N << " | " <<  (w_SS_all)/(N_SS+N_OS_SS) << " | " << D_SS_all/N << endl << endl;
    cout << "OS only  | " << N_OS/N << " | " <<  w_OS/N_OS << " | " << N_OS/N * D_OS/N_OS << endl;
    cout << "SS only  | " << N_SS/N << " | " <<  w_SS/N_SS << " | " << N_SS/N * D_SS/N_SS << endl;
    cout << "OS+SS    | " << N_OS_SS/N << " | " <<  w_OS_SS/N_OS_SS << " | " << N_OS_SS/N * D_OS_SS/N_OS_SS << endl;
    cout << "Combined | " << (N_OS+N_SS+N_OS_SS)/N << " | "<<  (w_OS+w_SS+w_OS_SS)/(N_OS+N_SS+N_OS_SS) << " | " << (N_OS+N_SS+N_OS_SS)/N * D_comb/(N_OS+N_SS+N_OS_SS) << endl << endl ;

    N_Run1 = N_Run1_t0 + N_Run1_t1;
    N_Run2 = N_Run2_t0 + N_Run2_t1;
    /// Generate toys 
    DalitzEventList toys;
    if(mode == "gen"){
        if(doSimFit) {
            toys.Add(t_pdf_Run1.generateToys(N_scale_toys * N_Run1,1,0));
            toys.Add(t_pdf_Run2.generateToys(N_scale_toys *N_Run2,2,0));
        }
        else  toys.Add(t_pdf.generateToys(N_scale_toys *N));
        
        t_pdf.saveEventListToFile(toys,((string)OutputDir+"toys_"+anythingToString((int)step)+".root").c_str());

        return;
    }

    /// Save results 
    TMatrixTSym<double> cov_full = mini.covMatrixFull();     
    //cov_full.Print();
    ofstream fitResults;
    fitResults.open(((string)OutputDir+"FitResults.txt").c_str(),std::ofstream::trunc);

    vector<string> prefix;
    if(doSimFit){
		prefix.push_back("_Run1");
    		prefix.push_back("_Run2");
    }
    else prefix.push_back("");

    for(int p = 0 ; p < prefix.size(); p++){
	
	fitResults << "\"" << "ConstrainMulti_Tagging" << prefix[p] << "\"" << "     " << "\"" ;
	for(int i = 0 ; i < mps->size(); i++){	
		if(mps->getParPtr(i)->iFixInit() == 0 && A_is_in_B(prefix[p],mps->getParPtr(i)->name()))
			fitResults << mps->getParPtr(i)->name() << " ";
	}	
	fitResults << "\"";
	fitResults<< "\n";
	fitResults<< "\n";
	fitResults << "\"" << "ConstrainMulti_Tagging" << prefix[p] << "_corr" << "\"" << "     " << "\"" ;
	for(int i= 0; i< cov_full.GetNcols(); i++){
			if(!(mps->getParPtr(i)->iFixInit() == 0 && A_is_in_B(prefix[p],mps->getParPtr(i)->name())))continue;
			for(int j= i; j< cov_full.GetNcols(); j++){
				if(!(mps->getParPtr(j)->iFixInit() == 0 && A_is_in_B(prefix[p],mps->getParPtr(j)->name())))continue;
				if(i == j) fitResults << 1. << " " ;
				else if(cov_full[i][i] == 0. || cov_full[j][j] == 0.) fitResults << 0. << " " ;
				//else if(cov_full[i][j]/sqrt(cov_full[i][i])/sqrt(cov_full[j][j]) < 0.00001) fitResults << 0. << " " ; 
				else fitResults << (cov_full[i][j]/sqrt(cov_full[i][i])/sqrt(cov_full[j][j])) << " ";
			}
	}
	fitResults << "\"";
	fitResults<< "\n";
	fitResults<< "\n";
    }

     for(int i = 0 ; i < mps->size(); i++)	    
		fitResults << "\"" << mps->getParPtr(i)->name() << "\"" << "    " << mps->getParPtr(i)->iFixInit() << "    " << mps->getParPtr(i)->mean() << "    " << mps->getParPtr(i)->err() <<  "\n";
     fitResults.close();

    
    /// Save pulls
    gDirectory->cd();
    string paraFileName = (string)OutputDir+"pull_"+ (string)doSystematic+ "_" + anythingToString((int)step)+".root";
    if(doAccSystematics) paraFileName = (string)OutputDir+"pullAcc_"+anythingToString((int)step)+".root";
    if(doAccSystematics && useCholDec) paraFileName = (string)OutputDir+"pullAccChol_"+anythingToString((int)step)+".root";
    cout << paraFileName << endl;


    TFile* paraFile = new TFile( paraFileName.c_str(), "RECREATE");
    paraFile->cd();
    TNtupleD* ntp=0;
    mps->fillNtp(paraFile, ntp);
    ntp->AutoSave();
    paraFile->Close();
    delete paraFile;

    /// Create tagging perfromance tables
    /*
    if(updateAnaNote){

	ofstream resultsfile;
	resultsfile.open(("../../../../../TD-AnaNote/latex/tables/timeFit/"+(string)OutputDir+"result.tex").c_str(),std::ofstream::trunc);
	resultsfile << "\\begin{table}[h]" << "\n";
	resultsfile << "\\centering" << "\n";
// 	resultsfile << "\\small" << "\n";
	resultsfile << "\\caption{Result of the phase-space integrated fit to "; 
	if((string)channel == "norm")resultsfile << "$B_s \\to D_s \\pi \\pi \\pi$";
	else if((string)channel == "signal")resultsfile << "$B_s \\to D_s K \\pi \\pi$";
	resultsfile << " data.}\n";
	resultsfile << "\\begin{tabular}{c c c}" << "\n";
	resultsfile << "\\hline" << "\n";
	resultsfile << "\\hline" << "\n";
	resultsfile << "& Fit parameter & Value \\\\" << "\n";
	resultsfile << "\\hline" << "\n";

	if((string)channel == "norm"){
		resultsfile << std::fixed << std::setprecision(4) 
		<< "Run-I & $p_{0}^{\\text{OS}}$ & " <<  mps->getParPtr("p0_os_Run1")->mean() << " $\\pm$ " << mps->getParPtr("p0_os_Run1")->err() << "\\\\" << "\n"
		<< "&$p_{1}^{\\text{OS}}$  & " <<  mps->getParPtr("p1_os_Run1")->mean() << " $\\pm$ " << mps->getParPtr("p1_os_Run1")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{0}^{\\text{OS}}$  & " <<  mps->getParPtr("delta_p0_os_Run1")->mean() << " $\\pm$ " << mps->getParPtr("delta_p0_os_Run1")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{1}^{\\text{OS}}$  & " <<  mps->getParPtr("delta_p1_os_Run1")->mean() << " $\\pm$ " << mps->getParPtr("delta_p1_os_Run1")->err() << "\\\\" << "\n"
		<< "&$\\epsilon_{tag}^{\\text{OS}}$  & " <<  mps->getParPtr("tageff_os_Run1")->mean() << " $\\pm$ " << mps->getParPtr("tageff_os_Run1")->err() << "\\\\" << "\n"
		<< "&$\\Delta\\epsilon_{tag}^{\\text{OS}}$  & " <<  mps->getParPtr("tageff_asym_os_Run1")->mean() << " $\\pm$ " << mps->getParPtr("tageff_asym_os_Run1")->err() << "\\\\" << "\n"
		
		<< "& $p_{0}^{\\text{SS}}$ & " <<  mps->getParPtr("p0_ss_Run1")->mean() << " $\\pm$ " << mps->getParPtr("p0_ss_Run1")->err() << "\\\\" << "\n"
		<< "&$p_{1}^{\\text{SS}}$  & " <<  mps->getParPtr("p1_ss_Run1")->mean() << " $\\pm$ " << mps->getParPtr("p1_ss_Run1")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{0}^{\\text{SS}}$  & " <<  mps->getParPtr("delta_p0_ss_Run1")->mean() << " $\\pm$ " << mps->getParPtr("delta_p0_ss_Run1")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{1}^{\\text{SS}}$  & " <<  mps->getParPtr("delta_p1_ss_Run1")->mean() << " $\\pm$ " << mps->getParPtr("delta_p1_ss_Run1")->err() << "\\\\" << "\n"
		<< "&$\\epsilon_{tag}^{\\text{SS}}$  & " <<  mps->getParPtr("tageff_ss_Run1")->mean() << " $\\pm$ " << mps->getParPtr("tageff_ss_Run1")->err() << "\\\\" << "\n"
		<< "&$\\Delta\\epsilon_{tag}^{\\text{SS}}$  & " <<  mps->getParPtr("tageff_asym_ss_Run1")->mean() << " $\\pm$ " << mps->getParPtr("tageff_asym_ss_Run1")->err() << "\\\\" << "\n"
	
		<< "&$A_{p}$ & " <<  mps->getParPtr("production_asym_Run1")->mean() << " $\\pm$ " << mps->getParPtr("production_asym_Run1")->err() << "\\\\" << "\n";
	
		resultsfile << "\\\\" << "\n" ;
		resultsfile << std::fixed << std::setprecision(4) 
		<< "Run-II & $p_{0}^{\\text{OS}}$  & " <<  mps->getParPtr("p0_os_Run2")->mean() << " $\\pm$ " << mps->getParPtr("p0_os_Run2")->err() << "\\\\" << "\n"
		<< "&$p_{1}^{\\text{OS}}$  & " <<  mps->getParPtr("p1_os_Run2")->mean() << " $\\pm$ " << mps->getParPtr("p1_os_Run2")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{0}^{\\text{OS}}$  & " <<  mps->getParPtr("delta_p0_os_Run2")->mean() << " $\\pm$ " << mps->getParPtr("delta_p0_os_Run2")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{1}^{\\text{OS}}$  & " <<  mps->getParPtr("delta_p1_os_Run2")->mean() << " $\\pm$ " << mps->getParPtr("delta_p1_os_Run2")->err() << "\\\\" << "\n"
		<< "&$\\epsilon_{tag}^{\\text{OS}}$  & " <<  mps->getParPtr("tageff_os_Run2")->mean() << " $\\pm$ " << mps->getParPtr("tageff_os_Run2")->err() << "\\\\" << "\n"
		<< "&$\\Delta\\epsilon_{tag}^{\\text{OS}}$  & " <<  mps->getParPtr("tageff_asym_os_Run2")->mean() << " $\\pm$ " << mps->getParPtr("tageff_asym_os_Run2")->err() << "\\\\" << "\n"
	
		<< "& $p_{0}^{\\text{SS}}$  & " <<  mps->getParPtr("p0_ss_Run2")->mean() << " $\\pm$ " << mps->getParPtr("p0_ss_Run2")->err() << "\\\\" << "\n"
		<< "&$p_{1}^{\\text{SS}}$  & " <<  mps->getParPtr("p1_ss_Run2")->mean() << " $\\pm$ " << mps->getParPtr("p1_ss_Run2")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{0}^{\\text{SS}}$  & " <<  mps->getParPtr("delta_p0_ss_Run2")->mean() << " $\\pm$ " << mps->getParPtr("delta_p0_ss_Run2")->err() << "\\\\" << "\n"
		<< "&$\\Delta p_{1}^{\\text{SS}}$  & " <<  mps->getParPtr("delta_p1_ss_Run2")->mean() << " $\\pm$ " << mps->getParPtr("delta_p1_ss_Run2")->err() << "\\\\" << "\n"
		<< "&$\\epsilon_{tag}^{\\text{SS}}$  & " <<  mps->getParPtr("tageff_ss_Run2")->mean() << " $\\pm$ " << mps->getParPtr("tageff_ss_Run2")->err() << "\\\\" << "\n"
		<< "&$\\Delta\\epsilon_{tag}^{\\text{SS}}$  & " <<  mps->getParPtr("tageff_asym_ss_Run2")->mean() << " $\\pm$ " << mps->getParPtr("tageff_asym_ss_Run2")->err() << "\\\\" << "\n"
	
		<< "&$A_{p}$ & " <<  mps->getParPtr("production_asym_Run2")->mean() << " $\\pm$ " << mps->getParPtr("production_asym_Run2")->err() << "\\\\" << "\n";
	
		resultsfile << "\\\\" << "\n" ;
		resultsfile << std::fixed << std::setprecision(4) 
		<< "&$\\Delta m_{s}$ & " 
		//<<   mps->getParPtr("dm")->mean() 
		<< " xx.xx " 
		<< " $\\pm$ " << mps->getParPtr("dm")->err() << "\\\\" << "\n";
	}

	else {
		resultsfile << std::fixed << std::setprecision(3) 
		<< "& $C$ & " 
	// 	<<   mps->getParPtr("C")->mean() 
		<< " xx.xx " 
		<< " $\\pm$ " << mps->getParPtr("C")->err() << "\\\\" << "\n";
		resultsfile << std::fixed << std::setprecision(3) 
		<< "&$D$ & " 
	// 	<<   mps->getParPtr("D")->mean() 
		<< " xx.xx " 
		<< " $\\pm$ " << mps->getParPtr("D")->err() << "\\\\" << "\n";
		resultsfile << std::fixed << std::setprecision(3) 
		<< "&$\\bar D$ & " 
	// 	<<   mps->getParPtr("D_bar")->mean() 
		<< " xx.xx " 
		<< " $\\pm$ " << mps->getParPtr("D_bar")->err() << "\\\\" << "\n";
		resultsfile << std::fixed << std::setprecision(3) 
		<< "& $S$ & " 
	// 	<<   mps->getParPtr("S")->mean() 
		<< " xx.xx " 
		<< " $\\pm$ " << mps->getParPtr("S")->err() << "\\\\" << "\n";
		resultsfile << std::fixed << std::setprecision(3) 
		<< "& $\\bar S$ & " 
	// 	<<   mps->getParPtr("S_bar")->mean() 
		<< " xx.xx " 
		<< " $\\pm$ " << mps->getParPtr("S_bar")->err() << "\\\\" << "\n";
	}

	resultsfile << "\\hline" << "\n";
	resultsfile << "\\hline" << "\n";
	resultsfile << "\\end{tabular}" << "\n";
	resultsfile << "\\label{table:timeFit_" << (string) channel << "}" << "\n";
	resultsfile << "\\end{table}";	
	resultsfile.close();

	for(int p = 0 ; p < prefix.size(); p++){

		vector<string> cov_params;
//     		if(!mps->getParPtr("Gamma")->iFixInit())cov_params.push_back("Gamma");
//     		if(!mps->getParPtr("dGamma")->iFixInit())cov_params.push_back("dGamma");
//     		if(!mps->getParPtr("dm")->iFixInit())cov_params.push_back("dm");

//     		if(!mps->getParPtr("offset_sigma_dt"+prefix[p])->iFixInit())cov_params.push_back("offset_sigma_dt"+prefix[p]);
//     		if(!mps->getParPtr("scale_sigma_dt"+prefix[p])->iFixInit())cov_params.push_back("scale_sigma_dt"+prefix[p]);

		if(!mps->getParPtr("p0_os"+prefix[p])->iFixInit())cov_params.push_back("p0_os"+prefix[p]);
    		if(!mps->getParPtr("p1_os"+prefix[p])->iFixInit())cov_params.push_back("p1_os"+prefix[p]);
    		if(!mps->getParPtr("delta_p0_os"+prefix[p])->iFixInit())cov_params.push_back("delta_p0_os"+prefix[p]);
    		if(!mps->getParPtr("delta_p1_os"+prefix[p])->iFixInit())cov_params.push_back("delta_p1_os"+prefix[p]);
//     		if(!mps->getParPtr("avg_eta_os"+prefix[p])->iFixInit())cov_params.push_back("avg_eta_os"+prefix[p]);
    		if(!mps->getParPtr("tageff_os"+prefix[p])->iFixInit())cov_params.push_back("tageff_os"+prefix[p]);
    		//if(!mps->getParPtr("tageff_asym_os"+prefix[p])->iFixInit())cov_params.push_back("tageff_asym_os"+prefix[p]);

		if(!mps->getParPtr("p0_ss"+prefix[p])->iFixInit())cov_params.push_back("p0_ss"+prefix[p]);
    		if(!mps->getParPtr("p1_ss"+prefix[p])->iFixInit())cov_params.push_back("p1_ss"+prefix[p]);
    		if(!mps->getParPtr("delta_p0_ss"+prefix[p])->iFixInit())cov_params.push_back("delta_p0_ss"+prefix[p]);
    		if(!mps->getParPtr("delta_p1_ss"+prefix[p])->iFixInit())cov_params.push_back("delta_p1_ss"+prefix[p]);
//     		if(!mps->getParPtr("avg_eta_ss"+prefix[p])->iFixInit())cov_params.push_back("avg_eta_ss"+prefix[p]);
    		if(!mps->getParPtr("tageff_ss"+prefix[p])->iFixInit())cov_params.push_back("tageff_ss"+prefix[p]);
    		//if(!mps->getParPtr("tageff_asym_ss"+prefix[p])->iFixInit())cov_params.push_back("tageff_asym_ss"+prefix[p]);

//     		if(!mps->getParPtr("production_asym"+prefix[p])->iFixInit())cov_params.push_back("production_asym"+prefix[p]);
//     		if(!mps->getParPtr("detection_asym"+prefix[p])->iFixInit())cov_params.push_back("detection_asym"+prefix[p]);

		vector<int> cov_params_id;
		RooArgList xvec, mu;
		

		for(int i = 0; i < cov_params.size(); i++){
			cov_params_id.push_back(mps->findParPtr(cov_params[i]));
			double mean = mps->getParPtr(cov_params[i])->mean();	
			double error = mps->getParPtr(cov_params[i])->err();	
		
			RooRealVar* x = new RooRealVar(("x_"+cov_params[i]).c_str(), ("x_"+cov_params[i]).c_str(),mean-10.*error,mean+10.*error);
			xvec.add(*x);
			mu.add(RooRealConstant::value(mean));
		}
		
		if(cov_params.size()>0){

			TMatrixTSym<double> cov(cov_params.size());	
			for(int i = 0; i < cov_params_id.size(); i++)for(int j = 0; j < cov_params_id.size(); j++) cov[i][j] = cov_full[cov_params_id[i]][cov_params_id[j]]; 
				
			cov.Print();
			xvec.Print();
			mu.Print();
			RooMultiVarGaussian gauss_cov("gauss_cov","gauss_cov",xvec, mu, cov);
	
			const int N_toys_cov = 1000; 
			RooDataSet* data_cov = gauss_cov.generate(xvec, N_toys_cov);
	
			double N_tot = 0;
			vector<double> v_N_OS(N_toys_cov,0.);
			vector<double> v_N_SS(N_toys_cov,0.);
			vector<double> v_N_OS_SS(N_toys_cov,0.);
			
			vector<double> v_w_OS(N_toys_cov,0.);
			vector<double> v_w_SS(N_toys_cov,0.);
			vector<double> v_w_OS_SS(N_toys_cov,0.);
				
			vector<double> v_D_OS(N_toys_cov,0.);
			vector<double> v_D_SS(N_toys_cov,0.);
			vector<double> v_D_OS_SS(N_toys_cov,0.);
			vector<double> v_D_comb(N_toys_cov,0.);
			
			for (unsigned int i=0; i<eventList.size(); i++) {
			
				int f_evt = eventList[i].getValueFromVector(2);
				int q1 = eventList[i].getValueFromVector(3);
				int q2 = eventList[i].getValueFromVector(5);   
				int q_eff = 0;
				double w_eff = 0.5;
				int run_evt = eventList[i].getValueFromVector(7);   
				int trigger_evt = eventList[i].getValueFromVector(8);   
	
				if(A_is_in_B("Run1",prefix[p]) && (run_evt == 2 || run_evt == 3)) continue;
				else if(A_is_in_B("Run2",prefix[p]) && run_evt == 1) continue;
	
				N_tot += eventList[i].getWeight();
	
				for(int j = 0 ; j < N_toys_cov; j++){
					RooArgSet* xvec_cov= (RooArgSet*)data_cov->get(j);
	
					double x_avg_eta_ss = xvec_cov->find(("x_avg_eta_ss"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_avg_eta_ss"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("avg_eta_ss"+prefix[p])->mean(); 
					double x_p0_ss = xvec_cov->find(("x_p0_ss"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_p0_ss"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("p0_ss"+prefix[p])->mean();
					double x_p1_ss = xvec_cov->find(("x_p1_ss"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_p1_ss"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("p1_ss"+prefix[p])->mean(); 
					double x_delta_p0_ss = xvec_cov->find(("x_delta_p0_ss"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_delta_p0_ss"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("delta_p0_ss"+prefix[p])->mean();
					double x_delta_p1_ss = xvec_cov->find(("x_delta_p1_ss"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_delta_p1_ss"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("delta_p1_ss"+prefix[p])->mean(); 
	
					double x_avg_eta_os = xvec_cov->find(("x_avg_eta_os"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_avg_eta_os"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("avg_eta_os"+prefix[p])->mean(); 
					double x_p0_os = xvec_cov->find(("x_p0_os"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_p0_os"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("p0_os"+prefix[p])->mean();
					double x_p1_os = xvec_cov->find(("x_p1_os"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_p1_os"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("p1_os"+prefix[p])->mean(); 
					double x_delta_p0_os = xvec_cov->find(("x_delta_p0_os"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_delta_p0_os"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("delta_p0_os"+prefix[p])->mean();
					double x_delta_p1_os = xvec_cov->find(("x_delta_p1_os"+prefix[p]).c_str()) ? ((RooRealVar*)xvec_cov->find(("x_delta_p1_os"+prefix[p]).c_str()))->getVal() :  mps->getParPtr("delta_p1_os"+prefix[p])->mean(); 
	
					std::pair<double, double> calibrated_mistag_os;
					std::pair<double, double> calibrated_mistag_ss;
		
					calibrated_mistag_ss = t_pdf.getCalibratedMistag_SS(eventList[i],x_avg_eta_ss,x_p0_ss,x_p1_ss,x_delta_p0_ss,x_delta_p1_ss);
					calibrated_mistag_os = t_pdf.getCalibratedMistag_OS(eventList[i],x_avg_eta_os,x_p0_os,x_p1_os,x_delta_p0_os,x_delta_p1_os);    
		
					double p = ( (1.-q1)/2. + q1 * (1.- calibrated_mistag_os.first )) * ( (1.-q2)/2. + q2 * (1.- calibrated_mistag_ss.first ));
					double p_bar = ( (1.+q1)/2. - q1 * (1.- calibrated_mistag_os.second )) * ( (1.+q2)/2. - q2 * (1.- calibrated_mistag_ss.second ));
					
					if( p/(p+p_bar) > 0.5 ){ 
						q_eff = 1;
						w_eff = 1-p/(p+p_bar);
					}
					else if( p/(p+p_bar) < 0.5 ){
						q_eff = -1;
						w_eff = p/(p+p_bar);
					}
				
					if(q1 != 0 && q2 != 0){
					if(q_eff != 0){
						v_N_OS_SS[j] += eventList[i].getWeight();
						v_w_OS_SS[j] += w_eff * eventList[i].getWeight();
						v_D_OS_SS[j] += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
					}
					}
					else if( q1 != 0){
						//q_eff = q1;  flip tag ???
						v_N_OS[j] += eventList[i].getWeight();
						v_w_OS[j] += w_eff * eventList[i].getWeight(); 
						v_D_OS[j] += pow(1.-2.*w_eff,2)* eventList[i].getWeight();
					}
					else if( q2 != 0){
						//q_eff = q2;
						v_N_SS[j] += eventList[i].getWeight();
						v_w_SS[j] += w_eff * eventList[i].getWeight(); 
						v_D_SS[j] += pow(1.-2.*w_eff,2)* eventList[i].getWeight(); 
					} 
				}
			}
	
			double eff_OS = 0;
			double eff_OS_err = 0;
			double w_OS = 0;
			double w_OS_err = 0;
			double e_OS = 0;
			double e_OS_err = 0;
	
			double eff_SS = 0;
			double eff_SS_err = 0;
			double w_SS = 0;
			double w_SS_err = 0;
			double e_SS = 0;
			double e_SS_err = 0;
	
			double eff_OS_SS = 0;
			double eff_OS_SS_err = 0;
			double w_OS_SS = 0;
			double w_OS_SS_err = 0;
			double e_OS_SS = 0;
			double e_OS_SS_err = 0;
		
			for(int j = 0 ; j < N_toys_cov; j++){
				eff_OS += (v_N_OS[j])/N_tot/N_toys_cov;
				w_OS += (v_w_OS[j])/(v_N_OS[j])/N_toys_cov;
				e_OS += (v_D_OS[j])/N_tot/N_toys_cov;
	
				eff_SS += (v_N_SS[j])/N_tot/N_toys_cov;
				w_SS += (v_w_SS[j])/(v_N_SS[j])/N_toys_cov;
				e_SS += (v_D_SS[j])/N_tot/N_toys_cov;
	
				eff_OS_SS += (v_N_OS_SS[j])/N_tot/N_toys_cov;
				w_OS_SS += (v_w_OS_SS[j])/(v_N_OS_SS[j])/N_toys_cov;
				e_OS_SS += (v_D_OS_SS[j])/N_tot/N_toys_cov;	
			}
			for(int j = 0 ; j < N_toys_cov; j++){
				eff_OS_err += pow(((v_N_OS[j])/N_tot-eff_OS),2)/(N_toys_cov-1.);
				w_OS_err += pow(((v_w_OS[j])/(v_N_OS[j]) - w_OS),2)/(N_toys_cov-1.);
				e_OS_err += pow((v_D_OS[j]/N_tot-e_OS),2)/(N_toys_cov-1.);
	
				eff_SS_err += pow(((v_N_SS[j])/N_tot-eff_SS),2)/(N_toys_cov-1.);
				w_SS_err += pow(((v_w_SS[j])/(v_N_SS[j]) - w_SS),2)/(N_toys_cov-1.);
				e_SS_err += pow((v_D_SS[j]/N_tot-e_SS),2)/(N_toys_cov-1.);
	
				eff_OS_SS_err += pow(((v_N_OS_SS[j])/N_tot-eff_OS_SS),2)/(N_toys_cov-1.);
				w_OS_SS_err += pow(((v_w_OS_SS[j])/(v_N_OS_SS[j]) - w_OS_SS),2)/(N_toys_cov-1.);
				e_OS_SS_err += pow((v_D_OS_SS[j]/N_tot-e_OS_SS),2)/(N_toys_cov-1.);
			}
			
			double rel_eff_err_OS = mps->getParPtr("tageff_os"+prefix[p])->iFixInit() ? 0. : mps->getParPtr("tageff_os"+prefix[p])->err() / mps->getParPtr("tageff_os"+prefix[p])->mean(); 
	
			eff_OS_err += pow( rel_eff_err_OS * eff_OS,2);
			w_OS_err += pow( rel_eff_err_OS * w_OS,2);
			e_OS_err += pow( rel_eff_err_OS * e_OS,2);
	
			double rel_eff_err_SS = mps->getParPtr("tageff_ss"+prefix[p])->iFixInit() ? 0. : mps->getParPtr("tageff_ss"+prefix[p])->err() / mps->getParPtr("tageff_ss"+prefix[p])->mean(); 
	
			eff_SS_err += pow( rel_eff_err_SS * eff_SS,2);
			w_SS_err += pow( rel_eff_err_SS * w_SS,2);
			e_SS_err += pow( rel_eff_err_SS * e_SS,2);
	
			eff_OS_SS_err += pow( rel_eff_err_OS * eff_OS_SS,2) + pow( rel_eff_err_SS * eff_OS_SS,2);
			w_OS_SS_err += pow( rel_eff_err_OS * w_OS_SS,2) + pow( rel_eff_err_SS * w_OS_SS,2);
			e_OS_SS_err += pow( rel_eff_err_OS * e_OS_SS,2) + pow( rel_eff_err_SS * e_OS_SS,2);
	
			double eff_tot = eff_OS + eff_SS + eff_OS_SS;
			double eff_tot_err = eff_OS_err + eff_SS_err + eff_OS_SS_err;
			double w_tot = (eff_OS * w_OS + eff_SS * w_SS + eff_OS_SS * w_OS_SS)/eff_tot;
			double w_tot_err = (eff_OS * w_OS_err + eff_SS * w_SS_err + eff_OS_SS * w_OS_SS_err)/eff_tot;
			double e_tot = e_OS + e_SS + e_OS_SS;
			double e_tot_err = e_OS_err + e_SS_err + e_OS_SS_err;
	
			ofstream datafile;
			datafile.open(("../../../../../TD-AnaNote/latex/tables/Tagging/"+(string)OutputDir + "tagPower"+ prefix[p] + ".tex").c_str(),std::ofstream::trunc);
			//datafile << "\\begin{table}[h]" << "\n";
			//datafile << "\\centering" << "\n";
			//datafile << "\\caption{The flavour tagging performances for only OS tagged, only SS tagged and both OS and SS tagged events";
			//if(A_is_in_B("Run1", prefix[p])) datafile << " for Run-I data"; 
			//else if(A_is_in_B("Run2", prefix[p])) datafile << " for Run-II data"; 
			//datafile << ".}\n";;
			datafile << "\\begin{tabular}{c c c c}" << "\n";
			datafile << "\\hline" << "\n";
			datafile << "\\hline" << "\n";
			if((string)channel == "norm")datafile << "$ B_s \\to D_s \\pi \\pi \\pi$";
			else if((string)channel == "signal")datafile << "$ B_s \\to D_s K \\pi \\pi$";
			datafile << " & $\\epsilon_{tag} [\\%]$ & $\\langle \\omega \\rangle [\\%] $ & $\\epsilon_{eff} [\\%]$ \\\\" << "\n";
			datafile << "\\hline" << "\n";
			datafile << std::fixed << std::setprecision(2) << "Only OS & " 
			<< eff_OS * 100. << " $\\pm$ " << sqrt(eff_OS_err) * 100. << " & " 
			<< w_OS * 100. << " $\\pm$ " << sqrt(w_OS_err) * 100. << " & "
			<< e_OS * 100.<< " $\\pm$ " << sqrt(e_OS_err) * 100. << "\\\\" << "\n";
	
			datafile << std::fixed << std::setprecision(2) << "Only SS & " 
			<< eff_SS * 100.<< " $\\pm$ " << sqrt(eff_SS_err) * 100. << " & " 
			<< w_SS * 100.<< " $\\pm$ " << sqrt(w_SS_err) * 100. << " & "
			<< e_SS * 100.<< " $\\pm$ " << sqrt(e_SS_err) * 100. << "\\\\" << "\n";
	
			datafile << std::fixed << std::setprecision(2) << "Both OS-SS & " 
			<< eff_OS_SS * 100.<< " $\\pm$ " << sqrt(eff_OS_SS_err) * 100. << " & " 
			<< w_OS_SS * 100.<< " $\\pm$ " << sqrt(w_OS_SS_err) * 100. << " & "
			<< e_OS_SS * 100.<< " $\\pm$ " << sqrt(e_OS_SS_err) * 100. << "\\\\" << "\n";
	
			datafile << "\\hline" << "\n" << std::fixed << std::setprecision(2) << "Combined & " 
			<< eff_tot * 100.<< " $\\pm$ " << sqrt(eff_tot_err) * 100. << " & " 
			<< w_tot * 100.<< " $\\pm$ " << sqrt(w_tot_err) * 100. << " & "
			<< e_tot * 100.<< " $\\pm$ " << sqrt(e_tot_err) * 100. << "\\\\" << "\n";
	
			datafile << "\\hline" << "\n";
			datafile << "\\hline" << "\n";
			datafile << "\\end{tabular}" << "\n";
			//datafile << "\\label{table:tagging" << prefix[p] << "}" << "\n";
			//datafile << "\\end{table}";	
		}
	}
    }
     */
//     nBinst *= 2;
    nBinsAsym *= scale_asym;

    TH1D* h_t_fit = new TH1D("h_t_fit",";t",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_p_fit = new TH1D("h_t_p",";t (ps);Events (a.u.) ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_t_m_fit = new TH1D("h_t_m",";t (ps);Events (a.u.) ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_t_N_fit = new TH1D("h_t_N_fit",";t (ps);Events (a.u.) ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_t_Nbar_fit = new TH1D("h_t_Nbar_fit",";t (ps);Events (a.u.) ",nBinsAsym,min_TAU,max_TAU);

    TH1D* h_t_mixed_fit = new TH1D("h_t_mixed_fit",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed_fit = new TH1D("h_t_unmixed_fit",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_untagegged_fit = new TH1D("h_t_untagegged_fit",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_t_OS_fit = new TH1D("h_t_OS_fit",";t",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_mixed_OS_fit = new TH1D("h_t_mixed_OS_fit",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed_OS_fit = new TH1D("h_t_unmixed_OS_fit",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_SS_fit = new TH1D("h_t_SS_fit",";t",nBinst,min_TAU,max_TAU);    
    TH1D* h_t_mixed_SS_fit = new TH1D("h_t_mixed_SS_fit",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_unmixed_SS_fit = new TH1D("h_t_unmixed_SS_fit",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_t_fit_mp = new TH1D("h_t_fit_mp",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_0p = new TH1D("h_t_fit_0p",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_pp = new TH1D("h_t_fit_pp",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_mm = new TH1D("h_t_fit_mm",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_0m = new TH1D("h_t_fit_0m",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_t_fit_pm = new TH1D("h_t_fit_pm",";t (ps);Events (a.u.) ",nBinst,min_TAU,max_TAU_ForMixingPlot);

    TH1D* h_dt_fit = new TH1D("h_dt_fit",";#sigma_{t} (ps);Events (a.u.) ",nBinst,0,0.15);
    TH1D* h_eta_OS_fit = new TH1D("h_eta_OS_fit",";#eta_{OS};Events (a.u.) ",nBinst,0,0.5);
    TH1D* h_eta_SS_fit = new TH1D("h_eta_SS_fit",";#eta_{SS};Events (a.u.) ",nBinst,0,0.5);
    TH1D* h_q_OS_fit = new TH1D("h_q_OS_fit",";q_{OS};Events (a.u.) ",3,-1.5,1.5);
    TH1D* h_q_SS_fit = new TH1D("h_q_SS_fit",";q_{SS};Events (a.u.) ",3,-1.5,1.5);
    TH1D* h_f_fit = new TH1D("h_f_fit",";q_{f};Events (a.u.) ",2,-2,2);

    TH1D* h_N_mixed_fit = new TH1D("h_N_mixed_fit",";t modulo (2#pi/#Deltam_{s}) (ps); A_{mix} ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_N_unmixed_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_unmixed_fit");

    TH1D* h_N_mixed_OS_fit = new TH1D("h_N_mixed_OS_fit",";t modulo (2#pi/#Deltam_{s}) (ps); A_{mix} ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_N_unmixed_OS_fit = (TH1D*) h_N_mixed_OS_fit->Clone("h_N_unmixed_OS_fit");
    TH1D* h_N_mixed_SS_fit = new TH1D("h_N_mixed_SS_fit",";t modulo (2#pi/#Deltam_{s}) (ps); A_{mix} ",nBinsAsym,min_TAU,max_TAU);
    TH1D* h_N_unmixed_SS_fit = (TH1D*) h_N_mixed_SS_fit->Clone("h_N_unmixed_SS_fit");

    TH1D* h_N_mixed_p_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_mixed_p_fit");
    TH1D* h_N_unmixed_p_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_unmixed_p_fit");
    TH1D* h_N_mixed_m_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_mixed_m_fit");
    TH1D* h_N_unmixed_m_fit = (TH1D*) h_N_mixed_fit->Clone("h_N_unmixed_m_fit");

    TH1D* h_N_mixed_p_fit_unfolded = new TH1D("h_N_mixed_p_fit_unfolded",";t/#tau;A_{CP}(t) ",nBinst,min_TAU,max_TAU_ForMixingPlot);
    TH1D* h_N_unmixed_p_fit_unfolded = (TH1D*) h_N_mixed_p_fit_unfolded->Clone("h_N_unmixed_p_fit");
    TH1D* h_N_mixed_m_fit_unfolded = (TH1D*) h_N_mixed_p_fit_unfolded->Clone("h_N_mixed_m_fit");
    TH1D* h_N_unmixed_m_fit_unfolded = (TH1D*) h_N_mixed_p_fit_unfolded->Clone("h_N_unmixed_m_fit");

    double N_OS_MC = 0;
    double N_SS_MC = 0;
    double N_OS_SS_MC = 0;
    double N_MC = 0;

    double w_OS_MC = 0;
    double w_SS_MC = 0;
    double w_OS_SS_MC = 0;
 	
    double D_OS_MC = 0;
    double D_SS_MC = 0;
    double D_OS_SS_MC = 0;
    double D_comb_MC = 0;

    double N_OS_all_MC = 0;
    double N_SS_all_MC = 0;
    double w_OS_all_MC = 0;
    double w_SS_all_MC = 0;
    double D_OS_all_MC = 0;
    double D_SS_all_MC = 0;
    
    double N_Run1_t0_MC = 0;
    double N_Run1_t1_MC = 0;
    double N_Run2_t0_MC = 0;
    double N_Run2_t1_MC = 0;

   if (doPlots){ 
	for(int n = 0; n < N_plot_it; n++){   /// Multiple iterations needed to release memory 
		int N_sample = 500000;
		DalitzEventList sampleEvents;
		if(doSimFit) {
			if(N_Run1>0)sampleEvents.Add(t_pdf_Run1.generateToys(N_sample * N_Run1/N,1,0));
			if(N_Run2>0)sampleEvents.Add(t_pdf_Run2.generateToys(N_sample *N_Run2/N,2,0));
		}
		else sampleEvents.Add(t_pdf_Run1.generateToys(N_sample));	

		for(int i = 0; i < sampleEvents.size(); i++){

			DalitzEvent evt = sampleEvents[i];
			double t_MC = evt.getValueFromVector(0) ;
			double dt_MC = evt.getValueFromVector(1) ;
			int f_MC = evt.getValueFromVector(2) ;
			int q_OS_MC = evt.getValueFromVector(3) ;
			double eta_OS_MC = evt.getValueFromVector(4) ;
			int q_SS_MC = evt.getValueFromVector(5);
			double eta_SS_MC = evt.getValueFromVector(6);
			int run_MC = evt.getValueFromVector(7) ;
			int trigger_MC = evt.getValueFromVector(8) ;

			double weight = 1; //t_pdf_Run1_t0.getVal(evt)/evt.getGeneratorPdfRelativeToPhaseSpace();
			N_MC += weight;
		
			h_t_fit->Fill(t_MC,weight);
			h_dt_fit->Fill(dt_MC,weight);
			if(q_OS_MC != 0)h_eta_OS_fit->Fill(eta_OS_MC,weight);
			if(q_SS_MC != 0)h_eta_SS_fit->Fill(eta_SS_MC,weight);
			
			int f_evt = f_MC;
			int q1 = q_OS_MC;
			int q2 = q_SS_MC;   
			int q_eff = 0;
			double w_eff = 0.5;
		
			h_q_OS_fit->Fill(q1,weight);
			h_q_SS_fit->Fill(q2,weight);
			h_f_fit->Fill(f_evt,weight);
			
			if( q1 != 0)h_t_OS_fit->Fill(t_MC,weight);
			if( q2 != 0)h_t_SS_fit->Fill(t_MC,weight);

			std::pair<double, double> calibrated_mistag_os;
			std::pair<double, double> calibrated_mistag_ss;
			if(doSimFit){
				if(run_MC==1){
					calibrated_mistag_os = t_pdf_Run1.getCalibratedMistag_OS(eta_OS_MC);
					calibrated_mistag_ss = t_pdf_Run1.getCalibratedMistag_SS(eta_SS_MC);
				}
				else{
					calibrated_mistag_os = t_pdf_Run2.getCalibratedMistag_OS(eta_OS_MC);
					calibrated_mistag_ss = t_pdf_Run2.getCalibratedMistag_SS(eta_SS_MC);                
				}
			}
			else{
				calibrated_mistag_os = t_pdf.getCalibratedMistag_OS(eta_OS_MC);
				calibrated_mistag_ss = t_pdf.getCalibratedMistag_SS(eta_SS_MC);        
			}
			
			double p = ( (1.-q1)/2. + q1 * (1.- calibrated_mistag_os.first )) * ( (1.-q2)/2. + q2 * (1.- calibrated_mistag_ss.first ));
			double p_bar = ( (1.+q1)/2. - q1 * (1.- calibrated_mistag_os.second )) * ( (1.+q2)/2. - q2 * (1.- calibrated_mistag_ss.second ));
			
			if( p/(p+p_bar) > 0.5 ){ 
				q_eff = 1;
				w_eff = 1-p/(p+p_bar);
			}
			else if( p/(p+p_bar) < 0.5 ){
				q_eff = -1;
				w_eff = p/(p+p_bar);
			}
			
			if(q1 != 0 && q2 != 0){
				if(q_eff != 0){
					N_OS_SS_MC += weight;
					w_OS_SS_MC += w_eff * weight;
					D_OS_SS_MC += pow(1.-2.*w_eff,2)* weight;
				}
			}
			else if( q1 != 0){
				//q_eff = q1;  flip tag ???
				N_OS_MC +=weight;
				w_OS_MC += w_eff * weight; 
				D_OS_MC += pow(1.-2.*w_eff,2)* weight;
			}
			else if( q2 != 0){
				//q_eff = q2;
				N_SS_MC += weight;
				w_SS_MC += w_eff * weight; 
				D_SS_MC += pow(1.-2.*w_eff,2)* weight; 
			} 
		
			D_comb_MC += pow(1.-2.*w_eff,2)* weight;
			double D_res = 1.;
			if(doSimFit){
				if(run_MC==1){
					D_res = exp(-pow(t_pdf_Run1.getCalibratedResolution(evt.getValueFromVector(1))*dm,2)/2.);
				}
				else if(run_MC==2){
					D_res = exp(-pow(t_pdf_Run2.getCalibratedResolution(evt.getValueFromVector(1))*dm,2)/2.);
				}
			}
			double D_tag = 0.;
			if(q_eff != 0) D_tag = (1.-2.*abs(w_eff));
			double D_tot = D_tag * D_res;
// 			double D_tot =  D_res;
			if(!dilutionWeight)D_tot = 1.;
// 			double D_tot = (1.-2.*abs(w_eff)) * exp(-pow(t_pdf.getCalibratedResolution(evt.getValueFromVector(1))*dm,2)/2.);

			if(q1 != 0) N_OS_all_MC += weight;
			if(q2 != 0)N_SS_all_MC += weight;
			
			if(q1>0){
					w_OS_all_MC +=  calibrated_mistag_os.first * weight;
					D_OS_all_MC +=  pow(1.-2.*calibrated_mistag_os.first,2)* weight;
			} 
			else if(q1<0){
					w_OS_all_MC +=  calibrated_mistag_os.second * weight;	
					D_OS_all_MC +=  pow(1.-2.*calibrated_mistag_os.second,2)* weight;
			}
		
			if(q2>0){
					w_SS_all_MC +=  calibrated_mistag_ss.first * weight;
					D_SS_all_MC +=  pow(1.-2.*calibrated_mistag_ss.first,2)* weight;
			} 
			else if(q2<0){
					w_SS_all_MC +=  calibrated_mistag_ss.second * weight;	
					D_SS_all_MC +=  pow(1.-2.*calibrated_mistag_ss.second,2)* weight;
			}
		
			if((string)channel=="signal"){

				if(f_evt == 1){ 
					h_t_p_fit->Fill(t_MC,weight);
				}
				else h_t_m_fit->Fill(t_MC,weight);
			
	   			if(q_eff == 1)h_t_N_fit->Fill(t_MC,weight*D_tot);
	    			else if(q_eff == -1)h_t_Nbar_fit->Fill(t_MC,weight*D_tot);

				if(q_eff==-1 && f_evt == 1){
					h_t_fit_mp->Fill(t_MC,weight*D_tot);
					if(w_eff<w_max){
						h_N_mixed_p_fit->Fill(t_MC,weight*D_tot);
						h_N_mixed_p_fit_unfolded->Fill(t_MC,weight*D_tot);
					}
				}
				else if(q_eff==0 && f_evt == 1)h_t_fit_0p->Fill(t_MC,weight*D_tot);
				else if(q_eff==1 && f_evt == 1){
					h_t_fit_pp->Fill(t_MC,weight*D_tot);
					if(w_eff<w_max){
						h_N_unmixed_p_fit->Fill(t_MC,weight*D_tot);
						h_N_unmixed_p_fit_unfolded->Fill(t_MC,weight*D_tot);
					}
				}
				else if(q_eff==-1 && f_evt == -1){
					h_t_fit_mm->Fill(t_MC,weight*D_tot);
					if(w_eff<w_max){
						h_N_unmixed_m_fit->Fill(t_MC,weight*D_tot);
						h_N_unmixed_m_fit_unfolded->Fill(t_MC,weight*D_tot);
					}	
				}
				else if(q_eff==0 && f_evt == -1)h_t_fit_0m->Fill(t_MC,weight*D_tot);
				else if(q_eff==1 && f_evt == -1){
					h_t_fit_pm->Fill(t_MC,weight*D_tot);
					if(w_eff<w_max){
						h_N_mixed_m_fit->Fill(t_MC,weight*D_tot);
						h_N_mixed_m_fit_unfolded->Fill(t_MC,weight*D_tot);
					}
				}
			}
//  			else {   
				if(q_eff == 0)h_t_untagegged_fit->Fill(t_MC,weight);
				else if(q_eff*f_evt > 0  ){
					h_t_mixed_fit->Fill(t_MC,weight*D_tot);
					//if(w_eff<w_max)
					h_N_mixed_fit->Fill(t_MC,weight*D_tot);

					if( q1 != 0)h_t_mixed_OS_fit->Fill(t_MC,weight*D_tot);
					if( q1 != 0)h_N_mixed_OS_fit->Fill(t_MC,weight*D_tot);

					if( q2 != 0)h_t_mixed_SS_fit->Fill(t_MC,weight*D_tot);
					if( q2 != 0)h_N_mixed_SS_fit->Fill(t_MC,weight*D_tot);
				}
				else{ 
					h_t_unmixed_fit->Fill(t_MC,weight*D_tot);
					//if(w_eff<w_max)
					h_N_unmixed_fit->Fill(t_MC,weight*D_tot);

					if( q1 != 0)h_t_unmixed_OS_fit->Fill(t_MC,weight*D_tot);
					if( q1 != 0)h_N_unmixed_OS_fit->Fill(t_MC,weight*D_tot);

					if( q2 != 0)h_t_unmixed_SS_fit->Fill(t_MC,weight*D_tot);
					if( q2 != 0)h_N_unmixed_SS_fit->Fill(t_MC,weight*D_tot);
				}
 			}
// 		}
	}
	
	cout << "N_MC = " << N_MC << endl << endl;
	cout << "N_OS = " << N_OS_all_MC << endl;
	cout << "N_SS = " << N_SS_all_MC << endl << endl;
	
	cout << "eff_OS = " <<(N_OS_all_MC)/N_MC << endl;
	cout << "eff_SS = " <<(N_SS_all_MC)/N_MC << endl;
		
	cout << "Tagging perfromance " << endl << endl;        
	cout << "Tagger | eff_tag | <w> | e_eff " <<  endl;
	
	cout << "OS  | " << (N_OS_MC+N_OS_SS_MC)/N_MC << " | " <<  (w_OS_all_MC)/(N_OS_MC+N_OS_SS_MC) << " | " << D_OS_all_MC/N_MC << endl;
	cout << "SS  | " << (N_SS_MC+N_OS_SS_MC)/N_MC << " | " <<  (w_SS_all_MC)/(N_SS_MC+N_OS_SS_MC) << " | " << D_SS_all_MC/N_MC << endl << endl;
	
	cout << "OS only  | " << N_OS_MC/N_MC << " | " <<  w_OS_MC/N_OS_MC << " | " << N_OS_MC/N_MC * D_OS_MC/N_OS_MC << endl;
	cout << "SS only  | " << N_SS_MC/N_MC << " | " <<  w_SS_MC/N_SS_MC << " | " << N_SS_MC/N_MC * D_SS_MC/N_SS_MC << endl;
	cout << "OS+SS    | " << N_OS_SS_MC/N_MC << " | " <<  w_OS_SS_MC/N_OS_SS_MC << " | " << N_OS_SS_MC/N_MC * D_OS_SS_MC/N_OS_SS_MC << endl;
	cout << "Combined | " << (N_OS_MC+N_SS_MC+N_OS_SS_MC)/N_MC << " | "<<  (w_OS_MC+w_SS_MC+w_OS_SS_MC)/(N_OS_MC+N_SS_MC+N_OS_SS_MC) << " | " << (N_OS_MC+N_SS_MC+N_OS_SS_MC)/N_MC * D_comb_MC/(N_OS_MC+N_SS_MC+N_OS_SS_MC) << endl << endl ;
	
	/// Plots
	TGraph* graph = new TGraph(2);
	graph->SetPoint(1,min_TAU,0);
	graph->SetPoint(2,max_TAU,0);
	graph->SetLineStyle(kDashed);

	h_t->SetMinimum(0.1);    
	h_t->SetLineColor(kBlack);
	h_t->DrawNormalized("e",1);
		
	h_t_fit->SetLineColor(kBlue);
	h_t_fit->SetLineWidth(3);
	h_t_fit->SetMarkerColor(kBlue); 
	h_t_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_t.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_t.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"h_t_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_t_log.pdf").c_str());
	gPad->SetLogy(0);
	
	h_t_OS->SetMinimum(0.1);    
	h_t_OS->SetLineColor(kBlack);
	h_t_OS->DrawNormalized("e",1);
	h_t_OS_fit->SetLineColor(kBlue);
	h_t_OS_fit->SetLineWidth(3);
	h_t_OS_fit->SetMarkerColor(kBlue); 
	h_t_OS_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_t_OS.eps").c_str());

	h_t_SS->SetMinimum(0.1);    
	h_t_SS->SetLineColor(kBlack);
	h_t_SS->DrawNormalized("e",1);
	h_t_SS_fit->SetLineColor(kBlue);
	h_t_SS_fit->SetLineWidth(3);
	h_t_SS_fit->SetMarkerColor(kBlue); 
	h_t_SS_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_t_SS.eps").c_str());

	h_t_untagegged->SetMinimum(0.1);    
	h_t_untagegged->SetLineColor(kBlack);
	h_t_untagegged->DrawNormalized("e",1);
	h_t_untagegged_fit->SetLineColor(kBlue);
	h_t_untagegged_fit->SetLineWidth(3);
	h_t_untagegged_fit->SetMarkerColor(kBlue); 
	h_t_untagegged_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_t_untagged.eps").c_str());

	h_dt->SetMinimum(0);        
	h_dt->SetLineColor(kBlack);
	h_dt->DrawNormalized("e1",1);
	h_dt_fit->SetLineColor(kBlue);
	h_dt_fit->SetLineWidth(3);
	h_dt_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_dt.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_dt.pdf").c_str());
	
	h_eta_OS->SetMinimum(0);        
	h_eta_OS->SetLineColor(kBlack);
	h_eta_OS->DrawNormalized("e1",1);
	h_eta_OS_fit->SetLineColor(kBlue);
	h_eta_OS_fit->SetLineWidth(3);
	h_eta_OS_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_eta_OS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_eta_OS.pdf").c_str());
	
	h_eta_SS->SetMinimum(0);        
	h_eta_SS->SetLineColor(kBlack);
	h_eta_SS->DrawNormalized("e1",1);
	h_eta_SS_fit->SetLineColor(kBlue);
	h_eta_SS_fit->SetLineWidth(3);
	h_eta_SS_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_eta_SS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_eta_SS.pdf").c_str());
	
	h_q_OS->SetMinimum(0);        
	h_q_OS->SetLineColor(kBlack);
	h_q_OS->DrawNormalized("e1",1);
	h_q_OS_fit->SetLineColor(kBlue);
	h_q_OS_fit->SetLineWidth(3);
	h_q_OS_fit->DrawNormalized("histsame",1);
	c->Print(((string)OutputDir+"h_q_OS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_q_OS.pdf").c_str());
	
	h_q_SS->SetMinimum(0);        
	h_q_SS->SetLineColor(kBlack);
	h_q_SS->DrawNormalized("e1",1);
	h_q_SS_fit->SetLineColor(kBlue);
	h_q_SS_fit->SetLineWidth(3);
	h_q_SS_fit->DrawNormalized("histsame",1);
	c->Print(((string)OutputDir+"h_q_SS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_q_SS.pdf").c_str());
	
	h_f->SetMinimum(0);        
	h_f->SetLineColor(kBlack);
	h_f->DrawNormalized("e1",1);
	h_f_fit->SetLineColor(kBlue);
	h_f_fit->SetLineWidth(3);
	h_f_fit->DrawNormalized("histsame",1);
	c->Print(((string)OutputDir+"h_f.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_f.pdf").c_str());
	
	if((string)channel=="norm"){
	
		h_t_mixed->SetMarkerColor(kRed); 
		h_t_mixed->SetLineColor(kRed);
		h_t_mixed->DrawNormalized("e",1);
		
		h_t_unmixed->SetMarkerColor(kBlue); 
		h_t_unmixed->SetLineColor(kBlue);
		h_t_unmixed->DrawNormalized("esame",1);
		
		h_t_mixed_fit->SetMarkerColor(kRed); 
		h_t_mixed_fit->SetLineColor(kRed);
		h_t_mixed_fit->DrawNormalized("histcsame",1);
		
		h_t_unmixed_fit->SetMarkerColor(kBlue); 
		h_t_unmixed_fit->SetLineColor(kBlue);
		h_t_unmixed_fit->DrawNormalized("histcsame",1);
		
		h_t_untagegged->SetMarkerColor(kGreen); 
		h_t_untagegged->SetLineColor(kGreen);
		//h_t_untagegged->DrawNormalized("esame",1);
		
		h_t_untagegged_fit->SetMarkerColor(kGreen); 
		h_t_untagegged_fit->SetLineColor(kGreen);
		//h_t_untagegged_fit->DrawNormalized("histcsame",1);

		TLegend leg(0.6,0.6,0.9,0.9,"");
		leg.SetLineStyle(0);
		leg.SetLineColor(0);
		leg.SetFillColor(0);
		leg.SetTextFont(22);
		leg.SetTextColor(1);
		leg.SetTextSize(0.05);
		leg.SetTextAlign(12);

		//leg.AddEntry((TObject*)0,"#font[22]{LHCb}","");
		TLegendEntry* le = leg.AddEntry(h_t_mixed,"B^{0}_{s}(t)#rightarrow f","l");
		le->SetTextColor(kRed);
		le = leg.AddEntry(h_t_unmixed,"#bar{B^{0}_{s}}(t)#rightarrow f","l");
	        le->SetTextColor(kBlue);    
		leg.Draw();
	
		c->Print(((string)OutputDir+"h_t_mixed.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_t_mixed.pdf").c_str());
	

		TLegend leg2(0.3,0.75,0.6,0.92,"");
		leg2.SetLineStyle(0);
		leg2.SetLineColor(0);
		leg2.SetFillColor(0);
		leg2.SetTextFont(22);
		leg2.SetTextColor(1);
		leg2.SetTextSize(0.065);
		leg2.SetTextAlign(12);
		le = leg2.AddEntry((TObject*)0,"0.30 < #eta_{SS} < 0.40","");
		le->SetTextColor(kBlack);    				

		TH1D* h_asym = (TH1D*) h_N_mixed->GetAsymmetry(h_N_unmixed);	
		double max_asym = max(h_asym->GetMaximum(),fabs(h_asym->GetMinimum())) *1.25;
		h_asym->SetMaximum(max_asym);
		h_asym->SetMinimum(-max_asym);

 		h_asym->SetMaximum(1);
		h_asym->SetMinimum(-1);

		TH1D* h_asym_fit = (TH1D*) h_N_mixed_fit->GetAsymmetry(h_N_unmixed_fit);	
		h_asym_fit->SetLineColor(kRed);
		h_asym_fit->SetLineWidth(3);
		h_asym->Draw("e");
		h_asym_fit->Draw("histcsame");
// 		leg2.Draw();
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym.pdf").c_str());


		h_t_mixed_OS->SetMarkerColor(kRed); 
		h_t_mixed_OS->SetLineColor(kRed);
		h_t_mixed_OS->DrawNormalized("e",1);
		
		h_t_unmixed_OS->SetMarkerColor(kBlue); 
		h_t_unmixed_OS->SetLineColor(kBlue);
		h_t_unmixed_OS->DrawNormalized("esame",1);
		
		h_t_mixed_OS_fit->SetMarkerColor(kRed); 
		h_t_mixed_OS_fit->SetLineColor(kRed);
		h_t_mixed_OS_fit->DrawNormalized("histcsame",1);
		
		h_t_unmixed_OS_fit->SetMarkerColor(kBlue); 
		h_t_unmixed_OS_fit->SetLineColor(kBlue);
		h_t_unmixed_OS_fit->DrawNormalized("histcsame",1);
		c->Print(((string)OutputDir+"h_t_mixed_OS.eps").c_str());

		h_t_mixed_SS->SetMarkerColor(kRed); 
		h_t_mixed_SS->SetLineColor(kRed);
		h_t_mixed_SS->DrawNormalized("e",1);
		
		h_t_unmixed_SS->SetMarkerColor(kBlue); 
		h_t_unmixed_SS->SetLineColor(kBlue);
		h_t_unmixed_SS->DrawNormalized("esame",1);
		
		h_t_mixed_SS_fit->SetMarkerColor(kRed); 
		h_t_mixed_SS_fit->SetLineColor(kRed);
		h_t_mixed_SS_fit->DrawNormalized("histcsame",1);
		
		h_t_unmixed_SS_fit->SetMarkerColor(kBlue); 
		h_t_unmixed_SS_fit->SetLineColor(kBlue);
		h_t_unmixed_SS_fit->DrawNormalized("histcsame",1);
		c->Print(((string)OutputDir+"h_t_mixed_SS.eps").c_str());
	}
	
	else{
	
		h_t_mp->Scale(1./h_t_mp->Integral());
		double maxY= h_t_mp->GetMaximum()*1.3;        
		h_t_mp->SetMinimum(0.);  
		h_t_mp->SetMaximum(maxY);        
		h_t_mp->SetMarkerColor(kBlue); 
		h_t_mp->SetLineColor(kBlue);
		h_t_mp->DrawNormalized("e",1);
	
		h_t_pp->SetMarkerColor(kRed); 
		h_t_pp->SetLineColor(kRed);
		h_t_pp->DrawNormalized("esame",1);
		
		h_t_fit_mp->SetMarkerColor(kBlue); 
		h_t_fit_mp->SetLineColor(kBlue);
		h_t_fit_mp->DrawNormalized("histcsame",1);
		
		h_t_fit_pp->SetMarkerColor(kRed); 
		h_t_fit_pp->SetLineColor(kRed);
		h_t_fit_pp->DrawNormalized("histcsame",1);
		
		c->Print(((string)OutputDir+"h_t_mixed_p.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_t_mixed_p.pdf").c_str());
	
		h_t_pm->Scale(1./h_t_pm->Integral());
		h_t_pm->SetMinimum(0.);     
		h_t_pm->SetMaximum(maxY);        
		h_t_pm->SetMarkerColor(kBlue); 
		h_t_pm->SetLineColor(kBlue);
		h_t_pm->DrawNormalized("e",1);
		
		h_t_mm->SetMarkerColor(kRed); 
		h_t_mm->SetLineColor(kRed);
		h_t_mm->DrawNormalized("esame",1);
		
		h_t_fit_pm->SetMarkerColor(kBlue); 
		h_t_fit_pm->SetLineColor(kBlue);
		h_t_fit_pm->DrawNormalized("histcsame",1);
		
		h_t_fit_mm->SetMarkerColor(kRed); 
		h_t_fit_mm->SetLineColor(kRed);
		h_t_fit_mm->DrawNormalized("histcsame",1);
		
		c->Print(((string)OutputDir+"h_t_mixed_m.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_t_mixed_m.pdf").c_str());
	
		cout << h_N_unmixed_p->GetEntries() << endl;
		cout << h_N_mixed_p->GetEntries() << endl;
	
		TH1D* h_asym_p = (TH1D*) h_N_unmixed_p->GetAsymmetry(h_N_mixed_p);	
		//h_asym_p->SetMinimum(-20);
		//h_asym_p->SetMaximum(20);
		TH1D* h_asym_p_fit = (TH1D*) h_N_unmixed_p_fit->GetAsymmetry(h_N_mixed_p_fit);	
		h_asym_p_fit->SetLineColor(kRed);
		h_asym_p->Draw("e");
		h_asym_p_fit->Draw("histcsame");
		c->Print(((string)OutputDir+"h_asym_p.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym_p.pdf").c_str());

		TH1D* h_asym_m = (TH1D*) h_N_unmixed_m->GetAsymmetry(h_N_mixed_m);	
		//h_asym_m->SetMinimum(-20);
		//h_asym_m->SetMaximum(20);
		TH1D* h_asym_m_fit = (TH1D*) h_N_unmixed_m_fit->GetAsymmetry(h_N_mixed_m_fit);	
		h_asym_m_fit->SetLineColor(kRed);
		h_asym_m->Draw("e");
		h_asym_m_fit->Draw("histcsame");
		c->Print(((string)OutputDir+"h_asym_m.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym_m.pdf").c_str());
	

		double max_asym = max(max(h_asym_p->GetMaximum(),h_asym_m->GetMaximum()),fabs(min(h_asym_p->GetMinimum(),h_asym_m->GetMinimum()))) *1.25;
		h_asym_p->SetMaximum(1);
		h_asym_p->SetMinimum(-1);	
// 		h_asym_p->SetMarkerColor(kRed);
// 		h_asym_p->SetLineColor(kRed);
		h_asym_p->Draw("e");
		h_asym_m->SetLineColor(kGreen+3);
		h_asym_m->SetMarkerColor(kGreen+3);
		h_asym_m->Draw("esame");
		h_asym_p_fit->SetLineColor(kBlack);
		h_asym_p_fit->SetLineWidth(5);
		h_asym_m_fit->SetLineWidth(5);
		h_asym_p_fit->Draw("histcsame");
		h_asym_m_fit->SetLineColor(kGreen+3);
		h_asym_m_fit->SetLineStyle(kDashed);
		h_asym_m_fit->Draw("histcsame");

		TLegend leg2(0.4,0.75,0.6,0.92,"");
		leg2.SetLineStyle(0);
		leg2.SetLineColor(0);
		leg2.SetFillColor(0);
		leg2.SetTextFont(22);
		leg2.SetTextColor(1);
		leg2.SetTextSize(0.06);
		leg2.SetTextAlign(12);
		TLegendEntry* le = leg2.AddEntry(h_asym_p,"f=D_{s}^{-}K^{+}#pi^{+}#pi^{-}","l");
		le->SetTextColor(kBlack);                    
		le = leg2.AddEntry(h_asym_m,"#bar{f}=D_{s}^{+}K^{-}#pi^{-}#pi^{+}","l");
		le->SetTextColor(kGreen+3); 
		leg2.Draw();                   
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym.eps").c_str());	
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym.pdf").c_str());
	
		TH1D* h_asym_p_fit_unfolded = (TH1D*) h_N_unmixed_p_fit_unfolded->GetAsymmetry(h_N_mixed_p_fit_unfolded);	
		h_asym_p_fit_unfolded->SetLineWidth(3);
		h_asym_p_fit_unfolded->SetLineColor(kRed);
		h_asym_p_fit_unfolded->Draw("histc");
	
		TH1D* h_asym_m_fit_unfolded = (TH1D*) h_N_unmixed_m_fit_unfolded->GetAsymmetry(h_N_mixed_m_fit_unfolded);
		h_asym_m_fit_unfolded->SetLineWidth(3);
		h_asym_m_fit_unfolded->SetLineColor(kBlue);
		h_asym_m_fit_unfolded->Draw("histcsame");
		c->Print(((string)OutputDir+"h_asym_unfolded.eps").c_str());	
		c->Print(((string)OutputDir+"h_asym_unfolded.pdf").c_str());	
	
		TH1D* h_asym_p_unfolded = (TH1D*) h_N_unmixed_p_unfolded->GetAsymmetry(h_N_mixed_p_unfolded);	
		h_asym_p_unfolded->SetLineWidth(3);
		h_asym_p_unfolded->SetLineColor(kRed);
		h_asym_p_unfolded->Draw("e");
	
		TH1D* h_asym_m_unfolded = (TH1D*) h_N_unmixed_m_unfolded->GetAsymmetry(h_N_mixed_m_unfolded);
		h_asym_m_unfolded->SetLineWidth(3);
		h_asym_m_unfolded->SetLineColor(kBlue);
		h_asym_m_unfolded->Draw("esame");
	
		h_asym_p_fit_unfolded->Draw("histcsame");
		h_asym_m_fit_unfolded->Draw("histcsame");
		c->Print(((string)OutputDir+"h_asym_unfolded_fit.eps").c_str());	


// 		TH1D* h_asym_CP_f = (TH1D*) h_N_mixed->GetAsymmetry(h_N_unmixed);
// 		h_asym_CP_f->SetTitle(";t modulo (2#pi/#Deltam_{s}) (ps); A_{CP}^{#LTf#GT} ");	
// 		h_asym_CP_f->GetYaxis()->SetTitleOffset(0.9);	
// 		//h_asym_p->SetMinimum(-20);
// 		//h_asym_p->SetMaximum(20);
// 		TH1D* h_asym_CP_f_fit = (TH1D*) h_N_mixed_fit->GetAsymmetry(h_N_unmixed_fit);	
// 		h_asym_CP_f_fit->SetLineColor(kBlue+1);
// 		h_asym_CP_f_fit->SetLineWidth(5);
// 		h_asym_CP_f->Draw("e");
// 		h_asym_CP_f_fit->Draw("histcsame");
// 		graph->Draw("same");
// 		c->Print(((string)OutputDir+"h_asym_CP_f.eps").c_str());
// 		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym_CP_f.pdf").c_str());

		TH1D* h_asym_CP_f = (TH1D*) h_t_N->GetAsymmetry(h_t_Nbar);
		h_asym_CP_f->SetTitle(";t modulo (2#pi/#Deltam_{s}) (ps); A_{CP}^{#LTf#GT} ");	
		h_asym_CP_f->GetYaxis()->SetTitleOffset(0.9);	
		//h_asym_p->SetMinimum(-20);
		//h_asym_p->SetMaximum(20);
		TH1D* h_asym_CP_f_fit = (TH1D*) h_t_N_fit->GetAsymmetry(h_t_Nbar_fit);	
		h_asym_CP_f_fit->SetLineColor(kBlue+1);
		h_asym_CP_f_fit->SetLineWidth(5);
		h_asym_CP_f->Draw("e");
		h_asym_CP_f_fit->Draw("histcsame");
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym_CP_f.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym_CP_f.pdf").c_str());


		TH1D* h_asym_CP_i = (TH1D*) h_t_p->GetAsymmetry(h_t_m);	
		h_asym_CP_i->SetTitle(";t (ps); A_{CP}^{#LTi#GT} ");	
		h_asym_CP_i->GetYaxis()->SetTitleOffset(0.9);		
		//h_asym_p->SetMinimum(-20);
		//h_asym_p->SetMaximum(20);
		TH1D* h_asym_CP_i_fit = (TH1D*) h_t_p_fit->GetAsymmetry(h_t_m_fit);	
		h_asym_CP_i_fit->SetLineColor(kBlue+1);
		h_asym_CP_i_fit->SetLineWidth(5);
		h_asym_CP_i->Draw("e");
		h_asym_CP_i_fit->Draw("histcsame");
		graph->Draw("same");
		c->Print(((string)OutputDir+"h_asym_CP_i.eps").c_str());
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/timeFit/"+(string)OutputDir +"h_asym_CP_i.pdf").c_str());
     }
   }
   return;
}

void produceMarginalPdfs(){
    
    NamedParameter<string> InputFileName("InputFileName", (std::string) "/auto/data/dargent/BsDsKpipi/Final/signal_tagged.root", (char*) 0);
    NamedParameter<string> weightName("weightName", (std::string) "N_Bs_sw", (char*) 0);

    TString prefix = "";
    NamedParameter<double> min_year("min_year", 11);
    NamedParameter<double> max_year("max_year", 17);
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> min_TAUERR("min_TAUERR", 0.);
    NamedParameter<double> max_TAUERR("max_TAUERR", 0.1);    

    /// Load files
    // Data
    Int_t q_OS,f,Ds_ID,q_SS;
    Double_t w_OS,w_SS;
    double sw;
    int run,year,Ds_finalState,trigger;
    double t,dt;
    
    TChain* tree_norm=new TChain("DecayTree");
    tree_norm->Add( ((string)InputFileName).c_str());
    tree_norm->SetBranchStatus("*",0);
    tree_norm->SetBranchStatus("*sw*",1);
    tree_norm->SetBranchStatus("year",1);
    tree_norm->SetBranchStatus("*DEC",1);
    tree_norm->SetBranchStatus("*PROB",1);
    tree_norm->SetBranchStatus("*OS",1);
    tree_norm->SetBranchStatus("*TAU*",1);
    tree_norm->SetBranchStatus("run",1);
    tree_norm->SetBranchStatus("TriggerCat",1);
    tree_norm->SetBranchStatus("Ds_ID",1);

    tree_norm->SetBranchAddress("OS_Combination_DEC",&q_OS);
    tree_norm->SetBranchAddress("OS_Combination_PROB",&w_OS);
    tree_norm->SetBranchAddress("SS_Combination_DEC",&q_SS);
    tree_norm->SetBranchAddress("SS_Combination_PROB",&w_SS);
    tree_norm->SetBranchAddress(((string)weightName).c_str(),&sw);
    tree_norm->SetBranchAddress("year",&year);
    tree_norm->SetBranchAddress("run",&run);
    tree_norm->SetBranchAddress("Ds_finalState",&Ds_finalState);
    tree_norm->SetBranchAddress("B_FullDTF_TAU",&t);
    tree_norm->SetBranchAddress("B_FullDTF_TAUERR",&dt);
    tree_norm->SetBranchAddress("TriggerCat",&trigger);
    tree_norm->SetBranchAddress("Ds_ID",&Ds_ID);

    ///Make histograms
    int bins = 60;
    TH1D* h_w_OS_norm = new TH1D("h_w_OS_norm_comb","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run1 = new TH1D("h_w_OS_norm_Run1","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2 = new TH1D("h_w_OS_norm_Run2","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run1_t0 = new TH1D("h_w_OS_norm_Run1_t0","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2_t0 = new TH1D("h_w_OS_norm_Run2_t0","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run1_t1 = new TH1D("h_w_OS_norm_Run1_t1","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2_t1 = new TH1D("h_w_OS_norm_Run2_t1","; #eta_{OS}",bins,0,0.5);
    
    TH1D* h_w_SS_norm = new TH1D("h_w_SS_norm_comb","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run1 = new TH1D("h_w_SS_norm_Run1","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2 = new TH1D("h_w_SS_norm_Run2","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run1_t0 = new TH1D("h_w_SS_norm_Run1_t0","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2_t0 = new TH1D("h_w_SS_norm_Run2_t0","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run1_t1 = new TH1D("h_w_SS_norm_Run1_t1","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2_t1 = new TH1D("h_w_SS_norm_Run2_t1","; #eta_{SS}",bins,0,0.5);
    
    TH1D* h_q_OS_norm = new TH1D("h_q_OS_norm_comb","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run1 = new TH1D("h_q_OS_norm_Run1","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2 = new TH1D("h_q_OS_norm_Run2","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run1_t0 = new TH1D("h_q_OS_norm_Run1_t0","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2_t0 = new TH1D("h_q_OS_norm_Run2_t0","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run1_t1 = new TH1D("h_q_OS_norm_Run1_t1","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2_t1 = new TH1D("h_q_OS_norm_Run2_t1","; q_{OS}",3,-1.5,1.5);
    
    TH1D* h_q_SS_norm = new TH1D("h_q_SS_norm_comb","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run1 = new TH1D("h_q_SS_norm_Run1","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2 = new TH1D("h_q_SS_norm_Run2","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run1_t0 = new TH1D("h_q_SS_norm_Run1_t0","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2_t0 = new TH1D("h_q_SS_norm_Run2_t0","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run1_t1 = new TH1D("h_q_SS_norm_Run1_t1","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2_t1 = new TH1D("h_q_SS_norm_Run2_t1","; q_{SS}",3,-1.5,1.5);

    TH1D* h_q_f_norm = new TH1D("h_q_f_norm_comb","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run1 = new TH1D("h_q_f_norm_Run1","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2 = new TH1D("h_q_f_norm_Run2","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run1_t0 = new TH1D("h_q_f_norm_Run1_t0","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2_t0 = new TH1D("h_q_f_norm_Run2_t0","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run1_t1 = new TH1D("h_q_f_norm_Run1_t1","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2_t1 = new TH1D("h_q_f_norm_Run2_t1","; q_{f}",2,-2,2);
    
    TH1D* h_t_norm = new TH1D("h_t_norm_comb",";t (ps);Events (a.u.) ",bins,min_TAU,max_TAU);
    TH1D* h_t_norm_Run1 = new TH1D("h_t_norm_Run1",";t (ps);Events (a.u.) ",bins,min_TAU,max_TAU);
    TH1D* h_t_norm_Run2 = new TH1D("h_t_norm_Run2",";t (ps);Events (a.u.) ",bins,min_TAU,max_TAU);

    TH1D* h_dt_norm = new TH1D("h_dt_norm_comb",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run1 = new TH1D("h_dt_norm_Run1",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2 = new TH1D("h_dt_norm_Run2",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run1_t0 = new TH1D("h_dt_norm_Run1_t0",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_t0 = new TH1D("h_dt_norm_Run2_t0",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run1_t1 = new TH1D("h_dt_norm_Run1_t1",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_t1 = new TH1D("h_dt_norm_Run2_t1",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);

    TH1D* h_w_OS_norm_Run2_17 = new TH1D("h_w_OS_norm_Run2_17","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2_17_t0 = new TH1D("h_w_OS_norm_Run2_17_t0","; #eta_{OS}",bins,0,0.5);
    TH1D* h_w_OS_norm_Run2_17_t1 = new TH1D("h_w_OS_norm_Run2_17_t1","; #eta_{OS}",bins,0,0.5);
    
    TH1D* h_w_SS_norm_Run2_17 = new TH1D("h_w_SS_norm_Run2_17","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2_17_t0 = new TH1D("h_w_SS_norm_Run2_17_t0","; #eta_{SS}",bins,0,0.5);
    TH1D* h_w_SS_norm_Run2_17_t1 = new TH1D("h_w_SS_norm_Run2_17_t1","; #eta_{SS}",bins,0,0.5);
    
    TH1D* h_q_OS_norm_Run2_17 = new TH1D("h_q_OS_norm_Run2_17","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2_17_t0 = new TH1D("h_q_OS_norm_Run2_17_t0","; q_{OS}",3,-1.5,1.5);
    TH1D* h_q_OS_norm_Run2_17_t1 = new TH1D("h_q_OS_norm_Run2_17_t1","; q_{OS}",3,-1.5,1.5);
    
    TH1D* h_q_SS_norm_Run2_17 = new TH1D("h_q_SS_norm_Run2_17","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2_17_t0 = new TH1D("h_q_SS_norm_Run2_17_t0","; q_{SS}",3,-1.5,1.5);
    TH1D* h_q_SS_norm_Run2_17_t1 = new TH1D("h_q_SS_norm_Run2_17_t1","; q_{SS}",3,-1.5,1.5);

    TH1D* h_q_f_norm_Run2_17 = new TH1D("h_q_f_norm_Run2_17","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2_17_t0 = new TH1D("h_q_f_norm_Run2_17_t0","; q_{f}",2,-2,2);
    TH1D* h_q_f_norm_Run2_17_t1 = new TH1D("h_q_f_norm_Run2_17_t1","; q_{f}",2,-2,2);
    
    TH1D* h_t_norm_Run2_17 = new TH1D("h_t_norm_Run2_17",";t (ps);Events (a.u.) ",bins,min_TAU,max_TAU);

    TH1D* h_dt_norm_Run2_17 = new TH1D("h_dt_norm_Run2_17",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_17_t0 = new TH1D("h_dt_norm_Run2_17_t0",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_17_t1 = new TH1D("h_dt_norm_Run2_17_t1",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    
    ///loop over data events
    for(int i=0; i< tree_norm->GetEntries(); i++)
    {    
        tree_norm->GetEntry(i);
        if(year < min_year || year > max_year) continue;
	if(Ds_ID>0)f = -1;
	else f = 1;

        h_t_norm->Fill(t,sw);
        h_dt_norm->Fill(dt,sw);
        h_q_OS_norm->Fill((double)q_OS,sw);
        h_q_SS_norm->Fill((double)q_SS,sw);
        h_q_f_norm->Fill((double)f,sw);
        if(q_OS != 0)h_w_OS_norm->Fill(w_OS,sw);
        if(q_SS != 0)h_w_SS_norm->Fill(w_SS,sw);
            
        if(run==1){
            h_t_norm_Run1->Fill(t,sw);
            h_dt_norm_Run1->Fill(dt,sw);
            h_q_OS_norm_Run1->Fill((double)q_OS,sw);
            h_q_SS_norm_Run1->Fill((double)q_SS,sw);
	    h_q_f_norm_Run1->Fill((double)f,sw);
            if(q_OS != 0)h_w_OS_norm_Run1->Fill(w_OS,sw);
            if(q_SS != 0)h_w_SS_norm_Run1->Fill(w_SS,sw);
	    if(trigger == 0){
		h_dt_norm_Run1_t0->Fill(dt,sw);
            	h_q_OS_norm_Run1_t0->Fill((double)q_OS,sw);
            	h_q_SS_norm_Run1_t0->Fill((double)q_SS,sw);
        	h_q_f_norm_Run1_t0->Fill((double)f,sw);
	        if(q_OS != 0)h_w_OS_norm_Run1_t0->Fill(w_OS,sw);
                if(q_SS != 0)h_w_SS_norm_Run1_t0->Fill(w_SS,sw);
	    }
	    else if(trigger == 1){
		h_dt_norm_Run1_t1->Fill(dt,sw);
            	h_q_OS_norm_Run1_t1->Fill((double)q_OS,sw);
            	h_q_SS_norm_Run1_t1->Fill((double)q_SS,sw);
       		h_q_f_norm_Run1_t1->Fill((double)f,sw);
	        if(q_OS != 0)h_w_OS_norm_Run1_t1->Fill(w_OS,sw);
                if(q_SS != 0)h_w_SS_norm_Run1_t1->Fill(w_SS,sw);
	    }
        }
        else if(run==2){
            h_t_norm_Run2->Fill(t,sw);
            h_dt_norm_Run2->Fill(dt,sw);
            h_q_OS_norm_Run2->Fill((double)q_OS,sw);
            h_q_SS_norm_Run2->Fill((double)q_SS,sw);
       	    h_q_f_norm_Run2->Fill((double)f,sw);
            if(q_OS != 0)h_w_OS_norm_Run2->Fill(w_OS,sw);
            if(q_SS != 0)h_w_SS_norm_Run2->Fill(w_SS,sw);
	    if(trigger == 0){
		if(year < 17){
			h_dt_norm_Run2_t0->Fill(dt,sw);
			h_q_OS_norm_Run2_t0->Fill((double)q_OS,sw);
			h_q_SS_norm_Run2_t0->Fill((double)q_SS,sw);
			h_q_f_norm_Run2_t0->Fill((double)f,sw);
			if(q_OS != 0)h_w_OS_norm_Run2_t0->Fill(w_OS,sw);
			if(q_SS != 0)h_w_SS_norm_Run2_t0->Fill(w_SS,sw);
		}
		else {
			h_dt_norm_Run2_17_t0->Fill(dt,sw);
			h_q_OS_norm_Run2_17_t0->Fill((double)q_OS,sw);
			h_q_SS_norm_Run2_17_t0->Fill((double)q_SS,sw);
			h_q_f_norm_Run2_17_t0->Fill((double)f,sw);
			if(q_OS != 0)h_w_OS_norm_Run2_17_t0->Fill(w_OS,sw);
			if(q_SS != 0)h_w_SS_norm_Run2_17_t0->Fill(w_SS,sw);
		}
	    }
	    else if(trigger == 1){
		if(year < 17){
			h_dt_norm_Run2_t1->Fill(dt,sw);
			h_q_OS_norm_Run2_t1->Fill((double)q_OS,sw);
			h_q_SS_norm_Run2_t1->Fill((double)q_SS,sw);
			h_q_f_norm_Run2_t1->Fill((double)f,sw);
			if(q_OS != 0)h_w_OS_norm_Run2_t1->Fill(w_OS,sw);
			if(q_SS != 0)h_w_SS_norm_Run2_t1->Fill(w_SS,sw);
		}
		else {
			h_dt_norm_Run2_17_t1->Fill(dt,sw);
			h_q_OS_norm_Run2_17_t1->Fill((double)q_OS,sw);
			h_q_SS_norm_Run2_17_t1->Fill((double)q_SS,sw);
			h_q_f_norm_Run2_17_t1->Fill((double)f,sw);
			if(q_OS != 0)h_w_OS_norm_Run2_17_t1->Fill(w_OS,sw);
			if(q_SS != 0)h_w_SS_norm_Run2_17_t1->Fill(w_SS,sw);
		}
	    }
        }
       
    }
    
    TFile* out = new TFile("Mistag_pdfs.root","RECREATE");
    h_t_norm->Write();
    h_dt_norm->Write();
    h_q_OS_norm->Write();
    h_w_OS_norm->Write();
    h_q_SS_norm->Write();
    h_w_SS_norm->Write();
    h_q_f_norm->Write();

    h_t_norm_Run1->Write();
    h_dt_norm_Run1->Write();
    h_q_OS_norm_Run1->Write();
    h_w_OS_norm_Run1->Write();
    h_q_SS_norm_Run1->Write();
    h_w_SS_norm_Run1->Write();
    h_q_f_norm_Run1->Write();

    h_dt_norm_Run1_t0->Write();
    h_q_OS_norm_Run1_t0->Write();
    h_w_OS_norm_Run1_t0->Write();
    h_q_SS_norm_Run1_t0->Write();
    h_w_SS_norm_Run1_t0->Write();
    h_q_f_norm_Run1_t0->Write();

    h_dt_norm_Run1_t1->Write();
    h_q_OS_norm_Run1_t1->Write();
    h_w_OS_norm_Run1_t1->Write();
    h_q_SS_norm_Run1_t1->Write();
    h_w_SS_norm_Run1_t1->Write();
    h_q_f_norm_Run1_t1->Write();

    h_t_norm_Run2->Write();
    h_dt_norm_Run2->Write();
    h_q_OS_norm_Run2->Write();
    h_w_OS_norm_Run2->Write();
    h_q_SS_norm_Run2->Write();
    h_w_SS_norm_Run2->Write();
    h_q_f_norm_Run2->Write();

    h_dt_norm_Run2_t0->Write();
    h_q_OS_norm_Run2_t0->Write();
    h_w_OS_norm_Run2_t0->Write();
    h_q_SS_norm_Run2_t0->Write();
    h_w_SS_norm_Run2_t0->Write();
    h_q_f_norm_Run2_t0->Write();

    h_dt_norm_Run2_t1->Write();
    h_q_OS_norm_Run2_t1->Write();
    h_w_OS_norm_Run2_t1->Write();
    h_q_SS_norm_Run2_t1->Write();
    h_w_SS_norm_Run2_t1->Write();
    h_q_f_norm_Run2_t1->Write();

    h_dt_norm_Run2_17_t0->Write();
    h_q_OS_norm_Run2_17_t0->Write();
    h_w_OS_norm_Run2_17_t0->Write();
    h_q_SS_norm_Run2_17_t0->Write();
    h_w_SS_norm_Run2_17_t0->Write();
    h_q_f_norm_Run2_17_t0->Write();

    h_dt_norm_Run2_17_t1->Write();
    h_q_OS_norm_Run2_17_t1->Write();
    h_w_OS_norm_Run2_17_t1->Write();
    h_q_SS_norm_Run2_17_t1->Write();
    h_w_SS_norm_Run2_17_t1->Write();
    h_q_f_norm_Run2_17_t1->Write();

    out->Write();
}

void calculateAverageReso(){

    NamedParameter<string> InputDir("InputDir", (std::string) "/auto/data/dargent/BsDsKpipi/Final/", (char*) 0);

    int run,year,Ds_finalState,trigger;
    double t,dt;
    double sw;
    double K_P;

    TChain* tree_norm=new TChain("DecayTree");
    tree_norm->Add( ((string)InputDir + "Data/norm_tagged.root").c_str());
    tree_norm->SetBranchStatus("*",0);
    tree_norm->SetBranchStatus("*TAU*",1);
    tree_norm->SetBranchStatus("N_Bs_sw",1);
    tree_norm->SetBranchStatus("year",1);
    tree_norm->SetBranchStatus("run",1);
    tree_norm->SetBranchStatus("TriggerCat",1);
    tree_norm->SetBranchStatus("K_minus_fromDs_P",1);
//     tree_norm->SetBranchStatus("K_plus_P",1);
    tree_norm->SetBranchStatus("*finalState*",1);

    tree_norm->SetBranchAddress("Bs_BsDTF_TAU",&t);
    tree_norm->SetBranchAddress("Bs_BsDTF_TAUERR",&dt);
    tree_norm->SetBranchAddress("run",&run);
    tree_norm->SetBranchAddress("TriggerCat",&trigger);
    tree_norm->SetBranchAddress("N_Bs_sw",&sw);
    tree_norm->SetBranchAddress("year",&year);
    tree_norm->SetBranchAddress("K_minus_fromDs_P",&K_P);
//     tree_norm->SetBranchAddress("K_plus_P",&K_P);
    tree_norm->SetBranchAddress("Ds_finalState",&Ds_finalState);

    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> min_TAUERR("min_TAUERR", 0.);
    NamedParameter<double> max_TAUERR("max_TAUERR", 0.1);    
    FitParameter  offset_sigma_dt_Run1("offset_sigma_dt_Run1",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run1("scale_sigma_dt_Run1",1,1.2,0.1);
    FitParameter  offset_sigma_dt_Run2("offset_sigma_dt_Run2",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run2("scale_sigma_dt_Run2",1,1.2,0.1);
    FitParameter  offset_sigma_dt_Run2_17("offset_sigma_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run2_17("scale_sigma_dt_Run2_17",1,1.2,0.1);

    int bins = 60;
    TH1D* h_dt_norm = new TH1D("h_dt_norm_comb",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run1 = new TH1D("h_dt_norm_Run1",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2 = new TH1D("h_dt_norm_Run2",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);
    TH1D* h_dt_norm_Run2_17 = new TH1D("h_dt_norm_Run2_17",";#sigma_{t} (ps);Events (a.u.) ",bins,min_TAUERR,max_TAUERR);

    TFile* f_det_asym_Run1 = new TFile("../Asymmetries/AsymmetryHistos/det_asym_Run1.root");
    TH1D* det_asym_Run1 = (TH1D*) f_det_asym_Run1->Get("det_asym_Run1");
    TH1D* det_asym_Run1_toy = (TH1D*) f_det_asym_Run1->Get("det_asym_Run1");
    TFile* f_det_asym_Run2 = new TFile("../Asymmetries/AsymmetryHistos/det_asym_Run2.root");
    TH1D* det_asym_Run2 = (TH1D*) f_det_asym_Run2->Get("det_asym_Run2");
    TH1D* det_asym_Run2_toy = (TH1D*) f_det_asym_Run2->Get("det_asym_Run2");

    double sigma_dt = 0.;
    double sigma_dt_Run1 = 0.;
    double sigma_dt_Run2 = 0.;
    double sigma_dt_Run2_17 = 0.;
    double N_Run1 = 0.;
    double N_Run2 = 0.;
    double N_Run2_17 = 0.;

    double det_asy_Run1= 0.;
    double det_asy_Run2= 0.;
    double det_asy_N_Run1= 0.;
    double det_asy_N_Run2= 0.;

    double n_4 = 0.;

    for(int i=0; i< tree_norm->GetEntries(); i++)
    {    
        tree_norm->GetEntry(i);
	if(t < min_TAU || t > max_TAU )continue;
	if( dt < min_TAUERR || dt > max_TAUERR )continue;

	if(run==1){
		sigma_dt += sw*(offset_sigma_dt_Run1 + scale_sigma_dt_Run1 * dt);
		sigma_dt_Run1 += sw*(offset_sigma_dt_Run1 + scale_sigma_dt_Run1 * dt);
		N_Run1 += sw;
		h_dt_norm_Run1->Fill(offset_sigma_dt_Run1 + scale_sigma_dt_Run1 * dt,sw);
	}
	else if(year < 17){
		sigma_dt += sw*(offset_sigma_dt_Run2 + scale_sigma_dt_Run2 * dt);
		sigma_dt_Run2 += sw*(offset_sigma_dt_Run2 + scale_sigma_dt_Run2 * dt);
		N_Run2 += sw;
		h_dt_norm_Run2->Fill(offset_sigma_dt_Run2 + scale_sigma_dt_Run2 * dt,sw);
	}
	else{
		sigma_dt += sw*(offset_sigma_dt_Run2_17 + scale_sigma_dt_Run2_17 * dt);
		sigma_dt_Run2_17 += sw*(offset_sigma_dt_Run2_17 + scale_sigma_dt_Run2_17 * dt);
		N_Run2_17 += sw;
		h_dt_norm_Run2_17->Fill(offset_sigma_dt_Run2_17 + scale_sigma_dt_Run2_17 * dt,sw);
	}

	if(Ds_finalState == 4)n_4 += sw;

	double det_asy= 0.;
	if(run==2){
			double kaon_p = min(69.99,max(3.01,K_P/1000.));
			det_asy = -0.0091;
// 			det_asy = (Ds_finalState == 4) ? 0. : det_asym_Run2->GetBinContent(det_asym_Run2->FindBin(kaon_p))/100.;
			det_asy = (Ds_finalState != 4) ? 0. : det_asym_Run2->GetBinContent(det_asym_Run2->FindBin(kaon_p))/100.;
			det_asy_Run2 += det_asy * sw;
			det_asy_N_Run2 += sw;
	}
	else{
			double kaon_p = min(99.99,max(5.01,K_P/1000.));
			det_asy = -0.01;
			det_asy = (Ds_finalState != 4) ? 0. : det_asym_Run1->GetBinContent(det_asym_Run1->FindBin(kaon_p))/100.;
// 			det_asy = (Ds_finalState == 4) ? 0. : det_asym_Run1->GetBinContent(det_asym_Run1->FindBin(kaon_p))/100.;
			det_asy_Run1 += det_asy * sw;
			det_asy_N_Run1 += sw;
	}
	

    }

    cout << "n4 = " << n_4/(N_Run1+N_Run2+N_Run2_17) <<  endl;

    cout << "det_asy_Run1 = " << det_asy_Run1/det_asy_N_Run1 <<  endl;
    cout << "det_asy_Run2 = " << det_asy_Run2/det_asy_N_Run2 <<  endl;

    cout << "Reso Run1 = " << sigma_dt_Run1/N_Run1 *1000.<< endl;
    cout << "Reso Run2 = " << sigma_dt_Run2/N_Run2 *1000.<< endl;
    cout << "Reso Run2_17 = " << sigma_dt_Run2_17/N_Run2_17 *1000.<< endl;
    cout << "Reso tot = " << sigma_dt/(N_Run1+N_Run2+N_Run2_17) *1000.<< endl;

return;

    TCanvas* c = new TCanvas();
    h_dt_norm_Run1->DrawNormalized("e1",1);
    h_dt_norm_Run2->SetMarkerColor(kRed);
    h_dt_norm_Run2->SetLineColor(kRed);
    h_dt_norm_Run2_17->SetMarkerColor(kBlue);
    h_dt_norm_Run2_17->SetLineColor(kBlue);
    h_dt_norm_Run2->DrawNormalized("histcsame",1);
    h_dt_norm_Run2_17->DrawNormalized("histcsame",1);
    c->Print("dt_calib.eps");

// 	return;

    int N_toys = 100;
    double rho = -0.9 ;    
    TMatrixDSym cov(2);
    cov(0,0) = 1.;
    cov(0,1) = cov(1,0) = rho;    
    cov(1,1) = 1.;    

    RooArgList x_Run1;
    RooArgList m_Run1;
    RooRealVar* off_Run1 = new  RooRealVar("off_Run1", "off_Run1", offset_sigma_dt_Run1);
    RooRealVar* scale_Run1 = new  RooRealVar("scale_Run1", "scale_Run1", scale_sigma_dt_Run1);
    x_Run1.add(*off_Run1); 
    x_Run1.add(*scale_Run1);
    m_Run1.add(RooRealConstant::value(offset_sigma_dt_Run1));
    m_Run1.add(RooRealConstant::value(scale_sigma_dt_Run1));
        
    vector<double> sigma_Run1;    
    sigma_Run1.push_back(1.5/1000.);
    sigma_Run1.push_back(0.042);    
    
    TMatrixDSym cov_Run1(cov);
    for(int i=0; i < cov.GetNcols(); i++){
        for(int j=0; j < cov.GetNcols(); j++){    
            cov_Run1(i,j) = cov(i,j) * sigma_Run1[i] * sigma_Run1[j];
        }
    }

    RooMultiVarGaussian* gauss_cov_Run1 = new RooMultiVarGaussian("gauss_cov_Run1","gauss_cov_Run1",x_Run1, m_Run1, cov_Run1);
    RooDataSet* data_Run1 = gauss_cov_Run1->generate(x_Run1,N_toys);


    RooArgList x_Run2;
    RooArgList m_Run2;
    RooRealVar* off_Run2 = new  RooRealVar("off_Run2", "off_Run2", offset_sigma_dt_Run2);
    RooRealVar* scale_Run2 = new  RooRealVar("scale_Run2", "scale_Run2", scale_sigma_dt_Run2);
    x_Run2.add(*off_Run2); 
    x_Run2.add(*scale_Run2);
    m_Run2.add(RooRealConstant::value(offset_sigma_dt_Run2));
    m_Run2.add(RooRealConstant::value(scale_sigma_dt_Run2));
        
    vector<double> sigma_Run2;    
    sigma_Run2.push_back(1.5/1000.);
    sigma_Run2.push_back(0.042);    
    
    TMatrixDSym cov_Run2(cov);
    for(int i=0; i < cov.GetNcols(); i++){
        for(int j=0; j < cov.GetNcols(); j++){    
            cov_Run2(i,j) = cov(i,j) * sigma_Run2[i] * sigma_Run2[j];
        }
    }

    RooMultiVarGaussian* gauss_cov_Run2 = new RooMultiVarGaussian("gauss_cov_Run2","gauss_cov_Run2",x_Run2, m_Run2, cov_Run2);
    RooDataSet* data_Run2 = gauss_cov_Run2->generate(x_Run2,N_toys);


    RooArgList x_Run2_17;
    RooArgList m_Run2_17;
    RooRealVar* off_Run2_17 = new  RooRealVar("off_Run2_17", "off_Run2_17", offset_sigma_dt_Run2_17);
    RooRealVar* scale_Run2_17 = new  RooRealVar("scale_Run2_17", "scale_Run2_17", scale_sigma_dt_Run2_17);
    x_Run2_17.add(*off_Run2_17); 
    x_Run2_17.add(*scale_Run2_17);
    m_Run2_17.add(RooRealConstant::value(offset_sigma_dt_Run2_17));
    m_Run2_17.add(RooRealConstant::value(scale_sigma_dt_Run2_17));
        
    vector<double> sigma_Run2_17;    
    sigma_Run2_17.push_back(1.5/1000.);
    sigma_Run2_17.push_back(0.042);    
    
    TMatrixDSym cov_Run2_17(cov);
    for(int i=0; i < cov.GetNcols(); i++){
        for(int j=0; j < cov.GetNcols(); j++){    
            cov_Run2_17(i,j) = cov(i,j) * sigma_Run2_17[i] * sigma_Run2_17[j];
        }
    }

    RooMultiVarGaussian* gauss_cov_Run2_17 = new RooMultiVarGaussian("gauss_cov_Run2_17","gauss_cov_Run2_17",x_Run2_17, m_Run2_17, cov_Run2_17);
    RooDataSet* data_Run2_17 = gauss_cov_Run2_17->generate(x_Run2_17,N_toys);

    
    vector<double> sigma_dt_vec,sigma_dt_Run1_vec,sigma_dt_Run2_vec,sigma_dt_Run2_17_vec;
    vector<double> det_asy_Run1_vec,det_asy_Run2_vec;
    sigma_dt_vec.push_back(sigma_dt/(N_Run1+N_Run2+N_Run2_17));       
    sigma_dt_Run1_vec.push_back(sigma_dt_Run1/N_Run1);       
    sigma_dt_Run2_vec.push_back(sigma_dt_Run2/N_Run2);       
    sigma_dt_Run2_17_vec.push_back(sigma_dt_Run2_17/N_Run2_17);   
    det_asy_Run1_vec.push_back(det_asy_Run1/det_asy_N_Run1);
    det_asy_Run2_vec.push_back(det_asy_Run2/det_asy_N_Run2);  

    for(Int_t j=0;j<N_toys;j++) 
    {
        RooArgList* l_Run1 = (RooArgList*) data_Run1->get(j);    
        double off_Run1_val = ((RooRealVar *) l_Run1->at(0))->getVal();
        double scale_Run1_val = ((RooRealVar *) l_Run1->at(1))->getVal();

        RooArgList* l_Run2 = (RooArgList*) data_Run2->get(j);    
        double off_Run2_val = ((RooRealVar *) l_Run2->at(0))->getVal();
        double scale_Run2_val = ((RooRealVar *) l_Run2->at(1))->getVal();

        RooArgList* l_Run2_17 = (RooArgList*) data_Run2_17->get(j);    
        double off_Run2_17_val = ((RooRealVar *) l_Run2_17->at(0))->getVal();
        double scale_Run2_17_val = ((RooRealVar *) l_Run2_17->at(1))->getVal();

	sigma_dt = 0.;
	sigma_dt_Run1 = 0.;
	sigma_dt_Run2 = 0.;
	sigma_dt_Run2_17 = 0.;

	det_asy_Run1= 0.;
    	det_asy_Run2= 0.;
    	det_asy_N_Run1= 0.;
    	det_asy_N_Run2= 0.;

	for(int n = 1; n <= det_asym_Run1->GetNbinsX(); n++){
		det_asym_Run1_toy->SetBinContent(n, det_asym_Run1->GetBinContent(n) + gRandom->Gaus(0.,det_asym_Run1->GetBinError(n)) );
	}

	for(int n = 1; n <= det_asym_Run2->GetNbinsX(); n++){
		det_asym_Run2_toy->SetBinContent(n, det_asym_Run2->GetBinContent(n) + gRandom->Gaus(0.,det_asym_Run2->GetBinError(n)) );
	}


	
	for(int i=0; i< tree_norm->GetEntries(); i++)
	{    
		tree_norm->GetEntry(i);
		if(t < min_TAU || t > max_TAU )continue;
		if( dt < min_TAUERR || dt > max_TAUERR )continue;
	
		if(run==1){
			sigma_dt += sw*(off_Run1_val + scale_Run1_val * dt);
			sigma_dt_Run1 += sw*(off_Run1_val + scale_Run1_val * dt);
		}
		else if(year < 17){
			sigma_dt += sw*(off_Run2_val + scale_Run2_val * dt);
			sigma_dt_Run2 += sw*(off_Run2_val + scale_Run2_val * dt);
		}
		else{
			sigma_dt += sw*(off_Run2_17_val + scale_Run2_17_val * dt);
			sigma_dt_Run2_17 += sw*(off_Run2_17_val + scale_Run2_17_val * dt);
		}

		double det_asy= 0.;
		if(run==2){
			double kaon_p = min(69.99,max(3.01,K_P/1000.));
			det_asy = -0.0091;
// 			det_asy = (Ds_finalState == 4) ? 0. : det_asym_Run2->GetBinContent(det_asym_Run2->FindBin(kaon_p))/100.;
			det_asy = (Ds_finalState != 4) ? 0. : det_asym_Run2->GetBinContent(det_asym_Run2->FindBin(kaon_p))/100.;
			if(Ds_finalState == 4)det_asy_Run2 += det_asy * sw;
			if(Ds_finalState == 4)det_asy_N_Run2 += sw;
		}
		else{
			double kaon_p = min(99.99,max(5.01,K_P/1000.));
			det_asy = -0.01;
// 			det_asy = (Ds_finalState == 4) ? 0. : det_asym_Run1->GetBinContent(det_asym_Run1->FindBin(kaon_p))/100.;
			det_asy = (Ds_finalState != 4) ? 0. : det_asym_Run1->GetBinContent(det_asym_Run1->FindBin(kaon_p))/100.;
			if(Ds_finalState == 4)det_asy_Run1 += det_asy * sw;
			if(Ds_finalState == 4)det_asy_N_Run1 += sw;
		}

    	}
	sigma_dt_vec.push_back(sigma_dt/(N_Run1+N_Run2+N_Run2_17));       
	sigma_dt_Run1_vec.push_back(sigma_dt_Run1/N_Run1);       
	sigma_dt_Run2_vec.push_back(sigma_dt_Run2/N_Run2);       
	sigma_dt_Run2_17_vec.push_back(sigma_dt_Run2_17/N_Run2_17);   
	det_asy_Run1_vec.push_back(det_asy_Run1/det_asy_N_Run1);
	det_asy_Run2_vec.push_back(det_asy_Run2/det_asy_N_Run2);     
    }
    
    double sigma_dt_var = 0;
    double sigma_dt_var_Run1 = 0;
    double sigma_dt_var_Run2 = 0;
    double sigma_dt_var_Run2_17 = 0;

    double sigma_det_asy_Run1 = 0;
    double sigma_det_asy_Run2 = 0;

    for(int j=0; j < N_toys; j++){
        for(int k=j+1; k < N_toys; k++){
            sigma_dt_var += pow(sigma_dt_vec[j] - sigma_dt_vec[k],2);
            sigma_dt_var_Run1 += pow(sigma_dt_Run1_vec[j] - sigma_dt_Run1_vec[k],2);
            sigma_dt_var_Run2 += pow(sigma_dt_Run2_vec[j] - sigma_dt_Run2_vec[k],2);
            sigma_dt_var_Run2_17 += pow(sigma_dt_Run2_17_vec[j] - sigma_dt_Run2_17_vec[k],2);

            sigma_det_asy_Run1 += pow(det_asy_Run1_vec[j] - det_asy_Run1_vec[k],2);
            sigma_det_asy_Run2 += pow(det_asy_Run2_vec[j] - det_asy_Run2_vec[k],2);
        }
    }
    cout << "Reso Run1 = " << sigma_dt_Run1/N_Run1 *1000.<< " +- " <<  sqrt(sigma_dt_var_Run1/((double)N_toys*((double)N_toys-.1))) *1000 << endl;
    cout << "Reso Run2 = " << sigma_dt_Run2/N_Run2 *1000.<< " +- " <<  sqrt(sigma_dt_var_Run2/((double)N_toys*((double)N_toys-.1))) *1000 << endl;
    cout << "Reso Run2_17 = " << sigma_dt_Run2_17/N_Run2_17 *1000.<< " +- " <<  sqrt(sigma_dt_var_Run2_17/((double)N_toys*((double)N_toys-.1))) *1000 << endl;
    cout << "Reso tot = " << sigma_dt/(N_Run1+N_Run2+N_Run2_17) *1000.<< " +- " <<  sqrt(sigma_dt_var/((double)N_toys*((double)N_toys-.1))) *1000 << endl;

    cout << "DA Run1 = " <<  det_asy_Run1_vec[0] *100.<< " +- " <<  sqrt(sigma_det_asy_Run1/((double)N_toys*((double)N_toys-.1))) *100 << endl;
    cout << "DA Run2 = " <<  det_asy_Run2_vec[0] *100.<< " +- " <<  sqrt(sigma_det_asy_Run2/((double)N_toys*((double)N_toys-.1))) *100 << endl;
}

void calculateSweightsForToys(int step = 0){

	/// Load file
        NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
	TString inFileName = ((string)OutputDir+"toys_"+anythingToString((int)step)+".root").c_str();
	TString outFileName = ((string)OutputDir+"sw_toys_"+anythingToString((int)step)+".root").c_str();

	TFile* f = new TFile(inFileName,"Open");
	TTree* tree = (TTree*) f->Get("DecayTree");
	tree->SetBranchStatus("N_Bs_sw",0);

	TFile* out = new TFile(outFileName,"RECREATE");
	TTree* out_tree = tree->CopyTree("");

 	double sw;
     	TBranch* b_w = out_tree->Branch("N_Bs_sw", &sw, "N_Bs_sw/D");

	TString channelString = "m(D_{s}K#pi#pi)" ;
	RooRealVar DTF_Bs_M("Bs_DTF_MM", channelString, 5200, 5700,"MeV/c^{2}");
	RooArgList list =  RooArgList(DTF_Bs_M);
        RooDataSet* data = new RooDataSet("data","data",list,Import(*out_tree));
	
	/// Signal pdf
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

	/// Total pdf
	RooRealVar n_sig("n_sig", "n_sig", data->numEntries()*0.75, 0., data->numEntries());
	RooRealVar n_bkg("n_bkg", "n_bkg", data->numEntries()*0.25, 0., data->numEntries());

	RooAddPdf* pdf = new RooAddPdf("pdf", "pdf", RooArgList(*signal, *bkg), RooArgList(n_sig, n_bkg));

	/// Fit
	RooFitResult* result = pdf->fitTo(*data,Save(kTRUE),NumCPU(6),Extended(kTRUE));
	result->Print();

	/// Calculate sweights
	SPlot sPlot("sPlot","sPlot",*data,pdf,RooArgList(n_sig,n_bkg)); 
	
	for(int n = 0; n < tree->GetEntries(); n++){
			sw = sPlot.GetSWeight(n,"n_sig_sw");
			b_w->Fill();
	}
	
	/// Plotting
	TCanvas* c = new TCanvas();

	RooPlot* frame= DTF_Bs_M.frame();
	frame->SetTitle("");
 	data->plotOn(frame,Name("data"),Binning(200));
	pdf->plotOn(frame,Name("pdf"),LineColor(kBlue));
	pdf->plotOn(frame,Name("pdf"),LineColor(kRed),LineStyle(kDashed),Components("signal"));
	pdf->plotOn(frame,Name("pdf"),LineColor(kRed),LineStyle(kDashed),Components("bkg"));
	frame->Draw();
	c->Print("toyFit.eps");

	out_tree->Write();
	out->Close();
	f->Close();
}

int main(int argc, char** argv){

  NamedParameter<int>  addBkgToToys("addBkgToToys", 0);

  time_t startTime = time(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();
  gROOT->ProcessLine(".x ../lhcbStyle.C");

  produceMarginalPdfs();
  fullTimeFit(atoi(argv[1]),(string)argv[2]);
  if((string)argv[2] == "gen" && addBkgToToys)calculateSweightsForToys(atoi(argv[1]));

  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
  
  return 0;
}
//
