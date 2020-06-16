// Full td amplitude fit
// author: Philippe d'Argent
#include "Mint/FitParameter.h"
#include "Mint/NamedParameter.h"
#include "Mint/Minimiser.h"
#include "Mint/Neg2LL.h"
#include "Mint/Neg2LLSum.h"
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
#include "Neg2LLMultiConstraint.h"
#include "Mint/FitFractionList.h"
#include "Mint/FitFraction.h"
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
#include "RooTrace.h"
#include "RooConstVar.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooAddPdf.h"
#include "RooRandom.h"
#include "RooUniform.h"
#include "RooExponential.h"
#include "RooGaussian.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "Mint/RooCubicSplineFun.h"
#include "Mint/RooCubicSplineKnot.h"
#include "Mint/RooGaussEfficiencyModel.h"
#include "Mint/DecRateCoeff_Bd.h"
#include "Mint/TimePdfMaster.h"
#include "Mint/FullTimePdf.h"
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
#include <TF2.h>
#include "TProfile.h"
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
#include "Mint/HyperHistogram.h"
#include "Mint/LASSO.h"
#include "Mint/LASSO_flexi.h"
#include "RooStats/SPlot.h"

using namespace std;
using namespace RooFit ;
using namespace RooStats;
using namespace MINT;

class AmpsPdfFlexiFast
: public DalitzPdfBaseFlexiFastInteg
{
protected:
    TRandom* _localRnd;
    SignalGenerator* _sgGen;
    FromFileGenerator* _fileGen;
    IEventGenerator<IDalitzEvent>* _chosenGen;
    NamedParameter<std::string> _integratorSource;
    std::string _integratorFileName;
public:
    double un_normalised_noPs(IDalitzEvent& evt){
        double ampSq =  _amps->RealVal(evt);
        return ampSq;// * getEvent()->phaseSpace();
    }
    
    std::complex<double> ComplexVal_un_normalised_noPs(IDalitzEvent& evt){
        return  _amps->ComplexVal(evt);
    }
    
    AmpsPdfFlexiFast(const DalitzEventPattern& pat
                     , IFastAmplitudeIntegrable* amps
                     , MinuitParameterSet* mps
                     , double precision=1.e-4
                     , std::string method="efficient"
                     , std::string fname =  "SignalIntegrationEvents.root", bool generateNew = false, bool genMoreEvents = false
                     )
    : DalitzPdfBaseFlexiFastInteg(pat, 0, amps, precision, mps)
    , _localRnd(0)
    , _sgGen(0)
    , _fileGen(0)
    , _chosenGen(0)
    , _integratorSource("IntegratorSource", (std::string) "new", (char*) 0)
    , _integratorFileName(fname)
    {
        cout << " AmpsPdfFlexiFast with integ method " << method << endl;
        bool nonFlat = "efficient" == method;
        //bool generateNew = ((string)_integratorSource == (string)"new");
        if(nonFlat){
            cout << "AmpsPdfFlexiFast uses nonFlat integration." << endl;
            if(generateNew){
                _sgGen =  new SignalGenerator(pat);
                _sgGen->setWeighted();
                //_sgGen->setSaveEvents(_integratorFileName);
                _chosenGen = _sgGen;
            }else{
                // here, SignalGenerator is used by FromFileGenerator, to fill
                // up missing events in case more are needed than found in the
                // file.  Since we don't know with which random seed the
                // events in the file were generated, we supply a random
                // number generator with randomised seed.
                _localRnd = new TRandom3(time(0));
                _sgGen =  new SignalGenerator(pat, _localRnd);
                _sgGen->setWeighted();
                _sgGen->dontSaveEvents();// saving events is done by FromFileGenerator
                if(genMoreEvents) _fileGen   = new FromFileGenerator(_integratorFileName, _sgGen);
                else{
                    _fileGen = new FromFileGenerator(_integratorFileName, 0, "OPEN");
                    cout << "not going to generate any more events" << endl;
                }
                _chosenGen = _fileGen;
            }
            this->setEventGenerator(_chosenGen);
        }else{
            cout << "AmpsPdfFlexiFast uses flat integration." << endl;
        }
    }
    
    IFastAmplitudeIntegrable* getAmpSum(){ return _amps;}
    
    ~AmpsPdfFlexiFast(){
        if(0 != _fileGen)  delete _fileGen;
        if(0 != _sgGen)    delete _sgGen;
        if(0 != _localRnd) delete _localRnd;
    }
};

TString latexName(TString s){
    if(s == "Bs0_K_1__1270_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_1(1270) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_K_1__1270_p_D___Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_1(1270)[D] \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_K_1__1270_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_1(1270) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_K_1__1270_p__rho_1450_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_1(1270) \\to K \\, \\rho(1450) )$";
    if(s == "Bs0_K_1__1270_p__K_0_s_1430_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_1(1270) \\to K^{*}_{0}(1430) \\, \\pi )$";
    if(s == "Bs0_K_1__1400_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_1(1400) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_K_1__1400_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_1(1400) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_K_2_s_1430_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_2^{*}(1430) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_K_2_s_1430_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_2^{*}(1430) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_K_1460_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K(1460) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_K_1460_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K(1460) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_K_1460_p__sigma10__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K(1460) \\to K \\, \\sigma )$";
    if(s == "Bs0_Ks_1410_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_Ks_1410_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K^{*}(1410) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_Ks_1680_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K^{*}(1680) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_Ks_1680_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K^{*}(1680) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_K_2__1770_p__Ks_892_0__Kppim_pip_Dsm") return "$B_s \\to D_s \\, ( K_2(1770) \\to K^{*}(892) \\, \\pi )$";
    if(s == "Bs0_K_2__1770_p__rho_770_0__pippim_Kp_Dsm") return "$B_s \\to D_s \\, ( K_2(1770) \\to K \\, \\rho(770) )$";
    if(s == "Bs0_NonResS0__Dsmpip_Ks_892_0__Kppim_") return "$B_s \\to ( D_s \\, \\pi)_{S} \\, \\, K^{*}(892)$";
    if(s == "Bs0_NonResV0__Dsmpip_Ks_892_0__Kppim_") return "$B_s \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892)$";
    if(s == "Bs0_P__NonResV0__Dsmpip_Ks_892_0__Kppim_") return "$B_s[P] \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892)$";
    if(s == "Bs0_D__NonResV0__Dsmpip_Ks_892_0__Kppim_") return "$B_s[D] \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}(892)$";
    if(s == "Bs0_NonResS0__DsmKp_rho_770_0__pippim_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, \\rho(770)$";
    if(s == "Bs0_NonResV0__DsmKp_rho_770_0__pippim_") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, \\rho(770)$";
    if(s == "Bs0_P__NonResV0__DsmKp_rho_770_0__pippim_") return "$B_s[P] \\to ( D_s \\, K)_{P} \\, \\, \\rho(770)$";
    if(s == "Bs0_D__NonResV0__DsmKp_rho_770_0__pippim_") return "$B_s[D] \\to ( D_s \\, K)_{P} \\, \\, \\rho(770)$";
    if(s == "Bs0_NonResS0__DsmKp_sigma10__pippim_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, \\sigma$";
    if(s == "Bs0_NonResV0__DsmKp_sigma10__pippim_") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, \\sigma$";
    if(s == "Bs0_NonResS0__DsmKp_f_0__980_0__pippim_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, f_0(980)$";
    if(s == "Bs0_f_2__1270_0__pippim_NonResS0__DsmKp_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, f_2(1270)$";
    if(s == "Bs0_f_2__1270_0__pippim_NonResV0__DsmKp_") return "$B_s \\to ( D_s \\, K)_{P} \\, \\, f_2(1270)$";
    if(s == "Bs0_f_0__1370_0__pippim_NonResS0__DsmKp_") return "$B_s \\to ( D_s \\, K)_{S} \\, \\, f_0(1370)$";
    if(s == "Bs0_K_0_s_1430_0__Kppim_NonResS0__Dsmpip_") return "$B_s \\to ( D_s \\, \\pi)_{P} \\, \\, K^{*}_{0}(1430)$";
    
    return s;
}

vector<string> split(string s, string delimiter){
		size_t pos = 0;
		string token;
		vector<string> list;
		while ((pos = s.find(delimiter)) != string::npos) {
			token = s.substr(0, pos);
			list.push_back(token);
			s.erase(0, pos + delimiter.length());
		}
		return list;
}

std::vector<double> coherenceFactor(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DalitzEventList& eventList){
    
    cout << "Calculating coherence factor ..." << endl << endl;
    fas.print();
    fas_bar.print();
    
    std::complex<double> valK(0,0);
    double val1 = 0;
    double val2 = 0;
    
    const complex<double> phase_diff = polar(r, delta/360.*2*pi);

    for(unsigned int i=0; i<eventList.size(); i++){
        const std::complex<double> amp = fas.getVal(eventList[i]) ;
        const std::complex<double> amp_bar = fas_bar.getVal(eventList[i])*phase_diff ;
        valK += amp_bar*conj(amp)*eventList[i].getWeight()/ eventList[i].getGeneratorPdfRelativeToPhaseSpace();
        val1 += norm(amp)*eventList[i].getWeight()/ eventList[i].getGeneratorPdfRelativeToPhaseSpace();
        val2 += norm(amp_bar)*eventList[i].getWeight()/ eventList[i].getGeneratorPdfRelativeToPhaseSpace();
    }

    std::vector<double> result;
    result.push_back(sqrt(val2/val1));
    result.push_back(std::abs(valK)/sqrt(val1)/sqrt(val2));
    result.push_back(std::arg(valK)/(2.*pi)*360.);

    result.push_back((1.-pow(result[0],2))/(1.+pow(result[0],2)));
    result.push_back(-2. * valK.real()/val1 /(1.+pow(result[0],2)) );
    result.push_back(2. * valK.imag()/val1 /(1.+pow(result[0],2)) );
    
    cout << "r = " << result[0] << endl;
    cout << "k = " << result[1] << endl;
    cout << "phase = " << result[2] << " [deg]" << endl << endl;

    cout << "C = " << result[3] << endl;
    cout << "D = " << result[4] << endl;
    cout << "S = " << result[5] << endl;

    return result;
}

std::vector<double> coherenceFactor(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DiskResidentEventList& eventList, DalitzEventList& eventListData){
    
    cout << "Calculating coherence factor ..." << endl << endl;
    fas.print();
    fas_bar.print();
    
    std::complex<double> valK(0,0);
    double val1 = 0;
    double val2 = 0;
    
    const complex<double> phase_diff = polar(r, delta/360.*2*pi);
    
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    int nBins = 20;
    vector<int> s123;
    s123.push_back(1);
    s123.push_back(2);
    s123.push_back(3);
       
    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);
        
    vector<int> s134;
    s134.push_back(1);
    s134.push_back(3);
    s134.push_back(4);
        
    vector<int> s124;
    s124.push_back(1);
    s124.push_back(2);
    s124.push_back(4);

    TH1D* m_Kpi = new TH1D("",";#left[m(K^{+} #pi^{-})#right] (GeV/c^{2});Events (a.u.) ",nBins,0.6,1.3);
    TH1D* r_Kpi = (TH1D*)m_Kpi->Clone();
    TH1D* k_Kpi = (TH1D*)m_Kpi->Clone();
    TH1D* rk_Kpi = (TH1D*)m_Kpi->Clone();

    TH1D* m_Kpipi = new TH1D("",";#left[m(K^{+} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (a.u.) ",nBins,1,2);
    TH1D* r_Kpipi = (TH1D*)m_Kpipi->Clone();
    TH1D* k_Kpipi = (TH1D*)m_Kpipi->Clone();
    TH1D* rk_Kpipi = (TH1D*)m_Kpipi->Clone();

    TH1D* m_pipi = new TH1D("",";#left[m(#pi^{+} #pi^{-})#right] (GeV/c^{2});Events (a.u.) ",nBins,0.2,1.3);
    TH1D* r_pipi = (TH1D*)m_pipi->Clone();
    TH1D* k_pipi = (TH1D*)m_pipi->Clone();
    TH1D* rk_pipi = (TH1D*)m_pipi->Clone();

    TH1D* m_Dspipi = new TH1D("",";#left[m(D_{s}^{-} #pi^{+} #pi^{-})#right] (GeV/c^{2});Events (a.u.) ",nBins,2,5.5);
    TH1D* r_Dspipi = (TH1D*)m_Dspipi->Clone();
    TH1D* k_Dspipi = (TH1D*)m_Dspipi->Clone();
    TH1D* rk_Dspipi = (TH1D*)m_Dspipi->Clone();

    TH1D* m_Dspi = new TH1D("",";#left[m(D_{s}^{-} #pi^{+})#right] (GeV/c^{2});Events (a.u.) ",nBins,1.5,5);
    TH1D* r_Dspi = (TH1D*)m_Dspi->Clone();
    TH1D* k_Dspi = (TH1D*)m_Dspi->Clone();
    TH1D* rk_Dspi = (TH1D*)m_Dspi->Clone();

    vector< complex<double> > valK_Kpipi(nBins+1,0);
    vector<double> val1_Kpipi(nBins+1,0);
    vector<double> val2_Kpipi(nBins+1,0);

    vector< complex<double> > valK_Kpi(nBins+1,0);
    vector<double> val1_Kpi(nBins+1,0);
    vector<double> val2_Kpi(nBins+1,0);

    vector< complex<double> > valK_pipi(nBins+1,0);
    vector<double> val1_pipi(nBins+1,0);
    vector<double> val2_pipi(nBins+1,0);

    vector< complex<double> > valK_Dspipi(nBins+1,0);
    vector<double> val1_Dspipi(nBins+1,0);
    vector<double> val2_Dspipi(nBins+1,0);

    vector< complex<double> > valK_Dspi(nBins+1,0);
    vector<double> val1_Dspi(nBins+1,0);
    vector<double> val2_Dspi(nBins+1,0);

    for(unsigned int i=0; i<eventList.size(); i++){
        DalitzEvent evt = eventList.getEvent(i);

        const std::complex<double> amp = fas.getVal(evt) ;
        const std::complex<double> amp_bar = fas_bar.getVal(evt)*phase_diff ;
        valK += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        val1 += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        val2 += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();

	int bin_Kpi = r_Kpi->FindBin(sqrt(evt.s(2,4)/(GeV*GeV))); 
	if(!(bin_Kpi < 0 || bin_Kpi > nBins)) {
	        valK_Kpi[bin_Kpi] += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        	val1_Kpi[bin_Kpi] += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        	val2_Kpi[bin_Kpi] += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	}
	int bin_Kpipi = r_Kpipi->FindBin(sqrt(evt.sij(s234)/(GeV*GeV))); 
	if(!(bin_Kpipi < 0 || bin_Kpipi > nBins)) {
	        valK_Kpipi[bin_Kpipi] += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        	val1_Kpipi[bin_Kpipi] += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	        val2_Kpipi[bin_Kpipi] += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	}
	int bin_pipi = r_pipi->FindBin(sqrt(evt.s(3,4)/(GeV*GeV))); 
	if(!(bin_pipi < 0 || bin_pipi > nBins)) {
	        valK_pipi[bin_pipi] += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        	val1_pipi[bin_pipi] += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	        val2_pipi[bin_pipi] += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	}
	int bin_Dspipi = r_Dspipi->FindBin(sqrt(evt.sij(s134)/(GeV*GeV))); 
	if(!(bin_Dspipi < 0 || bin_Dspipi > nBins)) {
	        valK_Dspipi[bin_Dspipi] += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        	val1_Dspipi[bin_Dspipi] += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	        val2_Dspipi[bin_Dspipi] += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	}
	int bin_Dspi = r_Dspi->FindBin(sqrt(evt.s(1,3)/(GeV*GeV))); 
	if(!(bin_Dspi < 0 || bin_Dspi > nBins)) {
	        valK_Dspi[bin_Dspi] += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        	val1_Dspi[bin_Dspi] += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	        val2_Dspi[bin_Dspi] += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
	}

    }

    for(unsigned int i=0; i<eventListData.size(); i++){
	m_Kpipi->Fill(sqrt(eventListData[i].sij(s234)/(GeV*GeV)),eventListData[i].getWeight()) ;
	m_Kpi->Fill(sqrt(eventListData[i].s(2,4)/(GeV*GeV)),eventListData[i].getWeight()) ;
	m_pipi->Fill(sqrt(eventListData[i].s(3,4)/(GeV*GeV)),eventListData[i].getWeight()) ;
	m_Dspipi->Fill(sqrt(eventListData[i].sij(s134)/(GeV*GeV)),eventListData[i].getWeight()) ;
	m_Dspi->Fill(sqrt(eventListData[i].s(1,3)/(GeV*GeV)),eventListData[i].getWeight()) ;
    }

   for(int i = 1; i <= nBins; i++){
	if(val1_Kpipi[i]>0){
		r_Kpipi->SetBinContent(i,sqrt(val2_Kpipi[i]/val1_Kpipi[i]));
		k_Kpipi->SetBinContent(i,abs(valK_Kpipi[i])/sqrt(val1_Kpipi[i])/sqrt(val2_Kpipi[i]));
		rk_Kpipi->SetBinContent(i,r_Kpipi->GetBinContent(i)*k_Kpipi->GetBinContent(i));
	}
	if(val1_Kpi[i]>0){
		r_Kpi->SetBinContent(i,sqrt(val2_Kpi[i]/val1_Kpi[i]));
		k_Kpi->SetBinContent(i,abs(valK_Kpi[i])/sqrt(val1_Kpi[i])/sqrt(val2_Kpi[i]));
		rk_Kpi->SetBinContent(i,r_Kpi->GetBinContent(i)*k_Kpi->GetBinContent(i));
	}
	if(val1_pipi[i]>0){
		r_pipi->SetBinContent(i,sqrt(val2_pipi[i]/val1_pipi[i]));
		k_pipi->SetBinContent(i,abs(valK_pipi[i])/sqrt(val1_pipi[i])/sqrt(val2_pipi[i]));
		rk_pipi->SetBinContent(i,r_pipi->GetBinContent(i)*k_pipi->GetBinContent(i));
	}
	if(val1_Dspipi[i]>0){
		r_Dspipi->SetBinContent(i,sqrt(val2_Dspipi[i]/val1_Dspipi[i]));
		k_Dspipi->SetBinContent(i,abs(valK_Dspipi[i])/sqrt(val1_Dspipi[i])/sqrt(val2_Dspipi[i]));
		rk_Dspipi->SetBinContent(i,r_Dspipi->GetBinContent(i)*k_Dspipi->GetBinContent(i));
	}	
	if(val1_Dspi[i]>0){
		r_Dspi->SetBinContent(i,sqrt(val2_Dspi[i]/val1_Dspi[i]));
		k_Dspi->SetBinContent(i,abs(valK_Dspi[i])/sqrt(val1_Dspi[i])/sqrt(val2_Dspi[i]));
		rk_Dspi->SetBinContent(i,r_Dspi->GetBinContent(i)*k_Dspi->GetBinContent(i));
	}
    }	
    
    TCanvas* c = new TCanvas();

    m_Kpipi->SetFillColor(kGray+1);
    m_Kpipi->SetLineColor(kGray+1);
    m_Kpipi->Scale(1./m_Kpipi->GetMaximum());
    m_Kpipi->SetMaximum(1.2);
    m_Kpipi->SetMinimum(0);
    m_Kpipi->Draw("hist");

    r_Kpipi->SetLineColor(kRed);
    r_Kpipi->Scale(1./r_Kpipi->GetMaximum());
    r_Kpipi->Draw("histcsame");

    k_Kpipi->SetLineColor(kBlue);
    k_Kpipi->Scale(1./k_Kpipi->GetMaximum());
    k_Kpipi->Draw("histcsame");

    rk_Kpipi->SetLineColor(kGreen+3);
    rk_Kpipi->Scale(1./rk_Kpipi->GetMaximum());
    rk_Kpipi->Draw("histcsame");

    c->Print(((string)OutputDir+"coherence_Kpipi.eps").c_str());

    m_Kpi->SetFillColor(kGray+1);
    m_Kpi->SetLineColor(kGray+1);
    m_Kpi->Scale(1./m_Kpi->GetMaximum());
    m_Kpi->SetMaximum(1.2);
    m_Kpi->SetMinimum(0);    
    m_Kpi->Draw("hist");

    r_Kpi->SetLineColor(kRed);
    r_Kpi->Scale(1./r_Kpi->GetMaximum());
    r_Kpi->Draw("histcsame");

    k_Kpi->SetLineColor(kBlue);
    k_Kpi->Scale(1./k_Kpi->GetMaximum());
    k_Kpi->Draw("histcsame");

    rk_Kpi->SetLineColor(kGreen+3);
    rk_Kpi->Scale(1./rk_Kpi->GetMaximum());
    rk_Kpi->Draw("histcsame");

    c->Print(((string)OutputDir+"coherence_Kpi.eps").c_str());

    m_pipi->SetFillColor(kGray+1);
    m_pipi->SetLineColor(kGray+1);
    m_pipi->Scale(1./m_pipi->GetMaximum());
    m_pipi->SetMaximum(1.2);
    m_pipi->SetMinimum(0);
    m_pipi->Draw("hist");

    r_pipi->SetLineColor(kRed);
    r_pipi->Scale(1./r_pipi->GetMaximum());
    r_pipi->Draw("histcsame");

    k_pipi->SetLineColor(kBlue);
    k_pipi->Scale(1./k_pipi->GetMaximum());
    k_pipi->Draw("histcsame");

    rk_pipi->SetLineColor(kGreen+3);
    rk_pipi->Scale(1./rk_pipi->GetMaximum());
    rk_pipi->Draw("histcsame");

    c->Print(((string)OutputDir+"coherence_pipi.eps").c_str());

    m_Dspipi->SetFillColor(kGray+1);
    m_Dspipi->SetLineColor(kGray+1);
    m_Dspipi->Scale(1./m_Dspipi->GetMaximum());
    m_Dspipi->SetMaximum(1.2);
    m_Dspipi->SetMinimum(0);
    m_Dspipi->Draw("hist");

    r_Dspipi->SetLineColor(kRed);
    r_Dspipi->Scale(1./r_Dspipi->GetMaximum());
    r_Dspipi->Draw("histcsame");

    k_Dspipi->SetLineColor(kBlue);
    k_Dspipi->Scale(1./k_Dspipi->GetMaximum());
    k_Dspipi->Draw("histcsame");

    rk_Dspipi->SetLineColor(kGreen+3);
    rk_Dspipi->Scale(1./rk_Dspipi->GetMaximum());
    rk_Dspipi->Draw("histcsame");

    c->Print(((string)OutputDir+"coherence_Dspipi.eps").c_str());

    m_Dspi->SetFillColor(kGray+1);
    m_Dspi->SetLineColor(kGray+1);
    m_Dspi->Scale(1./m_Dspi->GetMaximum());
    m_Dspi->SetMaximum(1.2);
    m_Dspi->SetMinimum(0);
    m_Dspi->Draw("hist");

    r_Dspi->SetLineColor(kRed);
    r_Dspi->Scale(1./r_Dspi->GetMaximum());
    r_Dspi->Draw("histcsame");

    k_Dspi->SetLineColor(kBlue);
    k_Dspi->Scale(1./k_Dspi->GetMaximum());
    k_Dspi->Draw("histcsame");

    rk_Dspi->SetLineColor(kGreen+3);
    rk_Dspi->Scale(1./rk_Dspi->GetMaximum());
    rk_Dspi->Draw("histcsame");

    c->Print(((string)OutputDir+"coherence_Dspi.eps").c_str());

    std::vector<double> result;
    result.push_back(sqrt(val2/val1));
    result.push_back(std::abs(valK)/sqrt(val1)/sqrt(val2));
    result.push_back(std::arg(valK)/(2.*pi)*360.);

    result.push_back((1.-pow(result[0],2))/(1.+pow(result[0],2)));
    result.push_back(-2. * valK.real()/val1 /(1.+pow(result[0],2)) );
    result.push_back(2. * valK.imag()/val1 /(1.+pow(result[0],2)) );
    
    cout << "r = " << result[0] << endl;
    cout << "k = " << result[1] << endl;
    cout << "d = " << result[2] << " [deg]" << endl << endl;

    cout << "C = " << result[3] << endl;
    cout << "D = " << result[4] << endl;
    cout << "S = " << result[5] << endl;
    
    return result;
}

std::vector<double> coherenceFactor(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DiskResidentEventList& eventList, int N = -1){
    
    cout << "Calculating coherence factor ..." << endl << endl;
    //fas.print();
    //fas_bar.print();
    
    std::complex<double> valK(0,0);
    double val1 = 0;
    double val2 = 0;
    
    const complex<double> phase_diff = polar(r, delta/360.*2*pi);
    
    N = (N > 0) ? N : eventList.size();
    for(unsigned int i=0; i< N; i++){
        DalitzEvent evt = eventList.getEvent(i);
       
        const std::complex<double> amp = fas.getVal(evt) ;
        const std::complex<double> amp_bar = fas_bar.getVal(evt)*phase_diff ;
        valK += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        val1 += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        val2 += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
    }
    
    std::vector<double> result;
    result.push_back(sqrt(val2/val1));
    result.push_back(std::abs(valK)/sqrt(val1)/sqrt(val2));
    result.push_back(std::arg(valK)/(2.*pi)*360.);
    
    result.push_back((1.-pow(result[0],2))/(1.+pow(result[0],2)));
    result.push_back(-2. * valK.real()/val1 /(1.+pow(result[0],2)) );
    result.push_back(-2. * valK.imag()/val1 /(1.+pow(result[0],2)) );
/*
    cout << "r = " << result[0] << endl;
    cout << "k = " << result[1] << endl;
    cout << "d = " << result[2] << " [deg]" << endl << endl;

    cout << "C = " << result[3] << endl;
    cout << "Dbar = " << result[4]<< endl;
    cout << "Sbar = " << result[5] << endl;
*/  
    return result;
}

std::vector<double> coherenceFactor_CP(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DiskResidentEventList& eventList){
    
    cout << "Calculating coherence factor ..." << endl << endl;
    fas.print();
    fas_bar.print();
    
    std::complex<double> valK(0,0);
    double val1 = 0;
    double val2 = 0;
    
    const complex<double> phase_diff = polar(r, delta/360.*2*pi);
    
    for(unsigned int i=0; i<eventList.size(); i++){
        DalitzEvent evt = eventList.getEvent(i);
        evt.CP_conjugateYourself();
        
        const std::complex<double> amp = fas.getVal(evt) ;
        const std::complex<double> amp_bar = fas_bar.getVal(evt)*phase_diff ;
        valK += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        val1 += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
        val2 += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
    }
    
    std::vector<double> result;
    result.push_back(sqrt(val2/val1));
    result.push_back(std::abs(valK)/sqrt(val1)/sqrt(val2));
    result.push_back(std::arg(valK)/(2.*pi)*360.);
    
    result.push_back((1.-pow(result[0],2))/(1.+pow(result[0],2)));
    result.push_back(-2. * valK.real()/val1 /(1.+pow(result[0],2)) );
    result.push_back(-2. * valK.imag()/val1 /(1.+pow(result[0],2)) );

    cout << "r = " << result[0] << endl;
    cout << "k = " << result[1] << endl;
    cout << "d = " << result[2] << " [deg]" << endl << endl;

    cout << "C = " << result[3] << endl;
    cout << "Dbar = " << result[4]<< endl;
    cout << "Sbar = " << result[5] << endl;
    
    return result;
}

std::vector<double> coherenceFactor_CP(FitAmpSum& fas, FitAmpSum& fas_bar, double r, double delta, DalitzEventList& eventList){
    
        cout << "Calculating coherence factor ..." << endl << endl;
        fas.print();
        fas_bar.print();
    
        std::complex<double> valK(0,0);
        double val1 = 0;
        double val2 = 0;
    
        const complex<double> phase_diff = polar(r, delta/360.*2*pi);
    
        for(unsigned int i=0; i<eventList.size(); i++){
                DalitzEvent evt(eventList[i]);
                evt.CP_conjugateYourself();
                const std::complex<double> amp = fas.getVal(evt) ;
                const std::complex<double> amp_bar = fas_bar.getVal(evt)*phase_diff ;
                valK += amp_bar*conj(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
                val1 += norm(amp)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
                val2 += norm(amp_bar)*evt.getWeight()/ evt.getGeneratorPdfRelativeToPhaseSpace();
            }
	
	std::vector<double> result;
	result.push_back(sqrt(val2/val1));
	result.push_back(std::abs(valK)/sqrt(val1)/sqrt(val2));
	result.push_back(std::arg(valK)/(2.*pi)*360.);
	
	result.push_back((1.-pow(result[0],2))/(1.+pow(result[0],2)));
	result.push_back(-2. * valK.real()/val1 /(1.+pow(result[0],2)) );
	result.push_back(-2. * valK.imag()/val1 /(1.+pow(result[0],2)) );
	
	cout << "r = " << result[0] << endl;
	cout << "k = " << result[1] << endl;
	cout << "d = " << result[2] << " [deg]" << endl << endl;
	
	cout << "C = " << result[3] << endl;
	cout << "Dbar = " << result[4]<< endl;
	cout << "Sbar = " << result[5] << endl;
    
        return result;
}

// Full time-dependent PDF
class FullAmpsPdfFlexiFastCPV : 
	public MINT::PdfBase<IDalitzEvent>, virtual public IDalitzPdf{
    
protected:
    AmpsPdfFlexiFast* _amps;
    AmpsPdfFlexiFast* _amps_bar;

    AmpsPdfFlexiFast* _amps_CP;
    AmpsPdfFlexiFast* _amps_bar_CP;

    AmpsPdfFlexiFast* _ampsSum;
    AmpsPdfFlexiFast* _ampsSum_CP;

    double _intA;
    double _intAbar;
    double _intA_CP;
    double _intAbar_CP;
    complex<double> _intAAbar;   
    complex<double> _intAAbar_CP;   
 
    // Fit parameters
    const MINT::FitParameter& _r;
    const MINT::FitParameter& _delta;
    const MINT::FitParameter& _gamma;

    const MINT::FitParameter& _xm;
    const MINT::FitParameter& _ym;

    const MINT::FitParameter& _xp;
    const MINT::FitParameter& _yp;

    const MINT::FitParameter& _Gamma;
    const MINT::FitParameter& _dGamma;
    const MINT::FitParameter& _dm;
    
    const MINT::FitParameter& _offset_mean_dt;
    const MINT::FitParameter& _scale_mean_dt;
    const MINT::FitParameter& _scale_mean_2_dt;
    
    const MINT::FitParameter& _offset_sigma_dt;    
    const MINT::FitParameter& _scale_sigma_dt;
    const MINT::FitParameter& _scale_sigma_2_dt;

    const MINT::FitParameter& _offset_sigma2_dt;    
    const MINT::FitParameter& _scale_sigma2_dt;
    const MINT::FitParameter& _scale_sigma2_2_dt;

    const MINT::FitParameter& _offset_sigma3_dt;    
    const MINT::FitParameter& _scale_sigma3_dt;
    const MINT::FitParameter& _scale_sigma3_2_dt;

    const MINT::FitParameter& _offset_f_dt;    
    const MINT::FitParameter& _scale_f_dt;
    const MINT::FitParameter& _scale_f_2_dt;

    const MINT::FitParameter& _offset_f2_dt;    
    const MINT::FitParameter& _scale_f2_dt;
    const MINT::FitParameter& _scale_f2_2_dt;

    const MINT::FitParameter& _c0;
    const MINT::FitParameter& _c1;
    const MINT::FitParameter& _c2;
    const MINT::FitParameter& _c3;
    const MINT::FitParameter& _c4;
    const MINT::FitParameter& _c5;
    const MINT::FitParameter& _c6;
    const MINT::FitParameter& _c7;
    const MINT::FitParameter& _c8;
    const MINT::FitParameter& _c9;
    
    const MINT::FitParameter& _p0_os;
    const MINT::FitParameter& _p1_os;
    const MINT::FitParameter& _delta_p0_os;
    const MINT::FitParameter& _delta_p1_os;
    const MINT::FitParameter& _avg_eta_os;
    const MINT::FitParameter& _tageff_os;
    const MINT::FitParameter& _tageff_asym_os;
    
    const MINT::FitParameter& _p0_ss;
    const MINT::FitParameter& _p1_ss;
    const MINT::FitParameter& _delta_p0_ss;
    const MINT::FitParameter& _delta_p1_ss;
    const MINT::FitParameter& _avg_eta_ss;
    const MINT::FitParameter& _tageff_ss;
    const MINT::FitParameter& _tageff_asym_ss;
    
    const MINT::FitParameter& _production_asym;
    const MINT::FitParameter& _detection_asym;
    
    string _marginalPdfsPrefix;
    
// Time pdf master
    TimePdfMaster* _timePdfMaster;
    NamedParameter<int> _useCartCoord;
    NamedParameter<int> _directCPV;
    DalitzEventPattern _pat;
    DalitzEventPattern _pat_CP;

    // Limits
    NamedParameter<double> _min_TAU;
    NamedParameter<double> _max_TAU;
    NamedParameter<double> _min_TAUERR;
    NamedParameter<double> _max_TAUERR;

    // Toy generation
    FitAmpIncoherentSum* _fasGen;
    SignalGenerator* _sg;
    FitAmpIncoherentSum* _fasGen_CP;
    SignalGenerator* _sg_CP;
    vector<int> _s234;
    TF2 *_f_bkg;


public:
    void parametersChanged(){
        _ampsSum->parametersChanged();
        _intA = (_ampsSum->ComplexIntegralForTags(1,1)).real();
        _intAbar = (_ampsSum->ComplexIntegralForTags(-1,-1)).real();        
        _intAAbar = _ampsSum->ComplexIntegralForTags(1,-1);

	if(_directCPV){
		_ampsSum_CP->parametersChanged();
		_intA_CP = (_ampsSum_CP->ComplexIntegralForTags(1,1)).real();
		_intAbar_CP = (_ampsSum_CP->ComplexIntegralForTags(-1,-1)).real();   
		_intAAbar_CP = _ampsSum_CP->ComplexIntegralForTags(1,-1);
	}
	else {
		_intA_CP = _intA;
		_intAbar_CP = _intAbar;   
		_intAAbar_CP = _intAAbar;
	}
    }
    void beginFit(){
        _ampsSum->beginFit();
        _ampsSum->redoIntegrator();

	if(_directCPV){      
		_ampsSum_CP->beginFit();
		_ampsSum_CP->redoIntegrator();
	}
	parametersChanged(); 
        printIntegralVals();
	_timePdfMaster->listFitParDependencies();
    }
    void endFit(){
        printIntegralVals();
        _ampsSum->endFit();
        if(_directCPV)_ampsSum_CP->endFit();
    }
    
    void printIntegralVals(){
        cout << "intSum = " << _ampsSum->getIntegralValue() << endl;
        cout << "intA = " << _intA << endl;
        cout << "intAbar = " << _intAbar << endl;
        cout << "mag intAAbar = " << std::abs(_intAAbar) << endl;
	cout << "phase intAAbar = " << std::arg(_intAAbar)/(2.*pi)*360. << endl;

	if(_directCPV){
		cout << "intSum_CP = " << _ampsSum_CP->getIntegralValue() << endl;
		cout << "intA_CP = " << _intA_CP << endl;
		cout << "intAbar_CP = " << _intAbar_CP << endl;
		cout << "mag intAAbar_CP = " << std::abs(_intAAbar_CP) << endl;
		cout << "phase intAAbar_CP = " << std::arg(_intAAbar_CP)/(2.*pi)*360. << endl;
	}
    }

    inline double un_normalised_noPs(IDalitzEvent& evt){

        const double f = (double) evt.getValueFromVector(2);
        _timePdfMaster->setAllObservablesAndFitParameters(evt);
                
	complex<double> phase_delta_0 = polar(1.,-std::arg(_intAAbar));
	complex<double> phase_delta_0_CP = polar(1.,-std::arg(_intAAbar_CP));

        complex<double> phase_diff = (_useCartCoord) ? complex<double>((double)_xm,(double)_ym) : polar((double)_r,((double) _delta -(double)_gamma)/360.*2*pi);
	phase_diff *= phase_delta_0;

        complex<double> phase_diff_CP = (_useCartCoord) ? complex<double>((double)_xp,(double)_yp) : polar((double)_r,((double) _delta +(double)_gamma)/360.*2*pi);
	phase_diff_CP *= phase_delta_0_CP;

        complex<double> amp(0,0);
        complex<double> amp_bar(0,0);
	double norm_amp = 0.;

        complex<double> amp_CP(0,0);
        complex<double> amp_bar_CP(0,0);
        double norm_amp_CP = 0.;

        if(f>0){
            amp = _amps->ComplexVal_un_normalised_noPs(evt)/sqrt(_intA) ;
            amp_bar = phase_diff * _amps_bar->ComplexVal_un_normalised_noPs(evt)/sqrt(_intAbar) ;
            norm_amp = (norm(amp) + norm(amp_bar));
        }
        else {
            amp_CP = _amps_CP->ComplexVal_un_normalised_noPs(evt)/sqrt(_intA_CP) ;
            amp_bar_CP = phase_diff_CP * _amps_bar_CP->ComplexVal_un_normalised_noPs(evt)/sqrt(_intAbar_CP) ;
            norm_amp_CP = (norm(amp_CP) + norm(amp_bar_CP));
        }

        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster->setCP_coeff(
            		norm_amp,
			norm_amp_CP,
        		 ( norm(amp) - norm(amp_bar) ) ,
        		 -( norm(amp_CP) - norm(amp_bar_CP) ),
        		 -2.* real(amp_bar*conj(amp)) ,
        		 -2.* real(amp_bar_CP*conj(amp_CP)) ,
        		 2. * imag(amp_bar*conj(amp)) ,   /// - sign included in DecRateCoeff !!!
        		 -2. * imag(amp_bar_CP*conj(amp_CP)) /// - sign included in DecRateCoeff !!!
			);
        
        const double val = 
             ( _timePdfMaster->get_cosh_term_Val(evt)
             +  _timePdfMaster->get_cos_term_Val(evt)
             +  _timePdfMaster->get_sinh_term_Val(evt)
             +  _timePdfMaster->get_sin_term_Val(evt)
             )  ; //* _timePdfMaster->get_marginalPdfs_product(evt); //* _timePdfMaster->get_marginalPdfs_Val(evt);

        return val;
    }
    
    virtual double getVal(IDalitzEvent& evt){
        
        const double val = un_normalised_noPs(evt);
    
	complex<double> phase_delta_0 = polar(1.,-std::arg(_intAAbar));
	complex<double> phase_delta_0_CP = polar(1.,-std::arg(_intAAbar_CP));

        complex<double> phase_diff = 
		(_useCartCoord) ? complex<double>((double)_xm,(double)_ym) : polar((double)_r,((double) _delta -(double)_gamma)/360.*2*pi);
	phase_diff *= phase_delta_0;

        complex<double> phase_diff_CP = 
		(_useCartCoord) ? complex<double>((double)_xp,(double)_yp) : polar((double)_r,((double) _delta +(double)_gamma)/360.*2*pi);
	phase_diff_CP *= phase_delta_0_CP;

        const complex<double> int_interference =  _intAAbar  * phase_diff/sqrt(_intA)/sqrt(_intAbar) ;
        const complex<double> int_interference_CP = _intAAbar_CP  * phase_diff_CP/sqrt(_intA_CP)/sqrt(_intAbar_CP) ;

        // C,Cbar,D,Dbar,S,Sbar
         _timePdfMaster->setCP_coeff(
			(1. +  norm(phase_diff)),
			(1. +  norm(phase_diff_CP)),
        		(1. -  norm(phase_diff)),
        		-(1. -  norm(phase_diff_CP)) , 
        		(- int_interference.real() ), /// *2 included in integration !!!
        		(- int_interference_CP.real() ), 
        		(int_interference.imag() ),  /// - sign included in DecRateCoeff !!!
        		(- int_interference_CP.imag() )   
		 ); 
        
        double norm =  // (_intA + r* r * _intAbar)  *
        (     _timePdfMaster->get_cosh_term_Integral(evt)
            +  _timePdfMaster->get_cos_term_Integral(evt)
            +  _timePdfMaster->get_sinh_term_Integral(evt)
            +  _timePdfMaster->get_sin_term_Integral(evt) );

        return val/norm;
    }
    
    inline double getVal_timeIntegrated(IDalitzEvent& evt){

        const double f = (double) evt.getValueFromVector(2);
        _timePdfMaster->setAllObservablesToMean(evt);
                
	complex<double> phase_delta_0 = polar(1.,-std::arg(_intAAbar));
	complex<double> phase_delta_0_CP = polar(1.,-std::arg(_intAAbar_CP));

        complex<double> phase_diff = 
		(_useCartCoord) ? complex<double>((double)_xm,(double)_ym) : polar((double)_r,((double) _delta -(double)_gamma)/360.*2*pi);
	phase_diff *= phase_delta_0;

        complex<double> phase_diff_CP = 
		(_useCartCoord) ? complex<double>((double)_xp,(double)_yp) : polar((double)_r,((double) _delta +(double)_gamma)/360.*2*pi);
	phase_diff_CP *= phase_delta_0_CP;

        complex<double> amp(0,0);
        complex<double> amp_bar(0,0);
	double norm_amp = 0.;

        complex<double> amp_CP(0,0);
        complex<double> amp_bar_CP(0,0);
        double norm_amp_CP = 0.;

        if(f>0){
            amp = _amps->ComplexVal_un_normalised_noPs(evt)/sqrt(_intA) ;
            amp_bar = phase_diff * _amps_bar->ComplexVal_un_normalised_noPs(evt)/sqrt(_intAbar) ;
            norm_amp = (norm(amp) + norm(amp_bar));
        }
        else {
            amp_CP = _amps_CP->ComplexVal_un_normalised_noPs(evt)/sqrt(_intA_CP) ;
            amp_bar_CP = phase_diff_CP * _amps_bar_CP->ComplexVal_un_normalised_noPs(evt)/sqrt(_intAbar_CP) ;
            norm_amp_CP = (norm(amp_CP) + norm(amp_bar_CP));
        }
       
        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster->setCP_coeff(
            		norm_amp,
			norm_amp_CP,
        		 ( norm(amp) - norm(amp_bar) ) ,
        		 -( norm(amp_CP) - norm(amp_bar_CP) ),
        		 -2.* real(amp_bar*conj(amp)) ,
        		 -2.* real(amp_bar_CP*conj(amp_CP)) ,
        		 2. * imag(amp_bar*conj(amp)) ,   /// - sign included in DecRateCoeff !!!
        		 -2. * imag(amp_bar_CP*conj(amp_CP)) /// - sign included in DecRateCoeff !!!
			);
        
        const double val = 
        (     _timePdfMaster->get_cosh_term_Integral(evt)
            +  _timePdfMaster->get_cos_term_Integral(evt)
            +  _timePdfMaster->get_sinh_term_Integral(evt)
            +  _timePdfMaster->get_sin_term_Integral(evt) );

        return val;
    }

    inline double getAmpVal_timeIntegrated(IDalitzEvent& evt, FitAmpSum& fas, FitAmpSum& fas_bar,FitAmpSum& fas_CP, FitAmpSum& fas_bar_CP, std::vector<std::string> ampNames, bool CC = false){

        const double f = (double) evt.getValueFromVector(2);
        _timePdfMaster->setAllObservablesToMean(evt);
                
	complex<double> phase_delta_0 = polar(1.,-std::arg(_intAAbar));
	complex<double> phase_delta_0_CP = polar(1.,-std::arg(_intAAbar_CP));

        complex<double> phase_diff = 
		(_useCartCoord) ? complex<double>((double)_xm,(double)_ym) : polar((double)_r,((double) _delta -(double)_gamma)/360.*2*pi);
	phase_diff *= phase_delta_0;

        complex<double> phase_diff_CP = 
		(_useCartCoord) ? complex<double>((double)_xp,(double)_yp) : polar((double)_r,((double) _delta +(double)_gamma)/360.*2*pi);
	phase_diff_CP *= phase_delta_0_CP;

        complex<double> amp(0,0);
        complex<double> amp_bar(0,0);
	double norm_amp = 0.;

        complex<double> amp_CP(0,0);
        complex<double> amp_bar_CP(0,0);
        double norm_amp_CP = 0.;

        if(f>0){
            amp = fas.getAmpComplex(evt, ampNames, CC)/sqrt(_intA) ;
            amp_bar = phase_diff * fas_bar.getAmpComplex(evt, ampNames, CC)/sqrt(_intAbar) ;
            norm_amp = (norm(amp) + norm(amp_bar));
        }
        else {
            amp_CP = fas_CP.getAmpComplex(evt, ampNames, CC)/sqrt(_intA_CP) ;
            amp_bar_CP = phase_diff_CP * fas_bar_CP.getAmpComplex(evt, ampNames, CC)/sqrt(_intAbar_CP) ;
            norm_amp_CP = (norm(amp_CP) + norm(amp_bar_CP));
        }
       
        // C,Cbar,D,Dbar,S,Sbar
        _timePdfMaster->setCP_coeff(
            		norm_amp,
			norm_amp_CP,
        		 ( norm(amp) - norm(amp_bar) ) ,
        		 -( norm(amp_CP) - norm(amp_bar_CP) ),
        		 -2.* real(amp_bar*conj(amp)) ,
        		 -2.* real(amp_bar_CP*conj(amp_CP)) ,
        		 2. * imag(amp_bar*conj(amp)) ,   /// - sign included in DecRateCoeff !!!
        		 -2. * imag(amp_bar_CP*conj(amp_CP)) /// - sign included in DecRateCoeff !!!
			);
        
        const double val = 
        (     _timePdfMaster->get_cosh_term_Integral(evt)
            +  _timePdfMaster->get_cos_term_Integral(evt)
            +  _timePdfMaster->get_sinh_term_Integral(evt)
            +  _timePdfMaster->get_sin_term_Integral(evt) );

        return val;
    }

    virtual double getVal_withPs(IDalitzEvent& evt){return getVal(evt);}
    virtual double getVal_noPs(IDalitzEvent& evt){return getVal(evt);}
    
    virtual double getVal(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal(*evt);
    }
    virtual double getVal_withPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_withPs(*evt);
    }
    virtual double getVal_noPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_noPs(*evt);
    }

    vector<double> calculateCP_coeff(bool print = true){

	vector<double> result;
	if(_intA == -1)this->parametersChanged();
    
	double delta_0 = std::arg(_intAAbar)/(2.*pi)*360.;
	double delta_0_CP = std::arg(_intAAbar_CP)/(2.*pi)*360.;

	complex<double> phase_delta_0 = polar(1.,-std::arg(_intAAbar));
	complex<double> phase_delta_0_CP = polar(1.,-std::arg(_intAAbar_CP));

        complex<double> phase_diff = 
		(_useCartCoord) ? complex<double>((double)_xm,(double)_ym) : polar((double)_r,((double) _delta -(double)_gamma)/360.*2*pi);
	phase_diff *= phase_delta_0;

        complex<double> phase_diff_CP = 
		(_useCartCoord) ? complex<double>((double)_xp,(double)_yp) : polar((double)_r,((double) _delta +(double)_gamma)/360.*2*pi);
	phase_diff_CP *= phase_delta_0_CP;

        const complex<double> int_interference =  phase_diff * _intAAbar/sqrt(_intA)/sqrt(_intAbar) ;
        const complex<double> int_interference_CP = phase_diff_CP * _intAAbar_CP/sqrt(_intA_CP)/sqrt(_intAbar_CP) ;

        result.push_back((1. -  norm(phase_diff))/(1. +  norm(phase_diff)));
        result.push_back(-(1. -  norm(phase_diff_CP))/(1. +  norm(phase_diff_CP)));
 
        result.push_back((- int_interference.real() )/(1. +  norm(phase_diff)));
        result.push_back((- int_interference_CP.real() )/(1. +  norm(phase_diff_CP))); 

        result.push_back((int_interference.imag())/(1. +  norm(phase_diff)));
        result.push_back((- int_interference_CP.imag() )/(1. +  norm(phase_diff_CP)));   

        result.push_back(abs(_intAAbar/sqrt(_intA)/sqrt(_intAbar))/2.);   

	if(print){
		cout << "CP coeffs:: (C,Cbar,D,Dbar,S,Sbar,k) " << endl << result << endl;
		cout <<  abs(phase_diff) << endl;
		cout <<  arg(phase_diff)/(2*pi)*360. << endl;
		cout <<  abs(phase_diff_CP) << endl;
		cout <<  arg(phase_diff_CP)/(2*pi)*360. << endl<< endl;
	}
	return result;
    }

    double getCPcoeffChi2(vector<double>& coeff,double prec){
	vector<double> val = calculateCP_coeff(false);
	double chi2  = 0.;
	for(int i= 0; i< coeff.size();i++){
		chi2+= pow((coeff[i]-val[i])/prec,2);
	}
	return chi2;
    }

    std::pair<double, double> getCalibratedMistag_OS(IDalitzEvent& evt){
        return _timePdfMaster->getCalibratedMistag_OS(evt);
    }
    
    std::pair<double, double> getCalibratedMistag_SS(IDalitzEvent& evt){
        return _timePdfMaster->getCalibratedMistag_SS(evt);
    }

    std::pair<double, double> getCalibratedMistag(double eta,double avg_eta,double p0,double p1,double delta_p0,double delta_p1 ){
        return _timePdfMaster->getCalibratedMistag(eta, avg_eta, p0, p1, delta_p0, delta_p1);
    }
    
    double getCalibratedResolution(double dt){
        return _timePdfMaster->getCalibratedResolution(dt);
    }

    RooDataSet* sampleEvents(int N = 10000){
	return _timePdfMaster->sampleEvents(N);
    }

    inline double getSampledPdfVal(IDalitzEvent& evt){
	return _timePdfMaster->getSamplingPdfVal(evt);
    }


    DalitzEventList generateBkgToys(int N, int run = -1, int trigger = -1,string input = "/auto/data/dargent/BsDsKpipi/BDT/Data/signal_SS.root"){

	DalitzEventList eventList;
	NamedParameter<int>  applyPhspCuts("generateBkgToys::applyPhspCuts", 0);

	if(_f_bkg==0){
		NamedParameter<int>  correlate("readBkgData::correlate", 1);
	
		double t_bkg,dt_bkg,m_bkg;
		
		TChain* tree;
		tree=new TChain("DecayTree");
		tree->Add(((string)input).c_str());
		tree->SetBranchStatus("*",0);
		tree->SetBranchStatus("*DTF*",1);
		
		tree->SetBranchAddress("Bs_DTF_MM",&m_bkg);
		tree->SetBranchAddress("Bs_DTF_TAU",&t_bkg);
		tree->SetBranchAddress("Bs_DTF_TAUERR",&dt_bkg);
		
		TH2D* h_bkg = new TH2D("h_bkg","h_bkg",100,5200,5700,100,_min_TAU,_max_TAU);

		for(int i=0; i< tree->GetEntries(); i++)
		{	
			if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << tree->GetEntries() << endl;
			tree->GetEntry(i);

			h_bkg->Fill(m_bkg,t_bkg);	
		}

		TCanvas* c = new TCanvas();
		h_bkg->SetMarkerSize(0.1);
		h_bkg->Draw();
		TProfile* h_bkg_prof = h_bkg->ProfileX("p",1,-1,"");
		h_bkg_prof->SetLineColor(kRed);
		h_bkg_prof->SetMarkerColor(kRed);
		h_bkg_prof->Sumw2();
		h_bkg_prof->Rebin(4);
		cout << "bkg correlation = " << h_bkg->GetCorrelationFactor()*100 << endl;
		h_bkg_prof->Draw("histsame");
		c->Print("h_bkg.eps");

		 if(correlate)_f_bkg = new TF2("f_bkg","[0]*(1-exp(-[3]*(y-0.4)))*exp(-x*[1])*exp(-y*[2]*(1+(5200-x)*[4]))",5200,5700,_min_TAU,_max_TAU);
		 else _f_bkg = new TF2("f_bkg","[0]*(1-exp(-[3]*(y-0.4)))*exp(-x*[1])*exp(-y*[2])",5200,5700,_min_TAU,_max_TAU);
		_f_bkg->SetNpx(100);
		_f_bkg->SetNpy(100);

		_f_bkg->SetParLimits(0,1,10000000);
		_f_bkg->SetParLimits(1,0,1);
		_f_bkg->SetParLimits(2,0,10);
		_f_bkg->SetParLimits(3,-200,200);
		if(correlate)_f_bkg->SetParLimits(4,-0.1,0.1);

		_f_bkg->SetParameter(0,1000);
		_f_bkg->SetParameter(1,0.);
		_f_bkg->SetParameter(2,1.);
		_f_bkg->SetParameter(3,0.);
		if(correlate)_f_bkg->FixParameter(4,0.);

		h_bkg->Fit("f_bkg");
		if(correlate){ 
			_f_bkg->ReleaseParameter(4);
			h_bkg->Fit("f_bkg");
		}
		
		TH2D* h_fit = (TH2D*)_f_bkg->CreateHistogram();	
		h_fit->SetMarkerSize(0.1);
		h_fit->Draw();
		cout << "bkg pdf correlation = " << h_fit->GetCorrelationFactor()*100 << endl;
		c->Print("f_bkg.eps");
	}

	while(eventList.size()<N)
	{	

		double t,m;
		_f_bkg->GetRandom2(m,t);

		vector<double> marginal_vals = _timePdfMaster->getRandom_marginalVals();
		double dt = marginal_vals[0] ;
		double eta_OS = marginal_vals[1] ;
		double eta_SS = marginal_vals[2] ;
		int f = (gRandom->Uniform() > 0.5) ? 1 : -1;

		DalitzEvent evt;
		if(f < 0)evt = DalitzEvent(_pat,gRandom);
		else evt = DalitzEvent(_pat_CP,gRandom);	
		
		if(!(evt.phaseSpace() > 0.))continue;
		if(applyPhspCuts)if(sqrt(evt.sij(_s234)/(GeV*GeV)) > 1.95 || sqrt(evt.s(2,4)/(GeV*GeV)) > 1.2 || sqrt(evt.s(3,4)/(GeV*GeV)) > 1.2) continue;	
	
		evt.setWeight(1.);
		evt.setValueInVector(0, t);
		evt.setValueInVector(1, dt);   
		if(f<0)evt.setValueInVector(2, 1);
		else if(f > 0)evt.setValueInVector(2, -1);
		evt.setValueInVector(3, TMath::Nint(gRandom->Uniform(-1,1)));
		evt.setValueInVector(4, eta_OS);
		evt.setValueInVector(5, TMath::Nint(gRandom->Uniform(-1,1)));
		evt.setValueInVector(6, eta_SS);		
		evt.setValueInVector(7, run);
		evt.setValueInVector(8, trigger);
		evt.setValueInVector(9, m);

		eventList.Add(evt);
	}
	//saveEventListToFile(eventList);
	return eventList;
    }

   DalitzEventList readBkgData(string input = "/auto/data/dargent/BsDsKpipi/BDT/Data/signal_SS.root"){

        NamedParameter<int>  correlate("readBkgData::correlate", 1);

	double t,dt,mB;
	int f;
	int q_OS;
	double Bs_ID,Ds_ID;
	Int_t q_SS;
	double eta_OS;
	Double_t eta_SS;
	double sw;
	int year,run,Ds_finalState,trigger;
	double K[4];
	double pip[4];
	double pim[4];
	double Ds_Kp[4],Ds_Km[4],Ds_pim[4],Ds[4];
	
	TChain* tree;

	tree=new TChain("DecayTree");
	tree->Add(((string)input).c_str());
	tree->SetBranchStatus("*",0);
	tree->SetBranchStatus("year",1);
	tree->SetBranchStatus("*TAU*",1);
	tree->SetBranchStatus("*ID*",1);
	tree->SetBranchStatus("weight",1);
	tree->SetBranchStatus("TriggerCat",1);
	tree->SetBranchStatus("run",1);
	tree->SetBranchStatus("BsDTF_*P*",1);
	tree->SetBranchStatus("Bs_DTF_MM",1);
	
	tree->SetBranchAddress("Bs_DTF_MM",&mB);
	tree->SetBranchAddress("Bs_DTF_TAU",&t);
	tree->SetBranchAddress("Bs_DTF_TAUERR",&dt);
	tree->SetBranchAddress("Ds_ID",&f);
	tree->SetBranchAddress("year",&year);
	tree->SetBranchAddress("run",&run);
	tree->SetBranchAddress("TriggerCat",&trigger);
	tree->SetBranchAddress("BsDTF_Kplus_PX",&K[0]);
	tree->SetBranchAddress("BsDTF_Kplus_PY",&K[1]);
	tree->SetBranchAddress("BsDTF_Kplus_PZ",&K[2]);
	tree->SetBranchAddress("BsDTF_Kplus_PE",&K[3]);
	tree->SetBranchAddress("BsDTF_piplus_PX",&pip[0]);
	tree->SetBranchAddress("BsDTF_piplus_PY",&pip[1]);
	tree->SetBranchAddress("BsDTF_piplus_PZ",&pip[2]);
	tree->SetBranchAddress("BsDTF_piplus_PE",&pip[3]);    
	tree->SetBranchAddress("BsDTF_piminus_PX",&pim[0]);
	tree->SetBranchAddress("BsDTF_piminus_PY",&pim[1]);
	tree->SetBranchAddress("BsDTF_piminus_PZ",&pim[2]);
	tree->SetBranchAddress("BsDTF_piminus_PE",&pim[3]);    
	tree->SetBranchAddress("BsDTF_Ds_Kplus_PX",&Ds_Kp[0]);
	tree->SetBranchAddress("BsDTF_Ds_Kplus_PY",&Ds_Kp[1]);
	tree->SetBranchAddress("BsDTF_Ds_Kplus_PZ",&Ds_Kp[2]);
	tree->SetBranchAddress("BsDTF_Ds_Kplus_PE",&Ds_Kp[3]);    
	tree->SetBranchAddress("BsDTF_Ds_Kminus_PX",&Ds_Km[0]);
	tree->SetBranchAddress("BsDTF_Ds_Kminus_PY",&Ds_Km[1]);
	tree->SetBranchAddress("BsDTF_Ds_Kminus_PZ",&Ds_Km[2]);
	tree->SetBranchAddress("BsDTF_Ds_Kminus_PE",&Ds_Km[3]);
	tree->SetBranchAddress("BsDTF_Ds_piminus_PX",&Ds_pim[0]);
	tree->SetBranchAddress("BsDTF_Ds_piminus_PY",&Ds_pim[1]);
	tree->SetBranchAddress("BsDTF_Ds_piminus_PZ",&Ds_pim[2]);
	tree->SetBranchAddress("BsDTF_Ds_piminus_PE",&Ds_pim[3]);
	
	DalitzEventList eventList;

	for(int i=0; i< tree->GetEntries()-1; i++)
	{	
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << tree->GetEntries() << endl;
		tree->GetEntry(i);

		vector<double> marginal_vals = _timePdfMaster->getRandom_marginalVals();
		dt = marginal_vals[0] ;
		double eta_OS = marginal_vals[1] ;
		double eta_SS = marginal_vals[2] ;
		f = (gRandom->Uniform() > 0.5) ? 1 : -1;		
		
		double sign = 1.;
		TLorentzVector K_p(sign*K[0],sign*K[1],sign*K[2],K[3]);
		TLorentzVector pip_p(sign*pip[0],sign*pip[1],sign*pip[2],pip[3]);
		TLorentzVector pim_p(sign*pim[0],sign*pim[1],sign*pim[2],pim[3]);
		TLorentzVector D_p;
		TLorentzVector D_Kp_p(sign*Ds_Kp[0],sign*Ds_Kp[1],sign*Ds_Kp[2],Ds_Kp[3]);
		TLorentzVector D_Km_p(sign*Ds_Km[0],sign*Ds_Km[1],sign*Ds_Km[2],Ds_Km[3]);
		TLorentzVector D_pim_p(sign*Ds_pim[0],sign*Ds_pim[1],sign*Ds_pim[2],Ds_pim[3]);
		D_p = D_Kp_p + D_Km_p + D_pim_p;
		
		TLorentzVector B_p = K_p + pip_p + pim_p + D_p;
		// array of vectors
		vector<TLorentzVector> vectorOfvectors;
	
		vectorOfvectors.push_back(B_p*MeV);
		vectorOfvectors.push_back(D_p*MeV);
		vectorOfvectors.push_back(K_p*MeV);
		vectorOfvectors.push_back(pip_p*MeV);
		vectorOfvectors.push_back(pim_p*MeV);
		DalitzEvent evt;
	
		if(f < 0)evt = DalitzEvent(_pat, vectorOfvectors);
		else evt = DalitzEvent(_pat_CP, vectorOfvectors);
	
		if(correlate==2){
			if(f < 0)evt = DalitzEvent(_pat,gRandom);
			else evt = DalitzEvent(_pat_CP,gRandom);
		}

		if(correlate==3) t = gRandom->Exp(1./_Gamma*0.75);
	
// 		t = gRandom->Exp(1./_Gamma);
// 		if(t < _min_TAU || t > _max_TAU )continue;
/*		if( dt < min_TAUERR || dt > max_TAUERR )continue;
		if(year < min_year || year > max_year) continue;
		if(eta_SS < w_min || eta_SS > w_max )continue;*/
	
		if(!(evt.phaseSpace() > 0.))continue;
			
		evt.setWeight(1.);
		evt.setValueInVector(0, t);
		evt.setValueInVector(1, dt);   
		if(f<0)evt.setValueInVector(2, 1);
		else if(f > 0)evt.setValueInVector(2, -1);
		else {
			cout << "ERROR:: Undefined final state " << f << endl;  
			throw "ERROR";
		}
		evt.setValueInVector(3, TMath::Nint(gRandom->Uniform(-1,1)));
		evt.setValueInVector(4, eta_OS);
		evt.setValueInVector(5, TMath::Nint(gRandom->Uniform(-1,1)));
		evt.setValueInVector(6, eta_SS);
		evt.setValueInVector(7, run);
		evt.setValueInVector(8, trigger);

		if(correlate==0)tree->GetEntry(i+1);
		mB = mB - 200;
		if(mB <= 5200 || mB >= 5700 )continue;

		evt.setValueInVector(9, mB);
		eventList.Add(evt);
        }
	return eventList;
    }

    DalitzEventList generateBkgToysBootstrap(int N,int run = -1 , int trigger = -1,string input = "/auto/data/dargent/BsDsKpipi/BDT/Data/signal_SS.root"){

	DalitzEventList eventList;
	DalitzEventList eventListData = readBkgData(input);
	int N_sample = eventListData.size();

	vector<int> b_indices;
	while( b_indices.size() < N )b_indices.push_back(TMath::Nint(gRandom->Uniform(0,N_sample-1)));
	sort(b_indices.begin(), b_indices.end());

	for(int i=0; i< N; i++)
	{	
		DalitzEvent evt = eventListData[b_indices[i]];

		evt.setValueInVector(7, run);
		evt.setValueInVector(8, trigger);
		//if(t < _min_TAU || t > _max_TAU )continue;

		eventList.Add(evt);
	}
	//saveEventListToFile(eventList);
	return eventList;
    }

    DalitzEvent generateWeightedEvent(){

	while(true){
		double t_MC = gRandom->Exp(1./_Gamma);
                if(t_MC > _max_TAU || t_MC < _min_TAU)continue;

		int f_MC = (gRandom->Uniform() > 0.5) ? 1 : -1;		

		counted_ptr<IDalitzEvent> evtPtr;

		if(f_MC > 0)evtPtr = _sg->newEvent();
		else evtPtr = _sg_CP->newEvent();

		DalitzEvent evt(evtPtr.get());
                if(!((sqrt(evt.sij(_s234)/(GeV*GeV)) < 1.95 && sqrt(evt.s(2,4)/(GeV*GeV)) < 1.2 && sqrt(evt.s(3,4)/(GeV*GeV)) < 1.2)))continue;

		vector<double> marginal_vals = _timePdfMaster->getRandom_marginalVals();
		double dt_MC = marginal_vals[0] ;
		double eta_OS_MC = marginal_vals[1] ;
		double eta_SS_MC = marginal_vals[2] ;
	
		// true flavor
		int q_MC = (gRandom->Uniform() > 0.5) ? 1 : -1;

 	        int q_SS_MC = (gRandom->Uniform() > 2./3.) ? 0 : q_MC ;
         	int q_OS_MC = (gRandom->Uniform() > 2./3.) ? 0 : q_MC ;
 
		q_OS_MC = (gRandom->Uniform() < 0.5) ? - q_OS_MC : q_OS_MC;
		q_SS_MC = (gRandom->Uniform() < 0.5) ? - q_SS_MC : q_SS_MC;
		
		eta_OS_MC = (q_OS_MC == 0) ? 0.5 : eta_OS_MC;
		eta_SS_MC = (q_SS_MC == 0) ? 0.5 : eta_SS_MC;

// 		if(f_MC<0)evt.CP_conjugateYourself();
		evt.setValueInVector(0, t_MC);
		evt.setValueInVector(1, dt_MC);
		evt.setValueInVector(2, f_MC);
		evt.setValueInVector(3, q_OS_MC);
		evt.setValueInVector(4, eta_OS_MC);
		evt.setValueInVector(5, q_SS_MC);
		evt.setValueInVector(6, eta_SS_MC);
		evt.setValueInVector(9, gRandom->Gaus(5367,20));

		evt.setGeneratorPdfRelativeToPhaseSpace(_Gamma * evt.getGeneratorPdfRelativeToPhaseSpace() * (exp(-t_MC*_Gamma) / ( ( exp(-_min_TAU*_Gamma) - exp(-_max_TAU*_Gamma) ))));
		return evt;
	}
    }

    DalitzEventList generateToys(int N = 10000, int run = -1 , int trigger = -1){

	time_t startTime = time(0);

	cout << "Generating " << N << " events" << endl;
	DalitzEventList eventList;
	if(_fasGen==0){
		        /// Simple amplitude model for importance sampling
   		        NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
        		DalitzEventPattern pat(EventPattern.getVector());
		        DalitzEventPattern pat_CP = pat.makeCPConjugate();     		

			_fasGen = new FitAmpIncoherentSum((DalitzEventPattern)pat);
        		_fasGen->print();
			{
				DalitzEventList eventListPhsp;
				eventListPhsp.generatePhaseSpaceEvents(200000,pat);
				_fasGen->normalizeAmps(eventListPhsp);
			}
	
			_fasGen_CP = new FitAmpIncoherentSum(*_fasGen);
			_fasGen_CP->CPConjugateSameFitParameters();

			{
				DalitzEventList eventListPhsp_CP;
				eventListPhsp_CP.generatePhaseSpaceEvents(2,pat_CP);
				_fasGen_CP->getVal(eventListPhsp_CP[0]);
			}

		        _sg = new SignalGenerator(pat,_fasGen);
		        _sg_CP = new SignalGenerator(pat_CP,_fasGen_CP);
	}

	/// Estimate max val
	vector<double> vals;	
	for(int i = 0; i < 100000; i++){
		DalitzEvent evt = generateWeightedEvent();
		double val = getVal(evt)/evt.getGeneratorPdfRelativeToPhaseSpace(); ///_timePdfMaster->get_marginalPdfs_product(evt);
		vals.push_back(val);	
	}

	cout << "Now calculating maximum val " << vals.size() << endl;
	double amax,pmax;
	generalisedPareto_estimateMaximum(vals,0.999,amax,pmax);
	
	double pdf_max = 1.;
	if(!TMath::IsNaN(pmax) && pmax > 0 && pmax < 100 * amax)pdf_max = pmax;
	else if(!TMath::IsNaN(amax))pdf_max = amax;
	if(amax > pmax)pdf_max = amax;
	// for safety
 	pdf_max *= 1.5;	

	cout << "pdf_max " << pdf_max << endl;

	int N_gen = 0;
	int N_tot = 0;
	while(true){
			DalitzEvent evt = generateWeightedEvent();
			double pdfVal = getVal(evt)/evt.getGeneratorPdfRelativeToPhaseSpace(); // /_timePdfMaster->get_marginalPdfs_product(evt) ;
			
			const double height = gRandom->Uniform(0,pdf_max);
			
			///Safety check on the maxmimum generated height
			if( pdfVal > pdf_max ){
				std::cout << "ERROR: PDF above determined maximum." << std::endl;
				std::cout << pdfVal << " > " << pdf_max << std::endl;
				pdf_max = pdf_max * 2.;
			}
			
			///Hit-and-miss
			if( height < pdfVal ) { 
				evt.setValueInVector(7, run);
				evt.setValueInVector(8, trigger);
				eventList.Add(evt);
				N_gen++;
				if (0ul == (N_gen % 500ul)) cout << "Generated event " << N_gen << "/" << N << endl;
			}		
			N_tot ++;
			if(N_gen == N)break;
	}

	cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
 	cout << "Generated " << N_gen << " events ! Efficiecy = " << (double)N_gen/(double)N_tot << endl;

	saveEventListToFile(eventList);
 	return eventList;
   }

    void saveEventListToFile(DalitzEventList& eventList, string name = "toys.root"){

	    NamedParameter<int>  addBkgToToys("addBkgToToys", 0);

	    TFile* out = new TFile(name.c_str(),"RECREATE");
	    TTree* tree = new TTree("DecayTree","DecayTree");
    		
	    double t,dt;
	    int Ds_ID,Ds_finalState;
	    int q_OS,q,q_SS;
	    double eta_OS;
	    double eta_SS;
	    double sw,w;
	    int run,trigger;
	    int year;
		
	    double K[4];
	    double pip[4];
	    double pim[4];
	    double Ds[4];
	    double mB,m_Kpipi,m_Kpi,m_pipi;

    	    TBranch* br_mB = tree->Branch( "Bs_DTF_MM", &mB, "Bs_DTF_MM/D" );
    	    TBranch* br_sw;
	    if(!addBkgToToys) br_sw = tree->Branch( "N_Bs_sw", &sw, "N_Bs_sw/D" );
    	    TBranch* br_w = tree->Branch( "weight", &w, "weight/D" );

    	    TBranch* br_t = tree->Branch( "Bs_BsDTF_TAU", &t, "Bs_BsDTF_TAU/D" );
    	    TBranch* br_dt = tree->Branch( "Bs_BsDTF_TAUERR", &dt, "Bs_BsDTF_TAUERR/D" );

	    TBranch* br_Ds_ID = tree->Branch("Ds_ID",&Ds_ID,"Ds_ID/I");
	    TBranch* br_Ds_finalState = tree->Branch("Ds_finalState",&Ds_finalState,"Ds_finalState/I");

	    TBranch* br_q_OS =tree->Branch("OS_Combination_DEC",&q_OS,"OS_Combination_DEC/I");
	    TBranch* br_eta_OS =  tree->Branch("OS_Combination_PROB",&eta_OS,"OS_Combination_PROB/D");
	    TBranch* br_q_SS = tree->Branch("SS_Kaon_DEC",&q_SS,"SS_Kaon_DEC/I");
	    TBranch* br_eta_SS = tree->Branch("SS_Kaon_PROB",&eta_SS,"SS_Kaon_PROB/D");

	    TBranch* br_run = tree->Branch("run",&run,"run/I");
    	    TBranch* br_year = tree->Branch( "year", &year, "year/I" );
	    TBranch* br_trigger = tree->Branch("TriggerCat",&trigger,"TriggerCat/I");
	    
	    TBranch* br_K0 = tree->Branch("BsDTF_Kplus_PX",&K[0],"BsDTF_Kplus_PX/D");
	    TBranch* br_K1 =tree->Branch("BsDTF_Kplus_PY",&K[1],"BsDTF_Kplus_PY/D");
	    TBranch* br_K2 =tree->Branch("BsDTF_Kplus_PZ",&K[2],"BsDTF_Kplus_PZ/D");
	    TBranch* br_K3 =tree->Branch("BsDTF_Kplus_PE",&K[3],"BsDTF_Kplus_PE/D");
	
            TBranch* br_pip0 =tree->Branch("BsDTF_piplus_PX",&pip[0],"BsDTF_piplus_PX/D");
	    TBranch* br_pip1 =tree->Branch("BsDTF_piplus_PY",&pip[1],"BsDTF_piplus_PY/D");
	    TBranch* br_pip2 =tree->Branch("BsDTF_piplus_PZ",&pip[2],"BsDTF_piplus_PZ/D");
	    TBranch* br_pip3 =tree->Branch("BsDTF_piplus_PE",&pip[3],"BsDTF_piplus_PE/D");
	
	    TBranch* br_pim0 =tree->Branch("BsDTF_piminus_PX",&pim[0],"BsDTF_piminus_PX/D");
	    TBranch* br_pim1 =tree->Branch("BsDTF_piminus_PY",&pim[1],"BsDTF_piminus_PY/D");
	    TBranch* br_pim2 =tree->Branch("BsDTF_piminus_PZ",&pim[2],"BsDTF_piminus_PZ/D");
	    TBranch* br_pim3 =tree->Branch("BsDTF_piminus_PE",&pim[3],"BsDTF_piminus_PE/D");
	
	    TBranch* br_Ds0 =tree->Branch("BsDTF_Ds_PX",&Ds[0],"BsDTF_Ds_PX/D");
	    TBranch* br_Ds1 =tree->Branch("BsDTF_Ds_PY",&Ds[1],"BsDTF_Ds_PY/D");
	    TBranch* br_Ds2 =tree->Branch("BsDTF_Ds_PZ",&Ds[2],"BsDTF_Ds_PZ/D");
	    TBranch* br_Ds3 =tree->Branch("BsDTF_Ds_PE",&Ds[3],"BsDTF_Ds_PE/D");

    	    TBranch* br_m_Kpipi = tree->Branch( "m_Kpipi", &m_Kpipi, "m_Kpipi/D" );
    	    TBranch* br_m_Kpi = tree->Branch( "m_Kpi", &m_Kpi, "m_Kpi/D" );
    	    TBranch* br_m_pipi = tree->Branch( "m_pipi", &m_pipi, "m_pipi/D" );

	    for(int i= 0; i < eventList.size(); i++){

		t = eventList[i].getValueFromVector(0);
		dt = eventList[i].getValueFromVector(1);

		Ds_ID = - eventList[i].getValueFromVector(2);
		Ds_finalState = 0;

		q_OS = eventList[i].getValueFromVector(3);
		eta_OS = eventList[i].getValueFromVector(4);

		q_SS = eventList[i].getValueFromVector(5);
		eta_SS = eventList[i].getValueFromVector(6);

		run = eventList[i].getValueFromVector(7);
		if(run == 2) year = 16;
		else if(run == 3) year = 17;
		else year = 12;
		trigger = eventList[i].getValueFromVector(8);

		mB = eventList[i].getValueFromVector(9);
		sw = eventList[i].getWeight();
		w = eventList[i].getWeight();

		Ds[0] = eventList[i].p(1).Px()/MeV;
		Ds[1] = eventList[i].p(1).Py()/MeV;
		Ds[2] = eventList[i].p(1).Pz()/MeV;
		Ds[3] = eventList[i].p(1).E()/MeV;

		K[0] = eventList[i].p(2).Px()/MeV;
		K[1] = eventList[i].p(2).Py()/MeV;
		K[2] = eventList[i].p(2).Pz()/MeV;
		K[3] = eventList[i].p(2).E()/MeV;

		pip[0] = eventList[i].p(3).Px()/MeV;
		pip[1] = eventList[i].p(3).Py()/MeV;
		pip[2] = eventList[i].p(3).Pz()/MeV;
		pip[3] = eventList[i].p(3).E()/MeV;

		pim[0] = eventList[i].p(4).Px()/MeV;
		pim[1] = eventList[i].p(4).Py()/MeV;
		pim[2] = eventList[i].p(4).Pz()/MeV;
		pim[3] = eventList[i].p(4).E()/MeV;

		m_Kpipi = sqrt(eventList[i].sij(_s234))/GeV;        
		m_Kpi = sqrt(eventList[i].s(2,4))/GeV;        
		m_pipi = sqrt(eventList[i].s(3,4))/GeV;        

		tree->Fill();
	     }

	     tree->Write();
	     out->Write();
	     out->Close();
    }
    
    virtual DalitzHistoSet histoSet(){return _ampsSum->histoSet();}
    
    void doFinalStatsAndSaveForAmp12(MINT::Minimiser* min=0,const std::string& fname = "FitAmpResults", const std::string& fnameROOT="fitFractions"){
        _amps->redoIntegrator();
        _amps_bar->redoIntegrator();
        _amps->doFinalStatsAndSave(min,((string)fname+".tex").c_str(),((string)fnameROOT+".root").c_str());
        _amps_bar->doFinalStatsAndSave(min,((string)fname+"_bar.tex").c_str(),((string)fnameROOT+"_Bar.root").c_str());        
    }
    
    FullAmpsPdfFlexiFastCPV(
		AmpsPdfFlexiFast* amps, AmpsPdfFlexiFast* amps_bar, 
		AmpsPdfFlexiFast* amps_CP, AmpsPdfFlexiFast* amps_bar_CP,
		AmpsPdfFlexiFast* ampsSum, AmpsPdfFlexiFast* ampsSum_CP, 
                const MINT::FitParameter& r,const MINT::FitParameter& delta,const MINT::FitParameter& gamma,
                const MINT::FitParameter& xm,const MINT::FitParameter& ym,const MINT::FitParameter& xp,const MINT::FitParameter& yp,
		const MINT::FitParameter& Gamma, const MINT::FitParameter& dGamma, const MINT::FitParameter& dm
		,const MINT::FitParameter& offset_mean_dt,const MINT::FitParameter& scale_mean_dt,const MINT::FitParameter& scale_mean_2_dt
                ,const MINT::FitParameter& offset_sigma_dt, const MINT::FitParameter& scale_sigma_dt, const MINT::FitParameter& scale_sigma_2_dt
                ,const MINT::FitParameter& offset_sigma2_dt, const MINT::FitParameter& scale_sigma2_dt, const MINT::FitParameter& scale_sigma2_2_dt
                ,const MINT::FitParameter& offset_sigma3_dt, const MINT::FitParameter& scale_sigma3_dt, const MINT::FitParameter& scale_sigma3_2_dt
                ,const MINT::FitParameter& offset_f_dt, const MINT::FitParameter& scale_f_dt, const MINT::FitParameter& scale_f_2_dt
                ,const MINT::FitParameter& offset_f2_dt, const MINT::FitParameter& scale_f2_dt, const MINT::FitParameter& scale_f2_2_dt
                ,const MINT::FitParameter& c0, const MINT::FitParameter& c1, const MINT::FitParameter& c2
                ,const MINT::FitParameter& c3, const MINT::FitParameter& c4, const MINT::FitParameter& c5
                ,const MINT::FitParameter& c6, const MINT::FitParameter& c7, const MINT::FitParameter& c8
                ,const MINT::FitParameter& c9,
                const MINT::FitParameter& p0_os, const MINT::FitParameter& p1_os, const MINT::FitParameter& delta_p0_os, const MINT::FitParameter& delta_p1_os, 
                const MINT::FitParameter& avg_eta_os, const MINT::FitParameter& tageff_os, const MINT::FitParameter& tageff_asym_os, 
                const MINT::FitParameter& p0_ss, const MINT::FitParameter& p1_ss, const MINT::FitParameter& delta_p0_ss, const MINT::FitParameter& delta_p1_ss, 
                const MINT::FitParameter& avg_eta_ss, const MINT::FitParameter& tageff_ss, const MINT::FitParameter& tageff_asym_ss, 
                const MINT::FitParameter& production_asym, const MINT::FitParameter& detection_asym, string marginalPdfsPrefix = ""
                ):
    _amps(amps),_amps_bar(amps_bar),_amps_CP(amps_CP),_amps_bar_CP(amps_bar_CP),_ampsSum(ampsSum),_ampsSum_CP(ampsSum_CP),
    _intA(-1),_intAbar(-1),_intAAbar(-1),_intA_CP(-1),_intAbar_CP(-1),_intAAbar_CP(-1),
    _r(r),_delta(delta),_gamma(gamma),
    _xm(xm),_ym(ym),_xp(xp),_yp(yp),
    _Gamma(Gamma),
    _dGamma(dGamma),
    _dm(dm),
    _offset_mean_dt(offset_mean_dt),
    _scale_mean_dt(scale_mean_dt),
    _scale_mean_2_dt(scale_mean_2_dt),
    _offset_sigma_dt(offset_sigma_dt),    
    _scale_sigma_dt(scale_sigma_dt),
    _scale_sigma_2_dt(scale_sigma_2_dt),
    _offset_sigma2_dt(offset_sigma2_dt),    
    _scale_sigma2_dt(scale_sigma2_dt),
    _scale_sigma2_2_dt(scale_sigma2_2_dt),
    _offset_sigma3_dt(offset_sigma3_dt),    
    _scale_sigma3_dt(scale_sigma3_dt),
    _scale_sigma3_2_dt(scale_sigma3_2_dt),
    _offset_f_dt(offset_f_dt),    
    _scale_f_dt(scale_f_dt),
    _scale_f_2_dt(scale_f_2_dt),
    _offset_f2_dt(offset_f2_dt),    
    _scale_f2_dt(scale_f2_dt),
    _scale_f2_2_dt(scale_f2_2_dt),
    _c0(c0),
    _c1(c1),
    _c2(c2),
    _c3(c3),
    _c4(c4),
    _c5(c5),
    _c6(c6),
    _c7(c7),
    _c8(c8),
    _c9(c9),
    _p0_os(p0_os),
    _p1_os(p1_os),
    _delta_p0_os(delta_p0_os),
    _delta_p1_os(delta_p1_os),
    _avg_eta_os(avg_eta_os),
    _tageff_os(tageff_os),
    _tageff_asym_os(tageff_asym_os),
    _p0_ss(p0_ss),
    _p1_ss(p1_ss),
    _delta_p0_ss(delta_p0_ss),
    _delta_p1_ss(delta_p1_ss),
    _avg_eta_ss(avg_eta_ss),
    _tageff_ss(tageff_ss),
    _tageff_asym_ss(tageff_asym_ss),
    _production_asym(production_asym),
    _detection_asym(detection_asym),
    _marginalPdfsPrefix(marginalPdfsPrefix),
    _useCartCoord("FullAmpsPdfFlexiFastCPV::useCartCoord",1),
    _directCPV("FullAmpsPdfFlexiFastCPV::directCPV",0),
    _min_TAU("min_TAU", 0.4),
    _max_TAU("max_TAU", 10.),
    _min_TAUERR("min_TAUERR", 0.),
    _max_TAUERR("max_TAUERR", 0.1),
    _fasGen(0),
    _sg(0),
    _fasGen_CP(0),
    _sg_CP(0)
    {
		        _timePdfMaster = new TimePdfMaster(_Gamma, _dGamma, _dm
					  ,_offset_mean_dt,_scale_mean_dt,_scale_mean_2_dt
                                          ,_offset_sigma_dt, _scale_sigma_dt, _scale_sigma_2_dt
                                          ,_offset_sigma2_dt, _scale_sigma2_dt, _scale_sigma2_2_dt
                                          ,_offset_sigma3_dt, _scale_sigma3_dt, _scale_sigma3_2_dt
                                          ,_offset_f_dt, _scale_f_dt, _scale_f_2_dt
                                          ,_offset_f2_dt, _scale_f2_dt, _scale_f2_2_dt
                                          ,_c0, _c1, _c2
                                          ,_c3, _c4, _c5
                                          ,_c6, _c7, _c8
                                          ,_c9,
                                          _p0_os, _p1_os, _delta_p0_os, _delta_p1_os, 
                                          _avg_eta_os, _tageff_os, _tageff_asym_os, 
                                          _p0_ss, p1_ss, _delta_p0_ss, _delta_p1_ss, 
                                          _avg_eta_ss, _tageff_ss, _tageff_asym_ss, 
                                          _production_asym, _detection_asym,_marginalPdfsPrefix);

			_timePdfMaster->setAllFitParameters();

        		_s234.push_back(2);
		        _s234.push_back(3);
        		_s234.push_back(4);

			NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
			_pat = DalitzEventPattern(EventPattern.getVector());
			_pat_CP = _pat.makeCPConjugate();
			
			_f_bkg = 0;
    }
};


// Bkg PDF
class FullAmpsPdfFlexiFastBkg : 
	public MINT::PdfBase<IDalitzEvent>, virtual public IDalitzPdf{
    
protected:
    AmpsPdfFlexiFast* _amps;
    AmpsPdfFlexiFast* _amps_CP;

public:
    void parametersChanged(){
        _amps->parametersChanged();
        _amps_CP->parametersChanged();
    }
    void beginFit(){
        _amps->beginFit();
        _amps_CP->beginFit();
    }
    void endFit(){
        _amps->endFit();
        _amps_CP->endFit();
    }
    
    inline double un_normalised_noPs(IDalitzEvent& evt){
        const double f = (double) evt.getValueFromVector(2);
        if(f>0) return _amps->un_normalised_noPs(evt);
	else return _amps_CP->un_normalised_noPs(evt);
    }
    
    virtual double getVal(IDalitzEvent& evt){
        const double f = (double) evt.getValueFromVector(2);
        if(f>0) return _amps->getVal(evt);
	else return _amps_CP->getVal(evt);
    }
    
    virtual double getVal_withPs(IDalitzEvent& evt){return getVal(evt);}
    virtual double getVal_noPs(IDalitzEvent& evt){return getVal(evt);}
    
    virtual double getVal(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal(*evt);
    }
    virtual double getVal_withPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_withPs(*evt);
    }
    virtual double getVal_noPs(IDalitzEvent* evt){
        if(0 == evt) return 0;
        return getVal_noPs(*evt);
    }

    virtual DalitzHistoSet histoSet(){return _amps->histoSet();}

    FullAmpsPdfFlexiFastBkg(AmpsPdfFlexiFast* amps, AmpsPdfFlexiFast* amps_CP ):
    _amps(amps),_amps_CP(amps_CP)
    { ;  }
};


class CPcoeffLL : public Minimisable{
  FullAmpsPdfFlexiFastCPV* _pdf;
  vector<double> coeff;
  double _prec;
  double _C,_Cbar,_S,_Sbar,_D,_Dbar;
  double _r,_k,_delta,_gamma;
public:
  CPcoeffLL(FullAmpsPdfFlexiFastCPV* pdf,double C,double Cbar,double S,double Sbar,double D,double Dbar, double prec) 
	: _pdf(pdf),_C(C),_Cbar(Cbar),_D(D),_Dbar(Dbar),_S(S),_Sbar(Sbar),_prec(prec){
	    _pdf->parametersChanged();//makes sure we are initialised
	    coeff.push_back(_C);
	    coeff.push_back(_Cbar);
	    coeff.push_back(_D);
	    coeff.push_back(_Dbar);
	    coeff.push_back(_S);
	    coeff.push_back(_Sbar);
  }
  double getVal(){
    _pdf->parametersChanged();
    return _pdf->getCPcoeffChi2(coeff,_prec);
  }
};

class FracLL : public Minimisable{
  FlexiFastAmplitudeIntegrator* _integ;
public:
  FracLL(FlexiFastAmplitudeIntegrator* integ) : _integ(integ){
    _integ->getVal();//makes sure we are initialised
  }
  double getVal(){
    return _integ->getFractionChi2();
  }
};

vector<double> getChi2(DalitzEventList& data, DiskResidentEventList& mc){
	
    double minBinWidth = 0.;
    const int dim = 5;
    
    NamedParameter<int> EventPattern("Event Pattern",  531, -431, 321, 211, -211);
    DalitzEventPattern pdg(EventPattern.getVector());
          
    NamedParameter<int> minEventsPerBin("minEventsPerBin", 25);
    HyperPointSet points( dim );
    HyperPoint min(pdg.sijMin(1,3),pdg.sijMin(2,4),pdg.sijMin(3,4),pdg.sijMin(1,2,4),pdg.sijMin(2,3,4));
    HyperPoint max(pdg.sijMax(1,3),pdg.sijMax(2,4),pdg.sijMax(3,4),pdg.sijMax(1,2,4),pdg.sijMax(2,3,4));
    //HyperPoint max(pdg.sijMax(1,3),pdg.sijMax(2,4),pdg.sijMax(3,4),pdg.sijMax(1,2,4),4 * GeV * GeV);

    HyperCuboid limits(min, max );

    vector<int> s124;
    s124.push_back(1);
    s124.push_back(2);
    s124.push_back(4);

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);    	
        
    for (int i = 0; i < data.size(); i++){
        DalitzEvent evt = data[i];
	HyperPoint point( dim );
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.addWeight(evt.getWeight());
      	points.push_back(point);
    }

    HyperHistogram dataHist(limits, points, 
                         /*** Name of the binning algorithm you want to use     */
                         HyperBinningAlgorithms::SMART_MULTI,
                         /***  The minimum number of events allowed in each bin */
                         /***  from the HyperPointSet provided (points1)        */
                         AlgOption::MinBinContent      (minEventsPerBin),                    
                         /*** This minimum bin width allowed. Can also pass a   */
                         /*** HyperPoint if you would like different min bin    */
                         /*** widths for each dimension                         */
                         AlgOption::MinBinWidth        (0.),                                                 
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

//    dataHist.save("histData.root");
//     HyperHistogram binningHist("histData.root",5);    
//     HyperHistogram dataHist( binningHist.getBinning() );
//     dataHist.fill(points); 

    HyperPointSet pointsMC( dim);
    for (int i = 0; i < mc.size(); i++){
     	DalitzEvent evt = mc.getEvent(i);
	HyperPoint point( dim);
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.addWeight(evt.getWeight());
	//if(!(evt.phaseSpace() > 0.))continue;
      	pointsMC.push_back(point);
    }

    HyperHistogram mcHist( dataHist.getBinning() );
    mcHist.fill(pointsMC); 
    //mcHist.save("histMC.root");
    //data.normalise(1);
    mcHist.normalise(dataHist.integral());

    double chi2 = dataHist.chi2(mcHist);
    int nBins   = dataHist.getNBins();

    cout << "chi2 = " << (double)chi2/(nBins-1.) << endl;
    //mcHist.divide(dataHist);
    //mcHist.save("histMC_Data.root");

    vector<double> vals;
    vals.push_back(chi2);
    vals.push_back((double)nBins);

    return vals;
}

vector<double> getChi2LL(DalitzEventList& data, DiskResidentEventList& mc){
	
    double minBinWidth = 0.;
    const int dim = 5;
    
    NamedParameter<int> EventPattern("Event Pattern",  531, -431, 321, 211, -211);
    DalitzEventPattern pdg(EventPattern.getVector());
          
    NamedParameter<int> minEventsPerBin("minEventsPerBin", 25);
    HyperPointSet points( dim );
    HyperPoint min(pdg.sijMin(1,3),pdg.sijMin(2,4),pdg.sijMin(3,4),pdg.sijMin(1,2,4),pdg.sijMin(2,3,4));
    HyperPoint max(pdg.sijMax(1,3),pdg.sijMax(2,4),pdg.sijMax(3,4),pdg.sijMax(1,2,4),pdg.sijMax(2,3,4));
    //HyperPoint max(pdg.sijMax(1,3),pdg.sijMax(2,4),pdg.sijMax(3,4),pdg.sijMax(1,2,4),4 * GeV * GeV);

    HyperCuboid limits(min, max );

    vector<int> s124;
    s124.push_back(1);
    s124.push_back(2);
    s124.push_back(4);

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);    	
        
    for (int i = 0; i < data.size(); i++){
        DalitzEvent evt = data[i];
	HyperPoint point( dim );
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.addWeight(evt.getWeight());
      	points.push_back(point);
    }

    HyperHistogram dataHist(limits, points, 
                         /*** Name of the binning algorithm you want to use     */
                         HyperBinningAlgorithms::LIKELIHOOD,
                         /***  The minimum number of events allowed in each bin */
                         /***  from the HyperPointSet provided (points1)        */
                         AlgOption::MinBinContent      (minEventsPerBin),                    
                         /*** This minimum bin width allowed. Can also pass a   */
                         /*** HyperPoint if you would like different min bin    */
                         /*** widths for each dimension                         */
                         AlgOption::MinBinWidth        (0.),                                                 
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

//    dataHist.save("histData.root");
//     HyperHistogram binningHist("histData.root",5);    
//     HyperHistogram dataHist( binningHist.getBinning() );
//     dataHist.fill(points); 

    HyperPointSet pointsMC( dim);
    for (int i = 0; i < mc.size(); i++){
     	DalitzEvent evt = mc.getEvent(i);
	HyperPoint point( dim);
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.addWeight(evt.getWeight());
	//if(!(evt.phaseSpace() > 0.))continue;
      	pointsMC.push_back(point);
    }

    HyperHistogram mcHist( dataHist.getBinning() );
    mcHist.fill(pointsMC); 
    //mcHist.save("histMC.root");
    //data.normalise(1);
    mcHist.normalise(dataHist.integral());

    double chi2 = dataHist.chi2(mcHist);
    int nBins   = dataHist.getNBins();

    cout << "chi2 = " << (double)chi2/(nBins-1.) << endl;
    //mcHist.divide(dataHist);
    //mcHist.save("histMC_Data.root");

    vector<double> vals;
    vals.push_back(chi2);
    vals.push_back((double)nBins);

    return vals;
}

double getChi2_6D(DalitzEventList& data, DiskResidentEventList& mc){
	
    double minBinWidth = 0.;
    const int dim = 6;
    
    NamedParameter<int> EventPattern("Event Pattern",  531, -431, 321, 211, -211);
    DalitzEventPattern pdg(EventPattern.getVector());
          
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);

    NamedParameter<int> minEventsPerBin("minEventsPerBin", 50);       
    HyperPointSet points( dim );
    HyperPoint min(pdg.sijMin(1,3),pdg.sijMin(2,4),pdg.sijMin(3,4),pdg.sijMin(1,2,4),pdg.sijMin(2,3,4),(double)min_TAU);
    HyperPoint max(pdg.sijMax(1,3),pdg.sijMax(2,4),pdg.sijMax(3,4),pdg.sijMax(1,2,4),pdg.sijMax(2,3,4),(double)max_TAU);
    //HyperPoint max(pdg.sijMax(1,3),pdg.sijMax(2,4),pdg.sijMax(3,4),pdg.sijMax(1,2,4),4 * GeV * GeV);

    HyperCuboid limits(min, max );

    vector<int> s124;
    s124.push_back(1);
    s124.push_back(2);
    s124.push_back(4);

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);    	
        
    for (int i = 0; i < data.size(); i++){
        DalitzEvent evt = data[i];
	HyperPoint point( dim );
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.at(5)= evt.getValueFromVector(0);
      	point.addWeight(evt.getWeight());
      	points.push_back(point);
    }

    HyperHistogram dataHist(limits, points, 
                         /*** Name of the binning algorithm you want to use     */
                         HyperBinningAlgorithms::SMART,  
                         /***  The minimum number of events allowed in each bin */
                         /***  from the HyperPointSet provided (points1)        */
                         AlgOption::MinBinContent      (minEventsPerBin),                    
                         /*** This minimum bin width allowed. Can also pass a   */
                         /*** HyperPoint if you would like different min bin    */
                         /*** widths for each dimension                         */
                         AlgOption::MinBinWidth        (0.),                                                 
                         /*** If you want to use the sum of weights rather than */
                         /*** the number of events, set this to true.           */    
                         AlgOption::UseWeights         (true),                         
                         /*** Some algorithms use a random number generator. Set*/
                         /*** the seed here                                     */
                         AlgOption::RandomSeed         (1),                         
                         /*** What dimesnion would you like to split first? Only*/
                         /*** applies to certain algortihms                     */
                         AlgOption::StartDimension     (4)                        
                         /*** What dimesnions would you like to bin in?         */
                         //AlgOption::BinningDimensions  (binningDims),                      
                         /*** Setting this option will make the agorithm draw   */
                         /*** the binning scheme at each iteration              */
                         //AlgOption::DrawAlgorithm("Algorithm")                 
                         );

//    dataHist.save("histData.root");
//     HyperHistogram binningHist("histData.root",5);    
//     HyperHistogram dataHist( binningHist.getBinning() );
//     dataHist.fill(points); 

    HyperPointSet pointsMC( dim);
    for (int i = 0; i < mc.size(); i++){
     	DalitzEvent evt = mc.getEvent(i);
	HyperPoint point( dim);
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.at(5)= evt.getValueFromVector(0);
	point.addWeight(evt.getWeight());
	//if(!(evt.phaseSpace() > 0.))continue;
      	pointsMC.push_back(point);
    }

    HyperHistogram mcHist( dataHist.getBinning() );
    mcHist.fill(pointsMC); 
    //mcHist.save("histMC.root");
    //data.normalise(1);
    mcHist.normalise(dataHist.integral());

    double chi2 = dataHist.chi2(mcHist);
    int nBins   = dataHist.getNBins();

    cout << "chi2 = " << (double)chi2/(nBins-1.) << endl;
    //mcHist.divide(dataHist);
    //mcHist.save("histMC_Data.root");

    return (double)chi2/(nBins-1.);
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


void ampFit(int step=0, string mode = "fit"){

    /// Generate list of amplitudes
    FitAmplitude::AutogenerateFitFile();

    /// Options
    NamedParameter<int> updateAnaNote("updateAnaNote", 0);
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    int seed = RandomSeed + step;
    ranLux.SetSeed((int)seed);
    gRandom = &ranLux;
    RooRandom::randomGenerator()->SetSeed(seed);
    TString prefix = "";
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    DalitzEventPattern pat_CP = pat.makeCPConjugate();

    NamedParameter<string> InputFileName("InputFileName", (std::string) "/auto/data/dargent/BsDsKpipi/Final/signal_tagged.root", (char*) 0);
    NamedParameter<string> channel("channel", (std::string) "norm", (char*) 0);
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);

    NamedParameter<string> weightName("weightName", (std::string) "N_Bs_sw", (char*) 0);
    NamedParameter<double> min_TAU("min_TAU", 0.4);
    NamedParameter<double> max_TAU("max_TAU", 10.);
    NamedParameter<double> max_TAU_ForMixingPlot("max_TAU_ForMixingPlot", 4.);
    NamedParameter<double> min_TAUERR("min_TAUERR", 0.);
    NamedParameter<double> max_TAUERR("max_TAUERR", 0.1);
    NamedParameter<double> w_max("w_max", 0.5);

    NamedParameter<double> min_year("min_year", 10);
    NamedParameter<double> max_year("max_year", 20);
    NamedParameter<double> min_trigger("min_trigger", -10);
    NamedParameter<double> max_trigger("max_trigger", 10);

    NamedParameter<int>  doPlots("doPlots", 1);
    NamedParameter<int>  nBins("nBins", 50);
    NamedParameter<int>  nBinst("nBinst", 50);
    NamedParameter<int>  nBinsAsym("nBinsAsym", 8);
    NamedParameter<int>  N_plot_it("N_plot_it",1);

    NamedParameter<int>  randomizeStartVals("randomizeStartVals", 0);
    NamedParameter<int>  doSimFit("doSimFit", 0);

    NamedParameter<int>  do2DScan("do2DScan", 0);
    NamedParameter<int>  doScanParam("doScanParam", 0);
    NamedParameter<string> scanParam("scanParam", (std::string) "", (char*) 0);
    NamedParameter<double>  scanParam_min("scanParam_min", 0.);
    NamedParameter<double>  scanParam_max("scanParam_max", 1.);
    NamedParameter<double>  scanParam_steps("scanParam_steps", 100.);

    NamedParameter<int>  useLASSO("useLASSO", 0);
    NamedParameter<double>  lambda("lambda", 1.);

    NamedParameter<int>  initCPcoeff("initCPcoeff", 0);
    NamedParameter<int>  doBootstrap("doBootstrap", 0);
    NamedParameter<int>  N_bootstrap("N_bootstrap", 10000);
    NamedParameter<int> doFractions("doFractions",1);
    NamedParameter<int> doFractionsErr("doFractionsErr",0);

    NamedParameter<int>  doToyStudy("doToyStudy", 0);
    NamedParameter<double> N_scale_toys("N_scale_toys", 1);
    NamedParameter<int>  addBkgToToys("addBkgToToys", 0);

    NamedParameter<int>  useGaussConstrainsTagging("useGaussConstrainsTagging", 0);
    NamedParameter<int>  usePerEventDetAsym("usePerEventDetAsym", 0);
    NamedParameter<int>  useGaussConstrainsBias("useGaussConstrainsBias", 0);

    NamedParameter<int>  doAccSystematics("doAccSystematics", 0);
    NamedParameter<int>  useCholDec("useCholDec", 0);
    NamedParameter<int>  varPerParChol("varPerParChol", 100);
    int chol_index = (step-1)/varPerParChol ;

    NamedParameter<int>  doRadiusBWSystematic("doRadiusBWSystematic", 0);
    NamedParameter<string> doSystematic("doSystematic", (std::string) "", (char*) 0);
    NamedParameter<string> doLineshapeSystematic("doLineshapeSystematic", (std::string) "", (char*) 0);
    NamedParameter<string> addAmpName("addAmpName", (std::string) "", (char*) 0);

    NamedParameter<double> integPrecision("IntegPrecision", 1.e-2);
    NamedParameter<std::string> integMethod("IntegMethod", (std::string)"efficient");
    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    NamedParameter<string> IntegratorEventFileCP("IntegratorEventFileCP", (std::string) "", (char*) 0);
    NamedParameter<string> IntegratorEventFileBkg("IntegratorEventFileBkg", (std::string) "", (char*) 0);
    TString integratorEventFile = (string) IntegratorEventFile;
    TString integratorEventFile_CP = (string) IntegratorEventFileCP;
    if(integratorEventFile_CP == ""){
        integratorEventFile_CP = integratorEventFile; 
        integratorEventFile_CP.ReplaceAll(".root","_CP.root");
    }
    TString integratorEventFileBkg = (string) IntegratorEventFileBkg;
    if(integratorEventFileBkg == ""){
        integratorEventFileBkg = integratorEventFile; 
    }
    TString integratorEventFileBkg_CP = integratorEventFileBkg;
    integratorEventFileBkg_CP.ReplaceAll(".root","_CP.root");
    
    /// Define PDF 

    /// Define amplitude model
    DalitzEventList eventListPhsp,eventListPhsp_CP;
    eventListPhsp.generatePhaseSpaceEvents(2,pat);
    eventListPhsp_CP.generatePhaseSpaceEvents(2,pat_CP);

    FitAmpSum fas((DalitzEventPattern)pat);
    FitAmpSum fas_bar((DalitzEventPattern)pat,"Bar_","");
    fas.getVal(eventListPhsp[0]);
    fas_bar.getVal(eventListPhsp[0]);
    fas.print();
    fas_bar.print();
    /// Define B -> f amplitude        
    fas.setTag(1);
    /// Define Bbar -> f amplitude
    fas_bar.setTag(-1);
    
    FitAmpIncoherentSum fasBkg(pat,"Bkg_","Bkg_"); 
    fasBkg.getVal(eventListPhsp[0]);
    fasBkg.print();

    /// Lineshape systematics
    MinuitParameterSet* mps = MinuitParameterSet::getDefaultSet();
    Neg2LLMultiConstraint constrains_lineshapeSys(MinuitParameterSet::getDefaultSet(),("_" + (string)doLineshapeSystematic).c_str());
    if((string)doLineshapeSystematic != "" && mode == "fit")constrains_lineshapeSys.smearInputValues();

    if(doRadiusBWSystematic && mode == "fit"){
				double val = gRandom->Uniform(0.,0.003);
				mps->getParPtr("BW_radius")->setCurrentFitVal(val);
				((FitParameter*)mps->getParPtr("BW_radius"))->setInit(val);
				cout << "Set parameter " << "BW_radius" << " to " << val << endl;  
    }
    /// Normalize amps
    {
    	DalitzEventList eventListNorm;
        TFile *file =  TFile::Open(((string)IntegratorEventFile).c_str());
        TTree* tree=dynamic_cast<TTree*>(file->Get("DalitzEventList"));
        eventListNorm.fromNtuple(tree,0.5);
        fas.normalizeAmps(eventListNorm);
        fas_bar.normalizeAmps(eventListNorm);
        // fasBkg.normalizeAmps(eventListNorm);
        file->Close();
    }
    
    /// CP conjugate amplitudes
    FitAmpSum fas_CP(fas);
    fas_CP.CPConjugateSameFitParameters();

    FitAmpSum fas_bar_CP(fas_bar);
    fas_bar_CP.CPConjugateSameFitParameters();

    FitAmpIncoherentSum fasBkg_CP(fasBkg);
    fasBkg_CP.CPConjugateSameFitParameters();

    /// Add amplitudes: A + r e^(i gamma) Abar
    counted_ptr<FitAmpList> sumList = fas.GetCloneSameFitParameters();
    FitAmpSum fas_sum(*sumList);
    fas_sum.addAsList(fas_bar,1.);
    fas_sum.getVal(eventListPhsp[0]);
    
    AmpsPdfFlexiFast ampsSig(pat, &fas, 0, integPrecision,integMethod, (std::string) integratorEventFile);
    AmpsPdfFlexiFast ampsSig_bar(pat, &fas_bar, 0, integPrecision,integMethod, (std::string) integratorEventFile);
    AmpsPdfFlexiFast ampsSum(pat, &fas_sum, 0, integPrecision,integMethod, (std::string) integratorEventFile);

    counted_ptr<FitAmpList> sumList_CP = fas_CP.GetCloneSameFitParameters();
    FitAmpSum fas_sum_CP(*sumList_CP);
    fas_sum_CP.addAsList(fas_bar_CP,1.);
    fas_sum_CP.getVal(eventListPhsp_CP[0]);
    
    AmpsPdfFlexiFast ampsSig_CP(pat_CP, &fas_CP, 0, integPrecision,integMethod, (std::string) integratorEventFile_CP);
    AmpsPdfFlexiFast ampsSig_bar_CP(pat_CP, &fas_bar_CP, 0, integPrecision,integMethod, (std::string) integratorEventFile_CP);
    AmpsPdfFlexiFast ampsSum_CP(pat_CP, &fas_sum_CP, 0, integPrecision,integMethod, (std::string) integratorEventFile_CP);

    /// Bkg
    FitParameter sigfraction("SigFraction",1,0.999999,0.01);
    AmpsPdfFlexiFast ampsBkg(pat, &fasBkg, 0, integPrecision*10,integMethod, (std::string) integratorEventFileBkg);
    AmpsPdfFlexiFast ampsBkg_CP(pat_CP, &fasBkg_CP, 0, integPrecision*10,integMethod, (std::string) integratorEventFileBkg_CP);
    FullAmpsPdfFlexiFastBkg pdf_bkg(&ampsBkg,&ampsBkg_CP);

    /// Fit parameters
    /// Randomize start vals
    if(randomizeStartVals){
        for(int i=0; i < mps->size(); i++){
            if(((FitParameter*)mps->getParPtr(i))->iFixInit()!=0)continue;
            if(A_is_in_B("_Amp",((FitParameter*)mps->getParPtr(i))->name()) || A_is_in_B("_Phase",((FitParameter*)mps->getParPtr(i))->name()) ){
                double val = A_is_in_B("_Amp",((FitParameter*)mps->getParPtr(i))->name() ) ? gRandom->Uniform(0.5,1.5) * mps->getParPtr(i)->mean() : gRandom->Uniform(-180,180);
                mps->getParPtr(i)->setCurrentFitVal(val);
                ((FitParameter*)mps->getParPtr(i))->setInit(val);    
                cout << "Setting " << ((FitParameter*)mps->getParPtr(i))->name() << "  to  " << val << endl; 
            }
        }
    } 

    FitParameter  r("r",1,0.,0.1);
    FitParameter  delta("delta",1,100.,1.);
    FitParameter  gamma("gamma",1,70,1.);

    FitParameter xm("xm",1,0,0.01);
    FitParameter ym("ym",1,0,0.01); 
    FitParameter xp("xp",1,0,0.01);
    FitParameter yp("yp",1,0,0.01); 

    if(doScanParam){
        for(int i=0; i < mps->size(); i++){
            if(A_is_in_B(scanParam,((FitParameter*)mps->getParPtr(i))->name())){
                double val = scanParam_min + step * (scanParam_max - scanParam_min)/scanParam_steps; 
                ((FitParameter*)mps->getParPtr(i))->setInit(val);
                ((FitParameter*)mps->getParPtr(i))->fix();		
                cout << "Setting " << ((FitParameter*)mps->getParPtr(i))->name() << "  to  " << val << endl; 
            }
        }
    } 

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
    FitParameter  p1_os("p1_os",1,0.,0.);
    FitParameter  delta_p0_os("delta_p0_os",1,0.,0.);
    FitParameter  delta_p1_os("delta_p1_os",1,0.,0.);
    FitParameter  avg_eta_os("avg_eta_os",1,0.,0.);
    FitParameter  tageff_os("tageff_os",1,0.,0.);
    FitParameter  tageff_asym_os("tageff_asym_os",1,0.,0.);
    FitParameter  p0_ss("p0_ss",1,0.,0.);
    FitParameter  p1_ss("p1_ss",1,0.,0.);
    FitParameter  delta_p0_ss("delta_p0_ss",1,0.,0.);
    FitParameter  delta_p1_ss("delta_p1_ss",1,0.,0.);
    FitParameter  avg_eta_ss("avg_eta_ss",1,0.,0.);
    FitParameter  tageff_ss("tageff_ss",1,0.,0.);
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

    FitParameter  offset_mean_dt_Run2_17("offset_mean_dt_Run2_17",1,0,0.1);
    FitParameter  scale_mean_dt_Run2_17("scale_mean_dt_Run2_17",1,0,0.1);
    FitParameter  scale_mean_2_dt_Run2_17("scale_mean_2_dt_Run2_17",1,0,0.1);
    FitParameter  offset_sigma_dt_Run2_17("offset_sigma_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_sigma_dt_Run2_17("scale_sigma_dt_Run2_17",1,1.,0.1);
    FitParameter  scale_sigma_2_dt_Run2_17("scale_sigma_2_dt_Run2_17",1,0.,0.1);
    FitParameter  offset_sigma2_dt_Run2_17("offset_sigma2_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_sigma2_dt_Run2_17("scale_sigma2_dt_Run2_17",1,1.,0.1);
    FitParameter  scale_sigma2_2_dt_Run2_17("scale_sigma2_2_dt_Run2_17",1,0.,0.1);
    FitParameter  offset_sigma3_dt_Run2_17("offset_sigma3_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_sigma3_dt_Run2_17("scale_sigma3_dt_Run2_17",1,1.,0.1);
    FitParameter  scale_sigma3_2_dt_Run2_17("scale_sigma3_2_dt_Run2_17",1,0.,0.1);
    FitParameter  offset_f_dt_Run2_17("offset_f_dt_Run2_17",1,1,0.1);
    FitParameter  scale_f_dt_Run2_17("scale_f_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_f_2_dt_Run2_17("scale_f_2_dt_Run2_17",1,0.,0.1);
    FitParameter  offset_f2_dt_Run2_17("offset_f2_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_f2_dt_Run2_17("scale_f2_dt_Run2_17",1,0.,0.1);
    FitParameter  scale_f2_2_dt_Run2_17("scale_f2_2_dt_Run2_17",1,0.,0.1);

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

    /// Make full time-dependent PDF
    string marginalPdfsPrefix = "comb";
    FullAmpsPdfFlexiFastCPV pdf(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP
		      ,r,delta,gamma
		      , xm, ym, xp, yp
		      ,Gamma, dGamma, dm
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

    // Simultaneous pdfs
    FullAmpsPdfFlexiFastCPV pdf_Run1(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP
		      ,r,delta,gamma, xm, ym, xp, yp,
		      Gamma, dGamma, dm,
		      offset_mean_dt_Run1,scale_mean_dt_Run1,scale_mean_2_dt_Run1
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

    FullAmpsPdfFlexiFastCPV pdf_Run2(&ampsSig,&ampsSig_bar,&ampsSig_CP,&ampsSig_bar_CP,&ampsSum,&ampsSum_CP
		      ,r,delta,gamma, xm, ym, xp, yp,
		      Gamma, dGamma, dm,
		      offset_mean_dt_Run2,scale_mean_dt_Run2,scale_mean_2_dt_Run2
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
    
    tree_norm->SetBranchAddress("B_DTF_MM",&mB);
    tree_norm->SetBranchAddress("B_FullDTF_TAU",&t);
    tree_norm->SetBranchAddress("B_FullDTF_TAUERR",&dt);
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
        
        if(t < min_TAU || t > max_TAU )continue;
        if( dt < min_TAUERR || dt > max_TAUERR )continue;
        if(year < min_year || year > max_year) continue;
        if(trigger < min_trigger || trigger > max_trigger) continue;
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

    /// Total model
    SumPdf<IDalitzEvent> tot_pdf_Run1(sigfraction,pdf_Run1,pdf_bkg);
    SumPdf<IDalitzEvent> tot_pdf_Run2(sigfraction,pdf_Run2,pdf_bkg);

    /// Likelihood
    Neg2LL neg2LL(pdf, eventList);   
    Neg2LL neg2LL_Run1(tot_pdf_Run1, eventList_Run1);    
    Neg2LL neg2LL_Run2(tot_pdf_Run2, eventList_Run2);    
 
    Neg2LLSum neg2LL_sim;
    if(eventList_Run1.size()>0)neg2LL_sim.add(&neg2LL_Run1);
    if(eventList_Run2.size()>0)neg2LL_sim.add(&neg2LL_Run2);

    /// Add tagging constraints
    Neg2LLMultiConstraint constrains_tagging_Run1(MinuitParameterSet::getDefaultSet(),"_Tagging_Run1");
    Neg2LLMultiConstraint constrains_tagging_Run2(MinuitParameterSet::getDefaultSet(),"_Tagging_Run2");
    if(useGaussConstrainsTagging){
 	    neg2LL_sim.add(&constrains_tagging_Run1);
 	    neg2LL_sim.add(&constrains_tagging_Run2);
    }

    /// Systematic studies
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

    /// LASSO
    double stepSize = 1;
    lambda = lambda + (step-1) * stepSize;
    LASSO_flexi lasso(&ampsSig,lambda);
    LASSO_flexi lasso_bar(&ampsSig_bar,lambda);
    
    Neg2LLSum neg2LL_lasso(&neg2LL,&lasso_bar);
    Neg2LLSum neg2LL_sim_lasso(&neg2LL_sim,&lasso_bar);    
    if(useLASSO == 1) { 
        neg2LL_lasso.add(&lasso);    
        neg2LL_sim_lasso.add(&lasso);
    }

    /// Fit
    Minimiser mini;
    if(useLASSO>0){
    	if(doSimFit)mini.attachFunction(&neg2LL_sim_lasso);
    	else mini.attachFunction(&neg2LL_lasso);    
    }
    else {
    	if(doSimFit)mini.attachFunction(&neg2LL_sim);
    	else mini.attachFunction(&neg2LL);    
    }
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
    
    TH1D* m_DKs = new TH1D("m_DKs","; m(DK_{s}) [GeV]; Yield",nBins,2,5.5);
    TH1D* m_Dpi = new TH1D("m_Dpi","; m(D#pi) [GeV]; Yield",nBins,1,5.5);
    TH1D* m_Kspi = new TH1D("m_Kspi","; m(K_{s}#pi) [GeV]; Yield",nBins,0,4);
    
    TH2D* m_Kspi_DKs = new TH2D("m_Kspi_DKs","; m(K_{s}#pi) [GeV]; m(DK_{s}) [GeV];  Yield",nBins,0,4,nBins,2,5.5);
    TH2D* m_Kspi_Dpi = new TH2D("m_Kspi_Dpi","; m(K_{s}#pi) [GeV]; m(D#pi) [GeV];  Yield",nBins,0,4,nBins,2,5.5);
    TH2D* m_Dpi_DKs = new TH2D("m_Dpi_DKs","; m(D#pi) [GeV]; m(DK_{s}) [GeV];  Yield",nBins,2,5.5,nBins,2,5.5);
    
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
        
        m_DKs->Fill(sqrt(eventList[i].s(1,2)/(GeV*GeV)),eventList[i].getWeight());
        m_Dpi->Fill(sqrt(eventList[i].s(1,3)/(GeV*GeV)),eventList[i].getWeight());
        m_Kspi->Fill(sqrt(eventList[i].s(2,3)/(GeV*GeV)),eventList[i].getWeight());
    
        N += eventList[i].getWeight();
        if(!doSimFit)sigma_t_eff += pdf.getCalibratedResolution(eventList[i].getValueFromVector(1)) * eventList[i].getWeight();
        
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
                calibrated_mistag_os = pdf_Run1.getCalibratedMistag_OS(eventList[i]);
                calibrated_mistag_ss = pdf_Run1.getCalibratedMistag_SS(eventList[i]);
            }
            else{
                calibrated_mistag_os = pdf_Run2.getCalibratedMistag_OS(eventList[i]);
                calibrated_mistag_ss = pdf_Run2.getCalibratedMistag_SS(eventList[i]);                
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
                D_res = exp(-pow(pdf_Run1.getCalibratedResolution(eventList[i].getValueFromVector(1))*dm,2)/2.);
            }
            else if(run_evt==2){
                D_res = exp(-pow(pdf_Run2.getCalibratedResolution(eventList[i].getValueFromVector(1))*dm,2)/2.);
            }
        }
        
        double D_tag = 0.;
        if(q_eff != 0) D_tag = (1.-2.*abs(w_eff));
        double D_tot = D_tag * D_res;
        
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
            toys.Add(pdf_Run1.generateToys(N_scale_toys * N_Run1,1,0));
            toys.Add(pdf_Run2.generateToys(N_scale_toys *N_Run2,2,0));
        }
        else toys.Add(pdf.generateToys(N_scale_toys *N));

        pdf.saveEventListToFile(toys,((string)OutputDir+"toys_"+anythingToString((int)step)+".root").c_str());
        return;
    }

    string outTableName = (string)OutputDir+"FitAmpResults";//_"+anythingToString((int)step);
    string fractionFileName = (string)OutputDir+"fitFractions";//_"+ (string)doSystematic+ "_" + anythingToString((int)step);

    double C_val,D_val,S_val;
    double Cbar_val,Dbar_val,Sbar_val;
    double k_val;
    vector<double> CP_coeff;    
    if(mode != "gen"){
        CP_coeff = pdf_Run1.calculateCP_coeff();
        C_val = CP_coeff[0];
        Cbar_val = CP_coeff[1];
        D_val = CP_coeff[2];
        Dbar_val = CP_coeff[3];
        S_val = CP_coeff[4];
        Sbar_val = CP_coeff[5];
        k_val = CP_coeff[6];
        if(mode == "fit")pdf.doFinalStatsAndSaveForAmp12(&mini,outTableName,fractionFileName);
    }

    /// Loop over MC
    if(doPlots){

	/// Fit histograms
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
        
        TH1D* m_DKs_fit = (TH1D*) m_DKs->Clone("m_DKs_fit");
        TH1D* m_Dpi_fit = (TH1D*) m_Dpi->Clone("m_Dpi_fit");
        TH1D* m_Kspi_fit = (TH1D*) m_Kspi->Clone("m_Kspi_fit");
        
        TH1D* m_DKs_fit_1 = (TH1D*) m_DKs->Clone("m_DKs_fit_1");
        TH1D* m_Dpi_fit_1 = (TH1D*) m_Dpi->Clone("m_Dpi_fit_1");
        TH1D* m_Kspi_fit_1 = (TH1D*) m_Kspi->Clone("m_Kspi_fit_1");
        
        TH1D* m_DKs_fit_2 = (TH1D*) m_DKs->Clone("m_DKs_fit_2");
        TH1D* m_Dpi_fit_2 = (TH1D*) m_Dpi->Clone("m_Dpi_fit_2");
        TH1D* m_Kspi_fit_2 = (TH1D*) m_Kspi->Clone("m_Kspi_fit_2");
        
        TH1D* m_DKs_fit_3 = (TH1D*) m_DKs->Clone("m_DKs_fit_3");
        TH1D* m_Dpi_fit_3 = (TH1D*) m_Dpi->Clone("m_Dpi_fit_3");
        TH1D* m_Kspi_fit_3 = (TH1D*) m_Kspi->Clone("m_Kspi_fit_3");
        
        TH1D* m_DKs_fit_4 = (TH1D*) m_DKs->Clone("m_DKs_fit_4");
        TH1D* m_Dpi_fit_4 = (TH1D*) m_Dpi->Clone("m_Dpi_fit_4");
        TH1D* m_Kspi_fit_4 = (TH1D*) m_Kspi->Clone("m_Kspi_fit_4");
        
        TH1D* m_DKs_fit_5 = (TH1D*) m_DKs->Clone("m_DKs_fit_5");
        TH1D* m_Dpi_fit_5 = (TH1D*) m_Dpi->Clone("m_Dpi_fit_5");
        TH1D* m_Kspi_fit_5 = (TH1D*) m_Kspi->Clone("m_Kspi_fit_5");
        
        TH1D* m_DKs_fit_A = (TH1D*) m_DKs->Clone("m_DKs_fit_A");
        TH1D* m_Dpi_fit_A = (TH1D*) m_Dpi->Clone("m_Dpi_fit_A");
        TH1D* m_Kspi_fit_A = (TH1D*) m_Kspi->Clone("m_Kspi_fit_A");
        
        TH1D* m_DKs_fit_Abar = (TH1D*) m_DKs->Clone("m_DKs_fit_Abar");
        TH1D* m_Dpi_fit_Abar = (TH1D*) m_Dpi->Clone("m_Dpi_fit_Abar");
        TH1D* m_Kspi_fit_Abar = (TH1D*) m_Kspi->Clone("m_Kspi_fit_Abar");
	
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

	DiskResidentEventList eventListMC(((string) IntegratorEventFile).c_str(),"OPEN");
	DiskResidentEventList eventListMC_rw(pat,("dummy_"+anythingToString(step)+".root").c_str(),"RECREATE");
	
        vector<string> ampNames1;
        ampNames1.push_back("K*(892)+");
        
        vector<string> ampNames2;
        ampNames2.push_back("K(0)*(1430)+");
        ampNames2.push_back("K(0)*(1430)+");
        ampNames2.push_back("K(2)*(1430)+");
        ampNames2.push_back("K*(1680)+");
        
        vector<string> ampNames3;
        ampNames3.push_back("D(0)*(2300)0");
        ampNames3.push_back("D(2)*(2460)0");
        ampNames3.push_back("D*(2600)0");
        ampNames3.push_back("D(3)*(2750)0");
        ampNames3.push_back("D(3000)0");
        ampNames3.push_back("X_S0");
        ampNames3.push_back("X_V0");
        
        vector<string> ampNames4;
        ampNames4.push_back("D(s2)(2573)");
        ampNames4.push_back("D(s1)(2700)");
        ampNames4.push_back("D(s1)*(2860)");
        ampNames4.push_back("D(s3)*(2860)");
        ampNames4.push_back("D(s2)(3040)-");
        ampNames4.push_back("X2_S");
        ampNames4.push_back("X2_V");
        
        vector<string> ampNames5;
        ampNames5.push_back("NonRes");

	///Dalitz plots 
	for(int i = 0; i < eventListMC.size(); i++){
	
			DalitzEvent evt(eventListMC.getEvent(i));
            evt.setWeight(evt.getWeight()*1000000);

			int f = gRandom->Uniform()>0.5 ? 1 : -1;
			if(f==-1)evt.CP_conjugateYourself();
			evt.setValueInVector(2, f);

			double pdfVal = 0;
			if(doSimFit) {
				pdfVal += pdf_Run1.getVal_timeIntegrated(evt) * N_Run1/N;
				pdfVal += pdf_Run2.getVal_timeIntegrated(evt) * N_Run2/N;
			}
			else pdfVal = pdf.getVal_timeIntegrated(evt);
			double weight = pdfVal*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();			
			
			double weight_A,weight_Abar;

			double weight1= pdf_Run2.getAmpVal_timeIntegrated(evt, fas, fas_bar,fas_CP, fas_bar_CP, ampNames1)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
			double weight2= pdf_Run2.getAmpVal_timeIntegrated(evt, fas, fas_bar,fas_CP, fas_bar_CP, ampNames2)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
			double weight3= pdf_Run2.getAmpVal_timeIntegrated(evt, fas, fas_bar,fas_CP, fas_bar_CP, ampNames3)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
			double weight4= pdf_Run2.getAmpVal_timeIntegrated(evt, fas, fas_bar,fas_CP, fas_bar_CP, ampNames4)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
			double weight5= pdf_Run2.getAmpVal_timeIntegrated(evt, fas, fas_bar,fas_CP, fas_bar_CP, ampNames5)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();

			if(f==1){ 
				weight_A = fas.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
				weight_Abar = fas_bar.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
			}
			else {		
				weight_A = fas_CP.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
				weight_Abar = fas_bar_CP.RealVal(evt)*evt.getWeight()/evt.getGeneratorPdfRelativeToPhaseSpace();
			}
        
            if(i < 20){
                cout << "evt " << i << endl;
                cout <<  pdfVal << endl;
                cout <<  weight << endl;
                cout <<   evt.getWeight() << endl;
                cout <<   evt.getGeneratorPdfRelativeToPhaseSpace() << endl;
                cout <<   weight1<< endl;
                cout <<   weight2<< endl;
                cout <<   weight3<< endl;
                cout <<   weight4<< endl;
                cout <<   weight5<< endl;
                cout <<   weight_A<< endl;
                cout <<   weight_Abar<< endl;
            }
        
            m_DKs_fit->Fill(sqrt(evt.s(1,2)/(GeV*GeV)),weight);
            m_Dpi_fit->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight);
            m_Kspi_fit->Fill(sqrt(evt.s(2,3)/(GeV*GeV)),weight);
            
            m_DKs_fit_1->Fill(sqrt(evt.s(1,2)/(GeV*GeV)),weight1);
            m_Dpi_fit_1->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight1);
            m_Kspi_fit_1->Fill(sqrt(evt.s(2,3)/(GeV*GeV)),weight1);
            
            m_DKs_fit_2->Fill(sqrt(evt.s(1,2)/(GeV*GeV)),weight2);
            m_Dpi_fit_2->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight2);
            m_Kspi_fit_2->Fill(sqrt(evt.s(2,3)/(GeV*GeV)),weight2);
            
            m_DKs_fit_3->Fill(sqrt(evt.s(1,2)/(GeV*GeV)),weight3);
            m_Dpi_fit_3->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight3);
            m_Kspi_fit_3->Fill(sqrt(evt.s(2,3)/(GeV*GeV)),weight3);
            
            m_DKs_fit_4->Fill(sqrt(evt.s(1,2)/(GeV*GeV)),weight4);
            m_Dpi_fit_4->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight4);
            m_Kspi_fit_4->Fill(sqrt(evt.s(2,3)/(GeV*GeV)),weight4);
            
            m_DKs_fit_5->Fill(sqrt(evt.s(1,2)/(GeV*GeV)),weight5);
            m_Dpi_fit_5->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight5);
            m_Kspi_fit_5->Fill(sqrt(evt.s(2,3)/(GeV*GeV)),weight5);
        
            m_DKs_fit_A->Fill(sqrt(evt.s(1,2)/(GeV*GeV)),weight_A);
            m_Dpi_fit_A->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight_A);
            m_Kspi_fit_A->Fill(sqrt(evt.s(2,3)/(GeV*GeV)),weight_A);
            
            m_DKs_fit_Abar->Fill(sqrt(evt.s(1,2)/(GeV*GeV)),weight_Abar);
            m_Dpi_fit_Abar->Fill(sqrt(evt.s(1,3)/(GeV*GeV)),weight_Abar);
            m_Kspi_fit_Abar->Fill(sqrt(evt.s(2,3)/(GeV*GeV)),weight_Abar);
			
			evt.setWeight(weight);
			eventListMC_rw.Add(evt);
	}

	/// Time plots	
	FitParameter  C("C",1,CP_coeff[0],0.1);
	FitParameter  D("D",1,CP_coeff[2],0.1);
	FitParameter  D_bar("D_bar",1,CP_coeff[3],0.1);
	FitParameter  S("S",1,CP_coeff[4],0.1);
	FitParameter  S_bar("S_bar",1,CP_coeff[5],0.1);
	FitParameter  k("k",1,1.,0.1);
	
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
    
	for(int n = 0; n < N_plot_it; n++){   /// Multiple iterations needed to release memory 
		int N_sample = 250000;
		DalitzEventList sampleEvents;
   		if(doSimFit) {
			sampleEvents.Add(t_pdf_Run1.generateToys(N_sample * N_Run1/N,1,0));
			sampleEvents.Add(t_pdf_Run2.generateToys(N_sample *N_Run2/N,2,0));
		}
		else sampleEvents.Add(t_pdf.generateToys(N_sample));	

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
					D_res = exp(-pow(pdf_Run1.getCalibratedResolution(evt.getValueFromVector(1))*dm,2)/2.);
				}
				else if(run_MC==2){
					D_res = exp(-pow(pdf_Run2.getCalibratedResolution(evt.getValueFromVector(1))*dm,2)/2.);
				}
			}
			double D_tag = 0.;
			if(q_eff != 0) D_tag = (1.-2.*abs(w_eff));
			double D_tot = D_tag * D_res;

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
				else if(q_eff==0 && f_evt == 1)h_t_fit_0p->Fill(t_MC,weight);
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
			else {   
				if(q_eff == 0)h_t_untagegged_fit->Fill(t_MC,weight);
				else if(q_eff*f_evt > 0  ){
					h_t_mixed_fit->Fill(t_MC,weight);
					if(w_eff<w_max)h_N_mixed_fit->Fill(t_MC,weight);
				}
				else{ 
					h_t_unmixed_fit->Fill(t_MC,weight);
					if(w_eff<w_max)h_N_unmixed_fit->Fill(t_MC,weight);
				}
			}
		}
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

	/// Plot
	TGraph* graph = new TGraph(2);
	graph->SetPoint(1,min_TAU,0);
	graph->SetPoint(2,max_TAU,0);
	graph->SetLineStyle(kDashed);

	TCanvas* c= new TCanvas("");

	TLegend leg(0.6,0.8,0.9,0.9,"");
	leg.SetLineStyle(0);
	leg.SetLineColor(0);
	leg.SetFillColor(0);
	leg.SetTextFont(22);
	leg.SetTextColor(1);
	leg.SetTextSize(0.05);
	leg.SetTextAlign(12);
	leg.AddEntry((TObject*)0,"#font[22]{LHCb unofficial}","");

	h_t->SetMinimum(0.1);    
	h_t->SetLineColor(kBlack);
	h_t->DrawNormalized("e",1);
	
	h_t_fit->SetLineColor(kBlue);
	h_t_fit->SetLineWidth(3);
	h_t_fit->SetMarkerColor(kBlue); 
	h_t_fit->DrawNormalized("histcsame",1);
	leg.Draw();

	c->Print(((string)OutputDir+"h_t.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t.pdf").c_str());
	gPad->SetLogy(1);
	c->Print(((string)OutputDir+"h_t_log.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_log.pdf").c_str());
	gPad->SetLogy(0);
	
	h_dt->SetMinimum(0);        
	h_dt->SetLineColor(kBlack);
	h_dt->DrawNormalized("e1",1);
	h_dt_fit->SetLineColor(kBlue);
	h_dt_fit->SetLineWidth(3);
	h_dt_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_dt.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_dt.pdf").c_str());
	
	h_eta_OS->SetMinimum(0);        
	h_eta_OS->SetLineColor(kBlack);
	h_eta_OS->DrawNormalized("e1",1);
	h_eta_OS_fit->SetLineColor(kBlue);
	h_eta_OS_fit->SetLineWidth(3);
	h_eta_OS_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_eta_OS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_eta_OS.pdf").c_str());
	
	h_eta_SS->SetMinimum(0);        
	h_eta_SS->SetLineColor(kBlack);
	h_eta_SS->DrawNormalized("e1",1);
	h_eta_SS_fit->SetLineColor(kBlue);
	h_eta_SS_fit->SetLineWidth(3);
	h_eta_SS_fit->DrawNormalized("histcsame",1);
	c->Print(((string)OutputDir+"h_eta_SS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_eta_SS.pdf").c_str());

	h_q_OS->SetMinimum(0);        
	h_q_OS->SetLineColor(kBlack);
	h_q_OS->DrawNormalized("e1",1);
	h_q_OS_fit->SetLineColor(kBlue);
	h_q_OS_fit->SetLineWidth(3);
	h_q_OS_fit->DrawNormalized("histsame",1);
	c->Print(((string)OutputDir+"h_q_OS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_q_OS.pdf").c_str());
	
	h_q_SS->SetMinimum(0);        
	h_q_SS->SetLineColor(kBlack);
	h_q_SS->DrawNormalized("e1",1);
	h_q_SS_fit->SetLineColor(kBlue);
	h_q_SS_fit->SetLineWidth(3);
	h_q_SS_fit->DrawNormalized("histsame",1);
	c->Print(((string)OutputDir+"h_q_SS.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_q_SS.pdf").c_str());
	
	h_f->SetMinimum(0);        
	h_f->SetLineColor(kBlack);
	h_f->DrawNormalized("e1",1);
	h_f_fit->SetLineColor(kBlue);
	h_f_fit->SetLineWidth(3);
	h_f_fit->DrawNormalized("histsame",1);
	c->Print(((string)OutputDir+"h_f.eps").c_str());
	if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_f.pdf").c_str());
	
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
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_mixed_p.pdf").c_str());
		
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
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_t_mixed_m.pdf").c_str());
		
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
		h_asym_p->SetMaximum(max_asym);
		h_asym_p->SetMinimum(-max_asym);	
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
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_asym.pdf").c_str());

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
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_asym_CP_f.pdf").c_str());


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
		if(updateAnaNote)c->Print(("../../../../../TD-AnaNote/latex/figs/fullFit/"+(string)OutputDir +"h_asym_CP_i.pdf").c_str());

	double r_val = (double)r ;	
	double norm_A = 1./(1.+r_val*r_val); 
	double norm_Abar = r_val*r_val/(1.+r_val*r_val); 		 
	cout << "ratio = " << sqrt(norm_Abar/norm_A) << endl;

        m_Kspi->SetMinimum(0.01);
        m_Kspi->SetLineColor(kBlack);
        m_Kspi->DrawNormalized("e1",1);
        m_Kspi_fit->SetLineColor(kBlue);
        m_Kspi_fit->SetLineWidth(3);
        m_Kspi_fit->DrawNormalized("histcsame",1);
    // 	m_Kspi_fit_A->SetLineColor(kRed+1);
        m_Kspi_fit_A->SetLineWidth(2);
        m_Kspi_fit_A->SetFillColor(kGray);
        //m_Kspi_fit_A->SetFillStyle(3353);
        m_Kspi_fit_A->DrawNormalized("histcsame",norm_A);
        m_Kspi_fit_Abar->SetLineColor(kGreen+3);
        m_Kspi_fit_Abar->SetLineWidth(2);
        m_Kspi_fit_Abar->SetFillColor(kGreen+3);
        m_Kspi_fit_Abar->SetFillStyle(3353);
        m_Kspi_fit_Abar->DrawNormalized("histcsame",norm_Abar);
        m_Kspi->DrawNormalized("e1same",1);
        c->Print(((string)OutputDir+"m_Kspi_mod.eps").c_str());
        gPad->SetLogy(1);
        c->Print(((string)OutputDir+"m_Kspi_mod_log.eps").c_str());
        gPad->SetLogy(0);

        m_DKs->SetMinimum(0.01);
        m_DKs->SetLineColor(kBlack);
        m_DKs->DrawNormalized("e1",1);
        m_DKs_fit->SetLineColor(kBlue);
        m_DKs_fit->SetLineWidth(3);
        m_DKs_fit->DrawNormalized("histcsame",1);
        //     m_DKs_fit_A->SetLineColor(kRed+1);
        m_DKs_fit_A->SetLineWidth(2);
        m_DKs_fit_A->SetFillColor(kGray);
        //m_DKs_fit_A->SetFillStyle(3353);
        m_DKs_fit_A->DrawNormalized("histcsame",norm_A);
        m_DKs_fit_Abar->SetLineColor(kGreen+3);
        m_DKs_fit_Abar->SetLineWidth(2);
        m_DKs_fit_Abar->SetFillColor(kGreen+3);
        m_DKs_fit_Abar->SetFillStyle(3353);
        m_DKs_fit_Abar->DrawNormalized("histcsame",norm_Abar);
        m_DKs->DrawNormalized("e1same",1);
        c->Print(((string)OutputDir+"m_DKs_mod.eps").c_str());
        gPad->SetLogy(1);
        c->Print(((string)OutputDir+"m_DKs_mod_log.eps").c_str());
        gPad->SetLogy(0);
        
        m_Dpi->SetMinimum(0.01);
        m_Dpi->SetLineColor(kBlack);
        m_Dpi->DrawNormalized("e1",1);
        m_Dpi_fit->SetLineColor(kBlue);
        m_Dpi_fit->SetLineWidth(3);
        m_Dpi_fit->DrawNormalized("histcsame",1);
        //     m_Dpi_fit_A->SetLineColor(kRed+1);
        m_Dpi_fit_A->SetLineWidth(2);
        m_Dpi_fit_A->SetFillColor(kGray);
        //m_Dpi_fit_A->SetFillStyle(3353);
        m_Dpi_fit_A->DrawNormalized("histcsame",norm_A);
        m_Dpi_fit_Abar->SetLineColor(kGreen+3);
        m_Dpi_fit_Abar->SetLineWidth(2);
        m_Dpi_fit_Abar->SetFillColor(kGreen+3);
        m_Dpi_fit_Abar->SetFillStyle(3353);
        m_Dpi_fit_Abar->DrawNormalized("histcsame",norm_Abar);
        m_Dpi->DrawNormalized("e1same",1);
        c->Print(((string)OutputDir+"m_Dpi_mod.eps").c_str());
        gPad->SetLogy(1);
        c->Print(((string)OutputDir+"m_Dpi_mod_log.eps").c_str());
        gPad->SetLogy(0);
        
        
        m_DKs->SetMinimum(0.01);
        m_DKs->SetLineColor(kBlack);
        m_DKs->DrawNormalized("e1",1);
        m_DKs_fit->SetLineColor(kBlue);
        m_DKs_fit->SetLineWidth(3);
        m_DKs_fit->DrawNormalized("histcsame",1);
        m_DKs_fit_5->SetLineColor(kGray+3);
        m_DKs_fit_5->SetLineWidth(2);
        m_DKs_fit_5->SetFillColor(kGray+3);
        m_DKs_fit_5->SetFillStyle(1001);
        m_DKs_fit_5->DrawNormalized("histcsame",m_DKs_fit_5->Integral()/m_DKs_fit->Integral());
        m_DKs_fit_1->SetLineColor(kRed+1);
        m_DKs_fit_1->SetLineWidth(2);
        //m_DKs_fit_1->SetFillColor(kRed+1);
        //m_DKs_fit_1->SetFillStyle(3353);
        m_DKs_fit_1->DrawNormalized("histcsame",m_DKs_fit_1->Integral()/m_DKs_fit->Integral());
        m_DKs_fit_2->SetLineColor(kGreen+3);
        m_DKs_fit_2->SetLineWidth(2);
        //m_DKs_fit_2->SetFillColor(kGreen+3);
        //m_DKs_fit_2->SetFillStyle(3353);
        m_DKs_fit_2->DrawNormalized("histcsame",m_DKs_fit_2->Integral()/m_DKs_fit->Integral());
        m_DKs_fit_3->SetLineColor(kMagenta+3);
        m_DKs_fit_3->SetLineWidth(2);
        //m_DKs_fit_3->SetFillColor(kMagenta+3);
        //m_DKs_fit_3->SetFillStyle(3353);
        m_DKs_fit_3->DrawNormalized("histcsame",m_DKs_fit_3->Integral()/m_DKs_fit->Integral());
        m_DKs_fit_4->SetLineColor(kBlack);
        m_DKs_fit_4->SetLineWidth(3);
        m_DKs_fit_4->SetLineStyle(kDashed);
        m_DKs_fit_4->DrawNormalized("histcsame",m_DKs_fit_4->Integral()/m_DKs_fit->Integral());
        m_DKs->DrawNormalized("e1same",1);
        c->Print(((string)OutputDir+"m_DKs.eps").c_str());
        gPad->SetLogy(1);
        c->Print(((string)OutputDir+"m_DKs_log.eps").c_str());
        gPad->SetLogy(0);
        
        m_Dpi->SetMinimum(0.01);
        m_Dpi->SetLineColor(kBlack);
        m_Dpi->DrawNormalized("e1",1);
        m_Dpi_fit->SetLineColor(kBlue);
        m_Dpi_fit->SetLineWidth(3);
        m_Dpi_fit->DrawNormalized("histcsame",1);
        m_Dpi_fit_5->SetLineColor(kGray+3);
        m_Dpi_fit_5->SetLineWidth(2);
        m_Dpi_fit_5->SetFillColor(kGray+3);
        m_Dpi_fit_5->SetFillStyle(1001);
        m_Dpi_fit_5->DrawNormalized("histcsame",m_Dpi_fit_5->Integral()/m_Dpi_fit->Integral());
        m_Dpi_fit_1->SetLineColor(kRed+1);
        m_Dpi_fit_1->SetLineWidth(2);
        //m_Dpi_fit_1->SetFillColor(kRed+1);
        //m_Dpi_fit_1->SetFillStyle(3353);
        m_Dpi_fit_1->DrawNormalized("histcsame",m_Dpi_fit_1->Integral()/m_Dpi_fit->Integral());
        m_Dpi_fit_2->SetLineColor(kGreen+3);
        m_Dpi_fit_2->SetLineWidth(2);
        //m_Dpi_fit_2->SetFillColor(kGreen+3);
        //m_Dpi_fit_2->SetFillStyle(3353);
        m_Dpi_fit_2->DrawNormalized("histcsame",m_Dpi_fit_2->Integral()/m_Dpi_fit->Integral());
        m_Dpi_fit_3->SetLineColor(kMagenta+3);
        m_Dpi_fit_3->SetLineWidth(2);
        //m_Dpi_fit_3->SetFillColor(kMagenta+3);
        //m_Dpi_fit_3->SetFillStyle(3353);
        m_Dpi_fit_3->DrawNormalized("histcsame",m_Dpi_fit_3->Integral()/m_Dpi_fit->Integral());
        m_Dpi_fit_4->SetLineColor(kBlack);
        m_Dpi_fit_4->SetLineWidth(3);
        m_Dpi_fit_4->SetLineStyle(kDashed);
        m_Dpi_fit_4->DrawNormalized("histcsame",m_Dpi_fit_4->Integral()/m_Dpi_fit->Integral());
        m_Dpi->DrawNormalized("e1same",1);
        c->Print(((string)OutputDir+"m_Dpi.eps").c_str());
        gPad->SetLogy(1);
        c->Print(((string)OutputDir+"m_Dpi_log.eps").c_str());
        gPad->SetLogy(0);
        
        cout << "int" << endl;
        cout << m_Kspi_fit->Integral() << endl;
        cout << m_Kspi_fit_1->Integral() << endl;
        cout << m_Kspi_fit_2->Integral() << endl;
        cout << m_Kspi_fit_3->Integral() << endl;
        cout << m_Kspi_fit_4->Integral() << endl;
        cout << m_Kspi_fit_5->Integral() << endl;
        cout << m_Kspi_fit_A->Integral() << endl;
        cout << m_Kspi_fit_Abar->Integral() << endl;
        
        m_Kspi->SetMinimum(0.0001);
        m_Kspi->SetLineColor(kBlack);
        m_Kspi->DrawNormalized("e1",1);
        m_Kspi_fit->SetLineColor(kBlue);
        m_Kspi_fit->SetLineWidth(3);
        m_Kspi_fit->DrawNormalized("histcsame",1);
        m_Kspi_fit_5->SetLineColor(kGray+3);
        m_Kspi_fit_5->SetLineWidth(2);
        m_Kspi_fit_5->SetFillColor(kGray+3);
        m_Kspi_fit_5->SetFillStyle(1001);
        m_Kspi_fit_5->DrawNormalized("histcsame",m_Kspi_fit_5->Integral()/m_Kspi_fit->Integral());
        m_Kspi_fit_1->SetLineColor(kRed+1);
        m_Kspi_fit_1->SetLineWidth(2);
        //m_Kspi_fit_1->SetFillColor(kRed+1);
        //m_Kspi_fit_1->SetFillStyle(3353);
        m_Kspi_fit_1->DrawNormalized("histcsame",m_Kspi_fit_1->Integral()/m_Kspi_fit->Integral());
        m_Kspi_fit_2->SetLineColor(kGreen+3);
        m_Kspi_fit_2->SetLineWidth(2);
        //m_Kspi_fit_2->SetFillColor(kGreen+3);
        //m_Kspi_fit_2->SetFillStyle(3353);
        m_Kspi_fit_2->DrawNormalized("histcsame",m_Kspi_fit_2->Integral()/m_Kspi_fit->Integral());
        m_Kspi_fit_3->SetLineColor(kMagenta+3);
        m_Kspi_fit_3->SetLineWidth(2);
        //m_Kspi_fit_3->SetFillColor(kMagenta+3);
        //m_Kspi_fit_3->SetFillStyle(3353);
        m_Kspi_fit_3->DrawNormalized("histcsame",m_Kspi_fit_3->Integral()/m_Kspi_fit->Integral());
        m_Kspi_fit_4->SetLineColor(kBlack);
        m_Kspi_fit_4->SetLineWidth(3);
        m_Kspi_fit_4->SetLineStyle(kDashed);
        m_Kspi_fit_4->DrawNormalized("histcsame",m_Kspi_fit_4->Integral()/m_Kspi_fit->Integral());
        m_Kspi->DrawNormalized("e1same",1);
        c->Print(((string)OutputDir+"m_Kspi.eps").c_str());
        gPad->SetLogy(1);
        c->Print(((string)OutputDir+"m_Kspi_log.eps").c_str());
        gPad->SetLogy(0);
        
        
        c->Clear();
        TLegend leg_mod(0.,0.,1,1,"");
        leg_mod.SetLineStyle(0);
        leg_mod.SetLineColor(0);
        leg_mod.SetFillColor(0);
        leg_mod.SetTextFont(22);
        leg_mod.SetTextColor(1);
        leg_mod.SetTextSize(0.075);
        leg_mod.SetTextAlign(12);
        leg_mod.AddEntry(m_Kspi_fit,"Full PDF","l");
        leg_mod.AddEntry(m_Kspi_fit_A,"|A^{c}(x)|^{2}","f");
        leg_mod.AddEntry(m_Kspi_fit_Abar,"r^{2}|A^{u}(x)|^{2}","f");
        leg_mod.Draw();
        c->Print(((string)OutputDir+"leg_mod.eps").c_str());


	int nPars = 0;
	for(int i=0; i < mps->size(); i++){
		if(((FitParameter*)mps->getParPtr(i))->iFixInit()!=0)continue;
		nPars++;
	}

	//vector<double> chi2 = getChi2(eventList,eventListMC_rw);
	//double chi2_val = chi2[0]/(chi2[1]-1.);
	//double nu_val = (double)nPars;
	//cout << "chi2 = " << chi2_val << endl;	
	//cout << "nPars = " << nPars << endl;	
   }

    /*
   if(useLASSO>0){
		paraFile->cd();

		double N_sig = 0;
       		for (int i=0; i<eventList.size(); i++) N_sig += eventList[i].getWeight();

		Double_t x[1],y[1];
		vector<double> thresholds;
		thresholds.push_back(0.001);
		thresholds.push_back(0.002);
		thresholds.push_back(0.003);
		thresholds.push_back(0.004);
		thresholds.push_back(0.005);
		thresholds.push_back(0.006);
		thresholds.push_back(0.007);
		thresholds.push_back(0.008);
		thresholds.push_back(0.009);
		thresholds.push_back(0.01);
		thresholds.push_back(0.02);
		thresholds.push_back(0.05);
		
		x[0]=lambda;
		for(int i = 0; i < thresholds.size() ; i++){
			if(doSimFit)y[0]=neg2LL_sim.getVal(); 
			else y[0]=neg2LL.getVal();
			if(useLASSO==1)y[0]+= 2. * lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			y[0]+= 2. * lasso_bar.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			TGraph* aic = new TGraph(1,x,y);
			aic->SetName( ("AIC_"+anythingToString((int) (thresholds[i]*1000))).c_str());
			aic->SetTitle("");
			aic->GetXaxis()->SetTitle("#lambda");
			aic->GetXaxis()->SetTitleOffset(0.65);
			aic->GetYaxis()->SetTitle("AIC");
			aic->Draw("A*");
			aic->Write();
		}
		
		for(int i = 0; i < thresholds.size() ; i++){
			if(doSimFit)y[0]=neg2LL_sim.getVal(); 
			else y[0]=neg2LL.getVal();
			if(useLASSO==1)y[0]+= lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]) * log(N_sig);
			y[0]+= lasso_bar.numberOfFitFractionsLargerThanThreshold(thresholds[i]) * log(N_sig);
			TGraph* bic = new TGraph(1,x,y);
			bic->SetName( ("BIC_"+anythingToString((int) (thresholds[i]*1000))).c_str());
			bic->SetTitle("");
			bic->GetXaxis()->SetTitle("#lambda");
			bic->GetXaxis()->SetTitleOffset(0.65);
			bic->GetYaxis()->SetTitle("BIC");
			bic->Draw("A*");
			bic->Write();
		}
		
		for(int i = 0; i < thresholds.size() ; i++){
			y[0]=lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			y[0]+= lasso_bar.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			TGraph* r = new TGraph(1,x,y);
			r->SetName( ("r_"+anythingToString((int) (thresholds[i]*1000))).c_str());
			r->SetTitle("");
			r->GetXaxis()->SetTitle("#lambda");
			r->GetXaxis()->SetTitleOffset(0.65);
			r->GetYaxis()->SetTitle("Number of fit fractions larger than threshold");
			r->Draw("A*");
			r->Write();
		}

		for(int i = 0; i < thresholds.size() ; i++){

			y[0]=lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			TGraph* r = new TGraph(1,x,y);
			r->SetName( ("r_A_"+anythingToString((int) (thresholds[i]*1000))).c_str());
			r->SetTitle("");
			r->GetXaxis()->SetTitle("#lambda");
			r->GetXaxis()->SetTitleOffset(0.65);
			r->GetYaxis()->SetTitle("Number of fit fractions larger than threshold");
			r->Draw("A*");
			r->Write();
		}

		for(int i = 0; i < thresholds.size() ; i++){
			y[0]=lasso_bar.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
			TGraph* r = new TGraph(1,x,y);
			r->SetName( ("r_Abar_"+anythingToString((int) (thresholds[i]*1000))).c_str());
			r->SetTitle("");
			r->GetXaxis()->SetTitle("#lambda");
			r->GetXaxis()->SetTitleOffset(0.65);
			r->GetYaxis()->SetTitle("Number of fit fractions larger than threshold");
			r->Draw("A*");
			r->Write();
		}

		if(doSimFit)y[0]=neg2LL_sim.getVal() ;
		else y[0]=neg2LL.getVal() ;
		TGraph* nll = new TGraph(1,x,y);
		nll->SetName("Neg2LL");
		nll->SetTitle("");
		nll->GetXaxis()->SetTitle("#lambda");
		nll->GetXaxis()->SetTitleOffset(0.65);
		nll->GetYaxis()->SetTitle("-2 logL");
		nll->Draw("A*");
		nll->Write();
		
		y[0]= lasso.sumOfFitFractions() ;
		TGraph* ff = new TGraph(1,x,y);
		ff->SetName("SumOfFitFractions A");
		ff->SetTitle("");
		ff->GetXaxis()->SetTitle("#lambda");
		ff->GetXaxis()->SetTitleOffset(0.65);
		ff->GetYaxis()->SetTitle("Total Fit Fraction");
		ff->Draw("A*");
		ff->Write();
		
		y[0]= lasso.absSumOfInterferenceFractions() ;
		TGraph* iff = new TGraph(1,x,y);
		iff->SetName("AbsSumOfInterferenceFractions A");
		iff->SetTitle("");
		iff->GetXaxis()->SetTitle("#lambda");
		iff->GetXaxis()->SetTitleOffset(0.65);
		iff->GetYaxis()->SetTitle("Sum of abs(Interference Fraction)");
		iff->Draw("A*");
		iff->Write();

		y[0]= lasso_bar.sumOfFitFractions() ;
		TGraph* ff_bar = new TGraph(1,x,y);
		ff_bar->SetName("SumOfFitFractions Abar");
		ff_bar->SetTitle("");
		ff_bar->GetXaxis()->SetTitle("#lambda");
		ff_bar->GetXaxis()->SetTitleOffset(0.65);
		ff_bar->GetYaxis()->SetTitle("Total Fit Fraction");
		ff_bar->Draw("A*");
		ff_bar->Write();
		
		y[0]= lasso_bar.absSumOfInterferenceFractions() ;
		TGraph* iff_bar = new TGraph(1,x,y);
		iff_bar->SetName("AbsSumOfInterferenceFractions Abar");
		iff_bar->SetTitle("");
		iff_bar->GetXaxis()->SetTitle("#lambda");
		iff_bar->GetXaxis()->SetTitleOffset(0.65);
		iff_bar->GetYaxis()->SetTitle("Sum of abs(Interference Fraction)");
		iff_bar->Draw("A*");
		iff_bar->Write();

		y[0]= chi2_val ;
		TGraph* chi2 = new TGraph(1,x,y);
		chi2->SetName("Chi2");
		chi2->SetTitle("");
		chi2->GetXaxis()->SetTitle("#lambda");
		chi2->GetXaxis()->SetTitleOffset(0.65);
		chi2->GetYaxis()->SetTitle("#Chi^{2}/Bin");
		chi2->Draw("A*");
		chi2->Write();
	}
    */
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


void produceIntegratorFile_CP(){
    
    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    TString integratorEventFile = (string) IntegratorEventFile;

    DiskResidentEventList eventListMC(((string) IntegratorEventFile).c_str(),"OPEN");
    DiskResidentEventList eventList(((string) integratorEventFile.ReplaceAll(".root","_CP.root")).c_str(),"RECREATE");

    for(int i = 0; i < eventListMC.size()/100; i++){
	DalitzEvent evt(eventListMC.getEvent(i));
	evt.CP_conjugateYourself();
	evt.P_conjugateYourself();
	eventList.Add(evt);
    }

    eventList.save();
    return;
}

void makeIntegratorFileForToys(int step = 0){
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    int seed = RandomSeed + step;
    ranLux.SetSeed((int)seed);
    gRandom = &ranLux;
    
    FitAmplitude::AutogenerateFitFile();
    
    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    NamedParameter<int>  IntegratorEvents("IntegratorEvents", 300000);    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    DalitzEventList eventListPhsp,eventList,eventList_cut,eventList_cut_CP;    
    eventListPhsp.generatePhaseSpaceEvents(100000,pat);
    
    FitAmpIncoherentSum fas((DalitzEventPattern)pat);
    fas.print();
    fas.getVal(eventListPhsp[0]);
    //fas.normalizeAmps(eventListPhsp);
    
    SignalGenerator sg(pat,&fas);
    
    sg.FillEventList(eventList, IntegratorEvents);
    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);
    
    for(int i = 0; i < eventList.size(); i++){
        if(sqrt(eventList[i].sij(s234)/(GeV*GeV)) < 1.95 && sqrt(eventList[i].s(2,4)/(GeV*GeV)) < 1.2 && sqrt(eventList[i].s(3,4)/(GeV*GeV)) < 1.2){
            eventList_cut.Add(eventList[i]);
            DalitzEvent evt(eventList[i]);
            evt.CP_conjugateYourself();
            eventList_cut_CP.Add(evt);
        }
    }
    
    cout << "Generated " << eventList_cut.size() << " events inside selected phasespace region" << endl;
    
    TString outputName = (string)IntegratorEventFile;
    if(step>0) outputName.ReplaceAll(".root",("_" + anythingToString(step) + ".root").c_str());
    eventList_cut.saveAsNtuple((string)outputName);
    eventList_cut_CP.saveAsNtuple((string)outputName.ReplaceAll(".root","_CP.root"));
    return;
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
	RooRealVar mean("mean", "mean", 5366.89,5350.,5390.); 
	RooRealVar sigma("sigma", "sigma", 20.,0.,80.); 
	RooGaussian* signal = new RooGaussian("signal","signal",DTF_Bs_M, mean,sigma);

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

  time_t startTime = time(0);

  TH1::SetDefaultSumw2();
  TH2::SetDefaultSumw2();

  gROOT->ProcessLine(".x ../lhcbStyle.C");

   //produceMarginalPdfs();
   //produceIntegratorFile_CP();
  //makeIntegratorFileForToys(atoi(argv[1]));
    
  NamedParameter<int>  addBkgToToys("addBkgToToys", 0);


  //for(int i = 0; i < 200; i++)ampFit(atoi(argv[1])+i);
     ampFit(atoi(argv[1]),(string)argv[2]);
     if((string)argv[2] == "gen" && addBkgToToys)calculateSweightsForToys(atoi(argv[1]));

//  	for(int i = 1; i <= 200; i++)calculateSweightsForToys(i);


//    animate(atoi(argv[1]));

  cout << "==============================================" << endl;
  cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
  cout << "==============================================" << endl;
  
  return 0;
}
//
