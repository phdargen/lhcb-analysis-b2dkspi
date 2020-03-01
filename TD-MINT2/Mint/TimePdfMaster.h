#ifndef TIMEPDFMASTER_HH
#define TIMEPDFMASTER_HH
// author: Philippe d'Argent

#include <complex>
#include <vector>
#include <iostream>

#include <TH1D.h>
#include <TCanvas.h>
#include <TRandom3.h>

#include "Mint/IReturnRealForEvent.h"
#include "Mint/IReturnComplexForEvent.h"
#include "Mint/DalitzEvent.h"
#include "Mint/FitParRef.h"
#include "Mint/FitParDependent.h"
#include "Mint/IFitParRegister.h"
//#include "Mint/CachedByEvent.h"

#include "RooConstVar.h"
#include "RooAbsReal.h"
#include "RooUniform.h"
#include "RooRealVar.h"
#include "RooRealConstant.h"
#include "RooCategory.h"
#include "RooDataHist.h"
#include "RooHistPdf.h"
#include "RooBDecay.h"
#include "RooEffProd.h"
#include "RooProdPdf.h"
#include "RooProduct.h"
#include "RooGaussModel.h"
#include "RooMCStudy.h"

#include "Mint/RooGaussEfficiencyModel.h"
#include "Mint/RooCubicSplineFun.h"
#include "Mint/RooCubicSplinePdf.h"
#include "Mint/DecRateCoeff_Bd.h"
#include "Mint/TimePdfIntegrator.h"


enum basisType { 
      noBasis=0  ,  expBasis= 3
    , sinBasis=13,  cosBasis=23
    , sinhBasis=63, coshBasis=53 };

using namespace std;
using namespace RooFit ;
using namespace MINT;

class TimePdfMaster
{
 protected:
    
    // Fit parameters
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

    /// Cast of MINT parameters to B2DX parameters
    RooRealVar* _r_t;
    RooRealVar* _r_dt;
    RooRealVar* _r_mean_scaled;
    RooRealVar* _r_dt_scaled;
    RooRealVar* _r_dt_scaled2;
    RooRealVar* _r_dt_scaled3;
    RooCategory* _r_q_OS;
    RooRealVar* _r_eta_OS;
    RooCategory* _r_q_SS;
    RooRealVar* _r_eta_SS;
    RooCategory* _r_f;
    RooCategory* _r_run;
    RooCategory* _r_trigger;

    RooRealVar* _r_tau;
    RooRealVar* _r_dGamma;
    RooRealVar* _r_dm;
    
    RooRealVar* _r_offset_mean_dt;
    RooRealVar* _r_scale_mean_dt;
    RooRealVar* _r_scale_mean_2_dt;  

    RooRealVar* _r_offset_sigma_dt;
    RooRealVar* _r_scale_sigma_dt;
    RooRealVar* _r_scale_sigma_2_dt;  

    RooRealVar* _r_offset_sigma2_dt;
    RooRealVar* _r_scale_sigma2_dt;
    RooRealVar* _r_scale_sigma2_2_dt;  

    RooRealVar* _r_offset_sigma3_dt;
    RooRealVar* _r_scale_sigma3_dt;
    RooRealVar* _r_scale_sigma3_2_dt;  

    RooRealVar* _r_offset_f_dt;
    RooRealVar* _r_scale_f_dt;
    RooRealVar* _r_scale_f_2_dt;  

    RooRealVar* _r_offset_f2_dt;
    RooRealVar* _r_scale_f2_dt;
    RooRealVar* _r_scale_f2_2_dt;  

    RooRealVar* _r_c0;
    RooRealVar* _r_c1;
    RooRealVar* _r_c2;
    RooRealVar* _r_c3;
    RooRealVar* _r_c4;
    RooRealVar* _r_c5;
    RooRealVar* _r_c6;
    RooRealVar* _r_c7;
    RooRealVar* _r_c8;
    RooRealVar* _r_c9;
    RooFormulaVar* _coeff_last;

    RooRealVar* _r_p0_os;
    RooRealVar* _r_p1_os;
    RooRealVar* _r_delta_p0_os;
    RooRealVar* _r_delta_p1_os;
    RooRealVar* _r_avg_eta_os;
    RooRealVar* _r_tageff_os;
    RooRealVar* _r_tageff_asym_os;
    RooRealVar* _r_p0_ss;
    RooRealVar* _r_p1_ss;
    RooRealVar* _r_delta_p0_ss;
    RooRealVar* _r_delta_p1_ss;
    RooRealVar* _r_avg_eta_ss;
    RooRealVar* _r_tageff_ss;
    RooRealVar* _r_tageff_asym_ss;
    
    RooRealVar* _r_production_asym;
    RooRealVar* _r_detection_asym;

    // Acceptance function
    RooCubicSplineFun* _spline;
    RooProduct* _splinePdf;
    RooGaussEfficiencyModel* _efficiency;
    RooGaussEfficiencyModel* _efficiency2;
    RooGaussEfficiencyModel* _efficiency3;
    
    // CP coefficients
    RooRealVar* _r_norm;
    RooRealVar* _r_norm_bar;
    RooRealVar* _r_C;
    RooRealVar* _r_C_bar;
    RooRealVar* _r_D;
    RooRealVar* _r_D_bar;
    RooRealVar* _r_S;
    RooRealVar* _r_S_bar;
    
    // Decay rate coefficients
    DecRateCoeff_Bd* _cos_coeff;
    DecRateCoeff_Bd* _cosh_coeff;
    DecRateCoeff_Bd* _sin_coeff;
    DecRateCoeff_Bd* _sinh_coeff;
    
    // Time pdf integrators
    TimePdfIntegrator* _cosh_term;
    TimePdfIntegrator* _sinh_term; 
    TimePdfIntegrator* _cos_term; 
    TimePdfIntegrator* _sin_term;
    
    // Marginal pdfs
    string _marginalPdfsPrefix;
    TFile* _f_pdfs;

    TH1D* _h_dt;
    RooDataHist* _r_h_dt;
    RooAbsPdf* _pdf_sigma_t;
    
    TH1D* _h_q_f;
    TH1D* _h_q_OS;
    TH1D* _h_q_SS;

    RooAbsPdf* _pdf_eta_OS;
    TH1D* _h_eta_OS;
    RooDataHist* _r_h_eta_OS;
    
    RooAbsPdf* _pdf_eta_SS;
    TH1D* _h_eta_SS;
    RooDataHist* _r_h_eta_SS;

    // Limits
    NamedParameter<double> _min_TAU;
    NamedParameter<double> _max_TAU;
    NamedParameter<double> _min_TAUERR;
    NamedParameter<double> _max_TAUERR;

    // 
    RooBDecay* _samplingPdf_t;
    RooBDecay* _fitPdf_t;
    RooGaussModel* _resmodel;
    RooEffProd* _samplingPdf_t_eff;
    RooProdPdf* _samplingPdf;
    RooProdPdf* _fitPdf;
    RooDataSet* _protoData; 
    RooMCStudy* _sampleGen;

 public:
    TimePdfMaster(const MINT::FitParameter& Gamma, const MINT::FitParameter& dGamma, const MINT::FitParameter& dm
                  ,const MINT::FitParameter& offset_mean_dt, const MINT::FitParameter& scale_mean_dt, const MINT::FitParameter& scale_mean_2_dt
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
                  const MINT::FitParameter& production_asym, const MINT::FitParameter& detection_asym, string marginalPdfsPrefix ): 
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
        _min_TAU("min_TAU", 0.4),
        _max_TAU("max_TAU", 10.),
        _min_TAUERR("min_TAUERR", 0.),
        _max_TAUERR("max_TAUERR", 0.1)
    {
        /// Init B2DX parameters
        
        // Observables
        _r_t = new RooRealVar("t", "t", _min_TAU, _max_TAU);
        _r_dt = new RooRealVar("dt", "dt",_min_TAUERR,_max_TAUERR);
        _r_f = new RooCategory("qf", "qf");
        _r_f->defineType("h+", +1);
        _r_f->defineType("h-", -1);
        _r_f->setRange("Range_p","h+");
        _r_f->setRange("Range_m","h-");
        _r_q_OS = new RooCategory("q_OS", "q_OS");
        _r_q_OS->defineType("B+", +1);
        _r_q_OS->defineType("B-", -1) ;
        _r_q_OS->defineType("untagged", 0);
        _r_eta_OS = new RooRealVar("eta_OS", "eta_OS",0.,0.5);
        _r_q_SS = new RooCategory("q_SS", "q_SS");
        _r_q_SS->defineType("B+", +1);
        _r_q_SS->defineType("B-", -1) ;
        _r_q_SS->defineType("untagged", 0);
        _r_eta_SS = new RooRealVar("eta_SS", "eta_SS",0.,0.5);
        _r_run = new RooCategory("run", "run");
        _r_run->defineType("Run1", 1);
        _r_run->defineType("Run2", 2);
        _r_trigger = new RooCategory("trigger", "trigger");
        _r_trigger->defineType("t0", 0);
        _r_trigger->defineType("t1", 1);

        // Fit parameters
        _r_tau = new RooRealVar("tau", "tau",1./_Gamma);
        _r_dGamma = new RooRealVar("dGamma", "dGamma",_dGamma);
        _r_dm = new RooRealVar("dm", "dm",_dm);
        
        _r_offset_mean_dt = new RooRealVar("offset_mean_dt", "offset_mean_dt", _offset_mean_dt);
        _r_scale_mean_dt = new RooRealVar("scale_mean_dt", "scale_mean_dt", _scale_mean_dt);
        _r_scale_mean_2_dt = new RooRealVar("scale_mean_2_dt", "scale_mean_2_dt", _scale_mean_2_dt);

        _r_offset_sigma_dt = new RooRealVar("offset_sigma_dt", "offset_sigma_dt", _offset_sigma_dt);
        _r_scale_sigma_dt = new RooRealVar("scale_sigma_dt", "scale_sigma_dt", _scale_sigma_dt);
	_r_scale_sigma_2_dt = new RooRealVar("scale_sigma_2_dt", "scale_sigma_2_dt", _scale_sigma_2_dt);      

        _r_offset_sigma2_dt = new RooRealVar("offset_sigma2_dt", "offset_sigma2_dt", _offset_sigma2_dt);
        _r_scale_sigma2_dt = new RooRealVar("scale_sigma2_dt", "scale_sigma2_dt", _scale_sigma2_dt);
	_r_scale_sigma2_2_dt = new RooRealVar("scale_sigma2_2_dt", "scale_sigma2_2_dt", _scale_sigma2_2_dt);      

        _r_offset_sigma3_dt = new RooRealVar("offset_sigma3_dt", "offset_sigma3_dt", _offset_sigma3_dt);
        _r_scale_sigma3_dt = new RooRealVar("scale_sigma3_dt", "scale_sigma3_dt", _scale_sigma3_dt);
	_r_scale_sigma3_2_dt = new RooRealVar("scale_sigma3_2_dt", "scale_sigma3_2_dt", _scale_sigma3_2_dt);      

        _r_offset_f_dt = new RooRealVar("offset_f_dt", "offset_f_dt", _offset_f_dt);
        _r_scale_f_dt = new RooRealVar("scale_f_dt", "scale_f_dt", _scale_f_dt);
	_r_scale_f_2_dt = new RooRealVar("scale_f_2_dt", "scale_f_2_dt", _scale_f_2_dt);      

        _r_offset_f2_dt = new RooRealVar("offset_f2_dt", "offset_f2_dt", _offset_f2_dt);
        _r_scale_f2_dt = new RooRealVar("scale_f2_dt", "scale_f2_dt", _scale_f2_dt);
	_r_scale_f2_2_dt = new RooRealVar("scale_f2_2_dt", "scale_f2_2_dt", _scale_f2_2_dt);      

        _r_c0 = new RooRealVar("coeff_0", "coeff_0",_c0);
        _r_c1 = new RooRealVar("coeff_1", "coeff_1",_c1);
        _r_c2 = new RooRealVar("coeff_2", "coeff_2",_c2);
        _r_c3 = new RooRealVar("coeff_3", "coeff_3",_c3);
        _r_c4 = new RooRealVar("coeff_4", "coeff_4",_c4);
        _r_c5 = new RooRealVar("coeff_5", "coeff_5",_c5);
        _r_c6 = new RooRealVar("coeff_6", "coeff_6",_c6);
        _r_c7 = new RooRealVar("coeff_7", "coeff_7",_c7);
        _r_c8 = new RooRealVar("coeff_8", "coeff_8",_c8);
        _r_c9 = new RooRealVar("coeff_9", "coeff_9",_c9);
        
        _r_p0_os = new RooRealVar("p0_os", "p0_os",_p0_os);
        _r_p1_os = new RooRealVar("p1_os", "p1_os",_p1_os);
        _r_delta_p0_os = new RooRealVar("delta_p0_os", "delta_p0_os",_delta_p0_os);
        _r_delta_p1_os = new RooRealVar("delta_p1_os", "delta_p1_os",_delta_p1_os);
        _r_avg_eta_os = new RooRealVar("avg_eta_os", "avg_eta_os",_avg_eta_os);
        _r_tageff_os = new RooRealVar("tageff_os", "tageff_os",_tageff_os);
        _r_tageff_asym_os = new RooRealVar("tageff_asym_os", "tageff_asym_os",_tageff_asym_os);
        _r_p0_ss = new RooRealVar("p0_ss", "p0_ss",_p0_ss);
        _r_p1_ss = new RooRealVar("p1_ss", "p1_ss",_p1_ss);
        _r_delta_p0_ss = new RooRealVar("delta_p0_ss", "delta_p0_ss",_delta_p0_ss);
        _r_delta_p1_ss = new RooRealVar("delta_p1_ss", "delta_p1_ss",_delta_p1_ss);
        _r_avg_eta_ss = new RooRealVar("avg_eta_ss", "avg_eta_ss",_avg_eta_ss);
        _r_tageff_ss = new RooRealVar("tageff_ss", "tageff_ss",_tageff_ss);
        _r_tageff_asym_ss= new RooRealVar("tageff_asym_ss", "tageff_asym_ss",_tageff_asym_ss);
        
        _r_production_asym = new RooRealVar("production_asym", "production_asym",_production_asym);
        _r_detection_asym = new RooRealVar("detection_asym", "detection_asym",_detection_asym);
        
        /// Acceptance
        vector<RooRealVar*> v_coeff;
        v_coeff.push_back(_r_c0);
        v_coeff.push_back(_r_c1);
        v_coeff.push_back(_r_c2);
        v_coeff.push_back(_r_c3);
        v_coeff.push_back(_r_c4);
        v_coeff.push_back(_r_c5);
        v_coeff.push_back(_r_c6);
        v_coeff.push_back(_r_c7);
        v_coeff.push_back(_r_c8);
        v_coeff.push_back(_r_c9);
 
        //SPLINE KNOTS
        NamedParameter<double> knot_positions("knot_positions", 0.5, 1., 1.5, 2., 3., 6., 9.5);
        vector<double> myBinning = knot_positions.getVector();
        
        RooArgList tacc_list;
        for(int i= 0; i<= myBinning.size(); i++){
            tacc_list.add(*v_coeff[i]);
        }
        
        _coeff_last = new RooFormulaVar(("coeff_"+anythingToString((int)myBinning.size()+1)).c_str(),("coeff_"+anythingToString((int)myBinning.size()+1)).c_str(), "@0 + ((@0-@1)/(@2-@3)) * (@4 - @2)", RooArgList(RooRealConstant::value(1.0), *tacc_list.find(("coeff_"+anythingToString((int)myBinning.size()-1)).c_str())  , RooRealConstant::value(myBinning[myBinning.size()-1]), RooRealConstant::value(myBinning[myBinning.size()-2]), RooRealConstant::value(_r_t->getMax()) ));
        
        tacc_list.add(*_coeff_last);	
        
        // Cubic spline function
        _spline = new RooCubicSplineFun("spline", "spline", *_r_t, myBinning, tacc_list);        
	_splinePdf = new RooProduct("splinePdf", "splinePdf", RooArgList(*_spline,RooRealConstant::value(0.1)));

        _r_mean_scaled = (RooRealVar*) new RooFormulaVar( "r_mean_scaled","r_mean_scaled", "@0+@1*@2+@3*@2*@2",RooArgList(*_r_offset_mean_dt,*_r_scale_mean_dt,*_r_dt,*_r_scale_mean_2_dt));
        _r_dt_scaled = (RooRealVar*) new RooFormulaVar( "r_dt_scaled","r_dt_scaled", "@0+@1*@2+@3*@2*@2",RooArgList(*_r_offset_sigma_dt,*_r_scale_sigma_dt,*_r_dt,*_r_scale_sigma_2_dt));
        _r_dt_scaled2 = (RooRealVar*) new RooFormulaVar( "r_dt_scaled2","r_dt_scaled2", "@0+@1*@2+@3*@2*@2",RooArgList(*_r_offset_sigma2_dt,*_r_scale_sigma2_dt,*_r_dt,*_r_scale_sigma2_2_dt));
        _r_dt_scaled3 = (RooRealVar*) new RooFormulaVar( "r_dt_scaled3","r_dt_scaled3", "@0+@1*@2+@3*@2*@2",RooArgList(*_r_offset_sigma3_dt,*_r_scale_sigma3_dt,*_r_dt,*_r_scale_sigma3_2_dt));

        _efficiency = new RooGaussEfficiencyModel("resmodel", "resmodel", *_r_t, *_spline, *_r_mean_scaled, *_r_dt_scaled, RooRealConstant::value(1.),RooRealConstant::value(1.) );
        _efficiency2 = new RooGaussEfficiencyModel("resmodel2", "resmodel2", *_r_t, *_spline, *_r_mean_scaled, *_r_dt_scaled2, RooRealConstant::value(1.),RooRealConstant::value(1.) );
        _efficiency3 = new RooGaussEfficiencyModel("resmodel3", "resmodel3", *_r_t, *_spline, *_r_mean_scaled, *_r_dt_scaled3, RooRealConstant::value(1.),RooRealConstant::value(1.) );

        // CP coefficients
        _r_norm = new RooRealVar("norm", "norm",1);
        _r_norm_bar = new RooRealVar("norm_bar", "norm_bar",1);
        _r_C = new RooRealVar("C", "C",1);
        _r_C_bar = new RooRealVar("Cbar", "Cbar",-1);
        //_r_C_bar= (RooRealVar*) new RooFormulaVar("Cbar","-1. * @0",RooArgList(*_r_C));
        _r_D= new RooRealVar("D", "D",0);
        _r_D_bar= new RooRealVar("Dbar", "Dbar",0);
        _r_S= new RooRealVar("S", "S",0);
        _r_S_bar= new RooRealVar("Sbar", "Sbar",0);
        
        /// Decay rate coefficients
        _cosh_coeff = new DecRateCoeff_Bd("cosh_coeff",
                                          "cosh_coeff",
                                          DecRateCoeff_Bd::kCosh ,
                                          *_r_f,
                                          *_r_norm,
                                          *_r_norm_bar,
                                          *_r_q_OS,
                                          *_r_eta_OS,
                                          *_r_p0_os,
                                          *_r_p1_os,
                                          *_r_delta_p0_os,
                                          *_r_delta_p1_os,
                                          *_r_avg_eta_os,
                                          *_r_tageff_os,
                                          *_r_tageff_asym_os,
                                          *_r_q_SS,
                                          *_r_eta_SS,
                                          *_r_p0_ss,
                                          *_r_p1_ss,
                                          *_r_delta_p0_ss,
                                          *_r_delta_p1_ss,
                                          *_r_avg_eta_ss,
                                          *_r_tageff_ss,
                                          *_r_tageff_asym_ss,
                                          *_r_production_asym,
                                          *_r_detection_asym );
        
        _cos_coeff = new DecRateCoeff_Bd("cos_coeff",
                                         "cos_coeff",
                                         DecRateCoeff_Bd::kCos ,
                                         *_r_f,
                                         *_r_C,
                                         *_r_C_bar,
                                         *_r_q_OS,
                                         *_r_eta_OS,
                                         *_r_p0_os,
                                         *_r_p1_os,
                                         *_r_delta_p0_os,
                                         *_r_delta_p1_os,
                                         *_r_avg_eta_os,
                                         *_r_tageff_os,
                                         *_r_tageff_asym_os,
                                         *_r_q_SS,
                                         *_r_eta_SS,
                                         *_r_p0_ss,
                                         *_r_p1_ss,
                                         *_r_delta_p0_ss,
                                         *_r_delta_p1_ss,
                                         *_r_avg_eta_ss,
                                         *_r_tageff_ss,
                                         *_r_tageff_asym_ss,
                                         *_r_production_asym,
                                         *_r_detection_asym );
        
        _sinh_coeff = new DecRateCoeff_Bd("sinh_coeff",
                                          "sinh_coeff",
                                          DecRateCoeff_Bd::kSinh ,
                                          *_r_f,
                                          *_r_D,
                                          *_r_D_bar,
                                          *_r_q_OS,
                                          *_r_eta_OS,
                                          *_r_p0_os,
                                          *_r_p1_os,
                                          *_r_delta_p0_os,
                                          *_r_delta_p1_os,
                                          *_r_avg_eta_os,
                                          *_r_tageff_os,
                                          *_r_tageff_asym_os,
                                          *_r_q_SS,
                                          *_r_eta_SS,
                                          *_r_p0_ss,
                                          *_r_p1_ss,
                                          *_r_delta_p0_ss,
                                          *_r_delta_p1_ss,
                                          *_r_avg_eta_ss,
                                          *_r_tageff_ss,
                                          *_r_tageff_asym_ss,
                                          *_r_production_asym,
                                          *_r_detection_asym );
        
        _sin_coeff = new DecRateCoeff_Bd("sin_coeff",
                                         "sin_coeff",
                                         DecRateCoeff_Bd::kSin,
                                         *_r_f,
                                         *_r_S,
                                         *_r_S_bar,
                                         *_r_q_OS,
                                         *_r_eta_OS,
                                         *_r_p0_os,
                                         *_r_p1_os,
                                         *_r_delta_p0_os,
                                         *_r_delta_p1_os,
                                         *_r_avg_eta_os,
                                         *_r_tageff_os,
                                         *_r_tageff_asym_os,
                                         *_r_q_SS,
                                         *_r_eta_SS,
                                         *_r_p0_ss,
                                         *_r_p1_ss,
                                         *_r_delta_p0_ss,
                                         *_r_delta_p1_ss,
                                         *_r_avg_eta_ss,
                                         *_r_tageff_ss,
                                         *_r_tageff_asym_ss,
                                         *_r_production_asym,
                                         *_r_detection_asym );
        
        /// Time pdf integrators
        _cosh_term = new TimePdfIntegrator(coshBasis,_efficiency,_efficiency2,_efficiency3,
                                           _Gamma,_dGamma,_dm
					   ,_offset_mean_dt,_scale_mean_dt,_scale_mean_2_dt
                                           ,_offset_sigma_dt, _scale_sigma_dt, _scale_sigma_2_dt
                                           ,_offset_sigma2_dt, _scale_sigma2_dt, _scale_sigma2_2_dt
                                           ,_offset_sigma3_dt, _scale_sigma3_dt, _scale_sigma3_2_dt
                                           ,_offset_f_dt, _scale_f_dt, _scale_f_2_dt
                                           ,_offset_f2_dt, _scale_f2_dt, _scale_f2_2_dt
                                           ,_c0, _c1, _c2
                                           ,_c3, _c4, _c5
                                           ,_c6, _c7, _c8
                                           ,_c9);
        _cos_term = new TimePdfIntegrator(cosBasis,_efficiency,_efficiency2,_efficiency3,
                                           _Gamma,_dGamma,_dm
					   ,_offset_mean_dt,_scale_mean_dt,_scale_mean_2_dt
                                           ,_offset_sigma_dt, _scale_sigma_dt, _scale_sigma_2_dt
                                           ,_offset_sigma2_dt, _scale_sigma2_dt, _scale_sigma2_2_dt
                                           ,_offset_sigma3_dt, _scale_sigma3_dt, _scale_sigma3_2_dt
                                           ,_offset_f_dt, _scale_f_dt, _scale_f_2_dt
                                           ,_offset_f2_dt, _scale_f2_dt, _scale_f2_2_dt
                                           ,_c0, _c1, _c2
                                           ,_c3, _c4, _c5
                                           ,_c6, _c7, _c8
                                           ,_c9);
        _sinh_term = new TimePdfIntegrator(sinhBasis,_efficiency,_efficiency2,_efficiency3,
                                           _Gamma,_dGamma,_dm
					   ,_offset_mean_dt,_scale_mean_dt,_scale_mean_2_dt
                                           ,_offset_sigma_dt, _scale_sigma_dt, _scale_sigma_2_dt
                                           ,_offset_sigma2_dt, _scale_sigma2_dt, _scale_sigma2_2_dt
                                           ,_offset_sigma3_dt, _scale_sigma3_dt, _scale_sigma3_2_dt
                                           ,_offset_f_dt, _scale_f_dt, _scale_f_2_dt
                                           ,_offset_f2_dt, _scale_f2_dt, _scale_f2_2_dt
                                           ,_c0, _c1, _c2
                                           ,_c3, _c4, _c5
                                           ,_c6, _c7, _c8
                                           ,_c9);
        _sin_term = new TimePdfIntegrator(sinBasis,_efficiency,_efficiency2,_efficiency3,
                                           _Gamma,_dGamma,_dm
					   ,_offset_mean_dt,_scale_mean_dt,_scale_mean_2_dt
                                           ,_offset_sigma_dt, _scale_sigma_dt, _scale_sigma_2_dt
                                           ,_offset_sigma2_dt, _scale_sigma2_dt, _scale_sigma2_2_dt
                                           ,_offset_sigma3_dt, _scale_sigma3_dt, _scale_sigma3_2_dt
                                           ,_offset_f_dt, _scale_f_dt, _scale_f_2_dt
                                           ,_offset_f2_dt, _scale_f2_dt, _scale_f2_2_dt
                                           ,_c0, _c1, _c2
                                           ,_c3, _c4, _c5
                                           ,_c6, _c7, _c8
                                           ,_c9);
        
        /// Marginal pdfs        
	cout << (string)_marginalPdfsPrefix << endl;
	if((string)_marginalPdfsPrefix == "Uniform"){
		_pdf_sigma_t = (RooAbsPdf*) (new RooUniform("pdf_sigma_t","pdf_sigma_t",*_r_dt));
		_pdf_eta_OS = (RooAbsPdf*) (new RooUniform("pdf_eta_OS","pdf_eta_OS",*_r_eta_OS));
		_pdf_eta_SS = (RooAbsPdf*) (new RooUniform("pdf_eta_SS","pdf_eta_SS",*_r_eta_SS));
	}
	else {
		_f_pdfs = new TFile("Mistag_pdfs.root","OPEN");
	
		_h_dt = new TH1D( *((TH1D*) _f_pdfs->Get(("h_dt_norm_"+(string)_marginalPdfsPrefix).c_str())));
		// RooHistPdf doesn't like negative or 0 bins (can happen due to sweights), set them to a small positive number
		_h_dt->Smooth();
		for(int i= 1 ; i<=_h_dt->GetNbinsX(); i++){
			if(_h_dt->GetBinContent(i) <= 0.)_h_dt->SetBinContent(i,0.000000001*_h_dt->GetMaximum());
		}
		_r_h_dt = new RooDataHist("r_h_dt","r_h_dt",*_r_dt,_h_dt);
		_pdf_sigma_t = (RooAbsPdf*) (new RooHistPdf("pdf_sigma_t","pdf_sigma_t",*_r_dt,*_r_h_dt));
		
		_h_eta_OS = new TH1D( *((TH1D*) _f_pdfs->Get(("h_w_OS_norm_"+(string)_marginalPdfsPrefix).c_str())));
		_h_eta_OS->Smooth();
		for(int i= 1 ; i<=_h_eta_OS->GetNbinsX(); i++){
			if(_h_eta_OS->GetBinContent(i) <= 0.)_h_eta_OS->SetBinContent(i,0.000000001*_h_eta_OS->GetMaximum());
		}
		_r_h_eta_OS = new RooDataHist("r_eta_OS","r_eta_OS",*_r_eta_OS,_h_eta_OS);
		_pdf_eta_OS = (RooAbsPdf*) (new RooHistPdf("pdf_eta_OS","pdf_eta_OS",*_r_eta_OS,*_r_h_eta_OS));
		
		_h_eta_SS = new TH1D( *((TH1D*) _f_pdfs->Get(("h_w_SS_norm_"+(string)_marginalPdfsPrefix).c_str())));
		_h_eta_SS->Smooth();
		for(int i= 1 ; i<=_h_eta_SS->GetNbinsX(); i++){
			if(_h_eta_SS->GetBinContent(i) <= 0.)_h_eta_SS->SetBinContent(i,0.000000001*_h_eta_SS->GetMaximum());
		}
		_r_h_eta_SS = new RooDataHist("r_eta_SS","r_eta_SS",*_r_eta_SS,_h_eta_SS);
		_pdf_eta_SS = (RooAbsPdf*) (new RooHistPdf("pdf_eta_SS","pdf_eta_SS",*_r_eta_SS,*_r_h_eta_SS));

// 		_h_q_f = new TH1D( *((TH1D*) _f_pdfs->Get(("h_q_f_norm_"+(string)_marginalPdfsPrefix).c_str())));
// 		_h_q_OS = new TH1D( *((TH1D*) _f_pdfs->Get(("h_q_OS_norm_"+(string)_marginalPdfsPrefix).c_str())));
// 		_h_q_SS = new TH1D( *((TH1D*) _f_pdfs->Get(("h_q_SS_norm_"+(string)_marginalPdfsPrefix).c_str())));
    	}

	_resmodel = new RooGaussModel("resmodel", "resmodel", *_r_t,  RooRealConstant::value(0.), *_r_dt_scaled);              

	_samplingPdf_t = new RooBDecay(("samplingPdf_t"+(string)_marginalPdfsPrefix).c_str(), ("samplingPdf_t"+(string)_marginalPdfsPrefix).c_str(),
				*_r_t,*_r_tau, *_r_dGamma, 
                                *_cosh_coeff,*_sinh_coeff,*_cos_coeff,*_sin_coeff,
				*_r_dm, *_resmodel, RooBDecay::SingleSided); 

   	_samplingPdf_t_eff = new RooEffProd(("samplingPdf_t_eff"+(string)_marginalPdfsPrefix).c_str(), ("samplingPdf_t_eff"+(string)_marginalPdfsPrefix).c_str(),
						 *_samplingPdf_t, *_splinePdf);
	
        _samplingPdf = new RooProdPdf(("samplingPdf"+(string)_marginalPdfsPrefix).c_str(),("samplingPdf"+(string)_marginalPdfsPrefix).c_str(),
				RooArgSet(*_pdf_sigma_t,*_pdf_eta_OS,*_pdf_eta_SS),Conditional(RooArgSet(*_samplingPdf_t_eff),RooArgSet(*_r_t,*_r_f,*_r_q_OS,*_r_q_SS)));


	_fitPdf_t = new RooBDecay(("fitPdf_t"+(string)_marginalPdfsPrefix).c_str(), ("fitPdf_t"+(string)_marginalPdfsPrefix).c_str(),
				*_r_t,*_r_tau, *_r_dGamma, 
                                *_cosh_coeff,*_sinh_coeff,*_cos_coeff,*_sin_coeff,
				*_r_dm, *_efficiency, RooBDecay::SingleSided); 

        _fitPdf = new RooProdPdf(("fitPdf"+(string)_marginalPdfsPrefix).c_str(),("fitPdf"+(string)_marginalPdfsPrefix).c_str(),
				RooArgSet(*_pdf_sigma_t,*_pdf_eta_OS,*_pdf_eta_SS),Conditional(RooArgSet(*_fitPdf_t),RooArgSet(*_r_t,*_r_f,*_r_q_OS,*_r_q_SS)));

	_sampleGen = 0;
    }
    
    RooDataSet* sampleEvents(int N = 10000){
        if(_sampleGen == 0){
    		_protoData = new RooDataSet("protoData","protoData",RooArgSet(*_r_dt,*_r_eta_OS,*_r_eta_SS));
		int N_proto = N > 100000 ? N : 100000;
		for(int i = 0 ; i < N; i++){
 		        _r_dt->setVal(_h_dt->GetRandom());        
         		_r_eta_OS->setVal(_h_eta_OS->GetRandom());
         		_r_eta_SS->setVal(_h_eta_SS->GetRandom());
		        _protoData->add(RooArgSet(*_r_dt,*_r_eta_OS,*_r_eta_SS));
		}
		_sampleGen = new RooMCStudy(*_samplingPdf,RooArgSet(*_r_t,*_r_q_OS,*_r_q_SS,*_r_f),ProtoData(*_protoData,kFALSE,kTRUE));
	}
	_sampleGen->generate(1,N,kTRUE);
	RooDataSet* data = (RooDataSet*)_sampleGen->genData(0);
	return data;
    }

    vector<double> getRandom_marginalVals(){
	vector<double> vals;
	vals.push_back(_h_dt->GetRandom());
	vals.push_back(_h_eta_OS->GetRandom());
	vals.push_back(_h_eta_SS->GetRandom());
	return vals;
    }

    void fillProtoData(double dt, int f, int q_OS, double eta_OS, int q_SS, double eta_SS){
        _r_dt->setVal(dt);        
        _r_q_SS->setIndex(q_SS);
        _r_q_OS->setIndex(q_OS);
        _r_f->setIndex(f);
        _r_eta_OS->setVal(eta_OS);
        _r_eta_SS->setVal(eta_SS);
        _protoData->add(RooArgSet(*_r_dt,*_r_f,*_r_q_OS,*_r_eta_OS,*_r_q_SS,*_r_eta_SS));
    }

    double getSamplingPdfVal(IDalitzEvent& evt){
	setAllObservablesAndFitParameters(evt);
	return _fitPdf->getVal(RooArgSet(*_r_t,*_r_dt,*_r_eta_OS,*_r_eta_SS,*_r_f,*_r_q_OS,*_r_q_SS));
    }

    double get_cosh_term_Val(IDalitzEvent& evt){
        return _cosh_coeff->evaluate() * _cosh_term->getVal(evt).real();
    }

    double get_cosh_term_Integral(IDalitzEvent& evt){
        return _cosh_coeff->analyticalIntegral(2) * _cosh_term->getVal(evt).imag();
    }

    double get_sinh_term_Val(IDalitzEvent& evt){
        return _sinh_coeff->evaluate() * _sinh_term->getVal(evt).real();
    }

    double get_sinh_term_Integral(IDalitzEvent& evt){
        return _sinh_coeff->analyticalIntegral(2) * _sinh_term->getVal(evt).imag();
    }
    
    double get_cos_term_Val(IDalitzEvent& evt){
        return _cos_coeff->evaluate() * _cos_term->getVal(evt).real();
    }
    
    double get_cos_term_Integral(IDalitzEvent& evt){
        return _cos_coeff->analyticalIntegral(2) * _cos_term->getVal(evt).imag();
    }
    
    double get_sin_term_Val(IDalitzEvent& evt){
        return _sin_coeff->evaluate() * _sin_term->getVal(evt).real();
    }
    
    double get_sin_term_Integral(IDalitzEvent& evt){
        return _sin_coeff->analyticalIntegral(2) * _sin_term->getVal(evt).imag();
    }

    double get_cosh_coeff_Val(IDalitzEvent& evt){
        return _cosh_coeff->evaluate();
    }

    double get_cos_coeff_Val(IDalitzEvent& evt){
        return _cos_coeff->evaluate();
    }

    double get_sinh_coeff_Val(IDalitzEvent& evt){
        return _sinh_coeff->evaluate();
    }

    double get_sin_coeff_Val(IDalitzEvent& evt){
        return _sin_coeff->evaluate();
    }
    
    double get_marginalPdfs_Val(IDalitzEvent& evt){
        const int q_OS = (int)evt.getValueFromVector(3);
        const int q_SS = (int)evt.getValueFromVector(5);
        
        /*
        // check norm 
        TRandom3 r(0);
        double sum1 = 0.;
        double sum2 = 0.;
        double sum3 = 0.;

        for (int i=0; i < 1000000; i++) {
            _r_eta_OS->setVal(r.Uniform(0,0.5));
            _r_eta_SS->setVal(r.Uniform(0,0.5));
            _r_dt->setVal(r.Uniform(0,0.1));
            sum1 += _pdf_eta_OS->getVal(RooArgSet(*_r_eta_OS));
            sum2 += _pdf_eta_SS->getVal(RooArgSet(*_r_eta_SS));
            sum3 += _pdf_sigma_t->getVal(RooArgSet(*_r_dt));
        }
        cout << "norm = " << 0.5 * sum1/1000000. << endl;
        cout << "norm = " << 0.5 * sum2/1000000. << endl;
        cout << "norm = " << 0.1 * sum3/1000000. << endl;
        throw "";
        */
        
        return  _pdf_sigma_t->getVal(RooArgSet(*_r_dt))
                * ( abs(q_OS) * _pdf_eta_OS->getVal(RooArgSet(*_r_eta_OS)) + ( 1. - abs(q_OS)) * 1./0.5 )
                * ( abs(q_SS) * _pdf_eta_SS->getVal(RooArgSet(*_r_eta_SS)) + ( 1. - abs(q_SS)) * 1./0.5 );
    }

    double get_marginalPdfs_product(IDalitzEvent& evt){
        return  _pdf_sigma_t->getVal(RooArgSet(*_r_dt)) * _pdf_eta_OS->getVal(RooArgSet(*_r_eta_OS)) * _pdf_eta_SS->getVal(RooArgSet(*_r_eta_SS));
    }

    double get_spline_Val(IDalitzEvent& evt){
        return  _spline->getVal();
    }
    
    double get_tau_Val(){
        return (double) 1./_Gamma;
    }
    
    double get_dm_Val(){
        return (double)_dm;
    }
    
    double get_dGamma_Val(){
        return (double)_dGamma;
    }
    
    std::pair<double, double> getCalibratedMistag_OS(IDalitzEvent& evt){
        return _cosh_coeff->calibrate(evt.getValueFromVector(4), _avg_eta_os, _p0_os, _p1_os, _delta_p0_os, _delta_p1_os);
    }

    std::pair<double, double> getCalibratedMistag_OS(double& eta_OS){
        return _cosh_coeff->calibrate(eta_OS, _avg_eta_os, _p0_os, _p1_os, _delta_p0_os, _delta_p1_os);
    }
    
    std::pair<double, double> getCalibratedMistag_OS(IDalitzEvent& evt,double& avg_eta_os,double& p0_os,double& p1_os,double& delta_p0_os,double& delta_p1_os ){
        return _cosh_coeff->calibrate(evt.getValueFromVector(4), avg_eta_os, p0_os, p1_os, delta_p0_os, delta_p1_os);
    }
    
    std::pair<double, double> getCalibratedMistag_SS(IDalitzEvent& evt){
        return _cosh_coeff->calibrate(evt.getValueFromVector(6), _avg_eta_ss, _p0_ss, _p1_ss, _delta_p0_ss, _delta_p1_ss);
    }

    std::pair<double, double> getCalibratedMistag_SS(IDalitzEvent& evt,double& avg_eta_ss,double& p0_ss,double& p1_ss,double& delta_p0_ss,double& delta_p1_ss ){
        return _cosh_coeff->calibrate(evt.getValueFromVector(6), avg_eta_ss, p0_ss, p1_ss, delta_p0_ss, delta_p1_ss);
    }

    std::pair<double, double> getCalibratedMistag(double eta,double avg_eta,double p0,double p1,double delta_p0,double delta_p1 ){
        return _cosh_coeff->calibrate(eta, avg_eta, p0, p1, delta_p0, delta_p1);
    }

    std::pair<double, double> getCalibratedMistag_SS(double& eta_SS){
        return _cosh_coeff->calibrate(eta_SS, _avg_eta_ss, _p0_ss, _p1_ss, _delta_p0_ss, _delta_p1_ss);
    }
        
    double getCalibratedResolution(double& dt){
        return _offset_sigma_dt + _scale_sigma_dt * dt + _scale_sigma_2_dt * dt *dt;
    }

    double getCalibratedResolution(double& dt,double& scale_sigma_dt,double& scale_sigma_2_dt){
        return _offset_sigma_dt + scale_sigma_dt * dt + scale_sigma_2_dt * dt *dt;
    }
    
   void listFitParDependencies(){
       std::cout << "TimePDfIntegrator depends on these fitParameters:" << std::endl;
       _cosh_term->listFitParDependencies(std::cout);
       //std::cout << "TimePDfMaster depends on these fitParameters (sinh):" << std::endl;
       //_sinh_term->listFitParDependencies(std::cout);
       //std::cout << "" << std::endl;
   }
   
   void setAllFitParameters(){       
       _r_tau->setVal(1./_Gamma);
       _r_dGamma->setVal(_dGamma);
       _r_dm->setVal(_dm);
       
       _r_offset_mean_dt->setVal(_offset_mean_dt);
       _r_scale_mean_dt->setVal(_scale_mean_dt);
       _r_scale_mean_2_dt->setVal(_scale_mean_2_dt);

       _r_offset_sigma_dt->setVal(_offset_sigma_dt);
       _r_scale_sigma_dt->setVal(_scale_sigma_dt);
       _r_scale_sigma_2_dt->setVal(_scale_sigma_2_dt);

       _r_offset_sigma2_dt->setVal(_offset_sigma2_dt);
       _r_scale_sigma2_dt->setVal(_scale_sigma2_dt);
       _r_scale_sigma2_2_dt->setVal(_scale_sigma2_2_dt);

       _r_offset_sigma3_dt->setVal(_offset_sigma3_dt);
       _r_scale_sigma3_dt->setVal(_scale_sigma3_dt);
       _r_scale_sigma3_2_dt->setVal(_scale_sigma3_2_dt);

       _r_offset_f_dt->setVal(_offset_f_dt);
       _r_scale_f_dt->setVal(_scale_f_dt);
       _r_scale_f_2_dt->setVal(_scale_f_2_dt);

       _r_offset_f2_dt->setVal(_offset_f2_dt);
       _r_scale_f2_dt->setVal(_scale_f2_dt);
       _r_scale_f2_2_dt->setVal(_scale_f2_2_dt);

       _r_c0->setVal(_c0);
       _r_c1->setVal(_c1);
       _r_c2->setVal(_c2);
       _r_c3->setVal(_c3);
       _r_c4->setVal(_c4);
       _r_c5->setVal(_c5);
       _r_c6->setVal(_c6);
       _r_c7->setVal(_c7);
       _r_c8->setVal(_c8);
       _r_c9->setVal(_c9);
       
       _r_p0_os->setVal(_p0_os);
       _r_p1_os->setVal(_p1_os);
       _r_delta_p0_os->setVal(_delta_p0_os);
       _r_delta_p1_os->setVal(_delta_p1_os);
       _r_avg_eta_os->setVal(_avg_eta_os);
       _r_tageff_os->setVal(_tageff_os);
       _r_tageff_asym_os->setVal(_tageff_asym_os);
       _r_p0_ss->setVal(_p0_ss);
       _r_p1_ss->setVal(_p1_ss);
       _r_delta_p0_ss->setVal(_delta_p0_ss);
       _r_delta_p1_ss->setVal(_delta_p1_ss);
       _r_avg_eta_ss->setVal(_avg_eta_ss);
       _r_tageff_ss->setVal(_tageff_ss);
       _r_tageff_asym_ss->setVal(_tageff_asym_ss);
       
       _r_production_asym->setVal(_production_asym);
       _r_detection_asym->setVal(_detection_asym);
   }
   
   void setAllObservables(IDalitzEvent& evt){
       const double t = (double) evt.getValueFromVector(0);
       const double dt = (double) evt.getValueFromVector(1);
       const int f = (int)evt.getValueFromVector(2);
       const int q_OS = (int)evt.getValueFromVector(3);
       const double eta_OS = (double) evt.getValueFromVector(4);
       const int q_SS = (int)evt.getValueFromVector(5);
       const double eta_SS = (double) evt.getValueFromVector(6);
       
       _r_t->setVal(t);
       _r_dt->setVal(dt);
       _r_f->setIndex(f);
       _r_q_OS->setIndex(q_OS);
       _r_eta_OS->setVal(eta_OS);
       _r_q_SS->setIndex(q_SS);
       _r_eta_SS->setVal(eta_SS);
    }
    
   void setAllObservablesToMean(IDalitzEvent& evt){
       _r_t->setVal(1./_Gamma);
       if((string)_marginalPdfsPrefix == "Uniform")_r_dt->setVal(1./_Gamma/100.);
       else _r_dt->setVal(_h_dt->GetMean());
       _r_f->setIndex((int)evt.getValueFromVector(2));
       _r_q_OS->setIndex(0);
       if((string)_marginalPdfsPrefix == "Uniform")_r_eta_OS->setVal(0.);
       else _r_eta_OS->setVal(_h_eta_OS->GetMean());
       _r_q_SS->setIndex(0);
       if((string)_marginalPdfsPrefix == "Uniform")_r_eta_SS->setVal(0.);
       else _r_eta_SS->setVal(_h_eta_SS->GetMean());
    }

    void setCP_coeff(double norm, double norm_bar,double C,double C_bar,double D,double D_bar,double S,double S_bar ){
        _r_norm->setVal(norm);
        _r_norm_bar->setVal(norm_bar);
        _r_C->setVal(C);
        _r_C_bar->setVal(C_bar);
        _r_D->setVal(D);
        _r_D_bar->setVal(D_bar);
        _r_S->setVal(S);
        _r_S_bar->setVal(S_bar);
    }

    vector<double> getCP_coeff(){

	vector<double> vals;
        vals.push_back(_r_norm->getVal());
        vals.push_back(_r_norm_bar->getVal());
        vals.push_back(_r_C->getVal());
        vals.push_back(_r_C_bar->getVal());
        vals.push_back(_r_D->getVal());
        vals.push_back(_r_D_bar->getVal());
        vals.push_back(_r_S->getVal());
        vals.push_back(_r_S_bar->getVal());
	
	return vals;
    }
    
    void setAllObservablesAndFitParameters(IDalitzEvent& evt){
        setAllObservables(evt);
        setAllFitParameters();
    }

   TH1D* plotSpline(){
       setAllFitParameters();
       TH1D *h_spline = new TH1D("", ";t (ps);Efficiency (a.u.)", 100, _min_TAU, _max_TAU);
       for (int i = 1; i<=h_spline->GetNbinsX(); i++) {
           _r_t->setVal(h_spline->GetXaxis()->GetBinCenter(i));
           h_spline->SetBinContent(i,_spline->getVal());
       }
       TCanvas* c = new TCanvas();
       h_spline->SetLineColor(kRed);
       h_spline->SetLineWidth(5);

       h_spline->SetMinimum(0.);
       h_spline->SetMaximum(1.2);

       h_spline->Draw("histc");
       c->Print("spline.eps");
       return h_spline;
   }

   virtual ~TimePdfMaster(){
       // Plot acceptance
       TH1D* h_spline = plotSpline();       
       TCanvas* c = new TCanvas();
       h_spline->SetLineColor(kRed);
       h_spline->Draw("histc");
       c->Print("spline.eps");
       //c->Print("spline.pdf");
   }

};

#endif
//
