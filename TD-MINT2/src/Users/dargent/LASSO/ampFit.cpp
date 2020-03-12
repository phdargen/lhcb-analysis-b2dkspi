// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk), Jack Benton (Jack.B.Benton@bristol.ac.uk)
// status:  Fri 28 Jun 2013 11:21:01 GMT
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
#include "TGraph.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TNtupleD.h"
#include "TTree.h"
#include "TFile.h"
#include <TStyle.h>
#include <TChain.h>
#include "TRandom2.h"
#include "TRandom3.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TROOT.h>
#include "Mint/HyperHistogram.h"
//#include "Mint/GofTests.h"
//#include "Mint/PermutationTest.h"
#include "Mint/LASSO.h"
#include "Mint/LASSO_flexi.h"

using namespace std;
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

double getChi2(DalitzEventList& data, DalitzEventList& mc){
	
    double minBinWidth = 0.;
    const int dim = 5;
    
    NamedParameter<int> EventPattern("Event Pattern",  531, -431, 321, 211, -211);
    DalitzEventPattern pdg(EventPattern.getVector());
    //cout << " got event pattern: " << pdg << endl;
          
    NamedParameter<int> minEventsPerBin("minEventsPerBin", 25);       
    HyperPointSet points( dim );
    HyperPoint min(pdg.sijMin(1,3),pdg.sijMin(2,4),pdg.sijMin(3,4),pdg.sijMin(1,2,4),pdg.sijMin(2,3,4));
    HyperPoint max(pdg.sijMax(1,3),pdg.sijMax(2,4),pdg.sijMax(3,4),pdg.sijMax(1,2,4),pdg.sijMax(2,3,4));
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

    //hist.save("histData.root");
    //HyperBinningHistogram binningHist("histData.root",5);    
    //HyperBinningHistogram dataHist( binningHist.getBinning() );
    //dataHist.fill(points); 

    HyperPointSet pointsMC( dim);
    for (int i = 0; i < mc.size(); i++){
     	DalitzEvent evt = mc[i];
	HyperPoint point( dim);
      	point.at(0)= evt.s(1,3);
      	point.at(1)= evt.s(2,4); 
      	point.at(2)= evt.s(3,4);
      	point.at(3)= evt.sij(s124);
      	point.at(4)= evt.sij(s234);
      	point.addWeight(evt.getWeight());
      	pointsMC.push_back(point);
    }

    HyperHistogram mcHist( dataHist.getBinning() );
    mcHist.fill(pointsMC); 
    //data.normalise(1);
    mcHist.normalise(dataHist.integral());

    double chi2 = dataHist.chi2(mcHist);
    int nBins   = dataHist.getNBins();

    cout << "chi2 = " << (double)chi2/(nBins-1.) << endl;

    return (double)chi2/(nBins-1.);
}

void AddScaledAmpsToList(FitAmpSum& fas_tmp, FitAmpSum& fas, std::string name, counted_ptr<IReturnComplex>& scale){
    counted_ptr<FitAmpList> List = fas_tmp.GetCloneOfSubsetSameFitParameters(name);
    FitAmpSum fas_2(*List);
    //fasCC_2.CPConjugateSameFitParameters();
    //fasCC_2.CConjugateFinalStateSameFitParameters();
    fas.multiply(scale);
    fas.addAsList(fas_2,1.);
}

void AddAmpsToList(FitAmpSum& fas_tmp, FitAmpSum& fas, std::string name){
    counted_ptr<FitAmpList> List = fas_tmp.GetCloneOfSubsetSameFitParameters(name);
    FitAmpSum fas_2(*List);
    //fasCC_2.CPConjugateSameFitParameters();
    //fasCC_2.CConjugateFinalStateSameFitParameters();
    fas.addAsList(fas_2,1.);
}

class AmpRatio : virtual public IReturnComplex{
    FitParameter& _re;
    FitParameter& _im;
    int _f;
public:
    AmpRatio(FitParameter& re, FitParameter& im,int f = 1)
    : _re(re), _im(im), _f(f) {}
    
    complex<double> ComplexVal(){
        std::complex<double> result(_re,static_cast<double>(_f) * _im); 
        return result;
    }
};

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

int ampFit(int step=0){
    TRandom3 ranLux;
    NamedParameter<int> RandomSeed("RandomSeed", 0);
    ranLux.SetSeed((int)RandomSeed);
    gRandom = &ranLux;
    
    FitAmplitude::AutogenerateFitFile();
    NamedParameter<string> channel("channel", (std::string) "norm", (char*) 0);
    NamedParameter<string> InputDir("InputDir", (std::string) "/auto/data/dargent/BsDsKpipi/Final/", (char*) 0);

    NamedParameter<int> updateAnaNote("updateAnaNote", 0);
    NamedParameter<string> InputFileName("InputFileName", (std::string) "");
    NamedParameter<string> InputTreeName("InputTreeName", (std::string) "DalitzEventList");
    std::string inputFile = InputFileName;
    std::string inputTreeName = InputTreeName;
    bool generateNew = (std::string) InputFileName == "";
    std::cout << "InputFileName: " << InputFileName << std::endl;
    
    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root" , (char*) 0);
    NamedParameter<string> OutputRootFile("OutputRootFile", (std::string) "OutputRootFile.root" , (char*) 0);
    
    NamedParameter<int>  Nevents("Nevents", 1000);
    NamedParameter<double> integPrecision("IntegPrecision", 1.e-2);
    NamedParameter<std::string> integMethod("IntegMethod", (std::string)"efficient");
    NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);
    
    NamedParameter<int>  useLASSO("useLASSO", 1);
    NamedParameter<double>  lambda("lambda", 1.);
    NamedParameter<int>  doBootstrap("doBootstrap", 0);
    NamedParameter<int>  N_bootstrap("N_bootstrap", 10000);
    NamedParameter<int>  doPlots("doPlots", 1);
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
	
    FitAmpSum fas(pat);
    FitAmpIncoherentSum fasBkg(pat);
    {
        fas.print();    
    	DalitzEventList eventListNorm;
        TFile *file =  TFile::Open(((string)IntegratorEventFile).c_str());
        TTree* tree=dynamic_cast<TTree*>(file->Get("DalitzEventList"));
        eventListNorm.fromNtuple(tree,1);
        fas.normalizeAmps(eventListNorm);
        fasBkg.normalizeAmps(eventListNorm);
    }
    
    FitParameter sigfraction("SigFraction",2,1.,0.01);
    DalitzEventList eventList;
    
    if(generateNew){
        SignalGenerator sg(pat,&fas);
        cout << "Generating " << Nevents << " MC events." << endl;
        sg.FillEventList(eventList, Nevents);
        eventList.saveAsNtuple(OutputRootFile);
    }
    
    if(!generateNew){
        /// Load data
        double sw;
        double Ks[4];
        double pi[4];
        double D[4];
        
        TChain* tree;
        tree=new TChain("DecayTree");
        tree->Add(((string)InputFileName).c_str());
        tree->SetBranchStatus("*",0);
        tree->SetBranchStatus("N_B_sw",1);
        tree->SetBranchStatus("weight",1);
        tree->SetBranchStatus("FullDTF*",1);

        tree->SetBranchAddress("N_B_sw",&sw);
        
        tree->SetBranchAddress("FullDTF_Ks_PX",&Ks[0]);
        tree->SetBranchAddress("FullDTF_Ks_PY",&Ks[1]);
        tree->SetBranchAddress("FullDTF_Ks_PZ",&Ks[2]);
        tree->SetBranchAddress("FullDTF_Ks_PE",&Ks[3]);
        
        tree->SetBranchAddress("FullDTF_D_PX",&D[0]);
        tree->SetBranchAddress("FullDTF_D_PY",&D[1]);
        tree->SetBranchAddress("FullDTF_D_PZ",&D[2]);
        tree->SetBranchAddress("FullDTF_D_PE",&D[3]);

        tree->SetBranchAddress("FullDTF_pi_PX",&pi[0]);
        tree->SetBranchAddress("FullDTF_pi_PY",&pi[1]);
        tree->SetBranchAddress("FullDTF_pi_PZ",&pi[2]);
        tree->SetBranchAddress("FullDTF_pi_PE",&pi[3]);

        int N_sample = tree->GetEntries();
        vector<int> b_indices;
        while( b_indices.size() < N_bootstrap )b_indices.push_back(TMath::Nint(ranLux.Uniform(0,N_sample)));
        sort(b_indices.begin(), b_indices.end());
        if(doBootstrap)N_sample = b_indices.size();

        TRandom3 rndm;
        int badEvents = 0;
        for(int i=0; i< N_sample; i++){
            if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << N_sample << endl;
            if(doBootstrap) tree->GetEntry(b_indices[i]);
            else tree->GetEntry(i);

            double sign = 1.;
            //if(f > 0) sign = -1.;
            TLorentzVector Ks_p(sign*Ks[0],sign*Ks[1],sign*Ks[2],Ks[3]);
            TLorentzVector D_p(sign*D[0],sign*D[1],sign*D[2],D[3]);
            TLorentzVector pi_p(sign*pi[0],sign*pi[1],sign*pi[2],pi[3]);
            TLorentzVector B_p = Ks_p + pi_p + D_p;
            // array of vectors
            vector<TLorentzVector> vectorOfvectors;        
            vectorOfvectors.push_back(B_p*MeV);
            vectorOfvectors.push_back(D_p*MeV);
            vectorOfvectors.push_back(Ks_p*MeV);
            vectorOfvectors.push_back(pi_p*MeV);
        
            DalitzEvent evt = DalitzEvent(pat, vectorOfvectors);
            if(!(evt.phaseSpace() > 0.)){
                badEvents++;
                //continue;
            }
            
            if(abs(sqrt(evt.s(2,3))-1968.30)<30)continue;
            
            evt.setWeight(sw);
            eventList.Add(evt);	
        }
    cout << endl << "bad events " << badEvents << " ( " << badEvents/(double) N_sample * 100. << " %)" << endl << endl;
    }
        
    AmpsPdfFlexiFast* ampsSig;
    if(useLASSO) ampsSig = new AmpsPdfFlexiFast(pat, &fas, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);
    else ampsSig = new AmpsPdfFlexiFast(pat, &fas, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);
    AmpsPdfFlexiFast ampsBkg(pat, &fasBkg, 0, integPrecision,integMethod, (std::string) IntegratorEventFile);
    DalitzSumPdf amps(sigfraction,*ampsSig,ampsBkg);
    
    Neg2LL neg2ll(amps, eventList);
    
    double stepSize = 1;
    lambda = lambda + (step-1) * stepSize;
    LASSO_flexi lasso(ampsSig,lambda);
    Neg2LLSum fcn(&neg2ll,&lasso);
    fcn.addConstraints(); 
    
    Minimiser mini;
    if(useLASSO)mini.attachFunction(&fcn);
    else mini.attachFunction(&neg2ll);
    mini.doFit();    
    mini.printResultVsInput();
        
	if(useLASSO)ampsSig->doFinalStats(&mini);
	else {
            string outTableName = (string)OutputDir+"FitAmpResults.tex";
        	ampsSig->doFinalStatsAndSave(&mini,outTableName,((string)OutputDir+"FitAmpResults.root").c_str());
        }
                
    if(useLASSO){
        
        double N_sig = 0;
        for (int i=0; i<eventList.size(); i++) N_sig += eventList[i].getWeight();


        TFile* out_LASSO = new TFile(((string)OutputDir+"LASSO_"+anythingToString(step)+".root").c_str(),"RECREATE");
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
            y[0]=neg2ll.getVal() + 2. * lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]);
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
            y[0]=neg2ll.getVal() +  lasso.numberOfFitFractionsLargerThanThreshold(thresholds[i]) * log(N_sig);
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
            TGraph* r = new TGraph(1,x,y);
            r->SetName( ("r_"+anythingToString((int) (thresholds[i]*1000))).c_str());
            r->SetTitle("");
            r->GetXaxis()->SetTitle("#lambda");
            r->GetXaxis()->SetTitleOffset(0.65);
            r->GetYaxis()->SetTitle("Number of fit fractions larger than threshold");
            r->Draw("A*");
            r->Write();
        }
        
        y[0]=neg2ll.getVal() ;
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
        ff->SetName("SumOfFitFractions");
        ff->SetTitle("");
        ff->GetXaxis()->SetTitle("#lambda");
        ff->GetXaxis()->SetTitleOffset(0.65);
        ff->GetYaxis()->SetTitle("Total Fit Fraction");
        ff->Draw("A*");
        ff->Write();
        
        y[0]= lasso.absSumOfInterferenceFractions() ;
        TGraph* iff = new TGraph(1,x,y);
        iff->SetName("AbsSumOfInterferenceFractions");
        iff->SetTitle("");
        iff->GetXaxis()->SetTitle("#lambda");
        iff->GetXaxis()->SetTitleOffset(0.65);
        iff->GetYaxis()->SetTitle("Sum of abs(Interference Fraction)");
        iff->Draw("A*");
        iff->Write();
        
        /*
        y[0]= getChi2( eventList, eventListMC) ;
        out_LASSO->cd();
        TGraph* chi2 = new TGraph(1,x,y);
        chi2->SetName("Chi2");
        chi2->SetTitle("");
        chi2->GetXaxis()->SetTitle("#lambda");
        chi2->GetXaxis()->SetTitleOffset(0.65);
        chi2->GetYaxis()->SetTitle("#Chi^{2}/Bin");
        chi2->Draw("A*");
        chi2->Write();
        
        for (int i = 0; i < eventListMC.size(); i++){
            //DalitzEvent evt = eventListMC[i];
            double ampVal=amps.getVal_noPs(eventListMC[i]);
            eventListMC[i].setWeight(ampVal*eventListMC[i].getWeight()/eventListMC[i].getGeneratorPdfRelativeToPhaseSpace());
            //weightedMC.Add(evt);
        }
        //weightedMC.save(((string)OutputDir+"mcEvents.root").c_str());
        //DalitzHistoSet mcH = eventListMC.weightedHistoSet();
        //datH.drawWithFitNorm(mcH, ((string)OutputDir+"mcFit_l_"+anythingToString(step)+"_").c_str(),"eps");
        
        y[0]= getChi2_pionSwitched( eventList, eventListMC) ;
        out_LASSO->cd();
        TGraph* chi2_pionSwitched = new TGraph(1,x,y);
        chi2_pionSwitched->SetName("Chi2_pionSwitched");
        chi2_pionSwitched->SetTitle("");
        chi2_pionSwitched->GetXaxis()->SetTitle("#lambda");
        chi2_pionSwitched->GetXaxis()->SetTitleOffset(0.65);
        chi2_pionSwitched->GetYaxis()->SetTitle("#Chi^{2}/Bin (pion switched)");
        chi2_pionSwitched->Draw("A*");
        chi2_pionSwitched->Write();
        */
        
        out_LASSO->Write();
        out_LASSO->Close(); 
    }
        
    if(doPlots){
        
            cout << "Now plotting:" << endl;
            double nBins = 60;
            TH1D* m_DKs = new TH1D("m_DKs","; m(DK_{s}) [GeV]; Yield",nBins,2,5.5);
            TH1D* m_Dpi = new TH1D("m_Dpi","; m(D#pi) [GeV]; Yield",nBins,1,5.5);
            TH1D* m_Kspi = new TH1D("m_Kspi","; m(K_{s}#pi) [GeV]; Yield",nBins,0,4);
            
            TH2D* m_Kspi_DKs = new TH2D("m_Kspi_DKs","; m(K_{s}#pi) [GeV]; m(DK_{s}) [GeV];  Yield",nBins,0,4,nBins,2,5.5);
            TH2D* m_Kspi_Dpi = new TH2D("m_Kspi_Dpi","; m(K_{s}#pi) [GeV]; m(D#pi) [GeV];  Yield",nBins,0,4,nBins,2,5.5);
            TH2D* m_Dpi_DKs = new TH2D("m_Dpi_DKs","; m(D#pi) [GeV]; m(DK_{s}) [GeV];  Yield",nBins,2,5.5,nBins,2,5.5);
        
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
        
            for (int i=0; i<eventList.size(); i++) {
                m_DKs->Fill(sqrt(eventList[i].s(1,2)/(GeV*GeV)),eventList[i].getWeight());
                m_Dpi->Fill(sqrt(eventList[i].s(1,3)/(GeV*GeV)),eventList[i].getWeight());
                m_Kspi->Fill(sqrt(eventList[i].s(2,3)/(GeV*GeV)),eventList[i].getWeight());
            }    
            
            DalitzEventList eventListMC;
            TFile *FileMC =  TFile::Open(((string) IntegratorEventFile).c_str());
            TTree* treeMC = dynamic_cast<TTree*>(FileMC->Get("DalitzEventList"));
            eventListMC.fromNtuple(treeMC,1);
            FileMC->Close();

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

            vector<string> ampNames4;
            ampNames4.push_back("D(s2)(2573)");
            ampNames4.push_back("D(s1)(2700)");
            ampNames4.push_back("D(s1)*(2860)");
            ampNames4.push_back("D(s3)*(2860)");
            ampNames4.push_back("D(s2)(3040)-");

            vector<string> ampNames5;
            ampNames5.push_back("NonRes");
                       
            for(int i = 0; i < eventListMC.size(); i++){
                                    
                    double weight = ampsSig->un_normalised_noPs(eventListMC[i])*eventListMC[i].getWeight()/eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
     
                    double weight1 = fas.getAmpSqr(eventListMC[i], ampNames1) 
                        * eventListMC[i].getWeight()/ eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
                    double weight2 = fas.getAmpSqr(eventListMC[i], ampNames2) 
                        * eventListMC[i].getWeight()/ eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
                    double weight3 = fas.getAmpSqr(eventListMC[i], ampNames3) 
                        * eventListMC[i].getWeight()/ eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
                    double weight4 = fas.getAmpSqr(eventListMC[i], ampNames4) 
                        * eventListMC[i].getWeight()/ eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
                    double weight5 = fas.getAmpSqr(eventListMC[i], ampNames5) 
                        * eventListMC[i].getWeight()/ eventListMC[i].getGeneratorPdfRelativeToPhaseSpace();
              
                    m_DKs_fit->Fill(sqrt(eventListMC[i].s(1,2)/(GeV*GeV)),weight);
                    m_Dpi_fit->Fill(sqrt(eventListMC[i].s(1,3)/(GeV*GeV)),weight);
                    m_Kspi_fit->Fill(sqrt(eventListMC[i].s(2,3)/(GeV*GeV)),weight);

                    m_DKs_fit_1->Fill(sqrt(eventListMC[i].s(1,2)/(GeV*GeV)),weight1);
                    m_Dpi_fit_1->Fill(sqrt(eventListMC[i].s(1,3)/(GeV*GeV)),weight1);
                    m_Kspi_fit_1->Fill(sqrt(eventListMC[i].s(2,3)/(GeV*GeV)),weight1);

                    m_DKs_fit_2->Fill(sqrt(eventListMC[i].s(1,2)/(GeV*GeV)),weight2);
                    m_Dpi_fit_2->Fill(sqrt(eventListMC[i].s(1,3)/(GeV*GeV)),weight2);
                    m_Kspi_fit_2->Fill(sqrt(eventListMC[i].s(2,3)/(GeV*GeV)),weight2);

                    m_DKs_fit_3->Fill(sqrt(eventListMC[i].s(1,2)/(GeV*GeV)),weight3);
                    m_Dpi_fit_3->Fill(sqrt(eventListMC[i].s(1,3)/(GeV*GeV)),weight3);
                    m_Kspi_fit_3->Fill(sqrt(eventListMC[i].s(2,3)/(GeV*GeV)),weight3);

                    m_DKs_fit_4->Fill(sqrt(eventListMC[i].s(1,2)/(GeV*GeV)),weight4);
                    m_Dpi_fit_4->Fill(sqrt(eventListMC[i].s(1,3)/(GeV*GeV)),weight4);
                    m_Kspi_fit_4->Fill(sqrt(eventListMC[i].s(2,3)/(GeV*GeV)),weight4);

                    m_DKs_fit_5->Fill(sqrt(eventListMC[i].s(1,2)/(GeV*GeV)),weight5);
                    m_Dpi_fit_5->Fill(sqrt(eventListMC[i].s(1,3)/(GeV*GeV)),weight5);
                    m_Kspi_fit_5->Fill(sqrt(eventListMC[i].s(2,3)/(GeV*GeV)),weight5);

                    eventListMC[i].setWeight(weight);
            }

            TCanvas* c = new TCanvas();
                     
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
            m_DKs_fit_1->SetFillColor(kRed+1);
            m_DKs_fit_1->SetFillStyle(3353);
            m_DKs_fit_1->DrawNormalized("histcsame",m_DKs_fit_1->Integral()/m_DKs_fit->Integral());
            m_DKs_fit_2->SetLineColor(kGreen+3);
            m_DKs_fit_2->SetLineWidth(2);
            m_DKs_fit_2->SetFillColor(kGreen+3);
            m_DKs_fit_2->SetFillStyle(3353);
            m_DKs_fit_2->DrawNormalized("histcsame",m_DKs_fit_2->Integral()/m_DKs_fit->Integral());
            m_DKs_fit_3->SetLineColor(kMagenta+3);
            m_DKs_fit_3->SetLineWidth(2);
            m_DKs_fit_3->SetFillColor(kMagenta+3);
            m_DKs_fit_3->SetFillStyle(3353);
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
             m_Dpi_fit_1->SetFillColor(kRed+1);
             m_Dpi_fit_1->SetFillStyle(3353);
             m_Dpi_fit_1->DrawNormalized("histcsame",m_Dpi_fit_1->Integral()/m_Dpi_fit->Integral());
             m_Dpi_fit_2->SetLineColor(kGreen+3);
             m_Dpi_fit_2->SetLineWidth(2);
             m_Dpi_fit_2->SetFillColor(kGreen+3);
             m_Dpi_fit_2->SetFillStyle(3353);
             m_Dpi_fit_2->DrawNormalized("histcsame",m_Dpi_fit_2->Integral()/m_Dpi_fit->Integral());
             m_Dpi_fit_3->SetLineColor(kMagenta+3);
             m_Dpi_fit_3->SetLineWidth(2);
             m_Dpi_fit_3->SetFillColor(kMagenta+3);
             m_Dpi_fit_3->SetFillStyle(3353);
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

            m_Kspi->SetMinimum(0.01);
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
             m_Kspi_fit_1->SetFillColor(kRed+1);
             m_Kspi_fit_1->SetFillStyle(3353);
             m_Kspi_fit_1->DrawNormalized("histcsame",m_Kspi_fit_1->Integral()/m_Kspi_fit->Integral());
             m_Kspi_fit_2->SetLineColor(kGreen+3);
             m_Kspi_fit_2->SetLineWidth(2);
             m_Kspi_fit_2->SetFillColor(kGreen+3);
             m_Kspi_fit_2->SetFillStyle(3353);
             m_Kspi_fit_2->DrawNormalized("histcsame",m_Kspi_fit_2->Integral()/m_Kspi_fit->Integral());
             m_Kspi_fit_3->SetLineColor(kMagenta+3);
             m_Kspi_fit_3->SetLineWidth(2);
             m_Kspi_fit_3->SetFillColor(kMagenta+3);
             m_Kspi_fit_3->SetFillStyle(3353);
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
    
            //getChi2(eventList,eventListMC);
    }
    return 0;
}

void makeIntegratorFile(){
    
    FitAmplitude::AutogenerateFitFile();

    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int>  IntegratorEvents("IntegratorEvents", 300000);
    
    DalitzEventList eventListPhsp,eventList,eventList_cut;
    
    eventListPhsp.generatePhaseSpaceEvents(100000,pat);
    
    FitAmpIncoherentSum fas((DalitzEventPattern)pat);
    fas.print();
    fas.getVal(eventListPhsp[0]);
    fas.normalizeAmps(eventListPhsp);
    
    SignalGenerator sg(pat,&fas);    
    sg.FillEventList(eventList, IntegratorEvents);

    for(int i = 0; i < eventList.size(); i++){
            if(abs(sqrt(eventList[i].s(2,3))-1968.30)<30)continue;
            eventList_cut.Add(eventList[i]);
    }

    cout << "Generated " << eventList_cut.size() << " events inside selected phasespace region" << endl;
    
    eventList_cut.saveAsNtuple(IntegratorEventFile);
    return;
}

void makeIntegratorFilePhsp(){
    
    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int>  IntegratorEvents("IntegratorEvents", 300000);
    
    DalitzEventList eventList;

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);

    int i = 0;
    while(i < IntegratorEvents){
	DalitzEvent evt(pat);
	evt.generateThisToPhaseSpace();
	if(sqrt(evt.sij(s234)/(GeV*GeV)) < 1.95 && sqrt(evt.s(2,4)/(GeV*GeV)) < 1.2 && sqrt(evt.s(3,4)/(GeV*GeV)) < 1.2){ 
		if(evt.phaseSpace() > 0.){
			eventList.Add(evt);
			i++;
		}
	}
    }

    cout << "Generated " << eventList.size() << " events inside selected phasespace region" << endl;
    
    eventList.saveAsNtuple(IntegratorEventFile);
    return;
}

void makeIntegratorFileEvtGen(){
    
    NamedParameter<int> EventPattern("Event Pattern", 521, 321, 211, -211, 443);
    DalitzEventPattern pat(EventPattern.getVector());
    cout << " got event pattern: " << pat << endl;
    
    NamedParameter<int>  IntegratorEvents("IntegratorEvents", 300000);
    
    DalitzEventList eventList_cut;
    DiskResidentEventList eventList(pat,"/auto/data/dargent/BsDsKpipi/MINT/SignalIntegrationEvents_EvtGen_new.root","OPEN");

    vector<int> s234;
    s234.push_back(2);
    s234.push_back(3);
    s234.push_back(4);

    for(int i = 0; i < IntegratorEvents; i++){

        DalitzEvent evt = eventList.getEvent(i);

	if(sqrt(evt.sij(s234)/(GeV*GeV)) < 1.95 && sqrt(evt.s(2,4)/(GeV*GeV)) < 1.2 && sqrt(evt.s(3,4)/(GeV*GeV)) < 1.2)eventList_cut.Add(evt);
    }

    cout << "Generated " << eventList_cut.size() << " events inside selected phasespace region" << endl;
    
    eventList_cut.saveAsNtuple("SignalIntegrationEvents_MINT_PhspCut.root");
    return;
}

int main(int argc, char** argv){
    
    time_t startTime = time(0);
    
    TH1::SetDefaultSumw2();
    TH2::SetDefaultSumw2();
    gROOT->ProcessLine(".x ../lhcbStyle.C");
    gStyle->SetPalette(1);

    NamedParameter<string> IntegratorEventFile("IntegratorEventFile", (std::string) "SignalIntegrationEvents.root", (char*) 0);
    if(! std::ifstream(((string)IntegratorEventFile).c_str()).good()) makeIntegratorFile();
  
    ampFit(atoi(argv[1]));
    
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
    return 0;
}
//
