#ifndef FLEXI_FAST_AMPLITUDE_INTEGRATOR_HH
#define FLEXI_FAST_AMPLITUDE_INTEGRATOR_HH
// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:02 GMT

#include "Mint/IDalitzIntegrator.h"

#include "Mint/DalitzEventPattern.h"

#include "Mint/DalitzEvent.h"
#include "Mint/DalitzEventList.h"

#include "TRandom.h"

#include "Mint/IEventGenerator.h"

#include "Mint/IFastAmplitudeIntegrable.h"
#include "Mint/IIntegrationCalculator.h"
#include "Mint/IntegCalculator.h"

#include "Mint/DalitzHistoSet.h"
#include "Mint/counted_ptr.h"

#include "Mint/DiskResidentEventList.h"

class FitAmpSum;
class IDalitzEvent;

namespace MINT{
  class Minimiser;
}
class FlexiFastAmplitudeIntegrator : virtual public IDalitzIntegrator{
  mutable double _mean;
  const static long int _minEvents=100;
  int _numEvents;
  int _Ncalls;
  int _NReEvaluations;

  double _sumT, _sumTRe, _sumTfast;
  bool _initialised;

  bool _db;
 protected:

  //DiskResidentEventList _fastFlexiEventList;
  DalitzEventList _fastFlexiEventList;

  DalitzEventPattern _pat;
  IFastAmplitudeIntegrable* _amps; // amplitudes
  //MINT::counted_ptr<IGetDalitzEvent> _weight; // for efficiencies etc
  //  FitAmpPairList _ampList;
  //MINT::counted_ptr<IIntegrationCalculator> _integCalc;
  //MINT::counted_ptr<IIntegrationCalculator> _integCalc_copyForDebug;
  //  MINT::counted_ptr<IntegCalculator> _integCalc;
  MINT::counted_ptr<FitAmpPairList> _integCalc;
  //  MINT::counted_ptr<IntegCalculator> _integCalc_copyForDebug;
  MINT::counted_ptr<FitAmpPairList> _integCalc_copyForDebug;

  TRandom* _rnd;
  MINT::counted_ptr<TRandom> _localRnd;
  double _precision;

  void forcedReIntegrate();
  void reIntegrate();

  double weight(IDalitzEvent* evtPtr);

  int updateEventSet(long int Nevents);
  double evaluateSum();
  int determineNumEvents();
  int generateEnoughEvents();

  void addEvent(IDalitzEvent& evt);
  //  void reAddEvent(DalitzEvent evt);
  void reAddEvent(DalitzEvent& evt);

  int addEvents(long int Nevents);

  TRandom* makeNewRnd(int seed=-9999);
  MINT::IEventGenerator<IDalitzEvent>* _generator;
 public:
  FlexiFastAmplitudeIntegrator(const std::string& fastFlexiEventFileName="fastFlexiEventFileName.root");

  FlexiFastAmplitudeIntegrator(const DalitzEventPattern& pattern
			       , IFastAmplitudeIntegrable* amps=0
			       , MINT::IEventGenerator<IDalitzEvent>* 
			       eventGenerator=0
			       , TRandom* rnd = 0
			       , double precision = 1.e-2
			       , const std::string& fastFlexiEventFileName="fastFlexiEventFileName.root"
			  );

  bool initialise(const DalitzEventPattern& pattern
		  , IFastAmplitudeIntegrable* amps=0
		  , MINT::IEventGenerator<IDalitzEvent>* eventGenerator=0
		  , TRandom* rnd = 0
		  , double precision = 1.e-2
		  );

  /* this is too risky until I save the fit parameter values with the file
  bool initialiseFromFile(const DalitzEventPattern& pattern
			  , IFastAmplitudeIntegrable* amps
			  , const std::string& commaSeparatedList
			  );
  */

  virtual bool add(const FlexiFastAmplitudeIntegrator& other);

  bool setValues(const DalitzEventPattern& pattern
		 , IFastAmplitudeIntegrable* amps=0
		 //, MINT::counted_ptr<IGetDalitzEvent> weight=0
		 , MINT::IEventGenerator<IDalitzEvent>* eventGenerator=0
		 , TRandom* rnd = 0
		 , double precision = 1.e-2
		 );

  void setPrecision(double prec){
    _precision=prec;
    if(initialised()){
      generateEnoughEvents();
    }
  }

  // warning - this will not re-evaluate the #events needed for
  // the precision etc.

  bool initialised() const{
    return _initialised;
  }

  double getVal();

  double RealVal(){
    return getVal();
  }
    
  double integralForMatchingPatterns(bool match,int pattern_sign){
        return _integCalc->integralForMatchingPatterns(match, pattern_sign);
  }
    
  std::complex<double> ComplexIntegralForTags(int tag1, int tag2){
        return _integCalc->ComplexIntegralForTags(tag1,tag2);
  }

  std::complex<double> ComplexSumForMatchingPatterns(bool match){
        return _integCalc->ComplexSumForMatchingPatterns(match);
  }

  double sumOfSqrtFitFractions() {
        return _integCalc->sumOfSqrtFitFractions();
  }
  
  double sumOfFitFractions() {
        return _integCalc->sumOfFitFractions();
  }  
    
  double absSumOfSqrtInterferenceFractions() {
        return _integCalc->absSumOfSqrtInterferenceFractions();
  }   
    
  double absSumOfInterferenceFractions() {
        return _integCalc->absSumOfInterferenceFractions();
  }   
    
  int numberOfFitFractionsLargerThanThreshold(double threshold){
        return _integCalc->numberOfFitFractionsLargerThanThreshold(threshold);
  } 

  FitFractionList getFractions() const{return _integCalc->getFractions();}
  FitFractionList getInterferenceTerms() const{return _integCalc->getInterferenceTerms();}
  bool doFractions() {return _integCalc->doFractions();}

  double variance() const;

  DalitzHistoSet histoSet() const;
  void saveEachAmpsHistograms(const std::string& prefix) const;
  std::vector<DalitzHistoSet> GetEachAmpsHistograms();

  DalitzHistoSet interferenceHistoSet() const;
  void saveInterferenceHistograms(const std::string& prefix) const;
  std::vector<DalitzHistoSet> GetInterferenceHistograms();

  virtual void doFinalStats(MINT::Minimiser* mini=0);
  void doFinalStatsAndSave(MINT::Minimiser* min=0,const std::string& fname = "FitAmpResults.txt", const std::string& fnameROOT="fitFractions.root");

  double getFractionChi2()const;

  bool ensureFreshEvents();



  //  bool save(const std::string& fname)const;

};
#endif
//
