#ifndef FIT_AMP_PAIR_LIST_HH
#define FIT_AMP_PAIR_LIST_HH
// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:03 GMT

#include "Mint/FitAmpPair.h"
#include "Mint/DalitzHistoSet.h"
#include "Mint/FitAmpPairCovariance.h"
#include "Mint/FitFractionList.h"

#include "Mint/IReturnRealForEvent.h"
#include "Mint/IIntegrationCalculator.h"

#include "Mint/IDalitzEvent.h"

#include "Mint/counted_ptr.h"
//#include <vector>

#include "Mint/PolymorphVector.h"

class FitAmplitude;
class IDalitzEvent;
namespace MINT{
  class Minimiser;
}

class FitAmpPairList 
: public MINT::PolymorphVector<FitAmpPair> // std::vector<FitAmpPair>
, virtual public IIntegrationCalculator
{

  MINT::NamedParameter<std::string> HistoOption;
  void applyHistoOption();

  int _Nevents;
  double _sum;
  double _sumsq;

  double _psSum;
  double _psSumSq;

  bool _slow;

  mutable FitAmpPairCovariance _cov;

  MINT::counted_ptr<MINT::IReturnRealForEvent<IDalitzEvent> > _efficiency;

  FitFractionList _singleAmpFractions, _interferenceFractions;

  double phaseSpaceIntegral()const; // only for debug

  std::string dirName() const;
  bool makeDirectory(const std::string& asSubdirOf=".")const;
  virtual double oldVariance() const;

  bool reset();
 public:
  FitAmpPairList();
  FitAmpPairList(const FitAmpPairList& other);
  virtual ~FitAmpPairList(){}

  MINT::counted_ptr<IIntegrationCalculator> 
    clone_IIntegrationCalculator() const;

  virtual void addAmps(FitAmplitude* a1, FitAmplitude* a2); // for backwards compatibility, will be removed
  virtual void addAmps(FitAmplitude& a1, FitAmplitude& a2);
  virtual void addEvent(IDalitzEvent* evtPtr, double weight=1);
  virtual void addEvent(IDalitzEvent& evt, double weight=1);
  virtual void addEvent(MINT::counted_ptr<IDalitzEvent> evtPtr
			, double weight=1);
  virtual void reAddEvent(IDalitzEvent& evt, double weight=1);
  //  virtual void reAddEvent(MINT::counted_ptr<IDalitzEvent> evtPtr
  //			, double weight=1);

  bool isCompatibleWith(const FitAmpPairList& other) const;
  virtual bool add(const FitAmpPairList& otherList);
  virtual bool add(const FitAmpPairList* otherListPtr);
  virtual bool add(MINT::const_counted_ptr<FitAmpPairList> otherListPtr);

  virtual bool append(const FitAmpPairList& otherListPtr);
  virtual bool append(const FitAmpPairList* otherListPtr);
  virtual bool append(MINT::const_counted_ptr<FitAmpPairList> otherListPtr);

  virtual int numEvents() const;
  virtual double integral() const;
  
  std::complex<double> ComplexIntegralForTags(int tag1, int tag2) const;

  double integralForMatchingPatterns(bool match,int pattern_sign) const;
  std::complex<double> ComplexSumForMatchingPatterns(bool match) const;

  virtual double variance() const;
  virtual double sumOfVariances() const;
  virtual std::complex<double> ComplexSum() const;
  virtual void Gradient(MINT::MinuitParameterSet* mps, std::vector<double>& grad);
  virtual void GradientForLasso(MINT::MinuitParameterSet* mps, std::vector<double>& grad);

  double sumOfSqrtFitFractions();
  double absSumOfInterferenceFractions();
  double absSumOfSqrtInterferenceFractions();
  double sumOfFitFractions();
  int numberOfFitFractionsLargerThanThreshold(double threshold);
    
  FitFractionList getFractions() const{return _singleAmpFractions;}
  FitFractionList getInterferenceTerms() const{return _interferenceFractions;}
  bool doFractions();

  void setEfficiency(MINT::counted_ptr<MINT::IReturnRealForEvent<IDalitzEvent> > eff);
  void unsetEfficiency();
  double efficiency(IDalitzEvent* evtPtr);
  bool haveEfficiency() const{return 0 != _efficiency;}

  virtual bool makeAndStoreFractions(MINT::Minimiser* mini=0){
    return makeAndStoreFractions("FitAmpResults.txt","fitFractions.root",mini);}
  virtual bool makeAndStoreFractions(const std::string& fname, const std::string& fnameROOT, MINT::Minimiser* min=0);
  virtual double getFractionChi2() const;
  
  virtual DalitzHistoSet histoSet() const;
  virtual DalitzHistoSet un_normalised_histoSetRe() const;
  virtual DalitzHistoSet un_normalised_histoSetIm() const;
  void saveEachAmpsHistograms(const std::string& prefix) const;
  std::vector<DalitzHistoSet> GetEachAmpsHistograms();

  DalitzHistoSet interferenceHistoSet() const;
  void saveInterferenceHistograms(const std::string& prefix) const;
  std::vector<DalitzHistoSet> GetInterferenceHistograms();

  virtual void doFinalStats(MINT::Minimiser* min=0);
  virtual void doFinalStatsAndSave(MINT::Minimiser* min=0,const std::string& fname = "FitAmpResults.txt", const std::string& fnameROOT="fitFractions.root"){
      makeAndStoreFractions(fname,fnameROOT,min);
  }

    
  virtual bool save(const std::string& asSubdirOf=".") const;
  virtual bool retrieve(const std::string& asSubdirOf=".") ;

  virtual void print(std::ostream& os=std::cout) const;

  bool needToReIntegrate() const;
  void startIntegration();
  void startReIntegration();
  void startForcedReIntegration();
  void endIntegration();


  void setSlow();
  void setFast();

  bool slow() const{return _slow;}
  bool fast() const{return !_slow;}

  FitAmpPairList& operator+=(const FitAmpPairList& other);
  FitAmpPairList operator+(const FitAmpPairList& other) const;
};

std::ostream& operator<<(std::ostream& os, const FitAmpPairList& fap);

#endif
//
