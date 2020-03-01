// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:03 GMT
#include "Mint/FitAmpPairList.h"

#include "Mint/IDalitzEvent.h"
#include "Mint/FitAmplitude.h"
#include "Mint/FitAmpPairFitCovariance.h"

#include "Mint/Minimiser.h"

#include "Mint/Utils.h"

#include <sys/types.h>
#include <sys/stat.h>

#include <algorithm>
#include <vector>

#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;
using namespace MINT;

//MINT::NamedParameter<std::string> FitAmpPairList::HistoOption("FitAmpPairList::HistoOption"
//							  , (std::string) "default");

/*
FitAmpPairList::HistoOption can be set to either fast (no
histograms will be created or retrieved by the integrator) or slow
(histograms will be created or retrieved). The default option is
"default" and the default behaviour with FastAmplitudeIntegrator used
by DalitzPdfFastInteg is "slow", i.e. with histograms (the behaviour
is a bit more complex for the case of the
FlexiFastAmplitudeIntegrator, you're better off not using this option
there at all). Omitting the histograms saves a bit of time and
memory. Note that if you use this option when generating sgIntegrator
files (not recommended), you will also need to use it when reading
them. Recommended use is to use it with DalitzPdfFastInteg, when
reading in pre-saved integrator directories ("sgIntegrator"), in cases
where you don't make any plots.
*/

void FitAmpPairList::applyHistoOption(){
  if((string)HistoOption == (string)"fast") setFast();
  if((string)HistoOption == (string)"slow") setSlow();
}

FitAmpPairList::FitAmpPairList()
  : MINT::PolymorphVector<FitAmpPair>()
  , _Nevents(0)
  , _sum(0)
  , _sumsq(0)
  , _psSum(0)
  , _psSumSq(0)
  , _slow(true)
  , _cov(this)
  , _efficiency(0)
  , HistoOption("FitAmpPairList::HistoOption", (std::string) "default")
{
  applyHistoOption();
}
FitAmpPairList::FitAmpPairList(const FitAmpPairList& other)
  : IIntegrationCalculator()
  , MINT::PolymorphVector<FitAmpPair>(other)
  , _Nevents(other._Nevents)
  , _sum(other._sum)
  , _sumsq(other._sumsq)
  , _psSum(other._psSum)
  , _psSumSq(other._psSumSq)
  , _slow(other._slow)
  , _cov(other._cov, this)
  , _efficiency(other._efficiency)
  , HistoOption("FitAmpPairList::HistoOption", (std::string) "default")
{
  applyHistoOption();
}
counted_ptr<IIntegrationCalculator> 
FitAmpPairList::clone_IIntegrationCalculator() const{
  counted_ptr<IIntegrationCalculator> ptr(new FitAmpPairList(*this));
  return ptr;
}

void FitAmpPairList::addAmps(FitAmplitude* a1, FitAmplitude* a2){
  if(0 == a1 || 0 == a2){
    cout << "Error in FitAmpPairList::addAmps: " << a1 << ", " << a2 << endl;
    throw "zero pointer";
  }
  FitAmpPair amps(*a1, *a2);
  this->push_back(amps);
}
void FitAmpPairList::addAmps(FitAmplitude& a1, FitAmplitude& a2){
  FitAmpPair amps(a1, a2);
  this->push_back(amps);
}

void FitAmpPairList::addEvent(IDalitzEvent* evtPtr, double weight){
  if(0 == evtPtr) return;
  addEvent(*evtPtr, weight);
}

void FitAmpPairList::addEvent(IDalitzEvent& evt, double weight){
  bool dbThis=false;
  static unsigned int nZeroPs=0;
  _Nevents++;
  double x=0;
  double ps = evt.phaseSpace();
  if(ps <= 0.0){
	  if (!(nZeroPs%1000))
	  {
		cout << "WARNING in FitAmpPairList::addToHistograms"
		 << " " << ++nZeroPs << " th"
		 << " event with phase space = " << ps << endl;
		return; // should not happen.
	  }
  }

  for(unsigned int i=0; i< this->size(); i++){
    //    x += this->at(i).add(evtPtr, weight*efficiency(evtPtr));
    x += this->at(i).add(evt, weight, efficiency(&evt));
  }
  if(dbThis) cout << "FitAmpPairList::addEvent(): adding event" << endl;
  //if(_cov.isValid()) 
  _cov.addLastEventFromList();

  _sum   += x;
  _sumsq += x*x;
  
  //if(0 == evtPtr) return;
  //double ps = evtPtr->phaseSpace();
  _psSum   += ps;
  _psSumSq += ps*ps;

  return;
}
void FitAmpPairList::addEvent(counted_ptr<IDalitzEvent> evtPtr, double weight){
  return addEvent(evtPtr.get(), weight);
}

void FitAmpPairList::reAddEvent(IDalitzEvent& evt, double weight){
  if(evt.phaseSpace() < 0) return;

  for(unsigned int i=0; i< this->size(); i++){
    if( this->at(i).acceptEvents()){
      this->at(i).reAdd(evt, weight, efficiency(&evt));
    }
  }

  //  _cov.addLastEventFromList(); (need to do something about it)

  return;
}
//void FitAmpPairList::reAddEvent(counted_ptr<IDalitzEvent> evtPtr, double weight){
//  return reAddEvent(evtPtr.get(), weight);
//}

bool FitAmpPairList::isCompatibleWith(const FitAmpPairList& otherList) const{
  if(this->size() != otherList.size()) return false;
  for(unsigned int i=0; i < this->size(); i++){
    if(this->at(i).name() != otherList[i].name()) return false;
  }
  return true;
}
bool FitAmpPairList::add(const FitAmpPairList* otherListPtr){
  if(0 != otherListPtr){
    return add(*otherListPtr);
  }
  return true;
}
bool FitAmpPairList::add(const_counted_ptr<FitAmpPairList> otherListPtr){
  if(0 != otherListPtr){
    return add(*otherListPtr);
  }
  return true;
}
bool FitAmpPairList::add(const FitAmpPairList& otherList){
  if(! isCompatibleWith(otherList)){
    cout << "ERROR in FitAmpPairList::add "
	 << "trying to add incompatible lists"
	 << endl;
    return false;
  }

  _Nevents  += otherList._Nevents;
  _sum      += otherList._sum;
  _sumsq    += otherList._sumsq;
  _psSum    += otherList._psSum;
  _psSumSq  += otherList._psSumSq;

  _cov      += otherList._cov;

  for(unsigned int i=0; i< this->size(); i++){
    this->at(i) += otherList[i];
  }
  
  return true;

}

bool FitAmpPairList::append(MINT::const_counted_ptr<FitAmpPairList> otherListPtr){
  return append(*otherListPtr);
}
bool FitAmpPairList::append(const FitAmpPairList* otherListPtr){
  return append(*otherListPtr);
}
bool FitAmpPairList::append(const FitAmpPairList& otherList){
  // WARNING: this will reset everything to zero.

  for(unsigned int i=0; i< otherList.size(); i++){
    this->push_back(otherList[i]);
  }
  this->reset();

  return true;

}
bool FitAmpPairList::reset(){
  _Nevents=0;
  _sum=0;
  _sumsq=0;
  _psSum=0;
  _psSumSq=0;
  _cov.reset();
  _singleAmpFractions.clear();
  _interferenceFractions.clear();

  return true;
}
int FitAmpPairList::numEvents() const{
  return _Nevents;
}

double FitAmpPairList::integral() const{
  double sum = 0;
  for(unsigned int i=0; i< this->size(); i++){
    sum += this->at(i).integral();
  }
  return sum;
}

std::complex<double> FitAmpPairList::ComplexIntegralForTags(int tag1, int tag2) const{
    std::complex<double> sum = 0;

    for(unsigned int i=0; i< this->size(); i++){
        if((this->at(i).fitAmp1().getTag() == tag1) && ( this->at(i).fitAmp2().getTag() == tag2)){
            sum += conj(this->at(i).complexVal());       
        }
        else if((this->at(i).fitAmp1().getTag() == tag2) && ( this->at(i).fitAmp2().getTag() == tag1)){
            sum += this->at(i).complexVal();
        }
    }
    return sum;
}

double FitAmpPairList::integralForMatchingPatterns(bool match, int pattern_sign) const{
    double sum = 0;
    for(unsigned int i=0; i< this->size(); i++){
        if(this->at(i).hasMatchingPattern()==match){
            if(match){
                if(static_cast<double>(pattern_sign) * this->at(i).fitAmp1().theBareDecay().getVal() < 0 ) continue;
            }
            sum += this->at(i).integral();
        }
    }
    return sum;
}

std::complex<double> FitAmpPairList::ComplexSumForMatchingPatterns(bool match) const{
    std::complex<double> sum = 0;
    for(unsigned int i=0; i< this->size(); i++){
        if(this->at(i).hasMatchingPattern()==match){
            std::complex<double> val = this->at(i).complexVal();
            if(this->at(i).fitAmp1().theBareDecay().getVal() >  this->at(i).fitAmp2().theBareDecay().getVal() ) val = conj(val);
            sum += val;
        }
    }
    return sum;
}

std::complex<double> FitAmpPairList::ComplexSum() const{ //laurenPsuedo
  bool dbThis = false;
  std::complex<double> sum = 0;
  for(unsigned int i=0; i< this->size(); i++){
    sum += this->at(i).complexVal();
  }
  return sum;
  (void)dbThis;
}

void FitAmpPairList::Gradient(MinuitParameterSet* mps, std::vector<double>& grad){
    
  for (unsigned int i=0; i<mps->size(); i++) {
    if(mps->getParPtr(i)->hidden())continue;
    
    if(i+1 >= grad.size()){
      cout << "WARNING in FitAmpPairList::Gradient:"
	   << " have to increase size of grad to avoid memory issues" << endl;
      grad.resize(i+2);
    }
 
    grad.at(i)=0; grad.at(i+1)=0;

    string name_i= mps->getParPtr(i)->name();
    if(name_i.find("_Re")!=std::string::npos){
      if(mps->getParPtr(i)->iFixInit() && mps->getParPtr(i+1)->iFixInit()){
	i++;
	continue;
      }
      name_i.replace(name_i.find("_Re"),3,"");
      complex<double> sum(0);
      for(unsigned int j=0; j< this->size(); j++){
	if(!A_is_in_B(name_i,this->at(j).name())) continue;
	if(A_is_in_B("Inco", name_i) != A_is_in_B("Inco",this->at(j).name()) ) continue;
	double singleAmpCorrection= 1.;
	if(this->at(j).isSingleAmp()) singleAmpCorrection = 2.;
	// 2 a_j^* A_i A_j^*
	if(A_is_in_B(name_i,this->at(j).fitAmp1().name())){
	  sum += singleAmpCorrection* this->at(j).valNoFitPars()* conj(this->at(j).fitAmp2().AmpPhase());
	}
	else {
	  sum += singleAmpCorrection* conj( this->at(j).valNoFitPars()* this->at(j).fitAmp1().AmpPhase() );
	}
        
      }
      grad.at(i)= sum.real();
      grad.at(i+1)= -sum.imag();
      i++;
    } 
    // Doesn't work. Don't use! 
    /* 
       else if(name_i.find("_Amp")!=std::string::npos){
       name_i.replace(name_i.find("_Amp"),4,""); 
       complex<double> sumAmp(0);
       complex<double> sumPhase(0);
       
       for(unsigned int j=0; j< this->size(); j++){
       if(!A_is_in_B(name_i,this->at(j).name())) continue;
       double singleAmpCorrection= 1.;
       if(this->at(j).isSingleAmp()) singleAmpCorrection = 2.;
       // 2 a_j^* A_i A_j^*
       if(A_is_in_B(name_i,this->at(j).fitAmp1().name())){
       sumAmp += singleAmpCorrection*this->at(j).complexVal()/std::abs(this->at(j).fitAmp1().AmpPhase());
       if(!this->at(j).isSingleAmp())sumPhase += this->at(j).complexVal()*std::arg(this->at(j).fitAmp1().AmpPhase());
       }
       else {
       sumAmp += singleAmpCorrection*conj(this->at(j).complexVal())/std::abs(this->at(j).fitAmp2().AmpPhase());
       if(!this->at(j).isSingleAmp())sumPhase += conj(this->at(j).complexVal())*std::arg(this->at(j).fitAmp2().AmpPhase());
       }
       
       }
       grad[i]= sumAmp.real();
       grad[i+1]= -sumPhase.imag();
       i++;
       }
    */
    else if(mps->getParPtr(i)->iFixInit())continue;
    else {
	  std::cout << "FitAmpPairList::Gradient() called. Sorry, I don't know how to calculate the derivative with respect to the fit parameter " << mps->getParPtr(i)->name() << " ! Please implement me or set useAnalytic Gradient to 0 in your options file. I'll crash now. " << std::endl;
            throw "crash";
        }
    }
    
}

void FitAmpPairList::GradientForLasso(MinuitParameterSet* mps
				      , std::vector<double>& grad){
    
    for (unsigned int i=0; i<mps->size(); i++) {
        if(mps->getParPtr(i)->hidden())continue;

        grad[i]=0; grad[i+1]=0;
        string name_i= mps->getParPtr(i)->name();
        if(name_i.find("_Re")!=std::string::npos){
            if(mps->getParPtr(i)->iFixInit() && mps->getParPtr(i+1)->iFixInit())continue;
            name_i.replace(name_i.find("_Re"),3,"");
            for(unsigned int j=0; j< this->size(); j++){
                if(!A_is_in_B(name_i,this->at(j).name())) continue;
                if(A_is_in_B("Inco", name_i) != A_is_in_B("Inco",this->at(j).name()) ) continue;
                if(!this->at(j).isSingleAmp()) continue;

		if(i+1 >= grad.size()){
		  cout << "WARNING in FitAmpPairList::GradientForLasso:"
		       << " have to increase size of grad to avoid memory issues" << endl;
		  grad.resize(i+2);
	}

                double integral = this->at(j).valNoFitPars().real()/sqrt(this->at(j).integral());
                // Re(a_j) A_j A_j^* dphi
                grad[i]= integral * this->at(j).fitAmp1().AmpPhase().real();
                // Im(a_j) A_j A_j^* dphi
                grad[i+1]= integral * this->at(j).fitAmp1().AmpPhase().imag();
                i++;
                break;
            }
        }
        else if(mps->getParPtr(i)->iFixInit())continue;
        else {
            std::cout << "Sorry, I don't know how to calculate the derivative with respect to the fit parameter " << mps->getParPtr(i)->name() << " ! Please implement me or set useAnalytic Gradient to 0 in your options file. I'll crash now. " << std::endl;
            throw "crash";
        }
    }
    
}

double FitAmpPairList::sumOfSqrtFitFractions() {
    double sum = 0;
    for(unsigned int i=0; i< this->size(); i++){
        if(this->at(i).isSingleAmp())sum += sqrt(this->at(i).integral());
    }
    return sum;
}

double FitAmpPairList::sumOfFitFractions() {
    double sum = 0;
    for(unsigned int i=0; i< this->size(); i++){
        if(this->at(i).isSingleAmp())sum += this->at(i).integral();
    }
    return sum/integral();
}

double FitAmpPairList::absSumOfSqrtInterferenceFractions() {
    double sum = 0;
    for(unsigned int i=0; i< this->size(); i++){
        if(!this->at(i).isSingleAmp())sum += sqrt(abs(this->at(i).integral()));
    }
    return sum;
}

double FitAmpPairList::absSumOfInterferenceFractions() {
    double sum = 0;
    for(unsigned int i=0; i< this->size(); i++){
        if(!this->at(i).isSingleAmp())sum += abs(this->at(i).integral());
    }
    return sum/integral();
}

int FitAmpPairList::numberOfFitFractionsLargerThanThreshold(double threshold){
    int n=0;
    double norm = integral();
    for(unsigned int i=0; i< this->size(); i++){
        if(this->at(i).isSingleAmp()){
            if(this->at(i).integral()/norm > threshold)n++ ;
        }
    }
    return n;
}


double FitAmpPairList::phaseSpaceIntegral()const{
  // only for debug
  if(_Nevents <= 0) return -9999;
  return _psSum/((double) _Nevents);
}

double FitAmpPairList::variance() const{
  bool dbThis=false;
  static int printouts=0;
  static bool printedWarning=false;
  if(_Nevents <=0)return 0;

  if(_cov.isValid()){
    if(dbThis && printouts < 100){
      printouts++;
      cout << "FitAmpPairList::variance() with "
	   << _Nevents << "(" << _cov.Nevents() << ") events. Compare: "
	   << " old: " << oldVariance()
	   << " sum: " << sumOfVariances()
	   << ", new: " << _cov.getIntegralVariance()
	   << endl;
    }
    if(_Nevents < 10000) return oldVariance();
    double var = _cov.getIntegralVariance();
    if(var < 1.e-6 * oldVariance()){
      if(dbThis) cout << "returning OLD variance " << var << endl;
      // (if var is tiny, it usually means something went wrong)
      return oldVariance(); 
    }
    if(dbThis) cout << "returning new variance " << var << endl;
    return var;
  }else{
    if(! printedWarning){
      cout << "WARNING in FitAmpPairList::variance()"
	   << " ignoring correlations because _cov is invalid."
	   << " (I'll print this only once)"
	   << endl;
      printedWarning=true;
    }
    return sumOfVariances();
  }
}

double FitAmpPairList::sumOfVariances() const{
  double sum=0;
  for(unsigned int i=0; i < this->size(); i++){    
    sum += this->at(i).variance();
  }
  return sum;
}
double FitAmpPairList::oldVariance() const{
  // this does not (yet) take into account
  // that the variance changes with the
  // fit paramters.
  bool dbThis=false;
  if(0 == _Nevents) return -9999.0;

  double dN = (double) _Nevents;

  double mean   = _sum/dN;
  double meansq = _sumsq/dN;
  double v = (meansq - mean*mean)/dN;
  double sf=1.0;
  if(0 != integral()){
    sf = integral()/mean;
    v        *= sf*sf;
    // this construct allows for a scaling factor (such as _weightSum)
    // in the integral.
  }
  if(dbThis){
    cout << "FitAmpPairList::oldVariance() "
	 << " N = " << dN
	 << " sum = " << _sum
	 << " sumsq = " << _sumsq
	 << " sf " << sf
	 << " v " << v << endl;
  }
  return v;
}
double FitAmpPairList::efficiency(IDalitzEvent* evtPtr){
  if(0 == _efficiency) return 1;
  if(0 == evtPtr) return 0;
  return _efficiency->RealVal(*evtPtr);
}
void FitAmpPairList::setEfficiency(counted_ptr<IReturnRealForEvent<IDalitzEvent> > eff){
  _efficiency=eff;
}
void FitAmpPairList::unsetEfficiency(){
  _efficiency=counted_ptr<IReturnRealForEvent<IDalitzEvent> >(0);
}
DalitzHistoSet FitAmpPairList::histoSet() const{
  return un_normalised_histoSetRe()/integral();
}
DalitzHistoSet FitAmpPairList::un_normalised_histoSetRe() const{
  DalitzHistoSet sum;
  for(unsigned int i=0; i< this->size(); i++){
    sum += this->at(i).histoSetRe();
  }
  return sum;
}
DalitzHistoSet FitAmpPairList::un_normalised_histoSetIm() const{
  DalitzHistoSet sum;
  for(unsigned int i=0; i< this->size(); i++){
    sum += this->at(i).histoSetIm();
  }
  return sum;
}

void FitAmpPairList::saveEachAmpsHistograms(const std::string& prefix) const{
  DalitzHistoSet sum;
  int counter=0;
  for(unsigned int i=0; i< this->size(); i++){
    if(this->at(i).isSingleAmp()){
      counter++;
      std::string name = prefix + "_" + anythingToString(counter) + ".root";
      DalitzHistoSet hs(this->at(i).histoSet());
      std::string title = this->at(i).fitAmp1().name();
      hs.setTitle(title);
      cout << "FitAmpPairList::saveEachAmpsHistograms: "
	   << "saving " << title << " as " << name << endl;
      hs.save(name);
    }
  }
  return;
}

std::vector<DalitzHistoSet> FitAmpPairList::GetEachAmpsHistograms(){
	std::vector<DalitzHistoSet> HistoSet;
	 DalitzHistoSet sum;
	  int counter=0;
	  for(unsigned int i=0; i< this->size(); i++){
	    if(this->at(i).isSingleAmp()){
	      counter++;
	      std::string name =  anythingToString(counter) + ".root";
              double frac = this->at(i).integral()/integral();  
	      DalitzHistoSet hs(this->at(i).histoSet());
	      std::string title = this->at(i).fitAmp1().name();
	      hs.setTitle(title);
	      cout << "FitAmpPairList::saveEachAmpsHistograms: "
		   << "saving " << title << " as " << name << endl;
              hs.scale(frac/hs.integral());  
	      HistoSet.push_back(hs);
	    }
	  }
	  return HistoSet;
}

DalitzHistoSet FitAmpPairList::interferenceHistoSet() const{
    DalitzHistoSet sum;
    for(unsigned int i=0; i< this->size(); i++){
        if(!this->at(i).isSingleAmp())sum += this->at(i).histoSet();
    }
    sum /= integral();
    //  if(_Nevents > 0) sum /= (double) _Nevents;
    // above two lines normalise this to 1.
    
    return sum;
}

void FitAmpPairList::saveInterferenceHistograms(const std::string& prefix) const{
    DalitzHistoSet sum;
    int counter=0;
    for(unsigned int i=0; i< this->size(); i++){
        if(!this->at(i).isSingleAmp()){
            counter++;
            std::string name = prefix + "_" + anythingToString(counter) + ".root";
            DalitzHistoSet hs(this->at(i).histoSet());
            std::string title = this->at(i).name();
            hs.setTitle(title);
            cout << "FitAmpPairList::saveEachAmpsHistograms: "
            << "saving " << title << " as " << name << endl;
            hs.save(name);
        }
    }
    return;
}

std::vector<DalitzHistoSet> FitAmpPairList::GetInterferenceHistograms(){
    std::vector<DalitzHistoSet> HistoSet;
    DalitzHistoSet sum;
    int counter=0;
    for(unsigned int i=0; i< this->size(); i++){
	    if(!this->at(i).isSingleAmp()){
            counter++;
            std::string name =  anythingToString(counter) + ".root";
            double frac = this->at(i).integral()/integral();  
            DalitzHistoSet hs(this->at(i).histoSet());
            std::string title = this->at(i).name();
            hs.setTitle(title);
            cout << "FitAmpPairList::saveEachAmpsHistograms: "
            << "saving " << title << " as " << name << endl;
            hs.scale(frac/hs.integral());  
            HistoSet.push_back(hs);
	    }
    }
    return HistoSet;
}


//Add Plot of eveything on top of each other here
//Get EachAmpHistograms
// Can we from there get the ROOT histograms?

void FitAmpPairList::doFinalStats(Minimiser* mini){
  bool dbThis=true;
  if(dbThis) cout << "FitAmpPairList::doFinalStats() called" << endl;
  makeAndStoreFractions(mini);
}


bool FitAmpPairList::doFractions(){

  if(this->empty()) return 0;

   _singleAmpFractions = FitFractionList() ;
   _interferenceFractions = FitFractionList();

  cout << "\n============================================"
       << "============================================"
       << "\n        Amplitude Fractions";

  double norm = integral();

  if(norm <= 0){
    cout << "ERROR in FitAmpPairList::makeAndStoreFractions()"
	 << " integral = " << integral()
	 << " won't do fractions."
	 << endl;
    return false;
  }

  for(unsigned int i=0; i < this->size(); i++){    
    double frac = this->at(i).integral()/norm;
    string name;
    if(this->at(i).isSingleAmp()){
      name = this->at(i).fitAmp1().name();
    }else{
//       name = this->at(i).name();
      name = this->at(i).fitAmp1().name() + "_x_" + this->at(i).fitAmp2().name();
    }
    FitFraction f(name, frac);

    if(this->at(i).isSingleAmp()){
      this->at(i).fitAmp1().setFraction(frac);
      _singleAmpFractions.add(f);
    }else{
      _interferenceFractions.add(f);
    }
  }

  cout << "filled FitFractionLists" << endl;
  FitFractionList interferenceFractionsSorted(_interferenceFractions);
  interferenceFractionsSorted.sortByMagnitudeDecending();
  
  cout <<   "================================================="
       << "\n================================================="
       << "\n FRACTIONS:"
       << "\n ^^^^^^^^^^" << endl;
  cout << _singleAmpFractions << endl;
  cout <<   "================================================="
       << "\n================================================="
       << "\n Interference terms (sorted by size)"
       << "\n ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
  cout << interferenceFractionsSorted << endl;
  cout <<   "================================================="
       << "\n=================================================" << endl;

  cout << "\n\n\t X-check: total sum of all terms"
       << " (amplitude fractions + interference terms): "
       << _singleAmpFractions.sum().frac() 
    + _interferenceFractions.sum().frac();

  return true;
}


bool FitAmpPairList::makeAndStoreFractions(const std::string& fname
					   , const std::string&  fnameROOT
					   , Minimiser* mini
					   ){
  bool dbThis=true;

   _singleAmpFractions = FitFractionList() ;
   _interferenceFractions = FitFractionList();

  if(this->empty()) return 0;
  counted_ptr<FitAmpPairFitCovariance> fcov(0);
  if(0 != mini){
    if(mini->parSet() != this->at(0).fitAmp1().FitAmpPhase().p1().parSet()){
      cout << "ERROR in FitAmpPairList::makeAndStoreFractions"
	   << "\n Minuit parameter set and fit parameters incompatible"
	   << endl;
    }else{
      cout << "getting minuitCov" << endl;
      TMatrixTSym<double> minuitCov(mini->covMatrixFull());
      cout << "got minuitCov" << endl;
      counted_ptr<FitAmpPairFitCovariance> 
	fc(new FitAmpPairFitCovariance(this, minuitCov));
      fcov = fc;
      cout << "got fcov" << endl;
    }
  }

  bool printFitErrors = (fcov != 0 && fcov->isValid());

  cout << "\n============================================"
       << "============================================"
       << "\n        Amplitude Fractions";
  if(0 != fcov) cout << " +/- fit error ";
  cout << " +/- integration error "
       << endl;
  double norm = integral();
  if(norm <= 0){
    cout << "ERROR in FitAmpPairList::makeAndStoreFractions()"
	 << " integral = " << integral()
	 << " won't do fractions."
	 << endl;
    return false;
  }

  for(unsigned int i=0; i < this->size(); i++){    
    double frac = this->at(i).integral()/norm;
    string name;
    if(this->at(i).isSingleAmp()){
      name = this->at(i).fitAmp1().name();
    }else{
      name = this->at(i).name();
    }
    FitFraction f(name, frac);
    if(printFitErrors){
      f.sigmaFit() = fcov->getFractionError(i);
    }
    if(_cov.isValid()){ 
      f.sigmaInteg() = _cov.getFractionError(i);
    }

    if(this->at(i).isSingleAmp()){
      this->at(i).fitAmp1().setFraction(frac);
      _singleAmpFractions.add(f);
    }else{
      _interferenceFractions.add(f);
    }
  }

  if(dbThis)cout << "filled FitFractionLists" << endl;
  if(printFitErrors){
    _singleAmpFractions.setSumFitError(fcov->getFractionSumError());
    _interferenceFractions.setSumFitError(fcov->getInterferenceFracSumError());
  }
  if(_cov.isValid()){
    _singleAmpFractions.setSumIntegError(_cov.getFractionSumError());
  }

  if(dbThis) cout << "... now including errors" << endl;
  if(dbThis)cout  << " sorting" << endl;
  _interferenceFractions.sortByMagnitudeDecending();
  if(dbThis)cout << "sorted" << endl;

  cout <<   "================================================="
       << "\n================================================="
       << "\n FRACTIONS:"
       << "\n ^^^^^^^^^^" << endl;
  cout << _singleAmpFractions << endl;
  cout <<   "================================================="
       << "\n================================================="
       << "\n Interference terms (sorted by size)"
       << "\n ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
  cout << _interferenceFractions << endl;
  cout <<   "================================================="
       << "\n=================================================" << endl;

  cout << "\n\n\t X-check: total sum of all terms"
       << " (amplitude fractions + interference terms): "
       << _singleAmpFractions.sum().frac() 
    + _interferenceFractions.sum().frac();
  cout << endl;
  ofstream os(fname.c_str(),std::ofstream::trunc);
//   if(os){
//     os << _singleAmpFractions << endl;
//     os << _interferenceFractions << endl;
//     os.close();
//   }
  /// latex table
  if(!os)return true;
  os << "\\begin{tabular}{l r}" << "\n";
  os << "\\hline" << "\n";
  os << "\\hline" << "\n";
  os << "Decay channel & Fraction [$\\%$] \\\\" << "\n";
  os << "\\hline" << "\n";

  for(unsigned int i=0; i < this->size(); i++){    
    double frac = this->at(i).integral()/norm;
    TString name;
    if(this->at(i).isSingleAmp()) name = this->at(i).fitAmp1().name();
    else continue;
    FitFraction f((string)name, frac);
    double frac_err = fcov->getFractionError(i);
    
    name.ReplaceAll("->"," \\to ");
    name.ReplaceAll("Bs0","B_s");
    name.ReplaceAll("Ds","D_s");
    name.ReplaceAll("pi","\\pi");
    name.ReplaceAll("rho","\\rho");
    name.ReplaceAll("sigma10","\\sigma");
    name.ReplaceAll("GS","");
    name.ReplaceAll("SBW","");
    name.ReplaceAll("RhoOmega","");
    name.ReplaceAll("GLass","");
    name.ReplaceAll("Lass","");
    name.ReplaceAll("Bugg","");
    name.ReplaceAll("Flatte","");
    name.ReplaceAll("+","^+");
    name.ReplaceAll("-","^-");
    name.ReplaceAll("*","^*");
    name.ReplaceAll(")0",")^0");
    name.ReplaceAll(","," \\, ");
    name.ReplaceAll("NonResS0( \\to D_s^- \\, \\pi^+)","( D_s^- \\, \\pi^+)_{S}");
    name.ReplaceAll("NonResV0( \\to D_s^- \\, \\pi^+)","( D_s^- \\, \\pi^+)_{P}");
    name.ReplaceAll("NonResS0( \\to D_s^- \\, K^+)","( D_s^- \\, K^+)_{S}");
    name.ReplaceAll("NonResV0( \\to D_s^- \\, K^+)","( D_s^- \\, K^+)_{P}");

    os << std::fixed << std::setprecision(2) << "$" << name << "$ & " << frac * 100. << " $\\pm$ "  << frac_err * 100. <<  " \\\\" << "\n";
  }
  os <<  " \\hline" << "\n";
  os << " Sum & " << std::fixed << std::setprecision(2) << _singleAmpFractions.sum().frac() * 100. << " $\\pm$ "  << fcov->getFractionSumError() * 100. <<  " \\\\" << "\n";
  os << "\\hline" << "\n";
  os << "\\hline" << "\n";
  os << "\\end{tabular}" << "\n";
  os.close();

  TFile* paraFile = new TFile((fnameROOT).c_str(), "RECREATE");
  TTree* tree = new TTree("fractions","fractions");
    
  for(unsigned int i=0; i < this->size(); i++){    
	double frac = this->at(i).integral()/norm;
	TString name;
	if(this->at(i).isSingleAmp())name = TString(this->at(i).fitAmp1().name());
	else name = TString(this->at(i).name());
	FitFraction f((string)name, frac);	

	name.ReplaceAll("A(","");
	name.ReplaceAll("fit","");
	name.ReplaceAll("[","_");
	name.ReplaceAll("]","_");
	name.ReplaceAll("(","_");
	name.ReplaceAll(")","_");
	name.ReplaceAll("->","_");
	name.ReplaceAll("+","p");
	name.ReplaceAll("-","m");
	name.ReplaceAll("*","s");
	name.ReplaceAll(",","");
	name.ReplaceAll("GS","");
	name.ReplaceAll("GLass","");
	name.ReplaceAll("Lass","");
	name.ReplaceAll("RhoOmega","");
	name.ReplaceAll("Bugg","");
	name.ReplaceAll("Flatte","");

	if(0 != tree){
		double val = frac;
		if(!this->at(i).isSingleAmp()) name = "int_" + name;
		TBranch* Bra_val = tree->Branch( ((string)name).c_str(), &val, ((string)(name+"/D")).c_str());
		if(0 != Bra_val)  Bra_val->Fill();
	}
  }
  double val = _singleAmpFractions.sum().frac();
  TBranch* Bra_val = tree->Branch( "Sum", &val, "Sum/D");
  if(0 != Bra_val) Bra_val->Fill();

  if(0 != tree) tree->Fill();

  tree->Write();  
  paraFile->Close();
  delete paraFile;
  
  return true;
}



/*
bool FitAmpPairList::makeAndStoreFractions(const std::string& fname, const std::string& fnameROOT, Minimiser* mini  ){
  bool dbThis=true;
  //return 0;
  if(this->empty()) return 0;
    
  NamedParameter<string> OutputDir("OutputDir", (std::string) "", (char*) 0);  
    
  counted_ptr<FitAmpPairFitCovariance> fcov(0);
  if(0 != mini){
    if(mini->parSet() != this->at(0].fitAmp1().FitAmpPhase().p1().parSet()){
      cout << "ERROR in FitAmpPairList::makeAndStoreFractions"
	   << "\n Minuit parameter set and fit parameters incompatible"
	   << endl;
    }else{
      counted_ptr<FitAmpPairFitCovariance> 
	fc(new FitAmpPairFitCovariance(this, mini->covMatrixFull()));
      fcov = fc;
    }
  }
  
  // return 0;

  bool printFitErrors = (fcov != 0 && fcov->isValid());

  cout << "\n============================================"
       << "============================================"
       << "\n        Amplitude Fractions";
  if(0 != fcov) cout << " +/- fit error ";
  cout << " +/- integration error "
       << endl;
  double norm = integral();
  if(norm <= 0){
    cout << "ERROR in FitAmpPairList::makeAndStoreFractions()"
	 << " integral = " << integral()
	 << " won't do fractions."
	 << endl;
    return false;
  }

  TFile* paraFile = new TFile(((string)OutputDir + fnameROOT).c_str(), "RECREATE");
  //  if(0 != paraFile) paraFile->cd();
  TTree* tree = new TTree("fractions","fractions");
  int counter = 0;  
    
  for(unsigned int i=0; i < this->size(); i++){    
    double frac = this->at(i).integral()/norm;
    string name;
    if(this->at(i).isSingleAmp()){
      name = this->at(i).fitAmp1().name();
    }else{
      name = this->at(i).name();
    }
    FitFraction f(name, frac);
    if(printFitErrors){
      f.sigmaFit() = fcov->getFractionError(i);
    }
    if(_cov.isValid()){ 
      f.sigmaInteg() = _cov.getFractionError(i);
    }

    if(0 != tree && this->at(i).isSingleAmp()){
        counter++;
        double val = frac;
        double diff = (frac-this->at(i).fitAmp1().getFraction());
	double sigma = -9999;
	if(0 != fcov && fcov->isValid()) sigma = fcov->getFractionError(i);
	double pull = -9999;
        if(sigma > 0) pull = diff/sigma;
        std::string name_val = anythingToString((int)counter) + "_val";
        std::string name_pull = anythingToString((int)counter) + "_pull";
        std::string name_diff = anythingToString((int)counter) + "_diff";

        TBranch* Bra_val = tree->Branch( name_val.c_str(), &val, (name_val+"/D").c_str());
        TBranch* Bra_pull = tree->Branch( name_pull.c_str(), &pull, (name_pull+"/D").c_str());
        TBranch* Bra_diff = tree->Branch( name_diff.c_str(), &diff, (name_diff+"/D").c_str());
        if(0 != Bra_pull) Bra_pull->Fill();
        if(0 != Bra_diff) Bra_diff->Fill();
        if(0 != Bra_val)  Bra_val->Fill();
        // tree->Fill();
        //this->at(i).fitAmp1().setFraction(frac);
        _singleAmpFractions.add(f);
    }else{
      _interferenceFractions.add(f);
    }
  }

  if(0 != tree) tree->Fill();

  if(dbThis)cout << "filled FitFractionLists" << endl;
  if(printFitErrors){
    _singleAmpFractions.setSumFitError(fcov->getFractionSumError());
    _interferenceFractions.setSumFitError(fcov->getInterferenceFracSumError());
  }
  if(_cov.isValid()){
    _singleAmpFractions.setSumIntegError(_cov.getFractionSumError());
  }

  if(dbThis) cout << "... now including errors" << endl;
  if(dbThis)cout  << " sorting" << endl;
  _interferenceFractions.sortByMagnitudeDecending();
  if(dbThis)cout << "sorted" << endl;

  cout <<   "================================================="
       << "\n================================================="
       << "\n FRACTIONS:"
       << "\n ^^^^^^^^^^" << endl;
  cout << _singleAmpFractions << endl;
  cout <<   "================================================="
       << "\n================================================="
       << "\n Interference terms (sorted by size)"
       << "\n ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^" << endl;
  cout << _interferenceFractions << endl;
  cout <<   "================================================="
       << "\n=================================================" << endl;

  cout << "\n\n\t X-check: total sum of all terms"
       << " (amplitude fractions + interference terms): "
       << _singleAmpFractions.sum().frac() 
    + _interferenceFractions.sum().frac();
  cout << endl;
  ofstream os(((string) OutputDir + fname).c_str());
  if(os){
    os << _singleAmpFractions << endl;
    os << _interferenceFractions << endl;
    os.close();
  }
    
  tree->Write();  
  //paraFile->Close();
  //delete paraFile;
  
  return true;
}
*/


double FitAmpPairList::getFractionChi2() const{
  bool ErrorProportionalToSqrtTarget = true;
  double norm = integral();
  if(norm <= 0){
    cout << "ERROR in FitAmpPairList::getFractionChi2()"
	 << " integral = " << integral()
	 << " won't do fractions. Returning " << 9999.0e30
	 << endl;
    return 9999.0e30;
  }

  double canonicalError=0.001;

  int counter=0;
  double sum=0;
  for(unsigned int i=0; i < this->size(); i++){
    if(this->at(i).isSingleAmp()){
      counter++;
      double frac = this->at(i).integral()/norm;
      double target_f = this->at(i).fitAmp1().getFraction();
      //if(i < 10){
	//cout << "(" << frac << " - " << target_f << ")" << endl;
      //}
      double dfs=(frac - target_f)/canonicalError;
      if(ErrorProportionalToSqrtTarget && target_f > 0) dfs /= sqrt(0.01*target_f);
      sum += dfs*dfs;
    }
  }
  return sum;
}


std::string FitAmpPairList::dirName() const{
  return "FitAmpPairList";
}
bool FitAmpPairList::save(const std::string& asSubdirOf) const{
  bool sc=true;
  makeDirectory(asSubdirOf);
  ofstream os((asSubdirOf +"/"+ dirName() + "values.txt").c_str());
  os << "N " << _Nevents << endl;
  os << "sum " << setprecision(20) << _sum << endl;
  os << "sumsq " << setprecision(20) << _sumsq << endl;
  os << "psSum " << setprecision(20) << _psSum << endl;
  os << "psSumSq " << setprecision(20) << _psSumSq << endl;
  os.close();

  for(unsigned int i=0; i < this->size(); i++){
    sc &= this->at(i).save(asSubdirOf +"/"+ dirName());
  }
  if(_cov.isValid()){
    _cov.save(asSubdirOf + "/" + dirName());
  }else{
    cout << "WARNING IN FitAmpPairList::save(" << asSubdirOf
	 << "): invalid _cov" << endl;
  }

  return sc;
}
bool FitAmpPairList::retrieve(const std::string& asSubdirOf){
  bool sc=true;
  bool verbose=true;
  if(verbose){
    cout << "FitAmpPairList::retrieve: reading from "
	 << asSubdirOf << endl; 
  }

  struct stat buf;
  if(stat(asSubdirOf.c_str(), &buf) != 0){
    if(verbose) cout << "FitAmpPairList::retrieve: " << asSubdirOf
		     << " does not exist" << endl;
    return false;
  }

  ifstream is((asSubdirOf +"/"+ dirName() + "values.txt").c_str());
  std::string dummy;
  is >> dummy >> _Nevents;
  is >> dummy >> _sum;
  is >> dummy >> _sumsq;
  is >> dummy >> _psSum;
  is >> dummy >> _psSumSq;
  is.close();
  for(unsigned int i=0; i < this->size(); i++){
    sc &= this->at(i).retrieve(asSubdirOf +"/"+ dirName());
  }
  _cov.retrieve(asSubdirOf + "/" + dirName());
  return sc;
}

bool FitAmpPairList::makeDirectory(const std::string& asSubdirOf)const{
  /*
    A mode is created from or'd permission bit masks defined
     in <sys/stat.h>:
           #define S_IRWXU 0000700     RWX mask for owner 
           #define S_IRUSR 0000400     R for owner 
           #define S_IWUSR 0000200     W for owner 
           #define S_IXUSR 0000100     X for owner 

           #define S_IRWXG 0000070     RWX mask for group 
           #define S_IRGRP 0000040     R for group 
           #define S_IWGRP 0000020     W for group 
           #define S_IXGRP 0000010     X for group 

           #define S_IRWXO 0000007     RWX mask for other 
           #define S_IROTH 0000004     R for other 
           #define S_IWOTH 0000002     W for other 
           #define S_IXOTH 0000001     X for other 

           #define S_ISUID 0004000     set user id on execution 
           #define S_ISGID 0002000     set group id on execution 
           #define S_ISVTX 0001000     save swapped text even after use
   */

  mode_t mode = S_IRWXU | S_IRWXG | S_IRWXO 
              | S_ISUID | S_ISGID; 
  // see above for meaning. I want everybody to be allowed to read/write/exec.
  // Not sure about the last two bits.

  int zeroForSuccess = 0;
  zeroForSuccess |= mkdir( (asSubdirOf ).c_str(), mode );
  zeroForSuccess |= mkdir( (asSubdirOf + "/" + dirName() ).c_str(), mode );
  return (0 == zeroForSuccess);
}

void FitAmpPairList::print(std::ostream& os) const{
  os << "FitAmpPairList with " << this->size() << " FitAmpPairs:";
  for(unsigned int i=0; i < this->size(); i++){
    os << "\n\t" << i << ") " << this->at(i);
  }
}

bool FitAmpPairList::needToReIntegrate() const{
  for(unsigned int i=0; i < this->size(); i++){    
    if (this->at(i).needToReIntegrate() ) return true;
  }
  return false;
}
void FitAmpPairList::startIntegration(){
  for(unsigned int i=0; i < this->size(); i++){    
    this->at(i).startIntegration();
  }
}
void FitAmpPairList::startReIntegration(){
  for(unsigned int i=0; i < this->size(); i++){    
    if( this->at(i).needToReIntegrate() ) this->at(i).startReIntegration();
  }
}
void FitAmpPairList::startForcedReIntegration(){
  for(unsigned int i=0; i < this->size(); i++){    
    this->at(i).startReIntegration();
  }
}
void FitAmpPairList::endIntegration(){
  for(unsigned int i=0; i < this->size(); i++){    
    this->at(i).endIntegration();
  }
}

void FitAmpPairList::setSlow(){
  _slow=true;
  for(unsigned int i=0; i < this->size(); i++){    
    this->at(i).setSlow();
  }
}
void FitAmpPairList::setFast(){
  _slow=false;
  for(unsigned int i=0; i < this->size(); i++){    
    this->at(i).setFast();
  }
}

FitAmpPairList& FitAmpPairList::operator+=(const FitAmpPairList& other){
  this->add(other);
  return *this;
}
FitAmpPairList FitAmpPairList::operator+(const FitAmpPairList& other) const{
  FitAmpPairList returnVal(*this);
  returnVal.add(other);
  return returnVal;
}

std::ostream& operator<<(std::ostream& os, const FitAmpPairList& fap){
  fap.print(os);
  return os;
}

/*
double FitAmpPairList::variance() const{
  double sum = 0;
  for(unsigned int i=0; i< this->size(); i++){
    sum += this->at(i).variance();
  }
  return sum;
}
*/

/*
double FitAmpPairList::variance() const{
  double sumReSquared = 0;
  double sumImSquared = 0;
  double sumImRe      = 0;

  for(unsigned int i=0; i< this->size(); i++){
    sumReSquared += this->at(i).ReSquared();
    sumImSquared += this->at(i).ImSquared();
    sumImRe      += this->at(i).ImRe();
  }
  double sumsq = sumReSquared + sumImSquared - 2.0*sumImRe;

  double meansq = sumsq/((double)_Nevents);
  
  double mean = integral();

  double var = (meansq - mean*mean)/((double)_Nevents);

  return var;
}
*/

//
