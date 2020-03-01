// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:03 GMT
#include "Mint/FitAmpListBase.h"

#include "Mint/FitAmplitude.h"
#include "Mint/MinuitParameterSet.h"
#include "Mint/NamedDecayTreeList.h"
#include "Mint/FitAmplitude.h"
#include "Mint/FitAmpSum.h"
#include "Mint/FitAmpIncoherentSum.h"
#include "Mint/FitAmpPairList.h"
#include "TRandom3.h"
#include <iostream>

using namespace std;
using namespace MINT;

FitAmpListBase::FitAmpListBase()
  : _efficiency(0)
{
}

FitAmpListBase::FitAmpListBase(const FitAmpListBase& other)
  : _efficiency(other._efficiency)
{
  bool dbThis=false;
  this->deleteAll();
  if(dbThis)cout << "copy-ctor FitAmpListBase, done deleteAll()" << endl;
  /* 
     There'll be 'physical' copies of all Amplitudes, but the
     FitParameters remain the same (pointers to the same
     FitParameter Object).  This is useful for the CP-con coding
     as it is now, but perhaps a bit counter-intuitive.  Needs to
     be reviewed at some point. This behaviour is defined in the
     copy construcopy-ctor of the FitAmplitude class.
  */

  addCopyWithSameFitParameters(other);
}

FitAmpListBase::FitAmpListBase(const FitAmpListBase& other, string name)
: _efficiency(other._efficiency)
{
    bool dbThis=false;
    this->deleteAll();
    if(dbThis)cout << "copy-ctor FitAmpListBaseSubset, done deleteAll()" << endl;
    /* 
     There'll be 'physical' copies of all Amplitudes, but the
     FitParameters remain the same (pointers to the same
     FitParameter Object).  This is useful for the CP-con coding
     as it is now, but perhaps a bit counter-intuitive.  Needs to
     be reviewed at some point. This behaviour is defined in the
     copy construcopy-ctor of the FitAmplitude class.
     */
    
    addCopyOfSubsetWithSameFitParameters(other,(string) name);
}

FitAmpListBase& FitAmpListBase::operator=(const FitAmpListBase& other){
  if(&other == this) return *this;
  _efficiency    = other._efficiency;
  deleteAll();  
  addCopyWithSameFitParameters(other);
  return *this;
}
int FitAmpListBase::add(const FitAmpListBase& other, double factor){
  return addCopyWithSameFitParameters(other, factor);
}
int FitAmpListBase::addAsList(const FitAmpListBase& other, double factor){
  if(0 == other.size()) return size();
  // keeps them together
  //cout << "Hello from FitAmpListBase::addAsList" << endl;
  counted_ptr<FitAmpListBase> newList(other.GetCloneSameFitParameters());
  if(1.0 != factor) newList->multiply(factor);
  _fitAmpLists.push_back(newList);
  //cout << "all done in FitAmpListBase::addAsList" << endl;
  return size();
}

int FitAmpListBase::addCopyWithSameFitParameters(const FitAmpListBase& other
					     , double factor){
  for(unsigned int i=0; i < other._fitAmps.size(); i++){
    FitAmplitude* fa = other._fitAmps[i];
    FitAmplitude* newFa = new FitAmplitude(*fa);
    if(1.0 != factor) newFa->multiply(factor);
    //_fitAmps.push_back(newFa);
    this->addAmplitude(newFa);
  }
  for(unsigned int i=0; i < other._fitAmpLists.size(); i++){
    this->addAsList(*(other._fitAmpLists[i]->GetCloneSameFitParameters()));
  }
  return this->size();
}

int FitAmpListBase::addCopyOfSubsetWithSameFitParameters(const FitAmpListBase& other, string name, double factor){
    for(unsigned int i=0; i < other._fitAmps.size(); i++){
        if(! A_is_in_B(name,other._fitAmps[i]->name()))continue;
        FitAmplitude* fa = other._fitAmps[i];
        FitAmplitude* newFa = new FitAmplitude(*fa);
        if(1.0 != factor) newFa->multiply(factor);
        this->addAmplitude(newFa);
    }
    for(unsigned int i=0; i < other._fitAmpLists.size(); i++){
        this->addAsList(*(other._fitAmpLists[i]->GetCloneOfSubsetSameFitParameters(name)));
    }
    return this->size();
}

unsigned int FitAmpListBase::size() const{
  unsigned int sz(_fitAmps.size());
  for(unsigned int i=0; i < _fitAmpLists.size(); i++){
    sz += _fitAmpLists[i]->size();
  }
  return sz;
}

FitAmplitude* FitAmpListBase::getAmpPtr(unsigned int in){

  if(in < _fitAmps.size()) return _fitAmps[in];
  unsigned int sz = _fitAmps.size();
  for(unsigned int i = 0; i < _fitAmpLists.size(); i++){
    int index = in - sz;
    if(index >= 0 && static_cast<unsigned int>(index) < _fitAmpLists[i]->size()) return _fitAmpLists[i]->getAmpPtr(index);
    sz += _fitAmpLists[i]->size();
  }
  
  cout << " FitAmpListBase::getAmp index out of range"
       << endl;
  return 0;
}
const FitAmplitude* FitAmpListBase::getAmpPtr(unsigned int in) const{

  if(in < _fitAmps.size()) return _fitAmps[in];
  unsigned int sz = _fitAmps.size();
  for(unsigned int i = 0; i < _fitAmpLists.size(); i++){
    int index = in - sz;
    if(index >= 0 && static_cast<unsigned int>(index) < _fitAmpLists[i]->size()) return _fitAmpLists[i]->getAmpPtr(index);
    sz += _fitAmpLists[i]->size();
  }
  
  cout << " FitAmpListBase::getAmp index out of range"
       << endl;
  return 0;
}

bool FitAmpListBase::CPConjugateSameFitParameters(){
  bool dbThis=false;
  if(dbThis) cout << "FitAmpListBase::CPConjugateSameFitParameters()" << endl;

  bool success=true;
  for(unsigned int i=0; i< this->size(); i++){
    success &= getAmpPtr(i)->CPConjugateSameFitParameters();
  }
  return success;
}

bool FitAmpListBase::CConjugateFinalStateSameFitParameters(){
    bool dbThis=false;
    if(dbThis) cout << "FitAmpListBase::CConjugateFinalStateSameFitParameters()" << endl;
    
    bool success=true;
    for(unsigned int i=0; i< this->size(); i++){
        success &= getAmpPtr(i)->CConjugateFinalStateSameFitParameters();
    }
    return success;
}

bool FitAmpListBase::CConjugateInitialStateSameFitParameters(){
    bool dbThis=false;
    if(dbThis) cout << "FitAmpListBase::CConjugateFinalStateSameFitParameters()" << endl;
    
    bool success=true;
    for(unsigned int i=0; i< this->size(); i++){
        success &= getAmpPtr(i)->CConjugateInitialStateSameFitParameters();
    }
    return success;
}

bool FitAmpListBase::setLSameFitParameters(int L){
    bool dbThis=false;
    if(dbThis) cout << "FitAmpListBase::setLSameFitParameters()" << endl;
    
    bool success=true;
    for(unsigned int i=0; i< this->size(); i++){
        success &= getAmpPtr(i)->setLSameFitParameters(L);
    }
    return success;
}

counted_ptr<FitAmpListBase> FitAmpListBase::GetCloneSameFitParameters() const{ 
  bool dbThis=false;
  if(dbThis) cout << "FitAmpListBase::GetCloneSameFitParameters()" << endl;
  /* 
     There'll be 'physical' copies of all Amplitudes, but the
     FitParameters remain the same (pointers to the same
     FitParameter Object).  This is useful for the CP-con coding
     as it is now, but perhaps a bit counter-intuitive.  needs to
     be reviewed at some point. This behaviour is defined in the
     copy constructor of the FitAmplitude class.
  */

  counted_ptr<FitAmpListBase> newList(new FitAmpListBase(*this));
  return newList;
}

counted_ptr<FitAmpListBase> FitAmpListBase::GetCloneOfSubsetSameFitParameters(string name) const{ 
    bool dbThis=false;
    if(dbThis) cout << "FitAmpListBase::GetCloneSameFitParameters()" << endl;
    /* 
     There'll be 'physical' copies of all Amplitudes, but the
     FitParameters remain the same (pointers to the same
     FitParameter Object).  This is useful for the CP-con coding
     as it is now, but perhaps a bit counter-intuitive.  needs to
     be reviewed at some point. This behaviour is defined in the
     copy constructor of the FitAmplitude class.
     */
    
    counted_ptr<FitAmpListBase> newList(new FitAmpListBase(*this,name));
    return newList;
}

counted_ptr<FitAmpListBase> FitAmpListBase::GetCPConjugateSameFitParameters() const{
  bool dbThis=false;
  if(dbThis) cout << "FitAmpListBase::GetCPConjugateSameFitParameters()" << endl;

  counted_ptr<FitAmpListBase> newList = this->GetCloneSameFitParameters();
  newList->CPConjugateSameFitParameters();
  return newList;
}

counted_ptr<FitAmpListBase> FitAmpListBase::GetCConjugateFinalStateSameFitParameters() const{
    bool dbThis=false;
    if(dbThis) cout << "FitAmpListBase::GetCConjugateFinalStateSameFitParameters()" << endl;
    
    counted_ptr<FitAmpListBase> newList = this->GetCloneSameFitParameters();
    newList->CConjugateFinalStateSameFitParameters();
    return newList;
}

counted_ptr<FitAmpListBase> FitAmpListBase::GetCConjugateInitialStateSameFitParameters() const{
    bool dbThis=false;
    if(dbThis) cout << "FitAmpListBase::GetCConjugateFinalStateSameFitParameters()" << endl;
    
    counted_ptr<FitAmpListBase> newList = this->GetCloneSameFitParameters();
    newList->CConjugateInitialStateSameFitParameters();
    return newList;
}

counted_ptr<FitAmpListBase> FitAmpListBase::GetDifferentLSameFitParameters(int L) const{
    bool dbThis=false;
    if(dbThis) cout << "FitAmpListBase::GetDifferentLSameFitParameters()" << endl;
    
    counted_ptr<FitAmpListBase> newList = this->GetCloneSameFitParameters();
    newList->setLSameFitParameters(L);
    return newList;
}

bool FitAmpListBase::addAmplitude(FitAmplitude* fa){
  bool dbThis=false;
  if(0 == fa) return false;
  if(dbThis) cout << "check init values: " << *fa << endl;
  _fitAmps.push_back(fa);
  return true;
}
void FitAmpListBase::printLargestAmp(IDalitzEvent& evt, std::ostream& os){
  bool dbthis=false;
  if(size()==0){
    os << "FitAmpListBase::printLargestAmp: list is empty" << endl;
    return;
  }

  double largestValue = -9999;
  std::string largestName = "none";

  for(unsigned int i=0; i< this->size(); i++){
    if(dbthis){
      cout << "FitAmpListBase::printLargestAmp()"
	   << "\n     > for " << getAmpPtr(i)->theBareDecay().oneLiner()
	   << "\n     > I get " << getAmpPtr(i)->getVal(evt)
	   << endl;
    }
    double val = norm(getAmpPtr(i)->getVal(evt));
    if(val > largestValue){
      largestValue = val;
      largestName = getAmpPtr(i)->name();
    }
  }
  os << "largest amp for event " << evt
     << "\n is " << largestName 
     << " with value " << largestValue 
     << endl;
}


void FitAmpListBase::printAllAmps(IDalitzEvent& evt, std::ostream& os){
  bool dbThis=false;
  if(size()==0){
    os << "FitAmpListBase::printAllAmps: list is empty" << endl;
    return;
  }

  std::string largestName = "none";
  if(dbThis) cout << "Debug mode for  FitAmpListBase::printAllAmps" << endl;

  os << "FitAmpListBase::printAllAmps()\n====================";

  for(unsigned int i=0; i< this->size(); i++){
    os << "\n\t" << getAmpPtr(i)->theBareDecay().oneLiner()
       << " \t" << getAmpPtr(i)->getVal(evt)
       << endl;
  }
}
void FitAmpListBase::printAllAmps(std::ostream& os)const{
  bool dbThis=false;
  if(size() == 0){
    os << "FitAmpListBase::printAllAmps: list is empty" << endl;
    return;
  }

  std::string largestName = "none";
  if(dbThis) cout << "Debug mode for  FitAmpListBase::printAllAmps" << endl;

  os << "FitAmpListBase::printAllAmps()\n====================";

  for(unsigned int i=0; i< this->size(); i++){
    os << "\n\t" << getAmpPtr(i)->theBareDecay().oneLiner()
       << endl;
  }
}

void FitAmpListBase::printNonZeroWithValue(IDalitzEvent& evt, std::ostream& os){
  bool dbThis=false;
  if(size() == 0){
    os << "FitAmpListBase::printNonZeroWithValue: list is empty" << endl;
    return;
  }

  std::string largestName = "none";
  if(dbThis) cout << "Debug mode for  FitAmpListBase::printAllAmps" << endl;
  
  os << "FitAmpListBase::printNonZeroWithValue\n====================\n";

  for(unsigned int i=0; i < this->size(); i++){
    if(getAmpPtr(i)->isZero()) continue;
    os << "\t" << getAmpPtr(i)->theBareDecay().oneLiner()
       << " \t" << getAmpPtr(i)->getVal(evt)
       << endl;
  }
}
void FitAmpListBase::print(std::ostream& os) const{
   os << "FitAmpListBase::print\n====================";

  for(unsigned int i=0; i< this->size(); i++){
    os << "\n\t" << getAmpPtr(i)->theBareDecay().oneLiner()
       << endl;
  }
}
void FitAmpListBase::printNonZero(std::ostream& os) const{
   os << "FitAmpListBase::printNonZero\n====================";

  for(unsigned int i=0; i< this->size(); i++){
    if(getAmpPtr(i)->isZero()) continue;
    os << "\n\t" << getAmpPtr(i)->theBareDecay().oneLiner()
       << endl;
  }
}

void FitAmpListBase::setAllAmpsTo(std::complex<double> z){
  for(unsigned int i=0; i< this->size(); i++){
    getAmpPtr(i)->FitAmpPhase().set(z);
  }
}

DalitzBoxSet FitAmpListBase::makeBoxes(const DalitzEventPattern& pat
				   , IReturnRealForEvent<IDalitzEvent>* pdf
				   , double nSigma){
  DalitzBoxSet boxes;
  DalitzBox phaseSpaceBox(pat);
  boxes.add(phaseSpaceBox);

  for(unsigned int i=0; i< this->size(); i++){
    boxes.add( getAmpPtr(i)->MakeBoxes(pat, nSigma) );
  }
  boxes.setPDF(pdf);
  return boxes;
}

DalitzBWBoxSet FitAmpListBase::makeBWBoxes(const DalitzEventPattern& pat
				       , IReturnRealForEvent<IDalitzEvent>* pdf
				       , TRandom* rnd
				      ){
  DalitzBWBoxSet boxes(pdf, rnd);
  //  DalitzBox phaseSpaceBox(pat);
  //  boxes.add(phaseSpaceBox);

  for(unsigned int i=0; i< this->size(); i++){
    DalitzBWBoxSet oneAmpsBox(getAmpPtr(i)->MakeBWBoxes(pat, rnd));
    boxes.add(oneAmpsBox);
  }
  return boxes;
}


void FitAmpListBase::multiply(double r){
  for(unsigned int i=0; i< this->size(); i++){
    if(0 != getAmpPtr(i)){
      getAmpPtr(i)->multiply(r);
    }
  }
}
void FitAmpListBase::multiply(const std::complex<double>& z){
  for(unsigned int i=0; i< this->size(); i++){
    if(0 != getAmpPtr(i)){
      getAmpPtr(i)->multiply(z);
    }
  }
}
void FitAmpListBase::multiply(const MINT::counted_ptr<MINT::IReturnComplex> irc){
  for(unsigned int i=0; i< this->size(); i++){
    if(0 != getAmpPtr(i)){
      getAmpPtr(i)->multiply(irc);
    }
  }
}
void FitAmpListBase::multiply(MINT::counted_ptr<MINT::IReturnComplexForEvent<IDalitzEvent> > irc){
  for(unsigned int i=0; i< this->size(); i++){
    if(0 != getAmpPtr(i)){
      getAmpPtr(i)->multiply(irc);
    }
  }
}

void FitAmpListBase::deleteAll(){
  bool dbThis=false;
  if(dbThis) cout << "Hello from FitAmpListBase::deleteAll()" << endl;
  for(unsigned int i=0; i<_fitAmps.size(); i++){
    if(0 != (_fitAmps[i])){
      delete (_fitAmps[i]);
      _fitAmps[i]=0;
    }
  }
  for(unsigned int i=0; i < _fitAmpLists.size(); i++){
    if(0 != _fitAmpLists[i]){
      _fitAmpLists[i]->deleteAll();
    }
  }
 
  _fitAmps.clear();
  _fitAmpLists.clear();
  if(dbThis) cout << "all fine in FitAmpListBase::deleteAll()" << endl;
}

void FitAmpListBase::setEfficiency(const counted_ptr<IReturnRealForEvent<IDalitzEvent> >& eff){
  _efficiency=eff;
}
double FitAmpListBase::efficiency(IDalitzEvent& evt){
  if(0 == _efficiency) return 1.0;
  double eff = _efficiency->RealVal(evt);
  if(eff < 0) return 0;
  return eff;
}

void FitAmpListBase::normalizeAmps(DalitzEventList& evtList){
    
    for(unsigned int i=0; i< this->size(); i++){
      if(0 == getAmpPtr(i))continue;
      double integral=0.;
      double weight_sum=0.;
      for (unsigned int j=0; j<evtList.size(); j++) {
	double weight = evtList[j].getWeight()/evtList[j].getGeneratorPdfRelativeToPhaseSpace();
	weight_sum += weight;
	integral += weight * std::norm(getAmpPtr(i)->getValWithoutFitParameters(evtList[j]));
      }
      if(weight_sum==0)weight_sum = evtList.size(); 
      if(integral>0)(getAmpPtr(i))->multiply(sqrt(weight_sum/integral));
    }
}

vector<double> FitAmpListBase::normFactors(DalitzEventList& evtList){
    
    	vector<double> tmp;
        	for(unsigned int i=0; i< this->size(); i++){
              		if(0 == getAmpPtr(i))continue;
              		double integral=0.;
              		double weight_sum=0.;
              		for (unsigned int j=0; j<evtList.size(); j++) {
            			double weight = evtList[j].getWeight()/evtList[j].getGeneratorPdfRelativeToPhaseSpace();
            			weight_sum += weight;
            			integral += weight * std::norm(getAmpPtr(i)->amp().getNewVal(evtList[j]));
                  		}
              		if(weight_sum==0)weight_sum = evtList.size(); 
                 	if(integral>0)tmp.push_back(sqrt(weight_sum/integral));
            	}
        	return tmp;
}

void FitAmpListBase::randomizeStartVals(int seed){
        TRandom3* r = new TRandom3(seed);    
    
        for(unsigned int i=0; i< this->size(); i++){
              //if(0 == getAmpPtr(i))continue;
              double mag = r->Uniform(0.,1.);
              double re,im;
              r->Circle(re,im,mag);
              complex<double> c(re,im);
              //getAmpPtr(i)->multiply(c);
             if(!(getAmpPtr(i)->FitAmpPhase().isConstant()))getAmpPtr(i)->FitAmpPhase().set(c);
            }
}


void FitAmpListBase::randomizePhaseStartVals(int seed){
        TRandom3* r = new TRandom3(seed);    
    
        for(unsigned int i=0; i< this->size(); i++){
              //if(0 == getAmpPtr(i))continue;
              double phase = r->Uniform(-pi,pi);
              if(!(getAmpPtr(i)->FitAmpPhase().isConstant()))getAmpPtr(i)->FitAmpPhase().set( polar(getAmpPtr(i)->FitAmpPhase().getAmp(),phase)  );
            }
}

void FitAmpListBase::setTag(int tag){
    for(unsigned int i=0; i< this->size(); i++){
        this->getAmpPtr(i)->setTag(tag);
    }
}

FitAmpListBase::~FitAmpListBase(){
  deleteAll();
}

FitAmpListBase& FitAmpListBase::operator+=(const FitAmpListBase& other){
  add(other);
  return *this;
}
FitAmpListBase FitAmpListBase::operator+(const FitAmpListBase& rhs) const{
  FitAmpListBase fas(*this);
  fas.add(rhs);
  return fas;
}


FitAmpListBase& FitAmpListBase::operator*=(double r){
  multiply(r);
  return *this;
}
FitAmpListBase& FitAmpListBase::operator*=(const complex<double>& z){
  multiply(z);
  return *this;
}
FitAmpListBase& FitAmpListBase::operator*=(const counted_ptr<IReturnComplex>& irc){
  multiply(irc);
  return *this;
}

FitAmpListBase FitAmpListBase::operator*(double r) const{
  FitAmpListBase fas(*this);
  fas.multiply(r);
  return fas;
}
FitAmpListBase FitAmpListBase::operator*(const complex<double>& z) const{
  FitAmpListBase fas(*this);
  fas.multiply(z);
  return fas;
}
FitAmpListBase FitAmpListBase::operator*(const counted_ptr<IReturnComplex>& irc) const{
  FitAmpListBase fas(*this);
  fas.multiply(irc);
  return fas;
}


FitAmpListBase operator*(double r, const FitAmpListBase& rhs){
  FitAmpListBase fas(rhs);
  fas.multiply(r);
  return fas;
}
FitAmpListBase operator*(const complex<double>& z, const FitAmpListBase& rhs){
  FitAmpListBase fas(rhs);
  fas.multiply(z);
  return fas;
}
FitAmpListBase operator*(const counted_ptr<IReturnComplex>& irc
		     , const FitAmpListBase& rhs){
  FitAmpListBase fas(rhs);
  fas.multiply(irc);
  return fas;
}

std::ostream& operator<<(std::ostream& os, const FitAmpListBase& fal){
  fal.print(os);
  return os;
}

//
