#ifndef IDALITZ_EVENT_HH
#define IDALITZ_EVENT_HH
// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:00 GMT

#include "TLorentzVector.h"
#include "TVector3.h"

#include "Mint/IWeightedEvent.h"
#include "Mint/DalitzEventPattern.h"

#include <vector>
#include <iostream>
#include <complex>

class IDalitzEvent : public virtual MINT::IWeightedEvent{
 protected:
  IDalitzEvent(){}
 public:
  virtual void setAValue(double aValue)=0;
  virtual double getAValue()const=0;

  virtual void setWeight(double w)=0;
  virtual double getWeight()const=0;

  virtual void   setGeneratorPdfRelativeToPhaseSpace(double gpdf)=0;
  virtual double getGeneratorPdfRelativeToPhaseSpace()const=0;
    
  virtual const std::vector<double>& getVectorOfValues() const=0;
  virtual std::vector<double>& getVectorOfValues()=0;
  virtual const std::vector<double>& getVectorOfWeights() const=0;
  virtual std::vector<double>& getVectorOfWeights()=0;
  virtual void setValueInVector(unsigned int i, double value)=0;
  virtual void setWeightInVector(unsigned int i, double weight)=0;
  virtual double getValueFromVector(unsigned int i) const=0;
  virtual double getWeightFromVector(unsigned int i) const=0;

  virtual const DalitzEventPattern& eventPattern() const=0;
  virtual const TLorentzVector& p(unsigned int i) const= 0;

  virtual void setMothers3Momentum(const TVector3& mp3)=0;

  virtual double s(unsigned int i, unsigned int j) const= 0;
  virtual double sij(const MINT::PolymorphVector<int>& indices) const= 0;
  virtual double t(unsigned int i, unsigned int j) const= 0;
  
  virtual double sijMin(const MINT::PolymorphVector<int>& indices) const=0;
  virtual double sijMax(const MINT::PolymorphVector<int>& indices) const=0;

  virtual double phaseSpace() const= 0;

  virtual void print(std::ostream& os = std::cout) const=0;

  virtual bool retrieveValue(int i, std::complex<double>& value, long int configNumber)=0;
  virtual void setValue(int i, const std::complex<double>& value, long int configNumber)=0;

  virtual bool retrieveValue(int i, double value, long int configNumber)=0;
  virtual void setValue(int i, double value, long int configNumber)=0;

  virtual int numPermutations() const=0;
  virtual void setPermutationIndex(int i)=0;
  virtual int permutationIndex() const=0;  

  virtual IDalitzEvent* clone() const=0;
  virtual ~IDalitzEvent(){}
};

bool EqualEvent(const IDalitzEvent* a, const IDalitzEvent* b);

std::ostream& operator<<(std::ostream& os, const IDalitzEvent& de);

#endif
//
