// author: Jonas Rademacker (Jonas.Rademacker@bristol.ac.uk)
// status:  Mon 9 Feb 2009 19:18:13 GMT

#include "Mint/SpinFactors.h"

#include "Mint/DecayTree.h"
#include "Mint/SpinSumV.h"
#include "Mint/SpinSumT.h"
#include "Mint/Utils.h"
#include "Mint/CLHEPSystemOfUnits.h"

#include "Mint/ZTspin1.h"
#include "Mint/ZTspin2.h"

#include <iostream>

using namespace std;
using namespace MINT;

DecayTree* SpinFactor3::_exampleDecay=0;

const DecayTree& SpinFactor3::getExampleDecay(){
  if(0==_exampleDecay){
    _exampleDecay = new DecayTree(421);
    _exampleDecay->addDgtr(310, 113)->addDgtr(211, -211);// Ks rho - why not
    _exampleDecay->addDgtr(-211, 200321)->addDgtr(321, -10323);
  }
  return *_exampleDecay;
}

const DecayTree& SpinFactor3::exampleDecay(){
  return getExampleDecay();
}

bool SpinFactor3::setSpin(){
  if(0 == R){
    cout << " SpinFactor3::setSpin can't set spin, R==0"
	 << " probably need to parse tree, first."
	 << endl;
    return false;
  }
  if(R->getVal().J() == "0"){
    _spin = 0;
  }else if(R->getVal().J() == "1"){
    _spin = 1;
  }else if(R->getVal().J() == "2"){
    _spin = 2;
  }else if(R->getVal().J() == "3"){
      _spin = 3;
  }else{    
    cout << "SpinFactor3::setSpin() Don't know how to"
	 << " handle resonance with spin = " << R->getVal().J()
	 << "\n    > looking at decay\n" << theDecay()
	 << "\n    > Sorry, will crash now!" << endl;
    throw "sorry";
  }
  return true;
}

bool SpinFactor3::parseTree(const DalitzEventPattern& pat){
  if(fsPS.size() < 3) fsPS.reserve(3);

  if(theDecay(pat).nDgtr() == (int)(theDecay(pat).finalState().size())){
    _nonResonant = true;
    return true;
  }

  for(int i=0; i< theDecay(pat).nDgtr(); i++){
    const_counted_ptr<AssociatedDecayTree> dgtr= theDecay(pat).getDgtrTreePtr(i);
    if   (! dgtr->isFinalState())      R = dgtr;
    else                         fsPS[0] = dgtr;
  }
  if(0==R || 0==fsPS[0]){
    cout << "SpinFactor3::parseTree"
	 << " Didn't find R or P1 " << R << ", " << fsPS[0] << endl;
    return false;
  }
  if(R->nDgtr() != 2){
    cout << "ERROR in SpinFactor3::parseTree"
	 << " Resonance should have 2 daughters, but it says it has "
	 << R->nDgtr() << "."
	 << endl;
    return false;
  }
  fsPS[1] = R->getDgtrTreePtr(0);
  fsPS[2] = R->getDgtrTreePtr(1);

  //normalOrder(fsPS[1], fsPS[2]);

//  printYourself();

  setSpin();
  return true;

}

void SpinFactor3::printYourself(std::ostream& os)const{
  os << "INFO from SpinFactor3::printYourself():\n"
     << "    > parsed the following tree "
     << theDecay().oneLiner() << " like this:\n"
     << "    > R =  " << R->getVal().name()
     << ", fsPS[0]=C= " << fsPS[0]->getVal().name()
     << ", fsPS[1]=A= " << fsPS[1]->getVal().name()
     << ", fsPS[2]=B= " << fsPS[2]->getVal().name()
     << endl;
}

double SpinFactor3::getVal(IDalitzEvent& evt){
  if(_nonResonant){
    return nonResVal();
  }
  
  if(! ( fsPS[0] && fsPS[1] && fsPS[2])) parseTree(evt.eventPattern());
  
  if(0 == R){
    cout << "ERROR in SpinFactor3::getVal(): 0 == R"
	 << endl;
    throw "such things shouldn't happen";
  }

  if(_spin == 0){
    return spinZeroVal();
  }else if(_spin == 1){
    return spinOneFromZemach(evt);
  }else if(_spin == 2){
    return spinTwoFromZemach(evt);
  }else if(_spin == 3){
      return spinThreeFromZemach(evt);
  }else{    
    cout << "SpinFactor3::getVal() Don't know how to"
	 << " handle resonance with spin = " << R->getVal().J()
	 << "\n    > looking at decay\n" << theDecay()
	 << "\n    > Sorry, will crash now!" << endl;
    throw "sorry";
  }
  return -9999;
}

double SpinFactor3::GSSpinFactor(IDalitzEvent& evt){

  double m2AC = (p(0, evt) + p(1, evt)).M2();
  double m2BC = (p(0, evt) + p(2, evt)).M2();

  //     cout << "GS spin factor " << m2AC - m2BC << endl;
  return (m2AC - m2BC)/(GeV*GeV);
}

double SpinFactor3::spinOneVal(IDalitzEvent& evt){
  // parsed as:
  // D -> V P0; V->P1, P2;
  // fsPS[0] && fsPS[1] && fsPS[2]

  bool dbThis=false;

  TLorentzVector pV = p(1, evt) + p(2, evt);
  TLorentzVector pD = pV + p(0, evt);
  double mr = mRes(R, evt);
  SpinSumV spin_sum(pV, mr);

  TLorentzVector lhs = pD + p(0, evt);
  TLorentzVector rhs = p(1, evt) - p(2, evt);

  if(dbThis){
    cout << " spinOneVal for " << theDecay().oneLiner()
	 << "\n   > compare: Sandwich: "  
	 << -spin_sum.Sandwich(lhs, rhs)
	 << "\n   > from masses             "  
	 << spinOneFromMasses(evt) 
	 << "\n   > ratio sw/m         "  
	 << -spin_sum.Sandwich(lhs, rhs)/spinOneFromMasses(evt) 
	 << "\n   > from Zemach " << spinOneFromZemach(evt)
	 << "\n   > ratios sw/Zemach "
	 << spin_sum.Sandwich(lhs, rhs)/spinOneFromZemach(evt)
	 << endl;
  }

  // note '-' sign. This makes my decay-tree
  // sorting algorithm and this spin factor
  // comparible with BaBar/BELLE spin
  // factor conventions.
  return -spin_sum.Sandwich(lhs, rhs)/(GeV*GeV);
}

double SpinFactor3::spinOneFromZemach(IDalitzEvent& evt){
  TLorentzVector pR = p(1, evt) + p(2, evt);
  TLorentzVector qR = p(1, evt) - p(2, evt);
  TLorentzVector pD = pR + p(0, evt);
  TLorentzVector qD = pR - p(0, evt);

  double mR = mRes(R, evt);
  double mD = pD.M();

  ZTspin1 LR(qR, pR, mR);
  ZTspin1 LD(qD, pD, mD);
  return LD.Contract(LR)/(GeV*GeV);
}

double SpinFactor3::spinTwoFromZemach(IDalitzEvent& evt){
    TLorentzVector pR = p(1, evt) + p(2, evt);
    TLorentzVector qR = p(1, evt) - p(2, evt);
    TLorentzVector pD = pR + p(0, evt);
    TLorentzVector qD = pR - p(0, evt);
    
    double mR = mRes(R, evt);
    double mD = pD.M();
    
    ZTspin2 LR(qR, pR, mR);
    ZTspin2 LD(qD, pD, mD);
    return LD.Contract_2(LR)/(GeV*GeV*GeV*GeV);
}

double SpinFactor3::spinThreeFromZemach(IDalitzEvent& evt){
    TLorentzVector pR = p(1, evt) + p(2, evt);
    TLorentzVector qR = p(1, evt) - p(2, evt);
    TLorentzVector pD = pR + p(0, evt);
    TLorentzVector qD = pR - p(0, evt);
    
    double mR = mRes(R, evt);
    double mD = pD.M();
    
    ZTspin1 LR(qR, pR, mR);
    ZTspin1 LD(qD, pD, mD);
    
    SpinSumV PR(pR,mR);  
    SpinSumV PD(pD,mD);  

    double val = pow(LD.Contract(LR),3);
    val -= 3./5. * pR.M2() * LD.Contract(LR) * PR.Sandwich(LD,LD);
    val -= 3./5. * pD.M2() * LD.Contract(LR) * PD.Sandwich(LR,LR);
    val += 3./25. * pR.M2() * pD.M2() * LD.Contract(LR) * (2. + pow(pR.Dot(pD),2)/(pR.M2() * pD.M2()));
    val += 6./25. * pR.M2() * pD.M2() * (PD.Dot(LR)).Dot((PR.Dot(LD)));    

    return val/(GeV*GeV*GeV*GeV*GeV*GeV);
}



double SpinFactor3::spinOneFromMasses(IDalitzEvent& evt){
  // parsed as:
  // D -> V P0; V->P1, P2;
  // fsPS[0] && fsPS[1] && fsPS[2]
  
  double MV = mRes(R, evt);
  double mA = fsPS[1]->getVal().mass();
  double mB = fsPS[2]->getVal().mass();

  double mC = fsPS[0]->getVal().mass();

  double mD = (p(0, evt)+p(1, evt)+p(2, evt)).M();

  double m2AC = (p(0, evt) + p(1, evt)).M2();
  double m2BC = (p(0, evt) + p(2, evt)).M2();

  return (m2AC - m2BC + (mD*mD - mC*mC)*(mB*mB-mA*mA)/(MV*MV))/(GeV*GeV);

}

double SpinFactor3::spinTwoFromMasses(IDalitzEvent& evt){
  // parsed as:
  // D -> V P0; V->P1, P2;
  // fsPS[0] && fsPS[1] && fsPS[2]
  
  double MV = mRes(R, evt);
  double mA = fsPS[1]->getVal().mass();
  double mB = fsPS[2]->getVal().mass();

  double mC = fsPS[0]->getVal().mass();

  double mD = (p(0, evt)+p(1, evt)+p(2, evt)).M();

  double m2AB = (p(1, evt) + p(2, evt)).M2();
  double m2AC = (p(0, evt) + p(1, evt)).M2();
  double m2BC = (p(0, evt) + p(2, evt)).M2();
  
  Double_t term1 = m2BC-m2AC+(mD*mD-mC*mC)*(mA*mA-mB*mB)/(MV*MV);
  Double_t term2 = m2AB-2.0*mD*mD-2.0*mC*mC+pow(mD*mD-mC*mC,2)/(MV*MV);
  Double_t term3 = m2AB-2.0*mA*mA-2.0*mB*mB+pow(mA*mA-mB*mB,2)/(MV*MV);
  
  return (pow(term1,2)-(1.0/3.0)*term2*term3)/(GeV*GeV*GeV*GeV);
}
//
double SpinFactor3::spinTwoVal(IDalitzEvent& evt){
  // parsed as:
  // D -> T P0; T->P1, P2;
  bool dbThis=false;

  TLorentzVector pT = p(1, evt) + p(2, evt);
  TLorentzVector pD = pT + p(0, evt);
  double mr = mRes(R, evt);
  SpinSumT spin_sum(pT, mr);

  TLorentzVector lhs = pD + p(0, evt);
  TLorentzVector rhs = p(1, evt) - p(2, evt);

  double returnVal = spin_sum.Sandwich(lhs, lhs, rhs, rhs)/(GeV*GeV*GeV*GeV);
  if(dbThis){
    cout << "spin-2 spin factor got called for "
	 << theDecay().oneLiner()
	 << "\n   >  returning " << returnVal
	 << "\n   >  compare " << spinTwoFromMasses(evt)
	 << "\n   >  ratio " << returnVal/spinTwoFromMasses(evt)
	 << endl;
  }

  return returnVal;
}

//
