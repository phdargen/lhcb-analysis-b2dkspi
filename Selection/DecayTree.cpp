#define DecayTree_cxx
#include "DecayTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <ctime>
#include <math.h>
#include <sstream>
#include "TLorentzVector.h"

using namespace std; 

TTree* DecayTree::GetInputTree(){
    
    TString tupleName;
    if(_data==DataType::data) tupleName = "DecayTree";
    else tupleName = "B2DKspi_DD_Tuple/DecayTree";
    TChain* chain = new TChain(tupleName);
    
    TString fileName = _inFileLoc;
    stringstream ss_year;
    ss_year << _year;
    string str_year = ss_year.str();
    
    if(_data==DataType::data) {
        //fileName+= "Data/";
        fileName+= str_year; 
    }
    else {
        //fileName+= "MC/";
        fileName+= "bkg/";
        //fileName+= str_year; 
        //fileName+= "U";
    }
    //fileName+= "/b2dkspi*.root"; 
    fileName+= "/bkg_bs2dstarkspi.root"; 
    
    cout << "Using the files: " << endl;
    cout << fileName << endl << endl;
    chain->Add(fileName);

    if(chain->GetEntries()==0){
        cout << "ERROR: No events found" << endl;
        throw "ERROR";
    }
    
    return (TTree*)chain;
}

DecayTree::DecayTree(Decay::Type decay, Year::Type year, DataType::Type dataType, TString polarity, TString inFileLoc, TString outFileLoc ) : 
fChain(0), _decay(decay), _year(year), _data(dataType), _polarity(polarity), _inFileLoc(inFileLoc), _outFileLoc(outFileLoc)
{    
    cout << "Requested to process files with options: " << endl << endl;
    
    TString s1,s2,s3,s4;
    if(_data==DataType::data)s1="Data";
    else s1 = "MC";
    if(_decay==Decay::LL)s2="LL";
    else s2 = "DD";
    stringstream ss_year;
    ss_year << _year;
    string str_year = ss_year.str();
    
    cout << "DataType: " << s1 << endl;
    cout << "Decay: " << s2 << endl;
    cout << "Year: " << _year << endl;
    cout << "Polarity: " << _polarity  << endl;
    
    _outFileName = _outFileLoc;   
    _outFileName += s1;  
    _outFileName += "_";  
    //_outFileName += "b2dkspi_";
    _outFileName += "bkg_bs2dstarkspi_";
    _outFileName += s2;   
    _outFileName += "_";  
    _outFileName += str_year;
    if(_polarity == "Up") _outFileName += "_up";
    if(_polarity == "Down") _outFileName += "_down";
    _outFileName += ".root";   
    
    cout << "Producing new file " << _outFileName << endl;    
}

inline Bool_t DecayTree::TriggerCuts(Long64_t i){
    
    b_B_L0Global_TIS->GetEntry(i);
    b_B_L0HadronDecision_TOS->GetEntry(i);
    if( (!B_L0Global_TIS) && (!B_L0HadronDecision_TOS)) return false;
    
    if(_year == 15 || _year == 16 || _year == 17 || _year == 18){
        b_B_Hlt1TrackMVADecision_TOS->GetEntry(i);
        b_B_Hlt1TwoTrackMVADecision_TOS->GetEntry(i);
        if((!B_Hlt1TrackMVADecision_TOS) && (!B_Hlt1TwoTrackMVADecision_TOS) ) return false;
        
        b_B_Hlt2Topo2BodyDecision_TOS->GetEntry(i);
        b_B_Hlt2Topo3BodyDecision_TOS->GetEntry(i);
        b_B_Hlt2Topo4BodyDecision_TOS->GetEntry(i);
        if(_year == 15){
            b_B_Hlt2IncPhiDecision_TOS->GetEntry(i);
            if((!B_Hlt2Topo2BodyDecision_TOS) &&  (!B_Hlt2Topo3BodyDecision_TOS) && (!B_Hlt2Topo4BodyDecision_TOS) && (!B_Hlt2IncPhiDecision_TOS) ) return false;    
        }
        else {
            b_B_Hlt2PhiIncPhiDecision_TOS->GetEntry(i);
            if((!B_Hlt2Topo2BodyDecision_TOS) &&  (!B_Hlt2Topo3BodyDecision_TOS) && (!B_Hlt2Topo4BodyDecision_TOS) && (!B_Hlt2PhiIncPhiDecision_TOS) ) return false;
        }
    }
  /*  
    else if(_year == 11 || _year == 12){
        b_B_Hlt1TrackAllL0Decision_TOS->GetEntry(i);        
        if(!B_Hlt1TrackAllL0Decision_TOS) return false;
        
        b_B_Hlt2Topo2BodyBBDTDecision_TOS->GetEntry(i);
        b_B_Hlt2Topo3BodyBBDTDecision_TOS->GetEntry(i);
        b_B_Hlt2Topo4BodyBBDTDecision_TOS->GetEntry(i);
        b_B_Hlt2IncPhiDecision_TOS->GetEntry(i);
        if((!B_Hlt2Topo2BodyBBDTDecision_TOS) &&  (!B_Hlt2Topo3BodyBBDTDecision_TOS) && (!B_Hlt2Topo4BodyBBDTDecision_TOS) 
           && (!B_Hlt2IncPhiDecision_TOS)) return false;
    }
    */
    return true;
}

inline Bool_t DecayTree::LooseCuts(Long64_t i){
    
    b_B_DIRA_OWNPV->GetEntry(i);
    if(B_DIRA_OWNPV<0.99994) return false;
    
    b_B_IPCHI2_OWNPV->GetEntry(i);
    if(B_IPCHI2_OWNPV>20) return false;
    
    b_B_FDCHI2_OWNPV->GetEntry(i);
    if(B_FDCHI2_OWNPV<100) return false;
    
    b_B_ENDVERTEX_CHI2->GetEntry(i);
    b_B_ENDVERTEX_NDOF->GetEntry(i);
    if((B_ENDVERTEX_CHI2/B_ENDVERTEX_NDOF)> 8) return false;
    
    //     b_B_TAU->GetEntry(i);
    //     if(B_TAU < 0.0002) return false;
    
    b_B_MM->GetEntry(i);
    if(B_MM < 4800. || B_MM > 6000.) return false;
    
    b_B_DTF_M->GetEntry(i);
    if(B_DTF_M[0] < 4800. || B_DTF_M[0] > 6000.) return false;
        
    b_B_PV_M->GetEntry(i);
    if (B_PV_M[0] < 4800. || B_PV_M[0] > 6000.) return false;
    
    b_D_ENDVERTEX_Z->GetEntry(i);
    b_B_ENDVERTEX_Z->GetEntry(i);
    if((D_ENDVERTEX_Z - B_ENDVERTEX_Z) < -1) return false;
    
    //b_D_MM->GetEntry(i);
    //if(fabs(D_MM-massDminus) > 20 ) return false;
    
    //b_Ks_MM->GetEntry(i);
    //if(fabs(Ks_MM-massKs) > 20 ) return false;
    
    b_D_FDCHI2_ORIVX->GetEntry(i);
    if(D_FDCHI2_ORIVX < 0) return false; 
    
    b_K_D_ProbNNk->GetEntry(i);
    if(K_D_ProbNNk < 0.15) return false;

    b_pi1_D_ProbNNpi->GetEntry(i);
    if(pi1_D_ProbNNpi < 0.1) return false;

    b_pi2_D_ProbNNpi->GetEntry(i);
    if(pi2_D_ProbNNpi < 0.1) return false;
    
    b_pi_ProbNNpi->GetEntry(i);
    if(pi_ProbNNpi < 0.15) return false;
    
    b_pip_Ks_ProbNNpi->GetEntry(i);
    if(pip_Ks_ProbNNpi < 0.1) return false;
    
    b_pim_Ks_ProbNNpi->GetEntry(i);
    if(pim_Ks_ProbNNpi < 0.1) return false;

    b_pi_isMuon->GetEntry(i);
    if(pi_isMuon == true) return false;
    
    b_pi_hasRich->GetEntry(i);
    if(pi_hasRich == false) return false;
    
    b_K_D_hasRich->GetEntry(i);
    if(K_D_hasRich == false) return false;

    b_D_DIRA_OWNPV->GetEntry(i);
    if(D_DIRA_OWNPV<0.) return false;
    
    b_Ks_DIRA_OWNPV->GetEntry(i);
    if(Ks_DIRA_OWNPV<0.) return false;

    b_Ks_PT->GetEntry(i);
    if(Ks_PT< 500.) return false;

    
    return true;
}


inline void DecayTree::set_LorentzVectors(){
    
        pi.SetXYZM(pi_PX,pi_PY,pi_PZ,massPion);        
       
        pip_Ks.SetXYZM(pip_Ks_PX,pip_Ks_PY,pip_Ks_PZ,massPion);
        pim_Ks.SetXYZM(pim_Ks_PX,pim_Ks_PY,pim_Ks_PZ,massPion);
        
        pi1_D.SetXYZM(pi1_D_PX,pi1_D_PY,pi1_D_PZ,massPion);
        pi2_D.SetXYZM(pi2_D_PX,pi2_D_PY,pi2_D_PZ,massPion);
        K_D.SetXYZM(K_D_PX,K_D_PY,K_D_PZ,massKaon);

        Ks = pip_Ks + pim_Ks;
        D = pi1_D + pi2_D + K_D;
    
FullDTF_pi1_D=TLorentzVector(B_FullDTF_Dplus_piplus_PX[0],B_FullDTF_Dplus_piplus_PY[0],B_FullDTF_Dplus_piplus_PZ[0],B_FullDTF_Dplus_piplus_PE[0]);  
FullDTF_pi2_D=TLorentzVector(B_FullDTF_Dplus_piplus_0_PX[0],B_FullDTF_Dplus_piplus_0_PY[0],B_FullDTF_Dplus_piplus_0_PZ[0],B_FullDTF_Dplus_piplus_0_PE[0]);  
FullDTF_K_D=TLorentzVector(B_FullDTF_Dplus_Kplus_PX[0],B_FullDTF_Dplus_Kplus_PY[0],B_FullDTF_Dplus_Kplus_PZ[0],B_FullDTF_Dplus_Kplus_PE[0]);  

FullDTF_pi=TLorentzVector(B_FullDTF_piplus_PX[0],B_FullDTF_piplus_PY[0],B_FullDTF_piplus_PZ[0],B_FullDTF_piplus_PE[0]);  

FullDTF_pip_Ks=TLorentzVector(B_FullDTF_KS0_piplus_PX[0],B_FullDTF_KS0_piplus_PY[0],B_FullDTF_KS0_piplus_PZ[0],B_FullDTF_KS0_piplus_PE[0]);  
FullDTF_pim_Ks=TLorentzVector(B_FullDTF_KS0_piplus_0_PX[0],B_FullDTF_KS0_piplus_0_PY[0],B_FullDTF_KS0_piplus_0_PZ[0],B_FullDTF_KS0_piplus_0_PE[0]);  

        FullDTF_D = FullDTF_pi1_D + FullDTF_K_D + FullDTF_pi2_D;
        FullDTF_Ks = FullDTF_pip_Ks + FullDTF_pim_Ks;
    
        K_fromD_asP_MissID.SetXYZM(K_D_PX,K_D_PY,K_D_PZ, massProton);
        K_fromD_asPi_MissID.SetXYZM(K_D_PX,K_D_PY,K_D_PZ,massPion);
    
        pi_asP_MissID.SetXYZM(pi_PX,pi_PY,pi_PZ,massProton);        
        pi_asK_MissID.SetXYZM(pi_PX,pi_PY,pi_PZ,massKaon);        
    
        pi1_fromD_asP_MissID.SetXYZM(pi_PX,pi_PY,pi_PZ,massProton);        
        pi1_fromD_asK_MissID.SetXYZM(pi_PX,pi_PY,pi_PZ,massKaon);        
    
        pi2_fromD_asP_MissID.SetXYZM(pi_PX,pi_PY,pi_PZ,massProton);        
        pi2_fromD_asK_MissID.SetXYZM(pi_PX,pi_PY,pi_PZ,massKaon);        
}
    
    
void DecayTree::Loop()
{
    
    time_t startTime = time(0);
    
    Init();
    if (fChain == 0) return;
    
    fChain->SetBranchStatus("*",0);  // disable all branches
    
    // activate branchname
    fChain->SetBranchStatus("*ID",1);  
    
    fChain->SetBranchStatus("B_*_T*S",1);  
    fChain->SetBranchStatus("B_*Muon*_T*S",0);  
    fChain->SetBranchStatus("B_*Hlt*Phys*_T*S",0);  
    fChain->SetBranchStatus("B_*Hlt*Global*_T*S",0);  
    fChain->SetBranchStatus("B_*DTF*",1);  
    fChain->SetBranchStatus("B_*PV*",1);  
    
    fChain->SetBranchStatus("B*ENDVERTEX*",1);  
    fChain->SetBranchStatus("B*OWNPV*",1);  
    fChain->SetBranchStatus("B*ENDVERTEX_COV*",0);  
    fChain->SetBranchStatus("B*OWNPV_COV*",0);  
    
    fChain->SetBranchStatus("D_M",1);         
    fChain->SetBranchStatus("D*ENDVERTEX*",1);  
    fChain->SetBranchStatus("D*OWNPV*",1);  
    fChain->SetBranchStatus("D*ENDVERTEX_COV*",0);  
    fChain->SetBranchStatus("D*OWNPV_COV*",0);  
    fChain->SetBranchStatus("D*ORIVX*",1);  
    fChain->SetBranchStatus("D*ORIVX_COV*",0);  

    fChain->SetBranchStatus("Ks_M",1);         
    fChain->SetBranchStatus("Ks*ENDVERTEX*",1);  
    fChain->SetBranchStatus("Ks*OWNPV*",1);  
    fChain->SetBranchStatus("Ks*ENDVERTEX_COV*",0);  
    fChain->SetBranchStatus("Ks*OWNPV_COV*",0);  
    fChain->SetBranchStatus("Ks*ORIVX*",1);  
    fChain->SetBranchStatus("Ks*ORIVX_COV*",0);  

    fChain->SetBranchStatus("*IP*",1);  
    fChain->SetBranchStatus("*IPCHI2*",1);  
    fChain->SetBranchStatus("*FD*",1);  
    fChain->SetBranchStatus("*FDCHI2*",1);  
    fChain->SetBranchStatus("*P",1);  
    fChain->SetBranchStatus("*PT",1);  
    fChain->SetBranchStatus("*PE",1);  
    fChain->SetBranchStatus("*PX",1);  
    fChain->SetBranchStatus("*PY",1);  
    fChain->SetBranchStatus("*PZ",1);  
    fChain->SetBranchStatus("*ETA",1);  
    fChain->SetBranchStatus("*MM*",1);  
    fChain->SetBranchStatus("*TAU*",1);  
    fChain->SetBranchStatus("*ptasy_1.00",1);  
    
    fChain->SetBranchStatus("*DIRA*",1);  
    fChain->SetBranchStatus("D_DOCA*",1);  
    
    fChain->SetBranchStatus("*PID*",1);  
    fChain->SetBranchStatus("*PIDe*",0);  
    fChain->SetBranchStatus("*ProbNN*",1);  
    fChain->SetBranchStatus("*ProbNNe*",0);  
    fChain->SetBranchStatus("*TRACK_Ghost*",1);  
    fChain->SetBranchStatus("*TRACK_CHI2*",1);  
    fChain->SetBranchStatus("*isMuon*",1);  
    fChain->SetBranchStatus("*hasRich",1);  
    
    fChain->SetBranchStatus("nCandidate",1) ;
    fChain->SetBranchStatus("nTracks",1) ;
    fChain->SetBranchStatus("nPV",1) ;
    fChain->SetBranchStatus("eventNumber",1) ;
    fChain->SetBranchStatus("runNumber",1) ;
    fChain->SetBranchStatus("EventInSequence",1) ;
    fChain->SetBranchStatus("totCandidates",1) ;
    fChain->SetBranchStatus("Polarity",1) ;

    fChain->SetBranchStatus("*TOP*",0) ;
    fChain->SetBranchStatus("*All*",0) ;
    fChain->SetBranchStatus("*TuneV*",0) ;
    
    if(!_data){
        fChain->SetBranchStatus("*TRUE*",1) ;
        fChain->SetBranchStatus("*BKG*",1) ;
    }
    
    TFile* output = new TFile(_outFileName,"RECREATE");
    fChain->SetBranchStatus("*ptasy_1.00",0);  
    fChain->SetBranchStatus("*DTF*",0);  
    fChain->SetBranchStatus("*PV*",0);  
    fChain->SetBranchStatus("*OWNPV*",1);  
    TTree* summary_tree = fChain->CloneTree(0);
    fChain->SetBranchStatus("*ptasy_1.00",1);  
    fChain->SetBranchStatus("*DTF*",1);  
    fChain->SetBranchStatus("*PV*",1);  

    // Dalitz stuff
    double m_D_Kpi,m_D_Kpi2,m_D_pipi;
    summary_tree->Branch("m_D_Kpi", &m_D_Kpi, "m_D_Kpi/D");
    summary_tree->Branch("m_D_Kpi2", &m_D_Kpi2, "m_D_Kpi2/D");
    summary_tree->Branch("m_D_pipi", &m_D_pipi, "m_D_pipi/D");
    
    double bkg_pipipi_as_D,bkg_ppipi_as_D;
    summary_tree->Branch("bkg_pipipi_as_D", &bkg_pipipi_as_D, "bkg_pipipi_as_D/D");
    summary_tree->Branch("bkg_ppipi_as_D", &bkg_ppipi_as_D, "bkg_ppipi_as_D/D");

    double bkg_KKpi_as_D,bkg_Kppi_as_D,bkg_KKpi2_as_D,bkg_Kppi2_as_D;
    summary_tree->Branch("bkg_KKpi_as_D", &bkg_KKpi_as_D, "bkg_KKpi_as_D/D");
    summary_tree->Branch("bkg_Kppi_as_D", &bkg_Kppi_as_D, "bkg_Kppi_as_D/D");
    summary_tree->Branch("bkg_KKpi2_as_D", &bkg_KKpi2_as_D, "bkg_KKpi2_as_D/D");
    summary_tree->Branch("bkg_Kppi2_as_D", &bkg_Kppi2_as_D, "bkg_Kppi2_as_D/D");
 
    double bkg_DKsK_as_B,bkg_DKsp_as_B;
    summary_tree->Branch("bkg_DKsK_as_B", &bkg_DKsK_as_B, "bkg_DKsK_as_B/D");
    summary_tree->Branch("bkg_DKsp_as_B", &bkg_DKsp_as_B, "bkg_DKsp_as_B/D");

    double m_DKs,m_Dpi,m_Kspi;
    summary_tree->Branch("m_DKs", &m_DKs, "m_DKs/D");
    summary_tree->Branch("m_Dpi", &m_Dpi, "m_Dpi/D");
    summary_tree->Branch("m_Kspi", &m_Kspi, "m_Kspi/D");
    
    double FullDTF_m_DKs,FullDTF_m_Dpi,FullDTF_m_Kspi;
    summary_tree->Branch("FullDTF_m_DKs", &FullDTF_m_DKs, "FullDTF_m_DKs/D");
    summary_tree->Branch("FullDTF_m_Dpi", &FullDTF_m_Dpi, "FullDTF_m_Dpi/D");
    summary_tree->Branch("FullDTF_m_Kspi", &FullDTF_m_Kspi, "FullDTF_m_Kspi/D");
    
    int B_BKGCAT = 0;
    if(_data)summary_tree->Branch("B_BKGCAT", &B_BKGCAT, "B_BKGCAT/I");

    int KsCat, TriggerCat,run;
    int year = (int)_year;
    if(year == 11 || year == 12) run = 1;
    else run = 2;
    if(_decay==Decay::LL) KsCat = 0;
    else KsCat = 1;
    summary_tree->Branch("KsCat",&KsCat,"KsCat/I");    
    summary_tree->Branch("year", &year, "year/I");
    summary_tree->Branch("run", &run, "run/I");
    summary_tree->Branch("TriggerCat", &TriggerCat, "TriggerCat/I");
    
    double weight = 1.;
    summary_tree->Branch("weight",&weight,"weight/D"); 
    
    // Variables for BDT 
    double DDaughters_min_IPCHI2 = 0;
    double KsDaughters_min_IPCHI2 = 0;
    double DDaughters_max_IPCHI2 = 0;
    double KsDaughters_max_IPCHI2 = 0;
    double DDaughters_min_PT = 0;
    double KsDaughters_min_PT = 0;
    double B_max_DOCA = 0;
    double Ks_max_DOCA = 0;
    double D_max_DOCA = 0;
    double Ks_max_ghostProb = 0;
    double D_max_ghostProb = 0;
    double max_ghostProb = 0;
    double track_min_PT,track_min_P,track_min_IPCHI2,track_max_IPCHI2,Ks_ptasy,D_ptasy,B_ptasy;    
    
    summary_tree->Branch("DDaughters_min_IPCHI2",&DDaughters_min_IPCHI2,"DDaughters_min_IPCHI2/D");
    summary_tree->Branch("KsDaughters_min_IPCHI2",&KsDaughters_min_IPCHI2,"KsDaughters_min_IPCHI2/D");
    summary_tree->Branch("DDaughters_max_IPCHI2",&DDaughters_max_IPCHI2,"DDaughters_max_IPCHI2/D");
    summary_tree->Branch("KsDaughters_max_IPCHI2",&KsDaughters_max_IPCHI2,"KsDaughters_max_IPCHI2/D");
    summary_tree->Branch("track_min_IPCHI2",&track_min_IPCHI2,"track_min_IPCHI2/D");
    summary_tree->Branch("DDaughters_min_PT",&DDaughters_min_PT,"DDaughters_min_PT/D");
    summary_tree->Branch("KsDaughters_min_PT",&KsDaughters_min_PT,"KsDaughters_min_PT/D");
    summary_tree->Branch("track_min_PT",&track_min_PT,"track_min_PT/D");
    summary_tree->Branch("track_min_P",&track_min_P,"track_min_P/D");

    summary_tree->Branch("D_ptasy",&D_ptasy,"D_ptasy/D");
    summary_tree->Branch("Ks_ptasy",&Ks_ptasy,"Ks_ptasy/D");
    summary_tree->Branch("B_ptasy",&B_ptasy,"B_ptasy/D");
    
    summary_tree->Branch("Ks_max_DOCA",&Ks_max_DOCA,"Ks_max_DOCA/D");
    summary_tree->Branch("D_max_DOCA",&D_max_DOCA,"D_max_DOCA/D");
    summary_tree->Branch("Ks_max_ghostProb",&Ks_max_ghostProb,"Ks_max_ghostProb/D");
    summary_tree->Branch("D_max_ghostProb",&D_max_ghostProb,"D_max_ghostProb/D");
    summary_tree->Branch("max_ghostProb",&max_ghostProb,"max_ghostProb/D");
    
    double angKs,angPi,angPiKs, maxCos, maxCos2;
    summary_tree->Branch("angKs", &angKs, "angKs/D");
    summary_tree->Branch("angPi", &angPi, "angPi/D");
    summary_tree->Branch("angPiKs", &angPiKs, "angPiKs/D");
    summary_tree->Branch("maxCos", &maxCos, "maxCos/D");
    summary_tree->Branch("maxCos2", &maxCos2, "maxCos2/D");
    
    double B_RFD,D_RFD,D_FDsig,D_z,Ks_RFD,Ks_FDsig,Ks_z;
    summary_tree->Branch("B_RFD", &B_RFD, "B_RFD/D");
    summary_tree->Branch("D_RFD", &D_RFD, "D_RFD/D");
    summary_tree->Branch("D_FDsig", &D_FDsig, "D_FDsig/D");
    summary_tree->Branch("D_z", &D_z, "D_z/D");
    summary_tree->Branch("Ks_RFD", &Ks_RFD, "Ks_RFD/D");
    summary_tree->Branch("Ks_FDsig", &Ks_FDsig, "Ks_FDsig/D");
    summary_tree->Branch("Ks_z", &Ks_z, "Ks_z/D");
    
    // DTF stuff
    double PV_status, DTF_status, FullDTF_status;
    double PV_CHI2NDOF, DTF_CHI2NDOF, FullDTF_CHI2NDOF;
    summary_tree->Branch("PV_status", &PV_status, "PV_status/D");
    summary_tree->Branch("DTF_status", &DTF_status, "DTF_status/D");
    summary_tree->Branch("FullDTF_status", &FullDTF_status, "FullDTF_status/D");
    summary_tree->Branch("PV_CHI2NDOF", &PV_CHI2NDOF, "PV_CHI2NDOF/D");
    summary_tree->Branch("DTF_CHI2NDOF", &DTF_CHI2NDOF, "DTF_CHI2NDOF/D");
    summary_tree->Branch("FullDTF_CHI2NDOF", &FullDTF_CHI2NDOF, "FullDTF_CHI2NDOF/D");
    
    double B_PV_TAU,B_PV_TAUERR;
    double B_DTF_TAU,B_DTF_TAUERR;
    double B_FullDTF_TAU, B_FullDTF_TAUERR;
    summary_tree->Branch("B_PV_TAU", &B_PV_TAU, "B_PV_TAU/D");
    summary_tree->Branch("B_PV_TAUERR", &B_PV_TAUERR, "B_PV_TAUERR/D");
    summary_tree->Branch("B_DTF_TAU", &B_DTF_TAU, "B_DTF_TAU/D");
    summary_tree->Branch("B_DTF_TAUERR", &B_DTF_TAUERR, "B_DTF_TAUERR/D");
    summary_tree->Branch("B_FullDTF_TAU", &B_FullDTF_TAU, "B_FullDTF_TAU/D");
    summary_tree->Branch("B_FullDTF_TAUERR", &B_FullDTF_TAUERR, "B_FullDTF_TAUERR/D");
    
    double B_DTF_MM,B_DTF_MMERR;
    summary_tree->Branch("B_DTF_MM", &B_DTF_MM, "B_DTF_MM/D");
    summary_tree->Branch("B_DTF_MMERR", &B_DTF_MMERR, "B_DTF_MMERR/D");
    double B_PV_MM,B_PV_MMERR;
    summary_tree->Branch("B_PV_MM", &B_PV_MM, "B_PV_MM/D");
    summary_tree->Branch("B_PV_MMERR", &B_PV_MMERR, "B_PV_MMERR/D");    
    double D_PV_MM,D_PV_MMERR;
    summary_tree->Branch("D_PV_MM", &D_PV_MM, "D_PV_MM/D");
    summary_tree->Branch("D_PV_MMERR", &D_PV_MMERR, "D_PV_MMERR/D");    
    double Ks_PV_MM,Ks_PV_MMERR;
    summary_tree->Branch("Ks_PV_MM", &Ks_PV_MM, "Ks_PV_MM/D");
    summary_tree->Branch("Ks_PV_MMERR", &Ks_PV_MMERR, "Ks_PV_MMERR/D");    
    
    double FullDTF_D_PX,FullDTF_D_PY,FullDTF_D_PZ,FullDTF_D_PE;
    summary_tree->Branch("FullDTF_D_PX", &FullDTF_D_PX, "FullDTF_D_PX/D");
    summary_tree->Branch("FullDTF_D_PY", &FullDTF_D_PY, "FullDTF_D_PY/D");
    summary_tree->Branch("FullDTF_D_PZ", &FullDTF_D_PZ, "FullDTF_D_PZ/D");
    summary_tree->Branch("FullDTF_D_PE", &FullDTF_D_PE, "FullDTF_D_PE/D");   

    double FullDTF_Ks_PX,FullDTF_Ks_PY,FullDTF_Ks_PZ,FullDTF_Ks_PE;
    summary_tree->Branch("FullDTF_Ks_PX", &FullDTF_Ks_PX, "FullDTF_Ks_PX/D");
    summary_tree->Branch("FullDTF_Ks_PY", &FullDTF_Ks_PY, "FullDTF_Ks_PY/D");
    summary_tree->Branch("FullDTF_Ks_PZ", &FullDTF_Ks_PZ, "FullDTF_Ks_PZ/D");
    summary_tree->Branch("FullDTF_Ks_PE", &FullDTF_Ks_PE, "FullDTF_Ks_PE/D");   

    double FullDTF_pi_PX,FullDTF_pi_PY,FullDTF_pi_PZ,FullDTF_pi_PE;
    summary_tree->Branch("FullDTF_pi_PX", &FullDTF_pi_PX, "FullDTF_pi_PX/D");
    summary_tree->Branch("FullDTF_pi_PY", &FullDTF_pi_PY, "FullDTF_pi_PY/D");
    summary_tree->Branch("FullDTF_pi_PZ", &FullDTF_pi_PZ, "FullDTF_pi_PZ/D");
    summary_tree->Branch("FullDTF_pi_PE", &FullDTF_pi_PE, "FullDTF_pi_PE/D");   

    
    Long64_t nentries = fChain->GetEntries();
    cout << "Have " << nentries << " events" <<  endl << endl;
    
    for (Long64_t i=0; i<nentries;i++) {
        
        if(0ul == (i % 100000ul)) cout << "Read event " << i << "/" << nentries <<
            "  ( " << i/(double)nentries * 100. << " % )" << endl;
        
        // Read from individual branches rather than whole tree,
        // messy and prone to errors but benefical to performance
        // fChain->GetEntry(i);   
        
        Long64_t j = LoadTree(i);
        if (j < 0) break;
        
        if(!TriggerCuts(j)) continue;
        else if(!LooseCuts(j)) continue;
        
        fChain->GetEntry(i);   
        set_LorentzVectors();
    
        if(B_L0Global_TIS)TriggerCat = 2;
        else TriggerCat = 3;
        if(B_L0HadronDecision_TOS)TriggerCat = 0;
        else TriggerCat = 1;
        
        B_TAU = B_TAU*1000.;
        B_TAUERR = B_TAUERR*1000.;
        B_PV_TAU = B_PV_ctau[0] * 3.33564095; // 1e-3*1/(2.99792458*1e8)*1e12
        B_PV_TAUERR = B_PV_ctauErr[0] * 3.33564095;
        B_DTF_TAU = B_DTF_ctau[0] * 3.33564095;
        B_DTF_TAUERR = B_DTF_ctauErr[0] * 3.33564095;
        B_FullDTF_TAU = B_FullDTF_ctau[0] * 3.33564095;
        B_FullDTF_TAUERR = B_FullDTF_ctauErr[0] * 3.33564095;
        
        //if(B_FullDTF_TAU < 0.4 || B_FullDTF_TAU > 10.) continue;
        if(B_FullDTF_TAUERR < 0. || B_FullDTF_TAUERR > 0.5) continue;

        B_RFD = sqrt(pow(B_ENDVERTEX_X-B_OWNPV_X,2)+pow(B_ENDVERTEX_Y-B_OWNPV_Y,2));
        D_RFD = sqrt(pow(D_ENDVERTEX_X-D_OWNPV_X,2)+pow(D_ENDVERTEX_Y-D_OWNPV_Y,2));
        D_FDsig = (D_ENDVERTEX_Z-D_ORIVX_Z)/sqrt(pow(D_ENDVERTEX_ZERR,2)+pow(D_ORIVX_ZERR,2));
        D_z = D_ENDVERTEX_Z - B_ENDVERTEX_Z;
      
        Ks_RFD = sqrt(pow(Ks_ENDVERTEX_X-Ks_OWNPV_X,2)+pow(Ks_ENDVERTEX_Y-Ks_OWNPV_Y,2));
        Ks_FDsig = (Ks_ENDVERTEX_Z-Ks_ORIVX_Z)/sqrt(pow(Ks_ENDVERTEX_ZERR,2)+pow(Ks_ORIVX_ZERR,2));
        Ks_z = Ks_ENDVERTEX_Z - B_ENDVERTEX_Z;

        PV_status = B_PV_status[0];
        DTF_status = B_DTF_status[0];
        FullDTF_status = B_FullDTF_status[0];
        PV_CHI2NDOF = B_PV_chi2[0]/B_PV_nDOF[0];
        DTF_CHI2NDOF = B_DTF_chi2[0]/B_DTF_nDOF[0];
        FullDTF_CHI2NDOF = B_FullDTF_chi2[0]/B_FullDTF_nDOF[0];
        
        B_PV_MM = B_PV_M[0];
        B_DTF_MM = B_DTF_M[0];
        D_PV_MM = B_PV_Dplus_M[0];
        Ks_PV_MM = B_PV_KS0_M[0];
        
        B_PV_MMERR = B_PV_MERR[0];
        B_DTF_MMERR = B_DTF_MERR[0];
        D_PV_MMERR = B_PV_Dplus_MERR[0];
        Ks_PV_MMERR = B_PV_KS0_MERR[0];
                
        FullDTF_D_PX = FullDTF_D.Px() ;
        FullDTF_D_PY = FullDTF_D.Py() ;
        FullDTF_D_PZ = FullDTF_D.Pz() ;
        FullDTF_D_PE = FullDTF_D.E() ;
        
        FullDTF_Ks_PX = FullDTF_Ks.Px() ;
        FullDTF_Ks_PY = FullDTF_Ks.Py() ;
        FullDTF_Ks_PZ = FullDTF_Ks.Pz() ;
        FullDTF_Ks_PE = FullDTF_Ks.E() ;
        
        FullDTF_pi_PX = FullDTF_pi.Px() ;
        FullDTF_pi_PY = FullDTF_pi.Py() ;
        FullDTF_pi_PZ = FullDTF_pi.Pz() ;
        FullDTF_pi_PE = FullDTF_pi.E() ;

        m_DKs = (D+Ks).M();
        m_Dpi = (D+pi).M();
        m_Kspi = (pi+Ks).M();
        
        FullDTF_m_DKs = (FullDTF_D+FullDTF_Ks).M();
        FullDTF_m_Dpi = (FullDTF_D+FullDTF_pi).M();
        FullDTF_m_Kspi = (FullDTF_pi+FullDTF_Ks).M();
            
        m_D_Kpi = (K_D+pi1_D).M();
        m_D_Kpi2 = (K_D+pi2_D).M();
        m_D_pipi = (pi2_D+pi1_D).M();
                
        bkg_DKsK_as_B = (D+Ks+pi_asK_MissID).M();
        bkg_DKsp_as_B = (D+Ks+pi_asP_MissID).M();
        
        bkg_pipipi_as_D = (K_fromD_asPi_MissID + pi1_D + pi2_D).M();
        bkg_ppipi_as_D = (K_fromD_asP_MissID + pi1_D + pi2_D).M();
        bkg_KKpi_as_D = (K_D + pi1_fromD_asK_MissID +pi2_D).M();
        bkg_Kppi_as_D = (K_D + pi1_fromD_asP_MissID +pi2_D).M();
        bkg_KKpi2_as_D = (K_D + pi2_fromD_asK_MissID +pi1_D).M();
        bkg_Kppi2_as_D = (K_D + pi2_fromD_asP_MissID +pi1_D).M();
        
        TVector3 v_D(D_PX,D_PY,0.);
        TVector3 v_Ks(Ks_PX,Ks_PY,0.);
        TVector3 v_pi(pi_PX,pi_PY,0.);
        angKs= v_D.Angle(v_Ks);
        angPi= v_D.Angle(v_pi);
        angPiKs= v_Ks.Angle(v_pi);
        maxCos = cos(max(angKs,angPi));
        maxCos2 = cos(max(max(angPiKs,angPi),angKs));
        
        if(maxCos2<-0.95)continue;
        
        KsDaughters_min_IPCHI2 = min(pip_Ks_IPCHI2_OWNPV,pim_Ks_IPCHI2_OWNPV);            
        KsDaughters_max_IPCHI2 = max(pip_Ks_IPCHI2_OWNPV,pim_Ks_IPCHI2_OWNPV);            
        KsDaughters_min_PT = min(pip_Ks_PT,pim_Ks_PT);          
        //Ks_max_DOCA = max(Ks_DOCA1,Ks_DOCA2);
        Ks_max_ghostProb = max(pip_Ks_TRACK_GhostProb,pim_Ks_TRACK_GhostProb);   
        Ks_ptasy = Ks_ptasy_1_00;

        DDaughters_min_IPCHI2 = min(pi1_D_IPCHI2_OWNPV,min(pi2_D_IPCHI2_OWNPV,K_D_IPCHI2_OWNPV));            
        DDaughters_max_IPCHI2 = max(pi1_D_IPCHI2_OWNPV,max(pi2_D_IPCHI2_OWNPV,K_D_IPCHI2_OWNPV));            
        DDaughters_min_PT = min(pi1_D_PT,min(pi2_D_PT,K_D_PT));            
        //D_max_DOCA = max(D_DOCA1,max(D_DOCA2,D_DOCA3));
        D_max_ghostProb = max(pi1_D_TRACK_GhostProb,max(K_D_TRACK_GhostProb,pi2_D_TRACK_GhostProb));  
        D_ptasy = D_ptasy_1_00;

        track_min_IPCHI2 = min(pi_IPCHI2_OWNPV,min(KsDaughters_min_IPCHI2,DDaughters_min_IPCHI2));            
        track_max_IPCHI2 = max(pi_IPCHI2_OWNPV,max(KsDaughters_max_IPCHI2,DDaughters_max_IPCHI2));            
        track_min_PT = min(pi_PT,min(KsDaughters_min_PT,DDaughters_min_PT));          
        //B_max_DOCA = max(B_DOCA1,max(B_DOCA2,B_DOCA3));
        max_ghostProb = max(pi_TRACK_GhostProb,max(Ks_max_ghostProb,D_max_ghostProb));   
        B_ptasy = B_ptasy_1_00;
        
        summary_tree->Fill();
        if(0ul == (i % 100000ul))summary_tree->AutoSave();
    }
    
    cout << "Selected " << summary_tree->GetEntries() << " events" <<  endl;
    cout << "Efficiency = " << summary_tree->GetEntries()/(double)nentries * 100. << " %" <<  endl;
    
    summary_tree->Write();
    output->Close();
    
    cout << "Created new file " << _outFileName << endl;
    
    cout << "==============================================" << endl;
    cout << " Done. " << " Total time since start " << (time(0) - startTime)/60.0 << " min." << endl;
    cout << "==============================================" << endl;
    
}
