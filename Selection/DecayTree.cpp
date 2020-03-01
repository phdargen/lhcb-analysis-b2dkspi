#define DecayTree_cxx
#include "DecayTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <ctime>
#include <math.h>
#include <sstream>

using namespace std; 

TTree* DecayTree::GetInputTree(){
    
    TString tupleName = "B2DKspi_Tuple/DecayTree";
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
        fileName+= "MC/";
        fileName+= str_year; 
        //fileName+= "U";
    }
    fileName+= "/b2dkspi*.root"; 
    
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
    _outFileName += "/";  
    _outFileName += "b2dkspi_";
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
    //if(fabs(D_MM-massDminus) > 60 ) return false;
    
    b_D_FDCHI2_ORIVX->GetEntry(i);
    if(D_FDCHI2_ORIVX < 0) return false; 
    
     /*
    b_K_D_ProbNNk->GetEntry(i);
    if(K_D_ProbNNk < 0.2) return false;

    b_pi1_D_ProbNNpi->GetEntry(i);
    if(pi1_D_ProbNNpi < 0.2) return false;

    b_pi2_D_ProbNNpi->GetEntry(i);
    if(pi2_D_ProbNNpi < 0.2) return false;
    
    b_pi_ProbNNpi->GetEntry(i);
    if(pi_ProbNNpi < 0.2) return false;
    
    b_pip_Ks_ProbNNpi->GetEntry(i);
    if(pip_Ks_ProbNNpi < 0.2) return false;
    
    b_pim_Ks_ProbNNpi->GetEntry(i);
    if(pim_Ks_ProbNNpi < 0.2) return false;
    */
    
    return true;
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
    
    if(!_data){
        fChain->SetBranchStatus("*TRUE*",1) ;
        fChain->SetBranchStatus("*BKG*",1) ;
    }
    
    TFile* output = new TFile(_outFileName,"RECREATE");
    TTree* summary_tree = fChain->CloneTree(0);
    
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
        
        //if(!TriggerCuts(j)) continue;
        else if(!LooseCuts(j)) continue;
        
        fChain->GetEntry(i);   
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
