#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>

using namespace std;

void mergeTrees(string old_fileName , string tagging_fileName_Run1, string tagging_fileName_Run2, string tagging_fileNameSS_Run1, string tagging_fileNameSS_Run2 ){

	cout << "Adding OS/SS combo branch to file" << old_fileName << ".root" << endl;

	string outputName = old_fileName;
	outputName.append("_tagged.root"); 

	string oldFileName = old_fileName.append(".root");
	TFile* old_file= new TFile(oldFileName.c_str());
	TTree* old_tree = (TTree*) old_file->Get("DecayTree");
	old_tree->SetBranchStatus("OS_Combination_*",0);
        old_tree->SetBranchStatus("SS_Combination_*",0);
	
	string taggingFileName_Run1 = tagging_fileName_Run1.append(".root");
	TFile* tagging_file_Run1= new TFile(taggingFileName_Run1.c_str());
	TTree* tagging_tree_Run1 = (TTree*) tagging_file_Run1->Get("TaggingTree");

	string taggingFileName_Run2 = tagging_fileName_Run2.append(".root");
	TFile* tagging_file_Run2= new TFile(taggingFileName_Run2.c_str());
	TTree* tagging_tree_Run2 = (TTree*) tagging_file_Run2->Get("TaggingTree");

        string taggingFileNameSS_Run1 = tagging_fileNameSS_Run1.append(".root");
        TFile* tagging_fileSS_Run1= new TFile(taggingFileNameSS_Run1.c_str());
        TTree* tagging_treeSS_Run1 = (TTree*) tagging_fileSS_Run1->Get("TaggingTree");

        string taggingFileNameSS_Run2 = tagging_fileNameSS_Run2.append(".root");
        TFile* tagging_fileSS_Run2= new TFile(taggingFileNameSS_Run2.c_str());
        TTree* tagging_treeSS_Run2 = (TTree*) tagging_fileSS_Run2->Get("TaggingTree");

	Int_t OS_Combination_DEC;
	Double_t OS_Combination_ETA;
	Int_t SS_Combination_DEC;
	Double_t SS_Combination_ETA;
	Int_t OS_Muon_DEC;
	Double_t OS_Muon_ETA;
	Int_t SS_Proton_DEC;
	Double_t SS_Proton_ETA;

	tagging_tree_Run1->SetBranchAddress("OS_Combination_DEC" , &OS_Combination_DEC );
	tagging_tree_Run1->SetBranchAddress("OS_Combination_ETA" , &OS_Combination_ETA );
	tagging_tree_Run1->SetBranchAddress("OS_Muon_DEC" , &OS_Muon_DEC );
	tagging_tree_Run1->SetBranchAddress("OS_Muon_ETA" , &OS_Muon_ETA );
       
	tagging_tree_Run2->SetBranchAddress("OS_Combination_DEC" , &OS_Combination_DEC );
	tagging_tree_Run2->SetBranchAddress("OS_Combination_ETA" , &OS_Combination_ETA );
	tagging_tree_Run2->SetBranchAddress("OS_Muon_DEC" , &OS_Muon_DEC );
	tagging_tree_Run2->SetBranchAddress("OS_Muon_ETA" , &OS_Muon_ETA );

        tagging_treeSS_Run1->SetBranchAddress("SS_Combination_DEC" , &SS_Combination_DEC );
        tagging_treeSS_Run1->SetBranchAddress("SS_Combination_ETA" , &SS_Combination_ETA );
        tagging_treeSS_Run1->SetBranchAddress("SS_Proton_DEC" , &SS_Proton_DEC );
        tagging_treeSS_Run1->SetBranchAddress("SS_Proton_ETA" , &SS_Proton_ETA );

        tagging_treeSS_Run2->SetBranchAddress("SS_Combination_DEC" , &SS_Combination_DEC );
        tagging_treeSS_Run2->SetBranchAddress("SS_Combination_ETA" , &SS_Combination_ETA );
        tagging_treeSS_Run2->SetBranchAddress("SS_Proton_DEC" , &SS_Proton_DEC );
        tagging_treeSS_Run2->SetBranchAddress("SS_Proton_ETA" , &SS_Proton_ETA );

	TFile* output=new TFile(outputName.c_str(),"RECREATE");
	TTree* new_tree_Run1 = old_tree->CopyTree("run == 1");
	TTree* new_tree_Run2 = old_tree->CopyTree("run == 2");
	
	Int_t Bs_OS_Muon_DEC;
	Double_t Bs_OS_Muon_ETA;
        Int_t Bs_SS_Proton_DEC;
        Double_t Bs_SS_Proton_ETA;

	TBranch* b_OS_Muon_DEC_Run1; 
	TBranch* b_OS_Muon_PROB_Run1;
	TBranch* b_OS_Muon_DEC_Run2;
	TBranch* b_OS_Muon_PROB_Run2;

        TBranch* b_SS_Proton_DEC_Run1;
        TBranch* b_SS_Proton_PROB_Run1;
        TBranch* b_SS_Proton_DEC_Run2;
        TBranch* b_SS_Proton_PROB_Run2;

	new_tree_Run1->SetBranchAddress("B_OS_Muon_TAGDEC" , &Bs_OS_Muon_DEC, &b_OS_Muon_DEC_Run1 );
	new_tree_Run1->SetBranchAddress("B_OS_Muon_TAGETA" , &Bs_OS_Muon_ETA, &b_OS_Muon_PROB_Run1);
	new_tree_Run2->SetBranchAddress("B_OS_Muon_TAGDEC" , &Bs_OS_Muon_DEC, &b_OS_Muon_DEC_Run2);
	new_tree_Run2->SetBranchAddress("B_OS_Muon_TAGETA" , &Bs_OS_Muon_ETA, &b_OS_Muon_PROB_Run2);

        new_tree_Run1->SetBranchAddress("B_SS_Proton_TAGDEC" , &Bs_SS_Proton_DEC, &b_SS_Proton_DEC_Run1 );
        new_tree_Run1->SetBranchAddress("B_SS_Proton_TAGETA" , &Bs_SS_Proton_ETA, &b_SS_Proton_PROB_Run1);
        new_tree_Run2->SetBranchAddress("B_SS_Proton_TAGDEC" , &Bs_SS_Proton_DEC, &b_SS_Proton_DEC_Run2);
        new_tree_Run2->SetBranchAddress("B_SS_Proton_TAGETA" , &Bs_SS_Proton_ETA, &b_SS_Proton_PROB_Run2);

	Int_t Bs_OS_Combination_DEC;
	Double_t Bs_OS_Combination_PROB;
	Int_t Bs_SS_Combination_DEC;
	Double_t Bs_SS_Combination_PROB;

	TBranch* OS_Comb_DEC_Branch_Run1 = new_tree_Run1->Branch("OS_Combination_DEC", &Bs_OS_Combination_DEC, "OS_Combination_DEC/I");
	TBranch* OS_Comb_PROB_Branch_Run1 = new_tree_Run1->Branch("OS_Combination_PROB", &Bs_OS_Combination_PROB, "OS_Combination_PROB/D");
	TBranch* OS_Comb_DEC_Branch_Run2 = new_tree_Run2->Branch("OS_Combination_DEC", &Bs_OS_Combination_DEC, "OS_Combination_DEC/I");
	TBranch* OS_Comb_PROB_Branch_Run2 = new_tree_Run2->Branch("OS_Combination_PROB", &Bs_OS_Combination_PROB, "OS_Combination_PROB/D");

        TBranch* SS_Comb_DEC_Branch_Run1 = new_tree_Run1->Branch("SS_Combination_DEC", &Bs_SS_Combination_DEC, "SS_Combination_DEC/I");
        TBranch* SS_Comb_PROB_Branch_Run1 = new_tree_Run1->Branch("SS_Combination_PROB", &Bs_SS_Combination_PROB, "SS_Combination_PROB/D");
        TBranch* SS_Comb_DEC_Branch_Run2 = new_tree_Run2->Branch("SS_Combination_DEC", &Bs_SS_Combination_DEC, "SS_Combination_DEC/I");
        TBranch* SS_Comb_PROB_Branch_Run2 = new_tree_Run2->Branch("SS_Combination_PROB", &Bs_SS_Combination_PROB, "SS_Combination_PROB/D");

	int numEvents_Run1 = tagging_tree_Run1->GetEntries();
	if(new_tree_Run1->GetEntries() != numEvents_Run1){
	      cout << "ERROR::Number of events in Run1 data tree and tagging tree are different" << endl <<  new_tree_Run1->GetEntries() << endl << numEvents_Run1 << endl;
		throw "ERROR";
	}

        if(tagging_treeSS_Run1->GetEntries() != numEvents_Run1){
	  cout << "ERROR::Number of events in Run1 tagging trees are different" << endl <<  tagging_treeSS_Run1->GetEntries() << endl << numEvents_Run1 << endl;
	  throw "ERROR";
        }

	for(int i=0; i< numEvents_Run1; i++){
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_Run1 << endl;
		tagging_tree_Run1->GetEntry(i);
		tagging_treeSS_Run1->GetEntry(i);
	
		b_OS_Muon_DEC_Run1->GetEntry(i);
		b_OS_Muon_PROB_Run1->GetEntry(i);
		b_SS_Proton_DEC_Run1->GetEntry(i);
                b_SS_Proton_PROB_Run1->GetEntry(i);
		if(Bs_OS_Muon_ETA != OS_Muon_ETA){
			cout << "ERROR::Inconsistent events" << endl;
			throw "ERROR";
		}
		if(Bs_SS_Proton_ETA != SS_Proton_ETA){
		  cout << "ERROR::Inconsistent events" << endl;
		  throw "ERROR";
                }

		Bs_OS_Combination_DEC = OS_Combination_DEC;
		Bs_OS_Combination_PROB = OS_Combination_ETA;
		Bs_SS_Combination_DEC = SS_Combination_DEC;
		Bs_SS_Combination_PROB = SS_Combination_ETA;
	
		OS_Comb_DEC_Branch_Run1->Fill();
		OS_Comb_PROB_Branch_Run1->Fill();
		SS_Comb_DEC_Branch_Run1->Fill();
		SS_Comb_PROB_Branch_Run1->Fill();
	}
	
	int numEvents_Run2 = tagging_tree_Run2->GetEntries();
	if(new_tree_Run2->GetEntries() != numEvents_Run2){
		 cout << "ERROR: Number of events in Run2 data tree and tagging tree are different" << endl <<  new_tree_Run2->GetEntries() << endl << numEvents_Run2 << endl;
		throw "ERROR";
	}
        if(tagging_treeSS_Run2->GetEntries() != numEvents_Run2){
	  cout << "ERROR: Number of events in Run2 tagging trees are different" << endl <<  tagging_treeSS_Run2->GetEntries() << endl << numEvents_Run2 << endl;
	  throw "ERROR";
        }

	for(int i=0; i< numEvents_Run2; i++){
		if (0ul == (i % 10000ul)) cout << "Read event " << i << "/" << numEvents_Run2 << endl;
		tagging_tree_Run2->GetEntry(i);
                tagging_treeSS_Run2->GetEntry(i);
	
		b_OS_Muon_DEC_Run2->GetEntry(i);
		b_OS_Muon_PROB_Run2->GetEntry(i);
		b_SS_Proton_DEC_Run2->GetEntry(i);
                b_SS_Proton_PROB_Run2->GetEntry(i);
		if(Bs_OS_Muon_ETA != OS_Muon_ETA){
			cout << "ERROR::Inconsistent events" << endl;
			cout << Bs_OS_Muon_DEC << endl;
			cout << OS_Muon_DEC  << endl;
			cout << Bs_OS_Muon_ETA << endl;
                        cout << OS_Muon_ETA << endl;
			throw "ERROR";
		}
                if(Bs_SS_Proton_ETA != SS_Proton_ETA){
		  cout << "ERROR::Inconsistent events" << endl;
		  cout << Bs_SS_Proton_DEC  << endl;
                  cout << SS_Proton_DEC  << endl;
                  cout << Bs_SS_Proton_ETA  << endl;
                  cout <<  SS_Proton_ETA  << endl;
		  throw "ERROR";
                }

		Bs_OS_Combination_DEC = OS_Combination_DEC;
		Bs_OS_Combination_PROB = OS_Combination_ETA;
		Bs_SS_Combination_DEC = SS_Combination_DEC;
		Bs_SS_Combination_PROB = SS_Combination_ETA;
	
		OS_Comb_DEC_Branch_Run2->Fill();
		OS_Comb_PROB_Branch_Run2->Fill();
		SS_Comb_DEC_Branch_Run2->Fill();
		SS_Comb_PROB_Branch_Run2->Fill();
	}

  	TList *list = new TList; 
  	list->Add(new_tree_Run1); 
  	list->Add(new_tree_Run2); 
  
  	TTree *new_tree = TTree::MergeTrees(list); 
  	new_tree->Write();

	output->Close();	
	old_file->Close();
	tagging_file_Run1->Close();
	tagging_file_Run2->Close();
        tagging_fileSS_Run1->Close();
        tagging_fileSS_Run2->Close();
	cout << "Wrote output file " << outputName << endl;
}

int main(int argc, char** argv){
        mergeTrees("b2dkspi_sw","OS_combo_Run1","OS_combo_Run2","SS_combo_Run1","SS_combo_Run2");
	return 0;
}
