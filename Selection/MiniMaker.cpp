#include "DecayTree.h"
#include <iostream>
#include <stdlib.h>
#include <cstdlib>

int main(int argc, char** argv){
        
    Decay::Type decay;
    Year::Type year;
    DataType::Type dataType;
    McEventType::Type mcEventType = McEventType::BdDKspi;

    TString polarity = "Both";
    
    if((string)argv[1] == "B2DKspi_LL") decay = Decay::B2DKspi_LL ;
    else if((string)argv[1] == "B2DKspi_DD") decay = Decay::B2DKspi_DD ;
    else if((string)argv[1] == "B2DKsK_LL") decay = Decay::B2DKsK_LL ;
    else if((string)argv[1] == "B2DKsK_DD") decay = Decay::B2DKsK_DD ;
    else {
		  cout << "Missing or wrong config! I'll crash now, see commands.txt how to run me " << endl;
		  throw "ERROR" ;
	  }

    if((string)argv[2] == "Data") dataType = DataType::data ;
    else if((string)argv[2] == "MC") dataType = DataType::mc ;
    else {
		  cout << "Missing or wrong config! I'll crash now, see commands.txt how to run me " << endl;
		  throw "ERROR" ;
	  }

    if(atoi(argv[3]) == 11) year = Year::y11;
    else if(atoi(argv[3]) == 12) year = Year::y12;
    else if(atoi(argv[3]) == 15) year = Year::y15;
    else if(atoi(argv[3]) == 16) year = Year::y16;
    else if(atoi(argv[3]) == 17) year = Year::y17;
    else if(atoi(argv[3]) == 18) year = Year::y18;

    else {
		  cout << "Missing or wrong config! I'll crash now, see commands.txt how to run me " << endl;
		  throw "ERROR" ;
	  }    

    if(argv[4] != 0){
	    //enum  Type { BdDKspi, BsDKsK, BsDstKsK, BdDstKspi, BsDstKspi };
      if((string)argv[4] == "BsDKsK") mcEventType = McEventType::BsDKsK;
      else if((string)argv[4] == "BsDstKsK") mcEventType = McEventType::BsDstKsK;
      else if((string)argv[4] == "BdDstKspi") mcEventType = McEventType::BdDstKspi;
      else if((string)argv[4] == "BsDstKspi") mcEventType = McEventType::BsDstKspi;
    }

    /// bools : bkg, ltu, ss
    DecayTree d(decay, year, dataType, polarity,"./","/eos/lhcb/user/p/phdargen/B2DKspi/Preselected/",mcEventType);
    d.Loop();

    return 0;
}
