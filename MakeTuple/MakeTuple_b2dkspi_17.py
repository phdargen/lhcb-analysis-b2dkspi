from Gaudi.Configuration import *
from Configurables import DaVinci
from Configurables import CombineParticles

############# Global settings
year = "2017"
data = True
down = False
stream = "AllStreams"
if (data):
    stream = "Bhadron"

# Event filter
from PhysSelPython.Wrappers import AutomaticData, Selection, SelectionSequence
from Configurables import FilterDesktop
line = 'B02DKsPiDDD2HHHCFPIDBeauty2CharmLine'
inputs = 'Phys/{0}/Particles'.format(line)
reqsel = AutomaticData(Location = inputs)
Bs2DsXSel = FilterDesktop("Bs2DsXSel", Code = "(INTREE((ABSID=='B0')&(M>5))) & (INTREE((ABSID=='B0')&(M<60000)))")
                          #(INTREE((ABSID=='D+')&(M>1888))) & (INTREE((ABSID=='D+')&(M<2048))) & (INTREE((ABSID=='B0')&(M>0))) & (INTREE((ABSID=='B0')&(M<60000))) & (INTREE((ABSID=='K_1(1270)+')&(M<3000)))
#                          ")
MyFilterSel = Selection("MyFilterSel", Algorithm = Bs2DsXSel, RequiredSelections = [reqsel])

from Configurables import CheckPV
checkPVs = CheckPV("checkPVs")
checkPVs.MinPVs = 1
checkPVs.MaxPVs = -1


triggerlines_Run1 = [
		"L0HadronDecision", 
                "L0MuonDecision",
		"L0DiMuonDecision",
		"L0ElectronDecision",
		"L0PhotonDecision",
                "Hlt1TrackAllL0Decision", 
                'Hlt2IncPhiDecision', 
                'Hlt2Topo2BodyBBDTDecision', 
                'Hlt2Topo3BodyBBDTDecision', 
                'Hlt2Topo4BodyBBDTDecision' 
]

triggerlines_Run2 = [
		"L0HadronDecision", 
                "L0MuonDecision",
		"L0DiMuonDecision",
		"L0ElectronDecision",
		"L0PhotonDecision",
                "Hlt1TrackMVADecision",
                "Hlt1TwoTrackMVADecision",
		'Hlt2IncPhiDecision', 
		"Hlt2PhiIncPhiDecision",
                'Hlt2Topo2BodyDecision',
                'Hlt2Topo3BodyDecision',
                'Hlt2Topo4BodyDecision'
]

triggerlines = triggerlines_Run2


###############   Pre Filter, does not really do much except choose only candidates passing the Stripping line, maybe beneficial to performance
from Configurables import LoKi__HDRFilter as StripFilter
stripFilter = StripFilter( 'stripPassFilter',\
                           Code = "HLT_PASS('StrippingB02DKsPiDDD2HHHCFPIDBeauty2CharmLineDecision')",\
                           Location= "/Event/Strip/Phys/DecReports")

############# DecayTreeTuple
from DecayTreeTuple.Configuration import *
from Configurables import TupleToolTISTOS
from Configurables import PrintDecayTree, PrintDecayTreeTool
##subpid stuff
#from Configurables import SubPIDMMFilter
from Configurables import SubstitutePID
from Configurables import TupleToolDecayTreeFitter, TupleToolTrackIsolation, TupleToolTagging, TupleToolRecoStats, TupleToolKinematic, TupleToolGeometry, TupleToolVtxIsoln
from Configurables import LoKi__Hybrid__TupleTool

#B0 -> (D- -> K K pi) (K_1(1270)+ -> K+ pi+ pi-)
b2dkpipiTuple = DecayTreeTuple("B2DKspi_Tuple")
b2dkpipiTuple.Decay = "[[B0]cc -> ^(D- -> ^K+ ^pi- ^pi-) ^(KS0 -> ^pi+ ^pi-) ^pi+]CC"
b2dkpipiTuple.Branches= {
"B" : "^([[B0]cc -> (D- -> K+ pi- pi-) (KS0 -> pi+ pi-) pi+]CC)" ,
"Ks" : "[[B0]cc -> (D- -> K+ pi- pi-) ^(KS0 -> pi+ pi-) pi+]CC",
"D" : "[[B0]cc -> ^(D- -> K+ pi- pi-) (KS0 -> pi+ pi-) pi+]CC",
"pip_Ks" : "[[B0]cc -> (D- -> K+ pi- pi-) (KS0 -> ^pi+ pi-) pi+]CC",
"pim_Ks" : "[[B0]cc -> (D- -> K+ pi- pi-) (KS0 -> pi+ ^pi-) pi+]CC",
"K_D" : "[[B0]cc -> (D- -> ^K+ pi- pi-) (KS0 -> pi+ pi-) pi+]CC",
"pi1_D" : "[[B0]cc -> (D- -> K+ ^pi- pi-) (KS0 -> pi+ pi-) pi+]CC",
"pi2_D" : "[[B0]cc -> (D- -> K+ pi- ^pi-) (KS0 -> pi+ pi-) pi+]CC",
"pi" : "[[B0]cc -> (D- -> K+ pi- pi-) (KS0 -> pi+ pi-) ^pi+]CC",
}
b2dkpipiTuple.ReFitPVs = True

#config tools
b2dkpipiTuple.ToolList +=  ["TupleToolGeometry", \
                            "TupleToolKinematic", \
                            "TupleToolPrimaries", \
                           # "TupleToolEventInfo", \
                            "TupleToolTrackInfo", \
                            "TupleToolRecoStats", \
                            #"TupleToolAngles", \
                            "TupleToolPid", \
                            "TupleToolTrackIsolation",
                           # "TupleToolVtxIsoln",
                            "TupleToolTagging" 
			    ]

if (data==False):
    b2dkpipiTuple.ToolList +=  [
                            "TupleToolMCTruth", \
                            "TupleToolMCBackgroundInfo",
			     "TupleToolPhotonInfo"
				]
    MCTruth = TupleToolMCTruth()
    MCTruth.ToolList =  [
         "MCTupleToolHierarchy"
        , "MCTupleToolKinematic"
        , "MCTupleToolReconstructed"
        ]
    b2dkpipiTuple.addTool(MCTruth) 

b2dkpipiTuple.addTool(TupleToolTrackIsolation, name="TupleToolTrackIsolation")
b2dkpipiTuple.TupleToolTrackIsolation.FillAsymmetry = True
b2dkpipiTuple.TupleToolTrackIsolation.FillDeltaAngles = False
b2dkpipiTuple.TupleToolTrackIsolation.MinConeAngle = 1.0

b2dkpipiTuple.addTool( TupleToolKinematic,name = "TupleToolKinematic" )
b2dkpipiTuple.TupleToolKinematic.Verbose = False

b2dkpipiTuple.addTool( TupleToolGeometry, name = "TupleToolGeometry" )
b2dkpipiTuple.TupleToolGeometry.Verbose = True
b2dkpipiTuple.TupleToolGeometry.RefitPVs = True
b2dkpipiTuple.TupleToolGeometry.FillMultiPV = True

b2dkpipiTuple.addTool(TupleToolRecoStats, name="TupleToolRecoStats")
b2dkpipiTuple.TupleToolRecoStats.Verbose = True
b2dkpipiTuple.UseLabXSyntax = True                          
b2dkpipiTuple.RevertToPositiveID = False

b2dkpipiTuple.addTool(TupleToolDecay, name="B")
b2dkpipiTuple.addTool(TupleToolDecay, name="Ks")
b2dkpipiTuple.addTool(TupleToolDecay, name="D")

#b2dkpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("DTF"))
#b2dkpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/DTF" ]        
#b2dkpipiTuple.Bs.DTF.constrainToOriginVertex = True
#b2dkpipiTuple.Bs.DTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+'  }
#b2dkpipiTuple.Bs.DTF.daughtersToConstrain = ["D_s-", "D_s+", ]  
#b2dkpipiTuple.Bs.DTF.UpdateDaughters = True
#b2dkpipiTuple.Bs.DTF.Verbose = True

#b2dkpipiTuple.Bs.addTool(TupleToolDecayTreeFitter("B0DTF"))
#b2dkpipiTuple.Bs.ToolList +=  ["TupleToolDecayTreeFitter/B0DTF" ]
#b2dkpipiTuple.Bs.B0DTF.constrainToOriginVertex = True
#b2dkpipiTuple.Bs.B0DTF.Substitutions = { 'Beauty -> X+ ^(Charm & X-)' : 'D_s-', 'Beauty -> X- ^(Charm & X+)' : 'D_s+' }
#b2dkpipiTuple.Bs.B0DTF.daughtersToConstrain = ["D_s-", "B0", "D_s+", "B~0" ]  
#b2dkpipiTuple.Bs.B0DTF.UpdateDaughters = True
#b2dkpipiTuple.Bs.B0DTF.Verbose = True

b2dkpipiTuple.B.addTool(TupleToolDecayTreeFitter("DTF"))
b2dkpipiTuple.B.ToolList +=  ["TupleToolDecayTreeFitter/DTF" ]
b2dkpipiTuple.B.DTF.constrainToOriginVertex = True
b2dkpipiTuple.B.DTF.daughtersToConstrain = ["D-", "KS0"]
#b2dkpipiTuple.B.DTF.UpdateDaughters = True
b2dkpipiTuple.B.DTF.Verbose = True

b2dkpipiTuple.B.addTool(TupleToolDecayTreeFitter("BDTF"))
b2dkpipiTuple.B.ToolList +=  ["TupleToolDecayTreeFitter/BDTF" ]
b2dkpipiTuple.B.BDTF.constrainToOriginVertex = True
b2dkpipiTuple.B.BDTF.daughtersToConstrain = ["B0"]
#b2dkpipiTuple.B.BDTF.UpdateDaughters = True                                                                                                            
b2dkpipiTuple.B.BDTF.Verbose = True

b2dkpipiTuple.B.addTool(TupleToolDecayTreeFitter("FullDTF"))
b2dkpipiTuple.B.ToolList +=  ["TupleToolDecayTreeFitter/FullDTF" ]         
b2dkpipiTuple.B.FullDTF.constrainToOriginVertex = True
b2dkpipiTuple.B.FullDTF.daughtersToConstrain = ["D-", "B0", "KS0"]  
b2dkpipiTuple.B.FullDTF.UpdateDaughters = True
b2dkpipiTuple.B.FullDTF.Verbose = True

b2dkpipiTuple.B.addTool(TupleToolDecayTreeFitter("PV"))
b2dkpipiTuple.B.ToolList +=  ["TupleToolDecayTreeFitter/PV" ]         
b2dkpipiTuple.B.PV.constrainToOriginVertex = True 
b2dkpipiTuple.B.PV.UpdateDaughters = True
b2dkpipiTuple.B.PV.Verbose = True

LoKiTool = b2dkpipiTuple.addTupleTool("LoKi::Hybrid::TupleTool/LoKiTool")
LoKiTool.Variables = { "ETA" : "ETA" };

LoKiToolB = b2dkpipiTuple.B.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToolB")
LoKiToolB.Variables = { 
                       "ETA" : "ETA" , "DOCA1" : "DOCA(1,2)", "TAU" : "BPVLTIME()", "TAUERR" : "BPVLTERR()"
                    	};

LoKiToolKs = b2dkpipiTuple.Ks.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToolKs")
LoKiToolKs.Variables = {
                       "ETA" : "ETA" , "DOCA1" : "DOCA(1,2)", "TAU" : "BPVLTIME()", "TAUERR" : "BPVLTERR()"
                        };

LoKiToolD = b2dkpipiTuple.D.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToolD")
LoKiToolD.Variables = {
                       "ETA" : "ETA" , "DOCA1" : "DOCA(1,2)", "TAU" : "BPVLTIME()", "TAUERR" : "BPVLTERR()"
                        };

#b2dkpipiTuple.addTool(TupleToolDecay, name="Ds")
#LoKiToolDs = b2dkpipiTuple.Ds.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToolDs")
#LoKiToolDs.Variables = { "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)"
 #                      };

#b2dkpipiTuple.addTool(TupleToolDecay, name="K_1_1270_plus")
#LoKiToolK_1_1270_plus = b2dkpipiTuple.K_1_1270_plus.addTupleTool("LoKi::Hybrid::TupleTool/LoKiToolK_1_1270_plus")
#LoKiToolK_1_1270_plus.Variables = { "DOCA1" : "DOCA(1,2)" , "DOCA2" : "DOCA(1,3)" , "DOCA3" : "DOCA(2,3)"
 #                      };

##tagging config
#b2dkpipiTuple.addTool(TupleToolTagging, name="TupleToolTagging")
#b2dkpipiTuple.TupleToolTagging.Verbose = True
#b2dkpipiTuple.TupleToolTagging.StoreTaggersInfo = True
#tag=b2dkpipiTuple.Bs.addTupleTool( TupleToolTagging, name = "BsAll")
#configureTaggingTools(tag, "Bs")

#from Configurables import BTaggingTool
#tt_tagging = b2dkpipiTuple.addTupleTool("TupleToolTagging") 
#tt_tagging.Verbose = True
#tt_tagging.AddMVAFeatureInfo = False
#tt_tagging.AddTagPartsInfo = False 

#btagtool = tt_tagging.addTool(BTaggingTool , name = "MyBTaggingTool")
#from FlavourTagging.Tunings import applyTuning as applyFTTuning # pick the right tuning here ...
#applyFTTuning(btagtool , tuning_version="Summer2017Optimisation_v4_Run2")
#tt_tagging.TaggingToolName = btagtool.getFullName ()

TupleToolTagging = TupleToolTagging('TupleToolTagging', useFTonDST=True)
TupleToolTagging.Verbose = True
b2dkpipiTuple.addTool(TupleToolTagging)



#trigger config
b2dkpipitt = b2dkpipiTuple.B.addTupleTool(TupleToolTISTOS)
b2dkpipitt.TriggerList = triggerlines
b2dkpipitt.FillL0 = True
b2dkpipitt.FillHlt1 = True
b2dkpipitt.FillHlt2 = True
b2dkpipitt.Verbose = True
b2dkpipitt.VerboseL0 = True
b2dkpipitt.VerboseHlt1 = True
b2dkpipitt.VerboseHlt2 = True

#b2dkpipiTuple.Bs.addTool(TupleToolTISTOS,name="TisTosB")
#b2dkpipiTuple.Bs.TisTosB.Verbose=True
#b2dkpipiTuple.Bs.TisTosB.TriggerList= triggers

#printer
#b2dkpipiprinter = PrintDecayTree("PrintB2Dkpipi")
#b2dkpipiprinter.addTool( PrintDecayTreeTool, name = "PrintDecay" )
#b2dkpipiprinter.PrintDecay.Information = "Name M P Px Py Pz Pt chi2"
#b2dkpipiprinter.Inputs = [  "/Event/"+stream+"/Phys/B02DKPiPiD2HHHPIDBeauty2CharmLine/Particles"  ]

#main sequence
makeb2dkpipiseq = SelectionSequence("makeb2dkpipiseq", TopSelection = MyFilterSel)
b2dkpipiTuple.Inputs = [makeb2dkpipiseq.outputLocation()]
b2dkpipiseq = GaudiSequencer("B2dkpipiSeq")
#b2dkpipiseq.RootInTES = '/Event/{0}'.format(stream)
b2dkpipiseq.Members += [makeb2dkpipiseq.sequence(),b2dkpipiTuple]

#b2dkpipiTuple.Inputs = [ "/Event/"+stream+"/Phys/B02DKPiPiD2HHHPIDBeauty2CharmLine/Particles" ]
#b2dkpipiseq = GaudiSequencer("B2dkpipiSeq")
#b2dkpipiseq.Members += [b2dkpipiprinter,b2dkpipiTuple]


#

#
#
#
DaVinci().InputType = 'MDST'
DaVinci().RootInTES = '/Event/{0}'.format(stream) #to be used for micro-DSTs
DaVinci().EventPreFilters = [stripFilter]
DaVinci().UserAlgorithms += [ b2dkpipiseq]

DaVinci().DataType = year
if (data):
    DaVinci().Simulation = False
else:
    DaVinci().Simulation = True
    
DaVinci().EvtMax = -1
#DaVinci().EvtMax = 1000
DaVinci().SkipEvents = 0
DaVinci().PrintFreq = 50000
DaVinci().TupleFile = "b2dkspi.root"
DaVinci().Lumi = True

from Configurables import CondDB, CondDBAccessSvc

if (data):
	CondDB().LatestGlobalTagByDataType = year
    
else:  
    if(year == "2012"):
	    DaVinci().DDDBtag = "dddb-20130929-1"
    	    if (down):
            	DaVinci().CondDBtag = "sim-20141210-1-vc-md100"
    	    else:
        	DaVinci().CondDBtag = "sim-20141210-1-vc-mu100" 

    if(year == "2011"):
	    DaVinci().DDDBtag = "dddb-20130929"
    	    if (down):
            	DaVinci().CondDBtag = "sim-20141210-vc-md100"
    	    else:
        	DaVinci().CondDBtag = "sim-20141210-vc-mu100" 



## Use the local input data
#from GaudiConf import IOHelper
#IOHelper().inputFiles([
#    '00071501_00000026_1.bhadron.mdst'
#], clear=True)
