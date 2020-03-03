#j = Job(name = "17")
#try:
#    myApp = prepareGaudiExec("DaVinci","v42r7p2", myPath=".")
#except:
#    myApp = GaudiExec()
#    myApp.directory = "./DaVinciDev_v42r7p2"

#j.application = myApp
#j.application.options = ["MakeTuple_b2dkspi_17.py"]
#j.backend=Dirac()
#j.application.platform = "x86_64-slc6-gcc49-opt"
#datatmp=BKQuery("/LHCb/Collision17/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco17/Stripping29r2/90000000/BHADRON.MDST", dqflag=["OK"]).getDataset()
#datatmp2=BKQuery("/LHCb/Collision17/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco17/Stripping29r2/90000000/BHADRON.MDST", dqflag=["OK"]).getDataset()

#for f in datatmp2.files:
#    datatmp.append(f)

#j.inputdata = datatmp
#j.splitter = SplitByFiles( filesPerJob = 120)
#j.splitter.ignoremissing= True
#j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
#j.do_auto_resubmit = False
#j.parallel_submit = True
#j.submit()

j = Job(name = "17psi")
try:
    myApp = prepareGaudiExec("DaVinci","v42r7p2", myPath=".")
except:
    myApp = GaudiExec()
    myApp.directory = "./DaVinciDev_v42r7p2"

j.application = myApp
j.application.options = ["MakeTuple_b2psikpipi_17.py"]
j.backend=Dirac()
j.application.platform = "x86_64-slc6-gcc49-opt"
datatmp=BKQuery('/LHCb/Collision17/Beam6500GeV-VeloClosed-MagDown/Real Data/Reco17/Stripping29r2/90000000/LEPTONIC.MDST', dqflag=['OK']).getDataset()
#datatmp2=BKQuery('/LHCb/Collision17/Beam6500GeV-VeloClosed-MagUp/Real Data/Reco17/Stripping29r2/90000000/LEPTONIC.MDST', dqflag=['OK']).getDataset()

#for f in datatmp2.files:
#    datatmp.append(f)

j.inputdata = datatmp
j.splitter = SplitByFiles( filesPerJob = 100)
j.splitter.ignoremissing= True
j.outputfiles = [DiracFile('*.root'), LocalFile('stdout')]
j.do_auto_resubmit = False
j.parallel_submit = True
j.submit()
