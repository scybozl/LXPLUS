#! /usr/bin/python

import os,sys
import subprocess

ListRoutines = "RivetATLAS_2015_I1404878_custom.so,ATLAS_2015_I1404878_custom.yoda"

OutputFlag   = "V5"

ContList  = []
ContList.append(['user.lscyboz:user.lscyboz.100019.Herwig7_Matchbox_Internal_ttbar_H7_UE_NNPDF30__8TeV_Incl_V1_EXT0.133357661',   253.0])

for Cont in ContList:
    Container = Cont[0]
    Xsec      = Cont[1]

    Sample     = Container.split(".")
    NewCont    = Sample[2]+"."+Sample[3].replace("_EXT0/", "")+"I1404878"

    if "mc15" in Container:
        NewCont    = Sample[1]+"."+Sample[2].replace("_EXT0/", "")+""

    OutputFile = "Output_"+NewCont+".yoda"
     
     # now make JO script (python)
    submitFileNamePY = "JO_"+NewCont+".py"

    submitfile = open(submitFileNamePY, "w")
    codeLines = []
    codeLines.append("theApp.EvtMax = -1")
    codeLines.append("import AthenaPoolCnvSvc.ReadAthenaPool")
    codeLines.append("svcMgr.EventSelector.InputCollections =['bla']")
    codeLines.append("OutputYoda  = \""+OutputFile+"\"")
    codeLines.append("from AthenaCommon.AlgSequence import AlgSequence")
    codeLines.append("job = AlgSequence()")
    codeLines.append("from Rivet_i.Rivet_iConf import Rivet_i")
    codeLines.append("rivet = Rivet_i()")
    codeLines.append("rivet.AnalysisPath = os.environ['PWD']")
    
    codeLines.append("rivet.Analyses += [ 'ATLAS_2015_I1404878_custom' ] ")
    codeLines.append("rivet.CrossSection = "+str(Xsec))

    codeLines.append("rivet.HistoFile = OutputYoda")
    codeLines.append("job += rivet")
    codeLines.append("from GaudiSvc.GaudiSvcConf import THistSvc")
    codeLines.append("svcMgr += THistSvc()")
    

    for codeLine in codeLines:
        submitfile.write(codeLine+" \n")
        
    submitfile.close()

    os.system("pathena -c 'xs = 1.0 ' "+submitFileNamePY+" --extFile="+ListRoutines+" --nJobs=236 --nFilesPerJob=1 --long  --extOutFile="+OutputFile+" --inDS="+Container+" --outDS=user.lscyboz."+NewCont+"_"+OutputFlag+"/")
