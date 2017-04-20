#! /usr/bin/python

import os,sys

OutputFlag = "_8TeV_Incl_V1/"
CME        = "8000"

FullList = []
FullList.append(["MC15.100014.Herwig7_Matchbox_Internal_ttbar_H7_UE_NNPDF30.py", "20"])

ExtraJO = "MC15.100014.Herwig7_Matchbox_Internal_ttbar_H7_UE_NNPDF30.py"

for entry in FullList:
    ecmEnergy = CME
    jO        = entry[0]
    numJobs   = entry[1]
           
    DSID         = jO.replace("MC15.", "")
    DSID         = DSID[0:6]

    outputFolder = jO.replace("MC15.", "user.lscyboz.")
    outputFolder = outputFolder.replace(".py", "_"+OutputFlag)
          

    cmd = "pathena --site=INFN-LECCE --trf=\"Generate_tf.py --ecmEnergy="+ecmEnergy+" --runNumber="+DSID+" --firstEvent=1 --maxEvents=5000 --randomSeed=%RNDM:100 --jobConfig="+jO+" --outputEVNTFile=test.pool.root \"  --long --split "+numJobs+" --outDS "+outputFolder+" --extOutFile test.pool.root --extFile="+ExtraJO+","+jO

    #print cmd

    os.system(cmd)

