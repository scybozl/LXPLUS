theApp.EvtMax = -1

import os
import AthenaPoolCnvSvc.ReadAthenaPool
svcMgr.EventSelector.InputCollections = ["user.lscyboz.100013.Herwig7_Matchbox_Internal_ttbar_H7_UE_NNPDF30__13TeV_Incl_V1_EXT0.132871236/user.lscyboz.11201462.EXT0._000001.test.pool.root"]

from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()

from Rivet_i.Rivet_iConf import Rivet_i

rivet = Rivet_i()
rivet.AnalysisPath = os.environ['PWD']
rivet.Analyses += [ 'ATLAS_2015_I1404878_custom' ]
rivet.RunName = ""
rivet.HistoFile = "myanalysis"
rivet.CrossSection = 2.5300E+02
job += rivet
