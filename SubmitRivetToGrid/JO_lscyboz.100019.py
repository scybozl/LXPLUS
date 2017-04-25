theApp.EvtMax = -1 
import AthenaPoolCnvSvc.ReadAthenaPool 
svcMgr.EventSelector.InputCollections =['bla'] 
OutputYoda  = "Output_lscyboz.100019.yoda" 
from AthenaCommon.AlgSequence import AlgSequence 
job = AlgSequence() 
from Rivet_i.Rivet_iConf import Rivet_i 
rivet = Rivet_i() 
rivet.AnalysisPath = os.environ['PWD'] 
rivet.Analyses += [ 'ATLAS_2015_I1404878_custom' ]  
rivet.CrossSection = 253.0 
rivet.HistoFile = OutputYoda 
job += rivet 
from GaudiSvc.GaudiSvcConf import THistSvc 
svcMgr += THistSvc() 
