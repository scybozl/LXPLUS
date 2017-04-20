
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/
asetup 19.2.5.17.1,MCProd


Generate_tf.py --ecmEnergy=13000 --maxEvents=50 --runNumber=999999 --firstEvent=0 --randomSeed=123456 --outputEVNTFile=999999_1.pool.root --jobConfig=MC15.999999.Herwig7_Matchbox_Internal_ttbar_H7_UE_NNPDF30.py
