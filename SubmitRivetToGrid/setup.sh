#! /bin/bash

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
asetup 20.7.7.7,AtlasProduction,here
#asetup 19.2.1.3,slc6,64


localSetupPandaClient
localSetupRucioClients
voms-proxy-init --voms=atlas
