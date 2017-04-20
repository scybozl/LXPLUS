#! /bin/bash

export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/
source /cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase/user/atlasLocalSetup.sh
asetup 19.2.5.17.1,MCProd,here

localSetupPandaClient
localSetupRucioClients
voms-proxy-init --voms=atlas
