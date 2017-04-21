from Herwig7_i.Herwig7_iConf import Herwig7
from Herwig7_i.Herwig7ConfigMatchbox import Hw7ConfigMatchbox

genSeq += Herwig7()

## Provide config information
evgenConfig.generators += ["Herwig7"]
evgenConfig.description = "Herwig7 ttbar sample with NNPDF30 ME PDF"
evgenConfig.keywords    = ["SM","ttbar"]
evgenConfig.contact     = ["Ludovic Scyboz"]


## initialize generator configuration object
generator = Hw7ConfigMatchbox(genSeq, runArgs, run_name="HerwigMatchbox", beams="pp")

## configure generator
generator.me_pdf_commands(order="NLO", name="NNPDF30_nlo_as_0118")
generator.tune_commands()


generator.add_commands("""
##################################################
## Process selection
##################################################

## Model assumptions
read Matchbox/StandardModelLike.in

## Set the hard process
set /Herwig/MatrixElements/Matchbox/Factory:OrderInAlphaS 2
set /Herwig/MatrixElements/Matchbox/Factory:OrderInAlphaEW 0
do /Herwig/MatrixElements/Matchbox/Factory:Process p p -> t tbar

read Matchbox/OnShellTopProduction.in


##################################################
## Matrix element library selection
##################################################

# read Matchbox/MadGraph-GoSam.in
# read Matchbox/MadGraph-MadGraph.in
# read Matchbox/MadGraph-NJet.in
read Matchbox/MadGraph-OpenLoops.in
# read Matchbox/HJets.in
# read Matchbox/VBFNLO.in

#cd /Herwig/MatrixElements/Matchbox
#set Amplitudes/GoSam:SetupInFilename gosamtTNLO172.5.rc

cd /Herwig/MatrixElements/Matchbox
insert Factory:DiagramGenerator:ExcludeInternal 0 /Herwig/Particles/e-
insert Factory:DiagramGenerator:ExcludeInternal 0 /Herwig/Particles/nu_ebar
insert Factory:DiagramGenerator:ExcludeInternal 0 /Herwig/Particles/mu+
insert Factory:DiagramGenerator:ExcludeInternal 0 /Herwig/Particles/nu_mu
insert Factory:DiagramGenerator:ExcludeInternal 0 /Herwig/Particles/h0


##################################################
## Cut selection
## See the documentation for more options
##################################################

## cuts on additional jets

# read Matchbox/DefaultPPJets.in

# insert JetCuts:JetRegions 0 FirstJet
# insert JetCuts:JetRegions 1 SecondJet
# insert JetCuts:JetRegions 2 ThirdJet
# insert JetCuts:JetRegions 3 FourthJet

##################################################
## Scale choice
## See the documentation for more options
##################################################

cd /Herwig/MatrixElements/Matchbox
set Factory:ScaleChoice Scales/TopPairMTScale

##################################################
## Matching and shower selection
## Please also see flavour scheme settings
## towards the end of the input file.
##################################################

# read Matchbox/MCatNLO-DefaultShower.in
read Matchbox/Powheg-DefaultShower.in
## use for strict LO/NLO comparisons
# read Matchbox/MCatLO-DefaultShower.in
## use for improved LO showering
# read Matchbox/LO-DefaultShower.in

set /Herwig/Shower/GtoQQbarSplitFn:AngularOrdered Yes
set /Herwig/Shower/Evolver:MECorrMode 1
set /Herwig/Shower/PartnerFinder:PartnerMethod Random
set /Herwig/Shower/PartnerFinder:ScaleChoice Partner
set /Herwig/Shower/ShowerHandler:RestrictPhasespace On
set /Herwig/Shower/ShowerHandler:MaxPtIsMuF Yes
set /Herwig/Shower/GammatoQQbarSudakov:Alpha /Herwig/Shower/AlphaQED

# read Matchbox/MCatNLO-DipoleShower.in
# read Matchbox/Powheg-DipoleShower.in
## use for strict LO/NLO comparisons
# read Matchbox/MCatLO-DipoleShower.in
## use for improved LO showering
# read Matchbox/LO-DipoleShower.in

# read Matchbox/LO-NoShower.in
# read Matchbox/NLO-NoShower.in

cd /Herwig/Particles
set t:NominalMass 172.5*GeV
set t:HardProcessMass 172.5*GeV
set t:Width 1.3167*GeV

set W+:NominalMass 80.399*GeV
set W+:HardProcessMass 80.399*GeV
set W+:Width 2.09974*GeV

set Z0:NominalMass 91.1876*GeV
set Z0:HardProcessMass 91.1876*GeV
set Z0:Width 2.50966*GeV

set /Herwig/Model:EW/Scheme GMuScheme
set /Herwig/Model:EW/FermiConstant 1.16637e-05
set /Herwig/Model:EW/RecalculateEW On
set /Herwig/MatrixElements/Matchbox/Factory:FixedQEDCouplings On


#cd /Herwig/Generators
#insert EventGenerator:AnalysisHandlers 0 /Herwig/Analysis/HepMCFile
#set /Herwig/Analysis/HepMCFile:PrintEvent 1000000
#set /Herwig/Analysis/HepMCFile:Format GenEvent
#set /Herwig/Analysis/HepMCFile:Units GeV_mm
#set /Herwig/Analysis/HepMCFile:Filename 14.02.17.hepmc

cd /Herwig/Analysis
set Basics:CheckQuark No

##################################################
## PDF choice
##################################################

read Matchbox/FiveFlavourScheme.in
## required for dipole shower and fixed order in five flavour scheme
# read Matchbox/FiveFlavourNoBMassScheme.in
""")

## CAUTION: Extremely crude sampling for testing purposes
# generator.sampler_commands("CellGridSampler", 1000, 2, 1000, 1, 100)
generator.sampler_commands("MonacoSampler", 10000, 2, 1000, 1, 100)
# generator.sampler_commands("FlatBinSampler", 1000, 2, 1000, 1, 100)

## run generator
generator.run(cleanup_herwig_scratch=False)

