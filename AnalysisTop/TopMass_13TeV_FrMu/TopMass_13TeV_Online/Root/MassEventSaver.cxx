/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#include "TopMass_13TeV_Online/MassEventSaver.h"
#include "TopEvent/Event.h"
#include "TopEventSelectionTools/TreeManager.h"
#include "TopParticleLevel/ParticleLevelEvent.h"
#include "TopConfiguration/TopConfig.h"

#include <TRandom3.h>


namespace top{
  ///-- Construcutor --///
  MassEventSaver::MassEventSaver() :
    m_tma_nbjets(0.),
    m_tma_etdr(0.),
    m_tma_met(0.),
    m_tma_nleps(0.),
    m_tma_nelecs(0.),
    m_tma_nmuons(0.),
    m_tma_njets(0.),
    m_tma_mlb_minavg(0.),
    m_tma_mlb_minavglow(0.),
    m_tma_mlb_minavghigh(0.),
    m_tma_mlb_minmax(0.),
    m_tma_mlb_minmaxlow(0.),
    m_tma_mlb_minmaxavg(0.),
    m_tma_pTlb_1(0.),
    m_tma_pTlb_2(0.),
    m_tma_dRlb_1(0.),
    m_tma_dRlb_2(0.),
    m_tma_mll(0.),
    m_tma_pTll(0.),
    m_tma_dRll(0.),
    m_tma_mbb(0.),
    m_tma_pTbb(0.),
    m_tma_dRbb(0.),
    m_tma_Rbq_avgLJ(0.),
    m_tma_Rbq_leadLJ(0.)
  {
  }
  
  ///-- initialize - done once at the start of a job before the loop over events --///
  void MassEventSaver::initialize(std::shared_ptr<top::TopConfig> config, TFile* file, const std::vector<std::string>& extraBranches)
{
    ///-- Let the base class do all the hard work --///
    ///-- It will setup TTrees for each systematic with a standard set of variables --///
    top::EventSaverFlatNtuple::initialize(config, file, extraBranches);
    
    ///-- Loop over the systematic TTrees and add the custom variables --///
    for (auto systematicTree : treeManagers()) {

      // "tma" stands here for "top mass analysis"
      systematicTree->makeOutputVariable(m_tma_klfitter_mtop_param, "tma_klfitter_mtop_param");
      systematicTree->makeOutputVariable(m_tma_original_mw,         "tma_original_mw");
      systematicTree->makeOutputVariable(m_tma_original_rbq,        "tma_original_rbq");
      systematicTree->makeOutputVariable(m_tma_bjet_pt,             "tma_bjet_pt");
      systematicTree->makeOutputVariable(m_tma_bjet_eta,            "tma_bjet_eta");
      systematicTree->makeOutputVariable(m_tma_bjet_phi,            "tma_bjet_phi");
      systematicTree->makeOutputVariable(m_tma_bjet_e,              "tma_bjet_e");
      systematicTree->makeOutputVariable(m_tma_nbjets,              "tma_nbjets");


      // plots for dilepton channel: reco-level                                                                                                                                                                        
      systematicTree->makeOutputVariable(m_tma_etdr,           "tma_etdr");
      systematicTree->makeOutputVariable(m_tma_met,            "tma_met");
      systematicTree->makeOutputVariable(m_tma_nleps,          "tma_nleps");
      systematicTree->makeOutputVariable(m_tma_lep_pt,         "tma_lep_pt");
      systematicTree->makeOutputVariable(m_tma_lep_eta,        "tma_lep_eta");
      systematicTree->makeOutputVariable(m_tma_lep_phi,        "tma_lep_phi");
      systematicTree->makeOutputVariable(m_tma_lep_type,       "tma_lep_type");
      systematicTree->makeOutputVariable(m_tma_lep_charge,     "tma_lep_charge");
      systematicTree->makeOutputVariable(m_tma_lep_e,          "tma_lep_e");
      systematicTree->makeOutputVariable(m_tma_el_pt,          "tma_el_pt");
      systematicTree->makeOutputVariable(m_tma_el_eta,         "tma_el_eta");
      systematicTree->makeOutputVariable(m_tma_mu_pt,          "tma_mu_pt");
      systematicTree->makeOutputVariable(m_tma_mu_eta,         "tma_mu_eta");
      systematicTree->makeOutputVariable(m_tma_njets,          "tma_njets");
      systematicTree->makeOutputVariable(m_tma_jet_pt,         "tma_jet_pt");
      systematicTree->makeOutputVariable(m_tma_jet_eta,        "tma_jet_eta");
      systematicTree->makeOutputVariable(m_tma_jet_phi,        "tma_jet_phi");
      systematicTree->makeOutputVariable(m_tma_jet_e,          "tma_jet_e");
      systematicTree->makeOutputVariable(m_tma_mlb_minavg,     "tma_mlb_minavg");
      systematicTree->makeOutputVariable(m_tma_mlb_minavglow,  "tma_mlb_minavglow");
      systematicTree->makeOutputVariable(m_tma_mlb_minavghigh, "tma_mlb_minavghigh");
      systematicTree->makeOutputVariable(m_tma_mlb_minmax,     "tma_mlb_minmax");
      systematicTree->makeOutputVariable(m_tma_mlb_minmaxlow,  "tma_mlb_minmaxlow");
      systematicTree->makeOutputVariable(m_tma_mlb_minmaxavg,  "tma_mlb_minmaxavg");
      systematicTree->makeOutputVariable(m_tma_pTlb_1,         "tma_pTlb_1");
      systematicTree->makeOutputVariable(m_tma_pTlb_2,         "tma_pTlb_2");
      systematicTree->makeOutputVariable(m_tma_dRlb_1,         "tma_dRlb_1");
      systematicTree->makeOutputVariable(m_tma_dRlb_2,         "tma_dRlb_2");
      systematicTree->makeOutputVariable(m_tma_mll,            "tma_mll");
      systematicTree->makeOutputVariable(m_tma_pTll,           "tma_pTll");
      systematicTree->makeOutputVariable(m_tma_dRll,           "tma_dRll");
      systematicTree->makeOutputVariable(m_tma_mbb,            "tma_mbb");
      systematicTree->makeOutputVariable(m_tma_pTbb,           "tma_pTbb"); 
      systematicTree->makeOutputVariable(m_tma_dRbb,           "tma_dRbb");
      systematicTree->makeOutputVariable(m_tma_Rbq_avgLJ,      "tma_Rbq_avgLJ");
      systematicTree->makeOutputVariable(m_tma_Rbq_leadLJ,     "tma_Rbq_leadLJ");

    }

  if ( topConfig()->doTopParticleLevel() ){

      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_met,                  "tma_met");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_etdr,                 "tma_etdr");

      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_nbjets,               "tma_nbjets");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_bjet_pt,              "tma_bjet_pt");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_bjet_eta,             "tma_bjet_eta");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_bjet_phi,             "tma_bjet_phi");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_bjet_e,               "tma_bjet_e");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_njets,                "tma_njets");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_jet_pt,               "tma_jet_pt");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_jet_eta,              "tma_jet_eta");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_jet_phi,              "tma_jet_phi");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_jet_e,                "tma_jet_e");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_mlb_minavg,           "tma_mlb_minavg");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_mlb_minavglow,        "tma_mlb_minavglow");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_mlb_minavghigh,       "tma_mlb_minavghigh");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_mlb_minmax,           "tma_mlb_minmax");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_mlb_minmaxlow,        "tma_mlb_minmaxlow");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_mlb_minmaxavg,        "tma_mlb_minmaxavg");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_pTlb_1,               "tma_pTlb_1");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_pTlb_2,               "tma_pTlb_2");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_dRlb_1,               "tma_dRlb_1");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_dRlb_2,               "tma_dRlb_2");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_mll,                  "tma_mll");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_pTll,                 "tma_pTll");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_dRll,                 "tma_dRll");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_mbb,                  "tma_mbb");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_pTbb,                 "tma_pTbb");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_dRbb,                 "tma_dRbb");

      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_lep_pt,               "tma_lep_pt");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_lep_eta,              "tma_lep_eta");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_lep_type,             "tma_lep_type");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_lep_charge,           "tma_lep_charge");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_lep_e,                "tma_lep_e");

      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_nelecs,               "tma_nelecs");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_el_pt,                "tma_lep_pt");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_el_eta,               "tma_lep_eta");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_el_e,                 "tma_lep_e");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_nmuons,               "tma_nmuons");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_mu_pt,                "tma_mu_pt");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_mu_eta,               "tma_mu_eta");
      particleLevelTreeManager()->makeOutputVariable(m_tma_particle_mu_e,                 "tma_mu_e");

  }
}
  
  ///-- saveEvent - run for every systematic and every event --///
  void MassEventSaver::saveEvent(const top::Event& event) 
  {
    ///-- set our variables to zero --///
    m_tma_nbjets = 0;
    m_tma_njets  = 0;    
    m_tma_lep_pt.clear();
    m_tma_lep_eta.clear();
    m_tma_lep_type.clear();
    m_tma_lep_charge.clear();
    m_tma_lep_e.clear();
    m_tma_el_e.clear();
    m_tma_mu_e.clear();

    // calculate here now the flat variables neccessary for the top mass analysis, want to avoid to recalculate them every time offline!
    
    // get vector of b-tagged jets
    for (const auto* const jetPtr : event.m_jets) {
      
      //if(jetPtr->isAvailable<char>("isbtagged_FixedCutBEff_77"))

      m_tma_njets++;

      const xAOD::BTagging* btag(nullptr);
      btag = jetPtr->btagging();	
      double mvx = -999;
      if (btag) btag->MVx_discriminant("MV2c10", mvx);

      if(mvx > 0.64)
	m_tma_nbjets++;
      
    }

    m_tma_bjet_pt.resize(m_tma_nbjets);
    m_tma_bjet_eta.resize(m_tma_nbjets);
    m_tma_bjet_phi.resize(m_tma_nbjets);
    m_tma_bjet_e.resize(m_tma_nbjets);

    m_tma_jet_pt.resize(m_tma_njets);
    m_tma_jet_eta.resize(m_tma_njets);
    m_tma_jet_phi.resize(m_tma_njets);
    m_tma_jet_e.resize(m_tma_njets);

    int i = 0;
    int j = 0;

    std::vector<TLorentzVector> bJets;
    std::vector<TLorentzVector> LJets;
    std::vector<double> LJets_mvx;
    for (const auto* const jetPtr : event.m_jets) {

      m_tma_jet_pt[j]	= jetPtr->pt();
      m_tma_jet_eta[j]	= jetPtr->eta();
      m_tma_jet_phi[j]	= jetPtr->phi();
      m_tma_jet_e[j]	= jetPtr->e();

      ++j;

      const xAOD::BTagging* btag(nullptr);
      btag = jetPtr->btagging();
      double mvx = -999;
      if(btag) btag->MVx_discriminant("MV2c10", mvx);
      if(mvx > 0.64){
	  
	m_tma_bjet_pt[i]  = jetPtr->pt();
	m_tma_bjet_eta[i] = jetPtr->eta();
	m_tma_bjet_phi[i] = jetPtr->phi();
	m_tma_bjet_e[i]   = jetPtr->e();
	   
	TLorentzVector help;
	help.SetPtEtaPhiE(jetPtr->pt(), jetPtr->eta(), jetPtr->phi(), jetPtr->e());
 
	bJets.push_back(help);

	++i;
      
      }
      else{
	
	TLorentzVector help;
	help.SetPtEtaPhiE(jetPtr->pt(), jetPtr->eta(), jetPtr->phi(), jetPtr->e());
 
	LJets.push_back(help);	
	LJets_mvx.push_back(mvx);	
	
      }

    }
      
    TLorentzVector blep_original;;
    TLorentzVector bhad_original;;
    TLorentzVector lq1_original;
    TLorentzVector lq2_original;
    TLorentzVector Whad_original;

    int nPermutations  = 0; 
    bool validKLFitter = false;

    if (event.m_KLFitterResults != nullptr) {
      
      validKLFitter = true;
      nPermutations = event.m_KLFitterResults->size();
      
    }

    m_tma_klfitter_mtop_param.resize(nPermutations);
    m_tma_original_mw.resize(nPermutations);
    m_tma_original_rbq.resize(nPermutations);

    int iPerm = 0;

    if (validKLFitter) {

      for (const auto* const klPtr : *event.m_KLFitterResults) {


	int blep_index = klPtr->model_blep_jetIndex();
	int bhad_index = klPtr->model_bhad_jetIndex();
	int lq1_index  = klPtr->model_lq1_jetIndex();
	int lq2_index  = klPtr->model_lq2_jetIndex();

	blep_original.SetPtEtaPhiE(event.m_jets.at(blep_index) -> pt()/1000.0,
				   event.m_jets.at(blep_index) -> eta(),
				   event.m_jets.at(blep_index) -> phi(),
				   event.m_jets.at(blep_index) -> e()/1000.0);

	bhad_original.SetPtEtaPhiE(event.m_jets.at(bhad_index) -> pt()/1000.0,
				   event.m_jets.at(bhad_index) -> eta(),
                                   event.m_jets.at(bhad_index) -> phi(),
				   event.m_jets.at(bhad_index) -> e()/1000.0);

	lq1_original.SetPtEtaPhiE(event.m_jets.at(lq1_index) -> pt()/1000.0,
				   event.m_jets.at(lq1_index) -> eta(),
                                   event.m_jets.at(lq1_index) -> phi(),
				   event.m_jets.at(lq1_index) -> e()/1000.0);

	lq2_original.SetPtEtaPhiE(event.m_jets.at(lq2_index) -> pt()/1000.0,
				   event.m_jets.at(lq2_index) -> eta(),
                                   event.m_jets.at(lq2_index) -> phi(),
				   event.m_jets.at(lq2_index) -> e()/1000.0);

	
	Whad_original = lq1_original + lq2_original;

	m_tma_original_mw[iPerm] = Whad_original.M();
	
	if(m_tma_nbjets == 1){
	  m_tma_original_rbq[iPerm] = m_tma_bjet_pt[0]/1000.0/((lq1_original.Pt()+lq2_original.Pt())*0.5);
	}
	else if(m_tma_nbjets >= 2){
	  m_tma_original_rbq[iPerm] = (m_tma_bjet_pt[0] + m_tma_bjet_pt[1])/1000.0/(lq1_original.Pt()+lq2_original.Pt());
	}

	int nr_param = klPtr->parameters().size();
	m_tma_klfitter_mtop_param[iPerm] = klPtr->parameters().at(nr_param-1);
	
	// std::cout << m_tma_original_mw[iPerm] << "\t" << m_tma_original_rbq[iPerm] << "\t" << m_tma_klfitter_mtop_param[iPerm] << std::endl;  
		
	++iPerm;
	
      }
      
    }

    
    std::vector<TLorentzVector> goodLeptons;
    std::vector<TLorentzVector> goodElectrons;
    std::vector<TLorentzVector> goodMuons;

    for (const auto* const elPtr : event.m_electrons) {

      TLorentzVector help;
      help.SetPtEtaPhiE(elPtr->pt(), elPtr->eta(), elPtr->phi(), elPtr->e());

      m_tma_el_pt.push_back(elPtr->pt()/1000.);
      m_tma_el_eta.push_back(elPtr->eta());
      m_tma_el_e.push_back(elPtr->e()/1000.);

      m_tma_lep_pt.push_back(elPtr->pt()/1000.);
      m_tma_lep_eta.push_back(elPtr->eta());
      m_tma_lep_e.push_back(elPtr->e()/1000.);
      m_tma_lep_type.push_back(11);
      m_tma_lep_charge.push_back(elPtr->charge());

      goodLeptons.push_back(help);
      goodElectrons.push_back(help);

    }
    for (const auto* const muPtr : event.m_muons) {

      TLorentzVector help;
      help.SetPtEtaPhiE(muPtr->pt(), muPtr->eta(), muPtr->phi(), muPtr->e());

      m_tma_mu_pt.push_back(muPtr->pt()/1000.);
      m_tma_mu_eta.push_back(muPtr->eta());
      m_tma_mu_e.push_back(muPtr->e()/1000.);
      
      m_tma_lep_pt.push_back(muPtr->pt()/1000.);
      m_tma_lep_eta.push_back(muPtr->eta());
      m_tma_lep_e.push_back(muPtr->e()/1000.);
      m_tma_lep_type.push_back(13);
      m_tma_lep_charge.push_back(muPtr->charge());

      goodLeptons.push_back(help);
      goodMuons.push_back(help);

    }

    m_tma_nleps  = goodLeptons.size();
    m_tma_nelecs = goodElectrons.size();
    m_tma_nmuons = goodMuons.size();

    m_tma_met = event.m_met->met()/1000.;

    if(m_tma_nbjets >= 2){ 

      // Sort b-tagged jets by pT (*BP* or by MV2)  
      TLorentzVector b1;
      TLorentzVector b2;
      if(bJets[0].Pt() >= bJets[1].Pt()){
	b1 = bJets[0];
	b2 = bJets[1];
      }
      else{
	b1 = bJets[1];
	b2 = bJets[0];
      }
      
      m_tma_mbb  = (b1 + b2).M();
      m_tma_pTbb = (b1 + b2).Pt();
      m_tma_dRbb = b1.DeltaR(b2);

      // for dilepton channel only
      if(goodLeptons.size() == 2){

	// Sort leptons by pT  
	TLorentzVector L1;
	TLorentzVector L2;
	if(goodLeptons.at(0).Pt() >= goodLeptons.at(1).Pt()){
	  L1 = goodLeptons.at(0);
	  L2 = goodLeptons.at(1);
	}
	else{
	  L1 = goodLeptons.at(1);
	  L2 = goodLeptons.at(0);
	}

	//Pairing decision             
	TLorentzVector L1b1 = L1 + b1;
	TLorentzVector L2b2 = L2 + b2;
	TLorentzVector L1b2 = L1 + b2;
	TLorentzVector L2b1 = L2 + b1;
	double avgMass1 = (L1b1.M() + L2b2.M())/2 ;
	double avgMass2 = (L1b2.M() + L2b1.M())/2 ;
	
	TLorentzVector LBpair1;
	TLorentzVector LBpair2;
	TLorentzVector LBpair1_reject;
	TLorentzVector LBpair2_reject;
	double LBpair_avgMass = -1.0;
	double LBpair1_dR = -1.0;
	double LBpair2_dR = -1.0;
	
	if( avgMass1 <= avgMass2 ){
	  LBpair1 = L1b1;
	  LBpair2 = L2b2;
	  LBpair1_reject = L1b2;
	  LBpair2_reject = L2b1;
	  LBpair_avgMass = avgMass1;
	  LBpair1_dR = b1.DeltaR(L1);
	  LBpair2_dR = b2.DeltaR(L2);
	}
	else{
	  LBpair1 = L1b2;
	  LBpair2 = L2b1;
	  LBpair1_reject = L1b1;
	  LBpair2_reject = L2b2;
	  LBpair_avgMass = avgMass2;
	  LBpair1_dR = b2.DeltaR(L1);
	  LBpair2_dR = b1.DeltaR(L2);
	}

        m_tma_etdr = (L1.Et()*LBpair1_dR + L2.Et()*LBpair2_dR)/2.;
	m_tma_mlb_minavg     = LBpair_avgMass; 
	
	if(LBpair1.M() < LBpair2.M()){
	  m_tma_mlb_minavglow   = LBpair1.M();
	  m_tma_mlb_minavghigh  = LBpair2.M();
	}
	else{
	  m_tma_mlb_minavglow   = LBpair2.M();
	  m_tma_mlb_minavghigh  = LBpair1.M();
	}
	
	m_tma_pTlb_1 = LBpair1.Pt();
	m_tma_pTlb_2 = LBpair2.Pt();
	
	m_tma_dRlb_1 = LBpair1_dR;
	m_tma_dRlb_2 = LBpair2_dR;
	
	m_tma_mll    = (L1 + L2).M();
	m_tma_pTll   = (L1 + L2).Pt();
	m_tma_dRll   = L1.DeltaR(L2);

	double Rbq_avgLJ  = -1;
	double Rbq_leadLJ = -1;
	if(LJets.size() >= 1){
	  double bJets_sumPt = 0;
	  double LJets_sumPt = 0;
	  for(int i = 0; i < (int)bJets.size(); ++i)
	    bJets_sumPt += bJets.at(i).Pt();
	  for(int i = 0; i < (int)LJets.size(); ++i)
	    LJets_sumPt += LJets.at(i).Pt();
	  
	  Rbq_avgLJ  = (bJets_sumPt/bJets.size())/(LJets_sumPt/LJets.size());
	  Rbq_leadLJ = (bJets_sumPt/bJets.size())/LJets.at(0).Pt();
	}
	m_tma_Rbq_avgLJ  = Rbq_avgLJ;
	m_tma_Rbq_leadLJ = Rbq_leadLJ;


	//Alternate Pairing decision min(max)       
	TLorentzVector LBpair1_maxMass;
	TLorentzVector LBpair1_minMass;
	TLorentzVector LBpair2_maxMass;
	TLorentzVector LBpair2_minMass;
	if(L1b1.M() > L2b2.M()){
	  LBpair1_maxMass = L1b1;
	  LBpair1_minMass = L2b2;
	}
	else{
	  LBpair1_maxMass = L2b2;
	  LBpair1_minMass = L1b1;
	}
	if(L1b2.M() > L2b1.M()){
	  LBpair2_maxMass = L1b2;
	  LBpair2_minMass = L2b1;
	}
	else{
	  LBpair2_maxMass = L2b1;
	  LBpair2_minMass = L1b2;
	}	
	
	if(LBpair1_maxMass.M() < LBpair2_maxMass.M())
	  m_tma_mlb_minmax = LBpair1_maxMass.M();
	else
	  m_tma_mlb_minmax = LBpair2_maxMass.M();
	
	if(LBpair1_maxMass.M() < LBpair2_maxMass.M()){
	  m_tma_mlb_minmaxlow = LBpair1_minMass.M();
	}
	else{
	  m_tma_mlb_minmaxlow = LBpair2_minMass.M();
	}
	
	if(LBpair1_maxMass.M() < LBpair2_maxMass.M()){
	  m_tma_mlb_minmaxavg = (LBpair1_maxMass.M() + LBpair1_minMass.M())/2;
	}
	else{
	  m_tma_mlb_minmaxavg = (LBpair2_maxMass.M() + LBpair2_minMass.M())/2;
	}
	
	
      }
      
    }
    
    ///-- Let the base class do all the hard work --///
    top::EventSaverFlatNtuple::saveEvent(event);

  } 

void MassEventSaver::saveParticleLevelEvent(const top::ParticleLevelEvent& plEvent){

  if( !topConfig()->doTopParticleLevel() ){
	return;
  }

  m_tma_particle_nbjets = 0;
  m_tma_particle_bjet_pt.clear();
  m_tma_particle_bjet_eta.clear();
  m_tma_particle_bjet_phi.clear();
  m_tma_particle_bjet_e.clear();

  m_tma_particle_njets  = 0;
  m_tma_particle_jet_pt.clear();
  m_tma_particle_jet_eta.clear();
  m_tma_particle_jet_phi.clear();
  m_tma_particle_jet_e.clear();

  m_tma_particle_nelecs  = 0;
  m_tma_particle_el_pt.clear();
  m_tma_particle_el_eta.clear();
  m_tma_particle_el_e.clear();
  m_tma_particle_nmuons  = 0;
  m_tma_particle_mu_pt.clear();
  m_tma_particle_mu_eta.clear();
  m_tma_particle_mu_e.clear();

  std::vector<TLorentzVector> particle_bJets;
  std::vector<TLorentzVector> particle_goodLeptons;

  for (const auto & jetPtr : * plEvent.m_jets) {
    m_tma_particle_jet_pt.push_back(jetPtr->pt()/1000.);
    m_tma_particle_jet_eta.push_back(jetPtr->eta());
    m_tma_particle_jet_phi.push_back(jetPtr->phi());
    m_tma_particle_jet_e.push_back(jetPtr->e()/1000.);
    ++m_tma_particle_njets;

    if( jetPtr->auxdata<int>( "GhostBHadronsFinalCount" ) > 0){
      m_tma_particle_bjet_pt.push_back(jetPtr->pt()/1000.);
      m_tma_particle_bjet_eta.push_back(jetPtr->eta());
      m_tma_particle_bjet_phi.push_back(jetPtr->phi());
      m_tma_particle_bjet_e.push_back(jetPtr->e()/1000.);
      ++m_tma_particle_nbjets;

      TLorentzVector help;
      help.SetPtEtaPhiE(jetPtr->pt(), jetPtr->eta(), jetPtr->phi(), jetPtr->e());

      particle_bJets.push_back(help);
    }
  }

  for (const auto & elPtr : * plEvent.m_electrons) {

    m_tma_particle_el_pt.push_back(elPtr->pt()/1000.);
    m_tma_particle_el_eta.push_back(elPtr->eta());
    m_tma_particle_el_e.push_back(elPtr->e()/1000.);
    ++m_tma_particle_nelecs;

    m_tma_particle_lep_type.push_back(11);
    m_tma_particle_lep_pt.push_back(elPtr->pt()/1000.);
    m_tma_particle_lep_eta.push_back(elPtr->eta());
    m_tma_particle_lep_charge.push_back(elPtr->charge());
    ++m_tma_particle_nleps;

    TLorentzVector help;
    help.SetPtEtaPhiE(elPtr->pt(), elPtr->eta(), elPtr->phi(), elPtr->e());
    particle_goodLeptons.push_back(help);
  }

  for (const auto & muPtr : * plEvent.m_muons) {

    m_tma_particle_mu_pt.push_back(muPtr->pt()/1000.);
    m_tma_particle_mu_eta.push_back(muPtr->eta());
    m_tma_particle_mu_e.push_back(muPtr->e()/1000.);
    ++m_tma_particle_nmuons;

    m_tma_particle_lep_type.push_back(13);
    m_tma_particle_lep_pt.push_back(muPtr->pt()/1000.);
    m_tma_particle_lep_eta.push_back(muPtr->eta());
    m_tma_particle_lep_charge.push_back(muPtr->charge());
    ++m_tma_particle_nleps;

    TLorentzVector help;
    help.SetPtEtaPhiE(muPtr->pt(), muPtr->eta(), muPtr->phi(), muPtr->e());
    particle_goodLeptons.push_back(help);
  }

  m_tma_particle_met = plEvent.m_met->met()/1000.;

  if(m_tma_particle_nbjets >= 2){ 

      // Sort b-tagged jets by pT (*BP* or by MV2)  
      TLorentzVector pb1;
      TLorentzVector pb2;
      if(particle_bJets[0].Pt() >= particle_bJets[1].Pt()){
	pb1 = particle_bJets[0];
	pb2 = particle_bJets[1];
      }
      else{
	pb1 = particle_bJets[1];
	pb2 = particle_bJets[0];
      }
      
      m_tma_particle_mbb  = (pb1 + pb2).M();
      m_tma_particle_pTbb = (pb1 + pb2).Pt();
      m_tma_particle_dRbb = pb1.DeltaR(pb2);

      // for dilepton channel only
      if(particle_goodLeptons.size() == 2){

	// Sort leptons by pT  
	TLorentzVector pL1;
	TLorentzVector pL2;
	if(particle_goodLeptons.at(0).Pt() >= particle_goodLeptons.at(1).Pt()){
	  pL1 = particle_goodLeptons.at(0);
	  pL2 = particle_goodLeptons.at(1);
	}
	else{
	  pL1 = particle_goodLeptons.at(1);
	  pL2 = particle_goodLeptons.at(0);
	}

	//Pairing decision             
	TLorentzVector pL1b1 = pL1 + pb1;
	TLorentzVector pL2b2 = pL2 + pb2;
	TLorentzVector pL1b2 = pL1 + pb2;
	TLorentzVector pL2b1 = pL2 + pb1;
	double pavgMass1 = (pL1b1.M() + pL2b2.M())/2 ;
	double pavgMass2 = (pL1b2.M() + pL2b1.M())/2 ;
	
	TLorentzVector pLBpair1;
	TLorentzVector pLBpair2;
	TLorentzVector pLBpair1_reject;
	TLorentzVector pLBpair2_reject;
	double pLBpair_avgMass = -1.0;
	double pLBpair1_dR = -1.0;
	double pLBpair2_dR = -1.0;
	
	if( pavgMass1 <= pavgMass2 ){
	  pLBpair1 = pL1b1;
	  pLBpair2 = pL2b2;
	  pLBpair1_reject = pL1b2;
	  pLBpair2_reject = pL2b1;
	  pLBpair_avgMass = pavgMass1;
	  pLBpair1_dR = pb1.DeltaR(pL1);
	  pLBpair2_dR = pb2.DeltaR(pL2);
	}
	else{
	  pLBpair1 = pL1b2;
	  pLBpair2 = pL2b1;
	  pLBpair1_reject = pL1b1;
	  pLBpair2_reject = pL2b2;
	  pLBpair_avgMass = pavgMass2;
	  pLBpair1_dR = pb2.DeltaR(pL1);
	  pLBpair2_dR = pb1.DeltaR(pL2);
	}

        m_tma_particle_etdr = (pL1.Et()*pLBpair1_dR + pL2.Et()*pLBpair2_dR)/2.;
	m_tma_particle_mlb_minavg     = pLBpair_avgMass; 
	
	if(pLBpair1.M() < pLBpair2.M()){
	  m_tma_particle_mlb_minavglow   = pLBpair1.M();
	  m_tma_particle_mlb_minavghigh  = pLBpair2.M();
	}
	else{
	  m_tma_particle_mlb_minavglow   = pLBpair2.M();
	  m_tma_particle_mlb_minavghigh  = pLBpair1.M();
	}
	
	m_tma_particle_pTlb_1 = pLBpair1.Pt();
	m_tma_particle_pTlb_2 = pLBpair2.Pt();
	
	m_tma_particle_dRlb_1 = pLBpair1_dR;
	m_tma_particle_dRlb_2 = pLBpair2_dR;
	
	m_tma_particle_mll    = (pL1 + pL2).M();
	m_tma_particle_pTll   = (pL1 + pL2).Pt();
	m_tma_particle_dRll   = pL1.DeltaR(pL2);

     }
  }

  top::EventSaverFlatNtuple::saveParticleLevelEvent(plEvent);

}
}
