/*
  Copyright (C) 2002-2017 CERN for the benefit of the ATLAS collaboration
*/

#ifndef TOPMASS_13TEV_ONLINE_MASSEVENTSAVER_H
#define TOPMASS_13TEV_ONLINE_MASSEVENTSAVER_H

#include "TopAnalysis/EventSaverFlatNtuple.h"

/**
 * 
 * This class inherits from top::EventSaverFlatNtuple, which will be doing all the hard work 
 * 
 */

namespace top{
  class MassEventSaver : public top::EventSaverFlatNtuple {
    public:
      ///-- Default constrcutor with no arguments - needed for ROOT --///
      MassEventSaver();
      ///-- Destructor does nothing --///
      virtual ~MassEventSaver(){}
      
      ///-- initialize function for top::EventSaverFlatNtuple --///
      ///-- We will be setting up out custom variables here --///
      virtual void initialize(std::shared_ptr<top::TopConfig> config, TFile* file, const std::vector<std::string>& extraBranches) override;
      
      ///-- Keep the asg::AsgTool happy --///
      virtual StatusCode initialize() override {return StatusCode::SUCCESS;}      
      
      ///-- saveEvent function for top::EventSaverFlatNtuple --///
      ///-- We will be setting our custom variables on a per-event basis --///
      virtual void saveEvent(const top::Event& event) override;
      virtual void saveParticleLevelEvent(const top::ParticleLevelEvent& plEvent) override;
      
    private:
      ///-- Some additional custom variables for the output --///
      //float m_randomNumber;
      //float m_someOtherVariable;
      
      // "tma" always stands for "top mass analysis"

      std::vector<double> m_tma_klfitter_mtop_param;
      std::vector<double> m_tma_original_mw;
      std::vector<double> m_tma_original_rbq;

      int m_tma_nbjets;
      std::vector<double> m_tma_bjet_pt;
      std::vector<double> m_tma_bjet_eta;
      std::vector<double> m_tma_bjet_phi;
      std::vector<double> m_tma_bjet_e;


      // plots for dilepton channel
      double m_tma_etdr;
      double m_tma_met;

      int m_tma_nleps;
      std::vector<double> m_tma_lep_pt;
      std::vector<double> m_tma_lep_eta;
      std::vector<double> m_tma_lep_phi;
      std::vector<double> m_tma_lep_type;
      std::vector<double> m_tma_lep_charge;
      std::vector<double> m_tma_lep_e;
      int m_tma_nelecs;
      std::vector<double> m_tma_el_pt;
      std::vector<double> m_tma_el_eta;
      std::vector<double> m_tma_el_e;
      int m_tma_nmuons;
      std::vector<double> m_tma_mu_pt;
      std::vector<double> m_tma_mu_eta;
      std::vector<double> m_tma_mu_e;

      int m_tma_njets;
      std::vector<double> m_tma_jet_pt;
      std::vector<double> m_tma_jet_eta;
      std::vector<double> m_tma_jet_phi;
      std::vector<double> m_tma_jet_e;
      double m_tma_mlb_minavg;
      double m_tma_mlb_minavglow;
      double m_tma_mlb_minavghigh;
      double m_tma_mlb_minmax;
      double m_tma_mlb_minmaxlow;
      double m_tma_mlb_minmaxavg;
      double m_tma_pTlb_1;
      double m_tma_pTlb_2;
      double m_tma_dRlb_1;
      double m_tma_dRlb_2;
      double m_tma_mll;
      double m_tma_pTll;
      double m_tma_dRll;
      double m_tma_mbb;
      double m_tma_pTbb;
      double m_tma_dRbb;
      double m_tma_Rbq_avgLJ;
      double m_tma_Rbq_leadLJ;

      // particle level
      double m_tma_particle_etdr;
      double m_tma_particle_met;

      int m_tma_particle_nbjets;
      std::vector<double> m_tma_particle_bjet_pt;
      std::vector<double> m_tma_particle_bjet_eta;
      std::vector<double> m_tma_particle_bjet_phi;
      std::vector<double> m_tma_particle_bjet_e;
      int m_tma_particle_njets;
      std::vector<double> m_tma_particle_jet_pt;
      std::vector<double> m_tma_particle_jet_eta;
      std::vector<double> m_tma_particle_jet_phi;
      std::vector<double> m_tma_particle_jet_e;
      double m_tma_particle_mlb_minavg;
      double m_tma_particle_mlb_minavglow;
      double m_tma_particle_mlb_minavghigh;
      double m_tma_particle_mlb_minmax;
      double m_tma_particle_mlb_minmaxlow;
      double m_tma_particle_mlb_minmaxavg;
      double m_tma_particle_pTlb_1;
      double m_tma_particle_pTlb_2;
      double m_tma_particle_dRlb_1;
      double m_tma_particle_dRlb_2;
      double m_tma_particle_mll;
      double m_tma_particle_pTll;
      double m_tma_particle_dRll;
      double m_tma_particle_mbb;
      double m_tma_particle_pTbb;
      double m_tma_particle_dRbb;

      int m_tma_particle_nleps;
      std::vector<double> m_tma_particle_lep_pt;
      std::vector<double> m_tma_particle_lep_eta;
      std::vector<double> m_tma_particle_lep_type;
      std::vector<double> m_tma_particle_lep_charge;
      std::vector<double> m_tma_particle_lep_e;
      int m_tma_particle_nelecs;
      std::vector<double> m_tma_particle_el_pt;
      std::vector<double> m_tma_particle_el_eta;
      std::vector<double> m_tma_particle_el_e;
      int m_tma_particle_nmuons;
      std::vector<double> m_tma_particle_mu_pt;
      std::vector<double> m_tma_particle_mu_eta;
      std::vector<double> m_tma_particle_mu_e;

      ///-- Tell RootCore to build a dictionary (we need this) --///
      ClassDef(top::MassEventSaver, 0);
  };
}

#endif
