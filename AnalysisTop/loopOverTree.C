#include <TFile.h>
#include <TH1.h>
#include <TTree.h>

#include <iostream>
#include <string>
#include <vector>

void loopOverTree() {

  TCanvas *c1 = new TCanvas("c1", "c1",900,900);
  gStyle->SetOptStat(0);

  TFile *myFile = TFile::Open("output_mc.root");
  TTree *recoLevel = (TTree*)myFile->Get("nominal");
  TTree *partLevel = (TTree*)myFile->Get("particleLevel");

  TH2F *tMatrix_el_pt	= new TH2F("tMatrix_el_pt","",100,0,5.8e+05,100,0,5.8e+05);
  TH2F *tMatrix_el_eta	= new TH2F("tMatrix_el_eta","",100,-2.9,2.9,100,-2.9,2.9);
  TH2F *tMatrix_mu_pt	= new TH2F("tMatrix_mu_pt","",100,0,6.7e+05,100,0,6.7e+05);
  TH2F *tMatrix_mu_eta	= new TH2F("tMatrix_mu_eta","",100,-2.9,2.9,100,-2.9,2.9);
  TH2F *tMatrix_n_jets	= new TH2F("tMatrix_n_jets","",15,0,15,15,0,15);
/*  TH2F *tMatrix_jet_pt	= new TH2F("tMatrix_jet_pt","",100,0,1.04e+06,100,0,1.04e+06);
  TH2F *tMatrix_jet_eta	= new TH2F("tMatrix_jet_eta","",100,-2.9,2.9,100,-2.9,2.9);
  TH2F *tMatrix_met	= new TH2F("tMatrix_met","",100,0,6.5e+05,100,0,6.5e+05);
*/
  recoLevel->BuildIndex("runNumber", "eventNumber");
  partLevel->BuildIndex("runNumber", "eventNumber");

  UInt_t runNumber;
  ULong64_t eventNumber;

  recoLevel->SetBranchAddress("runNumber", &runNumber);
  recoLevel->SetBranchAddress("eventNumber", &eventNumber);

  vector<float> *bjet_pt_reco;
  recoLevel->SetBranchAddress("tma_bjet_pt", bjet_pt_reco);
  vector<float> *bjet_pt_part;
  partLevel->SetBranchAddress("tma_bjet_pt", bjet_pt_part);
  vector<float> *bjet_eta_reco;
  recoLevel->SetBranchAddress("tma_bjet_eta", bjet_eta_reco);
  vector<float> *bjet_eta_part;
  partLevel->SetBranchAddress("tma_bjet_eta", bjet_eta_part);
  vector<float> *bjet_phi_reco;
  recoLevel->SetBranchAddress("tma_bjet_phi", bjet_phi_reco);
  vector<float> *bjet_phi_part;
  partLevel->SetBranchAddress("tma_bjet_phi", bjet_phi_part);
  vector<float> *bjet_e_reco;
  recoLevel->SetBranchAddress("tma_bjet_e", bjet_e_reco);
  vector<float> *bjet_e_part;
  partLevel->SetBranchAddress("tma_bjet_e", bjet_e_part);
  vector<float> *nbjets_reco;
  recoLevel->SetBranchAddress("tma_nbjets", nbjets_reco);
  vector<float> *nbjets_part;
  partLevel->SetBranchAddress("tma_nbjets", nbjets_part);

  vector<float> *etdr_reco;
  recoLevel->SetBranchAddress("tma_etdr", etdr_reco);
  vector<float> *etdr_part;
  partLevel->SetBranchAddress("tma_etdr", etdr_part);
  vector<float> *met_reco;
  recoLevel->SetBranchAddress("tma_met", met_reco);
  vector<float> *met_part;
  partLevel->SetBranchAddress("tma_met", met_part);


  vector<float> *nleps_reco;
  recoLevel->SetBranchAddress("tma_nleps", nleps_reco);
  vector<float> *nleps_part;
  partLevel->SetBranchAddress("tma_nleps", nleps_part);
  vector<float> *lep_pt_reco;
  recoLevel->SetBranchAddress("tma_lep_pt", lep_pt_reco);
  vector<float> *lep_pt_part;
  partLevel->SetBranchAddress("tma_lep_pt", lep_pt_part);
  vector<float> *lep_eta_reco;
  recoLevel->SetBranchAddress("tma_lep_eta", lep_eta_reco);
  vector<float> *lep_eta_part;
  partLevel->SetBranchAddress("tma_lep_eta", lep_eta_part);
  vector<float> *lep_phi_reco;
  recoLevel->SetBranchAddress("tma_lep_phi", lep_phi_reco);
  vector<float> *lep_phi_part;
  partLevel->SetBranchAddress("tma_lep_phi", lep_phi_part);
  vector<float> *lep_type_reco;
  recoLevel->SetBranchAddress("tma_lep_type", lep_type_reco);
  vector<float> *lep_type_part;
  partLevel->SetBranchAddress("tma_lep_type", lep_type_part);
  vector<float> *lep_charge_reco;
  recoLevel->SetBranchAddress("tma_lep_charge", lep_charge_reco);
  vector<float> *lep_charge_part;
  partLevel->SetBranchAddress("tma_lep_charge", lep_charge_part);
  vector<float> *lep_e_reco;
  recoLevel->SetBranchAddress("tma_lep_e", lep_e_reco);
  vector<float> *lep_e_part;
  partLevel->SetBranchAddress("tma_lep_e", lep_e_part);

  vector<float> *el_pt_reco;
  recoLevel->SetBranchAddress("tma_el_pt", &el_pt_reco);
  vector<float> *el_pt_part;
  partLevel->SetBranchAddress("tma_el_pt", &el_pt_part);
  vector<float> *el_eta_reco;
  recoLevel->SetBranchAddress("tma_el_eta", &el_eta_reco);
  vector<float> *el_eta_part;
  partLevel->SetBranchAddress("tma_el_eta", &el_eta_part);
  vector<float> *mu_pt_reco;
  recoLevel->SetBranchAddress("tma_mu_pt", &mu_pt_reco);
  vector<float> *mu_pt_part;
  partLevel->SetBranchAddress("tma_mu_pt", &mu_pt_part);
  vector<float> *mu_eta_reco;
  recoLevel->SetBranchAddress("tma_mu_eta", &mu_eta_reco);
  vector<float> *mu_eta_part;
  partLevel->SetBranchAddress("tma_mu_eta", &mu_eta_part);

  vector<float> *jet_pt_reco;
  recoLevel->SetBranchAddress("tma_jet_pt", &jet_pt_reco);
  vector<float> *jet_pt_part;
  partLevel->SetBranchAddress("tma_jet_pt", &jet_pt_part);
  vector<float> *jet_eta_reco;
  recoLevel->SetBranchAddress("tma_jet_eta", &jet_eta_reco);
  vector<float> *jet_eta_part;
  partLevel->SetBranchAddress("tma_jet_eta", &jet_eta_part);
  vector<float> *jet_phi_reco;
  recoLevel->SetBranchAddress("tma_jet_phi", &jet_phi_reco);
  vector<float> *jet_phi_part;
  partLevel->SetBranchAddress("tma_jet_phi", &jet_phi_part);
  vector<float> *jet_e_reco;
  recoLevel->SetBranchAddress("tma_jet_e", &jet_e_reco);
  vector<float> *jet_e_part;
  partLevel->SetBranchAddress("tma_jet_e", &jet_e_part);
  vector<float> *njets_reco;
  recoLevel->SetBranchAddress("tma_njets", &njets_reco);
  vector<float> *njets_part;
  partLevel->SetBranchAddress("tma_njets", &njets_part);

  vector<float> *mlb_reco;
  recoLevel->SetBranchAddress("tma_mlb_minavg", &mlb_reco);
  vector<float> *mlb_part;
  partLevel->SetBranchAddress("tma_mlb_minavg", &mlb_part);
  vector<float> *pTlb1_reco;
  recoLevel->SetBranchAddress("tma_pTlb_1", &pTlb1_reco);
  vector<float> *pTlb1_part;
  partLevel->SetBranchAddress("tma_pTlb_1", &pTlb1_part);
  vector<float> *pTlb2_reco;
  recoLevel->SetBranchAddress("tma_pTlb_2", &pTlb2_reco);
  vector<float> *pTlb2_part;
  partLevel->SetBranchAddress("tma_pTlb_2", &pTlb2_part);
  vector<float> *dRlb1_reco;
  recoLevel->SetBranchAddress("tma_dRlb_1", &dRlb1_reco);
  vector<float> *dRlb1_part;
  partLevel->SetBranchAddress("tma_dRlb_1", &dRlb1_part);
  vector<float> *dRlb2_reco;
  recoLevel->SetBranchAddress("tma_dRlb_2", &dRlb2_reco);
  vector<float> *dRlb2_part;
  partLevel->SetBranchAddress("tma_dRlb_2", &dRlb2_part);

  vector<float> *mll_reco;
  recoLevel->SetBranchAddress("tma_mll", &mll_reco);
  vector<float> *mll_part;
  partLevel->SetBranchAddress("tma_mll", &mll_part);
  vector<float> *pTll_reco;
  recoLevel->SetBranchAddress("tma_pTll", &pTll_reco);
  vector<float> *pTll_part;
  partLevel->SetBranchAddress("tma_pTll", &pTll_part);

/*  vector<float> *jet_eta_reco;
  recoLevel->SetBranchAddress("jet_eta", &jet_eta_reco);
  vector<float> *jet_eta_part;
  partLevel->SetBranchAddress("jet_eta", &jet_eta_part);
  vector<float> *met_reco;
  recoLevel->SetBranchAddress("met_met", &met_reco);
  vector<float> *met_part;
  partLevel->SetBranchAddress("met_met", &met_part);
*/
  ULong64_t recoOnly;  // Fakes
  ULong64_t partOnly;  // Efficiency
  ULong64_t matchedFromReco;
  ULong64_t matchedFromPart; // Just to debug now
	
  for( Long64_t i=0; i < recoLevel->GetEntries(); ++i ) {

	recoLevel->GetEntry(i);
	// Is there a match in partLevel?
	if( partLevel->GetEntryWithIndex(runNumber, eventNumber) < 0 ) recoOnly++;
	else {
	  matchedFromReco++;
	  UInt_t size_el_pt_reco = el_pt_reco->size();
	  UInt_t size_el_pt_part = el_pt_part->size();
          UInt_t size_el_eta_reco = el_eta_reco->size();
          UInt_t size_el_eta_part = el_eta_part->size();
	  UInt_t size_mu_pt_reco = mu_pt_reco->size();
          UInt_t size_mu_pt_part = mu_pt_part->size();
          UInt_t size_mu_eta_reco = mu_eta_reco->size();
          UInt_t size_mu_eta_part = mu_eta_part->size();
	  UInt_t size_jet_pt_reco = jet_pt_reco->size();
          UInt_t size_jet_pt_part = jet_pt_part->size();

	  if( size_el_pt_reco == 1 && size_el_pt_part == 1 ) { 
	    Float_t el1_pt_reco = Float_t((*(el_pt_reco))[0]);
	    Float_t el1_pt_part = Float_t((*(el_pt_part))[0]);

	    tMatrix_el_pt->Fill(el1_pt_part,el1_pt_reco);
	  }
          if( size_el_eta_reco == 1 && size_el_eta_part == 1 ) {
            Float_t el1_eta_reco = Float_t((*(el_eta_reco))[0]);
            Float_t el1_eta_part = Float_t((*(el_eta_part))[0]);

            tMatrix_el_eta->Fill(el1_eta_part,el1_eta_reco);
          }

	  if( size_mu_pt_reco == 1 && size_mu_pt_part == 1 ) {
            Float_t mu1_pt_reco = Float_t((*(mu_pt_reco))[0]);
            Float_t mu1_pt_part = Float_t((*(mu_pt_part))[0]);

            tMatrix_mu_pt->Fill(mu1_pt_part,mu1_pt_reco);
          }
          if( size_mu_eta_reco == 1 && size_mu_eta_part == 1 ) {
            Float_t mu1_eta_reco = Float_t((*(mu_eta_reco))[0]);
            Float_t mu1_eta_part = Float_t((*(mu_eta_part))[0]);

            tMatrix_mu_eta->Fill(mu1_eta_part,mu1_eta_reco);
          }
	  Int_t n_jets_reco = 0;
	  Int_t n_jets_part = 0;
	  for( int j=0; j < size_jet_pt_reco; ++j) { ++n_jets_reco;}
	  for( int j=0; j < size_jet_pt_part; ++j) { ++n_jets_part;}
	  tMatrix_n_jets->Fill(n_jets_reco,n_jets_part);

//	  Long64_t entry_part = partLevel->GetEntryNumberWithIndex(runNumber, eventNumber);
//	  partLevel->GetEntry(i);
	}
  }

  // Normalization
  for( int xbins=0; xbins < 100; ++xbins) {
     Float_t sum = 0;
     for( int ybins=0; ybins < 100; ++ybins ) {  sum += tMatrix_el_pt->GetBinContent( tMatrix_el_pt->GetBin(xbins,ybins) ); }
     for( int ybins=0; ybins < 100 && sum != 0; ++ybins) { tMatrix_el_pt->SetBinContent( tMatrix_el_pt->GetBin(xbins,ybins), 
                                tMatrix_el_pt->GetBinContent( tMatrix_el_pt->GetBin(xbins,ybins) )/sum ); }
  }
  for( int xbins=0; xbins < 100; ++xbins) {
     Float_t sum = 0;
     for( int ybins=0; ybins < 100; ++ybins ) {  sum += tMatrix_el_eta->GetBinContent( tMatrix_el_eta->GetBin(xbins,ybins) ); }
     for( int ybins=0; ybins < 100 && sum != 0; ++ybins) { tMatrix_el_eta->SetBinContent( tMatrix_el_eta->GetBin(xbins,ybins),
                                tMatrix_el_eta->GetBinContent( tMatrix_el_eta->GetBin(xbins,ybins) )/sum ); }
  }
  for( int xbins=0; xbins < 100; ++xbins) {
     Float_t sum = 0;
     for( int ybins=0; ybins < 100; ++ybins ) {  sum += tMatrix_mu_pt->GetBinContent( tMatrix_mu_pt->GetBin(xbins,ybins) ); }
     for( int ybins=0; ybins < 100 && sum != 0; ++ybins) { tMatrix_mu_pt->SetBinContent( tMatrix_mu_pt->GetBin(xbins,ybins),
                                tMatrix_mu_pt->GetBinContent( tMatrix_mu_pt->GetBin(xbins,ybins) )/sum ); }
  }
  for( int xbins=0; xbins < 100; ++xbins) {
     Float_t sum = 0;
     for( int ybins=0; ybins < 100; ++ybins ) {  sum += tMatrix_mu_eta->GetBinContent( tMatrix_mu_eta->GetBin(xbins,ybins) ); }
     for( int ybins=0; ybins < 100 && sum != 0; ++ybins) { tMatrix_mu_eta->SetBinContent( tMatrix_mu_eta->GetBin(xbins,ybins),
                                tMatrix_mu_eta->GetBinContent( tMatrix_mu_eta->GetBin(xbins,ybins) )/sum ); }
  }
  for( int xbins=0; xbins < 15; ++xbins) {
     Float_t sum = 0;
     for( int ybins=0; ybins < 15; ++ybins ) {	sum += tMatrix_n_jets->GetBinContent( tMatrix_n_jets->GetBin(xbins,ybins) ); }
     for( int ybins=0; ybins < 15 && sum != 0; ++ybins) { tMatrix_n_jets->SetBinContent( tMatrix_n_jets->GetBin(xbins,ybins), 
				tMatrix_n_jets->GetBinContent( tMatrix_n_jets->GetBin(xbins,ybins) )/sum ); }
  }

  c1->cd();
  gStyle->SetPalette(1);
  tMatrix_el_pt->GetXaxis()->SetTitle("Electron p_{t}^{reco}");
  tMatrix_el_pt->GetYaxis()->SetTitle("Electron p_{t}^{part}");
  tMatrix_el_pt->Draw("COLZ");
  c1->SaveAs("el_pt.pdf");

  tMatrix_mu_pt->GetXaxis()->SetTitle("Muon p_{t}^{reco}");
  tMatrix_mu_pt->GetYaxis()->SetTitle("Muon p_{t}^{part}");
  tMatrix_mu_pt->Draw("COLZ");
  c1->SaveAs("mu_pt.pdf");

  tMatrix_el_eta->GetXaxis()->SetTitle("Electron \eta^{reco}");
  tMatrix_el_eta->GetYaxis()->SetTitle("Electron \eta^{part}");
  tMatrix_el_eta->Draw("COLZ");
  c1->SaveAs("el_eta.pdf");

  tMatrix_mu_eta->GetXaxis()->SetTitle("Muon \eta^{reco}");
  tMatrix_mu_eta->GetYaxis()->SetTitle("Muon \eta^{part}");
  tMatrix_mu_eta->Draw("COLZ");
  c1->SaveAs("mu_eta.pdf");

  tMatrix_n_jets->GetXaxis()->SetTitle("n_{jets}^{reco}");
  tMatrix_n_jets->GetYaxis()->SetTitle("n_{jets}^{part}");
  tMatrix_n_jets->Draw("COLZ");
  c1->SaveAs("n_jets.pdf");

/*
  recoLevel->ResetBranchAddresses();

  partLevel->SetBranchAddress("runNumber", &runNumber);
  partLevel->SetBranchAddress("eventNumber", &eventNumber);

  for( Long64_t i=0; i < partLevel->GetEntries(); ++i ) {

	partLevel->GetEntry(i);
	// Is there a match in recoLevel?
	if( recoLevel->GetEntryWithIndex(runNumber, eventNumber) < 0 ) partOnly++;
	else matchedFromPart++;
  }

  assert( (matchedFromReco + recoOnly) == recoLevel->GetEntries() );
  assert( (matchedFromPart + partOnly) == partLevel->GetEntries() );

  // Debug
//  std::cout << "Matched from Reco = " << matchedFromReco << "\n"
//	<< "Matched from Part = " << matchedFromPart << "\n"
//	<< "Only in Reco = " << recoOnly << "\n"
//	<< "Only in Part = " << partOnly << "\n"
//	<< "Total in Reco = " << recoLevel->GetEntries() << "\n"
//	<< "Total in Part = " << partLevel->GetEntries() << "\n";

*/
}
