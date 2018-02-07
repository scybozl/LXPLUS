#include <TFile.h>
#include <TH1.h>
#include <TTree.h>
#include <TGaxis.h>
#include <TStyle.h>

#include <iostream>
#include <string>
#include <vector>

void loopOverTree() {

  bool debug = false;

  TCanvas *c1 = new TCanvas("c1", "c1",900,900);
  gStyle->SetOptStat(0);

  TFile *myFile = TFile::Open("output_25K.root");
  TTree *recoLevel = (TTree*)myFile->Get("nominal");
  TTree *partLevel = (TTree*)myFile->Get("particleLevel");

  TH2F *tMatrix_el_pt	= new TH2F("tMatrix_el_pt","",100,0,580,100,0,580);
  TH1F *tPart_el_pt	= new TH1F("tPart_el_pt","",100,0,580);
  TH1F *tReco_el_pt	= new TH1F("tReco_el_pt","",100,0,580);
  TH1F *tEff_el_pt	= new TH1F("tEff_el_pt","",100,0,580);
  TH1F *tFacc_el_pt	= new TH1F("tFacc_el_pt","",100,0,580);

  TH2F *tMatrix_el_eta	= new TH2F("tMatrix_el_eta","",100,-2.9,2.9,100,-2.9,2.9);
  TH2F *tMatrix_mu_pt	= new TH2F("tMatrix_mu_pt","",100,0,6.7e+05,100,0,6.7e+05);
  TH2F *tMatrix_mu_eta	= new TH2F("tMatrix_mu_eta","",100,-2.9,2.9,100,-2.9,2.9);

  TH2F *tMatrix_njets	= new TH2F("tMatrix_njets","",10,0,10,10,0,10);
  TH1F *tPart_njets	= new TH1F("tPart_njets","",10,0,10);
  TH1F *tReco_njets	= new TH1F("tReco_njets","",10,0,10);
  TH1F *tEff_njets	= new TH1F("tEff_njets","",10,0,10);
  TH1F *tFacc_njets	= new TH1F("tFacc_njets","",10,0,10);

  TH2F *tMatrix_mlb	= new TH2F("tMatrix_mlb","",100,0,140e+03,100,0,140e+03);
  TH1F *tPart_matched	= new TH1F("tPart_matched","",10,0,10);
  TH1F *tReco_matched	= new TH1F("tReco_matched","",10,0,10);
  TH1F *tPart_mlb	= new TH1F("tPart_mlb","",100,0,140e+03);
  TH1F *tReco_mlb	= new TH1F("tReco_mlb","",100,0,140e+03);
  TH1F *tEff_mlb	= new TH1F("tEff_mlb","",100,0,140e+03);
  TH1F *tFacc_mlb	= new TH1F("tFacc_mlb","",100,0,140e+03);


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
  recoLevel->SetBranchAddress("tma_bjet_pt", &bjet_pt_reco);
  vector<float> *bjet_pt_part;
  partLevel->SetBranchAddress("tma_bjet_pt", &bjet_pt_part);
  vector<float> *bjet_eta_reco;
  recoLevel->SetBranchAddress("tma_bjet_eta", &bjet_eta_reco);
  vector<float> *bjet_eta_part;
  partLevel->SetBranchAddress("tma_bjet_eta", &bjet_eta_part);
  vector<float> *bjet_phi_reco;
  recoLevel->SetBranchAddress("tma_bjet_phi", &bjet_phi_reco);
  vector<float> *bjet_phi_part;
  partLevel->SetBranchAddress("tma_bjet_phi", &bjet_phi_part);
  vector<float> *bjet_e_reco;
  recoLevel->SetBranchAddress("tma_bjet_e", &bjet_e_reco);
  vector<float> *bjet_e_part;
  partLevel->SetBranchAddress("tma_bjet_e", &bjet_e_part);
  Int_t nbjets_reco;
  recoLevel->SetBranchAddress("tma_nbjets", &nbjets_reco);
  Int_t nbjets_part;
  partLevel->SetBranchAddress("tma_nbjets", &nbjets_part);

  Double_t etdr_reco;
  recoLevel->SetBranchAddress("tma_etdr", etdr_reco);
  Double_t etdr_part;
  partLevel->SetBranchAddress("tma_etdr", etdr_part);
  Double_t met_reco;
  recoLevel->SetBranchAddress("tma_met", met_reco);
  Double_t met_part;
  partLevel->SetBranchAddress("tma_met", met_part);


//  Int_t *nleps_reco;
//  recoLevel->SetBranchAddress("tma_nleps", nleps_reco);
//  Int_t *nleps_part;
//  partLevel->SetBranchAddress("tma_nleps", nleps_part);
  vector<float> *lep_pt_reco;
  recoLevel->SetBranchAddress("tma_lep_pt", &lep_pt_reco);
  vector<float> *lep_pt_part;
  partLevel->SetBranchAddress("tma_lep_pt", &lep_pt_part);
  vector<float> *lep_eta_reco;
  recoLevel->SetBranchAddress("tma_lep_eta", &lep_eta_reco);
  vector<float> *lep_eta_part;
  partLevel->SetBranchAddress("tma_lep_eta", &lep_eta_part);
//  vector<float> *lep_phi_reco;
//  recoLevel->SetBranchAddress("tma_lep_phi", &lep_phi_reco);
//  vector<float> *lep_phi_part;
//  partLevel->SetBranchAddress("tma_lep_phi", &lep_phi_part);
  vector<int> *lep_type_reco;
  recoLevel->SetBranchAddress("tma_lep_type", &lep_type_reco);
  vector<int> *lep_type_part;
  partLevel->SetBranchAddress("tma_lep_type", &lep_type_part);
  vector<int> *lep_charge_reco;
  recoLevel->SetBranchAddress("tma_lep_charge", &lep_charge_reco);
  vector<int> *lep_charge_part;
  partLevel->SetBranchAddress("tma_lep_charge", &lep_charge_part);
  vector<float> *lep_e_reco;
  recoLevel->SetBranchAddress("tma_lep_e", &lep_e_reco);
  vector<float> *lep_e_part;
  partLevel->SetBranchAddress("tma_lep_e", &lep_e_part);

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
  Int_t njets_reco;
  recoLevel->SetBranchAddress("tma_njets", &njets_reco);
  Int_t njets_part;
  partLevel->SetBranchAddress("tma_njets", &njets_part);

  Double_t mlb_reco;
  recoLevel->SetBranchAddress("tma_mlb_minavg", &mlb_reco);
  Double_t mlb_part;
  partLevel->SetBranchAddress("tma_mlb_minavg", &mlb_part);
  Double_t pTlb1_reco;
  recoLevel->SetBranchAddress("tma_pTlb_1", &pTlb1_reco);
  Double_t pTlb1_part;
  partLevel->SetBranchAddress("tma_pTlb_1", &pTlb1_part);
  Double_t pTlb2_reco;
  recoLevel->SetBranchAddress("tma_pTlb_2", &pTlb2_reco);
  Double_t pTlb2_part;
  partLevel->SetBranchAddress("tma_pTlb_2", &pTlb2_part);
  Double_t dRlb1_reco;
  recoLevel->SetBranchAddress("tma_dRlb_1", &dRlb1_reco);
  Double_t dRlb1_part;
  partLevel->SetBranchAddress("tma_dRlb_1", &dRlb1_part);
  Double_t dRlb2_reco;
  recoLevel->SetBranchAddress("tma_dRlb_2", &dRlb2_reco);
  Double_t dRlb2_part;
  partLevel->SetBranchAddress("tma_dRlb_2", &dRlb2_part);

  Double_t mll_reco;
  recoLevel->SetBranchAddress("tma_mll", &mll_reco);
  Double_t mll_part;
  partLevel->SetBranchAddress("tma_mll", &mll_part);
  Double_t pTll_reco;
  recoLevel->SetBranchAddress("tma_pTll", &pTll_reco);
  Double_t pTll_part;
  partLevel->SetBranchAddress("tma_pTll", &pTll_part);
  Double_t dRll_reco;
  recoLevel->SetBranchAddress("tma_dRll", &dRll_reco);
  Double_t dRll_part;
  partLevel->SetBranchAddress("tma_dRll", &dRll_part);
  Double_t mbb_reco;
  recoLevel->SetBranchAddress("tma_mbb", &mbb_reco);
  Double_t mbb_part;
  partLevel->SetBranchAddress("tma_mbb", &mbb_part);
  Double_t pTbb_reco;
  recoLevel->SetBranchAddress("tma_pTbb", &pTbb_reco);
  Double_t pTbb_part;
  partLevel->SetBranchAddress("tma_pTbb", &pTbb_part);
  Double_t dRbb_reco;
  recoLevel->SetBranchAddress("tma_dRbb", &dRbb_reco);
  Double_t dRbb_part;
  partLevel->SetBranchAddress("tma_dRbb", &dRbb_part);


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
	tReco_el_pt->Fill((*(el_pt_reco))[0]);
	tReco_mlb->Fill(mlb_reco);
	tReco_njets->Fill(njets_reco);
	// Is there a match in partLevel?
	if( partLevel->GetEntryWithIndex(runNumber, eventNumber) < 0 ) {
	  recoOnly++;
	}
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


	  if( size_el_pt_reco >= 1 && size_el_pt_part >= 1 ) { 
	    Float_t el1_pt_reco = Float_t((*(el_pt_reco))[0]);
	    Float_t el1_pt_part = Float_t((*(el_pt_part))[0]);

	    tEff_el_pt->Fill((*(el_pt_part))[0]);
	    tFacc_el_pt->Fill((*(el_pt_reco))[0]);
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
	  tEff_njets->Fill(njets_part);
	  tFacc_njets->Fill(njets_reco);
	  tMatrix_njets->Fill(njets_part,njets_reco);

	  tReco_matched->Fill(njets_reco);
	  tPart_matched->Fill(njets_part);
	  tEff_mlb->Fill(mlb_part);
	  tFacc_mlb->Fill(mlb_reco);
	  tMatrix_mlb->Fill(mlb_part,mlb_reco);

//	  Long64_t entry_part = partLevel->GetEntryNumberWithIndex(runNumber, eventNumber);
//	  partLevel->GetEntry(i);

	  if( i % 100 == 0) cout << "Processing event number " << i << endl;
	  if( debug == true && i >= 100) break;
	}
  }

  cout << "----- STARTING NORMALIZATION -----" << endl;

  // Normalization
  for( int xbins=0; xbins <= 100; ++xbins) {
     Float_t sum = 0;
     for( int ybins=0; ybins <= 100; ++ybins ) {  sum += tMatrix_el_pt->GetBinContent( tMatrix_el_pt->GetBin(xbins,ybins) ); }
     for( int ybins=0; ybins <= 100 && sum != 0; ++ybins) { tMatrix_el_pt->SetBinContent( tMatrix_el_pt->GetBin(xbins,ybins), 
                                tMatrix_el_pt->GetBinContent( tMatrix_el_pt->GetBin(xbins,ybins) )/Float_t(sum) ); }
  }
  for( int xbins=0; xbins <= 100; ++xbins) {
     Float_t sum = 0;
     for( int ybins=0; ybins <= 100; ++ybins ) {  sum += tMatrix_el_eta->GetBinContent( tMatrix_el_eta->GetBin(xbins,ybins) ); }
     for( int ybins=0; ybins <= 100 && sum != 0; ++ybins) { tMatrix_el_eta->SetBinContent( tMatrix_el_eta->GetBin(xbins,ybins),
                                tMatrix_el_eta->GetBinContent( tMatrix_el_eta->GetBin(xbins,ybins) )/sum ); }
  }
  for( int xbins=0; xbins <= 100; ++xbins) {
     Float_t sum = 0;
     for( int ybins=0; ybins <= 100; ++ybins ) {  sum += tMatrix_mu_pt->GetBinContent( tMatrix_mu_pt->GetBin(xbins,ybins) ); }
     for( int ybins=0; ybins <= 100 && sum != 0; ++ybins) { tMatrix_mu_pt->SetBinContent( tMatrix_mu_pt->GetBin(xbins,ybins),
                                tMatrix_mu_pt->GetBinContent( tMatrix_mu_pt->GetBin(xbins,ybins) )/sum ); }
  }
  for( int xbins=0; xbins <= 100; ++xbins) {
     Float_t sum = 0;
     for( int ybins=0; ybins <= 100; ++ybins ) {  sum += tMatrix_mu_eta->GetBinContent( tMatrix_mu_eta->GetBin(xbins,ybins) ); }
     for( int ybins=0; ybins <= 100 && sum != 0; ++ybins) { tMatrix_mu_eta->SetBinContent( tMatrix_mu_eta->GetBin(xbins,ybins),
                                tMatrix_mu_eta->GetBinContent( tMatrix_mu_eta->GetBin(xbins,ybins) )/sum ); }
  }
  for( int xbins=0; xbins <= 10; ++xbins) {
     Float_t sum = 0;
     for( int ybins=0; ybins <= 10; ++ybins ) {	sum += tMatrix_njets->GetBinContent( tMatrix_njets->GetBin(xbins,ybins) ); }
     for( int ybins=0; ybins <= 10 && sum != 0; ++ybins) { tMatrix_njets->SetBinContent( tMatrix_njets->GetBin(xbins,ybins), 
				tMatrix_njets->GetBinContent( tMatrix_njets->GetBin(xbins,ybins) )/Float_t(sum) ); }
  }
  for( int xbins=0; xbins <= 100; ++xbins) {
     Float_t sum = 0;
     for( int ybins=0; ybins <= 100; ++ybins ) {	sum += tMatrix_mlb->GetBinContent( tMatrix_mlb->GetBin(xbins,ybins) ); }
     for( int ybins=0; ybins <= 100 && sum != 0; ++ybins) { tMatrix_mlb->SetBinContent( tMatrix_mlb->GetBin(xbins,ybins), 
				tMatrix_mlb->GetBinContent( tMatrix_mlb->GetBin(xbins,ybins) )/Float_t(sum) ); }
  }
  cout << "----- NORMALIZATION FINISHED -----" << endl;

  c1->cd();
  gStyle->SetPalette(1);
  tMatrix_mlb->GetXaxis()->SetTitle("m_{lb}");
  tMatrix_mlb->GetYaxis()->SetTitle("m_{lb}");
  tMatrix_mlb->Draw("COLZ");
  c1->SaveAs("mlb.pdf");

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

  tMatrix_njets->GetXaxis()->SetTitle("n_{jets}^{reco}");
  tMatrix_njets->GetYaxis()->SetTitle("n_{jets}^{part}");
  tMatrix_njets->Draw("COLZ");
  c1->SaveAs("njets.pdf");

  cout << "----- MIGRATION MATRICES SAVED -----" << endl;

  recoLevel->ResetBranchAddresses();

  partLevel->SetBranchAddress("runNumber", &runNumber);
  partLevel->SetBranchAddress("eventNumber", &eventNumber);

  for( Long64_t i=0; i < partLevel->GetEntries(); ++i ) {

	partLevel->GetEntry(i);
	if (el_pt_part->size() >= 1) tPart_el_pt->Fill((*(el_pt_part))[0]);
	tPart_njets->Fill(njets_part);
	tPart_mlb->Fill(mlb_part);
//	tPart_el_pt->Fill(Float_t((*(el_pt_part))[0]));
//	tPart_el_eta->Fill(Float_t((*(el_eta_part))[0]));
	// Is there a match in recoLevel?
	if( recoLevel->GetEntryWithIndex(runNumber, eventNumber) < 0 ) {
	    partOnly++;
	}
	else matchedFromPart++;
  }

  cout << "----- PARTICLE LEVEL BRANCHES LOOPED THROUGH -----" << endl;

  // Debug

//  std::cout << "Total eff mlb entries: " << tEff_mlb->GetEntries() << "\n"
//	<< "Total facc mlb entries: " << tFacc_mlb->GetEntries() << "\n"
//	<< "Total part mlb entries: " << tPart_mlb->GetEntries() << endl;

//  tEff_el_pt->Divide(tPart_el_pt);
//  tFacc_el_pt->Divide(tReco_el_pt);
//  tEff_el_eta->Divide(tPart_el_eta);
//  tFacc_el_eta->Divide(tReco_el_eta);
  tEff_el_pt->Divide(tPart_el_pt);
  tFacc_el_pt->Divide(tReco_el_pt);
  tEff_njets->Divide(tPart_njets);
  tFacc_njets->Divide(tReco_njets);
  tEff_mlb->Divide(tPart_mlb);
  tFacc_mlb->Divide(tReco_mlb);

  tEff_el_pt->GetXaxis()->SetTitle("el p_T");
  tEff_el_pt->GetYaxis()->SetTitle("\epsilon_{eff} el p_T");
  tEff_el_pt->Draw();
  c1->SaveAs("efficiency_el_pt.pdf");

  tFacc_el_pt->GetXaxis()->SetTitle("el p_T");
  tFacc_el_pt->GetYaxis()->SetTitle("f_{acc} el p_T");
  tFacc_el_pt->Draw();
  c1->SaveAs("facc_el_pt.pdf");
/*
  tEff_el_eta->GetXaxis()->SetTitle("el \eta");
  tEff_el_eta->GetYaxis()->SetTitle("\epsilon_{eff} el \eta");
  tEff_el_eta->Draw();
  c1->SaveAs("efficiency_el_eta.pdf");

  tFacc_el_eta->GetXaxis()->SetTitle("el \eta");
  tFacc_el_eta->GetYaxis()->SetTitle("f_{acc} el \eta");
  tFacc_el_eta->Draw();
  c1->SaveAs("facc_el_eta.pdf");
*/

  tEff_njets->GetXaxis()->SetTitle("n_{jets}");
  tEff_njets->GetYaxis()->SetTitle("\epsilon_{eff} n_{jets}");
  tEff_njets->Draw();
  c1->SaveAs("efficiency_njets.pdf");

  tFacc_njets->GetXaxis()->SetTitle("n_{jets}");
  tFacc_njets->GetYaxis()->SetTitle("f_{acc} n_{jets}");
  tFacc_njets->Draw();
  c1->SaveAs("facc_njets.pdf");

  tEff_mlb->GetXaxis()->SetTitle("m_{lb}");
  tEff_mlb->GetYaxis()->SetTitle("\epsilon_{eff} m_{lb}");
  tEff_mlb->Draw();
  c1->SaveAs("efficiency_mlb.pdf");

  tFacc_mlb->GetXaxis()->SetTitle("m_{lb}");
  tFacc_mlb->GetYaxis()->SetTitle("f_{acc} m_{lb}");
  tFacc_mlb->Draw();
  c1->SaveAs("facc_mlb.pdf");

  cout << "----- EFFICIENCIES & FAKES CALCULATED -----" << endl;

  assert( (matchedFromReco + recoOnly) == recoLevel->GetEntries() );
  assert( (matchedFromPart + partOnly) == partLevel->GetEntries() );

  // Check that the matched events are correctly reconstructed

  TH1F *tReco_matched_computed       = new TH1F("tReco_matched_computed","",10,0,10);
  for ( int xbins = 0; xbins <= 10; ++xbins ) {
	float sum = 0;
	for (int j = 0; j <= 10; ++j) {sum += tMatrix_njets->GetBinContent(tMatrix_njets->GetBin(j, xbins))*tPart_matched->GetBinContent(j); }
	tReco_matched_computed->SetBinContent(xbins, sum);
	std::cout << "Reco reconstr: " << sum << "\nReco truth: " << tReco_matched->GetBinContent(xbins) << "\n";
  }


   // Define the Canvas
   TCanvas *c = new TCanvas("c", "canvas", 800, 800);

  // Plot both the original reco and the computed one
  // Upper plot will be in pad1
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   pad1->SetBottomMargin(0); // Upper and lower plot are joined
   pad1->SetGridx();         // Vertical grid
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();               // pad1 becomes the current pad
   tReco_matched_computed->SetStats(0);          // No statistics on upper plot
   tReco_matched_computed->Draw();               // Draw tReco_clos_mlb
   tReco_matched->Draw("same");         // Draw tReco_mlb on top of tReco_clos_mlb
   
   // Do not draw the Y axis label on the upper plot and redraw a small
   // axis instead, in order to avoid the first label (0) to be clipped.
   tReco_matched_computed->GetYaxis()->SetLabelSize(0.);
   TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   axis->SetLabelSize(10);
   axis->Draw();

   // lower plot will be in pad
   c->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd();       // pad2 becomes the current pad

   // Define the ratio plot
   TH1F *h3 = (TH1F*)tReco_matched_computed->Clone("h3"); 
   h3->SetLineColor(kBlack);
   h3->SetMinimum(0.8);  // Define Y ..
   h3->SetMaximum(1.35); // .. range
   h3->Sumw2();
   h3->SetStats(0);      // No statistics on lower plot
   h3->Divide(tReco_matched);
   h3->SetMarkerStyle(21);
   h3->Draw("ep");       // Draw the ratio plot
   
   // tReco_clos_mlb settings
   tReco_matched_computed->SetLineColor(kBlue+1);
   tReco_matched_computed->SetLineWidth(2);

   // Y axis tReco_clos_mlb plot settings
   tReco_matched_computed->GetYaxis()->SetTitleSize(20);
   tReco_matched_computed->GetYaxis()->SetTitleFont(43);
   tReco_matched_computed->GetYaxis()->SetTitleOffset(1.55);

   // tReco_mlb settings
   tReco_matched->SetLineColor(kRed);
   tReco_matched->SetLineWidth(2);

   // Ratio plot (h3) settings
   h3->SetTitle(""); // Remove the ratio title

   // Y axis ratio plot settings
   h3->GetYaxis()->SetTitle("ratio tReco_clos_mlb/tReco_mlb ");
   h3->GetYaxis()->SetNdivisions(505);
   h3->GetYaxis()->SetTitleSize(20);
   h3->GetYaxis()->SetTitleFont(43);
   h3->GetYaxis()->SetTitleOffset(1.55);
   h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h3->GetYaxis()->SetLabelSize(10);

   // X axis ratio plot settings
   h3->GetXaxis()->SetTitleSize(20);
   h3->GetXaxis()->SetTitleFont(43);
   h3->GetXaxis()->SetTitleOffset(4.);
   h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h3->GetXaxis()->SetLabelSize(10);

   c->SaveAs("closure_matched_njets.pdf");

   cout << "----- MATCHING CLOSURE TEST FINISHED -----" << endl;

  // Debug
//  std::cout << "Matched from Reco = " << matchedFromReco << "\n"
//	<< "Matched from Part = " << matchedFromPart << "\n"
//	<< "Only in Reco = " << recoOnly << "\n"
//	<< "Only in Part = " << partOnly << "\n"
//	<< "Total in Reco = " << recoLevel->GetEntries() << "\n"
//	<< "Total in Part = " << partLevel->GetEntries() << "\n";

  // Closure test: reco_i = (A_ij*part_j)*eff_i + (sum_j part_j)*facc_i

  TH1F *tReco_clos_njets       = new TH1F("tReco_clos_njets","",10,0,10);
  for ( int xbins = 0; xbins <=10; ++xbins ) {
	float sum = 0;
	for (int j = 0; j <= 10; ++j) sum += tMatrix_njets->GetBinContent(tMatrix_njets->GetBin(j, xbins))*tPart_njets->GetBinContent(j)*tEff_njets->GetBinContent(j);
//	if (xbins==50) float temp = sum;
	if (tFacc_njets->GetBinContent(xbins) != 0 ) sum *= 1./(tFacc_njets->GetBinContent(xbins));
	tReco_clos_njets->SetBinContent(xbins, sum);

  }

	// Debug
//	std::cout << "Part level mlb for bin 50: " << tPart_mlb->GetBinContent(50) << "\n"
//		<< "Reco level after applying A_ij: " << temp << "\n"
//		<< "Efficiency for bin 50: " << tEff_mlb->GetBinContent(50) << ", reco = " << temp*tEff_mlb->GetBinContent(50) << "\n"
//		<< "Fake rate for bin 50: " << 1./(tFacc_mlb->GetBinContent(50)) << ", reco = " << "\n"
//		<< "Reconstructed: " << tReco_clos_mlb->GetBinContent(50) << "\n"
//		<< "Truth: " << tReco_mlb->GetBinContent(50) << std::endl;

   // Define the Canvas
   TCanvas *c = new TCanvas("c", "canvas", 800, 800);

  // Plot both the original reco and the computed one
  // Upper plot will be in pad1
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   pad1->SetBottomMargin(0); // Upper and lower plot are joined
   pad1->SetGridx();         // Vertical grid
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();               // pad1 becomes the current pad
   tReco_clos_njets->SetStats(0);          // No statistics on upper plot
   tReco_clos_njets->Draw();               // Draw tReco_clos_mlb
   tReco_njets->Draw("same");         // Draw tReco_mlb on top of tReco_clos_mlb
   
   // Do not draw the Y axis label on the upper plot and redraw a small
   // axis instead, in order to avoid the first label (0) to be clipped.
   tReco_clos_njets->GetYaxis()->SetLabelSize(0.);
   TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   axis->SetLabelSize(10);
   axis->Draw();

   // lower plot will be in pad
   c->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd();       // pad2 becomes the current pad

   // Define the ratio plot
   TH1F *h3 = (TH1F*)tReco_clos_njets->Clone("h3"); 
   h3->SetLineColor(kBlack);
   h3->SetMinimum(0.5);  // Define Y ..
   h3->SetMaximum(1.5); // .. range
   h3->Sumw2();
   h3->SetStats(0);      // No statistics on lower plot
   h3->Divide(tReco_njets);
   h3->SetMarkerStyle(21);
   h3->Draw("ep");       // Draw the ratio plot
   
   // tReco_clos_mlb settings
   tReco_clos_njets->SetLineColor(kBlue+1);
   tReco_clos_njets->SetLineWidth(2);

   // Y axis tReco_clos_mlb plot settings
   tReco_clos_njets->GetYaxis()->SetTitleSize(20);
   tReco_clos_njets->GetYaxis()->SetTitleFont(43);
   tReco_clos_njets->GetYaxis()->SetTitleOffset(1.55);

   // tReco_mlb settings
   tReco_njets->SetLineColor(kRed);
   tReco_njets->SetLineWidth(2);

   // Ratio plot (h3) settings
   h3->SetTitle(""); // Remove the ratio title

   // Y axis ratio plot settings
   h3->GetYaxis()->SetTitle("ratio tReco_clos_mlb/tReco_mlb ");
   h3->GetYaxis()->SetNdivisions(505);
   h3->GetYaxis()->SetTitleSize(20);
   h3->GetYaxis()->SetTitleFont(43);
   h3->GetYaxis()->SetTitleOffset(1.55);
   h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h3->GetYaxis()->SetLabelSize(10);

   // X axis ratio plot settings
   h3->GetXaxis()->SetTitleSize(20);
   h3->GetXaxis()->SetTitleFont(43);
   h3->GetXaxis()->SetTitleOffset(4.);
   h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h3->GetXaxis()->SetLabelSize(10);

   c->SaveAs("closure_ratio_njets.pdf");

   cout << "----- CLOSURE TEST FOR NJETS FINISHED -----" << endl;

  TH1F *tReco_clos_mlb       = new TH1F("tReco_clos_mlb","",100,0,140e+03);
  for ( int xbins = 0; xbins <=100; ++xbins ) {
        float sum = 0;
        for (int j = 0; j <= 100; ++j) sum += tMatrix_mlb->GetBinContent(tMatrix_mlb->GetBin(j, xbins))*tPart_mlb->GetBinContent(j)*tEff_mlb->GetBinContent(j);
//      if (xbins==50) float temp = sum;
        if (tFacc_mlb->GetBinContent(xbins) != 0 ) sum *= 1./(tFacc_mlb->GetBinContent(xbins));
        tReco_clos_mlb->SetBinContent(xbins, sum);

  }

   // Define the Canvas
   TCanvas *c = new TCanvas("c", "canvas", 800, 800);

  // Plot both the original reco and the computed one
  // Upper plot will be in pad1
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   pad1->SetBottomMargin(0); // Upper and lower plot are joined
   pad1->SetGridx();         // Vertical grid
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();               // pad1 becomes the current pad
   tReco_clos_mlb->SetStats(0);          // No statistics on upper plot
   tReco_clos_mlb->Draw();               // Draw tReco_clos_mlb
   tReco_mlb->Draw("same");         // Draw tReco_mlb on top of tReco_clos_mlb

   // Do not draw the Y axis label on the upper plot and redraw a small
   // axis instead, in order to avoid the first label (0) to be clipped.
   tReco_clos_mlb->GetYaxis()->SetLabelSize(0.);
   TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   axis->SetLabelSize(10);
   axis->Draw();

   // lower plot will be in pad
   c->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd();       // pad2 becomes the current pad

   // Define the ratio plot
   TH1F *h3 = (TH1F*)tReco_clos_mlb->Clone("h3");
   h3->SetLineColor(kBlack);
   h3->SetMinimum(0.5);  // Define Y ..
   h3->SetMaximum(1.5); // .. range
   h3->Sumw2();
   h3->SetStats(0);      // No statistics on lower plot
   h3->Divide(tReco_mlb);
   h3->SetMarkerStyle(21);
   h3->Draw("ep");       // Draw the ratio plot

   // tReco_clos_mlb settings
   tReco_clos_mlb->SetLineColor(kBlue+1);
   tReco_clos_mlb->SetLineWidth(2);

   // Y axis tReco_clos_mlb plot settings
   tReco_clos_mlb->GetYaxis()->SetTitleSize(20);
   tReco_clos_mlb->GetYaxis()->SetTitleFont(43);
   tReco_clos_mlb->GetYaxis()->SetTitleOffset(1.55);

   // tReco_mlb settings
   tReco_mlb->SetLineColor(kRed);
   tReco_mlb->SetLineWidth(2);

// Ratio plot (h3) settings
   h3->SetTitle(""); // Remove the ratio title

   // Y axis ratio plot settings
   h3->GetYaxis()->SetTitle("ratio tReco_clos_mlb/tReco_mlb ");
   h3->GetYaxis()->SetNdivisions(505);
   h3->GetYaxis()->SetTitleSize(20);
   h3->GetYaxis()->SetTitleFont(43);
   h3->GetYaxis()->SetTitleOffset(1.55);
   h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h3->GetYaxis()->SetLabelSize(10);

   // X axis ratio plot settings
   h3->GetXaxis()->SetTitleSize(20);
   h3->GetXaxis()->SetTitleFont(43);
   h3->GetXaxis()->SetTitleOffset(4.);
   h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h3->GetXaxis()->SetLabelSize(10);

   c->SaveAs("closure_ratio_mlb.pdf");
   cout << "----- CLOSURE TEST FOR MLB FINISHED -----" << endl;

   cout << "----- CLOSURE TEST WITH STAT. INDEP. SAMPLE -----" << endl;

   TFile *myFile = TFile::Open("output_25K.root");
   TTree *recoLevel = (TTree*)myFile->Get("nominal");
   TTree *partLevel = (TTree*)myFile->Get("particleLevel");   
   TH1F *mlb	 = new TH1F("mlb","",100,0,140e+03);

   Double_t mlb_ind;
   partLevel->SetBranchAddress("tma_mlb_minavg", &mlb_ind);

   for( Long64_t i=0; i < partLevel->GetEntries(); ++i ) {

        partLevel->GetEntry(i);
        mlb->Fill(mlb_ind);
   }

   TH1F *tclos_mlb       = new TH1F("tclos_mlb","",100,0,140e+03);
  for ( int xbins = 0; xbins <=100; ++xbins ) {
        float sum = 0;
        for (int j = 0; j <= 100; ++j) sum += tMatrix_mlb->GetBinContent(tMatrix_mlb->GetBin(j, xbins))*mlb->GetBinContent(j)*tEff_mlb->GetBinContent(j);
//      if (xbins==50) float temp = sum;
        if (tFacc_mlb->GetBinContent(xbins) != 0 ) sum *= 1./(tFacc_mlb->GetBinContent(xbins));
        tclos_mlb->SetBinContent(xbins, sum);

  }

   // Define the Canvas
   TCanvas *c = new TCanvas("c", "canvas", 800, 800);

  // Plot both the original reco and the computed one
  // Upper plot will be in pad1
   TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
   pad1->SetBottomMargin(0); // Upper and lower plot are joined
   pad1->SetGridx();         // Vertical grid
   pad1->Draw();             // Draw the upper pad: pad1
   pad1->cd();               // pad1 becomes the current pad
   tclos_mlb->SetStats(0);          // No statistics on upper plot
   tclos_mlb->Draw();               // Draw tReco_clos_mlb
   mlb->Draw("same");         // Draw tReco_mlb on top of tReco_clos_mlb

   // Do not draw the Y axis label on the upper plot and redraw a small
   // axis instead, in order to avoid the first label (0) to be clipped.
   tclos_mlb->GetYaxis()->SetLabelSize(0.);
   TGaxis *axis = new TGaxis( -5, 20, -5, 220, 20,220,510,"");
   axis->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   axis->SetLabelSize(10);
   axis->Draw();

   // lower plot will be in pad
   c->cd();          // Go back to the main canvas before defining pad2
   TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.3);
   pad2->SetTopMargin(0);
   pad2->SetBottomMargin(0.2);
   pad2->SetGridx(); // vertical grid
   pad2->Draw();
   pad2->cd();       // pad2 becomes the current pad

   // Define the ratio plot
   TH1F *h3 = (TH1F*)tclos_mlb->Clone("h3");
   h3->SetLineColor(kBlack);
   h3->SetMinimum(0.5);  // Define Y ..
   h3->SetMaximum(1.5); // .. range
   h3->Sumw2();
   h3->SetStats(0);      // No statistics on lower plot
   h3->Divide(mlb);
   h3->SetMarkerStyle(21);
   h3->Draw("ep");       // Draw the ratio plot

   // tReco_clos_mlb settings
   tclos_mlb->SetLineColor(kBlue+1);
   tclos_mlb->SetLineWidth(2);

   // Ratio plot (h3) settings
   h3->SetTitle(""); // Remove the ratio title

   // Y axis ratio plot settings
   h3->GetYaxis()->SetTitle("ratio tReco_clos_mlb/tReco_mlb ");
   h3->GetYaxis()->SetNdivisions(505);
   h3->GetYaxis()->SetTitleSize(20);
   h3->GetYaxis()->SetTitleFont(43);
   h3->GetYaxis()->SetTitleOffset(1.55);
   h3->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h3->GetYaxis()->SetLabelSize(10);

   // X axis ratio plot settings
   h3->GetXaxis()->SetTitleSize(20);
   h3->GetXaxis()->SetTitleFont(43);
   h3->GetXaxis()->SetTitleOffset(4.);
   h3->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
   h3->GetXaxis()->SetLabelSize(10);

   c->SaveAs("closure_independent_mlb.pdf");

}

