#include "TCanvas.h"
#include "TAxis.h"
#include "TCut.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TImage.h"
#include "TROOT.h"
#include "fstream"
#include "string"
#include "sstream"
#include "iostream"
#include "iomanip"
#include "TCut.h"
#include "vector"
#include "TLegend.h"
#include "TMath.h"
#include "TBranch.h"
#include "TFile.h"
#include "TTree.h" 
#include "THStack.h"
#include "TLatex.h"
#include "Drawing.h"
#include "TF1.h"
#include "Corrections.h"
//#include "array.h"



// get the distribution for the Z pt
TH1F* Z_Pt()
{
	TFile *inclusive = new TFile( "/afs/cern.ch/work/b/bbachu/private/Z_nunu/boosted-combo-pfmetraw-fj200.root");
	TTree *tree = new TTree();
	tree = (TTree*) inclusive->FindObjectAny( "Znunu_signal");
	Int_t event;
	Float_t genVpt;
	Float_t weight;
	Float_t jet1pt;
	Float_t mvamet;
	Int_t nentries = (Int_t) tree->GetEntries();
	//read all entries and fill the hist
	TH1F *hx = new TH1F( "Znunu_signal_genVpt", "Z Pt" , 75, 250, 1000);
	hx->Sumw2();
	tree->SetBranchAddress("weight", &weight);
	tree->SetBranchAddress("genVpt", &genVpt);
	tree->SetBranchAddress("mvamet", &mvamet);
	tree->SetBranchAddress("jet1pt", &jet1pt);
	for ( Int_t i = 0 ; i < nentries ; i++)
	{
		tree->GetEntry(i);
		if ( mvamet < 250) continue;
		if ( jet1pt < 200) continue;
		hx->Fill(genVpt , weight);
	}
return hx;
}

TH1F* Photon_Pt()
{
	TFile *inclusive = new TFile( "/afs/cern.ch/work/b/bbachu/private/Z_nunu/boosted-combo-pfmetraw-fj200.root");
	TTree *tree = new TTree();
	tree = (TTree*) inclusive->FindObjectAny( "Photon_photon_control");
	Int_t event;
	Float_t ptpho;
	Float_t mvamet;
	Float_t jet1pt;
	Float_t weight;
	Int_t nentries = (Int_t) tree->GetEntries();
	//read all entries and fill the hist
	TH1F *hx = new TH1F( "Photon_photon_control_genVpt", "#gamma Pt " , 75 , 250 , 1000);
	hx->Sumw2();
	tree->SetBranchAddress("weight", &weight);
	tree->SetBranchAddress("ptpho", &ptpho);
	tree->SetBranchAddress("mvamet", &mvamet);
	tree->SetBranchAddress("jet1pt", &jet1pt);
	for ( Int_t i = 0 ; i < nentries ; i++)
	{
		//for boosted we need mvamet > 250 and jet1pt >200
		tree->GetEntry(i);
		if ( mvamet < 250) continue;
		if ( jet1pt < 200) continue;
		hx->Fill(ptpho , weight  /* apply efficiency factor */  /* *0.971*/ );
	}
return hx;
}

TH1F* Rebinned_hist(TH1F* h)
{
	Double_t xbins[7] = {250, 300, 350, 400, 450, 500, 1000};
	TH1F* h_rebinned = (TH1F*) h->Rebin( 6  ,"H" , xbins);
	//plot the bin density inserad
	for (Int_t bin_no = 1 ; bin_no < 6 ; bin_no++)
	{
		Double_t Bin_Density = h_rebinned->GetBinContent(bin_no) / h_rebinned->GetBinWidth(bin_no) ;
		Double_t Error_Density = h_rebinned->GetBinError(bin_no) / h_rebinned->GetBinWidth(bin_no) ;
		h_rebinned->SetBinContent(bin_no , Bin_Density);
		h_rebinned->SetBinError(bin_no , Error_Density);
	}
return h_rebinned;
}

TH1F* Ratio_hist(TH1F* h1, TH1F* h2)
{
	Float_t xbins[7] = {250, 300, 350, 400, 450, 500, 1000};
	TH1F* h_ratio = new TH1F("Transfer_ratios", "Transfer_ratios" , 6 , xbins);
	h_ratio->Divide(h1, h2);
return h_ratio;
}

void transfer_test()
{
	cout <<"1"<< endl;
	TH1F* h_Z_pt = (TH1F*) Z_Pt();
	cout <<"2"<< endl;
	TH1F* h_Photon_Pt = (TH1F*) Photon_Pt();
	//rebin both histograms
	cout <<"1"<< endl;
	TH1F* h_ZPT_rebinned = (TH1F*) Rebinned_hist(h_Z_pt);
	cout <<"3"<< endl;
	TH1F* h_Photon_Pt_rebinned = (TH1F*) Rebinned_hist(h_Photon_Pt);
	//make the tranfer fucntion
	cout <<"4"<< endl;
	TH1F* transfer_ratios = (TH1F*) Ratio_hist(h_ZPT_rebinned , h_Photon_Pt_rebinned);
	transfer_ratios->SetMaximum(1.0); transfer_ratios->SetMinimum(0.0);
	cout <<"5"<< endl;
	TCanvas *c = new TCanvas();
	transfer_ratios->Draw();
}
