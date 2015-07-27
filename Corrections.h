#ifndef CORRECTIONS_H
#define CORRECTIONS_H
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
#include "TF1.h"


using namespace std;


//Given any two histograms, make another histogram representing the ratio of the histograms
// h1 / h2  MAKE SURE BOTH HISTOGRAMS HAVE THE SAME BINNING
TH1F* h_Ratio(TH1F* h1 , TH1F* h2 )
{
	//TCanvas *c = new TCanvas();
	cout << "Getting the ratio of " << h1->GetName() << " to " << h2->GetName() << endl;
	//first check to make sure both histograms have the same number of bins
	Int_t h1_Nbins = h1->GetNbinsX();
	Int_t h2_Nbins = h2->GetNbinsX();
	Double_t h_xmin = h1->GetXaxis()->GetXmin();
	Double_t h_xmax = h1->GetXaxis()->GetXmax();
	TH1F *h_ratio = new TH1F("hRatio", "Ratio of " + TString(h1->GetTitle())  + " to "  + TString(h2->GetTitle()) , h1_Nbins, h_xmin, h_xmax);
	h_ratio->Sumw2();
	if(!(h1_Nbins == h2_Nbins))
	{
		cout << "THE NUMBER OF BINS DO NOT MATCH" << endl;
	}
	else
	{
		//declare new histogram to display the bin by bin ratio
		h_ratio->Divide( h1, h2,1,1);
		// TCanvas *c = new TCanvas();
		// h_ratio->Draw();
	}
return h_ratio;
}

//Transfer factor to correct for the Pt difference
Double_t R(Float_t photon_Pt, TH1F *hratios )
{
	Int_t N = (photon_Pt/ 10);
	Int_t Bin_no = N+1 ;
	Double_t Pt_factor = hratios->GetBinContent(Bin_no);
	//cout << "Photon pt: " << photon_Pt << ". Scaled with factor: " << Pt_factor << endl;
return Pt_factor;
}

// TH1F* Apply_transfer_function_event_by_event( TString root_directory ,const char* filename , TH1F* h_transfer , Float_t &photon_Pt)
// {
// 	cout << "Applying corrections based on Photon Pt" << endl;
// 	// Start by looping over the Photon + Jet 
// 	TFile *inclusive = new TFile(root_directory);
// 	TTree *tree = new TTree();
// 	tree = (TTree*) inclusive->FindObjectAny( filename );
// 	Float_t mvamet;
// 	//Float_t photon_Pt;
// 	Int_t event; 
// 	if (TString(filename) == "Photon_photon_control" )
// 	{
// 		tree->SetBranchAddress("genVpt" , &photon_Pt);
// 	}
// 	else
// 	{
// 		tree->SetBranchAddress("ptpho" , &photon_Pt);
// 	}
// 	// tree->SetBranchAddress("mvamet" , &mvamet);
// 	//create hist to store the corrected final met
// 	TH1F *hfinal_MET = new TH1F( "mvamet","Met Estimation from " + TString(filename), 100, 0, 1000);
	
// 	//read all entries and fill the hist
// 	Long64_t nentries = (Int_t) tree->GetEntries();
	
// 	cout << "Starting the event by event correction of "<< filename << endl;
// 	for (Long64_t entry = 0; entry<nentries; entry++ )
// 	{
// 		tree->GetEntry(entry);
// 		//correct for difference in Z pt and Photon pt
// 		Double_t R_factor = R(photon_Pt, h_transfer);
// 		//cout  <<  entry << " : " << R_factor << " x "<< photon_Pt << endl;
// 		//fill the final historgam
// 		hfinal_MET->Fill(mvamet, R_factor);
// 	}
// 	cout << "Finished with met correction" << endl;
// 	//TCanvas *c = new TCanvas();
// 	//hfinal_MET->Draw();
// return hfinal_MET;
// }

TH1F* Apply_Pt_Correction(TString root_directory, TString category , const char* filename, TH1F* Z_pt_Photon_Pt_Ratio_hist, Float_t &photon_Pt)
{
	cout << "Applying corrections based on Photon Pt" << endl;
	// Start by looping over the Photon + Jet 
	TFile *inclusive = new TFile(root_directory);
	TTree *tree = new TTree();
	tree = (TTree*) inclusive->FindObjectAny( filename );
	Float_t mvamet;
	//Float_t photon_Pt;
	Int_t event; 
	// tree->SetBranchAddress("genVpt" , &photon_Pt);
	tree->SetBranchAddress("mvamet" , &mvamet);
	if (TString(filename) == "Photon_photon_control" )
	{
		tree->SetBranchAddress("genVpt" , &photon_Pt);
	}
	else if ( TString(filename) == "data_photon_control")
	{
		tree->SetBranchAddress("ptpho" , &photon_Pt);
	}
	else 
	{
		cout<< "THIS DOES NOT MAKE SENSE" << endl;
	}
	//create hist to store the corrected final met
	Int_t nbins ; Double_t xMin ; Double_t xMax ;
	if (category == "monojet")
	{
		nbins = 80 ; xMin = 200 ; xMax = 1000;
	}
	else
	{
		nbins = 75 ; xMin = 250 ; xMax = 1000;
	}
	TH1F *hfinal_MET = new TH1F( "mvamet","Met Estimation", nbins, xMin, xMax);
	hfinal_MET->Sumw2();
	//read all entries and fill the hist
	Long64_t nentries = (Int_t) tree->GetEntries();
	cout << "Starting the event by event correction of "<< filename << endl;
	for (Long64_t entry = 0; entry<nentries; entry++ )
	{
		tree->GetEntry(entry);
		//correct for difference in Z pt and Photon pt
		Double_t R_factor = R(photon_Pt, Z_pt_Photon_Pt_Ratio_hist);
		//fill the final historgam
		hfinal_MET->Fill(mvamet, R_factor);
	}
	cout << "Finished with MET Estimation from the Photon + Jet control region" << endl;
return hfinal_MET;
}



void Bin_by_Bin_comparison(TH1F* h1, TH1F* h2)
{
	for(Int_t bin_no =1 ; bin_no < 100+1 ; bin_no++) 
	{
		//cout << "The value of h1 at " << bin_no << " = " << h1->GetBinContent(bin_no) << " compared to " << h2->GetBinContent(bin_no) << endl;
	}
}

//Apply bin by bin correction from a ratio hist
TH1F* Bin_by_Bin_Correction(TH1F* h , TH1F* h_ratio)
{
//	TCanvas *cbb = new TCanvas();
	//check to make sure the h_ratio and h histograms match in x range and nbins
	cout << h_ratio->GetNbinsX() << " " << h_ratio->GetXaxis()->GetXmax() << " " << h_ratio->GetXaxis()->GetXmin() << endl;
	cout << h->GetNbinsX() << " " << h->GetXaxis()->GetXmax() << " " << h->GetXaxis()->GetXmin() << endl;
	TH1F* h_Corrected = new TH1F("h_new", "h_new", h_ratio->GetNbinsX(), h_ratio->GetXaxis()->GetXmin(), h_ratio->GetXaxis()->GetXmax());
	if ( !(h->GetXaxis()->GetXmin() == h_ratio->GetXaxis()->GetXmin()) || !(h->GetXaxis()->GetXmax() == h_ratio->GetXaxis()->GetXmax()) || !(h->GetNbinsX() == h_ratio->GetNbinsX()) )
	{
		cout << "THE RATIOS AND THE HISTOGRAM DO NOT MAKE SENSE FOR CORRECTION "<< endl;
	}
	else
	{
		//make new histogram to fill with the corrected values
		Int_t N_bins = h->GetNbinsX();
		Double_t h_value;
		Double_t h_ratio_value ;
		Double_t corrected_value;
		//loop over the 2 histograms and use the values to do corrections
		for ( Int_t bin_no = 1 ; bin_no < N_bins+1 ; bin_no++)
		{
			corrected_value = h->GetBinContent(bin_no) * h_ratio->GetBinContent(bin_no) ;
			//cout << " The value at " << bin_no << " = " << h->GetBinContent(bin_no) << ". Correcting with factor " << h_ratio->GetBinContent(bin_no) << " = " << corrected_value<< endl;
			h_Corrected->SetBinContent(bin_no, corrected_value);
		}
	}

//	h_Corrected->Draw();
return h_Corrected;
}

#endif
