
#ifndef BACKGROUNDS_H
#define BACKGROUNDS_H
#include "TCanvas.h"
#include "TAxis.h"
#include "TCut.h"
#include "TSystem.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TChain.h"
#include "TCanvas.h"
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
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooMinuit.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"

using namespace std;


// =============================================================================
// Get the Histograms for Backgounds in MET that will be used in the Likelihood
// =============================================================================
std::vector< std::vector<TH1F*> > Backgrounds_MET(TString file_directory , TString category)
{
	cout << "a" << endl;
	TString rootfile = file_directory + category + "-combo-pfmetraw-fj200.root" ;

	Double_t bin_edges_1[19] = {200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 350, 380, 430, 500, 1000};
	Double_t bin_edges_2[7] = { 250, 300, 350, 400, 450, 500, 1000};
	Int_t n_edges;
	if (category == "monojet"){ n_edges = sizeof(bin_edges_1) / sizeof(bin_edges_1[1]) ;} else {n_edges = sizeof(bin_edges_2) / sizeof(bin_edges_2[1]) ;}
	Double_t bin_edges[n_edges];
	if (category == "monojet"){	for (Int_t i = 0; i < 19 ; ++i)	{bin_edges[i] = bin_edges_1[i];	}}else {for (Int_t i = 0; i < 7; ++i){bin_edges[i] = bin_edges_2[i];}}

	Int_t nbins = 0;
	Double_t xMin = 0;
	Double_t xMax = 0;

	if (category == "monojet"){		nbins = 80; xMin = 200 ; xMax = 1000;	}	else	{		nbins = 75 ; xMin = 250 ; xMax = 1000;	}

	std::vector< std::vector<TH1F*> > v_All_Backgrounds;
	std::vector<TH1F*> v_Photon_Backgrounds;
	std::vector<TH1F*> v_Dimuon_Backgrounds;
	std::vector<TH1F*> v_Single_Muon_Backgrounds;

	// **** Photon+Jet Control Region ****
	// Calculated using purity
	TH1F* h_Photon_Photon_control_background = (TH1F*) Make_my_hist(rootfile, category, "Photon_photon_control", "mvamet", "background", nbins, xMin, xMax, xMin);

	TH1F* h_Photon_Photon_control_background_rebinned = (TH1F*) Rebin_with_density(h_Photon_Photon_control_background, bin_edges, n_edges) ;		v_Photon_Backgrounds.push_back(h_Photon_Photon_control_background_rebinned) ;
	cout << "d" << endl;
	
	// **** Di-muon Control Region ****
	// Top and Dibosons
	TH1F* h_WW_di_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "WW_di_muon_control", "mvamet", "background" , nbins , xMin, xMax, xMin) ; 
	TH1F* h_WW_di_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_WW_di_muon_control_weighted, bin_edges, n_edges) ;	v_Dimuon_Backgrounds.push_back(h_WW_di_muon_control_weighted_rebinned) ;
	TH1F* h_WZ_di_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "WZ_di_muon_control", "mvamet", "background" , nbins , xMin, xMax, xMin) ; 
	TH1F* h_WZ_di_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_WZ_di_muon_control_weighted, bin_edges, n_edges) ;	v_Dimuon_Backgrounds.push_back(h_WZ_di_muon_control_weighted_rebinned) ; 
	TH1F* h_ZZ_di_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "ZZ_di_muon_control", "mvamet", "background" , nbins , xMin, xMax, xMin) ; 
	TH1F* h_ZZ_di_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_ZZ_di_muon_control_weighted , bin_edges, n_edges) ;	v_Dimuon_Backgrounds.push_back(h_ZZ_di_muon_control_weighted_rebinned) ;
	TH1F* h_SingleTop_di_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "SingleTop_di_muon_control", "mvamet", "background" , nbins , xMin, xMax, xMin) ; 
	TH1F* h_SingleTop_di_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_SingleTop_di_muon_control_weighted, bin_edges, n_edges) ;	v_Dimuon_Backgrounds.push_back(h_SingleTop_di_muon_control_weighted_rebinned) ;
	TH1F* h_ttbar_di_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "ttbar_di_muon_control", "mvamet", "background" , nbins , xMin, xMax, xMin); 
	TH1F* h_ttbar_di_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_ttbar_di_muon_control_weighted, bin_edges, n_edges) ;	v_Dimuon_Backgrounds.push_back(h_ttbar_di_muon_control_weighted_rebinned) ;

	// **** Single Muon Control Region ***
	// Top, Dibosons, Z->ll , QCD
	TH1F* h_QCD_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "QCD_single_muon_control", "mvamet", "background" , nbins , xMin, xMax, xMin) ; 
	TH1F* h_QCD_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_QCD_single_muon_control_weighted, bin_edges, n_edges) ;		v_Single_Muon_Backgrounds.push_back(h_QCD_single_muon_control_weighted_rebinned) ;
	TH1F* h_Zll_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "Zll_single_muon_control", "mvamet", "background" , nbins , xMin, xMax, xMin) ; 
	TH1F* h_Zll_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_Zll_single_muon_control_weighted, bin_edges, n_edges) ;		v_Single_Muon_Backgrounds.push_back(h_Zll_single_muon_control_weighted_rebinned) ;
	TH1F* h_WW_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "WW_single_muon_control", "mvamet", "background" , nbins , xMin, xMax, xMin) ; 
	TH1F* h_WW_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_WW_single_muon_control_weighted, bin_edges, n_edges) ;		v_Single_Muon_Backgrounds.push_back(h_WW_single_muon_control_weighted_rebinned) ;
	TH1F* h_WZ_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "WZ_single_muon_control", "mvamet", "background" , nbins , xMin, xMax, xMin) ; 
	TH1F* h_WZ_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_WZ_single_muon_control_weighted, bin_edges, n_edges) ;	v_Single_Muon_Backgrounds.push_back(h_WW_single_muon_control_weighted_rebinned) ;
	TH1F* h_ZZ_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "ZZ_single_muon_control", "mvamet", "background" , nbins , xMin, xMax, xMin) ; 
	TH1F* h_ZZ_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_ZZ_single_muon_control_weighted, bin_edges, n_edges) ;		v_Single_Muon_Backgrounds.push_back(h_ZZ_single_muon_control_weighted_rebinned) ;
	TH1F* h_ttbar_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "ttbar_single_muon_control", "mvamet", "background" , nbins , xMin, xMax, xMin) ; 
	TH1F* h_ttbar_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_ttbar_single_muon_control_weighted, bin_edges, n_edges) ;		v_Single_Muon_Backgrounds.push_back(h_ttbar_single_muon_control_weighted_rebinned) ;
	TH1F* h_SingleTop_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "SingleTop_single_muon_control", "mvamet", "background" , nbins , xMin, xMax, xMin) ; 
	TH1F* h_SingleTop_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_SingleTop_single_muon_control_weighted, bin_edges, n_edges) ;	v_Single_Muon_Backgrounds.push_back(h_SingleTop_single_muon_control_weighted_rebinned) ;


	// Add the backgrounds from the different control regions to the main vector
	v_All_Backgrounds.push_back(v_Photon_Backgrounds) ;
	v_All_Backgrounds.push_back(v_Dimuon_Backgrounds) ;
	v_All_Backgrounds.push_back(v_Single_Muon_Backgrounds) ;

return v_All_Backgrounds;
}


// Plot these to compare to figure 37
// void Backgrounds_Single_Muon(TString file_directory , TString category , TString variable)
// {
// 	// I need to plot distributions for 3 variables: mass, lead muon pt, N jets
// 	TString variables[3] = { }

// 	// **** Single Muon Control Region ***
// 	// Top, Dibosons, Z->ll , QCD
// 	TH1F* h_QCD_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "QCD_single_muon_control", "mvamet", "weight" , nbins , xMin, xMax, xMin) ; 
// 	TH1F* h_QCD_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_QCD_single_muon_control_weighted, bin_edges, n_edges) ;		v_Single_Muon_Backgrounds.push_back(h_QCD_single_muon_control_weighted_rebinned) ;
// 	TH1F* h_Zll_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "Zll_single_muon_control", "mvamet", "weight" , nbins , xMin, xMax, xMin) ; 
// 	TH1F* h_Zll_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_Zll_single_muon_control_weighted, bin_edges, n_edges) ;		v_Single_Muon_Backgrounds.push_back(h_Zll_single_muon_control_weighted_rebinned) ;
// 	TH1F* h_WW_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "WW_single_muon_control", "mvamet", "weight" , nbins , xMin, xMax, xMin) ; 
// 	TH1F* h_WW_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_WW_single_muon_control_weighted, bin_edges, n_edges) ;		v_Single_Muon_Backgrounds.push_back(h_WW_single_muon_control_weighted_rebinned) ;
// 	TH1F* h_WZ_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "WZ_single_muon_control", "mvamet", "weight" , nbins , xMin, xMax, xMin) ; 
// 	TH1F* h_WZ_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_WZ_single_muon_control_weighted, bin_edges, n_edges) ;	v_Single_Muon_Backgrounds.push_back(h_WW_single_muon_control_weighted_rebinned) ;
// 	TH1F* h_ZZ_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "ZZ_single_muon_control", "mvamet", "weight" , nbins , xMin, xMax, xMin) ; 
// 	TH1F* h_ZZ_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_ZZ_single_muon_control_weighted, bin_edges, n_edges) ;		v_Single_Muon_Backgrounds.push_back(h_ZZ_single_muon_control_weighted_rebinned) ;
// 	TH1F* h_ttbar_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "ttbar_single_muon_control", "mvamet", "weight" , nbins , xMin, xMax, xMin) ; 
// 	TH1F* h_ttbar_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_ttbar_single_muon_control_weighted, bin_edges, n_edges) ;		v_Single_Muon_Backgrounds.push_back(h_ttbar_single_muon_control_weighted_rebinned) ;
// 	TH1F* h_SingleTop_single_muon_control_weighted = (TH1F*) Make_my_hist(rootfile, category, "SingleTop_single_muon_control", "mvamet", "weight" , nbins , xMin, xMax, xMin) ; 
// 	TH1F* h_SingleTop_single_muon_control_weighted_rebinned = (TH1F*) Rebin_with_density(h_SingleTop_single_muon_control_weighted, bin_edges, n_edges) ;	v_Single_Muon_Backgrounds.push_back(h_SingleTop_single_muon_control_weighted_rebinned) ;

// }

// Summ all the hists in a vector of hists 
TH1F* Sum_all_hists_in_vector(std::vector<TH1F*> v_hist)
{
	//declare the base hist that will hold the sum of hists
	TH1F* h_all = new TH1F() ;
	if (v_hist.size() == 1)
	{
		h_all->Add(v_hist.at(0) , 1);
	}
	else
	{
		for (int i = 0; i < v_hist.size() ; ++i)
		{
			h_all->Add(v_hist.at(i) , 1) ;
		}
	}
return h_all;
}

#endif