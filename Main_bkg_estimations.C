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


using namespace std;

void Get_Estimate_from(TString file_directory , TString category)
{

	TString rootfile = file_directory + category + "-combo-pfmetraw-fj200.root" ;

	Double_t bin_edges_1[19] = {200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 350, 380, 430, 500, 1000};
	Double_t bin_edges_2[7] = { 250, 300, 350, 400, 450, 500, 1000};
	Int_t n_edges;
	if (category == "monojet"){ n_edges = sizeof(bin_edges_1) / sizeof(bin_edges_1[1]) ;} else {n_edges = sizeof(bin_edges_2) / sizeof(bin_edges_2[1]) ;}
	Double_t bin_edges[n_edges];
	if (category == "monojet"){	for (Int_t i = 0; i < 19 ; ++i)	{bin_edges[i] = bin_edges_1[i];	}}else {for (Int_t i = 0; i < 7; ++i){bin_edges[i] = bin_edges_2[i];}}

	//make a root file to save all my histograms
	//TFile *all_hist = new TFile("/afs/cern.ch/work/b/bbachu/private/Z_nunu/Z_nunu_Estimate_hists","RECREATE");

	Int_t nbins = 0;
	Double_t xMin = 0;
	Double_t xMax = 0;

	if (category == "monojet")
	{
		nbins = 80; xMin = 200 ; xMax = 1000;
	}
	else
	{
		nbins = 75 ; xMin = 250 ; xMax = 1000;
	}


	// // in MC unweighted
	// TH1F* h_Znunu_Signal_genVpt_unweighted = (TH1F*) Make_my_hist(rootfile, category ,"Znunu_signal", "genVpt","weight", 100, 0, 1000, 170);

	// TH1F* h_Znunu_Signal_mvamet_unweighted = (TH1F*) Make_my_hist(rootfile, category , "Znunu_signal", "mvamet","weight" , 80 , 200, 1000 , 0);
	// TH1F* h_Photon_Photon_Control_mvamet_unweighted = (TH1F*) Make_my_hist(rootfile, category, "Photon_photon_control", "mvamet","weight" , 80 , 200, 1000, 0);
	// TH1F* h_Photon_Photon_Control_genVpt_unweighted = (TH1F*) Make_my_hist(rootfile, category ,"Photon_photon_control", "genVpt","weight" , 100, 0, 1000 , 170);
	// TH1F* h_Zll_di_muon_control_mvamet_unweighted = (TH1F*) Make_my_hist(rootfile, category , "Zll_di_muon_control" , "mvamet", "weight" , 80 , 200, 1000 , 0);
	
	// //weighted
 // 	TH1F* h_Photon_Photon_Control_genVpt_weighted = (TH1F*) Make_my_hist(rootfile, category ,"Photon_photon_control", "genVpt", "weight" , 100, 0, 1000, 170 );
 // 	TH1F* h_Photon_Photon_Control_mvamet_weighted = (TH1F*) Make_my_hist(rootfile, category , "Photon_photon_control", "mvamet", "weight", 100, 0, 1000, 0);
 // 	TH1F* h_Znunu_Signal_genVpt_weighted = (TH1F*) Make_my_hist(rootfile, category , "Znunu_signal", "genVpt", "weight" ,100, 0, 1000 , 170);
	// TH1F* h_Znunu_Signal_mvamet_weighted = (TH1F*) Make_my_hist(rootfile, category ,"Znunu_signal", "mvamet", "weight", 80 , 200, 1000 , 0);
	// TH1F* h_Zll_di_muon_control_mvamet_weighted = (TH1F*) Make_my_hist(rootfile, category , "Zll_di_muon_control" , "mvamet" , "weight" , 80 , 200, 1000 ,0 );

	// // in Data
	// TH1F* h_data_photon_control_ptpho = (TH1F*) Make_my_hist(rootfile, category , "data_photon_control", "ptpho","data" , 100, 0, 1000, 0);
	// TH1F* h_data_photon_control_mvamet = (TH1F*) Make_my_hist(rootfile , category ,"data_photon_control", "mvamet", "data", 100, 0, 1000, 0);
	// TH1F* h_data_di_muon_control_mvamet = (TH1F*) Make_my_hist(rootfile, category ,"data_di_muon_control", "mvamet", "data" , 80 , 200, 1000,0);

	// in MC unweighted
	TH1F* h_Znunu_Signal_genVpt_unweighted = (TH1F*) Make_my_hist(rootfile, category ,"Znunu_signal", "genVpt","weight", 100, 0, 1000, 170);

	TH1F* h_Znunu_Signal_mvamet_unweighted = (TH1F*) Make_my_hist(rootfile, category , "Znunu_signal", "mvamet","weight" , nbins , xMin, xMax , xMin);
	TH1F* h_Photon_Photon_Control_mvamet_unweighted = (TH1F*) Make_my_hist(rootfile, category, "Photon_photon_control", "mvamet","weight" , nbins , xMin, xMax, xMin);
	TH1F* h_Photon_Photon_Control_genVpt_unweighted = (TH1F*) Make_my_hist(rootfile, category ,"Photon_photon_control", "genVpt","weight" , 100, 0, 1000 , 170);
	TH1F* h_Zll_di_muon_control_mvamet_unweighted = (TH1F*) Make_my_hist(rootfile, category , "Zll_di_muon_control" , "mvamet", "weight" , nbins , xMin, xMax , xMin);
	TH1F* h_Wjets_signal_mvamet_unweighted = (TH1F*) Make_my_hist(rootfile, category, "Wjets_signal", "mvamet", "weight", nbins , xMin , xMax, xMin);
	TH1F* h_Wjets_single_control_mvamet_unweighted = (TH1F*) Make_my_hist(rootfile, category, "Wjets_single_muon_control", "mvamet","weight", nbins, xMin , xMax,xMin);

	//weighted
 	TH1F* h_Photon_Photon_Control_genVpt_weighted = (TH1F*) Make_my_hist(rootfile, category ,"Photon_photon_control", "genVpt", "weight" , 100, 0, 1000, 170 );
 	TH1F* h_Photon_Photon_Control_mvamet_weighted = (TH1F*) Make_my_hist(rootfile, category , "Photon_photon_control", "mvamet", "weight", nbins, xMin, xMax, xMin);
 	TH1F* h_Znunu_Signal_genVpt_weighted = (TH1F*) Make_my_hist(rootfile, category , "Znunu_signal", "genVpt", "weight" ,100, 0, 1000 , 170);
	TH1F* h_Znunu_Signal_mvamet_weighted = (TH1F*) Make_my_hist(rootfile, category ,"Znunu_signal", "mvamet", "weight", nbins , xMin, xMax , 0);
	TH1F* h_Zll_di_muon_control_mvamet_weighted = (TH1F*) Make_my_hist(rootfile, category , "Zll_di_muon_control" , "mvamet" , "weight" , nbins , xMin, xMax ,xMin );
	TH1F* h_Wjets_signal_mvamet_weighted = (TH1F*) Make_my_hist(rootfile, category, "Wjets_signal", "mvamet", "weight", nbins , xMin , xMax, xMin);
	TH1F* h_Wjets_single_control_mvamet_weighted = (TH1F*) Make_my_hist(rootfile, category, "Wjets_single_muon_control", "mvamet", "weight", nbins, xMin , xMax,xMin);
	// in Data
	TH1F* h_data_photon_control_ptpho = (TH1F*) Make_my_hist(rootfile, category , "data_photon_control", "ptpho","data" , 100, 0, 1000, 0);
	TH1F* h_data_photon_control_mvamet = (TH1F*) Make_my_hist(rootfile , category ,"data_photon_control", "mvamet", "data", nbins, xMin, xMax, xMin);
	TH1F* h_data_di_muon_control_mvamet = (TH1F*) Make_my_hist(rootfile, category ,"data_di_muon_control", "mvamet", "data" , nbins , xMin, xMax,xMin);
	TH1F* h_data_single_muon_control_mvamet = (TH1F*) Make_my_hist(rootfile , category , "data_single_muon_control", "mvamet", "data", nbins, xMin , xMax, xMin);
	TCanvas *Sdfa = new TCanvas(); h_data_single_muon_control_mvamet->Draw();

	//compare the MET distributions of both control regions to the signal region
	//Draw_2_Hist(h_Photon_Photon_Control_mvamet_unweighted,h_Znunu_Signal_mvamet_unweighted , "Comparision of MET in Z#rightarrow#nu#nu Signal Region to #gamma + Jet Control Region " + category, "MET", "","#gamma + Jet MC" ,"Z#rightarrow#nu#nu MC" , 2, 1, "");
	
	//Draw_2_Hist(h_Zll_di_muon_control_mvamet_unweighted , h_Znunu_Signal_mvamet_unweighted, "Comparision of MET in Z#rightarrow#nu#nu Signal Region to Z#rightarrow#mu#mu Control Region " + category, "MET", "", "Z#rightarrow#mu#mu MC", "Z#rightarrow#mu#mu", 3, 1, "");
	
	// Match_all_integrals_of_hist(h_Znunu_Signal_mvamet_unweighted, h_Photon_Photon_Control_mvamet_unweighted, h_Zll_di_muon_control_mvamet_unweighted , "Comparison of MET with equal integrals " + category, "MET", "", "Z#rightarrow#nu#nu", "#gamma + Jet", "Z#rightarrow#mu#mu", 1, 2, 4, "" , "", "");

	TH1F* h_Znunu_Signal_mvamet_unweighted_rebinned = (TH1F*) Rebin_with_density(h_Znunu_Signal_mvamet_unweighted , bin_edges , n_edges);
	TH1F* h_Photon_Photon_Control_mvamet_unweighted_rebinned = (TH1F*) Rebin_with_density(h_Photon_Photon_Control_mvamet_unweighted , bin_edges , n_edges);
	TH1F* h_Zll_di_muon_control_mvamet_unweighted_rebinned = (TH1F*) Rebin_with_density(h_Zll_di_muon_control_mvamet_unweighted , bin_edges , n_edges);
	//Match_all_integrals_of_hist(h_Znunu_Signal_mvamet_unweighted_rebinned, h_Photon_Photon_Control_mvamet_unweighted_rebinned, h_Zll_di_muon_control_mvamet_unweighted_rebinned , "Comparison of MET with equal integrals " + category, "MET", "", "Z#rightarrow#nu#nu", "#gamma + Jet", "Z#rightarrow#mu#mu", 1, 2, 4, "" , "", "");

	//REBIN all the met distributions for comparison in MonoV paper pg 38
	TH1F* h_Znunu_Signal_mvamet_weighted_rebinned =(TH1F*) Rebin_with_density(h_Znunu_Signal_mvamet_weighted , bin_edges , n_edges);
	TH1F* h_data_photon_control_mvamet_rebinned = (TH1F*) Rebin_with_density(h_data_photon_control_mvamet, bin_edges , n_edges );
	TH1F* h_Photon_Photon_Control_mvamet_weighted_rebinned = (TH1F*) Rebin_with_density(h_Photon_Photon_Control_mvamet_weighted , bin_edges, n_edges);
	TH1F* h_data_di_muon_control_mvamet_rebinned = (TH1F*) Rebin_with_density(h_data_di_muon_control_mvamet , bin_edges , n_edges);
	TH1F* h_Zll_di_muon_control_mvamet_weighted_rebinned = (TH1F*) Rebin_with_density(h_Zll_di_muon_control_mvamet_weighted, bin_edges, n_edges);
	TH1F* h_Wjets_signal_mvamet_weighted_rebinned = (TH1F*) Rebin_with_density(h_Wjets_signal_mvamet_weighted, bin_edges, n_edges);
	TH1F* h_Wjets_single_control_mvamet_weighted_rebinned = (TH1F*) Rebin_with_density(h_Wjets_single_control_mvamet_weighted , bin_edges , n_edges);
	TH1F* h_data_single_muon_control_mvamet_rebinned = (TH1F*) Rebin_with_density(h_data_single_muon_control_mvamet , bin_edges , n_edges);

	//Draw_2_Hist(h_Photon_Photon_Control_mvamet_weighted_rebinned, h_data_photon_control_mvamet_rebinned, "Comparison of MET in #gamma Jet Data and MC", "mvamet" , "", "MC", "Data", 2, 1, "");
	// Draw_2_Hist(h_Zll_di_muon_control_mvamet_weighted_rebinned, h_data_di_muon_control_mvamet_rebinned, "Comparison of MET in Z#rightarrow#mu#my  Data and MC", "mvamet" , "", "MC", "Data", 2, 1, "");
	Draw_CMS_Preliminary_2h(category, bin_edges , n_edges ,  h_Photon_Photon_Control_mvamet_weighted_rebinned , h_data_photon_control_mvamet_rebinned,  "Expected (pre-fit)", "Data: #gamma + Jet", 4 , 1 , 1 , 20 , "Fake MET/GeV");
	Draw_CMS_Preliminary_2h(category, bin_edges , n_edges ,  h_Zll_di_muon_control_mvamet_weighted_rebinned , h_data_di_muon_control_mvamet_rebinned,  "Expected (pre-fit)" , "Data: Z#rightarrow#mu#mu", 4 , 1 , 1 , 20 , "Fake MET/GeV");
	Draw_CMS_Preliminary_2h(category, bin_edges , n_edges ,  h_Wjets_single_control_mvamet_weighted_rebinned , h_data_single_muon_control_mvamet_rebinned,  "Expected (pre-fit)" , "Data: W#rightarrow l#nu", 4 , 1 , 1 , 20 , "Fake MET/GeV");

// ===============================================
// ESTIMATE USING THE PHOTON + JET control region
// ===============================================

	// the first step is to construt the histogram of the ratios of the Z_pt to the photon_pt
 	TH1F* h_pho_transfer_function_unweighted = (TH1F*) h_Ratio(h_Znunu_Signal_genVpt_unweighted, h_Photon_Photon_Control_genVpt_unweighted);// TCanvas *e = new TCanvas(); h_pho_transfer_function_unweighted->Draw();
	Float_t genVpt ;
	Float_t ptpho;
	//plot the transfer function rebinned
	TH1F* h_Znunu_Signal_genVpt_unweighted_rebinned = (TH1F*) Rebin_with_density( h_Znunu_Signal_genVpt_unweighted , bin_edges , n_edges);
	TH1F* h_Photon_Photon_Control_genVpt_unweighted_rebinned = (TH1F*) Rebin_with_density( h_Photon_Photon_Control_genVpt_unweighted , bin_edges , n_edges);
	TH1F* h_transfer_Gamma_Jet_MC  = (TH1F*) Divide_with_variable_bins(h_Znunu_Signal_genVpt_unweighted_rebinned , h_Photon_Photon_Control_genVpt_unweighted_rebinned , category,  bin_edges ,  n_edges) ;
	//Apply this transfer function to the Photon+Jet control region in MC
	TH1F* h_Estimate_from_Pho_MC = (TH1F*) Apply_Pt_Correction(rootfile, category , "Photon_photon_control", h_pho_transfer_function_unweighted, genVpt);
	// TCanvas *c4 = new TCanvas();
	// h_Estimate_from_Pho_MC->Draw();
	// Draw_1_Hist( h_transfer_Gamma_Jet_MC , "Shape + Normalization Transfer Ratios (MC) " , "R^{#gamma}" ,"pT", 2 , "EP" , category , 0 , 1);
	//compare the estimate to the MC
	// Draw_2_Hist(h_Estimate_from_Pho_MC, h_Znunu_Signal_mvamet_unweighted, "Estimate of Z#rightarrow#nu#nu from MC " + category, "MET", "", "Estimate from #gamma + Jet MC", "Z#rightarrow#nu#nu MC", 2, 1, "");
	//now that this should work with MC, we can apply it to data
	//make the histogram that acts as the transfer function for the data
	TH1F* h_pho_transfer_function_weighted = (TH1F*) h_Ratio(h_Znunu_Signal_genVpt_weighted, h_Photon_Photon_Control_genVpt_weighted); //TCanvas *cnew = new TCanvas(); h_pho_transfer_function_weighted->Draw();
	//plot this transfer function for compasion in paper page 34
	//rebin the individual hists first
	TH1F* h_Znunu_Signal_genVpt_weighted_rebinned = (TH1F*) Rebin_with_density(h_Znunu_Signal_genVpt_weighted, bin_edges, n_edges); //h_Znunu_Signal_genVpt_weighted_rebinned->Draw();
	TH1F* h_Photon_Photon_Control_genVpt_weighted_rebinned = (TH1F*) Rebin_with_density(h_Photon_Photon_Control_genVpt_weighted, bin_edges, n_edges ); //h_Photon_Photon_Control_genVpt_weighted_rebinned->Draw();
	TH1F* h_transfer_Gamma_Jet = (TH1F*) Divide_with_variable_bins(h_Znunu_Signal_genVpt_weighted_rebinned , h_Photon_Photon_Control_genVpt_weighted_rebinned , category , bin_edges , n_edges);
	// TCanvas *transfer_data = new TCanvas(); h_transfer_Gamma_Jet->Draw();
	Draw_1_Hist( h_transfer_Gamma_Jet , "Shape + Normalization Transfer Ratios " , "R^{#gamma}" ,"pT", 2 , "EP" , category , 0 , 1);
	//apply this transfer function to the data
	TH1F* h_Estimate_from_Pho_Data = (TH1F*) Apply_Pt_Correction(rootfile, category, "data_photon_control", h_pho_transfer_function_weighted, ptpho);
	// Draw_2_Hist(h_Estimate_from_Pho_Data, h_Znunu_Signal_mvamet_weighted, "Estimate of Z#rightarrow#nu#nu from Data" + category, "MET", "", "Estimate from #gamma + Jet Data", "Z#rightarrow#nu#nu", 2 , 1, "Data");
	


//=================================================
// ESTIMATE USING THE DI-MUON cntrol region
// ===============================================
	// the first step is to construct the histogram of the ratios of the Z_nunu MC to Z_mumu MC
	TH1F* h_di_muon_transfer_function_unweighted = (TH1F*) h_Ratio(h_Znunu_Signal_mvamet_unweighted , h_Zll_di_muon_control_mvamet_unweighted);
	//now use this transfer function to correct the shape
	TH1F* h_Estimate_from_di_muon_MC = (TH1F*) Bin_by_Bin_Correction(h_Zll_di_muon_control_mvamet_unweighted , h_di_muon_transfer_function_unweighted);
	// Draw_2_Hist(h_Estimate_from_di_muon_MC, h_Znunu_Signal_mvamet_unweighted, "Estimate of Z#rightarrow#nu#nu from Di-muon MC" + category, "MET", "", "Estimate from Z#rightarrow#mu#mu MC", "Z#rightarrow#nu#nu ", 4, 1, "");
	//make the histogram so that we can do the estimate from data
	TH1F* h_di_muon_transfer_function_weighted = (TH1F*) h_Ratio( h_Znunu_Signal_mvamet_weighted , h_Zll_di_muon_control_mvamet_weighted );
	//plot this histogram to compare to MonoV paper , first rebin and then plot
	TH1F* h_transfer_rebinned_Zmumu = (TH1F*) Divide_with_variable_bins(h_Znunu_Signal_mvamet_weighted_rebinned , h_Zll_di_muon_control_mvamet_weighted_rebinned, category ,bin_edges , n_edges);
	Draw_1_Hist(h_transfer_rebinned_Zmumu, "Shape+Normalization Transfer Ratios " , "R^{Z}" ,"MET(GeV)", 4 , "EP" , category , 0, 14);
	//create the estimate from data
	TH1F* h_Estimate_from_di_muon_Data = (TH1F*) Bin_by_Bin_Correction(h_data_di_muon_control_mvamet , h_di_muon_transfer_function_weighted);
	// Draw_2_Hist(h_Estimate_from_di_muon_Data, h_Znunu_Signal_mvamet_weighted, "Estimate of Z#rightarrow#nu#nu from Dimuon Data" + category, "MET", "", "Estiamte from Z#rightarrow#mu#mu Data", "Z#rightarrow#nu#nu", 4, 1, "Data");
	

//=====================================================
// ESTIMATION FOR W+JETS signal region
//=====================================================

	// the first step is to construct the histogram of the ratios of the W_lnu MC to w_lnu MC
	TH1F* h_wjets_transfer_function_unweighted = (TH1F*) h_Ratio(h_Wjets_signal_mvamet_unweighted , h_Wjets_single_control_mvamet_unweighted );
	//now use this transfer function to correct the shape
	TH1F* h_Estimate_from_wjets_MC = (TH1F*) Bin_by_Bin_Correction(h_Wjets_single_control_mvamet_unweighted , h_wjets_transfer_function_unweighted);
	// Draw_2_Hist(h_Estimate_from_wjets_MC, h_Wjets_signal_mvamet_unweighted, "Estimate of W#rightarrow l#nu from W#rightarrow l#nu MC" + category, "MET", "", "Estimate from W#rightarrow l#nu MC", "W#rightarrow l#nu", 4, 1, "");
	
	//make the histogram so that we can do the estimate from data
	TH1F* h_wjets_transfer_function_weighted = (TH1F*) h_Ratio( h_Wjets_signal_mvamet_weighted , h_Wjets_single_control_mvamet_weighted );
	// plot this histogram to compare to MonoV paper , first rebin and then plot
	TH1F* h_wjets_transfer_rebinned = (TH1F*) Divide_with_variable_bins(h_Wjets_signal_mvamet_weighted_rebinned , h_Wjets_single_control_mvamet_weighted_rebinned, category ,bin_edges , n_edges);
	Draw_1_Hist(h_wjets_transfer_rebinned, "Shape+Normalization Transfer Ratios " , "R^{W}" ,"MET(GeV)", 4 , "EP" , category , 0, 4);
// 	//create the estimate from data
	TH1F* h_Estimate_from_single_muon_data = (TH1F*) Bin_by_Bin_Correction(h_data_single_muon_control_mvamet , h_wjets_transfer_function_weighted);
	Draw_2_Hist(h_Estimate_from_single_muon_data, h_Wjets_signal_mvamet_weighted, "Estimate of W#rightarrow l#nu from W#rightarrow l#nu Data" + category, "MET", "", "Estiamte from W#rightarrow l#nu Data", "W#rightarrow l#nu", 4, 1, "Data");
	


	// //REBIN ALL THE DATA ESTIMATES
	TH1F* h_Estimate_from_Pho_Data_rebinned = (TH1F*) Rebin_with_density(h_Estimate_from_Pho_Data , bin_edges , n_edges);
	TH1F* h_Estimate_from_di_muon_Data_rebinned = (TH1F*) Rebin_with_density(h_Estimate_from_di_muon_Data , bin_edges , n_edges);
	TH1F* h_Estimate_from_single_muon_data_rebinned = (TH1F*) Rebin_with_density(h_Estimate_from_single_muon_data , bin_edges, n_edges);
	// //Draw_3_Hist( h_Znunu_Signal_mvamet_weighted_rebinned , h_Estimate_from_Pho_Data_rebinned,h_Estimate_from_di_muon_Data_rebinned,  "Z#rightarrow#nu#nu MC", "Estimate from #gamma + Jet Data", "Estimate from Z#rightarrow#mu#mu Data" ,  1 , 2 , 4 );
	// //TCanvas* c5 = new TCanvas();
	Draw_CMS_Preliminary_2h(category ,bin_edges, n_edges,  h_Znunu_Signal_mvamet_weighted_rebinned , h_Estimate_from_Pho_Data_rebinned, "Z#rightarrow#nu#nu MC" , "Estiamte from #gamma + Jet Data", 1 ,  2 ,1, 20 , "MET/GeV");
	Draw_CMS_Preliminary_2h(category, bin_edges, n_edges , h_Znunu_Signal_mvamet_weighted_rebinned , h_Estimate_from_di_muon_Data_rebinned,  "Z#rightarrow#nu#nu MC" , "Estimate from Z#rightarrow#mu#mu Data" , 1 , 4, 1, 22 , "MET/GeV");
	Draw_CMS_Preliminary_2h(category, bin_edges, n_edges , h_Wjets_signal_mvamet_weighted_rebinned , h_Estimate_from_single_muon_data_rebinned,  "W#rightarrow l#nu MC" , "Estimate from W#rightarrow l#nu Data" , 1 , 6, 1, 34 , "MET/GeV");
	
	Draw_CMS_Preliminary_3h(category, bin_edges, n_edges,  h_Znunu_Signal_mvamet_weighted_rebinned , h_Estimate_from_Pho_Data_rebinned,h_Estimate_from_di_muon_Data_rebinned,  "Z#rightarrow#nu#nu MC", "Estimate from #gamma + Jet Data", "Estimate from Z#rightarrow#mu#mu Data" ,  1 , 2 , 4 , 20 , 22 , "MET/GeV");


}


void Main_bkg_estimations()
{
	Get_Estimate_from("/afs/cern.ch/work/b/bbachu/private/Z_nunu/" , "monojet");
	Get_Estimate_from("/afs/cern.ch/work/b/bbachu/private/Z_nunu/" , "boosted");
	Get_Estimate_from("/afs/cern.ch/work/b/bbachu/private/Z_nunu/" , "resolved");

}


