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

void Get_Estimate_from(TString file_directory , TString filename)
{
	TString rootfile = file_directory + filename + "-combo-pfmetraw-fj200.root" ;

	Double_t bin_edges_1[19] = {200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 350, 380, 430, 500, 1000};
	Double_t bin_edges_2[9] = {200, 250, 300, 350, 400, 450, 500, 550, 1000};
	Int_t n_edges;
	if (filename == "monojet"){ n_edges = sizeof(bin_edges_1) / sizeof(bin_edges_1[1]) ;} else {n_edges = sizeof(bin_edges_2) / sizeof(bin_edges_2[1]) ;}
	Double_t bin_edges[n_edges];
	if (filename == "monojet"){	for (Int_t i = 0; i < 19 ; ++i)	{bin_edges[i] = bin_edges_1[i];	}}else {for (Int_t i = 0; i < 9; ++i){bin_edges[i] = bin_edges_2[i];}}

	//make a root file to save all my histograms
	//TFile *all_hist = new TFile("/afs/cern.ch/work/b/bbachu/private/Z_nunu/Z_nunu_Estimate_hists","RECREATE");

	// in MC unweighted
	TH1F* h_Znunu_Signal_genVpt_unweighted = (TH1F*) Make_my_hist(rootfile,"Znunu_signal", "genVpt","weight", 100, 0, 1000, 170);
	TH1F* h_Znunu_Signal_mvamet_unweighted = (TH1F*) Make_my_hist(rootfile, "Znunu_signal", "mvamet","weight" , 80 , 200, 1000 , 0);
	TH1F* h_Photon_Photon_Control_mvamet_unweighted = (TH1F*) Make_my_hist(rootfile, "Photon_photon_control", "mvamet","weight" , 80 , 200, 1000, 0);
	TH1F* h_Photon_Photon_Control_genVpt_unweighted = (TH1F*) Make_my_hist(rootfile, "Photon_photon_control", "genVpt","weight" , 100, 0, 1000 , 170);
	TH1F* h_Zll_di_muon_control_mvamet_unweighted = (TH1F*) Make_my_hist(rootfile, "Zll_di_muon_control" , "mvamet", "weight" , 80 , 200, 1000 , 0);
	
	//weighted
 	TH1F* h_Photon_Photon_Control_genVpt_weighted = (TH1F*) Make_my_hist(rootfile, "Photon_photon_control", "genVpt", "weight" , 100, 0, 1000, 170 );
 	TH1F* h_Photon_Photon_Control_mvamet_weighted = (TH1F*) Make_my_hist(rootfile, "Photon_photon_control", "mvamet", "weight", 100, 0, 1000, 0);
 	TH1F* h_Znunu_Signal_genVpt_weighted = (TH1F*) Make_my_hist(rootfile, "Znunu_signal", "genVpt", "weight" ,100, 0, 1000 , 170);
	TH1F* h_Znunu_Signal_mvamet_weighted = (TH1F*) Make_my_hist(rootfile, "Znunu_signal", "mvamet", "weight", 80 , 200, 1000 , 0);
	TH1F* h_Zll_di_muon_control_mvamet_weighted = (TH1F*) Make_my_hist(rootfile, "Zll_di_muon_control" , "mvamet" , "weight" , 80 , 200, 1000 ,0 );

	// in Data
	TH1F* h_data_photon_control_ptpho = (TH1F*) Make_my_hist(rootfile, "data_photon_control", "ptpho","data" , 100, 0, 1000, 0);
	TH1F* h_data_photon_control_mvamet = (TH1F*) Make_my_hist(rootfile , "data_photon_control", "mvamet", "data", 100, 0, 1000, 0);
	TH1F* h_data_di_muon_control_mvamet = (TH1F*) Make_my_hist(rootfile, "data_di_muon_control", "mvamet", "data" , 80 , 200, 1000,0);


	//compare the MET distributions of both control regions to the signal region
	Draw_2_Hist(h_Photon_Photon_Control_mvamet_unweighted,h_Znunu_Signal_mvamet_unweighted , "Comparision of MET in Z#rightarrow#nu#nu Signal Region to #gamma + Jet Control Region " + filename, "MET", "","#gamma + Jet MC" ,"Z#rightarrow#nu#nu MC" , 2, 1, "");
	
	Draw_2_Hist(h_Zll_di_muon_control_mvamet_unweighted , h_Znunu_Signal_mvamet_unweighted, "Comparision of MET in Z#rightarrow#nu#nu Signal Region to Z#rightarrow#mu#mu Control Region " + filename, "MET", "", "Z#rightarrow#mu#mu MC", "Z#rightarrow#mu#mu", 3, 1, "");
	
	Match_all_integrals_of_hist(h_Znunu_Signal_mvamet_unweighted, h_Photon_Photon_Control_mvamet_unweighted, h_Zll_di_muon_control_mvamet_unweighted , "Comparison of MET with equal integrals " + filename, "MET", "", "Z#rightarrow#nu#nu", "#gamma + Jet", "Z#rightarrow#mu#mu", 1, 2, 4, "" , "", "");

	TH1F* h_Znunu_Signal_mvamet_unweighted_rebinned = (TH1F*) Rebin_with_density(h_Znunu_Signal_mvamet_unweighted , bin_edges , n_edges);
	TH1F* h_Photon_Photon_Control_mvamet_unweighted_rebinned = (TH1F*) Rebin_with_density(h_Photon_Photon_Control_mvamet_unweighted , bin_edges , n_edges);
	TH1F* h_Zll_di_muon_control_mvamet_unweighted_rebinned = (TH1F*) Rebin_with_density(h_Zll_di_muon_control_mvamet_unweighted , bin_edges , n_edges);
	//Match_all_integrals_of_hist(h_Znunu_Signal_mvamet_unweighted_rebinned, h_Photon_Photon_Control_mvamet_unweighted_rebinned, h_Zll_di_muon_control_mvamet_unweighted_rebinned , "Comparison of MET with equal integrals " + filename, "MET", "", "Z#rightarrow#nu#nu", "#gamma + Jet", "Z#rightarrow#mu#mu", 1, 2, 4, "" , "", "");
	//TCanvas *cc = new TCanvas();
	//Draw_CMS_Preliminary(momentum_bins_1 ,19, h_Znunu_Signal_mvamet_unweighted_rebinned , h_Photon_Photon_Control_mvamet_unweighted_rebinned, h_Zll_di_muon_control_mvamet_unweighted_rebinned , "Z#rightarrow#nu#nu Signal MC" ,  "Estimate from #gamma + Jet MC", "Estimate from Z#rightarrow#mu#mu MC" , 1, 2, 4);

	//REBIN all the met distributions for comparison in MonoV paper pg 38
	TH1F* h_Znunu_Signal_mvamet_weighted_rebinned =(TH1F*) Rebin_with_density(h_Znunu_Signal_mvamet_weighted , bin_edges , n_edges);
	TH1F* h_data_photon_control_mvamet_rebinned = (TH1F*) Rebin_with_density(h_data_photon_control_mvamet, bin_edges , n_edges );
	TH1F* h_Photon_Photon_Control_mvamet_weighted_rebinned = (TH1F*) Rebin_with_density(h_Photon_Photon_Control_mvamet_weighted , bin_edges, n_edges);
	TH1F* h_data_di_muon_control_mvamet_rebinned = (TH1F*) Rebin_with_density(h_data_di_muon_control_mvamet , bin_edges , n_edges);
	TH1F* h_Zll_di_muon_control_mvamet_weighted_rebinned = (TH1F*) Rebin_with_density(h_Zll_di_muon_control_mvamet_weighted, bin_edges, n_edges);
	//Draw_CMS_Preliminary_2h(filename, bin_edges , n_edges ,  h_Photon_Photon_Control_mvamet_weighted_rebinned , h_data_photon_control_mvamet_rebinned,  "Expected (pre-fit)", "Data: #gamma + Jet", 4 , 1 ,20);
	//Draw_CMS_Preliminary_2h(filename, bin_edges , n_edges ,  h_Zll_di_muon_control_mvamet_weighted_rebinned , h_data_di_muon_control_mvamet_rebinned,  "Expected (pre-fit)" , "Data: Z#rightarrow#mu#mu", 4 , 1 ,20);

	// ===============================================
	// ESTIMATE USING THE PHOTON + JET control region

	// the first step is to construt the histogram of the ratios of the Z_pt to the photon_pt
 	TH1F* h_pho_transfer_function_unweighted = (TH1F*) h_Ratio(h_Znunu_Signal_genVpt_unweighted, h_Photon_Photon_Control_genVpt_unweighted); TCanvas *e = new TCanvas(); h_pho_transfer_function_unweighted->Draw();
	Float_t genVpt ;
	Float_t ptpho;
	//Apply this transfer function to the Photon+Jet control region in MC
	TH1F* h_Estimate_from_Pho_MC = (TH1F*) Apply_Pt_Correction(rootfile, "Photon_photon_control", h_pho_transfer_function_unweighted, genVpt);
	// TCanvas *c4 = new TCanvas();
	// h_Estimate_from_Pho_MC->Draw();
	//compare the estimate to the MC
	Draw_2_Hist(h_Estimate_from_Pho_MC, h_Znunu_Signal_mvamet_unweighted, "Estimate of Z#rightarrow#nu#nu from MC " + filename, "MET", "", "Estimate from #gamma + Jet MC", "Z#rightarrow#nu#nu MC", 2, 1, "");
	//now that this should work with MC, we can apply it to data
	//make the histogram that acts as the transfer function for the data
	TH1F* h_pho_transfer_function_weighted = (TH1F*) h_Ratio(h_Znunu_Signal_genVpt_weighted, h_Photon_Photon_Control_genVpt_weighted);TCanvas *cnew = new TCanvas();h_pho_transfer_function_weighted->Draw();
	//plot this transfer function for compasion in paper page 34
	//rebin the individual hists first
	TH1F* h_Znunu_Signal_genVpt_weighted_rebinned = (TH1F*) Rebin_with_density(h_Znunu_Signal_genVpt_weighted, bin_edges, n_edges); h_Znunu_Signal_genVpt_weighted_rebinned->Draw();
	TH1F* h_Photon_Photon_Control_genVpt_weighted_rebinned = (TH1F*) Rebin_with_density(h_Photon_Photon_Control_genVpt_weighted, bin_edges, n_edges ); h_Photon_Photon_Control_genVpt_weighted_rebinned->Draw();
	TH1F* h_transfer_Gamma_Jet = (TH1F*) Divide_with_variable_bins(h_Znunu_Signal_genVpt_weighted_rebinned , h_Photon_Photon_Control_genVpt_weighted_rebinned , bin_edges , n_edges);
	Draw_1_Hist(h_transfer_Gamma_Jet, "Shape+Normalization" , "R^{#gamma}" ,"MET(GeV)", 4 , "EP" , filename , 0, 1);
	//apply this transfer function to the data
	TH1F* h_Estimate_from_Pho_Data = (TH1F*) Apply_Pt_Correction(rootfile, "data_photon_control", h_pho_transfer_function_weighted, ptpho);
	Draw_2_Hist(h_Estimate_from_Pho_Data, h_Znunu_Signal_mvamet_weighted, "Estimate of Z#rightarrow#nu#nu from Data" + filename, "MET", "", "Estimate from #gamma + Jet Data", "Z#rightarrow#nu#nu", 2 , 1, "Data");
	

//==========================================
	
	//figue out what is going wrong with the photon pt 

	//plot the photon pt and Z pt distributions in MC
	// Draw_2_Hist(h_Photon_Photon_Control_genVpt_unweighted, h_Znunu_Signal_genVpt_weighted, "Comparison of the Z and #gamma Pt in " + filename, "Pt", "", "#gamma" , "Z", 2, 1, "");
	// Draw_2_Hist(h_PhotonPt_scaled_to_ZPt, h_Znunu_Signal_genVpt_unweighted, "#gamma pt scaled to Z pt in " + filename, "Pt", "", "#gamma", "Z", 2, 1, "");
	// TCanvas* a = new TCanvas();
	// h_ratio_equal_integrals->Draw();

	// TH1F* h_ZPt_scaled_to_PhotonPt = (TH1F*) Scale_h1_to_h2(h_Znunu_Signal_genVpt_unweighted, h_Photon_Photon_Control_genVpt_unweighted);
	// Draw_2_Hist(h_Photon_Photon_Control_genVpt_unweighted, h_ZPt_scaled_to_PhotonPt, "Z pt scaled to #gamma pt in " + filename , "Pt",  "", "#gamma" , "Z" , 2 , 1 , "");
	// TH1F* h_ratio_equal_integrals_2 = (TH1F*) h_Ratio(h_ZPt_scaled_to_PhotonPt , h_Photon_Photon_Control_genVpt_unweighted);
	// //Draw_2_Hist(TH1F* h_1, TH1F* h_2, const Char_t *Name, const Char_t *XAxis, const Char_t *YAxis, const Char_t *label1, const Char_t * label2, Color_t lcolor1, Color_t lcolor2, const Char_t* h1source)
	// TCanvas* b = new TCanvas();
	// h_ratio_equal_integrals_2->Draw();
	

// //=================================================


	// ===============================================
	// ESTIMATE USING THE DI-MUON cntrol region
	// the first step is to construct the histogram of the ratios of the Z_nunu MC to Z_mumu MC
	TH1F* h_di_muon_transfer_function_unweighted = (TH1F*) h_Ratio(h_Znunu_Signal_mvamet_unweighted , h_Zll_di_muon_control_mvamet_unweighted);
	//now use this transfer function to correct the shape
	TH1F* h_Estimate_from_di_muon_MC = (TH1F*) Bin_by_Bin_Correction(h_Zll_di_muon_control_mvamet_unweighted , h_di_muon_transfer_function_unweighted);
	Draw_2_Hist(h_Estimate_from_di_muon_MC, h_Znunu_Signal_mvamet_unweighted, "Estimate of Z#rightarrow#nu#nu from Di-muon MC" + filename, "MET", "", "Estimate from Z#rightarrow#mu#mu MC", "Z#rightarrow#nu#nu ", 4, 1, "");
	//make the histogram so that we can do the estimate from data
	TH1F* h_di_muon_transfer_function_weighted = (TH1F*) h_Ratio( h_Znunu_Signal_mvamet_weighted , h_Zll_di_muon_control_mvamet_weighted );
	//plot this histogram to compare to MonoV paper , first rebin and then plot
	TH1F* h_transfer_rebinned_Zmumu = (TH1F*) Divide_with_variable_bins(h_Znunu_Signal_mvamet_weighted_rebinned , h_Zll_di_muon_control_mvamet_weighted_rebinned ,bin_edges , n_edges);
	Draw_1_Hist(h_transfer_rebinned_Zmumu, "Shape+Normalization" , "R^{Z}" ,"MET(GeV)", 4 , "EP" , filename , 0, 14);
	//create the estimate from data
	TH1F* h_Estimate_from_di_muon_Data = (TH1F*) Bin_by_Bin_Correction(h_data_di_muon_control_mvamet , h_di_muon_transfer_function_weighted);
	Draw_2_Hist(h_Estimate_from_di_muon_Data, h_Znunu_Signal_mvamet_weighted, "Estimate of Z#rightarrow#nu#nu from Dimuon Data" + filename, "MET", "", "Estiamte from Z#rightarrow#mu#mu Data", "Z#rightarrow#nu#nu", 4, 1, "Data");


	

	//REBIN ALL THE DATA ESTIMATES
	TH1F* h_Estimate_from_Pho_Data_rebinned = (TH1F*) Rebin_with_density(h_Estimate_from_Pho_Data , bin_edges , n_edges);
	TH1F* h_Estimate_from_di_muon_Data_rebinned = (TH1F*) Rebin_with_density(h_Estimate_from_di_muon_Data , bin_edges , n_edges);
	//Draw_3_Hist( h_Znunu_Signal_mvamet_weighted_rebinned , h_Estimate_from_Pho_Data_rebinned,h_Estimate_from_di_muon_Data_rebinned,  "Z#rightarrow#nu#nu MC", "Estimate from #gamma + Jet Data", "Estimate from Z#rightarrow#mu#mu Data" ,  1 , 2 , 4 );
	TCanvas* c5 = new TCanvas();
	Draw_CMS_Preliminary_2h(filename ,bin_edges, n_edges,  h_Znunu_Signal_mvamet_weighted_rebinned , h_Estimate_from_Pho_Data_rebinned, "Z#rightarrow#nu#nu MC" , "Estiamte from #gamma + Jet Data", 1 ,  2 , 20);
	Draw_CMS_Preliminary_2h(filename, bin_edges, n_edges , h_Znunu_Signal_mvamet_weighted_rebinned , h_Estimate_from_di_muon_Data_rebinned,  "Z#rightarrow#nu#nu MC" , "Estimate from Z#rightarrow#mu#mu Data" , 1 , 4, 22);
	Draw_CMS_Preliminary_3h(filename, bin_edges, n_edges,  h_Znunu_Signal_mvamet_weighted_rebinned , h_Estimate_from_Pho_Data_rebinned,h_Estimate_from_di_muon_Data_rebinned,  "Z#rightarrow#nu#nu MC", "Estimate from #gamma + Jet Data", "Estimate from Z#rightarrow#mu#mu Data" ,  1 , 2 , 4 , 20 , 22 );



}


void Z_nunu_Sig_Estimate_both_Cr()
{
	Get_Estimate_from("/afs/cern.ch/work/b/bbachu/private/Z_nunu/" , "monojet");
	Get_Estimate_from("/afs/cern.ch/work/b/bbachu/private/Z_nunu/" , "boosted");
	Get_Estimate_from("/afs/cern.ch/work/b/bbachu/private/Z_nunu/" , "resolved");

}


