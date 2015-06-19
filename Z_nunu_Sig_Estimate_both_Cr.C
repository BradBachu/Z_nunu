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


using namespace std;


void Z_nunu_Sig_Estimate_both_Cr()
{

	// in MC unweighted
	TH1F* h_Znunu_Signal_genVpt_unweighted = (TH1F*) Make_my_hist("Znunu_signal", "genVpt","", 100, 0, 1000);
	TH1F* h_Znunu_Signal_mvamet_unweighted = (TH1F*) Make_my_hist("Znunu_signal", "mvamet","" , 80 , 200, 1000);
	TH1F* h_Photon_Photon_Control_mvamet_unweighted = (TH1F*) Make_my_hist("Photon_photon_control", "mvamet","" , 80 , 200, 1000);
	TH1F* h_Photon_Photon_Control_genVpt_unweighted = (TH1F*) Make_my_hist("Photon_photon_control", "genVpt","" , 100, 0, 1000);
	TH1F* h_Zll_di_muon_control_mvamet_unweighted = (TH1F*) Make_my_hist("Zll_di_muon_control" , "mvamet", "" , 80 , 200, 1000);
	
	//weighted
	TH1F* h_Photon_Photon_Control_genVpt_weighted = (TH1F*) Make_my_hist("Photon_photon_control", "genVpt", "weight" , 100, 0, 1000);
	TH1F* h_Znunu_Signal_genVpt_weighted = (TH1F*) Make_my_hist("Znunu_signal", "genVpt", "weight" ,100, 0, 1000 );
	TH1F* h_Znunu_Signal_mvamet_weighted = (TH1F*) Make_my_hist("Znunu_signal", "mvamet", "weight", 80 , 200, 1000);
	TH1F* h_Zll_di_muon_control_mvamet_weighted = (TH1F*) Make_my_hist("Zll_di_muon_control" , "mvamet" , "weight" , 80 , 200, 1000);

	// in Data
	TH1F* h_data_photon_control_ptpho = (TH1F*) Make_my_hist("data_photon_control", "ptpho","" , 100, 0, 1000);
	TH1F* h_data_di_muon_control_mvamet = (TH1F*) Make_my_hist("data_di_muon_control", "mvamet", "" , 80 , 200, 1000);
	Double_t momentum_bins[19] = {200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310, 320, 330, 350, 380, 430, 500, 1000};

	//compare the MET distributions of both control regions to the signal region
	Draw_2_Hist(h_Photon_Photon_Control_mvamet_unweighted,h_Znunu_Signal_mvamet_unweighted , "Comparision of MET in Z#rightarrow#nu#nu Signal Region to #gamma + Jet Control Region", "MET", "","#gamma + Jet MC" ,"Z#rightarrow#nu#nu MC" , 2, 1, "");
	Draw_2_Hist(h_Zll_di_muon_control_mvamet_unweighted , h_Znunu_Signal_mvamet_unweighted, "Comparision of MET in Z#rightarrow#nu#nu Signal Region to Z#rightarrow#mu#mu Control Region", "MET", "", "Z#rightarrow#mu#mu MC", "Z#rightarrow#mu#mu", 3, 1, "");
	Match_all_integrals_of_hist(h_Znunu_Signal_mvamet_unweighted, h_Photon_Photon_Control_mvamet_unweighted, h_Zll_di_muon_control_mvamet_unweighted , "Comparison of MET with equal integrals", "MET", "", "Z#rightarrow#nu#nu", "#gamma + Jet", "Z#rightarrow#mu#mu", 1, 2, 4, "" , "", "");

	TH1F* h_Znunu_Signal_mvamet_unweighted_rebinned = (TH1F*) Rebin_with_density(h_Znunu_Signal_mvamet_unweighted , momentum_bins , 19);
	TH1F* h_Photon_Photon_Control_mvamet_unweighted_rebinned = (TH1F*) Rebin_with_density(h_Photon_Photon_Control_mvamet_unweighted , momentum_bins , 19);
	TH1F* h_Zll_di_muon_control_mvamet_unweighted_rebinned = (TH1F*) Rebin_with_density(h_Zll_di_muon_control_mvamet_unweighted , momentum_bins , 19);
	Match_all_integrals_of_hist(h_Znunu_Signal_mvamet_unweighted_rebinned, h_Photon_Photon_Control_mvamet_unweighted_rebinned, h_Zll_di_muon_control_mvamet_unweighted_rebinned , "Comparison of MET with equal integrals", "MET", "", "Z#rightarrow#nu#nu", "#gamma + Jet", "Z#rightarrow#mu#mu", 1, 2, 4, "" , "", "");
//TCanvas *cc = new TCanvas();
	//Draw_CMS_Preliminary(momentum_bins ,19, h_Znunu_Signal_mvamet_unweighted_rebinned , h_Photon_Photon_Control_mvamet_unweighted_rebinned, h_Zll_di_muon_control_mvamet_unweighted_rebinned , "Z#rightarrow#nu#nu Signal MC" ,  "Estimate from #gamma + Jet MC", "Estimate from Z#rightarrow#mu#mu MC" , 1, 2, 4);

	// ===============================================
	// ESTIMATE USING THE PHOTON + JET control region

	// the first step is to construt the histogram of the ratios of the Z_pt to the photon_pt
 	TH1F* h_pho_transfer_function_unweighted = (TH1F*) h_Ratio(h_Znunu_Signal_genVpt_unweighted, h_Photon_Photon_Control_genVpt_unweighted);
	Float_t genVpt ;
	Float_t ptpho;
	//Apply this transfer function to the Photon+Jet control region in MC
	TH1F* h_Estimate_from_Pho_MC = (TH1F*) Apply_transfer_function_event_by_event("Photon_photon_control", h_pho_transfer_function_unweighted, genVpt);
	TH1F* h_Estimate_from_MC = (TH1F*) Apply_Pt_Correction("Photon_photon_control", h_pho_transfer_function_unweighted, genVpt);
	//TCanvas *c4 = new TCanvas();
	//h_Estimate_from_MC->Draw();
	//compare the estimate to the MC
	Draw_2_Hist(h_Estimate_from_MC, h_Znunu_Signal_mvamet_unweighted, "Estimate of Z#rightarrow#nu#nu from MC", "MET", "", "Estimate from #gamma + Jet MC", "Z#rightarrow#nu#nu MC", 2, 1, "");
	//now that this should work with MC, we can apply it to data
	//make the histogram that acts as the transfer function for the data
	TH1F* h_pho_transfer_function_weighted = (TH1F*) h_Ratio(h_Znunu_Signal_genVpt_weighted, h_Photon_Photon_Control_genVpt_weighted);
	//apply this transfer function to the data
	TH1F* h_Estimate_from_Pho_Data = (TH1F*) Apply_Pt_Correction("data_photon_control", h_pho_transfer_function_weighted, ptpho);
	Draw_2_Hist(h_Estimate_from_Pho_Data, h_Znunu_Signal_mvamet_weighted, "Estimate of Z#rightarrow#nu#nu from Data", "MET", "", "Estimate from #gamma + Jet Data", "Z#rightarrow#nu#nu", 2 , 1, "Data");

	// ===============================================
	// ESTIMATE USING THE DI-MUON cntrol region
	// the first step is to construct the histogram of the ratios of the Z_nunu MC to Z_mumu MC
	TH1F* h_di_muon_transfer_function_unweighted = (TH1F*) h_Ratio(h_Znunu_Signal_mvamet_unweighted , h_Zll_di_muon_control_mvamet_unweighted);
	//now use this transfer function to correct the shape
	TH1F* h_Estimate_from_di_muon_MC = (TH1F*) Bin_by_Bin_Correction(h_Zll_di_muon_control_mvamet_unweighted , h_di_muon_transfer_function_unweighted);
	Draw_2_Hist(h_Estimate_from_di_muon_MC, h_Znunu_Signal_mvamet_unweighted, "Estimate of Z#rightarrow#nu#nu from Di-muon MC", "MET", "", "Estimate from Z#rightarrow#mu#mu MC", "Z#rightarrow#nu#nu ", 4, 1, "");
	//make the histogram so that we can do the estimate from data
	TH1F* h_di_muon_transfer_function_weighted = (TH1F*) h_Ratio( h_Znunu_Signal_mvamet_weighted , h_Zll_di_muon_control_mvamet_weighted );
	//create the estimate from data
	TH1F* h_Estimate_from_di_muon_Data = (TH1F*) Bin_by_Bin_Correction(h_data_di_muon_control_mvamet , h_di_muon_transfer_function_weighted);
	Draw_2_Hist(h_Estimate_from_di_muon_Data, h_Znunu_Signal_mvamet_weighted, "Estimate of Z#rightarrow#nu#nu from Dimuon Data", "MET", "", "Estiamte from Z#rightarrow#mu#mu Data", "Z#rightarrow#nu#nu", 4, 1, "Data");


	//REBIN ALL THE DATA ESTIMATES
	TH1F* h_Estimate_from_Pho_Data_rebinned = (TH1F*) Rebin_with_density(h_Estimate_from_Pho_Data , momentum_bins , 19);
	TH1F* h_Estimate_from_di_muon_Data_rebinned = (TH1F*) Rebin_with_density(h_Estimate_from_di_muon_Data , momentum_bins , 19);
	TH1F* h_Znunu_Signal_mvamet_weighted_rebinned =(TH1F*) Rebin_with_density(h_Znunu_Signal_mvamet_weighted , momentum_bins , 19);
	//	Draw_3_Hist( h_Znunu_Signal_mvamet_weighted_rebinned , h_Estimate_from_Pho_Data_rebinned,h_Estimate_from_di_muon_Data_rebinned,  "Z#rightarrow#nu#nu MC", "Estimate from #gamma + Jet Data", "Estimate from Z#rightarrow#mu#mu Data" ,  1 , 2 , 4 );
	TCanvas *c5 = new TCanvas();
	Draw_CMS_Preliminary(momentum_bins, 19,  h_Znunu_Signal_mvamet_weighted_rebinned , h_Estimate_from_Pho_Data_rebinned,h_Estimate_from_di_muon_Data_rebinned,  "Z#rightarrow#nu#nu MC", "Estimate from #gamma + Jet Data", "Estimate from Z#rightarrow#mu#mu Data" ,  1 , 2 , 4 );


	// Bin_by_Bin_comparison(h_Estimate_from_di_muon_MC, h_Znunu_Signal_mvamet_weighted);

}


