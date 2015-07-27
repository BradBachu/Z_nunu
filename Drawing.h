#ifndef DRAWING_H
#define DRAWING_H

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


//need to turn the spaces into _ to save my png's

TString Transform_string(TString string) 
{
	TString newstring = string.Copy();
	Int_t size = string.Sizeof();
	//loop through string and replace  space  + or # with underscore
	for ( Int_t i = 0 ; i < size ; i++)
	{
		if (string(i) == char(32) || string(i) == char(35) || string(i) == char(43) || string(i) == char(46) || string(i) == char(58) )
		{
			newstring.Replace( i , 1 , "_" , 1);
		}
		else
		{
			newstring(i) = string(i) ;
		}
	}
return newstring;
}

Float_t Photon_Efficiency_factor(Float_t genVpt)
{
	//Float_t genVpt;

	Float_t matching_factor;
	Float_t scale_factors[5] = {0.9425 , 0.959 , 0.972 , 0.973 , 0.974 } ;
	Float_t bin_edges[6] = { 0 , 20 ,30 , 40 , 50 , 1000} ;
	if ( genVpt < bin_edges[1])
	{
		matching_factor = scale_factors[0];
	}
	else if ( genVpt > bin_edges[1] && genVpt < bin_edges[2] )
	{
		matching_factor = scale_factors[1];
	}
	else if ( genVpt > bin_edges[2]  && genVpt < bin_edges[3] )
	{
		matching_factor = scale_factors[2];
	}
	else if (genVpt > bin_edges[3] && genVpt < bin_edges[4])
	{
		matching_factor = scale_factors[3];
	}
	else
	{
		matching_factor = scale_factors[4];
	}
return matching_factor;
}

//get the efficiency factor associated with the pt of a muon
Float_t Muon_Efficiency_factor( Float_t ptll)
{

	Float_t matching_factor;

	Float_t scale_factors[6] = { 0.958, 0.96, 0.95 , 0.97, 0.985 ,0.99};
	Float_t bin_edges[7] = {0, 20, 40, 65, 100, 150, 300};
	if (ptll < 20)
	{
		matching_factor = scale_factors[0] ;
	}
	else if (ptll > bin_edges[1] && ptll < bin_edges[2])
	{
		matching_factor = scale_factors[1];
	}
	else if (ptll > bin_edges[2] && ptll < bin_edges[3])
	{
		matching_factor = scale_factors[2];
	}
	else if (ptll >bin_edges[3] && ptll <bin_edges[4] )
	{
		matching_factor = scale_factors[3];
	}
	else if (ptll > bin_edges[5] && ptll < bin_edges[6])
	{
		matching_factor = scale_factors[4];
	}
	else 
	{
		matching_factor = scale_factors[5];
	}
return matching_factor;
}

Int_t Category_Cut(TString category , TString variable , Float_t x , Float_t mvamet, Float_t jet1pt)
{
	if ( variable == "mvamet" )
	{
		mvamet = x ;
	}
	if (variable == "jet1pt")
	{
		jet1pt = x ;
	}

	Int_t pass_or_fail_cut; //passing a cut means it satisfied the criteria and we give it value 1
	if (category == "boosted")
	{
		if ( mvamet > 250 && jet1pt > 200 )
		{
			pass_or_fail_cut = 1 ;
		}
		else
		{
			pass_or_fail_cut = 0 ;
		}
		// cout << "Category: " << category << ". mvamet = " << mvamet << ". jet1pt = " << jet1pt << endl;
	}
	else if (category == "resolved")
	{
		//each jet has to have a pt greater than 30 and abs eta < 2.5. 

		if (jet1pt > 30 && mvamet > 250 )
		{
			pass_or_fail_cut = 1;
		}
		else
		{
			pass_or_fail_cut = 0;
		}
		// cout << "Category: " << category << ". mvamet = " << mvamet << ". jet1pt = " << jet1pt << endl;
	}
	else if (category == "monojet")
	{
		if( jet1pt > 30)
		{
			pass_or_fail_cut = 1;
		}
		else
		{
			pass_or_fail_cut = 0;
		}
		// cout << "Category: " << category << ". mvamet = " << mvamet << ". jet1pt = " << jet1pt << endl;
	}
	else
	{
		cout << "MY INPUT DOES NOT MAKE SENSE" << endl;
		cout << "Category: " << category << ". mvamet = " << mvamet << ". jet1pt = " << jet1pt << endl;
	}
	//cout << "Result: " << pass_or_fail_cut << endl;
return pass_or_fail_cut;
}


//Use to import a ntuple from a directory
TH1F* Make_my_hist(TString root_directory , TString category , const Char_t* filename , TString variable, const Char_t *option  , Double_t nbins, Double_t xMin, Double_t xMax , Double_t Min_cut )
{
	TFile *inclusive = new TFile(root_directory);
	TTree *tree = new TTree();
	tree = (TTree*) inclusive->FindObjectAny( filename);
	Float_t x;
	Int_t event;
	//tree->SetBranchAddress(variable , &x);
	Float_t weight;
	tree->SetBranchAddress("weight", &weight);
	//read all entries and fill the hist
	Int_t nentries = (Int_t) tree->GetEntries();
	TH1F *hx = new TH1F( variable,  variable + " distribution  from " + filename + " with option " + option + ". Category: " + category, nbins, xMin, xMax);
	hx->Sumw2();
	Float_t mvamet;
	Float_t jet1pt;

	if ( !(variable == "mvamet" ) && !(variable == "jet1pt" ) )
	{
		tree->SetBranchAddress(variable , &x);
		tree->SetBranchAddress("mvamet", &mvamet);
		tree->SetBranchAddress("jet1pt", &jet1pt);
	}
	else if ( variable == "mvamet")
	{
		tree->SetBranchAddress(variable, &x);
		tree->SetBranchAddress("jet1pt" , &jet1pt);
	}
	else if (variable == "jet1pt")
	{
		tree->SetBranchAddress(variable, &x);
		tree->SetBranchAddress("mvamet", &mvamet);
	}
	else
	{
		cout << "THE VARIABLES / LOGIC DOES NOT MAKE SENSE" << endl;
	}

	//now fill the histogram according to the option
	if (option == "data")
	{
		for (Int_t i = 0 ; i < nentries ; i++)
		{
			tree->GetEntry(i);
			if (x > Min_cut) 
				{
					hx->Fill(x  , Category_Cut( category , variable, x , mvamet, jet1pt) ) ;
				}
		}	
	}
	else if ( option == "weight" )
	{
		for (Int_t i = 0 ; i <nentries ; i++)
		{
			tree->GetEntry(i);
			if (x > Min_cut) 
			{
				hx->Fill(x, weight  * Category_Cut(category , variable ,x, mvamet, jet1pt));
			}
		}
	}
	else if ( option == "unweighted" )
	{
		for (Int_t i = 0 ; i <nentries ; i++)
		{
			tree->GetEntry(i);
			if (x > Min_cut) 
				{
					hx->Fill(x, weight * Category_Cut(category , variable ,x,  mvamet , jet1pt)  );
				}
		}
	}

	else if (option == "weight+efficiency")
	{
		if (filename = "Photon_photon_control")
		{
			//we need to make this messy because the genvpt affects itself
			for (Int_t i = 0; i < nentries ; i++)
			{
				tree->GetEntry(i);
				if (x> Min_cut){hx->Fill(x, weight* Photon_Efficiency_factor(x) * Category_Cut(category , variable, x ,mvamet, jet1pt) );}
			}

		}
		else /*if (filename = "Zll_di_muon_control")*/
		{
			Float_t lep1pt;
			Float_t lep2pt;
			tree->SetBranchAddress("lep1pt", &lep1pt);
			tree->SetBranchAddress("lep2pt", &lep2pt);
					for (Int_t i = 0 ; i< nentries ; i++)
			{
				tree->GetEntry(i);
				if (x>Min_cut){hx->Fill(x , weight * Muon_Efficiency_factor(lep1pt) * Muon_Efficiency_factor(lep2pt)  * Category_Cut(category ,variable , x ,mvamet, jet1pt) );}
			}
		}
	}
	// Background should only be used with the Photon + jet CR
	else if (option == "background")
	{
		if (filename = "Photon_photon_control")
		{
			//we need to make this messy because the genvpt affects itself
			for (Int_t i = 0; i < nentries ; i++)
			{
				tree->GetEntry(i);
				if (x> Min_cut){hx->Fill(x,  (1-Photon_Efficiency_factor(x)) * Category_Cut(category , variable, x ,mvamet, jet1pt) );}
			}

		}
		else if (filename = "Zll_di_muon_control")
		{
			cout << "DONT USE BACKGROUND HERE FOR BACKGROUN"<< endl;
		}
	}
	else 
	{
		cout<< "THE FUNCTION INPUT DOES NOT MAKE SENSE " << filename << " " << variable << endl;
	}

 	// TCanvas* c1 = new TCanvas();
	// hx->Draw();
 	// c1->SaveAs(variable + " distribution  from " + filename + " with option " + option + ". Category: " + category+ ".png");
//	cout <<"end" << endl;
	
return hx;
}

// //rebin a histogram given an array
// TH1F* Rebin_Hist(TH1F* h, Double_t xbins[] , Int_t size)
// {
// 	cout<< "Rebinning histogram" << endl;
// 	TH1F* h2 = h1->Rebin(size,"hnew", xbins);
// 	return h2;
// }

// TH1F* Divide_with_variable_bins(TH1F* h1 , TH1F* h2 , TString category, Double_t xbins[] , Int_t size)
// {
// 	Double_t ngroup = size -1 ;
// 	Int_t nbins ; Double_t xMin ; Double_t xMax;
// 	if (category =="monojet")
// 	{
// 		nbins = 80 ; xMin = 200 ; xMax = 1000 ;
// 	}
// 	else
// 	{
// 		nbins = 75 ; xMin = 250 ; xMax = 1000 ;
// 	}

// 	TH1F* hnew = new TH1F("hnew", "hnew", nbins, xMin , xMax);
// 	hnew->Sumw2();
// 	TH1F* hnew_rebinned = (TH1F*) hnew->Rebin(ngroup , "Division of " + TString(h1->GetName()) +" by " + TString(h2->GetName()) , xbins);
// 	hnew_rebinned->Sumw2();
// 	for (Int_t bin_no = 1 ; bin_no < ngroup +1 ; bin_no++)
// 	{
// 		hnew_rebinned->SetBinContent(bin_no , h1->GetBinContent(bin_no) / h2->GetBinContent(bin_no)  );
// 	}
// 	hnew_rebinned->SetStats(0);
// 	hnew_rebinned->SetMaximum(1.0);
// 	hnew_rebinned->SetMinimum(0.0);
// 	return hnew_rebinned;
// }

TH1F* Divide_with_variable_bins(TH1F* h1 , TH1F* h2 , TString category, Double_t xbins[] , Int_t size)
{
	Double_t ngroup = size -1 ;
	Int_t nbins ; Double_t xMin ; Double_t xMax;
	if (category =="monojet")
	{
		nbins = 80 ; xMin = 200 ; xMax = 1000 ;
	}
	else
	{
		nbins = 75 ; xMin = 250 ; xMax = 1000 ;
	}

	TH1F* hnew = new TH1F("hnew", "hnew", nbins, xMin , xMax);
	hnew->Sumw2();
	TH1F* hnew_rebinned = (TH1F*) hnew->Rebin(ngroup , "Division of " + TString(h1->GetName()) +" by " + TString(h2->GetName()) , xbins);
	hnew_rebinned->Sumw2();
	hnew_rebinned->Divide( h1 , h2 );
	hnew_rebinned->SetStats(0);
	hnew_rebinned->SetMaximum(1.0);
	hnew_rebinned->SetMinimum(0.0);
	return hnew_rebinned;
}

//take a histogram, redraw with the bin contnet divided by bin width
TH1F* Rebin_with_density(TH1F *h , Double_t xbins[] , Int_t size)
{
	//TCanvas *c1 = new TCanvas();
	Double_t ngroup = size - 1 ;
	// FIGURE OUT HOW TO PULL THE HIST NAME IN FROM THE H
	TH1F *hnew = (TH1F*) h->Rebin( ngroup , h->GetTitle() , xbins);
	hnew->Sumw2();
	//loop through the bins and rebinn according to bin density
		// hnew->SetBinContent(0,0.00001);
		// hnew->SetBinContent(20, 0.0001);
	for( Int_t bin_no = 1 ; bin_no < ngroup + 1; bin_no++ )
	{
		//cout<< "Bin #" <<  bin_no << " has content = " << hnew->GetBinContent(bin_no) << " with width " << hnew-> GetBinWidth(bin_no) << endl;
		Double_t Bin_Density = hnew->GetBinContent( bin_no) / hnew->GetBinWidth(bin_no);
		Double_t Error_density = hnew->GetBinError(bin_no) / hnew->GetBinWidth(bin_no);
		hnew->SetBinContent( bin_no , Bin_Density);
		hnew->SetBinError(bin_no , Error_density);
		//cout <<"Bin Density = "<< Bin_Density << " = " <<  hnew->GetBinContent(bin_no) <<endl;
	}
//	hnew->SetMinimum(0.00001);
	hnew->SetStats(0);
	//hnew->Draw();



return hnew;
}

TH1F* Create_Ratio_Hist(TH1F* h1, TH1F* h2 , TString category ,  Double_t xbins[] , Int_t size, Color_t lcolor , Size_t MarkerStyle)
{
	Int_t nbins ; Double_t xMin ; Double_t xMax ;
	if (category == "monojet")
	{
		nbins = 80 ; xMin = 200 ; xMax = 1000;
	}
	else
	{
		nbins = 75 ; xMin = 250 ; xMax = 1000;
	}
	TH1F* h_ratio = new TH1F("h_ratio","Ratio", nbins, xMin , xMax );
	h_ratio->Sumw2();
	TH1F* h_ratio_rebinned = (TH1F*)h_ratio->Rebin(size -1, "h_ratio_rebinned", xbins);
	h_ratio_rebinned->Sumw2();
	// for(Int_t bin_no = 1 ; bin_no < size ; bin_no++)
	// {
	// 	//calculate ratio
	// 	// cout << "Bin Content of h1 = " << h1->GetBinContent(bin_no) << " and Bin Content of h2 = " << h2->GetBinContent(bin_no) << endl;
	// 	Double_t ratio = h2->GetBinContent(bin_no) / h1->GetBinContent(bin_no) ;
	// 	h_ratio_rebinned->SetBinContent(bin_no,ratio);
	// 	//cout << "Bin Content at bin "<< bin_no << " = " << ratio << endl ;
	// }
	h_ratio_rebinned->Divide(h2,h1 );

	h_ratio_rebinned->SetMaximum(1.5); h_ratio->SetMinimum(0.5);
	h_ratio_rebinned->GetYaxis()->SetTitle("Data/bkg");
	h_ratio_rebinned->SetMarkerColor(1);
	h_ratio_rebinned->SetMarkerStyle(MarkerStyle);
	h_ratio_rebinned->SetLineColor(lcolor);
	h_ratio_rebinned->SetMinimum(0.5);
	h_ratio_rebinned->SetMaximum(1.5);
	//TCanvas *c = new TCanvas();
	//h_ratio_rebinned->Draw("e");
return h_ratio_rebinned;
}

//Function to get the y max in a histogram to the max y value between any histogram
Double_t Get_Max_Y(TH1F* h1, TH1F* h2)
{
	Int_t Bin_Location_1;
	Int_t Bin_Location_2;
	Double_t Bin_Content_1; 
	Double_t Bin_Content_2;
	Double_t Max_Y = 0;
	
	Bin_Location_1 = h1->GetMaximumBin();
	Bin_Location_2 = h2->GetMaximumBin();

	Bin_Content_1 = h1->GetBinContent(Bin_Location_1);
	Bin_Content_2 = h2->GetBinContent(Bin_Location_2);

	if (Bin_Content_1 > Bin_Content_2) 
	{
		Max_Y = Bin_Content_1;
	}
	else
	{
		Max_Y = Bin_Content_2;
	}

return (Max_Y * 1.1);
}

Double_t Get_Max_Y_3(TH1F* h1, TH1F* h2 , TH1F* h3){

	Int_t Bin_Location_1;
	Int_t Bin_Location_2;
	Int_t Bin_Location_3;
	Double_t Bin_Content_1; 
	Double_t Bin_Content_2;
	Double_t Bin_Content_3;
	Double_t Max_Y = 0;
	
	Bin_Location_1 = h1->GetMaximumBin();
	Bin_Location_2 = h2->GetMaximumBin();
	Bin_Location_3 = h3->GetMaximumBin();

	Bin_Content_1 = h1->GetBinContent(Bin_Location_1);
	Bin_Content_2 = h2->GetBinContent(Bin_Location_2);
	Bin_Content_3 = h3->GetBinContent(Bin_Location_3);

	if (Bin_Content_1 > Bin_Content_2 && Bin_Content_1 > Bin_Content_3) 
	{
		Max_Y = Bin_Content_1;
	}
	else if (Bin_Content_2> Bin_Content_1 && Bin_Content_2 > Bin_Content_1)
	{
		Max_Y = Bin_Content_2;
	}
	else if (Bin_Content_3> Bin_Content_1 && Bin_Content_3 > Bin_Content_1)
	{
		Max_Y = Bin_Content_3;
	}
	else 
	{
		Max_Y = 0;
		cout <<"YOU ARE DOING RUBBISH!!!" << endl;
	}

return (Max_Y * 1.1);
}

//we also need to use another function that was in zeyneps code
void plot_cms(TString lumi , TCanvas *c4 , TString category)
{

	TLatex *latex2 = new TLatex();
	latex2->SetNDC();
	latex2->SetTextSize(0.035);
	latex2->SetTextAlign(31);  //align right
	latex2->DrawLatex(0.87,0.95, lumi + "fb^{-1}(8TeV)") ;

	TLatex *latex3 = new TLatex();
	latex3->SetNDC();
	latex3->SetTextSize(0.75 * c4->GetTopMargin());
	latex3->SetTextFont(62);
	latex3->SetTextAlign(11); //align right
	latex3->DrawLatex(0.22, 0.85, "CMS");
	latex3->SetTextSize(0.5 * c4->GetTopMargin());
	latex3->SetTextAlign(11);
	latex3->DrawLatex(0.20, 0.80, "Preliminary");

	TLatex *latex1 = new TLatex();
	latex1->SetNDC();
	latex1->SetTextSize(0.035);
	latex1->SetTextAlign(13);
	latex1->SetTextFont(12);
	latex1->DrawLatex(0.1 , 0.95, category);
}

//using this to make it look like how Zeynep wants it
void Draw_CMS_Preliminary_3Stack(TString category, Double_t xbins[] , Int_t size,  
								TH1F *h1 , TH1F *h2, TH1F *h3,  
								const Char_t* label1, const Char_t* label2, const Char_t* label3 ,  
								Color_t lcolor1 , Color_t lcolor2 , Color_t lcolor3 , 
								Size_t MarkerStyle_1, Size_t MarkerStyle_2 , Size_t MarkerStyle_3 ,
								TString draw_option_1, TString draw_option_2, TString draw_option_3,
								TString XAxis, TString YAxis )
{

	const char* s1 = h1->GetName(); const char* s2 = h2->GetName(); const char* s3 = h3->GetName();
	//TCanvas *cc = new TCanvas();
	//from here I am using Zeyeps code
	h2->SetMarkerStyle(MarkerStyle_2); h2->SetMarkerColor(lcolor2); h2->SetMarkerSize(1);
	h3->SetMarkerStyle(MarkerStyle_3); h3->SetMarkerColor(lcolor3); h2->SetMarkerSize(1);

	TLegend *legend = new TLegend(.60, .60, .92,.92);
	legend->AddEntry(h1 , label1 ,"l" );
	legend->AddEntry(h2, label2 , "l");
	legend->AddEntry(h3, label3 , "p"); 
	
	TCanvas *c4 = new TCanvas( Transform_string(TString(s1)) + "_" + Transform_string(TString(s2)) + "_"  + Transform_string(TString(s3)) , Transform_string(TString(s1)) + "_" + Transform_string(TString(s2)) + "_"  + Transform_string(TString(s3)) + "_"+ category + "THSTACK",  900, 1000);
	c4->SetBottomMargin(0.3);
	c4->SetRightMargin(0.06);

	c4->SetLogy();

	h1->SetTitle("");
	h1->Draw("HIST");
	h1->GetYaxis()->SetTitle("Events/GeV");
	h1->GetYaxis()->CenterTitle();
	h1->GetXaxis()->SetTitle("fake MET GeV");
	h1->GetXaxis()->SetLabelSize(0);
	h1->GetXaxis()->SetTitle("");
	h1->SetMaximum(Get_Max_Y_3(h1, h2 , h3));

	h2->SetTitle("");
	h2->Draw("EsameHIST");
	h3->SetTitle("");
	h3->Draw("EsameHIST");

	THStack *hs = new THStack("hs","Stacked 1D histograms");
	//draw background first
	hs->Add(h1);
	hs->Add(h2);
	// draw same for data
	hs->Add(h3 , "same");

	legend->SetShadowColor(0);
	legend->SetFillColor(0);
	legend->SetLineColor(1);

	legend->Draw("same");
	plot_cms("19.7", c4 , category);

	TPad* Pad = new TPad("pad", "pad", 0.0,0.0,1.0,1.0);
	Pad->SetTopMargin(0.7);
	Pad->SetFillColor(0);
	Pad->SetGridy(1);
	Pad->SetFillStyle(0);
	Pad->Draw();
	Pad->cd(0);
	Pad->SetRightMargin(0.06);

	TH1F* ratio_2 = (TH1F*) Create_Ratio_Hist(h1 , h2 , category , xbins, size , lcolor2 , MarkerStyle_2);
	TH1F* ratio_3 = (TH1F*) Create_Ratio_Hist(h1 , h3 , category , xbins , size ,  lcolor3 , MarkerStyle_3);
	ratio_2->SetStats(0);
	ratio_2->SetTitle("");
	ratio_2->SetLineColor(lcolor2);
	ratio_2->GetXaxis()->SetTitle(XAxis);
	ratio_2->SetMarkerColor(lcolor2);
	ratio_2->SetMarkerStyle(MarkerStyle_2);
	ratio_2->Draw("E");	
	ratio_3->SetStats(0);
	ratio_3->SetTitle("");
	ratio_3->SetLineColor(lcolor3);
	ratio_3->SetMarkerColor(lcolor3);
	ratio_3->SetMarkerStyle(MarkerStyle_3);
	ratio_3->GetXaxis()->SetTitle(XAxis);
	ratio_3->Draw("Esame");


	TF1* f1 = new TF1("f1", "1", -5000, 5000);
	f1->SetLineColor(4);
	f1->SetLineStyle(2);
	f1->SetLineWidth(2);
	f1->Draw("same");
	c4->SaveAs("/afs/cern.ch/work/b/bbachu/private/Z_nunu/Hists/" + Transform_string(TString(s1)) + "_" + Transform_string(TString(s2)) + "_" + Transform_string(TString(s3)) + "_" + category +"STACK.png");
	//TCanvas *c = new TCanvas();
}

//using this to make it look like how Zeynep wants it
void Draw_CMS_Preliminary_3h(TString category, Double_t xbins[] , Int_t size,  TH1F *h1 , TH1F *h2, TH1F *h3,  const Char_t* label1, const Char_t* label2, const Char_t* label3 ,  Color_t lcolor1 , Color_t lcolor2 , Color_t lcolor3 , Size_t MarkerStyle_2 , Size_t MarkerStyle_3 , TString XAxis )
{

	const char* s1 = h1->GetName(); const char* s2 = h2->GetName(); const char* s3 = h3->GetName();
	//TCanvas *cc = new TCanvas();
	//from here I am using Zeyeps code
	h2->SetMarkerStyle(MarkerStyle_2); h2->SetMarkerColor(lcolor2); h2->SetMarkerSize(1);
	h3->SetMarkerStyle(MarkerStyle_3); h3->SetMarkerColor(lcolor3); h2->SetMarkerSize(1);

	TLegend *legend = new TLegend(.60, .60, .92,.92);
	legend->AddEntry(h1 , label1 ,"l" );
	legend->AddEntry(h2, label2 , "p");
	legend->AddEntry(h3, label3 , "p"); 
	
	TCanvas *c4 = new TCanvas( Transform_string(TString(s1)) + "_" + Transform_string(TString(s2)) + "_"  + Transform_string(TString(s3)) , Transform_string(TString(s1)) + "_" + Transform_string(TString(s2)) + "_"  + Transform_string(TString(s3)) + "_"+ category,  900, 1000);
	c4->SetBottomMargin(0.3);
	c4->SetRightMargin(0.06);

	c4->SetLogy();

	h1->SetTitle("");
	h1->Draw("HIST");
	h1->GetYaxis()->SetTitle("Events/GeV");
	h1->GetYaxis()->CenterTitle();
	h1->GetXaxis()->SetTitle("fake MET GeV");
	h1->GetXaxis()->SetLabelSize(0);
	h1->GetXaxis()->SetTitle("");
	h1->SetMaximum(Get_Max_Y_3(h1, h2 , h3));

	h2->SetTitle("");
	h2->Draw("EsameHIST");
	h3->SetTitle("");
	h3->Draw("EsameHIST");

	legend->SetShadowColor(0);
	legend->SetFillColor(0);
	legend->SetLineColor(1);

	legend->Draw("same");
	plot_cms("19.7", c4 , category);

	TPad* Pad = new TPad("pad", "pad", 0.0,0.0,1.0,1.0);
	Pad->SetTopMargin(0.7);
	Pad->SetFillColor(0);
	Pad->SetGridy(1);
	Pad->SetFillStyle(0);
	Pad->Draw();
	Pad->cd(0);
	Pad->SetRightMargin(0.06);

	TH1F* ratio_2 = (TH1F*) Create_Ratio_Hist(h1 , h2 , category , xbins, size , lcolor2 , MarkerStyle_2);
	TH1F* ratio_3 = (TH1F*) Create_Ratio_Hist(h1 , h3 , category , xbins , size ,  lcolor3 , MarkerStyle_3);
	ratio_2->SetStats(0);
	ratio_2->SetTitle("");
	ratio_2->SetLineColor(lcolor2);
	ratio_2->GetXaxis()->SetTitle(XAxis);
	ratio_2->SetMarkerColor(lcolor2);
	ratio_2->SetMarkerStyle(MarkerStyle_2);
	ratio_2->Draw("E");	
	ratio_3->SetStats(0);
	ratio_3->SetTitle("");
	ratio_3->SetLineColor(lcolor3);
	ratio_3->SetMarkerColor(lcolor3);
	ratio_3->SetMarkerStyle(MarkerStyle_3);
	ratio_3->GetXaxis()->SetTitle(XAxis);
	ratio_3->Draw("Esame");


	TF1* f1 = new TF1("f1", "1", -5000, 5000);
	f1->SetLineColor(4);
	f1->SetLineStyle(2);
	f1->SetLineWidth(2);
	f1->Draw("same");
	c4->SaveAs("/afs/cern.ch/work/b/bbachu/private/Z_nunu/Hists/" + Transform_string(TString(s1)) + "_" + Transform_string(TString(s2)) + "_" + Transform_string(TString(s3)) + "_" + category +".png");
	//TCanvas *c = new TCanvas();
}

//using this to make it look like how Zeynep wants it
void Draw_CMS_Preliminary_4h(TString category, Double_t xbins[] , Int_t size,  TH1F *h1 , TH1F *h2, TH1F *h3, TH1F *h4 ,  const Char_t* label1, const Char_t* label2, const Char_t* label3 , const Char_t* label4, Color_t lcolor1 , Color_t lcolor2 , Color_t lcolor3 , Color_t lcolor4 ,  Size_t MarkerStyle_2 , Size_t MarkerStyle_3 , Size_t MarkerStyle_4 , TString XAxis )
{

	const char* s1 = h1->GetName(); const char* s2 = h2->GetName(); const char* s3 = h3->GetName(); const char* s4 = h4->GetName();
	//TCanvas *cc = new TCanvas();
	//from here I am using Zeyeps code
	h2->SetMarkerStyle(MarkerStyle_2); h2->SetMarkerColor(lcolor2); h2->SetMarkerSize(1);
	h3->SetMarkerStyle(MarkerStyle_3); h3->SetMarkerColor(lcolor3); h3->SetMarkerSize(1);
	h4->SetMarkerStyle(MarkerStyle_4); h3->SetMarkerColor(lcolor4); h4->SetMarkerSize(1);

	TLegend *legend = new TLegend(.60, .60, .92,.92);
	legend->AddEntry(h1 , label1 ,"l" );
	legend->AddEntry(h2, label2 , "p");
	legend->AddEntry(h3, label3 , "p"); 
	legend->AddEntry(h4, label4 , "p"); 
	
	TCanvas *c4 = new TCanvas( Transform_string(TString(s1)) + "_" + Transform_string(TString(s2)) + "_"  + Transform_string(TString(s3)) , Transform_string(TString(s1)) + "_" + Transform_string(TString(s2)) + "_"  + Transform_string(TString(s3)) + "_"+ Transform_string(TString(s2))+ "_"+ category,  900, 1000);
	c4->SetBottomMargin(0.3);
	c4->SetRightMargin(0.06);

	c4->SetLogy();

	h1->SetTitle("");
	h1->Draw("HIST");
	h1->GetYaxis()->SetTitle("Events/GeV");
	h1->GetYaxis()->CenterTitle();
	h1->GetXaxis()->SetTitle("fake MET GeV");
	h1->GetXaxis()->SetLabelSize(0);
	h1->GetXaxis()->SetTitle("");
	h1->SetMaximum(Get_Max_Y_3(h1, h2 , h3));

	h2->SetTitle("");
	h2->Draw("EsameHIST");
	h3->SetTitle("");
	h3->Draw("EsameHIST");
	h4->SetTitle("");
	h4->Draw("EsameHIST");

	legend->SetShadowColor(0);
	legend->SetFillColor(0);
	legend->SetLineColor(1);

	legend->Draw("same");
	plot_cms("19.7", c4 , category);

	TPad* Pad = new TPad("pad", "pad", 0.0,0.0,1.0,1.0);
	Pad->SetTopMargin(0.7);
	Pad->SetFillColor(0);
	Pad->SetGridy(1);
	Pad->SetFillStyle(0);
	Pad->Draw();
	Pad->cd(0);
	Pad->SetRightMargin(0.06);

	TH1F* ratio_2 = (TH1F*) Create_Ratio_Hist(h1 , h2 , category , xbins, size , lcolor2 , MarkerStyle_2);
	TH1F* ratio_3 = (TH1F*) Create_Ratio_Hist(h1 , h3 , category , xbins , size ,  lcolor3 , MarkerStyle_3);
	TH1F* ratio_4 = (TH1F*) Create_Ratio_Hist(h1 , h4 , category , xbins , size , lcolor4 , MarkerStyle_4);
	ratio_2->SetStats(0);
	ratio_2->SetTitle("");
	ratio_2->SetLineColor(lcolor2);
	ratio_2->GetXaxis()->SetTitle(XAxis);
	ratio_2->SetMarkerColor(lcolor2);
	ratio_2->SetMarkerStyle(MarkerStyle_2);
	ratio_2->Draw("E");	
	ratio_3->SetStats(0);
	ratio_3->SetTitle("");
	ratio_3->SetLineColor(lcolor3);
	ratio_3->SetMarkerColor(lcolor3);
	ratio_3->SetMarkerStyle(MarkerStyle_3);
	ratio_3->GetXaxis()->SetTitle(XAxis);
	ratio_3->Draw("Esame");
	ratio_4->SetStats(0);
	ratio_4->SetTitle("");
	ratio_4->SetLineColor(lcolor4);
	ratio_4->SetMarkerColor(lcolor4);
	ratio_4->SetMarkerStyle(MarkerStyle_4);
	ratio_4->GetXaxis()->SetTitle(XAxis);
	ratio_4->Draw("Esame");


	TF1* f1 = new TF1("f1", "1", -5000, 5000);
	f1->SetLineColor(4);
	f1->SetLineStyle(2);
	f1->SetLineWidth(2);
	f1->Draw("same");
	c4->SaveAs("/afs/cern.ch/work/b/bbachu/private/Z_nunu/Hists/" + Transform_string(TString(s1)) +"_" + Transform_string(TString(s2)) + "_" + Transform_string(TString(s3)) +"_" + Transform_string(TString(s4)) + "_" + category +".png");
	//TCanvas *c = new TCanvas();
}

//using this to make it look like how Zeynep wants it
void Draw_CMS_Preliminary_2h(TString category, Double_t xbins[] , Int_t size,  TH1F *h1 , TH1F *h2,  const Char_t* label1, const Char_t* label2, Color_t lcolor1 , Color_t lcolor2 , Size_t MarkerStyle_1 , Size_t MarkerStyle_2 , TString XAxis)
{

	const char* s1 = h1->GetName(); const char* s2 = h2->GetName(); TString s3 = category ; 
	//TCanvas *c = new TCanvas();
	//from here I am using Zeyeps code
	h1->SetMarkerStyle(MarkerStyle_1); h1->SetMarkerColor(lcolor1); h1->SetMarkerSize(1);
	h2->SetMarkerStyle(MarkerStyle_2); h2->SetMarkerColor(lcolor2); h2->SetMarkerSize(1);

	TLegend *legend = new TLegend(.60, .60, .92,.92);
	legend->AddEntry(h1 , label1 ,"l" );
	legend->AddEntry(h2, label2 , "p");

	TCanvas *c4 = new TCanvas( Transform_string(TString(s1)) + "_" + Transform_string(TString(s2)) + "_" + category , Transform_string(TString(s1)) + "_" + Transform_string(TString(s2)) + "_" +category ,  900, 1000);
	c4->SetBottomMargin(0.3);
	c4->SetRightMargin(0.06);

	c4->SetLogy();

	h1->SetTitle("");
	h1->Draw("HIST");
	h1->SetLineColor(lcolor1);
	h1->GetYaxis()->SetTitle("Events/GeV");
	h1->GetYaxis()->CenterTitle();
	h1->GetXaxis()->SetTitle("fake MET GeV");
	h1->GetXaxis()->SetLabelSize(0);
	h1->GetXaxis()->SetTitle("");
	h1->SetMaximum( Get_Max_Y(h1,h2) );
	h2->SetTitle("");
	h2->SetLineColor(lcolor2);
	h2->Draw("EsameHIST");

	legend->SetShadowColor(0);
	legend->SetFillColor(0);
	legend->SetLineColor(1);

	legend->Draw("same");
	plot_cms(" 19.7 ", c4, category);

	TPad* Pad = new TPad("pad", "pad", 0.0,0.0,1.0,1.0);
	Pad->SetTopMargin(0.7);
	Pad->SetFillColor(0);
	Pad->SetGridy(1);
	Pad->SetFillStyle(0);
	Pad->Draw();
	Pad->cd(0);
	Pad->SetRightMargin(0.06);

	TH1F* ratio_2 = (TH1F*) Create_Ratio_Hist(h1,h2 , category, xbins, size , lcolor2 , MarkerStyle_2);
	ratio_2->SetStats(0);
	ratio_2->SetTitle("");
	ratio_2->SetLineColor(lcolor2);
	ratio_2->SetMarkerStyle(MarkerStyle_2); ratio_2->SetMarkerColor(lcolor2); h2->SetMarkerSize(1);
	ratio_2->Draw("E");	
	ratio_2->GetXaxis()->SetTitle(XAxis);

	TF1* f1 = new TF1("f1", "1", -5000, 5000);
	f1->SetLineColor(4);
	f1->SetLineStyle(2);
	f1->SetLineWidth(2);
	f1->Draw("same");
	//cout << "!!! " << "/afs/cern.ch/work/b/bbachu/private/Z_nunu/" << Transform_string(TString(h1->GetTitle())) << TString("_") << Transform_string(h2->GetTitle()) << TString("_") << Transform_string(filename.Data()) << TString(".png") << endl;
	//c4->SaveAs("/afs/cern.ch/work/b/bbachu/private/Z_nunu/" + Transform_string(TString(h1->GetTitle())) +"_" + Transform_string(TString(h2->GetTitle())) + "_"+ Transform_string(category) + ".png");
	//printf("%s\n",s1);
	cout << s1 << endl;
	cout << s2 << endl;
	c4->SaveAs("/afs/cern.ch/work/b/bbachu/private/Z_nunu/Hists/" + Transform_string(TString(s1))  + Transform_string(TString(s2)) + "_" + category +".png");
	//TCanvas *cc = new TCanvas();

}


//Draw 2 legends in a (2)histogram
// DO NOT USE THIS!. It is used as a sub function
void Draw2Legend(TH1F *histo1,TH1F *histo2,const Char_t *label1,const Char_t *label2)
{
// example for adding a legend to the plot of two histos
  TLegend *legend=new TLegend(0.7,0.75,0.89,0.89);  //geometry of legend
  //legend->SetHeader(header);  //leg. header (name of histogram or sth. else)
  legend->AddEntry(histo1,label1,"L");  // see *http://root.cern.ch/root/html/TLegend.html *  for details
  legend->AddEntry(histo2,label2,"L");
  legend->SetBorderSize(1);  //no border for legend
  legend->SetFillColor(0);  //fill color is white
  legend->Draw();
}

// DO NOT USE THIS!. It is used as a sub function
void Draw3Legend(TH1F *histo1,TH1F *histo2, TH1F *histo3,  const Char_t *label1,const Char_t *label2 , const Char_t* label3)
{
// example for adding a legend to the plot of two histos
  TLegend *legend=new TLegend(0.7,0.75,0.89,0.89);  //geometry of legend
  //legend->SetHeader(header);  //leg. header (name of histogram or sth. else)
  legend->AddEntry(histo1,label1,"L");  // see *http://root.cern.ch/root/html/TLegend.html *  for details
  legend->AddEntry(histo2,label2,"L");
  legend->AddEntry(histo3, label3, "L");
  legend->SetBorderSize(1);  //no border for legend
  legend->SetFillColor(0);  //fill color is white
  legend->Draw();
}





void Draw_1_Hist(TH1F* h , TString title , TString YAxis ,TString XAxis, Color_t color , TString option ,TString category , Double_t Ymin , Double_t Ymax)
{
	const char* s1 = h->GetName();
	TCanvas *c = new TCanvas();
	TH1F* h1 = (TH1F*) h->Clone("h1");
	h1->SetTitle(category + ": " + title);
	h1->GetXaxis()->SetTitle(XAxis);
	h1->GetYaxis()->SetTitle(YAxis);
	h1->SetLineColor(color);
	h1->SetMinimum(Ymin);
	h1->SetMaximum(Ymax);
	// h1->SetMarkerStyle(MarkerStyle); 
	// h1->SetMarkerColor(color); 
	// h1->SetMarkerSize( markersize);
	h1->SetMarkerStyle(20); h1->SetMarkerColor(color); h1->SetMarkerSize(1);
	h1->Draw("E");

	c->SaveAs("/afs/cern.ch/work/b/bbachu/private/Z_nunu/Hists/" + Transform_string(TString(s1)) + Transform_string(TString(h->GetYaxis()->GetTitle() )) + Transform_string(category)+".png");

}

//MAIN DRAW 2 HIST FUNCTION
//Draw 2 histograms with 2 legend
void Draw_2_Hist(TH1F* h_1, TH1F* h_2, const Char_t *Name, const Char_t *XAxis, const Char_t *YAxis, const Char_t *label1, const Char_t * label2, Color_t lcolor1, Color_t lcolor2, const Char_t* h1source)
{
	const char* s1 = h_1->GetName() ; const char* s2 = h_2->GetName();
	TCanvas* tcanvas = new TCanvas( TString(Name) , TString(Name) , 700, 500);
	//tcanvas->SetLogy();
	TH1F *h1 = (TH1F*) h_1->Clone("h1");
	TH1F *h2 = (TH1F*) h_2->Clone("h2");
	h1->GetXaxis()->SetTitle(XAxis); 
	h2->GetXaxis()->SetTitle(XAxis);
	h1->GetYaxis()->SetTitle(YAxis);
	h2->GetYaxis()->SetTitle(YAxis);
	h1->SetTitle(Name);
	h2->SetTitle(Name);
	h1->SetLineColor(lcolor1);
	h2->SetLineColor(lcolor2);
	Double_t Max = Get_Max_Y(h1,h2);
	h1->SetMaximum(Max);
	h2->SetMaximum(Max);
	Double_t Min = 0;
	h1->SetMinimum(Min);
	h2->SetMinimum(Min);
	h1->SetStats(0);
	h2->SetStats(0);
	// h1->SetFillColorAlpha(lcolor1,0.35);
	// h2->SetFillColorAlpha(lcolor2,0.35);
	h2->Draw();
	if ( h1source == "Data")
	{
		h1->Draw("P E1 same");
	}
	else
	{
		h1->Draw("same");
	}	
	Draw2Legend(h1, h2, label1, label2);
	tcanvas->SaveAs("/afs/cern.ch/work/b/bbachu/private/Z_nunu/Hists/" + Transform_string(TString(Name)) + Transform_string(TString(s1)) + "_" + Transform_string(TString(s2)) + ".png");


}


//MAIN DRAW 2 HIST FUNCTION
//Draw 2 histograms with 2 legend
void Draw_3_Hist(TH1F* h_1, TH1F* h_2, TH1F* h_3 , const Char_t *Name, const Char_t *XAxis, const Char_t *YAxis, const Char_t *label1, const Char_t * label2, const Char_t *label3, Color_t lcolor1, Color_t lcolor2, Color_t lcolor3, const Char_t* h1source , const Char_t* h2source, const Char_t* h3source)
{
	const char* s1 = h_1->GetName() ; const char* s2 = h_2->GetName() ; const char* s3 = h_3->GetName();
	TCanvas* tcanvas = new TCanvas(TString(Name) , TString(Name), 700 , 500  );
	tcanvas->SetLogy();
	TH1F *h1 = (TH1F*) h_1->Clone("h1");
	TH1F *h2 = (TH1F*) h_2->Clone("h2");
	TH1F *h3 = (TH1F*) h_3->Clone("h3");
	h1->GetXaxis()->SetTitle(XAxis); 
	h2->GetXaxis()->SetTitle(XAxis);
	h3->GetXaxis()->SetTitle(XAxis);
	h1->GetYaxis()->SetTitle(YAxis);
	h2->GetYaxis()->SetTitle(YAxis);
	h3->GetYaxis()->SetTitle(YAxis);
	h1->SetTitle(Name);
	h2->SetTitle(Name);
	h3->SetTitle(Name);
	h1->SetLineColor(lcolor1);
	h2->SetLineColor(lcolor2);
	h3->SetLineColor(lcolor3);
	Double_t Max = Get_Max_Y_3(h1,h2,h3);
	h1->SetMaximum(Max);
	h2->SetMaximum(Max);
	h3->SetMaximum(Max);
	Double_t Min = 1;
	h1->SetMinimum(Min);
	h2->SetMinimum(Min);
	h3->SetMinimum(Min);
	//h1->SetFillColorAlpha(lcolor1,0.35);
	//h2->SetFillColorAlpha(lcolor2,0.35);
	h1->SetStats(0);
	h2->SetStats(0);
	h3->SetStats(0);

	if ( h1source == "Data")
	{
		h1->Draw("P E1 ");
	}
	else
	{
		h1->Draw();
	}	
	if ( h2source == "Data")
	{
		h2->Draw("P E1 same");
	}
	else
	{
		h2->Draw("same");
	}	

	if ( h3source == "Data")
	{
		h3->Draw("P E1 same");
	}
	else
	{
		h3->Draw("same");
	}	
	//h2->Draw("same");
	Draw3Legend(h1, h2, h3 ,label1, label2, label3);
	tcanvas->SaveAs("/afs/cern.ch/work/b/bbachu/private/Z_nunu/Hists" + Transform_string( TString(s1) )+"_" + Transform_string(TString(s2)) + "_" + Transform_string(TString(s3)) + ".png");
}




// Method to make Histograms match by integral
TH1F* Scale_h1_to_h2(TH1F *h1, TH1F *h2)
{
	//clone both histograms 
	TH1F *h1a = (TH1F*) h1->Clone("h1a");
	Double_t transfer_const = (h2->Integral()) / (h1a->Integral());
	h1a->Scale(transfer_const);
	return h1a;
}

//give me 3 histograms and i will replot them to match the integrals of the first one
void Match_all_integrals_of_hist(TH1F* h1, TH1F* h2, TH1F* h3 , const Char_t *Name, const Char_t *XAxis, const Char_t *YAxis, const Char_t *label1, const Char_t * label2, const Char_t *label3, Color_t lcolor1, Color_t lcolor2, Color_t lcolor3, const Char_t* h1source , const Char_t* h2source, const Char_t* h3source)
{
	TH1F* h2_scaled_to_h1 = (TH1F*) Scale_h1_to_h2(h2 , h1);
	TH1F* h3_scaled_to_h1 = (TH1F*) Scale_h1_to_h2(h3, h1);
	Draw_3_Hist(h1, h2_scaled_to_h1, h3_scaled_to_h1, Name, XAxis, YAxis, label1, label2, label3, lcolor1,  lcolor2,  lcolor3,  h1source ,h2source,  h3source);
}

#endif