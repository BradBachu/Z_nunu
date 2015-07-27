gSystem->AddIncludePath(TString("-l/cvmfs/cms.cern.ch/")+TString(gSystem->Getenv("SCRAM_ARCH"))+
			       TString("/lcg/roofit/5.32.00/include/"));
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooExponential.h"
#include "TCanvas.h"
#include "RooPlot.h"
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
#include "RooAddPdf.h"
#include "vector.h"
#include "My_Functions.h"



typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> > LorentzVector;


using namespace RooFit ;

// Example of fitting Gaussian PDF to data
void Fit1()
{
	cout<<"Start"<<endl;
	
	//Gauss_PDF();
	
	TCanvas *c1 = new TCanvas();
	//Build Gaussian PDF
	RooRealVar x("x","x", -10, 10);
	RooRealVar mean("mean", "mean of gaussian", 0, -10, 10);
	RooRealVar sigma("sigma", "width of gaussian", 3);
	// Define the PDF for a gaussian 
	RooGaussian gauss("gauss", "gaussian PDF", x, mean, sigma);
	
	//plot PDF
	RooPlot *xframe = x.frame();
	gauss.plotOn(xframe);

	//change sigma and plot on same frame
	sigma = 2;
	gauss.plotOn(xframe, LineColor(7));
	xframe->Draw();
	cout << "Created Gauss PDF" << endl ;
	TImage *img = TImage::Create();
	img->FromPad(c1);
	img->WriteImage("canvas1.png");
	
	cout << "Starting to create random binned data" << endl;
	TCanvas *c = new TCanvas();
	cout << "creating histrogram" << endl;
	TH1D *h = new TH1D("ahisto", "ahisto", 100, -10, 10);
	cout << "Filling histogram" << endl;
	h->FillRandom("gaus",100);
	h->Draw();
	
	cout << "Processing events" << endl;
	gSystem->ProcessEvents();
	TImage *img = TImage::Create();
	
	img->FromPad(c);
	img->WriteImage("canvas.png");

	cout << "Created random histogram"<< endl;

	//TH1* hh = (TH1*) gDirectory->Get("ahisto")
	//cout  << "Beginning the import of data"<< endl
	TCanvas *c2 = new TCanvas();
	cout << "Create new tcanvas" << endl;
	//You can import the contents of any ROOT histogram in a RooDataHist object
	//RooRealVar z("z", "z", -10, 10);
	cout <<"real var created"<< endl;
	RooDataHist data("data", "dataset with x", x, h, 1.0);
	cout <<"defined var and datahist" << endl;
	//A RooDataHist always associates the histogram with a RooFit variable object of type RooRealVar
	//In this way it always knows what kind of data is stored in the histogram
	
	//A RooDataHist can be visualized in the same way as a function can be visualized
	RooPlot *xframe2 = x.frame();
	cout << "roo plot "<< endl;
	data.plotOn(xframe2);
	xframe2->Draw();
	TImage *img = TImage::Create();
	img->FromPad(c2);
	img->WriteImage("canvas2.png");

	gauss.fitTo(data);
	mean.Print();
	sigma.Print();
	
	TCanvas *c4 = new TCanvas();
	RooPlot* xframe3 = x.frame();
	data.plotOn(xframe3);
	gauss.plotOn(xframe3);
	xframe3->Draw();

}


//Read the N-Tuple and print on screen its content
TH1F* Read_NTuple(const char* filename , const char* variable)
{
	TCanvas *c1 = new TCanvas();
	TFile *inclusive = new TFile("/home/bbachu/cms/cmssw/032/CMSSW_5_3_14/src/MitMonoJet/macros/Fits/inclusive.root");
	//inclusive->ls();
	//TTree *tree =  (TTree*) inclusive->FindObjectAny("Photon_photon_control");
	TTree *tree = new TTree();
	tree = (TTree*) inclusive->FindObjectAny(filename);
	tree->Draw(variable);
	TH1F *htemp = gROOT->FindObject("htemp");
	Int_t xbins = htemp->GetNbinsX();
	cout << "Nbins = " << xbins << endl;
	//htemp->Draw();
Draw_1_Hist( htemp,variable, variable, filename, 4, "MC");
return htemp;
}

TH1F* Get_Hist(const char* ttree_name)
{
	TCanvas* c1 = new TCanvas();
	TFile *inclusive = new TFile("/home/bbachu/cms/cmssw/032/CMSSW_5_3_14/src/MitMonoJet/macros/Fits/inclusive.root");
	//inclusive->ls();
	//TTree *tree =  (TTree*) inclusive->FindObjectAny("Photon_photon_control");
	TTree *tree = new TTree();
	tree = (TTree*) inclusive->FindObjectAny( ttree_name);
	Float_t mvamet;
	Int_t event;
	tree->SetBranchAddress("mvamet", &mvamet);
	//create hist
	TH1F *hmvamet = new TH1F("hmvamet", "distribution", 100, 0, 1000);
	//read all entries and fill the hist
	Int_t nentries = (Int_t) tree->GetEntries();
	for (Int_t i = 0; i<nentries; i++ )
	{
		tree->GetEntry(i);
		hmvamet->Fill(mvamet);
	}
	hmvamet->Draw();
	return hmvamet;
}


TH1F* Make_Hist(const char* filename , TString variable)
{
	TCanvas* c1 = new TCanvas();
	TFile *inclusive = new TFile("/home/bbachu/cms/cmssw/032/CMSSW_5_3_14/src/MitMonoJet/macros/Fits/inclusive.root");
	TTree *tree = new TTree();
	tree = (TTree*) inclusive->FindObjectAny( filename);
	Float_t x;
	Int_t event;
	tree->SetBranchAddress(variable , &x);
	//create hist
	TH1F *hx = new TH1F( variable,  variable + " distribution", 100, 0, 1000);
	//read all entries and fill the hist
	Int_t nentries = (Int_t) tree->GetEntries();
	for (Int_t i = 0; i<nentries; i++ )
	{
		tree->GetEntry(i);
		hx->Fill(x);
	}
	hx->Draw();
	return hx;
}

// Dont use
// First try: Fits exp nd gauss to data
void Fit2(const TH1F* histo = 0, const char* variable)
{
	TCanvas *c1 = new TCanvas();
	cout <<"Start the fitting of the photon+Jet MC for variable " << variable <<  endl;
	
	//To start we will use a exponential (although a double exponential is required
	//variable for pdf
	RooRealVar *mvamet =  new RooRealVar(variable, variable, 0, 1000);

	//create exponential function
	RooRealVar *a =  new RooRealVar("a","a", 0.003 ,-0.100, 0.100);
	RooExponential *exp_fit = new RooExponential("Exponential", "Exponential fit to mvamet", *mvamet , *a);

	//create gaussian
	RooRealVar *mean = new RooRealVar("mean", "mean of gaussian", 180, -300, 300);
	RooRealVar *sigma = new RooRealVar("sigma", "width of gaussian", 200, -1000, 1000);

	//If you want to fix a parameter in the fit, you just specify that as a property of the RooRewalVar parameter object
	mean->setConstant(kTRUE);
	RooGaussian *gauss = new RooGaussian("gauss", "gaussian PDF", *mvamet, *mean, *sigma);
	
	//Import the histrogram with the required variable
	TH1F *hh = histo ;
	//make histogram for Gaussian fit
	RooDataHist *data = (RooDataHist*) new RooDataHist("data", "dataset with mvamet", *mvamet , hh);
	//make histgoram for exponential fit
	TH1F *hexp = (TH1F*) hh->Clone("hexp");
	RooDataHist *data_exp = (RooDataHist*) new RooDataHist("data", "dataset with mvamet", *mvamet , hexp);
	
	//Fit Gaussian to all data
	RooFitResult *r_full_gauss = gauss->fitTo(*data, Save(kTRUE));
	//Fit exp to all data
	RooFitResult *r_full_exp = exp_fit->fitTo(*data, Save(kTRUE));

	//Define main signal region in data as [175, 1000]
	mvamet->setRange("signal", 175, 1000);

	// Fit gaussian pdf only to data in the "signal" region
	RooFitResult *r_sig_gauss = gauss->fitTo(*data, Save(kTRUE), Range("signal"));
	// Fit exponential pdf to data in the "signal" region 
	RooFitResult *r_sig_exp = exp_fit->fitTo(*data, Save(kTRUE), Range("signal"));

	//PLOT / Print results
	//----------------------------

	//make plost fame in mvamet and add data and fitted model
	RooPlot *frame = mvamet.frame(Title("Fitting a subrange"));
	data->plotOn(frame);
	//exp_fit->plotOn(frame, Range("Full"), LineStyle(kDashed), LineColor(kRed));//Add shape in full range dashed
	exp_fit->plotOn(frame, Range("signal"), LineColor(kRed));//by default only fitted range is shown
	gauss->plotOn(frame, Range("signal"), LineColor(kBlue));
	// //Print fit results
	cout<<"Result of fit on all data"<< endl;
	r_full_exp->Print();
	cout<< "Result of fit on signal "<< endl;
	r_sig_exp->Print();
	frame->Draw();
}

//first try to fit a gaussian and exponential to data
// did not use RooResult
void Fit_Pho_Jet(const TH1F* histo = 0, const char* variable)
{
	TCanvas *c1 = new TCanvas();
	cout <<"Start the fitting of the photon+Jet MC for variable " << variable <<  endl;
	
	//	observable for pdf
	RooRealVar *mvamet =  new RooRealVar(variable, variable, 0, 1000);

	//create exponential function
	RooRealVar *a =  new RooRealVar("a","a", 0.003 ,-0.100, 0.100);
	RooExponential *exp_fit = new RooExponential("Exponential", "Exponential fit to mvamet", *mvamet , *a);

	//create gaussian
	RooRealVar *mean = new RooRealVar("mean", "mean of gaussian", 180, -300, 300);
	RooRealVar *sigma = new RooRealVar("sigma", "width of gaussian", 200, -1000, 1000);
	//If you want to fix a parameter in the fit, you just specify that as a property of the RooRewalVar parameter object
	mean->setConstant(kTRUE);
	RooGaussian *gauss = new RooGaussian("gauss", "gaussian PDF", *mvamet, *mean, *sigma);
	
	//Import the histrogram with the required variable
	TH1F *hh = histo ;

	//make histogram for Gaussian fit
	RooDataHist *data_gauss = (RooDataHist*) new RooDataHist("data", "dataset with mvamet", *mvamet , hh);
	//make histgoram for exponential fit
	TH1F *hexp = (TH1F*) hh->Clone("hexp");
	RooDataHist *data_exp = (RooDataHist*) new RooDataHist("data", "dataset with mvamet", *mvamet , hexp);
	
	//Define main signal region in data as [175, 1000]
	mvamet->setRange("signal", 175, 700);

	// Fit gaussian pdf only to data in the "signal" region
	//RooFitResult *r_sig_gauss = gauss->fitTo(*data_gauss, Save(kTRUE), Range("signal"));
	// Fit exponential pdf to data in the "signal" region 
	//RooFitResult *r_sig_exp = exp_fit->fitTo(*data_exp, Save(kTRUE), Range("signal"));

	//Fit Gaussian to all data
	//RooFitResult *r_all_gauss = gauss->fitTo(*data_gauss, Save(kTRUE));
	//Fit exp to all data
	//RooFitResult *r_all_exp = exp_fit->fitTo(*data_exp, Save(kTRUE));

	//---------------------------------
	//fit gaussian model to data
	gauss->fitTo(*data_gauss, Range(170, 500));
	//fit exponential model to data
	exp_fit->fitTo(*data_exp, Range("signal"));
	//---------------------------------------


	//PLOT / Print results
	//----------------------------

	//make plost fame in mvamet and add data and fitted model
	// RooPlot *frame = mvamet.frame(Title("Fitting a subrange Exp"));
	// data_exp->plotOn(frame);
	// //exp_fit->plotOn(frame, Range("Full"), LineStyle(kDashed), LineColor(kRed));//Add shape in full range dashed
	// exp_fit->plotOn(frame, Range("signal"));//by default only fitted range is shown

	// //Print fit results
	// cout<<"Result of fit on all data"<< endl;
	// //r_all_exp->Print();
	// cout<< "Result of fit on signal "<< endl;
	// r_sig_exp->Print();
	// frame->Draw();

	//Plot PDF and data overlaid
	RooPlot *xframe = mvamet->frame();
	data_exp->plotOn(xframe);
	data_gauss->plotOn(xframe);
	exp_fit->plotOn(xframe, LineColor(2));
	gauss->plotOn(xframe, LineColor(4));
	xframe->Draw();

	//-----------------------------------------------------------
	//What are the variables of my model?
	//Given any composite RooFit value object, the getVariables() method returns you a RooArgSet 
	//with all the parameters of your model
	//-----------------------------------------------------------
	//Get parameters for gaussian fit
	RooArgSet *params_gauss = gauss->getVariables();
	params_gauss->Print("v");
	//Get parameters for exponential fit
	RooArgSet *params_exp = exp_fit->getVariables();
	params_exp->Print("v");
}


std::vector<Double_t> Fit_Photon_Photon_Control_model(const TH1F* histo = 0, const char* variable)
{
	TCanvas *c1 = new TCanvas();
	cout <<"Start the fitting of the photon+Jet MC for variable " << variable <<  endl;
	
	Double_t nEntries = histo->GetEntries();

	//Declare observable
	RooRealVar *mvamet =  new RooRealVar(variable, variable, 0, 1000);

	//create EXPONENTIAL PDF
	RooRealVar *a =  new RooRealVar("a","a", 0.0005 ,-0.100, 0.100);
	RooExponential *exp_fit = new RooExponential("Exponential", "Exponential fit to mvamet", *mvamet , *a);

	//create GAUSSIAN PDF
	RooRealVar *mean = new RooRealVar("mean", "mean of gaussian", 206, -300, 300);
	RooRealVar *sigma = new RooRealVar("sigma", "width of gaussian", 30, -10, 10);
	//If you want to fix a parameter in the fit, you just specify that as a property of the RooRewalVar parameter object
	mean->setConstant(kTRUE);
	sigma->setConstant(kTRUE);
	RooGaussian *gauss_fit = new RooGaussian("gauss", "gaussian PDF", *mvamet, *mean, *sigma);
	
	// Import the histrogram with the required variable
	TH1F *hh = histo ;
	hh->SetName("hist_data");
	//make histogram for fit
	RooDataHist *data = (RooDataHist*) new RooDataHist("data", "dataset with mvamet", *mvamet , hh);

	// Start creating first composite pdf
	// Clone original exponential to give us an exponential to make composite pdf
	RooExponential *exp_fit_comp_1 = exp_fit->clone("exp_fit_comp_1");
	// Clone origianl gaussian to give us an gaussian to make composite pdf
	RooGaussian *gauss_fit_comp_1 = gauss_fit->clone("gauss_fit_comp_1");
	// define the first ratio to build composite pdf
	 RooRealVar *f1_1 = new RooRealVar("Gauss_Frac", "gaussian fraction", 0.8, -1, 1.);
	// RooRealVar *f1_1 = new RooRealVar("Gauss_Frac", "gaussian fraction", 0, nEntries);
	// RooRealVar *f1_2 = new RooRealVar("Exp_Frac", "exponential fraction", 0, nEntries);
	// model1 = (gauss_fit_comp_1)*(c1) + (1 - c1)*(exp_fit_comp_1)
	RooAddPdf *model1 = (RooAddPdf*) new RooAddPdf("model1", "First Model", RooArgList(*gauss_fit_comp_1,*exp_fit_comp_1),*f1_1);

	// //Fit Gaussian to all data
	// RooFitResult *r_full_gauss = gauss_fit->fitTo(*data/*, Save(kTRUE)*/);
	// //Fit exp to all data
	// RooFitResult *r_full_exp = exp_fit->fitTo(*data/*, Save(kTRUE)*/);
	// //Fit Composite PDF to all data
	// RooFitResult *r_full_composite1 = model1->fitTo(*data /*, Save(kTRUE)*/);

	//Define main signal region in data as [175, 1000]
	mvamet->setRange("signal", 180, 1000);

	// // Fit gaussian pdf only to data in the "signal" region
	// RooFitResult *r_sig_gauss = gauss_fit->fitTo(*data /*, Save(kTRUE)*/, Range("signal"));
	// // Fit exponential pdf to data in the "signal" region 
	// RooFitResult *r_sig_exp = exp_fit->fitTo(*data /*, Save(kTRUE)*/, Range("signal"));
	// // Fit Composite model to data in the "signal" region
	// RooFitResult *r_model1 = model1->fitTo(*data/*, Save(kTRUE)*/, Range("signal"));

	model1->fitTo(*data, Range("signal"));
	// PLOT / Print results
	//----------------------------

	//make plot frame in mvamet and add data and fitted model
	RooPlot *frame = mvamet.frame(Title("Fitting model for MC Photon+Jet mvamet"));
	data->plotOn(frame);
	//exp_fit->plotOn(frame, Range("Full"), LineStyle(kDashed), LineColor(kRed));//Add shape in full range dashed
	exp_fit->plotOn(frame, Range("signal"), LineStyle(2) ,LineColor(2), Name("exp_fit"));//by default only fitted range is shown
	gauss_fit->plotOn(frame, Range("signal"), LineColor(3),LineStyle(2) , Name("gauss_fit"));
	model1->plotOn(frame, Range("signal"),LineColor(4), Name("model1"));
	//model2->plotOn(frame, Range("signal"), LineStyle(3), LineColor(6), Name("model2"));
	frame->Print("v");
	frame->Draw();

	// TCanvas *c3 = new TCanvas();
	// Print_chiSquare(frame, "h_data", "model2");
	// frame->Print("v");

	//take integral over the range signal to get a normalization factor
	RooAbsReal* integral_over_signal = model1->createIntegral(*mvamet, Range("signal"));
	Double_t normalization_constant = integral_over_signal->getVal();

	//create vector to store the integrals of the bins
	std::vector<Double_t> Integral_vector;
	for (Double_t i = 0; i < 82; i++)
	{
		mvamet->setRange("bin_width", 180 + (i*10) , 190 + (i*10) );
		RooAbsReal* pdf_integral_width = model1->createIntegral(*mvamet, /* NormSet(*mvamet),*/ Range("bin_width"));
		Double_t correction_factor = (pdf_integral_width->getVal()) / normalization_constant ;
		Double_t scaled_correction_factor = correction_factor;
		//Double_t factor = (pdf_integral_width->getVal())*(nEntries);
		//calculate integral in this range
		//cout << "Normalized Integral between " << 180 + (i*10) << " and " << 190 + (i*10) << "           = " << pdf_integral_width->getVal() << endl;
		//cout << "Corrected Integral between " << 180 + (i*10) << " and " << 190 + (i*10) << " = " << factor << endl;
		
		//add this value to the vector storing the integrals
		Integral_vector.push_back(scaled_correction_factor);
	}
return Integral_vector;
}

std::vector<Double_t> Fit_Data_Photon_Control_model(const TH1F* histo = 0, const char* variable)
{
	TCanvas *c1 = new TCanvas();
	cout <<"Start the fitting of the photon+Jet Data for variable " << variable <<  endl;
	
	Double_t nEntries = histo->GetEntries();

	//Declare observable
	RooRealVar *mvamet =  new RooRealVar(variable, variable, 0, 1000);

	//create EXPONENTIAL PDF
	RooRealVar *a =  new RooRealVar("a","a", 0.003 ,-0.100, 0.100);
	RooExponential *exp_fit = new RooExponential("Exponential", "Exponential fit to mvamet", *mvamet , *a);

	//create GAUSSIAN PDF
	RooRealVar *mean = new RooRealVar("mean", "mean of gaussian", 193, -300, 300);
	RooRealVar *sigma = new RooRealVar("sigma", "width of gaussian", 20, -10, 10);
	//If you want to fix a parameter in the fit, you just specify that as a property of the RooRewalVar parameter object
	mean->setConstant(kTRUE);
	sigma->setConstant(kTRUE);
	RooGaussian *gauss_fit = new RooGaussian("gauss", "gaussian PDF", *mvamet, *mean, *sigma);
	
	// Import the histrogram with the required variable
	// TH1F *hh = (TH1F*) new TH1F("hh", "hh")
	TH1F *hh = histo ;
	hh->SetName("hist_data");
	//make histogram for fit
	RooDataHist *data = (RooDataHist*) new RooDataHist("data", "dataset with mvamet", *mvamet , hh);

	// Start creating first composite pdf
	// Clone original exponential to give us an exponential to make composite pdf
	RooExponential *exp_fit_comp_1 = exp_fit->clone("exp_fit_comp_1");
	// Clone origianl gaussian to give us an gaussian to make composite pdf
	RooGaussian *gauss_fit_comp_1 = gauss_fit->clone("gauss_fit_comp_1");
	// define the first ratio to build composite pdf
	RooRealVar *f1_1 = new RooRealVar("Gauss_Frac", "gaussian fraction", 0.8, -3, 1.);
	// RooRealVar *f1_1 = new RooRealVar("Gauss_Frac", "gaussian fraction", 0, nEntries);
	// RooRealVar *f1_2 = new RooRealVar("Exp_Frac", "exponential fraction", 0, nEntries);
	// model1 = (gauss_fit_comp_1)*(c1) + (1 - c1)*(exp_fit_comp_1)
	RooAddPdf *model1 = (RooAddPdf*) new RooAddPdf("model1", "First Model", RooArgList(*gauss_fit_comp_1,*exp_fit_comp_1), *f1_1);

	// //Fit Gaussian to all data
	// RooFitResult *r_full_gauss = gauss_fit->fitTo(*data /*, Save(kTRUE)*/);
	// //Fit exp to all data
	// RooFitResult *r_full_exp = exp_fit->fitTo(*data/*, Save(kTRUE)*/);
	// //Fit Composite PDF to all data
	// RooFitResult *r_full_composite1 = model1->fitTo(*data/*, Save(kTRUE)*/);

	// Define main signal region in data as [175, 1000]
	mvamet->setRange("signal", 180, 1000);

	// // Fit gaussian pdf only to data in the "signal" region
	// RooFitResult *r_sig_gauss = gauss_fit->fitTo(*data /*, Save(kTRUE)*/, Range("signal"));
	// // Fit exponential pdf to data in the "signal" region 
	// RooFitResult *r_sig_exp = exp_fit->fitTo(*data/*, Save(kTRUE)*/, Range("signal"));
	// // Fit Composite model to data in the "signal" region
	// RooFitResult *r_model1 = model1->fitTo(*data/*, Save(kTRUE)*/, Range("signal"));
	gauss_fit->fitTo(*data, Range("signal"));
	model1->fitTo(*data, Range("signal"));
	// PLOT / Print results
	//----------------------------

	//make plot frame in mvamet and add data and fitted model
	RooPlot *frame = mvamet.frame(Title("Fitting model for Data Photon+Jet mvamet"));
	data->plotOn(frame);
	//exp_fit->plotOn(frame, Range("Full"), LineStyle(kDashed), LineColor(kRed));//Add shape in full range dashed
	exp_fit->plotOn(frame, Range("signal"), LineStyle(2) ,LineColor(2), Name("exp_fit"));//by default only fitted range is shown
	gauss_fit->plotOn(frame, Range("signal"), LineColor(3),LineStyle(2) , Name("gauss_fit"));
	model1->plotOn(frame, Range("signal"),LineColor(4), Name("model1"));
	//model2->plotOn(frame, Range("signal"), LineStyle(3), LineColor(6), Name("model2"));
	frame->Print("v");
	// //Print fit results
	//======================
	// cout<<"Result of fit on all data"<< endl;
	// r_full_exp->Print();
	//cout<< "Result of exponential 1 fit on signal "<< endl;
	//r_sig_exp->Print();
	frame->Draw();
	TCanvas *c3 = new TCanvas();
	//Calculate ChiSquared
	// Double_t CS_Composite_2 = frame->chiSquare("model2", "h_data", 0);
	// cout <<	model2->GetName() << endl;
	// cout << "ChiSquare " <<  frame->chiSquare("exp_fit", "h_data", 0) << endl;
	// cout << "ChiSquare " <<  frame->chiSquare("gauss_fit", "h_data", 0) << endl;
	// cout << "ChiSquare " <<  frame->chiSquare("model1", "h_data", 0) << endl;
	// cout << "ChiSquare of 2nd Composite model : " << CS_Composite_2 << endl;
	
	//Print_chiSquare(frame, "h_data", "gauss_fit");
	//Print_chiSquare(frame, "h_data", "exp_fit");	
	//Print_chiSquare(frame, "h_data", "model1");
	Print_chiSquare(frame, "h_data", "model1");
	frame->Print("v");

	RooAbsReal* integral_over_signal = model1->createIntegral(*mvamet, Range("signal"));
	Double_t normalization_constant = integral_over_signal->getVal();

	//create vector to store the integrals of the bins
	std::vector<Double_t> Integral_vector;
	//loop though each bin width and get the integral of the pdf
	for (Double_t i = 0; i < 82; i++)
	{
		mvamet->setRange("bin_width", 180 + (i*10) , 190 + (i*10) );
		//RooAbsReal* Mc_pdf_integral_width = MC_PDF->createIntegral(*mvamet , NormSet(*mvamet), Range("bin_width"));
		RooAbsReal* pdf_integral_width = model1->createIntegral(*mvamet/*, NormSet(*mvamet)*/, Range("bin_width"));
		Double_t correction_factor = (pdf_integral_width->getVal()) / normalization_constant;
		Double_t scaled_correction_factor = correction_factor ;
		//calculate ratio of the pdf integrals in this range
		//cout << "Normalized Integral between " << 180 + (i*10) << " and " << 190 + (i*10) << "           = " << pdf_integral_width->getVal() << endl;
		//cout << "Corrected Integral between " << 180 + (i*10) << " and " << 190 + (i*10) << " = " << factor << endl;
		//add the value of the integral to the vector storing the integral values
		Integral_vector.push_back(scaled_correction_factor);
		//Integral_vector.push_back(gaussian_integral->getVal());
	}
		//RooAbsReal* pdf_integral = model1->createIntegral(*mvamet, NormSet(*mvamet));
		//cout << "Integral = " << pdf_integral->getVal() << endl ;
		cout << "nEntries = "<< nEntries << endl;

return Integral_vector;
}

//calculate the ratio of the values of two vectors and return a vector
std::vector<Double_t> Vector_Ratio_v1_divide_v2(std::vector<Double_t> v1, std::vector<Double_t> v2)
{
	cout << "Calculating the ratio of the integrals of the pdfs" << endl;
	//create vector to store the ratio of the vectors
	std::vector<Double_t> Vector_Ratio_v1_divide_v2;
	//assuming both vectors have the same size
	size_t size_vector = v1.size();
	for (size_t i = 0; i < size_vector; i++)
	{
		//get ratio of vectors at position i
		Double_t ratio = v1.at(i) / v2.at(i) ;
		//cout << "The ratio between "<< 180 + (i*10) << " and " << 190 + (i*10) << " = " << ratio << endl;
		Vector_Ratio_v1_divide_v2.push_back( ratio);
	}
return Vector_Ratio_v1_divide_v2;
}

//this is used for the second stage of correction
TH1F* Cal_Pt_Ratio_Hist(TH1F *hZ, TH1F *hPhoton)
{
	cout << "Getting the histogram of the Z and Photon Pt ratio"<< endl;
	//declare new histogram to display the bin by bin ratio
	TH1F *hRatioPt = new TH1F("hRatioPt", "Z Pt / Photon Pt", 100, 0, 1000);
	// TH1F *hZPt_1 = (TH1F*) hZ->Clone("hZPt_1");
	// TH1F *hPho1Pt_1 = (TH1F*) hPhoton->Clone("hPho1Pt_1");
	hRatioPt->Divide( hZ, hPhoton,1,1);
	Draw_1_Hist(hRatioPt,"Ratio of Z Pt to Photon Pt", "pt", "Pt Ratio", 8, " ");
return hRatioPt;
}

//used for second stage of correction
Double_t Second_Correction_Factor(Float_t photon_Pt, TH1F *hratios )
{
	Int_t N = (photon_Pt/ 10);
	Int_t Bin_no = N+1 ;
	Double_t Pt_factor = hratios->GetBinContent(Bin_no);
return Pt_factor;
}

//apply correction to mvamet in MC using the integral ratios
Float_t First_Correction_Factor(Double_t mvamet, std::vector<Double_t> ratio_v) 
{
	Float_t ratio;
	if ((mvamet > 180) && ( mvamet < 1000))
	{
		//find what bin range the mvamet would fall under in the histogram
		size_t N = (mvamet/10);
		//what  position on the vector does this correspond to
		size_t position = N- 18 ;
		//see what value this corresponds to in the ratio of integrals
		ratio = ratio_v.at(position);
	}
	else
	{
		ratio = 0;
	}
return ratio;
}

// //fsecond correction
// // uses the ratio of the Z pt and Photon Pt
// Double_t Second_Correction_Factor( Float_t genVpt, TH1F* Pt_Ratio)
// {
// 	//take the value of the photon PT and look for the corresponding value in the ratio histogram
// 	Double_t factor = Pt_Factor(genVpt, Pt_Ratio);
// return factor;
// }

//look at the shape change only
TH1F* Integral_Correction(TH1F* h_mvamet_before, std::vector<Double_t> Correction_ratios)
{
	//loop through the file to get the photon pt and then reweight the mvamet
	TFile *inclusive = new TFile("/home/bbachu/cms/cmssw/032/CMSSW_5_3_14/src/MitMonoJet/macros/Fits/inclusive.root");
	TTree *tree = new TTree();
	tree = (TTree*) inclusive->FindObjectAny( "Photon_photon_control");
	Float_t mvamet;
	Float_t genVpt;
	Int_t event;
	tree->SetBranchAddress("genVpt" , &genVpt);
	tree->SetBranchAddress("mvamet" , &mvamet);
	//create a hist to look at how the pdf ratios change the shape the mvamet
	TH1F *hintermediate = new TH1F("mvamet_temp", "Corrected using pdf ratios", 100, 0, 1000);
	//read all entries and fill the hist
	Long64_t nentries = (Int_t) tree->GetEntries();
	cout << "Starting the event by event correction for shape distribution only" << endl;
	for (Long64_t entry = 0; entry<nentries; entry++ )
	{
		tree->GetEntry(entry);
		//correct for difference in data and MC mvamet
		Double_t Correction_1 = First_Correction_Factor(mvamet, Correction_ratios);
		hintermediate->Fill(mvamet, Correction_1);
	}
	TCanvas *c1 = new TCanvas();
	Draw_2_Hist(h_mvamet_before, hintermediate , "Shape Correction","mvamet" , " ", "Data: Photon+Jet","MC: Corrected Photon+Jet" , 2, 4, "Data");
	return hintermediate;
}

// Combine correction 1 and 2 doing an event by event reweighting
TH1F* Correct(TH1F* Z_pt_Photon_Pt_Ratio_hist, std::vector<Double_t> Correction_ratios)
{
	//start looping over the events
	//loop through the file to get the photon pt and then reweight the mvamet
	TFile *inclusive = new TFile("/home/bbachu/cms/cmssw/032/CMSSW_5_3_14/src/MitMonoJet/macros/Fits/inclusive.root");
	TTree *tree = new TTree();
	tree = (TTree*) inclusive->FindObjectAny( "Photon_photon_control");
	Float_t mvamet;
	Float_t genVpt;
	Int_t event;
	tree->SetBranchAddress("genVpt" , &genVpt);
	tree->SetBranchAddress("mvamet" , &mvamet);
	//create hist to store the corrected final met
	TH1F *hfinal_MET = new TH1F( "mvamet",  "Data Driven met distribution", 100, 0, 1000);
	//create a hist to look at how the pdf ratios scale the mvamet
	TH1F *hintermediate = new TH1F("mvamet_temp", "Corrected using pdf ratios", 100, 0, 1000);
	//read all entries and fill the hist
	Long64_t nentries = (Int_t) tree->GetEntries();
	cout << "Starting the event by event correction" << endl;
	for (Long64_t entry = 0; entry<nentries; entry++ )
	{
		tree->GetEntry(entry);
		//correct for difference in data and MC mvamet
		Double_t Correction_1 = First_Correction_Factor(mvamet, Correction_ratios);
		//correct for difference in Z pt and Photon pt
		Double_t Correction_2 = Second_Correction_Factor(genVpt, Z_pt_Photon_Pt_Ratio_hist);
		//fill the final historgam
		//hintermediate->Fill(mvamet, Correction_1);
		hfinal_MET->Fill(mvamet, Correction_1*Correction_2);
	}
	TCanvas *c1 = new TCanvas();
	//draw this hist to look at effect of first correction
	//Draw_1_Hist(hintermediate,"Effect of pdf ratio on mvamet", "mvamet", "Reweighted photon mvamet", 8, " ");
	//draw the 2 hist to look at effect of second correction on first
	TCanvas *c2 = new TCanvas();
	//Draw_2_Hist(hintermediate, hfinal_MET, "Effect of second correcton on scaled MC", "mvamet", "", "After pdf correction", "After pt correction", 8, 4, "");
return hfinal_MET;
}	

TH1F* Correct_Pt_Only(TH1F* Z_pt_Photon_Pt_Ratio_hist)
{
	//start looping over the events
	//loop through the file to get the photon pt and then reweight the mvamet
	TFile *inclusive = new TFile("/home/bbachu/cms/cmssw/032/CMSSW_5_3_14/src/MitMonoJet/macros/Fits/inclusive.root");
	TTree *tree = new TTree();
	tree = (TTree*) inclusive->FindObjectAny( "Photon_photon_control");
	Float_t mvamet;
	Float_t genVpt;
	Int_t event;
	tree->SetBranchAddress("genVpt" , &genVpt);
	tree->SetBranchAddress("mvamet" , &mvamet);
	//create hist to store the corrected final met
	TH1F *hfinal_MET = new TH1F( "mvamet",  "Data Driven met distribution", 100, 0, 1000);
	//read all entries and fill the hist
	Long64_t nentries = (Int_t) tree->GetEntries();
	cout << "Starting the event by event correction for Pt only in UNCORRECTED MC Photon+Jet" << endl;
	for (Long64_t entry = 0; entry<nentries; entry++ )
	{
		tree->GetEntry(entry);
		//correct for difference in Z pt and Photon pt
		Double_t Correction_2 = Second_Correction_Factor(genVpt, Z_pt_Photon_Pt_Ratio_hist);
		//fill the final historgam
		hfinal_MET->Fill(mvamet, Correction_2);
	}
return hfinal_MET;
}	

//Main Function
void Fits_GJets()
{
	cout<<"Start main function"<< endl;

	cout <<"Get the required histograms"<< endl;
	//Make histograms of the main variables 
	TH1F *hPhoton_Photon_Control_mvamet = (TH1F*) Make_Hist("Photon_photon_control", "mvamet");
	TH1F *hPhoton_Photon_Control_genVpt = (TH1F*) Make_Hist("Photon_photon_control", "genVpt");
	TH1F *hData_Photon_Control_mvamet = (TH1F*) Make_Hist("data_photon_control", "mvamet");
	TH1F *hZnunu_Signal_genVpt = (TH1F*) Make_Hist("Znunu_signal", "genVpt");
	TH1F *hZnunu_Signal_mvamet = (TH1F*) Make_Hist("Znunu_signal", "mvamet");

	//look at how the mvamet distributions are
	//compare the photon+jet data to photon+jet monte carlo for mvamet
	Draw_2_Hist(hData_Photon_Control_mvamet,hPhoton_Photon_Control_mvamet , "Comparison of Photon Data and MC mvamet", "mvamet", "", "Data: Photon+Jet", "MC: Photon+Jet", 2, 4,"Data" );
	//compare the Znunu montecarlo to the photon + Jet montecarlo for mvamet
	Draw_2_Hist(hZnunu_Signal_mvamet,hPhoton_Photon_Control_mvamet, "Comparison of Z->nunu and Photon+Jet mvamet", "mvamet", "", "Z->nunu", "MC: Photon+Jet", 3, 4, "");
	//lets also scale the mvamet distributions to show the difference in the shape
	TH1F* hZnunu_Signal_mvamet_scaled_to_hPhoton_Photon_Control_mvamet = (TH1F*) Scale_h1_to_h2(hZnunu_Signal_mvamet,hPhoton_Photon_Control_mvamet);
	Draw_2_Hist(hZnunu_Signal_mvamet_scaled_to_hPhoton_Photon_Control_mvamet, hPhoton_Photon_Control_mvamet, "Comparison with same Integrals", "mvamet", "", "Z->nunu", "MC: Photon+Jet", 3, 4, "");
	//look at the difference in distributions of the Zpt and Photon Pt
	Draw_2_Hist(hPhoton_Photon_Control_genVpt, hZnunu_Signal_genVpt, "Comparison of the Z and Photon Pt", "Pt", " ", "Photon Pt", "Z Pt", 4, 3, "");

	//make the histogram for comparison of the Z pt to Photon Pt ratio (shape and normalization)
	TH1F *h_Z_Photon_Pt_Ratio = (TH1F*) Cal_Pt_Ratio_Hist(hZnunu_Signal_genVpt, hPhoton_Photon_Control_genVpt);

	//scale photon pt to the Z pt
	TH1F *hPhoton_Photon_Control_genVpt_scaled = (TH1F*) Scale_h1_to_h2(hPhoton_Photon_Control_genVpt,  hZnunu_Signal_genVpt);
	//USE A DIFFERENT DRAW FUNCTION TO GET THE LOWER PIECES
	Draw_2_Hist(hZnunu_Signal_genVpt, hPhoton_Photon_Control_genVpt_scaled, "Z pt scaled to Photon Pt", "Pt", "", "Z Pt", "Photon Pt", 8, 4, "");
	//get the ratio histogram for the scaled pt (shape only)
	TH1F *h_Z_Photon_Pt_Ratio_scaled = (TH1F*) Cal_Pt_Ratio_Hist(hZnunu_Signal_genVpt,   hPhoton_Photon_Control_genVpt_scaled);
	Draw_1_Hist(h_Z_Photon_Pt_Ratio_scaled,"Ratio of scaled Pt", "Pt", " Scaled ratio", 6, "");
	
	//make the vectors storing the integrals
	std::vector<Double_t> Photon_photon_control_integral_vector = Fit_Photon_Photon_Control_model(hPhoton_Photon_Control_mvamet , "mvamet");
	std::vector<Double_t> Data_photon_control_integral_vector = Fit_Data_Photon_Control_model(hData_Photon_Control_mvamet , "mvamet");
	//make the vector storing the ratios
	std::vector<Double_t> Correction_ratios = Vector_Ratio_v1_divide_v2(Data_photon_control_integral_vector,Photon_photon_control_integral_vector);

	//make historgram showing the correction to the mvamet by looking at the pdf correction only
	TH1F *h_integral_correction = (TH1F*) Integral_Correction(hData_Photon_Control_mvamet, Correction_ratios);
	//compare the shape correction of the new MC to the old MC
	Draw_2_Hist(hPhoton_Photon_Control_mvamet, h_integral_correction, "Shape Correction of MC Photon+Jet", "mvamet", " " , "Before", "After", 4, 46, "");


	TH1F *Data_Driven_mvamet = (TH1F*) Correct(h_Z_Photon_Pt_Ratio, Correction_ratios);

	//draw final histogram
	//Draw_1_Hist(Data_Driven_mvamet,"Data Driven Met", "met", "Corrected Photon Met", 4, "");

	//make a comparison of this to the montecarlo Znunu met 
	Draw_2_Hist(hZnunu_Signal_mvamet, Data_Driven_mvamet, "Met Distributions", "mvamet", "" , "MC: Znunu", "Estimation", 3, 7," ");
	//Draw_2_Hist(TH1F* h_1, TH1F* h_2, const Char_t *Name, const Char_t *XAxis, const Char_t *YAxis, const Char_t *label1, const Char_t * label2, Color_t lcolor1, Color_t lcolor2, const Char_t* h1source)

	//output the vectors to take a look at the integrals
//	Output_Vector_Double_t("The value of data integral at ", Data_photon_control_integral_vector);
//	Output_Vector_Double_t("The value of montecarlo integral at ", Photon_photon_control_integral_vector);

	//cout << "Sum of  Data_photon_control_integral_vector = " << Sum_vector(Data_photon_control_integral_vector) << endl;
	//cout << "Sum of Photon_photon_control_integral_vector = " << Sum_vector(Photon_photon_control_integral_vector) << endl;
	//Output_Vector_Double_t(" The Correction Ratio at ", Correction_ratios);

	//apply the Z Pt, Photon Pt correction to the (UNCORRECTED) MC photon plus jet
	TH1F* h_MC_PT_correction_only = (TH1F*) Correct_Pt_Only(h_Z_Photon_Pt_Ratio);
	//Compare correction from pt only to the Znunu MC
	Draw_2_Hist(h_MC_PT_correction_only, hZnunu_Signal_mvamet, "Using only Pt Ratios", "mvamet","", "Estimate", "MC: Z->nunu", 40, 3, "");


	cout<<"end main function"<< endl ;
}
