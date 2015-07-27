// #ifndef __CINT__
// #include "RooGlobalFunc.h"
// #else
// // Refer to a class implemented in libRooFit to force its loading
// // via the autoloader.
// class Roo2DKeysPdf;
// #endif

TInterpreter *g = new TInterpreter();
g->AddIncludePath('/cvmfs/cms.cern.ch/'+ gSystem.Getenv("SCRAM_ARCH")+'/lcg/roofit/5.34.07-cms/include');

#include <iostream>
#include <iomanip>
#include <math.h>
#include <TH1.h>
#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TClonesArray.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TImage.h>
#include <TROOT.h>
#include <TCut.h>
#include <vector>
#include <TLegend.h>
#include <TMath.h>
#include <TBranch.h>
#include <THStack.h>
#include "TH2.h"
#include "TNtuple.h"
#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TFrame.h"
#include "TInterpreter.h"
#include "TVirtualHistPainter.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooRealVar.h"
#include "RooPlot.h"
#include " RooGaussian.h"
using namespace RooFit;

/* void Slide 8()//In RooFit all objects are self documented
{
	//objects representing a 'real' value
	RooRealVar mass("mass", "Invariant mass", 5.2, 5.3); //initial range
	RooRealVar width("width", "BO mass width", 0.00027, "GeV");//initial value, optional unit
	RooRealVar mb0("mb0", "B0 mass", 5.2794. "Gev");//initial value, optional unit
	
	//PDF Object
	RooGaussian b0sig("b0sig", "B0 sig PDF", mass ,mb0, width); //reference to variables
} */

/*  void Gauss_PDF()
{
	TCanvas *c1 = new TCanvas();
	//Build Gaussian PDF
	RooRealVar x("x","x", -10, 10);
	RooRealVar mean("mean", "mean of gaussian", 0, -10, 10);
	RooRealVar sigma("sigma", "width of gaussian", 3);
	
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
}  */



//Data in general comes in two flavours:
// 1) Unbinned data: represented in ROOT by class TTree and,
// 2) Binned data represented in ROOT by TH1, TH2, ...

/* void Binned_Data( TH1D *h1) //In RooFit, binned data is represented by the class RooDataHist
{
	TH1D* h2 = gDirectory->Get("");
	cout  << "Beginning the import of data"<< endl;
	TCanvas *c2 = new TCanvas();
	cout << "Create new tcanvas" << endl;
	//You can import the contents of any ROOT histogram in a RooDataHist object
	RooRealVar z("z", "z", -10, 10);
	cout <<"real var created"<< endl;
	RooDataHist data("data", "dataset with z", z, h2, 1.0);
	cout <<"defined var and datahist" << endl;
	//A RooDataHist always associates the histogram with a RooFit variable object of type RooRealVar
	//In this way it always knows what kind of data is stored in the histogram
	
	//A RooDataHist can be visualized in the same way as a function can be visualized
	RooPlot *xframe2 = z.frame();
	cout << "roo plot "<, endl;
	data.plotOn(xframe2);
	frame->Draw();
	TImage *img = TImage::Create();
	img->FromPad(c2);
	img->WriteImage("canvas2.png");
} */


TH1D* TH1_data_gauss()
{
	cout << "Starting to create random binned data" << endl;
	TCanvas *c = new TCanvas();
	cout << "creating histrogram" << endl;
	TH1D *h = new TH1D("gaus", "gaus", 100, -10, 10);
	cout << "Filling histogram" << endl;
	h->FillRandom("gaus",100);
	h->Draw();
	
	cout << "Processing events" << endl;
	gSystem->ProcessEvents();
	TImage *img = TImage::Create();
	
	img->FromPad(c);
	img->WriteImage("canvas.png");
	return h;
	
}

void example(const TH1* histo = 0)
{
	//Build Gaussian PDF
	RooRealVar x("x", "x", -10, 10);
	RooRealVar mean("mean", "mean of Gaussian", 0, -100 , 100);
	RooRealVar sigma("sigma", "width of Gaussian", 3.0, 0. , 10.);
	RooGaussian gauss("gauss","gaussian PDF" , x, mean, sigma);
	
	RooAbsData* data = 0;
	if (histo)
	{
		//If a histogram is given, import it into a RooDataHist - Binned dATA
		data = new RooDataHist("data", "data", x, histo);
	}
	else
	{
		//If no histogram is given, generate some toy data - Unbinned data
		data = gauss.generate(x,10000);
	}
	
	//Fit model to data//Note here that fiTo accepts both binned and unbinned data
	gauss.fitTo(*data);
	
	//Plot PDF and toy data overlaid
	RooPlot* xframe = x.frame();
	data->plotOn(xframe);
	gauss.plotOn(xframe);
	xframe->Draw();
	
	//Print final values of the parameters
	mean.Print();
	sigma.Print();
	
	//Delete the data
	delete data;
	
}

void example_my()
{
	cout<<"Start"<<endl;
	
	//Gauss_PDF();
	
	TCanvas *c1 = new TCanvas();
	//Build Gaussian PDF
	RooRealVar x("x","x", -10, 10);
	RooRealVar mean("mean", "mean of gaussian", 0, -10, 10);
	RooRealVar sigma("sigma", "width of gaussian", 3);
	
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
	
//	TH1D *h1 = TH1_data_gauss();

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
//	TH1* hh = (TH1*) gDirectory->Get("ahisto");
//	cout  << "Beginning the import of data"<< endl;
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

void Fit1()
{
	example();
}