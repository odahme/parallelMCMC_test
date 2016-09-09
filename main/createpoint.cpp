//This programm will create points of a MCMC and store them in a txt file

//general stuff
#include <fstream>
#include "string"

//Root stuff
#include "RooRealVar.h"
#include "RooAbsReal.h"
#include "RooDataSet.h"
#include "RooArgSet.h"
#include "RooGaussian.h"
#include "TRandom3.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "RooPlot.h"


using namespace std;


//Olivers stuff
#include "MCMC/RooMinuitMCMC.hpp"



int main(int argc,  char **argv) {

  //initialisation
  string inputfile = argv[1];
  int nthreads = atoi(argv[2]);
  TCanvas* c1 = new TCanvas("c1","c1",1,1,1080,1900);


  //filling DataSet
  ifstream input(inputfile);
  RooRealVar x("x","x",-10,10) ;
  RooDataSet data("data", "data",RooArgSet(x));
  int index;
  double point;
  string intro;
  getline(input,intro);
  while (input >> index >> point) {
    x.setVal(point);
    data.add(RooArgSet(x));
  }

  double ubmu = 4; //upper bound for mu
  double lbmu = 1.5; //lower bound for mu
  double ubsigma = 3; //upper bound for sigma
  double lbsigma = 0.5; //lower bound for sigma

  TRandom3 *rnd = new TRandom3(13);
  double first_mu = rnd->Uniform(lbmu,ubmu);
  double first_simga = rnd->Uniform(lbsigma,ubsigma);
  RooRealVar mean("mean","mean of gaussian",first_mu,lbmu,ubmu) ;
  RooRealVar sigma("sigma","width of gaussian",first_simga,lbsigma,ubsigma) ;
  RooGaussian* gauss = new RooGaussian("gauss","gaussian PDF",x,mean,sigma) ;


    TFile* file = new TFile("parallelMCMC_test_plots.root", "recreate");
    file->cd();
    RooAbsReal* nll = gauss->createNLL(data) ;
    RooMinuitMCMC* markov = new RooMinuitMCMC(*nll);
    markov->setSeed(35);
    markov->mcmc(2000,200,"gaus",nthreads);
    std::cout << "walk finished" << std::endl;
    TMultiGraph* grmean = markov->getWalkDis("mean",kFALSE);
    TMultiGraph* grsigma = markov->getWalkDis("sigma",kFALSE);
    grmean->Write();
    grsigma->Write();
    file->Close();
    std::cout << "I am finished" << std::endl;





// Construct unbinned likelihood of model w.r.t. data

  // markov->mcmcMulti(10000,200,"gaus",processnumber);





  // TFile* file = new TFile("alphastarTest.root", "recreate");
  // file->cd();
  // for (size_t i = 1; i < 100; i++) {
  //   RooAbsReal* nll = gauss->createNLL(data) ;
  //   RooMinuitMCMC* markov = new RooMinuitMCMC(*nll);
  //   markov->setSeed(35);
  //   Double_t alpha = (double) i/100.0;
  //   string grmName = "HistMean_at_"+to_string(alpha);
  //   string grsName = "HistSigma_at_"+to_string(alpha);
  //   markov->setAlphaStar(alpha);
  //   markov->mcmc(1000,100,"gaus");
  //   TMultiGraph* grmean = markov->getWalkDis("mean",kFALSE);
  //   grmean->SetName(grmName.c_str());
  //   TMultiGraph* grsigma = markov->getWalkDis("sigma",kFALSE);
  //   grsigma->SetName(grsName.c_str());
  //   grmean->Write();
  //   grsigma->Write();
  // }
  // file->Close();



  //c1->cd();
  // TFile *walkDisFile1 = new TFile("walkDisMean,o.root");
  // TH1F* hist1 = new TH1F();
  // hist1 = (TH1F*) walkDisFile1->get("Dismean");






  //   RooAbsReal* nll = gauss->createNLL(data) ;
  //   RooMinuitMCMC* markov = new RooMinuitMCMC(*nll);
  //   markov->setSeed(35);
  // c1->cd();
  // markov->mcmcMulti(1000,100,"gaus");
  //
  // TGraph* grmean = markov->getProfile("mean",kFALSE);
  // TMultiGraph* grsigma = markov->getWalkDis("sigma",kFALSE);
  //
  // grmean->Draw();
  // c2.cd();
  // grsigma->Draw();

  // c1->SaveAs("testplot.pdf");





  return 0;
}
