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


using namespace std;


//Olivers stuff
#include "MCMC/RooMinuitMCMC.hpp"



int main(int argc,  char **argv) {

  //initialisation
  string inputfile = argv[1];
  int processnumber = atoi(argv[2]);


  string outputfile = "/home/oliver/parallelMCMC_test/pointStorage/MCMCpoints/pointlist"+to_string(processnumber);

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

  TRandom3 *rnd = new TRandom3(13*5*processnumber);
  double first_mu = rnd->Uniform(lbmu,ubmu);
  double first_simga = rnd->Uniform(lbsigma,ubsigma);
  RooRealVar mean("mean","mean of gaussian",first_mu,lbmu,ubmu) ;
  RooRealVar sigma("sigma","width of gaussian",first_simga,lbsigma,ubsigma) ;
  RooGaussian* gauss = new RooGaussian("gauss","gaussian PDF",x,mean,sigma) ;


// Construct unbinned likelihood of model w.r.t. data
  RooAbsReal* nll = gauss->createNLL(data) ;
  RooMinuitMCMC* markov = new RooMinuitMCMC(*nll);
  markov->mcmc(10000,100);
  markov->saveCandidatesAs(outputfile.c_str());
  return 0;
}
