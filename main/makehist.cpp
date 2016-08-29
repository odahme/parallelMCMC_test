//this programm will collect MCMC points and plot them

//general stuff
#include <fstream>
#include "string"
#include "iostream"
#include "vector"

//Root stuff
#include "TCanvas.h"
#include "TH1.h"


using namespace std;


//Olivers stuff
#include "MCMC/RooMinuitMCMC.hpp"


int main(int argc,  char **argv) {

  string inputdirname = argv[1];
  size_t nofToys = atoi(argv[3]);
  size_t nofparams = atoi(argv[2]);

  TCanvas *c1 = new TCanvas("c1", "c1", 1200,1200);
  std::vector<TH1F*> histVec;
  histVec.reserve(nofparams);

  string inputintro = inputdirname;
  inputintro +="pointlist1";
  ifstream inputin(inputintro);
  for (size_t i = 1; i <= nofparams; i++) {
    string name;
    inputin >> name;
    string histname = "hist"+to_string(i);
    string histtitel = "Histogram of "+name;
    TH1F* hist = new TH1F(histname.c_str(),histtitel.c_str(),100,0,4);
    histVec.push_back(hist);
  }
  inputin.close();


  for (size_t i = 1; i <= nofToys; i++) {
    string inputfile = inputdirname;
    inputfile +="pointlist";
    inputfile +=to_string(i);
    ifstream input(inputfile);
    string intro2;
    getline(input,intro2);
    while (input) {
      for (size_t j = 0; j < nofparams; j++) {
        Double_t val;
        input >>  val;
        histVec[j]->Fill(val);
      }
    }
  }

  for (size_t i = 1; i <= nofparams; i++) {
    string outputhist = "/home/oliver/parallelMCMC_test/plots/hist"+to_string(i);
    c1->cd();
    histVec[i-1]->Draw();
    c1->SaveAs(outputhist.c_str());
  }




}
