//this programm simulates gauss points and writes them into a txt
#include "TRandom3.h"
#include <fstream>
#include "iostream"

using namespace std;

int main(int argc, char **argv) {

  //initialisation
  unsigned int np = atoi(argv[3]);
  double sigma = stod(argv[2]);
  double mu = stod(argv[1]);
  TRandom3 *rnd = new TRandom3(35);

  ofstream gausstream;
  gausstream.open("/home/oliver/parallelMCMC_test/pointStorage/simulatedgauspoints/gaus.txt");
  gausstream<<"These are"<<np<<" gaus distributed points with a mean of"<<mu<<
  "and a sigma of "<<sigma<<"\n";

  //simulation
  std::cout << "simulating "<<np<<" gaus distributed points" << std::endl;
  for (size_t i = 0; i < np; i++) {
     gausstream<<i<<"\t"<< rnd->Gaus(mu,sigma) <<"\n";
  }
  gausstream.close();
    return 0;
}
