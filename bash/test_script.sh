
#this script ist just for testing the parallalized MCMC



#simulate gaus distrubuted points
mean=3
sigma=1.5
np=100
./bin/simulategauspoints $mean $sigma $np

#running a MCMC to fit these points
inputfilename="/home/oliver/parallelMCMC_test/pointStorage/simulatedgauspoints/gaus.txt"
ncpu=4
for (( id = 1; id < ncpu+1; id++ )); do
  ./bin/createpoint $inputfilename $id &
  sleep 1
done
wait


#filling Histograms with Data and plot them
inputdirname="/home/oliver/parallelMCMC_test/pointStorage/MCMCpoints/"
toys=4
params=2
./bin/makehist $inputdirname $params $toys
