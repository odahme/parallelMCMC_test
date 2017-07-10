
#this script ist just for testing the parallalized MCMC



#simulate gaus distrubuted points
mean=3
sigma=1.5
np=1000
./bin/simulategauspoints $mean $sigma $np

#running a MCMC to fit these points
inputfilename="/home/oliver/parallelMCMC_test/pointStorage/simulatedgauspoints/gaus.txt"
ncpu=1

./bin/createpoint $inputfilename $ncpu
