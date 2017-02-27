//this program is a test of the Tthread class

//general stuff
#include <fstream>
#include "string"
#include <iostream>
#include <iterator>
#include <vector>
#ifndef _WIN32
#include <unistd.h>
#endif

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
#include "TThread.h"
#include "TThreadPool.h"
#include "Riostream.h"
#include "TFrame.h"
#include "TH1F.h"


using namespace std;


//Olivers stuff
#include "MCMC/RooMinuitMCMC.hpp"


//=============================================================================
const size_t g_sleeptime = 1; // in secs.
const size_t g_multTasks = 50;
//=============================================================================

// define a custom parameters type for task objects
enum EProc {start, clean};

// a class defining task objects
class TTestTask: public TThreadPoolTaskImp<TTestTask, EProc>
{
public:
   double result;
   bool runTask(EProc /*_param*/) {
      m_tid = TThread::SelfId();
      //TThread::Sleep(g_sleeptime, 0L);
      TRandom3 *rnd = new TRandom3(m_tid);

      result = 0;
      size_t count = 10;
      for (size_t i = 0; i < count; i++) {
        result += rnd->Uniform(0,1);
      }
      return true;
   }
   unsigned long threadID() const {
      return m_tid;
   }

private:
   unsigned long m_tid;
};

//=============================================================================
void threadPool(size_t _numThreads = 10, bool _needDbg = false)
{
   cout << "ThreadPool: starting..." << endl;
   // number of tasks to process
   size_t numTasks(_numThreads * g_multTasks);

   // create a thread pool object
   // _numThreads - a number of threads in the pool
   // _needDbg - defines whether to show debug messages
   TThreadPool<TTestTask, EProc> threadPool(_numThreads, _needDbg);

   // create a container of tasks
   vector <TTestTask> tasksList(numTasks);

   cout << "ThreadPool: getting tasks..." << endl;
   cout << "ThreadPool: processing tasks..." << endl;
   // push tasks to the ThreadPool
   // tasks can be also pushed asynchronously
   for (size_t i = 0; i < numTasks; ++i) {
      threadPool.PushTask(tasksList[i], start);
   }

   // Stop thread pool.
   // The parameter "true" requests the calling thread to wait,
   // until the thread pool task queue is drained.
   threadPool.Stop(true);
   cout << "ThreadPool: done" << endl;
   for (size_t i = 0; i < tasksList.size(); i++) {
     std::cout << "result of task"<<i<<" = "<<tasksList[i].result << std::endl;
   }
}


int main(int argc, char const *argv[]) {
  for (size_t i = 0; i < 10; i++) {
    threadPool();
  }
  return 0;
}
