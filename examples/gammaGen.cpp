#include "AmpGen/Particle.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/SumPDF.h"
#include "AmpGen/FitResult.h"
#include "AmpGen/Minimiser.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/Utilities.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/Generator.h"
#include "AmpGen/PolarisedSum.h"
#include "AmpGen/AddCPConjugate.h"
#include "AmpGen/QMI.h"
#ifdef _OPENMP
  #include <omp.h>
#endif
#if ENABLE_AVX
  #include "AmpGen/EventListSIMD.h"
  using EventList_type = AmpGen::EventListSIMD;
#else
  #include "AmpGen/EventList.h"
  using EventList_type = AmpGen::EventList; 
#endif

#include "TRandom3.h"
#include "TFile.h"
#include "TTree.h"



using namespace AmpGen;


real_t testSum(EventList_type& mcSig){
  real_t a =0;

  #pragma omp parallel for reduction(+:a)
  for(unsigned i=0;i<mcSig.size();++i) a+=1;
  return a; 
}





int main(int argc , char* argv[] ){
  OptionsParser::setArgs( argc, argv );
  MinuitParameterSet MPS;
  MPS.loadFromStream();

  TRandom3 rndm(0);

  EventType sigType(NamedParameter<std::string>("EventType", "", "Signal Type to generate"));
  size_t NInt(NamedParameter<size_t>("NInt", 1e7, "Number of events to calculate normalisation - should be large"));
  size_t seed(NamedParameter<size_t>("Seed", 0, "Random seed for generation"));
  size_t blockSize(NamedParameter<size_t>("blockSize", NInt, "Blocksize for generation"));
  size_t nEvents(NamedParameter<size_t>("nEvents", 1000, "number of events to generate, multiplies by the BR of each tag"));
  size_t nBins(NamedParameter<size_t>("nBins", 100, "Number of bins for projection histograms"));
  std::string outputFile(NamedParameter<std::string>("Output", "Generate_Output.root", "output file"));
  auto tags = NamedParameter<std::string>("BTagTypes").getVector();
  const size_t      nThreads = NamedParameter<size_t>     ("nCores"    , 8           , "Number of threads to use" );
  #ifdef _OPENMP
    omp_set_num_threads( nThreads );
    INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
    omp_set_dynamic( 0 );
  #endif

  auto  xy = QMI::xy(sigType, MPS); 
  EventList_type mc =  Generator<>(sigType, &rndm).generate(NInt);
  //real_t alpha(testSum(mcSig));
  //INFO("A = "<<alpha);
  //return 0;
  ProfileClock t_compilePoly, t_calcPoly, t_calcPolyManual;
  t_compilePoly.start();
  std::map<std::string, CompiledExpression<real_t(const real_t*, const real_t*)> > compiledPC(QMI::cPhaseCorrection(sigType, MPS));
  t_compilePoly.stop();
  INFO("Took "<<t_compilePoly.t_duration<<"ms to compile "<<NamedParameter<size_t>("PhaseCorrection::Order", 4)<<" polynomial");

  for (auto& p : compiledPC){
    INFO("Have "<<p.first<<"(0) = "<<p.second(mc[0].address())<<" x = "<<xy[0](mc[0].address())<<" y = "<<xy[1](mc[0].address()) ); 
  }

  MinuitParameterSet MPS_Kspipi;
  MPS_Kspipi.loadFromFile("Kspipi.opt");
  AddCPConjugate(MPS_Kspipi);
  CoherentSum A(sigType, MPS_Kspipi);

  CoherentSum Abar(sigType.conj(true), MPS_Kspipi);

  A.setEvents(mc);
  A.setMC(mc);
  A.prepare();
  Abar.setMC(mc);
  Abar.prepare();

  INFO("A(0) = "<<A.getValNoCache(mc[0]));
  INFO("A(0) = "<<A.getVal(mc[0]));
  INFO("NA = "<<A.norm());
  std::map<std::string, std::vector<real_t> > amps(QMI::AmpArrays(mc, A, Abar));
  TFile * tOutFile = TFile::Open(outputFile.c_str(), "RECREATE");
  tOutFile->cd();
  real_t C(QMI::cosTerm(amps, MPS, compiledPC, mc));
  real_t S(QMI::sinTerm(amps, MPS, compiledPC, mc));
  auto my_A_real = [&A](Event& evt){ return std::real(A.getValNoCache(evt)); };
  auto my_A_imag = [&A](Event& evt){ return std::imag(A.getValNoCache(evt)); };
  auto my_Abar_real = [&Abar](Event& evt){ return std::real(Abar.getValNoCache(evt)); };
  auto my_Abar_imag = [&Abar](Event& evt){ return std::imag(Abar.getValNoCache(evt)); };



  for (auto tag : tags){
    INFO("tag = "<<tag);
    auto a = split(tag, ' ');
    std::string BName(a[0]);
    EventList list(sigType);
    size_t nEventsB(std::stod(a[2]) * nEvents);
    Int_t gammaSign(std::stoi(a[1]));
    //real_t norm(QMI::ckm_norm(A, Abar, MPS, mc, compiledPC, gammaSign, amps));
    real_t norm = 0;
    complex_t zB(QMI::ckm_zB(MPS, gammaSign));
    if (gammaSign >0) norm = Abar.norm() + A.norm() * std::norm(zB) + 2 * (std::real(zB) * C - gammaSign * std::imag(zB) * S);
    if (gammaSign <0) norm = A.norm() + Abar.norm() * std::norm(zB) + 2 * (std::real(zB) * C - gammaSign * std::imag(zB) * S);
    auto psi = [&A, &Abar, &MPS, &compiledPC, &gammaSign, &amps, &norm](Event evt){
      return QMI::probCKM_unnorm(evt, A, Abar, MPS, compiledPC, gammaSign)/norm;
    };
    QMI::generateCKM(psi, list, nEventsB, seed, blockSize);
    if (NamedParameter<bool>("boostCKM", true)){
      EventList list2(list.eventType());
      QMI::boostCKMEvents(list, list2);
      QMI::writeEvents( list,"RF"+  a[0], nBins);
      QMI::writeEvents( list2,a[0], nBins);
    }
    else{
      QMI::writeEvents( list, a[0], nBins);
    }
    if (NamedParameter<bool>("writeAmpData", true)){
      QMI::writeValues(list, my_A_real, a[0] + "_A_real");
      QMI::writeValues(list, my_A_imag, a[0] + "_A_imag");
      QMI::writeValues(list, my_Abar_real, a[0] + "_Abar_real");
      QMI::writeValues(list, my_Abar_imag,a[0] + "_Abar_imag");
    }
 
  }
  auto my_dd = [&A, &Abar](Event& evt){
    return QMI::dd(evt, A, Abar);
  };
  auto my_corr = [&compiledPC, &MPS](Event& evt){
    return QMI::phaseCorrection(evt, compiledPC, MPS);
  };

  QMI::writeValues(mc, my_dd, "dd");
  QMI::writeValues(mc, my_corr, "corr");

  QMI::writeDalitz(mc);
  tOutFile->Close();

  return 0;
}
