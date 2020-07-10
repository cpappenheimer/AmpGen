#include "AmpGen/Psi3770.h"
#include "AmpGen/pCorrelatedSum.h"
#include "AmpGen/pCoherentSum.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/CorrelatedSum.h"
#include "AmpGen/CorrelatedLL.h"
#include "AmpGen/corrEventList.h"
//#include "AmpGen/SumLL.h"
#include "AmpGen/SimPDF.h"
#include "AmpGen/SumPDF.h"
//#include "AmpGen/CombCorrLL.h"
//#include "AmpGen/CombGamCorrLL.h"
#include "AmpGen/CombLL.h"
#include "AmpGen/MetaUtils.h"
#include <typeinfo>

//#include <boost/algorithm/string.hpp>
using namespace AmpGen;
using namespace std::complex_literals;




int main( int argc, char* argv[] )
{
  /* The user specified options must be loaded at the beginning of the programme, 
     and these can be specified either at the command line or in an options file. */   
  OptionsParser::setArgs( argc, argv );
  //OptionsParser::setArgs( argc, argv, "Toy simulation for Quantum Correlated Ψ(3770) decays");
  /* */
  //auto time_wall = std::chrono::high_resolution_clock::now();
  //auto time      = std::clock();


  /* Parameters that have been parsed can be accessed anywhere in the program 
     using the NamedParameter<T> class. The name of the parameter is the first option,
     then the default value, and then the help string that will be printed if --h is specified 
     as an option. */


  bool m_debug = NamedParameter<bool>("debug", false);
  int nThreads = NamedParameter<int>("nThreads", 12);

   #ifdef _OPENMP
  omp_set_num_threads( nThreads );
  if (m_debug) INFO( "Setting " << nThreads << " fixed threads for OpenMP" );
  omp_set_dynamic( 0 );
  #endif


   auto pNames = NamedParameter<std::string>("EventType" , ""    
      , "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(); 
 
  

  auto BTags   = NamedParameter<std::string>("BTagTypes" , std::string(), "").getVector();
 

  MinuitParameterSet MPS;
  MPS.loadFromStream();

  EventType eventType = EventType(pNames);

  std::vector<EventList> SigData;
  std::vector<EventList> SigInt;
  std::vector<EventType> SigType;
  std::vector<std::string> sumFactors;
  std::vector<int> gammaSigns;
  std::vector<int> useXYs;
  for (auto& BTag : BTags){

    INFO("B DecayType = "<<BTag);
 
 
    

    auto B_Name = split(BTag,' ')[0];
    auto B_Pref = split(BTag,' ')[1];
    int B_Conj = std::stoi(split(BTag,' ')[2]);
    int gammaSign = std::stoi(split(BTag,' ')[3]);
    bool useXY = std::stoi(split(BTag,' ')[4]);
 
    
    INFO("GammaSign = "<<gammaSign);
    if (B_Conj == 1){
      eventType = eventType.conj(true);
    }

    auto sig = pCoherentSum(eventType, MPS ,B_Pref, gammaSign, useXY);



    std::string DataFile = NamedParameter<std::string>("DataSample", "");
    std::string IntFile = NamedParameter<std::string>("IntegrationSample", "");

    std::stringstream DataSS;
    DataSS<<DataFile<<":"<<B_Name;
    std::string DataLoc = DataSS.str();


    std::stringstream IntSS;
    IntSS<<IntFile<<":"<<B_Name;
    std::string IntLoc = IntSS.str();

    EventList Data = EventList(DataLoc, eventType);
    EventList Int = EventList(IntLoc, eventType);
    
    sig.setEvents(Data);
    sig.setMC(Int);


    SigData.push_back(Data);
    SigInt.push_back(Int);
    SigType.push_back(eventType);
    sumFactors.push_back(B_Pref);
    gammaSigns.push_back(gammaSign);
    useXYs.push_back(useXY);


    if (NamedParameter<bool>("FitEach", false)){
    auto ll = make_likelihood(Data, sig);
    Minimiser mini = Minimiser(ll, &MPS);


 
  std::vector<std::string> params = {"pCoherentSum::x+", "pCoherentSum::y+", "D0{K*(892)bar-{K0S0,pi-},pi+}_Re"};
  std::vector<std::string> outputs = {"xp", "yp", "K892plus"};

  for (int i=0; i<params.size();i++){

    auto m_param = MPS[params[i]];
    auto mean = m_param->mean();
    auto err = m_param->err();
    auto min = -5 * err + mean;
    auto max = 5 * err + mean;
    auto step = 0.5 * err;
    int n = (max - min)/step;
    std::stringstream ss_file;
    ss_file<<outputs[i]<<".csv";
    auto fileName = ss_file.str();
    std::ofstream outfile;
    outfile.open(fileName);
    for (int i=0; i<n; i++){
      auto pVal = i * step + min;
      m_param->setCurrentFitVal(pVal);
      auto LLVal = mini.FCN();
      outfile<<pVal<<" "<<LLVal<<"\n";
    }
    outfile.close();
    MPS[params[i]]->setCurrentFitVal(mean);

  }
 


    mini.gradientTest();


   mini.doFit();

    std::ofstream fitOut;
    fitOut.open(NamedParameter<std::string>("fitOutput", "fit.csv"));
    fitOut<<"x+ "<<MPS["pCoherentSum::x+"]->mean()<<" "<<MPS["pCoherentSum::x+"]->err()<<"\n"
          <<"y+ "<<MPS["pCoherentSum::y+"]->mean()<<" "<<MPS["pCoherentSum::y+"]->err()<<"\n";
    fitOut.close();
    }


/*
    auto xp = MPS["pCoherentSum::x+"];
    auto yp = MPS["pCoherentSum::y+"];

    double stepxp = xp->err();
    double stepyp = yp->err();

    double xpMax = xp->mean() + 5*xp->err();
    double xpMin = -xpMax;
    double ypMax = yp->mean() + 5*yp->err();
    double ypMin = -ypMax;

    int Nxp = 10;
    int Nyp = 10;
    sig.prepare();
    std::ofstream NOut;
    NOut.open("norm.csv");
    INFO("Running through "<<Nxp<<" normalisations");
    INFO("Running through "<<Nyp<<" normalisations");
    for (int i=0; i <Nxp; i++){
      for (int j=0; j<Nyp; j++){
        INFO("At "<<i<<j); 
        MPS["pCoherentSum::x+"]->setCurrentFitVal(xpMin + i*stepxp);
        INFO("x+ = "<<xpMin + i * stepxp);
        MPS["pCoherentSum::y+"]->setCurrentFitVal(ypMin + j*stepyp);
        INFO("y+ = "<<ypMin + j * stepyp);

        double slowNorm = sig.norm();

        INFO("slowNorm = "<<slowNorm);
        double fastNorm = sig.getFastNorm();
        INFO("FastNorm = "<<fastNorm);

        INFO(xpMin + i * stepxp<<" "<<ypMin + j * stepyp<<" "<<slowNorm<<" "<<fastNorm);
        NOut<<xpMin + i * stepxp<<" "<<ypMin + j * stepyp<<" "<<slowNorm<<" "<<fastNorm<<"\n";
      }
    }
    NOut.close();
*/


  }

  auto LLC = CombLL(SigData, SigInt, SigType, MPS, sumFactors, gammaSigns, useXYs);
  Minimiser mini(LLC, &MPS);
  mini.gradientTest();
  mini.doFit();
  FitResult * fr = new FitResult(mini);
  fr->writeToFile("BFit.log");

  return 0;
}


std::vector<std::string> makeBranches(EventType Type, std::string prefix){
  auto n = Type.finalStates().size();
  std::vector<std::string> branches;
  std::vector<std::string> varNames = {"PX", "PY", "PZ", "E"};
  for (long unsigned int i=0; i<n; i++){
    auto part = replaceAll(Type.finalStates()[i], "+", "p");
    part = replaceAll(part, "-", "m");
    for (auto varName : varNames){
      std::ostringstream stringStream;
      stringStream<<prefix<<part<<"_"<<varName;
      branches.push_back(stringStream.str());
    }
  }
  return branches;
}


