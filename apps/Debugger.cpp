#include <TLorentzVector.h>
#include <cmath>
#include <complex>
#include <fstream>
#include <map>
#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "AmpGen/AmplitudeRules.h"
#include "AmpGen/CompiledExpression.h"
#include "AmpGen/EventList.h"
#include "AmpGen/MinuitParameter.h"
#include "AmpGen/MsgService.h"
#include "AmpGen/OptionsParser.h"
#include "AmpGen/Particle.h"
#include "AmpGen/ParticleProperties.h"
#include "AmpGen/ParticlePropertiesList.h"
#include "AmpGen/Utilities.h"
#include "TRandom3.h"
#include "AmpGen/AddCPConjugate.h"

#ifdef _OPENMP
#include <omp.h>
#endif

#include "AmpGen/EventType.h"
#include "AmpGen/CoherentSum.h"
#include "AmpGen/IncoherentSum.h"
#include "AmpGen/Generator.h"
#include "AmpGen/Kinematics.h"
#include "AmpGen/MinuitParameterSet.h"
#include "AmpGen/NamedParameter.h"
#include "AmpGen/PolarisedSum.h"


using namespace AmpGen;


template < class FCN > void printIntensity(FCN& sig, const EventList& accepted)
{
    // print(accepted[0], sig.matrixElements() , false);
    INFO("A(x) = " << sig.getValNoCache( accepted[0] ) );
    
    complex_t a = sig.getValNoCache( accepted[0] );
    real_t a_mag = (std::conj(a) * a).real();
    INFO("Intensity = " << a_mag);
    INFO(std::fixed << std::setprecision(10) << "1 / intensity = " << 1 / a_mag);
}


void invertParity( Event& event, const size_t& nParticles)
{
  for( size_t i = 0 ; i < nParticles; ++i )
  {
    event[4*i + 0] = -event[4*i+0];
    event[4*i + 1] = -event[4*i+1];
    event[4*i + 2] = -event[4*i+2];
  }
}

template < class FCN > void debug( FCN& sig, EventList& accepted){
  INFO("Debugging: ");
  unsigned eventToDebug = 0;
  sig.setEvents( accepted );
  INFO("Debugging this event:");
  accepted[eventToDebug].print();
  sig.prepare();
  sig.debug( accepted[eventToDebug] );
  printIntensity(sig, accepted);

  INFO("\n\nInverting parity");
  for( unsigned int i = 0 ; i != accepted.size(); ++i ) 
    invertParity(accepted[i], accepted.eventType().size() );
  accepted[eventToDebug].print();
  sig.reset();
  sig.setEvents(accepted);
  sig.prepare();
  sig.debug( accepted[eventToDebug] );
  printIntensity(sig, accepted);
}

int main( int argc, char** argv )
{
  OptionsParser::setArgs( argc, argv );

  int seed = NamedParameter<int>( "Seed", 156 );
  TRandom3* rndm = new TRandom3( seed );

  EventType eventType( NamedParameter<std::string>( "EventType" , "", "EventType to generate, in the format: \033[3m parent daughter1 daughter2 ... \033[0m" ).getVector(),
                       NamedParameter<bool>( "GenerateTimeDependent", false , "Flag to include possible time dependence of the amplitude") );

  bool verbose = NamedParameter<bool>("CoherentSum::Debug", 0 ) || 
                 NamedParameter<bool>("PolarisedSum::Debug", 0 );
  INFO("Using verbose mode: " << verbose );
  AmpGen::MinuitParameterSet MPS;
  MPS.loadFromStream();

  if ( NamedParameter<bool>( "conj", false ) == true ) 
  {
    eventType = eventType.conj();
    INFO( eventType );
    AddCPConjugate(MPS);
  }
  if( NamedParameter<bool>( "AddConj", false) == true && NamedParameter<bool>( "conj", false ) == false )
  {
    AddCPConjugate(MPS); 
  }
  INFO( "EventType = " << eventType );
  
  std::string infile = NamedParameter<std::string>("InputFile","");
  EventList accepted = infile == "" ? EventList( eventType ) : EventList( infile, eventType );
  
  std::string input_units = NamedParameter<std::string>("Units","GeV");
  if( input_units == "MeV" && infile != "") accepted.transform([](auto& event){ for( unsigned i = 0;i< event.size();++i) event[i]/=1000; } );
  if( infile == "" ){
    // for( unsigned i = 0 ; i != 16; ++i ){
    //   Event evt = PhaseSpace( eventType, rndm ).makeEvent();
    //   evt.setIndex(i);
    //   accepted.push_back(evt);
    // }

    // set event 4 mom here in evt_data
    // RS EventType D0 K+ pi- pi- pi+
    // WS EventType D0 K- pi+ pi+ pi-

    real_t my_evt_data[16] = {
      -451.77782853/1000.0, 224.31692554/1000.0, 36.67676649/1000.0, 706.7414343/1000.0,
      376.61949967/1000.0, 232.89066672/1000.0, -248.6031968/1000.0, 526.65329656/1000.0,
      -30.84046039/1000.0, -61.19955508/1000.0, 97.49327222/1000.0, 183.52463817/1000.0,
      105.99878925/1000.0, -396.00803718/1000.0, 114.43315809/1000.0, 447.9206309/1000.0
    };

    // real_t my_evt_data[16] = { 
    //   546.42464231,  364.22822821,   72.98156009,  824.79414425,
    //   -448.52703318, -410.75077151,   77.66495823,  628.81187701,
    //   66.46711795,   17.64945192, -108.36458439,  189.61038795,
    //   -164.36472708,   28.87309138, -42.28193393, 221.6235908
    // };

    Event my_evt(my_evt_data, 16);
    accepted.push_back(my_evt);
    INFO("Got event");
  }
  std::vector<double> event = NamedParameter<double>("Event",0).getVector();
  if( event.size() != 1 ) accepted[0].set( event.data() );
  INFO("Set event");

  std::string type = NamedParameter<std::string>("Type","CoherentSum");
  INFO("Type = " << type);
  if( type == "PolarisedSum")
  {
    PolarisedSum sig( eventType, MPS );  
    sig.setEvents( accepted );
    sig.prepare();
    debug( sig, accepted);
    sig.setMC(accepted);
    INFO("norm = " << sig.norm() );   
  }
  else if( type == "CoherentSum" )
  {
    INFO("Building sig ...");
    CoherentSum sig(eventType, MPS);
    INFO("Built sig");
    debug(sig, accepted);
  }
  else if( type == "IncoherentSum" )
  {
    IncoherentSum sig(eventType, MPS,"Inco"); 
    sig.setMC(accepted);
    sig.prepare();
    debug(sig, accepted);  
    INFO("norm = " << sig.norm() );   
  }
  else {
    ERROR( "Type: " << type << " is not recognised");
  }
}
